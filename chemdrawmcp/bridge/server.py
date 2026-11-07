"""Local HTTP bridge between ChemDraw's Agent Chat add-in and the user's Claude Code / Codex CLI.

Standard library only.  Binds 127.0.0.1 only; every request (except CORS preflight) must carry the bridge token
(``X-Bridge-Token`` header or ``?token=`` query, the latter for EventSource which cannot set headers).  The
token lives in ``<CHEMDRAWMCP_HOME>/bridge.json`` together with the port and pid, where the add-in reads it.

Endpoints (JSON unless noted)::

    GET  /health                 {ok, version, backends:{claude,codex}, default_backend, settings, runs_active}
    GET  /settings               settings snapshot
    POST /settings               {backend?, model?, effort?, default_backend?} -> settings snapshot
    POST /run                    {backend?, prompt, session_id?, new?, context?, allow_shell?, model?, effort?}
                                 -> {run_id, backend, session_id}   (409 while that backend is busy)
    GET  /events/<run_id>        text/event-stream: replay + live events, ``event: <kind>`` / ``data: <json>``,
                                 ``id: <n>`` (Last-Event-ID / ?after= resume), final ``event: end``
    GET  /runs/<run_id>          {run_id, backend, status, events_count, final_text, session_id, ...}
    POST /cancel/<run_id>        {ok, status}
    GET  /history?backend=       [ {prompt, final_text, session_id, backend, started, ...}, ... ]
    POST /clear                  {backend?} -> {ok, sessions}
    GET  /sessions               {claude: id|null, codex: id|null}
"""
from __future__ import annotations

import argparse
import hmac
import json
import os
import secrets
import signal
import sys
import threading
import time
import urllib.request
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, urlsplit

from . import agentcli, settings as cfg

try:
    from .. import __version__
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

DEFAULT_PORT = 47325
BIND_HOST = "127.0.0.1"
MAX_RUNS = 20
KEEPALIVE_SECONDS = 15.0
MAX_BODY = 1 << 20
EVENT_KINDS = ("session", "delta", "text", "tool", "tool_result", "stats", "done", "error")
_LOCAL_HOSTS = {"127.0.0.1", "localhost", "::1", "[::1]"}


def _env_port() -> int:
    try:
        return int(os.environ.get("CHEMDRAWMCP_BRIDGE_PORT") or DEFAULT_PORT)
    except ValueError:
        return DEFAULT_PORT


class BridgeError(Exception):
    def __init__(self, status: int, message: str, **extra):
        super().__init__(message)
        self.status = status
        self.extra = extra


def origin_allowed(origin: Optional[str]) -> bool:
    """Only *remote* web origins are refused.  ChemDraw's add-in web view runs from a file:// page and sends no
    Origin, ``null`` or ``file://``; those (and any custom scheme) are fine, as is a local http(s) page."""
    if not origin:
        return True
    o = origin.strip()
    if o.lower() in ("null", "file://"):
        return True
    low = o.lower()
    if low.startswith("http://") or low.startswith("https://"):
        try:
            host = urlsplit(o).hostname or ""
        except ValueError:
            return False
        return host.lower() in _LOCAL_HOSTS
    return True


# ---------------------------------------------------------------------------
# runs
# ---------------------------------------------------------------------------

class Run:
    """One turn: a CLI process plus the ordered list of events it produced (appended by the agent thread)."""

    def __init__(self, run_id: str, backend: str, prompt: str, session_id: Optional[str]):
        self.run_id = run_id
        self.backend = backend
        self.prompt = prompt
        self.session_id = session_id
        self.started = time.time()
        self.finished: Optional[float] = None
        self.status = "running"            # running | done | error | cancelled
        self.events: List[Tuple[str, Any]] = []
        self.cond = threading.Condition()
        self.final_text = ""
        self.error: Optional[str] = None
        self.cancelled = False
        self.agent: Optional[agentcli.AgentRun] = None

    def to_dict(self) -> dict:
        with self.cond:
            return {"run_id": self.run_id, "backend": self.backend, "status": self.status,
                    "events_count": len(self.events), "final_text": self.final_text,
                    "session_id": self.session_id, "prompt": self.prompt, "started": self.started,
                    "finished": self.finished, "error": self.error}

    def history_entry(self) -> dict:
        return {"run_id": self.run_id, "backend": self.backend, "prompt": self.prompt,
                "final_text": self.final_text, "session_id": self.session_id, "started": self.started,
                "finished": self.finished, "status": self.status, "error": self.error}


class Bridge:
    """State shared by all handler threads: token, runs, session store."""

    def __init__(self, token: str, store: Optional[cfg.SessionStore] = None, max_runs: int = MAX_RUNS,
                 cwd: Optional[str] = None):
        self.token = token
        self.store = store or cfg.SessionStore()
        self.max_runs = max_runs
        self.cwd = cwd
        self.lock = threading.RLock()
        self.runs: "OrderedDict[str, Run]" = OrderedDict()
        self.active: Dict[str, Run] = {}
        self.stopping = threading.Event()
        self._backends_cache: Tuple[float, dict] = (0.0, {})

    # -- backends / settings ------------------------------------------------
    def backends(self, refresh: bool = False) -> dict:
        now = time.time()
        ts, cached = self._backends_cache
        if cached and not refresh and now - ts < 30:
            return dict(cached)
        found = {b: agentcli.find_backend(b) for b in cfg.BACKENDS}
        self._backends_cache = (now, found)
        return dict(found)

    def default_backend(self) -> str:
        stored = cfg.stored_default_backend()
        avail = [b for b, p in self.backends().items() if p]
        if stored in avail:
            return stored
        return avail[0] if avail else (stored or "claude")

    def health(self) -> dict:
        with self.lock:
            active = len(self.active)
        return {"ok": True, "version": __version__, "app": agentcli.APP, "backends": self.backends(),
                "default_backend": self.default_backend(), "settings": cfg.settings_snapshot(self.default_backend()),
                "runs_active": active, "pid": os.getpid(), "cwd": self.cwd or cfg.workspace_dir(),
                "sessions": self.store.sessions()}

    def update_settings(self, body: dict) -> dict:
        backend = body.get("backend")
        if backend is not None and backend not in cfg.BACKENDS:
            raise BridgeError(400, "backend must be one of %s" % (cfg.BACKENDS,))
        model = body.get("model")
        effort = body.get("effort")
        if (model is not None or effort is not None) and backend is None:
            raise BridgeError(400, "model/effort need a backend")
        if effort is not None:
            effort = str(effort)
            if effort not in cfg.EFFORT_CHOICES[backend]:
                raise BridgeError(400, "effort for %s must be one of %s" % (backend, cfg.EFFORT_CHOICES[backend]))
        if model is not None:
            model = str(model)
        if backend is not None and (model is not None or effort is not None):
            cfg.set_backend_settings(backend, model, effort)
        default = body.get("default_backend")
        if default is not None:
            if default not in cfg.BACKENDS:
                raise BridgeError(400, "default_backend must be one of %s" % (cfg.BACKENDS,))
            cfg.set_default_backend(default)
        return cfg.settings_snapshot(self.default_backend())

    # -- runs ---------------------------------------------------------------
    def get_run(self, run_id: str) -> Run:
        with self.lock:
            run = self.runs.get(run_id)
        if run is None:
            raise BridgeError(404, "unknown run %s" % run_id)
        return run

    def start_run(self, body: dict) -> Run:
        prompt = body.get("prompt")
        if not isinstance(prompt, str) or not prompt.strip():
            raise BridgeError(400, "prompt (non-empty string) is required")
        backend = body.get("backend") or self.default_backend()
        if backend not in cfg.BACKENDS:
            raise BridgeError(400, "backend must be one of %s" % (cfg.BACKENDS,))
        exe = agentcli.find_backend(backend)
        if not exe:
            raise BridgeError(503, "%s CLI not found (install it or set %s)" % (
                cfg.BACKEND_LABELS[backend], "CHEMDRAWMCP_CLAUDE_EXE" if backend == "claude" else "CHEMDRAWMCP_CODEX_EXE"))
        sid = body.get("session_id")
        if body.get("new"):
            session_id = None
        elif isinstance(sid, str) and sid.strip():
            session_id = sid.strip()
        else:
            session_id = self.store.get_session(backend)
        model = body.get("model")
        effort = body.get("effort")
        if effort is not None and str(effort) not in cfg.EFFORT_CHOICES[backend]:
            raise BridgeError(400, "effort for %s must be one of %s" % (backend, cfg.EFFORT_CHOICES[backend]))
        context = body.get("context")
        if context is not None and not isinstance(context, str):
            context = json.dumps(context)
        with self.lock:
            busy = self.active.get(backend)
            if busy is not None:
                raise BridgeError(409, "%s is busy with run %s" % (backend, busy.run_id), run_id=busy.run_id)
            run = Run(secrets.token_hex(8), backend, prompt, session_id)
            try:
                run.agent = agentcli.AgentRun(
                    backend, prompt, callback=lambda kind, payload: self._on_event(run, kind, payload),
                    session_id=session_id, cwd=self.cwd, allow_shell=bool(body.get("allow_shell")), exe=exe,
                    model=None if model is None else str(model), effort=None if effort is None else str(effort),
                    extra_context=context)
            except Exception as e:
                raise BridgeError(500, "could not prepare %s run: %s" % (backend, e))
            self.active[backend] = run
            self.runs[run.run_id] = run
            self._prune()
            run.agent.start()
        return run

    def _prune(self):
        with self.lock:
            while len(self.runs) > self.max_runs:
                victim = None
                for rid, r in self.runs.items():
                    if r.status != "running":
                        victim = rid
                        break
                if victim is None:
                    break
                del self.runs[victim]

    def _on_event(self, run: Run, kind: str, payload: Any):
        """AgentRun callback (agent thread): record the event and, on 'done', finish the run."""
        if kind == "session" and payload:
            run.session_id = str(payload)
            try:
                self.store.set_session(run.backend, run.session_id)
            except Exception:
                pass
        finished = False
        with run.cond:
            if kind == "text" and payload:
                run.final_text = str(payload)
            elif kind == "error":
                run.error = str(payload)
            run.events.append((kind, payload))
            if kind == "done":
                if payload:
                    run.final_text = str(payload)
                run.finished = time.time()
                # The turn failed when the event right before 'done' is an error (CLI crash, is_error result,
                # cancellation); an error earlier in the turn (a failed tool call) that the agent recovered
                # from still counts as done.
                ended_on_error = len(run.events) >= 2 and run.events[-2][0] == "error"
                if run.cancelled:
                    run.status = "cancelled"
                elif ended_on_error:
                    run.status = "error"
                else:
                    run.status = "done"
                finished = True
            run.cond.notify_all()
        if finished:
            with self.lock:
                if self.active.get(run.backend) is run:
                    del self.active[run.backend]
            try:
                self.store.add_history(run.history_entry())
            except Exception:
                pass

    def cancel(self, run_id: str) -> dict:
        run = self.get_run(run_id)
        with run.cond:
            running = run.status == "running"
            if running:
                run.cancelled = True
        if running and run.agent is not None:
            run.agent.cancel()
        return {"ok": True, "run_id": run_id, "status": run.status, "was_running": running}

    def shutdown(self):
        self.stopping.set()
        with self.lock:
            runs = list(self.runs.values())
        for r in runs:
            if r.status == "running" and r.agent is not None:
                r.cancelled = True
                r.agent.cancel()
            with r.cond:
                r.cond.notify_all()


# ---------------------------------------------------------------------------
# HTTP handler
# ---------------------------------------------------------------------------

class BridgeHandler(BaseHTTPRequestHandler):
    server_version = "chemdraw-agent-bridge/%s" % __version__
    sys_version = ""
    quiet = False

    @property
    def bridge(self) -> Bridge:
        return self.server.bridge  # type: ignore[attr-defined]

    # -- plumbing -----------------------------------------------------------
    def log_message(self, fmt, *args):
        if self.quiet:
            return
        try:
            sys.stderr.write("[bridge %s] %s\n" % (time.strftime("%H:%M:%S"), fmt % args))
        except Exception:
            pass

    def _cors(self, preflight: bool = False):
        origin = self.headers.get("Origin")
        low = (origin or "").strip().lower()
        if low.startswith("http://") or low.startswith("https://"):
            self.send_header("Access-Control-Allow-Origin", origin.strip())
            self.send_header("Vary", "Origin")
        else:
            self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "X-Bridge-Token, Content-Type, Last-Event-ID")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Expose-Headers", "Content-Type")
        self.send_header("Cache-Control", "no-store")
        if preflight:
            self.send_header("Access-Control-Allow-Private-Network", "true")
            self.send_header("Access-Control-Max-Age", "600")

    def _json(self, status: int, obj: Any):
        data = json.dumps(obj).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self._cors()
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(data)

    def _read_body(self) -> dict:
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            raise BridgeError(400, "bad Content-Length")
        if n > MAX_BODY:
            raise BridgeError(413, "body too large")
        raw = self.rfile.read(n) if n else b""
        if not raw.strip():
            return {}
        try:
            body = json.loads(raw.decode("utf-8"))
        except Exception:
            raise BridgeError(400, "body must be JSON")
        if not isinstance(body, dict):
            raise BridgeError(400, "body must be a JSON object")
        return body

    def _authorized(self, query: dict) -> bool:
        token = self.headers.get("X-Bridge-Token")
        if not token:
            token = (query.get("token") or [None])[0]
        return bool(token) and hmac.compare_digest(str(token), self.bridge.token)

    def do_OPTIONS(self):
        self._dispatch("OPTIONS")

    def do_GET(self):
        self._dispatch("GET")

    def do_HEAD(self):
        self._dispatch("GET")

    def do_POST(self):
        self._dispatch("POST")

    def _dispatch(self, method: str):
        try:
            parts = urlsplit(self.path)
            path = parts.path.rstrip("/") or "/"
            query = parse_qs(parts.query)
            if not origin_allowed(self.headers.get("Origin")):
                return self._json(403, {"error": "origin not allowed: %s" % self.headers.get("Origin")})
            if method == "OPTIONS":
                self.send_response(204)
                self._cors(preflight=True)
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            if not self._authorized(query):
                return self._json(401, {"error": "missing or invalid bridge token (X-Bridge-Token header or ?token=)"})
            return self._route(method, path, query)
        except BridgeError as e:
            self._json(e.status, dict({"error": str(e)}, **e.extra))
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError):
            pass
        except Exception as e:  # pragma: no cover - defensive
            try:
                self._json(500, {"error": "%s: %s" % (type(e).__name__, e)})
            except Exception:
                pass

    # -- routes -------------------------------------------------------------
    def _route(self, method: str, path: str, query: dict):
        b = self.bridge
        if method == "GET" and path == "/":
            return self._json(200, {"ok": True, "service": "chemdraw-agent-bridge", "version": __version__,
                                    "endpoints": ["GET /health", "GET|POST /settings", "POST /run", "GET /events/<run_id>",
                                                  "GET /runs/<run_id>", "POST /cancel/<run_id>", "GET /history",
                                                  "POST /clear", "GET /sessions"]})
        if method == "GET" and path == "/health":
            return self._json(200, b.health())
        if path == "/settings":
            if method == "GET":
                return self._json(200, cfg.settings_snapshot(b.default_backend()))
            if method == "POST":
                return self._json(200, b.update_settings(self._read_body()))
        if method == "POST" and path == "/run":
            run = b.start_run(self._read_body())
            return self._json(200, {"run_id": run.run_id, "backend": run.backend, "session_id": run.session_id,
                                    "started": run.started})
        if method == "GET" and path.startswith("/events/"):
            run = b.get_run(path[len("/events/"):])
            after = self.headers.get("Last-Event-ID") or (query.get("after") or ["0"])[0]
            try:
                after_n = max(0, int(after))
            except ValueError:
                after_n = 0
            return self._stream(run, after_n)
        if method == "GET" and path.startswith("/runs/"):
            return self._json(200, b.get_run(path[len("/runs/"):]).to_dict())
        if method == "GET" and path == "/runs":
            with b.lock:
                runs = [r.to_dict() for r in b.runs.values()]
            return self._json(200, runs)
        if method == "POST" and path.startswith("/cancel/"):
            return self._json(200, b.cancel(path[len("/cancel/"):]))
        if method == "GET" and path == "/history":
            backend = (query.get("backend") or [None])[0]
            if backend and backend not in cfg.BACKENDS:
                raise BridgeError(400, "backend must be one of %s" % (cfg.BACKENDS,))
            limit = (query.get("limit") or [None])[0]
            try:
                limit_n = int(limit) if limit else None
            except ValueError:
                limit_n = None
            return self._json(200, b.store.history(backend, limit_n))
        if method == "POST" and path == "/clear":
            body = self._read_body()
            backend = body.get("backend")
            if backend and backend not in cfg.BACKENDS:
                raise BridgeError(400, "backend must be one of %s" % (cfg.BACKENDS,))
            sessions = b.store.clear(backend or None, history=bool(body.get("history")))
            return self._json(200, {"ok": True, "sessions": sessions})
        if method == "GET" and path == "/sessions":
            return self._json(200, b.store.sessions())
        raise BridgeError(404, "no such endpoint: %s %s" % (method, path))

    # -- SSE ----------------------------------------------------------------
    def _stream(self, run: Run, after: int):
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream; charset=utf-8")
        self.send_header("X-Accel-Buffering", "no")
        self.send_header("Connection", "close")
        self._cors()
        self.end_headers()
        if self.command == "HEAD":
            return
        w = self.wfile
        w.write(b"retry: 2000\n\n")
        w.flush()
        idx = after
        stopping = self.bridge.stopping
        while True:
            with run.cond:
                timed_out = False
                if idx >= len(run.events) and run.status == "running" and not stopping.is_set():
                    timed_out = not run.cond.wait(KEEPALIVE_SECONDS)
                batch = run.events[idx:]
                finished = run.status != "running"
            for kind, payload in batch:
                idx += 1
                w.write(("id: %d\nevent: %s\ndata: %s\n\n" % (idx, kind, json.dumps(payload))).encode("utf-8"))
            if batch:
                w.flush()
                if finished:
                    with run.cond:
                        if idx < len(run.events):
                            continue
            if finished and idx >= len(run.events):
                w.write(("id: %d\nevent: end\ndata: %s\n\n" % (idx, json.dumps({
                    "run_id": run.run_id, "status": run.status, "final_text": run.final_text,
                    "session_id": run.session_id, "error": run.error}))).encode("utf-8"))
                w.flush()
                return
            if stopping.is_set():
                w.write(b"event: end\ndata: {\"status\": \"bridge_stopping\"}\n\n")
                w.flush()
                return
            if timed_out and not batch:
                w.write(b": keep-alive\n\n")
                w.flush()


# ---------------------------------------------------------------------------
# server wrapper
# ---------------------------------------------------------------------------

def token_path() -> str:
    """The persistent bridge token file (``<CHEMDRAWMCP_HOME>/bridge_token``, mode 0600)."""
    return os.path.join(cfg.home_dir(), "bridge_token")


def ensure_token(new_token: bool = False) -> str:
    """Return the persistent bridge token, creating it when missing (used by the installer to build the add-in
    config and by the bridge at start, so a restarted bridge keeps the token the installed add-in knows)."""
    path = token_path()
    if not new_token:
        try:
            with open(path, encoding="utf-8") as fh:
                tok = fh.read().strip()
            if len(tok) >= 16:
                return tok
        except OSError:
            pass
        info = bridge_info()
        if isinstance(info, dict) and isinstance(info.get("token"), str) and len(info["token"]) >= 16:
            tok = info["token"]
        else:
            tok = secrets.token_urlsafe(24)
    else:
        tok = secrets.token_urlsafe(24)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(tok)
    return tok


def _load_token(new_token: bool) -> str:
    return ensure_token(new_token)


class BridgeServer:
    """Owns the ThreadingHTTPServer, its thread and ``bridge.json``."""

    def __init__(self, port: Optional[int] = None, host: str = BIND_HOST, new_token: bool = False,
                 quiet: bool = False, cwd: Optional[str] = None, write_info: bool = True, token: Optional[str] = None):
        self.port = _env_port() if port is None else int(port)
        self.host = host
        self.quiet = quiet
        self.write_info = write_info
        self.token = token or _load_token(new_token)
        self.bridge = Bridge(self.token, cwd=cwd)
        self.httpd: Optional[ThreadingHTTPServer] = None
        self.thread: Optional[threading.Thread] = None
        self.started: Optional[float] = None
        self.info_path = cfg.bridge_info_path()

    @property
    def url(self) -> str:
        return "http://%s:%d" % (self.host, self.port)

    def start(self) -> "BridgeServer":
        handler = type("BridgeHandler", (BridgeHandler,), {"quiet": self.quiet})
        self.httpd = ThreadingHTTPServer((self.host, self.port), handler)
        self.httpd.daemon_threads = True
        self.httpd.bridge = self.bridge  # type: ignore[attr-defined]
        self.port = self.httpd.server_address[1]
        self.started = time.time()
        self.thread = threading.Thread(target=self.httpd.serve_forever, kwargs={"poll_interval": 0.2},
                                       name="chemdrawmcp-bridge", daemon=True)
        self.thread.start()
        if self.write_info:
            cfg.write_json(self.info_path, self.info(), private=True)
        return self

    def info(self) -> dict:
        return {"port": self.port, "token": self.token, "pid": os.getpid(), "started": self.started,
                "url": self.url, "version": __version__}

    def stop(self, timeout: float = 5.0):
        self.bridge.shutdown()
        if self.httpd is not None:
            try:
                self.httpd.shutdown()
                self.httpd.server_close()
            except Exception:
                pass
        if self.thread is not None:
            self.thread.join(timeout)
        if self.write_info:
            info = cfg.read_json(self.info_path)
            if isinstance(info, dict) and info.get("pid") == os.getpid():
                try:
                    os.remove(self.info_path)
                except OSError:
                    pass


# ---------------------------------------------------------------------------
# helpers for the installer / MCP server
# ---------------------------------------------------------------------------

def bridge_info() -> Optional[dict]:
    """Contents of ``bridge.json`` ({port, token, pid, started, url}) or ``None`` when no bridge wrote one."""
    info = cfg.read_json(cfg.bridge_info_path())
    if isinstance(info, dict) and isinstance(info.get("port"), int) and isinstance(info.get("token"), str):
        return info
    return None


def bridge_health(info: Optional[dict] = None, timeout: float = 1.0) -> Optional[dict]:
    """GET /health of the bridge described by ``info`` (default: bridge.json); the JSON or ``None``."""
    info = info or bridge_info()
    if not info:
        return None
    req = urllib.request.Request("http://%s:%d/health" % (BIND_HOST, info["port"]),
                                 headers={"X-Bridge-Token": info["token"]})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            d = json.loads(r.read().decode("utf-8"))
        return d if isinstance(d, dict) and d.get("ok") else None
    except Exception:
        return None


def is_bridge_running(info: Optional[dict] = None, timeout: float = 1.0) -> bool:
    return bridge_health(info, timeout) is not None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="chemdraw-agent-bridge",
                                 description="Local HTTP bridge between ChemDraw's Agent Chat panel and the Claude Code / Codex CLI.")
    ap.add_argument("--port", type=int, default=None, help="listen port (default $CHEMDRAWMCP_BRIDGE_PORT or %d)" % DEFAULT_PORT)
    ap.add_argument("--new-token", action="store_true", help="generate a fresh token instead of reusing bridge.json's")
    ap.add_argument("--quiet", action="store_true", help="no request log / banner")
    ap.add_argument("--check", action="store_true", help="print the CLIs found and exit")
    ap.add_argument("--cwd", default=None, help="working directory for the CLIs (default $CHEMDRAWMCP_BRIDGE_CWD or <home>/workspace)")
    args = ap.parse_args(argv)

    if args.check:
        found = {b: agentcli.find_backend(b) for b in cfg.BACKENDS}
        for b in cfg.BACKENDS:
            print("%-6s %s" % (b + ":", found[b] or "not found"))
        print("default: %s" % agentcli.default_backend())
        print("home:    %s" % cfg.home_dir())
        info = bridge_info()
        if info:
            print("bridge:  %s (%s)" % (info.get("url") or info.get("port"), "running" if is_bridge_running(info) else "stale bridge.json"))
        return 0 if any(found.values()) else 1

    srv = BridgeServer(port=args.port, new_token=args.new_token, quiet=args.quiet, cwd=args.cwd)
    try:
        srv.start()
    except OSError as e:
        info = bridge_info()
        if info and info.get("port") == srv.port and is_bridge_running(info):
            print("chemdraw-agent-bridge already running at %s (pid %s)" % (info.get("url"), info.get("pid")))
            return 0
        print("chemdraw-agent-bridge: cannot listen on %s: %s" % (srv.url, e), file=sys.stderr)
        return 1
    if not args.quiet:
        print("chemdraw-agent-bridge %s listening on %s" % (__version__, srv.url))
        print("token file: %s" % srv.info_path)
        found = srv.bridge.backends()
        print("backends: %s" % ", ".join("%s=%s" % (b, found[b] or "-") for b in cfg.BACKENDS))
        sys.stdout.flush()

    stop = threading.Event()

    def _sig(signum, frame):
        stop.set()
    for name in ("SIGTERM", "SIGINT", "SIGBREAK"):
        s = getattr(signal, name, None)
        if s is not None:
            try:
                signal.signal(s, _sig)
            except Exception:
                pass
    try:
        while not stop.wait(0.5):
            pass
    except KeyboardInterrupt:
        pass
    finally:
        srv.stop()
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
