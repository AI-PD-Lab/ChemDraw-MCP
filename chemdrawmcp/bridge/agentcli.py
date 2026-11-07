"""Drive Claude Code or Codex CLI non-interactively for the ChemDraw Agent Chat panel.

Standard library only.  The CLI is the user's own installed Claude Code / Codex; it loads the ``chemdraw`` MCP
server from the user's MCP config (registered by the chemdrawmcp installer) and acts on the running ChemDraw
through its tools.  Adapted from molmcp's ``inapp/shared/agentcli.py`` (no process pinning: the chemdraw MCP
server always talks to the one running ChemDraw).

Events delivered to the callback as (kind, payload):
    ("session", id)          conversation id usable for follow-up turns
    ("delta", str)           streamed partial assistant text (Claude only)
    ("text", str)            complete assistant message
    ("tool", str)            tool call: "<server>/<tool> <json args>"
    ("tool_result", str)     tool result text (truncated)
    ("stats", dict)          cost/usage at the end of a turn
    ("done", str)            final answer (also delivered as text)
    ("error", str)
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import threading
from typing import Callable, List, Optional

from .settings import (BACKENDS, BACKEND_LABELS, EFFORT_CHOICES, MODEL_CHOICES, backend_settings,  # noqa: F401
                       home_dir, load_settings, save_settings, set_backend_settings, set_default_backend,
                       stored_default_backend, workspace_dir)

APP = "chemdraw"
APP_TITLE = "ChemDraw"

SYSTEM_CONTEXT = (
    "You are embedded in ChemDraw's Agent Chat panel through the chemdrawmcp bridge. The user is looking at the "
    "ChemDraw window right now and reads your answers in a small side panel. Use the 'chemdraw' MCP server tools "
    "to inspect and change what they see: status, documents, read_document, selection, draw_structures, "
    "draw_reaction, draw_scheme, draw_mechanism, insert, render, screenshot, run_command, properties, "
    "name_to_structure, convert, ... Prefer acting on the front document (document=None) and on the current "
    "selection when the user refers to 'this' structure; open a new document only when asked or when the front "
    "one would be clobbered. Never call quit_chemdraw or activate_chemdraw unless the user explicitly asks. "
    "Describe what you changed and keep answers concise (a few sentences; no long preambles)."
)


def _win() -> bool:
    return sys.platform.startswith("win")


def _candidates_from_env(*names: str) -> List[str]:
    out = []
    for n in names:
        v = os.environ.get(n)
        if v:
            out.append(v)
    return out


_validated: dict = {}


def _works(exe: str) -> bool:
    """A candidate counts only if `exe --version` succeeds (stale npm installs, wrong-arch builds...)."""
    if exe in _validated:
        return _validated[exe]
    ok = False
    try:
        kw = dict(stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, timeout=20,
                  env=clean_env())
        if _win():
            kw["creationflags"] = 0x08000000
        r = subprocess.run([exe, "--version"], **kw)
        ok = r.returncode == 0
    except Exception:
        ok = False
    _validated[exe] = ok
    return ok


def _first(paths: List[str]) -> Optional[str]:
    for p in paths:
        for hit in sorted(glob.glob(os.path.expanduser(os.path.expandvars(p))), reverse=True):
            if os.path.isfile(hit) and os.access(hit, os.X_OK) and _works(hit):
                return hit
    return None


def find_claude() -> Optional[str]:
    for v in _candidates_from_env("CHEMDRAWMCP_CLAUDE_EXE", "CLAUDE_CLI"):
        if os.path.exists(v):
            return v
    for name in ("claude", "claude.exe", "claude.cmd"):
        p = shutil.which(name)
        if p and _works(p):
            return p
    home = os.path.expanduser("~")
    pats = [
        os.path.join(home, ".claude", "local", "claude"),
        os.path.join(home, ".claude", "local", "node_modules", ".bin", "claude"),
        os.path.join(home, ".vscode*", "extensions", "anthropic.claude-code-*", "resources", "native-binary", "claude*"),
        os.path.join(home, ".cursor", "extensions", "anthropic.claude-code-*", "resources", "native-binary", "claude*"),
        os.path.join(home, ".npm-global", "bin", "claude"),
        "/usr/local/bin/claude", "/opt/homebrew/bin/claude",
    ]
    if _win():
        pats += [os.path.join(os.environ.get("APPDATA", ""), "npm", "claude.cmd"),
                 os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "claude*", "claude.exe")]
    return _first(pats)


def find_codex() -> Optional[str]:
    for v in _candidates_from_env("CHEMDRAWMCP_CODEX_EXE", "CODEX_CLI"):
        if os.path.exists(v):
            return v
    for name in ("codex", "codex.exe", "codex.cmd"):
        p = shutil.which(name)
        if p and _works(p):
            return p
    home = os.path.expanduser("~")
    pats = []
    if sys.platform == "darwin":
        pats += ["/Applications/ChatGPT.app/Contents/Resources/codex-cli/bin/codex",
                 os.path.join(home, "Applications", "ChatGPT.app", "Contents", "Resources", "codex-cli", "bin", "codex")]
    pats += [os.path.join(home, ".npm-global", "bin", "codex"), "/usr/local/bin/codex", "/opt/homebrew/bin/codex"]
    if _win():
        pats += [os.path.join(os.environ.get("APPDATA", ""), "npm", "codex.cmd"),
                 os.path.join(os.environ.get("LOCALAPPDATA", ""), "Programs", "ChatGPT", "resources", "codex-cli", "bin", "codex*")]
    return _first(pats)


def find_backend(backend: str) -> Optional[str]:
    return find_claude() if backend == "claude" else find_codex()


def available_backends() -> List[str]:
    """Backends whose CLI is installed and runs (so the UI can hide the missing one)."""
    return [b for b in BACKENDS if find_backend(b)]


def default_backend() -> str:
    d = stored_default_backend()
    avail = available_backends()
    if d in avail:
        return d
    return avail[0] if avail else (d or "claude")


# Variables a host application's launcher may set for its *own* bundled Python; inherited by the CLI (and by
# the MCP server it spawns from a different Python) they break interpreter start-up.
_STRIP_ENV_EXACT = {"PYTHONHOME", "PYTHONPATH", "PYTHONEXECUTABLE", "PYTHONNOUSERSITE", "PYTHONSTARTUP", "PYTHONSAFEPATH",
                    "__PYVENV_LAUNCHER__", "DYLD_LIBRARY_PATH", "DYLD_FRAMEWORK_PATH", "DYLD_FALLBACK_LIBRARY_PATH",
                    "LD_LIBRARY_PATH", "LD_PRELOAD", "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "QML2_IMPORT_PATH",
                    "QT_QPA_PLATFORM", "SSL_CERT_FILE", "REQUESTS_CA_BUNDLE", "FONTCONFIG_PATH", "FONTCONFIG_FILE",
                    "GSETTINGS_SCHEMA_DIR", "GIO_MODULE_DIR", "XDG_DATA_DIRS", "VIRTUAL_ENV"}
_STANDARD_PATH = ["/opt/homebrew/bin", "/usr/local/bin", "/usr/bin", "/bin", "/usr/sbin", "/sbin"]


def clean_env() -> dict:
    """A copy of the environment without private interpreter settings; PATH gets the usual user directories
    back (GUI apps / launch agents started from a dock, Finder or Start menu see a minimal PATH)."""
    env = {k: v for k, v in os.environ.items() if k not in _STRIP_ENV_EXACT}
    if not _win():
        home = os.path.expanduser("~")
        parts = [p for p in env.get("PATH", "").split(os.pathsep) if p]
        for extra in [os.path.join(home, ".local", "bin"), os.path.join(home, ".npm-global", "bin"),
                      os.path.join(home, ".claude", "local")] + _STANDARD_PATH:
            if extra not in parts:
                parts.append(extra)
        env["PATH"] = os.pathsep.join(parts)
    return env


def system_context(extra_context: Optional[str] = None) -> str:
    """The system prompt appended for every turn; ``extra_context`` is panel-supplied per-turn state (e.g. the
    SMILES of the current selection or the front document's name)."""
    ctx = SYSTEM_CONTEXT
    if extra_context and str(extra_context).strip():
        ctx += "\n\nCurrent ChemDraw context supplied by the panel:\n" + str(extra_context).strip()
    return ctx


def build_argv(backend: str, exe: str, prompt: str, session_id: Optional[str] = None, allow_shell: bool = False,
               model: str = "", effort: str = "", extra_context: Optional[str] = None) -> List[str]:
    context = system_context(extra_context)
    if backend == "claude":
        argv = [exe, "-p", prompt, "--output-format", "stream-json", "--verbose", "--include-partial-messages",
                "--allowedTools", "mcp__%s" % APP + (",Bash,Edit,Write" if allow_shell else ""),
                "--append-system-prompt", context]
        if model:
            argv += ["--model", model]
        if effort:
            argv += ["--effort", effort]
        if session_id:
            argv += ["--resume", session_id]
        return argv
    if backend == "codex":
        # Codex has no system-prompt flag: the context rides in the first turn's prompt; on resumed turns only the
        # per-turn panel context is prepended.
        if session_id:
            full_prompt = prompt
            if extra_context and str(extra_context).strip():
                full_prompt = "Current ChemDraw context supplied by the panel:\n%s\n\nUser request: %s" % (
                    str(extra_context).strip(), prompt)
        else:
            full_prompt = "%s\n\nUser request: %s" % (context, prompt)
        opts = ["--json", "--approve-for-me", "--skip-git-repo-check"]
        if model:
            opts += ["-m", model]
        if effort:
            opts += ["-c", 'model_reasoning_effort="%s"' % effort]
        if session_id:
            return [exe, "exec", "resume"] + opts + [session_id, full_prompt]
        return [exe, "exec"] + opts + [full_prompt]
    raise ValueError("unknown backend %r" % backend)


def _unwrap_result(obj) -> str:
    """MCP result envelopes ({"content":[{"type":"text","text":...}]}, {"result": "<json>"}) -> readable text."""
    for _ in range(3):
        if isinstance(obj, str):
            try:
                parsed = json.loads(obj)
            except Exception:
                return obj
            if isinstance(parsed, (dict, list)):
                obj = parsed
                continue
            return obj
        if isinstance(obj, dict):
            if isinstance(obj.get("content"), list):
                parts = []
                for c in obj["content"]:
                    if isinstance(c, dict) and c.get("type") == "text":
                        parts.append(str(c.get("text", "")))
                    elif isinstance(c, dict) and c.get("type") == "image":
                        parts.append("[image]")
                obj = "\n".join(parts)
                continue
            if set(obj) == {"result"} and isinstance(obj["result"], str):
                obj = obj["result"]
                continue
        break
    return obj if isinstance(obj, str) else json.dumps(obj)


def _tool_result_text(block) -> str:
    content = block.get("content")
    if isinstance(content, list):
        parts = []
        for c in content:
            if isinstance(c, dict) and c.get("type") == "text":
                parts.append(str(c.get("text", "")))
            elif isinstance(c, dict) and c.get("type") == "image":
                parts.append("[image]")
        text = "\n".join(parts)
    else:
        text = str(content or "")
    return _unwrap_result(text)


def _parse_claude(ev: dict, cb: Callable):
    t = ev.get("type")
    if t == "system" and ev.get("subtype") == "init" and ev.get("session_id"):
        cb("session", ev["session_id"])
    elif t == "stream_event":
        e = ev.get("event", {}) or {}
        if e.get("type") == "content_block_delta":
            d = e.get("delta", {}) or {}
            if d.get("type") == "text_delta" and d.get("text"):
                cb("delta", d["text"])
    elif t == "assistant":
        for c in ev.get("message", {}).get("content", []) or []:
            if c.get("type") == "text" and c.get("text"):
                cb("text", c["text"])
            elif c.get("type") == "tool_use":
                name = str(c.get("name", "")).replace("mcp__", "").replace("__", "/")
                args = json.dumps(c.get("input", {}))
                cb("tool", "%s %s" % (name, args[:300]))
    elif t == "user":
        for c in ev.get("message", {}).get("content", []) or []:
            if isinstance(c, dict) and c.get("type") == "tool_result":
                cb("tool_result", _tool_result_text(c)[:4000])
    elif t == "result":
        cb("stats", {k: ev.get(k) for k in ("duration_ms", "total_cost_usd", "num_turns", "usage") if k in ev})
        if ev.get("is_error") or str(ev.get("subtype", "")).startswith("error"):
            cb("error", str(ev.get("result") or ev.get("subtype")))
        cb("done", str(ev.get("result") or ""))


def _parse_codex(ev: dict, cb: Callable):
    t = ev.get("type", "")
    if t == "thread.started" and ev.get("thread_id"):
        cb("session", ev["thread_id"])
    elif t.startswith("item."):
        it = ev.get("item", {}) or {}
        kind = it.get("type")
        if kind == "agent_message" and t == "item.completed" and it.get("text"):
            cb("text", it["text"])
        elif kind == "mcp_tool_call":
            if t == "item.started":
                cb("tool", "%s/%s %s" % (it.get("server"), it.get("tool"), json.dumps(it.get("arguments", it.get("input", {})))[:300]))
            elif t == "item.completed":
                res = it.get("result")
                if res is None:
                    res = it.get("output")
                if res is not None:
                    cb("tool_result", _unwrap_result(res)[:4000])
                elif it.get("status") and it.get("status") != "completed":
                    cb("error", "%s/%s %s" % (it.get("server"), it.get("tool"), it.get("status")))
        elif kind == "command_execution" and t == "item.completed":
            cb("tool", "shell %s" % str(it.get("command", ""))[:200])
            if it.get("aggregated_output"):
                cb("tool_result", str(it["aggregated_output"])[:4000])
    elif t == "turn.completed":
        cb("stats", {"usage": ev.get("usage") or {}})
        cb("done", "")
    elif t in ("turn.failed", "error"):
        cb("error", json.dumps({k: v for k, v in ev.items() if k != "type"})[:500])


class AgentRun:
    """One prompt -> one CLI process, streamed on a background thread."""

    def __init__(self, backend: str, prompt: str, callback: Callable, session_id: Optional[str] = None,
                 cwd: Optional[str] = None, allow_shell: bool = False, exe: Optional[str] = None,
                 extra_env: Optional[dict] = None, model: Optional[str] = None, effort: Optional[str] = None,
                 extra_context: Optional[str] = None):
        if backend not in BACKENDS:
            raise ValueError("backend must be one of %s" % (BACKENDS,))
        self.backend = backend
        self.exe = exe or find_backend(backend)
        if not self.exe:
            raise RuntimeError("%s CLI not found. Install it or set %s." % (
                BACKEND_LABELS[backend], "CHEMDRAWMCP_CLAUDE_EXE" if backend == "claude" else "CHEMDRAWMCP_CODEX_EXE"))
        st = backend_settings(backend)
        self.model = st["model"] if model is None else model
        self.effort = st["effort"] if effort is None else effort
        self.prompt = prompt
        self.argv = build_argv(backend, self.exe, prompt, session_id, allow_shell, self.model, self.effort, extra_context)
        self.callback = callback
        self.cwd = cwd or workspace_dir()
        self.env = clean_env()
        self.env.setdefault("CHEMDRAWMCP_HOME", home_dir())
        if extra_env:
            self.env.update(extra_env)
        self.proc: Optional[subprocess.Popen] = None
        self.session_id: Optional[str] = session_id
        self.final_text = ""
        self._draft = ""
        self.thread: Optional[threading.Thread] = None
        self.returncode: Optional[int] = None
        self.cancelled = False

    def start(self):
        self.thread = threading.Thread(target=self._run, name="chemdrawmcp-agent-%s" % self.backend, daemon=True)
        self.thread.start()
        return self

    def cancel(self):
        self.cancelled = True
        if self.proc and self.proc.poll() is None:
            try:
                if _win():
                    # terminate() only ends the launcher (.cmd / node shim); take the whole process tree down
                    subprocess.run(["taskkill", "/F", "/T", "/PID", str(self.proc.pid)], capture_output=True, timeout=15)
                else:
                    self.proc.terminate()
            except Exception:
                try:
                    self.proc.kill()
                except Exception:
                    pass

    def wait(self, timeout: Optional[float] = None):
        if self.thread:
            self.thread.join(timeout)
        return self.final_text

    def _emit(self, kind, payload):
        if kind == "session":
            self.session_id = payload
        if kind == "done":
            self.final_text = payload or self.final_text
        if kind == "text":
            self.final_text = payload
        elif kind == "delta":
            self._draft += payload
        try:
            self.callback(kind, payload)
        except Exception:
            pass

    def _run(self):
        kw = dict(cwd=self.cwd, env=self.env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                  stderr=subprocess.PIPE, text=True, bufsize=1, encoding="utf-8", errors="replace")
        if _win():
            kw["creationflags"] = 0x08000000  # CREATE_NO_WINDOW
        try:
            self.proc = subprocess.Popen(self.argv, **kw)
        except Exception as e:
            self._emit("error", "could not start %s: %s" % (self.exe, e))
            self._emit("done", "")
            return
        if self.cancelled:  # cancelled between construction and process start
            self.cancel()
        parser = _parse_claude if self.backend == "claude" else _parse_codex
        err_lines = []

        def drain_err():
            try:
                for line in self.proc.stderr:
                    err_lines.append(line.rstrip())
            except Exception:
                pass
        threading.Thread(target=drain_err, daemon=True).start()
        done_seen = False
        try:
            for line in self.proc.stdout:
                line = line.strip()
                if not line:
                    continue
                try:
                    ev = json.loads(line)
                except Exception:
                    if line.startswith("Reading additional input"):
                        continue
                    self._emit("text", line)
                    continue
                if not isinstance(ev, dict):
                    continue
                parser(ev, self._emit)
                if ev.get("type") in ("result", "turn.completed"):
                    done_seen = True
        finally:
            self.returncode = self.proc.wait()
            if self.cancelled and not done_seen:
                self._emit("error", "cancelled")
            elif self.returncode not in (0, None) and not done_seen:
                tail = "\n".join(err_lines[-8:])
                self._emit("error", "%s exited with code %s%s" % (self.backend, self.returncode, (": " + tail) if tail else ""))
            if not done_seen:
                self._emit("done", self.final_text)
