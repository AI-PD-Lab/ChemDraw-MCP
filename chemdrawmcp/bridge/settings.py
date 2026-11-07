"""Persistent state of the ChemDraw agent bridge (standard library only).

Everything lives under ``$CHEMDRAWMCP_HOME`` (default ``~/.chemdrawmcp``):

* ``agent_settings.json``  - model/effort per backend and the default backend (``CHEMDRAWMCP_AGENT_SETTINGS``
  overrides the path);
* ``agent_sessions.json``  - the current conversation id per backend plus the turn history, so a bridge restart
  keeps conversations resumable;
* ``bridge.json``          - port/token/pid of the running bridge (written by :mod:`chemdrawmcp.bridge.server`).
"""
from __future__ import annotations

import json
import os
import threading
import time
from typing import Any, Dict, List, Optional

BACKENDS = ("claude", "codex")
BACKEND_LABELS = {"claude": "Claude Code", "codex": "Codex"}
# Model choices offered in the UI (free text is always allowed)
MODEL_CHOICES = {
    "claude": [("", "default"), ("fable", "fable"), ("opus", "opus"), ("sonnet", "sonnet"), ("haiku", "haiku")],
    "codex": [("", "default (config.toml)"), ("gpt-5-codex", "gpt-5-codex"), ("gpt-5", "gpt-5"), ("o3", "o3"),
              ("gpt-6-sol", "gpt-6-sol")],
}
EFFORT_CHOICES = {
    "claude": ["", "low", "medium", "high", "xhigh", "max"],
    "codex": ["", "minimal", "low", "medium", "high", "xhigh"],
}
HISTORY_LIMIT = 200


# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------

def home_dir() -> str:
    """``$CHEMDRAWMCP_HOME`` or ``~/.chemdrawmcp`` (created when possible)."""
    d = os.environ.get("CHEMDRAWMCP_HOME") or os.path.join(os.path.expanduser("~"), ".chemdrawmcp")
    d = os.path.abspath(os.path.expanduser(d))
    try:
        os.makedirs(d, exist_ok=True)
    except OSError:
        pass
    return d


def settings_path() -> str:
    return os.environ.get("CHEMDRAWMCP_AGENT_SETTINGS") or os.path.join(home_dir(), "agent_settings.json")


def sessions_path() -> str:
    return os.path.join(home_dir(), "agent_sessions.json")


def bridge_info_path() -> str:
    return os.path.join(home_dir(), "bridge.json")


def workspace_dir() -> str:
    """Working directory the CLIs run in (``CHEMDRAWMCP_BRIDGE_CWD`` or ``<home>/workspace``).  Claude Code keys
    its resumable sessions by working directory, so every turn must use the same one."""
    d = os.environ.get("CHEMDRAWMCP_BRIDGE_CWD") or os.path.join(home_dir(), "workspace")
    d = os.path.abspath(os.path.expanduser(d))
    try:
        os.makedirs(d, exist_ok=True)
    except OSError:
        return os.path.expanduser("~")
    return d


def read_json(path: str) -> Any:
    """Parsed JSON or ``None`` (missing/corrupt file)."""
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def write_json(path: str, obj: Any, private: bool = False) -> bool:
    """Atomic JSON write (temp file + rename); ``private`` restricts the file to the user (tokens)."""
    try:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        tmp = "%s.%d.tmp" % (path, os.getpid())
        flags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
        fd = os.open(tmp, flags, 0o600 if private else 0o644)
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=1)
        os.replace(tmp, path)
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# agent settings (model / effort / default backend)
# ---------------------------------------------------------------------------

def load_settings() -> dict:
    d = read_json(settings_path())
    return d if isinstance(d, dict) else {}


def save_settings(d: dict) -> None:
    write_json(settings_path(), d)


def backend_settings(backend: str) -> dict:
    """{"model": str, "effort": str} for a backend, persisted per user."""
    d = load_settings().get(backend, {})
    if not isinstance(d, dict):
        d = {}
    return {"model": str(d.get("model", "") or ""), "effort": str(d.get("effort", "") or "")}


def set_backend_settings(backend: str, model: Optional[str] = None, effort: Optional[str] = None) -> dict:
    d = load_settings()
    cur = d.get(backend, {})
    if not isinstance(cur, dict):
        cur = {}
    if model is not None:
        cur["model"] = model
    if effort is not None:
        cur["effort"] = effort
    d[backend] = cur
    d.setdefault("default_backend", backend)
    save_settings(d)
    return {"model": str(cur.get("model", "") or ""), "effort": str(cur.get("effort", "") or "")}


def stored_default_backend() -> Optional[str]:
    d = load_settings().get("default_backend")
    return d if d in BACKENDS else None


def set_default_backend(backend: str) -> None:
    d = load_settings()
    d["default_backend"] = backend
    save_settings(d)


def settings_snapshot(default: Optional[str] = None) -> dict:
    """Everything the panel needs to render its settings UI."""
    return {
        "claude": backend_settings("claude"),
        "codex": backend_settings("codex"),
        "default_backend": default or stored_default_backend() or "claude",
        "model_choices": {b: [c[0] for c in MODEL_CHOICES[b]] for b in BACKENDS},
        "effort_choices": EFFORT_CHOICES,
    }


# ---------------------------------------------------------------------------
# sessions + history
# ---------------------------------------------------------------------------

class SessionStore:
    """Conversation ids per backend and the turn history, persisted to ``agent_sessions.json``.

    Thread-safe; every mutation is written through immediately so a crash of the bridge loses nothing."""

    def __init__(self, path: Optional[str] = None, history_limit: int = HISTORY_LIMIT):
        self.path = path or sessions_path()
        self.history_limit = history_limit
        self._lock = threading.RLock()
        self._sessions: Dict[str, Optional[str]] = {b: None for b in BACKENDS}
        self._history: List[dict] = []
        self.load()

    def load(self) -> None:
        d = read_json(self.path)
        with self._lock:
            if isinstance(d, dict):
                s = d.get("sessions") or {}
                for b in BACKENDS:
                    v = s.get(b) if isinstance(s, dict) else None
                    self._sessions[b] = str(v) if v else None
                h = d.get("history") or []
                self._history = [e for e in h if isinstance(e, dict)][-self.history_limit:]

    def save(self) -> None:
        with self._lock:
            write_json(self.path, {"sessions": dict(self._sessions), "history": list(self._history),
                                   "saved": time.time()})

    def sessions(self) -> Dict[str, Optional[str]]:
        with self._lock:
            return dict(self._sessions)

    def get_session(self, backend: str) -> Optional[str]:
        with self._lock:
            return self._sessions.get(backend)

    def set_session(self, backend: str, session_id: Optional[str]) -> None:
        with self._lock:
            self._sessions[backend] = session_id or None
            self.save()

    def clear(self, backend: Optional[str] = None, history: bool = False) -> Dict[str, Optional[str]]:
        """Forget the conversation of one backend (or all); the next run starts a new one."""
        with self._lock:
            for b in ([backend] if backend else list(BACKENDS)):
                self._sessions[b] = None
                if history:
                    self._history = [e for e in self._history if e.get("backend") != b]
            self.save()
            return dict(self._sessions)

    def add_history(self, entry: dict) -> None:
        with self._lock:
            self._history.append(dict(entry))
            del self._history[:-self.history_limit]
            self.save()

    def history(self, backend: Optional[str] = None, limit: Optional[int] = None) -> List[dict]:
        with self._lock:
            out = [dict(e) for e in self._history if not backend or e.get("backend") == backend]
        if limit:
            out = out[-limit:]
        return out
