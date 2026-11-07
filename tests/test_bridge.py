"""Tests for the ChemDraw agent bridge (chemdrawmcp.bridge): HTTP/SSE server driving a *fake* Claude Code /
Codex CLI.  No real CLI, no ChemDraw; the whole file runs in a few seconds."""
from __future__ import annotations

import json
import os
import sys
import time
import urllib.error
import urllib.request

import pytest

from chemdrawmcp.bridge import agentcli, server, settings as cfg
from chemdrawmcp.bridge.server import BridgeServer, origin_allowed

# ---------------------------------------------------------------------------
# fake CLIs
# ---------------------------------------------------------------------------

FAKE_CLAUDE = r'''
import json, os, sys, time
argv = sys.argv[1:]
if "--version" in argv:
    print("9.9.9-fake"); sys.exit(0)
p = os.environ.get("FAKE_CLI_ARGV")
if p:
    with open(p, "w") as f:
        json.dump(argv, f)
sid = "sess-123"
if "--resume" in argv:
    sid = argv[argv.index("--resume") + 1]
def out(o):
    sys.stdout.write(json.dumps(o) + "\n"); sys.stdout.flush()
out({"type": "system", "subtype": "init", "session_id": sid, "model": "fake-model", "tools": ["mcp__chemdraw__status"]})
sleep = float(os.environ.get("FAKE_CLI_SLEEP", "0") or 0)
if sleep:
    time.sleep(sleep)
if os.environ.get("FAKE_CLI_FAIL"):
    sys.stderr.write("boom\n"); sys.stderr.flush(); sys.exit(2)
out({"type": "stream_event", "event": {"type": "content_block_delta", "index": 0, "delta": {"type": "text_delta", "text": "Hel"}}})
out({"type": "stream_event", "event": {"type": "content_block_delta", "index": 0, "delta": {"type": "text_delta", "text": "lo"}}})
out({"type": "assistant", "message": {"role": "assistant", "content": [
    {"type": "text", "text": "Hello"},
    {"type": "tool_use", "id": "tu1", "name": "mcp__chemdraw__status", "input": {}}]}, "session_id": sid})
out({"type": "user", "message": {"role": "user", "content": [
    {"type": "tool_result", "tool_use_id": "tu1", "content": [{"type": "text", "text": "{\"running\": true, \"documents\": 1}"}]}]}, "session_id": sid})
out({"type": "assistant", "message": {"role": "assistant", "content": [{"type": "text", "text": "ChemDraw is running."}]}, "session_id": sid})
out({"type": "result", "subtype": "success", "is_error": False, "result": "ChemDraw is running.", "duration_ms": 1234,
     "total_cost_usd": 0.01, "num_turns": 2, "usage": {"input_tokens": 10, "output_tokens": 5}, "session_id": sid})
'''

FAKE_CODEX = r'''
import json, os, sys
argv = sys.argv[1:]
if "--version" in argv:
    print("codex-cli 0.0.0-fake"); sys.exit(0)
p = os.environ.get("FAKE_CLI_ARGV")
if p:
    with open(p, "w") as f:
        json.dump(argv, f)
sid = argv[-2] if argv[:2] == ["exec", "resume"] else "thr-1"
def out(o):
    sys.stdout.write(json.dumps(o) + "\n"); sys.stdout.flush()
out({"type": "thread.started", "thread_id": sid})
out({"type": "turn.started"})
out({"type": "item.started", "item": {"id": "c1", "type": "mcp_tool_call", "server": "chemdraw", "tool": "status", "arguments": {}, "status": "in_progress"}})
out({"type": "item.completed", "item": {"id": "c1", "type": "mcp_tool_call", "server": "chemdraw", "tool": "status", "arguments": {},
     "status": "completed", "result": {"content": [{"type": "text", "text": "{\"running\": true}"}]}}})
out({"type": "item.completed", "item": {"id": "c2", "type": "command_execution", "command": "ls renders", "aggregated_output": "a.png\n", "exit_code": 0, "status": "completed"}})
out({"type": "item.completed", "item": {"id": "m1", "type": "agent_message", "text": "ChemDraw is running."}})
out({"type": "turn.completed", "usage": {"input_tokens": 12, "cached_input_tokens": 0, "output_tokens": 3}})
'''


def make_cli(tmp_path, name, body):
    """An executable (sh wrapper, or .cmd on Windows) running ``body`` with this interpreter."""
    script = tmp_path / (name + ".py")
    script.write_text(body, encoding="utf-8")
    if sys.platform.startswith("win"):
        exe = tmp_path / (name + ".cmd")
        exe.write_text('@"%s" -I "%s" %%*\r\n' % (sys.executable, script), encoding="utf-8")
    else:
        exe = tmp_path / name
        exe.write_text('#!/bin/sh\nexec "%s" -I "%s" "$@"\n' % (sys.executable, script), encoding="utf-8")
        exe.chmod(0o755)
    return str(exe)


# ---------------------------------------------------------------------------
# fixtures / client helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(h))
    for v in ("CHEMDRAWMCP_AGENT_SETTINGS", "CHEMDRAWMCP_CLAUDE_EXE", "CHEMDRAWMCP_CODEX_EXE", "CLAUDE_CLI",
              "CODEX_CLI", "CHEMDRAWMCP_BRIDGE_CWD", "FAKE_CLI_SLEEP", "FAKE_CLI_FAIL"):
        monkeypatch.delenv(v, raising=False)
    monkeypatch.setenv("FAKE_CLI_ARGV", str(tmp_path / "argv.json"))
    return h


@pytest.fixture
def fakes(tmp_path, monkeypatch, home):
    claude = make_cli(tmp_path, "fake_claude", FAKE_CLAUDE)
    codex = make_cli(tmp_path, "fake_codex", FAKE_CODEX)
    monkeypatch.setattr(agentcli, "find_claude", lambda: claude)
    monkeypatch.setattr(agentcli, "find_codex", lambda: codex)
    return {"claude": claude, "codex": codex, "argv": tmp_path / "argv.json"}


@pytest.fixture
def bridge(fakes):
    srv = BridgeServer(port=0, quiet=True).start()
    yield srv
    srv.stop()


_NO_TOKEN = object()


def call(srv, method, path, body=None, token=None, headers=None):
    """(status, json, headers) for one request; ``token=_NO_TOKEN`` sends none."""
    data = json.dumps(body).encode("utf-8") if isinstance(body, (dict, list)) else body
    req = urllib.request.Request(srv.url + path, data=data, method=method)
    if token is not _NO_TOKEN:
        req.add_header("X-Bridge-Token", token or srv.token)
    if data is not None:
        req.add_header("Content-Type", "application/json")
    for k, v in (headers or {}).items():
        req.add_header(k, v)
    try:
        with urllib.request.urlopen(req, timeout=5) as r:
            raw = r.read()
            return r.status, (json.loads(raw) if raw else None), r.headers
    except urllib.error.HTTPError as e:
        raw = e.read()
        return e.code, (json.loads(raw) if raw else None), e.headers


def read_sse(srv, run_id, after=None, timeout=8):
    """Consume /events/<run_id> (token in the query, like EventSource) until the ``end`` event.
    Returns [(kind, payload, id), ...] and the keep-alive comments seen."""
    req = urllib.request.Request("%s/events/%s?token=%s" % (srv.url, run_id, srv.token))
    if after is not None:
        req.add_header("Last-Event-ID", str(after))
    events, comments, cur = [], [], {}
    with urllib.request.urlopen(req, timeout=timeout) as r:
        assert r.status == 200
        assert r.headers["Content-Type"].startswith("text/event-stream")
        for raw in r:
            line = raw.decode("utf-8").rstrip("\r\n")
            if line == "":
                if "event" in cur:
                    events.append((cur["event"], json.loads(cur.get("data", "null")), cur.get("id")))
                    if cur["event"] == "end":
                        break
                cur = {}
            elif line.startswith(":"):
                comments.append(line)
            else:
                k, _, v = line.partition(":")
                cur[k] = v.lstrip()
    return events, comments


def run_to_end(srv, body):
    status, resp, _ = call(srv, "POST", "/run", body)
    assert status == 200, resp
    events, _ = read_sse(srv, resp["run_id"])
    return resp, events


def argv_of(fakes):
    return json.loads(fakes["argv"].read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def test_origin_policy():
    for ok in (None, "", "null", "file://", "chemdraw://addin", "app://local", "http://127.0.0.1",
               "http://127.0.0.1:5173", "http://localhost", "https://localhost:8443", "http://[::1]:3000"):
        assert origin_allowed(ok), ok
    for bad in ("https://evil.example", "http://evil.example:80", "https://127.0.0.1.evil.example",
                "http://localhost.evil.example"):
        assert not origin_allowed(bad), bad


def test_health_bridge_file_and_helpers(bridge, fakes, home):
    status, h, headers = call(bridge, "GET", "/health")
    assert status == 200 and h["ok"] is True
    assert h["backends"] == {"claude": fakes["claude"], "codex": fakes["codex"]}
    assert h["default_backend"] == "claude" and h["runs_active"] == 0
    assert set(h["settings"]) >= {"claude", "codex", "default_backend"}
    assert h["version"] and headers["Cache-Control"] == "no-store"
    assert headers["Access-Control-Allow-Origin"] == "*"

    info_path = home / "bridge.json"
    assert info_path.exists()
    info = json.loads(info_path.read_text(encoding="utf-8"))
    assert info["port"] == bridge.port and info["token"] == bridge.token and info["pid"] == os.getpid()
    assert isinstance(info["started"], float)
    if not sys.platform.startswith("win"):
        assert (info_path.stat().st_mode & 0o777) == 0o600
    assert server.bridge_info()["port"] == bridge.port
    assert server.is_bridge_running() is True
    assert server.bridge_health()["ok"] is True

    bridge.stop()
    assert not info_path.exists()
    assert server.bridge_info() is None
    assert server.is_bridge_running() is False


def test_token_reuse_and_new_token(fakes, home):
    a = BridgeServer(port=0, quiet=True).start()
    tok = a.token
    a.stop()
    # stop() removed bridge.json; simulate a stale file from a crashed bridge -> token is reused
    cfg.write_json(str(home / "bridge.json"), {"port": 1, "token": tok, "pid": 0, "started": 0})
    b = BridgeServer(port=0, quiet=True)
    assert b.token == tok
    c = BridgeServer(port=0, quiet=True, new_token=True)
    assert c.token != tok and len(c.token) >= 24


def test_claude_run_stream_sessions_history(bridge, fakes, home):
    status, resp, _ = call(bridge, "POST", "/run", {"backend": "claude", "prompt": "Draw caffeine",
                                                    "context": "Selection SMILES: CCO"})
    assert status == 200 and resp["backend"] == "claude" and resp["session_id"] is None
    run_id = resp["run_id"]
    events, comments = read_sse(bridge, run_id)
    kinds = [k for k, _, _ in events]
    assert kinds == ["session", "delta", "delta", "text", "tool", "tool_result", "text", "stats", "done", "end"]
    assert [e[2] for e in events[:-1]] == [str(i) for i in range(1, 10)]
    payload = {k: p for k, p, _ in events}
    assert events[0][1] == "sess-123"
    assert events[1][1] == "Hel" and events[2][1] == "lo"
    assert events[3][1] == "Hello"
    assert events[4][1].startswith("chemdraw/status {")
    assert "running" in events[5][1]
    assert payload["stats"]["total_cost_usd"] == 0.01 and payload["stats"]["usage"]["input_tokens"] == 10
    assert payload["done"] == "ChemDraw is running."
    assert payload["end"] == {"run_id": run_id, "status": "done", "final_text": "ChemDraw is running.",
                              "session_id": "sess-123", "error": None}
    assert comments == []  # no keep-alives needed for a fast run

    status, r, _ = call(bridge, "GET", "/runs/" + run_id)
    assert status == 200 and r["status"] == "done" and r["events_count"] == 9
    assert r["final_text"] == "ChemDraw is running." and r["session_id"] == "sess-123"
    assert call(bridge, "GET", "/sessions")[1] == {"claude": "sess-123", "codex": None}
    status, hist, _ = call(bridge, "GET", "/history?backend=claude")
    assert status == 200 and len(hist) == 1
    assert hist[0]["prompt"] == "Draw caffeine" and hist[0]["final_text"] == "ChemDraw is running."
    assert hist[0]["backend"] == "claude" and hist[0]["session_id"] == "sess-123" and hist[0]["started"] > 0
    assert call(bridge, "GET", "/history?backend=codex")[1] == []
    persisted = json.loads((home / "agent_sessions.json").read_text(encoding="utf-8"))
    assert persisted["sessions"]["claude"] == "sess-123" and persisted["history"][0]["prompt"] == "Draw caffeine"
    assert call(bridge, "GET", "/health")[1]["runs_active"] == 0

    argv = argv_of(fakes)
    assert argv[:2] == ["-p", "Draw caffeine"]
    assert argv[argv.index("--allowedTools") + 1] == "mcp__chemdraw"
    assert "--include-partial-messages" in argv and "--output-format" in argv
    sysprompt = argv[argv.index("--append-system-prompt") + 1]
    assert "Agent Chat" in sysprompt and "quit_chemdraw" in sysprompt and "Selection SMILES: CCO" in sysprompt
    assert "--resume" not in argv and "--model" not in argv and "--effort" not in argv

    # replay from the start and resume after an id (EventSource reconnect with Last-Event-ID)
    again, _ = read_sse(bridge, run_id)
    assert [k for k, _, _ in again] == kinds
    tail, _ = read_sse(bridge, run_id, after=8)
    assert [k for k, _, _ in tail] == ["done", "end"]

    # second turn resumes the stored session
    resp2, ev2 = run_to_end(bridge, {"backend": "claude", "prompt": "again", "allow_shell": True})
    assert resp2["session_id"] == "sess-123" and ev2[-1][1]["status"] == "done"
    argv = argv_of(fakes)
    assert argv[argv.index("--resume") + 1] == "sess-123"
    assert argv[argv.index("--allowedTools") + 1] == "mcp__chemdraw,Bash,Edit,Write"
    assert len(call(bridge, "GET", "/history")[1]) == 2

    # /clear forgets the session -> next run starts fresh; explicit session_id wins; "new" forces fresh
    status, c, _ = call(bridge, "POST", "/clear", {"backend": "claude"})
    assert status == 200 and c["sessions"]["claude"] is None
    assert call(bridge, "GET", "/sessions")[1]["claude"] is None
    resp3, _ = run_to_end(bridge, {"backend": "claude", "prompt": "fresh"})
    assert resp3["session_id"] is None and "--resume" not in argv_of(fakes)
    resp4, _ = run_to_end(bridge, {"backend": "claude", "prompt": "x", "session_id": "sess-custom"})
    assert resp4["session_id"] == "sess-custom" and argv_of(fakes)[-1] == "sess-custom"
    assert call(bridge, "GET", "/sessions")[1]["claude"] == "sess-custom"
    resp5, _ = run_to_end(bridge, {"backend": "claude", "prompt": "x", "new": True})
    assert resp5["session_id"] is None and "--resume" not in argv_of(fakes)
    assert len(call(bridge, "GET", "/runs")[1]) == 5


def test_codex_run_stream(bridge, fakes, home):
    resp, events = run_to_end(bridge, {"backend": "codex", "prompt": "status?"})
    kinds = [k for k, _, _ in events]
    assert kinds == ["session", "tool", "tool_result", "tool", "tool_result", "text", "stats", "done", "end"]
    assert events[0][1] == "thr-1"
    assert events[1][1] == "chemdraw/status {}"
    assert events[2][1] == '{"running": true}'
    assert events[3][1] == "shell ls renders" and events[4][1] == "a.png\n"
    assert events[5][1] == "ChemDraw is running."
    assert events[6][1] == {"usage": {"input_tokens": 12, "cached_input_tokens": 0, "output_tokens": 3}}
    assert events[7][1] == ""
    assert events[8][1]["status"] == "done" and events[8][1]["final_text"] == "ChemDraw is running."
    assert call(bridge, "GET", "/sessions")[1]["codex"] == "thr-1"
    argv = argv_of(fakes)
    assert argv[:4] == ["exec", "--json", "--approve-for-me", "--skip-git-repo-check"]
    assert "Agent Chat" in argv[-1] and argv[-1].endswith("User request: status?")

    resp2, ev2 = run_to_end(bridge, {"backend": "codex", "prompt": "again", "context": "front doc: Untitled-1"})
    assert resp2["session_id"] == "thr-1" and ev2[0][1] == "thr-1"
    argv = argv_of(fakes)
    assert argv[:2] == ["exec", "resume"] and argv[-2] == "thr-1"
    assert "Agent Chat" not in argv[-1] and "front doc: Untitled-1" in argv[-1] and argv[-1].endswith("User request: again")
    hist = call(bridge, "GET", "/history?backend=codex")[1]
    assert [h["prompt"] for h in hist] == ["status?", "again"]


def test_settings_round_trip(bridge, fakes, home):
    status, s, _ = call(bridge, "GET", "/settings")
    assert status == 200 and s["claude"] == {"model": "", "effort": ""} and s["default_backend"] == "claude"
    status, s, _ = call(bridge, "POST", "/settings", {"backend": "claude", "model": "opus", "effort": "high",
                                                      "default_backend": "codex"})
    assert status == 200
    assert s["claude"] == {"model": "opus", "effort": "high"} and s["default_backend"] == "codex"
    assert call(bridge, "GET", "/settings")[1]["claude"]["model"] == "opus"
    assert call(bridge, "GET", "/health")[1]["default_backend"] == "codex"
    on_disk = json.loads((home / "agent_settings.json").read_text(encoding="utf-8"))
    assert on_disk["claude"] == {"model": "opus", "effort": "high"} and on_disk["default_backend"] == "codex"
    assert call(bridge, "POST", "/settings", {"backend": "claude", "effort": "turbo"})[0] == 400
    assert call(bridge, "POST", "/settings", {"backend": "gemini", "model": "x"})[0] == 400
    assert call(bridge, "POST", "/settings", {"model": "x"})[0] == 400

    # settings reach the CLI; a run-level override wins; backend defaults to default_backend (codex now)
    run_to_end(bridge, {"backend": "claude", "prompt": "hi"})
    argv = argv_of(fakes)
    assert argv[argv.index("--model") + 1] == "opus" and argv[argv.index("--effort") + 1] == "high"
    run_to_end(bridge, {"backend": "claude", "prompt": "hi", "effort": "low", "model": "haiku"})
    argv = argv_of(fakes)
    assert argv[argv.index("--model") + 1] == "haiku" and argv[argv.index("--effort") + 1] == "low"
    resp, _ = run_to_end(bridge, {"prompt": "hi"})
    assert resp["backend"] == "codex"
    assert call(bridge, "POST", "/settings", {"backend": "codex", "effort": "minimal", "model": "gpt-5"})[0] == 200
    run_to_end(bridge, {"backend": "codex", "prompt": "hi"})
    argv = argv_of(fakes)
    assert argv[argv.index("-m") + 1] == "gpt-5" and 'model_reasoning_effort="minimal"' in argv


def test_auth_and_origin(bridge):
    assert call(bridge, "GET", "/health", token=_NO_TOKEN)[0] == 401
    assert call(bridge, "GET", "/health", token="wrong-token")[0] == 401
    assert call(bridge, "POST", "/run", {"prompt": "x"}, token=_NO_TOKEN)[0] == 401
    status, h, _ = call(bridge, "GET", "/health?token=" + bridge.token, token=_NO_TOKEN)
    assert status == 200 and h["ok"]

    assert call(bridge, "GET", "/health", headers={"Origin": "https://evil.example"})[0] == 403
    assert call(bridge, "GET", "/health", headers={"Origin": "http://evil.example:80"})[0] == 403
    assert call(bridge, "GET", "/health", headers={"Origin": "http://localhost.evil.example"})[0] == 403
    for origin, acao in (("null", "*"), ("file://", "*"), ("chemdraw://addin", "*"),
                         ("http://localhost:5173", "http://localhost:5173"), ("http://127.0.0.1", "http://127.0.0.1")):
        status, _, headers = call(bridge, "GET", "/health", headers={"Origin": origin})
        assert status == 200, origin
        assert headers["Access-Control-Allow-Origin"] == acao, origin
        assert headers["Cache-Control"] == "no-store"

    # preflight needs no token
    status, body, headers = call(bridge, "OPTIONS", "/run", token=_NO_TOKEN,
                                 headers={"Origin": "null", "Access-Control-Request-Method": "POST",
                                          "Access-Control-Request-Headers": "x-bridge-token,content-type"})
    assert status == 204 and body is None
    assert headers["Access-Control-Allow-Origin"] == "*"
    assert "X-Bridge-Token" in headers["Access-Control-Allow-Headers"]
    assert "Content-Type" in headers["Access-Control-Allow-Headers"]
    assert "POST" in headers["Access-Control-Allow-Methods"]
    assert headers["Access-Control-Allow-Private-Network"] == "true"
    assert headers["Cache-Control"] == "no-store"
    assert call(bridge, "OPTIONS", "/run", token=_NO_TOKEN, headers={"Origin": "https://evil.example"})[0] == 403


def test_busy_cancel_and_failure(bridge, fakes, monkeypatch):
    monkeypatch.setenv("FAKE_CLI_SLEEP", "1")
    t0 = time.time()
    status, a, _ = call(bridge, "POST", "/run", {"backend": "claude", "prompt": "slow"})
    assert status == 200
    status, b, _ = call(bridge, "POST", "/run", {"backend": "claude", "prompt": "second"})
    assert status == 409 and b["run_id"] == a["run_id"] and "busy" in b["error"]
    assert call(bridge, "GET", "/health")[1]["runs_active"] == 1
    assert call(bridge, "GET", "/runs/" + a["run_id"])[1]["status"] == "running"
    # the other backend is free
    monkeypatch.delenv("FAKE_CLI_SLEEP")
    assert call(bridge, "POST", "/run", {"backend": "codex", "prompt": "parallel"})[0] == 200

    # let the fake print its init line (it then sleeps 1 s) so the cancel hits a live process
    deadline = time.time() + 5
    while call(bridge, "GET", "/runs/" + a["run_id"])[1]["events_count"] < 1 and time.time() < deadline:
        time.sleep(0.02)
    status, c, _ = call(bridge, "POST", "/cancel/" + a["run_id"])
    assert status == 200 and c["ok"] and c["was_running"]
    events, _ = read_sse(bridge, a["run_id"])
    kinds = [k for k, _, _ in events]
    assert kinds == ["session", "error", "done", "end"]
    assert events[1][1] == "cancelled" and events[-1][1]["status"] == "cancelled"
    assert time.time() - t0 < 3.0  # cancelled well before the fake's 1 s sleep + overhead would matter
    assert call(bridge, "GET", "/runs/" + a["run_id"])[1]["status"] == "cancelled"
    assert call(bridge, "POST", "/cancel/" + a["run_id"])[1]["was_running"] is False
    assert call(bridge, "POST", "/cancel/nope")[0] == 404
    hist = call(bridge, "GET", "/history?backend=claude")[1]
    assert hist[-1]["status"] == "cancelled" and hist[-1]["prompt"] == "slow"
    # the session id of a cancelled turn is still stored (the CLI-side conversation exists)
    assert call(bridge, "GET", "/sessions")[1]["claude"] == "sess-123"

    # backend free again; a crashing CLI is reported as an error run
    monkeypatch.setenv("FAKE_CLI_FAIL", "1")
    resp, events = run_to_end(bridge, {"backend": "claude", "prompt": "crash"})
    kinds = [k for k, _, _ in events]
    assert kinds == ["session", "error", "done", "end"]
    assert "exited with code 2" in events[1][1] and "boom" in events[1][1]
    assert events[-1][1]["status"] == "error"
    assert call(bridge, "GET", "/runs/" + resp["run_id"])[1]["status"] == "error"
    monkeypatch.delenv("FAKE_CLI_FAIL")
    # wait for the parallel codex run to finish so the fixture teardown is clean
    read_sse(bridge, call(bridge, "GET", "/runs")[1][1]["run_id"])
    assert call(bridge, "GET", "/health")[1]["runs_active"] == 0


def test_errors_and_missing_backend(bridge, fakes, monkeypatch):
    assert call(bridge, "GET", "/events/nope?token=" + bridge.token, token=_NO_TOKEN)[0] == 404
    assert call(bridge, "GET", "/runs/nope")[0] == 404
    assert call(bridge, "GET", "/nope")[0] == 404
    assert call(bridge, "POST", "/run", b"{not json")[0] == 400
    assert call(bridge, "POST", "/run", {})[0] == 400
    assert call(bridge, "POST", "/run", {"prompt": "   "})[0] == 400
    assert call(bridge, "POST", "/run", {"backend": "gemini", "prompt": "x"})[0] == 400
    assert call(bridge, "POST", "/run", {"backend": "claude", "prompt": "x", "effort": "turbo"})[0] == 400
    assert call(bridge, "GET", "/history?backend=gemini")[0] == 400
    assert call(bridge, "POST", "/clear", {"backend": "gemini"})[0] == 400
    monkeypatch.setattr(agentcli, "find_codex", lambda: None)
    status, body, _ = call(bridge, "POST", "/run", {"backend": "codex", "prompt": "x"})
    assert status == 503 and "CHEMDRAWMCP_CODEX_EXE" in body["error"]
    assert call(bridge, "GET", "/health", headers={})[1]["backends"]["claude"] == fakes["claude"]


def test_env_exe_override_and_default_backend(tmp_path, monkeypatch, home):
    claude = make_cli(tmp_path, "fake_claude", FAKE_CLAUDE)
    monkeypatch.setenv("CHEMDRAWMCP_CLAUDE_EXE", claude)
    monkeypatch.setattr(agentcli, "find_codex", lambda: None)
    assert agentcli.find_claude() == claude
    assert agentcli.available_backends() == ["claude"]
    assert agentcli.default_backend() == "claude"
    assert server.bridge_info() is None and server.is_bridge_running() is False
    srv = BridgeServer(port=0, quiet=True).start()
    try:
        h = call(srv, "GET", "/health")[1]
        assert h["backends"] == {"claude": claude, "codex": None} and h["default_backend"] == "claude"
        resp, events = run_to_end(srv, {"prompt": "hello"})
        assert resp["backend"] == "claude" and events[-1][1]["status"] == "done"
    finally:
        srv.stop()


def test_store_survives_restart(fakes, home):
    a = BridgeServer(port=0, quiet=True).start()
    try:
        run_to_end(a, {"backend": "claude", "prompt": "first"})
    finally:
        a.stop()
    b = BridgeServer(port=0, quiet=True).start()
    try:
        assert call(b, "GET", "/sessions")[1]["claude"] == "sess-123"
        assert [h["prompt"] for h in call(b, "GET", "/history")[1]] == ["first"]
        resp, _ = run_to_end(b, {"backend": "claude", "prompt": "second"})
        assert resp["session_id"] == "sess-123"
        assert argv_of(fakes)[argv_of(fakes).index("--resume") + 1] == "sess-123"
    finally:
        b.stop()


def test_main_check_and_parse(fakes, capsys):
    assert server.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert fakes["claude"] in out and fakes["codex"] in out and "default: claude" in out
