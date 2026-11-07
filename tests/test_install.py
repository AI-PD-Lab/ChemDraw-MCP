"""Tests for the installer (MCP client registration) and ChemDraw discovery.

Everything that touches config files runs against a temporary HOME; the only
test that looks at the real machine is the macOS-only ChemDraw discovery one.
"""
import json
import os
import plistlib
import re
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemdrawmcp import discovery, install  # noqa: E402

try:
    import tomllib  # py3.11+
except ImportError:  # pragma: no cover - py3.10
    import tomli as tomllib


def load_toml(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


FOREIGN = {"command": "npx", "args": ["-y", "some-other-chemdraw-server"]}


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """Point every config location (and HOME/APPDATA) at tmp_path."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.setenv("APPDATA", str(home / "AppData" / "Roaming"))
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(home / ".chemdrawmcp"))
    monkeypatch.delenv("CHEMDRAWMCP_OUTPUT_DIR", raising=False)
    claude = str(home / ".claude.json")
    codex = str(home / ".codex" / "config.toml")
    desktop = str(home / "Claude" / "claude_desktop_config.json")
    monkeypatch.setattr(install, "claude_config_path", lambda: claude)
    monkeypatch.setattr(install, "codex_config_path", lambda: codex)
    monkeypatch.setattr(install, "claude_desktop_config_path", lambda: desktop if os.path.exists(desktop) else None)
    return {"home": str(home), "claude": claude, "codex": codex, "desktop": desktop}


# ---------------------------------------------------------------------------
# Claude Code (~/.claude.json)
# ---------------------------------------------------------------------------

def test_register_claude_into_empty_config(cfg):
    assert not os.path.exists(cfg["claude"])
    install.register_claude(dry=False)
    data = json.loads(read(cfg["claude"]))
    entry = data["mcpServers"]["chemdraw"]
    assert entry["type"] == "stdio"
    assert install._is_ours(entry)
    cmd, args = install.server_entry()
    assert entry["command"] == cmd and entry["args"] == args
    assert not os.path.exists(cfg["claude"] + install.BACKUP_SUFFIX)  # nothing to back up


def test_register_claude_preserves_other_servers(cfg):
    original = {"numStartups": 3, "mcpServers": {"foo": {"type": "stdio", "command": "foo-mcp", "args": []}},
                "projects": {"/x": {"allowedTools": []}}}
    write(cfg["claude"], json.dumps(original))
    install.register_claude(dry=False)
    data = json.loads(read(cfg["claude"]))
    assert data["numStartups"] == 3 and data["projects"] == original["projects"]
    assert data["mcpServers"]["foo"] == original["mcpServers"]["foo"]
    assert "chemdraw" in data["mcpServers"]
    assert json.loads(read(cfg["claude"] + install.BACKUP_SUFFIX)) == original
    # idempotent
    install.register_claude(dry=False)
    assert json.loads(read(cfg["claude"])) == data


def test_uninstall_claude_removes_only_ours(cfg):
    write(cfg["claude"], json.dumps({"mcpServers": {"foo": {"command": "foo-mcp", "args": []}}}))
    install.register_claude(dry=False)
    install.register_claude(dry=False, remove=True)
    data = json.loads(read(cfg["claude"]))
    assert "chemdraw" not in data["mcpServers"]
    assert data["mcpServers"]["foo"] == {"command": "foo-mcp", "args": []}
    # a foreign "chemdraw" entry is never removed
    write(cfg["claude"], json.dumps({"mcpServers": {"chemdraw": FOREIGN}}))
    install.register_claude(dry=False, remove=True)
    assert json.loads(read(cfg["claude"]))["mcpServers"]["chemdraw"] == FOREIGN


def test_uninstall_claude_without_config_creates_nothing(cfg):
    install.register_claude(dry=False, remove=True)
    assert not os.path.exists(cfg["claude"])


def test_is_ours():
    assert install._is_ours({"command": "/x/.venv/bin/chemdraw-mcp", "args": []})
    assert install._is_ours({"command": r"C:\x\.venv\Scripts\chemdraw-mcp.exe", "args": []})
    assert install._is_ours({"command": "/x/.venv/bin/python", "args": ["-m", "chemdrawmcp.server"]})
    assert not install._is_ours(FOREIGN)
    assert not install._is_ours({"command": "/x/.venv/bin/pymol-mcp", "args": []})
    assert not install._is_ours({})
    assert not install._is_ours(None)


# ---------------------------------------------------------------------------
# Codex CLI (~/.codex/config.toml)
# ---------------------------------------------------------------------------

CODEX_BASE = 'model = "o3"\napproval_policy = "never"\n\n[mcp_servers.other]\ncommand = "other-mcp"\nargs = ["--flag"]\n\n[projects."/tmp/x"]\ntrust_level = "trusted"\n'


def test_codex_add_replace_remove(cfg, monkeypatch):
    write(cfg["codex"], CODEX_BASE)
    install.register_codex(dry=False)
    data = load_toml(cfg["codex"])
    cmd, args = install.server_entry()
    assert data["model"] == "o3" and data["approval_policy"] == "never"
    assert data["mcp_servers"]["other"] == {"command": "other-mcp", "args": ["--flag"]}
    assert data["projects"]["/tmp/x"] == {"trust_level": "trusted"}
    assert data["mcp_servers"]["chemdraw"] == {"command": cmd, "args": args}
    assert read(cfg["codex"] + install.BACKUP_SUFFIX) == CODEX_BASE

    # replace: a different server command updates the block in place, nothing is duplicated
    monkeypatch.setattr(install, "server_entry", lambda: ("/new/.venv/bin/chemdraw-mcp", []))
    install.register_codex(dry=False)
    text = read(cfg["codex"])
    assert text.count("[mcp_servers.chemdraw]") == 1
    data = load_toml(cfg["codex"])
    assert data["mcp_servers"]["chemdraw"] == {"command": "/new/.venv/bin/chemdraw-mcp", "args": []}
    assert data["mcp_servers"]["other"]["command"] == "other-mcp"
    assert data["projects"]["/tmp/x"]["trust_level"] == "trusted"

    # remove: only our table goes, everything else survives and the file stays valid TOML
    install.register_codex(dry=False, remove=True)
    data = load_toml(cfg["codex"])
    assert "chemdraw" not in data["mcp_servers"]
    assert data["mcp_servers"]["other"] == {"command": "other-mcp", "args": ["--flag"]}
    assert data["model"] == "o3" and data["projects"]["/tmp/x"]["trust_level"] == "trusted"
    assert "[mcp_servers.chemdraw]" not in read(cfg["codex"])


def test_codex_block_in_the_middle_and_subtables(cfg):
    write(cfg["codex"], '[mcp_servers.chemdraw]\ncommand = "old"\nargs = []\n\n[mcp_servers.chemdraw.env]\nFOO = "1"\n\n'
                        '[mcp_servers.zzz]\ncommand = "zzz"\nextra = "x"\n')
    install.register_codex(dry=False)
    data = load_toml(cfg["codex"])
    assert "env" not in data["mcp_servers"]["chemdraw"]  # stale sub-table dropped together with the old block
    assert data["mcp_servers"]["chemdraw"]["command"] != "old"
    assert data["mcp_servers"]["zzz"] == {"command": "zzz", "extra": "x"}
    assert read(cfg["codex"]).count("[mcp_servers.chemdraw") == 1
    install.register_codex(dry=False, remove=True)
    data = load_toml(cfg["codex"])
    assert data == {"mcp_servers": {"zzz": {"command": "zzz", "extra": "x"}}}


def test_codex_creates_config_when_missing(cfg):
    install.register_codex(dry=False)
    data = load_toml(cfg["codex"])
    assert set(data) == {"mcp_servers"} and set(data["mcp_servers"]) == {"chemdraw"}
    assert not os.path.exists(cfg["codex"] + install.BACKUP_SUFFIX)


def test_codex_uninstall_without_config_creates_nothing(cfg):
    install.register_codex(dry=False, remove=True)
    assert not os.path.exists(cfg["codex"])


def test_toml_str_escapes_windows_paths():
    s = install._toml_str(r"C:\Users\me\.venv\Scripts\chemdraw-mcp.exe")
    assert tomllib.loads("x = %s" % s)["x"] == r"C:\Users\me\.venv\Scripts\chemdraw-mcp.exe"


# ---------------------------------------------------------------------------
# Claude Desktop
# ---------------------------------------------------------------------------

def test_claude_desktop_only_touched_when_config_exists(cfg):
    install.register_claude_desktop(dry=False)
    assert not os.path.exists(cfg["desktop"])
    write(cfg["desktop"], json.dumps({"mcpServers": {"foo": {"command": "foo", "args": []}}, "theme": "dark"}))
    install.register_claude_desktop(dry=False)
    data = json.loads(read(cfg["desktop"]))
    assert data["theme"] == "dark" and data["mcpServers"]["foo"] == {"command": "foo", "args": []}
    assert install._is_ours(data["mcpServers"]["chemdraw"])
    assert "type" not in data["mcpServers"]["chemdraw"]  # Desktop entries are just command/args
    assert os.path.exists(cfg["desktop"] + install.BACKUP_SUFFIX)


def test_claude_desktop_keeps_foreign_chemdraw_entry(cfg):
    write(cfg["desktop"], json.dumps({"mcpServers": {"chemdraw": FOREIGN}}))
    install.register_claude_desktop(dry=False)
    data = json.loads(read(cfg["desktop"]))
    assert data["mcpServers"]["chemdraw"] == FOREIGN            # untouched
    assert install._is_ours(data["mcpServers"][install.ALT_SERVER_NAME])  # ours under the renamed key
    install.register_claude_desktop(dry=False, remove=True)
    data = json.loads(read(cfg["desktop"]))
    assert data["mcpServers"] == {"chemdraw": FOREIGN}


# ---------------------------------------------------------------------------
# dry run / main
# ---------------------------------------------------------------------------

def _snapshot(root):
    out = {}
    for d, _dirs, files in os.walk(root):
        for fn in files:
            p = os.path.join(d, fn)
            with open(p, "rb") as f:
                out[p] = f.read()
    return out


def test_dry_run_writes_nothing(cfg, capsys):
    write(cfg["claude"], json.dumps({"mcpServers": {"foo": {"command": "foo", "args": []}}}))
    write(cfg["codex"], CODEX_BASE)
    write(cfg["desktop"], json.dumps({"mcpServers": {}}))
    before = _snapshot(cfg["home"])
    assert install.main(["--dry-run"]) == 0
    assert install.main(["--dry-run", "--uninstall"]) == 0
    assert install.main(["--dry-run", "--codex"]) == 0
    assert _snapshot(cfg["home"]) == before
    assert not any(p.endswith(install.BACKUP_SUFFIX) for p in _snapshot(cfg["home"]))
    out = capsys.readouterr().out
    assert "chemdraw ->" in out and "dry run" in out


def test_main_registers_everything_and_prints_verification(cfg, capsys):
    write(cfg["desktop"], json.dumps({"mcpServers": {}}))
    assert install.main([]) == 0
    out = capsys.readouterr().out
    assert "/mcp" in out and "codex mcp list" in out
    assert "chemdraw" in json.loads(read(cfg["claude"]))["mcpServers"]
    assert "chemdraw" in load_toml(cfg["codex"])["mcp_servers"]
    assert "chemdraw" in json.loads(read(cfg["desktop"]))["mcpServers"]
    assert install.main(["--uninstall"]) == 0
    assert "chemdraw" not in json.loads(read(cfg["claude"]))["mcpServers"]
    assert "chemdraw" not in load_toml(cfg["codex"]).get("mcp_servers", {})
    assert "chemdraw" not in json.loads(read(cfg["desktop"]))["mcpServers"]


def test_main_check_reports(cfg, capsys):
    assert install.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert "== chemdrawmcp " in out
    assert "python:" in out and "chemdraw:" in out and "rdkit:" in out and "server:" in out and "home:" in out
    assert not os.path.exists(cfg["claude"])  # --check never registers


def test_server_entry_is_ours():
    cmd, args = install.server_entry()
    assert os.path.isabs(cmd)
    assert install._is_ours({"command": cmd, "args": args})


# ---------------------------------------------------------------------------
# discovery
# ---------------------------------------------------------------------------

def test_output_dir_honours_env_var(tmp_path, monkeypatch):
    monkeypatch.setenv("CHEMDRAWMCP_OUTPUT_DIR", str(tmp_path / "out" / "renders"))
    assert discovery.output_dir() == str(tmp_path / "out" / "renders")
    assert os.path.isdir(tmp_path / "out" / "renders")
    monkeypatch.delenv("CHEMDRAWMCP_OUTPUT_DIR")
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path / "home"))
    assert discovery.home_dir() == str(tmp_path / "home")
    assert discovery.output_dir() == str(tmp_path / "home" / "renders")
    assert os.path.isdir(tmp_path / "home" / "renders")


def _fake_app(root, name, version, bundle_id="com.revvity.ChemDraw"):
    app = root / name
    (app / "Contents" / "MacOS").mkdir(parents=True)
    with open(app / "Contents" / "Info.plist", "wb") as f:
        plistlib.dump({"CFBundleIdentifier": bundle_id, "CFBundleShortVersionString": version,
                       "CFBundleExecutable": "ChemDraw"}, f)
    return str(app)


def test_find_chemdraw_picks_highest_version_bundle(tmp_path, monkeypatch):
    monkeypatch.delenv("CHEMDRAW_APP", raising=False)
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    monkeypatch.setattr(discovery, "MAC", True)
    monkeypatch.setattr(discovery, "WIN", False)
    monkeypatch.setattr(discovery, "MAC_APP_DIRS", [str(tmp_path)])
    monkeypatch.setattr(discovery, "_mdfind_apps", lambda: [])
    plain = _fake_app(tmp_path, "ChemDraw.app", "21.0.0.1")
    old = _fake_app(tmp_path, "ChemDraw 22.0.0.app", "22.0.0.5")
    new = _fake_app(tmp_path, "ChemDraw 23.1.2.app", "23.1.2.7")
    _fake_app(tmp_path, "ChemDraw Collections.app", "99.0", bundle_id="com.revvity.chemdrawcollections")
    assert discovery.find_chemdraw() == new
    assert discovery.chemdraw_version(new) == "23.1.2.7"
    assert discovery.chemdraw_version(os.path.join(new, "Contents", "MacOS", "ChemDraw")) == "23.1.2.7"
    assert discovery.macos_bundle_id(new) == "com.revvity.ChemDraw"
    assert discovery.version_tuple("ChemDraw 23.1.2.app") == (23, 1, 2)
    assert discovery.version_tuple("ChemDraw.app") == ()
    # a versionless name falls back to Info.plist: "ChemDraw Professional.app" at 30.0 beats 23.1.2
    pro = _fake_app(tmp_path, "ChemDraw Professional.app", "30.0.0.1")
    assert discovery.find_chemdraw() == pro
    # the environment variable wins over everything
    monkeypatch.setenv("CHEMDRAW_APP", old)
    assert discovery.find_chemdraw() == old
    monkeypatch.setenv("CHEMDRAW_APP", str(tmp_path / "does-not-exist.app"))
    assert discovery.find_chemdraw() == pro  # non-existent -> ignored
    assert plain != new


def test_mac_name_filter():
    ok = discovery._mac_name_ok
    assert ok("/Applications/ChemDraw.app")
    assert ok("/Applications/ChemDraw 23.1.2.app")
    assert ok("/Applications/ChemDraw Professional 22.0.app")
    assert not ok("/Applications/ChemDraw Collections.app")
    assert not ok("/Applications/ChemOffice.app")


def test_windows_version_from_path():
    assert discovery._version_from_path(r"C:\Program Files\PerkinElmerInformatics\ChemOffice 20.1\ChemDraw\ChemDraw.exe") == "20.1"
    assert discovery._version_from_path(r"C:\Program Files\Revvity Signals\ChemDraw 23.1\ChemDraw.exe") == "23.1"
    assert discovery._version_from_path(r"C:\Program Files\Revvity Signals\ChemDraw\ChemDraw.exe") is None


def test_find_chemdraw_linux_returns_none(tmp_path, monkeypatch):
    monkeypatch.delenv("CHEMDRAW_APP", raising=False)
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    monkeypatch.setattr(discovery, "MAC", False)
    monkeypatch.setattr(discovery, "WIN", False)
    assert discovery.find_chemdraw() is None
    exe = tmp_path / "ChemDraw"
    exe.write_text("")
    monkeypatch.setenv("CHEMDRAW_EXE", str(exe))
    assert discovery.find_chemdraw() == str(exe)


@pytest.mark.skipif(sys.platform != "darwin", reason="looks for the real ChemDraw.app on macOS")
def test_find_chemdraw_on_this_mac(monkeypatch):
    monkeypatch.delenv("CHEMDRAW_APP", raising=False)
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    app = discovery.find_chemdraw()
    if app is None:
        pytest.skip("ChemDraw is not installed on this Mac")
    assert "ChemDraw" in app and app.endswith(".app") and "Collections" not in app
    assert os.path.isdir(app)
    version = discovery.chemdraw_version(app)
    assert version and re.match(r"^\d+(\.\d+)+$", version)
    assert discovery.macos_bundle_id(app) in discovery.chemdraw_bundle_ids()


# ---------------------------------------------------------------------------
# regressions: duplicate Codex definitions, pre-write validation, TOML strings, backups, env overrides
# ---------------------------------------------------------------------------

CODEX_INLINE = 'model = "o3"\n\n[mcp_servers]\nchemdraw = { command = "old", args = [] }\n'
CODEX_DOTTED = 'model = "o3"\nmcp_servers.chemdraw.command = "old"\n'


@pytest.mark.parametrize("parser", [True, False], ids=["parser", "no-parser"])
@pytest.mark.parametrize("existing", [CODEX_INLINE, CODEX_DOTTED], ids=["inline-table", "dotted-key"])
def test_codex_refuses_to_duplicate_a_non_header_definition(cfg, capsys, monkeypatch, existing, parser):
    if not parser:
        monkeypatch.setattr(install, "_toml_module", lambda: None)  # python 3.10 without tomli
    write(cfg["codex"], existing)
    install.register_codex(dry=False)
    assert read(cfg["codex"]) == existing  # untouched and still valid TOML
    assert load_toml(cfg["codex"])["mcp_servers"]["chemdraw"]["command"] == "old"
    assert not os.path.exists(cfg["codex"] + install.BACKUP_SUFFIX)
    out = capsys.readouterr().out
    assert "already defined" in out and "saved" not in out


@pytest.mark.parametrize("parser", [True, False], ids=["parser", "no-parser"])
def test_codex_replaces_a_quoted_header_in_place(cfg, monkeypatch, parser):
    if not parser:
        monkeypatch.setattr(install, "_toml_module", lambda: None)
    write(cfg["codex"], '[mcp_servers."chemdraw"]\ncommand = "old"\nargs = []\n\n[mcp_servers.zzz]\ncommand = "zzz"\n')
    install.register_codex(dry=False)
    data = load_toml(cfg["codex"])
    assert install._is_ours(data["mcp_servers"]["chemdraw"])
    assert data["mcp_servers"]["zzz"] == {"command": "zzz"}
    assert read(cfg["codex"]).count("[mcp_servers.") == 2
    install.register_codex(dry=False, remove=True)
    assert load_toml(cfg["codex"]) == {"mcp_servers": {"zzz": {"command": "zzz"}}}


def test_codex_validates_before_writing(cfg, capsys, monkeypatch):
    write(cfg["codex"], CODEX_BASE)
    monkeypatch.setattr(install, "_toml_str", lambda s: s)  # unquoted value -> invalid TOML block
    install.register_codex(dry=False)
    assert read(cfg["codex"]) == CODEX_BASE
    assert not os.path.exists(cfg["codex"] + install.BACKUP_SUFFIX)
    assert not os.path.exists(cfg["codex"] + ".tmp")
    out = capsys.readouterr().out
    assert "left unchanged" in out and "saved" not in out


def test_codex_does_not_touch_an_unparsable_config(cfg, capsys):
    broken = "this is = not toml [\n"
    write(cfg["codex"], broken)
    install.register_codex(dry=False)
    install.register_codex(dry=False, remove=True)
    assert read(cfg["codex"]) == broken
    assert not os.path.exists(cfg["codex"] + install.BACKUP_SUFFIX)
    assert "not touching it" in capsys.readouterr().out


def test_toml_str_keeps_non_bmp_characters_and_escapes_controls():
    for s in ["/Users/me/Projects \U0001F9EA/chemdraw-mcp/.venv/bin/chemdraw-mcp",
              r"C:\Users\me\.venv\Scripts\chemdraw-mcp.exe", 'a"b\\c\td\x01e\x7ff', "caf\u00e9/\u5316\u5b66"]:
        assert tomllib.loads("x = %s" % install._toml_str(s))["x"] == s
    assert "\\ud83e" not in install._toml_str("\U0001F9EA")


def test_backup_survives_idempotent_rerun_and_uninstall(cfg):
    original = {"mcpServers": {"chemdraw": FOREIGN}}
    write(cfg["claude"], json.dumps(original))
    install.register_claude(dry=False)
    bak = cfg["claude"] + install.BACKUP_SUFFIX
    assert json.loads(read(bak)) == original
    after_first = read(cfg["claude"])
    install.register_claude(dry=False)  # idempotent re-run must not roll the backup forward
    assert read(cfg["claude"]) == after_first
    assert json.loads(read(bak)) == original
    install.register_claude(dry=False, remove=True)
    assert json.loads(read(bak)) == original
    # codex: the foreign table replaced on the first run is still in the backup after a re-run
    codex_original = '[mcp_servers.chemdraw]\ncommand = "old"\nargs = []\n'
    write(cfg["codex"], codex_original)
    install.register_codex(dry=False)
    install.register_codex(dry=False)
    install.register_codex(dry=False, remove=True)
    assert read(cfg["codex"] + install.BACKUP_SUFFIX) == codex_original


def test_rerun_does_not_rewrite_unchanged_configs(cfg, capsys):
    write(cfg["desktop"], json.dumps({"mcpServers": {}}))
    assert install.main([]) == 0
    before = _snapshot(cfg["home"])
    mtimes = {p: os.stat(p).st_mtime_ns for p in before}
    assert install.main([]) == 0
    assert _snapshot(cfg["home"]) == before
    assert {p: os.stat(p).st_mtime_ns for p in before} == mtimes
    assert "already up to date" in capsys.readouterr().out


def test_config_paths_honour_client_env_overrides(tmp_path, monkeypatch):
    home = str(tmp_path / "home")
    monkeypatch.setenv("HOME", home)
    monkeypatch.setenv("USERPROFILE", home)
    monkeypatch.delenv("CLAUDE_CONFIG_DIR", raising=False)
    monkeypatch.delenv("CODEX_HOME", raising=False)
    assert install.claude_config_path() == os.path.join(home, ".claude.json")
    assert install.codex_config_path() == os.path.join(home, ".codex", "config.toml")
    monkeypatch.setenv("CLAUDE_CONFIG_DIR", str(tmp_path / "claude-cfg"))
    monkeypatch.setenv("CODEX_HOME", str(tmp_path / "codex-home"))
    assert install.claude_config_path() == os.path.join(str(tmp_path / "claude-cfg"), ".claude.json")
    assert install.codex_config_path() == os.path.join(str(tmp_path / "codex-home"), "config.toml")
    monkeypatch.setenv("CLAUDE_CONFIG_DIR", "")  # empty counts as unset
    monkeypatch.setenv("CODEX_HOME", "")
    assert install.claude_config_path() == os.path.join(home, ".claude.json")
    assert install.codex_config_path() == os.path.join(home, ".codex", "config.toml")


def test_check_prints_resolved_config_paths(cfg, capsys):
    assert install.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert cfg["claude"] in out and cfg["codex"] in out
