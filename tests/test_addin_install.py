"""Tests for the add-in / Scripts-menu / launchd installer helpers (dry runs and pure functions only)."""
import os
import plistlib

import pytest

from chemdrawmcp import addin_install


def test_applescript_source_quotes_the_bridge_command(monkeypatch):
    monkeypatch.setattr(addin_install, "bridge_command", lambda: ["/path with space/chemdraw-agent-bridge"])
    src = addin_install._applescript_source()
    assert "do shell script" in src and "nohup" in src
    assert "'/path with space/chemdraw-agent-bridge' --quiet" in src
    assert "display notification" in src and "47325" in src


def test_sh_and_applescript_quoting():
    assert addin_install._sh_quote("a'b") == "'a'\\''b'"
    assert addin_install._as_quote('say "hi"') == '"say \\"hi\\""'


def test_dry_runs_write_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("CHEMDRAWMCP_ADDIN_DIR", str(tmp_path / "addins"))
    addin_install.install_addin(dry=True)
    addin_install.install_scripts_menu(dry=True)
    addin_install.install_launch_agent(dry=True)
    addin_install.uninstall_addin(dry=True)
    addin_install.uninstall_scripts_menu(dry=True)
    addin_install.uninstall_launch_agent(dry=True)
    assert not (tmp_path / "addins").exists()
    assert not list(tmp_path.rglob("*.scpt")) and not list(tmp_path.rglob("*.plist"))
    out = capsys.readouterr().out
    assert "would" in out


@pytest.mark.skipif(os.sys.platform != "darwin", reason="launchd / Scripts menu are macOS features")
def test_launch_agent_plist_shape(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path / "home"))
    monkeypatch.setattr(addin_install, "_bootout", lambda p: None)
    calls = []
    monkeypatch.setattr(addin_install.subprocess, "run", lambda *a, **k: calls.append(a) or type("R", (), {"returncode": 0, "stderr": "", "stdout": ""})())
    p = addin_install.install_launch_agent(dry=False)
    assert p and os.path.exists(p)
    with open(p, "rb") as f:
        d = plistlib.load(f)
    assert d["Label"] == addin_install.LAUNCH_AGENT_LABEL and d["KeepAlive"] is True and d["RunAtLoad"] is True
    assert d["ProgramArguments"][-1] == "--quiet"
    assert calls and calls[0][0][:2] == ["launchctl", "bootstrap"]
    addin_install.uninstall_launch_agent(dry=False)
    assert not os.path.exists(p)


@pytest.mark.skipif(os.sys.platform != "darwin", reason="osacompile is macOS only")
def test_scripts_menu_install_compiles(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    target = addin_install.install_scripts_menu(dry=False)
    assert target and os.path.exists(target) and target.endswith(".scpt")
    addin_install.uninstall_scripts_menu(dry=False)
    assert not os.path.exists(target)
