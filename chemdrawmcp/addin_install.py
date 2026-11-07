"""Install the "Agent Chat" ChemDraw add-in, the Scripts-menu starter and (optionally) a launchd agent for the bridge.

* The add-in is an HTML/JS panel (``chemdrawmcp/addin``) that ChemDraw shows from its Add-ins menu.  It talks to the
  local agent bridge (``chemdraw-agent-bridge``) which runs Claude Code / Codex.  The installer builds the add-in with
  the bridge port and the persistent bridge token baked into ``config.js`` and copies it into ChemDraw's user add-in
  folder (``~/Library/Application Support/com.revvity.ChemDraw/Add-ins`` on macOS); ChemDraw lists it after a restart.
* ``Start Agent Bridge`` in ChemDraw's *Scripts* menu (macOS: an AppleScript in
  ``~/Library/Application Support/com.revvity.ChemDraw/Scripts``) starts the bridge from inside ChemDraw.
* ``--launch-agent`` installs ``~/Library/LaunchAgents/com.chemdrawmcp.bridge.plist`` so the bridge is always running.
"""
from __future__ import annotations

import os
import plistlib
import shutil
import subprocess
import sys
from typing import Optional

ADDIN_NAME = "Agent Chat"
SCRIPT_NAME = "Start Agent Bridge"
LAUNCH_AGENT_LABEL = "com.chemdrawmcp.bridge"


def log(msg: str):
    print(msg, flush=True)


def bridge_executable() -> str:
    exe = os.path.join(os.path.dirname(os.path.abspath(sys.executable)), "chemdraw-agent-bridge" + (".exe" if sys.platform.startswith("win") else ""))
    return exe if os.path.exists(exe) else ""


def bridge_command() -> list:
    exe = bridge_executable()
    if exe:
        return [exe]
    return [sys.executable, "-m", "chemdrawmcp.bridge.server"]


def scripts_dir() -> Optional[str]:
    if sys.platform == "darwin":
        return os.path.expanduser("~/Library/Application Support/com.revvity.ChemDraw/Scripts")
    return None


def launch_agent_path() -> Optional[str]:
    if sys.platform == "darwin":
        return os.path.expanduser("~/Library/LaunchAgents/%s.plist" % LAUNCH_AGENT_LABEL)
    return None


# ---------------------------------------------------------------------------
# add-in
# ---------------------------------------------------------------------------

def install_addin(dry: bool, port: Optional[int] = None) -> Optional[str]:
    from . import addin_builder
    from .bridge import server as bridge_server
    log("\n== ChemDraw add-in '%s' ==" % ADDIN_NAME)
    dest = addin_builder.addin_install_dir()
    if not dest:
        log("no ChemDraw user add-in folder on this platform; skipping")
        return None
    port = port or bridge_server.DEFAULT_PORT
    if dry:
        log("would build the add-in with port %d and the bridge token into %s" % (port, os.path.join(dest, ADDIN_NAME)))
        return None
    token = bridge_server.ensure_token()
    path = addin_builder.build_addin(dest, port, token, name=ADDIN_NAME)
    pkg = addin_builder.package_addin(os.path.join(_home(), "build"), ADDIN_NAME, source_dir=path) if _has_kw(addin_builder.package_addin, "source_dir") else None
    log("installed %s" % path)
    if pkg:
        log("packaged %s (for ChemDraw's Add-in Manager on other machines)" % pkg)
    log("restart ChemDraw, then open  Add-ins > %s  (it needs the bridge: Scripts > %s, or run chemdraw-agent-bridge)" % (ADDIN_NAME, SCRIPT_NAME))
    return path


def uninstall_addin(dry: bool):
    from . import addin_builder
    log("\n== ChemDraw add-in '%s' (remove) ==" % ADDIN_NAME)
    path = addin_builder.installed_addin_dir(ADDIN_NAME)
    if path and os.path.isdir(path):
        if dry:
            log("would remove %s" % path)
        else:
            shutil.rmtree(path, ignore_errors=True)
            log("removed %s" % path)
    else:
        log("not installed")


# ---------------------------------------------------------------------------
# Scripts menu (macOS)
# ---------------------------------------------------------------------------

def _applescript_source() -> str:
    cmd = " ".join(_sh_quote(c) for c in bridge_command()) + " --quiet"
    return (
        '-- Starts the chemdrawmcp agent bridge (Claude Code / Codex runner for the Agent Chat add-in).\n'
        'set bridgeCommand to %s\n'
        'do shell script "nohup " & bridgeCommand & " >/dev/null 2>&1 &"\n'
        'delay 1\n'
        'display notification "Agent bridge started (http://127.0.0.1:%d)" with title "ChemDraw Agent Chat"\n'
    ) % (_as_quote(cmd), _port())


def _port() -> int:
    from .bridge import server as bridge_server
    return bridge_server.DEFAULT_PORT


def _sh_quote(s: str) -> str:
    return "'" + s.replace("'", "'\\''") + "'"


def _as_quote(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def install_scripts_menu(dry: bool) -> Optional[str]:
    d = scripts_dir()
    log("\n== ChemDraw Scripts menu: '%s' ==" % SCRIPT_NAME)
    if not d:
        log("not available on this platform; start the bridge with: chemdraw-agent-bridge")
        return None
    target = os.path.join(d, SCRIPT_NAME + ".scpt")
    if dry:
        log("would compile an AppleScript into %s" % target)
        return None
    os.makedirs(d, exist_ok=True)
    wanted = _applescript_source()
    if os.path.exists(target) and _decompiled(target) == _normalise(wanted):
        log("up to date: %s" % target)
        return target
    src = os.path.join(d, SCRIPT_NAME + ".applescript")
    with open(src, "w", encoding="utf-8") as f:
        f.write(wanted)
    try:
        subprocess.run(["osacompile", "-o", target, src], check=True, capture_output=True, text=True)
        os.remove(src)
        log("installed %s (ChemDraw shows it in its Scripts menu)" % target)
        return target
    except Exception as e:
        log("osacompile failed (%s); left the plain-text script at %s" % (e, src))
        return src


def _normalise(text: str) -> str:
    return "\n".join(line.strip() for line in text.strip().splitlines() if line.strip() and not line.strip().startswith("--"))


def _decompiled(target: str) -> str:
    try:
        r = subprocess.run(["osadecompile", target], capture_output=True, text=True, timeout=30)
        return _normalise(r.stdout) if r.returncode == 0 else ""
    except Exception:
        return ""


def uninstall_scripts_menu(dry: bool):
    d = scripts_dir()
    if not d:
        return
    log("\n== ChemDraw Scripts menu: '%s' (remove) ==" % SCRIPT_NAME)
    for ext in (".scpt", ".applescript"):
        p = os.path.join(d, SCRIPT_NAME + ext)
        if os.path.exists(p):
            if dry:
                log("would remove %s" % p)
            else:
                os.remove(p)
                log("removed %s" % p)


# ---------------------------------------------------------------------------
# launchd agent (macOS)
# ---------------------------------------------------------------------------

def install_launch_agent(dry: bool) -> Optional[str]:
    p = launch_agent_path()
    log("\n== launchd agent for the bridge ==")
    if not p:
        log("launchd agents exist only on macOS; run chemdraw-agent-bridge yourself (or add it to your startup items)")
        return None
    plist = {
        "Label": LAUNCH_AGENT_LABEL,
        "ProgramArguments": bridge_command() + ["--quiet"],
        "RunAtLoad": True,
        "KeepAlive": True,
        "StandardOutPath": os.path.join(_home(), "bridge.log"),
        "StandardErrorPath": os.path.join(_home(), "bridge.log"),
        "EnvironmentVariables": {k: v for k, v in os.environ.items() if k.startswith("CHEMDRAWMCP_") or k in ("PATH", "HOME")},
    }
    if dry:
        log("would write %s and load it with launchctl" % p)
        return None
    os.makedirs(os.path.dirname(p), exist_ok=True)
    _bootout(p)
    with open(p, "wb") as f:
        plistlib.dump(plist, f)
    r = subprocess.run(["launchctl", "bootstrap", "gui/%d" % os.getuid(), p], capture_output=True, text=True)
    if r.returncode == 0:
        log("installed and started %s" % p)
    else:
        log("wrote %s but launchctl bootstrap failed: %s" % (p, (r.stderr or r.stdout).strip()))
    return p


def uninstall_launch_agent(dry: bool):
    p = launch_agent_path()
    if not p or not os.path.exists(p):
        return
    log("\n== launchd agent (remove) ==")
    if dry:
        log("would unload and remove %s" % p)
        return
    _bootout(p)
    os.remove(p)
    log("removed %s" % p)


def _bootout(p: str):
    try:
        subprocess.run(["launchctl", "bootout", "gui/%d" % os.getuid(), p], capture_output=True, text=True, timeout=30)
    except Exception:
        pass


def _home() -> str:
    from . import discovery
    return discovery.home_dir()


def _has_kw(fn, name: str) -> bool:
    import inspect
    try:
        return name in inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False
