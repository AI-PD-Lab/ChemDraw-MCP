"""Build, package and locate the "Agent Chat" ChemDraw add-in (``chemdrawmcp/addin/``).

A ChemDraw add-in is a folder with ``chemdraw-addin-metadata.json`` plus HTML/JS/CSS that ChemDraw
(macOS and Windows) shows in a web view when the user picks it from the Add-ins menu.  Ours talks to
the local agent bridge (:mod:`chemdrawmcp.bridge.server`); the bridge port and token reach the page
through ``config.js``, which this module writes next to ``main.html``::

    window.CHEMDRAWMCP_BRIDGE = {port: 47325, token: "..."};

Usage (also wired into ``python3 install.py --addin``)::

    python -m chemdrawmcp.addin_builder                     # install into the user's Add-ins folder
    python -m chemdrawmcp.addin_builder --package dist/     # also zip it as "Agent Chat.chemdrawaddin"
    python -m chemdrawmcp.addin_builder --dest DIR --port 47325 --token SECRET

Port/token default to the running bridge's ``bridge.json`` (``<CHEMDRAWMCP_HOME>/bridge.json``),
then ``$CHEMDRAWMCP_BRIDGE_PORT`` / ``$CHEMDRAWMCP_BRIDGE_TOKEN``.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import zipfile
from typing import List, Optional

ADDIN_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "addin")
DEFAULT_NAME = "Agent Chat"
DEFAULT_PORT = 47325
METADATA_FILE = "chemdraw-addin-metadata.json"
CONFIG_FILE = "config.js"
PACKAGE_FILES = [METADATA_FILE, "main.html", "app.js", "style.css", "config.example.js", "README.md"]
PACKAGE_SUFFIX = ".chemdrawaddin"

MAC_ADDIN_DIR = os.path.join("~", "Library", "Application Support", "com.revvity.ChemDraw", "Add-ins")
WIN_ADDIN_DIR = os.path.join("RevvitySignalsSoftware", "ChemDrawApplications", "ChemDraw", "ChemDraw Items", "Add-ins")


# ---------------------------------------------------------------------------
# locations
# ---------------------------------------------------------------------------

def addin_install_dir() -> str:
    """The user's ChemDraw add-in folder for this platform (``$CHEMDRAWMCP_ADDIN_DIR`` overrides).

    macOS:   ~/Library/Application Support/com.revvity.ChemDraw/Add-ins
    Windows: %ProgramData%\\RevvitySignalsSoftware\\ChemDrawApplications\\ChemDraw\\ChemDraw Items\\Add-ins
    Linux (no ChemDraw): ~/.chemdrawmcp/addins, so packaging still works.
    """
    override = os.environ.get("CHEMDRAWMCP_ADDIN_DIR")
    if override:
        return os.path.abspath(os.path.expanduser(override))
    if sys.platform == "darwin":
        return os.path.expanduser(MAC_ADDIN_DIR)
    if sys.platform.startswith("win"):
        program_data = os.environ.get("PROGRAMDATA") or r"C:\ProgramData"
        return os.path.join(program_data, WIN_ADDIN_DIR)
    home = os.environ.get("CHEMDRAWMCP_HOME") or os.path.join(os.path.expanduser("~"), ".chemdrawmcp")
    return os.path.join(os.path.abspath(os.path.expanduser(home)), "addins")


def installed_addin_dir(name: str = DEFAULT_NAME) -> str:
    """Where ``build_addin`` puts (or would put) the add-in folder for ``name``."""
    return os.path.join(addin_install_dir(), name)


def is_installed(name: str = DEFAULT_NAME) -> bool:
    d = installed_addin_dir(name)
    return os.path.isfile(os.path.join(d, METADATA_FILE)) and os.path.isfile(os.path.join(d, "main.html"))


# ---------------------------------------------------------------------------
# building
# ---------------------------------------------------------------------------

def metadata() -> dict:
    with open(os.path.join(ADDIN_SRC, METADATA_FILE), encoding="utf-8") as f:
        return json.load(f)


def config_js(port: int, token: str) -> str:
    """The ``config.js`` text read by app.js (``window.CHEMDRAWMCP_BRIDGE``)."""
    port = int(port)
    if port <= 0 or port > 65535:
        raise ValueError("port must be between 1 and 65535, got %r" % (port,))
    body = json.dumps({"port": port, "token": str(token or "")}, indent=2)
    return ("// Written by chemdrawmcp (python3 install.py --addin). Do not edit by hand: rerun the installer\n"
            "// after `chemdraw-agent-bridge --new-token` or when the bridge port changes.\n"
            "window.CHEMDRAWMCP_BRIDGE = %s;\n" % body)


def build_addin(dest_dir: str, port: int, token: str, name: str = DEFAULT_NAME) -> str:
    """Copy the add-in package into ``dest_dir/<name>/`` and write ``config.js``; return that folder."""
    if not name or any(c in name for c in '/\\:*?"<>|') or name.strip() != name:
        raise ValueError("invalid add-in name %r" % (name,))
    target = os.path.join(os.path.abspath(os.path.expanduser(dest_dir)), name)
    os.makedirs(target, exist_ok=True)
    for fname in PACKAGE_FILES:
        src = os.path.join(ADDIN_SRC, fname)
        if not os.path.isfile(src):
            raise FileNotFoundError("add-in file missing from the package: %s" % src)
        shutil.copyfile(src, os.path.join(target, fname))
    meta = metadata()
    meta["name"] = name
    meta["menuItemText"] = name
    with open(os.path.join(target, METADATA_FILE), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
        f.write("\n")
    cfg_path = os.path.join(target, CONFIG_FILE)
    with open(cfg_path, "w", encoding="utf-8") as f:
        f.write(config_js(port, token))
    try:
        os.chmod(cfg_path, 0o600)   # the token is a secret
    except OSError:
        pass
    return target


def package_addin(dest_dir: str, name: str = DEFAULT_NAME) -> str:
    """Zip ``dest_dir/<name>/`` (built by :func:`build_addin`) as ``dest_dir/<name>.chemdrawaddin``.

    The files sit at the top level of the archive (``main.html``, not ``<name>/main.html``), which is
    what ChemDraw's Add-in Manager expects of a ``.chemdrawaddin`` package.
    """
    dest_dir = os.path.abspath(os.path.expanduser(dest_dir))
    folder = os.path.join(dest_dir, name)
    if not os.path.isfile(os.path.join(folder, METADATA_FILE)):
        raise FileNotFoundError("no built add-in in %s (run build_addin first)" % folder)
    out = os.path.join(dest_dir, name + PACKAGE_SUFFIX)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
        for root, dirs, files in os.walk(folder):
            dirs.sort()
            for fname in sorted(files):
                if fname.startswith(".") or fname.endswith(".pyc"):
                    continue
                full = os.path.join(root, fname)
                zf.write(full, os.path.relpath(full, folder).replace(os.sep, "/"))
    return out


def remove_addin(dest_dir: Optional[str] = None, name: str = DEFAULT_NAME) -> bool:
    """Delete ``dest_dir/<name>/`` (default: the installed add-in); True when something was removed."""
    folder = os.path.join(os.path.abspath(os.path.expanduser(dest_dir)), name) if dest_dir else installed_addin_dir(name)
    if not os.path.isdir(folder):
        return False
    if not os.path.isfile(os.path.join(folder, METADATA_FILE)):
        raise RuntimeError("refusing to delete %s: it is not a ChemDraw add-in folder" % folder)
    shutil.rmtree(folder)
    return True


# ---------------------------------------------------------------------------
# bridge connection details
# ---------------------------------------------------------------------------

def bridge_connection(port: Optional[int] = None, token: Optional[str] = None) -> dict:
    """Resolve port/token: explicit arguments, then the bridge's ``bridge.json``, then the environment."""
    info = None
    if port is None or not token:
        try:
            from chemdrawmcp.bridge import server as bridge_server  # local import: the bridge is optional
            info = bridge_server.bridge_info()
        except Exception:
            info = None
    if port is None:
        if info and isinstance(info.get("port"), int):
            port = info["port"]
        else:
            try:
                port = int(os.environ.get("CHEMDRAWMCP_BRIDGE_PORT") or DEFAULT_PORT)
            except ValueError:
                port = DEFAULT_PORT
    if not token:
        token = (info or {}).get("token") or os.environ.get("CHEMDRAWMCP_BRIDGE_TOKEN") or ""
    return {"port": int(port), "token": str(token), "source": "argument" if info is None and token else ("bridge.json" if info else "environment")}


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m chemdrawmcp.addin_builder", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", default=None, help="folder that receives <name>/ (default: the ChemDraw Add-ins folder)")
    ap.add_argument("--name", default=DEFAULT_NAME, help="add-in / menu item name (default: %(default)s)")
    ap.add_argument("--port", type=int, default=None, help="bridge port (default: bridge.json, $CHEMDRAWMCP_BRIDGE_PORT, 47325)")
    ap.add_argument("--token", default=None, help="bridge token (default: bridge.json, $CHEMDRAWMCP_BRIDGE_TOKEN)")
    ap.add_argument("--package", metavar="DIR", default=None, help="also write <name>.chemdrawaddin into DIR")
    ap.add_argument("--uninstall", action="store_true", help="remove the installed add-in folder")
    ap.add_argument("--check", action="store_true", help="only report locations")
    a = ap.parse_args(argv)
    dest = a.dest or addin_install_dir()
    if a.check:
        print("add-in folder: %s%s" % (installed_addin_dir(a.name), " (installed)" if is_installed(a.name) else " (not installed)"))
        conn = bridge_connection(a.port, a.token)
        print("bridge:        127.0.0.1:%d, token %s (%s)" % (conn["port"], "set" if conn["token"] else "MISSING", conn["source"]))
        return 0
    if a.uninstall:
        removed = remove_addin(a.dest, a.name)
        print("%s %s" % ("Removed" if removed else "Not installed:", os.path.join(dest, a.name)))
        return 0
    conn = bridge_connection(a.port, a.token)
    if not conn["token"]:
        print("warning: no bridge token found (start `chemdraw-agent-bridge` once, or pass --token); "
              "the add-in will report 'bridge not configured' until you rerun this with a token.", file=sys.stderr)
    folder = build_addin(dest, conn["port"], conn["token"], a.name)
    print("Installed add-in: %s (bridge 127.0.0.1:%d, token from %s)" % (folder, conn["port"], conn["source"]))
    if a.package:
        shutil.copytree(folder, os.path.join(a.package, a.name), dirs_exist_ok=True)
        print("Package:          %s" % package_addin(a.package, a.name))
    print("In ChemDraw: Add-ins menu > %s (restart ChemDraw if the menu does not list it yet)." % a.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
