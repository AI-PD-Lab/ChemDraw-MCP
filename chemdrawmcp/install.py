"""chemdrawmcp installer: register the ``chemdraw`` MCP server with Claude Code, Codex and Claude Desktop.

Run through the repository's ``install.py`` (which creates the virtualenv) or,
inside the virtualenv, ``chemdrawmcp-install``.  Works on macOS, Linux and Windows.
Nothing has to be installed inside ChemDraw itself: the server drives the app
through its automation interface (AppleScript on macOS, COM on Windows) and
falls back to the headless CDXML/RDKit engine when ChemDraw is not available.

Steps (all idempotent; choose with flags, default = everything that applies):
  --claude           register the server with Claude Code (user scope, ~/.claude.json,
                     or $CLAUDE_CONFIG_DIR/.claude.json when that variable is set)
  --codex            register with Codex CLI (~/.codex/config.toml, or $CODEX_HOME/config.toml)
  --claude-desktop   register with Claude Desktop (if its config file exists)
  --addin            install the 'Agent Chat' ChemDraw add-in + Scripts-menu bridge starter (default when ChemDraw is found)
  --launch-agent     keep the agent bridge running via launchd (macOS, opt-in)
  --uninstall        remove the registrations (and the add-in / launch agent)
  --dry-run          show what would be done
  --check            only report what was found
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
import re
import shutil
import sys
from typing import List, Optional, Tuple

from chemdrawmcp import __version__
from chemdrawmcp import discovery

SERVER_NAME = "chemdraw"
ALT_SERVER_NAME = "chemdraw-mcp"   # used when the user already has a foreign "chemdraw" server
SCRIPT_NAME = "chemdraw-mcp"
BACKUP_SUFFIX = ".chemdrawmcp-backup"


def log(msg: str):
    print(msg, flush=True)


def venv_bin_dir() -> str:
    return os.path.dirname(os.path.abspath(sys.executable))


def server_command() -> str:
    name = SCRIPT_NAME + (".exe" if sys.platform.startswith("win") else "")
    path = os.path.join(venv_bin_dir(), name)
    if not os.path.exists(path):
        # fall back to running the module with this interpreter
        return ""
    return path


def server_entry() -> Tuple[str, List[str]]:
    cmd = server_command()
    if cmd:
        return cmd, []
    return sys.executable, ["-m", "chemdrawmcp.server"]


def _is_ours(entry: dict) -> bool:
    cmd = ((entry or {}).get("command", "") or "").replace("\\", "/")
    base = cmd.rsplit("/", 1)[-1]
    args = " ".join(str(a) for a in ((entry or {}).get("args", []) or []))
    return base in ("chemdraw-mcp", "chemdraw-mcp.exe") or "chemdrawmcp." in args


# ---------------------------------------------------------------------------
# MCP client registration
# ---------------------------------------------------------------------------

def _backup(path: str, force: bool = False):
    """Copy ``path`` to ``path + BACKUP_SUFFIX``.

    An existing backup is kept (it holds the config as it was before chemdrawmcp first touched
    it) unless ``force`` is set, which callers use when they are about to replace something of
    the user's that they promise is kept in the backup.
    """
    if os.path.exists(path) and (force or not os.path.exists(path + BACKUP_SUFFIX)):
        shutil.copy(path, path + BACKUP_SUFFIX)


def _read_text(path: str) -> Optional[str]:
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write_text(path: str, text: str):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
    os.replace(tmp, path)


def _save(path: str, original: Optional[str], text: str, force_backup: bool = False) -> bool:
    """Back up ``path`` and write ``text`` to it, unless the file already holds exactly ``text``.

    ``original`` is the file's current contents (None when it does not exist yet).  Returns
    whether anything was written.
    """
    if original is not None and text == original:
        return False
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    _backup(path, force=force_backup)
    _write_text(path, text)
    return True


def _json_text(data: dict) -> str:
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def claude_config_path() -> str:
    # Claude Code keeps .claude.json inside $CLAUDE_CONFIG_DIR when that variable is set
    base = os.environ.get("CLAUDE_CONFIG_DIR") or os.path.expanduser("~")
    return os.path.join(base, ".claude.json")


def register_claude(dry: bool, remove: bool = False):
    path = claude_config_path()
    log("\n== Claude Code (user scope, %s) ==" % path)
    original = _read_text(path)
    data = {}
    if original is not None:
        try:
            data = json.loads(original)
        except Exception as e:
            log("cannot parse %s (%s); not touching it" % (path, e))
            return
    elif remove:
        log("no config file; nothing to remove")
        return
    servers = data.setdefault("mcpServers", {})
    replacing = False
    if remove:
        if SERVER_NAME in servers and _is_ours(servers[SERVER_NAME]):
            servers.pop(SERVER_NAME)
            log("removed %s" % SERVER_NAME)
        else:
            log("%s is not registered (or not ours); nothing to remove" % SERVER_NAME)
            return
    else:
        if SERVER_NAME in servers and not _is_ours(servers[SERVER_NAME]):
            replacing = True
            log("replacing an existing non-chemdrawmcp '%s' entry (kept in %s%s)" % (SERVER_NAME, path, BACKUP_SUFFIX))
        cmd, args = server_entry()
        servers[SERVER_NAME] = {"type": "stdio", "command": cmd, "args": args, "env": {}}
        log("%s -> %s %s" % (SERVER_NAME, cmd, " ".join(args)))
    if dry:
        return
    if _save(path, original, _json_text(data), force_backup=replacing):
        log("saved. In Claude Code run /mcp to see the server (restart running sessions).")
    else:
        log("already up to date; nothing written.")


def codex_config_path() -> str:
    # Codex keeps its config under $CODEX_HOME when that variable is set
    base = os.environ.get("CODEX_HOME") or os.path.join(os.path.expanduser("~"), ".codex")
    return os.path.join(base, "config.toml")


def _toml_str(s: str) -> str:
    """``s`` as a TOML basic string: raw UTF-8 with ``"``, ``\\`` and control characters escaped.

    ``ensure_ascii`` must be off: json would write non-BMP characters (emoji in a path) as a
    ``\\uXXXX`` surrogate pair, which TOML rejects.  json leaves U+007F alone but TOML wants it escaped.
    """
    return json.dumps(s, ensure_ascii=False).replace("\x7f", "\\u007F")


def _toml_key_pattern(name: str) -> str:
    """Regex for the TOML key ``name``, bare or quoted."""
    n = re.escape(name)
    return "(?:%s|\"%s\"|'%s')" % (n, n, n)


def _codex_block_pattern(name: str) -> "re.Pattern[str]":
    # the [mcp_servers.<name>] table (the key may be quoted) and any of its sub-tables
    # ([mcp_servers.<name>.env] ...), each running up to the next table header or the end of the file
    return re.compile(r"(?ms)^\[mcp_servers\.%s(?:\.[^\]\n]*)?\][^\n]*(?:\n|\Z).*?(?=^\[|\Z)" % _toml_key_pattern(name))


def _set_codex_block(text: str, name: str, block: Optional[str]) -> str:
    """Replace (or remove, when ``block`` is None) the ``[mcp_servers.<name>]`` table in ``text``."""
    pat = _codex_block_pattern(name)
    m = pat.search(text)
    if not m:
        if block is None:
            return text
        if text and not text.endswith("\n"):
            text += "\n"
        if text and not text.endswith("\n\n"):
            text += "\n"
        return text + block
    start = m.start()
    text = pat.sub("", text)
    if block is not None:
        text = text[:start] + block + text[start:]
    return text


def _toml_module():
    """The stdlib ``tomllib`` (python 3.11+) or ``tomli``; None when neither is importable."""
    for mod in ("tomllib", "tomli"):
        try:
            return importlib.import_module(mod)
        except ImportError:
            continue
    return None


def _validate_toml(text: str) -> Optional[str]:
    """Error text if ``text`` does not parse as TOML (None when fine or no parser is available)."""
    toml = _toml_module()
    if toml is None:
        return None
    try:
        toml.loads(text)
    except Exception as e:
        return str(e)
    return None


def _codex_entry(text: str, name: str):
    """The value of ``mcp_servers.<name>`` in ``text``, however it is spelled: a ``[mcp_servers.<name>]``
    table, an inline table under ``[mcp_servers]``, a dotted key ...  None when it is not defined.

    Without a TOML parser (python 3.10 without tomli) the common spellings are detected textually
    and ``{}`` stands in for the value.
    """
    toml = _toml_module()
    if toml is not None:
        try:
            servers = toml.loads(text).get("mcp_servers")
        except Exception:
            return None
        return servers.get(name) if isinstance(servers, dict) else None
    key = _toml_key_pattern(name)
    if re.search(r"(?m)^\s*mcp_servers\s*\.\s*%s\s*[.=]" % key, text):        # mcp_servers.<name>.command = ...
        return {}
    table = re.search(r"(?ms)^\s*\[\s*mcp_servers\s*\][^\n]*\n(.*?)(?=^\s*\[|\Z)", text)  # [mcp_servers] <name> = { ... }
    if table and re.search(r"(?m)^\s*%s\s*[.=]" % key, table.group(1)):
        return {}
    return {} if _codex_block_pattern(name).search(text) else None


def register_codex(dry: bool, remove: bool = False):
    path = codex_config_path()
    log("\n== Codex CLI (%s) ==" % path)
    original = _read_text(path)
    if original is None and remove:
        log("no config file; nothing to remove")
        return
    text = original or ""
    err = _validate_toml(text)
    if err:
        log("cannot parse %s (%s); not touching it" % (path, err))
        return
    has_block = _codex_block_pattern(SERVER_NAME).search(text) is not None
    replacing = False
    if remove:
        if not has_block:
            log("%s is not registered; nothing to remove" % SERVER_NAME)
            return
        text = _set_codex_block(text, SERVER_NAME, None)
        log("removed %s" % SERVER_NAME)
    else:
        existing = _codex_entry(text, SERVER_NAME)
        if existing is not None and not has_block:
            # e.g. `[mcp_servers]` + `chemdraw = { ... }` or `mcp_servers.chemdraw.command = ...`: appending
            # our table would declare mcp_servers.chemdraw twice and Codex would reject the whole file
            log("'%s' is already defined in %s outside a [mcp_servers.%s] table (inline table or dotted key); "
                "remove that definition and run the installer again" % (SERVER_NAME, path, SERVER_NAME))
            return
        if isinstance(existing, dict) and existing and not _is_ours(existing):
            replacing = True
            log("replacing an existing non-chemdrawmcp '%s' table (kept in %s%s)" % (SERVER_NAME, path, BACKUP_SUFFIX))
        cmd, args = server_entry()
        block = "[mcp_servers.%s]\ncommand = %s\nargs = [%s]\n\n" % (
            SERVER_NAME, _toml_str(cmd), ", ".join(_toml_str(a) for a in args))
        text = _set_codex_block(text, SERVER_NAME, block)
        log("%s -> %s %s" % (SERVER_NAME, cmd, " ".join(args)))
    text = text.rstrip("\n") + "\n" if text.strip() else ""
    err = _validate_toml(text)   # before anything is written: a broken config.toml makes Codex reject it whole
    if err:
        log("ERROR: the edited config would not parse as TOML (%s); %s left unchanged" % (err, path))
        return
    if dry:
        return
    if _save(path, original, text, force_backup=replacing):
        log("saved. In Codex run `codex mcp list` to verify.")
    else:
        log("already up to date; nothing written.")


def claude_desktop_config_path() -> Optional[str]:
    if sys.platform == "darwin":
        p = os.path.expanduser("~/Library/Application Support/Claude/claude_desktop_config.json")
    elif sys.platform.startswith("win"):
        p = os.path.join(os.environ.get("APPDATA", ""), "Claude", "claude_desktop_config.json")
    else:
        p = os.path.expanduser("~/.config/Claude/claude_desktop_config.json")
    return p if os.path.exists(p) else None


def register_claude_desktop(dry: bool, remove: bool = False):
    path = claude_desktop_config_path()
    if not path:
        return
    log("\n== Claude Desktop (%s) ==" % path)
    original = _read_text(path) or ""
    try:
        data = json.loads(original)
    except Exception as e:
        log("cannot parse config (%s); skipping" % e)
        return
    servers = data.setdefault("mcpServers", {})
    if remove:
        removed = False
        for n in (SERVER_NAME, ALT_SERVER_NAME):
            if n in servers and _is_ours(servers[n]):
                servers.pop(n)
                log("removed %s" % n)
                removed = True
        if not removed:
            log("nothing of ours registered; nothing to remove")
            return
    else:
        name = SERVER_NAME
        if name in servers and not _is_ours(servers[name]):
            name = ALT_SERVER_NAME   # keep the user's existing server untouched
            log("'%s' is already used by another server; registering as '%s'" % (SERVER_NAME, name))
        cmd, args = server_entry()
        servers[name] = {"command": cmd, "args": args}
        log("%s -> %s %s" % (name, cmd, " ".join(args)))
    if dry:
        return
    if _save(path, original, _json_text(data)):
        log("saved (restart Claude Desktop).")
    else:
        log("already up to date; nothing written.")


# ---------------------------------------------------------------------------

def _rdkit_version() -> str:
    try:
        rdkit = importlib.import_module("rdkit")
    except Exception:
        return "NOT INSTALLED (headless engine, naming and properties unavailable)"
    return getattr(rdkit, "__version__", "installed")


def check():
    log("== chemdrawmcp %s ==" % __version__)
    log("python:    %s" % sys.executable)
    log("platform:  %s" % discovery.platform_note())
    app = discovery.find_chemdraw()
    if app:
        details = []
        ver = discovery.chemdraw_version(app)
        if ver:
            details.append("version %s" % ver)
        if sys.platform == "darwin":
            bid = discovery.macos_bundle_id(app)
            if bid:
                details.append(bid)
        log("chemdraw:  %s%s" % (app, " (%s)" % ", ".join(details) if details else ""))
    else:
        log("chemdraw:  NOT FOUND (set CHEMDRAW_APP / CHEMDRAW_EXE)")
    log("rdkit:     %s" % _rdkit_version())
    cmd, args = server_entry()
    log("server:    %s %s" % (cmd, " ".join(args)))
    log("home:      %s" % discovery.home_dir())
    log("renders:   %s" % discovery.output_dir())
    log("configs:   Claude Code    %s" % claude_config_path())
    log("           Codex CLI      %s" % codex_config_path())
    log("           Claude Desktop %s" % (claude_desktop_config_path() or "(no config file found)"))
    try:
        from . import addin_builder, addin_install
        from .bridge import server as bridge_server
        inst = addin_builder.installed_addin_dir(addin_install.ADDIN_NAME)
        log("add-in:    %s" % (inst if inst and os.path.isdir(inst) else "not installed (python3 install.py --addin)"))
        info = bridge_server.bridge_info()
        running = bridge_server.is_bridge_running(info) if info else False
        log("bridge:    %s" % ((info.get("url") or "http://127.0.0.1:%s" % info.get("port")) + " (running)" if running
                               else "not running (Scripts > Start Agent Bridge, or chemdraw-agent-bridge)"))
    except Exception as e:  # pragma: no cover - optional pieces
        log("add-in:    (unavailable: %s)" % e)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="chemdrawmcp-install", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--claude", action="store_true")
    ap.add_argument("--codex", action="store_true")
    ap.add_argument("--claude-desktop", action="store_true")
    ap.add_argument("--addin", action="store_true", help="install the 'Agent Chat' ChemDraw add-in and the Scripts-menu bridge starter")
    ap.add_argument("--launch-agent", action="store_true", help="also install a launchd agent that keeps the bridge running (macOS)")
    ap.add_argument("--uninstall", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true", help="only report what was found")
    a = ap.parse_args(argv)
    check()
    if a.check:
        return 0
    chosen = any([a.claude, a.codex, a.claude_desktop, a.addin, a.launch_agent])
    addin_possible = bool(discovery.find_chemdraw()) and sys.platform in ("darwin", "win32")
    do = lambda flag: (flag or not chosen)  # noqa: E731
    dry = a.dry_run
    if dry:
        log("\n(dry run: nothing will be written)")
    if a.uninstall:
        if do(a.claude):
            register_claude(dry, remove=True)
        if do(a.codex):
            register_codex(dry, remove=True)
        if do(a.claude_desktop):
            register_claude_desktop(dry, remove=True)
        if a.addin or (not chosen and addin_possible):
            from . import addin_install
            addin_install.uninstall_addin(dry)
            addin_install.uninstall_scripts_menu(dry)
        if a.launch_agent or not chosen:
            from . import addin_install
            addin_install.uninstall_launch_agent(dry)
        log("\nDone. The chemdraw server registrations were removed%s." % (" (dry run)" if dry else ""))
        return 0
    if do(a.claude):
        register_claude(dry)
    if do(a.codex):
        register_codex(dry)
    if do(a.claude_desktop):
        register_claude_desktop(dry)
    if a.addin or (not chosen and addin_possible):
        from . import addin_install
        addin_install.install_addin(dry)
        addin_install.install_scripts_menu(dry)
    if a.launch_agent:
        from . import addin_install
        addin_install.install_launch_agent(dry)
    log("\nDone%s. Verify with:" % (" (dry run)" if dry else ""))
    if do(a.claude):
        log("  Claude Code:    /mcp  (restart running sessions)")
    if do(a.codex):
        log("  Codex CLI:      codex mcp list")
    if do(a.claude_desktop) and claude_desktop_config_path():
        log("  Claude Desktop: restart the app, then check Settings > Developer")
    if a.addin or (not chosen and addin_possible):
        log("  ChemDraw:       restart it, then Add-ins > Agent Chat (start the bridge via Scripts > Start Agent Bridge)")
    log("Then ask the agent for the chemdraw `status` tool%s." % (
        "" if discovery.find_chemdraw() or not sys.platform.startswith("linux") else " (headless mode on Linux)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
