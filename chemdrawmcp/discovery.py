"""Finding the ChemDraw application and the chemdrawmcp directories on any platform.

Order of precedence for the ChemDraw application (:func:`find_chemdraw`):

1. ``CHEMDRAW_APP`` / ``CHEMDRAW_EXE`` environment variables (a ``.app`` bundle on
   macOS, ``ChemDraw.exe`` on Windows).
2. macOS: the highest-versioned ``ChemDraw*.app`` in ``/Applications`` or
   ``~/Applications`` (``ChemDraw 23.1.2.app``, ``ChemDraw.app``,
   ``ChemDraw Professional*.app``; companions such as ``ChemDraw Collections.app``
   are ignored), then Spotlight (``mdfind``) by bundle identifier.
   Windows: the registry (Revvity Signals / PerkinElmer / CambridgeSoft / App Paths,
   including ``WOW6432Node``), then the usual ``Program Files`` locations; the newest
   executable wins.
   Linux: nothing - ChemDraw does not run on Linux; the server works headless there.
"""
from __future__ import annotations

import glob
import os
import plistlib
import re
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

MAC = sys.platform == "darwin"
WIN = sys.platform.startswith("win")

# bundle identifiers used by ChemDraw over the years (newest first)
BUNDLE_IDS = ["com.revvity.ChemDraw", "com.perkinelmer.ChemDraw", "com.cambridgesoft.ChemDraw"]

MAC_APP_DIRS = ["/Applications", "~/Applications"]
# companion apps whose name also starts with "ChemDraw" but which are not the drawing application
_MAC_REJECT = ("collections", "uninstall", "updater", "helper", "activation", "installer")

WIN_REGISTRY_KEYS = [
    r"SOFTWARE\Revvity Signals\ChemDraw",
    r"SOFTWARE\PerkinElmer\ChemDraw",
    r"SOFTWARE\CambridgeSoft\ChemDraw",
    r"SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\ChemDraw.exe",
]
WIN_GLOBS = [
    r"C:\Program Files*\Revvity Signals\ChemDraw*\ChemDraw.exe",
    r"C:\Program Files*\PerkinElmerInformatics\ChemOffice*\ChemDraw\ChemDraw.exe",
    r"C:\Program Files*\CambridgeSoft\ChemOffice*\ChemDraw\ChemDraw.exe",
]

_VERSION_RE = re.compile(r"\d+(?:\.\d+)*")


# ---------------------------------------------------------------------------
# directories
# ---------------------------------------------------------------------------

def home_dir() -> str:
    """``$CHEMDRAWMCP_HOME`` or ``~/.chemdrawmcp`` (created)."""
    d = os.environ.get("CHEMDRAWMCP_HOME") or os.path.join(os.path.expanduser("~"), ".chemdrawmcp")
    d = os.path.abspath(os.path.expanduser(d))
    os.makedirs(d, exist_ok=True)
    return d


def output_dir() -> str:
    """``$CHEMDRAWMCP_OUTPUT_DIR`` or ``<home>/renders`` (created); where renders/exports are written."""
    d = os.environ.get("CHEMDRAWMCP_OUTPUT_DIR") or os.path.join(home_dir(), "renders")
    d = os.path.abspath(os.path.expanduser(d))
    os.makedirs(d, exist_ok=True)
    return d


def platform_note() -> str:
    if MAC:
        return "macOS: live ChemDraw (AppleScript/JXA) + headless CDXML/RDKit engine"
    if WIN:
        return "Windows: live ChemDraw (COM automation) + headless CDXML/RDKit engine"
    return "Linux: headless mode only \u2014 CDXML/RDKit engine, no live ChemDraw"


# ---------------------------------------------------------------------------
# versions
# ---------------------------------------------------------------------------

def version_tuple(s: Optional[str]) -> Tuple[int, ...]:
    """First dotted number in ``s`` as a tuple of ints (``"ChemDraw 23.1.2.app"`` -> ``(23, 1, 2)``)."""
    m = _VERSION_RE.search(s or "")
    return tuple(int(x) for x in m.group(0).split(".")) if m else ()


def _version_from_path(path: str) -> Optional[str]:
    """``...\\ChemOffice 20.1\\ChemDraw\\ChemDraw.exe`` -> ``"20.1"``; ``ChemDraw 23.1.2.app`` -> ``"23.1.2"``."""
    for part in reversed(os.path.normpath(path).replace("\\", "/").split("/")):
        m = re.search(r"(?i)chem(?:draw|office)\D{0,3}(\d+(?:\.\d+)*)", part)
        if m:
            return m.group(1)
    return None


def _win_file_version(path: str) -> Optional[str]:
    """Version resource of a Windows executable via ``version.dll`` (ctypes only)."""
    if not WIN:
        return None
    try:
        import ctypes
        from ctypes import wintypes
        ver = ctypes.windll.version  # type: ignore[attr-defined]
        size = ver.GetFileVersionInfoSizeW(path, None)
        if not size:
            return None
        buf = ctypes.create_string_buffer(size)
        if not ver.GetFileVersionInfoW(path, 0, size, buf):
            return None
        ptr = ctypes.c_void_p()
        length = wintypes.UINT()
        if not ver.VerQueryValueW(buf, "\\", ctypes.byref(ptr), ctypes.byref(length)) or not ptr.value:
            return None
        fixed = ctypes.cast(ptr, ctypes.POINTER(ctypes.c_uint32 * 4)).contents  # VS_FIXEDFILEINFO head
        ms, ls = fixed[2], fixed[3]
        return "%d.%d.%d.%d" % (ms >> 16, ms & 0xFFFF, ls >> 16, ls & 0xFFFF)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# macOS bundles
# ---------------------------------------------------------------------------

def _app_root(path: str) -> str:
    """The ``.app`` bundle containing ``path`` (``path`` itself when it is the bundle)."""
    p = os.path.abspath(path)
    while p and p != os.path.dirname(p):
        if p.lower().endswith(".app"):
            return p
        p = os.path.dirname(p)
    return os.path.abspath(path)


def _read_info_plist(app: str) -> dict:
    try:
        with open(os.path.join(_app_root(app), "Contents", "Info.plist"), "rb") as f:
            return plistlib.load(f) or {}
    except Exception:
        return {}


def chemdraw_bundle_ids() -> List[str]:
    """Known ChemDraw bundle identifiers (current and historical vendors)."""
    return list(BUNDLE_IDS)


def macos_bundle_id(path: str) -> Optional[str]:
    """``CFBundleIdentifier`` of the ``.app`` at (or containing) ``path``."""
    v = _read_info_plist(path).get("CFBundleIdentifier")
    return str(v) if v else None


def chemdraw_version(path: Optional[str]) -> Optional[str]:
    """Version string of the ChemDraw install at ``path``.

    macOS: ``CFBundleShortVersionString`` from ``Contents/Info.plist``.  Windows: the
    executable's version resource, else whatever version appears in the path.
    """
    if not path:
        return None
    info = _read_info_plist(path)
    v = info.get("CFBundleShortVersionString") or info.get("CFBundleVersion")
    if v:
        return str(v)
    return _win_file_version(path) or _version_from_path(path)


def _mac_name_ok(app: str) -> bool:
    """Accept ``ChemDraw.app``, ``ChemDraw 23.1.2.app``, ``ChemDraw Professional 22.app``; reject companions."""
    name = os.path.basename(app)
    if name.lower().endswith(".app"):
        name = name[:-4]
    low = name.lower()
    if not low.startswith("chemdraw"):
        return False
    rest = low[len("chemdraw"):].strip(" -_")
    return not any(word in rest for word in _MAC_REJECT)


def _mac_app_key(app: str) -> Tuple[Tuple[int, ...], float]:
    ver = version_tuple(os.path.basename(app)) or version_tuple(chemdraw_version(app))
    try:
        mtime = os.path.getmtime(app)
    except OSError:
        mtime = 0.0
    return ver, mtime


def _best_mac_app(apps: List[str], check_name: bool = True) -> Optional[str]:
    """The highest-versioned ChemDraw bundle among ``apps`` (version from the name, else Info.plist)."""
    apps = [a for a in dict.fromkeys(os.path.abspath(a) for a in apps)
            if os.path.isdir(a) and (not check_name or _mac_name_ok(a))]
    if not apps:
        return None
    return max(apps, key=_mac_app_key)


def _mac_candidates() -> List[str]:
    out: List[str] = []
    for d in MAC_APP_DIRS:
        out += glob.glob(os.path.join(os.path.expanduser(d), "ChemDraw*.app"))
    return out


def _mdfind_apps() -> List[str]:
    """Spotlight lookup by bundle identifier (macOS only)."""
    mdfind = shutil.which("mdfind")
    if not mdfind:
        return []
    out: List[str] = []
    for bid in BUNDLE_IDS:
        try:
            r = subprocess.run([mdfind, "kMDItemCFBundleIdentifier == '%s'" % bid],
                               capture_output=True, text=True, timeout=15)
        except Exception:
            continue
        for line in (r.stdout or "").splitlines():
            line = line.strip()
            if (line.lower().endswith(".app") and os.path.isdir(line)
                    and "/.Trash/" not in line and "Backups.backupdb" not in line):
                out.append(line)
    return out


# ---------------------------------------------------------------------------
# Windows installs
# ---------------------------------------------------------------------------

def _dedupe(paths: List[str]) -> List[str]:
    seen, out = set(), []
    for p in paths:
        k = os.path.normcase(os.path.abspath(p))
        if k not in seen:
            seen.add(k)
            out.append(os.path.abspath(p))
    return out


def _newest(paths: List[str]) -> Optional[str]:
    paths = [p for p in _dedupe(paths) if os.path.isfile(p)]
    return max(paths, key=os.path.getmtime) if paths else None


def _exe_from_registry_value(value: str) -> Optional[str]:
    v = os.path.expandvars((value or "").strip().strip('"'))
    if not v:
        return None
    if v.lower().endswith("chemdraw.exe") and os.path.isfile(v):
        return v
    if os.path.isdir(v):
        for cand in (os.path.join(v, "ChemDraw.exe"), os.path.join(v, "ChemDraw", "ChemDraw.exe")):
            if os.path.isfile(cand):
                return cand
    return None


def _reg_exe_paths(winreg, root, key: str, depth: int) -> List[str]:
    found: List[str] = []
    try:
        h = winreg.OpenKey(root, key)
    except OSError:
        return found
    with h:
        values: List[str] = []
        try:
            values.append(winreg.QueryValue(h, None))  # (Default)
        except OSError:
            pass
        i = 0
        while True:
            try:
                _name, val, _typ = winreg.EnumValue(h, i)
            except OSError:
                break
            i += 1
            if isinstance(val, str):
                values.append(val)
        for v in values:
            exe = _exe_from_registry_value(v)
            if exe:
                found.append(exe)
        if depth > 0:  # per-version subkeys such as ...\ChemDraw\23.0
            j = 0
            while True:
                try:
                    sub = winreg.EnumKey(h, j)
                except OSError:
                    break
                j += 1
                found += _reg_exe_paths(winreg, root, key + "\\" + sub, depth - 1)
    return found


def _win_registry_candidates() -> List[str]:
    try:
        import winreg  # type: ignore[import-not-found]
    except ImportError:
        return []
    keys = list(WIN_REGISTRY_KEYS) + [k.replace("SOFTWARE\\", "SOFTWARE\\WOW6432Node\\", 1) for k in WIN_REGISTRY_KEYS]
    out: List[str] = []
    for root in (winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER):
        for key in keys:
            out += _reg_exe_paths(winreg, root, key, depth=2)
    return _dedupe(out)


def _win_glob_candidates() -> List[str]:
    drive = os.environ.get("SystemDrive", "C:")
    out: List[str] = []
    for pat in WIN_GLOBS:
        out += glob.glob(os.path.expandvars(pat.replace("C:", drive, 1)))
    return out


# ---------------------------------------------------------------------------
# the application
# ---------------------------------------------------------------------------

def find_chemdraw() -> Optional[str]:
    """Path of the ChemDraw application: the ``.app`` bundle on macOS, ``ChemDraw.exe`` on Windows, else None."""
    for var in ("CHEMDRAW_APP", "CHEMDRAW_EXE"):
        env = os.environ.get(var)
        if env:
            p = os.path.expanduser(env)
            if os.path.exists(p):
                return os.path.abspath(p)
    if MAC:
        return _best_mac_app(_mac_candidates()) or _best_mac_app(_mdfind_apps(), check_name=False)
    if WIN:
        return _newest(_win_registry_candidates()) or _newest(_win_glob_candidates())
    return None  # ChemDraw does not run on Linux


def chemdraw_info() -> Dict[str, Optional[str]]:
    """Everything a status report needs: path, version, bundle id (macOS) and the platform mode."""
    path = find_chemdraw()
    return {
        "path": path,
        "version": chemdraw_version(path) if path else None,
        "bundle_id": macos_bundle_id(path) if (path and MAC) else None,
        "platform": platform_note(),
    }


# ---------------------------------------------------------------------------
# launching
# ---------------------------------------------------------------------------

def app_launch_env() -> Dict[str, str]:
    env = dict(os.environ)
    env.setdefault("CHEMDRAWMCP_HOME", home_dir())
    return env


def popen_detached(argv: List[str], cwd: Optional[str] = None,
                   extra_env: Optional[Dict[str, str]] = None) -> subprocess.Popen:
    """Start a GUI application that outlives the MCP server process."""
    env = app_launch_env()
    if extra_env:
        env.update(extra_env)
    kw = dict(cwd=cwd, env=env, stdin=subprocess.DEVNULL,
              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if WIN:
        kw["creationflags"] = 0x00000008 | 0x00000200  # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    else:
        kw["start_new_session"] = True
    return subprocess.Popen(argv, **kw)
