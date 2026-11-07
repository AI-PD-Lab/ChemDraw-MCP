#!/usr/bin/env python3
"""Bootstrap installer for chemdrawmcp (the ChemDraw MCP server).

    python3 install.py           # create .venv, install chemdrawmcp, register with Claude Code, Codex & Claude Desktop
    python3 install.py --check   # only report what is found (ChemDraw app, RDKit, server command)
    python3 install.py --help    # all options

Works on macOS, Linux and Windows with Python 3.10+ (uses `uv` when available).
On Linux ChemDraw itself is unavailable; the server runs in headless (CDXML/RDKit) mode.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VENV = os.path.join(HERE, ".venv")
WIN = sys.platform.startswith("win")
VENV_PY = os.path.join(VENV, "Scripts" if WIN else "bin", "python.exe" if WIN else "python")


def run(cmd, **kw):
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ensure_venv():
    if os.path.exists(VENV_PY):
        return
    uv = shutil.which("uv")
    if uv:
        run([uv, "venv", "--python", ">=3.10", VENV])
        return
    if sys.version_info < (3, 10):
        sys.exit("Python 3.10+ is required (found %s). Install a newer Python or `pip install uv` and rerun." % sys.version.split()[0])
    run([sys.executable, "-m", "venv", VENV])


def install_package():
    target = HERE + ("[windows]" if WIN else "")   # pywin32 for COM automation on Windows
    uv = shutil.which("uv")
    if uv:
        run([uv, "pip", "install", "--python", VENV_PY, "-e", target])
    else:
        run([VENV_PY, "-m", "pip", "install", "--upgrade", "pip", "-q"])
        run([VENV_PY, "-m", "pip", "install", "-e", target])


def main():
    args = sys.argv[1:]
    if "--help" in args or "-h" in args:
        print(__doc__, flush=True)
        subprocess.run([VENV_PY, "-m", "chemdrawmcp.install", "--help"]) if os.path.exists(VENV_PY) else None
        return
    if "--check" in args and not os.path.exists(VENV_PY):
        # read-only report before anything is installed: the installer module is pure standard library
        env = dict(os.environ, PYTHONPATH=HERE + os.pathsep + os.environ.get("PYTHONPATH", ""))
        sys.exit(subprocess.call([sys.executable, "-m", "chemdrawmcp.install", "--check"], env=env, cwd=HERE))
    ensure_venv()
    if not ("--check" in args or "--uninstall" in args or "--dry-run" in args) or not _package_installed():
        install_package()
    sys.exit(subprocess.call([VENV_PY, "-m", "chemdrawmcp.install"] + args))


def _package_installed() -> bool:
    try:
        r = subprocess.run([VENV_PY, "-c", "import chemdrawmcp.install"], capture_output=True, timeout=60)
        return r.returncode == 0
    except Exception:
        return False


if __name__ == "__main__":
    main()
