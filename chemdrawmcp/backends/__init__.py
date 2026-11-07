"""ChemDraw automation backends (macOS AppleScript/JXA, Windows COM) and the selector."""
from __future__ import annotations

import os
import sys
from typing import Optional

from .base import (EXPORT_FORMATS, INSERT_KINDS, ChemDrawBackend, ChemDrawError, ChemDrawNotAvailable, CommandInfo,
                   DocumentInfo, NullBackend, normalize_format)

__all__ = ["get_backend", "ChemDrawBackend", "ChemDrawError", "ChemDrawNotAvailable", "CommandInfo", "DocumentInfo",
           "NullBackend", "EXPORT_FORMATS", "INSERT_KINDS", "normalize_format"]


def get_backend(kind: Optional[str] = None) -> ChemDrawBackend:
    """Pick the backend for this platform.  ``CHEMDRAWMCP_BACKEND=none`` forces headless mode."""
    kind = (kind or os.environ.get("CHEMDRAWMCP_BACKEND") or "auto").lower()
    if kind in ("none", "headless", "off"):
        return NullBackend("headless mode requested (CHEMDRAWMCP_BACKEND=none): no live ChemDraw")
    if kind in ("auto", "macos", "mac", "darwin") and sys.platform == "darwin":
        from .macos import MacBackend
        b = MacBackend()
        if b.app_path or kind != "auto":
            return b
        return NullBackend("ChemDraw.app was not found in /Applications (set CHEMDRAW_APP); running headless")
    if kind in ("auto", "windows", "win", "win32") and sys.platform.startswith("win"):
        from .windows import WindowsBackend
        return WindowsBackend()
    if kind == "auto":
        return NullBackend("ChemDraw only runs on macOS and Windows; this is %s, so only the headless CDXML/RDKit engine is available" % sys.platform)
    raise ValueError("unknown backend %r" % kind)
