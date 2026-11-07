"""Agent bridge: a localhost HTTP/SSE server that lets ChemDraw's Agent Chat add-in drive the user's Claude Code
or Codex CLI (which in turn acts on ChemDraw through the ``chemdraw`` MCP server).

Start it with ``chemdraw-agent-bridge`` (console script) or ``python -m chemdrawmcp.bridge.server``.  Other parts
of chemdrawmcp use :func:`bridge_info` / :func:`is_bridge_running` (lazily re-exported from
:mod:`chemdrawmcp.bridge.server`) to find a running bridge.
"""
from __future__ import annotations

DEFAULT_PORT = 47325

__all__ = ["DEFAULT_PORT", "BridgeServer", "bridge_health", "bridge_info", "is_bridge_running", "main"]


def __getattr__(name):  # lazy: keeps ``python -m chemdrawmcp.bridge.server`` free of the double-import warning
    if name in __all__:
        from . import server
        return getattr(server, name)
    raise AttributeError(name)
