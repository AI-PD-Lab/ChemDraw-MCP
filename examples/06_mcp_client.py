"""Talk to the chemdraw-mcp server over stdio with the official MCP Python client.

Starts `.venv/bin/chemdraw-mcp` from the repository root (fallback: `python -m chemdrawmcp.server`), lists
the tools, then calls `status`, `draw_reaction` (open_in_chemdraw=False) and `properties`.
Headless by default (CHEMDRAWMCP_BACKEND=none is passed to the server); CHEMDRAWMCP_EXAMPLES_LIVE=1 lets it use ChemDraw.
Run:   python examples/06_mcp_client.py            (works with mcp 1.x and 2.x)
"""
import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, get_default_environment, stdio_client

ROOT = Path(__file__).resolve().parents[1]
ESTERIFICATION = "acetic acid + ethanol --[H2SO4 (cat.)]--> ethyl acetate + O"


def server_params() -> StdioServerParameters:
    exe = ROOT / ".venv" / ("Scripts/chemdraw-mcp.exe" if os.name == "nt" else "bin/chemdraw-mcp")
    cmd, args = (str(exe), []) if exe.exists() else (sys.executable, ["-m", "chemdrawmcp.server"])
    env = get_default_environment()                        # the client passes only a minimal environment
    env.update({k: v for k, v in os.environ.items() if k.startswith(("CHEMDRAWMCP_", "CHEMDRAW_"))})
    if os.environ.get("CHEMDRAWMCP_EXAMPLES_LIVE") != "1":
        env["CHEMDRAWMCP_BACKEND"] = "none"                  # never touch a running ChemDraw
    env.setdefault("PYTHONPATH", str(ROOT))                  # lets the `python -m` fallback find the package
    print("starting:", cmd, *args)
    return StdioServerParameters(command=cmd, args=args, env=env, cwd=str(ROOT))


def field(obj, name):
    """mcp 2.x result fields are snake_case (server_info, is_error, mime_type); mcp 1.x used camelCase."""
    camel = name.split("_")[0] + "".join(p.title() for p in name.split("_")[1:])
    return getattr(obj, name) if hasattr(obj, name) else getattr(obj, camel, None)


def show(title, result, paths):
    """Print selected dotted keys (e.g. 'render.path') of a tool's JSON result; images are summarised."""
    print("\n== %s ==%s" % (title, "  (error)" if field(result, "is_error") else ""))
    for block in result.content:
        if block.type == "image":
            print("  [image %s, %d base64 chars]" % (field(block, "mime_type"), len(block.data)))
            continue
        try:
            data = json.loads(block.text)
        except (ValueError, AttributeError):
            print(" ", str(getattr(block, "text", block))[:600])
            continue
        if "error" in data:
            print("  error:", data["error"])
        for p in paths:
            v = data
            for k in p.split("."):
                v = v.get(k) if isinstance(v, dict) else None
            print("  %-24s %s" % (p, v))


async def main() -> int:
    async with stdio_client(server_params()) as (read, write):
        async with ClientSession(read, write) as session:
            init = await session.initialize()
            info = field(init, "server_info")
            print("connected to server %r, protocol %s" % (info.name, field(init, "protocol_version")))
            tools = (await session.list_tools()).tools
            print("%d tools: %s" % (len(tools), ", ".join(t.name for t in tools)))
            show("status", await session.call_tool("status", {}),
                 ["backend", "chemdraw_available", "headless_reason", "renders_dir", "default_style", "rdkit"])
            show("draw_reaction", await session.call_tool("draw_reaction", {
                "reaction": ESTERIFICATION, "conditions_below": ["reflux, 4 h"], "product_labels": ["85%", None],
                "arrow": "equilibrium", "style": "acs", "open_in_chemdraw": False}),
                 ["reaction_smiles", "cdxml_path", "opened", "render.path", "render.rendered_by", "render.analysis.summary"])
            show("properties", await session.call_tool("properties", {"structure": "aspirin"}),
                 ["rdkit.formula", "rdkit.molecular_weight", "rdkit.inchikey", "rdkit.logp", "rdkit.tpsa", "rdkit.elemental_analysis"])
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
