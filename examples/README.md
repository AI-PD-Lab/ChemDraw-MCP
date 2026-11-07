# chemdraw-mcp examples

Runnable scripts for the two halves of the project: the **library** (`chemdrawmcp.chem`, `layout`, `cdxml`,
`styles`: RDKit + the CDXML writer/reader, no ChemDraw required) and the **MCP server** (`chemdrawmcp.server`,
spoken to over stdio by Claude Code, Codex or any MCP client).  `prompts.md` collects prompts for agents.

## Setup

```bash
cd chemdraw-mcp
python install.py                 # or: python -m venv .venv && .venv/bin/pip install -e ".[dev]"
```

All scripts are short, stand-alone and exit 0 without ChemDraw.  Output (CDXML, PNG) goes to
`$CHEMDRAWMCP_OUTPUT_DIR`, by default `~/.chemdrawmcp/renders`; the examples print every path they write.

```bash
# headless: CDXML engine + RDKit only, no app, no network
export CHEMDRAWMCP_BACKEND=none CHEMDRAWMCP_OUTPUT_DIR=/tmp/chemdrawmcp-examples CHEMDRAWMCP_ONLINE=0
for s in examples/0*.py; do .venv/bin/python "$s" || echo "FAILED $s"; done
```

On Windows use `.venv\Scripts\python.exe` and `set VAR=value` (or `$env:VAR = "value"` in PowerShell).

## The examples

| script | shows |
|--------|-------|
| `01_structures_grid.py` | `chem.mol_from_input` (names from the offline dictionary, SMILES), `layout.build_structures` for a labelled 3 x 2 grid in ACS style, `doc.summary()` (every object with its bounding box), `cdxml.read_cdxml` to get SMILES / formula / MW back out of the file |
| `02_reaction.py` | one reaction with `layout.Step` / `Component`: reagents above the arrow, conditions below, a yield label, an equilibrium arrow, Nature style; the `<scheme>/<step>` record read back as reaction SMILES; an RDKit preview PNG via `chem.depict_png(chem.Reaction(...))` |
| `03_scheme.py` | `layout.build_scheme`: a four-step synthesis (benzene -> bromobenzene -> PhMgBr -> benzoic acid -> methyl benzoate) with a bold title, chained intermediates and rows wrapped at 400 pt; prints the reaction SMILES of each step |
| `04_read_cdxml.py` | the reader as a command-line tool: `python examples/04_read_cdxml.py file.cdxml` prints molecules, captions, arrows, graphics and reaction steps as tables (default input: the scheme from 03, also for files saved by ChemDraw itself) |
| `05_properties.py` | `chem.properties` (formula, MW, exact mass, InChIKey, logP, TPSA, HBD/HBA, Lipinski, QED), `chem.elemental_analysis`, `chem.describe`, and `chem.parse_reaction` + `reaction_properties` for balance, masses and atom economy, with a name resolved offline (`resolve_name("aspirin")`) |
| `06_mcp_client.py` | an MCP stdio client (`mcp.client.stdio.stdio_client` + `ClientSession`): starts the server, lists its 43 tools, calls `status`, `draw_reaction(open_in_chemdraw=False)` and `properties`, and prints the JSON results (the image content block of the render is summarised) |
| `07_custom_cdxml.py` | the low-level `CDXMLDocument` API: coloured fragments and atom colours, `add_plus`, `add_arrow`, a curved electron-pushing arrow (`add_curved_arrow`), captions with bold / italic / formula faces, `add_step`, a rounded-rectangle frame (`add_graphic`), all in ChemDraw's dark-mode style; every call is explained in a comment |

Each `.cdxml` the examples write opens directly in ChemDraw (File > Open, or drag it onto the app); the
exact settings of the stationery (bond length, fonts, line widths, page size) are in the file header.

## Headless vs. with ChemDraw

* **01 - 05, 07 never talk to ChemDraw.**  They use the library only, so they behave the same on macOS, Linux
  and Windows, in CI, and with ChemDraw open in another window.  `CHEMDRAWMCP_ONLINE=0` keeps name resolution
  offline (built-in dictionary + OPSIN if installed); set `CHEMDRAWMCP_ONLINE=1` to let unknown names fall back
  to PubChem.
* **06 talks to the server**, which by default is started with `CHEMDRAWMCP_BACKEND=none`, so `status` reports
  `chemdraw_available: false`, the reaction is laid out as CDXML and rendered with RDKit, and `properties` uses
  RDKit.  Set `CHEMDRAWMCP_EXAMPLES_LIVE=1` to pass your real environment through instead: the server then finds
  ChemDraw (macOS via AppleScript, Windows via COM), `status` lists the open documents, and chemical names that the
  built-in dictionary and OPSIN do not know are sent to ChemDraw's Name=Struct (this may launch the app unless
  `CHEMDRAWMCP_AUTOLAUNCH=0`).  The
  reaction is still only written to disk because the call passes `open_in_chemdraw: false`; change that to
  `true` (or drop it) to have it open as a new ChemDraw document and be rendered by ChemDraw itself.
* To use the server from Claude Code / Codex / Claude Desktop instead of a script, run `python install.py`
  once and then use the prompts in `prompts.md`.

## Adapting the examples

* Styles: any key from `chemdrawmcp.styles.STYLES` or the `styles` tool (`acs`, `nature`, `rsc`, `wiley`,
  `science_of_synthesis`, `slide`, `dark`, ...); aliases such as `jacs` work too.  `settings={...}` and
  `bond_length=` override single values.
* Inputs: everything `chem.mol_from_input` accepts: SMILES, InChI, molblock/SDF, CDXML, a file path or a name.
* Arrows: `full`, `equilibrium`, `retrosynthetic`, `double`/`resonance`, `dashed`, `bold`, `hollow`, `nogo`,
  `half_left`, `half_right`, `unbalanced_equilibrium` (see `cdxml.ARROW_KINDS`).
* The unit tests in `../tests` show more of the API (`test_cdxml.py` for the writer/reader and layouts,
  `test_server.py` for every tool in headless mode).
