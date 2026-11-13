# chemdraw-mcp in detail

Reference material that the [README](../README.md) leaves out: the complete tool table, environment variables,
the Agent Chat internals, how the two engines work, development and known limits.
Per-tool parameters and return values are in [tools.md](tools.md); agent workflows in [agent-guide.md](agent-guide.md).

**Contents:** [Installer steps](#installer-steps) · [Tool overview](#tool-overview) · [Styles](#styles) ·
[Environment variables](#environment-variables) · [Agent Chat in detail](#agent-chat-in-detail) ·
[How it works](#how-it-works) · [Development](#development) · [Notes and limits](#notes-and-limits)

## Installer steps

| step | effect |
|------|--------|
| Claude Code | adds `chemdraw` to `~/.claude.json` (user scope). Run `/mcp` in Claude Code to check. |
| Codex CLI | adds `[mcp_servers.chemdraw]` to `~/.codex/config.toml` (`codex mcp list`). |
| Claude Desktop | adds the server to `claude_desktop_config.json` if that file exists. |

Config files are backed up as `*.chemdrawmcp-backup`.  Manual registration for any MCP client: command
`<repo>/.venv/bin/chemdraw-mcp` (`.venv\Scripts\chemdraw-mcp.exe` on Windows), stdio transport, no arguments.

ChemDraw is found in `/Applications/ChemDraw*.app` (newest version) or via `CHEMDRAW_APP` / `CHEMDRAW_EXE`.
Python ≥ 3.10.  macOS: the first automation call asks for permission to control ChemDraw and System Events
(Automation / Accessibility) — grant it to your terminal or IDE.  `screenshot` also needs Screen Recording.

## Tool overview

| tool | purpose |
|------|---------|
| `status`, `launch`, `quit_chemdraw`, `activate_chemdraw` | discover / start ChemDraw; open documents; pending dialogs |
| *Agent Chat add-in* + `chemdraw-agent-bridge` | Claude / Codex chat panel inside ChemDraw (see below) |
| `documents`, `new_document(style)`, `open_document`, `close_document`, `save_document(path, format)` | document lifecycle; export to cdxml, cdx, mol / sdf / rxn (V2000) and molv3000 / sdfv3000 / rxnv3000, rdf, cml, ct, png, svg, pdf, eps, tiff, gif, jpeg, bmp, 3mf |
| `draw_structures(structures, labels, style, columns)` | grid of molecules from SMILES / InChI / molblock / CDXML / names |
| `draw_reaction(reaction | reactants+products, reagents_above, conditions_below, arrow, labels)` | one reaction with reagents, conditions, yields, atom mapping; arrows: full, equilibrium, retrosynthetic, double, dashed, bold, nogo, hollow |
| `draw_scheme(steps, title, max_width)` | multi-step synthesis, chained intermediates, wraps into rows |
| `draw_mechanism(stages, title)` | mechanism stages with curved electron-pushing arrows (`lp:0->0-3`, `1=2->2`, `h:3->3-4`, `:1e` fishhooks) and reaction/resonance arrows between stages |
| `named_reactions`, `apply_named_reaction(reaction, reactants)`, `suggest_reactions(reactants)` | library of 100 named-reaction templates (Suzuki, Grignard, Wittig, amide coupling, Boc, CuAAC, Diels–Alder, ...): browse, predict products and draw the atom-mapped reaction with typical reagents/conditions, or list what a set of reactants can do |
| `insert(data, kind)` | paste SMILES / name / CDXML / MOL / InChI / HELM / peptide / DNA / RNA / text through ChemDraw's importers |
| `add_caption`, `add_arrow` | extra text and arrows |
| `render(document, format, width)`, `quick_render(input)`, `screenshot` | images + numerical analysis + thumbnail |
| `read_document`, `read_file(path)` | molecules (SMILES, formula, MW, position), captions, arrows, reaction steps, settings — of an open document, or of a `.cdxml` / binary `.cdx` file on disk without opening it |
| `name_to_structure`, `structure_to_name` | ChemDraw Name=Struct / Struct=Name; names are resolved in this order: a built-in offline dictionary of ~3,100 reagents, solvents, drugs, amino acids, sugars, cofactors and natural products (7,800 aliases), OPSIN if installed, ChemDraw's Name=Struct, and PubChem only with `CHEMDRAWMCP_ONLINE=1` |
| `properties`, `analyze_reaction` | ChemDraw values (formula, MW, exact mass, elemental analysis) + RDKit (InChIKey, logP, TPSA, HBD/HBA, stereocentres, Lipinski, QED); reaction balance & atom economy |
| `convert(data, to)` | any format → any format (text results inline, files on disk) |
| `clean_structure`, `check_structure`, `predict_nmr(nucleus)` | ChemDraw clean-up, chemical warnings, 1H/13C NMR prediction (spectrum document + PNG) |
| `selection`, `select` | what is selected (SMILES, formula, MW, text, bounds) / select all, none, next, previous, invert |
| `commands(pattern)`, `run_command(command)` | ChemDraw's ~1600 named commands with descriptions and enabled state |
| `menus`, `do_menu(menu, item)` | menu items by title |
| `document_settings` | read/set bond length, line widths, fonts, sizes, magnification |
| `dismiss_dialogs`, `styles`, `analyze_image`, `compare_renders` | utilities |

Documents are addressed by name (window title) or window id; `None` means the front document.  ChemDraw renames
untitled documents after their first edit, so tools return and accept the stable window id.

## Styles

The `style` argument of the drawing tools uses the settings of the stationery shipped with ChemDraw 23:
`acs` (ACS Document 1996, default), `nature`, `rsc`, `rsc_1col`, `wiley`, `science_of_synthesis`, `j_het_chem`,
`j_chin_chem_soc`, `j_mol_mod_1col`, `j_mol_mod_2col`, `helv_chim_acta`, `adv_synth_catal`, `phytomedicine`,
`default`, `a4`, `slide`, `dark`, `idraw`, `terminal_carbons`.  `new_document(style)` opens the real stationery.

## Environment variables

| variable | meaning |
|----------|---------|
| `CHEMDRAW_APP` / `CHEMDRAW_EXE` | explicit ChemDraw application / executable |
| `CHEMDRAWMCP_HOME` | state directory (default `~/.chemdrawmcp`) |
| `CHEMDRAWMCP_OUTPUT_DIR` | where renders and CDXML files go (default `$CHEMDRAWMCP_HOME/renders`) |
| `CHEMDRAWMCP_AUTOLAUNCH=0` | never start ChemDraw from a tool call |
| `CHEMDRAWMCP_BACKEND=none` | force headless mode (CDXML + RDKit only) |
| `CHEMDRAWMCP_STYLE` | default style (default `acs`) |
| `CHEMDRAWMCP_ALLOW_ACTIVATE=1` | allow the backend to activate ChemDraw as a last resort (by default ChemDraw stays in the background and never takes keyboard focus) |
| `CHEMDRAWMCP_ONLINE=1` | allow PubChem name lookups when the offline dictionary, OPSIN and ChemDraw fail |
| `CHEMDRAWMCP_DEBUG=1` | log every JXA call and dialog poll with timings to stderr (macOS) |
| `CHEMDRAWMCP_LIVE=1` | enable the live ChemDraw tests (`pytest tests/test_live_macos.py -m live`) |
| `CHEMDRAWMCP_NO_FONTS=1`, `CHEMDRAWMCP_FONT_DIR` | force the heuristic text metrics / add a font directory for layout measurements |
| `CHEMDRAWMCP_CDX_PAIRS`, `CHEMDRAWMCP_SAMPLES` | directories with CDX/CDXML pairs and ChemDraw sample drawings for the reader tests |

## Agent Chat in detail

Like the PyMOL / ChimeraX panels in mol-mcp, ChemDraw gets a Claude-Code-style chat window of its own: an add-in
(**Add-ins ▸ Agent Chat**) that runs your installed Claude Code or Codex CLI through a small local bridge and streams
the agent's text, every `⏺ chemdraw/…` tool call with a collapsed `⎿` result preview, an animated status line and the
cost/tokens when a turn is done.  The agent uses the `chemdraw` MCP server, so it acts on the document you are looking
at; the panel also attaches your current selection (as SMILES) to each prompt.

```
python3 install.py            # on a Mac/PC with ChemDraw this also installs the add-in and a Scripts-menu starter
                              # (add --launch-agent to keep the bridge running permanently via launchd)
```

1. Restart ChemDraw once (add-ins are read at launch).
2. **Scripts ▸ Start Agent Bridge** (or run `chemdraw-agent-bridge` in a terminal; `--check` shows which CLIs it found).
3. **Add-ins ▸ Agent Chat**.  Enter sends, Shift+Enter newline, Esc interrupts, ↑/↓ recall prompts; `/model` picks
   provider · model · effort, `/claude` and `/codex` switch, `/clear` starts a new conversation, `/stop`, `/help`.

The bridge listens on `127.0.0.1:47325` only (`CHEMDRAWMCP_BRIDGE_PORT`), requires a per-install token
(`~/.chemdrawmcp/bridge_token`, baked into the add-in's `config.js`) and refuses requests from web origins, so other
websites cannot drive it.  Conversations resume across panel reopenings (`~/.chemdrawmcp/agent_sessions.json`).
CLI locations are found automatically (`CHEMDRAWMCP_CLAUDE_EXE` / `CHEMDRAWMCP_CODEX_EXE` override); the agent runs
with `~/.chemdrawmcp/workspace` as its working directory (`CHEMDRAWMCP_BRIDGE_CWD`).  Details: [agent-chat.md](agent-chat.md).

## How it works

* **macOS**: ChemDraw's AppleScript dictionary exposes documents (atoms, bonds, molecules, captions, arrows, the
  selection with SMILES / formula / MW / exact mass / elemental analysis, drawing settings), `open`, `save … as
  <type>` for every export format, menus and ~1600 named commands (`do`).  Data goes in through the pasteboard
  (`pasteAsSMILES`, `pasteAsCdxml`, `pasteAsStructureToName`, …) and out through `copyAsCdxml`/`copyAsSMILES`/…;
  the pasteboard is driven through NSPasteboard and your clipboard (every item and flavour, not just text) is
  snapshotted and restored inside the same script.  Modal alerts (e.g. "no reaction center found" when exporting
  RXN) are polled for every 2 s while a call runs, detected through Accessibility, dismissed, and their text is
  returned as the error.  Commands act on ChemDraw's front document window, which the backend reorders inside ChemDraw without
  activating the app, so ChemDraw stays in the background and never takes your keyboard focus.
* **Windows**: `chemdrawmcp/backends/windows.py` drives `ChemDraw.Application` through COM (pywin32 or comtypes)
  and can use the ChemScript SDK that ChemDraw installs.  It was written from the published object model without a
  Windows machine at hand — the candidate member names live in editable tables at the top of the file, and every
  error says which one failed.  Install with `python install.py` (adds pywin32).  See [windows.md](windows.md).
* **CDXML engine** (`chemdrawmcp/cdxml.py`, `layout.py`, `styles.py`): fragments from RDKit molecules (Kekulé bonds,
  wedges, charges, isotopes, radicals, atom maps), text with ChemDraw's formula face (automatic subscripts), straight
  and curved arrows of all kinds, plus signs, graphics, and `<scheme>/<step>` records so ChemDraw and the reader know
  what reacts to what.  Pages tile automatically so ChemDraw never asks to "expand the document".
* **Headless**: without ChemDraw, `draw_*` write CDXML and render with RDKit; `properties`, `convert`,
  `analyze_reaction`, `name_to_structure`, `named_reactions` / `apply_named_reaction` / `suggest_reactions` (100 reaction templates), `draw_mechanism` and
  `read_file` (own binary CDX reader) use RDKit and the offline name dictionary.

## Development

```bash
.venv/bin/python -m pytest -q          # ~2,070 unit tests, no ChemDraw needed (cdxml, cdx, layout, chem, names, reactions,
                                       # mechanism, server, installer, macOS + Windows backends against fakes)
CHEMDRAWMCP_LIVE=1 .venv/bin/python -m pytest tests/test_live_macos.py -m live -v   # 31 live tests against the real ChemDraw (macOS)
```

Layout: `chemdrawmcp/server.py` (tools), `backends/{base,macos,windows}.py`, `cdxml.py`, `cdx.py` (binary CDX reader),
`layout.py`, `textmetrics.py`, `abbreviations.py`, `mechanism.py` (electron-pushing arrows), `reactions.py` (named-reaction
templates), `styles.py`, `chem.py` (RDKit),
`imageanalysis.py`, `colornames.py`, `discovery.py`, `mcpcompat.py`, `install.py`.  Details, test fixtures and how to add a tool:
[development.md](development.md).

## Notes and limits

* ChemDraw's PNG export is always 300 dpi; `render(width=...)` resamples.  SVG/PDF/EPS are vector.
* Pasting into an existing document lets ChemDraw choose the position (the pasted objects stay selected so the user
  can move them); new documents get full layout control.
* Struct=Name, NMR prediction and some import formats depend on your ChemDraw edition/licence.
* ChemDraw is driven entirely in the background; only `activate_chemdraw` (or `CHEMDRAWMCP_ALLOW_ACTIVATE=1`) brings
  it forward.  Keep it that way while the live tests run: ChemDraw hotkeys relabel selected atoms if it ever has focus.
