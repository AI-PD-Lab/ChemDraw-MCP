# Tool reference

Every MCP tool exposed by `chemdrawmcp/server.py` (43 tools, server name `chemdraw`).  Verified against the
code by importing the module headless and listing `server.mcp._tool_manager.list_tools()`.

Contents: [Conventions](#conventions) · [Platform availability](#platform-availability) ·
[Lifecycle](#lifecycle) · [Documents](#documents) · [Drawing](#drawing) · [Rendering and looking](#rendering-and-looking) ·
[Chemistry](#chemistry) · [Selection, commands, menus](#selection-commands-menus) · [Utilities](#utilities)

## Conventions

**Results.**  Every tool returns a JSON string (`json.dumps(..., indent=1)`).  Tools that produce a raster
image return a two-element MCP result instead: the JSON text followed by a PNG thumbnail (long side at most
1024 px, `imageanalysis.thumbnail_png`) that the agent can look at.  Pass `return_image=false` to get the JSON
alone (the full-size file is still written to disk).

**Errors.**  Failures never raise; they come back as

```json
{"ok": false, "error": "no document named \"foo\"; open documents: [\"Untitled ACS Document 1996-1\"]",
 "error_type": "ChemDrawError"}
```

`error_type` is the Python exception class: `ChemDrawError` (anything ChemDraw-side, including the text of a
modal alert), `ChemDrawNotAvailable` (no live ChemDraw on this machine; the result then also carries a `hint`
saying which tools still work headless), `ChemError` (bad chemistry input), `KeyError` (unknown style), ...
Successful results do **not** always contain `"ok": true` (the `draw_*` tools, `status`, `read_document`,
`analyze_reaction` and others return plain data), so test for `"ok": false`.

**The `document` argument.**  Most live tools take `document`:

| value | meaning |
|-------|---------|
| omitted / `null` / `""` / `"front"` | the front document (first window; ChemDraw's auto-save documents are skipped) |
| integer, or a string like `"812"` / `"#812"` | a **window id** as returned in `window_id` by `documents`, `status` and every drawing tool |
| any other string | the document name as shown in the window title; on macOS matching is case-insensitive and also accepts the name without `.cdxml` or a unique prefix |

ChemDraw renames an untitled document after its first edit ("Untitled ACS Document 1996-1" becomes the
file name once saved, and stationery-based names change too), so address documents by `window_id` whenever
you got one.  On Windows `window_id` is the document's window handle when the COM object exposes one
(`Document.HWnd` or `Document.Window.HWnd`, see `docs/windows.md`); a number that matches no handle is taken
as a 1-based position in the document list, so use names there when `documents` reports `window_id: null`.

**Live vs headless.**  Tools marked *live* need the ChemDraw application.  They start it automatically when
it is installed but not running (`launch` is implicit) unless `CHEMDRAWMCP_AUTOLAUNCH=0`.  Tools marked
*live, no autolaunch* refuse with `"ChemDraw is not running; call `launch` first"` instead of starting it.
Tools marked *headless OK* work without ChemDraw (Linux, CI, `CHEMDRAWMCP_BACKEND=none`) using RDKit and the
CDXML engine.  Note that **chemical-name inputs that are not in the built-in dictionary (and not parsed by
OPSIN) are sent to ChemDraw's Name=Struct** whenever ChemDraw is installed, so e.g.
`properties(structure="2-chloro-N-(4-methoxyphenyl)acetamide")` can start ChemDraw on a Mac (dictionary names
such as `"aspirin"` never do); set `CHEMDRAWMCP_AUTOLAUNCH=0` to keep name resolution offline.

**Files.**  Everything the server writes goes to the renders directory reported by `status`
(`CHEMDRAWMCP_OUTPUT_DIR`, default `~/.chemdrawmcp/renders`) unless a tool takes an explicit `path`.
File names are `<tag>_<YYYYmmdd_HHMMSS_mmm>_<NN>.<ext>` where `<tag>` is a slug of the tool's `name`
argument (or a fixed tag such as `render`, `quick`, `screenshot`, `convert`, `nmr_1h`) and `<NN>` is a
per-process counter, so two calls in the same millisecond never collide.

**Render analysis.**  Wherever a result has an `analysis` key it is the dictionary produced by
`imageanalysis.analyze_render` on the PNG: `width`, `height`, `background` (`rgb`, `hex`, `name`,
`transparent`, `uniform_border`), `coverage` (fraction of non-background pixels), `empty`, `bbox` (content box:
`px` = `[x0, y0, x1, y1]` in pixels, `frac` = the same as fractions of the frame, `width_frac`, `height_frac`),
`centroid_frac`, `offset_from_centre_frac`, `margins_frac` (`left`/`right`/`top`/`bottom` as fractions of the
frame), `edge_coverage`, `clipped_edges`, `touching_edges`, `thirds`, `palette` (dominant colours as
`{name, hex, rgb, share, family, name_distance}`), `hue_shares`, `luminance`
(`fg_mean`, `fg_std`, `bg`, `contrast_ratio_fg_bg`), `edge_density`, `ascii` (20 rows x 48 columns density
map), `ascii_legend` and a one-paragraph `summary` such as *"Content covers 6% of the 1650x420 frame (bounding
box 92% x 61%). Centroid is near the centre. Free margins: left 4%, right 4%, top 20%, bottom 19%. Dominant
colours: black 100%. Background white; foreground/background contrast 21.0:1."*

**`DocumentInfo`.**  Wherever a document is described you get
`{"name", "path", "modified", "front", "atoms", "bonds", "molecules", "captions", "arrows", "window_id"}`.

## Platform availability

| tool | macOS (JXA) | Windows (COM) | headless |
|------|:-----------:|:-------------:|:--------:|
| `status`, `styles`, `analyze_image`, `compare_renders`, `analyze_reaction`, `named_reactions`, `suggest_reactions` | yes | yes | yes |
| `apply_named_reaction` | yes | yes | yes (prediction is RDKit; `draw=true` renders like `draw_reaction`) |
| `read_file` | yes | yes | yes (CDXML and binary CDX are parsed without ChemDraw) |
| `launch`, `quit_chemdraw`, `activate_chemdraw`, `documents`, `new_document`, `open_document`, `close_document`, `save_document`, `read_document`, `insert`, `add_caption`, `add_arrow`, `render`, `screenshot`, `selection`, `dismiss_dialogs` | yes | yes (best effort, see `docs/windows.md`) | no |
| `draw_structures`, `draw_reaction`, `draw_scheme`, `draw_mechanism`, `quick_render` | yes | yes | yes (CDXML file + RDKit PNG/SVG) |
| `name_to_structure`, `properties`, `convert` | yes | yes | yes (offline name dictionary, RDKit formats only) |
| `structure_to_name` | yes | yes (`Objects.Data("chemical/x-name")`) | no |
| `document_settings` | yes | read from the COM `Settings` object, else from the document's CDXML attributes | no |
| `clean_structure` | yes | `reaction=false` only | no |
| `select` | yes | `all` and `none` only | no |
| `commands`, `run_command` | all ~1600 named commands | curated list of ~19 commands (`COMMANDS` table) | no |
| `check_structure`, `predict_nmr`, `menus`, `do_menu` | yes | no (error) | no |

---

## Lifecycle

### `status`

ChemDraw availability and state.  Never launches ChemDraw.  *headless OK.*

No parameters.

Returns `server_version`, `platform` (`darwin`/`win32`/`linux`), `backend` (`macos`/`windows`/`none`),
`chemdraw_app` (path or null), `chemdraw_available`, `chemdraw_version` (from Info.plist / exe resource),
`renders_dir`, `home`, `rdkit` (bool), `styles` (list of style keys), `default_style`.  Headless adds
`headless_reason`.  With a backend: `running`; when running also `version` (as reported by the app),
`documents` (list of `DocumentInfo`, front first), `front_document` (name), and `pending_dialogs` + `hint`
when ChemDraw is showing an alert.  When installed but not running: `hint`.

```json
{}
```
```json
{"server_version": "0.1.0", "platform": "darwin", "backend": "macos",
 "chemdraw_app": "/Applications/ChemDraw 23.1.2.app", "chemdraw_available": true, "chemdraw_version": "23.1.2.7",
 "renders_dir": "/Users/me/.chemdrawmcp/renders", "home": "/Users/me/.chemdrawmcp", "rdkit": true,
 "styles": ["a4", "acs", "adv_synth_catal", "dark", "default", "..."], "default_style": "acs",
 "running": true, "version": "23.1.2.7",
 "documents": [{"name": "reaction_20261107_154217_512_03.cdxml", "path": "/Users/me/.chemdrawmcp/renders/reaction_20261107_154217_512_03.cdxml",
                "modified": false, "front": true, "atoms": 11, "molecules": 4, "captions": 3, "arrows": 1, "window_id": 812}],
 "front_document": "reaction_20261107_154217_512_03.cdxml"}
```

Call it first in a session: it tells you whether you are live or headless, where files go and which documents
are already open (and whether a dialog is blocking the app).  On macOS the document list comes from one
`snapshot()` round trip and omits `bonds`; `documents` has the full `DocumentInfo`.

### `launch`

Start ChemDraw if it is not running and wait until it answers.  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `wait_seconds` | float | 90 | how long to wait for the first successful automation call |

Returns (macOS) `{"running": true, "was_running": false, "version": "23.1.2.7", "app": "/Applications/ChemDraw 23.1.2.app", "seconds": 6.2}`;
(Windows) `{"already_running", "executable", "started_via" ("COM CreateObject(ChemDraw.Application)" or "exe"), "pid", "pids", "version", "prog_id", "seconds"}`.

On macOS the app is opened with `open -g` (in the background); the first automation call triggers the
Automation permission prompts for ChemDraw and System Events.  Other live tools call this implicitly.

### `quit_chemdraw`

Quit ChemDraw.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `save` | bool | false | `true` saves changed documents (`app.quit({saving: "yes"})`), `false` discards them |

Returns `{"ok": true}`.

### `activate_chemdraw`

Bring ChemDraw (and optionally one document) to the foreground and **leave it there**.  All other tools hand
keyboard focus back to the application the user was in; use this when the user wants to look at the drawing.
*live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document to raise (name or window id) |

Returns `{"ok": true}`.

---

## Documents

### `documents`

List open documents.  *live.*

No parameters.  Returns `{"documents": [DocumentInfo, ...]}` with the front document first; ChemDraw's hidden
auto-save documents (`ChemDraw Auto-Save ...`) are excluded.

### `new_document`

Create a document from a stationery/style.  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `style` | str | `"acs"` | a style key or alias (`acs`, `jacs`, `nature`, `rsc`, `rsc_1col`, `wiley`, `science_of_synthesis`, `j_het_chem`, `j_chin_chem_soc`, `j_mol_mod_1col`, `j_mol_mod_2col`, `helv_chim_acta`, `adv_synth_catal`, `phytomedicine`, `default`, `a4`, `slide`, `dark`, `idraw`, `terminal_carbons`; see `styles`), a stationery title (`"ACS Document 1996"`) or the path of a `.cds` file |

Returns `{"ok": true, ...DocumentInfo}`.

macOS opens `<ChemDraw.app>/Contents/Resources/SpecialPurpose/Stationery/<title>.cds` (ChemDraw then creates
"Untitled <title>-N"); Windows uses `Documents.Add(path)` with the file found in a `ChemDraw Items` folder.
The real stationery carries everything (fonts, colours, page setup, arrow defaults), whereas the `style`
argument of the drawing tools only reproduces the settings listed in `styles`.

```json
{"style": "nature"}
```
```json
{"ok": true, "name": "Untitled Nature Research Document-2", "path": null, "modified": false, "front": true,
 "atoms": 0, "bonds": 0, "molecules": 0, "captions": 0, "arrows": 0, "window_id": 815}
```

### `open_document`

Open a file in ChemDraw.  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `path` | str | required | any file ChemDraw reads: cdxml, cdx, mol, sdf, rxn, cml, skc, cif, jdx spectra, pdf, png, ... (`~` is expanded) |

Returns `{"ok": true, ...DocumentInfo}`.  If the same file is already open, ChemDraw just brings that window
forward and the tool returns its info after about half a second.  Waits up to 5 s for the new window.

### `close_document`

Close a document.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | front document by default |
| `save` | bool | false | save before closing |

Returns `{"ok": true, "remaining": ["name", ...]}`.

### `save_document`

Save or export a document.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `path` | str | required | target file; the matching extension is appended if missing (`out` -> `out.png`) |
| `document` | str | null | document to save |
| `format` | str | from the extension, else `cdxml` | `cdxml`, `cdx`, `mol` (MDL V2000), `molv3000`, `sdf` (V2000), `sdfv3000`, `rxn` (V2000), `rxnv3000`, `rdf`, `cml`, `ct`, `png`, `svg`, `pdf`, `eps`, `tiff`, `gif`, `jpeg`, `bmp`, `3mf` (aliases `jpg`, `tif`, `molfile`, `molv2000`, `sd`, `sdfv2000`, `rxnv2000`, `mol3000`, `v3000`, `chemdraw`, `xml`) |

Returns `{"ok": true, "path": "<actual path written>", "format": "png", "bytes": 48213}`.

Notes:
* An existing file at `path` is deleted first (so ChemDraw never asks "replace?").
* `rxn` needs a reaction: at least one arrow with reactants and products ChemDraw recognises.  Otherwise
  ChemDraw shows "No reaction center was found"; the alert is dismissed and returned as the error text.
* PNG is exported at ChemDraw's fixed 300 dpi.  SVG, PDF and EPS are vector.
* `mol` writes the first structure only as a V2000 molfile; use `sdf` for several and `molv3000` /
  `sdfv3000` / `rxnv3000` for V3000 output (`mol`, `sdf` and `rxn` are the V2000 variants).
* When ChemDraw appends its own extension the tool returns the file it actually wrote (newest matching file).

```json
{"path": "~/Desktop/scheme1.svg", "document": 812}
```

### `read_document`

Structured content of a document, obtained by exporting it as CDXML and parsing that (`cdxml.read_cdxml`).
Use it before editing a document the user made, and after your own edits to confirm what is on the page.
*live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document to read |
| `include_cdxml` | bool | false | also return the CDXML text (trimmed to 8000 characters: head and tail) |

Returns:

| key | content |
|-----|---------|
| `settings` | root attributes: `BondLength`, `LineWidth`, `BoldWidth`, `LabelSize`, `CaptionSize`, `BondSpacing`, `ChainAngle`, `HashSpacing`, `MarginWidth`, `CreationProgram`, `Name`, plus `LabelFont` / `CaptionFont` resolved to font names |
| `pages` | number of `<page>` elements |
| `molecules` | one entry per top-level `<fragment>` with atoms (nicknames expanded, groups flattened): `id`, `kind: "molecule"`, `bounds` (`x0,y0,x1,y1,width,height` in points), `atoms`, `bonds`, `smiles`, `formula`, `molecular_weight`, `heavy_atoms`; optionally `name` (from a ChemDraw chemical-property caption), `properties`, `labels` (atom labels that could not be interpreted), `warnings`; `error` when RDKit could not fully sanitize the drawing (valence error, unkekulizable ring) — the `smiles` is then a best-effort reading |
| `groups` | number of `<group>` elements |
| `other_objects`, `spectra`, `bracketed_groups` | present when the file holds things the reader passes through: counts of unhandled element types, number of `<spectrum>` objects, polymer/repeat brackets (`repeatCount`, `bracketUsage`, ...; repeat units are not expanded in the SMILES) |
| `texts` | captions (not atom labels): `id`, `kind: "text"`, `text`, `bounds`, optional `chemical_property` |
| `arrows` | `<arrow>` elements and arrow-type `<graphic>` lines: `id`, `bounds`, `head3D`, `tail3D`, `arrowheadHead`, `arrowheadTail`, `arrowheadType`, `lineType`, `noGo`, `angularSize` when present; a `<graphic>` carrying `SupersededBy` (the legacy twin ChemDraw writes next to every `<arrow>` for old readers) is skipped, so arrows are not counted twice |
| `graphics` | rectangles, ovals, lines: `type`, `rectangleType`, `ovalType`, `lineType` |
| `steps` | reaction steps from `<scheme>/<step>`: `reactants` / `products` (`{id, smiles, formula}`), `arrows` (ids), `above_arrow` / `below_arrow` (`{id, text}` or `{id, smiles, formula}`), `atom_map_pairs` |
| `reaction_smiles` | one `reactants>agents>products` string per step (only present when there are steps) |
| `warnings` | reader problems |
| `document` | `DocumentInfo` |
| `chemdraw_settings` | the live settings as returned by `document_settings` (absent if the backend cannot read them) |
| `cdxml` | the text when `include_cdxml=true` |

```json
{"document": 812}
```
```json
{"settings": {"BondLength": 14.4, "LineWidth": 0.6, "LabelFont": "Arial", "CaptionFont": "Arial", "Name": "reaction"},
 "pages": 1,
 "molecules": [{"id": "1", "kind": "molecule", "bounds": {"x0": 36, "y0": 54.8, "x1": 80.5, "y1": 89.2, "width": 44.5, "height": 34.4},
                "atoms": 4, "bonds": 3, "smiles": "CC(=O)O", "formula": "C2H4O2", "molecular_weight": 60.052, "heavy_atoms": 4}, "..."],
 "texts": [{"id": "9", "kind": "text", "text": "+", "bounds": {"...": 0}}, {"id": "24", "kind": "text", "text": "H2SO4", "bounds": {"...": 0}}],
 "arrows": [{"id": "23", "kind": "arrow", "bounds": {"...": 0}, "head3D": "190.9 72 0", "tail3D": "147.7 72 0", "arrowheadHead": "Full"}],
 "graphics": [],
 "steps": [{"id": "40", "reactants": [{"id": "1", "smiles": "CC(=O)O", "formula": "C2H4O2"}, {"id": "11", "smiles": "CCO", "formula": "C2H6O"}],
            "products": [{"id": "26", "smiles": "CCOC(C)=O", "formula": "C4H8O2"}, {"id": "35", "smiles": "O", "formula": "H2O"}],
            "arrows": ["23"], "above_arrow": [{"id": "24", "text": "H2SO4"}], "below_arrow": [], "atom_map_pairs": 0}],
 "reaction_smiles": ["CC(=O)O.CCO>>CCOC(C)=O.O"], "warnings": [],
 "document": {"name": "reaction_20261107_154217_512_03.cdxml", "window_id": 812, "atoms": 11, "molecules": 4, "captions": 3, "arrows": 1, "...": 0},
 "chemdraw_settings": {"name": "reaction_20261107_154217_512_03.cdxml", "bond_length": 14.4, "line_width": 0.6, "label_size": 10, "...": 0}}
```

Pitfalls: a user-drawn reaction without a `<scheme>` (ChemDraw records steps only after *Structure > Map
Reaction Atoms* / *Autonumber Reaction* or when the file was written by this server) has `steps: []` even
though arrows and molecules are present.  Nicknames ChemDraw has not expanded are looked up in the built-in
abbreviation table (`Ph`, `OMe`, `Boc`, `TBS`, `NBS`, ... in `chemdrawmcp/abbreviations.py`); unknown labels
become dummy atoms (`*` in the SMILES) and are listed under `labels`.

### `read_file`

Parse a ChemDraw file on disk without opening it in ChemDraw.  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `path` | str | required | a `.cdxml` file, or a binary `.cdx` file (recognised by the `VjCD0100` magic or the extension; converted with `chemdrawmcp.cdx`) |
| `include_cdxml` | bool | false | also return the CDXML text (for `.cdx` the converted text), trimmed to 8000 characters |

Returns the same structure as `read_document` (`settings`, `pages`, `groups`, `molecules`, `texts`, `arrows`,
`graphics`, `steps`, `reaction_smiles`, `warnings`, ...) plus `path`, but no `document` or
`chemdraw_settings`.  A missing file gives `"file not found: ..."`; a file that is neither CDXML nor CDX
gives an XML parse error.  Use it to inspect files the user points you to before deciding whether to open
them, and to check files you exported.

```json
{"path": "~/Desktop/scheme.cdx"}
```

### `document_settings`

Read or change the drawing settings of a document.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document |
| `settings` | object | null | keys to change; omitted = read only |

Keys (lengths in points): `bond_length`, `line_width`, `bold_width`, `margin_width`, `hash_spacing`,
`bond_spacing` (% of bond length), `chain_angle` (degrees), `label_size`, `caption_size`, `label_font`,
`caption_font`, `magnification` (%), `fractional_widths`, `rulers_visible`, `crosshair_visible`.
Returns `{"ok": true, "settings": {"name": "...", "bond_length": 14.4, ...}}` (the full set after the change).
macOS converts from ChemDraw's internal 1/20-point units.  Windows reads a COM `Settings` object when the
document has one (`"source": "com"`) and otherwise the attributes of the document's CDXML
(`"source": "cdxml"`); keys it could not read are listed under `unavailable`, and writing settings without a
`Settings` object only works when `Objects.SetData` is known to replace the document.

```json
{"document": 812, "settings": {"bond_length": 17, "label_size": 12}}
```

Changing settings does not rescale existing objects; it affects what is drawn next (and `Clean Up`).

---

## Drawing

The three `draw_*` tools share one pipeline: parse inputs with RDKit -> lay out a CDXML document in the
requested style (`layout.py`) -> write it to the renders directory -> open it in ChemDraw as a new document
(or paste it into `document`) -> export a PNG through ChemDraw -> analyse it.  They share these parameters:

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `style` | str | `"acs"` | style key / alias / stationery title (`styles`) |
| `bond_length` | float | style default | override the bond length in points |
| `document` | str | null | paste into this open document instead of creating a new one |
| `open_in_chemdraw` | bool | true | `false`: only write the CDXML (and an RDKit render) |
| `render` | bool | true | export a PNG and analyse it |
| `name` | str | tool-specific | tag for the output file names and the CDXML `Name` attribute |
| `return_image` | bool | true | attach the PNG thumbnail to the result |

Shared result keys:

| key | content |
|-----|---------|
| `cdxml_path` | the CDXML file written (always, also headless) |
| `layout` | `{"style": "ACS Document 1996", "bond_length": 14.4, "objects": [...], "steps": 1, "bounds": {...}}`; each object has `id`, `kind` and `bounds` plus extras: fragments `atom_ids`, `smiles`, `atoms`, `atom_positions` (page coordinates per RDKit atom index) and `labels` (the atom labels ChemDraw will draw, with their predicted boxes); texts `text`; straight arrows report their arrow style under `kind` (`"full"`, `"equilibrium"`, ...) with `tail`/`head`; curved arrows `kind: "curved"`, `angle`; graphics `kind: "graphic"`, `shape` |
| `layout_warnings` | only when `layout.layout_warnings` found something (overlapping fragment / text / arrow boxes, content wider than the page): a list of strings |
| `document`, `window_id`, `opened: true`, `atoms`, `molecules`, `captions`, `arrows` | when a new document was opened |
| `document`, `window_id`, `inserted: true`, `atoms_added`, `note` | when pasted into `document` (ChemDraw picks the position; the pasted objects stay selected so the user can drag them) |
| `opened: false`, `note` | headless or `open_in_chemdraw=false` |
| `render` | `{"path", "format": "png", "rendered_by": "ChemDraw", "bytes", "size_px", "dpi": [300, 300], "analysis"}` for a ChemDraw export.  For an RDKit fallback the shape depends on the tool: `{"path", "format": "png", "rendered_by": "RDKit...", "analysis"}` from `draw_reaction` / `quick_render`, `{"path", "rendered_by": "RDKit", "analysis"}` from `draw_structures`, and only `{"path", "rendered_by": "RDKit (components only ...)"}` / `{"path", "rendered_by": "RDKit (structures only ...)"}` with **no `analysis`** from `draw_scheme` / `draw_mechanism` |

**Note:** with `open_in_chemdraw=false` (or headless) the `render` is an **RDKit** picture of the molecules
you passed (`rendered_by` says so), not a rendering of the CDXML layout; open the CDXML in ChemDraw to see the
real page.  In paste mode (`document=...`) the whole target document is exported, so the PNG shows the
user's existing content too.

**Structure inputs** (everywhere a structure is accepted): SMILES, SMARTS (converted), InChI, a MOL/SDF block,
CDXML text, a path to a `.mol/.sdf/.cdxml/.smi` file, or a chemical name.  Names are resolved in this order:
the per-process cache -> the built-in dictionary (`chem.COMMON_NAMES`: ~3,100 structures under ~7,800 names
and aliases — solvents, reagents, drugs, amino acids, sugars, cofactors, natural products) -> OPSIN if installed
(`py2opsin`, an `opsin` executable or `$OPSIN_JAR`) -> ChemDraw's Name=Struct when ChemDraw is available (through a temporary document; this is the
only step that can launch ChemDraw) -> PubChem only when `CHEMDRAWMCP_ONLINE=1`.  Unresolvable names give
`"could not interpret 'X' as SMILES, InChI, molfile, CDXML or a chemical name"`; pass a SMILES instead.

### `draw_structures`

Molecules in a grid with optional captions.  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `structures` | list[str] | required | structure inputs (see above) |
| `labels` | list[str or null] | names for name inputs, else none | caption under each structure |
| `columns` | int | `min(4, ceil(sqrt(n)))` | grid columns |
| `colors` | list[str] | null | per-structure colour: `black`, `white`, `red`, `green`, `blue`, `orange`, `magenta`, `cyan`, `gray`, `purple`, `yellow`, `brown` or `#rrggbb` |
| + shared parameters (`name` defaults to `"structures"`) | | | |

Extra result key `molecules`: one summary string per input, e.g.
`"aspirin: C9H8O4, MW 180.159, SMILES CC(=O)Oc1ccccc1C(=O)O, neutral, no stereocentres"`.
Returns an image.

```json
{"structures": ["aspirin", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", "caffeine"],
 "labels": ["1", "(S)-2", "3"], "style": "nature", "columns": 3}
```
```json
{"cdxml_path": "/Users/me/.chemdrawmcp/renders/structures_20261107_160102_118_05.cdxml",
 "layout": {"style": "Nature Research Document", "bond_length": 10.8, "objects": ["..."], "steps": 0, "bounds": {"x0": 36, "y0": 36, "x1": 178.3, "y1": 83.5, "width": 142.3, "height": 47.5}},
 "document": "structures_20261107_160102_118_05.cdxml", "opened": true, "window_id": 819, "atoms": 41, "molecules": 3, "captions": 3, "arrows": 0,
 "molecules": ["aspirin: C9H8O4, MW 180.159, SMILES CC(=O)Oc1ccccc1C(=O)O, neutral, no stereocentres",
               "C13H18O2, MW 206.285, SMILES CC(C)Cc1ccc([C@H](C)C(=O)O)cc1, neutral, 1 stereocentre(s): 1 assigned, 0 unassigned",
               "caffeine: C8H10N4O2, MW 194.194, SMILES Cn1c(=O)c2c(ncn2C)n(C)c1=O, neutral, no stereocentres"],
 "render": {"path": "/Users/me/.chemdrawmcp/renders/structures_20261107_160103_402_06.png", "format": "png", "rendered_by": "ChemDraw",
            "bytes": 61230, "size_px": [1320, 470], "dpi": [300, 300], "analysis": {"summary": "Content covers 7% of the 1320x470 frame ...", "...": 0}}}
```

Headless the PNG is an RDKit `MolsToGridImage` (300x250 px per cell, 4 per row, labels as legends).
Stereo wedges follow RDKit's `WedgeMolBonds`; charges, isotopes, radicals and atom-map numbers are preserved.
Molecules with two or fewer heavy atoms (water, methanol, HCl) or a cumulated carbon (CO2, ketenes) get every
atom labelled (element, hydrogens and charge), as in `draw_reaction`.

### `draw_reaction`

One reaction on one line: reactants (with `+`), arrow with reagents above and conditions below, products.
*headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `reaction` | str | null | reaction SMILES `A.B>agents>C.D`, or arrow text `A + B -> C`, `A + B --[reagents]--> C` (see below) |
| `reactants`, `products` | list[str] | null | alternative to `reaction`: component lists (structure inputs, or dicts `{"structure": ..., "label": ..., "color": ..., "label_atoms": true}`) |
| `reagents_above` | list[str] | null | text lines above the arrow; lines that look like formulas (`NaBH4`, `H3O+`) get ChemDraw's automatic sub/superscripts |
| `conditions_below` | list[str] | null | text lines below the arrow |
| `arrow` | str | `"full"` | `full` (`forward`, `plain`), `equilibrium`, `unbalanced_equilibrium`, `double` / `resonance`, `retrosynthetic` / `retro`, `dashed`, `bold`, `hollow`, `nogo` / `cross`, `hash`, `line`, `half_left`, `half_right` |
| `reactant_labels`, `product_labels` | list[str or null] | null | captions under the components, e.g. `["85%", null]` |
| `atom_map` | bool | false | record atom mapping from map numbers in the SMILES (`[CH3:1]`) as `ReactionStepAtomMap` |
| + shared parameters (`name` defaults to `"reaction"`) | | | |

Reaction text syntax (`chem.parse_reaction` after names are resolved):
* arrows: `->`, `-->`, `>>`, `=>`, `→`, `⟶`, `⇒`, `⇌`, `⇄`, `↔`, ` = `, ` > `;
* `--[text]-->`: the bracket text is placed **above** the arrow (or appended below when `reagents_above` is
  also given);
* components separated by ` + ` (a `+` inside `[...]` such as `[NH4+]` is safe); a leading integer is a
  stoichiometric coefficient and duplicates the molecule (`2 H2O`);
* reaction SMILES agents (`A>agent>B`) are drawn as text above the arrow: the molecule's name when it came
  from a name, otherwise its formula;
* a single step only (`A -> B -> C` is an error: use `draw_scheme`).

Extra result key `reaction_smiles` (string, from the `<step>` written: `reactants>agents>products`).
Returns an image.

```json
{"reaction": "acetic acid + ethanol --[H2SO4 (cat.)]--> ethyl acetate + water",
 "conditions_below": ["reflux, 4 h"], "product_labels": ["85%", null], "style": "acs"}
```
```json
{"cdxml_path": "/Users/me/.chemdrawmcp/renders/reaction_20261107_154217_512_03.cdxml",
 "layout": {"style": "ACS Document 1996", "bond_length": 14.4, "objects": ["..."], "steps": 1, "bounds": {"...": 0}},
 "document": "reaction_20261107_154217_512_03.cdxml", "opened": true, "window_id": 812,
 "atoms": 11, "molecules": 4, "captions": 4, "arrows": 1,
 "reaction_smiles": "CC(=O)O.CCO>>CCOC(C)=O.O",
 "render": {"path": "/Users/me/.chemdrawmcp/renders/reaction_20261107_154218_061_04.png", "format": "png", "rendered_by": "ChemDraw",
            "bytes": 48213, "size_px": [1650, 420], "dpi": [300, 300], "analysis": {"summary": "...", "...": 0}}}
```

Equivalent with lists:

```json
{"reactants": ["CC(=O)O", "OCC"], "products": ["CC(=O)OCC", "O"],
 "reagents_above": ["H2SO4 (cat.)"], "conditions_below": ["reflux, 4 h"], "product_labels": ["85%", null]}
```

Notes: the arrow is at least 3 bond lengths long and grows to fit the widest text line; the `+` signs are
written with `InterpretChemically="no"` (otherwise ChemDraw 23 strips the formal charges of every atom in the
document); the `<scheme>/<step>` record makes ChemDraw treat it as a reaction (so `save_document(format="rxn")`
and `run_command("mapReactionAtoms")` work).  Tiny molecules that are unreadable as skeletal drawings (two or
fewer heavy atoms such as water, methanol, HCl, or cumulated carbons such as CO2 and ketenes) are drawn with
all atom labels automatically (element, hydrogens and charge, e.g. `OH-`, `NH3+`); a component dict can force this either way with `"label_atoms": true|false`.
Headless the PNG is RDKit's `DrawReaction` without the text.

### `draw_scheme`

Multi-step synthesis laid out left to right, chained through intermediates and wrapped into rows.
*headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `steps` | list[object] | required | one object per step: `{"reaction": "A --[Br2, FeBr3]--> B"}` or `{"reactants": [...], "products": [...], "above": [...], "below": [...], "arrow": "full", "reactant_labels": [...], "product_labels": [...], "atom_map": false}` (`reagents_above` / `conditions_below` are accepted as synonyms of `above` / `below`) |
| `title` | str | null | bold caption above the scheme |
| `max_width` | float | page width − 72 pt (min 300) | row width in points before wrapping |
| + shared parameters (`name` defaults to `"scheme"`) | | | |

Chaining: when the first reactant of step *i+1* is the same molecule (canonical SMILES) as the last product of
step *i*, the intermediate is drawn once and reused; further reactants of that step are added with `+`.  When
a step does not fit in the row, a new row starts and the intermediate is drawn again.  An empty `steps` list
is an error.

Extra result key `reaction_smiles`: a list with one string per step.  Returns an image (ChemDraw export when
live; headless only an RDKit grid of all components, `rendered_by` says *components only* and the `render`
has no `analysis`).

```json
{"title": "Scheme 1. Benzoic acid from benzene",
 "steps": [{"reaction": "benzene --[Br2, FeBr3]--> bromobenzene", "product_labels": ["2"]},
           {"reactants": ["bromobenzene"], "products": ["benzoic acid"], "above": ["1. Mg, THF", "2. CO2"], "below": ["3. H3O+"], "product_labels": ["3, 72%"]}],
 "style": "acs"}
```
```json
{"cdxml_path": "/Users/me/.chemdrawmcp/renders/scheme_20261107_160530_007_09.cdxml",
 "layout": {"style": "ACS Document 1996", "bond_length": 14.4, "objects": ["..."], "steps": 2, "bounds": {"...": 0}},
 "document": "scheme_20261107_160530_007_09.cdxml", "opened": true, "window_id": 821, "atoms": 22, "molecules": 3, "captions": 7, "arrows": 2,
 "reaction_smiles": ["c1ccccc1>>Brc1ccccc1", "Brc1ccccc1>>O=C(O)c1ccccc1"],
 "render": {"path": "/Users/me/.chemdrawmcp/renders/scheme_20261107_160531_220_10.png", "format": "png", "rendered_by": "ChemDraw", "...": 0}}
```

### `draw_mechanism`

A reaction mechanism: stages (starting materials, intermediates, transition states, resonance structures)
laid out left to right with curved electron-pushing arrows inside each stage and a reaction arrow to the
next stage, wrapping into rows.  *headless OK* (the CDXML carries the curved arrows; the headless PNG does not).

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `stages` | list[object] | required | one object per stage, see below |
| `title` | str | null | bold caption above the mechanism |
| `max_width` | float | page width − 72 pt (min 300) | row width in points before wrapping |
| + shared parameters (`name` defaults to `"mechanism"`) | | | |

A stage object: `{"structures": ["[OH-]", "CCl"]}` (or `"structure": "CC=O"` for one; structure inputs,
names resolved like everywhere else; a list entry may also be `{"structure": ..., "label": ...}`),
`"arrows": [...]` (electron-pushing arrows, see the grammar), `"label"` (caption below the stage), `"arrow"`
(kind of the reaction arrow to the next stage: `full` default, `equilibrium`, `resonance`, ... as in
`draw_reaction`), `"above"` / `"below"` (text on that arrow).  Several structures in one stage are drawn side
by side with `+`.

Arrow specs address atoms by **0-based RDKit index in SMILES order** (heavy atoms; explicit `[H]` atoms
written in the SMILES are kept so protons can be addressed), qualified as `fragment:atom` when a stage has
several structures:

| spec | meaning |
|------|---------|
| `lp:5->5-8` | lone pair of atom 5 forms the new bond 5–8 |
| `3-4->4` | bond 3–4 collapses onto atom 4 (lone pair / charge) |
| `3=4->4-7` | pi bond 3=4 becomes the new bond 4–7 |
| `3-4->3=4` | electrons move into the existing bond 3–4 (new pi bond) |
| `h:3->3-4` | the C–H bond of an implicit hydrogen on atom 3 as the source |
| `lp:5->5-8:1e` or `lp:5~>5-8` | one electron (fishhook head) |
| `lp:0:0->0:0-1:0` | lone pair of atom 0 of structure 0 attacks atom 0 of structure 1 |

Dict form: `{"from": "lp:5", "to": "5-8", "electrons": 1, "bulge": "up" | "down" | "auto", "label": "..."}`.
The arc side is chosen away from the molecule's centroid (`bulge: "auto"`) and the ends keep clear of atom
labels and neighbouring atoms.

Returns the shared drawing keys plus `mechanism`:
`{"name", "style", "bond_length", "title", "stages": [{"index", "label", "fragment_ids", "smiles", "row",
"electron_arrow_ids", "arrows": [{"id", "spec", "electrons", "kind", "source": {"kind", "atoms"}, "target": {"kind", "atoms"},
"tail", "head", "apex", "angle", "bulge", "label"}], "arrow_to_next", "above", "below", "reaction_arrow_id", "caption_id"}],
"electron_arrows": <count>, "reaction_arrows": [{"id", "kind", "from_stage", "to_stage", "tail", "head"}], "bounds", "warnings"?}`.
In `layout.objects` the curved arrows have `kind: "curved"`.  Returns an image (ChemDraw export when live;
headless only an RDKit grid of the structures, `rendered_by` says *structures only* and the `render` has no
`analysis`).

```json
{"title": "SN2: hydroxide + chloromethane",
 "stages": [{"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"], "label": "nucleophilic attack", "arrow": "full"},
            {"structures": ["CO", "[Cl-]"], "label": "products"}]}
```
```json
{"cdxml_path": "/Users/me/.chemdrawmcp/renders/mechanism_20261107_161849_153_02.cdxml",
 "layout": {"style": "ACS Document 1996", "bond_length": 14.4, "objects": ["..."], "steps": 0, "bounds": {"...": 0}},
 "document": "mechanism_20261107_161849_153_02.cdxml", "opened": true, "window_id": 830, "atoms": 7, "molecules": 4, "captions": 5, "arrows": 3,
 "mechanism": {"name": "mechanism", "style": "ACS Document 1996", "bond_length": 14.4, "title": "SN2: hydroxide + chloromethane",
               "stages": [{"index": 0, "label": "nucleophilic attack", "fragment_ids": [2, 5], "smiles": ["[OH-]", "CCl"], "row": 0,
                           "electron_arrow_ids": [9, 10],
                           "arrows": [{"id": 9, "spec": "lp:0:0->new:0:0-1:0", "electrons": 2, "kind": "full",
                                       "source": {"kind": "lp", "atoms": [[0, 0]]}, "target": {"kind": "new_bond", "atoms": [[0, 0], [1, 0]]},
                                       "tail": [57.06, 83.69], "head": [76.89, 84.65], "apex": [67.18, 79.88], "angle": -93.64, "bulge": "auto", "label": null},
                                      {"id": 10, "spec": "1:0-1:1->1:1", "...": 0}],
                           "arrow_to_next": "full", "above": [], "below": [], "reaction_arrow_id": 12, "caption_id": 11},
                          {"index": 1, "label": "products", "smiles": ["CO", "[Cl-]"], "arrows": [], "...": 0}],
               "electron_arrows": 2, "reaction_arrows": [{"id": 12, "kind": "full", "from_stage": 0, "to_stage": 1, "tail": [120.65, 84.65], "head": [163.85, 84.65]}],
               "bounds": {"...": 0}},
 "render": {"path": "...png", "format": "png", "rendered_by": "ChemDraw", "...": 0}}
```

Errors: `"stages is empty"`; a spec that points outside a structure (`"arrow lp:9->1 references atom 9 of
fragment 0 which has 3 atoms (0-based indices)"`); a bond spec for atoms that are not bonded; an
unresolvable name.  Check `mechanism.stages[*].arrows[*].tail/head` against the fragment's
`atom_positions` in `layout.objects` when an arrow looks misplaced in the thumbnail.

### `named_reactions`

Browse the library of 100 named-reaction templates (`chemdrawmcp/reactions.py`: Suzuki, Grignard, Wittig,
amide coupling, Boc protection/deprotection, CuAAC, Diels–Alder, aldol, oxidations, reductions, ...).
*headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `query` | str | null | free text matched against key, name, class, roles, aliases and description |
| `reaction_class` | str | null | class filter (substring): `C–C coupling`, `C–N coupling`, `acylation`, `addition`, `aromatic substitution`, `condensation`, `cycloaddition`, `deprotection`, `elimination`, `heterocycle synthesis`, `hydrolysis`, `multicomponent`, `olefination`, `oxidation`, `protection`, `rearrangement`, `reduction`, `substitution`, `C–C bond formation` |
| `name` | str | null | one template by key, display name or alias (fuzzy: `"suzuki"`, `"Suzuki-Miyaura"`, `"sm coupling"`) -> full description |

Returns `{"count", "classes": [...], "reactions": [{"key", "name", "class", "roles", "n_reactants", "equivalents",
"reagents", "conditions", "byproducts", "description", "aliases", "smarts", "example": {"reactants", "products"}}]}`,
or with `name` `{"template": {...same keys...}, "description": "<multi-line text>"}`.  An unknown name is an
error listing the closest candidates.

```json
{"name": "suzuki"}
```

### `apply_named_reaction`

Predict the product(s) of a named reaction for concrete reactants and, by default, draw the atom-mapped
reaction with the template's typical reagents and conditions.  *headless OK* (prediction is RDKit; the
drawing behaves like `draw_reaction`).

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `reaction` | str | required | template key / name / alias |
| `reactants` | list[str] | required | one structure input per template role, in any order (names are resolved) |
| `draw` | bool | true | also draw the reaction |
| `product_index` | int | 0 | which distinct outcome to draw when the template gives several (regio-isomers) |
| `style`, `document`, `open_in_chemdraw`, `render`, `return_image` | | as for `draw_reaction` | |
| `atom_map` | bool | true | write the atom mapping into the `<step>` (`ReactionStepAtomMap`) |
| `max_products` | int | 8 | maximum number of distinct outcomes returned |

Returns `{"template", "name", "outcomes": [{"template", "name", "class", "reactant_order", "reactants", "roles": {role: SMILES},
"products": [SMILES], "mapped": "<atom-mapped reaction SMILES>", "reagents", "conditions", "byproducts"}], "drawn"?: {...}}`
where `drawn` carries the usual drawing keys (`cdxml_path`, `layout`, `document`/`window_id`, `render`) and the
result includes the thumbnail.  By-products (H2O, HX, Ph3P=O, ...) are text only and are not drawn.
Errors say which role each reactant failed to match, or that the reactant count is wrong.

```json
{"reaction": "suzuki", "reactants": ["bromobenzene", "OB(O)c1ccccc1"]}
```
```json
{"template": "suzuki_coupling", "name": "Suzuki–Miyaura coupling",
 "outcomes": [{"template": "suzuki_coupling", "name": "Suzuki–Miyaura coupling", "class": "C–C coupling",
               "reactant_order": [0, 1], "reactants": ["Brc1ccccc1", "OB(O)c1ccccc1"],
               "roles": {"aryl/vinyl halide": "Brc1ccccc1", "boronic acid/ester": "OB(O)c1ccccc1"},
               "products": ["c1ccc(-c2ccccc2)cc1"],
               "mapped": "Br[c:1]1[cH:3][cH:5][cH:7][cH:6][cH:4]1.OB(O)[c:2]1[cH:8][cH:10][cH:12][cH:11][cH:9]1>>[c:1]1(-[c:2]2[cH:8][cH:10][cH:12][cH:11][cH:9]2)[cH:3][cH:5][cH:7][cH:6][cH:4]1",
               "reagents": ["Pd(PPh3)4 (5 mol%)", "K2CO3"], "conditions": ["dioxane/H2O, 80 °C"], "byproducts": ["B(OH)2X"]}],
 "drawn": {"cdxml_path": "...suzuki_coupling_<stamp>.cdxml", "document": "...", "opened": true, "window_id": 833, "render": {"...": 0}}}
```

### `suggest_reactions`

Which named reactions a set of reactants can undergo together.  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `reactants` | list[str] | required | structure inputs (names resolved); every template with the same number of roles is tried in every order |

Returns `{"reactants": [resolved SMILES], "suggestions": [{"template", "name", "class", "roles": {role: SMILES},
"products": [first outcome], "mapped", "n_outcomes", "reagents", "conditions"}]}` in registry order; an empty
list means no template matched.

```json
{"reactants": ["Brc1ccccc1", "OB(O)c1ccccc1"]}
```

### `insert`

Paste chemistry into an open document through ChemDraw's own importers.  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `data` | str | required | the text to paste |
| `kind` | str | `"auto"` | `smiles`, `inchi`, `mol`, `cdxml`, `name`, `helm`, `sequence`, `peptide`, `dna`, `rna`, `text` (plain caption); `auto` detects SMILES / InChI / molblock / SDF / CDXML / name and otherwise pastes as text |
| `document` | str | null | target document |

Returns (macOS) `{"document", "before": {"atoms", "molecules", "captions"}, "after": {...}, "selection": {"smiles", "bounds"} | null, "atoms_added", "kind", "window_id", "warning"?}`;
(Windows) `{"document", "kind", "mime", "mode" ("added" | "replaced+merged" | "nothing"), "atoms_before", "atoms_after", "atoms_added", "bonds_before", "bonds_after", "molecules", "captions", "alerts", "helm"?, "note"?, "warning"?, "fallbacks"?}`
— `replaced+merged` means `Objects.SetData` replaced the document and the previous content was merged back in
(the new objects are placed below the old ones).

The inserted objects stay selected (so `selection`, `clean_structure(select_all=false)`, `run_command` act on
them, and the user can drag them).  `auto` never picks the sequence kinds: pass `kind="peptide"` etc.
explicitly.  Reaction SMILES is rejected ("use draw_reaction for reactions").  On macOS the data goes through
the general pasteboard inside one JXA script (NSPasteboard: the user's clipboard is snapshotted with every item
and flavour -> the text is written -> `pasteAsSMILES` / `pasteAsCdxml` / `pasteAsStructureToName` ... -> the
snapshot is restored, whether or not the paste succeeded); `atoms_added <= 0` for a chemistry kind adds a `warning` (and any alert text).

```json
{"data": "Ala-Gly-Ser-Phe", "kind": "peptide", "document": 815}
```

### `add_caption`

Add a text caption (pasted as CDXML; ChemDraw chooses the position, the caption stays selected).  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `text` | str | required | caption text (`\n` for several lines) |
| `document` | str | null | target |
| `size` | float | style caption size | point size |
| `bold`, `italic` | bool | false | faces |
| `formula` | bool | false | ChemDraw formula face: automatic subscripts/superscripts in `H2SO4`, `Na+` |
| `color` | str | null | colour name or `#rrggbb` |

Returns `{"ok": true, "document": "...", "window_id": 815, "captions_after": 5}`.  The font comes from the default style
(`CHEMDRAWMCP_STYLE`, i.e. Arial for `acs`).

### `add_arrow`

Add a reaction arrow with optional text above/below (pasted; stays selected).  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `kind` | str | `"full"` | any arrow kind listed under `draw_reaction` |
| `length` | float | 3 x bond length (43.2 pt for `acs`) | arrow length in points |
| `document` | str | null | target |
| `text_above`, `text_below` | str | null | reagent / condition text (above uses the formula face) |

Returns `{"ok": true, "document": "...", "window_id": 815, "arrows_after": 2}`.

---

## Rendering and looking

### `render`

Export a document as an image and describe it numerically.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document |
| `format` | str | `"png"` | `png`, `svg`, `pdf`, `eps`, `tiff`, `gif`, `jpeg`, `bmp` (any `save_document` format is accepted) |
| `path` | str | renders dir | output file |
| `width` | int | null | resample the raster file to this width in pixels (in place, LANCZOS) |
| `return_image` | bool | true | attach a thumbnail (raster formats only) |

Returns `{"path", "format", "rendered_by": "ChemDraw", "bytes", "document": "<name>"}` and for raster formats
also `size_px`, `dpi`, `analysis`, `resized_to_width` (when `width` was given) — or `analysis_error`.
Returns an image for raster formats; vector formats give JSON only.

ChemDraw writes PNG at 300 dpi of the page content, so a one-line reaction is typically 1500-2500 px wide;
`width` shrinks the file afterwards (the dpi metadata is then lost).  SVG/PDF/EPS are true vector exports
and are what to hand to a user for a manuscript.

```json
{"document": 812, "format": "png", "width": 1200}
```

### `quick_render`

Render one structure or reaction to a file without leaving a document open.  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `input` | str | required | structure input, reaction SMILES, or `A + B -> C` text |
| `style` | str | `"acs"` | style |
| `format` | str | `"png"` | any `render` format when live; `png` or `svg` headless |
| `path` | str | renders dir | output file |
| `label` | str | null | caption under a single structure |
| `return_image` | bool | true | thumbnail |

Live: opens the generated CDXML as a temporary document, exports, closes it; returns `{"path", "format",
"rendered_by": "ChemDraw", "bytes", "size_px", "dpi", "analysis", "cdxml_path"}` plus an image for raster
formats (the CDXML stays in the renders dir).  Headless: `{"path", "format": "png", "rendered_by": "RDKit (ChemDraw not available)", "analysis"}`
plus image, or `{"path", "format": "svg", "rendered_by": "RDKit"}`; other formats error.

```json
{"input": "CC(=O)O + OCC -> CC(=O)OCC + O", "format": "svg", "path": "~/Desktop/ester.svg"}
```

### `screenshot`

Screenshot of a document window as the user sees it (window chrome included).  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | window to capture (brought to the front first) |
| `return_image` | bool | true | thumbnail |

Returns `{"path": ".../screenshot_<stamp>.png", "analysis": {...}}` plus the image.  macOS uses
`screencapture -l <window id>` and needs Screen Recording permission for the terminal/IDE running the server;
Windows grabs the main window rectangle (the window must not be minimised).  Prefer `render` for judging the
drawing itself; use `screenshot` to see the magnification, selection handles, rulers or a dialog.

### `analyze_image`

Numerical description of any PNG/JPEG on disk.  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `path` | str | required | image file |

Returns the full `analyze_render` dictionary (see Conventions) including `ascii` and `summary`.

### `compare_renders`

Pixel diff of two images (before/after a change).  *headless OK.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `path_a`, `path_b` | str | required | image files (`b` is resized to `a` when the sizes differ) |

Returns `{"width", "height", "resized_b", "changed_fraction", "changed_pixels", "mean_abs_diff", "bbox" (of the
changed region or null), "ascii", "summary"}`; identical images give `changed_fraction: 0` and
*"Images are identical within tolerance."* (tolerance 16/255 per channel).

---

## Chemistry

### `name_to_structure`

Chemical name -> structure.  *headless OK* (built-in dictionary / OPSIN / PubChem); names the offline
resolvers do not know are sent to ChemDraw's Name=Struct when ChemDraw is installed (same order as for
*Structure inputs* above: cache -> dictionary and OPSIN -> ChemDraw -> PubChem).

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `name` | str | required | the name |
| `document` | str | null | with `draw=true`: paste into this document instead of a new one |
| `draw` | bool | false | also draw it (default style, labelled with the name) |

Returns `{"name", "smiles", "properties": {basic RDKit properties: formula, molecular_weight, exact_mass,
canonical_smiles, isomeric_smiles, inchi, inchikey, heavy_atoms, atoms, bonds, charge, elemental_analysis,
rings, aromatic_rings, stereocenters, double_bond_stereo}, "source": "offline dictionary" | "OPSIN" | "ChemDraw Name=Struct" | "PubChem",
"drawn"?: {...}}`.  `source` is `"offline dictionary"` when the built-in dictionary answered, `"OPSIN"` when OPSIN did (ChemDraw
is not consulted then), `"ChemDraw Name=Struct"` when ChemDraw resolved the name, and `"PubChem"` when
PubChem (`CHEMDRAWMCP_ONLINE=1`) answered after ChemDraw failed or was unavailable.  `drawn` holds only the
document keys of a drawing result: `cdxml_path`, `layout`, `layout_warnings`? and either `document`,
`window_id`, `opened: true` (plus `atoms`, `molecules`, `captions`, `arrows` of the new document), or
`document`, `window_id`, `inserted: true`, `atoms_added`, `note`, or `opened: false`, `note` (headless) — no
`render` key and no thumbnail.

```json
{"name": "4-(dimethylamino)pyridine"}
```
```json
{"name": "4-(dimethylamino)pyridine", "smiles": "CN(C)c1ccncc1",
 "properties": {"formula": "C7H10N2", "molecular_weight": 122.171, "exact_mass": 122.08440, "inchikey": "VHYFNPMBLIVWCW-UHFFFAOYSA-N", "...": 0},
 "source": "offline dictionary"}
```

Results are cached per server process.  The built-in dictionary only knows common names; systematic IUPAC
names need OPSIN (offline) or ChemDraw's Name=Struct, which resolves most of them but takes a few seconds per
lookup.

### `structure_to_name`

IUPAC-style name from ChemDraw's Struct=Name.  *live.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `structure` | str | null | structure input; named in a temporary document that is closed afterwards |
| `document` | str | null | without `structure`: name the selection of this document (everything if nothing is selected) |

Returns `{"ok": <a name appeared>, "name", "smiles", "all_new_text", "alerts"}` and, in document mode,
`"note": "the name was added as a caption below the structure"`.  The tool runs the
`convertStructureToName` command and diffs the document's captions before and after, so in document mode the
caption **stays in the user's document** (undo with `run_command("undo")` if unwanted).  `ok: false` with
`alerts` means ChemDraw refused (unsupported structure, or Struct=Name not licensed).

```json
{"structure": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"}
```
```json
{"ok": true, "name": "(S)-2-(4-isobutylphenyl)propanoic acid", "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
 "all_new_text": ["(S)-2-(4-isobutylphenyl)propanoic acid"], "alerts": []}
```

### `properties`

Chemical properties from RDKit and/or ChemDraw.  *headless OK* with `structure`.

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `structure` | str | null | structure input -> RDKit properties |
| `document` | str | null | without `structure`: ChemDraw's values for the selection (whole document if nothing is selected) plus RDKit values from its SMILES |
| `extended` | bool | true | add `logp`, `tpsa`, `hbd`, `hba`, `rotatable_bonds`, `fraction_csp3`, `heavy_atom_count`, `lipinski_violations`, `qed` |
| `use_chemdraw` | bool | false | with `structure`: also compute ChemDraw's values in a temporary document (slower) |

Returns `{"rdkit": {...}, "chemdraw"?: {...}}` or `{"chemdraw": {...}, "rdkit"?: {...}, "rdkit_error"?}`.
`chemdraw` is the backend's selection info: macOS `{"document", "smiles", "formula", "molecular_weight",
"exact_mass", "elemental_analysis", "text", "bounds", "atoms", "bonds", "molecules", "captions", "arrows", "empty"}`
(e.g. `"elemental_analysis": "C, 60.00; H, 4.48; O, 35.52"`); Windows adds `selected`, `scope` and
`unavailable`.

```json
{"structure": "aspirin"}
```
```json
{"rdkit": {"formula": "C9H8O4", "molecular_weight": 180.159, "exact_mass": 180.04226, "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
           "isomeric_smiles": "CC(=O)Oc1ccccc1C(=O)O", "inchi": "InChI=1S/C9H8O4/...", "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
           "heavy_atoms": 13, "atoms": 21, "bonds": 13, "charge": 0, "elemental_analysis": "C, 60.00; H, 4.48; O, 35.52",
           "rings": 1, "aromatic_rings": 1, "stereocenters": {"assigned": 0, "unassigned": 0, "atoms": [], "labels": {}},
           "double_bond_stereo": {"assigned": 0, "unassigned": 0}, "name": "aspirin",
           "logp": 1.31, "tpsa": 63.6, "hbd": 1, "hba": 3, "rotatable_bonds": 2, "fraction_csp3": 0.111,
           "heavy_atom_count": 13, "lipinski_violations": 0, "qed": 0.55}}
```

### `analyze_reaction`

Balance check and bookkeeping for a single reaction.  *headless OK* (RDKit only).

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `reaction` | str | required | reaction SMILES or `A + B -> C` text (SMILES or names) |

Returns `{"reactants" / "agents" / "products": [{"formula", "molecular_weight", "smiles", "name"}],
"reaction_smiles", "conditions", "name", "element_balance": {"C": [4, 4], ...}, "unbalanced_elements",
"charge_balance": [0, 0], "balanced", "mass_reactants", "mass_products", "atom_economy"}` — atom economy is
MW(first product) / MW(all reactants) in percent and `null` when the reaction is not balanced.

```json
{"reaction": "acetic acid + ethanol -> ethyl acetate + water"}
```
```json
{"reactants": [{"formula": "C2H4O2", "molecular_weight": 60.052, "smiles": "CC(=O)O", "name": "acetic acid"}, {"formula": "C2H6O", "molecular_weight": 46.069, "smiles": "CCO", "name": "ethanol"}],
 "agents": [], "products": [{"formula": "C4H8O2", "molecular_weight": 88.106, "smiles": "CCOC(C)=O", "name": "ethyl acetate"}, {"formula": "H2O", "molecular_weight": 18.015, "smiles": "O", "name": "water"}],
 "reaction_smiles": "CC(=O)O.CCO>>CCOC(C)=O.O", "conditions": "", "name": "",
 "element_balance": {"C": [4, 4], "H": [10, 10], "O": [3, 3]}, "unbalanced_elements": [], "charge_balance": [0, 0], "balanced": true,
 "mass_reactants": 106.121, "mass_products": 106.121, "atom_economy": 83.0}
```

### `convert`

Convert between chemistry formats.  *headless OK* for the RDKit formats.

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `data` | str | required | SMILES, InChI, molblock/SDF, CDXML, a file path, or a name |
| `to` | str | required | `smiles`, `inchi`, `inchikey`, `mol`, `molv3000`, `sdf`, `cdxml`, `cdx`, `rxn`, `cml`, `helm`, `sln`, `name`, `png`, `svg`, `pdf`, `eps`, `tiff` (plus the other `save_document` formats) |
| `from_format` | str | auto-detected | force the input format (`smiles`, `inchi`, `molblock`, `sdf`, `cdxml`, `name`; aliases such as `smi`, `mol`, `sd` are accepted) |
| `path` | str | renders dir | where to write file formats (also honoured for text formats on the RDKit path) |
| `use_chemdraw` | bool | true | let ChemDraw convert when it is available; `false` forces RDKit |

Text results: `{"ok": true, "to": "inchikey", "result": "<text, trimmed to 8000 chars>", "converted_by": "ChemDraw" | "RDKit", "path"?}`.
File results: `{"ok": true, "to": "cdx", "path": "...", "bytes": 1234, "converted_by": "ChemDraw"}` (RDKit
file results have no `bytes`).  `to="name"` delegates to `structure_to_name`.

Which engine: with ChemDraw, the input is pasted into a temporary document (SMILES, InChI, molblock, CDXML
or name) and read back with `copyAs<Format>` (text formats `smiles`, `inchi`, `inchikey`, `mol`, `molv3000`,
`cdxml`, `helm`, `sln`) or exported with `save … as` (everything else).  Without ChemDraw RDKit handles
`smiles`, `inchi`, `inchikey`, `mol`, `molv3000`, `sdf`, `cdxml` (through the layout engine, default style),
`png`, `svg`; anything else fails with *"RDKit cannot write X; ChemDraw is needed for that format"*.
Reactions are not converted by this tool (reaction SMILES input goes to RDKit, which only handles single
structures): draw them with `draw_reaction` and `save_document(format="rxn")` instead.

```json
{"data": "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3", "to": "cdx", "path": "~/Desktop/caffeine.cdx"}
```

### `clean_structure`

ChemDraw's *Clean Up Structure* (or *Clean Up Reaction*).  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document |
| `reaction` | bool | false | run `cleanReaction` instead of `cleanStructure` |
| `select_all` | bool | true | select everything first; `false` cleans the current selection only |

Returns the `run_command` result: `{"ok", "command", "enabled_before", "document", "before": {"atoms", "molecules", "captions", "arrows"}, "after": {...}, "description", "alerts": []}`.

### `check_structure`

Chemical warnings.  *live, no autolaunch;* macOS only (the `checkStructure` command is not in the Windows
command table).

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document |

Returns `{"warnings": [{"object": "n" | "b" | "t" | "fragment", "warning": "..."}], "captions_not_chemical": <count>,
"check_structure": <alert text | error | "no dialog shown">, "ok": <no warnings>}`.  The warnings are the
`Warning="..."` attributes ChemDraw writes into the CDXML export (valence errors, unusual bonds, ...); captions
that ChemDraw merely could not interpret as chemistry are counted separately.  Then everything is selected and
ChemDraw's *Check Structure* command is run; its dialog is dismissed and reported.

### `predict_nmr`

1H or 13C NMR prediction with ChemDraw.  *live;* macOS only; needs a ChemDraw edition with NMR prediction.

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `structure` | str | null | structure input (temporary document) |
| `document` | str | null | without `structure`: predict for the selection of this document (everything if nothing is selected) |
| `nucleus` | str | `"1H"` | `"1H"` or `"13C"` (anything containing `13` means carbon) |
| `keep_open` | bool | true | leave the spectrum document open in ChemDraw |
| `return_image` | bool | true | thumbnail of the spectrum |

Runs `predict1H_NMR_Shifts` / `predict13C_NMR_Shifts`, waits for the new spectrum document, renders it to
PNG and collects its captions.  Returns `{"ok": true, "spectrum_document": DocumentInfo, "render": {...},
"text": [up to 200 caption strings: shift assignments, "Estimation quality ..." notes]}` plus the image, or
`{"ok": false, "error": "no spectrum document appeared", "command": {...}}`.

```json
{"structure": "CC(=O)OCC", "nucleus": "13C", "keep_open": false}
```

---

## Selection, commands, menus

### `selection`

What is selected.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `document` | str | null | document |
| `as_format` | str | `"info"` | `info`, or a conversion: `smiles`, `inchi`, `inchikey`, `mol`, `molv3000`, `cdxml`, `helm`, `sln` (macOS also `helm_natural`; Windows also `name`) |

`info` returns the selection-info dictionary shown under `properties`; a conversion returns
`{"format": "smiles", "result": "..."}`.  Conversions use ChemDraw's `copyAs…` commands, which are disabled
when nothing chemical is selected (error *"copyAsSMILES is not enabled: select one or more structures first"*).

### `select`

Change the selection.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `mode` | str | `"all"` | `all`, `none`, `next` (next structure), `previous`, `invert` |
| `document` | str | null | document |

Returns the selection info after the change.  `none` on macOS is implemented as select-all + invert; Windows
supports `all` and (if the COM document has a deselect member) `none`.

### `commands`

List ChemDraw's named commands.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `pattern` | str | null | case-insensitive regular expression on the command name (`"clean"`, `"copyAs"`, `"flip|rotate"`, `"^paste"`) |
| `enabled_only` | bool | false | only commands currently enabled for the front document |
| `limit` | int | 150 | maximum entries returned |

Returns `{"count": <matches>, "shown": <returned>, "commands": [{"name", "enabled", "short_description", "long_description"}]}`.
macOS reads `app.commands` (about 1600 entries); when descriptions or enabled flags are requested the match
list is capped at 400, so narrow the pattern.  Windows returns the curated `COMMANDS` table with `enabled: null`.

```json
{"pattern": "reaction|autonumber", "enabled_only": true}
```

### `run_command`

Execute a named command on a document.  *live, no autolaunch.*

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `command` | str | required | exact command name (case-sensitive on macOS), e.g. `cleanStructure`, `cleanReaction`, `mapReactionAtoms`, `autoNumberReaction`, `convertStructureToName`, `convertNameToStructure`, `toggleShowingAtomNumber`, `flipObjectsHorizontal`, `flipObjectsVertical`, `centerObjectOnPage`, `zoomToFit`, `analyzeStoichiometry`, `expandLabel`, `contractLabel`, `toggleKekulization`, `copyAsSMILES`, `selectAll`, `undo`, `redo` |
| `document` | str | null | document (made the key window first) |

Returns (macOS) `{"ok": true, "command", "enabled_before": true, "document", "before": {"atoms", "molecules", "captions", "arrows"}, "after": {...}, "description", "alerts": [...]}`
or `{"ok": false, "enabled_before": false, "error": "command is not enabled right now (<short>; <long>)", "alerts"}`;
(Windows) `{"ok": true, "command", "enabled_before": null, "alerts", "result"}`.

Commands act on the current selection of the key window: select first (`select`, or the objects left selected
by `insert`).  Any alert the command raises is dismissed with its preferred button and returned in `alerts`
(do not assume the command succeeded when `alerts` is non-empty).  Commands that open a dialog that is not a
simple alert (preferences, *Save As*) are dismissed too, so they achieve nothing.

```json
{"command": "mapReactionAtoms", "document": 812}
```

### `menus`

ChemDraw's menus and their item titles.  *live, no autolaunch;* macOS only.

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `menu` | str | null | an exact menu title returns only that menu's items; any other string lists the menu titles that start with it |

Returns `{"File": ["New Document", "Open…", ...], "Edit": [...], ...}`.  With `menu` set to an exact title
(case-sensitive, e.g. `"Structure"`) the result is `{"Structure": ["Clean Up Structure", ...]}`, that menu's
items.  With any other string the value is the list of **menu titles** that start with it, case-insensitively,
e.g. `menus(menu="struct")` -> `{"struct": ["Structure"]}`; call again with the exact title to get the items
(these are titles, not items, so do not pass them to `do_menu` as `item`).

### `do_menu`

Invoke a menu item by titles.  *live, no autolaunch;* macOS only.

| parameter | type | default | meaning |
|-----------|------|---------|---------|
| `menu` | str | required | menu title, e.g. `"Structure"` |
| `item` | str | required | item title exactly as listed by `menus`, e.g. `"Clean Up Structure"` |
| `document` | str | null | document to make key first |

Returns `{"ok": true, "menu", "item", "alerts"}` or `{"ok": false, "error": "menu item is disabled", "alerts"}`.
Items that open dialogs are dismissed like alerts.

### `dismiss_dialogs`

Close any modal alert or sheet ChemDraw is showing and return its text.  *live, no autolaunch.*

No parameters.  Returns `{"dismissed": ["No reaction center was found. [dismissed with 'OK']", ...]}` (empty
list when nothing was open).  macOS inspects ChemDraw's windows through System Events (sheets, `AXDialog`,
`AXSystemDialog`, unnamed floating windows) and clicks the first of *OK, Yes, Don't Save, Discard, Cancel,
Close, No*, else the first button; Windows walks `#32770` dialog windows.  Tools call this automatically when
a call times out or fails, so you rarely need it — `status` reports `pending_dialogs` when a dialog is
blocking the application (for example one the user opened).

---

## Utilities

### `styles`

The document styles the drawing tools know.  *headless OK.*

No parameters.  Returns `{"default": "acs", "styles": {"acs": {"name": "ACS Document 1996", "bond_length": 14.4,
"line_width": 0.6, "bold_width": 2, "label_font": "Arial", "label_size": 10, "caption_font": "Arial",
"caption_size": 10, "bond_spacing": 18, "page": [0, 0, 540, 720]}, "nature": {...}, ...}}`.

| key | stationery | bond | line | labels | captions | page (pt) |
|-----|-----------|------|------|--------|----------|-----------|
| `acs` (`jacs`, `acs1996`) | ACS Document 1996 | 14.4 | 0.6 | Arial 10 | Arial 10 | 540 x 720 |
| `nature` (`nature_research`) | Nature Research Document | 10.8 | 0.6 | Arial 6 | Arial 6 | 510 x 680 |
| `rsc` (`rsc_2col`) | RSC (2 Column) Document | 12.2 | 0.45 | Helvetica 7 | Helvetica 7 | 539 x 785 |
| `rsc_1col` | RSC (1 Column) Document | 12.2 | 0.45 | Helvetica 7 | Helvetica 7 | 235 x 646 |
| `wiley` | Wiley Document | 17 | 0.75 | Arial 12 | Arial 12 | 770 x 523 |
| `science_of_synthesis` (`science`) | Science of Synthesis | 17 | 0.8 | Arial 10 | Arial 10 | 561 x 770 |
| `j_het_chem` | J. Het. Chem | 10.8 | 0.6 | Times New Roman 7 | Times New Roman 7 | 540 x 720 |
| `j_chin_chem_soc` | J. Chin. Chem. Soc. Doc | 18 | 1 | Arial 12 | Arial 12 | 675 x 1170 |
| `j_mol_mod_1col` / `j_mol_mod_2col` | J. Mol. Mod. (1/2 Column) | 14.4 | 0.6 | Times New Roman 10 | Times New Roman 10 | 241 / 505 x 720 |
| `helv_chim_acta` (`hca`) | Verlag Helv. Chim. Acta Doc | 17 | 0.6 | Helvetica 10 | Helvetica 12 | 510 x 721 |
| `adv_synth_catal` (`advanced_synthesis`) | Adv. Synth. Catal. Document | 17 | 1 | Arial 10 | Times New Roman 12 | 369 x 721 |
| `phytomedicine` | Phytomedicine Document | 20 | 1 | Arial 12 | Arial 12 | 540 x 720 |
| `default` (`new`, `chemdraw`) | New Document | 30 | 1 | Helvetica 10 | Helvetica 12 | 540 x 720 |
| `a4` (`a4_document`) | New A4 Document | 30 | 1 | Helvetica 10 | Times New Roman 12 | 523 x 770 |
| `slide` (`presentation`) | New Slide | 30 | 1.6 | Helvetica 12 | Helvetica 16 | 1101 x 635, white on dark |
| `dark` (`dark_mode`) | Dark Mode | 14.4 | 1.61 | Arial 10 | Arial 10 | 1440 x 1080, light on black |
| `idraw` | I Draw Styles | 30 | 1 | Arial 10 | Times New Roman 12 | 540 x 720 |
| `terminal_carbons` | Show Terminal Carbon Labels | 30 | 1 | Helvetica 10 | Times New Roman 12 | 540 x 720, CH3 labels shown |

Style names, aliases and stationery titles are matched case-insensitively with spaces/dashes treated as
underscores.
