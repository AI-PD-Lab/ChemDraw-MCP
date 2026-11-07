# Architecture

## Two engines

```
                       ┌──────────────────────────────────────────────────────────┐
  MCP client           │ chemdrawmcp.server  (FastMCP "chemdraw", 43 tools)       │
  (Claude Code, Codex, │   parses inputs, picks the engine, writes files, returns  │
   Claude Desktop) ───►│   JSON (+ PNG thumbnail)                                 │
                       └───────┬───────────────────────────────┬──────────────────┘
                               │                               │
                 headless CDXML engine                  live ChemDraw backend
          ┌────────────────────┴──────────────┐   ┌────────────┴──────────────────────┐
          │ chem.py      RDKit: parse, 2D,    │   │ backends/base.py   abstract API    │
          │              properties, depict   │   │ backends/macos.py  JXA/AppleScript │
          │ layout.py    grids, reactions,    │   │ backends/windows.py COM (+ChemScript)│
          │              schemes, collisions  │   │ NullBackend        headless stub   │
          │ cdxml.py     CDXML writer/reader  │   └────────────┬──────────────────────┘
          │ cdx.py       binary CDX reader    │                │ documents, selection,
          │ textmetrics.py fonts, atom labels │                │ ~1600 commands, menus,
          │ abbreviations.py nickname table   │                │ save-as, pasteboard
          │ mechanism.py electron-pushing     │                ▼
          │ reactions.py named-reaction SMARTS│       ChemDraw application
          │ styles.py    stationery settings  │
          │ imageanalysis.py  render metrics  │
          └───────────────────────────────────┘
```

* **Live ChemDraw** does what only ChemDraw can do: journal-quality rendering, Name=Struct / Struct=Name,
  its chemical properties and elemental analysis, every export format (CDX, RXN, RDF, CML, CT, EPS, 3MF ...),
  NMR prediction, Clean Up, the named commands and menus.  The server never draws *inside* ChemDraw
  object-by-object; it always prepares a CDXML document and lets ChemDraw open or paste it.
* **The headless engine** never needs the application.  It turns SMILES/InChI/molblocks/names into RDKit
  molecules, lays out structures, reactions and schemes in a chosen style, writes CDXML that ChemDraw opens
  without complaint, reads CDXML back into SMILES/text/arrows/steps, and renders with RDKit when ChemDraw is
  absent.  It is what makes the server usable on Linux and in CI, and what makes every drawing deterministic.

`backends.get_backend()` chooses at first use: `CHEMDRAWMCP_BACKEND=none|headless|off` -> `NullBackend`;
`auto` -> `MacBackend` on Darwin (falls back to `NullBackend` if no app bundle is found), `WindowsBackend` on
Windows, `NullBackend` elsewhere; `macos`/`windows` force a backend.  `NullBackend` raises
`ChemDrawNotAvailable` (with the reason) from every method, which `server._err` turns into the JSON error plus
a hint about the headless tools.

## Data flow of a drawing call

```
input text ──► chem.detect_format ──► chem.mol_from_input ──► prepare_for_depiction
 (SMILES, InChI,   (smiles/inchi/molblock/   (RDKit Mol)        (sanitize, RemoveHs, 2D coords via CoordGen,
  molblock, CDXML,  sdf/cdxml/rxn/            name → _resolve_name  kekulize, wedge bonds from stereo)
  path, name)       reaction_smiles/name)      (dictionary → OPSIN   │
                                                → ChemDraw Name=Struct▼
                                                → PubChem)       layout.build_structures / build_reaction / build_scheme
                                                                  │  (cdxml.CDXMLDocument: fragments, text, arrows,
                                                                  │   plus signs, <scheme>/<step>, page tiling)
                                                                  ▼
                                               <renders>/<name>_<stamp>.cdxml   (always written)
                                                                  │
                     ┌────────────────────────────────────────────┴─────────────────────────┐
                     │ live                                                                 │ headless / open_in_chemdraw=false
                     ▼                                                                      ▼
   backend.open_document(path)  or  backend.insert(doc, xml, "cdxml")           chem.depict_png / Draw.MolsToGridImage
                     │                                                                      │
                     ▼                                                                      │
   backend.save_document(handle, <name>_<stamp>.png, "png")   (ChemDraw export, 300 dpi)    │
                     │                                                                      │
                     └──────────────────────────────► imageanalysis.analyze_render ◄────────┘
                                                      (coverage, bbox, margins, palette, contrast, ASCII map, summary)
                                                                  │
                                                                  ▼
                                   JSON (cdxml_path, layout, document/window_id, render{path, analysis}) + thumbnail
```

Reading goes the other way: `backend.document_cdxml()` (macOS: save-as CDXML into a temp dir; Windows:
`Objects.Data("chemical/x-cdxml")`) -> `cdxml.read_cdxml()` -> molecules as RDKit mols (`fragment_to_mol`)
-> SMILES/formula/MW, captions, arrows, graphics, steps -> reaction SMILES.  `read_file` feeds a `.cdxml`
from disk into the same reader, or first converts a binary `.cdx` with `cdx.cdx_to_cdxml()` (a tagged
little-endian stream parser verified against ChemDraw 23 files; object ids are identical in both formats).

Name resolution (`server._resolve_name`): per-process cache -> the built-in dictionary (`chem._COMMON`, ~3,100
structures under ~7,800 names and aliases, exposed as `chem.COMMON_NAMES`; `server._dictionary_lookup` tries the
case-sensitive table, then the normalised key with hydrate / purity / `L-`, `(S)-` qualifiers handled) -> OPSIN when
installed (`chem._opsin`: `py2opsin`, an `opsin` executable on PATH or `$OPSIN_JAR`) -> ChemDraw's Name=Struct when
ChemDraw is available (`insert(kind="name")` into a temporary document, then `selection_as("smiles")`; the only
step that can launch ChemDraw, a few seconds per lookup) -> PubChem PUG REST (`chem._pubchem`) only with
`CHEMDRAWMCP_ONLINE=1`.  `_NAME_SOURCE` remembers which engine answered (`"offline dictionary"`, `"OPSIN"`,
`"ChemDraw Name=Struct"`, `"PubChem"`); the library function `chem.resolve_name` used by the examples runs the
same chain without the ChemDraw step.

## Module map

| module | role |
|--------|------|
| `chemdrawmcp/server.py` | the MCP tools, argument normalisation (`_doc`, `_mol`, `_reaction_step`), output paths, materialisation in ChemDraw, rendering helpers, error envelope |
| `chemdrawmcp/mcpcompat.py` | imports `FastMCP`/`Image` from mcp 1.x (`mcp.server.fastmcp`) or mcp 2.x (`mcp.server.mcpserver.MCPServer`) |
| `chemdrawmcp/backends/base.py` | `ChemDrawBackend` ABC, `DocumentInfo`, `CommandInfo`, `EXPORT_FORMATS` (key -> extension, ChemDraw file-type name, MIME), `INSERT_KINDS`, `NullBackend`, exceptions |
| `chemdrawmcp/backends/macos.py` | `MacBackend`: JXA scripts over `osascript`, pasteboard round trips, Accessibility dialog handling, screenshots |
| `chemdrawmcp/backends/windows.py` | `WindowsBackend` (COM through pywin32/comtypes), `ChemScriptBridge`, Win32 dialog/screenshot helpers, editable candidate tables |
| `chemdrawmcp/cdxml.py` | `CDXMLDocument` writer (fragments with label-aware boxes, text, straight/curved arrows, plus signs, graphics, steps, page tiling) and reader (`read_cdxml`, `fragment_to_mol`, `cdxml_mols`, `document_text`) |
| `chemdrawmcp/cdx.py` | binary CDX reader: `parse_cdx` (object tree), `cdx_to_cdxml`, `read_cdx`, `dump_cdx` (debugging), `is_cdx`; the `OBJECT_TAGS` / `PROPERTY_TAGS` tables say which tags were verified on ChemDraw 23 files |
| `chemdrawmcp/textmetrics.py` | text widths/boxes from real fonts (Pillow) or an Arial advance table, formula-face runs (`formula_runs`), glyph centres, predicted atom labels and their boxes (`fragment_labels`, `atom_label`) |
| `chemdrawmcp/abbreviations.py` | `SUBSTITUENTS` (label -> SMILES with the attachment atom first: `OMe`, `Boc`, `TBS`, `Ph`, ...), `MOLECULES` (reagent/solvent abbreviations), `lookup`, `parse_atom_label` (`"NH3+"`, `"13C"`, `"D"`) — used by the CDXML reader for nicknames ChemDraw did not expand |
| `chemdrawmcp/mechanism.py` | electron-pushing arrows on placed fragments: `parse_arrow_spec` (`"lp:5->5-8"`, `"3=4->4-7"`, `"h:3->3-4"`, `:1e` fishhooks, `fragment:atom` addressing), `add_electron_arrows` (plans tail/head/arc around atom labels and neighbours), `Stage`, `build_mechanism` (stages joined by reaction/resonance arrows, wrapping rows, structures oriented towards their partner), `describe_mechanism`; exposed as the `draw_mechanism` tool |
| `chemdrawmcp/reactions.py` | registry of 100 named reactions (`REACTIONS`: `ReactionTemplate` with RDKit reaction SMARTS, class, reagents/conditions text, roles, aliases, by-products, worked example), `get_reaction` (fuzzy lookup), `list_reactions`, `apply_reaction` (products + atom-mapped reaction SMILES, any reactant order), `suggest_reactions`, `to_step` (a drawable `layout.Step`), `describe`; exposed as the `named_reactions`, `apply_named_reaction` and `suggest_reactions` tools |
| `chemdrawmcp/layout.py` | `Component`, `Step`, `build_structures`, `build_reaction`, `build_scheme`, `collisions`, `layout_warnings` |
| `chemdrawmcp/styles.py` | `STYLES` (19 stationery settings), `ALIASES`, `get_style`, `style_names` |
| `chemdrawmcp/chem.py` | RDKit layer: format detection, parsing (incl. reactions in SMILES/RXN/CDXML/arrow text), depiction preparation, properties, elemental analysis, reaction properties, PNG/SVG depiction, molblock/SDF writers, substructure/Tanimoto, name resolution |
| `chemdrawmcp/imageanalysis.py` | `analyze_render`, `summarize_render`, `diff_images`, `thumbnail_png`, ID-pass segmentation helpers |
| `chemdrawmcp/colornames.py` | nearest named colour for the palette report |
| `chemdrawmcp/discovery.py` | finding the ChemDraw app/exe and its version on macOS/Windows, `home_dir()`, `output_dir()` |
| `chemdrawmcp/install.py` | registration with Claude Code (`~/.claude.json`), Codex (`~/.codex/config.toml`) and Claude Desktop; `--check` |
| `install.py` (repo root) | bootstrap: venv (uv if present), `pip install -e .[windows]` on Windows, then `chemdrawmcp.install` |

`draw_mechanism` follows the same materialise/render path as the other drawing tools, with
`mechanism.build_mechanism` in place of `layout.build_*` and `describe_mechanism(doc)` added to the result.
`apply_named_reaction` runs the template's reaction SMARTS with RDKit (`reactions.apply_reaction`: every
reactant order is tried, products sanitised and deduplicated by canonical SMILES, grouped `(A.B)` product
templates split into fragments), then builds a `layout.Step` with `reactions.to_step` (atom-map numbers kept
so the `<step>` carries `ReactionStepAtomMap`; template reagents above, conditions below) and continues like
`draw_reaction`.

## CDXML conventions the engine relies on

All verified against files written by ChemDraw 23.1.

**Units and geometry.**  Points (1/72 inch); y grows downwards; atom positions `p="x y"`; every object has a
`BoundingBox="x0 y0 x1 y1"` and an increasing `Z`.  RDKit coordinates are scaled so the median bond length
equals the style's `bond_length` (and rescaled to 1.5 Å on the way back).  Fragment boxes are padded by
0.35 x bond length.

**Root element.**  `<CDXML CreationProgram="chemdrawmcp" Name=...>` with `BoundingBox` of the content,
`BondLength`, `LineWidth`, `BoldWidth`, `MarginWidth`, `HashSpacing`, `BondSpacing` (%), `ChainAngle`,
`LabelFont`/`LabelSize`/`LabelFace`, `CaptionFont`/`CaptionSize`/`CaptionFace`, `LabelJustification`,
`CaptionJustification`, `FractionalWidths`, `ShowTerminalCarbonLabels`, `HideImplicitHydrogens`,
`InterpretChemically="yes"`, `PrintMargins`, `color`/`bgcolor`, and autonumbering defaults
(`RxnAutonumberStyle="Roman" RxnAutonumberFormat="(#)"`).  `<fonttable>` ids start at 3 (`charset="iso-8859-1"`).

**Colours.**  Indices 0 (black) and 1 (white) are implicit; `<colortable>` entries are numbered from 2.  The
style tables keep ChemDraw's default palette (white, black, red, yellow, green, cyan, blue, magenta -> 2..9);
`CDXMLDocument.color_index` maps names/`#rrggbb`/RGB tuples onto existing entries or appends new ones.  The
`dark` and `slide` styles use `color`/`bgcolor` indices into their own tables (light text on dark pages).

**Atoms `<n>`.**  `Element` omitted for carbon; `NumHydrogens` written for heteroatoms and for any atom with
a charge, isotope or radical (ChemDraw otherwise recomputes valence); `Charge`, `Isotope`,
`Radical="Doublet|Singlet|Triplet"`, `AtomNumber` for atom-map numbers; dummy atoms become
`NodeType="GenericNickname" GenericNickname="R"`; with `label_atoms` every atom gets an explicit
`<t><s face="96">CH3</s></t>` label whose text includes the hydrogens and the formal charge (`OH-`, `NH3+`, `Fe3+`;
the formula face renders the trailing sign as a superscript).  `server._tiny` switches this on for molecules with at
most two heavy atoms or a cumulated carbon (water, methanol, HCl, CO2, ketenes), which are unreadable as skeletal
drawings.  Fragment boxes are **label-aware**: `textmetrics.fragment_labels`
predicts the label ChemDraw will draw for every labelled atom (`OH`, `H3N+`, `13C`; hydrogens first — `HO`,
`H2N` — when the neighbours lie to the right) and its box padded by the style's `margin_width`, and the
fragment's `BoundingBox` is the union of the atom box (padded by 0.35 L) and those label boxes.  The atom
positions and predicted labels are recorded on the fragment's `ObjRef` (`atom_positions`, `labels`) and
reported in the tools' `layout`; `mechanism.py` uses them to aim curved arrows.

**Bonds `<b>`.**  `B`/`E` node ids, `Order` 2 / 3 / 1.5 / dative / 4 / ionic (single omitted),
`Display="WedgeBegin|WedgedHashBegin|Wavy"` with the wedge starting at the stereocentre; with the style /
settings key `wavy_unspecified_double_bonds` a stereo-capable double bond of unspecified geometry is also
drawn wavy.  The reader also accepts `WedgeEnd`/`WedgedHashEnd` (flipped) and `Dash`/`Bold` are passed through.

**Text `<t>`.**  `p` is the baseline start (left edge for `Justification="Left"`, centre for `Center`),
`LineHeight="auto"`, runs `<s font= size= face= color=>`; face bits 1 bold, 2 italic, 4 underline, 32
subscript, 64 superscript, **96 = ChemDraw's formula face** (automatic sub/superscripts for digits and charges
— used for atom labels and for reagent lines that look like formulas).  Text widths come from
`textmetrics`: glyph advances measured once with Pillow when Arial, Helvetica, Liberation Sans, DejaVu Sans
or the requested font (e.g. Times New Roman) is installed, otherwise an Arial advance-width table
(`CHEMDRAWMCP_NO_FONTS=1` forces the table for reproducible output, `CHEMDRAWMCP_FONT_DIR` adds a font
directory); formula-face sub/superscript runs are measured at 60 % size; vertical metrics are 0.9 / 0.25 /
1.2 x size.  ChemDraw re-measures on open, so boxes stay approximate.
A caption consisting only of `+` is written with `InterpretChemically="no"`: ChemDraw interprets captions
chemically, and a lone `+` makes ChemDraw 23 strip the formal charges of every atom in the document.

**Arrows `<arrow>`.**  Straight: `Tail3D="x y 0"`, `Head3D="x y 0"`, `FillType="None"`,
`ArrowheadType="Solid"` (default), `HeadSize="1000"`, `ArrowheadCenterSize="875"`, `ArrowheadWidth="250"`,
then per kind: `ArrowheadHead="Full|HalfLeft|HalfRight"`, `ArrowheadTail`, `ArrowShaftSpacing="300"`
(equilibrium) / `"500"` (retrosynthetic, with `ArrowheadType="Angle"`), `ArrowheadType="Hollow"`,
`LineType="Dashed|Bold"`, `NoGo="Cross|Hash"`, `ArrowEquilibriumRatio`.  Curved (electron pushing):
`AngularSize` (signed sweep in degrees) plus `Center3D`, `MajorAxisEnd3D`, `MinorAxisEnd3D` of the circular
arc; `ArrowheadHead="HalfLeft"` for one-electron fishhooks.  Sign convention as verified in ChemDraw and used
by `mechanism.py`: a positive angle bulges to the *right* of the tail→head direction on screen (downwards for
a left-to-right arrow), a negative angle to the left (the older `add_curved_arrow` docstring says the
opposite; trust `mechanism.py`).  The reader also treats `<graphic GraphicType="Line" ArrowType=...>` as an
arrow (older files); ChemDraw itself writes such a legacy `<graphic>` twin next to every `<arrow>` it saves,
marked `SupersededBy="<arrow id>"`, and the reader skips those so arrows are not counted twice.

**Graphics `<graphic>`.**  `GraphicType="Rectangle|Oval|Line"` with `RectangleType` (`Plain`, `RoundEdge`,
`Dashed`, `Bold`, `Filled`, `Shaded`), `OvalType`, `LineType`, `CornerRadius`.

**Reaction semantics.**  `<scheme><step ReactionStepReactants="ids" ReactionStepProducts="ids"
ReactionStepArrows="ids" ReactionStepObjectsAboveArrow="ids" ReactionStepObjectsBelowArrow="ids"
ReactionStepPlusses="ids" ReactionStepAtomMap="r1 p1 r2 p2 ..."/></scheme>` at page level.  This is what makes
ChemDraw see a reaction (RXN export, autonumbering, stoichiometry) and what `read_document` turns into
reaction SMILES.

**Page tiling.**  `<page ... Width Height HeightPages WidthPages DrawingSpace="poster" PrintTrimMarks="yes">`.
`to_xml` computes `WidthPages = ceil((content.x1 + 36) / page_w)` and likewise for the height, so content that
is wider than the stationery page gets extra printer pages instead of triggering ChemDraw's "expand the
document?" dialog on open.  Content is laid out from the origin (36, 36) for grids and (36, 72) for reactions
and schemes.

**Layout constants** (`layout.py`, L = bond length): every component's label-aware box is centred on the row
baseline; `CAPTION_GAP` 0.35 L between a structure box and its caption; `PLUS_GAP` 0.6 L minimum clearance
on either side of a `+` (a caption 2 pt larger than the caption size, centred optically on the baseline using
the measured glyph centre); `ARROW_CLEARANCE` 0.4 L (plus half the shaft spacing for equilibrium /
retrosynthetic arrows) between the arrow line and reagent or condition text; arrow length max(3 L, widest
text line + `ARROW_TEXT_MARGIN` 0.8 L); `LINE_GAP` 0.2 x text size between stacked lines; grid gap 1.5 L,
reaction gap 1.0 L, scheme row gap 2.5 L; origins (36, 36) for grids and (36, 72) for reactions/schemes; row
heights account for captions, atom labels and arrow text; wrap at `max_width` (page width − 72, min 300);
chaining by canonical SMILES (memoised per build).  After every `build_*` call the document carries
`doc.warnings` from `layout.layout_warnings` (overlapping fragment / text / arrow-line boxes, content wider
than the page); `tests/test_layout_quality.py` keeps the layouts overlap-free and every drawing tool returns
the list as `layout_warnings` when it is non-empty (`server._materialize`).

**Reader.**  Walks pages and `<group>`s; each top-level `<fragment>` with atoms is one molecule.  Nested
fragments (nicknames, peptide residues) are expanded in place through their `ExternalConnectionPoint`s
(matched via `ConnectionOrder` / `BondOrdering`); nicknames without an inner fragment are expanded from
`abbreviations.py` when the label is known, otherwise they become dummy atoms with a `dummyLabel`; charge,
isotope, hydrogen count and radical fall back to the label text (`O-`, `NH3+`, `D`, `13C`); explicit hydrogen
nodes are folded into their heavy atoms after stereo perception; wedges and coordinates give tetrahedral and
E/Z stereo; `Order="1.5"` rings are kekulized when possible.  Problems that could not be repaired end up in
the molecule's `error` (best-effort SMILES still returned), softer notes in `warnings`.  Ids map step
references to molecules/texts; `Warning` attributes carry ChemDraw's chemical warnings; `<chemicalproperty>`
objects attach names/properties to the molecules they belong to (`BasisObjects`); bracketed groups,
spectra and unknown element types are counted rather than interpreted.

## macOS automation

Verified against ChemDraw 23.1.2, bundle id `com.revvity.ChemDraw` (older ids `com.perkinelmer.ChemDraw`,
`com.cambridgesoft.ChemDraw` are recognised).  The app is found by `CHEMDRAW_APP`/`CHEMDRAW_EXE`, then the
highest-versioned `/Applications/ChemDraw*.app` or `~/Applications/ChemDraw*.app` (companions such as *ChemDraw
Collections* rejected), then Spotlight by bundle id.

**JXA prelude.**  Every call runs `osascript -l JavaScript -` with a script made of a prelude plus the body of
`main()`.  The prelude defines `app = Application(BUNDLE)`, `ARGS` (the Python kwargs as JSON), `cmd(name)`
(`app.commands.whose({name})`), `docByName(name)` (null -> first non-auto-save window's document; number or
`#123` -> window id; string -> exact name, then case-insensitive / `.cdxml`-less / unique-prefix match),
`front(doc)`, `doCmd(name)` and `restoreFocus()`.  The script returns one JSON line `{ok, result}` or
`{ok: false, error}`; Python raises `ChemDrawError` for the latter.  Timeouts: 60 s default, 90 s for
open/insert, 120 s for save/commands/menus, 45 s for the status snapshot, 20 s for version.  While a script
runs, `_run_js` polls ChemDraw for modal alerts (first after 1.5 s, then every 2 s) and dismisses them; at
the timeout the `osascript` process is killed and the dismissed alert texts become part of the error.
`snapshot()` returns version, documents and the front document in one round trip for `status`.

**Front window without activation.**  Commands and pasteboard operations act on the front *document window inside
ChemDraw*; they work while the application itself stays in the background (verified: open, selectAll, copy/paste
commands and Struct=Name all succeed with another app frontmost).  `front(doc)` therefore only reorders windows:
`doc.window().index = 1`, then an Accessibility `AXRaise` through System Events if that did not take effect, and
`app.activate()` only as a last resort when `CHEMDRAWMCP_ALLOW_ACTIVATE=1` (then `restoreFocus()` gives focus back
afterwards).  The user keeps keyboard focus throughout; this matters because ChemDraw hotkeys relabel selected atoms
when stray keystrokes reach it.

**Documents.**  `app.documents` (name, modified, file, window, counts of atoms/bonds/molecules/captions/
arrows); the front document is the first `app.windows()` entry with a document; auto-save documents
(`ChemDraw Auto-Save …`) are filtered everywhere.  New documents: open the stationery
`<app>/Contents/Resources/SpecialPurpose/Stationery/<title>.cds` (ChemDraw creates "Untitled <title>-N"), or
`app.Document()` pushed into `app.documents` / the `newDocument` command.  `app.open(Path)` for files
(already-open files are just raised).  `doc.close({saving: "no"|"yes"})`.

**Save-as types.**  `app.save(doc, {in: Path, as: "<file type name>"})` with the names in
`base.EXPORT_FORMATS`: "ChemDraw XML", "ChemDraw", "MDL Molfile V2000" (`mol`), "MDL Molfile" (`molv3000`),
"MDL SDfile V2000" (`sdf`), "MDL SDfile" (`sdfv3000`), "Reaction Molfile V2000 (RXN)" (`rxn`),
"Reaction Molfile (RXN)" (`rxnv3000`), "MDL RDfile", "Chemical Markup Language (CML)", "Connection Table",
"PNG", "Scalable Vector Graphics (SVG)", "PDF", "Encapsulated PostScript", "TIFF", "GIF", "JPEG", "BMP",
"3D Manufacturing Format for Display (3MF)".  ChemDraw 23's unsuffixed Molfile / SDfile / RXN types write
V3000, which is why the V2000 variants back the plain `mol` / `sdf` / `rxn` keys.  The target is deleted
beforehand; if nothing appears, files with an extension ChemDraw appended are looked for, then dialogs are
dismissed and their text becomes the error.

**Settings.**  JXA document properties `fixedLength`, `lineWidth`, `boldWidth`, `marginWidth`,
`hashSpacing`, `labelSize`, `captionSize` are in 1/20 pt (divided by 20); `bondSpacing` (%), `chainAngle`
(°), `magnification` (%), `labelFont`, `captionFont`, `fractionalWidths`, `rulersVisible`, `crosshairVisible`
are plain.

**Pasteboard round trips.**  The general pasteboard is driven through NSPasteboard (`ObjC.import("AppKit")`
inside the JXA script, `JS_PASTEBOARD`), never through `pbcopy`/`pbpaste`, which would keep only plain text.
`pbSnapshot` records every item of the user's clipboard with every flavour it carries (images, Finder files,
rich text, ChemDraw's own CDX/CDXML flavours, as base64 per type) and `pbRestore` writes them back; snapshot,
paste/copy and restore all happen inside the same script, so nothing can slip in between and the clipboard is
restored even when the paste fails.  In: write the text as `public.utf8-plain-text`, run `pasteAsCdxml` /
`pasteAsSMILES` / `pasteAsInChi` / `pasteAsMol` / `pasteAsStructureToName` / `pasteAsHELM` / `pasteAsSequence` /
`pasteAsPeptide` / `pasteAsDNA` / `pasteAsRNA` / `paste` on the key window, wait (up to 2 s) until the object
counts change, report before/after counts and the new selection.  Out: clear the pasteboard, run `copyAsSMILES` /
`copyAsInChi` / `copyAsInChIKey` / `copyAsMol` / `copyAsMolV3000` / `copyAsCdxml` / `copyAsHELM` / `copyAsSLN` /
`copyAsHELMUsingNaturalAnalog`, wait until the pasteboard's `changeCount` moves (up to 3 s), read the plain-text
flavour.  The whole document is read by saving it as CDXML into a temporary directory.

**Selection.**  `doc.selection.smiles() / molecularFormula() / molecularWeight() / exactMass() /
elementalAnalysis() / characters() / bounds()` and the counts of its atoms/bonds/molecules/captions/arrows;
`empty` when none of atoms/captions/arrows.  Deselect = `selectAll` followed by `invertSelection`.

**Commands and menus.**  `app.commands` (name, enabled, shortDescription, longDescription; about 1600) run
with `app.do(command)`; listing with descriptions is capped at 400 matches.  `app.menus` / `menuItems` by
title, run with `app.do(item)`.  After every command or menu item the backend dismisses dialogs and returns
their text as `alerts`.

**Modal alerts.**  A modal alert blocks the Apple event.  `JS_UI` (a second osascript against System Events)
walks ChemDraw's windows: sheets of any window, windows with subrole `AXDialog` / `AXSystemDialog`, unnamed
`AXFloatingWindow`s; collects `staticTexts` and `buttons`, clicks the first of *OK, Yes, Don't Save, Discard,
Cancel, Close, No* (else the first button).  Run on timeouts, on "timed out"/"not available" errors, after
commands, and by `dismiss_dialogs`/`status`.  This needs the Accessibility/Automation permission for System
Events.

**Screenshots.**  `screencapture -l <window id> -x -o` (Screen Recording permission).

## Windows COM design

`backends/windows.py` was written on macOS from the published ChemDraw/ChemOffice automation examples and
the ChemScript documentation, without a Windows machine, so it is deliberately table-driven:

* `PROG_IDS` (`ChemDraw.Application`, then versioned ids) tried with `GetActiveObject` (running instance)
  and then `Dispatch`/`CreateObject` (starts ChemDraw; `connect_note` records when a second instance may
  have been started); the exe is found via `CHEMDRAW_EXE`, the registry (App Paths, the COM server's and the
  document class's `LocalServer32`, vendor keys incl. `WOW6432Node`) and `Program Files` globs.
* **One STA worker thread** (`_ComWorker`) owns the COM connection; every public method body is submitted to
  it (`@_com_method`).  The calling thread waits and meanwhile polls the ChemDraw process for `#32770` dialogs
  (`DIALOG_CLASSES`) every `dialog_poll` second, dismissing them with the preferred buttons
  (`DEFAULT_DIALOG_BUTTONS`: OK, Don't Save, No, Discard, Close, Cancel, Yes last), the default button,
  IDOK/IDCANCEL or WM_CLOSE — a modal dialog is what blocks a COM call forever.  A call that still has not
  returned after `call_timeout` (120 s) abandons the worker; the next call reconnects on a fresh thread.
* Every COM member is reached through `MEMBERS["<object>.<concept>"]`, a list of candidate names
  (`"doc.objects": ["Objects", "Content", "Contents"]`); `_get`/`_invoke`/`_set` try them in order, treat
  methods and properties alike (`_is_method`, falling back to `IDispatch::Invoke` for parameterised
  properties such as `Data(mime)`), detect 0- vs 1-based `Item()` and raise a `ChemDrawError` that names
  the key, the candidates and each failure.
* Data exchange uses `Objects.Data(mime)` / `Objects.SetData(mime, text)` — tried on `Document.Objects` and
  on the `Document` itself — with the `MIME` table (`chemical/x-cdxml`, `chemical/x-daylight-smiles`,
  `chemical/x-inchi`, `chemical/x-mdl-molfile`, `chemical/x-name`, `chemical/x-helm`, ...) and
  `MIME_ALTERNATIVES` spellings; sequences are converted to HELM first (`sequence_to_helm`), captions become
  a one-caption CDXML (`caption_cdxml`), and a clipboard `Paste()` is the last resort.  `insert` compares
  id-independent object signatures of the CDXML before and after (`cdxml_object_signatures`): when `SetData`
  turns out to *replace* the document (`INSERT_MODE="auto"`), the previous content is merged back with
  `merge_cdxml` (ids renumbered, fonts mapped, new objects placed below the old ones) and `setdata_replaces`
  is remembered.
* Exports use `Document.SaveAs(path)`, then `SaveAs(path, mime)`, `SaveAs(path, "<type name>")`,
  `SaveAs(path, "<ext>")`, each result sniffed (`_sniff_format`: magic bytes / markers) so a file written in
  the wrong format is not accepted, and finally `Objects.Data(mime)` written to disk.
* Document settings come from a COM `Settings` object when there is one, else from the CDXML root attributes
  (`cdxml_settings` / `cdxml_with_settings`, same keys as `styles.py`).
* Window ids are `Document.HWnd` / `Document.Window.HWnd` when exposed; `_doc()` resolves a number as such a
  handle and otherwise as a 1-based position in the document list.
* Commands are a curated table (`COMMANDS`, with `COMMAND_ALIASES` for the macOS spellings such as
  `copyAsMol`) implemented as `_cmd_*` handlers; menus are not available.
* `ChemScriptBridge` optionally imports the ChemScript SDK (`ChemScript25`...`ChemScript`, `CHEMSCRIPT_PATH`,
  `CHEMSCRIPT_MODULE`) and is used as a conversion fallback when `Data()` cannot produce a format.

**Verification checklist for the first Windows user** (the ten points in the module docstring, expanded in
[windows.md](windows.md)): connecting (ProgID, Running Object Table, same integrity level as ChemDraw);
`SetData` add-vs-replace semantics; which object holds `Data`/`SetData`; MIME spellings; `SaveAs` format
selection; the `Selection` collection; the property members; the optional `Settings` object; dialog
detection; ChemScript import.

## Headless fallbacks

| tool | without ChemDraw |
|------|------------------|
| `draw_structures` | CDXML written; RDKit `MolsToGridImage` PNG; `opened: false` |
| `draw_reaction` | CDXML written; RDKit `DrawReaction` PNG (no text) |
| `draw_scheme` | CDXML written; RDKit grid of all components |
| `draw_mechanism` | CDXML written (with the curved arrows); RDKit grid of the structures only |
| `named_reactions`, `suggest_reactions`, `apply_named_reaction` | unchanged (RDKit reaction SMARTS); `apply_named_reaction(draw=true)` falls back like `draw_reaction` |
| `quick_render` | RDKit PNG or SVG only |
| `name_to_structure` | dictionary -> OPSIN -> PubChem (online only on request) |
| `properties` | RDKit only (`structure` required) |
| `convert` | RDKit: smiles, inchi, inchikey, mol, molv3000, sdf, cdxml, png, svg |
| `read_file` | unchanged (CDXML parsed directly, CDX converted by `cdx.py`) |
| `analyze_reaction`, `analyze_image`, `compare_renders`, `styles`, `status` | unchanged |
| everything else | `ChemDrawNotAvailable` error with a hint |

The same fallbacks apply on macOS/Windows when `CHEMDRAWMCP_BACKEND=none` is set.

## Environment variables

| variable | effect |
|----------|--------|
| `CHEMDRAW_APP`, `CHEMDRAW_EXE` | explicit application bundle (macOS) or `ChemDraw.exe` (Windows); checked before any search |
| `CHEMDRAWMCP_BACKEND` | `auto` (default), `none`/`headless`/`off`, `macos`/`mac`/`darwin`, `windows`/`win`/`win32` |
| `CHEMDRAWMCP_HOME` | state directory, default `~/.chemdrawmcp` (created) |
| `CHEMDRAWMCP_OUTPUT_DIR` | where renders, CDXML files and conversions go, default `$CHEMDRAWMCP_HOME/renders` |
| `CHEMDRAWMCP_AUTOLAUNCH` | `0` -> never start ChemDraw from a tool call (also keeps name resolution offline) |
| `CHEMDRAWMCP_STYLE` | default style for `add_caption`, `add_arrow`, `name_to_structure(draw=true)`, `convert(to="cdxml")` and temporary documents (default `acs`; the `draw_*` tools have their own `style` argument defaulting to `acs`) |
| `CHEMDRAWMCP_ALLOW_ACTIVATE` | `1` -> allow activating ChemDraw as a last resort when a window cannot be raised otherwise (macOS; off by default) |
| `CHEMDRAWMCP_DEBUG` | `1` -> log every JXA call and every dialog poll with timings to stderr (macOS) |
| `CHEMDRAWMCP_ONLINE` | `1` -> allow PubChem lookups for unknown names |
| `OPSIN_JAR` | path of an OPSIN jar for offline IUPAC name parsing (`java` needed); `py2opsin` or an `opsin` executable are found automatically |
| `CHEMDRAWMCP_NO_FONTS` | `1` -> text metrics from the built-in Arial table instead of measuring installed fonts (reproducible layouts, e.g. in tests) |
| `CHEMDRAWMCP_FONT_DIR` | extra directory searched for `.ttf`/`.ttc` fonts by `textmetrics` |
| `CHEMDRAWMCP_LIVE` | `1` -> enable `tests/test_live_macos.py` (together with `pytest -m live`) |
| `CHEMDRAWMCP_CDX_PAIRS` | directory of `X.cdx` + `X.cdxml` pairs saved by ChemDraw for the CDX reader round-trip tests (default `~/.chemdrawmcp/testdata/cdxpairs`) |
| `CHEMDRAWMCP_SAMPLES` | folder of ChemDraw's sample drawings parsed by `tests/test_cdxml_robustness.py` (default: `<ChemDraw.app>/Contents/Resources/SpecialPurpose/Samples` of the discovered app; those tests skip without it) |
| `CHEMSCRIPT_PATH`, `CHEMSCRIPT_MODULE` | Windows: where ChemScript lives / its import name |

## File locations

| path | content |
|------|---------|
| `~/.chemdrawmcp/` | state directory (`CHEMDRAWMCP_HOME`) |
| `~/.chemdrawmcp/renders/` | every CDXML, PNG/SVG/..., screenshot, NMR render and conversion the server writes (`CHEMDRAWMCP_OUTPUT_DIR`); nothing is cleaned up automatically |
| system temp `chemdrawmcp_*/document.cdxml`, `chemdrawmcp_shot_*.png` | transient files for `read_document` and `screenshot` (deleted immediately) |
| `~/.claude.json` (`mcpServers.chemdraw`), `~/.codex/config.toml` (`[mcp_servers.chemdraw]`), `~/Library/Application Support/Claude/claude_desktop_config.json` or `%APPDATA%\Claude\claude_desktop_config.json` | MCP client registrations written by the installer; backups as `*.chemdrawmcp-backup` |
| `<repo>/.venv/bin/chemdraw-mcp` (`.venv\Scripts\chemdraw-mcp.exe`) | the server command (stdio) |
| `<ChemDraw.app>/Contents/Resources/SpecialPurpose/Stationery/*.cds` (macOS), `...\ChemDraw Items\*.cds` (Windows) | the stationery `new_document` opens |
