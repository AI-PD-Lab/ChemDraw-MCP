# Development

## Package layout

```
chemdraw-mcp/
├── install.py                  bootstrap: venv (uv if available), pip install -e .[windows] on Windows, registration
├── pyproject.toml              hatchling build; scripts chemdraw-mcp / chemdrawmcp-install; pytest (live marker) + ruff config
├── README.md                   entry point
├── CHANGELOG.md, LICENSE       release notes (Keep a Changelog), PolyForm Noncommercial licence
├── .github/workflows/ci.yml    GitHub Actions: pytest on Linux/macOS/Windows x Python 3.10/3.12, ruff
├── docs/                       tools.md, agent-guide.md, architecture.md, windows.md, development.md
├── chemdrawmcp/
│   ├── server.py               the 43 MCP tools (FastMCP instance `mcp`, `main()` runs stdio)
│   ├── mcpcompat.py            FastMCP/Image import shim for mcp 1.x and 2.x (MCP_MAJOR)
│   ├── backends/
│   │   ├── __init__.py         get_backend(): CHEMDRAWMCP_BACKEND / platform selection
│   │   ├── base.py             ChemDrawBackend ABC, DocumentInfo, CommandInfo, EXPORT_FORMATS, INSERT_KINDS, NullBackend
│   │   ├── macos.py            MacBackend (JXA via osascript, dialog polling, snapshot)
│   │   └── windows.py          WindowsBackend (COM on one STA worker thread), ChemScriptBridge, editable tables, merge_cdxml
│   ├── cdxml.py                CDXMLDocument writer (label-aware fragment boxes), read_cdxml / fragment_to_mol reader
│   ├── cdx.py                  binary CDX reader (parse_cdx, cdx_to_cdxml, read_cdx, dump_cdx)
│   ├── layout.py               build_structures / build_reaction / build_scheme, collisions, layout_warnings
│   ├── textmetrics.py          text widths from real fonts or an Arial table, formula runs, predicted atom labels
│   ├── abbreviations.py        nickname -> SMILES tables (substituents, reagents), parse_atom_label
│   ├── mechanism.py            electron-pushing arrows (parse_arrow_spec, add_electron_arrows, build_mechanism) -> draw_mechanism
│   ├── reactions.py            100 named-reaction templates (reaction SMARTS + reagents/conditions) -> named_reactions, apply_named_reaction, suggest_reactions
│   ├── styles.py               STYLES (settings of 19 ChemDraw stationery files), ALIASES, get_style
│   ├── chem.py                 RDKit utilities (parsing, depiction prep, properties, reactions, rendering, names)
│   ├── imageanalysis.py        numerical description of renders, diffs, thumbnails
│   ├── colornames.py           named-colour lookup for the palette report
│   ├── discovery.py            find ChemDraw (app/exe, version, bundle id), home/output dirs
│   └── install.py              MCP client registration (--claude --codex --claude-desktop --uninstall --dry-run --check)
└── tests/
    ├── test_server.py          every tool headless (CHEMDRAWMCP_BACKEND=none), return shapes, error envelopes
    ├── test_cdxml.py           writer/reader round trips, arrows, colours, page tiling, layout engines
    ├── test_cdxml_robustness.py reader edge cases: nicknames, residues, label-derived charges, wedges, brackets, ChemDraw's sample files
    ├── test_cdx.py             binary CDX reader: hand-assembled byte strings + round trips against ChemDraw-saved cdx/cdxml pairs
    ├── test_layout_quality.py  label-aware bounds, '+' and arrow placement, grids, text metrics, overlap detection
    ├── test_mechanism.py       arrow-spec parsing, endpoint planning, SN2 / carbonyl / E1 / resonance mechanisms
    ├── test_reactions.py       one worked example per named-reaction template, fuzzy lookup, apply/suggest/to_step
    ├── test_chem.py            format detection, parsing, depiction prep, properties, reactions, name resolution (offline)
    ├── test_windows_backend.py Windows backend against a fake COM application (runs on every OS)
    ├── test_macos_backend.py   macOS backend headless: the generated JXA runs under node against a fake ChemDraw object model, pasteboard helpers on a private named pasteboard
    ├── test_names.py           the built-in name dictionary and resolve_name normalisation, parametrised over the dictionary (~1,300 tests, offline)
    ├── test_install.py         registration against temporary config files; ChemDraw discovery on this Mac
    └── test_live_macos.py      live tests against the real ChemDraw (macOS, opt-in; marker `live`)
```

Adding a named reaction: append a `_t(key, name, cls, smarts, reagents, conditions, ...)` entry to
`_TEMPLATES` in `reactions.py` following the SMARTS conventions in the module docstring (only drawn products
on the product side, by-products as text, grouped `(A.B)` products for ring openings, repeated reactant
templates plus `stoichiometry` for > 1 equiv), give it `roles`, `aliases` and an `example`, and add the
worked example to `EXAMPLES` in `tests/test_reactions.py` — `test_registry_size_and_required_keys` and
`test_worked_example` fail otherwise.

## Running the tests

```bash
python install.py                      # once: creates .venv and installs the package (+ rdkit, numpy, pillow, mcp)
.venv/bin/pip install pytest           # or: uv pip install --python .venv/bin/python -e '.[dev]'
.venv/bin/python -m pytest -q          # ~2,070 unit tests; never touch ChemDraw or the network
```

`pyproject.toml` registers the `live` marker and sets `addopts = "-m 'not live'"`, so the default run
excludes the live tests.  The live suite drives the real application and is doubly opted in:

```bash
CHEMDRAWMCP_LIVE=1 .venv/bin/python -m pytest tests/test_live_macos.py -m live -v
```

Its 31 tests are skipped unless `sys.platform == "darwin"`, `CHEMDRAWMCP_LIVE=1` is set and `discovery.find_chemdraw()`
finds ChemDraw.app (`CHEMDRAW_APP` overrides the search).  It launches ChemDraw if needed (ChemDraw stays in the background: the backend never
activates it), and at the end closes *without saving* every document that was not open before the session.  What it covers: launch/version/capabilities, `js` passthrough, dialog handling when
idle, argument validation before any Apple event, documents and window ids surviving ChemDraw's rename,
insert (SMILES, name, text), select/deselect, every `selection_as` format, `save_document` for png/svg/cdxml/
pdf/mol (300 dpi PNG, vector SVG, round-tripped CDXML), extension handling, the RXN-without-reaction error,
CDXML round trips, ACS stationery settings and changing them, `cleanStructure` and `convertStructureToName`,
commands and menus listing, opening a reaction CDXML, screenshots, and the server tools `status`,
`draw_reaction` + `read_document`, `structure_to_name`, `quick_render` (no document left open) and
`new_document` + `insert` + `selection` + `render` with an image.  Expected values were verified on ChemDraw
23.1.2; assertions are tolerant (`in`, `approx`, canonical SMILES) where another version might differ.

How the unit tests stay offline and isolated:

* `tests/test_server.py` has a `server` fixture that sets `CHEMDRAWMCP_BACKEND=none`,
  `CHEMDRAWMCP_ONLINE=0`, `CHEMDRAWMCP_OUTPUT_DIR` and `CHEMDRAWMCP_HOME` to per-test directories, reloads
  `chemdrawmcp.server` and clears `_BACKEND` / `_NAME_CACHE`.  Tools are called as plain functions
  (`server.draw_structures(...)`), since with mcp 2 `@mcp.tool()` returns the undecorated function; results
  are JSON strings or `[json, Image]` (`unwrap()` helper).  The `NEEDS_CHEMDRAW` list asserts that every
  live-only tool returns the `ChemDrawNotAvailable` envelope headless.
* `tests/test_chem.py` has an autouse `offline` fixture that removes `py2opsin`, `OPSIN_JAR`, `shutil.which`
  hits and patches `urllib.request.urlopen` to raise, so name resolution exercises only the built-in
  dictionary.
* `tests/test_layout_quality.py` checks geometry (no overlapping boxes, `+` clearance, arrows longer than
  their text, captions under their structure, grid cells, wrapping) and the text metrics (monotonic, formula
  text narrower than plain, heuristic close to measured fonts).  Set `CHEMDRAWMCP_NO_FONTS=1` to make layout
  numbers identical on machines with different fonts.
* `tests/test_cdx.py` always runs its hand-assembled CDX byte strings; the round-trip half looks for
  `X.cdx` + `X.cdxml` pairs saved by ChemDraw in `CHEMDRAWMCP_CDX_PAIRS` (default
  `~/.chemdrawmcp/testdata/cdxpairs`) and skips when there are none — save a few documents in both formats
  there to exercise it.
* `tests/test_cdxml_robustness.py` parses ChemDraw's own sample drawings from
  `<ChemDraw.app>/Contents/Resources/SpecialPurpose/Samples` when that folder exists (macOS) or
  `CHEMDRAWMCP_SAMPLES` points at a copy of it, otherwise only the inline snippets run.
* `tests/test_reactions.py` carries one worked example per template (`EXAMPLES`), and
  `test_registry_size_and_required_keys` fails when a template is added without one.
* `tests/test_windows_backend.py` monkeypatches `windows._com_connect` to return a fake application; the
  pure helpers (`sequence_to_helm`, `caption_cdxml`, `cdxml_summary`, `cdxml_object_signatures`,
  `merge_cdxml`, `_ensure_extension`, `_describe_com_error`) are tested directly.
* `tests/test_macos_backend.py` never talks to ChemDraw: the JXA scripts the backend generates are executed
  with `node` against a fake object model (`FAKE_JS`; skipped when node is missing), `save_document` runs with
  `js` stubbed, and the `JS_PASTEBOARD` helpers run through `osascript` on a private *named* pasteboard, so the
  user's clipboard is never touched.
* `tests/test_names.py` walks the whole built-in dictionary and the `resolve_name` normalisation rules with
  `CHEMDRAWMCP_ONLINE=0`, so it never reaches PubChem.
* `tests/test_install.py` redirects `HOME`/`APPDATA` and the three config-path functions to a temp dir; the
  only test that looks at the real machine is the macOS-only `test_find_chemdraw_on_this_mac`.

Lint: `ruff check .` (configured in `pyproject.toml`: line length 130, rules E/F/W/I, with `E501` and `I001`
ignored globally and a few per-file ignores).

## Adding a tool

1. **Decide the engine.**  Headless-capable tools build on `chem`, `layout`, `cdxml`, `cdx`, `imageanalysis`
   and only touch the backend through `_has_live()` / `_materialize`.  Live tools get the backend with
   `b = _live()` (starts ChemDraw on demand) or `b = _live(autolaunch=False)` (for things that only make sense
   when the app is already running: reading, rendering, commands).  If the backend lacks a capability, add
   an abstract method to `backends/base.py`, implement it in `macos.py` **and** `windows.py` (or raise a
   clear `ChemDrawError` there), and add it to `NullBackend`'s `_raise` list.
2. **Write the tool** in the matching section of `server.py`:

   ```python
   @mcp.tool()
   def count_rings(structure: Optional[str] = None, document: Optional[str] = None) -> str:
       """One-paragraph description: this docstring is what the agent reads.  Say what the inputs accept,
       what comes back, and the pitfalls (needs a selection, macOS only, ...)."""
       try:
           if structure:
               mol = _mol(structure)                      # SMILES / InChI / molblock / CDXML / path / name
           else:
               b = _live(autolaunch=False)
               info = b.selection_info(_doc(document))    # _doc(): None/name/window id normalisation
               mol = chem.mol_from_input(info["smiles"], "smiles")
           from rdkit.Chem import rdMolDescriptors
           return dumps({"ok": True, "rings": rdMolDescriptors.CalcNumRings(mol)})
       except Exception as e:
           return _err(e)
   ```

   Rules: parameters must be JSON-schema friendly (`str`, `int`, `float`, `bool`, `Optional[...]`,
   `List[str]`, `Dict[str, Any]`); return `dumps(dict)` or `_with_image(dumps(dict), png_bytes)`; wrap the
   body in `try/except Exception` and return `_err(e)` so failures are JSON; write files with
   `_out_path(tag, ext, path)`; trim long text with `_trim`; address documents through `_doc()` and prefer
   returning `window_id` (`_doc_info`).  Keep the headless path working whenever the operation can be done
   with RDKit.
3. **Test it** in `tests/test_server.py`: a headless test of the success path (and the error envelope), and
   if it needs ChemDraw add it to `NEEDS_CHEMDRAW` so the `ChemDrawNotAvailable` envelope is checked; add a
   live test in `tests/test_live_macos.py` when ChemDraw behaviour matters.
4. **Document it**: the tools table in `README.md`, a section in `docs/tools.md`, and — if an agent should
   use it in a typical flow — `INSTRUCTIONS` in `server.py` and `docs/agent-guide.md`.

Backend changes on macOS are JXA snippets passed to `MacBackend.js(body, **args)`; the body is the inside of
`main()`, has `app`, `ARGS`, `cmd()`, `docByName()`, `front()`, `doCmd()` available and must `return` a JSON-
serialisable value.  Call `front(d)` before anything that acts on the key window.  On Windows add candidate
member names to `MEMBERS` and reach them through `self._get` / `self._invoke` / `self._set`; decorate public
methods with `@_com_method()` so they run on the COM worker thread.

## Writing live tests

Follow `tests/test_live_macos.py`: the module-level `pytestmark` applies the `live` marker and the skip
condition; the session-scoped `live` fixture (`LiveSession`) owns one `MacBackend(timeout=90)`, remembers the
window ids that were open before, and closes everything else at teardown.  Use `cd` for the backend and
`server` for `chemdrawmcp.server` wired to that backend with renders in `tmp_path`.

```python
def test_reaction_roundtrip(server, live, out_dir):
    res = server.draw_reaction(reaction="CC(=O)O + OCC --[H2SO4]--> CC(=O)OCC + O", return_image=False)
    info = json.loads(res[0] if isinstance(res, list) else res)
    wid = info["window_id"]
    try:
        assert info["opened"] and info["render"]["rendered_by"] == "ChemDraw"
        doc = json.loads(server.read_document(wid))
        assert len(doc["steps"]) == 1 and canon_rxn(doc["reaction_smiles"][0]) == canon_rxn(ESTERIFICATION_RXN)
    finally:
        live.close(wid)
```

Guidelines: create documents with `live.new_doc(style)` (returns `(window_id, name)` and registers the id) or
close the ones a tool opened with `live.close(wid)` in a `finally` block; address documents by window id;
write into `tmp_path`; do not assert exact pixel sizes (ChemDraw's export crops to content) but on
`analysis` fields and file signatures (`_check_png` requires 300 dpi); compare SMILES canonically; keep
assertions tolerant to ChemDraw versions.

## Regenerating `styles.py`

`STYLES` holds the *settings* of ChemDraw's stationery files, not the files themselves.  They were extracted
by opening each `.cds` in ChemDraw and saving it as CDXML (`.cds` is a binary CDX container; with
`chemdrawmcp.cdx` it could now be read directly, but the ChemDraw round trip is what was verified), then
copying the root attributes, font table, colour table and page box.  To add or refresh a style (for example
the two stationery files in ChemDraw 23.1 that are not yet covered, `Labeled Monomers.cds` and
`SYNTHESIS, SYNLETT Document.cds`, or a newer ChemDraw release):

1. In ChemDraw: *File > Open* `<ChemDraw.app>/Contents/Resources/SpecialPurpose/Stationery/<Title>.cds`
   (macOS) or the `ChemDraw Items` folder (Windows) — or with the server, `new_document(style="<Title>")`
   then `save_document("<Title>.cdxml", format="cdxml")`.
2. Turn the CDXML into a style entry:

   ```python
   import re, xml.etree.ElementTree as ET
   def style_from_cdxml(path, title):
       root = ET.fromstring(re.sub(r"<!DOCTYPE[^>]*>", "", open(path, encoding="utf-8").read()))
       fonts = {f.get("id"): f.get("name") for f in root.find("fonttable")}
       colors = [[float(c.get(k)) for k in "rgb"] for c in root.find("colortable")]
       page = [float(v) for v in root.find("page").get("BoundingBox").split()]
       num = lambda k, d=None: float(root.get(k)) if root.get(k) is not None else d
       return {
           "name": title, "stationery": title + ".cds",
           "bond_length": num("BondLength"), "line_width": num("LineWidth"), "bold_width": num("BoldWidth"),
           "margin_width": num("MarginWidth"), "hash_spacing": num("HashSpacing"), "bond_spacing": num("BondSpacing"),
           "chain_angle": num("ChainAngle"), "label_size": num("LabelSize"), "caption_size": num("CaptionSize"),
           "label_face": int(num("LabelFace", 96)),
           "print_margins": [float(v) for v in (root.get("PrintMargins") or "36 36 36 36").split()],
           "show_terminal_carbon_labels": root.get("ShowTerminalCarbonLabels") == "yes",
           "hide_implicit_hydrogens": root.get("HideImplicitHydrogens") == "yes",
           "label_justification": root.get("LabelJustification", "Auto"),
           "caption_justification": root.get("CaptionJustification", "Left"),
           "fractional_widths": root.get("FractionalWidths", "yes") == "yes",
           "label_font": fonts.get(root.get("LabelFont"), "Arial"),
           "caption_font": fonts.get(root.get("CaptionFont"), fonts.get(root.get("LabelFont"), "Arial")),
           "colortable": colors, "color": int(root.get("color", 0)), "bgcolor": int(root.get("bgcolor", 1)),
           "page": page,
       }
   ```

   Keys and units must match what `cdxml.CDXMLDocument.to_xml` reads (`bond_length`, `line_width`,
   `bold_width`, `margin_width`, `hash_spacing`, `bond_spacing`, `chain_angle`, `label_size`, `caption_size`,
   `label_face`, `caption_face` (optional), `print_margins`, `show_terminal_carbon_labels`,
   `hide_implicit_hydrogens`, `label_justification`, `caption_justification`, `fractional_widths`,
   `label_font`, `caption_font`, `colortable`, `color`, `bgcolor`, `page`), all in points.  Remember that
   colour indices 0/1 are implicit black/white and the table starts at index 2.
3. Add the entry to `STYLES` under a short lower-case key, an alias in `ALIASES` if useful, and extend the
   style tables in `README.md` and `docs/tools.md`.  `styles.get_style` also matches the stationery title, so
   `new_document(style="<Title>")` opens the real file as long as `stationery` is right.
4. `tests/test_cdxml.py::test_styles_lookup` pins a few values (ACS 14.4 pt, Nature 6 pt labels, RSC
   Helvetica); add one for the new style.

## Coding conventions

* **Python 3.10+, no hard dependency on the GUI.**  `rdkit`, `win32com`/`comtypes`, `winreg`, `ctypes.windll`
  and `PIL.ImageGrab` are imported lazily inside functions; the package must import on Linux without
  ChemDraw.  `chem._rd()` is the single RDKit import point and `chem._quiet()` silences RDKit's C++ log.
* **Errors are messages for an agent.**  Backends raise `ChemDrawError` (or `ChemDrawNotAvailable`) with
  text that says what failed and what to do (including dismissed alert text and, on Windows, the table entry
  to edit); the chemistry layer raises `chem.ChemError`; tools convert everything into the
  `{"ok": false, "error", "error_type"}` envelope.  Never let a traceback reach the MCP client.
* **Tools return JSON strings**, built with `server.dumps` (indent 1, UTF-8, `default=str`), optionally with a
  thumbnail through `_with_image`.  Long text is trimmed with `_trim` (8000 chars).  Docstrings are the tool
  descriptions the agent sees — write them as instructions, with examples of accepted inputs.
* **Documents by window id.**  Anything that returns a document should include `window_id`
  (`_doc_info`); anything that accepts one goes through `_doc()`.
* **Files only in the renders dir** (`_out_path`) or at a path the caller gave; temp files are cleaned up.
* **Headless parity.**  Every drawing tool writes the CDXML before touching ChemDraw and has an RDKit render
  fallback; `_has_live()` guards live-only extras; file readers (`read_file`, `cdx.py`) never need the app.
* **Layouts must not overlap.**  Measure with `textmetrics` / `layout._component_extent` before placing, keep
  `layout_warnings(doc)` empty, and add a case to `tests/test_layout_quality.py` for new placement code.
* **Backends stay thin and symmetric.**  macOS: one JXA script per operation, state in ChemDraw, nothing
  cached except `last_alerts`.  Windows: all names in tables, every COM access through `_get`/`_invoke`/
  `_set` on the worker thread, every uncertainty listed in the module docstring's checklist.
* **Empirical facts are cited.**  Comments and docstrings say which ChemDraw version a behaviour was verified
  on ("verified with ChemDraw 23"), and `cdx.py` marks each tag `(v)` verified or `(s)` from the spec.
* **Tests are offline and hermetic** (see above); new tests must not depend on a running ChemDraw unless
  marked `live`, and must not hit the network.
* **Style**: ruff (E/F/W/I, line length 130); `from __future__ import annotations`; `%`-formatting in the
  older modules and f-strings in `chem.py` are both accepted — match the file you are in.
