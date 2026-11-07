# Windows guide

The Windows backend (`chemdrawmcp/backends/windows.py`) drives ChemDraw through COM automation
(`ChemDraw.Application`), with an optional bridge to the ChemScript Python SDK.  It was written on macOS
from the published object model **without a Windows machine**, so the first Windows user has to verify a few
things.  The module is built for that: every uncertain name lives in editable tables at the top of the file,
every error says exactly which candidate failed and which table to edit, and nothing in the headless engine
depends on it.  The ten-point checklist in the module docstring is the authoritative list; this page walks
through it.

## 1. Install

Requirements: Windows 10/11, ChemDraw (any edition that registers `ChemDraw.Application`), Python 3.10+ on
`PATH`.

```bat
git clone <repo> chemdraw-mcp
cd chemdraw-mcp
python install.py
```

`install.py` creates `.venv` (with `uv` if present), installs the package as `chemdrawmcp[windows]` — the
extra pulls in **pywin32** — and registers the `chemdraw` server with Claude Code (`%USERPROFILE%\.claude.json`),
Codex (`%USERPROFILE%\.codex\config.toml`) and Claude Desktop (`%APPDATA%\Claude\claude_desktop_config.json`,
only if that file exists).  The server command is `.venv\Scripts\chemdraw-mcp.exe`.  `comtypes` is accepted
as an alternative COM library (`pip install comtypes` into the venv) if pywin32 cannot be installed.

```bat
python install.py --check
```

prints what was found: the ChemDraw executable (from `CHEMDRAW_EXE`, the registry — App Paths, the
`LocalServer32` of the application and document COM classes, `SOFTWARE\Revvity Signals|PerkinElmer|PerkinElmerInformatics|CambridgeSoft\ChemDraw|ChemOffice`
including `WOW6432Node` — then `Program Files*` globs), its version, RDKit, and the paths.  If ChemDraw is
installed somewhere unusual set `CHEMDRAW_EXE=C:\path\to\ChemDraw.exe` in the server's environment (the
installer writes an empty `env` block into `.claude.json` you can fill in).

Run the MCP server at the **same integrity level** as ChemDraw (both elevated or both not); COM refuses
cross-level connections.

## 2. What to verify first

Do these in the venv's Python with ChemDraw **running** and one document open, before trusting any tool.

### 2.1 Connecting: the ProgID and the Running Object Table

```python
import win32com.client
app = win32com.client.GetActiveObject("ChemDraw.Application")   # attach to the running instance
print(app.Version)
```

* `GetActiveObject` fails with *"Operation unavailable"* when ChemDraw does not register itself in the ROT.
  The backend then falls back to `Dispatch`, which may **start a second ChemDraw instance** (invisible until
  `Visible = True`); `launch()` / `capabilities()` report this as `connect_note`.  If that happens, run
  ChemDraw once as administrator to re-register its COM server, and tell us.
* If the ProgID is different (e.g. only `ChemDraw.Application.23` exists), put it first in `PROG_IDS`.
  The registry editor under `HKEY_CLASSES_ROOT` shows the registered ProgIDs (search for `ChemDraw.Application`).

### 2.2 The members the backend needs

```python
docs = app.Documents; print(docs.Count)
doc = docs.Item(1)                       # or docs(1) / docs[0]: the backend detects 0- vs 1-based
print(doc.Name, doc.FullName, doc.Saved)
objs = doc.Objects
print(objs.Data("chemical/x-cdxml")[:200])
print(objs.Data("chemical/x-daylight-smiles"))
print(objs.Formula, objs.MolecularWeight, objs.ExactMass, objs.ElementalAnalysis)
print(doc.Selection.Count)
```

`dir()` on a late-bound `CDispatch` object often shows nothing useful.  To see the real member names use the
type library: `win32com.client.gencache.EnsureDispatch("ChemDraw.Application")` then `dir(app)`, or
`python -m win32com.client.makepy` and pick the ChemDraw type library, or Microsoft's OLE/COM Object Viewer
on `ChemDraw.exe`.  The names you find go into the `MEMBERS` table (section 3).  COM name lookup is
case-insensitive for late-bound access, so only genuinely different names matter.

### 2.3 The backend end to end

The module docstring's one-liner, expanded:

```python
from chemdrawmcp.backends.windows import WindowsBackend
b = WindowsBackend()
print(b.capabilities())          # com_library, com_thread, prog_ids, connect_note, chemscript, insert_mode, setdata_replaces, commands, notes
print(b.launch())                # attaches (already_running=True) or starts ChemDraw; connect_note if a 2nd instance is possible
print(b.documents())             # DocumentInfo list: names, counts, window_id (= Document.HWnd / Window.HWnd when exposed)
name = b.new_document()          # Documents.Add()
r = b.insert(name, "CCO", "smiles"); print(r)      # mode: added | replaced+merged | nothing; atoms_added should be 3
print(b.selection_info(name))    # smiles / formula / molecular_weight / exact_mass / elemental_analysis; 'unavailable' lists the rest
print(b.selection_as(name, "inchikey"))
print(b.document_cdxml(name)[:300])
print(b.save_document(name, r"C:\Temp\ethanol.png", "png"))   # the written file is sniffed; see the error text if SaveAs ignored the extension
print(b.document_settings(name)) # source: com | cdxml
print(b.pending_dialogs())       # what the dialog enumerator sees (should be [])
b.close_document(name)
```

Every failure is a `ChemDrawError` whose text names the table entry to fix, e.g.

```
ChemDraw COM: none of the candidate members ['Objects', 'Content', 'Contents'] for 'doc.objects' worked on
CDispatch (Objects: not present; Content: not present; Contents: not present). If your ChemDraw version names
this member differently, edit MEMBERS['doc.objects'] in chemdrawmcp/backends/windows.py.
```

Then run the MCP tools through your client: `status`, `draw_structures(["aspirin"])`, `read_document`,
`render`, `structure_to_name("CCO")`, `save_document(... format="cdx")`.

### 2.4 Checklist (in the order most likely to need a fix)

1. **Connecting** — `PROG_IDS`, ROT registration, integrity level (2.1).
2. **Insert semantics** — does `Objects.SetData` *add to* or *replace* the document?  `insert` detects
   replacement (an id-independent signature of some old object is gone afterwards, `cdxml_object_signatures`)
   and merges the old content back (`merge_cdxml`); the result's `mode` says `added` or `replaced+merged`.  Force the behaviour with `INSERT_MODE = "add"` / `"replace"` and pick the
   collection with `INSERT_TARGET` (`"objects"` or `"selection"`).
3. **Data holders** — `Data`/`SetData` are tried on `Document.Objects` and on the `Document` itself
   (`MEMBERS["doc.objects"]`, `["objs.data"]`, `["objs.set_data"]`).
4. **MIME spellings** — `MIME` / `MIME_ALTERNATIVES`; `selection_as` and `insert` report which ones were
   refused (`fallbacks`).
5. **SaveAs** — `Document.SaveAs(path)` is expected to pick the format from the extension; the file is sniffed
   afterwards (`_sniff_format`) and `SaveAs(path, mime|name|ext)` plus `Objects.Data(mime)` are tried when the
   content does not match (`EXPORT_MIME`).
6. **Selection** — `Document.Selection` must be an Objects-like collection with `Count` (or at least `Data`);
   `MEMBERS["doc.selection"]`, `doc.select_all`, `doc.deselect`.
7. **Properties** — `Objects.Formula` / `MolecularWeight` / `ExactMass` / `ElementalAnalysis`
   (`MEMBERS["objs.*"]`); `selection_info` lists what was `unavailable`.
8. **Settings** — the `Document.Settings` object is optional: without it settings are read from the document
   CDXML (`CDXML_SETTINGS`) and written back only when `SetData` is known to replace.
9. **Dialogs** — alerts are expected to be `#32770` windows of `ChemDraw.exe` (`DIALOG_CLASSES`);
   `pending_dialogs()` shows what the enumerator sees.
10. **ChemScript** — the optional `ChemScriptNN` module only imports in the interpreter ChemDraw installed it
    for (section 4).

## 3. The tables at the top of `windows.py` and how to fix them

Everything uncertain about ChemDraw's COM API is data, not code:

| table | what it controls | how to fix |
|-------|------------------|------------|
| `PROG_IDS` | ProgIDs tried for `GetActiveObject` then `Dispatch` | put the ProgID that exists first |
| `FILE_PROG_IDS`, `FILE_EXTENSIONS` | document/ActiveX classes used only to locate `ChemDraw.exe` in the registry | add yours |
| `PROCESS_NAMES` | image names that count as "ChemDraw is running" (`psutil` if installed, else `tasklist`) | add your exe name in lower case |
| `MIME`, `MIME_ALTERNATIVES` | MIME strings for `Objects.Data` / `SetData`, and alternative spellings tried after a refusal | put the accepted spelling in front |
| `INSERT_MIME`, `SEQUENCE_KINDS`, `SEQUENCE_FALLBACK_MIME` | which MIME each `insert(kind=...)` uses; sequences are converted to HELM (`sequence_to_helm`), and if HELM is refused the raw sequence is offered with the fallback MIME names (guesses) | adjust after testing `insert(..., "peptide")` |
| `SELECTION_FORMATS` | MIME per `selection_as` format (incl. `name`) | — |
| `EXPORT_MIME` | second `SaveAs` argument / `Objects.Data` MIME per export format (defaults to `base.EXPORT_FORMATS[...]["mime"]`) | override single entries |
| `INSERT_MODE`, `INSERT_TARGET` | `"auto"` (detect add vs replace) / `"add"` / `"replace"`; paste into `"objects"` or `"selection"` | set once you know how your version behaves |
| `MEMBERS` | candidate member names, keyed `"<object>.<concept>"`: `app.*`, `coll.*`, `docs.*`, `doc.*` (incl. `doc.window`, `doc.hwnd`, `win.hwnd`), `objs.*` (Objects and Selection collections, also tried on the Document), `obj.text`, `rect.*` | put the real name first (or alone); keep the key |
| `SETTINGS_MEMBERS`, `CDXML_SETTINGS` | members of `Document.Settings` per settings key; the CDXML root attributes used instead when there is no `Settings` object | same |
| `REGISTRY_KEYS`, `APP_PATHS_KEY`, `REGISTRY_PATH_VALUES`, `EXE_GLOB_SUFFIXES` | where `find_chemdraw_windows()` looks for the exe | add your vendor/version path |
| `STATIONERY_GLOB_SUFFIXES` | where `new_document(style)` looks for `*.cds` (the `ChemDraw Items` folder) | add the folder |
| `CHEMSCRIPT_MODULES`, `CHEMSCRIPT_GLOB_SUFFIXES`, `CHEMSCRIPT_PYTHON_DIRS`, `CHEMSCRIPT_MEMBERS` | ChemScript import names, search paths and `StructureData` members (`ChemicalName()` is Struct=Name) | see section 4 |
| `COMMANDS`, `COMMAND_ALIASES` | the named commands `run_command` offers (each backed by a `_cmd_<name>` method) and the macOS spellings accepted as synonyms | add a handler to extend |
| `DIALOG_CLASSES`, `DEFAULT_DIALOG_BUTTONS` | window classes treated as modal dialogs; button labels pressed in order (`Yes` last, because answering "Save changes?" with Yes opens a file dialog) | add a class / reorder |

How a lookup works (`WindowsBackend._get` / `_invoke` / `_set`): for each candidate, `getattr`; if that is a
method (`_is_method`: Python method, pywin32 `MethodCaller`/`NamedProperty`, comtypes `method-wrapper`) it is
called with the given arguments, else the property value is returned; when the client library cannot tell
property from method (comtypes, parameterised properties like `Data(mime)`) the member is invoked through
`IDispatch::GetIDsOfNames`/`Invoke`; the first candidate that does not raise wins.  "Member not found"
HRESULTs (`DISP_E_UNKNOWNNAME`, `DISP_E_MEMBERNOTFOUND`) move on to the next candidate.  Collections are read
with `Count` (or `len`) and `Item(i)`; the first successful index base (1 or 0) is remembered.  COM errors are
described with HRESULT, description and source (`_describe_com_error`); the RPC "server unavailable"
HRESULTs make the backend reconnect on the next call.

**Threading and dialogs.**  All COM calls run on one dedicated STA thread (`_ComWorker`); the calling thread
waits and polls for modal dialogs every second (`dialog_poll`), dismissing them and restarting the time
budget.  A call that does not return within `call_timeout` (120 s) is abandoned ("the COM call is still
blocked inside ChemDraw ..."), and the next call reconnects on a fresh thread.  `WindowsBackend(com_thread=False)`
runs calls on the caller's thread instead (then `CoInitialize` yourself).

Because the tables are plain Python, the quickest loop is: run a tool, read the error, edit the table, rerun.
`tests/test_windows_backend.py` runs on any OS against a fake COM application (`FakeApp`/`FakeDocument`/
`FakeObjects` with `Objects.Data/SetData/Clean`, `Documents.Add/Open/Item`, `SaveAs/Save/Close/Activate/SelectAll`);
keep it green after edits and consider adding the members you discovered to the fakes.

## 4. ChemScript setup (optional)

ChemScript is the Python SDK ChemDraw/ChemOffice installs (`ChemScript<version>` package).  The backend uses
it only as a **conversion fallback** (`_convert_fallback`: when `Objects.Data(mime)` cannot produce a format,
the molfile/SMILES/CDXML of the selection is converted through `StructureData.LoadData` + `WriteData`;
RDKit is the last resort for InChI/InChIKey/molfile).  `ChemScriptBridge` also offers
`name_to_structure`, `structure_to_name` (`ChemicalName()`), `properties` for your own scripts.

* ChemScript is importable only from the Python interpreter it was installed for (typically `C:\Python39`),
  not from the project venv.  Either run the server with that interpreter (`<that python> -m pip install -e .`
  in the repo and point your MCP client at `<that python> -m chemdrawmcp.server`), or set
  `CHEMSCRIPT_PATH` to the directory holding the `ChemScriptNN` package and `CHEMSCRIPT_MODULE` to its name
  (the DLLs must still load in your interpreter's bitness).
* Check: `python -c "from chemdrawmcp.backends.windows import chemscript_available; print(chemscript_available())"`;
  `WindowsBackend().capabilities()["chemscript"]` and `["chemscript_module"]` show the same.
* The `StructureData` member names are candidates in `CHEMSCRIPT_MEMBERS`; fix them the same way.

## 5. Known unknowns

Things the code assumes and nobody has confirmed on a real Windows ChemDraw:

1. **ProgID and ROT registration** (2.1).  Without ROT registration every server start may create a second
   ChemDraw instance (`connect_note`).
2. **`Objects.SetData` semantics** (checklist 2).  `merge_cdxml` has only been exercised against the fake COM
   model; check that the merged document looks right (old content on top, new objects 36 pt below).
3. **`Document.SaveAs(path[, format])`** (checklist 5).  Vector formats (SVG/PDF/EPS) through `Objects.Data`
   are unverified; the format sniffing only knows signatures for the common formats.
4. **`Document.Settings`** (checklist 8).  Without it, `document_settings(settings=...)` can only write when
   `SetData` replaces the document.
5. **Window ids**: `documents()` reports `window_id` from `Document.HWnd` or `Document.Window.HWnd` when such
   members exist, and `_doc()` resolves them; a number that matches nothing is read as a 1-based index into
   the document list.  If your ChemDraw exposes no handle, `window_id` is `null` and tools fall back to
   names, which change when an untitled document is saved.
6. **Dialogs**: alerts are assumed to be standard `#32770` dialogs of the ChemDraw process; custom dialog
   classes will not be detected (add them to `DIALOG_CLASSES`); a blocked COM call then times out after
   `call_timeout`.
7. **Deselect**: `doc.deselect` candidates (`Deselect`, `DeselectAll`, `ClearSelection`, `UnselectAll`,
   `SelectNone`) are guesses; `select(mode="none")` may be unavailable.
8. **Struct=Name through COM**: `Objects.Data("chemical/x-name")` (then `Objects.ChemicalName`) is assumed to
   return the name; `structure_to_name` inserts the result as a caption.
9. **Commands**: only the curated `COMMANDS` table is available (`cleanStructure`, `convertStructureToName`,
   `copyAs*`, `chemicalProperties`, ...); `cleanReaction`, `checkStructure`, `mapReactionAtoms`,
   `predict1H_NMR_Shifts`, `selectNextStructure`, `undo` and the rest of ChemDraw's ~1600 named commands
   are **not** reachable, and `menus`/`do_menu` raise.  Tools that depend on them (`check_structure`,
   `predict_nmr`, `clean_structure(reaction=true)`, `select(next|previous|invert)`) error on Windows.
10. **Screenshots** use `PIL.ImageGrab` on the main window rectangle after `SetForegroundWindow`; multi-monitor
    DPI scaling may offset the capture.
11. **Bitness**: `ChemDraw.Application` is an out-of-process server, so a 64-bit Python should talk to a
    32-bit ChemDraw and vice versa; untested.
12. **Stationery location** (`ChemDraw Items`) and the `Documents.Add(path)` vs `Documents.Open(path)`
    behaviour for `.cds` files.

When you have answers, please record them in the module docstring of `windows.py` (the checklist) and in this
file, and send the output of `WindowsBackend().capabilities()` plus the exact error strings with any bug
report.
