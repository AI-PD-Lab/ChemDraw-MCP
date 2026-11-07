"""Windows backend: drives ChemDraw through COM automation (plus an optional ChemScript bridge).

How this module was written
---------------------------
It was authored on macOS without access to a Windows ChemDraw installation, from the published
ChemDraw/ChemOffice automation examples (``CreateObject("ChemDraw.Application")``, ``Documents.Add``,
``Objects.SetData "chemical/x-daylight-smiles", ...``, ``Objects.Clean``, ``Objects.Formula`` ...),
the ChemDraw ActiveX/Automation object model (Application > Documents > Document > Objects / Selection,
``Objects.Data(mime)`` / ``SetData(mime, data)``, ``Document.SaveAs``, ``Close(saveChanges)``) and
the ChemScript Python SDK documentation.  Everything is therefore split into two tiers:

* **Verified / platform independent** - the module imports and the class constructs on every OS;
  all Windows-only imports (``win32com``/``pythoncom``, ``comtypes``, ``winreg``, ``ctypes.windll``,
  ``PIL.ImageGrab``) are lazy and guarded.  Pure-Python helpers (MIME tables, HELM generation for
  sequences, caption CDXML, CDXML merging / settings, format sniffing, path handling, error
  formatting, the command table, the COM worker thread) are unit tested with a fake COM object
  model (``tests/test_windows_backend.py``).
* **Best effort** - the exact member names and call shapes of ChemDraw's COM object model.  Every
  COM member is looked up through the ``MEMBERS`` candidate table: the first candidate that exists on
  the live object wins; when a client library cannot decide whether a member is a property or a
  method (comtypes and parameterised properties such as ``Data(mime)``) the member is invoked
  through ``IDispatch::GetIDsOfNames``/``Invoke`` directly.  When nothing works the error message
  lists what was tried *and names the table entry to edit*.

Threading model
---------------
COM proxies belong to the apartment (thread) that created them; an MCP server may call tools from
the event-loop thread or from worker threads, so ``WindowsBackend`` owns **one dedicated STA thread**
(``_ComWorker``: ``CoInitialize`` once, then every public method body runs on it, serialised).  The
calling thread waits for the result and, while waiting, polls for modal dialogs (Win32 ``#32770``
windows of the ChemDraw process) and dismisses them - a modal dialog is what blocks a COM call
forever.  If a call still does not return within ``call_timeout`` seconds the worker is abandoned
and the next call reconnects on a fresh thread.  Pass ``com_thread=False`` to run calls on the
caller's thread instead (then CoInitialize per thread yourself).

Verification checklist for the first run on a real Windows machine
------------------------------------------------------------------
Run ``python -c "from chemdrawmcp.backends.windows import WindowsBackend as W; b=W(); print(b.launch());
print(b.new_document()); print(b.insert(None, 'CCO', 'smiles')); print(b.selection_info());
print(b.save_document(None, 'ethanol.cdxml', 'cdxml'))"`` and check, in this order (most likely to
need adjusting first):

1. **Connecting** - ``PROG_IDS``: ``GetActiveObject("ChemDraw.Application")`` must find the running
   instance (ChemDraw must register in the Running Object Table; if it does not, ``Dispatch`` starts
   a second ChemDraw - ``launch()`` reports ``connect_note`` when that risk exists).  Run the server
   at the same integrity level as ChemDraw (both elevated or both not), otherwise COM refuses.
2. **Insert semantics** - whether ``Objects.SetData`` *adds to* or *replaces* the document.  ``insert``
   detects replacement (object counts went down) and then re-inserts the merged CDXML (old + new,
   see :func:`merge_cdxml`); force the behaviour with ``INSERT_MODE`` = ``"add"``/``"replace"`` and
   choose the collection with ``INSERT_TARGET`` (``"objects"`` or ``"selection"``).
3. **Data holders** - ``Data``/``SetData`` are tried on ``Document.Objects`` and on the ``Document``
   itself (``MEMBERS["doc.objects"]``, ``MEMBERS["objs.data"]``, ``MEMBERS["objs.set_data"]``).
4. **MIME spellings** - ``MIME`` / ``MIME_ALTERNATIVES`` (e.g. ``chemical/x-mdl-molfile-v3000``,
   ``chemical/x-helm``, ``chemical/x-name``); ``selection_as`` reports which ones were refused.
5. **SaveAs** - ``Document.SaveAs(path)`` is expected to pick the format from the extension; the file
   is sniffed afterwards (:func:`_sniff_format`) and ``SaveAs(path, mime|name|ext)`` plus
   ``Objects.Data(mime)`` are tried when the content does not match (``EXPORT_MIME``).
6. **Selection** - ``Document.Selection`` must be an Objects-like collection with ``Count`` (or at
   least ``Data``); ``MEMBERS["doc.selection"]``, ``doc.select_all``, ``doc.deselect``.
7. **Properties** - ``Objects.Formula`` / ``MolecularWeight`` / ``ExactMass`` / ``ElementalAnalysis``
   (``MEMBERS["objs.*"]``); ``selection_info`` lists what was ``unavailable``.
8. **Settings** - the ``Document.Settings`` object is optional: without it settings are read from the
   document CDXML (``CDXML_SETTINGS``) and written back only when SetData is known to replace.
9. **Dialogs** - alerts are expected to be ``#32770`` windows of ChemDraw.exe (``DIALOG_CLASSES``);
   ``pending_dialogs()`` shows what the enumerator sees.
10. **ChemScript** - the optional ``ChemScriptNN`` module only imports in the interpreter ChemDraw
    installed it for (``CHEMSCRIPT_PATH`` / ``CHEMSCRIPT_MODULE``; ``CHEMSCRIPT_MODULES``).

Public API
----------
``WindowsBackend`` (implements :class:`chemdrawmcp.backends.base.ChemDrawBackend`),
``ChemScriptBridge`` / ``chemscript_available()``, ``find_chemdraw_windows()``, ``find_stationery()``,
``sequence_to_helm()``, ``caption_cdxml()``, ``cdxml_summary()``, ``cdxml_object_signatures()``, ``merge_cdxml()``,
``cdxml_settings()``, ``cdxml_with_settings()``; the editable tables ``PROG_IDS``, ``MIME``, ``MIME_ALTERNATIVES``, ``INSERT_MIME``,
``SELECTION_FORMATS``, ``EXPORT_MIME``, ``MEMBERS``, ``SETTINGS_MEMBERS``, ``CDXML_SETTINGS``, ``COMMANDS``,
``COMMAND_ALIASES``, ``CHEMSCRIPT_MODULES``, ``CHEMSCRIPT_MEMBERS``, ``INSERT_MODE``, ``INSERT_TARGET``,
``DIALOG_CLASSES``, ``DEFAULT_DIALOG_BUTTONS``.
"""
from __future__ import annotations

import csv
import ctypes
import functools
import glob
import importlib
import io
import os
import queue
import re
import subprocess
import sys
import tempfile
import threading
import time
import types
from concurrent.futures import Future
from concurrent.futures import TimeoutError as _FutureTimeout
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as _xml_escape

from .base import (
    EXPORT_FORMATS,
    INSERT_KINDS,
    ChemDrawBackend,
    ChemDrawError,
    ChemDrawNotAvailable,
    CommandInfo,
    DocumentInfo,
    normalize_format,
)

# =====================================================================================
# Editable tables - everything uncertain about ChemDraw's COM / ChemScript API lives here
# =====================================================================================

#: COM ProgIDs tried in order (``GetActiveObject`` first, then ``Dispatch``/``CreateObject``).
PROG_IDS: List[str] = ["ChemDraw.Application"] + ["ChemDraw.Application.%d" % v for v in range(25, 11, -1)]

#: ProgIDs / file classes used only to *locate* ChemDraw.exe in the registry (never connected to):
#: the document class behind .cdx/.cdxml files and the ActiveX control.
FILE_PROG_IDS: List[str] = ["ChemDraw.Document", "ChemDraw.Document.6.0", "ChemDraw.Document.1",
                            "ChemDrawControl.ChemDrawCtl", "ChemDrawCtl.ChemDrawCtrl", "ChemDrawCtl.ChemDrawCtl"]
FILE_EXTENSIONS: Tuple[str, ...] = (".cdxml", ".cdx")

#: Process image names (lower case) that count as "ChemDraw is running".
PROCESS_NAMES: Tuple[str, ...] = ("chemdraw.exe", "chemdraw64.exe", "chemdrawprofessional.exe", "chemdrawprime.exe")

#: MIME types ChemDraw understands in ``Objects.Data(mime)`` / ``Objects.SetData(mime, text)``.
MIME: Dict[str, str] = {
    "cdxml": "chemical/x-cdxml",
    "cdx": "chemical/x-cdx",
    "smiles": "chemical/x-daylight-smiles",
    "smarts": "chemical/x-daylight-smarts",
    "inchi": "chemical/x-inchi",
    "inchikey": "chemical/x-inchikey",
    "mol": "chemical/x-mdl-molfile",
    "molv3000": "chemical/x-mdl-molfile-v3000",
    "sdf": "chemical/x-mdl-sdfile",
    "sdfv3000": "chemical/x-mdl-sdfile-v3000",
    "rxn": "chemical/x-mdl-rxnfile",
    "rxnv3000": "chemical/x-mdl-rxnfile-v3000",
    "name": "chemical/x-name",
    "helm": "chemical/x-helm",
    "sln": "chemical/x-sln",
    "cml": "chemical/x-cml",
    "png": "image/png",
    "emf": "image/x-emf",
    "text": "text/plain",
}

#: Alternative spellings tried (in order) after the primary MIME type above was refused.  Keys are
#: the primary strings from ``MIME``; add the spelling your ChemDraw accepts to the front.
MIME_ALTERNATIVES: Dict[str, List[str]] = {
    MIME["molv3000"]: ["chemical/x-mdl-molfile-V3000", "chemical/x-mdl-molfile3000"],
    MIME["sdfv3000"]: ["chemical/x-mdl-sdfile-V3000", "chemical/x-mdl-sdfile3000", "chemical/x-mdl-sdf-v3000"],
    MIME["rxnv3000"]: ["chemical/x-mdl-rxnfile-V3000", "chemical/x-mdl-rxnfile3000"],
    MIME["inchikey"]: ["chemical/x-inchi-key", "chemical/x-inchikey-std"],
    MIME["helm"]: ["chemical/x-HELM", "chemical/x-helm-notation"],
    MIME["name"]: ["chemical/x-chemical-name", "chemical/x-iupac-name"],
    MIME["sln"]: ["chemical/x-tripos-sln"],
    MIME["emf"]: ["image/emf", "application/x-emf", "image/x-wmf"],
    MIME["sdf"]: ["chemical/x-mdl-sdf"],
    MIME["smarts"]: ["chemical/x-daylight-smarts-query"],
}

#: ``insert(kind=...)`` -> MIME type passed to ``SetData``.  Sequence kinds are converted to HELM
#: first (see :func:`sequence_to_helm`); ``text`` becomes a CDXML caption (:func:`caption_cdxml`).
INSERT_MIME: Dict[str, str] = {
    "cdxml": MIME["cdxml"],
    "smiles": MIME["smiles"],
    "inchi": MIME["inchi"],
    "mol": MIME["mol"],
    "name": MIME["name"],
    "helm": MIME["helm"],
    "sequence": MIME["helm"],
    "peptide": MIME["helm"],
    "dna": MIME["helm"],
    "rna": MIME["helm"],
    "text": MIME["cdxml"],
}
SEQUENCE_KINDS: Tuple[str, ...] = ("sequence", "peptide", "dna", "rna")
#: Raw-sequence MIME types tried when pasting HELM is refused (names are guesses).
SEQUENCE_FALLBACK_MIME: Dict[str, List[str]] = {
    "peptide": ["chemical/x-peptide-sequence", "chemical/x-sequence", "chemical/x-fasta"],
    "dna": ["chemical/x-dna-sequence", "chemical/x-sequence", "chemical/x-fasta"],
    "rna": ["chemical/x-rna-sequence", "chemical/x-sequence", "chemical/x-fasta"],
    "sequence": ["chemical/x-sequence", "chemical/x-fasta"],
}

#: ``selection_as(what=...)`` -> MIME type passed to ``Data``.
SELECTION_FORMATS: Dict[str, str] = {
    "smiles": MIME["smiles"],
    "inchi": MIME["inchi"],
    "inchikey": MIME["inchikey"],
    "mol": MIME["mol"],
    "molv3000": MIME["molv3000"],
    "cdxml": MIME["cdxml"],
    "helm": MIME["helm"],
    "sln": MIME["sln"],
    "name": MIME["name"],
}

#: ``save_document(fmt=...)`` -> MIME type tried as the second ``SaveAs`` argument and for the
#: ``Objects.Data`` fallback.  Defaults to ``base.EXPORT_FORMATS``; override single entries here.
EXPORT_MIME: Dict[str, str] = {fmt: spec["mime"] for fmt, spec in EXPORT_FORMATS.items()}
# base.EXPORT_FORMATS gives each V2000/V3000 pair one MIME type; the V3000 exports need their own spelling.
EXPORT_MIME.update({"molv3000": MIME["molv3000"], "sdfv3000": MIME["sdfv3000"], "rxnv3000": MIME["rxnv3000"]})

#: How ``Objects.SetData`` behaves: "auto" (detect: counts went down => it replaced the document and
#: the old content is merged back in), "add" (trust it to add), "replace" (always merge ourselves).
INSERT_MODE = "auto"
#: Which collection ``insert`` pastes into: "objects" (whole document) or "selection".
INSERT_TARGET = "objects"

#: Candidate COM member names, tried in order.  Keys are "<object>.<concept>".  COM name lookup is
#: case-insensitive for late-bound (IDispatch) access, so spelling variants only matter for
#: early-bound wrappers; *different* names are what the lists are for.
MEMBERS: Dict[str, List[str]] = {
    # Application
    "app.documents": ["Documents"],
    "app.active_document": ["ActiveDocument"],
    "app.version": ["Version", "AppVersion", "ApplicationVersion"],
    "app.visible": ["Visible"],
    "app.quit": ["Quit", "Exit"],
    "app.activate": ["Activate", "BringToFront"],
    "app.window_state": ["WindowState"],
    "app.path": ["FullName", "Path"],
    "app.hwnd": ["HWnd", "hWnd", "Hwnd", "WindowHandle"],
    # Documents collection / any collection
    "coll.count": ["Count", "Length"],
    "coll.item": ["Item", "item", "Items"],
    "docs.add": ["Add", "New", "AddNew"],
    "docs.open": ["Open", "OpenDocument"],
    # Document
    "doc.name": ["Name", "Title"],
    "doc.full_name": ["FullName", "FullPath", "FileName", "File"],   # the complete path of a saved document
    "doc.dir": ["Path"],                                              # Office semantics: the *folder* only
    "doc.saved": ["Saved"],
    "doc.modified": ["Modified", "Dirty", "IsModified"],
    "doc.save": ["Save"],
    "doc.save_as": ["SaveAs", "ExportAs", "Export"],
    "doc.close": ["Close"],
    "doc.activate": ["Activate"],
    "doc.objects": ["Objects", "Content", "Contents"],
    "doc.selection": ["Selection", "SelectedObjects"],
    "doc.select_all": ["SelectAll"],
    "doc.deselect": ["Deselect", "DeselectAll", "ClearSelection", "UnselectAll", "SelectNone"],
    "doc.paste": ["Paste"],
    "doc.settings": ["Settings", "DocumentSettings", "DrawingSettings", "Style"],
    "doc.window": ["Window", "ActiveWindow"],
    "doc.hwnd": ["HWnd", "hWnd", "Hwnd", "WindowHandle"],
    "win.hwnd": ["HWnd", "hWnd", "Hwnd", "Handle", "WindowHandle"],
    # Objects / Selection collections (ChemDrawObjects) - also tried on the Document itself
    "objs.data": ["Data", "GetData"],
    "objs.set_data": ["SetData", "PutData"],
    "objs.clean": ["Clean", "CleanUp", "CleanUpStructure"],
    "objs.select": ["Select", "SelectAll"],
    "objs.selected": ["Selected"],
    "objs.paste": ["Paste"],
    "objs.formula": ["Formula", "MolecularFormula"],
    "objs.molecular_weight": ["MolecularWeight", "MolWeight", "MW"],
    "objs.exact_mass": ["ExactMass", "MonoisotopicMass"],
    "objs.elemental_analysis": ["ElementalAnalysis", "Analysis"],
    "objs.name": ["ChemicalName", "Name", "IUPACName"],
    "objs.atoms": ["Atoms"],
    "objs.bonds": ["Bonds"],
    "objs.fragments": ["Fragments", "Molecules", "Structures"],
    "objs.captions": ["Captions", "Texts", "TextObjects"],
    "objs.arrows": ["Arrows", "Graphics"],
    "objs.bounds": ["BoundingBox", "Bounds", "BoundingRect", "Rect"],
    # individual objects
    "obj.text": ["Text", "Caption", "Value"],
    # bounds objects
    "rect.left": ["Left", "left", "x0", "X0"],
    "rect.top": ["Top", "top", "y0", "Y0"],
    "rect.right": ["Right", "right", "x1", "X1"],
    "rect.bottom": ["Bottom", "bottom", "y1", "Y1"],
}

#: Document settings keys (same names as ``chemdrawmcp.styles``) -> candidate members of the
#: ``Document.Settings`` object (only used when such an object exists).
SETTINGS_MEMBERS: Dict[str, List[str]] = {
    "bond_length": ["BondLength", "FixedLength"],
    "line_width": ["LineWidth"],
    "bold_width": ["BoldWidth"],
    "margin_width": ["MarginWidth"],
    "hash_spacing": ["HashSpacing"],
    "bond_spacing": ["BondSpacing"],
    "chain_angle": ["ChainAngle"],
    "label_font": ["LabelFont", "LabelFontName", "AtomLabelFont"],
    "label_size": ["LabelSize", "LabelFontSize", "AtomLabelSize"],
    "label_face": ["LabelFace", "LabelFontStyle"],
    "caption_font": ["CaptionFont", "CaptionFontName", "TextFont"],
    "caption_size": ["CaptionSize", "CaptionFontSize", "TextSize"],
    "caption_face": ["CaptionFace", "CaptionFontStyle"],
    "magnification": ["Magnification", "Zoom", "ZoomPercent"],
}
#: The same settings as attributes of the ``<CDXML>`` root element (points / percent / degrees -
#: the units ``document_settings`` returns on macOS too).  Font attributes hold font-table ids.
CDXML_SETTINGS: Dict[str, str] = {
    "bond_length": "BondLength", "line_width": "LineWidth", "bold_width": "BoldWidth",
    "margin_width": "MarginWidth", "hash_spacing": "HashSpacing", "bond_spacing": "BondSpacing",
    "chain_angle": "ChainAngle", "label_font": "LabelFont", "label_size": "LabelSize",
    "label_face": "LabelFace", "caption_font": "CaptionFont", "caption_size": "CaptionSize",
    "caption_face": "CaptionFace", "magnification": "Magnification",
}
_NUMERIC_SETTINGS = ("bond_length", "line_width", "bold_width", "margin_width", "hash_spacing", "bond_spacing",
                     "chain_angle", "label_size", "caption_size", "magnification", "label_face", "caption_face")

#: Registry keys (under HKLM and HKCU, 64- and 32-bit views) searched for the install directory.
REGISTRY_KEYS: List[str] = [
    r"SOFTWARE\Revvity Signals\ChemDraw",
    r"SOFTWARE\Revvity Signals\ChemOffice",
    r"SOFTWARE\PerkinElmer\ChemDraw",
    r"SOFTWARE\PerkinElmerInformatics\ChemDraw",
    r"SOFTWARE\PerkinElmerInformatics\ChemOffice",
    r"SOFTWARE\CambridgeSoft\ChemDraw",
    r"SOFTWARE\CambridgeSoft\ChemOffice",
]
APP_PATHS_KEY = r"SOFTWARE\Microsoft\Windows\CurrentVersion\App Paths\ChemDraw.exe"
REGISTRY_PATH_VALUES: List[str] = ["", "Path", "InstallPath", "InstallDir", "Install Directory", "ExePath",
                                   "ApplicationPath", "Location", "Dir", "Directory"]

#: Install-directory globs (relative to every Program Files variant) for ChemDraw.exe.
EXE_GLOB_SUFFIXES: List[str] = [
    r"Revvity Signals\ChemDraw*\ChemDraw.exe",
    r"Revvity Signals\ChemOffice*\ChemDraw\ChemDraw.exe",
    r"Revvity Signals\ChemOffice*\ChemDraw*\ChemDraw.exe",
    r"PerkinElmerInformatics\ChemDraw*\ChemDraw.exe",
    r"PerkinElmerInformatics\ChemOffice*\ChemDraw\ChemDraw.exe",
    r"PerkinElmerInformatics\ChemOffice*\ChemDraw*\ChemDraw.exe",
    r"PerkinElmer\ChemDraw*\ChemDraw.exe",
    r"PerkinElmer\ChemOffice*\ChemDraw\ChemDraw.exe",
    r"CambridgeSoft\ChemDraw*\ChemDraw.exe",
    r"CambridgeSoft\ChemOffice*\ChemDraw\ChemDraw.exe",
    r"CambridgeSoft\ChemOffice*\ChemDraw*\ChemDraw.exe",
    r"ChemDraw*\ChemDraw.exe",
]
#: Where stationery (``*.cds``) files live, relative to Program Files / ProgramData / APPDATA.
STATIONERY_GLOB_SUFFIXES: List[str] = [
    r"Revvity Signals\ChemDraw*\ChemDraw Items",
    r"Revvity Signals\ChemOffice*\ChemDraw\ChemDraw Items",
    r"PerkinElmerInformatics\ChemOffice*\ChemDraw\ChemDraw Items",
    r"PerkinElmerInformatics\ChemDraw*\ChemDraw Items",
    r"CambridgeSoft\ChemOffice*\ChemDraw\ChemDraw Items",
    r"Revvity Signals\ChemDraw\ChemDraw Items",
    r"PerkinElmer\ChemDraw\ChemDraw Items",
    r"ChemDraw*\ChemDraw Items",
]

#: ChemScript module names tried in order (``CHEMSCRIPT_MODULE`` env var is tried first).
CHEMSCRIPT_MODULES: List[str] = ["ChemScript%d" % v for v in range(25, 13, -1)] + ["ChemScript"]
#: Directories (globs) appended to ``sys.path`` before importing ChemScript (``CHEMSCRIPT_PATH`` first).
CHEMSCRIPT_GLOB_SUFFIXES: List[str] = [
    r"Revvity Signals\ChemDraw*\ChemScript\Python*",
    r"Revvity Signals\ChemOffice*\ChemScript\Python*",
    r"PerkinElmerInformatics\ChemOffice*\ChemScript\Python*",
    r"PerkinElmerInformatics\ChemDraw*\ChemScript\Python*",
    r"CambridgeSoft\ChemOffice*\ChemScript\Python*",
    r"ChemScript*\Python*",
]
CHEMSCRIPT_PYTHON_DIRS: List[str] = [r"C:\Python39\Lib\site-packages", r"C:\Python38\Lib\site-packages",
                                     r"C:\Python37\Lib\site-packages", r"C:\Python27\Lib\site-packages"]
#: Candidate ChemScript ``StructureData`` members (``ChemicalName()`` is Struct=Name; ``Name`` is the
#: molfile title, hence second).
CHEMSCRIPT_MEMBERS: Dict[str, List[str]] = {
    "load_data": ["LoadData"],
    "load_file": ["LoadFile"],
    "write_data": ["WriteData"],
    "write_file": ["WriteFile"],
    "name": ["ChemicalName", "Name", "GetName"],
    "clean": ["CleanUp", "Clean", "CleanUp2D"],
    "smiles": ["Smiles", "SMILES"],
    "inchi": ["InChI", "Inchi", "InChi"],
    "formula": ["Formula", "MolecularFormula"],
    "molecular_weight": ["MolecularWeight", "MolWeight"],
    "exact_mass": ["ExactMass", "MonoisotopicMass"],
}

#: Curated commands available through COM (``list_commands`` / ``run_command``).
COMMANDS: Dict[str, Tuple[str, str]] = {
    "newDocument": ("Create a new, empty document", "Documents.Add(); the new document becomes the front document."),
    "closeDocument": ("Close the document without saving", "Document.Close(False) - unsaved changes are discarded."),
    "saveDocument": ("Save the document to its current path", "Document.Save(); fails for untitled documents (use save_document)."),
    "selectAll": ("Select everything in the document", "Document.SelectAll() / Objects.Select() / Objects.Selected = True."),
    "deselectAll": ("Clear the selection", "Best effort: not every ChemDraw version exposes this through COM."),
    "cleanStructure": ("Clean up the selected structures", "Selection.Clean() (or Objects.Clean() when nothing is selected)."),
    "convertNameToStructure": ("Name > Structure for the selected caption text",
                               "Reads the selected caption(s) and pastes chemical/x-name; use insert(kind='name') when nothing is selected."),
    "convertStructureToName": ("Structure > Name for the selection",
                               "Objects.Data('chemical/x-name'); the name is returned and also inserted as a caption."),
    "copyAsSMILES": ("Selection as SMILES", "Returned in `result`; nothing is put on the clipboard."),
    "copyAsInChI": ("Selection as InChI", "Returned in `result`."),
    "copyAsInChIKey": ("Selection as InChIKey", "Returned in `result` (RDKit fallback when ChemDraw cannot)."),
    "copyAsMolfile": ("Selection as MDL molfile (V2000)", "Returned in `result`."),
    "copyAsMolfileV3000": ("Selection as MDL molfile (V3000)", "Returned in `result`."),
    "copyAsCDXML": ("Selection as CDXML", "Returned in `result`."),
    "copyAsHELM": ("Selection as HELM", "Returned in `result`."),
    "copyAsSLN": ("Selection as SLN", "Returned in `result`."),
    "chemicalProperties": ("Formula, MW, exact mass, elemental analysis of the selection", "Same data as selection_info."),
    "elementalAnalysis": ("Elemental analysis of the selection", "Objects.ElementalAnalysis."),
    "dismissDialogs": ("Close any modal dialog ChemDraw is showing", "Returns the dialog texts in `result`."),
}
#: macOS command names (see backends/macos.py) accepted as synonyms by ``run_command``.
COMMAND_ALIASES: Dict[str, str] = {
    "copyAsMol": "copyAsMolfile", "copyAsMolV3000": "copyAsMolfileV3000", "copyAsCdxml": "copyAsCDXML",
    "copyAsInChi": "copyAsInChI", "cleanUpStructure": "cleanStructure", "clean": "cleanStructure",
    "structureToName": "convertStructureToName", "nameToStructure": "convertNameToStructure",
    "close": "closeDocument", "save": "saveDocument", "new": "newDocument",
}

#: Win32 window classes treated as modal dialogs (plus any class containing "dialog").
DIALOG_CLASSES: Tuple[str, ...] = ("#32770",)
#: Button labels pressed (first match wins) when no preference is given; "Yes" last because answering
#: "Save changes?" with Yes opens a file dialog.
DEFAULT_DIALOG_BUTTONS: Tuple[str, ...] = ("OK", "Don't Save", "Don\u2019t Save", "No", "Discard", "Close", "Cancel", "Yes")

# Win32 constants used by the dialog / screenshot / clipboard helpers.
_WM_CLOSE = 0x0010
_WM_GETTEXT = 0x000D
_WM_GETTEXTLENGTH = 0x000E
_WM_COMMAND = 0x0111
_DM_GETDEFID = 0x0400
_DC_HASDEFID = 0x534B
_BM_CLICK = 0x00F5
_IDOK = 1
_IDCANCEL = 2
_SMTO_ABORTIFHUNG = 0x0002
_SW_RESTORE = 9
_GW_OWNER = 4
_CF_UNICODETEXT = 13
_GMEM_MOVEABLE = 0x0002
_CREATE_NO_WINDOW = 0x08000000
#: RPC_S_SERVER_UNAVAILABLE, RPC_S_CALL_FAILED, RPC_E_DISCONNECTED, RPC_S_CALL_FAILED_DNE: ChemDraw went away.
_DISCONNECTED_HRESULTS = {0x800706BA, 0x800706BE, 0x80010108, 0x800706BF}
#: DISP_E_UNKNOWNNAME, DISP_E_MEMBERNOTFOUND: the member does not exist (try the next candidate).
_NOT_FOUND_HRESULTS = {0x80020006, 0x80020003}
#: RPC_E_WRONG_THREAD: a proxy was used from the wrong apartment (should not happen with the worker).
_WRONG_THREAD_HRESULT = 0x8001010E

_MISSING = object()


# =====================================================================================
# Generic helpers (platform independent)
# =====================================================================================

def _is_windows() -> bool:
    return sys.platform == "win32"


def _is_com_error(exc: BaseException) -> bool:
    name = type(exc).__name__
    if name in ("com_error", "COMError"):
        return True
    for modname, attr in (("pywintypes", "com_error"), ("comtypes", "COMError")):
        mod = sys.modules.get(modname)
        cls = getattr(mod, attr, None) if mod else None
        if cls is not None and isinstance(exc, cls):
            return True
    return False


def _hresult(exc: BaseException) -> Optional[int]:
    h = getattr(exc, "hresult", None)
    if h is None:
        args = getattr(exc, "args", ())
        if args and isinstance(args[0], int):
            h = args[0]
    if h is None:
        return None
    return h & 0xFFFFFFFF


def _member_missing(exc: BaseException) -> bool:
    """True when an exception means "this member does not exist" (try the next candidate)."""
    if isinstance(exc, AttributeError):
        return True
    return _is_com_error(exc) and _hresult(exc) in _NOT_FOUND_HRESULTS


def _describe_com_error(exc: BaseException) -> str:
    """Human readable text for a pywin32 ``com_error`` or comtypes ``COMError``."""
    args = tuple(getattr(exc, "args", ()))
    hres = _hresult(exc)
    text = getattr(exc, "strerror", None) or (args[1] if len(args) > 1 and isinstance(args[1], str) else None)
    source = desc = None
    excepinfo = getattr(exc, "excepinfo", None)
    if excepinfo is None and len(args) > 2:
        excepinfo = args[2]
    if isinstance(excepinfo, (tuple, list)) and len(excepinfo) >= 3:      # pywin32
        source, desc = excepinfo[1], excepinfo[2]
    details = getattr(exc, "details", None)                                 # comtypes
    if isinstance(details, (tuple, list)) and details:
        desc = details[0] or desc
        if len(details) > 1:
            source = details[1] or source
    parts = []
    if hres is not None:
        parts.append("0x%08X" % hres)
    if desc:
        parts.append(str(desc).strip())
    if text and (not desc or str(text).strip() not in str(desc)):
        parts.append(str(text).strip())
    if not parts:
        parts.append(str(exc))
    msg = "COM error " + ": ".join(parts)
    if source:
        msg += " (source: %s)" % source
    if hres == _WRONG_THREAD_HRESULT:
        msg += " - the COM proxy was used from another thread; keep com_thread=True (the default)"
    return msg


def _exc_text(exc: BaseException) -> str:
    return _describe_com_error(exc) if _is_com_error(exc) else "%s: %s" % (type(exc).__name__, exc)


#: Type names of the callable wrappers the COM client libraries hand out for methods and for
#: parameterised properties (comtypes ``bound_named_property`` is ``Data[mime]`` *and* ``Data(mime)``).
_METHOD_TYPE_NAMES = frozenset((
    "MethodCaller", "_MethodCaller", "NamedProperty", "bound_named_property", "named_property",
    "method-wrapper", "_Method", "LazyMethod",
))


def _is_method(val: Any) -> bool:
    """True for a bound COM/Python method (vs. a property value, which may itself be a callable COM object)."""
    if isinstance(val, (types.MethodType, types.FunctionType, types.BuiltinMethodType, types.BuiltinFunctionType)):
        return True
    return type(val).__name__ in _METHOD_TYPE_NAMES


def _raw_invoke(obj: Any, name: str, args: Sequence[Any] = (), kind: str = "call") -> Any:
    """Invoke a member by name straight through ``IDispatch::GetIDsOfNames`` + ``Invoke``.

    Used when the client library's guess (property vs. method) fails, e.g. comtypes raising
    DISP_E_BADPARAMCOUNT for ``Data(mime)``.  ``kind`` is "call" (DISPATCH_METHOD|PROPERTYGET) or "put".
    Raises AttributeError for objects that are not COM dispatch wrappers (the test fakes).
    """
    ole = getattr(obj, "_oleobj_", None)                       # pywin32 (dynamic and makepy wrappers)
    if ole is not None:
        import pythoncom  # type: ignore

        try:
            dispid = ole.GetIDsOfNames(0, name)
        except TypeError:
            dispid = ole.GetIDsOfNames(name)
        if isinstance(dispid, (tuple, list)):
            dispid = dispid[0]
        if kind == "put":
            return ole.Invoke(dispid, 0, pythoncom.DISPATCH_PROPERTYPUT, 0, *args)
        return ole.Invoke(dispid, 0, pythoncom.DISPATCH_METHOD | pythoncom.DISPATCH_PROPERTYGET, 1, *args)
    com = getattr(obj, "_comobj", None)                        # comtypes dynamic Dispatch
    if com is None and hasattr(obj, "GetIDsOfNames") and hasattr(obj, "Invoke"):
        com = obj                                              # comtypes early-bound IDispatch pointer
    if com is not None:                                        # DISPATCH_METHOD=1, PROPERTYGET=2, PROPERTYPUT=4
        dispid = com.GetIDsOfNames(name)
        if isinstance(dispid, (tuple, list)):
            dispid = dispid[0]
        if kind == "put":
            return com.Invoke(dispid, *args, _invkind=4)
        return com.Invoke(dispid, *args, _invkind=1 | 2)
    raise AttributeError("%s is not a COM dispatch object (no IDispatch access for %r)" % (type(obj).__name__, name))


def _as_bytes(data: Any) -> bytes:
    if data is None:
        return b""
    if isinstance(data, (bytes, bytearray, memoryview)):
        return bytes(data)
    if isinstance(data, (tuple, list)) and data and all(isinstance(x, int) for x in data):
        return bytes(x & 0xFF for x in data)        # SAFEARRAY of VT_UI1 delivered as a tuple
    if isinstance(data, str):
        return data.encode("utf-8")
    return str(data).encode("utf-8")


def _as_text(data: Any) -> str:
    if data is None:
        return ""
    if isinstance(data, str):
        return data
    raw = _as_bytes(data)
    for enc in ("utf-8-sig", "utf-16"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1")


_FLOAT_RE = re.compile(r"[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?")


def _num(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, bool):
        return float(value)
    if isinstance(value, (int, float)):
        return float(value)
    m = _FLOAT_RE.search(str(value))
    return float(m.group(0).replace(",", ".")) if m else None


def _mime_candidates(mime: str) -> List[str]:
    out = [mime]
    for alt in MIME_ALTERNATIVES.get(mime, []):
        if alt not in out:
            out.append(alt)
    return out


def _ensure_extension(path: str, ext: str) -> str:
    """Absolute path whose extension is ``ext`` (replacing a known chemistry/image extension)."""
    path = os.path.abspath(os.path.expanduser(path))
    root, cur = os.path.splitext(path)
    if cur.lower() == ext.lower():
        return path
    known = {spec["ext"].lower() for spec in EXPORT_FORMATS.values()} | {".jpg", ".tif", ".xml"}
    if cur.lower() in known:
        return root + ext
    return path + ext


#: Leading bytes / text markers used to check that ``SaveAs`` really wrote the requested format.
_FORMAT_MAGIC: Dict[str, Tuple[bytes, ...]] = {
    "cdx": (b"VjCD0100",), "png": (b"\x89PNG",), "gif": (b"GIF8",), "jpeg": (b"\xff\xd8",), "bmp": (b"BM",),
    "tiff": (b"II*\x00", b"MM\x00*"), "pdf": (b"%PDF",), "eps": (b"%!PS", b"\xc5\xd0\xd3\xc6"), "3mf": (b"PK",),
}
#: The V2000 and V3000 molfile/SDfile/RXN variants are told apart on purpose (ChemDraw 23's plain "MDL
#: Molfile" writes V3000): a wrong variant makes ``save_document`` try the other ``SaveAs`` shapes and
#: ``Objects.Data``.
_FORMAT_MARKERS: Dict[str, Tuple[str, ...]] = {
    "cdxml": ("<CDXML",), "cml": ("<cml", "<molecule"), "svg": ("<svg",),
    "mol": ("V2000",), "molv3000": ("V3000",), "sdf": ("V2000",), "sdfv3000": ("V3000",),
    "rxn": ("$RXN\n", "$RXN\r\n", "$RXN\r"), "rxnv3000": ("$RXN V3000",), "rdf": ("$RDFILE",),
}


def _sniff_format(path: str, fmt: str) -> Optional[bool]:
    """True/False when the file's head does/does not look like ``fmt``; None when the format has no signature."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(4096)
    except OSError:
        return False
    if not head:
        return False
    if fmt in _FORMAT_MAGIC:
        return any(head.startswith(m) for m in _FORMAT_MAGIC[fmt])
    if fmt in _FORMAT_MARKERS:
        if head[:2] in (b"\xff\xfe", b"\xfe\xff") or (len(head) > 3 and head[1:2] == b"\x00"):
            text = head.decode("utf-16", errors="replace")
        else:
            text = head.decode("latin-1")
        return any(m in text for m in _FORMAT_MARKERS[fmt])
    return None


def _version_key(path: str) -> Tuple:
    nums = [int(x) for x in re.findall(r"\d+", path)]
    return tuple(nums), path


def _program_dirs() -> List[str]:
    dirs: List[str] = []
    for var in ("ProgramFiles", "ProgramW6432", "ProgramFiles(x86)", "ProgramData", "LOCALAPPDATA"):
        v = os.environ.get(var)
        if v:
            dirs.append(v)
    dirs += [r"C:\Program Files", r"C:\Program Files (x86)", r"C:\ProgramData"]
    seen: List[str] = []
    for d in dirs:
        if d.lower() not in [s.lower() for s in seen]:
            seen.append(d)
    return seen


def _glob_many(patterns: Iterable[str]) -> List[str]:
    out: List[str] = []
    for p in patterns:
        try:
            out.extend(glob.glob(p))
        except Exception:
            continue
    return sorted(set(out), key=_version_key, reverse=True)


def _exe_from_command(cmd: str) -> Optional[str]:
    cmd = (cmd or "").strip()
    if not cmd:
        return None
    m = re.match(r'^"([^"]+)"', cmd)
    if m:
        return m.group(1)
    m = re.match(r"^(.+?\.exe)", cmd, re.IGNORECASE)
    return m.group(1) if m else cmd


# =====================================================================================
# Chemistry text helpers (pure Python)
# =====================================================================================

_AMINO_ONE = set("ACDEFGHIKLMNPQRSTVWY")
_AMINO_THREE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _strip_fasta(seq: str) -> str:
    lines = [ln.strip() for ln in (seq or "").splitlines()]
    return "".join(ln for ln in lines if ln and not ln.startswith((">", ";")))


def sequence_to_helm(seq: str, kind: str = "sequence") -> str:
    """Convert a peptide / DNA / RNA sequence to HELM notation (what ChemDraw pastes as chemical/x-helm).

    ``kind`` is one of ``peptide``, ``dna``, ``rna`` or ``sequence`` (auto-detect: six or more
    letters all in ACGTUN are nucleic acid - RNA when a U occurs - otherwise a peptide).
    Peptides accept one-letter strings, three-letter codes separated by ``-``/space/comma, and
    ``[Xxx]`` bracketed non-standard monomers; FASTA headers are ignored.
    """
    kind = (kind or "sequence").lower()
    text = _strip_fasta(seq) if kind != "peptide" or "\n" in (seq or "") else (seq or "").strip()
    if not text:
        raise ChemDrawError("empty sequence")
    if kind == "sequence":
        letters = re.sub(r"[^A-Za-z]", "", text).upper()
        if len(letters) >= 6 and set(letters) <= set("ACGTUN"):
            kind = "rna" if "U" in letters else "dna"
        else:
            kind = "peptide"
    if kind == "peptide":
        stripped = text.strip()
        tokens: List[str] = []
        if re.search(r"[-\s,]", stripped):
            raw = [t for t in re.split(r"[-\s,]+", stripped) if t]
            if all(t.upper() in _AMINO_THREE or t.startswith("[") or len(t) == 1 for t in raw):
                tokens = [t if t.startswith("[") else (_AMINO_THREE[t.upper()] if len(t) == 3 else t.upper()) for t in raw]
            elif all(re.fullmatch(r"[A-Za-z]+", t) for t in raw) and set("".join(raw).upper()) <= _AMINO_ONE:
                tokens = list("".join(raw).upper())           # one-letter sequence written in blocks
            else:
                for t in raw:
                    if t.startswith("["):
                        tokens.append(t)
                    elif t.upper() in _AMINO_THREE:
                        tokens.append(_AMINO_THREE[t.upper()])
                    elif len(t) == 1:
                        tokens.append(t.upper())
                    else:
                        tokens.append("[%s]" % t)
        else:
            for tok in re.findall(r"\[[^\]]+\]|[A-Za-z]", stripped):
                tokens.append(tok if tok.startswith("[") else tok.upper())
        if not tokens:
            raise ChemDrawError("no residues in sequence %r" % seq)
        monomers = [t if t.startswith("[") or t in _AMINO_ONE else "[%s]" % t for t in tokens]
        return "PEPTIDE1{%s}$$$$" % ".".join(monomers)
    bases = re.sub(r"[^A-Za-z]", "", text).upper()
    if not bases:
        raise ChemDrawError("no bases in sequence %r" % seq)
    sugar = "[dR]" if kind == "dna" else "R"
    if kind == "dna":
        bases = bases.replace("U", "T")
    parts = ["%s(%s)" % (sugar, b) for b in bases]
    return "RNA1{%s}$$$$" % "P.".join(parts)


def _parse_cdxml(text: Any) -> Any:
    import xml.etree.ElementTree as ET

    raw = text.encode("utf-8") if isinstance(text, str) else bytes(text or b"")
    return ET.fromstring(raw)


def cdxml_summary(cdxml: str) -> Dict[str, Any]:
    """Counts, caption text and bounding box read from CDXML text (platform independent).

    Returns ``{"atoms", "bonds", "molecules", "captions", "arrows", "text", "bounds"}``; ``text`` joins
    the captions (``<t>`` elements that are not atom labels, i.e. not inside ``<n>``), ``bounds`` is the
    union of the ``BoundingBox`` attributes of the page's objects (CDXML points, y down) or None.
    """
    import xml.etree.ElementTree as ET

    out: Dict[str, Any] = {"atoms": 0, "bonds": 0, "molecules": 0, "captions": 0, "arrows": 0, "text": None, "bounds": None}
    try:
        root = _parse_cdxml(cdxml)
    except (ET.ParseError, ValueError, TypeError):
        return out
    parents = {child: parent for parent in root.iter() for child in parent}

    def inside_node(el: Any) -> bool:
        p = parents.get(el)
        while p is not None:
            if p.tag == "n":
                return True
            p = parents.get(p)
        return False

    texts: List[str] = []
    boxes: List[Tuple[float, float, float, float]] = []
    for el in root.iter():
        tag = el.tag
        if tag == "n":
            out["atoms"] += 1
        elif tag == "b":
            out["bonds"] += 1
        elif tag == "fragment" and not inside_node(el):
            out["molecules"] += 1
        elif tag == "t" and not inside_node(el):
            out["captions"] += 1
            t = "".join((s.text or "") for s in el.iter("s")) or (el.text or "")
            if t.strip():
                texts.append(t.strip())
        elif tag == "arrow" or (tag == "graphic" and el.get("ArrowType")):
            out["arrows"] += 1
        if tag in ("fragment", "t", "arrow", "graphic", "group") and not inside_node(el):
            bb = el.get("BoundingBox")
            if bb:
                nums = [_num(v) for v in bb.split()]
                if len(nums) == 4 and all(n is not None for n in nums):
                    boxes.append(tuple(nums))  # type: ignore[arg-type]
    if texts:
        out["text"] = "\n".join(texts)
    if boxes:
        l, t, r, b = (min(x[0] for x in boxes), min(x[1] for x in boxes), max(x[2] for x in boxes), max(x[3] for x in boxes))
        out["bounds"] = {"left": l, "top": t, "right": r, "bottom": b, "width": r - l, "height": b - t}
    return out


def _total_objects(counts: Dict[str, int]) -> int:
    return int(counts.get("atoms", 0)) + int(counts.get("captions", 0)) + int(counts.get("arrows", 0))


def cdxml_object_signatures(cdxml: str, with_box: bool = True) -> List[Tuple[Any, ...]]:
    """One id-independent signature per top-level page object: (tag, atoms, bonds, text[, rounded bounding box]).

    Used to tell whether ``SetData`` *added* to a document (every old signature is still present) or
    *replaced* it (some old object is gone), regardless of how many atoms the new content has.
    ``with_box=False`` ignores positions (ChemDraw may re-lay out or re-round boxes when pasting).
    """
    import xml.etree.ElementTree as ET

    try:
        root = _parse_cdxml(cdxml)
    except (ET.ParseError, ValueError, TypeError):
        return []
    out: List[Tuple[Any, ...]] = []
    for page in root.iter("page"):
        for el in page:
            atoms = sum(1 for _ in el.iter("n"))
            bonds = sum(1 for _ in el.iter("b"))
            text = "".join((s.text or "") for s in el.iter("s")).strip() if el.tag != "fragment" else ""
            bb = el.get("BoundingBox") if with_box else None
            box = None
            if bb:
                nums = [_num(v) for v in bb.split()]
                if len(nums) == 4 and None not in nums:
                    box = tuple(round(n, 1) for n in nums)  # type: ignore[arg-type]
            out.append((el.tag, atoms, bonds, text, box) if with_box else (el.tag, atoms, bonds, text))
    return out


def caption_cdxml(text: str, x: float = 100.0, y: float = 100.0, size: float = 10.0, font: str = "Arial") -> str:
    """Minimal CDXML document containing one caption (text object) for ``SetData('chemical/x-cdxml', ...)``."""
    body = _xml_escape(text or "")
    return (
        '<?xml version="1.0" encoding="UTF-8" ?>\n'
        '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n'
        '<CDXML CreationProgram="chemdrawmcp" BondLength="14.40" CaptionFont="3" CaptionSize="%g" '
        'LabelFont="3" LabelSize="%g">\n'
        '<fonttable><font id="3" charset="iso-8859-1" name="%s"/></fonttable>\n'
        '<page HeightPages="1" WidthPages="1">\n'
        '<t p="%g %g" LineHeight="auto"><s font="3" size="%g">%s</s></t>\n'
        '</page></CDXML>\n' % (size, size, _xml_escape(font), x, y, size, body)
    )


_CDXML_HEADER = ('<?xml version="1.0" encoding="UTF-8" ?>\n'
                 '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n')
#: CDXML attributes holding object ids (space separated) that must be renumbered when documents are merged.
_CDXML_ID_REF_ATTRS = ("B", "E", "BondCircularOrdering", "Objects", "BasisObjects", "AttachedAtoms", "ConnectionNums",
                       "SupersededBy", "ReactionStepReactants", "ReactionStepProducts", "ReactionStepArrows",
                       "ReactionStepPlusses", "ReactionStepObjectsAboveArrow", "ReactionStepObjectsBelowArrow",
                       "ReactionStepAtomMap", "ReactionStepAtomMapManual", "ReactionStepAtomMapAuto")
#: CDXML coordinate attributes: "x y" pairs, "l t r b" boxes and "x y z" triples (z untouched).
_CDXML_XY_ATTRS = ("p",)
_CDXML_BOX_ATTRS = ("BoundingBox",)
_CDXML_XYZ_ATTRS = ("xyz", "Center3D", "MajorAxisEnd3D", "MinorAxisEnd3D", "Head3D", "Tail3D")


def _shift_coords(el: Any, dx: float, dy: float) -> None:
    def fmt(v: float) -> str:
        return ("%.2f" % v).rstrip("0").rstrip(".") if abs(v) < 1e9 else str(v)

    for attr in _CDXML_XY_ATTRS:
        v = el.get(attr)
        if v:
            nums = [_num(x) for x in v.split()]
            if len(nums) == 2 and None not in nums:
                el.set(attr, "%s %s" % (fmt(nums[0] + dx), fmt(nums[1] + dy)))  # type: ignore[operator]
    for attr in _CDXML_BOX_ATTRS:
        v = el.get(attr)
        if v:
            nums = [_num(x) for x in v.split()]
            if len(nums) == 4 and None not in nums:
                el.set(attr, " ".join(fmt(n + (dx if i % 2 == 0 else dy)) for i, n in enumerate(nums)))  # type: ignore[operator]
    for attr in _CDXML_XYZ_ATTRS:
        v = el.get(attr)
        if v:
            nums = [_num(x) for x in v.split()]
            if len(nums) == 3 and None not in nums:
                el.set(attr, "%s %s %s" % (fmt(nums[0] + dx), fmt(nums[1] + dy), fmt(nums[2])))  # type: ignore[operator]


def merge_cdxml(base: str, addition: str, gap: float = 36.0) -> str:
    """Append the objects of CDXML document ``addition`` below those of ``base`` (one CDXML document out).

    Used when ``Objects.SetData`` turned out to *replace* the document: the previous content (``base``)
    is merged with what ChemDraw now holds (``addition``) and set again.  Object ids of the addition are
    renumbered above the base's highest id (bond ``B``/``E`` and the other id references follow), its
    fonts are mapped into the base font table, and its coordinates are shifted so that it starts ``gap``
    points below the base's bounding box, left-aligned with it.  Root attributes come from ``base``.
    """
    import xml.etree.ElementTree as ET

    if not (base or "").strip():
        return addition
    if not (addition or "").strip():
        return base
    broot = _parse_cdxml(base)
    aroot = _parse_cdxml(addition)
    apage = aroot.find("page")
    if apage is None or len(apage) == 0:
        return base
    bpage = broot.find("page")
    if bpage is None:
        bpage = ET.SubElement(broot, "page")
        bpage.set("HeightPages", "1")
        bpage.set("WidthPages", "1")
    # --- offset so the new objects do not overlap the old ones
    bsum, asum = cdxml_summary(base), cdxml_summary(addition)
    dx = dy = 0.0
    if bsum["bounds"] and asum["bounds"]:
        dx = bsum["bounds"]["left"] - asum["bounds"]["left"]
        dy = bsum["bounds"]["bottom"] + gap - asum["bounds"]["top"]
    # --- ids (font / color table entries carry ids of their own and are not objects)
    def objects_of(root: Any) -> List[Any]:
        return [el for el in root.iter() if el.tag not in ("fonttable", "font", "colortable", "color")]

    used = [int(el.get("id")) for el in objects_of(broot) if (el.get("id") or "").lstrip("-").isdigit()]
    next_id = (max(used) if used else 0) + 1
    mapping: Dict[str, str] = {}
    for el in objects_of(aroot):
        old = el.get("id")
        if old is not None and old not in mapping:
            mapping[old] = str(next_id)
            next_id += 1
    # --- fonts: map the addition's font ids onto the base font table (adding missing fonts)
    bft = broot.find("fonttable")
    aft = aroot.find("fonttable")
    font_map: Dict[str, str] = {}
    if aft is not None:
        if bft is None:
            bft = ET.Element("fonttable")
            broot.insert(0, bft)
        by_name = {(f.get("name") or "").lower(): f.get("id") for f in bft.findall("font")}
        font_ids = [int(f.get("id")) for f in bft.findall("font") if (f.get("id") or "").isdigit()]
        next_font = (max(font_ids) if font_ids else 2) + 1
        for f in aft.findall("font"):
            fid, name = f.get("id"), (f.get("name") or "")
            if fid is None:
                continue
            if name.lower() in by_name and by_name[name.lower()]:
                font_map[fid] = by_name[name.lower()]
            else:
                nf = ET.SubElement(bft, "font")
                for k, v in f.attrib.items():
                    nf.set(k, v)
                nf.set("id", str(next_font))
                by_name[name.lower()] = str(next_font)
                font_map[fid] = str(next_font)
                next_font += 1
    # --- rewrite the addition in place
    for el in objects_of(aroot):
        old = el.get("id")
        if old is not None:
            el.set("id", mapping[old])
        for attr in _CDXML_ID_REF_ATTRS:
            v = el.get(attr)
            if v:
                el.set(attr, " ".join(mapping.get(tok, tok) for tok in v.split()))
        if font_map:
            for attr in ("font", "LabelFont", "CaptionFont"):
                v = el.get(attr)
                if v in font_map:
                    el.set(attr, font_map[v])
        if dx or dy:
            _shift_coords(el, dx, dy)
    for child in list(apage):
        bpage.append(child)
    body = ET.tostring(broot, encoding="unicode")
    return _CDXML_HEADER + body + ("\n" if not body.endswith("\n") else "")


def cdxml_settings(cdxml: str) -> Dict[str, Any]:
    """Drawing settings from the ``<CDXML>`` root attributes (fonts resolved through the font table)."""
    import xml.etree.ElementTree as ET

    try:
        root = _parse_cdxml(cdxml)
    except (ET.ParseError, ValueError, TypeError):
        return {}
    if root.tag != "CDXML":
        return {}
    fonts = {f.get("id"): f.get("name") for f in root.iter("font")}
    out: Dict[str, Any] = {}
    for key, attr in CDXML_SETTINGS.items():
        v = root.get(attr)
        if v is None:
            continue
        if key in ("label_font", "caption_font"):
            out[key] = fonts.get(v, v)
        elif key in _NUMERIC_SETTINGS:
            n = _num(v)
            if n is not None:
                out[key] = int(n) if key.endswith("_face") else n
        else:
            out[key] = v
    return out


def cdxml_with_settings(cdxml: str, settings: Dict[str, Any]) -> str:
    """``cdxml`` with the given settings (keys of ``CDXML_SETTINGS``) written into the root element.

    Only the ``<CDXML ...>`` start tag (and, for new fonts, the font table) is rewritten, the rest of the
    text is kept byte for byte.  Unknown keys raise :class:`ChemDrawError`.
    """
    m = re.search(r"<CDXML\b([^>]*)>", cdxml)
    if not m:
        raise ChemDrawError("not a CDXML document (no <CDXML> root element)")
    attrs = m.group(1)
    fonts = {(f[1] or "").lower(): f[0] for f in re.findall(r'<font\s+id="(\d+)"[^>]*\bname="([^"]*)"', cdxml)}
    new_fonts: List[Tuple[str, str]] = []
    next_font = max([int(i) for i in fonts.values()] + [2]) + 1
    for key, value in settings.items():
        if key not in CDXML_SETTINGS:
            raise ChemDrawError("unknown setting %r; known: %s" % (key, ", ".join(CDXML_SETTINGS)))
        attr = CDXML_SETTINGS[key]
        if key in ("label_font", "caption_font"):
            name = str(value)
            if name.isdigit():
                sval = name
            elif name.lower() in fonts:
                sval = fonts[name.lower()]
            else:
                sval = str(next_font)
                fonts[name.lower()] = sval
                new_fonts.append((sval, name))
                next_font += 1
        elif isinstance(value, bool):
            sval = "yes" if value else "no"
        elif isinstance(value, (int, float)):
            sval = ("%.4f" % float(value)).rstrip("0").rstrip(".") if not float(value).is_integer() else str(int(value))
        else:
            sval = str(value)
        attrs, n = re.subn(r'\s%s="[^"]*"' % re.escape(attr), ' %s="%s"' % (attr, _xml_escape(sval)), attrs, count=1)
        if not n:
            attrs = attrs + ' %s="%s"' % (attr, _xml_escape(sval))
    out = cdxml[:m.start()] + "<CDXML" + attrs + ">" + cdxml[m.end():]
    if new_fonts:
        entries = "".join('<font id="%s" charset="iso-8859-1" name="%s"/>' % (i, _xml_escape(n)) for i, n in new_fonts)
        if "</fonttable>" in out:
            out = out.replace("</fonttable>", entries + "</fonttable>", 1)
        else:
            m2 = re.search(r"<CDXML\b[^>]*>", out)
            assert m2 is not None
            out = out[:m2.end()] + "\n<fonttable>" + entries + "</fonttable>" + out[m2.end():]
    return out


# =====================================================================================
# Process / executable discovery
# =====================================================================================

_PID_CACHE: Tuple[float, List[int]] = (0.0, [])


def _chemdraw_pids(max_age: float = 2.0) -> List[int]:
    """PIDs of running ChemDraw processes (psutil if present, else ``tasklist``); [] off Windows.

    Cached for ``max_age`` seconds: the dialog watcher asks every second during long COM calls.
    """
    global _PID_CACHE
    if not _is_windows():
        return []
    if max_age > 0 and time.monotonic() - _PID_CACHE[0] < max_age:
        return list(_PID_CACHE[1])
    pids = _scan_chemdraw_pids()
    _PID_CACHE = (time.monotonic(), pids)
    return list(pids)


def _scan_chemdraw_pids() -> List[int]:
    pids: List[int] = []
    try:
        import psutil  # type: ignore

        for p in psutil.process_iter(["pid", "name"]):
            try:
                if (p.info.get("name") or "").lower() in PROCESS_NAMES:
                    pids.append(int(p.info["pid"]))
            except Exception:
                continue
        return pids
    except ImportError:
        pass
    except Exception:
        pids = []
    for image in PROCESS_NAMES:
        try:
            out = subprocess.run(["tasklist", "/FI", "IMAGENAME eq %s" % image, "/FO", "CSV", "/NH"],
                                 capture_output=True, text=True, timeout=20, creationflags=_CREATE_NO_WINDOW)
        except Exception:
            continue
        for row in csv.reader(io.StringIO(out.stdout or "")):
            if len(row) >= 2 and row[0].lower() == image:
                try:
                    pids.append(int(row[1]))
                except ValueError:
                    continue
    return sorted(set(pids))


def _reg_read(hive: Any, path: str, value: str) -> Optional[str]:
    try:
        import winreg  # type: ignore
    except ImportError:
        return None
    for view in (getattr(winreg, "KEY_WOW64_64KEY", 0), getattr(winreg, "KEY_WOW64_32KEY", 0)):
        try:
            with winreg.OpenKey(hive, path, 0, winreg.KEY_READ | view) as key:
                data, _ = winreg.QueryValueEx(key, value)
                if isinstance(data, str) and data.strip():
                    return os.path.expandvars(data.strip())
        except OSError:
            continue
    return None


def _reg_walk_for_exe(hive: Any, path: str, depth: int = 3) -> List[str]:
    try:
        import winreg  # type: ignore
    except ImportError:
        return []
    found: List[str] = []
    for view in (getattr(winreg, "KEY_WOW64_64KEY", 0), getattr(winreg, "KEY_WOW64_32KEY", 0)):
        try:
            key = winreg.OpenKey(hive, path, 0, winreg.KEY_READ | view)
        except OSError:
            continue
        with key:
            i = 0
            while True:
                try:
                    name, data, _ = winreg.EnumValue(key, i)
                except OSError:
                    break
                i += 1
                if isinstance(data, str) and (name in REGISTRY_PATH_VALUES or "path" in name.lower() or "dir" in name.lower()):
                    data = os.path.expandvars(data.strip().strip('"'))
                    if data.lower().endswith(".exe"):
                        found.append(data)
                    elif data:
                        found.append(os.path.join(data, "ChemDraw.exe"))
                        found.append(os.path.join(data, "ChemDraw", "ChemDraw.exe"))
            if depth > 0:
                j = 0
                while True:
                    try:
                        sub = winreg.EnumKey(key, j)
                    except OSError:
                        break
                    j += 1
                    found.extend(_reg_walk_for_exe(hive, path + "\\" + sub, depth - 1))
    return found


def _registry_executable_candidates() -> List[str]:
    try:
        import winreg  # type: ignore
    except ImportError:
        return []
    out: List[str] = []
    hives = [winreg.HKEY_LOCAL_MACHINE, winreg.HKEY_CURRENT_USER]
    hkcr = winreg.HKEY_CLASSES_ROOT
    for hive in hives:                                             # App Paths
        v = _reg_read(hive, APP_PATHS_KEY, "")
        if v:
            out.append(_exe_from_command(v) or v)
        v = _reg_read(hive, APP_PATHS_KEY, "Path")
        if v:
            out.append(os.path.join(v, "ChemDraw.exe"))
    for prog in PROG_IDS:                                          # the automation server itself
        clsid = _reg_read(hkcr, prog + r"\CLSID", "")
        if clsid:
            for sub in (r"CLSID\%s\LocalServer32" % clsid, r"WOW6432Node\CLSID\%s\LocalServer32" % clsid):
                v = _reg_read(hkcr, sub, "")
                if v:
                    out.append(_exe_from_command(v) or v)
    file_classes = list(FILE_PROG_IDS)                             # .cdx / .cdxml file association
    for ext in FILE_EXTENSIONS:
        v = _reg_read(hkcr, ext, "")
        if v and v not in file_classes:
            file_classes.insert(0, v)
    for cls in file_classes:
        v = _reg_read(hkcr, cls + r"\shell\open\command", "")
        if v:
            out.append(_exe_from_command(v) or v)
        clsid = _reg_read(hkcr, cls + r"\CLSID", "")
        if clsid:
            for sub in (r"CLSID\%s\LocalServer32" % clsid, r"CLSID\%s\InprocServer32" % clsid):
                v = _reg_read(hkcr, sub, "")
                if v:
                    exe = _exe_from_command(v) or v
                    out.append(exe if exe.lower().endswith(".exe") else os.path.join(os.path.dirname(exe), "ChemDraw.exe"))
    v = _reg_read(hkcr, r"Applications\ChemDraw.exe\shell\open\command", "")
    if v:
        out.append(_exe_from_command(v) or v)
    for hive in hives:                                             # vendor keys
        for key in REGISTRY_KEYS:
            out.extend(_reg_walk_for_exe(hive, key))
    return out


def find_chemdraw_windows() -> Optional[str]:
    """Path of ChemDraw.exe: ``CHEMDRAW_EXE`` env var, registry (App Paths, COM server, file class, vendor keys), then globs."""
    env = os.environ.get("CHEMDRAW_EXE")
    if env and os.path.isfile(os.path.expandvars(os.path.expanduser(env))):
        return os.path.abspath(os.path.expandvars(os.path.expanduser(env)))
    if not _is_windows():
        return None
    for cand in _registry_executable_candidates():
        if cand and os.path.isfile(cand):
            return os.path.abspath(cand)
    patterns = [os.path.join(base, suffix) for base in _program_dirs() for suffix in EXE_GLOB_SUFFIXES]
    for match in _glob_many(patterns):
        if os.path.isfile(match):
            return os.path.abspath(match)
    return None


def find_stationery(name: str) -> Optional[str]:
    """Locate a ChemDraw stationery/style file (``*.cds``) by path, file name, style key or style name.

    A style key or alias (``acs``, ``jacs``, ``nature``, ``rsc``, ``wiley``, ``default``, ...) is mapped
    through :func:`chemdrawmcp.styles.get_style` to the title of the stationery ChemDraw ships; then the
    ``*.cds`` stems in the ``ChemDraw Items`` folders are matched exactly (ignoring case and spaces) and
    finally by prefix, like the macOS backend.
    """
    if not name:
        return None
    direct = os.path.expandvars(os.path.expanduser(name))
    if os.path.isfile(direct):
        return os.path.abspath(direct)
    wants: List[str] = []
    try:
        from ..styles import get_style

        title = get_style(name).get("stationery")
        if title:
            wants.append(os.path.splitext(title)[0].lower())
    except KeyError:
        pass
    want = os.path.splitext(os.path.basename(name))[0].lower().strip()
    if want and want not in wants:
        wants.append(want)
    bases = _program_dirs() + [os.environ.get("APPDATA") or "", os.path.expanduser("~")]
    dirs = _glob_many(os.path.join(b, s) for b in bases if b for s in STATIONERY_GLOB_SUFFIXES)
    stems: List[Tuple[str, str]] = [(os.path.splitext(os.path.basename(f))[0].lower(), f)
                                    for d in dirs for f in _glob_many([os.path.join(d, "*.cds")])]
    for w in wants:
        for stem, f in stems:
            if stem == w or stem.replace(" ", "") == w.replace(" ", ""):
                return f
    for w in wants:
        for stem, f in stems:
            if stem.startswith(w):
                return f
    return None


# =====================================================================================
# COM connection
# =====================================================================================

def _com_library() -> Tuple[str, Any]:
    """(name, module) of the available COM client library: pywin32's win32com.client or comtypes.client."""
    errors = []
    try:
        import win32com.client  # type: ignore

        return "win32com", win32com.client
    except ImportError as exc:
        errors.append("win32com: %s" % exc)
    try:
        import comtypes.client  # type: ignore

        return "comtypes", comtypes.client
    except ImportError as exc:
        errors.append("comtypes: %s" % exc)
    raise ChemDrawNotAvailable(
        "ChemDraw COM automation needs pywin32 (pip install pywin32) or comtypes (pip install comtypes) "
        "on Windows; neither could be imported (%s)" % "; ".join(errors))


def _com_initialize() -> Optional[str]:
    """``CoInitialize`` (single-threaded apartment) the calling thread; returns the library name or None."""
    try:
        libname, _ = _com_library()
    except ChemDrawNotAvailable:
        return None
    try:
        if libname == "win32com":
            import pythoncom  # type: ignore

            pythoncom.CoInitialize()
        else:
            import comtypes  # type: ignore

            comtypes.CoInitialize()
    except Exception:
        pass
    return libname


def _com_uninitialize(libname: Optional[str]) -> None:
    try:
        if libname == "win32com":
            import pythoncom  # type: ignore

            pythoncom.CoUninitialize()
        elif libname == "comtypes":
            import comtypes  # type: ignore

            comtypes.CoUninitialize()
    except Exception:
        pass


def _com_connect(prog_ids: Sequence[str] = PROG_IDS, create: bool = True) -> Tuple[Any, str, str]:
    """Connect to ChemDraw's COM server.  Returns ``(app, prog_id, "active"|"created"|"created-while-running")``.

    Tries ``GetActiveObject`` (running instance, needs the Running Object Table) for every ProgID,
    then (if ``create``) ``Dispatch``/``CreateObject`` which starts ChemDraw.  Must be called on the
    thread that will use the object (the backend's COM worker).  Tests monkeypatch this function.
    """
    libname, lib = _com_library()
    _com_initialize()
    errors: List[str] = []
    for pid in prog_ids:
        try:
            return lib.GetActiveObject(pid), pid, "active"
        except Exception as exc:
            errors.append("GetActiveObject(%s): %s" % (pid, _exc_text(exc)))
    if create:
        running = bool(_chemdraw_pids())
        for pid in prog_ids:
            try:
                app = lib.Dispatch(pid) if libname == "win32com" else lib.CreateObject(pid)
                return app, pid, ("created-while-running" if running else "created")
            except Exception as exc:
                errors.append("CreateObject(%s): %s" % (pid, _exc_text(exc)))
    raise ChemDrawNotAvailable("could not connect to ChemDraw through COM (%s). Is ChemDraw installed and "
                               "registered (run it once as administrator to re-register the COM server)?  Edit "
                               "PROG_IDS in chemdrawmcp/backends/windows.py if the ProgID differs."
                               % "; ".join(errors[:6]))


class _ComWorker:
    """One dedicated single-threaded-apartment thread that owns the COM connection and runs every call.

    COM proxies must be used from the apartment that created them; an MCP server calls tools from
    arbitrary threads, so the backend funnels everything through :meth:`submit`.  The thread is
    started lazily and is a daemon; :meth:`abandon` detaches a thread stuck inside ChemDraw so the
    next call starts a fresh one (and reconnects).
    """

    def __init__(self, name: str = "chemdrawmcp-com") -> None:
        self.name = name
        self._lock = threading.Lock()
        self._thread: Optional[threading.Thread] = None
        self._queue: Optional["queue.Queue[Any]"] = None
        self.libname: Optional[str] = None
        self.generation = 0

    def _ensure(self) -> "queue.Queue[Any]":
        with self._lock:
            if self._thread is None or not self._thread.is_alive() or self._queue is None:
                self._queue = queue.Queue()
                self.generation += 1
                self._thread = threading.Thread(target=self._loop, args=(self._queue,), name=self.name, daemon=True)
                self._thread.start()
            return self._queue

    def _loop(self, q: "queue.Queue[Any]") -> None:
        self.libname = _com_initialize()
        try:
            while True:
                item = q.get()
                if item is None:
                    return
                fn, fut = item
                if not fut.set_running_or_notify_cancel():
                    continue
                try:
                    fut.set_result(fn())
                except BaseException as exc:       # noqa: BLE001 - everything goes back to the caller
                    fut.set_exception(exc)
        finally:
            _com_uninitialize(self.libname)

    def on_worker(self) -> bool:
        t = self._thread
        return t is not None and threading.current_thread() is t

    def submit(self, fn: Callable[[], Any]) -> "Future[Any]":
        q = self._ensure()
        fut: "Future[Any]" = Future()
        q.put((fn, fut))
        return fut

    def abandon(self) -> None:
        with self._lock:
            q, self._thread, self._queue = self._queue, None, None
        if q is not None:
            q.put(None)            # the stuck thread exits once its blocked call returns

    close = abandon


# =====================================================================================
# Win32 window helpers (ctypes, no pywin32 needed)
# =====================================================================================

_USER32: Any = None
_KERNEL32: Any = None


def _user32() -> Any:
    """``ctypes.windll.user32`` with argtypes set for the functions we use; None off Windows."""
    global _USER32
    if _USER32 is not None:
        return _USER32 or None
    if not _is_windows() or not hasattr(ctypes, "windll"):
        _USER32 = False
        return None
    from ctypes import wintypes

    u = ctypes.windll.user32  # type: ignore[attr-defined]
    H = wintypes.HWND
    u.GetWindowTextLengthW.argtypes = [H]
    u.GetWindowTextW.argtypes = [H, wintypes.LPWSTR, ctypes.c_int]
    u.GetClassNameW.argtypes = [H, wintypes.LPWSTR, ctypes.c_int]
    u.GetWindowThreadProcessId.argtypes = [H, ctypes.POINTER(wintypes.DWORD)]
    u.GetWindowThreadProcessId.restype = wintypes.DWORD
    u.IsWindowVisible.argtypes = [H]
    u.IsWindow.argtypes = [H]
    u.IsIconic.argtypes = [H]
    u.GetWindowRect.argtypes = [H, ctypes.POINTER(wintypes.RECT)]
    u.GetWindow.argtypes = [H, wintypes.UINT]
    u.GetWindow.restype = H
    u.GetDlgCtrlID.argtypes = [H]
    u.PostMessageW.argtypes = [H, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM]
    u.SendMessageTimeoutW.argtypes = [H, wintypes.UINT, wintypes.WPARAM, wintypes.LPARAM, wintypes.UINT,
                                      wintypes.UINT, ctypes.POINTER(ctypes.c_size_t)]
    u.SendMessageTimeoutW.restype = ctypes.c_ssize_t
    u.SetForegroundWindow.argtypes = [H]
    u.ShowWindow.argtypes = [H, ctypes.c_int]
    u.OpenClipboard.argtypes = [H]
    u.CloseClipboard.argtypes = []
    u.EmptyClipboard.argtypes = []
    u.SetClipboardData.argtypes = [wintypes.UINT, wintypes.HANDLE]
    u.SetClipboardData.restype = wintypes.HANDLE
    u.GetClipboardData.argtypes = [wintypes.UINT]
    u.GetClipboardData.restype = wintypes.HANDLE
    u.IsClipboardFormatAvailable.argtypes = [wintypes.UINT]
    _USER32 = u
    return u


def _kernel32() -> Any:
    global _KERNEL32
    if _KERNEL32 is not None:
        return _KERNEL32 or None
    if not _is_windows() or not hasattr(ctypes, "windll"):
        _KERNEL32 = False
        return None
    from ctypes import wintypes

    k = ctypes.windll.kernel32  # type: ignore[attr-defined]
    k.GlobalAlloc.argtypes = [wintypes.UINT, ctypes.c_size_t]
    k.GlobalAlloc.restype = wintypes.HGLOBAL
    k.GlobalLock.argtypes = [wintypes.HGLOBAL]
    k.GlobalLock.restype = ctypes.c_void_p
    k.GlobalUnlock.argtypes = [wintypes.HGLOBAL]
    k.GlobalFree.argtypes = [wintypes.HGLOBAL]
    k.GlobalSize.argtypes = [wintypes.HGLOBAL]
    k.GlobalSize.restype = ctypes.c_size_t
    _KERNEL32 = k
    return k


def _enum_windows(parent: Optional[int] = None) -> List[int]:
    u = _user32()
    if u is None:
        return []
    from ctypes import wintypes

    proc_t = ctypes.WINFUNCTYPE(wintypes.BOOL, wintypes.HWND, wintypes.LPARAM)  # type: ignore[attr-defined]
    out: List[int] = []

    def _cb(hwnd: Any, _lp: Any) -> bool:
        out.append(int(hwnd or 0))
        return True

    cb = proc_t(_cb)
    if parent is None:
        u.EnumWindows.argtypes = [proc_t, wintypes.LPARAM]
        u.EnumWindows(cb, 0)
    else:
        u.EnumChildWindows.argtypes = [wintypes.HWND, proc_t, wintypes.LPARAM]
        u.EnumChildWindows(parent, cb, 0)
    return out


def _window_text(hwnd: int) -> str:
    u = _user32()
    if u is None:
        return ""
    n = int(u.GetWindowTextLengthW(hwnd) or 0)
    buf = ctypes.create_unicode_buffer(n + 2)
    u.GetWindowTextW(hwnd, buf, n + 2)
    text = buf.value
    if not text:   # controls of other processes: ask through WM_GETTEXT
        res = ctypes.c_size_t(0)
        if u.SendMessageTimeoutW(hwnd, _WM_GETTEXTLENGTH, 0, 0, _SMTO_ABORTIFHUNG, 500, ctypes.byref(res)) and res.value:
            buf = ctypes.create_unicode_buffer(int(res.value) + 2)
            u.SendMessageTimeoutW(hwnd, _WM_GETTEXT, len(buf), ctypes.addressof(buf),
                                  _SMTO_ABORTIFHUNG, 500, ctypes.byref(res))
            text = buf.value
    return text


def _window_class(hwnd: int) -> str:
    u = _user32()
    if u is None:
        return ""
    buf = ctypes.create_unicode_buffer(256)
    u.GetClassNameW(hwnd, buf, 256)
    return buf.value


def _window_pid(hwnd: int) -> int:
    u = _user32()
    if u is None:
        return 0
    from ctypes import wintypes

    pid = wintypes.DWORD(0)
    u.GetWindowThreadProcessId(hwnd, ctypes.byref(pid))
    return int(pid.value)


def _window_rect(hwnd: int) -> Tuple[int, int, int, int]:
    u = _user32()
    if u is None:
        return (0, 0, 0, 0)
    from ctypes import wintypes

    r = wintypes.RECT()
    u.GetWindowRect(hwnd, ctypes.byref(r))
    return int(r.left), int(r.top), int(r.right), int(r.bottom)


def _window_alive(hwnd: int) -> bool:
    u = _user32()
    return bool(u and u.IsWindow(hwnd) and u.IsWindowVisible(hwnd))


def _post(hwnd: int, msg: int, wparam: int = 0, lparam: int = 0) -> None:
    u = _user32()
    if u is not None:
        u.PostMessageW(hwnd, msg, wparam, lparam)


def _send(hwnd: int, msg: int, wparam: int = 0, lparam: int = 0, timeout_ms: int = 1000) -> Optional[int]:
    u = _user32()
    if u is None:
        return None
    res = ctypes.c_size_t(0)
    ok = u.SendMessageTimeoutW(hwnd, msg, wparam, lparam, _SMTO_ABORTIFHUNG, timeout_ms, ctypes.byref(res))
    return int(res.value) if ok else None


def _chemdraw_top_windows(pids: Sequence[int]) -> List[Tuple[int, str, str]]:
    """Visible top-level windows ``(hwnd, class, title)`` of the ChemDraw process(es)."""
    out = []
    for hwnd in _enum_windows():
        if not _window_alive(hwnd):
            continue
        title = _window_text(hwnd)
        if pids:
            if _window_pid(hwnd) not in pids:
                continue
        elif "chemdraw" not in title.lower():
            u = _user32()
            owner = int(u.GetWindow(hwnd, _GW_OWNER) or 0) if u else 0
            if not owner or "chemdraw" not in _window_text(owner).lower():
                continue
        out.append((hwnd, _window_class(hwnd), title))
    return out


def _is_dialog_class(cls: str) -> bool:
    return cls in DIALOG_CLASSES or "dialog" in (cls or "").lower()


def _dialog_texts(hwnd: int) -> List[str]:
    texts = []
    title = _window_text(hwnd)
    if title:
        texts.append(title)
    for child in _enum_windows(hwnd):
        cls = _window_class(child).lower()
        if cls in ("static", "button", "edit", "richedit20w", "richedit50w", "syslink") or "static" in cls or "edit" in cls:
            t = _window_text(child).strip()
            if t:
                texts.append(t.replace("&", ""))
    return texts


def _dialog_buttons(hwnd: int) -> List[Tuple[int, int, str]]:
    """``(button hwnd, control id, label)`` of the push buttons in a dialog."""
    u = _user32()
    out = []
    for child in _enum_windows(hwnd):
        if _window_class(child).lower() == "button":
            label = _window_text(child).replace("&", "").strip()
            out.append((child, int(u.GetDlgCtrlID(child) or 0), label))
    return out


def _find_dialogs(pids: Sequence[int]) -> List[Dict[str, Any]]:
    """Modal dialogs of the ChemDraw process(es): ``[{hwnd, class, title, text: [...], buttons: [(hwnd, id, label)]}]``.

    Tests replace this function with a fake enumerator; off Windows it returns [].
    """
    if _user32() is None:
        return []
    out = []
    for hwnd, cls, title in _chemdraw_top_windows(pids):
        if not _is_dialog_class(cls):
            continue
        out.append({"hwnd": hwnd, "class": cls, "title": title, "text": _dialog_texts(hwnd), "buttons": _dialog_buttons(hwnd)})
    return out


def _press_dialog(hwnd: int, prefer: Optional[Sequence[str]] = None,
                  buttons: Optional[Sequence[Tuple[int, int, str]]] = None) -> str:
    """Try to dismiss one dialog (preferred button, default button, IDOK/IDCANCEL, WM_CLOSE); returns what worked."""
    def _gone() -> bool:
        time.sleep(0.35)
        return not _window_alive(hwnd)

    if buttons is None:
        buttons = _dialog_buttons(hwnd)
    wanted = [p.lower().replace("&", "").replace("'", "").replace("\u2019", "").strip() for p in (prefer or DEFAULT_DIALOG_BUTTONS)]
    for w in wanted:
        for bh, cid, label in buttons:
            norm = label.lower().replace("'", "").replace("\u2019", "").strip()
            if norm == w or (w and w in norm):
                _post(hwnd, _WM_COMMAND, (0 << 16) | (cid & 0xFFFF), bh)   # BN_CLICKED
                if _gone():
                    return "button %r" % label
                _post(bh, _BM_CLICK, 0, 0)
                if _gone():
                    return "click %r" % label
    defid = _send(hwnd, _DM_GETDEFID)
    if defid and (defid >> 16) & 0xFFFF == _DC_HASDEFID:
        ctrl = defid & 0xFFFF
        _post(hwnd, _WM_COMMAND, ctrl, 0)
        if _gone():
            return "default button (id %d)" % ctrl
    for ctrl, what in ((_IDOK, "IDOK"), (_IDCANCEL, "IDCANCEL")):
        _post(hwnd, _WM_COMMAND, ctrl, 0)
        if _gone():
            return what
    _post(hwnd, _WM_CLOSE)
    if _gone():
        return "WM_CLOSE"
    return "still open"


# ------------------------------------------------------------------ clipboard (ctypes, CF_UNICODETEXT)

def _clipboard_open(u: Any) -> None:
    for _ in range(20):
        if u.OpenClipboard(None):
            return
        time.sleep(0.05)
    raise ChemDrawError("could not open the Windows clipboard (another application holds it)")


def _clipboard_get_text() -> Optional[str]:
    """Current clipboard text (CF_UNICODETEXT) or None when it holds no text / off Windows."""
    u, k = _user32(), _kernel32()
    if u is None or k is None:
        return None
    try:
        _clipboard_open(u)
    except ChemDrawError:
        return None
    try:
        if not u.IsClipboardFormatAvailable(_CF_UNICODETEXT):
            return None
        h = u.GetClipboardData(_CF_UNICODETEXT)
        if not h:
            return None
        p = k.GlobalLock(h)
        try:
            return ctypes.wstring_at(p) if p else None
        finally:
            k.GlobalUnlock(h)
    finally:
        u.CloseClipboard()


def _clipboard_set_text(text: str) -> None:
    """Put ``text`` on the clipboard as CF_UNICODETEXT (UTF-16LE, NUL terminated)."""
    u, k = _user32(), _kernel32()
    if u is None or k is None:
        raise ChemDrawError("the clipboard paste fallback needs Windows (user32/kernel32)")
    data = (text or "").encode("utf-16-le") + b"\x00\x00"
    _clipboard_open(u)
    try:
        u.EmptyClipboard()
        h = k.GlobalAlloc(_GMEM_MOVEABLE, len(data))
        if not h:
            raise ChemDrawError("GlobalAlloc failed for %d clipboard bytes" % len(data))
        p = k.GlobalLock(h)
        if not p:
            k.GlobalFree(h)
            raise ChemDrawError("GlobalLock failed")
        ctypes.memmove(p, data, len(data))
        k.GlobalUnlock(h)
        if not u.SetClipboardData(_CF_UNICODETEXT, h):
            k.GlobalFree(h)
            raise ChemDrawError("SetClipboardData(CF_UNICODETEXT) failed")
    finally:
        u.CloseClipboard()


# =====================================================================================
# ChemScript bridge (optional)
# =====================================================================================

class ChemScriptBridge:
    """Lazy wrapper around the ChemScript Python SDK that ChemDraw installs (optional).

    ChemScript is only importable from the Python interpreter ChemDraw installed it for (usually
    ``C:\\Python39``); set ``CHEMSCRIPT_PATH`` to the directory holding the ``ChemScriptNN``
    package and/or ``CHEMSCRIPT_MODULE`` to its import name when auto-detection fails.
    """

    def __init__(self, module_names: Optional[Sequence[str]] = None, extra_paths: Optional[Sequence[str]] = None):
        self.module_names = list(module_names or CHEMSCRIPT_MODULES)
        self.extra_paths = list(extra_paths or [])
        self.module_name: Optional[str] = None
        self._mod: Any = None
        self._tried = False
        self.error: Optional[str] = None

    # -- import -----------------------------------------------------------------------
    def _search_paths(self) -> List[str]:
        paths = [os.environ.get("CHEMSCRIPT_PATH") or ""] + self.extra_paths
        if _is_windows():
            paths += _glob_many(os.path.join(b, s) for b in _program_dirs() for s in CHEMSCRIPT_GLOB_SUFFIXES)
            paths += CHEMSCRIPT_PYTHON_DIRS
        return [p for p in paths if p and os.path.isdir(p)]

    def module(self) -> Any:
        if self._tried:
            return self._mod
        self._tried = True
        for p in self._search_paths():
            if p not in sys.path:
                sys.path.append(p)
        names = list(self.module_names)
        env = os.environ.get("CHEMSCRIPT_MODULE")
        if env:
            names.insert(0, env)
        errors = []
        for name in names:
            try:
                self._mod = importlib.import_module(name)
                self.module_name = name
                return self._mod
            except Exception as exc:   # ImportError, or a DLL load failure on the wrong interpreter
                errors.append("%s: %s" % (name, exc))
        self.error = "; ".join(errors[:4])
        return None

    def available(self) -> bool:
        return self.module() is not None

    def _require(self) -> Any:
        mod = self.module()
        if mod is None:
            raise ChemDrawNotAvailable(
                "ChemScript (ChemDraw's Python SDK) is not importable in this interpreter (%s). Run chemdrawmcp with "
                "the Python ChemDraw installed it for, or set CHEMSCRIPT_PATH / CHEMSCRIPT_MODULE (CHEMSCRIPT_MODULES "
                "in chemdrawmcp/backends/windows.py lists the names tried)." % (self.error or "not tried"))
        return mod

    # -- member probing -----------------------------------------------------------------
    @staticmethod
    def _member(obj: Any, key: str, *args: Any) -> Any:
        errors = []
        for name in CHEMSCRIPT_MEMBERS[key]:
            try:
                val = getattr(obj, name)
            except AttributeError:
                errors.append("%s: not present" % name)
                continue
            if callable(val):
                try:
                    return val(*args)
                except Exception as exc:
                    errors.append("%s: %s" % (name, _exc_text(exc)))
                    continue
            if args:
                errors.append("%s: is a property" % name)
                continue
            return val
        raise ChemDrawError("ChemScript: none of %s worked for %r (%s); edit CHEMSCRIPT_MEMBERS[%r] in "
                            "chemdrawmcp/backends/windows.py" % (CHEMSCRIPT_MEMBERS[key], key, "; ".join(errors), key))

    # -- API ----------------------------------------------------------------------------
    def load(self, data: str, mime: Optional[str] = None) -> Any:
        """``StructureData.LoadData(data[, mime])``; raises when ChemScript cannot parse the input."""
        mod = self._require()
        sd_cls = getattr(mod, "StructureData", None)
        if sd_cls is None:
            raise ChemDrawError("ChemScript module %s has no StructureData class" % self.module_name)
        sd = None
        errors: List[str] = []
        for m in (_mime_candidates(mime) if mime else [None]):
            try:
                sd = self._member(sd_cls, "load_data", data, m) if m else self._member(sd_cls, "load_data", data)
            except ChemDrawError as exc:
                errors.append(str(exc))
                continue
            if sd is not None:
                break
        if sd is None:
            raise ChemDrawError("ChemScript could not interpret the input as %s: %.80r%s"
                                % (mime or "chemistry", data, (" (%s)" % "; ".join(errors)) if errors else ""))
        return sd

    def convert(self, data: str, in_mime: str, out_mime: str, clean: bool = False) -> str:
        """Convert chemistry text between MIME types (e.g. chemical/x-name -> chemical/x-cdxml)."""
        if out_mime == MIME["name"]:
            return self.structure_to_name(data, in_mime)
        sd = self.load(data, in_mime)
        if clean:
            try:
                self._member(sd, "clean")
            except ChemDrawError:
                pass
        if out_mime == MIME["smiles"]:
            try:
                return _as_text(self._member(sd, "smiles"))
            except ChemDrawError:
                pass
        errors: List[str] = []
        for m in _mime_candidates(out_mime):
            try:
                return _as_text(self._member(sd, "write_data", m))
            except ChemDrawError as exc:
                errors.append(str(exc))
        raise ChemDrawError("ChemScript could not write %s (%s)" % (out_mime, "; ".join(errors)))

    def name_to_structure(self, name: str, out: str = "smiles") -> str:
        """Chemical name -> structure text (``out`` is a key of ``MIME`` or a MIME type)."""
        out_mime = MIME.get(out, out)
        return self.convert(name, MIME["name"], out_mime, clean=(out_mime == MIME["cdxml"]))

    def structure_to_name(self, data: str, in_mime: str = MIME["smiles"]) -> str:
        """Structure text -> IUPAC name (Struct=Name)."""
        sd = self.load(data, in_mime)
        try:
            return _as_text(self._member(sd, "name"))
        except ChemDrawError:
            return _as_text(self._member(sd, "write_data", MIME["name"]))

    def properties(self, data: str, in_mime: Optional[str] = None) -> Dict[str, Any]:
        sd = self.load(data, in_mime)
        out: Dict[str, Any] = {}
        for key in ("smiles", "inchi", "formula", "molecular_weight", "exact_mass"):
            try:
                val = self._member(sd, key)
                out[key] = _num(val) if key in ("molecular_weight", "exact_mass") else _as_text(val)
            except ChemDrawError:
                continue
        return out


_DEFAULT_BRIDGE = ChemScriptBridge()


def chemscript_available() -> bool:
    """True when ChemDraw's ChemScript Python SDK can be imported in this interpreter."""
    try:
        return _DEFAULT_BRIDGE.available()
    except Exception:
        return False


# =====================================================================================
# The backend
# =====================================================================================

def _com_method(extra_budget: Optional[Callable[..., float]] = None) -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Run the decorated method body on the backend's COM worker thread (see ``_ComWorker``).

    ``extra_budget(self, *args, **kw)`` may return extra seconds the method is allowed to take beyond
    ``call_timeout`` (``launch`` sleeps while ChemDraw starts).
    """
    def deco(method: Callable[..., Any]) -> Callable[..., Any]:
        @functools.wraps(method)
        def wrapper(self: "WindowsBackend", *args: Any, **kw: Any) -> Any:
            extra = float(extra_budget(self, *args, **kw)) if extra_budget else 0.0
            return self._run(lambda: method(self, *args, **kw), what=method.__name__, budget=self.call_timeout + extra)
        return wrapper
    return deco


class WindowsBackend(ChemDrawBackend):
    """ChemDraw automation on Windows through COM (``ChemDraw.Application``).

    ``app_factory`` (tests) returns the application object instead of connecting through COM;
    ``call_timeout`` is how long one tool call may block inside ChemDraw before the call is given up
    (modal dialogs found meanwhile are dismissed and the budget restarts); ``com_thread=False``
    disables the dedicated COM thread.
    """

    platform = "windows"

    def __init__(self, app_factory: Optional[Callable[[], Any]] = None, prog_ids: Optional[Sequence[str]] = None,
                 call_timeout: float = 120.0, executable_path: Optional[str] = None,
                 chemscript: Optional[ChemScriptBridge] = None, com_thread: bool = True, dialog_poll: float = 1.0):
        self.prog_ids = list(prog_ids or PROG_IDS)
        self._app_factory = app_factory
        self.call_timeout = float(call_timeout)
        self.dialog_poll = float(dialog_poll)
        self._exe = executable_path
        self._local = threading.local()
        self._worker: Optional[_ComWorker] = _ComWorker() if com_thread else None
        self.prog_id: Optional[str] = None
        self.connect_note: Optional[str] = None
        self.last_dialogs: List[str] = []
        self.chemscript = chemscript or _DEFAULT_BRIDGE
        self._index_base: Optional[int] = None
        self._prefer: Optional[Sequence[str]] = None
        #: Learned at run time: does Objects.SetData replace the document content?  None = not known yet.
        self.setdata_replaces: Optional[bool] = {"add": False, "replace": True}.get(INSERT_MODE)

    # ------------------------------------------------------------------ threading
    def close(self) -> None:
        """Stop the COM worker thread (the connection is dropped; the next call starts a new one)."""
        if self._worker is not None:
            self._worker.close()
        self._forget_app()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _run(self, fn: Callable[[], Any], what: str, budget: Optional[float] = None) -> Any:
        """Execute ``fn`` on the COM worker and wait, dismissing ChemDraw's modal dialogs meanwhile."""
        worker = self._worker
        if worker is None or worker.on_worker():
            return fn()
        budget = max(0.1, float(budget or self.call_timeout))
        fut = worker.submit(fn)
        deadline = time.monotonic() + budget
        dismissed: List[str] = []
        while True:
            remaining = deadline - time.monotonic()
            try:
                result = fut.result(timeout=max(0.0, min(self.dialog_poll, remaining)))
            except _FutureTimeout:
                pass
            else:
                return self._attach_alerts(result, dismissed)
            texts = self._safe_dismiss(self._prefer)
            if texts:
                dismissed.extend(texts)
                self.last_dialogs.extend(texts)
                deadline = time.monotonic() + budget       # ChemDraw was blocked by a dialog: fresh budget
                continue
            if time.monotonic() >= deadline:
                worker.abandon()
                msg = "ChemDraw did not answer %s within %.0fs" % (what, budget)
                if dismissed:
                    msg += "; dialogs dismissed meanwhile: " + " || ".join(dismissed)
                msg += (". The COM call is still blocked inside ChemDraw (look for a modal dialog or a hung "
                        "ChemDraw.exe; call dismiss_dialogs); the next call reconnects on a fresh COM thread.")
                raise ChemDrawError(msg)

    @staticmethod
    def _attach_alerts(result: Any, dismissed: List[str]) -> Any:
        """Dialogs dismissed by the waiting thread during a call belong in that call's ``alerts`` (the
        worker may have read ``last_dialogs`` before the waiter appended them)."""
        if dismissed and isinstance(result, dict) and isinstance(result.get("alerts"), list):
            for text in dismissed:
                if text not in result["alerts"]:
                    result["alerts"].append(text)
        return result

    # ------------------------------------------------------------------ COM plumbing
    def _app(self, reconnect: bool = False) -> Any:
        app = getattr(self._local, "app", None)
        if app is not None and not reconnect:
            return app
        if self._app_factory is not None:
            app = self._app_factory()
        else:
            app, self.prog_id, how = _com_connect(self.prog_ids)
            self.connect_note = None
            if how == "created-while-running":
                self.connect_note = ("ChemDraw.exe was already running but GetActiveObject found no instance in the "
                                     "Running Object Table, so CreateObject was used - this may have started a second "
                                     "ChemDraw. If so, let the MCP server launch ChemDraw itself next time.")
        self._local.app = app
        return app

    def _app_or_none(self) -> Any:
        try:
            return self._app()
        except ChemDrawError:
            return None

    def _forget_app(self) -> None:
        self._local.app = None

    def _safe_dismiss(self, prefer: Optional[Sequence[str]] = None) -> List[str]:
        try:
            return self._dismiss(prefer)
        except Exception:
            return []

    def _call(self, fn: Callable[..., Any], *args: Any, what: Optional[str] = None,
              prefer_buttons: Optional[Sequence[str]] = None, **kw: Any) -> Any:
        """Invoke one COM member, converting COM errors into ChemDrawError (with any dialog text appended).

        ``prefer_buttons`` tells the dialog watcher (``_run`` on the calling thread) which button to
        press should ChemDraw show a dialog during this call (e.g. "Don't Save" while closing).
        """
        label = what or getattr(fn, "__name__", None) or getattr(fn, "_name", None) or repr(fn)
        seen_before = len(self.last_dialogs)
        if prefer_buttons:
            self._prefer = prefer_buttons
        try:
            return fn(*args, **kw)
        except ChemDrawError:
            raise
        except Exception as exc:
            if not (_is_com_error(exc) or isinstance(exc, (TypeError, ValueError, OSError))):
                raise
            desc = _exc_text(exc)
            if _hresult(exc) in _DISCONNECTED_HRESULTS:
                self._forget_app()
                desc += " - ChemDraw seems to have quit; reconnecting on the next call"
            dialogs = list(self.last_dialogs[seen_before:]) + self._safe_dismiss(prefer_buttons or self._prefer)
            if len(dialogs) > len(self.last_dialogs[seen_before:]):
                self.last_dialogs.extend(dialogs[len(self.last_dialogs[seen_before:]):])
            msg = "ChemDraw COM call %s failed: %s" % (label, desc)
            if dialogs:
                msg += "; ChemDraw showed dialog(s): " + " || ".join(dialogs)
            raise ChemDrawError(msg) from exc
        finally:
            if prefer_buttons:
                self._prefer = None

    @staticmethod
    def _tried(key: str, obj: Any, errors: List[str]) -> str:
        names = MEMBERS.get(key) or []
        return ("ChemDraw COM: none of the candidate members %s for %r worked on %s (%s). If your ChemDraw version "
                "names this member differently, edit MEMBERS[%r] in chemdrawmcp/backends/windows.py."
                % (names, key, type(obj).__name__, "; ".join(errors) or "no candidates", key))

    def _get(self, obj: Any, key: str, default: Any = _MISSING) -> Any:
        """Read a property through the candidate table (a no-argument method is called)."""
        errors: List[str] = []
        if obj is None:
            errors.append("object is None")
        else:
            for name in MEMBERS[key]:
                try:
                    val = getattr(obj, name)
                except Exception as exc:
                    if _member_missing(exc):
                        errors.append("%s: not present" % name)
                        continue
                    # the client library could not read it as a plain property: ask IDispatch directly
                    try:
                        return self._call(_raw_invoke, obj, name, (), what="%s.%s (IDispatch)" % (key, name))
                    except (ChemDrawError, AttributeError) as exc2:
                        errors.append("%s: %s / %s" % (name, _exc_text(exc), exc2))
                        continue
                if _is_method(val):
                    try:
                        val = self._call(val, what="%s.%s" % (key, name))
                    except ChemDrawError as exc:
                        errors.append("%s: %s" % (name, exc))
                        continue
                return val
        if default is not _MISSING:
            return default
        raise ChemDrawError(self._tried(key, obj, errors))

    def _invoke(self, obj: Any, key: str, *args: Any, prefer_buttons: Optional[Sequence[str]] = None) -> Any:
        """Call a method through the candidate table (falling back to a raw IDispatch invoke)."""
        errors: List[str] = []
        if obj is None:
            raise ChemDrawError(self._tried(key, obj, ["object is None"]))
        for name in MEMBERS[key]:
            fn: Any = None
            note = ""
            try:
                fn = getattr(obj, name)
            except Exception as exc:
                if _member_missing(exc):
                    errors.append("%s: not present" % name)
                    continue
                note = _exc_text(exc)          # e.g. comtypes DISP_E_BADPARAMCOUNT for Data(mime): use IDispatch
            if fn is not None and not _is_method(fn):
                if not args:
                    return fn                  # a plain property read
                fn, note = None, "is a property, cannot be called with arguments"
            try:
                if fn is not None:
                    return self._call(fn, *args, what="%s.%s" % (key, name), prefer_buttons=prefer_buttons)
                return self._call(_raw_invoke, obj, name, tuple(args), what="%s.%s (IDispatch)" % (key, name),
                                  prefer_buttons=prefer_buttons)
            except AttributeError as exc:      # not a COM object (fakes) - nothing more to try for this name
                errors.append("%s: %s" % (name, note or exc))
            except ChemDrawError as exc:
                errors.append("%s: %s%s" % (name, ("%s / " % note) if note else "", exc))
        raise ChemDrawError(self._tried(key, obj, errors))

    def _set(self, obj: Any, key: str, value: Any, names: Optional[Sequence[str]] = None) -> str:
        errors: List[str] = []
        for name in names or MEMBERS[key]:
            try:
                setattr(obj, name, value)
                return name
            except Exception as exc:
                if _member_missing(exc):
                    errors.append("%s: not present" % name)
                    continue
                try:
                    self._call(_raw_invoke, obj, name, (value,), "put", what="%s.%s (IDispatch put)" % (key, name))
                    return name
                except (ChemDrawError, AttributeError) as exc2:
                    errors.append("%s: %s / %s" % (name, _exc_text(exc), exc2))
        raise ChemDrawError(self._tried(key, obj, errors) if key in MEMBERS else
                            "ChemDraw COM: could not set %s (%s); edit SETTINGS_MEMBERS in chemdrawmcp/backends/windows.py"
                            % (key, "; ".join(errors)))

    def _count(self, coll: Any, default: int = 0) -> int:
        if coll is None:
            return default
        val = self._get(coll, "coll.count", default=None)
        if val is None:
            try:
                return len(coll)
            except Exception:
                return default
        try:
            return int(val)
        except (TypeError, ValueError):
            return default

    def _item(self, coll: Any, index: int) -> Any:
        """``coll.Item(i)`` with 1-based indexing (0-based is detected and remembered)."""
        bases = [self._index_base] if self._index_base is not None else [1, 0]
        last: Optional[Exception] = None
        for base in bases:
            i = index - 1 + base
            for getter in ("Item", "item", "__call__", "__getitem__"):
                try:
                    fn = getattr(coll, getter)
                    val = fn(i)
                except (AttributeError, IndexError, KeyError, TypeError) as exc:
                    last = exc
                    continue
                except Exception as exc:
                    last = exc
                    if _is_com_error(exc):
                        continue
                    raise
                if val is not None:
                    self._index_base = base
                    return val
        raise ChemDrawError("ChemDraw COM: cannot read item %d of %s (%s); check MEMBERS['coll.item']"
                            % (index, type(coll).__name__, _exc_text(last) if last else "no Item member"))

    def _items(self, coll: Any) -> List[Any]:
        if coll is None:
            return []
        n = self._count(coll, default=-1)
        if n < 0:
            try:
                return list(coll)
            except Exception:
                return []
        return [self._item(coll, i) for i in range(1, n + 1)]

    # ------------------------------------------------------------------ documents plumbing
    def _documents_list(self) -> List[Any]:
        app = self._app()
        docs = self._get(app, "app.documents", default=None)
        if docs is None:
            return []
        out = []
        for d in self._items(docs):
            name = self._doc_name(d)
            if name.startswith("~") or name.lower().startswith("autosave"):
                continue
            out.append(d)
        return out

    def _doc_name(self, doc: Any) -> str:
        name = self._get(doc, "doc.name", default=None)
        if not name:
            full = self._get(doc, "doc.full_name", default=None)
            name = os.path.basename(str(full)) if full else "Untitled"
        return str(name)

    def _doc_path(self, doc: Any) -> Optional[str]:
        """Full path of a saved document; None for untitled ones (``FullName`` == ``Name``, or no folder)."""
        full = self._get(doc, "doc.full_name", default=None)
        if full:
            full = str(full)
            if os.path.isabs(full) or "\\" in full or "/" in full:
                return full
        folder = self._get(doc, "doc.dir", default=None)
        if folder and str(folder).strip():
            name = self._get(doc, "doc.name", default=None)
            if name:
                return os.path.join(str(folder), str(name))
        return None

    def _doc_modified(self, doc: Any) -> bool:
        saved = self._get(doc, "doc.saved", default=None)
        if saved is not None:
            return not bool(saved)
        return bool(self._get(doc, "doc.modified", default=False))

    def _active_doc(self) -> Any:
        app = self._app()
        try:
            doc = self._get(app, "app.active_document", default=None)
        except ChemDrawError:
            doc = None
        if doc is None:
            docs = self._documents_list()
            doc = docs[0] if docs else None
        return doc

    def _doc(self, name: Optional[Any]) -> Any:
        if name is None or not str(name).strip():
            doc = self._active_doc()
            if doc is None:
                raise ChemDrawError("no document is open in ChemDraw")
            return doc
        want = str(name).strip().lower()
        want_stem = os.path.splitext(os.path.basename(want))[0]
        names = []
        docs = self._documents_list()
        for d in docs:
            nm = self._doc_name(d)
            names.append(nm)
            cands = {nm.lower(), os.path.splitext(nm)[0].lower()}
            full = self._get(d, "doc.full_name", default=None)
            if full:
                cands.add(str(full).lower())
                cands.add(os.path.basename(str(full)).lower())
            if want in cands or want_stem in cands:
                return d
        if re.fullmatch(r"#?\d+", want):                 # a window id / index from an earlier call
            idx = int(want.lstrip("#"))
            for d in docs:
                if self._doc_window_id(d) == idx:
                    return d
            if 1 <= idx <= len(docs):
                return docs[idx - 1]
        raise ChemDrawError("no open document named %r; open documents: %s" % (name, ", ".join(names) or "none"))

    def _doc_window_id(self, doc: Any) -> Optional[int]:
        hwnd = self._get(doc, "doc.hwnd", default=None)
        if hwnd is None:
            win = self._get(doc, "doc.window", default=None)
            if win is not None:
                hwnd = self._get(win, "win.hwnd", default=None)
        return int(hwnd) if isinstance(hwnd, int) and hwnd and not isinstance(hwnd, bool) else None

    def _objects(self, doc: Any) -> Any:
        return self._get(doc, "doc.objects")

    # ------------------------------------------------------------------ Data / counts / selection
    def _data_raw(self, target: Any, mime: str, doc: Any = None) -> Any:
        """``target.Data(mime)`` trying the MIME alternatives, and the Document itself as a second holder."""
        errors: List[str] = []
        holders = [target] + ([doc] if doc is not None and doc is not target else [])
        for m in _mime_candidates(mime):
            for holder in holders:
                try:
                    return self._invoke(holder, "objs.data", m)
                except ChemDrawError as exc:
                    errors.append("%s on %s: %s" % (m, type(holder).__name__, exc))
        raise ChemDrawError("ChemDraw could not produce %s (%s). If ChemDraw spells this MIME type differently, edit "
                            "MIME / MIME_ALTERNATIVES in chemdrawmcp/backends/windows.py." % (mime, " / ".join(errors)))

    def _data(self, target: Any, mime: str, doc: Any = None) -> str:
        return _as_text(self._data_raw(target, mime, doc))

    def _cdxml_or_empty(self, target: Any, doc: Any = None) -> str:
        try:
            text = self._data(target, MIME["cdxml"], doc)
        except ChemDrawError:
            return ""
        return text if text.lstrip().startswith("<") else ""

    def _counts(self, objs: Any, cdxml: Optional[str] = None) -> Dict[str, int]:
        """atoms/bonds/molecules/captions/arrows of an Objects/Selection collection.

        Uses the ``Atoms``/``Bonds``/... sub-collections when the COM object has them, otherwise
        counts elements in its CDXML (``cdxml`` may be passed in to avoid a second ``Data`` call).
        """
        out: Dict[str, int] = {}
        missing = []
        for key, member in (("atoms", "objs.atoms"), ("bonds", "objs.bonds"), ("molecules", "objs.fragments"),
                            ("captions", "objs.captions"), ("arrows", "objs.arrows")):
            coll = self._get(objs, member, default=None) if objs is not None else None
            if coll is None:
                missing.append(key)
                out[key] = 0
            else:
                out[key] = self._count(coll)
        if missing and objs is not None:
            if cdxml is None:
                cdxml = self._cdxml_or_empty(objs)
            if cdxml and cdxml.lstrip().startswith("<"):
                summary = cdxml_summary(cdxml)
                for key in missing:
                    out[key] = int(summary.get(key) or 0)
        return out

    def _selection(self, doc: Any) -> Tuple[Any, bool]:
        """(Selection collection or None, something selected?)."""
        sel = self._get(doc, "doc.selection", default=None)
        if sel is None:
            return None, False
        n = self._count(sel, default=-1)
        if n > 0:
            return sel, True
        if n < 0:                               # no Count member: look at the content
            cdxml = self._cdxml_or_empty(sel)
            if cdxml and _total_objects(cdxml_summary(cdxml)) > 0:
                return sel, True
        return sel, False

    def _target(self, doc: Any) -> Tuple[Any, bool]:
        """(Selection if something is selected else Objects, selected?)."""
        sel, selected = self._selection(doc)
        if selected:
            return sel, True
        return self._objects(doc), False

    # ------------------------------------------------------------------ dialogs
    def _dismiss(self, prefer: Optional[Sequence[str]] = None) -> List[str]:
        texts: List[str] = []
        for dlg in _find_dialogs(_chemdraw_pids()):
            how = _press_dialog(dlg["hwnd"], prefer, dlg.get("buttons"))
            desc = " | ".join(dlg.get("text") or []) or dlg.get("title") or "<dialog>"
            texts.append("%s [%s]" % (desc, how))
        return texts

    def dismiss_dialogs(self) -> List[str]:
        """Close any modal dialog ChemDraw shows (runs on the calling thread so it works while a COM call is blocked)."""
        texts = self._dismiss()
        if texts:
            self.last_dialogs.extend(texts)
        return texts

    def pending_dialogs(self) -> List[Dict[str, Any]]:
        """Modal dialogs currently shown, without dismissing them: ``[{title, text, buttons, hwnd}]``."""
        out = []
        for dlg in _find_dialogs(_chemdraw_pids()):
            out.append({"title": dlg.get("title") or "", "text": list(dlg.get("text") or []),
                        "buttons": [b[2] for b in (dlg.get("buttons") or [])], "hwnd": dlg.get("hwnd")})
        return out

    # ------------------------------------------------------------------ lifecycle
    def executable(self) -> Optional[str]:
        if self._exe:
            return self._exe
        try:
            exe = find_chemdraw_windows()
        except Exception:
            exe = None
        if exe:
            self._exe = exe
        return exe

    @_com_method()
    def is_running(self) -> bool:
        if _chemdraw_pids():
            return True
        if self._app_factory is not None:
            try:
                app = self._app()
                return self._get(app, "app.version", default=None) is not None
            except Exception:
                return False
        if not _is_windows():
            return False
        try:
            _com_connect(self.prog_ids, create=False)
            return True
        except ChemDrawError:
            return False

    @_com_method(extra_budget=lambda self, wait_seconds=90.0, *a, **k: float(wait_seconds) + 10.0)
    def launch(self, wait_seconds: float = 90.0) -> dict:
        t0 = time.time()
        info: Dict[str, Any] = {"already_running": self.is_running(), "executable": self.executable()}
        if not info["already_running"] and self._app_factory is None:
            _com_library()          # no pywin32/comtypes: starting ChemDraw cannot help, so fail now (not after the wait)
            try:
                app, self.prog_id, how = _com_connect(self.prog_ids, create=True)
                self._local.app = app
                self.connect_note = None        # nothing was running: this instance is ours
                info["started_via"] = "COM CreateObject(%s)" % self.prog_id
            except ChemDrawNotAvailable as exc:
                exe = self.executable()
                if not exe:
                    raise ChemDrawNotAvailable("ChemDraw.exe was not found (set CHEMDRAW_EXE) and the COM server could "
                                               "not be started: %s" % exc) from exc
                try:
                    proc = subprocess.Popen([exe], close_fds=True)
                except OSError as e:
                    raise ChemDrawError("could not start %s: %s" % (exe, e)) from e
                info["started_via"] = "exe"
                info["pid"] = proc.pid
        deadline = t0 + max(1.0, wait_seconds)
        last_err: Optional[Exception] = None
        version: Optional[str] = None
        while time.time() < deadline:
            try:
                app = self._app()           # the proxy CreateObject returned above; otherwise (re)connect
                version = str(self._get(app, "app.version"))
                break
            except Exception as exc:
                last_err = exc
                self._forget_app()          # dropped: the next round connects afresh
                time.sleep(1.0)
        if version is None:
            raise ChemDrawError("ChemDraw did not answer through COM within %.0fs (%s). If ChemDraw is running but "
                                "GetActiveObject fails, check PROG_IDS in chemdrawmcp/backends/windows.py and that the "
                                "server runs at the same integrity level as ChemDraw." % (wait_seconds, last_err))
        try:
            self._set(self._app(), "app.visible", True)
        except ChemDrawError:
            pass
        pids = _chemdraw_pids()
        info.update({"pid": pids[0] if pids else info.get("pid"), "pids": pids, "version": version,
                     "prog_id": self.prog_id, "seconds": round(time.time() - t0, 1)})
        if self.connect_note:
            info["connect_note"] = self.connect_note
        return info

    @_com_method()
    def quit(self, save: bool = False) -> None:
        app = self._app_or_none()
        if app is None:
            return
        prefer = ("Save", "Yes") if save else ("Don't Save", "Don’t Save", "No", "Discard")
        for d in self._documents_list():
            try:
                if save and self._doc_modified(d) and self._doc_path(d):
                    self._invoke(d, "doc.save")
                if not save:
                    self._invoke(d, "doc.close", False, prefer_buttons=prefer)
            except ChemDrawError:
                continue
        try:
            self._invoke(app, "app.quit", prefer_buttons=prefer)
        finally:
            self._forget_app()
        time.sleep(0.5)
        self._safe_dismiss(prefer)

    @_com_method()
    def version(self) -> str:
        return str(self._get(self._app(), "app.version"))

    @_com_method()
    def activate(self) -> None:
        app = self._app()
        try:
            self._set(app, "app.visible", True)
        except ChemDrawError:
            pass
        try:
            self._invoke(app, "app.activate")
            return
        except ChemDrawError:
            pass
        hwnd = self._main_window(None)
        u = _user32()
        if hwnd and u:
            if u.IsIconic(hwnd):
                u.ShowWindow(hwnd, _SW_RESTORE)
            u.SetForegroundWindow(hwnd)

    # ------------------------------------------------------------------ documents
    @_com_method()
    def documents(self) -> List[DocumentInfo]:
        active = self._active_doc()
        active_key = (self._doc_name(active), self._doc_path(active)) if active is not None else None
        infos: List[DocumentInfo] = []
        for d in self._documents_list():
            name, path = self._doc_name(d), self._doc_path(d)
            objs = self._get(d, "doc.objects", default=None)
            counts = self._counts(objs)
            info = DocumentInfo(name=name, path=path, modified=self._doc_modified(d), front=((name, path) == active_key),
                                atoms=counts["atoms"], bonds=counts["bonds"], molecules=counts["molecules"],
                                captions=counts["captions"], arrows=counts["arrows"], window_id=self._doc_window_id(d))
            infos.append(info)
        if infos and not any(i.front for i in infos):
            infos[0].front = True
        infos.sort(key=lambda i: not i.front)
        return infos

    @_com_method()
    def front_document(self) -> Optional[str]:
        doc = self._active_doc()
        return self._doc_name(doc) if doc is not None else None

    @_com_method()
    def activate_document(self, name: Optional[str]) -> str:
        doc = self._doc(name)
        try:
            self._invoke(doc, "doc.activate")
        except ChemDrawError as exc:
            try:
                self._set(self._app(), "app.active_document", doc)
            except ChemDrawError:
                raise exc
        return self._doc_name(doc)

    @_com_method()
    def new_document(self, stationery: Optional[str] = None) -> str:
        app = self._app()
        docs = self._get(app, "app.documents")
        before = {id(d) for d in self._items(docs)} if stationery else set()
        if stationery:
            path = find_stationery(stationery)
            if not path:
                raise ChemDrawError("stationery %r not found (give the full path of a .cds file or a name from the "
                                    "'ChemDraw Items' folder)" % stationery)
            try:
                doc = self._invoke(docs, "docs.add", path)          # Documents.Add(Template) as in Office
            except ChemDrawError:
                doc = self._invoke(docs, "docs.open", path)         # opening a .cds creates an untitled document
        else:
            doc = self._invoke(docs, "docs.add")
        if doc is None or isinstance(doc, (bool, int, str)):
            doc = self._active_doc()
            if doc is None or (stationery and id(doc) in before):
                fresh = [d for d in self._items(docs) if id(d) not in before]
                doc = fresh[-1] if fresh else doc
        if doc is None:
            raise ChemDrawError("Documents.Add() returned nothing and no document is active; check MEMBERS['docs.add']")
        return self._doc_name(doc)

    @_com_method()
    def open_document(self, path: str) -> str:
        path = os.path.abspath(os.path.expanduser(path))
        if not os.path.isfile(path):
            raise ChemDrawError("file not found: %s" % path)
        docs = self._get(self._app(), "app.documents")
        doc = self._invoke(docs, "docs.open", path)
        if doc is None or isinstance(doc, (bool, int, str)):
            doc = self._active_doc()
        if doc is None:
            raise ChemDrawError("Documents.Open(%r) returned nothing and no document is active; check MEMBERS['docs.open']" % path)
        return self._doc_name(doc)

    @_com_method()
    def close_document(self, name: Optional[str] = None, save: bool = False) -> None:
        doc = self._doc(name)
        prefer = ("Save", "Yes") if save else ("Don't Save", "Don’t Save", "No", "Discard")
        if save:
            if not self._doc_path(doc):
                raise ChemDrawError("document %r has never been saved; use save_document(path) first, then close it"
                                    % self._doc_name(doc))
            self._invoke(doc, "doc.save")
        try:
            self._invoke(doc, "doc.close", bool(save), prefer_buttons=prefer)
        except ChemDrawError:
            self._invoke(doc, "doc.close", prefer_buttons=prefer)
        self._safe_dismiss(prefer)

    @_com_method()
    def save_document(self, name: Optional[str], path: str, fmt: str) -> str:
        fmt = normalize_format(fmt)
        spec = EXPORT_FORMATS[fmt]
        doc = self._doc(name)
        docname = self._doc_name(doc)
        path = _ensure_extension(path, spec["ext"])
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        mime = EXPORT_MIME.get(fmt, spec["mime"])
        errors: List[str] = []

        def _remove() -> None:
            if os.path.exists(path):
                try:
                    os.remove(path)           # avoids ChemDraw's "replace existing file?" prompt
                except OSError as exc:
                    raise ChemDrawError("cannot overwrite %s: %s" % (path, exc)) from exc

        for args in ((path,), (path, mime), (path, spec["name"]), (path, spec["ext"].lstrip("."))):
            _remove()
            try:
                self._invoke(doc, "doc.save_as", *args, prefer_buttons=("Yes", "OK", "Replace"))
            except ChemDrawError as exc:
                errors.append(str(exc))
                continue
            if os.path.isfile(path) and os.path.getsize(path) > 0:
                looks_right = _sniff_format(path, fmt)
                if looks_right is not False:
                    return path
                errors.append("SaveAs%r wrote a file whose content is not %s (ChemDraw ignored the extension?)"
                              % (args, spec["name"]))
                continue
            errors.append("SaveAs%r returned but wrote no file" % (args,))
        _remove()
        objs = self._get(doc, "doc.objects", default=None)
        if objs is not None:
            try:
                raw = _as_bytes(self._data_raw(objs, mime, doc))
                if raw:
                    with open(path, "wb") as fh:
                        fh.write(raw)
                    return path
                errors.append("Objects.Data(%s) returned nothing" % mime)
            except ChemDrawError as exc:
                errors.append(str(exc))
        raise ChemDrawError("could not save %r as %s to %s: %s. If SaveAs needs another format argument, edit "
                            "EXPORT_MIME / MEMBERS['doc.save_as'] in chemdrawmcp/backends/windows.py."
                            % (docname, fmt, path, " / ".join(errors)))

    @_com_method()
    def document_cdxml(self, name: Optional[str] = None) -> str:
        doc = self._doc(name)
        objs = self._get(doc, "doc.objects", default=None)
        err: Optional[str] = None
        if objs is not None:
            try:
                text = self._data(objs, MIME["cdxml"], doc)
                if text.lstrip().startswith("<"):
                    return text.replace("\r\n", "\n")
                err = "Objects.Data(chemical/x-cdxml) returned non-XML data"
            except ChemDrawError as exc:
                err = str(exc)
        tmpdir = tempfile.mkdtemp(prefix="chemdrawmcp-")
        tmp = os.path.join(tmpdir, "document.cdxml")
        try:
            self._invoke(doc, "doc.save_as", tmp)
            with open(tmp, "r", encoding="utf-8", errors="replace") as fh:
                return fh.read().replace("\r\n", "\n")
        except (OSError, ChemDrawError) as exc:
            raise ChemDrawError("could not read the document as CDXML (%s; SaveAs fallback: %s)" % (err, exc)) from exc
        finally:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
                os.rmdir(tmpdir)
            except OSError:
                pass

    # ------------------------------------------------------------------ settings
    def _settings_obj(self, doc: Any) -> Any:
        return self._get(doc, "doc.settings", default=None)

    def _settings_from_com(self, settings: Any) -> Dict[str, Any]:
        out: Dict[str, Any] = {}
        unavailable = []
        for key, names in SETTINGS_MEMBERS.items():
            val = _MISSING
            for n in names:
                try:
                    val = getattr(settings, n)
                    if _is_method(val):
                        val = val()
                    break
                except Exception:
                    continue
            if val is _MISSING:
                unavailable.append(key)
                continue
            out[key] = _num(val) if key in _NUMERIC_SETTINGS and key not in ("label_face", "caption_face") else val
        if not out:
            raise ChemDrawError("the Settings object exposes none of %s; edit SETTINGS_MEMBERS in "
                                "chemdrawmcp/backends/windows.py or set them in the CDXML instead" % list(SETTINGS_MEMBERS))
        if unavailable:
            out["unavailable"] = unavailable
        out["source"] = "com"
        return out

    @_com_method()
    def document_settings(self, name: Optional[str] = None) -> Dict[str, Any]:
        doc = self._doc(name)
        settings = self._settings_obj(doc)
        if settings is not None:
            out = self._settings_from_com(settings)
        else:
            cdxml = self._cdxml_or_empty(self._get(doc, "doc.objects", default=None), doc)
            out = cdxml_settings(cdxml)
            if not out:
                raise ChemDrawError("document settings are not available: this ChemDraw COM interface has no Settings "
                                    "object (MEMBERS['doc.settings']) and the document CDXML carries no drawing "
                                    "attributes (BondLength, LineWidth, LabelFont, ... on the <CDXML> element)")
            missing = [k for k in SETTINGS_MEMBERS if k not in out]
            if missing:
                out["unavailable"] = missing
            out["source"] = "cdxml"
        out["name"] = self._doc_name(doc)
        return out

    @_com_method()
    def set_document_settings(self, name: Optional[str], settings: Dict[str, Any]) -> Dict[str, Any]:
        doc = self._doc(name)
        errors = []
        camel = {re.sub(r"_(\w)", lambda m: m.group(1).upper(), k).lower(): k for k in SETTINGS_MEMBERS}
        wanted: Dict[str, Any] = {}
        for raw_key, value in (settings or {}).items():
            key = str(raw_key).strip()
            key = key if key in SETTINGS_MEMBERS else camel.get(key.lower(), re.sub(r"(?<!^)(?=[A-Z])", "_", key).lower())
            if key not in SETTINGS_MEMBERS:
                errors.append("%s: unknown setting (known: %s)" % (raw_key, ", ".join(SETTINGS_MEMBERS)))
                continue
            wanted[key] = value
        obj = self._settings_obj(doc)
        if obj is not None:
            for key, value in wanted.items():
                try:
                    self._set(obj, key, value, names=SETTINGS_MEMBERS[key])
                except ChemDrawError as exc:
                    errors.append("%s: %s" % (key, exc))
        elif wanted:
            # No Settings object: rewrite the <CDXML> root attributes and set the whole document again.
            # That only works when SetData replaces the content; if it adds, everything would be duplicated.
            if self.setdata_replaces is not True:
                raise ChemDrawError(
                    "cannot change document settings: this ChemDraw COM interface has no Settings object "
                    "(MEMBERS['doc.settings']) and rewriting the CDXML needs Objects.SetData to *replace* the document, "
                    "which is %s. Insert something first (insert() learns the behaviour), set INSERT_MODE = 'replace' "
                    "in chemdrawmcp/backends/windows.py, or generate the CDXML with these settings and open_document() it."
                    % ("not known yet" if self.setdata_replaces is None else "not the case (it adds)"))
            objs = self._objects(doc)
            cdxml = self._cdxml_or_empty(objs, doc)
            if not cdxml:
                raise ChemDrawError("cannot read the document as CDXML to change its settings")
            new = cdxml_with_settings(cdxml, wanted)
            self._invoke(objs, "objs.set_data", MIME["cdxml"], new)
        result = self.document_settings(name)
        if errors:
            result["errors"] = errors
        return result

    # ------------------------------------------------------------------ content
    def _insert_payloads(self, kind: str, data: str) -> List[Tuple[str, str]]:
        if kind == "text":
            return [(MIME["cdxml"], caption_cdxml(data))]
        if kind in SEQUENCE_KINDS:
            payloads = [(MIME["helm"], sequence_to_helm(data, kind))]
            payloads += [(m, data) for m in SEQUENCE_FALLBACK_MIME.get(kind, [])]
            return payloads
        return [(INSERT_MIME[kind], data)]

    def _paste_text(self, doc: Any, objs: Any, text: str) -> str:
        """Clipboard fallback: CF_UNICODETEXT on the clipboard, then ``Document.Paste()`` / ``Objects.Paste()``.

        Plain text pasted into ChemDraw becomes a caption unless ChemDraw recognises the format, so this
        is a last resort; the previous clipboard text is restored afterwards.
        """
        old = _clipboard_get_text()
        _clipboard_set_text(text)
        try:
            errors: List[str] = []
            for holder, key in ((doc, "doc.paste"), (objs, "objs.paste")):
                try:
                    self._invoke(holder, key)
                    return key
                except ChemDrawError as exc:
                    errors.append(str(exc))
            raise ChemDrawError("no Paste member found (%s)" % " / ".join(errors))
        finally:
            if old is not None:
                try:
                    _clipboard_set_text(old)
                except ChemDrawError:
                    pass

    @staticmethod
    def _insert_outcome(before: Dict[str, int], after: Dict[str, int], before_cdxml: str, after_cdxml: str) -> str:
        """"added" | "replaced" | "nothing": what SetData did to the document."""
        tb, ta = _total_objects(before), _total_objects(after)
        if tb == 0:
            return "added" if ta > 0 else "nothing"
        compared = False
        for with_box in (True, False):               # positions first; then ignoring re-layout / re-rounding
            old_sigs = cdxml_object_signatures(before_cdxml, with_box)
            if not old_sigs:
                break
            compared = True
            remaining = list(cdxml_object_signatures(after_cdxml, with_box))
            missing = False
            for sig in old_sigs:
                if sig in remaining:
                    remaining.remove(sig)
                else:
                    missing = True
                    break
            if not missing:
                return "added" if remaining else "nothing"
        if compared:
            return "replaced"                         # old objects are gone even ignoring positions
        if ta > tb:                                   # no CDXML to compare: fall back to the counts
            return "added"
        return "nothing" if ta == tb else "replaced"

    @_com_method()
    def insert(self, name: Optional[str], data: str, kind: str = "cdxml") -> dict:
        kind = (kind or "cdxml").strip().lower()
        if kind not in INSERT_KINDS:
            raise ChemDrawError("unknown insert kind %r; use one of %s" % (kind, ", ".join(INSERT_KINDS)))
        if not (data or "").strip():
            raise ChemDrawError("nothing to insert (empty data)")
        doc = self._doc(name)
        objs = self._objects(doc)
        target = objs
        if INSERT_TARGET == "selection":
            target = self._get(doc, "doc.selection", default=None) or objs
        before_cdxml = self._cdxml_or_empty(objs, doc)
        before = self._counts(objs, cdxml=before_cdxml)
        payloads = self._insert_payloads(kind, data)
        used: Optional[str] = None
        errors: List[str] = []
        n_dialogs = len(self.last_dialogs)        # dialogs the waiting thread dismisses during this call come after
        for mime, text in payloads:
            for m in _mime_candidates(mime):
                try:
                    self._invoke(target, "objs.set_data", m, text)
                    used = m
                    break
                except ChemDrawError as exc:
                    errors.append("%s: %s" % (m, exc))
            if used:
                break
        if used is None:
            try:
                used = "clipboard via " + self._paste_text(doc, objs, payloads[0][1])
            except ChemDrawError as exc:
                errors.append("clipboard: %s" % exc)
                raise ChemDrawError("insert of %s failed: %s. If ChemDraw wants another MIME type, edit INSERT_MIME / "
                                    "MIME_ALTERNATIVES in chemdrawmcp/backends/windows.py." % (kind, " / ".join(errors)))
        after_cdxml = self._cdxml_or_empty(objs, doc)
        after = self._counts(objs, cdxml=after_cdxml)
        outcome = self._insert_outcome(before, after, before_cdxml, after_cdxml)
        if INSERT_MODE == "add" and outcome == "replaced":
            outcome = "added"
        if INSERT_MODE == "replace" and _total_objects(before) > 0 and outcome != "nothing":
            outcome = "replaced"
        result: Dict[str, Any] = {"document": self._doc_name(doc), "kind": kind, "mime": used, "mode": outcome}
        if kind in SEQUENCE_KINDS:
            result["helm"] = payloads[0][1]
        if outcome == "replaced":
            # SetData replaced the document: set the old content merged with the new one.
            self.setdata_replaces = True
            try:
                merged = merge_cdxml(before_cdxml, after_cdxml)
                self._invoke(target, "objs.set_data", MIME["cdxml"], merged)
                after_cdxml = self._cdxml_or_empty(objs, doc)
                after = self._counts(objs, cdxml=after_cdxml)
                result["mode"] = "replaced+merged"
                result["note"] = ("Objects.SetData replaces the document on this ChemDraw; the previous content was "
                                  "merged back in (new objects placed below the old ones)")
            except ChemDrawError as exc:
                result["warning"] = ("Objects.SetData replaced the document content and merging the previous content "
                                     "back failed: %s" % exc)
        elif outcome == "added" and _total_objects(before) > 0 and self.setdata_replaces is None:
            self.setdata_replaces = False
        alerts = list(self.last_dialogs[n_dialogs:])
        if outcome == "nothing" and kind != "text":
            result["warning"] = "no objects were added" + ((": " + " | ".join(alerts)) if alerts else "")
        result.update({
            "atoms_before": before["atoms"], "atoms_after": after["atoms"], "atoms_added": after["atoms"] - before["atoms"],
            "bonds_before": before["bonds"], "bonds_after": after["bonds"],
            "molecules": after["molecules"], "captions": after["captions"], "alerts": alerts,
        })
        if errors:
            result["fallbacks"] = errors
        return result

    @_com_method()
    def select_all(self, name: Optional[str] = None) -> None:
        doc = self._doc(name)
        errors: List[str] = []
        try:
            self._invoke(doc, "doc.select_all")
            return
        except ChemDrawError as exc:
            errors.append(str(exc))
        objs = self._objects(doc)
        try:
            self._invoke(objs, "objs.select")
            return
        except ChemDrawError as exc:
            errors.append(str(exc))
        try:
            self._set(objs, "objs.selected", True)
            return
        except ChemDrawError as exc:
            errors.append(str(exc))
        raise ChemDrawError("select_all is not available through ChemDraw COM on this version (%s)" % " / ".join(errors))

    @_com_method()
    def deselect(self, name: Optional[str] = None) -> None:
        doc = self._doc(name)
        errors: List[str] = []
        try:
            self._invoke(doc, "doc.deselect")
            return
        except ChemDrawError as exc:
            errors.append(str(exc))
        try:
            self._set(self._objects(doc), "objs.selected", False)
            return
        except ChemDrawError as exc:
            errors.append(str(exc))
        raise ChemDrawError("deselect is not available through ChemDraw COM on this version (%s)" % " / ".join(errors))

    @_com_method()
    def selection_info(self, name: Optional[str] = None) -> Dict[str, Any]:
        doc = self._doc(name)
        target, selected = self._target(doc)
        info: Dict[str, Any] = {"document": self._doc_name(doc), "selected": selected,
                                "scope": "selection" if selected else "document"}
        cdxml = self._cdxml_or_empty(target, None if selected else doc)
        summary = cdxml_summary(cdxml) if cdxml else {}
        info.update(self._counts(target, cdxml=cdxml))
        info["empty"] = _total_objects(info) == 0
        unavailable: List[str] = []

        def grab(key: str, fn: Callable[[], Any]) -> None:
            try:
                val = fn()
            except Exception:
                unavailable.append(key)
                return
            if val is None or (isinstance(val, str) and not val.strip()):
                unavailable.append(key)
            else:
                info[key] = val

        holder = None if selected else doc
        grab("smiles", lambda: self._data(target, MIME["smiles"], holder).strip())
        grab("formula", lambda: _as_text(self._get(target, "objs.formula")).strip())
        grab("molecular_weight", lambda: _num(self._get(target, "objs.molecular_weight")))
        grab("exact_mass", lambda: _num(self._get(target, "objs.exact_mass")))
        grab("elemental_analysis", lambda: _as_text(self._get(target, "objs.elemental_analysis")).strip())
        grab("text", lambda: summary.get("text") or self._caption_text(target))
        grab("bounds", lambda: summary.get("bounds") or self._bounds(target))
        if unavailable:
            info["unavailable"] = unavailable
        return info

    def _caption_text(self, target: Any) -> Optional[str]:
        """Caption text through the ``Captions`` collection (CDXML parsing is preferred, see selection_info)."""
        caps = self._get(target, "objs.captions", default=None)
        if caps is None:
            return None
        texts: List[str] = []
        for it in self._items(caps):
            t = self._get(it, "obj.text", default=None)
            if isinstance(t, str) and t.strip():
                texts.append(t.strip())
        return "\n".join(texts) if texts else None

    def _bounds(self, target: Any) -> Optional[Dict[str, float]]:
        rect = self._get(target, "objs.bounds", default=None)
        if rect is None:
            return None
        if isinstance(rect, (tuple, list)) and len(rect) == 4:
            l, t, r, b = [float(v) for v in rect]
        else:
            vals = [self._get(rect, k, default=None) for k in ("rect.left", "rect.top", "rect.right", "rect.bottom")]
            if any(v is None for v in vals):
                return None
            l, t, r, b = [float(v) for v in vals]
        return {"left": l, "top": t, "right": r, "bottom": b, "width": r - l, "height": b - t}

    @_com_method()
    def selection_cdxml(self, name: Optional[str] = None) -> str:
        doc = self._doc(name)
        target, selected = self._target(doc)
        return self._data(target, MIME["cdxml"], None if selected else doc).replace("\r\n", "\n")

    @_com_method()
    def selection_as(self, name: Optional[str], what: str) -> str:
        key = (what or "").strip().lower().replace("-", "").replace("_", "").replace(" ", "")
        key = {"inchikey": "inchikey", "molfile": "mol", "v3000": "molv3000", "mol3000": "molv3000", "molfilev3000": "molv3000",
               "chemicalname": "name", "iupac": "name", "iupacname": "name"}.get(key, key)
        if key not in SELECTION_FORMATS:
            raise ChemDrawError("unknown conversion %r; use one of %s" % (what, ", ".join(SELECTION_FORMATS)))
        doc = self._doc(name)
        target, selected = self._target(doc)
        holder = None if selected else doc
        err: str
        try:
            text = self._data(target, SELECTION_FORMATS[key], holder).strip()
            if text:
                return text.replace("\r\n", "\n")
            err = "ChemDraw returned nothing for %s" % SELECTION_FORMATS[key]
        except ChemDrawError as exc:
            err = str(exc)
        if key == "name":
            chem_name = _as_text(self._get(target, "objs.name", default="")).strip()
            if chem_name:
                return chem_name
        fallback = self._convert_fallback(target, key, holder)
        if fallback:
            return fallback
        raise ChemDrawError("could not convert the selection to %s: %s" % (key, err))

    def _convert_fallback(self, target: Any, key: str, doc: Any = None) -> Optional[str]:
        """InChIKey / V3000 molfile etc. via ChemScript or RDKit, from whatever ChemDraw *can* produce."""
        try:
            from rdkit import Chem  # type: ignore
        except ImportError:
            Chem = None
        for in_mime in (MIME["mol"], MIME["smiles"], MIME["cdxml"]):
            try:
                text = self._data(target, in_mime, doc).strip()
            except ChemDrawError:
                continue
            if not text:
                continue
            if self.chemscript.available():
                try:
                    out = self.chemscript.convert(text, in_mime, SELECTION_FORMATS[key]).strip()
                    if out:
                        return out
                except ChemDrawError:
                    pass
            if Chem is None:
                continue
            mol = None
            if in_mime == MIME["mol"]:
                mol = Chem.MolFromMolBlock(text)
            elif in_mime == MIME["smiles"]:
                mol = Chem.MolFromSmiles(text)
            if mol is None:
                continue
            try:
                if key == "inchikey":
                    out = Chem.MolToInchiKey(mol)
                elif key == "inchi":
                    out = Chem.MolToInchi(mol)
                elif key == "smiles":
                    out = Chem.MolToSmiles(mol)
                elif key == "mol":
                    out = Chem.MolToMolBlock(mol)
                elif key == "molv3000":
                    out = Chem.MolToMolBlock(mol, forceV3000=True)
                else:
                    return None
            except Exception:
                continue
            if out:
                return out
        return None

    # ------------------------------------------------------------------ commands & menus
    @staticmethod
    def _norm_cmd(s: str) -> str:
        return re.sub(r"[^a-z0-9]", "", (s or "").lower())

    def _resolve_command(self, command: str) -> str:
        want = self._norm_cmd(command)
        for name in COMMANDS:
            if self._norm_cmd(name) == want:
                return name
        for alias, name in COMMAND_ALIASES.items():
            if self._norm_cmd(alias) == want:
                return name
        raise ChemDrawError("command %r is not available through ChemDraw COM on Windows; available: %s (edit COMMANDS "
                            "in chemdrawmcp/backends/windows.py to add one)" % (command, ", ".join(COMMANDS)))

    def list_commands(self, pattern: Optional[str] = None, enabled_only: bool = False,
                      with_descriptions: bool = True) -> List[CommandInfo]:
        rx = None
        if pattern:
            pat = pattern.strip()
            if any(ch in pat for ch in "*?[") and not pat.startswith("^"):
                import fnmatch
                rx = re.compile(fnmatch.translate(pat), re.IGNORECASE)
            else:
                try:
                    rx = re.compile(pat, re.IGNORECASE)
                except re.error:
                    rx = re.compile(re.escape(pat), re.IGNORECASE)
        out = []
        for name, (short, long_) in COMMANDS.items():
            if rx and not (rx.search(name) or rx.search(short)):
                continue
            out.append(CommandInfo(name=name, enabled=None, short_description=short if with_descriptions else "",
                                   long_description=long_ if with_descriptions else ""))
        return out

    @_com_method()
    def run_command(self, command: str, name: Optional[str] = None) -> dict:
        canonical = self._resolve_command(command)
        n_dialogs = len(self.last_dialogs)
        handler = getattr(self, "_cmd_" + canonical)
        result = handler(name)
        alerts = list(self.last_dialogs[n_dialogs:])
        return {"ok": True, "command": canonical, "enabled_before": None, "alerts": alerts, "result": result}

    # command handlers ------------------------------------------------------------------
    def _cmd_newDocument(self, name: Optional[str]) -> Any:
        return self.new_document()

    def _cmd_closeDocument(self, name: Optional[str]) -> Any:
        self.close_document(name, save=False)
        return "closed"

    def _cmd_saveDocument(self, name: Optional[str]) -> Any:
        doc = self._doc(name)
        path = self._doc_path(doc)
        if not path:
            raise ChemDrawError("document %r has never been saved; use save_document with a path" % self._doc_name(doc))
        self._invoke(doc, "doc.save")
        return path

    def _cmd_selectAll(self, name: Optional[str]) -> Any:
        self.select_all(name)
        return "selected"

    def _cmd_deselectAll(self, name: Optional[str]) -> Any:
        self.deselect(name)
        return "deselected"

    def _cmd_cleanStructure(self, name: Optional[str]) -> Any:
        doc = self._doc(name)
        target, selected = self._target(doc)
        self._invoke(target, "objs.clean")
        return "cleaned %s" % ("selection" if selected else "document")

    def _cmd_convertNameToStructure(self, name: Optional[str]) -> Any:
        doc = self._doc(name)
        target, selected = self._target(doc)
        text = None
        if selected:
            cdxml = self._cdxml_or_empty(target)
            text = (cdxml_summary(cdxml).get("text") if cdxml else None) or self._caption_text(target)
        if not text:
            raise ChemDrawError("select a caption containing a chemical name first (or use insert(kind='name'))")
        return self.insert(self._doc_name(doc), text, kind="name")

    def _cmd_convertStructureToName(self, name: Optional[str]) -> Any:
        doc = self._doc(name)
        chem_name = self.selection_as(self._doc_name(doc), "name")
        out: Dict[str, Any] = {"name": chem_name}
        try:
            self.insert(self._doc_name(doc), chem_name, kind="text")
            out["caption_inserted"] = True
        except ChemDrawError as exc:
            out["caption_inserted"] = False
            out["caption_error"] = str(exc)
        return out

    def _cmd_copyAsSMILES(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "smiles")

    def _cmd_copyAsInChI(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "inchi")

    def _cmd_copyAsInChIKey(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "inchikey")

    def _cmd_copyAsMolfile(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "mol")

    def _cmd_copyAsMolfileV3000(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "molv3000")

    def _cmd_copyAsCDXML(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "cdxml")

    def _cmd_copyAsHELM(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "helm")

    def _cmd_copyAsSLN(self, name: Optional[str]) -> Any:
        return self.selection_as(name, "sln")

    def _cmd_chemicalProperties(self, name: Optional[str]) -> Any:
        return self.selection_info(name)

    def _cmd_elementalAnalysis(self, name: Optional[str]) -> Any:
        doc = self._doc(name)
        target, _ = self._target(doc)
        return _as_text(self._get(target, "objs.elemental_analysis")).strip()

    def _cmd_dismissDialogs(self, name: Optional[str]) -> Any:
        return self.dismiss_dialogs()

    def menus(self) -> Dict[str, List[str]]:
        raise ChemDrawError("ChemDraw's menus are not scriptable through COM on Windows; use run_command "
                            "(list_commands shows what is available) or drive the CDXML directly")

    def do_menu(self, menu: str, item: str, name: Optional[str] = None) -> dict:
        raise ChemDrawError("menu item %r > %r cannot be invoked on Windows: ChemDraw exposes no menu automation "
                            "through COM; use run_command (see list_commands) instead" % (menu, item))

    # ------------------------------------------------------------------ pictures
    def _main_window(self, name: Optional[str]) -> Optional[int]:
        pids = _chemdraw_pids()
        wins = [(h, c, t) for h, c, t in _chemdraw_top_windows(pids) if not _is_dialog_class(c)]
        if not wins:
            return None
        if name:
            want = str(name).lower()
            for h, _c, t in wins:
                if want in t.lower():
                    return h
        for h, _c, t in wins:
            if "chemdraw" in t.lower():
                return h

        def area(h: int) -> int:
            l, t, r, b = _window_rect(h)
            return max(0, r - l) * max(0, b - t)

        return max((h for h, _c, _t in wins), key=area)

    @_com_method()
    def screenshot(self, name: Optional[str] = None) -> bytes:
        u = _user32()
        if u is None:
            raise ChemDrawNotAvailable("screenshots of ChemDraw windows are only possible on Windows")
        try:
            from PIL import ImageGrab  # type: ignore
        except ImportError as exc:
            raise ChemDrawError("screenshot needs Pillow (pip install pillow): %s" % exc) from exc
        if name:
            try:
                self.activate_document(name)
            except ChemDrawError:
                pass
        hwnd = self._main_window(name)
        if not hwnd:
            raise ChemDrawError("no visible ChemDraw window found (is ChemDraw running and not minimized to tray?)")
        try:
            u.SetProcessDPIAware()
        except Exception:
            pass
        if u.IsIconic(hwnd):
            u.ShowWindow(hwnd, _SW_RESTORE)
        u.SetForegroundWindow(hwnd)
        time.sleep(0.4)
        l, t, r, b = _window_rect(hwnd)
        if r <= l or b <= t:
            raise ChemDrawError("ChemDraw window has an empty rectangle %r" % ((l, t, r, b),))
        try:
            img = ImageGrab.grab(bbox=(l, t, r, b), all_screens=True)
        except TypeError:
            img = ImageGrab.grab(bbox=(l, t, r, b))
        if img.mode not in ("RGB", "RGBA"):
            img = img.convert("RGB")
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()

    # ------------------------------------------------------------------ info
    def capabilities(self) -> Dict[str, Any]:
        caps = super().capabilities()
        libname = None
        try:
            libname = _com_library()[0]
        except ChemDrawNotAvailable:
            pass
        caps.update({
            "com_library": libname,
            "com_thread": self._worker is not None,
            "prog_ids": list(self.prog_ids),
            "prog_id": self.prog_id,
            "connect_note": self.connect_note,
            "chemscript": self.chemscript.available(),
            "chemscript_module": self.chemscript.module_name,
            "running": self.is_running() if _is_windows() or self._app_factory else False,
            "menus": False,
            "commands": list(COMMANDS),
            "insert_mode": INSERT_MODE,
            "insert_target": INSERT_TARGET,
            "setdata_replaces": self.setdata_replaces,
            "notes": [
                "COM member names are looked up through MEMBERS in chemdrawmcp/backends/windows.py; error messages name the entry to edit",
                "all COM calls run on one STA worker thread; modal dialogs are dismissed while a call is pending",
                "insert() detects whether Objects.SetData adds or replaces (setdata_replaces) and merges the old content back when it replaces",
                "save_document sniffs the written file and falls back to Objects.Data(mime) when SaveAs ignored the extension",
                "menus()/do_menu() are not supported on Windows; run_command offers the curated COMMANDS",
            ],
        })
        return caps


__all__ = [
    "WindowsBackend", "ChemScriptBridge", "chemscript_available", "find_chemdraw_windows", "find_stationery",
    "sequence_to_helm", "caption_cdxml", "cdxml_summary", "cdxml_object_signatures", "merge_cdxml", "cdxml_settings",
    "cdxml_with_settings",
    "PROG_IDS", "MIME", "MIME_ALTERNATIVES", "INSERT_MIME", "SELECTION_FORMATS", "EXPORT_MIME", "MEMBERS",
    "SETTINGS_MEMBERS", "CDXML_SETTINGS", "COMMANDS", "COMMAND_ALIASES", "CHEMSCRIPT_MODULES", "CHEMSCRIPT_MEMBERS",
    "INSERT_MODE", "INSERT_TARGET", "DIALOG_CLASSES", "DEFAULT_DIALOG_BUTTONS",
]
