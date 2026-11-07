"""Platform-independent tests for the Windows (COM) backend.

Nothing here needs Windows, pywin32 or ChemDraw: the COM object model is replaced by a fake
application that models the most likely real one (``Documents`` with Add/Open/Item/Count,
``Document`` with Objects/Selection/SaveAs/Close/Name/FullName/Path/Saved, ``Objects`` with
Data/SetData/Count/Item/Clean and the property members) and records every member the backend calls.
``FakeApp(replace=True)`` models a ChemDraw whose ``Objects.SetData`` *replaces* the document.
"""
from __future__ import annotations

import inspect
import os
import sys
import threading
import time
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

import pytest

from chemdrawmcp.backends import windows
from chemdrawmcp.backends.base import ChemDrawBackend, ChemDrawError, ChemDrawNotAvailable, EXPORT_FORMATS, INSERT_KINDS
from chemdrawmcp.backends.windows import WindowsBackend

NOT_WINDOWS = sys.platform != "win32"

CDXML = "chemical/x-cdxml"
SMILES = "chemical/x-daylight-smiles"
NAME = "chemical/x-name"
MOL = "chemical/x-mdl-molfile"

ETHANOL_CDXML = (
    '<?xml version="1.0" encoding="UTF-8" ?>\n'
    '<CDXML BondLength="14.4" LabelFont="3" LabelSize="10" CaptionFont="3" CaptionSize="10" LineWidth="0.6">'
    '<fonttable><font id="3" charset="iso-8859-1" name="Arial"/></fonttable><page>'
    '<fragment id="10" BoundingBox="100 100 150 130" smiles="CCO">'
    '<n id="1" p="100 100"/><n id="2" p="120 110"/><n id="3" p="140 100" Element="8"><t><s>OH</s></t></n>'
    '<b id="4" B="1" E="2"/><b id="5" B="2" E="3"/>'
    '</fragment>'
    '<t id="20" p="100 160" BoundingBox="100 150 160 165"><s font="3" size="10">ethanol</s></t>'
    '</page></CDXML>'
)

#: smiles -> (formula, MW, exact mass, name, atoms, bonds)
KNOWN = {
    "CCO": ("C2H6O", "46.07", 46.0419, "ethanol", 3, 2),
    "c1ccccc1": ("C6H6", "78.11", 78.0470, "benzene", 6, 6),
    "CC(=O)O": ("C2H4O2", "60.05", 60.0211, "acetic acid", 4, 3),
}
NAMES = {v[3]: k for k, v in KNOWN.items()}


# ---------------------------------------------------------------------------- fake COM objects
class FakeItem:
    def __init__(self, kind, smiles=None, text=None, atoms=0, bonds=0):
        self.kind, self.smiles, self.text, self.atoms, self.bonds = kind, smiles, text, atoms, bonds

    @classmethod
    def fragment(cls, smiles):
        if smiles not in KNOWN:
            raise ValueError("Cannot interpret SMILES %r" % smiles)
        return cls("fragment", smiles=smiles, atoms=KNOWN[smiles][4], bonds=KNOWN[smiles][5])


DEFAULT_ATTRS = {"BondLength": "14.4", "LabelFont": "3", "LabelSize": "10", "CaptionFont": "3", "CaptionSize": "10",
                 "LineWidth": "0.6"}
DEFAULT_FONTS = {"3": "Arial"}


def render_cdxml(items, attrs=None, fonts=None):
    attrs = dict(DEFAULT_ATTRS if attrs is None else attrs)
    fonts = dict(DEFAULT_FONTS if fonts is None else fonts)
    parts, nid, y = [], 10, 100
    for it in items:
        if it.kind == "fragment":
            atoms = "".join('<n id="%d" p="%g %g"/>' % (nid + k + 1, 100 + 20 * k, y) for k in range(it.atoms))
            bonds = "".join('<b id="%d" B="%d" E="%d"/>' % (nid + it.atoms + k + 1, nid + k + 1, nid + (k + 1) % it.atoms + 1)
                            for k in range(it.bonds))          # chain, closing into a ring when bonds == atoms
            smi = ' smiles="%s"' % it.smiles if it.smiles else ""
            parts.append('<fragment id="%d" BoundingBox="100 %g %g %g"%s>%s%s</fragment>'
                         % (nid, y, 100 + 20 * max(it.atoms, 1), y + 30, smi, atoms, bonds))
            nid += it.atoms + it.bonds + 1
        else:
            parts.append('<t id="%d" p="100 %g" BoundingBox="100 %g 160 %g"><s font="3" size="10">%s</s></t>'
                         % (nid, y + 10, y, y + 15, escape(it.text or "")))
            nid += 1
        y += 60
    attr_text = " ".join('%s="%s"' % (k, escape(str(v))) for k, v in attrs.items())
    font_text = "".join('<font id="%s" charset="iso-8859-1" name="%s"/>' % (i, escape(n)) for i, n in fonts.items())
    return ('<?xml version="1.0" encoding="UTF-8" ?>\n'
            '<CDXML %s><fonttable>%s</fonttable><page>%s</page></CDXML>' % (attr_text, font_text, "".join(parts)))


def parse_cdxml(cdxml):
    """(items, root attributes, fonts) of a CDXML text."""
    root = ET.fromstring(cdxml.replace('<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >', ""))
    fonts = {f.get("id"): f.get("name") for f in root.iter("font")}
    return parse_items(cdxml), dict(root.attrib), fonts


def parse_items(cdxml):
    root = ET.fromstring(cdxml.replace('<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >', ""))
    page = root.find("page")
    out = []
    for el in (list(page) if page is not None else []):
        if el.tag == "fragment":
            out.append(FakeItem("fragment", smiles=el.get("smiles"), atoms=len(el.findall("n")), bonds=len(el.findall("b"))))
        elif el.tag == "t":
            out.append(FakeItem("caption", text="".join(s.text or "" for s in el.iter("s"))))
    return out


class FakeObjects:
    """ChemDraw ``Objects`` / ``Selection``: Count, Item(i), Data(mime), SetData(mime, data), Clean() + properties."""

    def __init__(self, doc, items=None, replace=False):
        self._doc = doc
        self.items = list(items or [])
        self.replace = replace
        self.calls = []
        self.threads = set()

    def __setattr__(self, name, value):      # like pywin32: assigning an unknown COM property raises
        if not name.startswith("_") and not hasattr(self, name) and name not in ("items", "replace", "calls", "threads"):
            raise AttributeError("fake ChemDraw Objects has no property %r" % name)
        object.__setattr__(self, name, value)

    @property
    def Count(self):
        return len(self.items)

    def Item(self, i):
        return self.items[i - 1]

    def _fragments(self):
        return [it for it in self.items if it.kind == "fragment" and it.smiles]

    def Data(self, mime):
        self.calls.append(("Data", mime))
        self.threads.add(threading.get_ident())
        if mime == CDXML:
            return render_cdxml(self.items, self._doc.attrs, self._doc.fonts)
        frags = self._fragments()
        if mime == SMILES:
            if not frags:
                raise ValueError("fake ChemDraw: no structure for SMILES")
            return ".".join(f.smiles for f in frags)
        if mime == NAME:
            if len(frags) != 1:
                raise ValueError("fake ChemDraw: Struct=Name needs one structure")
            return KNOWN[frags[0].smiles][3]
        if mime == MOL:
            if len(frags) != 1:
                raise ValueError("fake ChemDraw: one structure please")
            f = frags[0]
            return "\n  fake\n\n%3d%3d  0  0  0  0  0  0  0  0999 V2000\nM  END\n" % (f.atoms, f.bonds)
        if mime == "image/svg+xml":
            return '<svg xmlns="http://www.w3.org/2000/svg"/>'
        if mime == "application/pdf":
            return b"%PDF-1.4 fake"
        raise ValueError("fake ChemDraw cannot produce %s" % mime)

    def SetData(self, mime, data):
        self.calls.append(("SetData", mime, data))
        self.threads.add(threading.get_ident())
        if mime == SMILES:
            new = [FakeItem.fragment(data)]
        elif mime == NAME:
            if data not in NAMES:
                raise ValueError("Cannot interpret name %r" % data)
            new = [FakeItem.fragment(NAMES[data])]
        elif mime == CDXML:
            new, attrs, fonts = parse_cdxml(data)
            if self.replace:                 # a replaced document takes the settings of the new CDXML
                self._doc.attrs, self._doc.fonts = attrs, fonts
        else:
            raise ValueError("fake ChemDraw refuses %s" % mime)
        self.items = new if self.replace else self.items + new
        self._doc.Saved = False

    def Clean(self):
        self.calls.append(("Clean",))

    def _known(self):
        frags = self._fragments()
        if not frags:
            raise AttributeError("no structure")
        return KNOWN[frags[0].smiles]

    Formula = property(lambda self: self._known()[0])
    MolecularWeight = property(lambda self: self._known()[1])
    ExactMass = property(lambda self: self._known()[2])
    ElementalAnalysis = property(lambda self: "C, 52.14; H, 13.13; O, 34.73" if self._known()[0] == "C2H6O" else "n/a")


class FakeDocument:
    def __init__(self, app, name, path=None, replace=False):
        self._app = app
        self.Name = name
        self._path = path
        self.Saved = True
        self.attrs = dict(DEFAULT_ATTRS)
        self.fonts = dict(DEFAULT_FONTS)
        self.Objects = FakeObjects(self, replace=replace)
        self.Selection = FakeObjects(self)
        self.calls = []
        self.closed = False

    @property
    def FullName(self):
        return self._path or self.Name

    @property
    def Path(self):
        return os.path.dirname(self._path) if self._path else ""

    def SaveAs(self, path, *extra):
        self.calls.append(("SaveAs", path) + extra)
        ext = os.path.splitext(path)[1].lower()
        if ext == ".svg":
            raise ValueError("fake ChemDraw cannot save as SVG")
        if ext in (".cdxml", ".ct"):
            data = self.Objects.Data(CDXML).encode("utf-8")
        elif ext == ".mol":
            data = self.Objects.Data(MOL).encode("utf-8")
        elif ext == ".png":
            data = b"\x89PNG\r\n\x1a\nfake"
        elif ext == ".jpeg":
            data = b"\xff\xd8\xff\xe0fake"
        else:                                   # .cdx, and formats where ChemDraw "ignores" the extension
            data = b"VjCD0100fake"
        with open(path, "wb") as fh:
            fh.write(data)
        self._path = path
        self.Name = os.path.basename(path)
        self.Saved = True

    def Save(self):
        self.calls.append(("Save",))
        self.Saved = True

    def Close(self, save_changes=False):
        self.calls.append(("Close", save_changes))
        self.closed = True
        self._app.Documents._docs.remove(self)
        if self._app.ActiveDocument is self:
            self._app.ActiveDocument = self._app.Documents._docs[-1] if self._app.Documents._docs else None

    def Activate(self):
        self.calls.append(("Activate",))
        self._app.ActiveDocument = self

    def SelectAll(self):
        self.calls.append(("SelectAll",))
        self.Selection = FakeObjects(self, items=self.Objects.items)


class FakeDocuments:
    def __init__(self, app):
        self._app = app
        self._docs = []
        self.calls = []

    @property
    def Count(self):
        return len(self._docs)

    def Item(self, i):          # 1-based like real COM collections
        if i < 1 or i > len(self._docs):
            raise IndexError(i)
        return self._docs[i - 1]

    def Add(self):
        self.calls.append(("Add",))
        self._app.serial += 1
        doc = FakeDocument(self._app, "Untitled Document-%d" % self._app.serial, replace=self._app.replace)
        self._docs.append(doc)
        self._app.ActiveDocument = doc
        return doc

    def Open(self, path):
        self.calls.append(("Open", path))
        doc = FakeDocument(self._app, os.path.basename(path), path, replace=self._app.replace)
        with open(path, encoding="utf-8") as fh:
            doc.Objects.items = parse_items(fh.read())
        self._docs.append(doc)
        self._app.ActiveDocument = doc
        return doc


class FakeApp:
    def __init__(self, replace=False):
        self.replace = replace
        self.serial = 0
        self.Documents = FakeDocuments(self)
        self.ActiveDocument = None
        self.Visible = False
        self.quit_called = False
        self.version_gate = None          # a threading.Event: Version blocks until it is set (models a hung ChemDraw)

    @property
    def Version(self):
        gate = self.version_gate
        if gate is not None:
            gate.wait(30.0)
        return "23.1.2.7"

    def Quit(self):
        self.quit_called = True


class ComErrorLike(Exception):
    """Shaped like pywin32's ``com_error`` (type name and attributes are what the backend looks at)."""

    def __init__(self, hresult, text, source="ChemDraw"):
        super().__init__(hresult, "Exception occurred.", (0, source, text, None, 0, hresult), None)
        self.hresult, self.strerror, self.excepinfo, self.argerror = self.args


com_error = type("com_error", (ComErrorLike,), {})


def _connect(app):
    return lambda prog_ids=None, create=True: (app, "ChemDraw.Application", "fake")


@pytest.fixture
def fake_app(monkeypatch):
    app = FakeApp()
    # Replace the module's COM factory: the backend must never touch win32com/comtypes in tests.
    monkeypatch.setattr(windows, "_com_connect", _connect(app))
    return app


# A generous call budget: the fakes answer instantly, the budget only matters when a test deliberately blocks a
# call, and a loaded machine must never turn the dialog watcher's deadline into a spurious failure.
BUDGET = 60.0
POLL = 0.02


@pytest.fixture
def backend(fake_app):
    return WindowsBackend(call_timeout=BUDGET, dialog_poll=POLL)


@pytest.fixture
def replace_backend(monkeypatch):
    app = FakeApp(replace=True)
    monkeypatch.setattr(windows, "_com_connect", _connect(app))
    return WindowsBackend(call_timeout=BUDGET, dialog_poll=POLL), app


# ---------------------------------------------------------------------------- static checks
def test_module_imports_and_constructs():
    assert windows.WindowsBackend.platform == "windows"
    b = WindowsBackend()
    assert isinstance(b, ChemDrawBackend)
    assert not inspect.isabstract(WindowsBackend)
    assert b.setdata_replaces is None and b.capabilities()["com_thread"] is True


def test_every_abstract_method_is_implemented():
    missing = [m for m in ChemDrawBackend.__abstractmethods__ if m not in WindowsBackend.__dict__]
    assert not missing, "WindowsBackend does not implement: %s" % missing
    for name in ChemDrawBackend.__abstractmethods__:
        base_sig = inspect.signature(getattr(ChemDrawBackend, name))
        impl_sig = inspect.signature(getattr(WindowsBackend, name))
        assert list(base_sig.parameters) == list(impl_sig.parameters), name


def test_executable_does_not_crash(monkeypatch):
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    exe = WindowsBackend().executable()
    assert exe is None or isinstance(exe, str)


@pytest.mark.skipif(not NOT_WINDOWS, reason="ChemDraw may really be installed on Windows")
def test_find_chemdraw_windows_returns_none_off_windows(monkeypatch):
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    assert windows.find_chemdraw_windows() is None
    assert WindowsBackend().is_running() is False
    assert WindowsBackend().dismiss_dialogs() == []
    assert WindowsBackend().pending_dialogs() == []


def test_env_var_executable(tmp_path, monkeypatch):
    exe = tmp_path / "ChemDraw.exe"
    exe.write_bytes(b"MZ")
    monkeypatch.setenv("CHEMDRAW_EXE", str(exe))
    assert windows.find_chemdraw_windows() == str(exe)
    assert WindowsBackend().available() is True


def test_mime_tables_cover_insert_kinds_and_export_formats():
    assert set(windows.INSERT_MIME) == set(INSERT_KINDS)
    assert set(windows.EXPORT_MIME) == set(EXPORT_FORMATS)
    for kind in INSERT_KINDS:
        assert windows.INSERT_MIME[kind].startswith("chemical/")
    for fmt in EXPORT_FORMATS:
        assert "/" in windows.EXPORT_MIME[fmt]
    for what in ("smiles", "inchi", "inchikey", "mol", "molv3000", "cdxml", "helm", "sln", "name"):
        assert what in windows.SELECTION_FORMATS
    assert windows.SELECTION_FORMATS["inchikey"] == "chemical/x-inchikey"
    for primary, alts in windows.MIME_ALTERNATIVES.items():
        assert primary in windows.MIME.values() and alts and primary not in alts
    assert windows._mime_candidates(windows.MIME["molv3000"])[0] == "chemical/x-mdl-molfile-v3000"
    assert len(windows._mime_candidates(windows.MIME["helm"])) > 1


def test_members_table_is_well_formed():
    for key, names in windows.MEMBERS.items():
        assert "." in key and names and all(isinstance(n, str) for n in names), key
    for key, names in windows.SETTINGS_MEMBERS.items():
        assert names, key
    # selection helpers must never fall through to a destructive member
    for key in ("doc.select_all", "doc.deselect", "objs.select", "objs.selected", "objs.clean"):
        assert not {"Clear", "Delete", "Remove", "Cut"} & set(windows.MEMBERS[key]), key
    assert "Path" not in windows.MEMBERS["doc.full_name"]          # Path is the folder in Office-style models
    assert set(windows.CDXML_SETTINGS) == set(windows.SETTINGS_MEMBERS)
    assert all(v in windows.COMMANDS for v in windows.COMMAND_ALIASES.values())


def test_chemscript_bridge_is_optional():
    bridge = windows.ChemScriptBridge(module_names=["ChemScriptDoesNotExist99"])
    assert bridge.available() is False
    assert isinstance(windows.chemscript_available(), bool)
    with pytest.raises(ChemDrawNotAvailable, match="CHEMSCRIPT_PATH"):
        bridge.name_to_structure("ethanol")


def test_without_com_library_everything_raises_not_available(monkeypatch):
    if sys.platform == "win32":
        pytest.skip("pywin32 may be present")
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    b = WindowsBackend()
    with pytest.raises(ChemDrawNotAvailable):
        b.version()
    with pytest.raises(ChemDrawNotAvailable):
        b.front_document()
    with pytest.raises(ChemDrawNotAvailable):
        b.launch(wait_seconds=1)


# ---------------------------------------------------------------------------- pure helpers
def test_sequence_to_helm():
    assert windows.sequence_to_helm("ACDE", "peptide") == "PEPTIDE1{A.C.D.E}$$$$"
    assert windows.sequence_to_helm("Ala-Gly-[Aib]-Ser", "peptide") == "PEPTIDE1{A.G.[Aib].S}$$$$"
    assert windows.sequence_to_helm("ACDE FGHI", "peptide") == "PEPTIDE1{A.C.D.E.F.G.H.I}$$$$"
    assert windows.sequence_to_helm("ACGT", "dna") == "RNA1{[dR](A)P.[dR](C)P.[dR](G)P.[dR](T)}$$$$"
    assert windows.sequence_to_helm("acgu", "rna") == "RNA1{R(A)P.R(C)P.R(G)P.R(U)}$$$$"
    assert windows.sequence_to_helm(">seq\nACGTACGT", "sequence").startswith("RNA1{[dR](A)")
    assert windows.sequence_to_helm("ACGUACGU", "sequence").startswith("RNA1{R(A)")
    assert windows.sequence_to_helm("MKT", "sequence") == "PEPTIDE1{M.K.T}$$$$"
    with pytest.raises(ChemDrawError):
        windows.sequence_to_helm("", "peptide")


def test_caption_cdxml_and_summary():
    cdxml = windows.caption_cdxml("Hello <ChemDraw> & co")
    assert "&lt;ChemDraw&gt; &amp; co" in cdxml and "<t " in cdxml
    summary = windows.cdxml_summary(cdxml)
    assert summary["captions"] == 1 and summary["text"] == "Hello <ChemDraw> & co"
    summary = windows.cdxml_summary(ETHANOL_CDXML)
    assert (summary["atoms"], summary["bonds"], summary["molecules"], summary["captions"]) == (3, 2, 1, 1)
    assert summary["text"] == "ethanol"          # the "OH" atom label is not a caption
    assert summary["bounds"] == {"left": 100.0, "top": 100.0, "right": 160.0, "bottom": 165.0, "width": 60.0, "height": 65.0}
    assert windows.cdxml_summary("not xml")["atoms"] == 0


def test_ensure_extension(tmp_path):
    p = windows._ensure_extension(str(tmp_path / "out"), ".cdxml")
    assert p.endswith("out.cdxml")
    assert windows._ensure_extension(str(tmp_path / "out.jpg"), ".jpeg").endswith("out.jpeg")
    assert windows._ensure_extension(str(tmp_path / "out.txt"), ".png").endswith("out.txt.png")


def test_sniff_format(tmp_path):
    def w(name, data):
        p = tmp_path / name
        p.write_bytes(data)
        return str(p)

    assert windows._sniff_format(w("a.png", b"\x89PNG\r\n\x1a\n..."), "png") is True
    assert windows._sniff_format(w("b.png", b"VjCD0100"), "png") is False
    assert windows._sniff_format(w("c.cdx", b"VjCD0100"), "cdx") is True
    assert windows._sniff_format(w("d.cdxml", ETHANOL_CDXML.encode()), "cdxml") is True
    assert windows._sniff_format(w("e.cdxml", ETHANOL_CDXML.encode("utf-16")), "cdxml") is True
    assert windows._sniff_format(w("f.mol", b"\n  x\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"), "mol") is True
    assert windows._sniff_format(w("g.mol", b"<CDXML/>"), "mol") is False
    assert windows._sniff_format(w("h.ct", b"3 2\n"), "ct") is None          # no signature
    assert windows._sniff_format(str(tmp_path / "missing.png"), "png") is False
    assert windows._sniff_format(w("i.rxn", b"$RXN\n"), "rxn") is True and windows._sniff_format(w("j.pdf", b"%PDF-1.4"), "pdf") is True


def test_com_error_description():
    err = com_error(-2147352567, "Cannot interpret SMILES")
    assert windows._is_com_error(err)
    text = windows._describe_com_error(err)
    assert "0x80020009" in text and "Cannot interpret SMILES" in text and "ChemDraw" in text
    assert windows._member_missing(com_error(-2147352570, "Unknown name."))          # DISP_E_UNKNOWNNAME
    assert windows._member_missing(com_error(-2147352573, "Member not found."))      # DISP_E_MEMBERNOTFOUND
    assert not windows._member_missing(com_error(-2147352562, "Invalid number of parameters."))
    assert windows._member_missing(AttributeError("x"))
    assert "another thread" in windows._describe_com_error(com_error(-2147417842, "wrong thread"))


def test_merge_cdxml():
    base = render_cdxml([FakeItem.fragment("CCO"), FakeItem("caption", text="ethanol")])
    add = render_cdxml([FakeItem.fragment("c1ccccc1")])
    merged = windows.merge_cdxml(base, add)
    s = windows.cdxml_summary(merged)
    assert (s["atoms"], s["bonds"], s["molecules"], s["captions"]) == (9, 8, 2, 1)
    root = ET.fromstring(merged.split("\n", 2)[2])
    ids = [el.get("id") for el in root.iter() if el.get("id")]
    assert len(ids) == len(set(ids)), "ids must be unique after merging"
    frags = root.findall("page/fragment")
    assert len(frags) == 2 and frags[1].get("smiles") == "c1ccccc1"
    new_atoms = {n.get("id") for n in frags[1].findall("n")}
    for b in frags[1].findall("b"):
        assert b.get("B") in new_atoms and b.get("E") in new_atoms, "bond references follow the renumbering"
    # the addition sits below the base content, left aligned with it
    top = min(float(n.get("p").split()[1]) for n in frags[1].findall("n"))
    assert top >= s["bounds"]["top"] and top > 165.0
    assert float(frags[1].get("BoundingBox").split()[0]) == 100.0
    assert len(root.findall("fonttable/font")) == 1       # same font name -> mapped, not duplicated
    assert windows.merge_cdxml("", add) == add and windows.merge_cdxml(base, "") == base
    # a document without a font table / page still works and new fonts are appended
    other = add.replace('name="Arial"', 'name="Helvetica"')
    root2 = ET.fromstring(windows.merge_cdxml(base, other).split("\n", 2)[2])
    assert sorted(f.get("name") for f in root2.findall("fonttable/font")) == ["Arial", "Helvetica"]


def test_cdxml_settings_read_and_write():
    s = windows.cdxml_settings(ETHANOL_CDXML)
    assert s["bond_length"] == 14.4 and s["line_width"] == 0.6 and s["label_font"] == "Arial" and s["label_size"] == 10.0
    assert "magnification" not in s
    new = windows.cdxml_with_settings(ETHANOL_CDXML, {"bond_length": 20, "label_font": "Helvetica", "magnification": 150})
    s2 = windows.cdxml_settings(new)
    assert s2["bond_length"] == 20.0 and s2["label_font"] == "Helvetica" and s2["magnification"] == 150.0
    assert s2["caption_font"] == "Arial" and "ethanol" in new and 'name="Helvetica"' in new
    assert windows.cdxml_settings("<CDXML><page/></CDXML>") == {} and windows.cdxml_settings("junk") == {}
    with pytest.raises(ChemDrawError, match="unknown setting"):
        windows.cdxml_with_settings(ETHANOL_CDXML, {"colour": 1})
    with pytest.raises(ChemDrawError, match="not a CDXML"):
        windows.cdxml_with_settings("<svg/>", {"bond_length": 1})


def test_raw_invoke_fallback_for_parameterised_properties():
    """comtypes raises DISP_E_BADPARAMCOUNT when it reads ``Data`` as a plain property; the backend must
    then call it through IDispatch (GetIDsOfNames + Invoke) instead of giving up."""
    COMError = type("COMError", (ComErrorLike,), {})

    class RawDispatch:
        def __init__(self):
            self.calls = []

        def GetIDsOfNames(self, name):
            if name not in ("Data", "SetData", "Formula", "Count"):
                raise COMError(-2147352570, "Unknown name.")
            return [{"Data": 1, "SetData": 2, "Formula": 3, "Count": 4}[name]]

        def Invoke(self, dispid, *args, **kw):
            self.calls.append((dispid, args, kw.get("_invkind")))
            if dispid == 1:
                return "<CDXML/> for " + args[0]
            if dispid == 3:
                return "C2H6O"
            if dispid == 4:
                return 0
            return None

    class ComtypesLike:
        """``getattr`` fails like comtypes' dynamic Dispatch does for members needing arguments."""

        def __init__(self):
            self._comobj = RawDispatch()

        def __getattr__(self, name):
            if name.startswith("_"):
                raise AttributeError(name)
            if name in ("Data", "SetData"):
                raise COMError(-2147352562, "Invalid number of parameters.")
            if name == "Formula":
                raise COMError(-2147352562, "Invalid number of parameters.")
            raise COMError(-2147352570, "Unknown name.")

    obj = ComtypesLike()
    b = WindowsBackend(com_thread=False)
    assert b._invoke(obj, "objs.data", CDXML) == "<CDXML/> for " + CDXML
    assert b._get(obj, "objs.formula") == "C2H6O"
    assert b._count(obj) == 0
    with pytest.raises(ChemDrawError, match=r"MEMBERS\['objs.clean'\]"):
        b._invoke(obj, "objs.clean")
    assert obj._comobj.calls[0][1] == (CDXML,)
    with pytest.raises(AttributeError):
        windows._raw_invoke(object(), "Data", (CDXML,))


def test_press_dialog_button_matching(monkeypatch):
    monkeypatch.setattr(windows.time, "sleep", lambda s: None)
    buttons = [(11, 1, "Don't Save"), (12, 2, "Cancel"), (13, 3, "Save")]
    assert windows._press_dialog(99, ("Don\u2019t Save", "No"), buttons) == "button \"Don't Save\""
    assert windows._press_dialog(99, None, [(12, 2, "Cancel"), (14, 4, "OK")]) == "button 'OK'"
    assert windows._press_dialog(99, ("Yes",), [(12, 2, "Cancel")]) in ("IDOK", "IDCANCEL", "WM_CLOSE", "default button (id 0)", "still open")


# ---------------------------------------------------------------------------- fake COM round trips
def test_new_document_and_documents(backend, fake_app):
    assert backend.front_document() is None
    assert backend.documents() == []
    name = backend.new_document()
    assert name == "Untitled Document-1"
    assert fake_app.Documents.calls == [("Add",)]
    assert backend.front_document() == name
    assert backend.version() == "23.1.2.7"
    docs = backend.documents()
    assert len(docs) == 1 and docs[0].name == name and docs[0].front is True and docs[0].path is None
    assert docs[0].modified is False and docs[0].window_id is None
    assert docs[0].to_dict()["atoms"] == 0
    second = backend.new_document()
    assert backend.front_document() == second
    assert [d.name for d in backend.documents()] == [second, name]   # front document first
    assert backend.activate_document(name) == name
    assert fake_app.Documents.Item(1).calls[-1] == ("Activate",) and backend.front_document() == name
    assert backend.activate_document("#2") == second                 # 1-based index as a fallback id
    with pytest.raises(ChemDrawError, match="no open document named"):
        backend.activate_document("nope.cdxml")
    backend.insert(name, "CCO", kind="smiles")
    info = next(d for d in backend.documents() if d.name == name)
    assert info.modified is True and info.atoms == 3 and info.molecules == 1


def test_open_document_paths_and_lookup(backend, fake_app, tmp_path):
    src = tmp_path / "in.cdxml"
    src.write_text(ETHANOL_CDXML, encoding="utf-8")
    name = backend.open_document(str(src))
    assert name == "in.cdxml" and fake_app.Documents.calls == [("Open", str(src))]
    doc = backend.documents()[0]
    assert doc.path == str(src) and doc.front and doc.atoms == 3 and doc.captions == 1
    assert backend.activate_document("in") == "in.cdxml"            # stem
    assert backend.activate_document(str(src)) == "in.cdxml"        # full path
    assert backend.activate_document("IN.CDXML") == "in.cdxml"      # case-insensitive
    assert backend.document_cdxml(name).startswith("<?xml") and 'smiles="CCO"' in backend.document_cdxml(name)
    with pytest.raises(ChemDrawError, match="file not found"):
        backend.open_document(str(tmp_path / "missing.mol"))


def test_doc_path_from_folder_and_name():
    """Documents that only expose Office-style ``Path`` (folder) + ``Name`` still get a full path."""
    class Doc:
        Name = "aspirin.cdx"
        Path = r"C:\work"

    class Untitled:
        Name = "Untitled Document-1"
        Path = ""

    b = WindowsBackend(com_thread=False)
    assert b._doc_path(Doc()) == os.path.join(r"C:\work", "aspirin.cdx")
    assert b._doc_path(Untitled()) is None


def test_insert_smiles(backend, fake_app):
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    result = backend.insert(name, "CCO", kind="smiles")
    assert ("SetData", SMILES, "CCO") in doc.Objects.calls
    assert result["mime"] == SMILES and result["mode"] == "added"
    assert result["atoms_before"] == 0 and result["atoms_after"] == 3 and result["atoms_added"] == 3
    assert result["molecules"] == 1 and "warning" not in result
    assert result["document"] == name
    assert backend.setdata_replaces is None                 # an empty document cannot tell add from replace
    result = backend.insert(name, "c1ccccc1", kind="smiles")
    assert result["mode"] == "added" and result["atoms_after"] == 9 and result["molecules"] == 2
    assert backend.setdata_replaces is False and backend.capabilities()["setdata_replaces"] is False
    with pytest.raises(ChemDrawError, match="Cannot interpret SMILES"):
        backend.insert(name, "not-smiles", kind="smiles")
    with pytest.raises(ChemDrawError, match="nothing to insert"):
        backend.insert(name, "   ", kind="smiles")


def test_insert_cdxml_name_and_text(backend, fake_app):
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    res = backend.insert(name, ETHANOL_CDXML, kind="cdxml")
    assert res["mime"] == CDXML and res["atoms_added"] == 3 and res["captions"] == 1
    res = backend.insert(name, "benzene", kind="name")
    assert ("SetData", NAME, "benzene") in doc.Objects.calls and res["atoms_added"] == 6
    res = backend.insert(name, "a caption", kind="text")
    sets = [c for c in doc.Objects.calls if c[0] == "SetData"]
    assert sets[-1][1] == CDXML and "a caption" in sets[-1][2]
    assert res["captions"] == 2 and res["mime"] == CDXML and res["atoms_added"] == 0 and "warning" not in res
    with pytest.raises(ChemDrawError, match="unknown insert kind"):
        backend.insert(name, "x", kind="pdb")


def test_insert_sequence_fallbacks(backend, fake_app):
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    with pytest.raises(ChemDrawError, match="insert of peptide failed") as err:
        backend.insert(name, "ACDE", kind="peptide")      # fake refuses HELM and has no clipboard off Windows
    tried = [c[1] for c in doc.Objects.calls if c[0] == "SetData"]
    assert tried[0] == "chemical/x-helm" and "chemical/x-peptide-sequence" in tried
    assert "MIME_ALTERNATIVES" in str(err.value) and "clipboard" in str(err.value)
    assert ("SetData", "chemical/x-helm", "PEPTIDE1{A.C.D.E}$$$$") in doc.Objects.calls


def test_insert_when_setdata_replaces(replace_backend):
    backend, app = replace_backend
    name = backend.new_document()
    doc = app.ActiveDocument
    first = backend.insert(name, "CCO", kind="smiles")
    assert first["mode"] == "added" and first["atoms_after"] == 3
    second = backend.insert(name, "c1ccccc1", kind="smiles")
    assert second["mode"] == "replaced+merged" and "note" in second
    assert second["atoms_before"] == 3 and second["atoms_after"] == 9 and second["molecules"] == 2
    assert backend.setdata_replaces is True
    sets = [c for c in doc.Objects.calls if c[0] == "SetData"]
    assert sets[-1][1] == CDXML and "c1ccccc1" in sets[-1][2] and 'smiles="CCO"' in sets[-1][2]
    assert backend.selection_as(name, "smiles") == "CCO.c1ccccc1"
    third = backend.insert(name, "label", kind="text")
    assert third["mode"] == "replaced+merged" and third["captions"] == 1 and third["atoms_after"] == 9
    docs = backend.documents()
    assert docs[0].atoms == 9 and docs[0].molecules == 2 and docs[0].captions == 1


def test_insert_outcome_ignores_relayout():
    """Objects that only moved (ChemDraw re-laid out or re-rounded their boxes) still count as "added"."""
    before = render_cdxml([FakeItem.fragment("CCO")])
    moved = before.replace('BoundingBox="100 100 160 130"', 'BoundingBox="100.04 100 160 130.02"')
    after_added = windows.merge_cdxml(moved, render_cdxml([FakeItem.fragment("c1ccccc1")]))
    counts = lambda c: windows.cdxml_summary(c)   # noqa: E731
    outcome = WindowsBackend._insert_outcome
    assert outcome(counts(before), counts(after_added), before, after_added) == "added"
    assert outcome(counts(before), counts(moved), before, moved) == "nothing"
    replaced = render_cdxml([FakeItem.fragment("c1ccccc1"), FakeItem.fragment("CC(=O)O")])
    assert outcome(counts(before), counts(replaced), before, replaced) == "replaced"
    assert outcome(counts(""), counts(before), "", before) == "added"
    assert len(windows.cdxml_object_signatures(after_added)) == 2 and windows.cdxml_object_signatures("junk") == []


def test_insert_nothing_added_warning(backend, fake_app):
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    doc = fake_app.ActiveDocument
    doc.Objects.SetData = lambda mime, data: None            # ChemDraw silently ignores the data
    res = backend.insert(name, "c1ccccc1", kind="smiles")
    assert res["mode"] == "nothing" and res["warning"].startswith("no objects were added")
    assert res["atoms_added"] == 0


def test_selection_info_and_conversions(backend, fake_app):
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    backend.insert(name, "ethanol", kind="text")
    info = backend.selection_info(name)       # nothing selected -> whole document
    assert info["selected"] is False and info["scope"] == "document" and info["empty"] is False
    assert info["smiles"] == "CCO" and info["formula"] == "C2H6O"
    assert info["molecular_weight"] == pytest.approx(46.07) and info["exact_mass"] == pytest.approx(46.0419)
    assert info["elemental_analysis"].startswith("C, 52.14")
    assert info["atoms"] == 3 and info["bonds"] == 2 and info["molecules"] == 1 and info["captions"] == 1
    assert info["text"] == "ethanol" and info["bounds"]["left"] == 100.0 and info["bounds"]["height"] > 0
    backend.select_all(name)
    assert fake_app.ActiveDocument.calls[-1] == ("SelectAll",)
    info = backend.selection_info(name)
    assert info["selected"] is True and info["scope"] == "selection"
    assert backend.selection_cdxml(name).startswith("<?xml")
    assert backend.selection_as(name, "smiles") == "CCO"
    assert backend.selection_as(name, "cdxml") == render_cdxml(fake_app.ActiveDocument.Objects.items)
    assert backend.selection_as(name, "name") == "ethanol"
    assert "V2000" in backend.selection_as(name, "molfile")
    with pytest.raises(ChemDrawError, match="unknown conversion"):
        backend.selection_as(name, "pdb")
    with pytest.raises(ChemDrawError, match="not available"):
        backend.deselect(name)                                # the fake has no Deselect member


def test_selection_info_on_empty_document(backend, fake_app):
    name = backend.new_document()
    info = backend.selection_info(name)
    assert info["empty"] is True and info["atoms"] == 0 and "smiles" in info["unavailable"] and "formula" in info["unavailable"]
    with pytest.raises(ChemDrawError, match="could not convert"):
        backend.selection_as(name, "smiles")


def test_selection_without_count_member(backend, fake_app):
    """A Selection that only offers Data() is still recognised as non-empty from its CDXML content."""
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    doc = fake_app.ActiveDocument

    class DataOnlySelection:
        def __init__(self, items):
            self.items = items

        def Data(self, mime):
            if mime != CDXML:
                raise ValueError("only cdxml")
            return render_cdxml(self.items)

    doc.Selection = DataOnlySelection([])
    assert backend.selection_info(name)["selected"] is False
    doc.Selection = DataOnlySelection(doc.Objects.items)
    info = backend.selection_info(name)
    assert info["selected"] is True and info["atoms"] == 3 and "smiles" in info["unavailable"]


def test_selection_as_rdkit_fallback(backend, fake_app):
    pytest.importorskip("rdkit")
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    assert backend.selection_as(name, "inchikey") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert backend.selection_as(name, "inchi").startswith("InChI=1S/C2H6O")
    assert "V3000" in backend.selection_as(name, "molv3000")


def test_save_document(backend, fake_app, tmp_path):
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    doc = fake_app.ActiveDocument
    out = backend.save_document(name, str(tmp_path / "ethanol"), "cdxml")
    assert out == str(tmp_path / "ethanol.cdxml") and os.path.isfile(out)
    assert doc.calls[-1] == ("SaveAs", out)
    assert open(out, encoding="utf-8").read() == render_cdxml(doc.Objects.items)
    assert backend.documents()[0].path == out and backend.front_document() == "ethanol.cdxml"
    out2 = backend.save_document("ethanol.cdxml", str(tmp_path / "ethanol.jpg"), "jpg")
    assert out2.endswith("ethanol.jpeg") and doc.calls[-1] == ("SaveAs", out2)
    assert open(out2, "rb").read().startswith(b"\xff\xd8")
    # SaveAs ignores the extension for PDF (writes CDX): the sniff rejects it, every SaveAs shape is tried,
    # then Objects.Data("application/pdf") is written by hand
    n = len(doc.calls)
    out3 = backend.save_document(None, str(tmp_path / "ethanol.pdf"), "pdf")
    assert open(out3, "rb").read().startswith(b"%PDF") and len([c for c in doc.calls[n:] if c[0] == "SaveAs"]) == 4
    assert ("Data", "application/pdf") in doc.Objects.calls
    # SaveAs raises for SVG: Objects.Data("image/svg+xml") fallback
    out4 = backend.save_document(None, str(tmp_path / "ethanol.svg"), "svg")
    assert open(out4, encoding="utf-8").read().startswith("<svg")
    out5 = backend.save_document(None, str(tmp_path / "ethanol.ct"), "ct")     # no signature: trusted
    assert os.path.getsize(out5) > 0
    with pytest.raises(ChemDrawError, match="unknown export format"):
        backend.save_document(name, str(tmp_path / "x"), "docx")
    doc.Objects.Data = lambda mime: (_ for _ in ()).throw(ValueError("no"))
    doc.SaveAs = lambda path, *extra: None
    with pytest.raises(ChemDrawError, match="EXPORT_MIME"):
        backend.save_document(None, str(tmp_path / "nothing.png"), "png")


def test_close_document(backend, fake_app, tmp_path):
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    with pytest.raises(ChemDrawError, match="never been saved"):
        backend.close_document(name, save=True)
    backend.close_document(name, save=False)
    assert doc.calls[-1] == ("Close", False) and doc.closed and backend.documents() == []
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    backend.save_document(name, str(tmp_path / "a.cdxml"), "cdxml")
    backend.close_document("a.cdxml", save=True)
    assert ("Save",) in doc.calls and doc.calls[-1] == ("Close", True)


def test_commands(backend, fake_app):
    names = [c.name for c in backend.list_commands()]
    assert "cleanStructure" in names and "convertNameToStructure" in names and "copyAsSMILES" in names
    assert [c.name for c in backend.list_commands("copy*")] == [n for n in names if n.startswith("copyAs")]
    assert backend.list_commands("clean")[0].short_description
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    res = backend.run_command("clean-structure", name)
    assert res["ok"] and res["command"] == "cleanStructure" and res["alerts"] == []
    assert fake_app.ActiveDocument.Objects.calls[-1] == ("Clean",)
    assert backend.run_command("copyAsSMILES", name)["result"] == "CCO"
    assert backend.run_command("copyAsMol", name)["command"] == "copyAsMolfile"      # macOS alias
    assert backend.run_command("chemicalProperties", name)["result"]["formula"] == "C2H6O"
    res = backend.run_command("convertStructureToName", name)["result"]
    assert res["name"] == "ethanol" and res["caption_inserted"] is True
    assert backend.selection_info(name)["text"] == "ethanol"
    assert backend.run_command("selectAll", name)["result"] == "selected"
    assert backend.run_command("convertNameToStructure", name)["result"]["atoms_added"] == 3   # caption "ethanol" -> structure
    assert backend.run_command("dismissDialogs")["result"] == []
    with pytest.raises(ChemDrawError, match="not available through ChemDraw COM"):
        backend.run_command("flipHorizontal")
    with pytest.raises(ChemDrawError):
        backend.menus()
    with pytest.raises(ChemDrawError):
        backend.do_menu("Structure", "Clean Up Structure")


def test_document_settings_from_cdxml(backend, fake_app):
    name = backend.new_document()
    s = backend.document_settings(name)                   # no Settings object: read from the document CDXML
    assert s["source"] == "cdxml" and s["bond_length"] == 14.4 and s["label_font"] == "Arial" and s["name"] == name
    assert "magnification" in s["unavailable"]
    with pytest.raises(ChemDrawError, match="not known yet"):
        backend.set_document_settings(name, {"bond_length": 20})
    fake_app.ActiveDocument.Objects.Data = lambda mime: "<CDXML><page/></CDXML>"
    with pytest.raises(ChemDrawError, match="no Settings"):
        backend.document_settings(name)


def test_document_settings_written_through_cdxml_when_setdata_replaces(replace_backend):
    backend, app = replace_backend
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    backend.insert(name, "c1ccccc1", kind="smiles")
    assert backend.setdata_replaces is True
    s = backend.set_document_settings(name, {"bond_length": 20, "LabelFont": "Helvetica"})
    assert s["bond_length"] == 20.0 and s["label_font"] == "Helvetica" and "errors" not in s
    assert app.ActiveDocument.Objects.Count == 2           # content kept, settings changed
    s = backend.set_document_settings(name, {"colour": "red"})
    assert "colour" in s["errors"][0]


def test_document_settings_with_settings_object(backend, fake_app):
    name = backend.new_document()

    class Settings:
        BondLength = 14.4
        LineWidth = 0.6
        LabelFont = "Arial"
        LabelSize = 10

    fake_app.ActiveDocument.Settings = Settings()
    s = backend.document_settings(name)
    assert s["source"] == "com" and s["bond_length"] == 14.4 and s["label_font"] == "Arial" and "magnification" in s["unavailable"]
    s = backend.set_document_settings(name, {"bond_length": 20, "LabelSize": 12})
    assert s["bond_length"] == 20 and s["label_size"] == 12 and "errors" not in s


def test_com_error_is_converted_and_app_quit(backend, fake_app):
    name = backend.new_document()
    doc = fake_app.ActiveDocument

    def bad_set(mime, text):
        raise com_error(-2147352567, "Cannot interpret InChI")

    doc.Objects.SetData = bad_set
    with pytest.raises(ChemDrawError) as err:
        backend.insert(name, "CCO", kind="smiles")
    assert "Cannot interpret InChI" in str(err.value) and "0x80020009" in str(err.value)
    backend.quit()
    assert fake_app.quit_called and doc.calls[-1] == ("Close", False)


def test_dismiss_dialogs_with_fake_enumerator(backend, fake_app, monkeypatch):
    dialog = {"hwnd": 4242, "class": "#32770", "title": "ChemDraw",
              "text": ["ChemDraw", "No reaction center was found."], "buttons": [(11, 1, "OK"), (12, 2, "Cancel")]}
    shown = [dialog]
    pressed = []

    def find(pids):
        return list(shown)

    def press(hwnd, prefer=None, buttons=None):
        pressed.append((hwnd, tuple(prefer) if prefer else None, buttons))
        shown.clear()
        return "button 'OK'"

    monkeypatch.setattr(windows, "_find_dialogs", find)
    monkeypatch.setattr(windows, "_press_dialog", press)
    assert backend.pending_dialogs() == [{"title": "ChemDraw", "text": ["ChemDraw", "No reaction center was found."],
                                          "buttons": ["OK", "Cancel"], "hwnd": 4242}]
    texts = backend.dismiss_dialogs()
    assert texts == ["ChemDraw | No reaction center was found. [button 'OK']"]
    assert pressed == [(4242, None, dialog["buttons"])] and backend.last_dialogs == texts
    assert backend.dismiss_dialogs() == [] and backend.pending_dialogs() == []
    shown.append(dialog)
    assert backend.run_command("dismissDialogs")["result"] == texts


def test_dialog_blocking_a_com_call_is_dismissed(backend, fake_app, monkeypatch):
    """A modal alert during SetData blocks the COM call; the calling thread sees the dialog, presses
    its button (which unblocks ChemDraw) and the tool call completes with the alert text attached.

    Driven by Event handshakes: the fake dialog exists only while the fake SetData is blocked, and the
    fake SetData only returns once the button was pressed, so the ordering is fixed regardless of load.
    """
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    blocked, unblock = threading.Event(), threading.Event()
    order = []
    dialog = {"hwnd": 1, "class": "#32770", "title": "ChemDraw", "text": ["Cannot interpret the name."], "buttons": [(2, 1, "OK")]}
    original = doc.Objects.SetData

    def blocking_set(mime, data):
        order.append("blocked")
        blocked.set()                                   # the "dialog" appears now
        assert unblock.wait(30.0), "the dialog was never dismissed"
        order.append("resumed")
        return original(mime, data)

    doc.Objects.SetData = blocking_set
    monkeypatch.setattr(windows, "_find_dialogs", lambda pids: [dialog] if blocked.is_set() and not unblock.is_set() else [])

    def press(hwnd, prefer=None, buttons=None):
        order.append("pressed")
        unblock.set()                                   # pressing OK lets ChemDraw continue
        return "button 'OK'"

    monkeypatch.setattr(windows, "_press_dialog", press)
    res = backend.insert(name, "CCO", kind="smiles")
    assert order == ["blocked", "pressed", "resumed"]
    assert res["atoms_added"] == 3 and res["alerts"] == ["Cannot interpret the name. [button 'OK']"]
    assert backend.last_dialogs == res["alerts"]
    assert backend.insert(name, "c1ccccc1", kind="smiles")["alerts"] == []     # earlier dialogs are not repeated
    with pytest.raises(ChemDrawError, match="deselect is not available"):
        backend.deselect(name)


def test_call_timeout_abandons_worker_and_reconnects(fake_app):
    """A call that stays blocked inside ChemDraw past call_timeout raises; the blocked worker thread is
    abandoned (it exits once ChemDraw answers) and the next call runs on a fresh thread."""
    backend = WindowsBackend(call_timeout=0.3, dialog_poll=POLL)
    backend._worker._ensure()                            # start the worker so the stuck thread can be observed
    stuck = backend._worker._thread
    fake_app.version_gate = threading.Event()           # Version blocks until released
    with pytest.raises(ChemDrawError, match="did not answer version within"):
        backend.version()
    assert backend._worker._thread is None, "the abandoned worker must be detached from the backend"
    assert stuck.is_alive(), "the stuck thread is still inside the blocked call"
    fake_app.version_gate.set()                          # ChemDraw answers at last
    fake_app.version_gate = None
    assert backend.version() == "23.1.2.7"               # a fresh worker thread reconnected
    fresh = backend._worker._thread
    assert fresh is not None and fresh is not stuck and fresh.is_alive()
    assert backend.front_document() is None
    stuck.join(timeout=30.0)                             # the old worker finished its call and exited
    assert not stuck.is_alive()


def test_all_com_calls_run_on_one_worker_thread(backend, fake_app):
    name = backend.new_document()
    doc = fake_app.ActiveDocument
    errors = []

    def work(smiles):
        try:
            backend.insert(name, smiles, kind="smiles")
        except Exception as exc:           # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=work, args=(s,)) for s in ("CCO", "c1ccccc1", "CC(=O)O")]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    backend.selection_info(name)
    assert not errors and doc.Objects.Count == 3
    assert len(doc.Objects.threads) == 1 and threading.get_ident() not in doc.Objects.threads
    direct = WindowsBackend(com_thread=False)
    direct.insert(name, "CCO", kind="smiles")
    assert threading.get_ident() in doc.Objects.threads


def test_close_stops_worker_and_next_call_restarts_it(backend, fake_app):
    assert backend.version() == "23.1.2.7"
    worker = backend._worker
    t = worker._thread
    backend.close()
    t.join(timeout=30.0)
    assert not t.is_alive() and worker._thread is None
    assert backend.new_document() == "Untitled Document-1" and worker._thread is not None and worker._thread is not t


def test_capabilities(backend):
    caps = backend.capabilities()
    assert caps["platform"] == "windows" and caps["menus"] is False and "cleanStructure" in caps["commands"]
    assert isinstance(caps["chemscript"], bool) and caps["com_thread"] is True
    assert caps["insert_mode"] == "auto" and caps["setdata_replaces"] is None and caps["notes"]


# ---------------------------------------------------------------------------- regression: review findings
STATIONERY_TITLES = ("ACS Document 1996", "Nature Research Document", "RSC (1 Column) Document",
                     "RSC (2 Column) Document", "Wiley Document", "New Document", "New A4 Document")


@pytest.fixture
def chemdraw_items(tmp_path, monkeypatch):
    """A simulated ``<Program Files>\\ChemDraw 23\\ChemDraw Items`` folder holding ChemDraw's stationery files."""
    folder = tmp_path / "ChemDraw 23" / "ChemDraw Items"
    folder.mkdir(parents=True)
    for title in STATIONERY_TITLES:
        (folder / (title + ".cds")).write_text(ETHANOL_CDXML, encoding="utf-8")
    monkeypatch.setattr(windows, "_program_dirs", lambda: [str(tmp_path)])
    monkeypatch.setattr(windows, "STATIONERY_GLOB_SUFFIXES", [os.path.join("ChemDraw*", "ChemDraw Items")])
    monkeypatch.setenv("APPDATA", str(tmp_path / "nowhere"))
    return folder


def test_find_stationery_resolves_style_keys(chemdraw_items):
    """Style keys / aliases are mapped through styles.get_style; titles match exactly, then by prefix."""
    def found(name):
        p = windows.find_stationery(name)
        return os.path.basename(p) if p else None

    assert found("acs") == "ACS Document 1996.cds"             # the server's DEFAULT_STYLE
    assert found("jacs") == "ACS Document 1996.cds"            # alias
    assert found("nature") == "Nature Research Document.cds"
    assert found("rsc") == "RSC (2 Column) Document.cds"
    assert found("rsc_1col") == "RSC (1 Column) Document.cds"
    assert found("wiley") == "Wiley Document.cds"
    assert found("default") == "New Document.cds"
    assert found("a4") == "New A4 Document.cds"
    assert found("ACS Document 1996") == "ACS Document 1996.cds"
    assert found("ACS Document 1996.cds") == "ACS Document 1996.cds"
    assert found("acsdocument1996") == "ACS Document 1996.cds"
    assert found("Nature Research") == "Nature Research Document.cds"      # prefix, like the macOS backend
    assert found("RSC (1 Column)") == "RSC (1 Column) Document.cds"
    assert found("no such style") is None and found("") is None
    direct = str(chemdraw_items / "Wiley Document.cds")
    assert windows.find_stationery(direct) == direct


def test_new_document_with_style_key(backend, fake_app, chemdraw_items):
    """new_document('acs') (what server._temp_document / the new_document tool pass) opens the stationery."""
    name = backend.new_document("acs")
    assert name == "ACS Document 1996.cds"
    assert fake_app.Documents.calls[-1] == ("Open", str(chemdraw_items / "ACS Document 1996.cds"))
    assert backend.front_document() == name
    with pytest.raises(ChemDrawError, match="stationery 'nope' not found"):
        backend.new_document("nope")


V2000_MOL = b"\n  fake\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\nM  END\n"
V3000_MOL = b"\n  fake\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 END CTAB\nM  END\n"


def test_export_mime_and_sniff_tell_v2000_from_v3000(tmp_path):
    assert windows.EXPORT_MIME["molv3000"] == windows.MIME["molv3000"] != windows.EXPORT_MIME["mol"]
    assert windows.EXPORT_MIME["sdfv3000"] == windows.MIME["sdfv3000"] != windows.EXPORT_MIME["sdf"]
    assert windows.EXPORT_MIME["rxnv3000"] == windows.MIME["rxnv3000"] != windows.EXPORT_MIME["rxn"]
    assert windows._mime_candidates(windows.EXPORT_MIME["sdfv3000"])[0] == "chemical/x-mdl-sdfile-v3000"

    def w(name, data):
        p = tmp_path / name
        p.write_bytes(data)
        return str(p)

    v2, v3 = w("v2.mol", V2000_MOL), w("v3.mol", V3000_MOL)
    assert windows._sniff_format(v2, "mol") is True and windows._sniff_format(v3, "mol") is False
    assert windows._sniff_format(v3, "molv3000") is True and windows._sniff_format(v2, "molv3000") is False
    assert windows._sniff_format(w("v2.sdf", V2000_MOL + b"$$$$\n"), "sdf") is True
    assert windows._sniff_format(w("v3.sdf", V3000_MOL + b"$$$$\n"), "sdf") is False
    assert windows._sniff_format(w("v3b.sdf", V3000_MOL + b"$$$$\n"), "sdfv3000") is True
    assert windows._sniff_format(w("v2b.sdf", V2000_MOL + b"$$$$\n"), "sdfv3000") is False
    r2, r2crlf, r3 = (w("r2.rxn", b"$RXN\n\n  fake\n\n  1  1\n$MOL\n" + V2000_MOL),
                      w("r2crlf.rxn", b"$RXN\r\n\r\n  fake\r\n\r\n  1  1\r\n$MOL\r\n"),
                      w("r3.rxn", b"$RXN V3000\n\n  fake\n\nM  V30 COUNTS 1 1\n"))
    assert windows._sniff_format(r2, "rxn") is True and windows._sniff_format(r2crlf, "rxn") is True
    assert windows._sniff_format(r3, "rxn") is False
    assert windows._sniff_format(r3, "rxnv3000") is True and windows._sniff_format(r2, "rxnv3000") is False


def test_save_document_retries_when_the_other_molfile_variant_was_written(backend, fake_app, tmp_path):
    """ChemDraw 23's SaveAs(x.mol) writes V3000: save_document(fmt='mol') must not return it as the V2000 export
    (and fmt='molv3000' must not return a V2000 file) but fall through to Objects.Data with the right MIME."""
    name = backend.new_document()
    backend.insert(name, "CCO", kind="smiles")
    doc = fake_app.ActiveDocument
    orig_save_as, orig_data = doc.SaveAs, doc.Objects.Data

    def save_as_v3000(path, *extra):                 # the plain "MDL Molfile" type writes V3000
        doc.calls.append(("SaveAs", path) + extra)
        with open(path, "wb") as fh:
            fh.write(V3000_MOL)

    def data(mime):
        doc.Objects.calls.append(("Data", mime))
        if mime == windows.MIME["molv3000"]:
            return V3000_MOL.decode()
        return orig_data(mime)

    doc.SaveAs, doc.Objects.Data = save_as_v3000, data
    n = len(doc.calls)
    out = backend.save_document(name, str(tmp_path / "ethanol.mol"), "mol")
    text = open(out, "rb").read()
    assert b"V2000" in text and b"V3000" not in text
    assert len([c for c in doc.calls[n:] if c[0] == "SaveAs"]) == 4            # every SaveAs shape was tried
    assert ("Data", MOL) in doc.Objects.calls
    assert ("SaveAs", out, "MDL Molfile V2000") in doc.calls[n:]

    doc.SaveAs = orig_save_as                        # the fake's default SaveAs(x.mol) writes V2000 (via Data(MOL))
    n = len(doc.calls)
    out3 = backend.save_document(name, str(tmp_path / "ethanol_v3.mol"), "molv3000")
    text3 = open(out3, "rb").read()
    assert b"V3000" in text3 and b"V2000" not in text3
    assert len([c for c in doc.calls[n:] if c[0] == "SaveAs"]) == 4
    assert ("SaveAs", out3, "chemical/x-mdl-molfile-v3000") in doc.calls[n:]      # EXPORT_MIME tells them apart
    assert doc.Objects.calls[-1] == ("Data", "chemical/x-mdl-molfile-v3000")      # the fallback asked for V3000


def test_launch_fails_fast_without_com_library(monkeypatch, tmp_path):
    """Without pywin32/comtypes launch() must raise ChemDrawNotAvailable at once instead of starting ChemDraw.exe
    and retrying the same permanent error until wait_seconds is over (90 s by default)."""
    if sys.platform == "win32":
        pytest.skip("pywin32 may be present")
    started = []
    monkeypatch.setattr(windows.subprocess, "Popen", lambda *a, **k: started.append(a) or pytest.fail("Popen called"))
    exe = tmp_path / "ChemDraw.exe"
    exe.write_text("stub")
    b = WindowsBackend(executable_path=str(exe), dialog_poll=POLL)
    t0 = time.time()
    with pytest.raises(ChemDrawNotAvailable, match="pywin32"):
        b.launch(wait_seconds=10)
    assert time.time() - t0 < 5.0 and started == []


class FakeComLib:
    """A COM client library whose GetActiveObject never finds ChemDraw (no ROT entry while it starts up)."""

    def __init__(self):
        self.apps = []

    def GetActiveObject(self, prog_id):
        raise OSError(-2147221021, "Operation unavailable")

    def CreateObject(self, prog_id):
        self.apps.append(FakeApp())
        return self.apps[-1]

    Dispatch = CreateObject


def test_launch_reuses_the_proxy_create_object_returned(monkeypatch):
    """After launch() started ChemDraw through CreateObject it must talk to that proxy, not call CreateObject
    again (a second instance, and a bogus 'may have started a second ChemDraw' note once the pid shows up)."""
    lib = FakeComLib()
    monkeypatch.setattr(windows, "_com_library", lambda: ("fake", lib))
    monkeypatch.setattr(windows, "_chemdraw_pids", lambda max_age=2.0: [4242] if lib.apps else [])
    b = WindowsBackend(call_timeout=BUDGET, dialog_poll=POLL)
    info = b.launch(wait_seconds=5)
    assert len(lib.apps) == 1, "CreateObject must run once"
    assert info["started_via"] == "COM CreateObject(%s)" % windows.PROG_IDS[0] and info["version"] == "23.1.2.7"
    assert info["pid"] == 4242 and "connect_note" not in info and b.connect_note is None
    assert lib.apps[0].Visible is True
    assert b.version() == "23.1.2.7" and len(lib.apps) == 1     # later calls keep using the same proxy
