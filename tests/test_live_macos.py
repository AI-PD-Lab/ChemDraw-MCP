"""Live tests that drive the real ChemDraw application on macOS (``MacBackend`` and the server tools).

They are skipped unless all of the following hold:

* ``sys.platform == "darwin"``;
* ``CHEMDRAWMCP_LIVE=1`` is set in the environment;
* :func:`chemdrawmcp.discovery.find_chemdraw` finds ChemDraw.app (``CHEMDRAW_APP`` overrides the search).

Every test carries the ``live`` marker and the project's pytest ``addopts`` deselects that marker by default, so
run them with ``CHEMDRAWMCP_LIVE=1 python -m pytest tests/test_live_macos.py -m live -v``.

One ``MacBackend(timeout=90)`` is shared by the whole session (``live`` fixture).  It launches ChemDraw if
needed, records which documents were open before the session and, at teardown, closes *without saving* every
document that was not open before (the tests' own documents plus anything a failing test left behind).
Documents are addressed by *window id* because ChemDraw renames an untitled document after its first edit
("Untitled ACS Document 1996-1" becomes "ChemDraw Document 1"), so the id is resolved right after creation.
Each test closes the documents it created in a ``finally`` block; files go to ``tmp_path``.

The expected values were verified on ChemDraw 23.1.2.  Assertions stay tolerant (``in``, ``approx``, canonical
SMILES comparison) where another ChemDraw version might differ in detail.
"""
from __future__ import annotations

import io
import json
import os
import re
import sys
from typing import Optional, Set, Tuple

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

PILImage = pytest.importorskip("PIL.Image")

from chemdrawmcp import cdxml, chem, layout, styles  # noqa: E402
from chemdrawmcp.backends import ChemDrawError, CommandInfo, DocumentInfo  # noqa: E402


# ---------------------------------------------------------------------------
# skip unless a live ChemDraw is wanted and available
# ---------------------------------------------------------------------------

def _live_reason() -> Optional[str]:
    """Why the live tests cannot run here (None = they can).  Never launches or talks to ChemDraw."""
    if sys.platform != "darwin":
        return "live ChemDraw tests only run on macOS"
    if os.environ.get("CHEMDRAWMCP_LIVE") != "1":
        return "set CHEMDRAWMCP_LIVE=1 to run the live ChemDraw tests"
    try:
        from chemdrawmcp import discovery
        if not discovery.find_chemdraw():
            return "ChemDraw.app not found (set CHEMDRAW_APP to its path)"
    except Exception as e:  # pragma: no cover - discovery is defensive already
        return "ChemDraw discovery failed: %s" % e
    return None


LIVE_REASON = _live_reason()
pytestmark = [pytest.mark.live, pytest.mark.skipif(LIVE_REASON is not None, reason=LIVE_REASON or "")]


# ---------------------------------------------------------------------------
# constants & helpers
# ---------------------------------------------------------------------------

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
PHENOL = "c1ccccc1O"
PHENOL_INCHIKEY = "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
ESTERIFICATION = "CC(=O)O + OCC --[H2SO4]--> CC(=O)OCC + O"
ESTERIFICATION_RXN = "CC(=O)O.CCO>>CCOC(C)=O.O"
SUBSCRIPT_DIGITS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def canon(smiles: str) -> str:
    """Canonical SMILES (atom-map numbers dropped) so ChemDraw's and RDKit's spellings can be compared."""
    mol = Chem.MolFromSmiles(smiles.strip())
    assert mol is not None, "ChemDraw returned something RDKit cannot read as SMILES: %r" % smiles
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def canon_side(side: str):
    return sorted(canon(s) for s in side.strip().split(".") if s)


def canon_rxn(reaction_smiles: str):
    """(reactants, products) of a reaction SMILES, canonicalised and sorted; agents are ignored."""
    parts = reaction_smiles.strip().split(">")
    assert len(parts) == 3, "not a reaction SMILES: %r" % reaction_smiles
    return canon_side(parts[0]), canon_side(parts[2])


def formula_text(formula) -> str:
    """ChemDraw's formula without whitespace or subscript digits ("C6H6O" however it was typeset)."""
    return re.sub(r"\s+", "", str(formula)).translate(SUBSCRIPT_DIGITS)


def read(path: str, binary: bool = False):
    with open(path, "rb" if binary else "r", encoding=None if binary else "utf-8", errors=None if binary else "replace") as fh:
        return fh.read()


def is_under(path: str, root) -> bool:
    path, root = os.path.realpath(path), os.path.realpath(str(root))
    return os.path.commonpath([path, root]) == root


def assert_png(path: str, root=None) -> None:
    assert path.lower().endswith(".png"), path
    assert os.path.isfile(path), path
    assert read(path, binary=True)[:8] == PNG_MAGIC, path
    if root is not None:
        assert is_under(path, root), (path, root)


def loads(result) -> dict:
    """A tool's JSON string -> dict (asserting it is a plain string, not ``[text, Image]``)."""
    assert isinstance(result, str), "expected a JSON string, got %s" % type(result).__name__
    return json.loads(result)


def reaction_cdxml(path) -> str:
    """Write the esterification (with reagent/condition text) as CDXML via layout.build_reaction; returns the path."""
    def comp(smi):
        return layout.Component(chem.mol_from_input(smi, "smiles"))
    step = layout.Step([comp("CC(=O)O"), comp("OCC")], [comp("CC(=O)OCC"), comp("O")], above=["H2SO4"], below=["reflux"])
    doc = layout.build_reaction(step, style="acs", name="esterification")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc.to_xml())
    return str(path)


class LiveSession:
    """The shared ChemDraw session: the backend, what was open before the tests, what the tests created."""

    def __init__(self, backend, launch_info: dict):
        self.backend = backend
        self.launch_info = launch_info
        self.before: Set[int] = self.window_ids()
        self.created: Set[int] = set()

    def window_ids(self) -> Set[int]:
        return {d.window_id for d in self.backend.documents() if d.window_id is not None}

    def find(self, wid: int) -> Optional[DocumentInfo]:
        for d in self.backend.documents():
            if d.window_id == wid:
                return d
        return None

    def new_doc(self, style: str = "acs") -> Tuple[int, str]:
        """Create a document and resolve its window id at once (the name changes after the first edit)."""
        known = self.window_ids()
        name = self.backend.new_document(style)
        wid = None
        for d in self.backend.documents():
            if d.name == name and d.window_id is not None:
                wid = d.window_id
                break
        if wid is None:  # name lookup failed: fall back to "the one window that is new"
            fresh = self.window_ids() - known
            if len(fresh) == 1:
                wid = fresh.pop()
        if wid is None:
            pytest.fail("could not resolve the window id of the new document %r" % name)
        self.created.add(wid)
        return wid, name

    def close(self, wid) -> None:
        """Close a document without saving; never raises (meant for ``finally`` blocks)."""
        try:
            self.backend.close_document(wid, save=False)
        except ChemDrawError:
            pass
        self.created.discard(wid)

    def cleanup(self):
        """Dismiss stray alerts and close (without saving) everything that was not open before the session."""
        closed = []
        try:
            self.backend.dismiss_dialogs()
            docs = self.backend.documents()
        except ChemDrawError:
            return closed
        for d in docs:
            if d.window_id is None:
                continue
            if d.window_id in self.created or d.window_id not in self.before:
                try:
                    self.backend.close_document(d.window_id, save=False)
                    closed.append(d.name)
                except ChemDrawError:
                    pass
        self.created.clear()
        return closed


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def live():
    """One MacBackend for the whole session; ChemDraw launched if needed; test documents closed at the end."""
    from chemdrawmcp.backends.macos import MacBackend
    backend = MacBackend(timeout=90)
    info = backend.launch(wait_seconds=120)
    session = LiveSession(backend, info)
    yield session
    session.cleanup()


@pytest.fixture
def cd(live):
    return live.backend


@pytest.fixture
def out_dir(tmp_path):
    return tmp_path / "renders"


@pytest.fixture
def server(live, out_dir, monkeypatch):
    """chemdrawmcp.server wired to the session's live backend, writing renders into ``out_dir``."""
    monkeypatch.setenv("CHEMDRAWMCP_OUTPUT_DIR", str(out_dir))
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "0")
    monkeypatch.delenv("CHEMDRAWMCP_BACKEND", raising=False)
    monkeypatch.delenv("CHEMDRAWMCP_AUTOLAUNCH", raising=False)
    import chemdrawmcp.server as srv
    srv._BACKEND = live.backend
    srv._NAME_CACHE.clear()
    yield srv
    srv._BACKEND = None
    srv._NAME_CACHE.clear()


# ---------------------------------------------------------------------------
# application lifecycle
# ---------------------------------------------------------------------------

def test_launch_reports_running_and_version(live, cd):
    info = live.launch_info
    assert info["running"] is True
    assert isinstance(info["was_running"], bool)
    assert re.match(r"^\d+\.\d+", str(info["version"])), info
    assert info["app"] and os.path.exists(info["app"])
    v = cd.version()
    assert isinstance(v, str) and re.match(r"^\d+\.\d+", v), v
    assert v == str(info["version"])
    assert cd.is_running() is True
    assert cd.ping() is True
    # launching again is a no-op
    again = cd.launch(wait_seconds=30)
    assert again["running"] is True and again["was_running"] is True


def test_capabilities(cd):
    caps = cd.capabilities()
    assert caps["platform"] == "macos"
    assert caps["available"] is True and caps["running"] is True
    assert caps["executable"] == cd.app_path
    assert "chemdraw" in str(caps["bundle_id"]).lower()


def test_js_passthrough(cd):
    assert cd.js("return ARGS.a + ARGS.b;", a=40, b=2) == 42
    assert "chemdraw" in str(cd.js("return app.name();")).lower()
    with pytest.raises(ChemDrawError, match="boom"):
        cd.js("throw new Error('boom');")


def test_dismiss_and_pending_dialogs_when_idle(cd):
    assert cd.dismiss_dialogs() == []
    assert cd.last_alerts == []
    assert cd.pending_dialogs() == []


def test_argument_validation_happens_before_talking_to_chemdraw(cd):
    with pytest.raises(ChemDrawError, match="kind must be one of"):
        cd.insert(None, "CCO", "bogus")
    with pytest.raises(ChemDrawError, match="what must be one of"):
        cd.selection_as(None, "bogus")
    with pytest.raises(ChemDrawError, match="unknown export format"):
        cd.save_document(None, "/nonexistent/x.docx", "docx")
    with pytest.raises(ChemDrawError, match="file not found"):
        cd.open_document("/nonexistent/no_such_file.cdxml")
    with pytest.raises(ChemDrawError, match="stationery"):
        cd.new_document("no_such_style_xyz")


# ---------------------------------------------------------------------------
# documents
# ---------------------------------------------------------------------------

def test_documents_lists_new_document_with_window_id(live, cd):
    wid, name = live.new_doc("acs")
    try:
        docs = cd.documents()
        assert docs and all(isinstance(d, DocumentInfo) for d in docs)
        assert all(d.window_id is None or isinstance(d.window_id, int) for d in docs), [d.to_dict() for d in docs]
        mine = [d for d in docs if d.window_id == wid]
        assert len(mine) == 1, [d.to_dict() for d in docs]
        d = mine[0]
        assert isinstance(d.window_id, int)
        assert d.name == name
        assert not d.path                              # untitled
        assert d.atoms == 0 and d.bonds == 0 and d.molecules == 0 and d.captions == 0 and d.arrows == 0
        # the new document is the (only) front document and is listed first
        assert [x.window_id for x in docs if x.front] == [wid]
        assert docs[0].window_id == wid
        assert cd.front_document() == name
        assert cd.activate_document(wid) == name
        as_dict = d.to_dict()
        assert as_dict["window_id"] == wid and as_dict["front"] is True and "extra" not in as_dict
    finally:
        live.close(wid)
    assert wid not in live.window_ids()


def test_new_document_name_and_window_id_survives_rename(live, cd):
    wid, name = live.new_doc("acs")
    try:
        assert name.startswith("Untitled"), name
        assert isinstance(wid, int)
        r = cd.insert(wid, "CCO", "smiles")
        assert r["atoms_added"] == 3
        # ChemDraw may rename the untitled document after its first edit ("ChemDraw Document N");
        # the window id keeps addressing it
        d = live.find(wid)
        assert d is not None, "document vanished from documents() after its first edit"
        assert d.atoms == 3 and d.molecules == 1
        assert r["document"] in (name, d.name)
        assert cd.document_settings(wid)["name"] == d.name
        assert cd.activate_document(wid) == d.name
    finally:
        live.close(wid)


def test_close_document_by_window_id_then_unknown_id_raises(live, cd):
    wid, _ = live.new_doc("acs")
    cd.close_document(wid, save=False)
    live.created.discard(wid)
    assert wid not in live.window_ids()
    with pytest.raises(ChemDrawError, match="no document window with id"):
        cd.close_document(wid, save=False)


# ---------------------------------------------------------------------------
# content: insert, selection
# ---------------------------------------------------------------------------

def test_insert_smiles_and_selection_info(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        r = cd.insert(wid, PHENOL, "smiles")
        assert r["kind"] == "smiles"
        assert r["atoms_added"] == 7, r
        assert r["after"]["molecules"] == r["before"]["molecules"] + 1
        assert "warning" not in r, r
        if r.get("selection") and r["selection"].get("smiles"):   # the pasted structure is left selected
            assert canon(r["selection"]["smiles"]) == canon("Oc1ccccc1")
        info = cd.selection_info(wid)
        assert info["empty"] is False
        assert info["atoms"] == 7 and info["bonds"] == 7 and info["molecules"] == 1
        assert canon(info["smiles"]) == canon("Oc1ccccc1")
        assert formula_text(info["formula"]) == "C6H6O", info["formula"]
        assert float(info["molecular_weight"]) == pytest.approx(94.11, abs=0.02)
        if info.get("exact_mass") is not None:
            assert float(info["exact_mass"]) == pytest.approx(94.042, abs=0.01)
        assert "C, 76.57" in str(info["elemental_analysis"]), info["elemental_analysis"]
        assert info["bounds"] is not None
        assert isinstance(info["document"], str)
    finally:
        live.close(wid)


def test_insert_name_caffeine(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        r = cd.insert(wid, "caffeine", "name")
        assert r["kind"] == "name"
        assert r["atoms_added"] == 14, r
        info = cd.selection_info(wid)
        assert info["atoms"] == 14 and info["molecules"] == 1
        assert canon(info["smiles"]) == canon(CAFFEINE)
        assert formula_text(info["formula"]) == "C8H10N4O2", info["formula"]
    finally:
        live.close(wid)


def test_insert_text_adds_a_caption(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        r = cd.insert(wid, "Scheme 1", "text")
        assert r["kind"] == "text" and r["atoms_added"] == 0
        assert "warning" not in r, r
        assert r["after"]["captions"] >= r["before"]["captions"] + 1, r
        assert "Scheme 1" in cdxml.document_text(cd.document_cdxml(wid))
    finally:
        live.close(wid)


def test_select_all_and_deselect(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        assert cd.insert(wid, PHENOL, "smiles")["atoms_added"] == 7
        assert cd.insert(wid, "CCO", "smiles")["atoms_added"] == 3
        cd.deselect(wid)
        assert cd.selection_info(wid)["empty"] is True
        cd.select_all(wid)
        info = cd.selection_info(wid)
        assert info["empty"] is False
        assert info["atoms"] == 10 and info["molecules"] == 2, info
        assert canon_side(info["smiles"]) == canon_side("Oc1ccccc1.CCO"), info["smiles"]
    finally:
        live.close(wid)


def test_selection_as_formats(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        assert cd.selection_as(wid, "inchikey").strip() == PHENOL_INCHIKEY
        inchi = cd.selection_as(wid, "inchi").strip()
        assert inchi.startswith("InChI=1S/C6H6O/"), inchi
        assert canon(cd.selection_as(wid, "smiles")) == canon("Oc1ccccc1")
        mol = cd.selection_as(wid, "mol")
        assert "V2000" in mol and "M  END" in mol
        assert "V3000" in cd.selection_as(wid, "molv3000")
        xml = cd.selection_as(wid, "cdxml")
        assert "<fragment" in xml
        back = cdxml.read_cdxml(xml)
        assert len(back["molecules"]) == 1, back["molecules"]
        assert canon(back["molecules"][0]["smiles"]) == canon("Oc1ccccc1")
        assert "<fragment" in cd.selection_cdxml(wid)
    finally:
        live.close(wid)


# ---------------------------------------------------------------------------
# save / export
# ---------------------------------------------------------------------------

def _check_png(path: str) -> None:
    im = PILImage.open(path)
    assert im.format == "PNG"
    assert im.size[0] > 50 and im.size[1] > 50, im.size
    dpi = im.info.get("dpi")
    assert dpi, "PNG carries no dpi information"
    assert 280 <= float(dpi[0]) <= 320, dpi  # ChemDraw writes ~300 dpi (observed 301.9 x 308.0)


def _check_svg(path: str) -> None:
    text = read(path).lstrip("﻿ \t\r\n")
    assert text.startswith("<?xml") or text.startswith("<svg"), text[:120]
    assert "<svg" in text


def _check_cdxml(path: str) -> None:
    text = read(path)
    assert "<CDXML" in text
    back = cdxml.read_cdxml(text)
    assert len(back["molecules"]) == 1, back["molecules"]
    assert canon(back["molecules"][0]["smiles"]) == canon("Oc1ccccc1")
    assert back["settings"].get("BondLength") == pytest.approx(14.4, abs=0.05), back["settings"]


def _check_pdf(path: str) -> None:
    assert read(path, binary=True)[:4] == b"%PDF"


def _check_mol(path: str) -> None:
    text = read(path)
    assert "V2000" in text and "M  END" in text


FORMAT_CHECKS = {"png": _check_png, "svg": _check_svg, "cdxml": _check_cdxml, "pdf": _check_pdf, "mol": _check_mol}


@pytest.mark.parametrize("fmt", sorted(FORMAT_CHECKS))
def test_save_document_formats(live, cd, tmp_path, fmt):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        target = tmp_path / ("phenol." + fmt)
        written = cd.save_document(wid, str(target), fmt)
        assert os.path.isfile(written) and os.path.getsize(written) > 0, written
        assert os.path.realpath(os.path.dirname(written)) == os.path.realpath(str(tmp_path))
        assert os.path.basename(written).startswith("phenol"), written
        FORMAT_CHECKS[fmt](written)
    finally:
        live.close(wid)


def test_save_document_appends_extension_and_overwrites(live, cd, tmp_path):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        written = cd.save_document(wid, str(tmp_path / "noext"), "png")
        assert written.endswith(".png"), written
        assert_png(written, tmp_path)
        # jpg alias -> jpeg writer; the given extension is kept
        jpg = cd.save_document(wid, str(tmp_path / "phenol.jpg"), "jpg")
        assert jpg.endswith(".jpg") and os.path.getsize(jpg) > 0
        assert PILImage.open(jpg).format == "JPEG"
        # saving again to the same path replaces the file
        again = cd.save_document(wid, written, "png")
        assert again == written and os.path.getsize(again) > 0
    finally:
        live.close(wid)


def test_save_rxn_without_reaction_arrow_raises(live, cd, tmp_path):
    """ChemDraw refuses to write RXN for a document without a reaction ("No reaction center was found"): the
    alert blocks the Apple event, the backend dismisses it and reports its text.  This can take up to the
    backend's save timeout (120 s) when the alert blocks until the timeout fires."""
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        target = tmp_path / "phenol.rxn"
        with pytest.raises(ChemDrawError) as excinfo:
            cd.save_document(wid, str(target), "rxn")
        assert "reaction" in str(excinfo.value).lower(), str(excinfo.value)
        leftover = cd.dismiss_dialogs()
        assert leftover == [], "the backend should already have dismissed the alert: %r" % leftover
        assert cd.pending_dialogs() == []
    finally:
        cd.dismiss_dialogs()
        live.close(wid)


def test_document_cdxml_roundtrip(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        xml = cd.document_cdxml(wid)
        assert xml.startswith("<?xml") and "<CDXML" in xml and "\r\n" not in xml
        back = cdxml.read_cdxml(xml)
        assert [canon(m["smiles"]) for m in back["molecules"]] == [canon("Oc1ccccc1")]
        assert back["pages"] >= 1
    finally:
        live.close(wid)


# ---------------------------------------------------------------------------
# document settings
# ---------------------------------------------------------------------------

def test_document_settings_acs_stationery_and_change(live, cd):
    wid, name = live.new_doc("acs")
    try:
        s = cd.document_settings(wid)
        assert s["name"] == name
        acs = styles.get_style("acs")
        assert s["bond_length"] == pytest.approx(14.4, abs=0.01), s
        for key in ("label_size", "caption_size"):
            assert s.get(key) == pytest.approx(acs[key], abs=0.05), (key, s.get(key), acs[key])
        for key in ("line_width", "bold_width", "margin_width", "hash_spacing", "bond_spacing", "chain_angle",
                    "label_font", "caption_font"):
            assert key in s, "setting %s missing from %r" % (key, s)
        new = cd.set_document_settings(wid, {"bond_length": 20})
        assert new["bond_length"] == pytest.approx(20, abs=0.01), new
        assert cd.document_settings(wid)["bond_length"] == pytest.approx(20, abs=0.01)
        assert cd.set_document_settings(wid, {"bond_length": 14.4})["bond_length"] == pytest.approx(14.4, abs=0.01)
        with pytest.raises(ChemDrawError, match="unknown setting"):
            cd.set_document_settings(wid, {"no_such_setting": 1})
    finally:
        live.close(wid)


# ---------------------------------------------------------------------------
# commands & menus
# ---------------------------------------------------------------------------

def test_run_command_clean_structure_and_struct_to_name(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        cd.select_all(wid)
        r = cd.run_command("cleanStructure", wid)
        assert r["ok"] is True and r["enabled_before"] is True and r["command"] == "cleanStructure", r
        assert r["after"]["atoms"] == 7 and r["alerts"] == []
        assert isinstance(r.get("description"), str)

        cd.select_all(wid)
        r2 = cd.run_command("convertStructureToName", wid)
        assert r2["ok"] is True, r2
        assert r2["after"]["captions"] > r2["before"]["captions"], r2
        texts = cdxml.document_text(cd.document_cdxml(wid))
        assert any("phenol" in t.lower() for t in texts), texts
        d = live.find(wid)
        assert d is not None and d.captions >= 1

        with pytest.raises(ChemDrawError, match="no ChemDraw command named"):
            cd.run_command("noSuchCommandXyz", wid)
    finally:
        live.close(wid)


def test_list_commands_and_menus(cd):
    rows = cd.list_commands("^copyAs")
    assert rows and all(isinstance(r, CommandInfo) for r in rows)
    names = {r.name for r in rows}
    assert all(n.startswith("copyAs") for n in names), names
    assert {"copyAsSMILES", "copyAsCdxml", "copyAsMol", "copyAsInChIKey"} <= names, names
    smi = next(r for r in rows if r.name == "copyAsSMILES")
    assert smi.short_description or smi.long_description, smi  # some commands only carry a long description
    assert isinstance(smi.enabled, bool)
    assert smi.to_dict()["name"] == "copyAsSMILES"

    names_only = cd.list_commands("^copyAs", with_descriptions=False)
    assert {r.name for r in names_only} == names
    assert all(r.enabled is None and r.short_description == "" for r in names_only)

    everything = cd.list_commands(None, with_descriptions=False)
    assert len(everything) > 500, len(everything)          # ~1600 on ChemDraw 23
    assert {"selectAll", "cleanStructure", "convertStructureToName", "pasteAsSMILES"} <= {r.name for r in everything}

    menus = cd.menus()
    assert isinstance(menus, dict) and "Structure" in menus and "File" in menus and "Edit" in menus, sorted(menus)
    assert any(item.startswith("Clean Up Structure") for item in menus["Structure"]), menus["Structure"]


# ---------------------------------------------------------------------------
# opening files: a reaction scheme
# ---------------------------------------------------------------------------

def test_open_reaction_cdxml_document(live, cd, tmp_path):
    path = reaction_cdxml(tmp_path / "esterification.cdxml")
    before = live.window_ids()
    name = cd.open_document(path)
    fresh = live.window_ids() - before
    assert len(fresh) == 1, "expected exactly one new window, got %r" % fresh
    wid = fresh.pop()
    live.created.add(wid)
    try:
        d = live.find(wid)
        assert d is not None and d.name == name
        assert name.startswith("esterification"), name
        assert d.path and os.path.basename(d.path) == "esterification.cdxml", d.path
        assert d.arrows == 1, d.to_dict()
        assert d.molecules >= 4 and d.atoms >= 14, d.to_dict()   # AcOH (4) + EtOH (3) + ester (6) + water (1)
        assert d.captions >= 2                                      # "H2SO4", "reflux" (and the plus signs)
        cd.select_all(wid)
        info = cd.selection_info(wid)
        smiles = str(info["smiles"]).strip()
        assert ">" in smiles, smiles       # ChemDraw interprets the scheme as a reaction -> reaction SMILES
        assert smiles.count(">") == 2, smiles
        assert canon_rxn(smiles) == canon_rxn(ESTERIFICATION_RXN), smiles
        assert info["arrows"] == 1 and info["molecules"] >= 4, info
        back = cdxml.read_cdxml(cd.document_cdxml(wid))
        assert len(back["arrows"]) == 1
        assert len(back["steps"]) == 1, back["steps"]
        assert canon_rxn(back["reaction_smiles"][0]) == canon_rxn(ESTERIFICATION_RXN)
    finally:
        live.close(wid)
    assert wid not in live.window_ids()


# ---------------------------------------------------------------------------
# screenshot
# ---------------------------------------------------------------------------

def test_screenshot_returns_png(live, cd):
    wid, _ = live.new_doc("acs")
    try:
        cd.insert(wid, PHENOL, "smiles")
        try:
            png = cd.screenshot(wid)
        except ChemDrawError as e:
            pytest.skip("screenshot not possible here (Screen Recording permission?): %s" % e)
        assert isinstance(png, (bytes, bytearray)) and png[:8] == PNG_MAGIC
        im = PILImage.open(io.BytesIO(png))
        assert im.size[0] > 100 and im.size[1] > 100, im.size
    finally:
        live.close(wid)


# ---------------------------------------------------------------------------
# server-level tools against the live application
# ---------------------------------------------------------------------------

def test_server_status_live(server, live):
    st = loads(server.status())
    assert st["backend"] == "macos" and st["chemdraw_available"] is True and st["running"] is True, st
    assert re.match(r"^\d+\.\d+", str(st.get("version"))), st.get("version")
    assert isinstance(st.get("documents"), list)
    assert "pending_dialogs" not in st, st.get("pending_dialogs")
    assert os.path.isdir(st["renders_dir"])


def test_server_draw_reaction_read_document_and_close(server, live, out_dir):
    before = live.window_ids()
    result = server.draw_reaction(reaction=ESTERIFICATION, return_image=False)
    info = loads(result)                       # return_image=False -> plain JSON, no [text, Image]
    assert info.get("ok", True) is not False, info
    assert info.get("opened") is True, info
    wid = info["window_id"]
    live.created.add(wid)
    try:
        assert isinstance(wid, int) and wid not in before
        assert info.get("arrows") == 1, info
        assert is_under(info["cdxml_path"], out_dir)
        assert canon_rxn(info["reaction_smiles"]) == canon_rxn(ESTERIFICATION_RXN)
        render = info["render"]
        assert render["rendered_by"] == "ChemDraw" and render["format"] == "png", render
        assert_png(render["path"], out_dir)
        assert render.get("size_px") and render["size_px"][0] > 50, render
        assert "analysis" in render, render

        rd = loads(server.read_document(document=str(wid)))
        assert "error" not in rd, rd
        assert len(rd["steps"]) == 1, rd["steps"]
        assert len(rd["arrows"]) == 1, rd["arrows"]
        assert len(rd["molecules"]) >= 4, rd["molecules"]
        assert canon_rxn(rd["reaction_smiles"][0]) == canon_rxn(ESTERIFICATION_RXN), rd["reaction_smiles"]
        assert rd["document"]["window_id"] == wid, rd["document"]
        assert "bond_length" in rd.get("chemdraw_settings", {}), rd.get("chemdraw_settings")
        # ChemDraw chemically interprets the reagent caption: "H2SO4" comes back either as a caption or as a
        # fragment (molecule) carrying the label text / sulfuric acid SMILES
        texts = [t["text"] for t in rd["texts"]]
        texts += [str(m.get("text") or "") for m in rd["molecules"]] + [str(m.get("smiles") or "") for m in rd["molecules"]]
        texts += [lab.get("text", "") if isinstance(lab, dict) else str(lab) for m in rd["molecules"] for lab in (m.get("labels") or [])]
        norm = [t.replace("₂", "2").replace("₄", "4") for t in texts]
        assert any("H2SO4" in t or t in ("OS(=O)(=O)O", "O=S(=O)(O)O") for t in norm), norm
    finally:
        closed = loads(server.close_document(document=str(wid)))
    assert closed.get("ok") is True, closed
    live.created.discard(wid)
    assert wid not in live.window_ids()
    assert live.window_ids() == before


def test_server_structure_to_name(server, live):
    before = live.window_ids()
    out = loads(server.structure_to_name(structure="CC(=O)OCC"))
    assert out.get("ok") is True, out
    assert str(out["name"]).strip().lower() == "ethyl acetate", out
    assert canon(out["smiles"]) == canon("CCOC(C)=O")
    assert live.window_ids() == before          # the temporary document was closed


def test_server_quick_render_leaves_no_document_open(server, live, out_dir):
    before = live.window_ids()
    result = server.quick_render("caffeine", return_image=False)
    info = loads(result)
    assert "error" not in info, info
    assert info["rendered_by"] == "ChemDraw" and info["format"] == "png", info
    assert_png(info["path"], out_dir)
    assert os.path.isfile(info["cdxml_path"]) and is_under(info["cdxml_path"], out_dir)
    assert live.window_ids() == before, "quick_render left a document open"


def test_server_new_insert_selection_render_with_image(server, live, out_dir):
    from chemdrawmcp.mcpcompat import Image
    created = loads(server.new_document(style="acs"))
    assert created.get("ok") is True, created
    wid = created["window_id"]
    live.created.add(wid)
    try:
        assert isinstance(wid, int)
        ins = loads(server.insert(data=PHENOL, kind="smiles", document=str(wid)))
        assert ins.get("atoms_added") == 7 and ins.get("window_id") == wid, ins
        sel = loads(server.selection(document=str(wid)))
        assert canon(sel["smiles"]) == canon("Oc1ccccc1")
        assert loads(server.selection(document=str(wid), as_format="inchikey"))["result"].strip() == PHENOL_INCHIKEY

        result = server.render(document=str(wid), return_image=True)
        assert isinstance(result, list) and len(result) == 2, type(result)
        text, image = result
        info = loads(text)
        assert info["rendered_by"] == "ChemDraw" and info["format"] == "png", info
        assert_png(info["path"], out_dir)
        assert isinstance(image, Image)
        data = getattr(image, "data", None)
        if data is not None:
            assert bytes(data)[:8] == PNG_MAGIC

        plain = server.render(document=str(wid), return_image=False)
        assert isinstance(plain, str)
    finally:
        closed = loads(server.close_document(document=str(wid)))
    assert closed.get("ok") is True, closed
    live.created.discard(wid)
    assert wid not in live.window_ids()
