"""Offline tests for chemdrawmcp.server: every MCP tool exercised in headless mode (no ChemDraw).

With mcp >= 2 ``@mcp.tool()`` returns the undecorated function, so the tools are called directly.  The
``server`` fixture forces ``CHEMDRAWMCP_BACKEND=none`` (NullBackend), points ``CHEMDRAWMCP_OUTPUT_DIR`` at a
per-test directory and disables the online name resolver before (re)importing the module.
"""
from __future__ import annotations

import importlib
import json
import os
import re
import time
import types

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from chemdrawmcp import cdxml, chem  # noqa: E402
from chemdrawmcp.mcpcompat import Image  # noqa: E402

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ESTERIFICATION = "CC(=O)O + OCC --[H2SO4]--> CC(=O)OCC + O"
ESTERIFICATION_RXN = "CC(=O)O.CCO>>CCOC(C)=O.O"


# ---------------------------------------------------------------------------
# fixtures & helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def out_dir(tmp_path):
    return tmp_path / "renders"


@pytest.fixture
def server(tmp_path, out_dir, monkeypatch):
    """chemdrawmcp.server imported headless, writing into ``out_dir``; never touches ChemDraw or the network."""
    monkeypatch.setenv("CHEMDRAWMCP_BACKEND", "none")
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "0")
    monkeypatch.setenv("CHEMDRAWMCP_OUTPUT_DIR", str(out_dir))
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path / "home"))
    monkeypatch.delenv("CHEMDRAWMCP_STYLE", raising=False)
    monkeypatch.delenv("CHEMDRAWMCP_AUTOLAUNCH", raising=False)
    import chemdrawmcp.server as srv
    srv = importlib.reload(srv)
    srv._BACKEND = None
    srv._NAME_CACHE.clear()
    yield srv
    srv._BACKEND = None
    srv._NAME_CACHE.clear()


def load(result):
    """A tool's JSON string -> dict (asserting it is a plain string, not [text, Image])."""
    assert isinstance(result, str), "expected a JSON string, got %s" % type(result).__name__
    return json.loads(result)


def unwrap(result):
    """JSON text of a tool result that may be either a string or ``[text, Image]``."""
    if isinstance(result, list):
        assert isinstance(result[0], str)
        return json.loads(result[0])
    return load(result)


def canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, "not SMILES: %r" % smiles
    return Chem.MolToSmiles(mol)


def canon_side(side: str):
    return sorted(canon(s) for s in side.split(".") if s)


def canon_rxn(reaction_smiles: str):
    """(reactants, products) of a reaction SMILES, each canonicalised and sorted (agents ignored)."""
    parts = reaction_smiles.split(">")
    assert len(parts) == 3, reaction_smiles
    return canon_side(parts[0]), canon_side(parts[2])


def is_under(path: str, root) -> bool:
    path, root = os.path.realpath(path), os.path.realpath(str(root))
    return os.path.commonpath([path, root]) == root


def read(path: str, binary: bool = False):
    with open(path, "rb" if binary else "r", encoding=None if binary else "utf-8") as fh:
        return fh.read()


def assert_png(path: str, root=None):
    assert path.lower().endswith(".png")
    assert os.path.isfile(path)
    assert read(path, binary=True)[:8] == PNG_MAGIC
    if root is not None:
        assert is_under(path, root)


def texts_of(xml: str):
    return [t["text"] for t in cdxml.read_cdxml(xml, want_mols=False)["texts"]]


# ---------------------------------------------------------------------------
# module shape / status
# ---------------------------------------------------------------------------

def test_tools_are_plain_functions(server):
    for name in ("status", "draw_structures", "draw_reaction", "draw_scheme", "quick_render", "properties",
                 "convert", "analyze_reaction", "name_to_structure", "styles", "documents", "screenshot"):
        assert isinstance(getattr(server, name), types.FunctionType), name


def test_status_headless(server, out_dir):
    s = load(server.status())
    assert s["backend"] == "none"
    assert s["chemdraw_available"] is False
    assert s["chemdraw_app"] is None
    assert "headless_reason" in s and "headless" in s["headless_reason"].lower()
    assert "acs" in s["styles"]
    assert s["default_style"] == "acs"
    assert s["rdkit"] is True
    assert os.path.realpath(s["renders_dir"]) == os.path.realpath(str(out_dir))
    assert "running" not in s  # NullBackend short-circuits before any liveness check


# ---------------------------------------------------------------------------
# draw_structures
# ---------------------------------------------------------------------------

def test_draw_structures_headless(server, out_dir):
    res = server.draw_structures(["CC(=O)O", "aspirin", "c1ccccc1"], labels=["acetic acid", None, "benzene"],
                                 open_in_chemdraw=False, return_image=False)
    info = load(res)
    assert info["opened"] is False
    assert "ChemDraw not available" in info["note"]
    assert os.path.isfile(info["cdxml_path"]) and info["cdxml_path"].endswith(".cdxml")
    assert is_under(info["cdxml_path"], out_dir)
    assert len(info["molecules"]) == 3
    assert info["layout"]["style"] == "ACS Document 1996"

    render = info["render"]
    assert "RDKit" in render["rendered_by"]
    assert_png(render["path"], out_dir)
    assert "summary" in render["analysis"]

    doc = cdxml.read_cdxml(read(info["cdxml_path"]))
    assert len(doc["molecules"]) == 3
    assert [canon(m["smiles"]) for m in doc["molecules"]] == [canon("CC(=O)O"), canon(ASPIRIN), canon("c1ccccc1")]
    assert [t["text"] for t in doc["texts"]] == ["acetic acid", "benzene"]
    assert doc["steps"] == []


def test_draw_structures_return_image(server):
    res = server.draw_structures(["CCO"], open_in_chemdraw=False, return_image=True)
    assert isinstance(res, list) and len(res) == 2
    assert isinstance(res[0], str)
    assert isinstance(res[1], Image)
    assert os.path.isfile(json.loads(res[0])["cdxml_path"])


def test_draw_structures_names_become_labels_and_render_can_be_skipped(server):
    info = load(server.draw_structures(["aspirin"], open_in_chemdraw=False, render=False, return_image=False))
    assert "render" not in info
    assert texts_of(read(info["cdxml_path"])) == ["aspirin"]


def test_draw_structures_bad_input_is_error_json(server):
    err = load(server.draw_structures(["this is not a molecule"], open_in_chemdraw=False, return_image=False))
    assert err["ok"] is False
    assert err["error_type"] == "ChemError"


# ---------------------------------------------------------------------------
# draw_reaction
# ---------------------------------------------------------------------------

def test_draw_reaction_arrow_syntax(server, out_dir):
    info = load(server.draw_reaction(reaction=ESTERIFICATION, conditions_below=["reflux"],
                                     open_in_chemdraw=False, return_image=False))
    assert info["opened"] is False
    assert canon_rxn(info["reaction_smiles"]) == canon_rxn(ESTERIFICATION_RXN)
    assert is_under(info["cdxml_path"], out_dir)

    xml = read(info["cdxml_path"])
    assert "<scheme" in xml and "<step" in xml
    assert "H2SO4" in xml and "reflux" in xml
    doc = cdxml.read_cdxml(xml)
    assert len(doc["steps"]) == 1
    step = doc["steps"][0]
    assert len(step["reactants"]) == 2 and len(step["products"]) == 2
    assert "H2SO4" in [t.get("text") for t in step["above_arrow"]]
    assert "reflux" in [t.get("text") for t in step["below_arrow"]]
    assert canon_rxn(doc["reaction_smiles"][0]) == canon_rxn(ESTERIFICATION_RXN)

    assert "RDKit" in info["render"]["rendered_by"]
    assert_png(info["render"]["path"], out_dir)


def test_draw_reaction_product_labels_and_equilibrium_arrow(server):
    info = load(server.draw_reaction(reaction=ESTERIFICATION, product_labels=["85%", None], arrow="equilibrium",
                                     open_in_chemdraw=False, return_image=False))
    xml = read(info["cdxml_path"])
    assert 'ArrowheadHead="HalfLeft"' in xml
    doc = cdxml.read_cdxml(xml)
    assert len(doc["arrows"]) == 1
    assert doc["arrows"][0]["arrowheadHead"] == "HalfLeft"
    assert "85%" in texts_of(xml)
    assert canon_rxn(info["reaction_smiles"]) == canon_rxn(ESTERIFICATION_RXN)


def test_draw_reaction_from_component_lists(server):
    info = load(server.draw_reaction(reactants=["CC(=O)O", "OCC"], products=["CC(=O)OCC", "O"],
                                     reagents_above=["H2SO4"], conditions_below=["reflux"],
                                     open_in_chemdraw=False, return_image=False))
    assert canon_rxn(info["reaction_smiles"]) == canon_rxn(ESTERIFICATION_RXN)
    xml = read(info["cdxml_path"])
    assert "<scheme" in xml
    texts = texts_of(xml)
    assert "H2SO4" in texts and "reflux" in texts


def test_draw_reaction_without_input_is_error_json(server):
    err = load(server.draw_reaction(open_in_chemdraw=False, return_image=False))
    assert err["ok"] is False
    assert "reactants" in err["error"] and "products" in err["error"]


# ---------------------------------------------------------------------------
# draw_scheme
# ---------------------------------------------------------------------------

def test_draw_scheme(server, out_dir):
    steps = [
        {"reaction": "c1ccccc1 --[Br2, FeBr3]--> c1ccccc1Br"},
        {"reactants": ["c1ccccc1Br"], "products": ["c1ccccc1C(=O)O"], "above": ["1. Mg", "2. CO2"], "below": ["3. H3O+"]},
    ]
    info = load(server.draw_scheme(steps=steps, title="Scheme 1", open_in_chemdraw=False, return_image=False))
    assert info["opened"] is False
    assert is_under(info["cdxml_path"], out_dir)
    assert isinstance(info["reaction_smiles"], list) and len(info["reaction_smiles"]) == 2
    assert canon_rxn(info["reaction_smiles"][0]) == canon_rxn("c1ccccc1>>Brc1ccccc1")
    assert canon_rxn(info["reaction_smiles"][1]) == canon_rxn("Brc1ccccc1>>OC(=O)c1ccccc1")

    xml = read(info["cdxml_path"])
    doc = cdxml.read_cdxml(xml)
    assert len(doc["steps"]) == 2
    texts = texts_of(xml)
    assert "Scheme 1" in texts
    for reagent in ("Br2, FeBr3", "1. Mg", "2. CO2", "3. H3O+"):
        assert reagent in texts
    # chained: bromobenzene is drawn once and shared by both steps
    assert len(doc["molecules"]) == 3

    assert "RDKit" in info["render"]["rendered_by"]
    assert_png(info["render"]["path"], out_dir)


# ---------------------------------------------------------------------------
# quick_render
# ---------------------------------------------------------------------------

def test_quick_render_png(server, out_dir):
    # see test_quick_render_png_honours_return_image_false for the return type
    info = unwrap(server.quick_render("caffeine", return_image=False))
    assert info["format"] == "png"
    assert "RDKit" in info["rendered_by"]
    assert_png(info["path"], out_dir)
    assert "summary" in info["analysis"]


def test_quick_render_png_honours_return_image_false(server):
    assert isinstance(server.quick_render("caffeine", return_image=False), str)


def test_quick_render_svg(server, out_dir):
    info = load(server.quick_render("caffeine", format="svg", return_image=False))
    assert info["format"] == "svg"
    assert "RDKit" in info["rendered_by"]
    assert info["path"].endswith(".svg") and is_under(info["path"], out_dir)
    assert "<svg" in read(info["path"])


def test_quick_render_reaction_svg_to_given_path(server, tmp_path):
    target = tmp_path / "custom" / "ester.svg"
    info = load(server.quick_render("CC(=O)O + OCC -> CC(=O)OCC + O", format="svg", path=str(target), return_image=False))
    assert os.path.realpath(info["path"]) == os.path.realpath(str(target))
    assert "<svg" in read(str(target))


# ---------------------------------------------------------------------------
# chemistry tools
# ---------------------------------------------------------------------------

def test_properties_aspirin(server):
    out = load(server.properties(structure="aspirin"))
    rd = out["rdkit"]
    assert rd["formula"] == "C9H8O4"
    assert rd["molecular_weight"] == pytest.approx(180.159, abs=0.01)
    assert rd["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert canon(rd["canonical_smiles"]) == canon(ASPIRIN)
    assert "logp" in rd and "tpsa" in rd  # extended=True by default
    assert "chemdraw" not in out


def test_properties_without_structure_is_error(server):
    err = load(server.properties())
    assert err["ok"] is False
    assert "structure" in err["error"]
    assert err["error_type"] == "ChemError"


def test_convert_inchikey(server):
    out = load(server.convert("c1ccccc1O", to="inchikey"))
    assert out["ok"] is True
    assert out["result"] == "ISWSIDIOOBJBQZ-UHFFFAOYSA-N"
    assert out["converted_by"] == "RDKit"


def test_convert_smiles_is_canonically_phenol(server):
    out = load(server.convert("c1ccccc1O", to="smiles"))
    assert out["ok"] is True
    assert canon(out["result"]) == canon("Oc1ccccc1")


def test_convert_mol(server):
    out = load(server.convert("c1ccccc1O", to="mol"))
    assert out["ok"] is True and "V2000" in out["result"]
    assert Chem.MolFromMolBlock(out["result"]).GetNumAtoms() == 7


def test_convert_cdxml(server):
    out = load(server.convert("c1ccccc1O", to="cdxml"))
    assert out["ok"] is True and "<fragment" in out["result"]
    doc = cdxml.read_cdxml(out["result"])
    assert len(doc["molecules"]) == 1
    assert canon(doc["molecules"][0]["smiles"]) == canon("Oc1ccccc1")


def test_convert_png_to_path(server, tmp_path):
    target = tmp_path / "phenol.png"
    out = load(server.convert("c1ccccc1O", to="png", path=str(target)))
    assert out["ok"] is True
    assert os.path.realpath(out["path"]) == os.path.realpath(str(target))
    assert out["converted_by"] == "RDKit"
    assert_png(str(target))


def test_convert_cdx_needs_chemdraw(server):
    err = load(server.convert("c1ccccc1O", to="cdx"))
    assert err["ok"] is False
    assert "ChemDraw" in err["error"]


def test_analyze_reaction_balanced(server):
    out = load(server.analyze_reaction("CC(=O)O + OCC -> CC(=O)OCC + O"))
    assert out["balanced"] is True
    assert out["unbalanced_elements"] == []
    assert canon_rxn(out["reaction_smiles"]) == canon_rxn(ESTERIFICATION_RXN)
    assert out["element_balance"]["C"] == [4, 4]
    assert out["atom_economy"] == pytest.approx(83.0, abs=0.1)


def test_analyze_reaction_unbalanced(server):
    out = load(server.analyze_reaction("CC(=O)O -> CCO"))
    assert out["balanced"] is False
    assert set(out["unbalanced_elements"]) == {"H", "O"}


def test_analyze_reaction_resolves_names_offline(server):
    out = load(server.analyze_reaction("acetic acid + ethanol -> ethyl acetate + water"))
    assert out["balanced"] is True


def test_name_to_structure_aspirin(server):
    out = load(server.name_to_structure("aspirin"))
    assert canon(out["smiles"]) == canon(ASPIRIN)
    assert out["source"] in ("offline resolver", "offline dictionary")
    assert out["properties"]["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert "drawn" not in out


def test_name_to_structure_unknown_name(server):
    err = load(server.name_to_structure("notarealchemicalxyz"))
    assert err["ok"] is False
    assert "notarealchemicalxyz" in err["error"]


def test_styles(server):
    out = load(server.styles())
    assert out["default"] == "acs"
    assert out["styles"]["acs"]["bond_length"] == 14.4
    assert out["styles"]["acs"]["name"] == "ACS Document 1996"
    assert {"acs", "nature", "rsc"} <= set(out["styles"])
    assert "page" in out["styles"]["acs"]


# ---------------------------------------------------------------------------
# tools that need the live application
# ---------------------------------------------------------------------------

NEEDS_CHEMDRAW = [
    ("documents", (), {}),
    ("new_document", (), {}),
    ("read_document", (), {}),
    ("render", (), {}),
    ("structure_to_name", ("CCO",), {}),
    ("run_command", ("selectAll",), {}),
    ("commands", (), {}),
    ("selection", (), {}),
    ("insert", ("CCO",), {}),
    ("screenshot", (), {}),
    ("clean_structure", (), {}),
    ("predict_nmr", ("CCO",), {}),
    ("check_structure", (), {}),
    ("launch", (), {}),
    ("quit_chemdraw", (), {}),
    ("activate_chemdraw", (), {}),
    ("open_document", ("/nonexistent/x.cdxml",), {}),
    ("close_document", (), {}),
    ("save_document", ("/nonexistent/x.png",), {}),
    ("document_settings", (), {}),
    ("add_caption", ("hello",), {}),
    ("add_arrow", (), {}),
    ("select", (), {}),
    ("menus", (), {}),
    ("do_menu", ("File", "New"), {}),
    ("dismiss_dialogs", (), {}),
    ("convert", ("CCO",), {"to": "name"}),
]


@pytest.mark.parametrize("tool,args,kwargs", NEEDS_CHEMDRAW, ids=[t[0] for t in NEEDS_CHEMDRAW])
def test_tools_needing_chemdraw_return_not_available_error(server, tool, args, kwargs):
    res = getattr(server, tool)(*args, **kwargs)
    assert isinstance(res, str), "%s returned %s, not a JSON string" % (tool, type(res).__name__)
    err = json.loads(res)
    assert err["ok"] is False
    assert err["error_type"] == "ChemDrawNotAvailable"
    assert "headless" in err["error"]
    assert "ChemDraw is not available" in err["hint"]


def test_headless_tools_write_nothing_for_chemdraw_only_calls(server, out_dir):
    server.screenshot()
    server.render()
    assert not out_dir.exists() or os.listdir(str(out_dir)) == []


# ---------------------------------------------------------------------------
# image tools
# ---------------------------------------------------------------------------

def test_analyze_image_and_compare_renders(server):
    info = load(server.draw_structures(["CC(=O)O", "aspirin", "c1ccccc1"], open_in_chemdraw=False, return_image=False))
    png = info["render"]["path"]

    a = load(server.analyze_image(png))
    assert isinstance(a["summary"], str) and a["summary"]
    assert a["width"] > 0 and a["height"] > 0
    assert a["empty"] is False
    assert a["background"]["name"] == "white"

    d = load(server.compare_renders(png, png))
    assert d["changed_fraction"] == 0
    assert d["changed_pixels"] == 0
    assert d["bbox"] is None
    assert "identical" in d["summary"]


def test_analyze_image_missing_file_is_error(server, tmp_path):
    err = load(server.analyze_image(str(tmp_path / "missing.png")))
    assert err["ok"] is False


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def test_doc_normalisation(server):
    assert server._doc(None) is None
    assert server._doc("") is None
    assert server._doc("front") is None
    assert server._doc("#12") == 12
    assert server._doc("12") == 12
    assert server._doc(7) == 7
    assert server._doc("aspirin.mol") == "aspirin.mol"
    assert server._doc("  Untitled ACS Document 1996-1 ") == "Untitled ACS Document 1996-1"


def test_slug(server):
    assert server._slug("Scheme 1: test") == "scheme_1_test"
    assert server._slug("") == "item"
    assert server._slug("!!!") == "item"
    assert len(server._slug("x" * 100)) == 40


def test_resolve_names_in_reaction(server):
    text = server._resolve_names_in_reaction("acetic acid + ethanol -> ethyl acetate")
    for name in ("acetic acid", "ethanol", "ethyl acetate"):
        assert name not in text
    left, right = re.split(r"\s*->\s*", text)
    comps = [c.strip() for side in (left, right) for c in side.split(" + ")]
    assert len(comps) == 3
    for c in comps:
        assert Chem.MolFromSmiles(c) is not None, c
        assert chem.detect_format(c) == "smiles", c
    assert canon(comps[0]) == canon("CC(=O)O")
    assert canon(comps[1]) == canon("CCO")
    assert canon(comps[2]) == canon("CCOC(C)=O")


def test_resolve_names_in_reaction_passthrough(server):
    assert server._resolve_names_in_reaction("CCO>>CC=O") == "CCO>>CC=O"
    assert server._resolve_names_in_reaction("CC(=O)O + OCC -> CC(=O)OCC") == "CC(=O)O + OCC -> CC(=O)OCC"
    assert server._resolve_names_in_reaction("just some text") == "just some text"


def test_err_marks_chemdraw_not_available(server):
    from chemdrawmcp.backends import ChemDrawNotAvailable
    err = json.loads(server._err(ChemDrawNotAvailable("nope"), extra=1))
    assert err == {"ok": False, "error": "nope", "error_type": "ChemDrawNotAvailable", "hint": err["hint"], "extra": 1}
    plain = json.loads(server._err(ValueError("bad")))
    assert plain == {"ok": False, "error": "bad", "error_type": "ValueError"}


def test_output_paths_are_unique_within_one_second(server):
    for _ in range(3):
        second = time.strftime("%Y%m%d_%H%M%S")
        a = load(server.draw_structures(["CCO"], open_in_chemdraw=False, return_image=False))
        b = load(server.draw_structures(["CCN"], open_in_chemdraw=False, return_image=False))
        if time.strftime("%Y%m%d_%H%M%S") == second:
            break
    else:  # pragma: no cover - three straddled second boundaries in a row
        pytest.skip("could not run two calls inside one second")
    assert a["cdxml_path"] != b["cdxml_path"]
    assert a["render"]["path"] != b["render"]["path"]


# ---------------------------------------------------------------------------
# regression tests: reaction coefficients, " = " arrows, name-resolution order, paste counts, `~` paths
# ---------------------------------------------------------------------------

def test_resolve_names_in_reaction_keeps_coefficients(server):
    text = server._resolve_names_in_reaction("2 hydrogen + oxygen -> 2 water")
    left, right = re.split(r"\s*->\s*", text)
    assert [c.split(" ", 1)[0] for c in left.split(" + ")] == ["2", "O=O"]
    assert canon(left.split(" + ")[0].split(" ", 1)[1]) == canon("[H][H]")
    assert right.startswith("2 ") and canon(right[2:]) == canon("O")


def test_analyze_reaction_with_coefficients_and_names_is_balanced(server):
    out = load(server.analyze_reaction("2 hydrogen + oxygen -> 2 water"))
    assert out["balanced"] is True
    assert out["balanced"] == chem.reaction_properties(chem.parse_reaction("2 hydrogen + oxygen -> 2 water"))["balanced"]


def test_resolve_names_in_reaction_accepts_parser_arrows(server):
    for arrow in (" = ", " --[H2SO4]--> ", " > ", " => "):
        text = server._resolve_names_in_reaction("acetic acid + ethanol" + arrow + "ethyl acetate + water")
        assert arrow in text, text
        for name in ("acetic acid", "ethanol", "ethyl acetate", "water"):
            assert name not in text
    assert server._resolve_names_in_reaction("A -> B -> C") == "A -> B -> C"   # multi-step: left to parse_reaction


def test_looks_like_reaction(server):
    assert server._looks_like_reaction("acetic acid + ethanol = ethyl acetate + water")
    assert server._looks_like_reaction("CCO>>CC=O")
    assert server._looks_like_reaction("ethanol > acetaldehyde")
    assert not server._looks_like_reaction("[NH3]->[Cu]")   # dative bond, not an arrow
    assert not server._looks_like_reaction("caffeine")
    assert not server._looks_like_reaction("C=C")
    assert not server._looks_like_reaction(Chem.MolToMolBlock(Chem.MolFromSmiles("CCO")) + ">  <prop>\nvalue\n\n$$$$\n")


def test_quick_render_equals_arrow_reaction(server, out_dir):
    info = unwrap(server.quick_render("acetic acid + ethanol = ethyl acetate + water", return_image=False))
    assert info.get("ok", True), info
    assert info["format"] == "png" and "RDKit" in info["rendered_by"]
    assert_png(info["path"], out_dir)


def test_name_to_structure_reports_opsin_as_source(server, monkeypatch):
    calls = []

    def fake_opsin(name):
        calls.append(name)
        return "CCCCCCCCO" if name == "octan-1-ol" else None
    monkeypatch.setattr(chem, "_opsin", fake_opsin)
    out = load(server.name_to_structure("octan-1-ol"))
    assert canon(out["smiles"]) == canon("CCCCCCCCO")
    assert out["source"] == "OPSIN"
    assert calls == ["octan-1-ol"]
    # dictionary names never reach OPSIN
    assert load(server.name_to_structure("aspirin"))["source"] == "offline dictionary"
    assert calls == ["octan-1-ol"]


def test_name_to_structure_pubchem_only_when_online(server, monkeypatch):
    calls = []
    monkeypatch.setattr(chem, "_opsin", lambda name: None)
    monkeypatch.setattr(chem, "_pubchem", lambda name: calls.append(name) or "CCCCCCCCCCO")
    err = load(server.name_to_structure("decan-1-ol"))
    assert err["ok"] is False and calls == []
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "1")
    out = load(server.name_to_structure("decan-1-ol"))
    assert out["source"] == "PubChem" and calls == ["decan-1-ol"]


def test_analyze_image_and_compare_renders_expand_tilde(server, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    (tmp_path / "r").mkdir()
    png = load(server.convert("CCO", "png", path=str(tmp_path / "r" / "x.png")))["path"]
    assert os.path.isfile(png)
    a = load(server.analyze_image("~/r/x.png"))
    assert a.get("ok", True) and a["width"] > 0
    assert os.path.realpath(a["path"]) == os.path.realpath(png)
    d = load(server.compare_renders("~/r/x.png", png))
    assert d.get("ok", True) and d["changed_fraction"] == 0.0
    missing = load(server.analyze_image("~/r/nope.png"))
    assert missing["ok"] is False and missing["error_type"] == "ChemError"
    assert "nope.png" in missing["error"] and "~" not in missing["error"]


def test_convert_text_result_to_path_creates_directory(server, tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    out = load(server.convert("CCO", "mol", path="~/newdir/deeper/x.mol", use_chemdraw=False))
    assert out["ok"] is True
    assert os.path.realpath(out["path"]) == os.path.realpath(str(tmp_path / "newdir" / "deeper" / "x.mol"))
    assert "M  END" in read(out["path"])
    assert load(server.convert("CCO", "smiles", use_chemdraw=False))["path"] is None


# ---------------------------------------------------------------------------
# a fake live backend: _live()'s auto-launch gating, status(), paste counts, name-resolution order
# ---------------------------------------------------------------------------

from chemdrawmcp.backends import ChemDrawError  # noqa: E402
from chemdrawmcp.backends.base import ChemDrawBackend, CommandInfo, DocumentInfo  # noqa: E402

NOT_RUNNING = "ChemDraw is not running; call `launch` first"


class FakeBackend(ChemDrawBackend):
    """In-memory stand-in for a live ChemDraw (installed, initially not running): records every call."""

    platform = "fake"

    def __init__(self):
        self.running = False
        self.calls = []
        self.docs = {}            # window id -> DocumentInfo
        self.name_structs = {}    # what "Name=Struct" knows: name -> SMILES
        self._next = 100
        self._last_smiles = ""

    def _find(self, name):
        if name is None:
            return next((d for d in self.docs.values() if d.front), None) or next(iter(self.docs.values()))
        for d in self.docs.values():
            if d.name == name or d.window_id == name or (isinstance(name, str) and name.isdigit() and d.window_id == int(name)):
                return d
        raise ChemDrawError("no document %r" % (name,))

    def executable(self):
        return "/Applications/Fake ChemDraw.app"

    def is_running(self):
        return self.running

    def launch(self, wait_seconds=90.0):
        self.calls.append("launch")
        self.running = True
        return {"launched": True}

    def quit(self, save=False):
        self.calls.append("quit")
        self.running = False

    def version(self):
        self.calls.append("version")
        return "fake 1.0"

    def activate(self):
        self.calls.append("activate")

    def documents(self):
        self.calls.append("documents")
        return list(self.docs.values())

    def front_document(self):
        self.calls.append("front_document")
        d = self._find(None)
        return d.name if d else None

    def activate_document(self, name):
        self.calls.append("activate_document")
        return self._find(name).name

    def new_document(self, stationery=None):
        self.calls.append("new_document")
        for d in self.docs.values():
            d.front = False
        wid = self._next
        self._next += 1
        self.docs[wid] = DocumentInfo(name="Untitled-%d" % wid, window_id=wid, front=True)
        return self.docs[wid].name

    def open_document(self, path):
        self.calls.append("open_document")
        name = self.new_document()
        self._find(name).path = path
        return name

    def close_document(self, name=None, save=False):
        self.calls.append("close_document")
        d = self._find(name)
        if d is not None:
            del self.docs[d.window_id]

    def save_document(self, name, path, fmt):
        self.calls.append("save_document")
        with open(path, "wb") as fh:
            fh.write(chem.depict_png(chem.mol_from_input("CCO", "smiles")))
        return path

    def document_cdxml(self, name=None):
        self.calls.append("document_cdxml")
        return cdxml.CDXMLDocument("acs").to_xml()

    def document_settings(self, name=None):
        self.calls.append("document_settings")
        return {}

    def set_document_settings(self, name, settings):
        self.calls.append("set_document_settings")
        return dict(settings)

    def insert(self, name, data, kind="cdxml"):
        self.calls.append("insert")
        d = self._find(name)
        before = {"atoms": d.atoms, "molecules": d.molecules, "captions": d.captions}
        if kind == "name":
            smi = self.name_structs.get(data)
            if smi:
                d.atoms += Chem.MolFromSmiles(smi).GetNumAtoms()
                d.molecules += 1
                self._last_smiles = smi
        elif kind == "cdxml":
            d.captions += data.count("<t ")
            d.arrows += data.count("<arrow ")
        after = {"atoms": d.atoms, "molecules": d.molecules, "captions": d.captions}   # as the macOS backend: no arrows
        return {"document": d.name, "before": before, "after": after, "kind": kind,
                "atoms_added": after["atoms"] - before["atoms"]}

    def select_all(self, name=None):
        self.calls.append("select_all")

    def deselect(self, name=None):
        self.calls.append("deselect")

    def selection_info(self, name=None):
        self.calls.append("selection_info")
        return {"empty": True}

    def selection_cdxml(self, name=None):
        self.calls.append("selection_cdxml")
        return ""

    def selection_as(self, name, what):
        self.calls.append("selection_as")
        return self._last_smiles if what == "smiles" else ""

    def list_commands(self, pattern=None, enabled_only=False, with_descriptions=True):
        self.calls.append("list_commands")
        return [CommandInfo("selectAll", True)]

    def run_command(self, command, name=None):
        self.calls.append("run_command")
        return {"ok": True, "command": command, "enabled": True}

    def menus(self):
        self.calls.append("menus")
        return {"File": ["New"]}

    def do_menu(self, menu, item, name=None):
        self.calls.append("do_menu")
        return {"ok": True}

    def screenshot(self, name=None):
        self.calls.append("screenshot")
        return chem.depict_png(chem.mol_from_input("CCO", "smiles"))

    def dismiss_dialogs(self):
        self.calls.append("dismiss_dialogs")
        return []


@pytest.fixture
def fake(server):
    """``server`` with a FakeBackend installed as the live backend: ChemDraw "installed" but not running."""
    b = FakeBackend()
    server._BACKEND = b
    server._NAME_SOURCE.clear()
    yield b
    server._BACKEND = None


def test_live_refuses_to_launch_when_autolaunch_disabled(server, fake, monkeypatch):
    monkeypatch.setenv("CHEMDRAWMCP_AUTOLAUNCH", "0")
    err = load(server.documents())
    assert err == {"ok": False, "error": NOT_RUNNING, "error_type": "ChemDrawError"}
    assert fake.calls == [] and fake.running is False


def test_live_launches_once_by_default(server, fake, monkeypatch):
    monkeypatch.delenv("CHEMDRAWMCP_AUTOLAUNCH", raising=False)
    out = load(server.documents())
    assert out["documents"] == [] and fake.running is True
    assert fake.calls == ["launch", "documents"]
    load(server.documents())
    assert fake.calls.count("launch") == 1


NO_AUTOLAUNCH = [
    ("render", (), {}),
    ("screenshot", (), {}),
    ("close_document", (), {}),
    ("read_document", (), {}),
    ("selection", (), {}),
    ("select", (), {}),
    ("run_command", ("selectAll",), {}),
    ("commands", (), {}),
    ("menus", (), {}),
    ("do_menu", ("File", "New"), {}),
    ("dismiss_dialogs", (), {}),
    ("clean_structure", (), {}),
    ("check_structure", (), {}),
    ("document_settings", (), {}),
    ("save_document", ("/nonexistent/x.png",), {}),
    ("quit_chemdraw", (), {}),
]


@pytest.mark.parametrize("tool,args,kwargs", NO_AUTOLAUNCH, ids=[t[0] for t in NO_AUTOLAUNCH])
def test_read_only_tools_never_start_chemdraw(server, fake, monkeypatch, tool, args, kwargs):
    monkeypatch.delenv("CHEMDRAWMCP_AUTOLAUNCH", raising=False)
    res = getattr(server, tool)(*args, **kwargs)
    err = json.loads(res if isinstance(res, str) else res[0])
    assert err["ok"] is False and err["error"] == NOT_RUNNING and err["error_type"] == "ChemDrawError"
    assert fake.calls == [] and fake.running is False


def test_status_with_running_backend(server, fake):
    fake.running = True
    fake.new_document()
    s = load(server.status())
    assert s["backend"] == "fake" and s["chemdraw_available"] is True
    assert s["running"] is True and s["version"] == "fake 1.0"
    assert [d["name"] for d in s["documents"]] == ["Untitled-100"]
    assert s["front_document"] == "Untitled-100"
    assert "headless_reason" not in s and "pending_dialogs" not in s


def test_add_caption_and_add_arrow_report_real_counts(server, fake):
    fake.running = True
    name = fake.new_document()
    c = load(server.add_caption("hello", document=name))
    assert c["ok"] is True and c["document"] == name and c["window_id"] == 100
    assert c["captions_after"] == 1
    a = load(server.add_arrow(text_above="H2SO4", document=name))
    assert a["ok"] is True and a["window_id"] == 100
    assert a["arrows_after"] == 1
    assert load(server.add_arrow(document=name))["arrows_after"] == 2
    assert load(server.add_caption("again", document=name))["captions_after"] == 3   # 1 + the arrow's text + this one


def test_name_resolution_order_dictionary_opsin_chemdraw(server, fake, monkeypatch):
    opsin_calls = []
    monkeypatch.setattr(chem, "_opsin", lambda name: opsin_calls.append(name) or None)
    fake.name_structs["octan-1-ol"] = "CCCCCCCCO"
    # dictionary names never touch ChemDraw
    out = load(server.name_to_structure("aspirin"))
    assert out["source"] == "offline dictionary" and fake.calls == [] and opsin_calls == []
    # other names try OPSIN (offline, instant) first; when it has no answer ChemDraw's Name=Struct is used
    # (the scratch document is closed again)
    out = load(server.name_to_structure("octan-1-ol"))
    assert canon(out["smiles"]) == canon("CCCCCCCCO")
    assert out["source"] == "ChemDraw Name=Struct"
    assert "insert" in fake.calls and fake.calls.count("close_document") == 1 and fake.docs == {}
    assert opsin_calls == ["octan-1-ol"]
    # what neither OPSIN nor ChemDraw knows is unresolved (PubChem only with CHEMDRAWMCP_ONLINE=1)
    assert load(server.name_to_structure("nonan-1-ol"))["ok"] is False
    assert opsin_calls == ["octan-1-ol", "nonan-1-ol"]
