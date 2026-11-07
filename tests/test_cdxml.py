"""Tests for the CDXML writer/reader and the layout engines (no ChemDraw needed)."""
import re
import xml.etree.ElementTree as ET

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdDepictor  # noqa: E402

from chemdrawmcp import cdxml, layout, styles  # noqa: E402


def prep(smi):
    m = Chem.MolFromSmiles(smi)
    Chem.Kekulize(m, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    Chem.WedgeMolBonds(m, m.GetConformer())
    return m


def test_styles_lookup():
    assert styles.get_style("acs")["bond_length"] == 14.4
    assert styles.get_style("ACS Document 1996")["label_font"] == "Arial"
    assert styles.get_style("nature")["label_size"] == 6
    assert styles.get_style("rsc")["label_font"] == "Helvetica"
    with pytest.raises(KeyError):
        styles.get_style("no such style")
    assert "acs" in styles.style_names()


def test_fragment_roundtrip_and_header():
    doc = cdxml.CDXMLDocument("acs", name="t")
    ref = doc.add_fragment(prep("CC(=O)Oc1ccccc1C(=O)O"), 100, 100)
    assert ref.kind == "fragment" and ref.extra["atoms"] == 13
    xml = doc.to_xml()
    root = ET.fromstring(re.sub(r"<!DOCTYPE[^>]*>", "", xml))
    assert root.tag == "CDXML" and root.get("BondLength") == "14.4"
    frag = root.find(".//fragment")
    assert len(frag.findall("n")) == 13 and len(frag.findall("b")) == 13
    # carbons carry no Element attribute, oxygens do with hydrogens counted
    elements = [n.get("Element") for n in frag.findall("n")]
    assert elements.count(None) == 9 and elements.count("8") == 4
    assert any(n.get("NumHydrogens") == "1" for n in frag.findall("n") if n.get("Element") == "8")
    # median bond length equals the style's bond length
    pts = {n.get("id"): tuple(float(v) for v in n.get("p").split()) for n in frag.findall("n")}
    lens = [((pts[b.get("B")][0] - pts[b.get("E")][0]) ** 2 + (pts[b.get("B")][1] - pts[b.get("E")][1]) ** 2) ** 0.5 for b in frag.findall("b")]
    assert abs(sorted(lens)[len(lens) // 2] - 14.4) < 0.5
    back = cdxml.read_cdxml(xml)
    assert back["molecules"][0]["smiles"] == "CC(=O)Oc1ccccc1C(=O)O"
    assert back["molecules"][0]["formula"] == "C9H8O4"
    assert abs(back["molecules"][0]["molecular_weight"] - 180.159) < 0.01


def test_stereo_wedges_survive():
    doc = cdxml.CDXMLDocument("acs")
    doc.add_fragment(prep("N[C@@H](C)C(=O)O"), 50, 50)
    xml = doc.to_xml()
    assert "WedgeBegin" in xml or "WedgedHashBegin" in xml
    back = cdxml.read_cdxml(xml)
    assert "@" in back["molecules"][0]["smiles"]
    assert Chem.CanonSmiles(back["molecules"][0]["smiles"]) == Chem.CanonSmiles("N[C@@H](C)C(=O)O")


def test_charges_isotopes_and_radicals():
    doc = cdxml.CDXMLDocument("acs")
    doc.add_fragment(prep("[NH3+]CC([O-])=O"), 50, 50)
    doc.add_fragment(prep("[13CH3][Mg+]"), 150, 50)
    doc.add_fragment(prep("C[CH2]"), 250, 50)
    xml = doc.to_xml()
    assert 'Charge="1"' in xml and 'Charge="-1"' in xml and 'Isotope="13"' in xml and 'Radical="Doublet"' in xml
    back = cdxml.read_cdxml(xml)
    smiles = [m["smiles"] for m in back["molecules"]]
    assert "[NH3+]CC(=O)[O-]" in smiles and "[13CH3][Mg+]" in smiles


def test_text_arrow_plus_and_step():
    doc = cdxml.CDXMLDocument("acs")
    a = doc.add_fragment(prep("CCO"), 50, 100)
    p = doc.add_plus(a.box.x1 + 10, 100)
    b = doc.add_fragment(prep("CC(=O)O"), p.box.x1 + 30, 100)
    ar = doc.add_arrow(b.box.x1 + 10, 100, b.box.x1 + 60, 100, "equilibrium")
    t = doc.add_text("H2SO4", ar.box.cx, 90, formula=True)
    c = doc.add_fragment(prep("CC(=O)OCC"), ar.box.x1 + 40, 100)
    doc.add_step([a, b], [c], [ar], [t], [], [p])
    xml = doc.to_xml()
    assert 'ArrowheadHead="HalfLeft"' in xml and 'face="96"' in xml and "<scheme" in xml
    assert re.search(r'ReactionStepReactants="%d %d"' % (a.id, b.id), xml)
    assert re.search(r'ReactionStepProducts="%d"' % c.id, xml)
    back = cdxml.read_cdxml(xml)
    assert len(back["steps"]) == 1
    assert back["steps"][0]["above_arrow"][0]["text"] == "H2SO4"
    assert back["reaction_smiles"][0] == "CCO.CC(=O)O>>CCOC(C)=O"
    assert back["arrows"][0]["arrowheadHead"] == "HalfLeft"
    assert "+" in [t["text"] for t in back["texts"]]


def test_unknown_arrow_kind_and_colors():
    doc = cdxml.CDXMLDocument("acs")
    with pytest.raises(ValueError):
        doc.add_arrow(0, 0, 10, 0, "spiral")
    assert doc.color_index("black") == 0 and doc.color_index("white") == 1
    red = doc.color_index("red")
    assert red >= 2 and doc.color_index("#ff0000") == red           # de-duplicated
    doc.add_fragment(prep("c1ccccc1"), 50, 50, color="red")
    assert 'color="%d"' % red in doc.to_xml()


def test_curved_arrow_geometry():
    doc = cdxml.CDXMLDocument("acs")
    ref = doc.add_curved_arrow(0, 0, 60, 0, 90)
    xml = doc.to_xml()
    m = re.search(r'Center3D="([-\d.]+) ([-\d.]+) 0"', xml)
    cx, cy = float(m.group(1)), float(m.group(2))
    r = ((cx - 0) ** 2 + (cy - 0) ** 2) ** 0.5
    assert abs(r - 60 / (2 * (2 ** 0.5) / 2)) < 0.1      # chord / (2 sin 45°)
    assert 'AngularSize="90"' in xml and ref.kind == "arrow"


def test_page_tiles_when_content_is_wide():
    doc = cdxml.CDXMLDocument("acs")
    for i in range(6):
        doc.add_fragment(prep("c1ccccc1"), 100 + i * 150, 100)
    xml = doc.to_xml()
    page = re.search(r"<page[^>]*>", xml).group(0)
    assert 'WidthPages="2"' in page and 'HeightPages="1"' in page
    assert re.search(r'BoundingBox="0 0 1080 720"', page)


def test_multiline_text_and_faces():
    doc = cdxml.CDXMLDocument("acs")
    t = doc.add_text("line one\nline two", 10, 10, face=["bold", "italic"], justification="Left")
    xml = doc.to_xml()
    assert 'face="3"' in xml and t.box.height > doc.caption_size * 1.5
    assert cdxml.document_text(xml) == ["line one\nline two"]


def test_read_chemdraw_style_file_with_nested_label():
    # a fragment exported by ChemDraw: atom labels are nested <t> elements and must be ignored by the reader
    text = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >
<CDXML BondLength="14.40" LabelFont="3" CaptionFont="3"><fonttable><font id="3" charset="iso-8859-1" name="Arial"/></fonttable>
<page id="9"><fragment id="1" BoundingBox="0 0 40 40"><n id="2" p="0 0"/><n id="3" p="14 0" Element="8" NumHydrogens="1"><t p="14 3"><s font="3" size="10" face="96">OH</s></t></n>
<b id="4" B="2" E="3"/></fragment><t id="5" p="0 60"><s font="3" size="10">methanol</s><chemicalproperty Name="name"/></t>
<arrow id="6" BoundingBox="0 0 10 10" ArrowheadHead="Full" Head3D="10 0 0" Tail3D="0 0 0"/></page></CDXML>"""
    back = cdxml.read_cdxml(text)
    assert back["molecules"][0]["smiles"] == "CO" and back["settings"]["LabelFont"] == "Arial"
    assert [t["text"] for t in back["texts"]] == ["methanol"]
    assert back["texts"][0]["chemical_property"]["Name"] == "name"
    assert back["arrows"][0]["head3D"] == "10 0 0"
    mols = cdxml.cdxml_mols(text)
    assert len(mols) == 1 and mols[0].GetNumAtoms() == 2


def test_nickname_expansion():
    text = """<CDXML BondLength="14.4"><page id="1"><fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" NodeType="Nickname">
<t p="14 0"><s font="3" size="10">OMe</s></t><fragment id="5"><n id="6" p="0 0" Element="8"/><n id="7" p="14 0"/><n id="8" p="-14 0" NodeType="ExternalConnectionPoint"/>
<b id="9" B="6" E="7"/><b id="10" B="8" E="6"/></fragment></n><b id="11" B="3" E="4"/></fragment></page></CDXML>"""
    back = cdxml.read_cdxml(text)
    assert back["molecules"][0]["smiles"] == "COC"


# ---------------------------------------------------------------------------
# layouts
# ---------------------------------------------------------------------------

def test_build_structures_grid():
    mols = [prep(s) for s in ("CCO", "c1ccccc1", "CC(=O)O", "C1CCCCC1", "O")]
    doc = layout.build_structures(mols, labels=["ethanol", "benzene", "acetic acid", None, "water"], columns=3)
    xml = doc.to_xml()
    back = cdxml.read_cdxml(xml)
    assert len(back["molecules"]) == 5
    assert [t["text"] for t in back["texts"]] == ["ethanol", "benzene", "acetic acid", "water"]
    boxes = [r.box for r in doc.refs if r.kind == "fragment"]
    # three columns: first three share a row, no overlaps
    assert abs(boxes[0].cy - boxes[1].cy) < 1e-6 and boxes[3].cy > boxes[0].cy
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            assert a.x1 <= b.x0 or b.x1 <= a.x0 or a.y1 <= b.y0 or b.y1 <= a.y0


def test_build_reaction_layout_and_semantics():
    step = layout.Step([layout.Component(prep("CC(=O)O")), layout.Component(prep("OCC"))],
                       [layout.Component(prep("CC(=O)OCC"), label="85%"), layout.Component(prep("O"))],
                       above=["H2SO4 (cat.)"], below=["reflux, 4 h"], arrow="equilibrium")
    doc = layout.build_reaction(step, style="acs")
    xml = doc.to_xml()
    back = cdxml.read_cdxml(xml)
    assert back["reaction_smiles"] == ["CC(=O)O.CCO>>CCOC(C)=O.O"]
    st = back["steps"][0]
    assert [t["text"] for t in st["above_arrow"]] == ["H2SO4 (cat.)"]
    assert [t["text"] for t in st["below_arrow"]] == ["reflux, 4 h"]
    frags = [r for r in doc.refs if r.kind == "fragment"]
    arrow = [r for r in doc.refs if r.kind == "arrow"][0]
    # reactants left of the arrow, products right of it, all on one baseline
    assert all(f.box.x1 <= arrow.box.x0 + 1 for f in frags[:2]) and all(f.box.x0 >= arrow.box.x1 - 1 for f in frags[2:])
    assert max(abs(f.box.cy - frags[0].box.cy) for f in frags) < 1e-6
    # the arrow is long enough for the reagent text
    above = [r for r in doc.refs if r.kind == "text" and r.extra["text"] == "H2SO4 (cat.)"][0]
    assert above.box.width <= arrow.box.width
    assert "+" in [t["text"] for t in back["texts"]] and "85%" in [t["text"] for t in back["texts"]]


def test_atom_mapping_recorded():
    r = prep("[CH3:1][OH:2]")
    p = prep("[CH3:1][Cl:3]")
    doc = layout.build_reaction(layout.Step([layout.Component(r)], [layout.Component(p)]))
    xml = doc.to_xml()
    m = re.search(r'ReactionStepAtomMap="(\d+) (\d+)"', xml)
    assert m, xml
    frag_refs = [x for x in doc.refs if x.kind == "fragment"]
    assert int(m.group(1)) == frag_refs[0].extra["atom_ids"][0] and int(m.group(2)) == frag_refs[1].extra["atom_ids"][0]
    assert cdxml.read_cdxml(xml)["steps"][0]["atom_map_pairs"] == 1


def test_build_scheme_chains_and_wraps():
    A, B, C, D = (prep(s) for s in ("c1ccccc1", "c1ccccc1Br", "c1ccccc1C(=O)O", "c1ccccc1C(=O)OC"))
    steps = [layout.Step([layout.Component(A)], [layout.Component(B)], above=["Br2"], below=["FeBr3"]),
             layout.Step([layout.Component(B)], [layout.Component(C)], above=["1. Mg", "2. CO2"], below=["3. H3O+"]),
             layout.Step([layout.Component(C), layout.Component(prep("CO"))], [layout.Component(D)], above=["H+"])]
    doc = layout.build_scheme(steps, style="acs", max_width=320, title="Scheme 1")
    xml = doc.to_xml()
    back = cdxml.read_cdxml(xml)
    assert len(back["steps"]) == 3
    assert back["reaction_smiles"][0] == "c1ccccc1>>Brc1ccccc1"
    assert back["reaction_smiles"][2].startswith("O=C(O)c1ccccc1.CO>>")
    frag_refs = [x for x in doc.refs if x.kind == "fragment"]
    rows = sorted({round(f.box.cy) for f in frag_refs})
    assert len(rows) >= 2, "narrow max_width must wrap the scheme into rows"
    # the intermediate B is drawn once in row 1 (shared) — benzene bromide appears at most twice (once per row at most)
    smiles = [f.extra["smiles"] for f in frag_refs]
    assert smiles.count("Brc1ccccc1") <= 2
    assert "Scheme 1" in [t["text"] for t in back["texts"]]


def _arrow_centre(xml):
    m = re.search(r'Center3D="([-\d.]+) ([-\d.]+) 0"', xml)
    return float(m.group(1)), float(m.group(2))


def _sweep_about_centre(tail, head, centre, angle):
    """Signed rotation (degrees, screen sense: positive = clockwise with y down) of the |angle| arc from tail to head."""
    import math
    t0 = math.atan2(tail[1] - centre[1], tail[0] - centre[0])
    t1 = math.atan2(head[1] - centre[1], head[0] - centre[0])
    d = (t1 - t0 + math.pi) % (2 * math.pi) - math.pi
    if abs(angle) > 180:
        d = d - 2 * math.pi if d > 0 else d + 2 * math.pi
    return math.degrees(d)


def test_curved_arrow_sign_convention_and_major_arcs():
    """Positive angles bulge to the right of tail->head (downwards for a left-to-right arrow); the written
    Center3D must stay consistent with that rotation sense when the sweep exceeds 180 degrees."""
    import math
    tail, head = (0.0, 0.0), (60.0, 0.0)
    for angle, centre_side in ((90, -1), (-90, 1), (270, 1), (-270, -1), (181, 1), (179, -1)):
        doc = cdxml.CDXMLDocument("acs")
        doc.add_curved_arrow(tail[0], tail[1], head[0], head[1], angle)
        cx, cy = _arrow_centre(doc.to_xml())
        # a minor arc has its centre on the left (y < 0) for angle > 0, a major arc on the right (y > 0)
        assert cx == pytest.approx(30.0) and (cy < 0) == (centre_side < 0), (angle, cx, cy)
        # the arc through tail and head about that centre has the same rotation sense for every magnitude
        sweep = _sweep_about_centre(tail, head, (cx, cy), angle)
        assert sweep == pytest.approx(-angle, abs=0.05), (angle, sweep)      # Center3D is written to 0.01 pt
        # the apex (far point of the arc) lies on the written circle and on the documented side
        r = math.hypot(cx - tail[0], cy - tail[1])
        t0 = math.atan2(tail[1] - cy, tail[0] - cx)
        apex_y = cy + r * math.sin(t0 + math.radians(sweep) / 2.0)
        assert (apex_y > 0) == (angle > 0), (angle, apex_y)
        assert abs(apex_y) == pytest.approx(r + abs(cy) if abs(angle) > 180 else r - abs(cy), abs=0.05)


def test_radicals_beyond_two_electrons_are_not_written_as_triplet():
    """CDXML's Radical enumeration has no value for three unpaired electrons: Triplet means a diradical."""
    for smi, expected in (("[CH]", None), ("[N]", None), ("[C]", None), ("[CH2]", "Singlet"), ("C[CH2]", "Doublet")):
        doc = cdxml.CDXMLDocument("acs")
        doc.add_fragment(prep(smi), 50, 50)
        xml = doc.to_xml()
        found = re.findall(r'Radical="([^"]*)"', xml)
        assert found == ([expected] if expected else []), (smi, found)
        assert "Triplet" not in xml and "Quartet" not in xml
    # the reader keeps ChemDraw's meaning of Triplet (two unpaired electrons) and ignores values ChemDraw lacks
    doc = cdxml.CDXMLDocument("acs")
    doc.add_fragment(prep("[CH2]"), 50, 50)
    xml = doc.to_xml().replace('Radical="Singlet"', 'Radical="Triplet"')
    mol = Chem.MolFromSmiles(cdxml.read_cdxml(xml)["molecules"][0]["smiles"])
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 2
    xml = xml.replace('Radical="Triplet"', 'Radical="Quartet"')
    mol = Chem.MolFromSmiles(cdxml.read_cdxml(xml)["molecules"][0]["smiles"])
    assert all(a.GetNumRadicalElectrons() == 0 for a in mol.GetAtoms())


def test_scheme_atom_map_follows_the_reused_intermediate_atom_order():
    # step 2 writes the intermediate in a different atom order (same canonical SMILES, so the fragment is reused):
    # the reaction map must be built from the atom order of the fragment that was actually drawn
    steps = [layout.Step([layout.Component(prep("[CH3:1][CH2:2]Br"))], [layout.Component(prep("[CH3:1][CH2:2][OH:3]"))], above=["NaOH"]),
             layout.Step([layout.Component(prep("[OH:3][CH2:2][CH3:1]"))], [layout.Component(prep("[CH3:1][CH:2]=[O:3]"))], above=["PCC"])]
    doc = layout.build_scheme(steps, style="acs")
    frag_refs = [x for x in doc.refs if x.kind == "fragment"]
    assert len(frag_refs) == 3, "the intermediate must be drawn once"
    xml = doc.to_xml()
    element = {int(n.get("id")): int(n.get("Element", "6")) for n in ET.fromstring(xml).iter("n")}
    ids = doc.steps[1]["ReactionStepAtomMap"]
    assert 'ReactionStepAtomMap="%s"' % " ".join(str(i) for i in ids) in xml
    pairs = list(zip(ids[::2], ids[1::2]))
    assert len(pairs) == 3
    assert all(element[a] == element[b] for a, b in pairs), pairs
    assert {a for a, _ in pairs} == set(frag_refs[1].extra["atom_ids"])
    assert {b for _, b in pairs} == set(frag_refs[2].extra["atom_ids"])
    assert cdxml.read_cdxml(xml)["steps"][1]["atom_map_pairs"] == 3
