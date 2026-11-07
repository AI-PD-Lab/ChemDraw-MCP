"""Tests for electron-pushing arrows and mechanism layouts (no ChemDraw needed)."""
import math
import re
import xml.etree.ElementTree as ET

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402

from chemdrawmcp import cdxml, chem, mechanism  # noqa: E402
from chemdrawmcp.mechanism import AtomRef, ElectronArrow, parse_arrow_spec  # noqa: E402


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def dist(a, b):
    return math.hypot(a[0] - b[0], a[1] - b[1])


def box_dist(p, box):
    dx = max(box[0] - p[0], 0.0, p[0] - box[2])
    dy = max(box[1] - p[1], 0.0, p[1] - box[3])
    return math.hypot(dx, dy)


def electron_arrows(doc):
    return [r for r in doc.refs if r.kind == "arrow" and r.extra.get("mechanism")]


def fragments(doc, stage=None):
    out = [r for r in doc.refs if r.kind == "fragment"]
    if stage is not None:
        out = [r for r in out if r.extra.get("stage") == stage]
    return out


def label_box(ref, atom):
    for lab in ref.extra.get("labels") or []:
        if lab.get("atom") == atom:
            return lab["box"]
    return None


def assert_near_atom(point, ref, atom, L, lo=0.25, hi=0.6):
    """An arrow end sits 0.25-0.6 L from the referenced atom, or just outside the atom's label ("OH-", "Cl")."""
    p = mechanism.fragment_atom_positions(_DOC[0], ref)[atom]
    d = dist(point, p)
    box = label_box(ref, atom)
    ok = lo * L <= d <= hi * L
    if not ok and box is not None:
        ok = box_dist(point, box) <= 0.15 * L and d <= 1.2 * L
    assert ok, "end %s is %.2f L from atom %d (label box %s)" % (point, d / L, atom, box)


def assert_near_bond(point, ref, i, j, L):
    pos = mechanism.fragment_atom_positions(_DOC[0], ref)
    mid = ((pos[i][0] + pos[j][0]) / 2.0, (pos[i][1] + pos[j][1]) / 2.0)
    d = dist(point, mid)
    assert 0.2 * L <= d <= 0.5 * L, "end %s is %.2f L from the midpoint of bond %d-%d" % (point, d / L, i, j)
    assert dist(point, pos[i]) <= 0.85 * L and dist(point, pos[j]) <= 0.85 * L


def assert_bulges_away(ref, L):
    """Sweep 90-150 degrees, a clearly visible arc, apex further from the molecule centroid than the chord midpoint
    (when the centroid is clearly off the chord; otherwise obstacles decide the side)."""
    tail, head, apex, c = ref.extra["tail"], ref.extra["head"], ref.extra["apex"], ref.extra["centroid"]
    mid = ((tail[0] + head[0]) / 2.0, (tail[1] + head[1]) / 2.0)
    assert abs(ref.extra["angle"]) >= 90 - 1e-6 and abs(ref.extra["angle"]) <= 150 + 1e-6
    chord = dist(tail, head)
    assert dist(apex, mid) >= 0.2 * chord                       # sagitta of a >= 90 degree arc
    ux, uy = (head[0] - tail[0]) / chord, (head[1] - tail[1]) / chord
    vx, vy = mid[0] - c[0], mid[1] - c[1]
    off = abs(vx * (-uy) + vy * ux)          # distance of the centroid from the chord line
    if off > 0.5 * L:
        assert dist(apex, c) > dist(mid, c), "arc %s bulges towards the centroid" % ref.extra["spec"]


def arc_clearance(ref, point, L):
    """Smallest distance between the inner part of the arc and a point."""
    tail, head = tuple(ref.extra["tail"]), tuple(ref.extra["head"])
    pts = [p for p in mechanism.arc_points(tail, head, ref.extra["angle"], 40) if dist(p, tail) > 0.35 * L and dist(p, head) > 0.35 * L]
    return min(dist(p, point) for p in pts)


_DOC = [None]


def build(stages, **kw):
    doc = mechanism.build_mechanism(stages, **kw)
    _DOC[0] = doc
    return doc


# ---------------------------------------------------------------------------
# specs
# ---------------------------------------------------------------------------

def test_parse_compact_specs():
    a = parse_arrow_spec("lp:5->5-8")
    assert a.source == ("lp", (0, 5)) and a.target == ("bond", (0, 5), (0, 8)) and a.electrons == 2
    a = parse_arrow_spec("3-4->4")
    assert a.source == ("bond", (0, 3), (0, 4)) and a.target == ("atom", (0, 4))
    a = parse_arrow_spec("3=4->4-7")
    assert a.source_kind == "bond" and a.target_kind == "bond" and a.target_atoms == ((0, 4), (0, 7))
    assert parse_arrow_spec("lp:5->5-8:1e").electrons == 1
    assert parse_arrow_spec("lp:5~>5-8").electrons == 1
    assert parse_arrow_spec("lp:5 -> 5-8 (2e)").electrons == 2
    a = parse_arrow_spec("lp:0:0->0:0-1:0")
    assert a.source == ("lp", (0, 0)) and a.target == ("bond", (0, 0), (1, 0)) and a.fragments == [0, 1]
    a = parse_arrow_spec("h:3->3-4")
    assert a.source == ("h", (0, 3))
    a = parse_arrow_spec("3=4->new:4-7")
    assert a.target == ("new_bond", (0, 4), (0, 7))
    a = parse_arrow_spec("1:0-1:1->1:1")
    assert a.source == ("bond", (1, 0), (1, 1)) and a.target == ("atom", (1, 1)) and a.to_spec() == "1:0-1:1->1:1"


def test_parse_dict_and_object_specs():
    a = parse_arrow_spec({"from": "lp:5", "to": "5-8", "electrons": 1})
    assert a.source == ("lp", (0, 5)) and a.electrons == 1 and a.arrow_kind == "half_left"
    a = parse_arrow_spec({"from": {"lone_pair": 5}, "to": {"new_bond": [5, 8]}, "bulge": "up", "label": "attack"})
    assert a.target == ("new_bond", (0, 5), (0, 8)) and a.bulge == "up" and a.label == "attack"
    a = parse_arrow_spec({"source": [3, 4], "target": 4})
    assert a.source == ("bond", (0, 3), (0, 4)) and a.target == ("atom", (0, 4))
    a = parse_arrow_spec({"spec": "3-4->4", "electrons": 1})
    assert a.source_kind == "bond" and a.electrons == 1
    a = ElectronArrow(("lp", 5), ("new_bond", 5, 8))
    assert a.source == ("lp", (0, 5)) and a.target == ("new_bond", (0, 5), (0, 8))
    assert ElectronArrow((3, 4), 4).target == ("atom", (0, 4))
    assert ElectronArrow(("lp", AtomRef(1, 0)), ("new_bond", AtomRef(0, 0), AtomRef(1, 0))).fragments == [0, 1]
    assert ElectronArrow("lp:5", "5-8") == parse_arrow_spec("lp:5->5-8")
    assert parse_arrow_spec(("lp:5", "5-8", 1)).electrons == 1
    assert parse_arrow_spec(a) is a
    # round trip through the compact form
    for spec in ("lp:5->5-8", "3-4->4", "3=4->new:4-7:1e", "lp:0:0->0:0-1:0", "h:3->3-4"):
        b = parse_arrow_spec(spec)
        assert parse_arrow_spec(b.to_spec()) == b


def test_parse_spec_validation_against_molecule():
    mol = chem.mol_from_input("CC(=O)C")            # atoms: 0 C, 1 C, 2 O, 3 C
    a = parse_arrow_spec("lp:2->1-2", mol)
    assert a.target_kind == "bond"                   # C1-O2 exists -> electrons go into the bond
    a = parse_arrow_spec("lp:2->0-2", mol)
    assert a.target_kind == "new_bond"               # C0-O2 does not exist
    with pytest.raises(ValueError, match="atom 7"):
        parse_arrow_spec("lp:7->7-1", mol)
    with pytest.raises(ValueError, match="does not exist"):
        parse_arrow_spec("0-2->2", mol)              # no bond between C0 and O2
    with pytest.raises(ValueError, match="fragment 1"):
        parse_arrow_spec("lp:1:0->1:0-0:1", mol)
    mols = [chem.mol_from_input("[OH-]"), chem.mol_from_input("CCl")]
    a = parse_arrow_spec("lp:0:0->0:0-1:0", mols)
    assert a.target_kind == "new_bond"


@pytest.mark.parametrize("bad", ["lp:5", "5->", "->5", "lp:5->5-8:3e", "foo:5->5", "5-5->5", "lp:a->5"])
def test_bad_specs_raise(bad):
    with pytest.raises(ValueError):
        parse_arrow_spec(bad)


def test_bad_arrow_options_raise():
    with pytest.raises(ValueError):
        ElectronArrow("lp:5", "5-8", electrons=3)
    with pytest.raises(ValueError):
        ElectronArrow("lp:5", "5-8", bulge="left")
    with pytest.raises(ValueError):
        ElectronArrow("5-8", "lp:5")                 # a lone pair cannot be a target
    with pytest.raises(ValueError):
        ElectronArrow("new:5-8", "5")                # a new bond cannot be a source
    with pytest.raises(ValueError):
        ElectronArrow("lp:5", "5-8", angle=5)


# ---------------------------------------------------------------------------
# geometry on a single fragment
# ---------------------------------------------------------------------------

def test_fragment_positions_from_xml_fallback():
    doc = cdxml.CDXMLDocument("acs")
    mol = chem.mol_from_input("CC(=O)Oc1ccccc1C(=O)O")
    ref = doc.add_fragment(mol, 150, 120)
    recorded = ref.extra.pop("atom_positions", None)
    pos = mechanism.fragment_atom_positions(doc, ref)
    assert len(pos) == mol.GetNumAtoms()
    frag = ET.fromstring(next(o for o in doc.objects if o.startswith('<fragment id="%d"' % ref.id)))
    by_id = {n.get("id"): tuple(float(v) for v in n.get("p").split()) for n in frag.findall("n")}
    for i, nid in enumerate(ref.extra["atom_ids"]):
        assert dist(pos[i], by_id[str(nid)]) < 1e-6
    if recorded:
        for i, p in enumerate(recorded):
            assert dist(pos[i], p) < 0.02
    geom = mechanism.fragment_geometry(doc, ref)
    assert geom.has_bond(1, 2) and geom.bond_order(1, 2) == "2" and not geom.has_bond(0, 2)
    assert geom.elements[2] == 8 and geom.labelled[2] and not geom.labelled[0]
    # connectivity from the molecule when the XML is gone (e.g. a foreign document)
    ref2 = doc.add_fragment(mol, 300, 120)
    doc.objects = [o for o in doc.objects if not o.startswith('<fragment id="%d"' % ref2.id)]
    geom2 = mechanism.fragment_geometry(doc, ref2, mol=mol)
    assert geom2.has_bond(1, 2) and len(geom2.positions) == mol.GetNumAtoms()


def test_add_electron_arrows_single_fragment():
    doc = cdxml.CDXMLDocument("acs")
    _DOC[0] = doc
    L = doc.bond_length
    mol = chem.mol_from_input("C=C([O-])C")          # enolate: 0 C, 1 C, 2 O-, 3 C
    ref = doc.add_fragment(mol, 200, 200)
    refs = mechanism.add_electron_arrows(doc, ref, ["lp:2->1=2", ElectronArrow((0, 1), 0)], mol)
    assert len(refs) == 2 and all(r.kind == "arrow" and r.extra["mechanism"] for r in refs)
    lp, pi = refs
    assert lp.extra["target"]["kind"] == "bond" and lp.extra["source"]["kind"] == "lp"
    assert_near_atom(lp.extra["tail"], ref, 2, L)
    assert_near_bond(lp.extra["head"], ref, 1, 2, L)
    assert_near_bond(pi.extra["tail"], ref, 0, 1, L)
    assert_near_atom(pi.extra["head"], ref, 0, L)
    for r in refs:
        assert_bulges_away(r, L)
        assert r.box.width < 2.5 * L and r.box.height < 2.5 * L      # tight box, not the whole circle
        # the arc apex and the circle centre written to the XML lie on opposite sides of the chord
        tail, head, apex = r.extra["tail"], r.extra["head"], r.extra["apex"]
        mid = ((tail[0] + head[0]) / 2.0, (tail[1] + head[1]) / 2.0)
        el = ET.fromstring(next(o for o in doc.objects if o.startswith('<arrow id="%d"' % r.id)))
        cx, cy = (float(v) for v in el.get("Center3D").split()[:2])
        assert (apex[0] - mid[0]) * (cx - mid[0]) + (apex[1] - mid[1]) * (cy - mid[1]) < 0
        assert el.get("AngularSize") == cdxml._f(r.extra["angle"])
    xml = doc.to_xml()
    assert len(re.findall(r"<arrow[^>]*AngularSize=", xml)) == 2
    assert cdxml.read_cdxml(xml)["molecules"][0]["smiles"] == Chem.MolToSmiles(Chem.MolFromSmiles("C=C([O-])C"))
    # arrows keep their ends >= 0.3 L away from every other atom
    pos = mechanism.fragment_atom_positions(doc, ref)
    for r in refs:
        own = {tuple(a) for a in r.extra["source"]["atoms"] + r.extra["target"]["atoms"]}
        for end in (r.extra["tail"], r.extra["head"]):
            for i, p in enumerate(pos):
                if (0, i) not in own:
                    assert dist(end, p) >= 0.3 * L - 1e-6


def test_fishhook_and_bulge_override():
    doc = cdxml.CDXMLDocument("acs")
    mol = chem.mol_from_input("CCl")
    ref = doc.add_fragment(mol, 100, 100)
    up, down, fish = mechanism.add_electron_arrows(doc, ref, [
        ElectronArrow((0, 1), 1, bulge="up"), ElectronArrow((0, 1), 1, bulge="down"),
        ElectronArrow((0, 1), 0, electrons=1)], mol)
    assert up.extra["apex"][1] < min(up.extra["tail"][1], up.extra["head"][1]) + 1e-6 or up.extra["angle"] * (
        up.extra["head"][0] - up.extra["tail"][0]) < 0
    mid_up = (up.extra["tail"][1] + up.extra["head"][1]) / 2.0
    mid_down = (down.extra["tail"][1] + down.extra["head"][1]) / 2.0
    assert up.extra["apex"][1] < mid_up and down.extra["apex"][1] > mid_down
    assert fish.extra["arrow_kind"] == "half_left" and fish.extra["electrons"] == 1
    xml = doc.to_xml()
    el = ET.fromstring(next(o for o in doc.objects if o.startswith('<arrow id="%d"' % fish.id)))
    assert el.get("ArrowheadHead") == "HalfLeft" and el.get("AngularSize")
    assert len(re.findall(r'ArrowheadHead="HalfLeft"', xml)) == 1
    forced = mechanism.add_electron_arrows(doc, ref, [ElectronArrow((0, 1), 1, angle=-75)], mol)[0]
    assert forced.extra["angle"] == -75


def test_arrow_label_and_obstacles():
    doc = cdxml.CDXMLDocument("acs")
    mol = chem.mol_from_input("CC(=O)C")
    ref = doc.add_fragment(mol, 100, 100)
    plus = doc.add_plus(140, 100)
    refs = mechanism.add_electron_arrows(doc, ref, [{"from": "1=2", "to": "2", "label": "pi -> O"}], mol, obstacles=[plus])
    assert refs[0].extra["label"] == "pi -> O" and refs[0].extra.get("label_id")
    assert any(r.kind == "text" and r.extra.get("text") == "pi -> O" for r in doc.refs)


# ---------------------------------------------------------------------------
# mechanisms
# ---------------------------------------------------------------------------

def test_sn2_two_fragments_per_stage():
    doc = build([
        {"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"], "label": "SN2"},
        {"structures": ["CO", "[Cl-]"], "label": "products"},
    ], title="SN2")
    L = doc.bond_length
    hydroxide, chloromethane = fragments(doc, 0)
    assert hydroxide.extra["smiles"] == "[OH-]" and chloromethane.extra["smiles"] == "CCl"
    arrows = electron_arrows(doc)
    assert [a.extra["spec"] for a in arrows] == ["lp:0:0->new:0:0-1:0", "1:0-1:1->1:1"]
    lp, cl = arrows
    # lone pair tail beside the hydroxide oxygen (~0.5 L out, not on the charge superscript at the label's top right),
    # head close to the carbon being attacked (new bond), just outside its "CH3" label
    assert_near_atom(lp.extra["tail"], hydroxide, 0, L)
    o_pos = mechanism.fragment_atom_positions(doc, hydroxide)[0]
    assert dist(lp.extra["tail"], o_pos) >= 0.45 * L
    ob = label_box(hydroxide, 0)
    assert not (lp.extra["tail"][0] > ob[2] - 0.2 * L and lp.extra["tail"][1] < (ob[1] + ob[3]) / 2.0), "tail covers the charge"
    assert_near_atom(lp.extra["head"], chloromethane, 0, L, hi=0.65)
    c_pos = mechanism.fragment_atom_positions(doc, chloromethane)[0]
    assert lp.extra["head"][0] < c_pos[0]                        # approaches from the hydroxide side
    assert lp.extra["target"]["kind"] == "new_bond" and lp.extra["fragment_ids"] == [hydroxide.id, chloromethane.id]
    # C-Cl bond tail beside the bond (>= 0.3 L off the bond line), head beside the chlorine, clearly visible hook
    assert_near_bond(cl.extra["tail"], chloromethane, 0, 1, L)
    assert abs(cl.extra["tail"][1] - mechanism.fragment_atom_positions(doc, chloromethane)[0][1]) >= 0.3 * L
    assert_near_atom(cl.extra["head"], chloromethane, 1, L)
    assert dist(cl.extra["tail"], cl.extra["head"]) >= 0.8 * L
    for a in arrows:
        assert_bulges_away(a, L)
        mid = ((a.extra["tail"][0] + a.extra["head"][0]) / 2.0, (a.extra["tail"][1] + a.extra["head"][1]) / 2.0)
        assert dist(a.extra["apex"], mid) >= 0.25 * L
    # the lone-pair arc keeps clear of the '+' sign instead of passing through it
    plus = next(r for r in doc.refs if r.kind == "text" and r.extra.get("text") == "+")
    assert arc_clearance(lp, (plus.box.cx, plus.box.cy), L) >= 0.3 * L
    # document: 2 curved arrows + 1 reaction arrow, readable, one reaction step
    xml = doc.to_xml()
    assert xml.count("<arrow ") == 3
    assert len(re.findall(r"<arrow[^>]*AngularSize=", xml)) == 2
    back = cdxml.read_cdxml(xml)
    assert [m["smiles"] for m in back["molecules"]] == ["[OH-]", "CCl", "CO", "[Cl-]"]
    assert len(back["arrows"]) == 3 and len(back["steps"]) == 1
    assert back["reaction_smiles"] == ["[OH-].CCl>>CO.[Cl-]"]
    assert "SN2" in cdxml.document_text(xml) and "products" in cdxml.document_text(xml)
    info = mechanism.describe_mechanism(doc)
    assert info["title"] == "SN2" and info["electron_arrows"] == 2 and len(info["stages"]) == 2
    assert info["stages"][0]["arrows"][0]["spec"] == "lp:0:0->new:0:0-1:0"
    assert info["stages"][0]["arrow_to_next"] == "full" and info["stages"][1]["arrow_to_next"] is None
    assert info["reaction_arrows"][0]["from_stage"] == 0 and info["reaction_arrows"][0]["to_stage"] == 1
    assert getattr(doc, "warnings", []) == []


def test_carbonyl_addition():
    doc = build([
        {"structures": ["[OH-]", "CC(=O)C"], "arrows": ["lp:0:0->0:0-1:1", "1:1=1:2->1:2"], "arrow": "equilibrium"},
        {"structure": "CC(C)(O)[O-]", "label": "tetrahedral intermediate"},
    ])
    L = doc.bond_length
    hydroxide, acetone = fragments(doc, 0)
    lp, pi = electron_arrows(doc)
    assert_near_atom(lp.extra["tail"], hydroxide, 0, L)
    c1 = mechanism.fragment_atom_positions(doc, acetone)[1]
    assert 0.35 * L <= dist(lp.extra["head"], c1) <= 0.65 * L
    assert pi.extra["source"]["kind"] == "bond" and pi.extra["target"]["kind"] == "atom"
    assert_near_bond(pi.extra["tail"], acetone, 1, 2, L)
    assert_near_atom(pi.extra["head"], acetone, 2, L)
    assert_bulges_away(pi, L)
    assert_bulges_away(lp, L)
    xml = doc.to_xml()
    assert len(re.findall(r"<arrow[^>]*AngularSize=", xml)) == 2
    assert 'ArrowShaftSpacing="300"' in xml                       # equilibrium arrow between the stages
    back = cdxml.read_cdxml(xml)
    assert len(back["molecules"]) == 3 and len(back["arrows"]) == 3


def test_e1_cation_with_explicit_hydrogens():
    doc = build([
        {"structures": ["[H]C([H])([H])[C+](C)C", "O"], "arrows": ["lp:1:0->1:0-0:0", "0:0-0:1->0:1-0:4"], "label": "E1"},
        {"structures": ["C=C(C)C", "[OH3+]"]},
    ])
    L = doc.bond_length
    cation, water = fragments(doc, 0)
    assert cation.extra["atoms"] == 7                                  # explicit H atoms are kept
    frag_xml = next(o for o in doc.objects if o.startswith('<fragment id="%d"' % cation.id))
    assert frag_xml.count('Element="1"') == 3
    lp, ch = electron_arrows(doc)
    assert_near_atom(lp.extra["tail"], water, 0, L)
    h0 = mechanism.fragment_atom_positions(doc, cation)[0]
    assert dist(lp.extra["head"], h0) <= 0.65 * L or box_dist(lp.extra["head"], label_box(cation, 0)) <= 0.15 * L
    assert_near_bond(ch.extra["tail"], cation, 0, 1, L)
    assert ch.extra["target"]["kind"] == "bond"                        # C1-C4 exists: it becomes the pi bond
    assert_near_bond(ch.extra["head"], cation, 1, 4, L)
    pos = mechanism.fragment_atom_positions(doc, cation)
    # the abstracted hydrogen faces the base after orientation: it is on the water's side of the cation
    assert h0[0] >= sum(p[0] for p in pos) / len(pos)
    for a in (lp, ch):
        assert_bulges_away(a, L)
    xml = doc.to_xml()
    back = cdxml.read_cdxml(xml)
    assert len(back["molecules"]) == 4 and len([a for a in back["arrows"] if a.get("angularSize")]) == 2


def test_resonance_wrapping_title_and_fishhook_stage():
    stages = [
        {"structure": "C=C([O-])C", "arrows": ["lp:2->1=2", "0=1->0"], "arrow": "resonance", "label": "enolate"},
        {"structure": "[CH2-]C(=O)C", "label": "carbanion"},
        {"structures": ["c1ccccc1", "[N+](=O)=O"], "arrows": ["0:0=0:1->0:0-1:0", "1:0=1:1->1:1"], "above": ["H2SO4"]},
        {"structure": "[O-][N+](=O)C1C=CC=C[CH+]1", "arrows": ["h:3->3-4"]},
        {"structure": "C[CH2]", "arrows": ["lp:1->1-0:1e", "0-1->0:1e"], "label": "radical"},
    ]
    doc = build(stages, max_width=420, title="assorted")
    L = doc.bond_length
    rows = sorted({round(r.box.cy) for r in fragments(doc)})
    assert len(rows) >= 2                                               # wrapped into at least two rows
    assert rows[-1] - rows[0] > 2 * L
    arrows = electron_arrows(doc)
    assert len(arrows) == 7
    for a in arrows:
        assert_bulges_away(a, L)
        assert 90 - 1e-6 <= abs(a.extra["angle"]) <= 150 + 1e-6
    assert [a.extra["electrons"] for a in arrows[-2:]] == [1, 1]
    xml = doc.to_xml()
    assert len(re.findall(r"<arrow[^>]*AngularSize=", xml)) == 7
    assert len(re.findall(r'ArrowheadHead="HalfLeft"', xml)) == 2
    # resonance arrow: double headed, no reaction step recorded for it; 3 other transitions -> 3 steps
    assert len(re.findall(r'<arrow[^>]*ArrowheadTail="Full"', xml)) == 1
    assert xml.count("<step ") == 3
    back = cdxml.read_cdxml(xml)
    assert len(back["molecules"]) == 6 and len(back["arrows"]) == 11
    assert "assorted" in cdxml.document_text(xml) and "H2SO4" in cdxml.document_text(xml)
    info = mechanism.describe_mechanism(doc)
    assert [s["row"] for s in info["stages"]][0] == 0 and info["stages"][-1]["row"] >= 1
    assert info["stages"][0]["arrow_to_next"] == "resonance"
    assert info["stages"][2]["above"] == ["H2SO4"]
    assert info["stages"][-1]["arrows"][0]["kind"] == "half_left"


def test_single_atom_reagents_have_finite_coordinates():
    """[H-] / [H+] come out of the depictor with NaN coordinates: they must still lay out and take arrows."""
    for smi in ("[H-]", "[H+]", "[OH-]", "[Cl-]", "[Na+]"):
        m = mechanism._coerce_mol(smi)
        p = m.GetConformer().GetAtomPosition(0)
        assert m.GetNumAtoms() == 1 and math.isfinite(p.x) and math.isfinite(p.y)
    doc = build([
        {"structures": ["CC(=O)C", "[H-]"], "arrows": ["lp:1:0->1:0-0:1", "0:1=0:2->0:2"], "arrow": "full"},
        {"structure": "CC([O-])C", "label": "alkoxide", "arrow": "full", "above": ["H3O+"]},
        {"structure": "CC(O)C", "label": "isopropanol"},
    ], title="Hydride addition to acetone")
    L = doc.bond_length
    xml = doc.to_xml()
    assert re.search(r"(?i)(?<![a-z])nan(?![a-z])", xml) is None
    acetone, hydride = fragments(doc, 0)
    lp, pi = electron_arrows(doc)
    assert_near_atom(lp.extra["tail"], hydride, 0, L)
    assert_near_atom(lp.extra["head"], acetone, 1, L, hi=0.65)
    assert lp.extra["centroid"] != list(mechanism.fragment_atom_positions(doc, hydride)[0])   # union centroid
    assert_near_bond(pi.extra["tail"], acetone, 1, 2, L)
    assert_near_atom(pi.extra["head"], acetone, 2, L)
    for a in (lp, pi):
        assert_bulges_away(a, L)
    back = cdxml.read_cdxml(xml)
    assert [m["smiles"] for m in back["molecules"]][:2] == ["CC(C)=O", "[H-]"] and len(back["arrows"]) == 4
    # single atoms as source and as target, partner on either side, two single atoms, lone single-atom stages
    cases = [
        [{"structures": ["[Cl-]", "C[C+](C)C"], "arrows": ["lp:0:0->0:0-1:1"]}, {"structure": "CC(C)(C)Cl"}],
        [{"structures": ["CC([O-])C", "[H+]"], "arrows": ["lp:0:2->0:2-1:0"]}, {"structure": "CC(O)C"}],
        [{"structures": ["[H+]", "[OH-]"], "arrows": ["lp:1:0->1:0-0:0"]}, {"structure": "O"}],
        [{"structures": ["[Na+]", "[Cl-]"], "arrows": ["lp:1:0->1:0-0:0:1e", "lp:0:0->0:0-1:0"]}],
        [{"structures": ["[H-]", "[H+]"], "arrows": ["lp:0:0->0:0-1:0"]}, {"structure": "[H][H]"}],
        [{"structure": "[Cl-]", "label": "lone ion"}],
        [{"structure": "[Cl-]", "arrows": ["lp:0->0"]}],
    ]
    for stages in cases:
        doc = build(stages)
        xml = doc.to_xml()
        assert re.search(r"(?i)(?<![a-z])nan(?![a-z])", xml) is None
        n_arrows = sum(len(st.get("arrows", [])) for st in stages)
        assert len(electron_arrows(doc)) == n_arrows
        for a in electron_arrows(doc):
            assert all(math.isfinite(v) for v in a.extra["tail"] + a.extra["head"] + a.extra["apex"])
            assert dist(a.extra["tail"], a.extra["head"]) >= 0.3 * doc.bond_length
            for key in a.extra["source"]["atoms"] + a.extra["target"]["atoms"]:
                frag = next(r for r in fragments(doc, 0) if r.extra["fragment_index"] == key[0])
                if len(a.extra["source"]["atoms"]) == 1 and key == a.extra["source"]["atoms"][0]:
                    assert_near_atom(a.extra["tail"], frag, key[1], doc.bond_length)
        back = cdxml.read_cdxml(xml)
        assert len(back["molecules"]) == sum(len(st["structures"]) if "structures" in st else 1 for st in stages)


def test_tiny_structures_get_atom_labels():
    doc = build([
        {"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"]},
        {"structures": ["CO", "[Cl-]"]},
        {"structures": ["O=C=O", "CC(=O)C"]},
        {"structures": [{"smiles": "CCl", "label_atoms": False}], "label_atoms": True},
    ])
    def labelled_atoms(ref):
        frag = next(o for o in doc.objects if o.startswith('<fragment id="%d"' % ref.id))
        return len(re.findall(r"<n [^>]*><t ", frag))
    hydroxide, chloromethane = fragments(doc, 0)
    methanol, chloride = fragments(doc, 1)
    co2, acetone = fragments(doc, 2)
    (plain,) = fragments(doc, 3)
    assert labelled_atoms(chloromethane) == 2 and labelled_atoms(methanol) == 2      # CH3-Cl, CH3-OH
    assert labelled_atoms(co2) == 3                                                   # cumulated carbon
    assert labelled_atoms(acetone) == 0                                               # normal skeletal drawing
    assert labelled_atoms(hydroxide) == 1 and labelled_atoms(chloride) == 1            # single ions are tiny too
    assert labelled_atoms(plain) == 0                                                 # explicit override wins
    assert [lab["text"] for lab in chloromethane.extra["labels"]] == ["CH3", "Cl"]
    # atom indices are unaffected and arrow ends stay outside the (wider) labels
    lp, cl = electron_arrows(doc)[:2]
    assert lp.extra["target"]["atoms"] == [[0, 0], [1, 0]]
    assert box_dist(lp.extra["head"], label_box(chloromethane, 0)) >= -1e-6
    assert box_dist(cl.extra["head"], label_box(chloromethane, 1)) >= -1e-6
    assert cdxml.read_cdxml(doc.to_xml())["molecules"][1]["smiles"] == "CCl"


def test_stage_inputs_and_errors():
    mol = chem.mol_from_input("CC(=O)C")
    st = mechanism.Stage.coerce({"structure": mol, "arrows": ["1=2->2"], "label": "x"})
    assert len(st.mols) == 1 and st.arrows[0].target_kind == "atom"
    st = mechanism.Stage.coerce("c1ccccc1")
    assert st.arrows == [] and st.mols[0].GetNumAtoms() == 6
    st = mechanism.Stage.coerce({"structures": [{"smiles": "[OH-]", "label": "base"}, "CCl"], "arrows": []})
    assert st.labels == ["base", None]
    with pytest.raises(ValueError, match="atom 6"):
        mechanism.build_mechanism([{"structure": "c1ccccc1", "arrows": ["0=1->0-6"]}])
    with pytest.raises(ValueError, match="structure"):
        mechanism.build_mechanism([{"arrows": ["0=1->0"]}])
    with pytest.raises(ValueError):
        mechanism.build_mechanism([])
    doc = mechanism.build_mechanism([{"structure": "CC(=O)C", "arrows": ["lp:2->1-2"]}], style="nature", bond_length=20)
    assert doc.bond_length == 20 and len(electron_arrows(doc)) == 1
    assert cdxml.read_cdxml(doc.to_xml())["molecules"][0]["smiles"] == "CC(C)=O"


def test_describe_without_build_metadata():
    doc = cdxml.CDXMLDocument("acs")
    mol = chem.mol_from_input("CC(=O)C")
    ref = doc.add_fragment(mol, 100, 100)
    mechanism.add_electron_arrows(doc, ref, ["1=2->2"], mol)
    info = mechanism.describe_mechanism(doc)
    assert info["electron_arrows"] == 1 and len(info["stages"]) == 1
    assert info["stages"][0]["fragment_ids"] == [ref.id] and info["stages"][0]["arrows"][0]["spec"] == "1-2->2"
    assert info["reaction_arrows"] == [] and info["bounds"]


def test_same_bond_source_and_target_rejected_with_pointer_to_working_forms():
    """'3-4->3=4' (the bond's own electrons turning it into a double bond) is meaningless: the error says so and
    names the forms that work, which all draw a hook landing on the existing bond."""
    mol = chem.mol_from_input("CCC=CC")
    with pytest.raises(ValueError, match=r"same bond 3-4.*lp:3->3=4.*h:3->3=4"):
        mechanism.build_mechanism([{"structure": "CCC=CC", "arrows": ["3-4->3=4"]}])
    with pytest.raises(ValueError, match="same bond"):
        parse_arrow_spec("3-4->3=4", mol)
    with pytest.raises(ValueError, match="same bond"):
        parse_arrow_spec("4-3->new:3-4", mol)
    assert parse_arrow_spec("3-4->3=4").target_kind == "bond"          # parsing alone (no molecule) is lenient
    for smi, spec in (("CC([O-])C", "lp:2->1=2"), ("CCC=CC", "h:1->1=2"), ("CCC=CC", "0-1->1=2")):
        doc = build([{"structure": smi, "arrows": [spec]}])
        (frag,) = fragments(doc, 0)
        (arrow,) = electron_arrows(doc)
        assert arrow.extra["target"] == {"kind": "bond", "atoms": [[0, 1], [0, 2]]}
        assert_near_bond(arrow.extra["head"], frag, 1, 2, doc.bond_length)
        assert doc.warnings == []


def test_rows_and_title_clear_reagent_stacks_and_captions():
    """Row baselines and the title leave room for reagent lines stacked above an arrow and for captions."""
    def texts(doc, prefix):
        return [r for r in doc.refs if r.kind == "text" and str(r.extra.get("text", "")).startswith(prefix)]

    williamson = [{"structure": "CCO", "arrows": ["lp:2->2-1"], "above": ["1. NaH, THF", "2. MeI", "3. aq. NH4Cl"]},
                  {"structure": "CCOC"}]
    doc = build(williamson, title="Williamson ether synthesis")
    assert doc.warnings == []
    (title,) = texts(doc, "Williamson")
    reagents = texts(doc, "1.") + texts(doc, "2.") + texts(doc, "3.")
    assert len(reagents) == 3 and title.box.y0 >= 18.0
    assert title.box.y1 < min(t.box.y0 for t in reagents)
    # without a title the tallest stack is kept below the top margin instead of running off the page
    doc = build([{"structure": "CCO", "arrows": ["lp:2->2-1"], "above": ["reagent %d" % i for i in range(1, 6)]},
                 {"structure": "CCOC"}])
    assert doc.warnings == [] and min(r.box.y0 for r in doc.refs) >= 18.0
    # a stage with a tall reagent stack opening the second row keeps clear of the captions of the first row
    steroid = "CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O"
    doc = build([
        {"structure": steroid, "label": "steroid 1 with a long caption"},
        {"structure": steroid, "label": "steroid 2 with a long caption"},
        {"structure": "CCO", "arrows": ["lp:2->2-1"], "above": ["reagent %d" % i for i in range(1, 6)]},
        {"structure": "CCOC"}], max_width=330)
    assert doc.warnings == []
    assert [m["row"] for m in doc.mechanism_stages] == [0, 0, 1, 1]
    assert max(c.box.y1 for c in texts(doc, "steroid")) < min(t.box.y0 for t in texts(doc, "reagent"))
    # ordinary stages are unaffected: the first row still sits on the requested baseline
    doc = build([{"structure": "CCO", "arrows": ["lp:2->2-1"]}, {"structure": "CCOC"}])
    assert [round(r.box.cy, 6) for r in fragments(doc)] == [72.0, 72.0]
