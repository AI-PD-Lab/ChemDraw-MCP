"""Layout quality: label-aware bounds, '+' and arrow placement, grids, text metrics, overlap detection."""
import os
import timeit

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import rdDepictor  # noqa: E402

from chemdrawmcp import cdxml, layout, textmetrics  # noqa: E402
from chemdrawmcp.cdxml import Box  # noqa: E402


def prep(smi):
    m = Chem.MolFromSmiles(smi)
    Chem.Kekulize(m, clearAromaticFlags=True)
    rdDepictor.Compute2DCoords(m)
    return m


def comp(smi, label=None, **kw):
    return layout.Component(prep(smi), label, **kw)


def frags(doc):
    return [r for r in doc.refs if r.kind == "fragment"]


def texts(doc, exclude_plus=True):
    return [r for r in doc.refs if r.kind == "text" and not (exclude_plus and r.extra["text"] == "+")]


def plusses(doc):
    return [r for r in doc.refs if r.kind == "text" and r.extra["text"] == "+"]


def overlap(a: Box, b: Box, tol=0.01) -> bool:
    return a.x0 < b.x1 - tol and b.x0 < a.x1 - tol and a.y0 < b.y1 - tol and b.y0 < a.y1 - tol


# ---------------------------------------------------------------------------
# reactions
# ---------------------------------------------------------------------------

def test_reaction_with_terminal_heteroatoms_has_no_overlaps():
    step = layout.Step([comp("[NH3+]CC([O-])=O"), comp("CC(=O)O")], [comp("CCO")],
                       above=["NaBH4, MeOH"], below=["0 °C to rt, 2 h"])
    doc = layout.build_reaction(step, style="acs")
    assert layout.collisions(doc) == []
    assert doc.warnings == []
    # the ammonium label "H3N+" extends to the left of the first atom: the fragment box must include it
    gly = frags(doc)[0]
    labels = gly.extra["labels"]
    amm = [lab for lab in labels if lab["atom"] == 0][0]
    assert amm["text"].startswith("H3N") and amm["side"] == "left"
    lab_box = Box(*amm["box"])
    assert gly.box.x0 <= lab_box.x0 + 1e-6 and gly.box.x1 >= lab_box.x1 - 1e-6
    assert lab_box.x0 < gly.extra["atom_positions"][0][0] - 0.5 * doc.label_size   # sticks out to the left


def test_plus_sits_between_components_with_clearance():
    step = layout.Step([comp("[NH3+]CC([O-])=O"), comp("O=C(O)c1ccccc1"), comp("OCCO")],
                       [comp("CCOC(C)=O", label="92%"), comp("O")])
    doc = layout.build_reaction(step, style="acs")
    L = doc.bond_length
    y = 72.0
    fr = frags(doc)
    pl = plusses(doc)
    assert len(pl) == 3
    arrow = [r for r in doc.refs if r.kind == "arrow"][0]
    # every plus is strictly between its neighbours with >= 0.5 L clearance, and vertically centred on the baseline
    for p in pl:
        left = max((f for f in fr if f.box.x1 <= p.box.x0), key=lambda f: f.box.x1)
        right = min((f for f in fr if f.box.x0 >= p.box.x1), key=lambda f: f.box.x0)
        assert p.box.x0 - left.box.x1 >= 0.5 * L - 1e-6
        assert right.box.x0 - p.box.x1 >= 0.5 * L - 1e-6
        assert abs(p.box.cy - y) < 0.15 * (doc.caption_size + 2)
        assert not overlap(p.box, arrow.box)
    assert all(abs(f.box.cy - y) < 1e-6 for f in fr)
    assert layout.collisions(doc) == []


def test_arrow_is_longer_than_its_text_and_text_clears_the_line():
    step = layout.Step([comp("c1ccccc1Br")], [comp("c1ccccc1C(=O)O")],
                       above=["1. Mg, THF, reflux", "2. CO2 (dry ice)", "3. H3O+"],
                       below=["Pd(PPh3)4 (5 mol%), K2CO3", "dioxane / H2O, 100 °C, 12 h"])
    doc = layout.build_reaction(step, style="acs")
    L = doc.bond_length
    arrow = [r for r in doc.refs if r.kind == "arrow"][0]
    tail, head = arrow.extra["tail"], arrow.extra["head"]
    length = head[0] - tail[0]
    y = tail[1]
    above = [t for t in texts(doc) if t.box.y1 <= y]
    below = [t for t in texts(doc) if t.box.y0 >= y]
    assert len(above) == 3 and len(below) == 2
    widest = max(t.box.width for t in above + below)
    assert length >= widest + 0.8 * L - 0.02           # tail/head are rounded to 1/100 pt
    assert arrow.box.width - 0.5 * L >= widest + 0.8 * L - 1e-6   # unrounded: box = line padded by 0.25 L
    for t in above + below:
        assert t.box.x0 >= tail[0] - 1e-6 and t.box.x1 <= head[0] + 1e-6
        assert abs(t.box.cx - (tail[0] + head[0]) / 2.0) < 0.02      # tail/head are rounded to 1/100 pt
        assert min(abs(t.box.y0 - y), abs(t.box.y1 - y)) >= 0.35 * L - 1e-6
    # stacking order: first reagent line on top, lines do not overlap, conditions go downward
    above.sort(key=lambda t: t.box.y0)
    assert [t.extra["text"] for t in above] == ["1. Mg, THF, reflux", "2. CO2 (dry ice)", "3. H3O+"]
    for a, b in zip(above, above[1:]):
        assert a.box.y1 <= b.box.y0 + 1e-6
    below.sort(key=lambda t: t.box.y0)
    assert [t.extra["text"] for t in below] == ["Pd(PPh3)4 (5 mol%), K2CO3", "dioxane / H2O, 100 °C, 12 h"]
    assert below[0].box.y1 <= below[1].box.y0 + 1e-6
    assert layout.collisions(doc) == []


def test_reaction_with_captions_keeps_captions_under_their_structure():
    step = layout.Step([comp("O", label="water (a very long caption indeed)"), comp("CC(=O)O", label="2")],
                       [comp("CCO", label="3 (85%)")])
    doc = layout.build_reaction(step, style="acs")
    for f in frags(doc):
        cap = [t for t in texts(doc) if abs(t.box.cx - f.box.cx) < 1e-6 and t.box.y0 >= f.box.y1]
        assert len(cap) == 1
    assert layout.collisions(doc) == []


# ---------------------------------------------------------------------------
# grid and scheme
# ---------------------------------------------------------------------------

def test_grid_cells_do_not_overlap_and_are_centred():
    smiles = ["CCO", "[NH3+]CC([O-])=O", "c1ccccc1", "O", "CC(C)(C)OC(=O)N[C@@H](C)C(=O)O", "C1CCCCC1", "OS(=O)(=O)O"]
    labels = ["ethanol", "glycine zwitterion (1)", "benzene", "water", "Boc-Ala-OH, a long caption", None, "sulfuric acid"]
    doc = layout.build_structures([prep(s) for s in smiles], labels=labels, columns=3, style="acs")
    assert layout.collisions(doc) == [] and doc.warnings == []
    fr = frags(doc)
    assert len(fr) == 7
    # equal column widths: fragments in one column share the centre x; rows share the centre y
    for i in range(7):
        r, c = divmod(i, 3)
        assert abs(fr[i].box.cx - fr[c].box.cx) < 1e-6
        assert abs(fr[i].box.cy - fr[r * 3].box.cy) < 1e-6
    # labels centred under their structure: the nearest fragment above each caption shares its centre x
    caps = texts(doc)
    assert len(caps) == 6
    for t in caps:
        same_col = [f for f in fr if abs(f.box.cx - t.box.cx) < 1e-6 and f.box.y1 <= t.box.y0 + 1e-6]
        assert same_col, "caption %r is not centred under any structure" % t.extra["text"]
        owner = max(same_col, key=lambda f: f.box.y1)
        assert t.box.y0 - owner.box.y1 == pytest.approx(layout.CAPTION_GAP * doc.bond_length)
    # a long caption must not run into the neighbouring column
    for t in caps:
        for f in fr:
            assert not overlap(t.box, f.box)


def test_scheme_rows_wrap_without_overlaps_using_label_widths():
    steps = [layout.Step([comp("c1ccccc1", label="1"), comp("BrBr")], [comp("c1ccccc1Br", label="2 (bromobenzene, 95%)")],
                         above=["FeBr3"], below=["CH2Cl2, 0 °C"]),
             layout.Step([comp("c1ccccc1Br")], [comp("c1ccccc1C(=O)O", label="3")], above=["1. Mg, THF", "2. CO2", "3. H3O+"]),
             layout.Step([comp("c1ccccc1C(=O)O"), comp("CO")], [comp("c1ccccc1C(=O)OC", label="4 (methyl benzoate)")],
                         above=["H2SO4 (cat.)"], below=["reflux, 4 h"]),
             layout.Step([comp("c1ccccc1C(=O)OC")], [comp("[NH3+]CC([O-])=O", label="5")], above=["NH3, MeOH"])]
    doc = layout.build_scheme(steps, style="acs", max_width=330, title="Scheme 1")
    assert layout.collisions(doc) == []
    assert not [w for w in doc.warnings if "overlaps" in w]
    rows = sorted({round(f.box.cy) for f in frags(doc)})
    assert len(rows) >= 2
    # rows are separated: every object of a lower row starts below everything of the row above
    objs = [r for r in doc.refs if r.kind in ("fragment", "text")]
    for a in objs:
        for b in objs:
            if a.box.cy < b.box.cy - 3 * doc.bond_length:
                assert a.box.y1 <= b.box.y0 + 1e-6
    # the title sits above everything
    title = [t for t in texts(doc) if t.extra["text"] == "Scheme 1"][0]
    assert all(title.box.y1 <= o.box.y0 + 1e-6 for o in objs if o is not title)


def test_atom_positions_and_labels_are_recorded():
    doc = cdxml.CDXMLDocument("acs")
    mol = prep("[NH3+]CC([O-])=O")
    ref = doc.add_fragment(mol, 100, 100)
    pos = ref.extra["atom_positions"]
    assert len(pos) == mol.GetNumAtoms() == ref.extra["atoms"]
    for x, y in pos:
        assert ref.box.x0 <= x <= ref.box.x1 and ref.box.y0 <= y <= ref.box.y1
    labs = ref.extra["labels"]
    assert {lab["atom"] for lab in labs} == {0, 3, 4}
    for lab in labs:
        b = Box(*lab["box"])
        ax, ay = pos[lab["atom"]]
        assert b.x0 <= ax <= b.x1 and b.y0 <= ay <= b.y1
        tol = 0.01                                   # label boxes are rounded to 1/100 pt
        assert ref.box.x0 <= b.x0 + tol and b.x1 <= ref.box.x1 + tol and ref.box.y0 <= b.y0 + tol and b.y1 <= ref.box.y1 + tol
    by_atom = {lab["atom"]: lab for lab in labs}
    assert by_atom[3]["text"] in ("O-", "-O") and by_atom[4]["text"] == "O"
    # the atom-position list is JSON friendly (plain floats) and in page coordinates
    assert all(isinstance(v, float) for xy in pos for v in xy)
    # the summary stays serialisable
    import json
    json.dumps(doc.summary())


def test_label_side_follows_neighbour_direction():
    # methanol drawn with the O on the right extends right ("OH"); mirrored it extends left ("HO")
    doc = cdxml.CDXMLDocument("acs")
    m = Chem.MolFromSmiles("CO")
    conf = Chem.Conformer(2)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (1.5, 0.0, 0.0))
    m.AddConformer(conf)
    right = doc.add_fragment(m, 100, 100).extra["labels"][0]
    conf.SetAtomPosition(1, (-1.5, 0.0, 0.0))
    m.RemoveAllConformers()
    m.AddConformer(conf)
    left = doc.add_fragment(m, 100, 100).extra["labels"][0]
    assert right["side"] == "right" and right["text"] == "OH"
    assert left["side"] == "left" and left["text"] == "HO"
    # isolated atoms are centred
    water = doc.add_fragment(prep("O"), 100, 100).extra["labels"][0]
    assert water["side"] == "center"
    b = Box(*water["box"])
    assert abs((b.x0 + b.x1) / 2.0 - 100) < 1e-6


# ---------------------------------------------------------------------------
# text metrics
# ---------------------------------------------------------------------------

def test_text_metrics_are_monotonic_and_scale():
    assert textmetrics.text_width("ab", 10) > textmetrics.text_width("a", 10) > 0
    assert textmetrics.text_width("W", 10) > textmetrics.text_width("i", 10)
    assert textmetrics.text_width("reflux, 4 h", 10) > textmetrics.text_width("reflux", 10)
    assert abs(textmetrics.text_width("NaBH4", 20) - 2 * textmetrics.text_width("NaBH4", 10)) < 1e-6
    w, asc, desc = cdxml.text_box("one\ntwo", 10)
    assert asc == pytest.approx(9.0) and desc == pytest.approx(2.5 + 12.0)
    assert w == pytest.approx(max(textmetrics.text_width("one", 10), textmetrics.text_width("two", 10)))
    assert cdxml.text_width("", 10) == 0.0


def test_formula_text_is_narrower_than_plain_text():
    plain = textmetrics.text_width("NaBH4", 10)
    formula = textmetrics.text_width("NaBH4", 10, formula=True)
    assert formula < plain
    assert textmetrics.formula_runs("NaBH4") == (("NaBH", "normal"), ("4", "sub"))
    assert textmetrics.formula_runs("H2SO4 (cat.)") == (("H", "normal"), ("2", "sub"), ("SO", "normal"), ("4", "sub"), (" (cat.)", "normal"))
    assert textmetrics.formula_runs("3. H3O+") == (("3. H", "normal"), ("3", "sub"), ("O", "normal"), ("+", "sup"))
    assert textmetrics.formula_runs("1. Mg") == (("1. Mg", "normal"),)
    assert textmetrics.formula_runs("n-BuLi, -78 °C") == (("n-BuLi, -78 °C", "normal"),)
    # add_text passes the flag through: a formula caption gets a narrower box
    doc = cdxml.CDXMLDocument("acs")
    a = doc.add_text("NaBH4", 100, 100)
    b = doc.add_text("NaBH4", 100, 130, formula=True)
    assert b.box.width < a.box.width


# faces whose advance widths match the Arial/Helvetica table in textmetrics.ARIAL_WIDTHS; the sans fallback
# can also resolve to DejaVu Sans / Nimbus / FreeSans, which are ~10% wider and must not be held to the 2% check
ARIAL_METRIC_FACES = ("arial.ttf", "arialmt.ttf", "helvetica.ttc", "helvetica.ttf", "liberationsans-regular.ttf")


def test_heuristic_matches_measured_font():
    path = textmetrics.font_path("Arial")
    if path is None:
        pytest.skip("no font on this machine")
    face = os.path.basename(path).lower()
    if face not in ARIAL_METRIC_FACES:
        pytest.skip("measured font %s is not Arial-metric" % face)
    samples = ["H2SO4 (cat.)", "reflux, 4 h", "Pd(PPh3)4, K2CO3", "ethanol", "W", "i", "+", "acetic acid"]
    measured = {s: textmetrics.text_width(s, 10) for s in samples}
    old = os.environ.get("CHEMDRAWMCP_NO_FONTS")
    os.environ["CHEMDRAWMCP_NO_FONTS"] = "1"
    textmetrics.reset_cache()
    try:
        assert textmetrics.font_path("Arial") is None
        for s in samples:
            assert textmetrics.text_width(s, 10) == pytest.approx(measured[s], rel=0.02)
    finally:
        if old is None:
            del os.environ["CHEMDRAWMCP_NO_FONTS"]
        else:
            os.environ["CHEMDRAWMCP_NO_FONTS"] = old
        textmetrics.reset_cache()


def test_plus_glyph_centre():
    c = textmetrics.glyph_center("+", 12)
    assert 0.2 * 12 < c < 0.4 * 12


# ---------------------------------------------------------------------------
# overlap detection, determinism, speed
# ---------------------------------------------------------------------------

def test_collisions_detects_real_overlaps_but_not_arrow_text():
    doc = cdxml.CDXMLDocument("acs")
    a = doc.add_fragment(prep("c1ccccc1"), 100, 100)
    b = doc.add_fragment(prep("c1ccccc1"), 110, 100)
    pairs = layout.collisions(doc)
    assert len(pairs) == 1 and {pairs[0][0].id, pairs[0][1].id} == {a.id, b.id}
    assert any("overlaps" in w for w in layout.layout_warnings(doc))
    doc2 = cdxml.CDXMLDocument("acs")
    ar = doc2.add_arrow(100, 100, 180, 100)
    t = doc2.add_text("H2SO4", 140, 100 - 0.4 * doc2.bond_length, anchor="baseline", formula=True)
    assert overlap(t.box, ar.box)                 # inside the padded arrow box ...
    assert layout.collisions(doc2) == []          # ... but not touching the line
    bad = doc2.add_text("too close", 140, 100, anchor="center")
    assert any({p[0].id, p[1].id} == {ar.id, bad.id} for p in layout.collisions(doc2))


def test_layout_is_deterministic():
    def build():
        step = layout.Step([comp("[NH3+]CC([O-])=O"), comp("CC(=O)O")], [comp("CCO", label="85%")], above=["NaBH4"], below=["MeOH"])
        return layout.build_reaction(step, style="acs").to_xml()
    assert build() == build()


def test_scheme_layout_is_fast():
    smiles = ["c1ccccc1", "c1ccccc1Br", "c1ccccc1C(=O)O", "c1ccccc1C(=O)OC", "CC(=O)Oc1ccccc1C(=O)O",
              "[NH3+]CC([O-])=O", "OCC(O)CO", "c1ccc2ccccc2c1", "CC(C)(C)OC(=O)N[C@@H](C)C(=O)O", "O=C1CCCCC1"]
    mols = [prep(s) for s in smiles]
    steps = [layout.Step([layout.Component(mols[i], label=str(i + 1))], [layout.Component(mols[i + 1], label="%d%%" % (60 + i))],
                         above=["Pd(PPh3)4, K2CO3"], below=["DMF, 80 °C"]) for i in range(9)]
    layout.build_scheme(steps, style="acs")                      # warm up fonts / caches
    doc = layout.build_scheme(steps, style="acs", title="Scheme 1")
    rows = {round(f.box.cy) for f in frags(doc)}
    # chained: each intermediate drawn once, plus one repeat at the start of every wrapped row
    assert 10 <= len(frags(doc)) <= 10 + len(rows) - 1
    assert layout.collisions(doc) == []
    # best of several runs: a single wall-clock sample is at the mercy of the scheduler on a loaded CI runner
    dt = min(timeit.repeat(lambda: layout.build_scheme(steps, style="acs", title="Scheme 1"), number=1, repeat=5))
    assert dt < 0.2, "10-molecule scheme took %.0f ms" % (dt * 1e3)   # ~10 ms on a laptop; generous for CI


def test_scheme_wrapped_row_clears_the_tallest_step_of_the_row():
    S = "C1CCC2C(C1)CCC1C2CCC2CCCC12"      # gonane: tall structures in the first row
    steps = [layout.Step([comp(S), comp(S)], [comp(S, label="72% (2 steps)")], above=["X"]),
             layout.Step([comp(S)], [comp("CCCCCC")], above=["Y"]),
             layout.Step([comp("CCCCCC")], [comp("CC=O")], above=["reagent %d" % i for i in range(1, 7)])]
    doc = layout.build_scheme(steps, style="acs", max_width=450)
    fr = frags(doc)
    rows = sorted({round(f.box.cy) for f in fr})
    assert len(rows) == 2, rows
    # steps 2 and 3 share the wrapped row; the six-line reagent stack of step 3 (not the step that triggered the
    # wrap) is what has to clear the first row
    assert round(fr[-1].box.cy) == rows[1] and round(fr[-3].box.cy) == rows[1]
    assert layout.collisions(doc) == []
    assert not [w for w in doc.warnings if "overlaps" in w]
    first_of_row2 = next(f for f in fr if round(f.box.cy) == rows[1])
    k = doc.refs.index(first_of_row2)
    row1_bottom = max(r.box.y1 for r in doc.refs[:k])
    row2_top = min(r.box.y0 for r in doc.refs[k:])
    assert row2_top >= row1_bottom + 2.5 * doc.bond_length - 0.5


def test_scheme_steps_that_do_not_chain_start_a_new_row():
    steps = [layout.Step([comp("CCBr")], [comp("CCO")], above=["NaOH"]),
             layout.Step([comp("c1ccccc1Br")], [comp("c1ccccc1O")], above=["KOH"])]
    doc = layout.build_scheme(steps, style="acs")
    fr = frags(doc)
    assert len(fr) == 4 and len(doc.steps) == 2
    assert abs(fr[0].box.cy - fr[1].box.cy) < 1e-6 and abs(fr[2].box.cy - fr[3].box.cy) < 1e-6
    # the second reaction does not continue from the first: its own row, at the left margin, a row gap below
    assert fr[2].box.x0 == pytest.approx(fr[0].box.x0)
    row1_bottom = max(r.box.y1 for r in doc.refs[:doc.refs.index(fr[2])])
    assert fr[2].box.y0 >= row1_bottom + 2.5 * doc.bond_length - 0.5
    # chained steps still share a row; chain=False puts every step on its own row
    chained = [layout.Step([comp("c1ccccc1")], [comp("c1ccccc1Br")], above=["Br2"]),
               layout.Step([comp("c1ccccc1Br")], [comp("c1ccccc1C(=O)O")], above=["CO2"])]
    assert len({round(f.box.cy) for f in frags(layout.build_scheme(chained, style="acs"))}) == 1
    assert len({round(f.box.cy) for f in frags(layout.build_scheme(chained, style="acs", chain=False))}) == 2
