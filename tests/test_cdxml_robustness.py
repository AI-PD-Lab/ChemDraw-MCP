"""Robustness tests for the CDXML reader and the fragment writer chemistry.

Inline CDXML snippets cover nicknames (nested fragments and abbreviation tables), multi-attachment residues,
charges/isotopes derived from label text, explicit hydrogens, wedge/wavy/aromatic/dative bonds, document-level
bookkeeping (groups, pages, spectra, brackets, chemical properties) and a SMILES round-trip fuzz through the writer.
The ChemDraw sample drawings are parsed too when the discovered ``ChemDraw*.app`` bundle (see
:func:`chemdrawmcp.discovery.find_chemdraw`, ``CHEMDRAW_APP`` overrides the search) carries a Samples folder, or when
``CHEMDRAWMCP_SAMPLES`` points at a copy of that folder; otherwise those tests are skipped.
"""
import glob
import os
import re

import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem, rdBase  # noqa: E402

from chemdrawmcp import abbreviations, cdxml, chem, discovery  # noqa: E402


def _samples_dir() -> str:
    """ChemDraw's sample-drawing folder: ``CHEMDRAWMCP_SAMPLES``, else derived from the discovered ``.app`` bundle, else ""."""
    override = os.environ.get("CHEMDRAWMCP_SAMPLES")
    if override:
        return os.path.expanduser(override)
    app = discovery.find_chemdraw()
    if app and app.lower().endswith(".app"):
        return os.path.join(app, "Contents", "Resources", "SpecialPurpose", "Samples")
    return ""


SAMPLES = _samples_dir()


def doc(body: str, head: str = "") -> str:
    return '<?xml version="1.0" encoding="UTF-8" ?>\n<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >\n' \
           '<CDXML BondLength="14.4" %s><page id="1">%s</page></CDXML>' % (head, body)


def canon(smi: str) -> str:
    with rdBase.BlockLogs():
        return Chem.CanonSmiles(smi)


def first_mol(text: str) -> dict:
    back = cdxml.read_cdxml(text)
    assert back["molecules"], back
    return back["molecules"][0]


# ---------------------------------------------------------------------------
# abbreviation tables and label parsing
# ---------------------------------------------------------------------------

def test_abbreviation_tables_are_valid_smiles():
    with rdBase.BlockLogs():
        for label, smi in abbreviations.SUBSTITUENTS.items():
            m = Chem.MolFromSmiles(smi)
            assert m is not None, label
        for label, smi in abbreviations.MOLECULES.items():
            if smi:
                assert Chem.MolFromSmiles(smi) is not None, label
    for label in ("Me", "Et", "Pr", "iPr", "Bu", "tBu", "Ph", "Bn", "Bz", "Ac", "Boc", "Fmoc", "Cbz", "Ts", "Ms", "Tf", "TBS",
                  "TBDMS", "TIPS", "TMS", "PMB", "MOM", "THP", "OMe", "OEt", "OAc", "OTf", "NHBoc", "NMe2", "CF3", "CO2H",
                  "CO2Me", "CO2Et", "CN", "NO2", "SO3H", "PPh3", "N3"):
        assert abbreviations.lookup(label) is not None, label
    assert abbreviations.lookup("n-Bu") == ("CCCC", True)
    assert abbreviations.lookup("THF") == ("C1CCOC1", False)
    assert abbreviations.lookup("no such thing") is None
    assert abbreviations.lookup("hν") is None


def test_parse_atom_label():
    p = abbreviations.parse_atom_label
    assert p("NH3+") == {"symbol": "N", "hcount": 3, "charge": 1}
    assert p("O-") == {"symbol": "O", "charge": -1}
    assert p("O−") == {"symbol": "O", "charge": -1}
    assert p("Fe3+") == {"symbol": "Fe", "charge": 3} and p("Fe+3") == {"symbol": "Fe", "charge": 3}
    assert p("S2-") == {"symbol": "S", "charge": -2}
    assert p("13C") == {"symbol": "C", "isotope": 13}
    assert p("D") == {"symbol": "H", "isotope": 2} and p("T") == {"symbol": "H", "isotope": 3}
    assert p("OH") == {"symbol": "O", "hcount": 1}
    assert p("Ph") is None and p("OMe") is None and p("") is None


# ---------------------------------------------------------------------------
# nicknames and nested fragments
# ---------------------------------------------------------------------------

def test_nickname_without_inner_fragment_expands_from_table():
    text = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" NodeType="Nickname"><t p="14 3"><s font="3" size="10">OMe</s></t></n>'
               '<n id="5" p="-14 0" NodeType="Nickname"><t p="-14 3"><s font="3" size="10">Boc</s></t></n>'
               '<b id="6" B="3" E="4"/><b id="7" B="3" E="5"/></fragment>')
    m = first_mol(text)
    assert canon(m["smiles"]) == canon("COCC(=O)OC(C)(C)C")
    assert "labels" not in m and "error" not in m


def test_bare_substituent_label_keeps_attachment_point_and_reagent_is_molecule():
    text = doc('<fragment id="2"><n id="3" p="0 0" NodeType="Nickname"><t p="0 3"><s font="3" size="10">Ph</s></t></n></fragment>'
               '<fragment id="4"><n id="5" p="50 0" NodeType="Unspecified"><t p="50 3"><s font="3" size="10">NBS</s></t></n></fragment>'
               '<fragment id="6"><n id="7" p="90 0" NodeType="Fragment"><t p="90 3"><s font="3" size="10">1. Br2, H2O</s></t></n></fragment>')
    back = cdxml.read_cdxml(text)
    smiles = [canon(m["smiles"]) for m in back["molecules"]]
    assert smiles[0] == canon("*c1ccccc1")
    assert smiles[1] == canon("O=C1CCC(=O)N1Br")
    assert smiles[2] == "*" and back["molecules"][2]["labels"] == ["1. Br2, H2O"] and back["molecules"][2]["text"] == "1. Br2, H2O"


def test_nested_fragment_with_external_connection_point():
    text = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" NodeType="Nickname"><t p="14 0"><s font="3" size="10">OMe</s></t>'
               '<fragment id="5"><n id="6" p="0 0" Element="8"/><n id="7" p="14 0"/><n id="8" p="-14 0" NodeType="ExternalConnectionPoint"/>'
               '<b id="9" B="6" E="7"/><b id="10" B="8" E="6"/></fragment></n><b id="11" B="3" E="4"/></fragment>')
    assert first_mol(text)["smiles"] == "COC"


def test_multi_attachment_residues_use_connection_order_and_bond_ordering():
    # Gly-Ala dipeptide drawn as two residue nodes; each residue has an N-terminal and a C-terminal connection point,
    # and the inner fragments live in their own coordinate frame far from the sequence
    gly = ('<n id="10" p="0 0" NodeType="Fragment" BondOrdering="0 100"><fragment id="11" ConnectionOrder="12 13">'
           '<n id="12" p="500 500" NodeType="ExternalConnectionPoint"/><n id="14" p="520 500" Element="7" NumHydrogens="2"/>'
           '<n id="15" p="540 500"/><n id="16" p="560 500"/><n id="17" p="560 480" Element="8" NumHydrogens="0"/>'
           '<n id="13" p="580 500" NodeType="ExternalConnectionPoint"/>'
           '<b id="18" B="12" E="14"/><b id="19" B="14" E="15"/><b id="20" B="15" E="16"/><b id="21" B="16" E="17" Order="2"/><b id="22" B="16" E="13"/>'
           '</fragment></n>')
    ala = ('<n id="30" p="30 0" NodeType="Fragment" BondOrdering="100 0"><fragment id="31" ConnectionOrder="32 33">'
           '<n id="32" p="800 800" NodeType="ExternalConnectionPoint"/><n id="34" p="820 800" Element="7" NumHydrogens="1"/>'
           '<n id="35" p="840 800"/><n id="36" p="840 780"/><n id="37" p="860 800"/><n id="38" p="860 780" Element="8" NumHydrogens="0"/>'
           '<n id="39" p="880 800" Element="8" NumHydrogens="1"/><n id="33" p="820 820" NodeType="ExternalConnectionPoint"/>'
           '<b id="40" B="32" E="34"/><b id="41" B="34" E="35"/><b id="42" B="35" E="36" Display="WedgeBegin"/><b id="43" B="35" E="37"/>'
           '<b id="44" B="37" E="38" Order="2"/><b id="45" B="37" E="39"/><b id="46" B="34" E="33"/>'
           '</fragment></n>')
    # the chain bond is listed second in Gly's BondOrdering (C-terminus) and first in Ala's (N-terminus)
    text = doc('<fragment id="2">%s%s<b id="100" B="10" E="30"/></fragment>' % (gly, ala))
    m = first_mol(text)
    assert "error" not in m
    assert canon(m["smiles"]) in (canon("NCC(=O)N[C@@H](C)C(=O)O"), canon("NCC(=O)N[C@H](C)C(=O)O"))
    assert "@" in m["smiles"], "stereo inside the nested residue must survive"
    assert m["formula"] == "C5H10N2O3"


def test_nested_fragment_in_other_coordinate_frame_keeps_outer_stereo():
    """A stereocentre bonded to a nickname whose definition is drawn elsewhere at a different scale."""
    plain = ('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="-14 0" Element="17"/><n id="5" p="7 12"/>'
             '<n id="6" p="7 -12" Element="8"/><n id="7" p="21 -12"/>'
             '<b id="8" B="3" E="4" Display="WedgeBegin"/><b id="9" B="3" E="5"/><b id="10" B="3" E="6"/><b id="11" B="6" E="7"/></fragment>')
    nick = ('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="-14 0" Element="17"/><n id="5" p="7 12"/>'
            '<n id="6" p="7 -12" NodeType="Nickname"><t p="7 -12"><s font="3" size="10">OMe</s></t>'
            '<fragment id="20" ConnectionOrder="23"><n id="21" p="500 500" Element="8"/><n id="22" p="540 500"/>'
            '<n id="23" p="460 500" NodeType="ExternalConnectionPoint"/><b id="24" B="21" E="22"/><b id="25" B="23" E="21"/></fragment></n>'
            '<b id="8" B="3" E="4" Display="WedgeBegin"/><b id="9" B="3" E="5"/><b id="10" B="3" E="6"/></fragment>')
    a, b = first_mol(doc(plain))["smiles"], first_mol(doc(nick))["smiles"]
    assert "@" in a and canon(a) == canon(b)


def test_generic_nickname_element_list_and_ar_label_become_dummies():
    text = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" NodeType="GenericNickname" GenericNickname="R1"/>'
               '<n id="5" p="-14 0" NodeType="ElementList" ElementList="7 8"/>'
               '<n id="6" p="0 14" Element="18"><t p="0 17"><s font="3" size="10">Ar</s></t></n>'
               '<b id="7" B="3" E="4"/><b id="8" B="3" E="5"/><b id="9" B="3" E="6"/></fragment>')
    m = first_mol(text)
    assert m["smiles"].count("*") == 3 and "error" not in m
    assert sorted(m["labels"]) == ["Ar", "R1", "[N,O]"]


# ---------------------------------------------------------------------------
# atom attributes
# ---------------------------------------------------------------------------

def test_charge_hydrogens_and_isotope_from_label_text():
    text = doc('<fragment id="2"><n id="3" p="0 0" Element="7"><t p="0 3"><s font="3" size="10">NH3+</s></t></n><n id="4" p="14 0"/>'
               '<n id="5" p="28 0"/><n id="6" p="28 -14" Element="8"/><n id="7" p="42 0" Element="8"><t p="42 3"><s font="3" size="10">O-</s></t></n>'
               '<n id="8" p="14 14" Element="1"><t p="14 17"><s font="3" size="10">D</s></t></n>'
               '<b id="9" B="3" E="4"/><b id="10" B="4" E="5"/><b id="11" B="5" E="6" Order="2"/><b id="12" B="5" E="7"/><b id="13" B="4" E="8"/></fragment>')
    m = first_mol(text)
    assert canon(m["smiles"]) == canon("[NH3+]C([2H])C(=O)[O-]")


def test_charge_attribute_wins_over_label_and_radicals():
    text = doc('<fragment id="2"><n id="3" p="0 0" Element="8" Charge="-1" NumHydrogens="0"><t p="0 3"><s font="3" size="10">O</s></t></n>'
               '<n id="4" p="14 0" Radical="Doublet"/><b id="5" B="3" E="4"/></fragment>')
    m = first_mol(text)
    assert canon(m["smiles"]) == canon("[CH2][O-]")


def test_explicit_hydrogen_nodes_are_folded_but_isotopes_and_h2_kept():
    text = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" Element="8"/><n id="5" p="-14 0" Element="1" NumHydrogens="0"/>'
               '<n id="6" p="0 14" Element="1" Isotope="2" NumHydrogens="0"/>'
               '<b id="7" B="3" E="4" Order="2"/><b id="8" B="3" E="5"/><b id="9" B="3" E="6"/></fragment>'
               '<fragment id="10"><n id="11" p="60 0" Element="1" NumHydrogens="0"/><n id="12" p="74 0" Element="1" NumHydrogens="0"/><b id="13" B="11" E="12"/></fragment>'
               '<fragment id="14"><n id="15" p="100 0" Element="1" Charge="1" NumHydrogens="0"/></fragment>')
    back = cdxml.read_cdxml(text)
    smiles = [canon(m["smiles"]) for m in back["molecules"]]
    assert smiles == [canon("[2H]C=O"), "[H][H]", "[H+]"]
    assert back["molecules"][0]["formula"] == "CH2O"


# ---------------------------------------------------------------------------
# bonds and stereo
# ---------------------------------------------------------------------------

def _written(smi: str, **settings) -> str:
    d = cdxml.CDXMLDocument("acs", settings=settings or None)
    d.add_fragment(chem.mol_from_input(smi), 100, 100)
    return d.to_xml()


def test_wedge_end_forms_equal_wedge_begin_forms():
    xml = _written("N[C@@H](C)C(=O)O")
    ref = canon(first_mol(xml)["smiles"])
    assert "@" in ref

    def swap(m):
        b, e = m.group(2), m.group(3)
        return '%sB="%s" E="%s"%s' % (m.group(1), e, b, m.group(4))
    flipped = re.sub(r'(<b [^>]*?)B="(\d+)" E="(\d+)"([^>]*Display="Wedge(?:dHash)?Begin"[^>]*/>)', swap, xml)
    flipped = flipped.replace('Display="WedgeBegin"', 'Display="WedgeEnd"').replace('Display="WedgedHashBegin"', 'Display="WedgedHashEnd"')
    assert flipped != xml
    assert canon(first_mol(flipped)["smiles"]) == ref


def test_ez_from_coordinates_and_wavy_bonds():
    trans = doc('<fragment id="2"><n id="3" p="-7 -12"/><n id="4" p="0 0"/><n id="5" p="14 0"/><n id="6" p="21 12"/>'
                '<b id="7" B="3" E="4"/><b id="8" B="4" E="5" Order="2"/><b id="9" B="5" E="6"/></fragment>')
    cis = trans.replace('p="21 12"', 'p="21 -12"')
    assert first_mol(trans)["smiles"] == "C/C=C/C"
    assert first_mol(cis)["smiles"] == "C/C=C\\C"
    wavy_single = trans.replace('<b id="9" B="5" E="6"/>', '<b id="9" B="5" E="6" Display="Wavy"/>')
    assert first_mol(wavy_single)["smiles"] == "CC=CC"
    wavy_double = trans.replace('Order="2"', 'Order="2" Display="Wavy"')
    assert first_mol(wavy_double)["smiles"] == "CC=CC"
    # wavy single bond at a tetrahedral centre: unknown configuration
    chiral = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="-14 0" Element="17"/><n id="5" p="7 12" Element="35"/><n id="6" p="7 -12" Element="9"/>'
                 '<b id="7" B="3" E="4" Display="WedgeBegin"/><b id="8" B="3" E="5"/><b id="9" B="3" E="6"/></fragment>')
    assert "@" in first_mol(chiral)["smiles"]
    assert "@" not in first_mol(chiral.replace('Display="WedgeBegin"', 'Display="Wavy"'))["smiles"]


def test_aromatic_order_rings_kekulize_and_chain_aromatic_bonds_downgrade():
    ring = doc('<fragment id="2">' + "".join('<n id="%d" p="%d %d"/>' % (10 + i, [0, 12, 12, 0, -12, -12][i], [-14, -7, 7, 14, 7, -7][i]) for i in range(6)) +
               "".join('<b id="%d" B="%d" E="%d" Order="1.5"/>' % (20 + i, 10 + i, 10 + (i + 1) % 6) for i in range(6)) + '</fragment>')
    m = first_mol(ring)
    assert m["smiles"] == "c1ccccc1" and "error" not in m
    chain = doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" Element="8"/><b id="5" B="3" E="4" Order="1.5"/></fragment>')
    m = first_mol(chain)
    assert m["smiles"] == "CO" and "error" not in m and any("outside a ring" in w for w in m["warnings"])
    # cyclopentadienyl without a charge cannot be kekulized: reported, aromatic SMILES still returned
    cp = doc('<fragment id="2">' + "".join('<n id="%d" p="%d %d"/>' % (10 + i, [0, 13, 8, -8, -13][i], [-14, -4, 11, 11, -4][i]) for i in range(5)) +
             "".join('<b id="%d" B="%d" E="%d" Order="1.5"/>' % (20 + i, 10 + i, 10 + (i + 1) % 5) for i in range(5)) + '</fragment>')
    m = first_mol(cp)
    assert "error" in m and "kekulize" in m["error"].lower() and m["smiles"] == "c1cccc1"


def test_dative_duplicate_and_decorative_bond_displays():
    text = doc('<fragment id="2"><n id="3" p="0 0" Element="7" NumHydrogens="3"/><n id="4" p="14 0" Element="29" Charge="2"/>'
               '<b id="5" B="3" E="4" Order="dative"/><b id="6" B="4" E="3" Order="dative"/></fragment>'
               '<fragment id="7"><n id="8" p="50 0"/><n id="9" p="64 0"/><n id="10" p="78 0" Element="8"/>'
               '<b id="11" B="8" E="9" Display="Bold"/><b id="12" B="9" E="10" Display="Hash"/><b id="13" B="8" E="9" Display="Dash"/></fragment>')
    back = cdxml.read_cdxml(text)
    m1, m2 = back["molecules"]
    assert m1["smiles"] == "[NH3]->[Cu+2]" and any("duplicate" in w for w in m1["warnings"])
    assert m2["smiles"] == "CCO" and any("duplicate" in w for w in m2["warnings"])
    mols = cdxml.cdxml_mols(text)
    assert mols[0].GetBondWithIdx(0).GetBondType() == Chem.BondType.DATIVE


def test_invalid_valence_gives_error_entry_with_best_effort_smiles():
    text = doc('<fragment id="2"><n id="3" p="0 0"/>' + "".join('<n id="%d" p="%d %d"/>' % (10 + i, 14 * (i % 3) - 14, 14 * (i // 3) + 14) for i in range(5)) +
               "".join('<b id="%d" B="3" E="%d"/>' % (20 + i, 10 + i) for i in range(5)) + '</fragment>')
    back = cdxml.read_cdxml(text)
    m = back["molecules"][0]
    assert "error" in m and "valence" in m["error"].lower() and m["smiles"]
    assert any("fragment 2" in w for w in back["warnings"])


# ---------------------------------------------------------------------------
# document-level bookkeeping
# ---------------------------------------------------------------------------

def test_document_counts_and_chemical_property_association():
    text = doc('<group id="5"><fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" Element="8"/><b id="6" B="3" E="4"/></fragment>'
               '<t id="7" p="0 30"><s font="3" size="10">methanol</s></t></group>'
               '<chemicalproperty id="8" ChemicalPropertyDisplayID="7" ChemicalPropertyType="1" BasisObjects="3 4"/>'
               '<t id="9" p="0 60"><s font="3" size="10">Residue 1</s></t>'
               '<chemicalproperty id="10" Name="Residue" ChemicalPropertyType="0" BasisObjects="2"/>'
               '<spectrum id="11" XSpacing="1" XLow="0" XType="MS" YType="Abundance"/>'
               '<bracketedgroup id="12" BracketedObjectIDs="3 4" RepeatCount="2"><bracketattachment id="13" GraphicID="14"/></bracketedgroup>'
               '<table id="15"><page id="16"><fragment id="17"><n id="18" p="100 0" Element="17" Charge="-1" NumHydrogens="0"/><n id="19" p="114 0" Element="11" Charge="1" NumHydrogens="0"/></fragment></page></table>'
               '<fragment id="20"><graphic id="21" GraphicType="Oval"/></fragment>')
    back = cdxml.read_cdxml(text)
    assert back["pages"] == 1 and back["groups"] == 1 and back["spectra"] == 1
    assert back["bracketed_groups"] == [{"id": "12", "objects": 2, "repeatCount": "2"}]
    assert any("bracketed" in w for w in back["warnings"])
    assert back["other_objects"]["fragment (no atoms)"] == 1 and back["other_objects"]["spectrum"] == 1
    mols = {m["id"]: m for m in back["molecules"]}
    assert set(mols) == {"2", "17"}                       # the table cell's fragment is found, the empty one is not a molecule
    assert mols["2"]["name"] == "methanol"
    assert [p["text"] for p in mols["2"]["properties"]] == ["methanol", "Residue 1"]
    assert mols["2"]["properties"][1]["name"] == "Residue"
    texts = {t["id"]: t for t in back["texts"]}
    assert texts["7"]["chemical_property"]["ChemicalPropertyType"] == "1"
    assert texts["9"]["chemical_property"]["Name"] == "Residue"
    assert canon(mols["17"]["smiles"]) == canon("[Na+].[Cl-]")
    assert sorted(cdxml.document_text(text)) == ["Residue 1", "methanol"]


def test_steps_expand_groups_and_skip_text_only_agents():
    text = doc('<group id="30"><fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0" Element="8"/><b id="5" B="3" E="4"/></fragment>'
               '<fragment id="6"><n id="7" p="0 40" NodeType="Unspecified"><t p="0 43"><s font="3" size="10">cat. H+</s></t></n></fragment></group>'
               '<fragment id="8"><n id="9" p="100 0"/><n id="10" p="114 0"/><b id="11" B="9" E="10"/></fragment>'
               '<arrow id="12" Head3D="80 0 0" Tail3D="40 0 0"/><t id="13" p="60 -10"><s font="3" size="10">heat</s></t>'
               '<scheme id="14"><step id="15" ReactionStepReactants="30" ReactionStepProducts="8" ReactionStepArrows="12" ReactionStepObjectsAboveArrow="13 6"/></scheme>')
    back = cdxml.read_cdxml(text)
    st = back["steps"][0]
    assert [r.get("smiles") or r.get("text") for r in st["reactants"]] == ["CO", "cat. H+"]
    assert [a.get("text") for a in st["above_arrow"]] == ["heat", "cat. H+"]
    assert back["reaction_smiles"] == ["CO>>CC"]


def test_parse_accepts_bom_and_bytes():
    text = "﻿" + doc('<fragment id="2"><n id="3" p="0 0"/><n id="4" p="14 0"/><b id="5" B="3" E="4"/></fragment>')
    assert cdxml.read_cdxml(text)["molecules"][0]["smiles"] == "CC"
    assert cdxml.read_cdxml(text.encode("utf-8"))["molecules"][0]["smiles"] == "CC"


# ---------------------------------------------------------------------------
# writer round trip
# ---------------------------------------------------------------------------

ROUND_TRIP = [
    "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "COc1ccc2[nH]cc(CCN(C)C)c2c1",
    "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "N[C@@H](Cc1ccccc1)C(=O)O", "C[C@@H](O)[C@@H](N)C(=O)O",
    "CC[C@H](C)[C@@H](N)C(=O)O", "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O",
    "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C",
    "C/C=C/C(=O)O", "C/C=C\\C(=O)O", "F/C=C/Cl", "F/C=C\\Cl", "C/C=C/C=C/C", "C/C=C\\C=C/C", "CC(C)=CCC/C(C)=C/CO",
    "C=C=C", "CC=C=CC", "C[C@H](F)C=C=C", "C1CCCCCCC/C=C/1", "C1CCCCCCCCCCCC(=O)OCC1", "C[C@]12CC[C@H](C1)C2",
    "O[C@H]1CC[C@@H](O)CC1", "O[C@H]1CC[C@H](O)CC1", "C[C@@H]1CCCC[C@H]1C", "F[C@H](Cl)Br", "F[C@@H](Cl)Br",
    "[Na+].[Cl-]", "[NH4+].[Cl-]", "[K+].[K+].[O-]C(=O)[O-]", "C[N+](C)(C)C.[Cl-]", "[NH3+][C@@H](C)C(=O)[O-]",
    "[O-][N+](=O)c1ccccc1", "CN=[N+]=[N-]", "[C-]#[O+]", "[13CH3]C(=O)O", "[2H]C([2H])([2H])C(=O)O", "C[13C@H](N)C(=O)O",
    "C[CH2]", "[CH3]", "CC1(C)CCCC(C)(C)N1[O]", "[O][O]", "C[Mg+]", "C[Mg]Br", "[Li]CCCC", "C[Si](C)(C)Cl", "O=[Os](=O)(=O)=O",
    "[Fe+2].[cH-]1cccc1.[cH-]1cccc1", "F[B-](F)(F)F.[K+]", "[Cu+2]<-N",
    "CC(C)(C)C(C(C)(C)C)(C(C)(C)C)C(C)(C)C", "c1ccc2c(c1)ccc1ccccc12", "O=c1cc[nH]c(=O)[nH]1", "c1ccc(cc1)P(c1ccccc1)c1ccccc1",
    "[CH3:1][OH:2]", "[NH2:1][CH2:2][C:3](=[O:4])[OH:5]", "C[N+](C)(C)[O-]", "O=S1(=O)CCCC1", "C1CC1C1CC1",
    "CN1CCN(CC1)c1ccc(cc1)Nc1ncc(c(n1)-c1cnn(C)c1)Cl", "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O",
    "N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N.[K+].[K+].[K+].[K+]",
]


@pytest.mark.parametrize("smi", ROUND_TRIP)
def test_round_trip_through_writer_and_reader(smi):
    m = first_mol(_written(smi))
    assert "error" not in m, m
    assert canon(m["smiles"]) == canon(smi)


def test_unspecified_double_bond_wavy_option():
    # by default an unspecified alkene is drawn flat and reads back with the drawn geometry (E)
    assert first_mol(_written("CC=CC"))["smiles"] == "C/C=C/C"
    xml = _written("CC=CC", wavy_unspecified_double_bonds=True)
    assert 'Display="Wavy"' in xml and first_mol(xml)["smiles"] == "CC=CC"
    # specified geometry is never made wavy, nor are ring double bonds in small rings
    xml = _written("C/C=C/C.C1=CCCC1", wavy_unspecified_double_bonds=True)
    assert "Wavy" not in xml


def test_writer_dummy_atoms_have_no_element_attribute():
    xml = _written("*c1ccccc1")
    assert 'Element="0"' not in xml and 'GenericNickname="R"' in xml
    assert first_mol(xml)["smiles"] == "*c1ccccc1"
    m = chem.mol_from_input("*c1ccccc1")
    m.GetAtomWithIdx(0).SetProp("dummyLabel", "R1")
    d = cdxml.CDXMLDocument("acs")
    d.add_fragment(m, 100, 100)
    xml = d.to_xml()
    assert 'GenericNickname="R1"' in xml and first_mol(xml)["labels"] == ["R1"]


# ---------------------------------------------------------------------------
# ChemDraw's own sample drawings (skipped when the application folder is missing)
# ---------------------------------------------------------------------------

def test_samples_dir_follows_discovery_and_override(tmp_path, monkeypatch):
    """The sample folder is not pinned to one ChemDraw version: it follows find_chemdraw() or CHEMDRAWMCP_SAMPLES."""
    monkeypatch.delenv("CHEMDRAWMCP_SAMPLES", raising=False)
    monkeypatch.delenv("CHEMDRAW_EXE", raising=False)
    app = tmp_path / "ChemDraw 24.0.1.app"
    (app / "Contents").mkdir(parents=True)
    monkeypatch.setenv("CHEMDRAW_APP", str(app))
    assert _samples_dir() == str(app / "Contents" / "Resources" / "SpecialPurpose" / "Samples")
    exe = tmp_path / "ChemDraw.exe"
    exe.write_bytes(b"")
    monkeypatch.delenv("CHEMDRAW_APP")
    monkeypatch.setenv("CHEMDRAW_EXE", str(exe))
    assert _samples_dir() == ""                                        # Windows: no bundle, no samples
    monkeypatch.setattr(discovery, "find_chemdraw", lambda: None)
    assert _samples_dir() == ""                                        # Linux / nothing installed
    monkeypatch.setenv("CHEMDRAWMCP_SAMPLES", str(tmp_path / "samples"))
    assert _samples_dir() == str(tmp_path / "samples")                 # explicit override wins


sample_files = sorted(glob.glob(os.path.join(SAMPLES, "*.cdxml"))) if SAMPLES else []
needs_samples = pytest.mark.skipif(not sample_files, reason="ChemDraw sample folder not available")


def _sample(name: str) -> str:
    with open(os.path.join(SAMPLES, name), encoding="utf-8", errors="replace") as fh:
        return fh.read()


@needs_samples
@pytest.mark.parametrize("path", sample_files, ids=[os.path.basename(p) for p in sample_files])
def test_every_sample_file_reads_cleanly(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    back = cdxml.read_cdxml(text)
    mols = cdxml.cdxml_mols(text)
    assert back["pages"] >= 1 and isinstance(back["groups"], int)
    assert len(mols) == len(back["molecules"])
    for m in back["molecules"]:
        assert m.get("smiles") is not None, m
        if "error" in m:
            # only drawings ChemDraw itself flags as chemically invalid (orbital sketches, hypervalent atoms)
            assert os.path.basename(path).startswith(("Carbon Hybridization", "Hotkeys Cheat Sheet", "Stereoisomerism")), m
    for mol in mols:
        assert mol.GetNumConformers() == 1


@needs_samples
def test_sample_peptides_expand_residues():
    cyc = cdxml.read_cdxml(_sample("Cyclic Peptide.cdxml"))["molecules"]
    assert cyc[0]["formula"] == "C37H38N6O6" and cyc[0]["smiles"].count("@") >= 4
    assert any(p.get("name") == "Residue" for m in cyc for p in m.get("properties", []))
    vas = cdxml.read_cdxml(_sample("Vasopressin Peptides.cdxml"))
    assert vas["pages"] == 1 and vas["other_objects"].get("table") == 1
    assert vas["molecules"][0]["formula"] == "C46H65N15O12S2"          # arginine vasopressin
    assert "SS" in vas["molecules"][0]["smiles"]                       # the disulfide bridge via the third connection point
    hormones = cdxml.read_cdxml(_sample("Pituitary Hormones.cdxml"))["molecules"]
    assert hormones[0]["formula"] == "C43H67N11O13S2"                  # oxytocin, free acid as drawn


@needs_samples
def test_sample_names_stereo_and_steps():
    citric = cdxml.read_cdxml(_sample("Citric Acid Cycle.cdxml"))
    assert "Oxaloacetate" in [m.get("name") for m in citric["molecules"]]
    assert citric["steps"] and citric["reaction_smiles"]
    stereo = cdxml.read_cdxml(_sample("Stereoisomerism.cdxml"))
    assert sum("/" in m["smiles"] for m in stereo["molecules"]) >= 4
    # 2D wedge drawings (halomethanes, sugars) read; perspective sketches (chairs, wedged ring bonds) cannot
    assert sum("@" in m["smiles"] for m in stereo["molecules"]) >= 6
    alkenes = cdxml.read_cdxml(_sample("Alkenes.cdxml"))
    agents = [x.get("text") for s in alkenes["steps"] for x in s["above_arrow"] + s["below_arrow"]]
    assert "1. mCPBA" in agents
    assert cdxml.read_cdxml(_sample("b-Cypermethrin (stereochemistry).cdxml"))["bracketed_groups"][0]["repeatCount"] == "2"
    assert cdxml.read_cdxml(_sample("Hotkeys Cheat Sheet.cdxml"))["groups"] > 50
