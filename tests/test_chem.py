"""Offline tests for chemdrawmcp.chem (every public function)."""
from __future__ import annotations

import io
import sys
import urllib.request

import numpy as np
import pytest

from chemdrawmcp import chem
from chemdrawmcp.chem import ChemError, Reaction

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem, rdBase  # noqa: E402
from rdkit.Chem import AllChem, rdChemReactions  # noqa: E402
from rdkit.Geometry import Point3D  # noqa: E402

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
S_ALANINE = "C[C@H](N)C(=O)O"
R_IBUPROFEN = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
ESTERIFICATION = "CC(=O)O.OCC>[H+]>CC(=O)OCC.O"

CDXML_TWO = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML CreationProgram="test" BondLength="30">
<page id="2">
<fragment id="10">
<n id="11" p="100 100" Element="6"/>
<n id="12" p="126 115" Element="8"/>
<b id="13" B="11" E="12" Order="1"/>
</fragment>
<fragment id="20">
<n id="21" p="200 100" Element="6"/>
<n id="22" p="226 115" Element="7"/>
<n id="23" p="252 100" Element="6"/>
<b id="24" B="21" E="22" Order="1"/>
<b id="25" B="22" E="23" Order="1"/>
</fragment>
</page>
</CDXML>
"""

CDXML_SCHEME = CDXML_TWO.replace(
    "</page>",
    '<scheme id="30"><step id="31" ReactionStepReactants="10" '
    'ReactionStepProducts="20"/></scheme>\n</page>',
)


def canon(smiles: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def canon_mol(mol) -> str:
    m = Chem.Mol(mol)
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


def wedged_from(mol, atom_idx: int):
    """Bonds starting at atom_idx carrying a wedge/dash direction."""
    dirs = (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    return [
        b for b in mol.GetBonds()
        if b.GetBondDir() in dirs and b.GetBeginAtomIdx() == atom_idx
    ]


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    """Keep every test offline and independent of OPSIN installs."""
    monkeypatch.delenv("CHEMDRAWMCP_ONLINE", raising=False)
    monkeypatch.setitem(sys.modules, "py2opsin", None)
    monkeypatch.delenv("OPSIN_JAR", raising=False)
    monkeypatch.setattr(chem.shutil, "which", lambda *a, **k: None)

    def no_network(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("network access attempted while offline")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    yield


@pytest.fixture
def aspirin():
    return chem.mol_from_input(ASPIRIN)


@pytest.fixture
def sdf_two(tmp_path):
    a = Chem.MolFromSmiles("CCO")
    a.SetProp("_Name", "ethanol")
    b = Chem.MolFromSmiles("CC(=O)O")
    b.SetProp("_Name", "acetic acid")
    text = chem.mol_to_sdf([a, b])
    path = tmp_path / "two.sdf"
    path.write_text(text)
    return text, path


# --------------------------------------------------------------------------- #
# rdkit_available / detect_format
# --------------------------------------------------------------------------- #
def test_rdkit_available():
    assert chem.rdkit_available() is True


@pytest.mark.parametrize(
    "data, expected",
    [
        ("CCO", "smiles"),
        ("CCO ethanol", "smiles"),
        ("c1ccccc1 benzene\nCCO ethanol", "smiles"),
        ("[#6]", "smiles"),
        ("[CH3:1]C(=O)[OH:2]", "smiles"),
        ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "inchi"),
        (ESTERIFICATION, "reaction_smiles"),
        ("CCO -> CC=O", "reaction_smiles"),
        ("CCO >> CC=O", "reaction_smiles"),
        ("acetic acid", "name"),
        ("aspirin", "name"),
        ("2-(acetyloxy)benzoic acid", "name"),
        ("$RXN\n\n  RDKit\n\n  2  1\n", "rxn"),
        (CDXML_TWO, "cdxml"),
    ],
)
def test_detect_format_strings(data, expected):
    assert chem.detect_format(data) == expected


def test_detect_format_molblock_and_sdf(sdf_two):
    molblock = Chem.MolToMolBlock(Chem.MolFromSmiles(ASPIRIN))
    assert chem.detect_format(molblock) == "molblock"
    text, _ = sdf_two
    assert chem.detect_format(text) == "sdf"
    v3000 = Chem.MolToMolBlock(Chem.MolFromSmiles(ASPIRIN), forceV3000=True)
    assert chem.detect_format(v3000) == "molblock"


def test_detect_format_paths(tmp_path, sdf_two):
    mol_path = tmp_path / "asp.mol"
    mol_path.write_text(Chem.MolToMolBlock(Chem.MolFromSmiles(ASPIRIN)))
    assert chem.detect_format(str(mol_path)) == "molblock"
    _, sdf_path = sdf_two
    assert chem.detect_format(str(sdf_path)) == "sdf"
    cdxml_path = tmp_path / "frag.cdxml"
    cdxml_path.write_text(CDXML_TWO)
    assert chem.detect_format(str(cdxml_path)) == "cdxml"
    rxn_path = tmp_path / "r.rxn"
    rxn_path.write_text("$RXN\n")
    assert chem.detect_format(str(rxn_path)) == "rxn"
    smi_path = tmp_path / "list.smi"
    smi_path.write_text("CCO ethanol\n")
    assert chem.detect_format(str(smi_path)) == "smiles"


def test_detect_format_empty():
    with pytest.raises(ChemError):
        chem.detect_format("   ")


# --------------------------------------------------------------------------- #
# mol_from_input
# --------------------------------------------------------------------------- #
def test_mol_from_smiles_prepared(aspirin):
    assert aspirin.GetNumAtoms() == 13
    assert aspirin.GetNumConformers() == 1
    assert not any(b.GetIsAromatic() for b in aspirin.GetBonds())
    assert not any(a.GetIsAromatic() for a in aspirin.GetAtoms())
    assert canon_mol(aspirin) == canon(ASPIRIN)
    # implicit hydrogens are available without further sanitization
    assert aspirin.GetAtomWithIdx(0).GetTotalNumHs() == 3


def test_mol_from_smiles_unprepared():
    m = chem.mol_from_input(ASPIRIN, prepare=False)
    assert any(b.GetIsAromatic() for b in m.GetBonds())
    assert m.GetNumConformers() == 0


def test_mol_from_smiles_with_name():
    m = chem.mol_from_input("CCO ethanol")
    assert m.GetProp("_Name") == "ethanol"
    assert canon_mol(m) == "CCO"


def test_mol_from_inchi():
    m = chem.mol_from_input("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    assert canon_mol(m) == "CCO"


def test_mol_from_molblock_keeps_2d_coords():
    ref = Chem.MolFromSmiles(ASPIRIN)
    AllChem.Compute2DCoords(ref)
    conf = ref.GetConformer()
    for i in range(ref.GetNumAtoms()):  # translate so a recomputation is distinguishable
        x, y, _ = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(x + 5.0, y + 3.0, 0.0))
    block = Chem.MolToMolBlock(ref)
    m = chem.mol_from_input(block)
    ref_pos = ref.GetConformer().GetPositions()
    assert np.allclose(m.GetConformer().GetPositions(), ref_pos, atol=1e-3)  # molblock has 4 decimals
    m2 = chem.prepare_for_depiction(Chem.MolFromMolBlock(block), keep_coords=False)
    assert m2.GetNumConformers() == 1
    assert not np.allclose(m2.GetConformer().GetPositions(), ref_pos, atol=1e-3)


def test_mol_from_3d_molblock_recomputes_2d_and_keeps_chirality():
    m3 = Chem.AddHs(Chem.MolFromSmiles(S_ALANINE))
    assert AllChem.EmbedMolecule(m3, randomSeed=11) == 0
    block = Chem.MolToMolBlock(m3)
    m = chem.mol_from_input(block)
    pos = m.GetConformer().GetPositions()
    assert np.all(np.abs(pos[:, 2]) < 1e-6)
    assert m.GetNumAtoms() == 6  # hydrogens removed
    assert canon_mol(m) == canon(S_ALANINE)
    assert chem.properties(m)["stereocenters"]["labels"] == {1: "S"}


def test_mol_from_cdxml_text_and_file(tmp_path):
    m = chem.mol_from_input(CDXML_TWO)
    assert canon_mol(m) == "CO"
    path = tmp_path / "frag.cdxml"
    path.write_text(CDXML_TWO)
    assert canon_mol(chem.mol_from_input(str(path))) == "CO"


def test_mol_from_mol_file_path(tmp_path):
    path = tmp_path / "asp.mol"
    path.write_text(Chem.MolToMolBlock(Chem.MolFromSmiles(ASPIRIN)))
    m = chem.mol_from_input(str(path))
    assert canon_mol(m) == canon(ASPIRIN)


def test_mol_from_name_offline():
    m = chem.mol_from_input("aspirin", fmt="name")
    assert canon_mol(m) == canon(ASPIRIN)
    assert m.GetProp("_Name") == "aspirin"
    m2 = chem.mol_from_input("Acetic acid")  # auto-detected as a name
    assert canon_mol(m2) == "CC(=O)O"


def test_mol_from_unknown_name_raises():
    with pytest.raises(ChemError, match="unobtainium"):
        chem.mol_from_input("unobtainium")


def test_mol_from_bad_smiles_raises():
    with pytest.raises(ChemError):
        chem.mol_from_input("C1CC(", fmt="smiles")


def test_mol_from_reaction_raises():
    with pytest.raises(ChemError, match="reaction"):
        chem.mol_from_input(ESTERIFICATION)


def test_mol_from_input_special_species():
    mapped = chem.mol_from_input("[CH3:1]C(=O)[OH:2]")
    assert [a.GetAtomMapNum() for a in mapped.GetAtoms()] == [1, 0, 0, 2]
    charged = chem.mol_from_input("[NH4+]")
    assert Chem.GetFormalCharge(charged) == 1
    iso = chem.mol_from_input("[13CH4]")
    assert iso.GetAtomWithIdx(0).GetIsotope() == 13
    salt = chem.mol_from_input("[Na+].[Cl-]")
    assert salt.GetNumAtoms() == 2 and len(Chem.GetMolFrags(salt)) == 2
    xs = salt.GetConformer().GetPositions()[:, 0]
    assert abs(xs[0] - xs[1]) > 1.0  # fragments are laid out apart
    smarts = chem.mol_from_input("[CX4][OX2H]", fmt="smarts")  # query SMARTS -> plain structure
    assert canon_mol(smarts) == "CO"
    assert canon_mol(chem.mol_from_input("[#6][#8]")) == "[C][O]"  # RDKit SMILES extension, as-is


def test_unknown_format_raises():
    with pytest.raises(ChemError, match="Unknown format"):
        chem.mol_from_input("CCO", fmt="xyz")


# --------------------------------------------------------------------------- #
# mols_from_input
# --------------------------------------------------------------------------- #
def test_mols_from_sdf_two_records(sdf_two):
    text, path = sdf_two
    mols = chem.mols_from_input(text)
    assert [canon_mol(m) for m in mols] == ["CCO", "CC(=O)O"]
    assert [m.GetProp("_Name") for m in mols] == ["ethanol", "acetic acid"]
    assert all(m.GetNumConformers() == 1 for m in mols)
    mols2 = chem.mols_from_input(str(path))
    assert len(mols2) == 2


def test_mols_from_smiles_list_with_names():
    text = "CCO ethanol\nCC(=O)O acetic acid\nc1ccccc1 benzene\n"
    mols = chem.mols_from_input(text)
    assert len(mols) == 3
    assert [m.GetProp("_Name") for m in mols] == ["ethanol", "acetic acid", "benzene"]
    assert canon_mol(mols[2]) == "c1ccccc1"


def test_mols_from_whitespace_separated_smiles():
    mols = chem.mols_from_input("CCO CCN CCC")
    assert [canon_mol(m) for m in mols] == ["CCO", "CCN", "CCC"]


def test_mols_from_smiles_list_reports_bad_line():
    with pytest.raises(ChemError, match="line 2"):
        chem.mols_from_input("CCO\nC1CC(\nCCC")


def test_mols_from_cdxml_all_fragments():
    mols = chem.mols_from_input(CDXML_TWO)
    assert [canon_mol(m) for m in mols] == ["CO", "CNC"]


def test_mols_from_reaction_components():
    mols = chem.mols_from_input(ESTERIFICATION)
    assert [canon_mol(m) for m in mols] == ["CC(=O)O", "CCO", "[H+]", "CCOC(C)=O", "O"]


def test_mols_from_names():
    mols = chem.mols_from_input("water\nethanol", fmt="name")
    assert [canon_mol(m) for m in mols] == ["O", "CCO"]


# --------------------------------------------------------------------------- #
# prepare_for_depiction
# --------------------------------------------------------------------------- #
def test_prepare_wedge_s_alanine():
    m = chem.mol_from_input(S_ALANINE)
    bonds = wedged_from(m, 1)
    assert len(bonds) == 1
    assert bonds[0].GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)
    assert canon_mol(m) == canon(S_ALANINE)
    assert chem.properties(m)["stereocenters"]["labels"] == {1: "S"}


def test_prepare_wedge_r_ibuprofen():
    m = chem.mol_from_input(R_IBUPROFEN)
    centre = [i for i, lab in Chem.FindMolChiralCenters(m, useLegacyImplementation=False)]
    assert len(centre) == 1
    assert len(wedged_from(m, centre[0])) == 1
    assert chem.properties(m)["stereocenters"] == {
        "assigned": 1, "unassigned": 0, "atoms": centre, "labels": {centre[0]: "R"}
    }
    # a wedge only ever starts at the stereocentre
    for b in m.GetBonds():
        if b.GetBondDir() != Chem.BondDir.NONE:
            assert b.GetBeginAtomIdx() == centre[0]


def test_prepare_no_wedge_without_stereo(aspirin):
    assert all(b.GetBondDir() == Chem.BondDir.NONE for b in aspirin.GetBonds())


def test_prepare_options():
    m = chem.prepare_for_depiction(Chem.MolFromSmiles(ASPIRIN), kekulize=False)
    assert any(b.GetIsAromatic() for b in m.GetBonds())
    m = chem.prepare_for_depiction(Chem.MolFromSmiles(S_ALANINE), wedge=False)
    assert all(b.GetBondDir() == Chem.BondDir.NONE for b in m.GetBonds())


def test_prepare_keeps_double_bond_stereo():
    m = chem.mol_from_input("C/C=C/C(=O)O")
    assert canon_mol(m) == canon("C/C=C/C(=O)O")
    assert chem.properties(m)["double_bond_stereo"] == {"assigned": 1, "unassigned": 0}


def test_prepare_hydrogen_handling():
    assert chem.mol_from_input("[H]C([H])([H])O").GetNumAtoms() == 2
    assert chem.mol_from_input("[2H]C").GetNumAtoms() == 2  # isotopic H kept
    assert chem.mol_from_input("[H+]").GetNumAtoms() == 1


def test_prepare_does_not_modify_input():
    src = Chem.MolFromSmiles(ASPIRIN)
    chem.prepare_for_depiction(src)
    assert any(b.GetIsAromatic() for b in src.GetBonds())
    assert src.GetNumConformers() == 0


def test_prepare_sanitize_error():
    bad = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    with pytest.raises(ChemError, match="valence"):
        chem.prepare_for_depiction(bad)


# --------------------------------------------------------------------------- #
# properties / elemental_analysis / describe
# --------------------------------------------------------------------------- #
def test_properties_aspirin(aspirin):
    p = chem.properties(aspirin)
    assert p["formula"] == "C9H8O4"
    assert abs(p["molecular_weight"] - 180.159) < 0.01
    assert abs(p["exact_mass"] - 180.04226) < 1e-4
    assert p["elemental_analysis"] == "C, 60.00; H, 4.48; O, 35.52"
    assert p["canonical_smiles"] == canon(ASPIRIN)
    assert p["isomeric_smiles"] == canon(ASPIRIN)
    assert p["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert p["inchi"].startswith("InChI=1S/C9H8O4/")
    assert p["heavy_atoms"] == 13 and p["atoms"] == 21 and p["bonds"] == 13
    assert p["charge"] == 0
    assert p["rings"] == 1 and p["aromatic_rings"] == 1
    assert p["stereocenters"] == {"assigned": 0, "unassigned": 0, "atoms": [], "labels": {}}
    assert p["hbd"] == 1 and p["hba"] == 3
    assert p["rotatable_bonds"] == 2
    assert p["lipinski_violations"] == 0
    assert abs(p["tpsa"] - 63.6) < 0.05
    assert 0 < p["qed"] < 1
    assert 1.0 < p["logp"] < 1.5
    assert p["heavy_atom_count"] == 13
    assert 0 < p["fraction_csp3"] < 0.2


def test_properties_basic_only(aspirin):
    p = chem.properties(aspirin, extended=False)
    assert "logp" not in p and "formula" in p


def test_properties_charge_and_stereo():
    assert chem.properties(chem.mol_from_input("[NH4+]"))["charge"] == 1
    assert chem.properties(chem.mol_from_input("[NH4+]"))["formula"] == "H4N+"
    unassigned = chem.properties(chem.mol_from_input("CC(N)C(=O)O"))["stereocenters"]
    assert unassigned["assigned"] == 0 and unassigned["unassigned"] == 1 and unassigned["atoms"] == [1]


def test_properties_accepts_string():
    assert chem.properties("CCO")["formula"] == "C2H6O"


def test_elemental_analysis(aspirin):
    assert chem.elemental_analysis(aspirin) == "C, 60.00; H, 4.48; O, 35.52"
    assert chem.elemental_analysis(chem.mol_from_input("O")) == "H, 11.19; O, 88.81"
    assert chem.elemental_analysis(chem.mol_from_input("[Na+].[Cl-]")) == "Cl, 60.66; Na, 39.34"


def test_describe(aspirin):
    text = chem.describe(aspirin)
    assert "\n" not in text
    assert "C9H8O4" in text and "180.159" in text and canon(ASPIRIN) in text
    assert "no stereocentres" in text
    ala = chem.describe(chem.mol_from_input(S_ALANINE))
    assert "1 assigned" in ala
    named = chem.mol_from_input("CCO ethanol")
    assert chem.describe(named).startswith("ethanol: ")


# --------------------------------------------------------------------------- #
# reactions
# --------------------------------------------------------------------------- #
def test_parse_reaction_smiles():
    rxn = chem.parse_reaction(ESTERIFICATION)
    assert isinstance(rxn, Reaction)
    assert [canon_mol(m) for m in rxn.reactants] == ["CC(=O)O", "CCO"]
    assert [canon_mol(m) for m in rxn.agents] == ["[H+]"]
    assert [canon_mol(m) for m in rxn.products] == ["CCOC(C)=O", "O"]
    assert rxn.conditions == "" and rxn.name == ""
    for m in rxn.components:
        assert m.GetNumConformers() == 1
        assert not any(b.GetIsAromatic() for b in m.GetBonds())
    assert chem.reaction_smiles(rxn) == "CC(=O)O.CCO>[H+]>CCOC(C)=O.O"


def test_parse_reaction_smiles_no_agents_and_grouping():
    rxn = chem.parse_reaction("CCO>>CC=O")
    assert len(rxn.reactants) == 1 and not rxn.agents and len(rxn.products) == 1
    grouped = chem.parse_reaction("([Na+].[OH-]).CC(=O)O>>CC(=O)[O-].[Na+].O")
    assert len(grouped.reactants) == 2
    assert len(Chem.GetMolFrags(grouped.reactants[0])) == 2


def test_parse_reaction_arrow_syntax():
    rxn = chem.parse_reaction("CC(=O)O + OCC -> CC(=O)OCC + O")
    assert chem.reaction_smiles(rxn) == "CC(=O)O.CCO>>CCOC(C)=O.O"
    assert not rxn.agents
    for arrow in ("-->", "→", ">>", "=", "=>"):
        r = chem.parse_reaction(f"CCO {arrow} CC=O")
        assert chem.reaction_smiles(r) == "CCO>>CC=O"
    cond = chem.parse_reaction("CC(=O)O + OCC --[H2SO4, reflux]--> CC(=O)OCC + O")
    assert cond.conditions == "H2SO4, reflux"
    assert chem.reaction_smiles(cond) == "CC(=O)O.CCO>>CCOC(C)=O.O"
    charged = chem.parse_reaction("[NH4+] + [OH-] -> N + O")  # '+' inside brackets survives
    assert len(charged.reactants) == 2


def test_parse_reaction_with_names_and_coefficients():
    rxn = chem.parse_reaction("acetic acid + ethanol -> ethyl acetate + water")
    assert chem.reaction_smiles(rxn) == "CC(=O)O.CCO>>CCOC(C)=O.O"
    assert [m.GetProp("_Name") for m in rxn.components] == [
        "acetic acid", "ethanol", "ethyl acetate", "water"
    ]
    burn = chem.parse_reaction("2 H2 + O2 -> 2 H2O")
    assert len(burn.reactants) == 3 and len(burn.products) == 2
    assert chem.reaction_properties(burn)["balanced"] is True


def test_parse_reaction_failures():
    with pytest.raises(ChemError) as exc:
        chem.parse_reaction("foo + CCO -> bar")
    assert "'foo'" in str(exc.value) and "'bar'" in str(exc.value)
    with pytest.raises(ChemError, match="arrow"):
        chem.parse_reaction("CCO CC=O")
    with pytest.raises(ChemError, match="Multi-step"):
        chem.parse_reaction("CCO -> CC=O -> CC(=O)O")
    with pytest.raises(ChemError):
        chem.parse_reaction("   ")


def test_parse_reaction_rxn_block_roundtrip(tmp_path):
    src = rdChemReactions.ReactionFromSmarts("CC(=O)O.OCC>>CC(=O)OCC.O", useSmiles=True)
    block = rdChemReactions.ReactionToRxnBlock(src)
    assert chem.detect_format(block) == "rxn"
    rxn = chem.parse_reaction(block)
    assert chem.reaction_smiles(rxn) == "CC(=O)O.CCO>>CCOC(C)=O.O"
    assert all(m.GetNumConformers() == 1 for m in rxn.components)
    assert all(a.GetTotalNumHs() >= 0 for m in rxn.components for a in m.GetAtoms())
    lines = block.splitlines()
    lines[1] = "esterification"
    path = tmp_path / "ester.rxn"
    path.write_text("\n".join(lines) + "\n")
    from_file = chem.parse_reaction(str(path))
    assert from_file.name == "esterification"
    assert chem.reaction_smiles(from_file) == chem.reaction_smiles(rxn)
    assert chem.parse_reaction("\n".join(lines), fmt="rxn").name == "esterification"


def test_parse_reaction_cdxml_scheme():
    if not hasattr(rdChemReactions, "ReactionsFromCDXMLFile"):
        pytest.skip("RDKit build cannot read CDXML reactions")
    rxn = chem.parse_reaction(CDXML_SCHEME)
    assert chem.reaction_smiles(rxn) == "CO>>CNC"


def test_reaction_properties_balanced():
    rxn = chem.parse_reaction(ESTERIFICATION)
    rp = chem.reaction_properties(rxn)
    assert rp["balanced"] is True
    assert rp["element_balance"] == {"C": (4, 4), "H": (10, 10), "O": (3, 3)}
    assert rp["unbalanced_elements"] == []
    assert abs(rp["mass_reactants"] - 106.121) < 0.01
    assert abs(rp["mass_products"] - 106.121) < 0.01
    assert abs(rp["atom_economy"] - 83.0) < 0.2  # ethyl acetate / (acetic acid + ethanol)
    assert [c["formula"] for c in rp["reactants"]] == ["C2H4O2", "C2H6O"]
    assert [c["formula"] for c in rp["products"]] == ["C4H8O2", "H2O"]
    assert rp["agents"][0]["formula"] == "H+"
    assert rp["reaction_smiles"] == chem.reaction_smiles(rxn)


def test_reaction_properties_unbalanced():
    rp = chem.reaction_properties(chem.parse_reaction("CCO>>CC=O"))
    assert rp["balanced"] is False
    assert rp["element_balance"]["H"] == (6, 4)
    assert rp["unbalanced_elements"] == ["H"]
    assert rp["atom_economy"] is None
    assert rp["mass_reactants"] > rp["mass_products"] > 0


def test_reaction_dataclass_defaults():
    r = Reaction()
    assert r.reactants == [] and r.agents == [] and r.products == []
    assert r.name == "" and r.conditions == "" and r.components == []


# --------------------------------------------------------------------------- #
# depiction and export
# --------------------------------------------------------------------------- #
def test_depict_png(aspirin):
    png = chem.depict_png(aspirin, width=300, height=200, legend="aspirin")
    assert isinstance(png, bytes) and png.startswith(b"\x89PNG\r\n\x1a\n")
    png2 = chem.depict_png(aspirin, width=300, height=200, bond_length_px=20)
    assert png2.startswith(b"\x89PNG")
    rxn_png = chem.depict_png(chem.parse_reaction(ESTERIFICATION), width=600, height=200)
    assert rxn_png.startswith(b"\x89PNG")
    assert chem.depict_png("CCO", width=100, height=100).startswith(b"\x89PNG")


def test_depict_svg(aspirin):
    svg = chem.depict_svg(aspirin, width=300, height=200, legend="aspirin")
    assert "<svg" in svg and "</svg>" in svg
    assert len(svg) > len(chem.depict_svg(aspirin, width=300, height=200))  # legend glyphs drawn
    assert "<svg" in chem.depict_svg(chem.parse_reaction("CCO>>CC=O"))
    assert "<svg" in chem.depict_svg("CCO -> CC=O")


def test_mol_to_molblock_roundtrip():
    m = chem.mol_from_input(S_ALANINE)
    block = chem.mol_to_molblock(m)
    assert "M  END" in block and "V2000" in block
    back = Chem.MolFromMolBlock(block)
    assert Chem.MolToSmiles(back) == canon(S_ALANINE)
    v3 = chem.mol_to_molblock(m, v3000=True)
    assert "V3000" in v3 and Chem.MolToSmiles(Chem.MolFromMolBlock(v3)) == canon(S_ALANINE)
    assert "M  END" in chem.mol_to_molblock(Chem.MolFromSmiles("CCO"))  # no coords yet


def test_mol_to_sdf():
    a = chem.mol_from_input("CCO ethanol")
    b = chem.mol_from_input(ASPIRIN)
    b.SetProp("ID", "42")
    text = chem.mol_to_sdf([a, b])
    assert text.count("$$$$") == 2
    assert text.startswith("ethanol")
    assert "<ID>" in text and "42" in text
    back = chem.mols_from_input(text)
    assert [canon_mol(m) for m in back] == ["CCO", canon(ASPIRIN)]
    assert chem.mol_to_sdf(a).count("$$$$") == 1


def test_substructure_matches(aspirin):
    assert chem.substructure_matches(aspirin, "c1ccccc1") == [[4, 5, 6, 7, 8, 9]]
    assert chem.substructure_matches(aspirin, "C(=O)O") == [[1, 2, 3], [10, 11, 12]]
    assert chem.substructure_matches(aspirin, "[Cl]") == []
    with pytest.raises(ChemError, match="SMARTS"):
        chem.substructure_matches(aspirin, "C1CC(")


def test_tanimoto(aspirin):
    assert chem.tanimoto(aspirin, aspirin) == pytest.approx(1.0)
    sim = chem.tanimoto("CCO", "CCCO")
    assert 0.0 < sim < 1.0
    assert chem.tanimoto("CCO", "CCCO", radius=1) >= sim
    assert chem.tanimoto(aspirin, "c1ccccc1") < 0.5


# --------------------------------------------------------------------------- #
# resolve_name
# --------------------------------------------------------------------------- #
def test_resolve_name_builtin():
    assert canon(chem.resolve_name("aspirin")) == canon(ASPIRIN)
    assert canon(chem.resolve_name("  Aspirin. ")) == canon(ASPIRIN)
    assert chem.resolve_name("sodium borohydride") == "[Na+].[BH4-]"
    assert canon(chem.resolve_name("L-Alanine")) == canon(S_ALANINE)
    assert canon(chem.resolve_name("acetaminophen")) == canon(chem.resolve_name("paracetamol"))
    assert chem.resolve_name("Hünig's base") == chem.resolve_name("DIPEA")
    assert chem.resolve_name("") is None


def test_builtin_dictionary_is_valid():
    assert len(chem.COMMON_NAMES) >= 60
    with rdBase.BlockLogs():  # [Na+].[H-] logs a harmless RemoveHs warning
        for name, smiles in chem.COMMON_NAMES.items():
            assert Chem.MolFromSmiles(smiles) is not None, name
    for name, label in [("alanine", "S"), ("cysteine", "R"), ("(s)-ibuprofen", "S"), ("(r)-ibuprofen", "R")]:
        m = Chem.MolFromSmiles(chem.COMMON_NAMES[name])
        assert [lab for _, lab in Chem.FindMolChiralCenters(m, useLegacyImplementation=False)] == [label]


def test_resolve_name_unknown_offline_no_network(monkeypatch):
    calls = []

    def fake_urlopen(*args, **kwargs):  # pragma: no cover - must not run
        calls.append(args)
        raise AssertionError("network must not be used offline")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    assert chem.resolve_name("unobtainium") is None
    assert chem.resolve_name("unobtainium", online=False) is None
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "0")
    assert chem.resolve_name("unobtainium") is None
    assert calls == []


def test_resolve_name_online_uses_pubchem(monkeypatch):
    seen = []

    class FakeResponse(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def fake_urlopen(url, timeout=None):
        seen.append((url, timeout))
        return FakeResponse(b"CC(=O)Oc1ccccc1C(=O)O\n")

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    assert canon(chem.resolve_name("acetylsalicylic acid online", online=True)) == canon(ASPIRIN)
    assert len(seen) == 1
    assert "pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" in seen[0][0]
    assert "acetylsalicylic%20acid%20online" in seen[0][0]
    assert seen[0][1] == 10
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "1")
    assert chem.resolve_name("acetylsalicylic acid online") is not None
    assert len(seen) == 2


def test_resolve_name_online_network_error_returns_none(monkeypatch):
    def failing_urlopen(url, timeout=None):
        raise OSError("no network")

    monkeypatch.setattr(urllib.request, "urlopen", failing_urlopen)
    assert chem.resolve_name("unobtainium", online=True) is None


def test_resolve_name_opsin_cli(monkeypatch):
    """An ``opsin`` executable on PATH is used before PubChem."""
    monkeypatch.setattr(chem.shutil, "which", lambda name, *a, **k: "/usr/bin/opsin" if name == "opsin" else None)

    class Proc:
        stdout = "CCO\n"
        returncode = 0

    monkeypatch.setattr(chem.subprocess, "run", lambda *a, **k: Proc())
    assert chem.resolve_name("totally-unknown-name") == "CCO"


# --------------------------------------------------------------------------- #
# regressions: all-hydrogen species and dative bonds
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("smiles", ["[H][H]", "[H+]", "[H-]", "[2H][2H]", "[H]"])
def test_all_hydrogen_species_get_finite_2d_coordinates(smiles):
    m = chem.mol_from_input(smiles)
    assert m.GetNumConformers() == 1
    pos = np.asarray(m.GetConformer().GetPositions())
    assert pos.shape[0] == m.GetNumAtoms()
    assert np.all(np.isfinite(pos))
    assert not m.GetConformer().Is3D()
    if m.GetNumAtoms() == 1:
        assert np.allclose(pos[0], [0.0, 0.0, 0.0])
    else:
        assert m.GetNumAtoms() == 2
        assert np.allclose(pos[0], [0.0, 0.0, 0.0])
        assert pos[1][0] == pytest.approx(1.5) and pos[1][1] == pytest.approx(0.0)
    prepared = chem.prepare_for_depiction(m)  # idempotent, still finite
    assert np.all(np.isfinite(np.asarray(prepared.GetConformer().GetPositions())))
    assert "<svg" in chem.depict_svg(m)
    assert chem.depict_png(m)[:8] == b"\x89PNG\r\n\x1a\n"


def test_fallback_coords_keep_fragments_apart():
    m = chem.mol_from_input("[H+].[H-]")
    pos = np.asarray(m.GetConformer().GetPositions())
    assert np.all(np.isfinite(pos))
    assert abs(pos[0][0] - pos[1][0]) > 1.0


@pytest.mark.parametrize("smiles", ["[NH3]->[Pt]", "C1=CC=CC=C1->[Fe]", "[Fe]<-[NH3]", "[Cu]<-N(C)C->[Cu]"])
def test_dative_bonds_are_smiles_not_reactions(smiles):
    assert chem.detect_format(smiles) == "smiles"
    assert not chem._is_reaction_smiles(smiles)
    m = chem.mol_from_input(smiles)
    assert m.GetNumAtoms() == Chem.MolFromSmiles(smiles).GetNumAtoms()
    assert np.all(np.isfinite(np.asarray(m.GetConformer().GetPositions())))


@pytest.mark.parametrize(
    "text, expected",
    [
        ("CCO>>CC=O", "reaction_smiles"),
        ("CC->[Cu]>>CC", "reaction_smiles"),
        ("[NH3]->[Pt]>>[Pt]", "reaction_smiles"),
        ("CCO -> CC=O", "reaction_smiles"),
        ("CCO --> CC=O", "reaction_smiles"),
    ],
)
def test_reaction_detection_still_works_with_dative_tokens(text, expected):
    assert chem.detect_format(text) == expected


def test_reaction_smiles_with_dative_bond_component():
    rxn = chem.parse_reaction("CC->[Cu]>>CC")
    assert len(rxn.reactants) == 1 and len(rxn.products) == 1 and rxn.agents == []
    assert rxn.reactants[0].GetNumAtoms() == 3
    rxn2 = chem.parse_reaction("[NH3]->[Pt]>[Cu]<-N>[Pt]")
    assert len(rxn2.reactants) == 1 and len(rxn2.agents) == 1 and len(rxn2.products) == 1
    assert rxn2.agents[0].GetNumAtoms() == 2


# --- regression tests added after the review workflow -----------------------------------------------

def test_cxsmiles_extension_block_is_not_a_name():
    mols = chem.mols_from_input("CC(=O)O |c:1| acetic\nCCO ethanol")
    assert [m.GetProp("_Name") for m in mols] == ["acetic", "ethanol"]
    assert mols[0].GetNumAtoms() == 4


def test_malformed_reaction_smiles_raises_chemerror():
    import pytest
    for bad in ("CC->", "CC>>CC>>O", ">>"):
        with pytest.raises(ChemError):
            chem.parse_reaction(bad)


def test_cdx_file_path_is_read_through_the_cdx_reader(tmp_path):
    import os
    pairs = os.environ.get("CHEMDRAWMCP_CDX_PAIRS") or os.path.expanduser("~/.chemdrawmcp/testdata/cdxpairs")
    src = os.path.join(pairs, "Catecholamines.cdx")
    if not os.path.isfile(src):
        import pytest
        pytest.skip("no CDX sample available")
    mols = chem.mols_from_input(src)
    assert len(mols) >= 5 and all(m.GetNumAtoms() > 0 for m in mols)
    assert chem.detect_format(src) == "cdxml"
