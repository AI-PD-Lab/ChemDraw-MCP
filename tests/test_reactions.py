"""Offline tests for chemdrawmcp.reactions (named-reaction templates)."""
from __future__ import annotations

import json

import pytest

from chemdrawmcp import chem, reactions
from chemdrawmcp.reactions import REACTIONS, ReactionError, ReactionTemplate

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402


def canon(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


# --------------------------------------------------------------------------- #
# One worked example per template: reactants (role order) -> expected products
# (canonical isomeric SMILES).  Every key in REACTIONS must appear here.
# --------------------------------------------------------------------------- #
EXAMPLES = {
    "fischer_esterification": (["CC(=O)O", "CCO"], ["CCOC(C)=O"]),
    "steglich_esterification": (["OC(=O)c1ccccc1", "CC(C)(C)O"], ["CC(C)(C)OC(=O)c1ccccc1"]),
    "acid_chloride_alcohol": (["CC(=O)Cl", "OCc1ccccc1"], ["CC(=O)OCc1ccccc1"]),
    "amide_coupling": (["OC(=O)c1ccccc1", "CN"], ["CNC(=O)c1ccccc1"]),
    "acid_chloride_amine": (["CC(=O)Cl", "Nc1ccccc1"], ["CC(=O)Nc1ccccc1"]),
    "acid_chloride_formation": (["OC(=O)c1ccccc1"], ["O=C(Cl)c1ccccc1"]),
    "sulfonamide_formation": (["Cc1ccc(S(=O)(=O)Cl)cc1", "Nc1ccccc1"], ["Cc1ccc(S(=O)(=O)Nc2ccccc2)cc1"]),
    "urea_formation": (["O=C=Nc1ccccc1", "CN"], ["CNC(=O)Nc1ccccc1"]),
    "suzuki_coupling": (["Brc1ccccc1", "OB(O)c1ccccc1"], ["c1ccc(-c2ccccc2)cc1"]),
    "heck_reaction": (["Brc1ccccc1", "C=CC(=O)OC"], ["COC(=O)/C=C/c1ccccc1"]),
    "sonogashira_coupling": (["Ic1ccccc1", "C#Cc1ccccc1"], ["C(#Cc1ccccc1)c1ccccc1"]),
    "negishi_coupling": (["Brc1ccccc1", "CC[Zn]Cl"], ["CCc1ccccc1"]),
    "stille_coupling": (["Ic1ccccc1", "C=C[Sn](CCCC)(CCCC)CCCC"], ["C=Cc1ccccc1"]),
    "grignard_aldehyde": (["O=Cc1ccccc1", "C[Mg]Br"], ["CC(O)c1ccccc1"]),
    "grignard_ketone": (["CC(=O)c1ccccc1", "C[Mg]Br"], ["CC(C)(O)c1ccccc1"]),
    "grignard_ester": (["CCOC(=O)c1ccccc1", "C[Mg]Br"], ["CC(C)(O)c1ccccc1"]),
    "wittig_reaction": (["C=P(c1ccccc1)(c1ccccc1)c1ccccc1", "O=Cc1ccccc1"], ["C=Cc1ccccc1"]),
    "horner_wadsworth_emmons": (["CCOC(=O)CP(=O)(OCC)OCC", "O=Cc1ccccc1"], ["CCOC(=O)/C=C/c1ccccc1"]),
    "cross_metathesis": (["C=CCc1ccccc1", "C=CC(=O)OC"], ["COC(=O)C=CCc1ccccc1"]),
    "ring_closing_metathesis": (["C=CCC(CC=C)(C(=O)OCC)C(=O)OCC"], ["CCOC(=O)C1(C(=O)OCC)CC=CC1"]),
    "aldol_addition": (["CC(C)=O", "O=Cc1ccccc1"], ["CC(=O)CC(O)c1ccccc1"]),
    "aldol_condensation": (["CC(C)=O", "O=Cc1ccccc1"], ["CC(=O)C=Cc1ccccc1"]),
    "claisen_condensation": (["CCOC(C)=O", "CCOC(C)=O"], ["CCOC(=O)CC(C)=O"]),
    "knoevenagel_condensation": (["O=Cc1ccccc1", "N#CCC#N"], ["N#CC(C#N)=Cc1ccccc1"]),
    "michael_addition": (["CCOC(=O)CC(=O)OCC", "C=CC(C)=O"], ["CCOC(=O)C(CCC(C)=O)C(=O)OCC"]),
    "mannich_reaction": (["CC(=O)c1ccccc1", "C=O", "CNC"], ["CN(C)CCC(=O)c1ccccc1"]),
    "robinson_annulation": (["O=C1CCCCC1", "CC(=O)C=C"], ["O=C1C=C2CCCCC2CC1"]),
    "decarboxylation": (["CC(=O)CC(=O)O"], ["CC(C)=O"]),
    "diels_alder": (["C1=CCC=C1", "O=C1OC(=O)C=C1"], ["O=C1OC(=O)C2C3C=CC(C3)C12"]),
    "cuaac_click": (["[N-]=[N+]=NCc1ccccc1", "C#Cc1ccccc1"], ["c1ccc(Cn2cc(-c3ccccc3)nn2)cc1"]),
    "spaac_click": (["[N-]=[N+]=NCc1ccccc1", "C1CCCC#CCC1"], ["c1ccc(Cn2nnc3c2CCCCCC3)cc1"]),
    "sn2_amine": (["BrCc1ccccc1", "CN"], ["CNCc1ccccc1"]),
    "sn2_alkoxide": (["CCBr", "[Na+].CC[O-]"], ["CCOCC"]),
    "sn2_azide": (["BrCc1ccccc1", "[Na+].[N-]=[N+]=[N-]"], ["[N-]=[N+]=NCc1ccccc1"]),
    "sn2_cyanide": (["BrCc1ccccc1", "[Na+].[C-]#N"], ["N#CCc1ccccc1"]),
    "sn2_thiol": (["BrCc1ccccc1", "Sc1ccccc1"], ["c1ccc(CSc2ccccc2)cc1"]),
    "williamson_ether_synthesis": (["Oc1ccccc1", "CI"], ["COc1ccccc1"]),
    "mitsunobu_reaction": (["C[C@H](O)CC", "OC(=O)c1ccccc1"], ["CC[C@@H](C)OC(=O)c1ccccc1"]),
    "appel_reaction": (["OCc1ccccc1"], ["BrCc1ccccc1"]),
    "tosylation": (["CCO"], ["CCOS(=O)(=O)c1ccc(C)cc1"]),
    "mesylation": (["CCO"], ["CCOS(C)(=O)=O"]),
    "e2_elimination": (["CC(C)(C)Br"], ["C=C(C)C"]),
    "friedel_crafts_acylation": (["c1ccccc1", "CC(=O)Cl"], ["CC(=O)c1ccccc1"]),
    "friedel_crafts_alkylation": (["c1ccccc1", "CC(C)Cl"], ["CC(C)c1ccccc1"]),
    "aromatic_nitration": (["c1ccccc1"], ["O=[N+]([O-])c1ccccc1"]),
    "aromatic_bromination": (["c1ccccc1"], ["Brc1ccccc1"]),
    "nitro_reduction": (["O=[N+]([O-])c1ccccc1"], ["Nc1ccccc1"]),
    "alkene_hydrogenation": (["C=Cc1ccccc1"], ["CCc1ccccc1"]),
    "alkyne_hydrogenation": (["c1ccc(C#Cc2ccccc2)cc1"], ["c1ccc(CCc2ccccc2)cc1"]),
    "lindlar_hydrogenation": (["CC#Cc1ccccc1"], ["C/C=C\\c1ccccc1"]),
    "carbonyl_reduction": (["CC(=O)c1ccccc1"], ["CC(O)c1ccccc1"]),
    "ester_reduction": (["COC(=O)c1ccccc1"], ["OCc1ccccc1", "CO"]),
    "amide_reduction": (["CNC(=O)c1ccccc1"], ["CNCc1ccccc1"]),
    "nitrile_reduction": (["N#Cc1ccccc1"], ["NCc1ccccc1"]),
    "staudinger_reduction": (["[N-]=[N+]=NCc1ccccc1"], ["NCc1ccccc1"]),
    "wolff_kishner_reduction": (["CC(=O)c1ccccc1"], ["CCc1ccccc1"]),
    "clemmensen_reduction": (["CC(=O)c1ccccc1"], ["CCc1ccccc1"]),
    "birch_reduction": (["COc1ccccc1"], ["COC1=CCC=CC1"]),
    "swern_oxidation": (["OCc1ccccc1"], ["O=Cc1ccccc1"]),
    "dess_martin_oxidation": (["CC(O)c1ccccc1"], ["CC(=O)c1ccccc1"]),
    "secondary_alcohol_oxidation": (["OC1CCCCC1"], ["O=C1CCCCC1"]),
    "jones_oxidation": (["OCc1ccccc1"], ["O=C(O)c1ccccc1"]),
    "pinnick_oxidation": (["O=Cc1ccccc1"], ["O=C(O)c1ccccc1"]),
    "baeyer_villiger_oxidation": (["O=C1CCCCC1"], ["O=C1CCCCCO1"]),
    "ozonolysis_reductive": (["CC1=CCCCC1"], ["CC(=O)CCCCC=O"]),
    "epoxidation": (["C1CCC=CC1"], ["C1CCC2OC2C1"]),
    "sharpless_dihydroxylation": (["C=Cc1ccccc1"], ["OCC(O)c1ccccc1"]),
    "alkene_dibromination": (["C1CCC=CC1"], ["BrC1CCCCC1Br"]),
    "hydroboration_oxidation": (["C=Cc1ccccc1"], ["OCCc1ccccc1"]),
    "markovnikov_hydration": (["C=Cc1ccccc1"], ["CC(O)c1ccccc1"]),
    "epoxide_opening_amine": (["C1OC1c1ccccc1", "CN"], ["CNCC(O)c1ccccc1"]),
    "buchwald_hartwig_amination": (["Brc1ccccc1", "C1COCCN1"], ["c1ccc(N2CCOCC2)cc1"]),
    "reductive_amination": (["O=Cc1ccccc1", "CN"], ["CNCc1ccccc1"]),
    "imine_formation": (["O=Cc1ccccc1", "Nc1ccccc1"], ["C(=Nc1ccccc1)c1ccccc1"]),
    "oxime_formation": (["CC(C)=O", "NO"], ["CC(C)=NO"]),
    "boc_protection": (["Nc1ccccc1"], ["CC(C)(C)OC(=O)Nc1ccccc1"]),
    "boc_deprotection": (["CC(C)(C)OC(=O)Nc1ccccc1"], ["Nc1ccccc1"]),
    "fmoc_protection": (["NCC(=O)O"], ["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc21"]),
    "fmoc_deprotection": (["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc12"], ["NCC(=O)O"]),
    "tbs_protection": (["OCc1ccccc1"], ["CC(C)(C)[Si](C)(C)OCc1ccccc1"]),
    "tbs_deprotection": (["C[Si](C)(OCc1ccccc1)C(C)(C)C"], ["OCc1ccccc1"]),
    "benzyl_protection": (["OCCO"], ["OCCOCc1ccccc1"]),
    "benzyl_deprotection": (["COc1ccc(OCc2ccccc2)cc1"], ["COc1ccc(O)cc1"]),
    "acetylation_alcohol": (["OCc1ccccc1"], ["CC(=O)OCc1ccccc1"]),
    "acetal_formation": (["O=C1CCCCC1", "OCCO"], ["C1CCC2(CC1)OCCO2"]),
    "acetal_formation_dialkyl": (["O=Cc1ccccc1", "CO"], ["COC(OC)c1ccccc1"]),
    "acetal_hydrolysis": (["COC(OC)c1ccccc1"], ["O=Cc1ccccc1"]),
    "ester_hydrolysis": (["CCOC(C)=O"], ["CC(=O)O", "CCO"]),
    "amide_hydrolysis": (["CC(=O)Nc1ccccc1"], ["CC(=O)O", "Nc1ccccc1"]),
    "nitrile_hydrolysis": (["N#Cc1ccccc1"], ["O=C(O)c1ccccc1"]),
    "curtius_rearrangement": (["[N-]=[N+]=NC(=O)c1ccccc1"], ["O=C=Nc1ccccc1"]),
    "isocyanate_hydrolysis": (["O=C=Nc1ccccc1"], ["Nc1ccccc1"]),
    "curtius_to_amine": (["OC(=O)c1ccccc1"], ["Nc1ccccc1"]),
    "hofmann_rearrangement": (["NC(=O)c1ccccc1"], ["Nc1ccccc1"]),
    "fischer_indole_synthesis": (["NNc1ccccc1", "O=C1CCCCC1"], ["c1ccc2c3c([nH]c2c1)CCCC3"]),
    "paal_knorr_pyrrole": (["CC(=O)CCC(C)=O", "NCc1ccccc1"], ["Cc1ccc(C)n1Cc1ccccc1"]),
    "paal_knorr_furan": (["CC(=O)CCC(C)=O"], ["Cc1ccc(C)o1"]),
    "biginelli_reaction": (["O=Cc1ccccc1", "CCOC(=O)CC(C)=O", "NC(N)=O"], ["CCOC(=O)C1=C(C)NC(=O)NC1c1ccccc1"]),
    "ugi_reaction": (["O=Cc1ccccc1", "CN", "CC(=O)O", "CC(C)(C)[N+]#[C-]"],
                     ["CC(=O)N(C)C(C(=O)NC(C)(C)C)c1ccccc1"]),
    "passerini_reaction": (["O=Cc1ccccc1", "CC(=O)O", "CC(C)(C)[N+]#[C-]"], ["CC(=O)OC(C(=O)NC(C)(C)C)c1ccccc1"]),
}

# templates the module must provide (task list)
REQUIRED = [
    "fischer_esterification", "amide_coupling", "acid_chloride_amine", "suzuki_coupling", "heck_reaction",
    "sonogashira_coupling", "buchwald_hartwig_amination", "negishi_coupling", "stille_coupling",
    "grignard_aldehyde", "grignard_ketone", "grignard_ester", "wittig_reaction", "horner_wadsworth_emmons",
    "aldol_addition", "claisen_condensation", "diels_alder", "sn2_amine", "sn2_alkoxide", "sn2_azide",
    "sn2_cyanide", "sn2_thiol", "e2_elimination", "williamson_ether_synthesis", "friedel_crafts_acylation",
    "friedel_crafts_alkylation", "aromatic_nitration", "aromatic_bromination", "nitro_reduction",
    "alkene_hydrogenation", "alkyne_hydrogenation", "lindlar_hydrogenation", "swern_oxidation",
    "dess_martin_oxidation", "jones_oxidation", "secondary_alcohol_oxidation", "baeyer_villiger_oxidation",
    "ozonolysis_reductive", "hydroboration_oxidation", "epoxidation", "epoxide_opening_amine",
    "reductive_amination", "boc_protection", "boc_deprotection", "fmoc_deprotection", "tbs_protection",
    "tbs_deprotection", "acetal_formation", "ester_hydrolysis", "amide_hydrolysis", "nitrile_hydrolysis",
    "nitrile_reduction", "curtius_rearrangement", "hofmann_rearrangement", "mitsunobu_reaction",
    "cuaac_click", "spaac_click", "staudinger_reduction", "michael_addition", "knoevenagel_condensation",
    "mannich_reaction", "robinson_annulation", "fischer_indole_synthesis", "paal_knorr_pyrrole",
    "biginelli_reaction", "ugi_reaction", "passerini_reaction", "cross_metathesis",
    "ring_closing_metathesis", "sharpless_dihydroxylation", "wolff_kishner_reduction",
    "clemmensen_reduction", "birch_reduction",
]

KEYS = list(REACTIONS)


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
def test_registry_size_and_required_keys():
    assert len(REACTIONS) >= 45
    missing = [k for k in REQUIRED if k not in REACTIONS]
    assert not missing, missing
    assert set(EXAMPLES) == set(REACTIONS), (set(EXAMPLES) ^ set(REACTIONS))


def test_registry_metadata_complete():
    norm_labels = {}
    for key, t in REACTIONS.items():
        assert isinstance(t, ReactionTemplate)
        assert t.key == key and key == key.lower() and " " not in key and "-" not in key
        assert t.name and t.cls and t.description and t.roles and t.smarts
        assert all(isinstance(x, str) and x for x in [*t.reagents, *t.conditions, *t.aliases, *t.byproducts])
        assert len(t.equivalents) == len(t.roles) and all(e >= 1 for e in t.equivalents)
        assert len(t.example_reactants) == t.n_reactants and t.example_products
        # key / name / alias labels never collide between templates
        for label in (key, t.name, *t.aliases):
            n = reactions._norm(label)
            assert norm_labels.setdefault(n, key) == key, (label, key, norm_labels[n])


def test_classes_and_as_dict_json():
    cls = reactions.classes()
    assert "oxidation" in cls and "reduction" in cls and "C–C coupling" in cls
    for t in REACTIONS.values():
        d = t.as_dict()
        assert d["class"] == t.cls and d["key"] == t.key and d["n_reactants"] == t.n_reactants
        json.dumps(d)


@pytest.mark.parametrize("key", KEYS)
def test_smarts_parses(key):
    t = REACTIONS[key]
    rxn = t.rxn
    assert rxn is not None
    assert rxn.GetNumReactantTemplates() == sum(t.equivalents)
    assert rxn.GetNumProductTemplates() >= 1
    for role in range(t.n_reactants):
        assert t.role_template(role).GetNumAtoms() >= 1
    # the worked example fits every role in the stated order
    for role, smi in enumerate(t.example_reactants):
        assert t.matches_role(Chem.MolFromSmiles(smi), role), (key, role, smi)


@pytest.mark.parametrize("key", KEYS)
def test_worked_example(key):
    reactants, expected = EXAMPLES[key]
    t = REACTIONS[key]
    # the module's own example agrees with the hard-coded table
    assert list(t.example_reactants) == reactants
    assert sorted(canon(s) for s in t.example_products) == sorted(canon(s) for s in expected)
    results = reactions.apply_reaction(key, reactants)
    assert results, key
    want = sorted(canon(s) for s in expected)
    hits = [r for r in results if sorted(r["products"]) == want]
    assert hits, (key, [r["products"] for r in results])
    r = hits[0]
    assert r["template"] == key and r["name"] == t.name and r["class"] == t.cls
    assert r["reagents"] == list(t.reagents) and r["conditions"] == list(t.conditions)
    assert len(r["reactants"]) == t.n_reactants and len(r["roles"]) == t.n_reactants
    assert sorted(r["reactant_order"]) == list(range(t.n_reactants))
    # mapped reaction SMILES: parses, has sum(equivalents) reactants, map numbers consistent
    mapped = r["mapped"]
    assert ">>" in mapped
    rxn = AllChem.ReactionFromSmarts(mapped, useSmiles=True)
    assert rxn.GetNumReactantTemplates() == sum(t.equivalents)
    assert rxn.GetNumProductTemplates() == len(r["products"])
    r_maps = [a.GetAtomMapNum() for m in rxn.GetReactants() for a in m.GetAtoms() if a.GetAtomMapNum()]
    p_maps = [a.GetAtomMapNum() for m in rxn.GetProducts() for a in m.GetAtoms() if a.GetAtomMapNum()]
    assert p_maps and len(set(p_maps)) == len(p_maps) and len(set(r_maps)) == len(r_maps)
    assert set(p_maps) == set(r_maps)
    # product SMILES in the mapped string match the reported products
    unmapped = []
    for m in rxn.GetProducts():
        mm = Chem.Mol(m)
        for a in mm.GetAtoms():
            a.SetAtomMapNum(0)
        Chem.SanitizeMol(mm)
        unmapped.append(Chem.MolToSmiles(mm))
    assert sorted(unmapped) == want


# --------------------------------------------------------------------------- #
# Lookup
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("query, key", [
    ("suzuki_coupling", "suzuki_coupling"),
    ("Suzuki–Miyaura coupling", "suzuki_coupling"),
    ("Suzuki-Miyaura", "suzuki_coupling"),
    ("suzuki miyaura", "suzuki_coupling"),
    ("SUZUKI", "suzuki_coupling"),
    ("  Suzuki Coupling ", "suzuki_coupling"),
    ("suzuk", "suzuki_coupling"),                      # unique substring
    ("Buchwald–Hartwig", "buchwald_hartwig_amination"),
    ("buchwald", "buchwald_hartwig_amination"),
    ("click", "cuaac_click"),
    ("CuAAC", "cuaac_click"),
    ("HWE", "horner_wadsworth_emmons"),
    ("Horner-Wadsworth-Emmons", "horner_wadsworth_emmons"),
    ("fischer", "fischer_esterification"),
    ("Fischer indole", "fischer_indole_synthesis"),
    ("Boc", "boc_protection"),
    ("boc deprotection", "boc_deprotection"),
    ("Diels-Alder", "diels_alder"),
    ("[4+2]", "diels_alder"),
    ("Wittig", "wittig_reaction"),
    ("RCM", "ring_closing_metathesis"),
    ("saponification", "ester_hydrolysis"),
    ("TBAF", "tbs_deprotection"),
    ("SN2 amine", "sn2_amine"),
    ("sn2-azide", "sn2_azide"),
    ("E2", "e2_elimination"),
    ("nitration", "aromatic_nitration"),
    ("Dess Martin", "dess_martin_oxidation"),
    ("DMP", "dess_martin_oxidation"),
    ("Ugi", "ugi_reaction"),
    ("Staudinger", "staudinger_reduction"),
])
def test_get_reaction_fuzzy(query, key):
    assert reactions.get_reaction(query).key == key


def test_get_reaction_accepts_template_and_rejects_unknown():
    t = REACTIONS["heck_reaction"]
    assert reactions.get_reaction(t) is t
    with pytest.raises(ReactionError, match="Unknown reaction"):
        reactions.get_reaction("xyzzy plugh")
    with pytest.raises(ReactionError, match="did you mean"):
        reactions.get_reaction("oxidation")         # ambiguous: several oxidations listed
    with pytest.raises(ReactionError):
        reactions.get_reaction("")
    with pytest.raises(ReactionError):
        reactions.get_reaction("!!!")


def test_get_reaction_generic_alias_does_not_hijack_query():
    """A short alias contained *in* the query must not win the substring step.

    Regression: ``"bromination"`` (alias of aromatic_bromination) used to match
    ``"alkene bromination"`` because labels were also searched as substrings of
    the query, so the ring was brominated instead of the C=C bond.
    """
    assert reactions.get_reaction("alkene bromination").key == "alkene_dibromination"
    assert reactions.apply_reaction("alkene bromination", ["C=Cc1ccccc1"])[0]["products"] == ["BrCC(Br)c1ccccc1"]
    # queries that merely contain a generic alias plus other words fail loudly
    # with a hint instead of silently picking the aliased template
    for query in ("acetylation of amine", "ketone alpha bromination", "benzylic bromination"):
        with pytest.raises(ReactionError, match="Unknown reaction.*did you mean"):
            reactions.get_reaction(query)
    with pytest.raises(ReactionError, match="acid_chloride_amine"):   # N-acylation is suggested
        reactions.get_reaction("acetylation of amine")
    # the forward direction (query contained in a label) still works
    assert reactions.get_reaction("suzuk").key == "suzuki_coupling"
    assert reactions.get_reaction("bromination").key == "aromatic_bromination"   # exact alias, unchanged


def test_list_reactions_filters():
    everything = reactions.list_reactions()
    assert len(everything) == len(REACTIONS)
    assert all({"key", "name", "class", "roles", "reagents", "conditions", "smarts"} <= set(d) for d in everything)
    ox = reactions.list_reactions(cls="oxidation")
    assert ox and all(d["class"] == "oxidation" for d in ox)
    assert "swern_oxidation" in {d["key"] for d in ox}
    cc = reactions.list_reactions(cls="C-C coupling")   # ASCII dash matches the en dash
    assert {"suzuki_coupling", "heck_reaction", "sonogashira_coupling"} <= {d["key"] for d in cc}
    pd = {d["key"] for d in reactions.list_reactions(query="Pd")}
    assert "suzuki_coupling" in pd and "fischer_esterification" not in pd
    assert [d["key"] for d in reactions.list_reactions(query="grubbs")] == [
        "cross_metathesis", "ring_closing_metathesis"]
    assert reactions.list_reactions(query="no such reaction anywhere") == []
    json.dumps(everything)


# --------------------------------------------------------------------------- #
# apply_reaction behaviour
# --------------------------------------------------------------------------- #
def test_apply_reaction_any_order_names_and_mols():
    forward = reactions.apply_reaction("suzuki", ["Brc1ccccc1", "OB(O)c1ccccc1"])
    reverse = reactions.apply_reaction("suzuki", ["OB(O)c1ccccc1", "Brc1ccccc1"])
    assert forward[0]["products"] == reverse[0]["products"] == ["c1ccc(-c2ccccc2)cc1"]
    assert forward[0]["reactant_order"] == [0, 1] and reverse[0]["reactant_order"] == [1, 0]
    assert forward[0]["roles"] == {"aryl/vinyl halide": "Brc1ccccc1", "boronic acid/ester": "OB(O)c1ccccc1"}
    # names resolve through chem.resolve_name (offline dictionary)
    assert reactions.apply_reaction("fischer", ["acetic acid", "ethanol"])[0]["products"] == ["CCOC(C)=O"]
    # raw and depiction-prepared (Kekulé) RDKit mols are accepted; input map numbers are ignored
    raw = Chem.MolFromSmiles("Brc1ccccc1")
    prepared = chem.mol_from_input("OB(O)[c:7]1ccccc1")
    assert reactions.apply_reaction("suzuki", [raw, prepared])[0]["products"] == ["c1ccc(-c2ccccc2)cc1"]
    assert reactions.apply_reaction("suzuki", [raw, prepared])[0]["reactants"] == ["Brc1ccccc1", "OB(O)c1ccccc1"]


def test_apply_reaction_isomers_dedupe_and_limit():
    res = reactions.apply_reaction("nitration", ["Cc1ccccc1"])
    prods = sorted(r["products"][0] for r in res)
    assert prods == sorted(canon(s) for s in ["Cc1ccccc1[N+](=O)[O-]", "Cc1cccc([N+](=O)[O-])c1",
                                              "Cc1ccc([N+](=O)[O-])cc1"])
    assert len(reactions.apply_reaction("nitration", ["Cc1ccccc1"], max_products=1)) == 1
    # symmetric substrates collapse to one outcome
    assert len(reactions.apply_reaction("nitration", ["c1ccccc1"])) == 1
    assert len(reactions.apply_reaction("claisen", ["CCOC(C)=O", "CCOC(C)=O"])) == 1
    # both E2 regioisomers are reported
    e2 = {r["products"][0] for r in reactions.apply_reaction("e2", ["CCC(C)Br"])}
    assert e2 == {"C=CCC", "CC=CC"}


def test_apply_reaction_ring_opening_gives_single_molecule():
    assert reactions.apply_reaction("saponification", ["O=C1CCCCO1"])[0]["products"] == ["O=C(O)CCCCO"]
    assert reactions.apply_reaction("amide hydrolysis", ["O=C1CCCCN1"])[0]["products"] == ["NCCCCC(=O)O"]
    assert reactions.apply_reaction("ozonolysis", ["C1CCC=CC1"])[0]["products"] == ["O=CCCCCC=O"]
    assert sorted(reactions.apply_reaction("saponification", ["CCOC(C)=O"])[0]["products"]) == ["CC(=O)O", "CCO"]


def test_apply_reaction_stereochemistry():
    assert reactions.apply_reaction("lindlar", ["CCC#CCC"])[0]["products"] == [canon("CC/C=C\\CC")]
    assert reactions.apply_reaction("heck", ["Ic1ccccc1", "C=Cc1ccccc1"])[0]["products"] == [
        canon("c1ccccc1/C=C/c1ccccc1")]
    # SN2 and Mitsunobu invert: (R)-2-bromobutane -> (S)-2-azidobutane, (S)-2-butanol -> (R)-ester
    r_br = "C[C@@H](Br)CC"
    assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(r_br)) == [(1, "R")]
    azide = reactions.apply_reaction("sn2 azide", [r_br, "[N-]=[N+]=[N-]"])[0]["products"][0]
    assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(azide))[0][1] == "S"
    ester = reactions.apply_reaction("mitsunobu", ["C[C@H](O)CC", "OC(=O)c1ccccc1"])[0]["products"][0]
    assert Chem.FindMolChiralCenters(Chem.MolFromSmiles(ester))[0][1] == "R"


def test_apply_reaction_stoichiometry_expands_reactant():
    r = reactions.apply_reaction("grignard ester", ["CCOC(=O)c1ccccc1", "C[Mg]Br"])[0]
    assert r["products"] == ["CC(C)(O)c1ccccc1"]
    assert r["mapped"].count(">>") == 1 and r["mapped"].split(">>")[0].count("[Mg]") == 2
    assert len(r["reactants"]) == 2
    r = reactions.apply_reaction("dimethyl acetal", ["O=Cc1ccccc1", "CO"])[0]
    assert r["products"] == ["COC(OC)c1ccccc1"]


def test_apply_reaction_variants():
    # boronate ester + triflate, trifluoroborate, charge-separated ylide, both azide resonance forms
    assert reactions.apply_reaction("suzuki", ["CC1(C)OB(OC1(C)C)c1ccccc1", "FC(F)(F)S(=O)(=O)Oc1ccncc1"]
                                    )[0]["products"] == [canon("c1ccc(-c2ccncc2)cc1")]
    assert reactions.apply_reaction("suzuki", ["F[B-](F)(F)c1ccccc1", "Brc1ccccc1"])[0]["products"] == [
        "c1ccc(-c2ccccc2)cc1"]
    assert reactions.apply_reaction("wittig", ["[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1", "CC(C)=O"]
                                    )[0]["products"] == ["C=C(C)C"]
    for az in ("[N-]=[N+]=NCc1ccccc1", "[N-]([N+]#N)Cc1ccccc1"):
        assert reactions.apply_reaction("cuaac", [az, "C#Cc1ccccc1"])[0]["products"] == [
            "c1ccc(Cn2cc(-c3ccccc3)nn2)cc1"]
    # Baeyer–Villiger migrates the non-methyl group; hydroboration is anti-Markovnikov
    assert reactions.apply_reaction("baeyer villiger", ["CC(=O)c1ccccc1"])[0]["products"] == ["CC(=O)Oc1ccccc1"]
    assert reactions.apply_reaction("hydroboration", ["CC1=CCCCC1"])[0]["products"] == ["CC1CCCCC1O"]
    # secondary amines / ammonia
    assert reactions.apply_reaction("boc", ["C1CCNCC1"])[0]["products"] == ["CC(C)(C)OC(=O)N1CCCCC1"]
    assert reactions.apply_reaction("amide coupling", ["CC(=O)O", "N"])[0]["products"] == ["CC(N)=O"]
    assert reactions.apply_reaction("michael", ["C[N+](=O)[O-]", "C=CC(C)=O"])[0]["products"] == [
        canon("CC(=O)CCC[N+](=O)[O-]")]


def test_apply_reaction_error_messages():
    with pytest.raises(ReactionError, match=r"takes 2 reactant\(s\).*aryl/vinyl halide.*got 1"):
        reactions.apply_reaction("suzuki", ["Brc1ccccc1"])
    with pytest.raises(ReactionError) as exc:
        reactions.apply_reaction("suzuki", ["CCO", "CCO"])
    msg = str(exc.value)
    assert "do not fit the template roles" in msg
    assert "role 1 'aryl/vinyl halide'" in msg and "matched by none of the reactants" in msg
    assert "reactant 1 'CCO' matches no role" in msg
    with pytest.raises(ReactionError) as exc:
        reactions.apply_reaction("fischer", ["CCO", "CCO"])
    msg = str(exc.value)
    assert "role 1 'carboxylic acid'" in msg and "role 2 'alcohol'" in msg
    assert "reactant 1 'CCO', reactant 2 'CCO'" in msg      # both fit the alcohol role, none the acid
    with pytest.raises(ReactionError, match="Reactant 2 .*not a molecule"):
        reactions.apply_reaction("suzuki", ["Brc1ccccc1", "not a molecule at all"])
    with pytest.raises(ReactionError, match="Unknown reaction"):
        reactions.apply_reaction("no such reaction", ["CCO"])
    with pytest.raises(ReactionError, match="list of"):
        reactions.apply_reaction("suzuki", "Brc1ccccc1.OB(O)c1ccccc1")  # type: ignore[arg-type]
    # substrates that do not carry the functional group are refused, not silently returned
    for name, rs in [("wolff kishner", ["CC(=O)O"]), ("alkene hydrogenation", ["c1ccccc1"]),
                     ("tbs", ["CC(=O)O"]), ("benzyl deprotection", ["OCc1ccccc1"])]:
        with pytest.raises(ReactionError):
            reactions.apply_reaction(name, rs)


# --------------------------------------------------------------------------- #
# suggest_reactions
# --------------------------------------------------------------------------- #
def test_suggest_reactions():
    hits = reactions.suggest_reactions(["c1ccccc1Br", "OB(O)c1ccccc1"])
    keys = [h["template"] for h in hits]
    assert "suzuki_coupling" in keys
    suzuki = next(h for h in hits if h["template"] == "suzuki_coupling")
    assert suzuki["products"] == ["c1ccc(-c2ccccc2)cc1"]
    assert suzuki["roles"]["aryl/vinyl halide"] == "Brc1ccccc1"
    assert suzuki["n_outcomes"] == 1 and ">>" in suzuki["mapped"]
    assert "Pd(PPh3)4 (5 mol%)" in suzuki["reagents"]
    assert "heck_reaction" not in keys and "fischer_esterification" not in keys
    keys = [h["template"] for h in reactions.suggest_reactions(["O=Cc1ccccc1", "CN"])]
    assert "reductive_amination" in keys and "imine_formation" in keys and "suzuki_coupling" not in keys
    keys = [h["template"] for h in reactions.suggest_reactions(["OCc1ccccc1"])]
    assert {"swern_oxidation", "jones_oxidation", "tbs_protection", "appel_reaction"} <= set(keys)
    assert "benzyl_deprotection" not in keys and "suzuki_coupling" not in keys
    assert reactions.suggest_reactions(["[Na+].[Cl-]", "[K+].[Br-]", "O", "O"]) == []
    json.dumps(reactions.suggest_reactions(["CC(=O)O", "CN"]))


# --------------------------------------------------------------------------- #
# to_step / describe
# --------------------------------------------------------------------------- #
def test_to_step_builds_layout_step():
    from chemdrawmcp import layout

    step = reactions.to_step("suzuki", ["Brc1ccccc1", "OB(O)c1ccccc1"])
    assert isinstance(step, layout.Step)
    assert len(step.reactants) == 2 and len(step.products) == 1
    assert step.above == ["Pd(PPh3)4 (5 mol%)", "K2CO3"] and step.below == ["dioxane/H2O, 80 °C"]
    assert step.arrow == "full" and step.atom_map is True
    for comp in [*step.reactants, *step.products]:
        assert isinstance(comp, layout.Component)
        assert comp.mol.GetNumConformers() == 1                       # depiction-ready
        assert not any(b.GetIsAromatic() for b in comp.mol.GetBonds())  # kekulised
    prod_maps = {a.GetAtomMapNum() for a in step.products[0].mol.GetAtoms()} - {0}
    react_maps = {a.GetAtomMapNum() for c in step.reactants for a in c.mol.GetAtoms()} - {0}
    assert prod_maps and prod_maps <= react_maps
    assert Chem.MolToSmiles(chem._aromatic_copy(step.products[0].mol)).count(":") == 12
    # it renders through the CDXML layout engine (atom map pairs included)
    doc = layout.build_reaction(step, name="suzuki")
    xml = doc.to_xml() if hasattr(doc, "to_xml") else str(doc)
    assert "<step" in xml and "ReactionStepAtomMap" in xml

    bare = reactions.to_step("grignard ester", ["CCOC(=O)c1ccccc1", "C[Mg]Br"], include_reagents=False,
                             atom_map=False, labels=True)
    assert bare.above == [] and bare.below == [] and bare.atom_map is False
    assert len(bare.reactants) == 2                                   # one component per role, not per equivalent
    assert [c.label for c in bare.reactants] == ["ester", "Grignard reagent (2 equiv)"]
    assert all(a.GetAtomMapNum() == 0 for c in [*bare.reactants, *bare.products] for a in c.mol.GetAtoms())

    second = reactions.to_step("nitration", ["Cc1ccccc1"], product_index=2)
    assert len(second.products) == 1
    with pytest.raises(ReactionError, match="out of range"):
        reactions.to_step("nitration", ["Cc1ccccc1"], product_index=3)
    with pytest.raises(ReactionError):
        reactions.to_step("suzuki", ["CCO", "CCO"])


def test_describe():
    text = reactions.describe("Suzuki")
    assert "Suzuki–Miyaura coupling" in text and "suzuki_coupling" in text
    assert "class: C–C coupling" in text
    assert "1. aryl/vinyl halide" in text and "2. boronic acid/ester" in text
    assert "Pd(PPh3)4" in text and "dioxane/H2O, 80 °C" in text
    assert "SMARTS: " + REACTIONS["suzuki_coupling"].smarts in text
    assert "Brc1ccccc1 + OB(O)c1ccccc1 -> c1ccc(-c2ccccc2)cc1" in text
    assert "(2 equiv)" in reactions.describe("grignard_ester")
    with pytest.raises(ReactionError):
        reactions.describe("nope nope nope")


def test_mapped_smiles_round_trips_through_chem():
    r = reactions.apply_reaction("cuaac", ["[N-]=[N+]=NCc1ccccc1", "C#Cc1ccccc1"])[0]
    rxn = chem.parse_reaction(r["mapped"])
    assert len(rxn.reactants) == 2 and len(rxn.products) == 1
    maps = sorted(a.GetAtomMapNum() for a in rxn.products[0].GetAtoms() if a.GetAtomMapNum())
    assert maps == list(range(1, len(maps) + 1))
    props = chem.reaction_properties(rxn)
    assert props["balanced"] is True                       # CuAAC is atom-economic
    assert chem.reaction_smiles(rxn).count(">") == 2
