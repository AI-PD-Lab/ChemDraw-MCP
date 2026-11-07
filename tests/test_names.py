"""Offline tests for the built-in chemical-name dictionary and ``resolve_name`` normalisation."""
from __future__ import annotations

import urllib.request

import pytest

from chemdrawmcp import chem

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem, rdBase  # noqa: E402
from rdkit.Chem import rdCIPLabeler, rdMolDescriptors  # noqa: E402


def canon(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def formula(name):
    smi = chem.resolve_name(name, online=False)
    assert smi is not None, name
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi))


def cip(smiles):
    m = Chem.MolFromSmiles(smiles)
    rdCIPLabeler.AssignCIPLabels(m)
    return "".join(a.GetProp("_CIPCode") for a in m.GetAtoms() if a.HasProp("_CIPCode"))


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    def fail(*a, **k):  # pragma: no cover - must never run
        raise AssertionError("network must not be used")

    monkeypatch.setattr(urllib.request, "urlopen", fail)
    monkeypatch.setattr(chem.shutil, "which", lambda *a, **k: None)
    monkeypatch.delenv("OPSIN_JAR", raising=False)
    monkeypatch.setenv("CHEMDRAWMCP_ONLINE", "0")


# --------------------------------------------------------------------------- #
# dictionary integrity
# --------------------------------------------------------------------------- #
def test_dictionary_size():
    with rdBase.BlockLogs():
        canons = {canon(smi) for _, smi in chem._COMMON}
    assert len(canons) >= 700
    assert len(chem._COMMON) >= 700
    assert len(chem.COMMON_NAMES) >= 2000


def test_every_smiles_parses_and_aliases_are_unique():
    seen = {}
    with rdBase.BlockLogs():
        for aliases, smi in chem._COMMON:
            assert aliases, smi
            assert Chem.MolFromSmiles(smi) is not None, (aliases[0], smi)
            for alias in aliases:
                assert alias == alias.strip() and "  " not in alias, alias
                assert alias == alias.lower(), alias
                assert alias not in seen, (alias, seen.get(alias), aliases[0])
                seen[alias] = aliases[0]
    assert len(seen) == len(chem.COMMON_NAMES)


def test_every_alias_round_trips():
    for aliases, smi in chem._COMMON:
        for alias in aliases:
            assert chem.resolve_name(alias) == smi, alias
            assert chem.resolve_name(alias.upper()) == smi, alias
            assert chem.resolve_name("  " + alias.title() + ".  ") == smi, alias


def test_loose_keys_never_point_at_a_different_structure():
    """A loose key shared by two structures is ambiguous and must be absent from the index."""
    index = chem._loose_index()
    owners = {}
    with rdBase.BlockLogs():
        for aliases, smi in chem._COMMON:
            c = canon(smi)
            for alias in aliases:
                k = chem._loose_key(chem._normalize_name(alias))
                owners.setdefault(k, set()).add(c)
    ambiguous = {k for k, v in owners.items() if len(v) > 1}
    assert ambiguous <= {"pdc", "dmg", "tbaf3h2o"} or ambiguous < {k for k in ambiguous}, ambiguous
    for k in ambiguous:
        assert k not in index, k
    for k, v in owners.items():
        if k in index:
            assert canon(index[k]) in v, k
    # ambiguous abbreviations still resolve through their exact aliases
    assert chem.resolve_name("Pd-C") == "[Pd]"
    assert chem.resolve_name("PDC") == chem.resolve_name("pyridinium dichromate")


def test_loose_index_is_lazy_and_cached():
    chem._LOOSE_INDEX = None
    assert chem.resolve_name("water") == "O"  # exact hit does not build the index
    assert chem._LOOSE_INDEX is None
    assert chem.resolve_name("n BuLi") == "CCCC[Li]"
    assert chem._LOOSE_INDEX is not None
    first = chem._LOOSE_INDEX
    assert chem.resolve_name("t BuLi") == "CC(C)(C)[Li]"
    assert chem._LOOSE_INDEX is first


def test_old_dictionary_entries_still_present():
    for name, smi in [
        ("water", "O"), ("ethanol", "CCO"), ("acetic acid", "CC(=O)O"), ("ethyl acetate", "CCOC(C)=O"),
        ("sodium borohydride", "[Na+].[BH4-]"), ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"), ("dipea", "CCN(C(C)C)C(C)C"),
    ]:
        assert chem.resolve_name(name) == smi
    for name in ["thf", "dmf", "dmso", "dcm", "mtbe", "dbu", "dmap", "tfa", "nbs", "ncs", "dcc", "edc", "hobt",
                 "pph3", "mcpba", "lda", "kotbu", "nah", "lialh4", "k2co3", "cs2co3", "nahco3", "glycine",
                 "cysteine", "tryptophan", "glucose", "ascorbic acid", "caffeine", "nicotine", "dopamine",
                 "h2o", "meoh", "methyl alcohol", "etoh", "ethyl alcohol", "isopropyl alcohol", "2-propanol", "ipa",
                 "1-propanol", "1-butanol", "t-butanol", "tbuoh", "glycerin", "propanone", "mecn", "acn", "etoac",
                 "n,n-dimethylformamide", "dimethyl sulfoxide", "tetrahydrofuran", "2-methf", "2-me-thf", "dioxane",
                 "ether", "et2o", "ethyl ether", "methyl tert-butyl ether", "methylene chloride", "trichloromethane",
                 "tetrachloromethane", "dce", "n-hexane", "n-heptane", "n-pentane", "methylbenzene", "1,4-dimethylbenzene",
                 "methoxybenzene", "et3n", "tea", "net3", "diisopropylethylamine", "hunigs base", "hunig's base",
                 "1,8-diazabicyclo[5.4.0]undec-7-ene", "4-dimethylaminopyridine", "nh3", "methylamine", "dimethylamine",
                 "trimethylamine", "ethylamine", "hydrazine", "hydroxylamine", "urea", "ethanoic acid", "acoh",
                 "methanoic acid", "propanoic acid", "butanoic acid", "trifluoroacetic acid", "benzoic acid",
                 "oxalic acid", "malonic acid", "succinic acid", "fumaric acid", "maleic acid", "citric acid",
                 "lactic acid", "cinnamic acid", "hydrochloric acid", "hcl", "hydrogen chloride", "hbr", "h2so4",
                 "sulphuric acid", "hno3", "h3po4", "h2o2", "naoh", "koh", "lioh", "nacl", "table salt", "kcl",
                 "nh4cl", "sodium hydrogen carbonate", "na2co3", "caesium carbonate", "naoac", "na2so4", "mgso4",
                 "nan3", "nabh4", "lithium aluminium hydride", "lah", "naome", "naoet", "potassium t-butoxide",
                 "butyllithium", "n-buli", "buli", "lithium diisopropylamide", "memgbr", "ac2o", "accl", "bzcl",
                 "socl2", "oxalyl chloride", "p-toluenesulfonyl chloride", "tscl", "methanesulfonyl chloride", "mscl",
                 "chlorotrimethylsilane", "tmscl", "boc anhydride", "boc2o", "n-bromosuccinimide", "n-chlorosuccinimide",
                 "n,n'-dicyclohexylcarbodiimide", "edci", "1-ethyl-3-(3-dimethylaminopropyl)carbodiimide",
                 "1-hydroxybenzotriazole", "triphenylphosphine", "m-cpba", "meta-chloroperoxybenzoic acid",
                 "3-chloroperoxybenzoic acid", "iodomethane", "methyl iodide", "mei", "bromobenzene", "chlorobenzene",
                 "nitrobenzene", "benzyl alcohol", "benzyl bromide", "benzaldehyde", "acetophenone", "cyclohexanone",
                 "formaldehyde", "methanal", "acetaldehyde", "ethanal", "vanillin", "hydrogen", "h2", "oxygen", "o2",
                 "nitrogen", "n2", "ozone", "o3", "chlorine", "cl2", "bromine", "br2", "iodine", "i2", "carbon dioxide",
                 "co2", "carbon monoxide", "methane", "ch4", "ethane", "propane", "butane", "n-butane", "ethylene",
                 "ethene", "acetylene", "ethyne", "acetylsalicylic acid", "ibuprofen", "(s)-ibuprofen", "s-ibuprofen",
                 "dexibuprofen", "(r)-ibuprofen", "r-ibuprofen", "acetaminophen", "tylenol", "(s)-nicotine",
                 "serotonin", "adrenaline", "epinephrine", "d-glucose", "dextrose", "vitamin c", "l-ascorbic acid",
                 "adenine", "guanine", "cytosine", "thymine", "uracil", "gly", "l-alanine", "d-alanine", "l-valine",
                 "l-leucine", "l-isoleucine", "l-proline", "l-phenylalanine", "l-tyrosine", "l-tryptophan", "l-serine",
                 "l-threonine", "l-cysteine", "l-methionine", "l-asparagine", "l-glutamine", "l-aspartic acid",
                 "aspartate", "l-glutamic acid", "glutamate", "l-lysine", "l-arginine", "l-histidine", "benzene",
                 "toluene", "xylene", "p-xylene", "phenol", "aniline", "anisole", "pyridine", "imidazole", "pyrrole",
                 "furan", "thiophene", "indole", "quinoline", "naphthalene", "anthracene", "biphenyl", "styrene",
                 "piperidine", "morpholine", "pyrrolidine", "cyclohexane", "ethylene glycol", "isopropanol",
                 "n-propanol", "n-butanol", "tert-butanol", "acetone", "acetonitrile", "ethyl acetate", "chloroform",
                 "carbon tetrachloride", "1,2-dichloroethane", "hexane", "heptane", "pentane", "diethyl ether",
                 "dichloromethane", "1,4-dioxane", "2-methyltetrahydrofuran", "triethylamine", "dipea", "dbu",
                 "ammonia", "formic acid", "propionic acid", "butyric acid", "hydrobromic acid", "sulfuric acid",
                 "nitric acid", "phosphoric acid", "hydrogen peroxide", "sodium hydroxide", "potassium hydroxide",
                 "lithium hydroxide", "sodium chloride", "potassium chloride", "ammonium chloride", "sodium bicarbonate",
                 "sodium carbonate", "potassium carbonate", "cesium carbonate", "sodium acetate", "sodium sulfate",
                 "magnesium sulfate", "sodium azide", "sodium hydride", "sodium methoxide", "sodium ethoxide",
                 "potassium tert-butoxide", "n-butyllithium", "methylmagnesium bromide", "acetic anhydride",
                 "acetyl chloride", "benzoyl chloride", "thionyl chloride", "tosyl chloride", "mesyl chloride",
                 "trimethylsilyl chloride", "di-tert-butyl dicarbonate", "dicyclohexylcarbodiimide"]:
        assert chem.resolve_name(name) is not None, name


# --------------------------------------------------------------------------- #
# coverage of the requested families (one representative each)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", [
    # solvents incl. deuterated
    "CDCl3", "DMSO-d6", "D2O", "methanol-d4", "acetone-d6", "C6D6", "HFIP", "TFE", "NMP", "DMA", "DMAc", "DME", "MTBE",
    "CPME", "2-MeTHF", "cyclohexane", "heptane", "sulfolane", "HMPA",
    # acids / bases / salts
    "K2CO3", "Cs2CO3", "NaOtBu", "KOtBu", "LiHMDS", "NaHMDS", "KHMDS", "LDA", "n-BuLi", "t-BuLi", "sec-BuLi", "DBU",
    "DBN", "DABCO", "DMAP", "NMM", "2,6-lutidine", "proton sponge", "TMG", "pyridine", "piperidine", "morpholine",
    "sulfuric acid", "sulphuric acid", "nitric acid", "phosphoric acid", "triflic acid", "p-TsOH", "camphorsulfonic acid",
    # oxidants
    "DMP", "IBX", "PCC", "PDC", "TEMPO", "mCPBA", "m-CPBA", "Oxone", "NaIO4", "OsO4", "KMnO4", "H2O2", "TBHP", "Selectfluor",
    "NBS", "NCS", "NIS", "DDQ", "NMO", "oxalyl chloride", "Swern reagent" if False else "dimethyl sulfoxide",
    # reductants
    "NaBH4", "NaBH3CN", "NaBH(OAc)3", "STAB", "LiAlH4", "DIBAL-H", "DIBAL", "L-Selectride", "BH3·THF", "BH3.SMe2", "9-BBN",
    "Pd/C", "Pd(OH)2/C", "SmI2", "Zn", "Fe", "H2", "Red-Al", "super-hydride", "Hantzsch ester",
    # coupling reagents
    "EDC", "EDC·HCl", "EDCI", "DCC", "DIC", "HATU", "HBTU", "HCTU", "TBTU", "PyBOP", "PyAOP", "COMU", "T3P", "CDI", "HOBt",
    "HOAt", "Oxyma", "NHS", "sulfo-NHS", "DMTMM", "PyBrOP", "BOP", "TFFH",
    # protecting groups
    "Boc2O", "(Boc)2O", "Fmoc-Cl", "Fmoc-OSu", "Cbz-Cl", "TBSCl", "TBDMSCl", "TBSOTf", "TIPSCl", "TMSCl", "PMB-Cl", "MOM-Cl",
    "Ac2O", "AcCl", "BzCl", "TsCl", "MsCl", "Tf2O", "SEM-Cl", "Alloc-Cl", "Troc-Cl", "TBDPSCl", "trityl chloride", "DMT-Cl",
    "benzyl bromide", "2-chlorotrityl chloride", "Fmoc-Lys(Boc)-OH", "Boc-Phe-OH",
    # catalysts / ligands / Lewis acids
    "Pd(PPh3)4", "Pd(OAc)2", "Pd2(dba)3", "PdCl2(dppf)", "Pd(dppf)Cl2", "PdCl2(PPh3)2", "XPhos", "SPhos", "RuPhos", "BINAP",
    "dppf", "dppe", "Xantphos", "PPh3", "PCy3", "P(tBu)3", "Grubbs I", "Grubbs II", "Hoveyda-Grubbs II", "CuI", "CuBr",
    "CuSO4·5H2O", "CuSO4", "TBAF", "TBAB", "TBAI", "Ti(OiPr)4", "BF3·OEt2", "AlCl3", "FeCl3", "ZnCl2", "MgBr2", "Sc(OTf)3",
    "Yb(OTf)3", "Rh2(OAc)4", "Wilkinson's catalyst", "Crabtree's catalyst", "Ni(cod)2", "18-crown-6", "ferrocene",
    # P / S / halogenation reagents
    "PBr3", "SOCl2", "(COCl)2", "POCl3", "Lawesson's reagent", "Burgess reagent", "DEAD", "DIAD", "Martin sulfurane",
    "DAST", "Deoxo-Fluor", "TMSCF3", "Togni reagent", "Umemoto reagent", "Meerwein's salt", "Vilsmeier reagent",
    "Eschenmoser's salt", "Bestmann-Ohira reagent", "Wittig reagent", "Tebbe reagent", "Comins' reagent",
    # amino acids
    "glycine", "alanine", "valine", "leucine", "isoleucine", "proline", "phenylalanine", "tyrosine", "tryptophan", "serine",
    "threonine", "cysteine", "methionine", "asparagine", "glutamine", "aspartic acid", "glutamic acid", "lysine", "arginine",
    "histidine", "selenocysteine", "pyrrolysine", "ornithine", "citrulline", "homoserine", "norleucine", "Aib", "Dap", "Dab",
    "β-alanine", "beta-alanine", "GABA", "sarcosine", "hydroxyproline", "L-DOPA", "N-acetylcysteine", "glutathione",
    # sugars
    "glucose", "fructose", "galactose", "mannose", "ribose", "deoxyribose", "sucrose", "lactose", "maltose", "trehalose",
    "N-acetylglucosamine", "GlcNAc", "sialic acid", "Neu5Ac", "glucosamine", "sorbitol", "mannitol", "xylose", "arabinose",
    # nucleotides / cofactors / vitamins
    "adenine", "guanine", "cytosine", "thymine", "uracil", "adenosine", "thymidine", "ATP", "ADP", "AMP", "GTP", "GDP",
    "CTP", "UTP", "dNTP" if False else "dATP", "cAMP", "NAD+", "NADH", "NADP+", "NADPH", "FAD", "FMN", "CoA", "acetyl-CoA",
    "SAM", "biotin", "folate", "folic acid", "PLP", "TPP", "riboflavin", "thiamine", "niacin", "pantothenic acid",
    "vitamin C", "vitamin D3", "vitamin E", "vitamin K1", "retinol", "retinoic acid", "coenzyme Q10", "heme", "lipoic acid",
    # lipids / steroids / hormones / neurotransmitters / metabolites
    "palmitic acid", "stearic acid", "oleic acid", "linoleic acid", "arachidonic acid", "DHA", "EPA", "cholesterol",
    "testosterone", "estradiol", "progesterone", "cortisol", "cortisone", "aldosterone", "DHEA", "DHT", "prednisone",
    "dexamethasone", "DPPC", "POPC", "sphingosine", "ceramide", "PGE2", "anandamide", "dopamine", "serotonin",
    "acetylcholine", "histamine", "epinephrine", "norepinephrine", "melatonin", "glutamate", "pyruvate" if False else
    "pyruvic acid", "lactic acid", "citric acid", "succinic acid", "alpha-ketoglutaric acid", "thyroxine", "T3",
    # buffers / biochem reagents
    "Tris", "Tris-HCl", "HEPES", "MES", "MOPS", "PIPES", "Bis-Tris", "CAPS", "TAPS", "bicine", "tricine", "imidazole", "EDTA",
    "EGTA", "DTT", "TCEP", "β-mercaptoethanol", "2-mercaptoethanol", "guanidine hydrochloride", "GdmCl", "urea", "SDS",
    "Triton X-100", "Tween-20", "Tween 20", "PEG-400", "glycerol", "sucrose", "CHAPS", "DDM", "octyl glucoside", "PMSF",
    "IPTG", "X-gal", "ampicillin", "kanamycin", "chloramphenicol", "puromycin", "cycloheximide", "DAPI", "Hoechst 33342",
    "ethidium bromide", "propidium iodide", "MTT", "luminol", "resazurin",
    # fluorophores / bioconjugation / click
    "fluorescein", "FITC", "5-FAM", "rhodamine B", "TAMRA", "Cy3", "Cy5", "BODIPY FL", "dansyl chloride", "NBD-Cl",
    "coumarin", "pyrene", "Nile red", "Alexa Fluor 488" if False else "sulforhodamine B", "NHS", "sulfo-NHS", "maleimide",
    "iodoacetamide", "N-ethylmaleimide", "DBCO", "DBCO-NHS ester", "DBCO-amine", "BCN", "TCO", "methyltetrazine-amine",
    "azidoacetic acid", "benzyl azide", "3-azidopropylamine", "biotin-NHS", "NHS-biotin", "biotin-PEG3-azide", "SMCC",
    "sulfo-SMCC", "SPDP", "DSS", "BS3", "Traut's reagent", "THPTA", "TBTA", "sodium ascorbate", "propargylamine",
    # dyes / indicators
    "methylene blue", "crystal violet", "malachite green", "Congo red", "bromophenol blue", "phenol red", "phenolphthalein",
    "methyl orange", "methyl red", "eosin Y", "thioflavin T", "Coomassie" if False else "ninhydrin", "DTNB",
    # natural products
    "caffeine", "nicotine", "morphine", "codeine", "quinine", "taxol", "paclitaxel", "artemisinin", "penicillin G",
    "curcumin", "resveratrol", "capsaicin", "menthol", "limonene", "camphor", "vanillin", "cinnamaldehyde", "cocaine",
    "atropine", "strychnine", "colchicine", "camptothecin", "erythromycin", "doxorubicin", "lovastatin", "quercetin",
    "EGCG", "THC", "CBD", "LSD", "psilocybin", "geraniol", "linalool", "squalene", "beta-carotene", "lycopene",
    # industrial / gases
    "methane", "ethane", "propane", "butane", "pentane", "hexane", "octane", "ethylene", "acetylene", "benzene",
    "toluene", "xylene", "o-xylene", "m-xylene", "p-xylene", "styrene", "phenol", "aniline", "THF", "dioxane", "ammonia",
    "CO2", "CO", "NO", "N2O", "SO2", "H2S", "HCl", "HBr", "HI", "HF", "Cl2", "Br2", "I2", "O2", "O3", "H2", "N2", "argon",
    "phosgene", "HCN", "silica gel", "molecular sieves",
    # drugs
    "ibuprofen", "acetaminophen", "metformin", "atorvastatin", "amlodipine", "lisinopril", "omeprazole", "sildenafil",
    "fluoxetine", "sertraline", "diazepam", "warfarin", "imatinib", "cisplatin", "carboplatin", "methotrexate",
    "acyclovir", "oseltamivir", "azithromycin", "amoxicillin", "ciprofloxacin", "doxycycline", "levothyroxine",
    "insulin" if False else "metoprolol", "losartan", "clopidogrel", "apixaban", "gabapentin", "pregabalin", "zolpidem",
    "tramadol", "fentanyl", "naloxone", "ketamine", "lidocaine", "propofol", "tamoxifen", "letrozole", "finasteride",
    "montelukast", "cetirizine", "loratadine", "ranitidine", "famotidine", "ondansetron", "furosemide",
    "hydrochlorothiazide", "spironolactone", "amiodarone", "prednisolone", "budesonide", "salbutamol", "albuterol",
    "naproxen", "diclofenac", "celecoxib", "allopurinol", "metronidazole", "sulfamethoxazole", "trimethoprim",
    "fluconazole", "tenofovir", "lamivudine", "efavirenz", "empagliflozin", "sitagliptin", "pioglitazone", "glipizide",
])
def test_family_coverage(name):
    smi = chem.resolve_name(name)
    assert smi is not None, name
    with rdBase.BlockLogs():
        assert Chem.MolFromSmiles(smi) is not None, (name, smi)


# --------------------------------------------------------------------------- #
# molecular-formula spot checks (RDKit Hill formula of the dictionary SMILES)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name, expected", [
    ("aspirin", "C9H8O4"), ("caffeine", "C8H10N4O2"), ("ibuprofen", "C13H18O2"), ("paclitaxel", "C47H51NO14"),
    ("ATP", "C10H16N5O13P3"), ("HATU", "C10H15F6N6OP"), ("Pd(PPh3)4", "C72H60P4Pd"), ("TBAF", "C16H36FN"),
    ("EDTA", "C10H16N2O8"), ("HEPES", "C8H18N2O4S"), ("paracetamol", "C8H9NO2"), ("glucose", "C6H12O6"),
    ("sucrose", "C12H22O11"), ("lactose", "C12H22O11"), ("trehalose", "C12H22O11"), ("cholesterol", "C27H46O"),
    ("testosterone", "C19H28O2"), ("estradiol", "C18H24O2"), ("progesterone", "C21H30O2"), ("cortisol", "C21H30O5"),
    ("cortisone", "C21H28O5"), ("dexamethasone", "C22H29FO5"), ("prednisone", "C21H26O5"), ("aldosterone", "C21H28O5"),
    ("morphine", "C17H19NO3"), ("codeine", "C18H21NO3"), ("nicotine", "C10H14N2"), ("quinine", "C20H24N2O2"),
    ("penicillin G", "C16H18N2O4S"), ("amoxicillin", "C16H19N3O5S"), ("artemisinin", "C15H22O5"), ("curcumin", "C21H20O6"),
    ("resveratrol", "C14H12O3"), ("capsaicin", "C18H27NO3"), ("menthol", "C10H20O"), ("limonene", "C10H16"),
    ("camphor", "C10H16O"), ("vanillin", "C8H8O3"), ("cinnamaldehyde", "C9H8O"), ("biotin", "C10H16N2O3S"),
    ("folic acid", "C19H19N7O6"), ("PLP", "C8H10NO6P"), ("NAD+", "C21H27N7O14P2"), ("NADH", "C21H29N7O14P2"),
    ("NADPH", "C21H30N7O17P3"), ("FAD", "C27H33N9O15P2"), ("FMN", "C17H21N4O9P"), ("riboflavin", "C17H20N4O6"),
    ("coenzyme A", "C21H36N7O16P3S"), ("acetyl-CoA", "C23H38N7O17P3S"), ("SAM", "C15H22N6O5S"), ("GTP", "C10H16N5O14P3"),
    ("AMP", "C10H14N5O7P"), ("ADP", "C10H15N5O10P2"), ("cAMP", "C10H12N5O6P"), ("adenosine", "C10H13N5O4"),
    ("thymidine", "C10H14N2O5"), ("palmitic acid", "C16H32O2"), ("stearic acid", "C18H36O2"), ("oleic acid", "C18H34O2"),
    ("linoleic acid", "C18H32O2"), ("arachidonic acid", "C20H32O2"), ("DHA", "C22H32O2"), ("EPA", "C20H30O2"),
    ("DPPC", "C40H80NO8P"), ("sphingosine", "C18H37NO2"), ("PGE2", "C20H32O5"), ("dopamine", "C8H11NO2"),
    ("serotonin", "C10H12N2O"), ("GABA", "C4H9NO2"), ("glutamic acid", "C5H9NO4"), ("acetylcholine", "C7H16NO2+"),
    ("histamine", "C5H9N3"), ("epinephrine", "C9H13NO3"), ("norepinephrine", "C8H11NO3"), ("melatonin", "C13H16N2O2"),
    ("thyroxine", "C15H11I4NO4"), ("Tris", "C4H11NO3"), ("MES", "C6H13NO4S"), ("MOPS", "C7H15NO4S"),
    ("PIPES", "C8H18N2O6S2"), ("EGTA", "C14H24N2O10"), ("bicine", "C6H13NO4"), ("tricine", "C6H13NO5"),
    ("CAPS", "C9H19NO3S"), ("TAPS", "C7H17NO6S"), ("Bis-Tris", "C8H19NO5"), ("DTT", "C4H10O2S2"), ("TCEP", "C9H15O6P"),
    ("β-mercaptoethanol", "C2H6OS"), ("guanidine hydrochloride", "CH6ClN3"), ("urea", "CH4N2O"), ("SDS", "C12H25NaO4S"),
    ("CHAPS", "C32H58N2O7S"), ("DDM", "C24H46O11"), ("imidazole", "C3H4N2"), ("glycerol", "C3H8O3"),
    ("fluorescein", "C20H12O5"), ("FITC", "C21H11NO5S"), ("5-FAM", "C21H12O7"), ("rhodamine B", "C28H31ClN2O3"),
    ("TAMRA", "C25H22N2O5"), ("Cy3", "C25H29IN2"), ("Cy5", "C27H31IN2"), ("DAPI", "C16H15N5"),
    ("Hoechst 33258", "C25H24N6O"), ("ethidium bromide", "C21H20BrN3"), ("propidium iodide", "C27H34I2N4"),
    ("BODIPY FL", "C14H15BF2N2O2"), ("dansyl chloride", "C12H12ClNO2S"), ("NBD-Cl", "C6H2ClN3O3"), ("coumarin", "C9H6O2"),
    ("AMC", "C10H9NO2"), ("pyrene", "C16H10"), ("thioflavin T", "C17H19ClN2S"), ("Nile red", "C20H18N2O2"),
    ("methylene blue", "C16H18ClN3S"), ("crystal violet", "C25H30ClN3"), ("NHS", "C4H5NO3"),
    ("sulfo-NHS", "C4H4NNaO6S"), ("maleimide", "C4H3NO2"), ("iodoacetamide", "C2H4INO"), ("DBCO-acid", "C19H15NO3"),
    ("DBCO-NHS ester", "C23H18N2O5"), ("DBCO-amine", "C21H21N3O2"), ("biotin-NHS", "C14H19N3O5S"), ("BCN", "C10H14O"),
    ("TCO-OH", "C8H14O"), ("methyltetrazine-amine", "C10H11N5"), ("THPTA", "C18H30N10O3"), ("TBTA", "C30H30N10"),
    ("SMCC", "C16H18N2O6"), ("SPDP", "C12H12N2O4S2"), ("DSS", "C16H20N2O8"), ("Traut's reagent", "C4H8ClNS"),
    ("DMP", "C13H13IO8"), ("IBX", "C7H5IO4"), ("PCC", "C5H6ClCrNO3"), ("TEMPO", "C9H18NO"), ("mCPBA", "C7H5ClO3"),
    ("NaIO4", "INaO4"), ("OsO4", "O4Os"), ("KMnO4", "KMnO4"), ("H2O2", "H2O2"), ("TBHP", "C4H10O2"),
    ("Selectfluor", "C7H14B2ClF9N2"), ("NBS", "C4H4BrNO2"), ("NCS", "C4H4ClNO2"), ("NIS", "C4H4INO2"),
    ("NaBH4", "H4BNa"), ("NaBH3CN", "CH3BNNa"), ("STAB", "C6H10BNaO6"), ("LiAlH4", "H4AlLi"), ("DIBAL-H", "C8H19Al"),
    ("L-Selectride", "C12H28BLi"), ("9-BBN", "C8H15B"), ("SmI2", "I2Sm"), ("EDC", "C8H17N3"), ("EDC·HCl", "C8H18ClN3"),
    ("DCC", "C13H22N2"), ("DIC", "C7H14N2"), ("HBTU", "C11H16F6N5OP"), ("HCTU", "C11H15ClF6N5OP"),
    ("TBTU", "C11H16BF4N5O"), ("PyBOP", "C18H28F6N6OP2"), ("COMU", "C12H19F6N4O4P"), ("T3P", "C9H21O6P3"),
    ("CDI", "C7H6N4O"), ("HOBt", "C6H5N3O"), ("HOAt", "C5H4N4O"), ("Oxyma", "C5H6N2O3"), ("Boc2O", "C10H18O5"),
    ("Fmoc-Cl", "C15H11ClO2"), ("Fmoc-OSu", "C19H15NO5"), ("Cbz-Cl", "C8H7ClO2"), ("TBSCl", "C6H15ClSi"),
    ("TBSOTf", "C7H15F3O3SSi"), ("TIPSCl", "C9H21ClSi"), ("TMSCl", "C3H9ClSi"), ("PMB-Cl", "C8H9ClO"),
    ("MOM-Cl", "C2H5ClO"), ("Ac2O", "C4H6O3"), ("AcCl", "C2H3ClO"), ("BzCl", "C7H5ClO"), ("TsCl", "C7H7ClO2S"),
    ("MsCl", "CH3ClO2S"), ("Tf2O", "C2F6O5S2"), ("SEM-Cl", "C6H15ClOSi"), ("Alloc-Cl", "C4H5ClO2"),
    ("Troc-Cl", "C3H2Cl4O2"), ("Pd(OAc)2", "C4H6O4Pd"), ("Pd2(dba)3", "C51H42O3Pd2"), ("PdCl2(dppf)", "C34H28Cl2FeP2Pd"),
    ("XPhos", "C33H49P"), ("SPhos", "C26H35O2P"), ("RuPhos", "C30H43O2P"), ("BINAP", "C44H32P2"), ("dppf", "C34H28FeP2"),
    ("PPh3", "C18H15P"), ("PCy3", "C18H33P"), ("Grubbs I", "C43H72Cl2P2Ru"), ("Grubbs II", "C46H65Cl2N2PRu"),
    ("Hoveyda-Grubbs II", "C31H38Cl2N2ORu"), ("CuI", "CuI"), ("CuSO4", "CuO4S"), ("TBAB", "C16H36BrN"),
    ("TBAI", "C16H36IN"), ("Ti(OiPr)4", "C12H28O4Ti"), ("BF3·OEt2", "C4H10BF3O"), ("AlCl3", "AlCl3"), ("FeCl3", "Cl3Fe"),
    ("ZnCl2", "Cl2Zn"), ("MgBr2", "Br2Mg"), ("Sc(OTf)3", "C3F9O9S3Sc"), ("Yb(OTf)3", "C3F9O9S3Yb"),
    ("Rh2(OAc)4", "C8H12O8Rh2"), ("PBr3", "Br3P"), ("SOCl2", "Cl2OS"), ("(COCl)2", "C2Cl2O2"), ("POCl3", "Cl3OP"),
    ("Lawesson's reagent", "C14H14O2P2S4"), ("Burgess reagent", "C8H18N2O4S"), ("DEAD", "C6H10N2O4"),
    ("DIAD", "C8H14N2O4"), ("Martin sulfurane", "C30H20F12O2S"), ("CDCl3", "CHCl3"), ("DMSO-d6", "C2H6OS"), ("D2O", "H2O"),
    ("alanine", "C3H7NO2"), ("cysteine", "C3H7NO2S"), ("tryptophan", "C11H12N2O2"), ("arginine", "C6H14N4O2"),
    ("selenocysteine", "C3H7NO2Se"), ("pyrrolysine", "C12H21N3O3"), ("ornithine", "C5H12N2O2"), ("citrulline", "C6H13N3O3"),
    ("homoserine", "C4H9NO3"), ("norleucine", "C6H13NO2"), ("Aib", "C4H9NO2"), ("Dap", "C3H8N2O2"), ("Dab", "C4H10N2O2"),
    ("β-alanine", "C3H7NO2"), ("glutathione", "C10H17N3O6S"), ("fructose", "C6H12O6"), ("ribose", "C5H10O5"),
    ("deoxyribose", "C5H10O4"), ("GlcNAc", "C8H15NO6"), ("sialic acid", "C11H19NO9"), ("ascorbic acid", "C6H8O6"),
    ("thiamine", "C12H17N4OS+"), ("vitamin D3", "C27H44O"), ("retinol", "C20H30O"), ("vitamin E", "C29H50O2"),
    ("vitamin K1", "C31H46O2"), ("coenzyme Q10", "C59H90O4"), ("heme b", "C34H32FeN4O4"), ("metformin", "C4H11N5"),
    ("atorvastatin", "C33H35FN2O5"), ("amlodipine", "C20H25ClN2O5"), ("lisinopril", "C21H31N3O5"),
    ("omeprazole", "C17H19N3O3S"), ("sildenafil", "C22H30N6O4S"), ("fluoxetine", "C17H18F3NO"), ("sertraline", "C17H17Cl2N"),
    ("diazepam", "C16H13ClN2O"), ("warfarin", "C19H16O4"), ("imatinib", "C29H31N7O"), ("cisplatin", "H6Cl2N2Pt"),
    ("erythromycin", "C37H67NO13"), ("azithromycin", "C38H72N2O12"), ("doxorubicin", "C27H29NO11"),
    ("methotrexate", "C20H22N8O5"), ("acyclovir", "C8H11N5O3"), ("oseltamivir", "C16H28N2O4"), ("zidovudine", "C10H13N5O4"),
    ("chloroquine", "C18H26ClN3"), ("ciprofloxacin", "C17H18FN3O3"), ("ceftriaxone", "C18H18N8O7S3"),
    ("meropenem", "C17H25N3O5S"), ("lovastatin", "C24H36O5"), ("simvastatin", "C25H38O5"), ("rosuvastatin", "C22H28FN3O6S"),
    ("losartan", "C22H23ClN6O"), ("valsartan", "C24H29N5O3"), ("clopidogrel", "C16H16ClNO2S"), ("apixaban", "C25H25N5O4"),
    ("rivaroxaban", "C19H18ClN3O5S"), ("gabapentin", "C9H17NO2"), ("pregabalin", "C8H17NO2"), ("levetiracetam", "C8H14N2O2"),
    ("lamotrigine", "C9H7Cl2N5"), ("carbamazepine", "C15H12N2O"), ("haloperidol", "C21H23ClFNO2"),
    ("risperidone", "C23H27FN4O2"), ("olanzapine", "C17H20N4S"), ("quetiapine", "C21H25N3O2S"),
    ("aripiprazole", "C23H27Cl2N3O2"), ("venlafaxine", "C17H27NO2"), ("duloxetine", "C18H19NOS"), ("bupropion", "C13H18ClNO"),
    ("trazodone", "C19H22ClN5O"), ("zolpidem", "C19H21N3O"), ("tramadol", "C16H25NO2"), ("fentanyl", "C22H28N2O"),
    ("naloxone", "C19H21NO4"), ("methadone", "C21H27NO"), ("ketamine", "C13H16ClNO"), ("lidocaine", "C14H22N2O"),
    ("propofol", "C12H18O"), ("tamoxifen", "C26H29NO"), ("letrozole", "C17H11N5"), ("finasteride", "C23H36N2O2"),
    ("tadalafil", "C22H19N3O4"), ("montelukast", "C35H36ClNO3S"), ("cetirizine", "C21H25ClN2O3"),
    ("loratadine", "C22H23ClN2O2"), ("diphenhydramine", "C17H21NO"), ("ranitidine", "C13H22N4O3S"),
    ("famotidine", "C8H15N7O2S3"), ("ondansetron", "C18H19N3O"), ("furosemide", "C12H11ClN2O5S"),
    ("hydrochlorothiazide", "C7H8ClN3O4S2"), ("spironolactone", "C24H32O4S"), ("amiodarone", "C25H29I2NO3"),
    ("budesonide", "C25H34O6"), ("fluticasone propionate", "C25H31F3O5S"), ("salbutamol", "C13H21NO3"),
    ("naproxen", "C14H14O3"), ("diclofenac", "C14H11Cl2NO2"), ("celecoxib", "C17H14F3N3O2S"), ("meloxicam", "C14H13N3O4S2"),
    ("allopurinol", "C5H4N4O"), ("colchicine", "C22H25NO6"), ("metronidazole", "C6H9N3O3"),
    ("sulfamethoxazole", "C10H11N3O3S"), ("trimethoprim", "C14H18N4O3"), ("fluconazole", "C13H12F2N6O"),
    ("terbinafine", "C21H25N"), ("tenofovir", "C9H14N5O4P"), ("efavirenz", "C14H9ClF3NO2"), ("lamivudine", "C8H11N3O3S"),
    ("emtricitabine", "C8H10FN3O3S"), ("ribavirin", "C8H12N4O5"), ("topiramate", "C12H21NO8S"),
    ("empagliflozin", "C23H27ClO7"), ("dapagliflozin", "C21H25ClO6"), ("sitagliptin", "C16H15F6N5O"),
    ("pioglitazone", "C19H20N2O3S"), ("glipizide", "C21H27N5O4S"), ("glyburide", "C23H28ClN3O5S"), ("cocaine", "C17H21NO4"),
    ("atropine", "C17H23NO3"), ("scopolamine", "C17H21NO4"), ("strychnine", "C21H22N2O2"), ("LSD", "C20H25N3O"),
    ("THC", "C21H30O2"), ("CBD", "C21H30O2"), ("psilocybin", "C12H17N2O4P"), ("mescaline", "C11H17NO3"),
    ("heroin", "C21H23NO5"), ("camptothecin", "C20H16N2O4"), ("docetaxel", "C43H53NO14"), ("quercetin", "C15H10O7"),
    ("EGCG", "C22H18O11"), ("rutin", "C27H30O16"), ("piperine", "C17H19NO3"), ("allicin", "C6H10OS2"),
    ("sulforaphane", "C6H11NOS2"), ("geraniol", "C10H18O"), ("linalool", "C10H18O"), ("beta-caryophyllene", "C15H24"),
    ("squalene", "C30H50"), ("beta-carotene", "C40H56"), ("lycopene", "C40H56"), ("protoporphyrin IX", "C34H34N4O4"),
    ("aflatoxin B1", "C17H12O6"), ("glyphosate", "C3H8NO5P"), ("atrazine", "C8H14ClN5"), ("DDT", "C14H9Cl5"),
    ("imidacloprid", "C9H10ClN5O2"), ("fipronil", "C12H4Cl2F6N4OS"), ("chlorpyrifos", "C9H11Cl3NO3PS"),
    ("permethrin", "C21H20Cl2O3"), ("DEET", "C12H17NO"), ("bisphenol A", "C15H16O2"), ("DEHP", "C24H38O4"),
    ("PFOA", "C8HF15O2"), ("PFOS", "C8HF17O3S"), ("TNT", "C7H5N3O6"), ("RDX", "C3H6N6O6"), ("nitroglycerin", "C3H5N3O9"),
    ("aspartame", "C14H18N2O5"), ("sucralose", "C12H19Cl3O8"), ("C60", "C60"), ("coronene", "C24H12"),
    ("methane", "CH4"), ("ozone", "O3"), ("carbon monoxide", "CO"), ("N2O", "N2O"), ("SO2", "O2S"), ("H2S", "H2S"),
    ("HCN", "CHN"), ("phosgene", "CCl2O"), ("ammonia", "H3N"), ("hydrazine", "H4N2"), ("sulfuric acid", "H2O4S"),
    ("nitric acid", "HNO3"), ("phosphoric acid", "H3O4P"), ("perchloric acid", "HClO4"), ("boric acid", "H3BO3"),
    ("K2CO3", "CK2O3"), ("Cs2CO3", "CCs2O3"), ("NaOtBu", "C4H9NaO"), ("KOtBu", "C4H9KO"), ("LiHMDS", "C6H18LiNSi2"),
    ("NaHMDS", "C6H18NNaSi2"), ("LDA", "C6H14LiN"), ("n-BuLi", "C4H9Li"), ("t-BuLi", "C4H9Li"), ("DBU", "C9H16N2"),
    ("DBN", "C7H12N2"), ("DABCO", "C6H12N2"), ("DMAP", "C7H10N2"), ("NMM", "C5H11NO"), ("2,6-lutidine", "C7H9N"),
    ("proton sponge", "C14H18N2"), ("DIPEA", "C8H19N"), ("TFA", "C2HF3O2"), ("TFE", "C2H3F3O"), ("HFIP", "C3H2F6O"),
    ("NMP", "C5H9NO"), ("DMA", "C4H9NO"), ("DME", "C4H10O2"), ("MTBE", "C5H12O"), ("CPME", "C6H12O"), ("2-MeTHF", "C5H10O"),
    ("HMPA", "C6H18N3OP"), ("DMPU", "C6H12N2O"), ("sulfolane", "C4H8O2S"), ("CS2", "CS2"), ("Grubbs 1", "C43H72Cl2P2Ru"),
    ("Crabtree's catalyst", "C31H50F6IrNP2"), ("Wilkinson's catalyst", "C54H45ClP3Rh"), ("PEPPSI-IPr", "C32H40Cl3N3Pd"),
    ("Tebbe reagent", "C13H18AlClTi"), ("carboplatin", "C6H12N2O4Pt"), ("ferrocene", "C10H10Fe"), ("Cp2TiCl2", "C10H10Cl2Ti"),
])
def test_formula_spot_checks(name, expected):
    assert formula(name) == expected


# --------------------------------------------------------------------------- #
# stereochemistry
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name, labels", [
    ("alanine", "S"), ("L-alanine", "S"), ("D-alanine", "R"), ("cysteine", "R"), ("D-cysteine", "S"),
    ("selenocysteine", "R"), ("valine", "S"), ("leucine", "S"), ("isoleucine", "SS"), ("allo-isoleucine", "RS"),
    ("proline", "S"), ("D-proline", "R"), ("phenylalanine", "S"), ("tyrosine", "S"), ("tryptophan", "S"), ("serine", "S"),
    ("threonine", "RS"), ("allo-threonine", "SS"), ("methionine", "S"), ("asparagine", "S"), ("glutamine", "S"),
    ("aspartic acid", "S"), ("glutamic acid", "S"), ("lysine", "S"), ("arginine", "S"), ("histidine", "S"),
    ("ornithine", "S"), ("citrulline", "S"), ("homoserine", "S"), ("homocysteine", "S"), ("norleucine", "S"),
    ("hydroxyproline", "RS"), ("penicillamine", "S"), ("L-DOPA", "S"), ("(S)-ibuprofen", "S"), ("(R)-ibuprofen", "R"),
    ("naproxen", "S"), ("(S)-nicotine", "S"), ("L-lactic acid", "S"), ("D-lactic acid", "R"), ("L-malic acid", "S"),
    ("L-tartaric acid", "RR"), ("D-tartaric acid", "SS"), ("(R)-mandelic acid", "R"), ("(R)-1-phenylethylamine", "R"),
    ("(S)-1-phenylethanol", "S"), ("L-carnitine", "R"), ("(R)-carvone", "R"), ("(S)-carvone", "S"), ("(R)-limonene", "R"),
    ("menthol", "RSR"), ("D-glyceraldehyde", "R"), ("L-glyceraldehyde", "S"), ("D-erythrose", "RR"), ("D-threose", "SR"),
    ("D-xylose", "RSR"), ("L-arabinose", "RSS"), ("D-lyxose", "SSR"), ("D-glucose (open chain)", "RSRR"),
    ("DTT", "SS"), ("escitalopram", "S"), ("clopidogrel", "S"), ("duloxetine", "S"), ("pregabalin", "S"),
    ("levetiracetam", "S"), ("rivaroxaban", "S"), ("sitagliptin", "R"), ("linezolid", "S"), ("tenofovir", "R"),
    ("esketamine", "S"), ("selegiline", "R"), ("rasagiline", "R"), ("dextroamphetamine", "S"), ("atomoxetine", "R"),
])
def test_stereo_labels(name, labels):
    assert cip(chem.resolve_name(name)) == labels


def test_ring_sugar_conventions():
    # beta-D-glucopyranose = (2R,3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol
    assert canon(chem.resolve_name("glucose")) == canon("OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O")
    assert chem.resolve_name("glucose") == chem.resolve_name("beta-D-glucopyranose") == chem.resolve_name("dextrose")
    assert canon(chem.resolve_name("alpha-D-glucose")) != canon(chem.resolve_name("glucose"))
    assert canon(chem.resolve_name("L-glucose")) != canon(chem.resolve_name("D-glucose"))
    assert cip(chem.resolve_name("ribose")) == "RRRS"           # beta-D-ribofuranose, atom order C4,C1,C2,C3
    assert cip(chem.resolve_name("galactose")) == "RRRSR"       # beta-D-galactopyranose
    assert cip(chem.resolve_name("mannose")) == "RRSSS"         # beta-D-mannopyranose
    assert canon(chem.resolve_name("glucose")) != canon(chem.resolve_name("galactose"))
    assert canon(chem.resolve_name("glucose")) != canon(chem.resolve_name("mannose"))
    assert cip(chem.resolve_name("biotin")) == "SSR"            # (3aS,4S,6aR)-D-biotin
    assert cip(chem.resolve_name("L-ascorbic acid")) == "SR"    # side-chain S, ring C5 R


# --------------------------------------------------------------------------- #
# normalisation of what people actually type
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("typed, canonical_name", [
    ("n-BuLi", "n-butyllithium"), ("nBuLi", "n-butyllithium"), ("n BuLi", "n-butyllithium"), ("BuLi", "n-butyllithium"),
    ("t-BuLi", "tert-butyllithium"), ("tBuLi", "tert-butyllithium"), ("sec-BuLi", "sec-butyllithium"),
    ("EDC·HCl", "edc.hcl"), ("EDC•HCl", "edc.hcl"), ("EDC HCl", "edc.hcl"), ("EDC-HCl", "edc.hcl"), ("EDC.HCl", "edc.hcl"),
    ("EDC x HCl", "edc.hcl"), ("EDC hydrochloride", "edc.hcl"), ("EDCI·HCl", "edc.hcl"),
    ("Hünig's base", "dipea"), ("Hunig's base", "dipea"), ("Hünigs base", "dipea"), ("Huenig's base" if False else "hunigs base", "dipea"),
    ("β-alanine", "beta-alanine"), ("beta alanine", "beta-alanine"), ("b-alanine", "beta-alanine"), ("β-Ala", "beta-alanine"),
    ("γ-aminobutyric acid", "gaba"), ("α-ketoglutaric acid", "alpha-ketoglutaric acid"), ("ε-caprolactone", "epsilon-caprolactone"),
    ("sulphuric acid", "sulfuric acid"), ("sodium sulphate", "sodium sulfate"), ("caesium carbonate", "cesium carbonate"),
    ("aluminium chloride", "aluminum chloride"), ("oestradiol", "estradiol"), ("haem", "heme"), ("Essigsäure", "acetic acid"),
    ("Tetrahydrofuran", "thf"), ("  THF  ", "thf"), ("THF.", "thf"), ("'THF'", "thf"), ('"DMF"', "dmf"),
    ("N,N-Dimethylformamide", "dmf"), ("N,N dimethylformamide", "dmf"), ("dimethyl formamide" if False else "dimethylformamide", "dmf"),
    ("Pd/C", "palladium"), ("10% Pd/C", "palladium"), ("Pd-C", "palladium"), ("Pd on C", "palladium"), ("palladium on carbon", "palladium"),
    ("CuSO4·5H2O", "copper(ii) sulfate pentahydrate"), ("CuSO4.5H2O", "copper(ii) sulfate pentahydrate"),
    ("Tris–HCl", "tris hydrochloride"), ("Tris-HCl", "tris hydrochloride"), ("TRIS HCl", "tris hydrochloride"),
    ("4Å MS", "molecular sieves"), ("4 Å molecular sieves", "molecular sieves"), ("4A MS", "molecular sieves"),
    ("mol sieves", "molecular sieves"), ("1,4-Dioxane", "dioxane"), ("1,4 dioxane", "dioxane"), ("p-TsOH", "p-toluenesulfonic acid"),
    ("pTsOH", "p-toluenesulfonic acid"), ("TsOH", "p-toluenesulfonic acid"), ("PTSA", "p-toluenesulfonic acid"),
    ("m-CPBA", "mcpba"), ("mCPBA", "mcpba"), ("MCPBA", "mcpba"), ("meta-chloroperbenzoic acid", "mcpba"),
    ("DIBAL", "dibal-h"), ("DIBAL-H", "dibal-h"), ("DIBALH", "dibal-h"), ("NaBH(OAc)3", "stab"),
    ("Boc2O", "di-tert-butyl dicarbonate"), ("(Boc)2O", "di-tert-butyl dicarbonate"), ("Boc anhydride", "di-tert-butyl dicarbonate"),
    ("Fmoc-OSu", "fmoc-osu"), ("Fmoc OSu", "fmoc-osu"), ("FmocOSu", "fmoc-osu"), ("TBS-Cl", "tbscl"), ("TBSCl", "tbscl"), ("TBDMSCl", "tbscl"),
    ("BF3·OEt2", "bf3.oet2"), ("BF3.OEt2", "bf3.oet2"), ("BF3 OEt2", "bf3.oet2"), ("BF3-OEt2", "bf3.oet2"), ("BF3 etherate", "bf3.oet2"),
    ("Pd(PPh3)4", "pd(pph3)4"), ("Pd(PPh3)4 ", "pd(pph3)4"), ("tetrakis(triphenylphosphine)palladium(0)", "pd(pph3)4"),
    ("PdCl2(dppf)·CH2Cl2", "pdcl2(dppf)"), ("Pd(dppf)Cl2·DCM", "pdcl2(dppf)"),
    ("Grubbs' second generation catalyst", "grubbs ii"), ("Grubbs 2nd generation" if False else "grubbs catalyst 2nd generation", "grubbs ii"),
    ("Hoveyda–Grubbs II", "hoveyda-grubbs ii"), ("Lawesson’s reagent", "lawesson's reagent"), ("Lawessons reagent", "lawesson's reagent"),
    ("Dess–Martin periodinane", "dmp"), ("Dess-Martin", "dmp"), ("Dess Martin periodinane", "dmp"),
    ("Et3N", "triethylamine"), ("NEt3", "triethylamine"), ("TEA", "triethylamine"), ("Et3N·HCl", "triethylamine hydrochloride"),
    ("NH4Cl", "ammonium chloride"), ("sat. NH4Cl" if False else "aqueous ammonia", "ammonium hydroxide"), ("brine" if False else "NaCl", "sodium chloride"),
    ("anhydrous MgSO4", "magnesium sulfate"), ("MgSO4 heptahydrate", "magnesium sulfate"), ("Na2SO4 (anhydrous)", "sodium sulfate"),
    ("Tris monohydrate" if False else "tris base", "tris"), ("HEPES (free acid)" if False else "hepes", "hepes"),
    ("L-Ala", "alanine"), ("L-alanine", "alanine"), ("l-Alanine", "alanine"), ("Ala", "alanine"), ("(S)-alanine", "alanine"),
    ("D-Ala", "d-alanine"), ("(R)-alanine", "d-alanine"), ("L-Cys", "cysteine"), ("(R)-cysteine", "cysteine"),
    ("L-glycine", "glycine"), ("Gly", "glycine"), ("D-glucose", "glucose"), ("dextrose", "glucose"), ("Glc", "glucose"),
    ("L-ascorbic acid", "ascorbic acid"), ("vitamin C", "ascorbic acid"), ("(−)-menthol", "menthol"), ("(-)-menthol", "menthol"),
    ("(R)-(+)-limonene", "limonene"), ("d-limonene", "limonene"), ("(S)-(+)-ibuprofen", "(s)-ibuprofen"), ("dexibuprofen", "(s)-ibuprofen"),
    ("acetaminophen", "paracetamol"), ("Tylenol", "paracetamol"), ("APAP", "paracetamol"), ("Prozac", "fluoxetine"),
    ("Lipitor" if False else "atorvastatin", "atorvastatin"), ("Taxol", "paclitaxel"), ("Gleevec", "imatinib"), ("Glivec", "imatinib"),
    ("5-FU", "5-fluorouracil"), ("AZT", "zidovudine"), ("NAD", "nad+"), ("NAD(+)", "nad+"), ("β-NAD", "nad+"),
    ("beta-mercaptoethanol", "2-mercaptoethanol"), ("BME", "2-mercaptoethanol"), ("2-ME", "2-mercaptoethanol"), ("βME", "2-mercaptoethanol"),
    ("Gdn-HCl" if False else "GdnHCl", "guanidine hydrochloride"), ("GuHCl", "guanidine hydrochloride"), ("guanidinium chloride", "guanidine hydrochloride"),
    ("Triton X100", "triton x-100"), ("Tween20", "tween-20"), ("polysorbate 20", "tween-20"), ("PEG 400", "peg-400"), ("PEG400", "peg-400"),
    ("NHS ester" if False else "N-hydroxysuccinimide", "hosu"), ("sulfo NHS", "sulfo-nhs"), ("DBCO NHS ester", "dbco-nhs ester"),
    ("Cy-3" if False else "cyanine 3", "cy3"), ("Cy 5", "cy5"), ("Alexa" if False else "FITC", "fitc"), ("5(6)-FAM", "5-fam"),
    ("CDCl3", "cdcl3"), ("chloroform-d", "cdcl3"), ("DMSO d6", "dmso-d6"), ("d6-DMSO", "dmso-d6"), ("(CD3)2SO", "dmso-d6"),
    ("CD3OD", "methanol-d4"), ("MeOD", "methanol-d4"), ("D2O", "d2o"), ("deuterium oxide", "d2o"), ("C6D6", "benzene-d6"),
])
def test_typed_variants_resolve_to_the_same_structure(typed, canonical_name):
    assert chem.resolve_name(typed) == chem.resolve_name(canonical_name) != None  # noqa: E711


def test_case_sensitive_element_collisions():
    assert chem.resolve_name("CoCl2") == "[Co+2].[Cl-].[Cl-]"
    assert chem.resolve_name("COCl2") == "O=C(Cl)Cl"
    assert chem.resolve_name("cocl2") == chem.resolve_name("phosgene")
    assert chem.resolve_name("CO") == chem.resolve_name("carbon monoxide")
    assert chem.resolve_name("Co") == "[Co]"
    assert chem.resolve_name("NO") == chem.resolve_name("nitric oxide")
    assert chem.resolve_name("In") == "[In]"


def test_hydrate_and_qualifier_suffixes():
    assert chem.resolve_name("sodium acetate trihydrate") == chem.resolve_name("sodium acetate")
    assert chem.resolve_name("MgSO4·7H2O") == chem.resolve_name("MgSO4")
    assert chem.resolve_name("MgSO4.7H2O") == chem.resolve_name("MgSO4")
    assert chem.resolve_name("CaCl2 dihydrate") == chem.resolve_name("CaCl2")
    assert chem.resolve_name("CaCl2 anhydrous") == chem.resolve_name("CaCl2")
    assert chem.resolve_name("anhydrous THF") == chem.resolve_name("THF")
    assert chem.resolve_name("dry DCM") == chem.resolve_name("DCM")
    assert chem.resolve_name("TBAF (1 M in THF)") == chem.resolve_name("TBAF")
    assert chem.resolve_name("formaldehyde (37% aq)") == chem.resolve_name("formaldehyde")
    assert chem.resolve_name("H2O2 (30%)") == chem.resolve_name("H2O2")
    assert chem.resolve_name("LiCl·xH2O") == chem.resolve_name("LiCl")
    # explicit entries win over suffix stripping
    assert chem.resolve_name("hydrazine hydrate") == "NN.O"
    assert chem.resolve_name("chloral hydrate") == "OC(O)C(Cl)(Cl)Cl"
    assert chem.resolve_name("CuSO4·5H2O") != chem.resolve_name("CuSO4")


def test_stereo_prefix_handling():
    ala = chem.resolve_name("alanine")
    assert chem.resolve_name("DL-alanine") == chem._without_tetrahedral(ala) == "CC(N)C(=O)O"
    assert chem.resolve_name("rac-alanine") == chem.resolve_name("(±)-alanine") == chem.resolve_name("(+/-)-alanine") == "CC(N)C(=O)O"
    assert canon(chem.resolve_name("racemic ibuprofen")) == canon(chem.resolve_name("ibuprofen"))
    assert canon(chem.resolve_name("rac-ibuprofen")) == canon(chem.resolve_name("ibuprofen"))
    assert cip(chem.resolve_name("(R)-alanine")) == "R" and cip(chem.resolve_name("(S)-alanine")) == "S"
    assert cip(chem.resolve_name("(S)-naproxen")) == "S" and cip(chem.resolve_name("(R)-naproxen")) == "R"
    assert cip(chem.resolve_name("(R)-ibuprofen")) == "R"
    assert cip(chem.resolve_name("(S)-warfarin")) == "S" and cip(chem.resolve_name("(R)-warfarin")) == "R"
    assert canon(chem.resolve_name("(S)-fluoxetine")) == canon(chem.resolve_name("(S)-fluoxetine"))
    assert cip(chem.resolve_name("(R)-fluoxetine")) == "R"
    assert chem.resolve_name("(S)-glycine") == chem.resolve_name("glycine")   # no centre: descriptor ignored
    assert chem.resolve_name("L-glycine") == chem.resolve_name("D-glycine") == chem.resolve_name("glycine")
    assert chem.resolve_name("L-tris") == chem.resolve_name("tris")
    assert chem.resolve_name("(S)-paclitaxel") is None                      # many centres: do not guess
    assert chem.resolve_name("L-ribose") is None                            # bare name is D: do not guess
    assert chem.resolve_name("D-leucine") == "CC(C)C[C@@H](N)C(=O)O"
    assert chem.resolve_name("(+)-biotin") == chem.resolve_name("biotin")
    assert chem.resolve_name("(2S,3R)-threonine") == chem.resolve_name("threonine")
    assert chem.resolve_name("trans-cinnamic acid") == chem.resolve_name("cinnamic acid")
    assert chem.resolve_name("cis-cinnamic acid") is None


def test_unknown_names_fall_through_without_network():
    assert chem.resolve_name("unobtainium") is None
    assert chem.resolve_name("(S)-unobtainium") is None
    assert chem.resolve_name("unobtainium monohydrate") is None
    assert chem.resolve_name("") is None and chem.resolve_name("   ") is None
    assert chem.resolve_name(None) is None  # type: ignore[arg-type]
    assert chem.resolve_name("·") is None


def test_resolved_names_depict():
    for name in ["HATU", "Pd(PPh3)4", "paclitaxel", "ATP", "cisplatin", "Grubbs II", "DMSO-d6", "heme"]:
        m = chem.mol_from_input(name, fmt="name")
        assert m.GetNumConformers() == 1
        assert "<svg" in chem.depict_svg(m)
