"""RDKit properties, elemental analysis and reaction balance, entirely offline.

Inputs may be SMILES, InChI, molblocks or chemical names; names are resolved by the built-in dictionary
(then OPSIN if installed; PubChem only with CHEMDRAWMCP_ONLINE=1).  Reactions accept reaction SMILES or
"A + B --[reagents]--> C" text with SMILES or names.
Run:   CHEMDRAWMCP_ONLINE=0 python examples/05_properties.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chemdrawmcp import chem  # noqa: E402

STRUCTURES = [
    "aspirin",                                       # a name, resolved offline
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",                     # ibuprofen as SMILES
    "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",   # caffeine as InChI
]
REACTIONS = [
    "acetic acid + ethanol --[H2SO4]--> ethyl acetate + O",   # names and SMILES mixed; bracket = conditions
    "CC(=O)O.OCC>[H+]>CC(=O)OCC.O",                            # reaction SMILES with an agent
    "benzene + BrBr -> bromobenzene",                           # not balanced: HBr is missing
]
KEYS = ("formula", "molecular_weight", "exact_mass", "inchikey", "logp", "tpsa", "hbd", "hba",
        "rotatable_bonds", "lipinski_violations", "qed")


def main() -> int:
    print("name resolution: aspirin ->", chem.resolve_name("aspirin"), "| detect_format:",
          chem.detect_format("aspirin"), "/", chem.detect_format(STRUCTURES[2])[:5])
    for src in STRUCTURES:
        mol = chem.mol_from_input(src)
        p = chem.properties(mol, extended=True)
        print("\n" + chem.describe(mol))
        for k in KEYS:
            print("  %-20s %s" % (k, p[k]))
        print("  %-20s %s" % ("elemental analysis", chem.elemental_analysis(mol)))
        st = p["stereocenters"]
        print("  %-20s %d assigned, %d unassigned" % ("stereocentres", st["assigned"], st["unassigned"]))

    for text in REACTIONS:
        rxn = chem.parse_reaction(text)
        r = chem.reaction_properties(rxn)
        print("\n%s\n  reaction SMILES: %s" % (text, r["reaction_smiles"]))
        if r["conditions"]:
            print("  conditions:      %s" % r["conditions"])
        print("  components:      %s -> %s" % (" + ".join(x["formula"] for x in r["reactants"]),
                                                " + ".join(x["formula"] for x in r["products"])))
        print("  balanced:        %s (reactants %.3f g/mol, products %.3f g/mol)" % (
            r["balanced"], r["mass_reactants"], r["mass_products"]))
        if r["unbalanced_elements"]:
            print("  unbalanced:      %s" % {el: r["element_balance"][el] for el in r["unbalanced_elements"]})
        if r["atom_economy"] is not None:
            print("  atom economy:    %.1f %% (first product / all reactants)" % r["atom_economy"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
