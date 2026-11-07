"""Four-step synthesis of methyl benzoate as a multi-step scheme with a title, wrapped at 400 pt.

benzene -> bromobenzene -> phenylmagnesium bromide -> benzoic acid -> methyl benzoate

Consecutive steps that start with the previous product are chained: the intermediate is drawn once.  When a
row would exceed max_width the scheme wraps and the intermediate is repeated at the start of the next row.
Run:   CHEMDRAWMCP_BACKEND=none python examples/03_scheme.py      (writes <output dir>/03_scheme.cdxml)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chemdrawmcp import cdxml, chem, discovery, layout  # noqa: E402
from chemdrawmcp.layout import Component, Step  # noqa: E402


def step(reactant: str, product: str, above, below, label=None, arrow="full") -> Step:
    """One A -> B step; names are resolved by the offline dictionary, everything else is SMILES."""
    return Step(reactants=[Component(chem.mol_from_input(reactant))],
                products=[Component(chem.mol_from_input(product), label=label)],
                above=list(above), below=list(below), arrow=arrow)


STEPS = [
    step("benzene", "bromobenzene", ["Br2"], ["FeBr3, 0 °C"], label="1  (92%)"),
    step("bromobenzene", "Br[Mg]c1ccccc1", ["Mg"], ["THF, reflux"], label="2"),
    step("Br[Mg]c1ccccc1", "benzoic acid", ["1. CO2", "2. H3O+"], [], label="3  (78%)"),
    step("benzoic acid", "COC(=O)c1ccccc1", ["MeOH, H2SO4 (cat.)"], ["reflux, 3 h"], label="4  (90%)"),
]


def main() -> int:
    doc = layout.build_scheme(STEPS, style="acs", title="Scheme 1.  Methyl benzoate from benzene", max_width=400,
                              name="methyl benzoate")
    out = Path(discovery.output_dir()) / "03_scheme.cdxml"
    out.write_text(doc.to_xml(), encoding="utf-8")
    print("wrote", out)

    back = cdxml.read_cdxml(out.read_text(encoding="utf-8"))
    rows = sorted({round(a["bounds"]["y0"]) for a in back["arrows"]})
    print("%d molecules drawn for %d steps (intermediates shared), %d rows, %d pages" % (
        len(back["molecules"]), len(back["steps"]), len(rows), back["pages"]))
    print("title:", repr(back["texts"][0]["text"]))
    for i, (st, smi) in enumerate(zip(back["steps"], back["reaction_smiles"]), 1):
        print("\nstep %d: %s -> %s" % (i, " + ".join(r["formula"] for r in st["reactants"]),
                                       " + ".join(p["formula"] for p in st["products"])))
        print("   above: %-24s below: %s" % ([t["text"] for t in st["above_arrow"]], [t["text"] for t in st["below_arrow"]]))
        print("   reaction SMILES:", smi)
    return 0


if __name__ == "__main__":
    sys.exit(main())
