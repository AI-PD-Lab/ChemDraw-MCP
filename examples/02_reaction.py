"""Fischer esterification as one reaction in Nature style: reagents above the arrow, conditions below,
a yield under the ester and an equilibrium arrow.

Run:   CHEMDRAWMCP_BACKEND=none python examples/02_reaction.py
Writes <output dir>/02_reaction.cdxml (the ChemDraw document, with a <scheme>/<step> record so ChemDraw knows
what reacts to what) and 02_reaction.png, an RDKit preview rendered with chem.depict_png.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chemdrawmcp import cdxml, chem, discovery, layout  # noqa: E402
from chemdrawmcp.layout import Component, Step  # noqa: E402


def main() -> int:
    acid, alcohol = chem.mol_from_input("acetic acid"), chem.mol_from_input("ethanol")
    ester, water = chem.mol_from_input("ethyl acetate"), chem.mol_from_input("O")

    step = Step(
        reactants=[Component(acid), Component(alcohol)],
        products=[Component(ester, label="85%"), Component(water)],   # label = caption under the structure
        above=["H2SO4 (cat.)"],        # "H2SO4" looks like a formula, so ChemDraw's formula face is used (H₂SO₄)
        below=["reflux, 4 h"],
        arrow="equilibrium",           # full | equilibrium | retrosynthetic | double | dashed | bold | nogo | hollow
    )
    doc = layout.build_reaction(step, style="nature", name="esterification")

    out_dir = Path(discovery.output_dir())
    cdxml_path = out_dir / "02_reaction.cdxml"
    cdxml_path.write_text(doc.to_xml(), encoding="utf-8")
    print("wrote", cdxml_path)

    back = cdxml.read_cdxml(cdxml_path.read_text(encoding="utf-8"))
    st = back["steps"][0]
    arrow = back["arrows"][0]
    print("style:          %s (bond length %s pt)" % (doc.style["name"], back["settings"]["BondLength"]))
    print("reaction SMILES:", back["reaction_smiles"][0])
    print("reactants:      ", [r["formula"] for r in st["reactants"]])
    print("products:       ", [p["formula"] for p in st["products"]])
    print("above / below:  ", [t["text"] for t in st["above_arrow"]], [t["text"] for t in st["below_arrow"]])
    print("arrow heads:     head=%s tail=%s (equilibrium)" % (arrow.get("arrowheadHead"), arrow.get("arrowheadTail")))
    print("captions:       ", [t["text"] for t in back["texts"]])

    # RDKit preview: a chem.Reaction can be drawn directly (agents are not drawn by RDKit)
    rxn = chem.Reaction(reactants=[acid, alcohol], products=[ester, water])
    png_path = out_dir / "02_reaction.png"
    png_path.write_bytes(chem.depict_png(rxn, width=900, height=300))
    print("preview PNG:    ", png_path, "(%d bytes)" % png_path.stat().st_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
