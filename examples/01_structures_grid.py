"""Labelled grid of six drug molecules in ACS style, built with the library API (no ChemDraw needed).

Run:   CHEMDRAWMCP_BACKEND=none python examples/01_structures_grid.py
Writes <output dir>/01_structures.cdxml, prints the layout summary and reads the file back with read_cdxml.
The output dir is $CHEMDRAWMCP_OUTPUT_DIR (default ~/.chemdrawmcp/renders).  Open the .cdxml in ChemDraw to see it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # run from a checkout without installing
from chemdrawmcp import cdxml, chem, discovery, layout  # noqa: E402

# names come from the offline dictionary (no network); anything RDKit reads works too (SMILES, InChI, molblock)
DRUGS = [
    ("aspirin", "aspirin"),
    ("ibuprofen", "ibuprofen"),
    ("paracetamol", "paracetamol"),
    ("caffeine", "caffeine"),
    ("naproxen", "COc1ccc2cc(ccc2c1)[C@H](C)C(=O)O"),
    ("diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
]


def main() -> int:
    mols = [chem.mol_from_input(src) for _, src in DRUGS]          # sanitised, 2D coords, Kekulé, wedges
    labels = ["%d  %s" % (i + 1, name) for i, (name, _) in enumerate(DRUGS)]
    doc = layout.build_structures(mols, labels=labels, style="acs", columns=3, name="drug grid")

    out = Path(discovery.output_dir()) / "01_structures.cdxml"
    out.write_text(doc.to_xml(), encoding="utf-8")
    print("wrote", out)

    s = doc.summary()
    b = s["bounds"]
    print("\nlayout: style %r, bond length %.1f pt, %d objects, content %.0f x %.0f pt" % (
        s["style"], s["bond_length"], len(s["objects"]), b["width"], b["height"]))
    for o in s["objects"]:
        bb = o["bounds"]
        what = o.get("smiles") or o.get("text") or ""
        print("  %-8s id=%-3d (%6.1f, %6.1f)  %5.1f x %5.1f  %s" % (
            o["kind"], o["id"], bb["x0"], bb["y0"], bb["width"], bb["height"], what))

    back = cdxml.read_cdxml(out.read_text(encoding="utf-8"))
    print("\nread back: %d molecules, %d captions, bond length %s, font %s" % (
        len(back["molecules"]), len(back["texts"]), back["settings"].get("BondLength"), back["settings"].get("LabelFont")))
    for m, cap in zip(back["molecules"], back["texts"]):
        print("  %-14s %-10s MW %8.3f  %s" % (cap["text"], m["formula"], m["molecular_weight"], m["smiles"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
