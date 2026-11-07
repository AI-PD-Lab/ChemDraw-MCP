"""Parse a CDXML file and print molecules, captions, arrows and reaction steps as tables.

usage: python examples/04_read_cdxml.py [FILE.cdxml]
Default: the scheme written by 03_scheme.py (it is generated first if missing).  Works on files saved by
ChemDraw as well as on files written by this package: cdxml.read_cdxml converts every top-level <fragment>
back into an RDKit molecule (wedges -> stereo), collects <t> captions, <arrow>/<graphic> objects and the
<scheme>/<step> records, and derives reaction SMILES per step.
"""
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from chemdrawmcp import cdxml, discovery  # noqa: E402


def main(argv) -> int:
    if len(argv) > 1:
        path = Path(argv[1]).expanduser()
    else:
        path = Path(discovery.output_dir()) / "03_scheme.cdxml"
        if not path.exists():
            runpy.run_path(str(HERE / "03_scheme.py"), run_name="__main__")
    text = path.read_text(encoding="utf-8")
    doc = cdxml.read_cdxml(text)
    s = doc["settings"]
    print("%s: %s page(s), created by %s, bond length %s pt, fonts %s / %s" % (
        path.name, doc["pages"], s.get("CreationProgram", "?"), s.get("BondLength"), s.get("LabelFont"), s.get("CaptionFont")))

    print("\nmolecules (%d)" % len(doc["molecules"]))
    print("  %-4s %-12s %9s %5s  %s" % ("id", "formula", "MW", "heavy", "SMILES"))
    for m in doc["molecules"]:
        if "error" in m:
            print("  %-4s %s" % (m["id"], m["error"]))
            continue
        print("  %-4s %-12s %9.3f %5d  %s" % (m["id"], m["formula"], m["molecular_weight"], m["heavy_atoms"], m["smiles"]))

    print("\ncaptions (%d)" % len(doc["texts"]))
    for t in doc["texts"]:
        b = t["bounds"] or {}
        print("  %-4s at (%6.1f, %6.1f)  %s" % (t["id"], b.get("x0", 0), b.get("y0", 0), t["text"].replace("\n", " / ")))

    print("\narrows (%d)" % len(doc["arrows"]))
    for a in doc["arrows"]:
        kind = ", ".join("%s=%s" % (k, v) for k, v in a.items() if k not in ("id", "kind", "bounds", "head3D", "tail3D"))
        print("  %-4s %s -> %s  %s" % (a["id"], a.get("tail3D", "?")[:-2], a.get("head3D", "?")[:-2], kind))
    if doc["graphics"]:
        print("\ngraphics (%d): %s" % (len(doc["graphics"]), ", ".join("%s %s" % (g["id"], g["type"]) for g in doc["graphics"])))

    print("\nsteps (%d)" % len(doc["steps"]))
    for i, st in enumerate(doc["steps"]):
        print("  %-4s %s  ->  %s" % (st["id"], " + ".join(r.get("formula") or r["id"] for r in st["reactants"]),
                                      " + ".join(p.get("formula") or p["id"] for p in st["products"])))
        print("       above %-28s below %s" % ([t.get("text") for t in st["above_arrow"]], [t.get("text") for t in st["below_arrow"]]))
        print("       rxn SMILES:", doc["reaction_smiles"][i])
    if doc["warnings"]:
        print("\nwarnings:", *doc["warnings"], sep="\n  ")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
