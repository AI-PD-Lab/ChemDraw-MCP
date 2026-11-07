"""Low-level CDXMLDocument use: coloured fragments, an electron-pushing (curved) arrow, captions with
bold / italic / formula faces, a rounded frame and ChemDraw's dark-mode style.

Coordinates are absolute page positions in points (1/72 inch); y grows downwards.  Every add_* call returns
an ObjRef with the object's id and bounding box, so later objects can be placed relative to earlier ones.
Run:   CHEMDRAWMCP_BACKEND=none python examples/07_custom_cdxml.py   (writes <output dir>/07_custom.cdxml)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chemdrawmcp import cdxml, chem, discovery  # noqa: E402
from chemdrawmcp.cdxml import CDXMLDocument  # noqa: E402


def main() -> int:
    # "dark" = ChemDraw's Dark Mode stationery: dark page background, near-white foreground, 1.6 pt lines.
    # Any style key from styles.py works (acs, nature, rsc, wiley, slide ...); `settings` overrides single values.
    doc = CDXMLDocument("dark", name="custom document", settings={"caption_size": 9})
    L = doc.bond_length                                        # 14.4 pt in this style
    x, y = 72.0, 150.0                                         # left edge and vertical centre of the first row

    # 1. Fragments.  anchor="left" puts the left edge at x and the centre line at y.  `color` colours the whole
    #    fragment; `atom_colors` maps RDKit atom indices to colours ("red", "#ff8800", (r, g, b) or an index).
    hydroxide = doc.add_fragment(chem.mol_from_input("[OH-]"), x, y, anchor="left", color="cyan")
    plus = doc.add_plus(hydroxide.box.x1 + 0.6 * L, y)         # the "+" between reaction components
    acetone = doc.add_fragment(chem.mol_from_input("CC(C)=O"), plus.box.x1 + 0.6 * L, y, anchor="left",
                               atom_colors={1: "orange", 3: "red"})   # 1 = carbonyl C, 3 = O (RDKit indices)

    # 2. A straight reaction arrow: tail (x0, y0) -> head (x1, y1).  Kinds: full, equilibrium, retrosynthetic,
    #    double, dashed, bold, hollow, nogo, half_left, ...  The arrow's box is used to centre the condition text.
    arrow = doc.add_arrow(acetone.box.x1 + L, y, acetone.box.x1 + 4.5 * L, y, kind="full")
    conditions = doc.add_text("NaOH, H2O, 25 °C", arrow.box.cx, arrow.box.y0 - 0.2 * L, anchor="baseline",
                              formula=True)                    # formula face: digits become subscripts (H₂O)
    product = doc.add_fragment(chem.mol_from_input("CC(C)(O)[O-]"), arrow.box.x1 + L, y, anchor="left",
                               atom_colors={4: "cyan", 3: "red"})

    # 3. An electron-pushing arrow: a circular arc from the hydroxide oxygen towards the carbonyl carbon.
    #    `angle` is the sweep in degrees; positive bends to the left of the tail->head direction (counter-
    #    clockwise on screen), negative to the right.  Tail/head are taken from the fragment boxes here.
    doc.add_curved_arrow(hydroxide.box.cx + 0.2 * L, hydroxide.box.y0, acetone.box.x0 + 1.2 * L,
                         acetone.box.cy - 0.6 * L, angle=-70, color="magenta")

    # 4. Captions.  face = "bold" | "italic" | ["bold", "italic"] | an int bit mask (1 bold, 2 italic, 4 underline,
    #    96 formula).  anchor="top" positions the top of the text block; justification Left/Center/Right.
    title = doc.add_text("Nucleophilic addition of hydroxide to acetone", x, y - 5 * L, size=12, face="bold",
                         justification="Left", anchor="top", color="yellow")
    doc.add_text("the lone pair of the nucleophile attacks the carbonyl carbon; the C=O pi bond moves onto oxygen",
                 x, title.box.y1 + 0.3 * L, size=8, face="italic", justification="Left", anchor="top")
    doc.add_text("tetrahedral alkoxide", product.box.cx, product.box.y1 + 0.4 * L, anchor="top", face=["bold", "italic"])
    doc.add_text("Keq ~ 10-3", product.box.cx, product.box.y0 - 0.4 * L, anchor="baseline", formula=True, color="gray")

    # 5. Reaction semantics: ChemDraw (and read_cdxml) learn what reacts to what from a <scheme>/<step> record.
    doc.add_step([hydroxide, acetone], [product], arrows=[arrow], above=[conditions], plusses=[plus])

    # 6. A rounded rectangle around everything drawn so far.  content_box() is the union of all object boxes.
    #    shape: rectangle | rounded | oval | circle | line; style: Plain | Dashed | Bold | Filled | Shaded.
    frame = doc.content_box().padded(1.2 * L)
    doc.add_graphic(frame.x0, frame.y0, frame.x1, frame.y1, shape="rounded", style="Plain", color="gray", line_width=1.2)

    out = Path(discovery.output_dir()) / "07_custom.cdxml"
    out.write_text(doc.to_xml(), encoding="utf-8")             # to_xml tiles pages so ChemDraw never asks to expand
    print("wrote", out)
    back = cdxml.read_cdxml(out.read_text(encoding="utf-8"))
    print("read back: %d molecules, %d captions, %d arrows, %d graphics, %d colours in the colour table" % (
        len(back["molecules"]), len(back["texts"]), len(back["arrows"]), len(back["graphics"]), len(doc.colors)))
    print("step:", back["reaction_smiles"][0], "| conditions:", [t["text"] for t in back["steps"][0]["above_arrow"]])
    print("curved arrow sweep:", [a.get("angularSize") for a in back["arrows"] if a.get("angularSize")], "degrees")
    return 0


if __name__ == "__main__":
    sys.exit(main())
