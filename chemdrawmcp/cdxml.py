"""CDXML writer and reader.

CDXML is ChemDraw's XML document format.  This module builds complete CDXML
documents (molecules, reaction arrows, captions, plus signs, reaction-step
records) from RDKit molecules, and parses CDXML back into a structured
description.  Everything is in points (1/72 inch), y grows downwards.

Writer essentials (verified against files written by ChemDraw 23):

* ``<n>`` atoms: ``Element`` (omitted = carbon), ``NumHydrogens``, ``Charge``,
  ``Isotope``, ``Radical``; ``<b>`` bonds: ``B``/``E`` atom ids, ``Order``
  (1/2/3/1.5/dative), ``Display`` (WedgeBegin, WedgedHashBegin, Wavy, Dash, Bold...).
* ``<t>`` text: ``p`` = baseline start, ``Justification``; style runs ``<s font= size= face=>``;
  face bits 1 bold, 2 italic, 4 underline, 32 subscript, 64 superscript, 96 = "formula"
  (automatic sub/superscripts for digits/charges).
* ``<arrow>`` with ``Tail3D``/``Head3D`` for straight reaction arrows; ``ArrowheadHead``
  / ``ArrowheadTail`` Full | HalfLeft | HalfRight; ``ArrowheadType`` Solid | Hollow | Angle;
  ``ArrowShaftSpacing`` for equilibrium/retro double shafts; ``LineType`` Dashed | Bold;
  ``NoGo`` Cross | Hash; ``AngularSize`` for arcs.
* colour indices: 0 = black and 1 = white are implicit; ``<colortable>`` entries are 2, 3, ...
* ``<scheme><step ReactionStepReactants= ReactionStepProducts= ReactionStepArrows=
  ReactionStepObjectsAboveArrow= ReactionStepObjectsBelowArrow= ReactionStepPlusses=
  ReactionStepAtomMap=/></scheme>`` records reaction semantics.
"""
from __future__ import annotations

import math
import re
import statistics
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape, quoteattr

from . import textmetrics as _tm
from .styles import get_style

PT_PER_INCH = 72.0

# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------


@dataclass
class Box:
    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def width(self) -> float:
        return self.x1 - self.x0

    @property
    def height(self) -> float:
        return self.y1 - self.y0

    @property
    def cx(self) -> float:
        return (self.x0 + self.x1) / 2.0

    @property
    def cy(self) -> float:
        return (self.y0 + self.y1) / 2.0

    def translated(self, dx: float, dy: float) -> "Box":
        return Box(self.x0 + dx, self.y0 + dy, self.x1 + dx, self.y1 + dy)

    def union(self, other: Optional["Box"]) -> "Box":
        if other is None:
            return Box(self.x0, self.y0, self.x1, self.y1)
        return Box(min(self.x0, other.x0), min(self.y0, other.y0), max(self.x1, other.x1), max(self.y1, other.y1))

    def padded(self, p: float) -> "Box":
        return Box(self.x0 - p, self.y0 - p, self.x1 + p, self.y1 + p)

    def attr(self) -> str:
        return "%s %s %s %s" % (_f(self.x0), _f(self.y0), _f(self.x1), _f(self.y1))

    def to_dict(self) -> dict:
        return {"x0": round(self.x0, 2), "y0": round(self.y0, 2), "x1": round(self.x1, 2), "y1": round(self.y1, 2),
                "width": round(self.width, 2), "height": round(self.height, 2)}


def _f(v: float) -> str:
    """Format a float compactly for XML attributes."""
    if abs(v - round(v)) < 1e-6:
        return str(int(round(v)))
    return ("%.2f" % v).rstrip("0").rstrip(".")


def union_boxes(boxes: Iterable[Optional[Box]]) -> Optional[Box]:
    out: Optional[Box] = None
    for b in boxes:
        if b is None:
            continue
        out = b.union(out) if out else Box(b.x0, b.y0, b.x1, b.y1)
    return out


# ---------------------------------------------------------------------------
# text metrics (approximate: ChemDraw re-measures text when it opens the file)
# ---------------------------------------------------------------------------

# Arial advance widths (em) used when no font can be measured; see textmetrics.ARIAL_WIDTHS
_CHAR_W = _tm.ARIAL_WIDTHS


def text_width(text: str, size: float, font: Optional[str] = None, formula: bool = False) -> float:
    """Width (points) of a single line: measured with a real font when one is found, Arial heuristic otherwise."""
    return _tm.text_width(text, size, font, formula)


def text_box(text: str, size: float, font: Optional[str] = None, formula: bool = False) -> Tuple[float, float, float]:
    """(width, ascent, descent) in points for a (possibly multi-line) text (``formula``: ChemDraw formula face)."""
    return _tm.text_box(text, size, font, formula)


# ---------------------------------------------------------------------------
# element data
# ---------------------------------------------------------------------------

SYMBOLS = ["*", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
           "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",
           "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
           "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
           "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
           "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl",
           "Mc", "Lv", "Ts", "Og"]
ATOMIC_NUMBER = {s: i for i, s in enumerate(SYMBOLS)}

ARROW_KINDS = {
    # name: attributes on the <arrow> element
    "full": {"ArrowheadHead": "Full"},
    "forward": {"ArrowheadHead": "Full"},
    "plain": {"ArrowheadHead": "Full"},
    "equilibrium": {"ArrowheadHead": "HalfLeft", "ArrowheadTail": "HalfLeft", "ArrowShaftSpacing": "300"},
    "double": {"ArrowheadHead": "Full", "ArrowheadTail": "Full"},
    "resonance": {"ArrowheadHead": "Full", "ArrowheadTail": "Full"},
    "retrosynthetic": {"ArrowheadHead": "Full", "ArrowheadType": "Angle", "ArrowShaftSpacing": "500"},
    "retro": {"ArrowheadHead": "Full", "ArrowheadType": "Angle", "ArrowShaftSpacing": "500"},
    "dashed": {"ArrowheadHead": "Full", "LineType": "Dashed"},
    "bold": {"ArrowheadHead": "Full", "LineType": "Bold"},
    "hollow": {"ArrowheadHead": "Full", "ArrowheadType": "Hollow"},
    "nogo": {"ArrowheadHead": "Full", "NoGo": "Cross"},
    "cross": {"ArrowheadHead": "Full", "NoGo": "Cross"},
    "hash": {"ArrowheadHead": "Full", "NoGo": "Hash"},
    "line": {},
    "half_left": {"ArrowheadHead": "HalfLeft"},
    "half_right": {"ArrowheadHead": "HalfRight"},
    "unbalanced_equilibrium": {"ArrowheadHead": "HalfLeft", "ArrowheadTail": "HalfLeft", "ArrowShaftSpacing": "300",
                               "ArrowEquilibriumRatio": "0.4"},
}

FACE_BITS = {"plain": 0, "bold": 1, "italic": 2, "underline": 4, "outline": 8, "shadow": 16, "subscript": 32,
             "superscript": 64, "formula": 96}


# ---------------------------------------------------------------------------
# object references returned by the document builder
# ---------------------------------------------------------------------------

@dataclass
class ObjRef:
    id: int
    kind: str
    box: Box
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = {"id": self.id, "kind": self.kind, "bounds": self.box.to_dict()}
        d.update(self.extra)
        return d


# ---------------------------------------------------------------------------
# molecule -> fragment
# ---------------------------------------------------------------------------

def _mol_points(mol, bond_length: float) -> List[Tuple[float, float]]:
    """2D coordinates of all atoms scaled so that the median bond is ``bond_length`` points, y down."""
    if mol.GetNumConformers() == 0:
        from rdkit.Chem import rdDepictor
        rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    pts = [(conf.GetAtomPosition(i).x, -conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())]
    lens = []
    for b in mol.GetBonds():
        a, c = pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]
        lens.append(math.hypot(a[0] - c[0], a[1] - c[1]))
    median = statistics.median(lens) if lens else 1.5
    if median <= 1e-6:
        median = 1.5
    s = bond_length / median
    return [(x * s, y * s) for x, y in pts]


def _bond_display(bond, wavy_unspecified: bool = False) -> Optional[str]:
    """CDXML ``Display`` value for a bond: wedges from the RDKit bond direction, ``Wavy`` for unknown stereo.

    ``wavy_unspecified=True`` also draws a stereo-capable double bond whose geometry is unspecified as wavy, so the
    drawing does not silently read back as E (see ``CDXMLDocument.add_fragment``).
    """
    from rdkit import Chem
    d = bond.GetBondDir()
    if d == Chem.BondDir.BEGINWEDGE:
        return "WedgeBegin"
    if d == Chem.BondDir.BEGINDASH:
        return "WedgedHashBegin"
    if d == Chem.BondDir.UNKNOWN:
        return "Wavy"
    if bond.GetBondType() == Chem.BondType.DOUBLE:
        if bond.GetStereo() == Chem.BondStereo.STEREOANY:
            return "Wavy"
        if wavy_unspecified and bond.GetStereo() == Chem.BondStereo.STEREONONE:
            return "Wavy"
    return None


def _unspecified_double_bonds(mol) -> set:
    """Indices of stereo-capable double bonds with no E/Z assignment."""
    from rdkit import Chem
    out = set()
    try:
        for si in Chem.FindPotentialStereo(mol):
            if si.type == Chem.StereoType.Bond_Double and si.specified == Chem.StereoSpecified.Unspecified:
                out.add(si.centeredOn)
    except Exception:
        pass
    return out


def _bond_order(bond) -> Optional[str]:
    from rdkit import Chem
    bt = bond.GetBondType()
    if bt == Chem.BondType.SINGLE:
        return None
    if bt == Chem.BondType.DOUBLE:
        return "2"
    if bt == Chem.BondType.TRIPLE:
        return "3"
    if bt == Chem.BondType.AROMATIC:
        return "1.5"
    if bt == Chem.BondType.DATIVE:
        return "dative"
    if bt == Chem.BondType.QUADRUPLE:
        return "4"
    if bt == Chem.BondType.ZERO:
        return "ionic"
    return None


def _radical(atom) -> Optional[str]:
    """CDXML ``Radical`` value for an atom, or None.

    CDXML's enumeration (None, Singlet, Doublet, Triplet) can only express one unpaired electron (Doublet) or
    two (Singlet / Triplet, the spin state of a diradical); ChemDraw has no value for three or more, so those
    atoms are written without a Radical attribute (their NumHydrogens is still stated) rather than mislabelled.
    """
    n = atom.GetNumRadicalElectrons()
    if n == 1:
        return "Doublet"
    if n == 2:
        return "Singlet"
    return None


# ---------------------------------------------------------------------------
# document builder
# ---------------------------------------------------------------------------

class CDXMLDocument:
    """Incrementally build a CDXML document.

    Coordinates passed to ``add_*`` are absolute page coordinates in points.
    Call :meth:`to_xml` for the finished file.
    """

    def __init__(self, style: Optional[str] = None, name: str = "chemdrawmcp document",
                 bond_length: Optional[float] = None, settings: Optional[Dict[str, Any]] = None):
        self.style = get_style(style)
        if settings:
            self.style.update({k: v for k, v in settings.items() if v is not None})
        if bond_length:
            self.style["bond_length"] = float(bond_length)
        self.name = name
        self._next_id = 1
        self.objects: List[str] = []          # XML snippets inside <page>
        self.refs: List[ObjRef] = []
        self.steps: List[Dict[str, List[int]]] = []
        self.fonts: Dict[str, int] = {}       # font name -> id
        self.colors: List[Tuple[float, float, float]] = [tuple(c) for c in self.style.get("colortable", [])]
        self.font_id(self.style.get("label_font", "Arial"))
        self.font_id(self.style.get("caption_font", self.style.get("label_font", "Arial")))
        self._z = 1

    # -- bookkeeping ------------------------------------------------------
    @property
    def bond_length(self) -> float:
        return float(self.style["bond_length"])

    @property
    def label_size(self) -> float:
        return float(self.style.get("label_size", 10))

    @property
    def caption_size(self) -> float:
        return float(self.style.get("caption_size", self.label_size))

    def new_id(self) -> int:
        i = self._next_id
        self._next_id += 1
        return i

    def _next_z(self) -> int:
        self._z += 1
        return self._z

    def font_id(self, name: str) -> int:
        if name not in self.fonts:
            self.fonts[name] = 3 + len(self.fonts)
        return self.fonts[name]

    def color_index(self, color: Any) -> Optional[int]:
        """Resolve a colour spec (index, '#rrggbb', (r,g,b) 0-1 or 0-255, or a few names) to a CDXML colour index."""
        if color is None:
            return None
        if isinstance(color, int):
            return color
        names = {"black": (0, 0, 0), "white": (1, 1, 1), "red": (1, 0, 0), "green": (0, 0.6, 0), "blue": (0, 0, 1),
                 "orange": (1, 0.5, 0), "magenta": (1, 0, 1), "cyan": (0, 0.7, 0.9), "gray": (0.5, 0.5, 0.5),
                 "grey": (0.5, 0.5, 0.5), "purple": (0.5, 0, 0.5), "yellow": (0.9, 0.75, 0), "brown": (0.55, 0.27, 0.07)}
        if isinstance(color, str):
            c = color.strip().lower()
            if c == "black":
                return 0
            if c == "white":
                return 1
            if c in names:
                rgb = names[c]
            elif c.startswith("#") and len(c) == 7:
                rgb = tuple(int(c[i:i + 2], 16) / 255.0 for i in (1, 3, 5))
            else:
                raise ValueError("unknown colour %r" % color)
        else:
            rgb = tuple(float(v) for v in color)
            if max(rgb) > 1.0:
                rgb = tuple(v / 255.0 for v in rgb)
        rgb = tuple(round(v, 4) for v in rgb)
        for i, c in enumerate(self.colors):
            if all(abs(c[k] - rgb[k]) < 1e-3 for k in range(3)):
                return i + 2
        self.colors.append(rgb)
        return len(self.colors) + 1

    # -- fragments --------------------------------------------------------
    def add_fragment(self, mol, x: float, y: float, anchor: str = "center", color: Any = None,
                     atom_colors: Optional[Dict[int, Any]] = None, label_atoms: bool = False) -> ObjRef:
        """Add an RDKit molecule (2D conformer, Kekulé) centred (or anchored) at (x, y).

        anchor: "center" | "left" (x = left edge, y = vertical centre) | "topleft".
        Returns an ObjRef whose ``extra`` has ``atom_ids`` (node id per RDKit atom index) and ``smiles``.
        """
        from rdkit import Chem
        pts = _mol_points(mol, self.bond_length)
        if not pts:
            raise ValueError("molecule has no atoms")
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        raw = Box(min(xs), min(ys), max(xs), max(ys))
        pad = 0.35 * self.bond_length
        box0 = raw.padded(pad)
        # label-aware bounds: ChemDraw draws atom labels ("OH", "H3N+") outward from the atom position
        label_est = _tm.fragment_labels(mol, pts, self.label_size, label_atoms, self.style.get("label_font", "Arial"),
                                        float(self.style.get("margin_width", 1.6)),
                                        show_h=not self.style.get("hide_implicit_hydrogens"),
                                        show_terminal_carbons=bool(self.style.get("show_terminal_carbon_labels")))
        for lab in label_est:
            box0 = box0.union(Box(*lab["box"]))
        if anchor == "center":
            dx, dy = x - box0.cx, y - box0.cy
        elif anchor == "left":
            dx, dy = x - box0.x0, y - box0.cy
        elif anchor == "topleft":
            dx, dy = x - box0.x0, y - box0.y0
        else:
            raise ValueError("anchor must be center, left or topleft")
        pts = [(px + dx, py + dy) for px, py in pts]
        box = box0.translated(dx, dy)
        frag_id = self.new_id()
        atom_ids = [self.new_id() for _ in range(mol.GetNumAtoms())]
        frag_color = self.color_index(color)
        parts = ["<fragment id=\"%d\" BoundingBox=\"%s\" Z=\"%d\"%s>" % (
            frag_id, box.attr(), self._next_z(), (" color=\"%d\"" % frag_color) if frag_color is not None else "")]
        label_font = self.font_id(self.style.get("label_font", "Arial"))
        for atom in mol.GetAtoms():
            i = atom.GetIdx()
            attrs = ["id=\"%d\"" % atom_ids[i], "p=\"%s %s\"" % (_f(pts[i][0]), _f(pts[i][1]))]
            z = atom.GetAtomicNum()
            if z != 6 and z != 0:
                attrs.append("Element=\"%d\"" % z)
            show_label = bool(z != 6 or atom.GetDegree() == 0 or atom.GetFormalCharge() or atom.GetIsotope()
                              or atom.GetNumRadicalElectrons())
            if show_label or label_atoms:
                attrs.append("NumHydrogens=\"%d\"" % atom.GetTotalNumHs())
            if atom.GetFormalCharge():
                attrs.append("Charge=\"%d\"" % atom.GetFormalCharge())
            if atom.GetIsotope():
                attrs.append("Isotope=\"%d\"" % atom.GetIsotope())
            rad = _radical(atom)
            if rad:
                attrs.append("Radical=\"%s\"" % rad)
            if atom.GetAtomMapNum():
                attrs.append("AtomNumber=\"%d\"" % atom.GetAtomMapNum())
            if z == 0:
                # query / R-group atom: use its label if any
                lbl = atom.GetProp("dummyLabel") if atom.HasProp("dummyLabel") else (
                    atom.GetProp("_MolFileRLabel") if atom.HasProp("_MolFileRLabel") else "R")
                if not str(lbl).strip() or str(lbl).strip() == "*":
                    lbl = "R"      # RDKit labels bare "*" atoms "*"; ChemDraw shows R groups as R
                attrs.append("NodeType=\"GenericNickname\" GenericNickname=%s" % quoteattr(str(lbl)))
            acol = self.color_index(atom_colors.get(i)) if atom_colors and i in atom_colors else None
            if acol is not None:
                attrs.append("color=\"%d\"" % acol)
            if (label_atoms or z == 0) and z != 1:
                # explicit text label so ChemDraw shows e.g. "CH3" for carbons when asked
                txt = _atom_label_text(atom)
                parts.append("<n %s><t p=\"%s %s\" Justification=\"Center\"><s font=\"%d\" size=\"%s\" face=\"96\">%s</s></t></n>" % (
                    " ".join(attrs), _f(pts[i][0]), _f(pts[i][1] + self.label_size * 0.35), label_font, _f(self.label_size), escape(txt)))
            else:
                parts.append("<n %s/>" % " ".join(attrs))
        wavy_unspec = _unspecified_double_bonds(mol) if self.style.get("wavy_unspecified_double_bonds") else set()
        for bond in mol.GetBonds():
            attrs = ["id=\"%d\"" % self.new_id(), "B=\"%d\"" % atom_ids[bond.GetBeginAtomIdx()],
                     "E=\"%d\"" % atom_ids[bond.GetEndAtomIdx()]]
            order = _bond_order(bond)
            if order:
                attrs.append("Order=\"%s\"" % order)
            disp = _bond_display(bond, wavy_unspecified=bond.GetIdx() in wavy_unspec)
            if disp:
                attrs.append("Display=\"%s\"" % disp)
            parts.append("<b %s/>" % " ".join(attrs))
        parts.append("</fragment>")
        self.objects.append("".join(parts))
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = ""
        for lab in label_est:
            lab["box"] = [round(lab["box"][0] + dx, 2), round(lab["box"][1] + dy, 2), round(lab["box"][2] + dx, 2), round(lab["box"][3] + dy, 2)]
        ref = ObjRef(frag_id, "fragment", box, {"atom_ids": atom_ids, "smiles": smiles, "atoms": mol.GetNumAtoms(),
                                                "atom_positions": [[round(px, 2), round(py, 2)] for px, py in pts],
                                                "labels": label_est})
        self.refs.append(ref)
        return ref

    # -- text ---------------------------------------------------------------
    def add_text(self, text: str, x: float, y: float, size: Optional[float] = None, font: Optional[str] = None,
                 face: Any = "plain", color: Any = None, justification: str = "Center", anchor: str = "baseline",
                 formula: bool = False, line_height: Optional[float] = None, interpret: bool = True) -> ObjRef:
        """Add a caption.

        anchor="baseline": (x, y) is the baseline (x = left edge for Left, centre for Center, right edge for Right).
        anchor="top": y is the top of the text block; anchor="center": y is the vertical centre.
        face: "plain" | "bold" | "italic" | "formula" | int bit mask | list of names.  ``formula=True`` turns on
        ChemDraw's automatic sub/superscripts (face 96) — use it for reagents like "NaBH4", "H2SO4".
        """
        size = float(size or self.caption_size)
        font_name = font or self.style.get("caption_font", "Arial")
        fid = self.font_id(font_name)
        if isinstance(face, int):
            face_bits = face
        elif isinstance(face, (list, tuple)):
            face_bits = sum(FACE_BITS[f] for f in face)
        else:
            face_bits = FACE_BITS.get(str(face), 0)
        if formula:
            face_bits |= 96
        lines = text.split("\n")
        w, asc, desc = text_box(text, size, font_name, formula)
        if anchor == "top":
            y = y + asc
        elif anchor == "center":
            y = y + (asc - desc) / 2.0
        just = justification.capitalize()
        if just == "Center":
            x0 = x - w / 2.0
        elif just == "Right":
            x0 = x - w
        else:
            x0 = x
        box = Box(x0, y - asc, x0 + w, y + desc)
        tid = self.new_id()
        col = self.color_index(color)
        col_attr = (" color=\"%d\"" % col) if col is not None else ""
        lh = " LineHeight=\"%s\"" % _f(line_height) if line_height else " LineHeight=\"auto\""
        runs = "".join("<s font=\"%d\" size=\"%s\" face=\"%d\"%s>%s</s>" % (
            fid, _f(size), face_bits, col_attr, escape(ln + ("\n" if k < len(lines) - 1 else ""))) for k, ln in enumerate(lines))
        # ChemDraw chemically interprets captions; a caption that is exactly "+" makes it strip the formal
        # charges of every atom in the document (verified with ChemDraw 23), so plus signs opt out.
        interp = "" if interpret else " InterpretChemically=\"no\""
        xml = "<t id=\"%d\" p=\"%s %s\" BoundingBox=\"%s\" Z=\"%d\" Justification=\"%s\"%s%s%s>%s</t>" % (
            tid, _f(x), _f(y), box.attr(), self._next_z(), just, lh, col_attr, interp, runs)
        self.objects.append(xml)
        ref = ObjRef(tid, "text", box, {"text": text})
        self.refs.append(ref)
        return ref

    # -- arrows ---------------------------------------------------------------
    def add_arrow(self, x0: float, y0: float, x1: float, y1: float, kind: str = "full", color: Any = None,
                  line_width: Optional[float] = None, head_size: int = 1000) -> ObjRef:
        """Straight arrow from tail (x0, y0) to head (x1, y1).  kind: see ARROW_KINDS."""
        k = kind.lower().replace("-", "_").replace(" ", "_")
        if k not in ARROW_KINDS:
            raise ValueError("unknown arrow kind %r; known: %s" % (kind, ", ".join(sorted(ARROW_KINDS))))
        attrs = dict(ARROW_KINDS[k])
        aid = self.new_id()
        pad = 0.25 * self.bond_length
        box = Box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)).padded(pad)
        extra = ["%s=\"%s\"" % (a, v) for a, v in attrs.items()]
        col = self.color_index(color)
        if col is not None:
            extra.append("color=\"%d\"" % col)
        if line_width:
            extra.append("LineWidth=\"%s\"" % _f(line_width))
        xml = ("<arrow id=\"%d\" Z=\"%d\" BoundingBox=\"%s\" FillType=\"None\" ArrowheadType=\"Solid\" HeadSize=\"%d\" "
               "ArrowheadCenterSize=\"%d\" ArrowheadWidth=\"250\" Head3D=\"%s %s 0\" Tail3D=\"%s %s 0\" %s/>") % (
            aid, self._next_z(), box.attr(), head_size, int(head_size * 0.875), _f(x1), _f(y1), _f(x0), _f(y0), " ".join(extra))
        if "ArrowheadType" in attrs:   # kind overrides the default solid head
            xml = xml.replace("ArrowheadType=\"Solid\" ", "", 1)
        self.objects.append(xml)
        ref = ObjRef(aid, "arrow", box, {"kind": k, "tail": [round(x0, 2), round(y0, 2)], "head": [round(x1, 2), round(y1, 2)]})
        self.refs.append(ref)
        return ref

    def add_curved_arrow(self, x0: float, y0: float, x1: float, y1: float, angle: float = 90.0,
                         kind: str = "full", color: Any = None) -> ObjRef:
        """Circular-arc arrow (electron-pushing style) from tail to head sweeping ``angle`` degrees.

        Sign convention (matches ChemDraw's own files): a positive angle bulges to the RIGHT of the
        tail->head direction on screen (downwards for a left-to-right arrow; the sweep runs
        counter-clockwise on screen), a negative angle bulges to the LEFT (upwards for a
        left-to-right arrow).  |angle| > 180 draws the major arc on that same side.
        """
        k = kind.lower()
        attrs = dict(ARROW_KINDS.get(k, ARROW_KINDS["full"]))
        if "half" in k:
            attrs = {"ArrowheadHead": "HalfLeft"}
        chord = math.hypot(x1 - x0, y1 - y0)
        if chord < 1e-6:
            raise ValueError("curved arrow needs distinct tail and head")
        a = math.radians(max(1.0, min(359.0, abs(angle))))
        r = chord / (2.0 * math.sin(a / 2.0))
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        ux, uy = (x1 - x0) / chord, (y1 - y0) / chord
        # perpendicular (screen coordinates, y down): left of direction is (uy, -ux).  ChemDraw reads the sign of
        # AngularSize as the rotation sense from tail to head, so the centre of a minor arc (|angle| <= 180) lies on
        # the left for positive angles and the centre of a major arc on the right (and vice versa for negative).
        left = (angle > 0) == (abs(angle) <= 180.0)
        px, py = (uy, -ux) if left else (-uy, ux)
        h = math.sqrt(max(r * r - (chord / 2.0) ** 2, 0.0))
        cx, cy = mx + px * h, my + py * h
        aid = self.new_id()
        box = Box(min(x0, x1, cx - r), min(y0, y1, cy - r), max(x0, x1, cx + r), max(y0, y1, cy + r))
        extra = ["%s=\"%s\"" % (n, v) for n, v in attrs.items()]
        col = self.color_index(color)
        if col is not None:
            extra.append("color=\"%d\"" % col)
        xml = ("<arrow id=\"%d\" Z=\"%d\" BoundingBox=\"%s\" FillType=\"None\" ArrowheadType=\"Solid\" HeadSize=\"1000\" "
               "ArrowheadCenterSize=\"875\" ArrowheadWidth=\"250\" AngularSize=\"%s\" Head3D=\"%s %s 0\" Tail3D=\"%s %s 0\" "
               "Center3D=\"%s %s 0\" MajorAxisEnd3D=\"%s %s 0\" MinorAxisEnd3D=\"%s %s 0\" %s/>") % (
            aid, self._next_z(), box.attr(), _f(angle), _f(x1), _f(y1), _f(x0), _f(y0), _f(cx), _f(cy),
            _f(cx + r), _f(cy), _f(cx), _f(cy + r), " ".join(extra))
        self.objects.append(xml)
        ref = ObjRef(aid, "arrow", box, {"kind": "curved", "angle": angle})
        self.refs.append(ref)
        return ref

    # -- symbols / graphics ---------------------------------------------------
    def add_plus(self, x: float, y: float, size: Optional[float] = None) -> ObjRef:
        """A '+' between reaction components, centred at (x, y)."""
        size = size or self.caption_size + 2
        return self.add_text("+", x, y, size=size, anchor="center", justification="Center", interpret=False)

    def add_graphic(self, x0: float, y0: float, x1: float, y1: float, shape: str = "rectangle",
                    style: str = "Plain", color: Any = None, line_width: Optional[float] = None) -> ObjRef:
        """Rectangle / oval / line.  shape: rectangle | rounded | oval | circle | line; style: Plain | Dashed | Bold | Filled | Shaded."""
        gid = self.new_id()
        box = Box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
        shape = shape.lower()
        if shape in ("rectangle", "rect"):
            attrs = "GraphicType=\"Rectangle\" RectangleType=\"%s\"" % style
        elif shape in ("rounded", "roundrect"):
            attrs = "GraphicType=\"Rectangle\" RectangleType=\"%s\" CornerRadius=\"600\"" % (
                "RoundEdge" if style == "Plain" else "RoundEdge " + style)
        elif shape == "oval":
            attrs = "GraphicType=\"Oval\" OvalType=\"%s\"" % ("" if style == "Plain" else style)
        elif shape == "circle":
            attrs = "GraphicType=\"Oval\" OvalType=\"Circle%s\"" % ("" if style == "Plain" else " " + style)
        elif shape == "line":
            attrs = "GraphicType=\"Line\" LineType=\"%s\"" % ("Solid" if style == "Plain" else style)
            box = Box(x0, y0, x1, y1)   # lines keep their direction
        else:
            raise ValueError("unknown shape %r" % shape)
        col = self.color_index(color)
        extra = (" color=\"%d\"" % col) if col is not None else ""
        if line_width:
            extra += " LineWidth=\"%s\"" % _f(line_width)
        self.objects.append("<graphic id=\"%d\" Z=\"%d\" BoundingBox=\"%s\" %s%s/>" % (
            gid, self._next_z(), box.attr(), attrs, extra))
        ref = ObjRef(gid, "graphic", Box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)), {"shape": shape})
        self.refs.append(ref)
        return ref

    def add_raw(self, xml: str, box: Optional[Box] = None) -> None:
        """Append a raw XML snippet at page level (advanced)."""
        self.objects.append(xml)
        if box:
            self.refs.append(ObjRef(self.new_id(), "raw", box))

    # -- reaction semantics ---------------------------------------------------
    def add_step(self, reactants: Sequence[ObjRef], products: Sequence[ObjRef], arrows: Sequence[ObjRef] = (),
                 above: Sequence[ObjRef] = (), below: Sequence[ObjRef] = (), plusses: Sequence[ObjRef] = (),
                 atom_map: Optional[Sequence[Tuple[int, int]]] = None) -> None:
        """Record a reaction step so ChemDraw (and readers) know what reacts to what."""
        self.steps.append({
            "ReactionStepReactants": [r.id for r in reactants],
            "ReactionStepProducts": [p.id for p in products],
            "ReactionStepArrows": [a.id for a in arrows],
            "ReactionStepObjectsAboveArrow": [t.id for t in above],
            "ReactionStepObjectsBelowArrow": [t.id for t in below],
            "ReactionStepPlusses": [p.id for p in plusses],
            "ReactionStepAtomMap": [i for pair in (atom_map or []) for i in pair],
        })

    # -- output -----------------------------------------------------------------
    def content_box(self) -> Optional[Box]:
        return union_boxes(r.box for r in self.refs)

    def to_xml(self, fit_page: bool = True, margin: float = 36.0) -> str:
        st = self.style
        page = list(st.get("page", [0, 0, 540, 720]))
        page_w, page_h = max(page[2] - page[0], 72.0), max(page[3] - page[1], 72.0)
        pages_w = pages_h = 1
        content = self.content_box()
        if fit_page and content is not None:
            # ChemDraw asks "expand the document?" when objects fall outside the drawing space: tile printer pages
            pages_w = max(1, int(math.ceil((content.x1 + margin) / page_w)))
            pages_h = max(1, int(math.ceil((content.y1 + margin) / page_h)))
        page = [0, 0, page_w * pages_w, page_h * pages_h]
        label_font = self.font_id(st.get("label_font", "Arial"))
        caption_font = self.font_id(st.get("caption_font", st.get("label_font", "Arial")))
        head = [
            '<?xml version="1.0" encoding="UTF-8" ?>',
            '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
            '<CDXML CreationProgram="chemdrawmcp" Name=%s' % quoteattr(self.name),
            ' BoundingBox="%s"' % (content.attr() if content else "0 0 0 0"),
            ' WindowPosition="0 0" WindowSize="0 0"',
            ' FractionalWidths="%s"' % ("yes" if st.get("fractional_widths", True) else "no"),
            ' InterpretChemically="yes" ShowAtomQuery="yes" ShowAtomStereo="no" ShowAtomEnhancedStereo="yes"',
            ' ShowAtomNumber="no" ShowResidueID="no" ShowBondQuery="yes" ShowBondRxn="yes" ShowBondStereo="no"',
            ' ShowTerminalCarbonLabels="%s" ShowNonTerminalCarbonLabels="no" HideImplicitHydrogens="%s"' % (
                "yes" if st.get("show_terminal_carbon_labels") else "no", "yes" if st.get("hide_implicit_hydrogens") else "no"),
            ' LabelFont="%d" LabelSize="%s" LabelFace="%d"' % (label_font, _f(self.label_size), int(st.get("label_face", 96))),
            ' CaptionFont="%d" CaptionSize="%s" CaptionFace="%d"' % (caption_font, _f(self.caption_size), int(st.get("caption_face", 0))),
            ' HashSpacing="%s" MarginWidth="%s" LineWidth="%s" BoldWidth="%s" BondLength="%s" BondSpacing="%s" ChainAngle="%s"' % (
                _f(st.get("hash_spacing", 2.5)), _f(st.get("margin_width", 1.6)), _f(st.get("line_width", 0.6)),
                _f(st.get("bold_width", 2)), _f(self.bond_length), _f(st.get("bond_spacing", 18)), _f(st.get("chain_angle", 120))),
            ' LabelJustification="%s" CaptionJustification="%s"' % (st.get("label_justification", "Auto"), st.get("caption_justification", "Left")),
            ' AminoAcidTermini="HOH" ShowSequenceTermini="yes" ShowSequenceBonds="yes" ResidueWrapCount="40" ResidueBlockCount="10"',
            ' PrintMargins="%s"' % " ".join(_f(v) for v in st.get("print_margins", [36, 36, 36, 36])),
            ' color="%d" bgcolor="%d"' % (int(st.get("color", 0)), int(st.get("bgcolor", 1))),
            ' RxnAutonumberStart="1" RxnAutonumberConditions="no" RxnAutonumberStyle="Roman" RxnAutonumberFormat="(#)">',
        ]
        colortable = ["<colortable>"] + ["<color r=\"%s\" g=\"%s\" b=\"%s\"/>" % (_f(c[0]), _f(c[1]), _f(c[2])) for c in self.colors] + ["</colortable>"]
        fonttable = ["<fonttable>"] + ["<font id=\"%d\" charset=\"iso-8859-1\" name=%s/>" % (fid, quoteattr(name))
                                       for name, fid in sorted(self.fonts.items(), key=lambda kv: kv[1])] + ["</fonttable>"]
        page_xml = ['<page id="%d" BoundingBox="%s" Width="%s" Height="%s" HeaderPosition="36" FooterPosition="36" PageOverlap="0" '
                    'PrintTrimMarks="yes" HeightPages="%d" WidthPages="%d" DrawingSpace="poster">' % (
                        self.new_id(), " ".join(_f(v) for v in page), _f(page[2]), _f(page[3]), pages_h, pages_w)]
        body = list(self.objects)
        if self.steps:
            steps = []
            for s in self.steps:
                attrs = " ".join("%s=\"%s\"" % (k, " ".join(str(i) for i in v)) for k, v in s.items() if v)
                steps.append("<step id=\"%d\" %s/>" % (self.new_id(), attrs))
            body.append("<scheme id=\"%d\">%s</scheme>" % (self.new_id(), "".join(steps)))
        return "\n".join(head + colortable + fonttable + page_xml + body + ["</page>", "</CDXML>", ""])

    def summary(self) -> dict:
        content = self.content_box()
        return {"style": self.style.get("name"), "bond_length": self.bond_length, "objects": [r.to_dict() for r in self.refs],
                "steps": len(self.steps), "bounds": content.to_dict() if content else None}


def _atom_label_text(atom) -> str:
    sym = atom.GetSymbol() if atom.GetAtomicNum() else (atom.GetProp("dummyLabel") if atom.HasProp("dummyLabel") else "R")
    nh = atom.GetTotalNumHs()
    txt = sym
    if nh == 1:
        txt += "H"
    elif nh > 1:
        txt += "H%d" % nh
    chg = atom.GetFormalCharge()
    if chg:
        # ChemDraw's formula face renders a trailing sign as a superscript charge ("OH-", "NH3+", "Fe3+")
        txt += ("" if abs(chg) == 1 else str(abs(chg))) + ("+" if chg > 0 else "-")
    return txt


# ---------------------------------------------------------------------------
# reader
# ---------------------------------------------------------------------------

def _strip_doctype(text: str) -> str:
    if isinstance(text, bytes):
        text = text.decode("utf-8", "replace")
    text = text.replace("\r", "").lstrip("\ufeff \n\t")
    return re.sub(r"<!DOCTYPE[^>]*>", "", text)


def parse_cdxml(text: str) -> ET.Element:
    """Parse CDXML text (str or UTF-8 bytes) into an ElementTree root (DOCTYPE and BOM removed)."""
    return ET.fromstring(_strip_doctype(text))


def _text_of(t: ET.Element) -> str:
    return "".join((s.text or "") for s in t.findall("s")) or (t.text or "")


def _box_of(el: ET.Element) -> Optional[Box]:
    bb = el.get("BoundingBox")
    if not bb:
        return None
    try:
        v = [float(x) for x in bb.split()]
        return Box(min(v[0], v[2]), min(v[1], v[3]), max(v[0], v[2]), max(v[1], v[3]))
    except (ValueError, IndexError):
        return None


def fragment_to_mol(frag: ET.Element, sanitize: bool = True):
    """Convert one ``<fragment>`` element to an RDKit molecule (2D conformer, wedge- and coordinate-derived stereo).

    * Nested fragments (``NodeType="Nickname"``/``"Fragment"`` nodes with an inner ``<fragment>``) are expanded in
      place.  Their ``ExternalConnectionPoint`` nodes are matched to the outer bonds through the inner fragment's
      ``ConnectionOrder`` and the node's ``BondOrdering`` (peptide residues have two or three attachment points).
    * Nicknames without an inner fragment are expanded from :mod:`chemdrawmcp.abbreviations` when the label is
      known ("OMe", "Boc", "Ph", "NBS", ...); otherwise (R groups, element lists, uninterpretable labels) they become
      dummy atoms carrying a ``dummyLabel`` property.
    * Charge, isotope, hydrogen count and radical come from the attributes, falling back to the label text
      ("O-", "NH3+", "D", "13C").  Explicit hydrogen nodes are folded into their heavy atoms after stereo perception.
    * ``WedgeBegin``/``WedgedHashBegin`` (and the ``*End`` forms) give tetrahedral stereo, E/Z comes from the
      coordinates, ``Wavy`` marks unknown stereo, ``Order="1.5"`` rings are kekulized when possible and non-ring
      "aromatic" bonds are downgraded to single bonds.

    Chemistry problems that could not be repaired (valence errors, unkekulizable rings) are stored on the molecule as
    the ``_cdxml_errors`` property (one per line) and a best-effort molecule is still returned; softer notes go to
    ``_cdxml_warnings`` and uninterpreted labels to ``_cdxml_labels``.  Returns (mol, node_id_to_atom_index).
    """
    from rdkit import Chem, rdBase
    from rdkit.Geometry import Point3D
    from . import abbreviations as _ab

    rw = Chem.RWMol()
    coords: List[Tuple[float, float]] = []
    id_to_idx: Dict[str, int] = {}
    attach: Dict[str, List[int]] = {}                  # nested node id -> attachment atoms, one per connection point
    attach_bonds: Dict[str, List[str]] = {}            # nested node id -> BondOrdering (bond ids matching the above)
    attach_extra: Dict[str, List[Tuple[Optional[str], Optional[str]]]] = {}   # per connection point: (Order, Display)
    attach_used: Dict[str, set] = {}
    warnings: List[str] = []
    errors: List[str] = []
    labels: List[str] = []
    wavy_double: List[Tuple[int, int]] = []
    dup_bonds = 0
    DUMMY_TYPES = ("GenericNickname", "ElementList", "MultiAttachment", "VariableAttachment", "LinkNode", "Monomer",
                   "AnonymousAlternativeGroup", "NamedAlternativeGroup", "ExternalConnectionPoint")
    LABEL_TYPES = ("Nickname", "Fragment", "Unspecified", "Formula")

    def node_pos(node: ET.Element) -> Tuple[float, float]:
        p = (node.get("p") or "0 0").split()
        try:
            return (float(p[0]), -float(p[1])) if len(p) >= 2 else (0.0, 0.0)
        except ValueError:
            return (0.0, 0.0)

    def node_label(node: ET.Element) -> str:
        t = node.find("t")
        return _text_of(t).strip() if t is not None else ""

    def bonded_ids(f: ET.Element) -> set:
        s = set()
        for b in f.findall("b"):
            s.add(b.get("B", ""))
            s.add(b.get("E", ""))
        return s

    unit_cache: Dict[int, float] = {}

    def unit_length(f: ET.Element) -> float:
        """Median bond length (points) of a fragment's own nodes, for placing expanded abbreviations."""
        if id(f) in unit_cache:
            return unit_cache[id(f)]
        pos = {n.get("id"): node_pos(n) for n in f.findall("n")}
        lens = []
        for b in f.findall("b"):
            a, c = pos.get(b.get("B")), pos.get(b.get("E"))
            if a and c:
                lens.append(math.hypot(a[0] - c[0], a[1] - c[1]))
        lens = [v for v in lens if v > 1e-6]
        unit_cache[id(f)] = statistics.median(lens) if lens else 14.4
        return unit_cache[id(f)]

    def new_atom(z: int, x: float, y: float) -> int:
        idx = rw.AddAtom(Chem.Atom(z))
        coords.append((x, y))
        return idx

    def apply_attrs(idx: int, node: ET.Element, label: str) -> None:
        a = rw.GetAtomWithIdx(idx)
        parsed = _ab.parse_atom_label(label) if label else None
        if parsed and parsed.get("symbol") != a.GetSymbol():
            parsed = None
        charge = node.get("Charge")
        if charge:
            try:
                a.SetFormalCharge(int(float(charge)))
            except ValueError:
                pass
        elif parsed and parsed.get("charge"):
            a.SetFormalCharge(parsed["charge"])
        iso = node.get("Isotope")
        if iso:
            try:
                a.SetIsotope(int(float(iso)))
            except ValueError:
                pass
        elif parsed and parsed.get("isotope"):
            a.SetIsotope(parsed["isotope"])
        nh = node.get("NumHydrogens")
        if a.GetAtomicNum() != 6:
            # ChemDraw states the hydrogen count of labelled atoms; for carbon RDKit's implicit model is more reliable
            hcount = None
            if nh is not None:
                try:
                    hcount = int(float(nh))
                except ValueError:
                    hcount = None
            elif parsed and "hcount" in parsed:
                hcount = parsed["hcount"]
            if hcount is not None and hcount >= 0:
                a.SetNumExplicitHs(hcount)
                a.SetNoImplicit(True)
        rad = node.get("Radical")
        if rad:
            a.SetNumRadicalElectrons({"Doublet": 1, "Singlet": 2, "Triplet": 2}.get(rad, 0))
        elif parsed and parsed.get("radical"):
            a.SetNumRadicalElectrons(parsed["radical"])
        if node.get("AtomNumber"):
            try:
                a.SetAtomMapNum(int(float(node.get("AtomNumber"))))
            except ValueError:
                pass

    def add_dummy(node: ET.Element, label: str, x: float, y: float) -> int:
        idx = new_atom(0, x, y)
        a = rw.GetAtomWithIdx(idx)
        if label:
            a.SetProp("dummyLabel", label)
            labels.append(label)
        a.SetNoImplicit(True)
        if node.get("Charge"):
            try:
                a.SetFormalCharge(int(float(node.get("Charge"))))
            except ValueError:
                pass
        return idx

    def add_smiles(smi: str, x: float, y: float, unit: float) -> Optional[int]:
        """Add a SMILES expansion of an abbreviation with its first atom at (x, y); returns that atom's index."""
        with rdBase.BlockLogs():
            m = Chem.MolFromSmiles(smi)
        if m is None or m.GetNumAtoms() == 0:
            return None
        from rdkit.Chem import rdDepictor
        with rdBase.BlockLogs():
            rdDepictor.Compute2DCoords(m)
        c = m.GetConformer()
        p0 = c.GetAtomPosition(0)
        s = unit / 1.5
        base = rw.GetNumAtoms()
        for a in m.GetAtoms():
            pos = c.GetAtomPosition(a.GetIdx())
            idx = new_atom(a.GetAtomicNum(), x + (pos.x - p0.x) * s, y + (pos.y - p0.y) * s)
            na = rw.GetAtomWithIdx(idx)
            na.SetFormalCharge(a.GetFormalCharge())
            na.SetIsotope(a.GetIsotope())
            na.SetIsAromatic(a.GetIsAromatic())
            na.SetNumRadicalElectrons(a.GetNumRadicalElectrons())
            if a.GetNoImplicit() or a.GetNumExplicitHs():
                na.SetNumExplicitHs(a.GetNumExplicitHs())
                na.SetNoImplicit(True)
        for b in m.GetBonds():
            rw.AddBond(base + b.GetBeginAtomIdx(), base + b.GetEndAtomIdx(), b.GetBondType())
            nb = rw.GetBondBetweenAtoms(base + b.GetBeginAtomIdx(), base + b.GetEndAtomIdx())
            nb.SetIsAromatic(b.GetIsAromatic())
        return base

    def resolve(nid: str, bond_id: str) -> Optional[int]:
        """Atom index a bond end refers to; nested nodes pick the connection point matching the bond."""
        if nid in attach:
            lst = attach[nid]
            if not lst:
                return None
            if len(lst) == 1:
                return lst[0]
            used = attach_used.setdefault(nid, set())
            order = attach_bonds.get(nid) or []
            if bond_id in order:
                k = order.index(bond_id)
                if k < len(lst):
                    used.add(k)
                    return lst[k]
            for k, a in enumerate(lst):
                if k not in used:
                    used.add(k)
                    return a
            return lst[0]
        return id_to_idx.get(nid)

    def attach_hint(nid: str, bond_id: str) -> Tuple[Optional[str], Optional[str]]:
        """(Order, Display) of the inner bond behind a connection point, for the outer bond that replaces it."""
        if nid not in attach_extra:
            return None, None
        lst = attach[nid]
        if len(lst) == 1:
            return attach_extra[nid][0]
        order = attach_bonds.get(nid) or []
        if bond_id in order and order.index(bond_id) < len(lst):
            return attach_extra[nid][order.index(bond_id)]
        return None, None

    BOND_TYPES = {"1": Chem.BondType.SINGLE, "2": Chem.BondType.DOUBLE, "3": Chem.BondType.TRIPLE,
                  "1.5": Chem.BondType.AROMATIC, "dative": Chem.BondType.DATIVE, "4": Chem.BondType.QUADRUPLE,
                  "ionic": Chem.BondType.ZERO, "hydrogen": Chem.BondType.ZERO, "0.5": Chem.BondType.SINGLE,
                  "2.5": Chem.BondType.DOUBLE, "3.5": Chem.BondType.TRIPLE, "any": Chem.BondType.SINGLE,
                  "singleordouble": Chem.BondType.SINGLE, "singleoraromatic": Chem.BondType.SINGLE,
                  "doubleoraromatic": Chem.BondType.DOUBLE}

    def add_bond(b: ET.Element, order: Optional[str] = None, display: Optional[str] = None) -> None:
        nonlocal dup_bonds
        bid = b.get("id", "")
        bn, en = b.get("B", ""), b.get("E", "")
        bi, ei = resolve(bn, bid), resolve(en, bid)
        if bi is None or ei is None or bi == ei:
            return
        ho_b, hd_b = attach_hint(bn, bid)
        ho_e, hd_e = attach_hint(en, bid)
        order = b.get("Order") or order or ho_b or ho_e
        display = b.get("Display") or display or hd_b or (hd_e and _flip_wedge(hd_e)) or None
        key = (order or "1").lower()
        bt = BOND_TYPES.get(key)
        if bt is None:
            warnings.append("bond %s: unknown order %r treated as single" % (bid, order))
            bt = Chem.BondType.SINGLE
        if key in ("0.5", "2.5", "3.5"):
            warnings.append("bond %s: partial order %s rounded" % (bid, key))
        if rw.GetBondBetweenAtoms(bi, ei) is not None:
            dup_bonds += 1
            return
        disp = display or ""
        if disp in ("WedgeEnd", "WedgedHashEnd", "HollowWedgeEnd"):
            bi, ei = ei, bi
            disp = {"WedgeEnd": "WedgeBegin", "WedgedHashEnd": "WedgedHashBegin", "HollowWedgeEnd": "HollowWedgeBegin"}[disp]
        rw.AddBond(bi, ei, bt)
        bond = rw.GetBondBetweenAtoms(bi, ei)
        if bt == Chem.BondType.AROMATIC:
            bond.SetIsAromatic(True)
        if disp in ("WedgeBegin", "HollowWedgeBegin"):
            bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
        elif disp == "WedgedHashBegin":
            bond.SetBondDir(Chem.BondDir.BEGINDASH)
        elif disp in ("Wavy", "WavyWedgeBegin", "WavyWedgeEnd"):
            if bt == Chem.BondType.DOUBLE:
                wavy_double.append((bi, ei))
            bond.SetBondDir(Chem.BondDir.UNKNOWN)

    def add_nested(node: ET.Element, sub: ET.Element, x: float, y: float) -> None:
        """Expand an inner fragment; record attachment atoms for the node."""
        nid = node.get("id", "")
        start = rw.GetNumAtoms()
        ecps: Dict[str, Tuple[float, float]] = {}
        ecp_num: Dict[str, float] = {}
        for sn in sub.findall("n"):
            if sn.get("NodeType") == "ExternalConnectionPoint":
                ecps[sn.get("id", "")] = node_pos(sn)
                try:
                    ecp_num[sn.get("id", "")] = float(sn.get("ExternalConnectionNum") or 0)
                except ValueError:
                    pass
                continue
            add_node(sn, sub)
        ecp_attach: Dict[str, Tuple[int, Optional[str], Optional[str]]] = {}
        for b in sub.findall("b"):
            bn, en = b.get("B", ""), b.get("E", "")
            if bn in ecps or en in ecps:
                ecp, other = (bn, en) if bn in ecps else (en, bn)
                oidx = resolve(other, b.get("id", ""))
                if oidx is None:
                    continue
                disp = b.get("Display")
                if disp and bn == ecp:
                    disp = _flip_wedge(disp)
                ecp_attach[ecp] = (oidx, b.get("Order"), disp)
                continue
            add_bond(b)
        order = [e for e in (sub.get("ConnectionOrder") or "").split() if e in ecps]
        rest = [e for e in ecps if e not in order]
        rest.sort(key=lambda e: ecp_num.get(e, 0))
        order += rest
        atoms: List[int] = []
        extra: List[Tuple[Optional[str], Optional[str]]] = []
        for e in order:
            if e in ecp_attach:
                atoms.append(ecp_attach[e][0])
                extra.append((ecp_attach[e][1], ecp_attach[e][2]))
        inner_count = rw.GetNumAtoms() - start
        if inner_count == 0:
            attach[nid] = []
            return
        if not atoms:
            # no connection point: anchor on the inner atom closest to the node position, else the first atom
            best, bestd = start, None
            for i in range(start, rw.GetNumAtoms()):
                d = math.hypot(coords[i][0] - x, coords[i][1] - y)
                if bestd is None or d < bestd:
                    best, bestd = i, d
            atoms = [best]
            extra = [(None, None)]
        attach[nid] = atoms
        attach_extra[nid] = extra
        attach_bonds[nid] = (node.get("BondOrdering") or "").split()
        # move the expansion so its first attachment atom sits at the node position (ChemDraw keeps nickname
        # definitions in their own coordinate frame)
        ax, ay = coords[atoms[0]]
        dx, dy = x - ax, y - ay
        if abs(dx) > 1e-6 or abs(dy) > 1e-6:
            for i in range(start, rw.GetNumAtoms()):
                coords[i] = (coords[i][0] + dx, coords[i][1] + dy)
        id_to_idx[nid] = atoms[0]

    def add_node(node: ET.Element, parent: ET.Element) -> None:
        nid = node.get("id", "")
        x, y = node_pos(node)
        ntype = node.get("NodeType", "Element")
        label = node_label(node)
        sub = node.find("fragment")
        has_inner = sub is not None and bool(sub.findall("n"))
        if ntype in ("Fragment", "Nickname") and has_inner:
            add_nested(node, sub, x, y)
            return
        if ntype == "ExternalConnectionPoint":
            # a connection point of a nickname definition at top level: attachment marker
            id_to_idx[nid] = add_dummy(node, "", x, y)
            return
        if ntype in LABEL_TYPES or ntype in DUMMY_TYPES:
            text = label or node.get("GenericNickname", "")
            if ntype in LABEL_TYPES and not (ntype == "Unspecified" and len(text.strip()) == 1):
                hit = _ab.lookup(text)
                if hit is not None:
                    smi, is_sub = hit
                    idx = add_smiles(smi, x, y, unit_length(parent))
                    if idx is not None:
                        id_to_idx[nid] = idx
                        if node.get("Charge"):
                            apply_attrs(idx, node, "")
                        if is_sub and nid not in bonded_ids(parent):
                            # a bare substituent label ("Ph", "Boc"): keep its open valence as an attachment point
                            d = new_atom(0, x - unit_length(parent), y)
                            rw.GetAtomWithIdx(d).SetNoImplicit(True)
                            rw.AddBond(d, idx, Chem.BondType.SINGLE)
                        return
            if ntype == "ElementList" and node.get("ElementList"):
                try:
                    text = "[" + ",".join(SYMBOLS[int(v)] for v in node.get("ElementList").split()) + "]"
                except (ValueError, IndexError):
                    pass
            if not text:
                text = {"GenericNickname": "R", "MultiAttachment": "MultiAttachment", "VariableAttachment": "VariableAttachment",
                        "LinkNode": "LinkNode", "Monomer": "Monomer"}.get(ntype, ntype)
            if ntype == "MultiAttachment":
                warnings.append("node %s: multi-centre (haptic) attachment represented as a dummy atom" % nid)
            id_to_idx[nid] = add_dummy(node, text, x, y)
            return
        try:
            z = int(float(node.get("Element", "6")))
        except ValueError:
            z = 6
        if z == 18 and label == "Ar" and nid in bonded_ids(parent):
            # "Ar" on a bonded atom is an aryl group, not argon
            id_to_idx[nid] = add_dummy(node, "Ar", x, y)
            return
        if z <= 0 or z >= len(SYMBOLS):
            id_to_idx[nid] = add_dummy(node, label or node.get("GenericNickname", "") or "*", x, y)
            return
        idx = new_atom(z, x, y)
        apply_attrs(idx, node, label)
        id_to_idx[nid] = idx

    for node in frag.findall("n"):
        add_node(node, frag)
    for b in frag.findall("b"):
        add_bond(b)
    if dup_bonds:
        warnings.append("%d duplicate bond(s) ignored" % dup_bonds)

    mol = rw.GetMol()
    n = mol.GetNumAtoms()
    conf = Chem.Conformer(n)
    lens = [math.hypot(coords[b.GetBeginAtomIdx()][0] - coords[b.GetEndAtomIdx()][0],
                       coords[b.GetBeginAtomIdx()][1] - coords[b.GetEndAtomIdx()][1]) for b in mol.GetBonds()]
    lens = [v for v in lens if v > 1e-6]
    scale = 1.5 / statistics.median(lens) if lens else 1.5 / 14.4
    for i in range(n):
        conf.SetAtomPosition(i, Point3D(coords[i][0] * scale, coords[i][1] * scale, 0.0))
    conf.Set3D(False)
    mol.AddConformer(conf, assignId=True)

    if sanitize and n:
        with rdBase.BlockLogs():
            # aromatic (order 1.5) bonds: flag ring atoms aromatic so kekulization can run; chains cannot be aromatic
            arom = [b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.AROMATIC]
            if arom:
                Chem.FastFindRings(mol)
                ri = mol.GetRingInfo()
                for b in arom:
                    if ri.NumBondRings(b.GetIdx()) == 0:
                        b.SetBondType(Chem.BondType.SINGLE)
                        b.SetIsAromatic(False)
                        warnings.append("bond %d-%d drawn with order 1.5 outside a ring treated as single" % (
                            b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
                    else:
                        b.GetBeginAtom().SetIsAromatic(True)
                        b.GetEndAtom().SetIsAromatic(True)
                # RWMol.AddBond flags both atoms of an aromatic bond; clear atoms that lost all their aromatic bonds
                for a in mol.GetAtoms():
                    if a.GetIsAromatic() and not any(b.GetIsAromatic() for b in a.GetBonds()):
                        a.SetIsAromatic(False)
            ok = False
            try:
                Chem.SanitizeMol(mol)
                ok = True
            except Exception as e:
                mol.UpdatePropertyCache(strict=False)
                probs = []
                try:
                    probs = [p.Message() for p in Chem.DetectChemistryProblems(mol)]
                except Exception:
                    pass
                msg = str(e).strip().split("\n")[0]
                errors.extend(probs or [msg or "sanitization failed"])
                # best effort: skip only the failing step (valence check or kekulization) so a SMILES can still be written
                no_valence = Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES
                no_kekule = Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY
                both = no_valence ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY
                kek = "kekulize" in (msg + " ".join(probs)).lower()
                for ops in ((no_kekule, no_valence, both) if kek else (no_valence, no_kekule, both)):
                    try:
                        mol.UpdatePropertyCache(strict=False)
                        Chem.SanitizeMol(mol, sanitizeOps=ops)
                        break
                    except Exception:
                        continue
                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass
            try:
                Chem.AssignChiralTypesFromBondDirs(mol)
                Chem.DetectBondStereochemistry(mol)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            except Exception as e:
                warnings.append("stereo perception failed: %s" % e)
            for bi, ei in wavy_double:
                b = mol.GetBondBetweenAtoms(bi, ei)
                if b is not None and b.GetBondType() == Chem.BondType.DOUBLE:
                    b.SetStereo(Chem.BondStereo.STEREOANY)
            # fold drawn hydrogens into their heavy atoms (isotopes, charged and lone H atoms are kept)
            if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
                for a in mol.GetAtoms():
                    a.SetIntProp("_cdxml_idx", a.GetIdx())
                try:
                    heavy_of = {a.GetIdx(): (a.GetNeighbors()[0].GetIdx() if a.GetDegree() else None)
                                for a in mol.GetAtoms() if a.GetAtomicNum() == 1}
                    stripped = Chem.RemoveHs(mol, sanitize=ok)
                    remap = {a.GetIntProp("_cdxml_idx"): a.GetIdx() for a in stripped.GetAtoms()}
                    for k, v in list(id_to_idx.items()):
                        if v in remap:
                            id_to_idx[k] = remap[v]
                        else:
                            h = heavy_of.get(v)
                            id_to_idx[k] = remap[h] if h in remap else remap.get(v, 0)
                    mol = stripped
                except Exception:
                    pass
                for a in mol.GetAtoms():
                    if a.HasProp("_cdxml_idx"):
                        a.ClearProp("_cdxml_idx")
    if errors:
        mol.SetProp("_cdxml_errors", "\n".join(errors))
    if warnings:
        mol.SetProp("_cdxml_warnings", "\n".join(warnings))
    if labels:
        mol.SetProp("_cdxml_labels", "\n".join(labels))
    return mol, id_to_idx


def _flip_wedge(display: str) -> str:
    """Swap the Begin/End sense of a wedge display name (used when a bond is traversed backwards)."""
    if display.endswith("Begin"):
        return display[:-5] + "End"
    if display.endswith("End"):
        return display[:-3] + "Begin"
    return display


def read_cdxml(text: str, want_mols: bool = True) -> Dict[str, Any]:
    """Structured description of a CDXML document.

    Returns {settings, pages, groups, molecules:[{id, smiles, formula, molecular_weight, heavy_atoms, atoms, bonds,
    bounds, name?, properties?, labels?, text?, warnings?, error?}], texts:[{id, text, bounds, chemical_property?}],
    arrows:[...], graphics:[...], steps:[{reactants, products, arrows, above_arrow, below_arrow, atom_map_pairs}],
    reaction_smiles?, bracketed_groups?, spectra?, other_objects?, warnings:[...]}.

    * ``pages`` counts top-level pages (pages nested in tables are not counted); ``groups`` counts ``<group>`` elements.
    * ``labels`` lists labels that became dummy atoms (R groups, uninterpretable text); a fragment that is only a
      caption typed as an atom label has ``smiles == "*"`` and its ``text``.
    * ``name`` / ``properties`` come from ``<chemicalproperty>`` records (which usually follow the ``<t>`` that
      displays them as siblings); the text also gets the record as ``chemical_property``.
    * ``error`` on a molecule means RDKit could not fully sanitize the drawing (invalid valence, unkekulizable ring);
      the SMILES given with it is a best-effort interpretation.  Reaction steps that reference a ``<group>`` are
      expanded to the group's members.
    """
    root = parse_cdxml(text)
    top_pages = root.findall("page")
    out: Dict[str, Any] = {"settings": {}, "molecules": [], "texts": [], "arrows": [], "graphics": [], "steps": [],
                           "warnings": [], "pages": len(top_pages) if top_pages else len(root.findall(".//page")),
                           "groups": len(root.findall(".//group"))}
    for k in ("BondLength", "LineWidth", "BoldWidth", "LabelSize", "CaptionSize", "BondSpacing", "ChainAngle",
              "HashSpacing", "MarginWidth", "CreationProgram", "Name"):
        if root.get(k) is not None:
            v = root.get(k)
            try:
                out["settings"][k] = float(v) if re.match(r"^-?[\d.]+$", v) else v
            except ValueError:
                out["settings"][k] = v
    fonts = {f.get("id"): f.get("name") for f in root.findall(".//fonttable/font")}
    lf, cf = root.get("LabelFont"), root.get("CaptionFont")
    if lf in fonts:
        out["settings"]["LabelFont"] = fonts[lf]
    if cf in fonts:
        out["settings"]["CaptionFont"] = fonts[cf]

    by_id: Dict[str, Dict[str, Any]] = {}
    node_owner: Dict[str, Dict[str, Any]] = {}     # atom/nested-fragment id -> molecule entry
    other: Dict[str, int] = {}
    groups: List[ET.Element] = []
    # top-level fragments only (nested ones belong to abbreviations); groups are flattened
    for page in root.findall(".//page"):
        stack: List[ET.Element] = [page]
        while stack:
            parent = stack.pop()
            for el in list(parent):
                tag = el.tag
                if tag == "group":
                    stack.append(el)
                    groups.append(el)
                    continue
                if tag == "fragment":
                    if not el.findall("n"):
                        # container for graphics only (orbitals, curves) - not a molecule
                        other["fragment (no atoms)"] = other.get("fragment (no atoms)", 0) + 1
                        continue
                    box = _box_of(el)
                    entry: Dict[str, Any] = {"id": el.get("id"), "kind": "molecule", "bounds": (box.to_dict() if box else None),
                                             "atoms": len(el.findall("n")), "bonds": len(el.findall("b"))}
                    for sub in el.iter():
                        if sub is not el and sub.get("id"):
                            node_owner[sub.get("id")] = entry
                    if want_mols:
                        try:
                            from rdkit import Chem, rdBase
                            from rdkit.Chem import Descriptors, rdMolDescriptors
                            mol, _ = fragment_to_mol(el)
                            with rdBase.BlockLogs():
                                entry["smiles"] = Chem.MolToSmiles(mol)
                                try:
                                    entry["formula"] = rdMolDescriptors.CalcMolFormula(mol)
                                    entry["molecular_weight"] = round(Descriptors.MolWt(mol), 3)
                                except Exception:
                                    pass
                            entry["heavy_atoms"] = mol.GetNumHeavyAtoms()
                            if mol.HasProp("_cdxml_labels"):
                                entry["labels"] = mol.GetProp("_cdxml_labels").split("\n")
                            if mol.HasProp("_cdxml_warnings"):
                                entry["warnings"] = mol.GetProp("_cdxml_warnings").split("\n")
                            if mol.HasProp("_cdxml_errors"):
                                entry["error"] = "could not fully interpret fragment (%s); SMILES is a best-effort reading" % (
                                    "; ".join(mol.GetProp("_cdxml_errors").split("\n")))
                                out["warnings"].append("fragment %s: %s" % (el.get("id"), entry["error"]))
                            if mol.GetNumAtoms() == 0:
                                entry["error"] = "fragment has no atoms"
                            elif mol.GetNumAtoms() == 1 and mol.GetAtomWithIdx(0).GetAtomicNum() == 0 and entry.get("labels"):
                                # a caption typed as an atom label ("1. Br2, H2O"): keep its text for reaction steps
                                entry["text"] = entry["labels"][0]
                        except Exception as e:  # pragma: no cover - depends on exotic content
                            entry["error"] = "could not interpret fragment: %s" % e
                            out["warnings"].append("fragment %s: %s" % (el.get("id"), e))
                    out["molecules"].append(entry)
                    by_id[el.get("id", "")] = entry
                elif tag == "t":
                    txt = _text_of(el)
                    box = _box_of(el)
                    entry = {"id": el.get("id"), "kind": "text", "text": txt, "bounds": (box.to_dict() if box else None)}
                    cp = el.find("chemicalproperty")
                    if cp is not None:
                        entry["chemical_property"] = {k: v for k, v in cp.attrib.items()}
                    out["texts"].append(entry)
                    by_id[el.get("id", "")] = entry
                elif tag == "arrow":
                    box = _box_of(el)
                    entry = {"id": el.get("id"), "kind": "arrow", "bounds": (box.to_dict() if box else None)}
                    for k in ("Head3D", "Tail3D", "ArrowheadHead", "ArrowheadTail", "ArrowheadType", "LineType", "NoGo", "AngularSize"):
                        if el.get(k):
                            entry[k[0].lower() + k[1:]] = el.get(k)
                    out["arrows"].append(entry)
                    by_id[el.get("id", "")] = entry
                elif tag == "graphic":
                    if el.get("SupersededBy"):
                        continue   # legacy twin of an <arrow> element written by ChemDraw for old readers
                    gt = el.get("GraphicType", "")
                    box = _box_of(el)
                    entry = {"id": el.get("id"), "kind": "graphic", "type": gt, "bounds": (box.to_dict() if box else None)}
                    for k in ("ArrowType", "SymbolType", "RectangleType", "OvalType", "LineType"):
                        if el.get(k):
                            entry[k[0].lower() + k[1:]] = el.get(k)
                    if gt == "Line" and el.get("ArrowType"):
                        out["arrows"].append(entry)
                    else:
                        out["graphics"].append(entry)
                    by_id[el.get("id", "")] = entry
                elif tag in ("chemicalproperty", "scheme", "step", "page"):
                    continue
                else:
                    other[tag] = other.get(tag, 0) + 1
    for g in groups:
        # reaction steps may reference a group; it stands for its member objects
        members = [sub.get("id") for sub in g.iter() if sub is not g and sub.get("id") in by_id]
        by_id.setdefault(g.get("id", ""), {"id": g.get("id"), "kind": "group", "members": members})
    if other:
        out["other_objects"] = other
    spectra = root.findall(".//spectrum")
    if spectra:
        out["spectra"] = len(spectra)
    brackets = root.findall(".//bracketedgroup")
    if brackets:
        out["bracketed_groups"] = []
        for bg in brackets:
            info = {"id": bg.get("id"), "objects": len((bg.get("BracketedObjectIDs") or "").split())}
            for k in ("RepeatCount", "BracketUsage", "PolymerRepeatPattern", "PolymerFlipType", "SRULabel", "ComponentOrder"):
                if bg.get(k):
                    info[k[0].lower() + k[1:]] = bg.get(k)
            out["bracketed_groups"].append(info)
        out["warnings"].append("%d bracketed group(s) (polymer / repeat brackets): repeat units are not expanded in the SMILES" % len(brackets))
    # chemical properties (names, residue labels ...) are usually siblings of the <t> that displays them
    for parent in root.iter():
        kids = list(parent)
        for i, cp in enumerate(kids):
            if cp.tag != "chemicalproperty" or parent.tag == "t":
                continue
            attrs = {k: v for k, v in cp.attrib.items()}
            target = by_id.get(cp.get("ChemicalPropertyDisplayID") or "")
            if target is None or target.get("kind") != "text":
                prev = next((kids[j] for j in range(i - 1, -1, -1) if kids[j].tag == "t"), None)
                target = by_id.get(prev.get("id", "")) if prev is not None else None
            if target is not None and target.get("kind") == "text":
                target.setdefault("chemical_property", attrs)
                text_value = target.get("text", "")
            else:
                text_value = ""
            owners: List[Dict[str, Any]] = []
            for oid in (cp.get("BasisObjects") or "").split():
                mol_entry = by_id.get(oid) if by_id.get(oid, {}).get("kind") == "molecule" else node_owner.get(oid)
                if mol_entry is not None and mol_entry not in owners:
                    owners.append(mol_entry)
            for mol_entry in owners:
                prop = {"text": text_value, "type": (cp.get("ChemicalPropertyType") or "").strip()}
                if cp.get("Name"):
                    prop["name"] = cp.get("Name")
                mol_entry.setdefault("properties", []).append(prop)
                if text_value and prop["type"] == "1" and "name" not in mol_entry:
                    mol_entry["name"] = text_value
    for step in root.findall(".//step"):
        def ids(key):
            return [i for i in (step.get(key) or "").split()]

        def describe(idlist):
            res = []
            for i in idlist:
                e = by_id.get(i)
                if e is None:
                    res.append({"id": i})
                elif e["kind"] == "group":
                    res.extend(describe(e.get("members", [])))
                elif e["kind"] == "molecule":
                    if e.get("text") and e.get("smiles") == "*":
                        res.append({"id": i, "text": e["text"]})
                    else:
                        res.append({"id": i, "smiles": e.get("smiles"), "formula": e.get("formula")})
                elif e["kind"] == "text":
                    res.append({"id": i, "text": e.get("text")})
                else:
                    res.append({"id": i, "kind": e["kind"]})
            return res
        out["steps"].append({
            "id": step.get("id"),
            "reactants": describe(ids("ReactionStepReactants")),
            "products": describe(ids("ReactionStepProducts")),
            "arrows": ids("ReactionStepArrows"),
            "above_arrow": describe(ids("ReactionStepObjectsAboveArrow")),
            "below_arrow": describe(ids("ReactionStepObjectsBelowArrow")),
            "atom_map_pairs": len(ids("ReactionStepAtomMap")) // 2,
        })
    if out["steps"]:
        out["reaction_smiles"] = [step_reaction_smiles(s) for s in out["steps"]]
    return out


def step_reaction_smiles(step: Dict[str, Any]) -> str:
    """``reactants>agents>products`` for one step of :func:`read_cdxml`; uninterpretable reactants/products are ``*``,
    agents are limited to real molecules (text labels above/below the arrow are not included)."""
    r = ".".join(x.get("smiles") or "*" for x in step.get("reactants", []) if not (x.get("text") and not x.get("smiles")))
    p = ".".join(x.get("smiles") or "*" for x in step.get("products", []) if not (x.get("text") and not x.get("smiles")))
    agents = ".".join(x.get("smiles") for x in step.get("above_arrow", []) + step.get("below_arrow", [])
                      if x.get("smiles") and x.get("smiles") != "*")
    return "%s>%s>%s" % (r, agents, p)


def cdxml_mols(text: str) -> list:
    """All top-level molecules of a CDXML document as RDKit mols (with 2D conformers).

    Fragments that cannot be interpreted at all are skipped; partially interpreted ones carry the
    ``_cdxml_errors`` / ``_cdxml_warnings`` properties (see :func:`fragment_to_mol`).
    """
    root = parse_cdxml(text)
    mols = []
    for page in root.findall(".//page"):
        stack = [page]
        while stack:
            parent = stack.pop()
            for el in list(parent):
                if el.tag == "group":
                    stack.append(el)
                elif el.tag == "fragment" and el.findall("n"):
                    try:
                        mol = fragment_to_mol(el)[0]
                    except Exception:
                        continue
                    if mol.GetNumAtoms():
                        mol.SetProp("_cdxml_id", el.get("id", ""))
                        mols.append(mol)
    return mols


def document_text(text: str) -> List[str]:
    """All caption strings in a document (not atom labels)."""
    return [t["text"] for t in read_cdxml(text, want_mols=False)["texts"]]
