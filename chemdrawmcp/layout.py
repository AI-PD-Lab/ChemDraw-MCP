"""Layout engines on top of :mod:`chemdrawmcp.cdxml`.

* :func:`build_structures` - one or more molecules in a grid with optional labels.
* :func:`build_reaction`   - one reaction: reactants + reactants -> products with reagents above / conditions
                             below the arrow, optional yields, labels and atom mapping.
* :func:`build_scheme`     - multi-step scheme: a chain of intermediates and reagents laid out left to right,
                             wrapping into rows (intermediate repeated at the start of a new row).
* :func:`collisions` / :func:`layout_warnings` - cheap overlap check; ``build_*`` attach ``doc.warnings``.

All functions take RDKit molecules prepared for depiction (2D conformer, Kekulé bonds, wedge directions)
and return a :class:`cdxml.CDXMLDocument`.

Geometry conventions (points, y down): a reaction row is centred on a baseline ``y``; every component's
label-aware fragment box (see :func:`textmetrics.fragment_labels`) is centred on it, the "+" glyph is centred
on it visually and the arrow lies on it.  Captions hang below the structures; reagent lines stack upward
from the arrow and condition lines downward, never closer than 0.35 bond lengths to the arrow line.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import textmetrics as _tm
from .cdxml import Box, CDXMLDocument, ObjRef, _mol_points, text_box, union_boxes

CAPTION_GAP = 0.35      # bond lengths between a structure's label-aware box and its caption
PLUS_GAP = 0.6          # minimum bond lengths between a component and the "+" on either side
ARROW_CLEARANCE = 0.4   # bond lengths between the arrow line and reagent / condition text (spec: >= 0.35)
ARROW_TEXT_MARGIN = 0.8  # the arrow is at least this many bond lengths longer than its widest text line
LINE_GAP = 0.2          # fraction of the text size between stacked reagent / condition lines
ARROW_LINE_PAD = 2.0    # points: half the arrowhead width, used for the collision box of straight arrows
_SHAFT_HALF = {"equilibrium": 1.5, "unbalanced_equilibrium": 1.5, "retrosynthetic": 2.5, "retro": 2.5}


@dataclass
class Component:
    """One molecule in a reaction with its presentation options."""
    mol: Any
    label: Optional[str] = None        # caption drawn below the structure (name, number, yield...)
    color: Any = None
    atom_colors: Optional[Dict[int, Any]] = None
    label_atoms: bool = False

    @classmethod
    def coerce(cls, obj) -> "Component":
        if isinstance(obj, Component):
            return obj
        if isinstance(obj, dict):
            return cls(obj["mol"], obj.get("label"), obj.get("color"), obj.get("atom_colors"), bool(obj.get("label_atoms")))
        return cls(obj)


@dataclass
class Step:
    """One reaction step for :func:`build_scheme` / :func:`build_reaction`."""
    reactants: List[Component]
    products: List[Component]
    above: List[str] = field(default_factory=list)      # reagent lines above the arrow
    below: List[str] = field(default_factory=list)      # condition lines below the arrow
    arrow: str = "full"
    atom_map: bool = True                               # emit ReactionStepAtomMap from atom map numbers


# ---------------------------------------------------------------------------
# measurement helpers (no document mutation)
# ---------------------------------------------------------------------------

def _caption_font(doc: CDXMLDocument) -> str:
    return doc.style.get("caption_font", doc.style.get("label_font", "Arial"))


def _cache(doc: CDXMLDocument) -> Dict[Any, Any]:
    """Per-build memo (fragment boxes, canonical SMILES) keyed by id(mol); reset by every ``build_*`` call so
    that ids of molecules freed between calls can never alias."""
    c = doc.__dict__.get("_layout_cache")
    if c is None:
        c = doc.__dict__["_layout_cache"] = {}
    return c


def _begin(doc: CDXMLDocument) -> CDXMLDocument:
    doc.__dict__["_layout_cache"] = {}
    return doc


def _fragment_box(doc: CDXMLDocument, mol, label_atoms: bool = False) -> Box:
    """Label-aware box of a fragment in its own frame - the same computation as ``CDXMLDocument.add_fragment``."""
    cache = _cache(doc)
    key = ("box", id(mol), bool(label_atoms))
    hit = cache.get(key)
    if hit is not None:
        return hit
    pts = _mol_points(mol, doc.bond_length)
    if not pts:
        raise ValueError("molecule has no atoms")
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    box = Box(min(xs), min(ys), max(xs), max(ys)).padded(0.35 * doc.bond_length)
    st = doc.style
    for lab in _tm.fragment_labels(mol, pts, doc.label_size, label_atoms, st.get("label_font", "Arial"),
                                   float(st.get("margin_width", 1.6)), show_h=not st.get("hide_implicit_hydrogens"),
                                   show_terminal_carbons=bool(st.get("show_terminal_carbon_labels"))):
        box = box.union(Box(*lab["box"]))
    cache[key] = box
    return box


def _measure_fragment(doc: CDXMLDocument, mol, label_atoms: bool = False) -> Tuple[float, float]:
    """Width/height a fragment will occupy (label-aware, same as add_fragment)."""
    b = _fragment_box(doc, mol, label_atoms)
    return b.width, b.height


def _component_extent(doc: CDXMLDocument, comp: Component, label_size: float) -> Tuple[float, float, float, Box]:
    """(width, extent above the centre line, extent below it incl. the caption, fragment box) of a component."""
    fb = _fragment_box(doc, comp.mol, comp.label_atoms)
    w = fb.width
    up = down = fb.height / 2.0
    if comp.label:
        cw, asc, desc = text_box(comp.label, label_size, _caption_font(doc))
        w = max(w, cw)
        down += CAPTION_GAP * doc.bond_length + asc + desc
    return w, up, down, fb


def _plus_metrics(doc: CDXMLDocument, gap: float) -> Tuple[float, float, float]:
    """(width of the '+', clearance on each side, vertical centre offset to pass to add_plus) for a reaction row."""
    size = doc.caption_size + 2          # add_plus default
    font = _caption_font(doc)
    w, asc, desc = text_box("+", size, font)
    clearance = max(PLUS_GAP * doc.bond_length, 0.5 * gap)
    # add_plus centres the text *box* on y (baseline = y + (asc - desc) / 2); shift so the glyph centre sits on y
    yoff = _tm.glyph_center("+", size, font) - (asc - desc) / 2.0
    return w, clearance, yoff


def _looks_like_formula(text: str) -> bool:
    """Turn on ChemDraw's formula face for reagent-like text (digits after element symbols, charges)."""
    import re
    if re.search(r"[A-Z][a-z]?\d", text) or re.search(r"\d[+-]", text):
        return True
    return False


def _arrow_clearance(doc: CDXMLDocument, kind: str) -> float:
    k = (kind or "full").lower().replace("-", "_").replace(" ", "_")
    return ARROW_CLEARANCE * doc.bond_length + _SHAFT_HALF.get(k, 0.0)


def _text_stack(doc: CDXMLDocument, lines: Sequence[str], text_size: float) -> Tuple[float, float]:
    """(widest line, total height) of stacked reagent / condition lines."""
    widest = height = 0.0
    font = _caption_font(doc)
    for k, ln in enumerate(lines):
        w, asc, desc = text_box(ln, text_size, font, _looks_like_formula(ln))
        widest = max(widest, w)
        height += asc + desc + (LINE_GAP * text_size if k else 0.0)
    return widest, height


def _arrow_length(doc: CDXMLDocument, step: Step, min_length: float, text_size: float) -> float:
    widest, _ = _text_stack(doc, list(step.above) + list(step.below), text_size)
    return max(min_length, widest + ARROW_TEXT_MARGIN * doc.bond_length)


def _step_extents(doc: CDXMLDocument, step: Step, text_size: float, label_size: float) -> Tuple[float, float]:
    """(extent above the baseline, extent below it) of a whole step: structures, captions, arrow text."""
    up = down = 0.0
    for c in list(step.reactants) + list(step.products):
        _, u, d, _ = _component_extent(doc, c, label_size)
        up, down = max(up, u), max(down, d)
    clear = _arrow_clearance(doc, step.arrow)
    if step.above:
        up = max(up, clear + _text_stack(doc, step.above, text_size)[1])
    if step.below:
        down = max(down, clear + _text_stack(doc, step.below, text_size)[1])
    return up, down


# ---------------------------------------------------------------------------
# placement helpers
# ---------------------------------------------------------------------------

def _place_component(doc: CDXMLDocument, comp: Component, x_left: float, y_center: float, label_size: float
                     ) -> Tuple[ObjRef, List[ObjRef], Box]:
    """Place one component with its left edge (structure or caption, whichever is wider) at x_left.

    Returns (fragment ref, caption refs, box of structure + caption)."""
    w, _, _, fb = _component_extent(doc, comp, label_size)
    ref = doc.add_fragment(comp.mol, x_left + (w - fb.width) / 2.0, y_center, anchor="left", color=comp.color,
                           atom_colors=comp.atom_colors, label_atoms=comp.label_atoms)
    labels: List[ObjRef] = []
    box = ref.box
    if comp.label:
        t = doc.add_text(comp.label, ref.box.cx, ref.box.y1 + CAPTION_GAP * doc.bond_length, size=label_size, anchor="top")
        labels.append(t)
        box = box.union(t.box)
    return ref, labels, box


def _add_plus(doc: CDXMLDocument, x_prev_right: float, y: float, gap: float) -> ObjRef:
    """A '+' to the right of a component box with the required clearance, visually centred on y."""
    w, clearance, yoff = _plus_metrics(doc, gap)
    return doc.add_plus(x_prev_right + clearance + w / 2.0, y + yoff)


def _place_side(doc: CDXMLDocument, comps: Sequence[Component], x: float, y: float, gap: float, label_size: float,
                lead: Optional[Box] = None) -> Tuple[List[ObjRef], List[ObjRef], List[ObjRef], List[Box], float]:
    """Lay out components left to right with '+' between them.

    ``lead``: box of an already placed component that precedes these (a '+' is inserted after it).
    Returns (fragment refs, plus refs, caption refs, component boxes, x after the last component)."""
    refs: List[ObjRef] = []
    plusses: List[ObjRef] = []
    labels: List[ObjRef] = []
    boxes: List[Box] = []
    prev = lead
    for comp in comps:
        if prev is not None:
            plus = _add_plus(doc, prev.x1, y, gap)
            plusses.append(plus)
            x = plus.box.x1 + _plus_metrics(doc, gap)[1]
        ref, labs, cbox = _place_component(doc, comp, x, y, label_size)
        refs.append(ref)
        labels.extend(labs)
        boxes.append(cbox)
        prev = cbox
        x = cbox.x1
    return refs, plusses, labels, boxes, x


def _arrow_with_text(doc: CDXMLDocument, x: float, y: float, above: Sequence[str], below: Sequence[str],
                     kind: str, min_length: float, text_size: float) -> Tuple[ObjRef, List[ObjRef], List[ObjRef], float]:
    """Arrow from x along the baseline y, long enough for its text; reagent lines stacked above (first line on
    top), condition lines below, each line its own text object.  Returns (arrow, above refs, below refs, x end)."""
    L = doc.bond_length
    font = _caption_font(doc)
    above = [ln for ln in above if ln is not None and str(ln) != ""]
    below = [ln for ln in below if ln is not None and str(ln) != ""]
    widest, _ = _text_stack(doc, above + below, text_size)
    length = max(min_length, widest + ARROW_TEXT_MARGIN * L)
    arrow = doc.add_arrow(x, y, x + length, y, kind)
    cx = x + length / 2.0
    clear = _arrow_clearance(doc, kind)
    line_gap = LINE_GAP * text_size
    above_refs: List[ObjRef] = []
    bottom = y - clear
    for ln in reversed(above):
        f = _looks_like_formula(ln)
        _, _, desc = text_box(ln, text_size, font, f)
        ref = doc.add_text(ln, cx, bottom - desc, size=text_size, anchor="baseline", formula=f)
        above_refs.insert(0, ref)
        bottom = ref.box.y0 - line_gap
    below_refs: List[ObjRef] = []
    top = y + clear
    for ln in below:
        f = _looks_like_formula(ln)
        ref = doc.add_text(ln, cx, top, size=text_size, anchor="top", formula=f)
        below_refs.append(ref)
        top = ref.box.y1 + line_gap
    return arrow, above_refs, below_refs, x + length


# ---------------------------------------------------------------------------
# overlap check
# ---------------------------------------------------------------------------

def _collision_box(ref: ObjRef) -> Optional[Box]:
    """Box used for the overlap check: fragments and texts as they are, straight arrows by their line (the padded
    arrow box legitimately contains its reagent text).  Curved arrows and graphics overlap things on purpose."""
    if ref.kind == "arrow":
        if ref.extra.get("kind") == "curved":
            return None
        tail, head = ref.extra.get("tail"), ref.extra.get("head")
        if tail and head:
            return Box(min(tail[0], head[0]), min(tail[1], head[1]), max(tail[0], head[0]), max(tail[1], head[1])).padded(ARROW_LINE_PAD)
        return ref.box
    if ref.kind in ("fragment", "text"):
        return ref.box
    return None


def collisions(doc: CDXMLDocument, tolerance: float = 0.01) -> List[Tuple[ObjRef, ObjRef]]:
    """Pairs of document objects whose boxes overlap (sweep over x; O(n log n) for well-spread layouts)."""
    items = []
    for r in doc.refs:
        b = _collision_box(r)
        if b is not None:
            items.append((b, r))
    items.sort(key=lambda t: t[0].x0)
    out: List[Tuple[ObjRef, ObjRef]] = []
    for i, (bi, ri) in enumerate(items):
        for bj, rj in items[i + 1:]:
            if bj.x0 >= bi.x1 - tolerance:
                break
            if bj.y0 < bi.y1 - tolerance and bi.y0 < bj.y1 - tolerance:
                out.append((ri, rj))
    return out


def _describe(ref: ObjRef) -> str:
    if ref.kind == "fragment":
        return "fragment %d (%s)" % (ref.id, ref.extra.get("smiles", ""))
    if ref.kind == "text":
        return "text %d (%r)" % (ref.id, ref.extra.get("text", ""))
    return "%s %d" % (ref.kind, ref.id)


def layout_warnings(doc: CDXMLDocument) -> List[str]:
    """Human-readable layout problems: overlapping objects and content wider than the page."""
    out: List[str] = []
    for a, b in collisions(doc):
        ba, bb = _collision_box(a), _collision_box(b)
        ox = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
        oy = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
        out.append("%s overlaps %s by %.1f x %.1f pt" % (_describe(a), _describe(b), ox, oy))
    content = doc.content_box()
    if content is not None:
        page = doc.style.get("page", [0, 0, 540, 720])
        page_w = page[2] - page[0]
        if content.x1 > page_w:
            out.append("content is %.0f pt wide, wider than the %.0f pt page: ChemDraw tiles extra pages" % (content.x1, page_w))
    return out


def _finish(doc: CDXMLDocument) -> CDXMLDocument:
    doc.warnings = layout_warnings(doc)
    return doc


# ---------------------------------------------------------------------------

def build_structures(mols: Sequence[Any], labels: Optional[Sequence[Optional[str]]] = None, style: Optional[str] = None,
                     columns: Optional[int] = None, bond_length: Optional[float] = None, gap: Optional[float] = None,
                     colors: Optional[Sequence[Any]] = None, label_size: Optional[float] = None,
                     settings: Optional[Dict[str, Any]] = None, name: str = "structures",
                     origin: Tuple[float, float] = (36.0, 36.0), doc: Optional[CDXMLDocument] = None,
                     label_atoms: Optional[Sequence[bool]] = None) -> CDXMLDocument:
    """Molecules in a grid.  ``columns`` defaults to a roughly square arrangement (max 4).  ``origin`` = top-left.

    Every column is as wide as its widest structure or caption; every molecule is centred in its cell horizontally
    and vertically, its caption centred below it."""
    doc = _begin(doc or CDXMLDocument(style, name=name, bond_length=bond_length, settings=settings))
    L = doc.bond_length
    n = len(mols)
    if n == 0:
        return _finish(doc)
    columns = columns or min(4, max(1, int(math.ceil(math.sqrt(n)))))
    gap = gap if gap is not None else 1.5 * L
    label_size = label_size or doc.caption_size
    comps = [Component(mols[i], labels[i] if labels and i < len(labels) else None,
                       colors[i] if colors and i < len(colors) else None) for i in range(n)]
    ext = [_component_extent(doc, c, label_size) for c in comps]
    rows = int(math.ceil(n / columns))
    col_w = [0.0] * columns
    row_h = [0.0] * rows          # tallest structure per row
    row_cap = [0.0] * rows        # caption space per row
    for i, (w, up, down, fb) in enumerate(ext):
        r, c = divmod(i, columns)
        col_w[c] = max(col_w[c], w)
        row_h[r] = max(row_h[r], fb.height)
        row_cap[r] = max(row_cap[r], down - fb.height / 2.0)
    x0, y0 = origin
    y = y0
    for r in range(rows):
        x = x0
        for c in range(columns):
            i = r * columns + c
            if i >= n:
                break
            comp = comps[i]
            cx = x + col_w[c] / 2.0
            cy = y + row_h[r] / 2.0
            ref = doc.add_fragment(comp.mol, cx, cy, color=comp.color, label_atoms=bool(label_atoms[i]) if label_atoms and i < len(label_atoms) else False)
            if comp.label:
                doc.add_text(comp.label, ref.box.cx, ref.box.y1 + CAPTION_GAP * L, size=label_size, anchor="top")
            x += col_w[c] + gap
        y += row_h[r] + row_cap[r] + gap
    return _finish(doc)


def build_reaction(step: Step, style: Optional[str] = None, bond_length: Optional[float] = None,
                   settings: Optional[Dict[str, Any]] = None, name: str = "reaction", origin: Tuple[float, float] = (36.0, 72.0),
                   doc: Optional[CDXMLDocument] = None, text_size: Optional[float] = None, min_arrow: Optional[float] = None
                   ) -> CDXMLDocument:
    """One reaction on one line.  Appends to ``doc`` when given (origin = left x, baseline y)."""
    doc = _begin(doc or CDXMLDocument(style, name=name, bond_length=bond_length, settings=settings))
    L = doc.bond_length
    gap = 1.0 * L
    text_size = text_size or doc.caption_size
    x, y = origin
    # every component's label-aware box is centred on the baseline y; so are the plusses and the arrow
    reactant_refs, plusses, labels, _, x = _place_side(doc, step.reactants, x, y, gap, text_size)
    arrow, above, below, x = _arrow_with_text(doc, x + gap, y, step.above, step.below, step.arrow, min_arrow or 3.0 * L, text_size)
    product_refs, plusses2, labels2, _, x = _place_side(doc, step.products, x + gap, y, gap, text_size)
    atom_map = _atom_map_pairs(reactant_refs, product_refs, step) if step.atom_map else None
    doc.add_step(reactant_refs, product_refs, [arrow], above, below, plusses + plusses2, atom_map)
    return _finish(doc)


def _atom_map_pairs(reactant_refs: Sequence[ObjRef], product_refs: Sequence[ObjRef], step: Step,
                    reactants: Optional[Sequence[Component]] = None) -> List[Tuple[int, int]]:
    """Pairs of CDXML node ids for atoms that share an atom-map number on both sides.

    ``reactants``: the components the reactant fragments were drawn from when they are not ``step.reactants``
    (a scheme reuses the previous product's fragment, whose ``atom_ids`` follow *that* molecule's atom order)."""
    def collect(refs, comps):
        out = {}
        for ref, comp in zip(refs, comps):
            mol = comp.mol
            for i in range(mol.GetNumAtoms()):
                m = mol.GetAtomWithIdx(i).GetAtomMapNum()
                if m:
                    out[m] = ref.extra["atom_ids"][i]
        return out
    r = collect(reactant_refs, step.reactants if reactants is None else reactants)
    p = collect(product_refs, step.products)
    return [(r[m], p[m]) for m in sorted(r) if m in p]


def build_scheme(steps: Sequence[Step], style: Optional[str] = None, bond_length: Optional[float] = None,
                 settings: Optional[Dict[str, Any]] = None, name: str = "scheme", max_width: Optional[float] = None,
                 text_size: Optional[float] = None, row_gap: Optional[float] = None, title: Optional[str] = None,
                 chain: bool = True, origin: Tuple[float, float] = (36.0, 72.0), doc: Optional[CDXMLDocument] = None
                 ) -> CDXMLDocument:
    """Multi-step scheme.

    chain=True: consecutive steps share the intermediate (product of step i is drawn once and reused as the
    first reactant of step i+1 when it is the same molecule).  A step that does not continue from the previous
    product (every step when chain=False) starts a new row: steps are never drawn side by side on one baseline,
    which would read as a single reaction.  Rows wrap at ``max_width`` points (default: the style's page width
    minus margins); a new row restarts with the intermediate drawn again.  Row heights account for captions,
    atom labels and arrow text of every step in the row.
    """
    doc = _begin(doc or CDXMLDocument(style, name=name, bond_length=bond_length, settings=settings))
    L = doc.bond_length
    gap = 1.0 * L
    text_size = text_size or doc.caption_size
    page_w = doc.style.get("page", [0, 0, 540, 720])[2]
    max_width = max_width or max(300.0, page_w - 72.0)
    row_gap = row_gap if row_gap is not None else 2.5 * L
    x0, y = origin
    extents = [_step_extents(doc, st, text_size, text_size) for st in steps]
    up_all = max([u for u, _ in extents] or [L])
    if title:
        t = doc.add_text(title, x0, y - up_all - 1.5 * L, size=text_size + 1, face="bold",
                         justification="Left", anchor="top")
        y = t.box.y1 + 2.0 * L + up_all
    # pass 1: chaining and greedy row assignment (the width estimate mirrors the placement code below).  A step
    # that does not continue from the previous product gets its own row; a chained step wraps when it does not fit.
    rows: List[List[int]] = []
    reuse_flags: List[bool] = []
    x = x0
    for si, step in enumerate(steps):
        reuse = bool(chain and si > 0 and len(step.reactants) >= 1 and len(steps[si - 1].products) >= 1
                     and _same_mol(step.reactants[0].mol, steps[si - 1].products[-1].mol, doc))
        est = _estimate_step_width(doc, step, gap, text_size, skip_first_reactant=reuse)
        if reuse and x + est <= x0 + max_width:
            rows[-1].append(si)
            x += est                      # the step starts at the reused product's right edge (x - gap)
        else:
            if reuse:                     # does not fit: wrap and draw the intermediate again
                reuse = False
                est = _estimate_step_width(doc, step, gap, text_size)
            rows.append([si])
            x = x0 + est + gap
        reuse_flags.append(reuse)
    # pass 2: placement
    row_bottom = y
    last_product_refs: List[ObjRef] = []
    last_products: List[Component] = []
    last_boxes: List[Box] = []
    for ri, row in enumerate(rows):
        if ri > 0:
            # the baseline clears the tallest member of the whole row (structures, captions, reagent stacks)
            y = row_bottom + row_gap + max(extents[j][0] for j in row)
        x = x0
        for si in row:
            step = steps[si]
            if reuse_flags[si]:
                first = last_product_refs[-1]
                lead = last_boxes[-1] if last_boxes else first.box
                reactant_refs = [first]
                # the reused fragment's atom_ids follow the previous product's atom order, not this step's reactant
                map_reactants = [last_products[-1]] + list(step.reactants[1:])
                if len(step.reactants) > 1:
                    refs2, plusses, labels, _, x = _place_side(doc, step.reactants[1:], lead.x1, y, gap, text_size, lead=lead)
                    reactant_refs += refs2
                else:
                    plusses, labels = [], []
                    x = lead.x1
            else:
                reactant_refs, plusses, labels, _, x = _place_side(doc, step.reactants, x, y, gap, text_size)
                map_reactants = list(step.reactants)
            arrow, above, below, x = _arrow_with_text(doc, x + gap, y, step.above, step.below, step.arrow, 3.0 * L, text_size)
            product_refs, plusses2, labels2, boxes2, x = _place_side(doc, step.products, x + gap, y, gap, text_size)
            atom_map = _atom_map_pairs(reactant_refs, product_refs, step, map_reactants) if step.atom_map else None
            doc.add_step(reactant_refs, product_refs, [arrow], above, below, plusses + plusses2, atom_map)
            last_product_refs, last_products, last_boxes = product_refs, list(step.products), boxes2
            bottom = union_boxes([r.box for r in reactant_refs + product_refs + labels + labels2 + below + plusses + plusses2])
            row_bottom = max(row_bottom, bottom.y1 if bottom else y)
            x += gap
    return _finish(doc)


def _max_height(doc: CDXMLDocument, steps: Sequence[Step]) -> float:
    """Tallest label-aware structure height among the steps' components (at least two bond lengths)."""
    h = doc.bond_length * 2
    for st in steps:
        for c in list(st.reactants) + list(st.products):
            h = max(h, _measure_fragment(doc, c.mol, c.label_atoms)[1])
    return h


def _canonical(doc: CDXMLDocument, mol) -> str:
    cache = _cache(doc)
    key = ("smiles", id(mol))
    smi = cache.get(key)
    if smi is None:
        try:
            from rdkit import Chem
            smi = Chem.MolToSmiles(mol)
        except Exception:
            smi = ""
        cache[key] = smi
    return smi


def _same_mol(a, b, doc: Optional[CDXMLDocument] = None) -> bool:
    if doc is not None:
        sa, sb = _canonical(doc, a), _canonical(doc, b)
        return bool(sa) and sa == sb
    try:
        from rdkit import Chem
        return Chem.MolToSmiles(a) == Chem.MolToSmiles(b)
    except Exception:
        return False


def _estimate_step_width(doc: CDXMLDocument, step: Step, gap: float, text_size: float, skip_first_reactant: bool = False) -> float:
    """Width a step will take when placed by build_reaction / build_scheme (mirrors the placement code)."""
    L = doc.bond_length
    pw, pclear, _ = _plus_metrics(doc, gap)
    plus_w = pw + 2 * pclear
    comps = list(step.reactants[1:] if skip_first_reactant else step.reactants)
    w = sum(_component_extent(doc, c, text_size)[0] for c in comps)
    n_plus = max(0, len(comps) - 1) + (1 if skip_first_reactant and comps else 0)
    w += n_plus * plus_w
    w += gap + _arrow_length(doc, step, 3.0 * L, text_size) + gap
    w += sum(_component_extent(doc, c, text_size)[0] for c in step.products)
    w += max(0, len(step.products) - 1) * plus_w
    return w
