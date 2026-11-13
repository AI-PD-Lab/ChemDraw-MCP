"""Electron-pushing (curved "mechanism") arrows and mechanism layouts.

Builds on :mod:`chemdrawmcp.cdxml`: for a fragment placed with
:meth:`CDXMLDocument.add_fragment` this module works out where the tail and
the head of a curved arrow belong (lone pairs, bonds, bonds that are being
formed), to which side and how far the arc should bulge, and draws it with
:meth:`CDXMLDocument.add_curved_arrow`.  :func:`build_mechanism` lays several
such stages out left to right with reaction / resonance arrows in between,
wrapping into rows like :func:`chemdrawmcp.layout.build_scheme`.

Atom indices are 0-based RDKit indices of the molecule handed to
``add_fragment`` (SMILES atom order for heavy atoms).  Explicit ``[H]`` atoms
written in a SMILES are kept when a stage is built from it, so protons that
take part in a mechanism can be addressed.  In a stage with several
structures an atom is addressed as ``fragment:atom`` (``"1:0"`` = atom 0 of
the second structure).

Compact arrow specs understood by :func:`parse_arrow_spec`::

    lp:5->5-8         lone pair on atom 5 to the new bond 5-8
    3-4->4            bond 3-4 collapses onto atom 4 (becomes a lone pair / charge)
    3=4->4-7          pi bond 3=4 to the new bond 4-7
    2-3->3=4          bond 2-3 moves into the existing bond 3-4 (new pi bond); the electrons that
                      turn a bond into a double bond must come from elsewhere (a lone pair,
                      a C-H bond or a neighbouring bond), never from that bond itself
    h:3->3-4          the C-H bond on atom 3 (implicit hydrogen) to the bond 3-4
    lp:5->5-8:1e      one electron (fishhook head); "lp:5~>5-8" is the same
    lp:0:0->0:0-1:0   fragment 0 atom 0 attacks fragment 1 atom 0

Arc sign convention (verified in ChemDraw): :meth:`add_curved_arrow` with a
positive angle bulges to the right of the tail->head direction on screen
(downwards for a left-to-right arrow), a negative angle to the left.
``bulge="auto"`` picks the side away from the molecule's centroid.
"""
from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Sequence, Tuple, Union

from .cdxml import Box, CDXMLDocument, ObjRef, union_boxes

__all__ = [
    "AtomRef", "ElectronArrow", "Stage", "FragmentGeometry", "add_electron_arrows", "parse_arrow_spec",
    "build_mechanism", "describe_mechanism", "fragment_geometry", "fragment_atom_positions", "arc_points", "arc_apex",
]

Point = Tuple[float, float]
AtomKey = Tuple[int, int]          # (fragment index, atom index)
Obstacle = Union[Box, Tuple[Point, float]]   # a box, or (centre, radius)

# distances in units of the bond length L
ATOM_OFFSET = 0.35          # arrow head outside a bare (carbon) atom
LABEL_OFFSET = 0.42         # arrow head outside a labelled atom (O, N, Cl, charged C ...)
LONE_PAIR_OFFSET = 0.5      # arrow tail leaving a lone pair / anion
MAX_ATOM_OFFSET = 0.6       # normally never further than this from the referenced atom ...
MAX_LABEL_OFFSET = 1.0      # ... unless a wide label ("OH-", "H3N+") has to be cleared
BOND_OFFSET = 0.22          # arrow end beside a bond midpoint (ChemDraw starts bond arrows on the bond itself)
MIN_CLEARANCE = 0.3         # keep arrow ends this far from other atoms
MIN_SWEEP, MAX_SWEEP = 90.0, 150.0   # arc sweep range (degrees): clearly visible hooks at ACS size


class AtomRef(NamedTuple):
    """Fragment-qualified atom reference (``AtomRef(1, 0)`` == ``"1:0"``)."""
    fragment: int
    atom: int


# ---------------------------------------------------------------------------
# small vector helpers (page coordinates: points, y down)
# ---------------------------------------------------------------------------

def _sub(a: Point, b: Point) -> Point:
    return (a[0] - b[0], a[1] - b[1])


def _add(a: Point, b: Point) -> Point:
    return (a[0] + b[0], a[1] + b[1])


def _mul(a: Point, s: float) -> Point:
    return (a[0] * s, a[1] * s)


def _dot(a: Point, b: Point) -> float:
    return a[0] * b[0] + a[1] * b[1]


def _norm(a: Point) -> float:
    return math.hypot(a[0], a[1])


def _dist(a: Point, b: Point) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def _unit(a: Point, fallback: Point = (1.0, 0.0)) -> Point:
    n = _norm(a)
    if n < 1e-9:
        return fallback
    return (a[0] / n, a[1] / n)


def _mid(a: Point, b: Point) -> Point:
    return ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)


def _perp(u: Point) -> Point:
    """Perpendicular of a unit vector (rotated 90 degrees clockwise on screen)."""
    return (-u[1], u[0])


def _mean(points: Sequence[Point]) -> Point:
    if not points:
        return (0.0, 0.0)
    return (sum(p[0] for p in points) / len(points), sum(p[1] for p in points) / len(points))


def _clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def _box_dist(p: Point, box: Box) -> float:
    """Distance from a point to a box (0 inside)."""
    dx = max(box.x0 - p[0], 0.0, p[0] - box.x1)
    dy = max(box.y0 - p[1], 0.0, p[1] - box.y1)
    return math.hypot(dx, dy)


def _clearance(p: Point, ob: Obstacle) -> float:
    if isinstance(ob, Box):
        return _box_dist(p, ob)
    centre, r = ob
    return _dist(p, centre) - r


def _exit_distance(p: Point, d: Point, box: Optional[Box]) -> float:
    """Length of the ray p + t*d inside ``box`` (0 when p is outside)."""
    if box is None or _box_dist(p, box) > 1e-9:
        return 0.0
    t = float("inf")
    if d[0] > 1e-9:
        t = min(t, (box.x1 - p[0]) / d[0])
    elif d[0] < -1e-9:
        t = min(t, (box.x0 - p[0]) / d[0])
    if d[1] > 1e-9:
        t = min(t, (box.y1 - p[1]) / d[1])
    elif d[1] < -1e-9:
        t = min(t, (box.y0 - p[1]) / d[1])
    return 0.0 if t == float("inf") else max(t, 0.0)


# ---------------------------------------------------------------------------
# arc geometry (mirrors CDXMLDocument.add_curved_arrow)
# ---------------------------------------------------------------------------

def _arc_params(tail: Point, head: Point, angle: float) -> Tuple[Point, float, Point, float]:
    """(centre, radius, bulge direction, sagitta) of the arc add_curved_arrow draws for ``angle``."""
    chord = _dist(tail, head)
    if chord < 1e-9:
        raise ValueError("curved arrow needs distinct tail and head")
    a = math.radians(_clamp(abs(angle), 1.0, 359.0))
    r = chord / (2.0 * math.sin(a / 2.0))
    m = _mid(tail, head)
    u = _unit(_sub(head, tail))
    p = (u[1], -u[0]) if angle > 0 else (-u[1], u[0])      # centre side (same formula as add_curved_arrow)
    h = math.sqrt(max(r * r - (chord / 2.0) ** 2, 0.0))
    centre = _add(m, _mul(p, h))
    bulge_dir = (-p[0], -p[1])
    sagitta = r - h if a <= math.pi else r + h
    return centre, r, bulge_dir, sagitta


def arc_apex(tail: Point, head: Point, angle: float) -> Point:
    """Apex (point of maximal bulge) of the arc ``add_curved_arrow(tail, head, angle)`` draws."""
    _centre, _r, bulge_dir, sagitta = _arc_params(tail, head, angle)
    return _add(_mid(tail, head), _mul(bulge_dir, sagitta))


def arc_points(tail: Point, head: Point, angle: float, n: int = 16) -> List[Point]:
    """``n`` points along the arc from tail to head (inclusive)."""
    centre, r, _bulge_dir, _ = _arc_params(tail, head, angle)
    t0 = math.atan2(tail[1] - centre[1], tail[0] - centre[0])
    t1 = math.atan2(head[1] - centre[1], head[0] - centre[0])
    d = t1 - t0
    while d <= -math.pi:
        d += 2 * math.pi
    while d > math.pi:
        d -= 2 * math.pi
    if abs(angle) > 180.0:          # major arc: go the long way round
        d = d - 2 * math.pi if d > 0 else d + 2 * math.pi
    n = max(2, n)
    return [(centre[0] + r * math.cos(t0 + d * k / (n - 1)), centre[1] + r * math.sin(t0 + d * k / (n - 1))) for k in range(n)]


# ---------------------------------------------------------------------------
# endpoints and arrow specs
# ---------------------------------------------------------------------------

_SOURCE_KINDS = ("lp", "atom", "bond", "h")
_TARGET_KINDS = ("atom", "bond", "new_bond")
_TAGS = {
    "lp": "lp", "lone_pair": "lp", "lonepair": "lp", "lone": "lp", "pair": "lp",
    "atom": "atom", "a": "atom",
    "h": "h", "hydrogen": "h", "ch": "h",
    "bond": "bond", "sigma": "bond", "pi": "bond",
    "new_bond": "new_bond", "new": "new_bond", "newbond": "new_bond",
}
_ATOM_TOKEN = r"\d+(?:\s*:\s*\d+)?"
_ATOM_RE = re.compile(r"^\s*(\d+)\s*(?::\s*(\d+))?\s*$")
_LP_RE = re.compile(r"^(?:lp|lone|lonepair|lone_pair|pair)\s*:\s*(%s)$" % _ATOM_TOKEN, re.I)
_H_RE = re.compile(r"^(?:h|ch|hydrogen)\s*:\s*(%s)$" % _ATOM_TOKEN, re.I)
_BOND_RE = re.compile(r"^(?:(bond|sigma|pi|new|new_bond|newbond)\s*:\s*)?(%s)\s*([-=#])\s*(%s)$" % (_ATOM_TOKEN, _ATOM_TOKEN), re.I)
_ATOM_ONLY_RE = re.compile(r"^(?:(?:atom|a)\s*:\s*)?(%s)$" % _ATOM_TOKEN, re.I)
_ARROW_SPLIT_RE = re.compile(r"\s*(-->|->|~>|=>|→|⇒|⟶)\s*")
_ELECTRON_SUFFIX_RE = re.compile(r"[\s:,(\[]\s*([12])\s*e[-−]?\s*[)\]]?\s*$", re.I)


def _atom_key(value: Any, default_fragment: int = 0) -> AtomKey:
    """Coerce ``5`` / ``"5"`` / ``"1:5"`` / ``AtomRef(1, 5)`` / ``(1, 5)`` / ``{"fragment": 1, "atom": 5}`` to (fragment, atom)."""
    if isinstance(value, AtomRef):
        return (int(value.fragment), int(value.atom))
    if isinstance(value, bool):
        raise ValueError("atom reference must be an index, 'fragment:atom' or (fragment, atom), not %r" % (value,))
    if isinstance(value, int):
        if value < 0:
            raise ValueError("atom index must be >= 0, got %d" % value)
        return (default_fragment, value)
    if isinstance(value, float) and value.is_integer():
        return _atom_key(int(value), default_fragment)
    if isinstance(value, str):
        m = _ATOM_RE.match(value)
        if not m:
            raise ValueError("bad atom reference %r (expected 'atom' or 'fragment:atom', 0-based)" % value)
        if m.group(2) is None:
            return (default_fragment, int(m.group(1)))
        return (int(m.group(1)), int(m.group(2)))
    if isinstance(value, dict):
        if "atom" not in value:
            raise ValueError("atom reference dict needs an 'atom' key: %r" % (value,))
        return (int(value.get("fragment", default_fragment)), int(value["atom"]))
    if isinstance(value, (tuple, list)) and len(value) == 2 and all(isinstance(v, int) and not isinstance(v, bool) for v in value):
        return (int(value[0]), int(value[1]))
    raise ValueError("bad atom reference %r" % (value,))


def _parse_endpoint(text: str, role: str) -> tuple:
    s = text.strip()
    if not s:
        raise ValueError("empty %s in arrow spec" % role)
    m = _LP_RE.match(s)
    if m:
        return ("lp", _atom_key(m.group(1)))
    m = _H_RE.match(s)
    if m:
        return ("h", _atom_key(m.group(1)))
    m = _BOND_RE.match(s)
    if m:
        prefix, a, _sep, b = m.groups()
        kind = "new_bond" if prefix and prefix.lower() in ("new", "new_bond", "newbond") else "bond"
        return (kind, _atom_key(a), _atom_key(b))
    m = _ATOM_ONLY_RE.match(s)
    if m:
        return ("atom", _atom_key(m.group(1)))
    raise ValueError("cannot parse arrow %s %r; use e.g. 'lp:5', '3-4', '3=4', 'h:3', '4' or '1:0-1:1' (fragment:atom)" % (role, text))


def _endpoint_from_dict(d: dict, role: str) -> tuple:
    for key in ("lone_pair", "lonepair", "lp"):
        if key in d:
            return ("lp", _atom_key(d[key]))
    for key in ("hydrogen", "h"):
        if key in d:
            return ("h", _atom_key(d[key]))
    if "new_bond" in d:
        a, b = d["new_bond"]
        return ("new_bond", _atom_key(a), _atom_key(b))
    if "bond" in d:
        a, b = d["bond"]
        return ("bond", _atom_key(a), _atom_key(b))
    if "atoms" in d and "kind" in d:
        kind = _TAGS.get(str(d["kind"]).lower())
        atoms = [_atom_key(a) for a in d["atoms"]]
        if kind is None:
            raise ValueError("unknown endpoint kind %r" % d["kind"])
        return (kind, *atoms)
    if "atom" in d:
        return ("atom", _atom_key(d))
    if "spec" in d:
        return _parse_endpoint(str(d["spec"]), role)
    raise ValueError("endpoint dict for %s must have one of lone_pair/lp, bond, new_bond, atom, hydrogen: %r" % (role, d))


def _normalize_endpoint(value: Any, role: str) -> tuple:
    """Canonical endpoint tuple: ("atom"|"lp"|"h", key) or ("bond"|"new_bond", key, key)."""
    if isinstance(value, str):
        ep = _parse_endpoint(value, role)
    elif isinstance(value, dict):
        ep = _endpoint_from_dict(value, role)
    elif isinstance(value, AtomRef) or (isinstance(value, int) and not isinstance(value, bool)):
        ep = ("atom", _atom_key(value))
    elif isinstance(value, (tuple, list)):
        items = list(value)
        tag = _TAGS.get(items[0].strip().lower()) if items and isinstance(items[0], str) else None
        if tag is not None:
            rest = items[1:]
            if len(rest) == 1 and isinstance(rest[0], (tuple, list)) and len(rest[0]) == 2 and tag in ("bond", "new_bond"):
                rest = list(rest[0])
            if tag in ("lp", "atom", "h"):
                if len(rest) != 1:
                    raise ValueError("%r endpoint takes one atom: %r" % (tag, value))
                ep = (tag, _atom_key(rest[0]))
            else:
                if len(rest) != 2:
                    raise ValueError("%r endpoint takes two atoms: %r" % (tag, value))
                ep = (tag, _atom_key(rest[0]), _atom_key(rest[1]))
        elif len(items) == 2:
            ep = ("bond", _atom_key(items[0]), _atom_key(items[1]))
        else:
            raise ValueError("cannot interpret %s %r" % (role, value))
    else:
        raise ValueError("cannot interpret %s %r" % (role, value))
    allowed = _SOURCE_KINDS if role == "source" else _TARGET_KINDS
    if ep[0] not in allowed:
        raise ValueError("%s cannot be %r (allowed: %s)" % (role, ep[0], ", ".join(allowed)))
    if ep[0] in ("bond", "new_bond") and ep[1] == ep[2]:
        raise ValueError("bond needs two different atoms: %r" % (value,))
    return ep


def _fmt_key(key: AtomKey, qualified: bool) -> str:
    return "%d:%d" % key if qualified else "%d" % key[1]


@dataclass
class ElectronArrow:
    """One electron-pushing arrow.

    ``source``: atom index (electrons of that atom), a bond ``(i, j)``, ``("lp", atom)`` for a lone pair or
    ``("h", atom)`` for the C-H bond of an implicit hydrogen.  ``target``: atom index, an existing bond
    ``(i, j)`` (electrons move into it) or ``("new_bond", i, j)``.  Atoms may be fragment-qualified
    (``"1:0"``, ``AtomRef(1, 0)``).  After construction ``source``/``target`` hold the canonical tuples
    ``("lp", (fragment, atom))`` / ``("new_bond", (f, i), (f, j))`` ...  A plain bond target is resolved to
    ``bond`` or ``new_bond`` when the molecule is known (``parse_arrow_spec(spec, mol)``, ``add_electron_arrows``).

    ``electrons``: 2 (full head) or 1 (fishhook).  ``bulge``: "auto" (away from the molecule's centroid) | "up" |
    "down".  ``angle``: optional explicit signed sweep in degrees (overrides the automatic choice).
    """
    source: Any
    target: Any
    electrons: int = 2
    bulge: str = "auto"
    label: Optional[str] = None
    angle: Optional[float] = None

    def __post_init__(self):
        self.source = _normalize_endpoint(self.source, "source")
        self.target = _normalize_endpoint(self.target, "target")
        try:
            e = int(self.electrons)
        except (TypeError, ValueError):
            raise ValueError("electrons must be 1 or 2, got %r" % (self.electrons,))
        if e not in (1, 2):
            raise ValueError("electrons must be 1 (fishhook) or 2 (full head), got %d" % e)
        self.electrons = e
        b = str(self.bulge or "auto").strip().lower()
        if b not in ("auto", "up", "down"):
            raise ValueError("bulge must be 'auto', 'up' or 'down', got %r" % (self.bulge,))
        self.bulge = b
        if self.angle is not None:
            self.angle = float(self.angle)
            if not 10.0 <= abs(self.angle) <= 180.0:
                raise ValueError("angle must be between 10 and 180 degrees in magnitude")
        if self.label is not None:
            self.label = str(self.label)

    # -- accessors ------------------------------------------------------------
    @property
    def source_kind(self) -> str:
        return self.source[0]

    @property
    def source_atoms(self) -> Tuple[AtomKey, ...]:
        return tuple(self.source[1:])

    @property
    def target_kind(self) -> str:
        return self.target[0]

    @property
    def target_atoms(self) -> Tuple[AtomKey, ...]:
        return tuple(self.target[1:])

    @property
    def fragments(self) -> List[int]:
        return sorted({k[0] for k in self.source_atoms + self.target_atoms})

    @property
    def arrow_kind(self) -> str:
        return "full" if self.electrons == 2 else "half_left"

    def to_spec(self) -> str:
        q = any(f != 0 for f in self.fragments)
        src = {"lp": "lp:%s", "h": "h:%s", "atom": "%s"}.get(self.source_kind)
        if src:
            s = src % _fmt_key(self.source_atoms[0], q)
        else:
            s = "%s-%s" % tuple(_fmt_key(k, q) for k in self.source_atoms)
        if self.target_kind == "atom":
            t = _fmt_key(self.target_atoms[0], q)
        else:
            t = "%s-%s" % tuple(_fmt_key(k, q) for k in self.target_atoms)
            if self.target_kind == "new_bond":
                t = "new:" + t
        return "%s->%s%s" % (s, t, ":1e" if self.electrons == 1 else "")

    def describe(self) -> dict:
        return {"spec": self.to_spec(), "source": {"kind": self.source_kind, "atoms": [list(k) for k in self.source_atoms]},
                "target": {"kind": self.target_kind, "atoms": [list(k) for k in self.target_atoms]},
                "electrons": self.electrons, "bulge": self.bulge, "label": self.label, "angle": self.angle}

    # -- validation / resolution -------------------------------------------------
    def validate(self, n_atoms: Sequence[int], has_bond: Optional[Callable[[int, int, int], bool]] = None) -> None:
        """Check fragment / atom indices against ``n_atoms`` (atom count per fragment) and that a source bond exists."""
        for key in self.source_atoms + self.target_atoms:
            f, a = key
            if f < 0 or f >= len(n_atoms):
                raise ValueError("arrow %s references fragment %d but the stage has %d structure(s)" % (self.to_spec(), f, len(n_atoms)))
            if a < 0 or a >= n_atoms[f]:
                raise ValueError("arrow %s references atom %d of fragment %d which has %d atoms (0-based indices)" % (
                    self.to_spec(), a, f, n_atoms[f]))
        if has_bond is not None and self.source_kind == "bond":
            (f1, a), (f2, b) = self.source_atoms
            if f1 != f2 or not has_bond(f1, a, b):
                raise ValueError("arrow %s: source bond %s-%s does not exist" % (self.to_spec(), _fmt_key((f1, a), True), _fmt_key((f2, b), True)))
        if self.source_kind == "bond" and self.target_kind in ("bond", "new_bond") and set(self.source_atoms) == set(self.target_atoms):
            q = any(f != 0 for f in self.fragments)
            a, b = (_fmt_key(k, q) for k in self.source_atoms)
            raise ValueError("arrow %s: source and target are the same bond %s-%s; the electrons that turn it into a double bond "
                             "must come from elsewhere: a lone pair ('lp:%s->%s=%s'), a C-H bond ('h:%s->%s=%s') or a "
                             "neighbouring bond ('<neighbour>-%s->%s=%s')" % (self.to_spec(), a, b, a, a, b, a, a, b, a, a, b))

    def resolved(self, has_bond: Callable[[int, int, int], bool]) -> "ElectronArrow":
        """Copy with a bond target classified as ``bond`` (exists) or ``new_bond`` (does not exist yet)."""
        if self.target_kind not in ("bond", "new_bond"):
            return self
        (f1, a), (f2, b) = self.target_atoms
        exists = f1 == f2 and has_bond(f1, a, b)
        kind = "bond" if exists else "new_bond"
        if kind == self.target_kind:
            return self
        return replace(self, target=(kind, (f1, a), (f2, b)))


def _mol_list(mol: Any) -> List[Any]:
    if mol is None:
        return []
    if hasattr(mol, "GetNumAtoms"):
        return [mol]
    return list(mol)


def _has_bond_in_mols(mols: Sequence[Any]) -> Callable[[int, int, int], bool]:
    def has_bond(f: int, a: int, b: int) -> bool:
        m = mols[f]
        return m is not None and m.GetBondBetweenAtoms(a, b) is not None
    return has_bond


def parse_arrow_spec(spec: Union[str, dict, ElectronArrow, Sequence[Any]], mol: Any = None) -> ElectronArrow:
    """Parse a compact string (``"lp:5->5-8"``, ``"3-4->4"``, ``"3=4->4-7:1e"``, ``"lp:0:0->0:0-1:0"``), a dict
    (``{"from": "lp:5", "to": "5-8", "electrons": 1, "bulge": "up", "label": "..."}``) or a ``(source, target[, electrons])``
    tuple.  With ``mol`` (an RDKit mol or a list of them for multi-structure stages) indices are validated and bond
    targets are classified as existing bond / new bond."""
    if isinstance(spec, ElectronArrow):
        arrow = spec
    elif isinstance(spec, str):
        arrow = _parse_spec_string(spec)
    elif isinstance(spec, dict):
        arrow = _parse_spec_dict(spec)
    elif isinstance(spec, (tuple, list)) and 2 <= len(spec) <= 4 and not (isinstance(spec[0], str) and spec[0].lower() in _TAGS):
        arrow = ElectronArrow(*spec)
    else:
        raise ValueError("arrow spec must be a string like 'lp:5->5-8', a dict with from/to or an ElectronArrow, got %r" % (spec,))
    mols = _mol_list(mol)
    if mols:
        arrow.validate([m.GetNumAtoms() for m in mols], _has_bond_in_mols(mols))
        arrow = arrow.resolved(_has_bond_in_mols(mols))
    return arrow


def _parse_spec_string(spec: str) -> ElectronArrow:
    s = spec.strip()
    electrons = 2
    m = _ELECTRON_SUFFIX_RE.search(s)
    if m:
        electrons = int(m.group(1))
        s = s[:m.start()].strip()
    parts = _ARROW_SPLIT_RE.split(s, maxsplit=1)
    if len(parts) != 3 or not parts[0].strip() or not parts[2].strip():
        raise ValueError("arrow spec %r needs the form 'source->target' (e.g. 'lp:5->5-8', '3-4->4')" % spec)
    src, sep, tgt = parts
    if sep == "~>":
        electrons = 1
    return ElectronArrow(src, tgt, electrons)


def _parse_spec_dict(d: dict) -> ElectronArrow:
    base: Optional[ElectronArrow] = None
    if d.get("spec"):
        base = _parse_spec_string(str(d["spec"]))
    source = d.get("from", d.get("source", d.get("tail")))
    target = d.get("to", d.get("target", d.get("head")))
    if base is not None:
        source = base.source if source is None else source
        target = base.target if target is None else target
    if source is None or target is None:
        raise ValueError("arrow dict needs 'from' and 'to' (or 'spec'): %r" % (d,))
    electrons = d.get("electrons", d.get("e"))
    if electrons is None:
        electrons = 1 if d.get("fishhook") or d.get("single") or d.get("radical") else (base.electrons if base else 2)
    return ElectronArrow(source, target, electrons, d.get("bulge") or "auto", d.get("label"), d.get("angle"))


# ---------------------------------------------------------------------------
# fragment geometry (page coordinates)
# ---------------------------------------------------------------------------

@dataclass
class FragmentGeometry:
    """Atom positions and connectivity of a placed fragment, in page coordinates."""
    ref: ObjRef
    index: int = 0
    positions: List[Point] = field(default_factory=list)
    elements: List[int] = field(default_factory=list)
    hydrogens: List[int] = field(default_factory=list)       # hydrogens shown in the atom label
    labelled: List[bool] = field(default_factory=list)       # True when ChemDraw draws a text label for the atom
    charges: List[int] = field(default_factory=list)
    neighbors: List[List[int]] = field(default_factory=list)
    bond_orders: Dict[Tuple[int, int], str] = field(default_factory=dict)
    label_boxes: Dict[int, Box] = field(default_factory=dict)  # estimated label box per labelled atom (when recorded)
    label_sides: Dict[int, str] = field(default_factory=dict)  # 'left' (hydrogens first, 'H3C'), 'right' ('OH') or 'center'

    @property
    def n_atoms(self) -> int:
        return len(self.positions)

    @property
    def centroid(self) -> Point:
        return _mean(self.positions)

    @property
    def box(self) -> Box:
        return self.ref.box

    def has_bond(self, i: int, j: int) -> bool:
        return (min(i, j), max(i, j)) in self.bond_orders

    def bond_order(self, i: int, j: int) -> Optional[str]:
        return self.bond_orders.get((min(i, j), max(i, j)))


def _fragment_element(doc: CDXMLDocument, ref: ObjRef) -> Optional[ET.Element]:
    prefix = '<fragment id="%d"' % ref.id
    for snippet in doc.objects:
        if snippet.startswith(prefix):
            try:
                return ET.fromstring(snippet)
            except ET.ParseError:
                return None
    return None


def fragment_geometry(doc: CDXMLDocument, ref: ObjRef, index: int = 0, mol: Any = None) -> FragmentGeometry:
    """Positions (``ref.extra["atom_positions"]`` when recorded, else the ``p`` attributes of the fragment's ``<n>``
    nodes in the generated XML) plus connectivity (from the XML ``<b>`` elements, or from ``mol``)."""
    if ref.kind != "fragment":
        raise ValueError("expected a fragment ObjRef, got kind %r (id %d)" % (ref.kind, ref.id))
    el = _fragment_element(doc, ref)
    atom_ids = ref.extra.get("atom_ids")
    nodes_by_id: Dict[str, ET.Element] = {}
    node_list: List[ET.Element] = []
    if el is not None:
        node_list = el.findall("n")
        nodes_by_id = {nd.get("id", ""): nd for nd in node_list}
    n = int(ref.extra.get("atoms") or (mol.GetNumAtoms() if mol is not None else len(node_list)))
    ordered: List[ET.Element] = []
    if atom_ids and len(atom_ids) == n and nodes_by_id and all(str(i) in nodes_by_id for i in atom_ids):
        ordered = [nodes_by_id[str(i)] for i in atom_ids]
    elif len(node_list) == n:
        ordered = list(node_list)
    recorded = ref.extra.get("atom_positions")
    positions: List[Point] = []
    if recorded and len(recorded) == n:
        positions = [(float(p[0]), float(p[1])) for p in recorded]
    elif ordered:
        for nd in ordered:
            p = (nd.get("p") or "0 0").split()
            positions.append((float(p[0]), float(p[1])))
    elif mol is not None and mol.GetNumConformers():
        # last resort: scale the conformer like cdxml._mol_points and centre it in the fragment box
        from .cdxml import _mol_points
        pts = _mol_points(mol, doc.bond_length)
        raw = union_boxes([Box(x, y, x, y) for x, y in pts])
        dx, dy = ref.box.cx - raw.cx, ref.box.cy - raw.cy
        positions = [(x + dx, y + dy) for x, y in pts]
    else:
        raise ValueError("fragment %d: no atom positions recorded and its XML was not found in the document" % ref.id)
    geom = FragmentGeometry(ref, index, positions)
    geom.neighbors = [[] for _ in range(n)]
    if mol is not None and mol.GetNumAtoms() == n:
        for atom in mol.GetAtoms():
            z = atom.GetAtomicNum()
            geom.elements.append(z)
            geom.hydrogens.append(atom.GetTotalNumHs())
            geom.charges.append(atom.GetFormalCharge())
            geom.labelled.append(bool(z != 6 or atom.GetDegree() == 0 or atom.GetFormalCharge() or atom.GetIsotope()
                                      or atom.GetNumRadicalElectrons()))
    elif ordered:
        for nd in ordered:
            geom.elements.append(int(float(nd.get("Element", "6"))))
            nh = nd.get("NumHydrogens")
            geom.hydrogens.append(int(float(nh)) if nh is not None else 0)
            geom.labelled.append(nh is not None or nd.get("NodeType") is not None)
            geom.charges.append(int(float(nd.get("Charge", "0"))))
    else:
        geom.elements = [6] * n
        geom.hydrogens = [0] * n
        geom.labelled = [False] * n
        geom.charges = [0] * n
    if el is not None and ordered:
        id_to_idx = {nd.get("id", ""): i for i, nd in enumerate(ordered)}
        for b in el.findall("b"):
            i, j = id_to_idx.get(b.get("B", "")), id_to_idx.get(b.get("E", ""))
            if i is None or j is None or i == j:
                continue
            geom.bond_orders[(min(i, j), max(i, j))] = b.get("Order", "1")
            geom.neighbors[i].append(j)
            geom.neighbors[j].append(i)
    elif mol is not None:
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            geom.bond_orders[(min(i, j), max(i, j))] = ("%g" % b.GetBondTypeAsDouble())
            geom.neighbors[i].append(j)
            geom.neighbors[j].append(i)
    for lab in ref.extra.get("labels") or []:
        try:
            i = int(lab["atom"])
            geom.label_boxes[i] = Box(*[float(v) for v in lab["box"]])
            geom.label_sides[i] = str(lab.get("side") or "center")
            if 0 <= i < n:
                geom.labelled[i] = True
        except (KeyError, TypeError, ValueError):
            continue
    return geom


def fragment_atom_positions(doc: CDXMLDocument, ref: ObjRef, mol: Any = None) -> List[Point]:
    """Page coordinates (points, y down) of the fragment's atoms, indexed like the RDKit molecule."""
    return list(fragment_geometry(doc, ref, 0, mol).positions)


# ---------------------------------------------------------------------------
# arrow planning
# ---------------------------------------------------------------------------

class _Stage:
    """All fragments of one mechanism stage plus obstacles that arrows should avoid."""

    def __init__(self, frags: Sequence[FragmentGeometry], bond_length: float, obstacles: Optional[Sequence[Any]] = None):
        self.frags = list(frags)
        self.L = bond_length
        self.obstacles: List[Obstacle] = []
        for ob in obstacles or []:
            self.add_obstacle(ob)

    def add_obstacle(self, ob: Any) -> None:
        if isinstance(ob, ObjRef):
            if ob.kind == "text" and str(ob.extra.get("text", "")).strip() == "+":
                self.obstacles.append(((ob.box.cx, ob.box.cy), 0.4 * self.L))
                return
            if ob.kind == "arrow" and ob.extra.get("tail") and ob.extra.get("head") and not ob.extra.get("mechanism"):
                t, h = ob.extra["tail"], ob.extra["head"]
                self.obstacles.append(Box(min(t[0], h[0]), min(t[1], h[1]), max(t[0], h[0]), max(t[1], h[1])).padded(0.1 * self.L))
                return
            ob = ob.box
        if isinstance(ob, Box):
            self.obstacles.append(ob)
        elif isinstance(ob, (tuple, list)) and len(ob) == 3:
            self.obstacles.append(((float(ob[0]), float(ob[1])), float(ob[2])))
        elif isinstance(ob, (tuple, list)) and len(ob) == 2:
            self.obstacles.append(((float(ob[0]), float(ob[1])), 0.25 * self.L))
        else:
            raise ValueError("obstacle must be an ObjRef, Box or (x, y, radius), got %r" % (ob,))

    def frag(self, f: int) -> FragmentGeometry:
        if f < 0 or f >= len(self.frags):
            raise ValueError("fragment index %d out of range (stage has %d structure(s))" % (f, len(self.frags)))
        return self.frags[f]

    def point(self, key: AtomKey) -> Point:
        fr = self.frag(key[0])
        if key[1] < 0 or key[1] >= fr.n_atoms:
            raise ValueError("atom %d out of range for fragment %d (%d atoms)" % (key[1], key[0], fr.n_atoms))
        return fr.positions[key[1]]

    def label_side(self, key: AtomKey) -> str:
        return self.frag(key[0]).label_sides.get(key[1], "center")

    def label_box(self, key: AtomKey) -> Optional[Box]:
        return self.frag(key[0]).label_boxes.get(key[1])

    def n_atoms(self) -> List[int]:
        return [fr.n_atoms for fr in self.frags]

    def has_bond(self, f: int, a: int, b: int) -> bool:
        return self.frag(f).has_bond(a, b)

    def atoms_except(self, keys: Sequence[AtomKey]) -> List[Point]:
        skip = set(keys)
        return [fr.positions[i] for fr in self.frags for i in range(fr.n_atoms) if (fr.index, i) not in skip]

    def label_boxes_except(self, keys: Sequence[AtomKey]) -> List[Box]:
        skip = set(keys)
        return [b for fr in self.frags for i, b in fr.label_boxes.items() if (fr.index, i) not in skip]

    def centroid_of(self, frag_ids: Sequence[int]) -> Point:
        return _mean([p for f in frag_ids for p in self.frag(f).positions])


@dataclass
class _Plan:
    tail: Point
    head: Point
    angle: float
    apex: Point
    centroid: Point
    samples: List[Point]
    source_kind: str
    target_kind: str


def _hydrogen_direction(frag: FragmentGeometry, i: int, bond_dirs: Sequence[Point]) -> Optional[Point]:
    """Where ChemDraw writes the H's of a labelled atom (left/right/up/down), or None when there are none."""
    if not frag.labelled[i] or frag.hydrogens[i] <= 0:
        return None
    if not bond_dirs:
        return (1.0, 0.0)
    m = (-sum(b[0] for b in bond_dirs), -sum(b[1] for b in bond_dirs))
    if _norm(m) < 0.2:
        return None
    u = _unit(m)
    if abs(u[0]) >= abs(u[1]):
        return (1.0 if u[0] > 0 else -1.0, 0.0)
    return (0.0, 1.0 if u[1] > 0 else -1.0)


def _charge_corner(p: Point, box: Box) -> Point:
    """Where ChemDraw puts the charge superscript of a label: top right of the text (top left for labels drawn to
    the left of the atom)."""
    extends_right = (box.x1 - p[0]) >= (p[0] - box.x0) - 1e-6
    return (box.x1, box.y0) if extends_right else (box.x0, box.y0)


def _atom_endpoint(stage: _Stage, key: AtomKey, other: Optional[Point], centroid: Point, base: Optional[float] = None,
                   role: str = "tail") -> Point:
    """A point beside the atom (tails ~0.5 L out, heads ~0.35 L / 0.42 L for labelled atoms, always outside the label
    box when known), on a side free of bonds, H labels and the charge superscript.  Tails prefer the side facing
    ``other`` (the opposite end of the arrow); heads prefer the far side so the hook stays visible."""
    L = stage.L
    frag = stage.frag(key[0])
    i = key[1]
    p = stage.point(key)
    if base is None:
        if role == "tail":
            base = LONE_PAIR_OFFSET * L
        else:
            base = (LABEL_OFFSET if frag.labelled[i] else ATOM_OFFSET) * L
    lbox = stage.label_box(key)
    charge_dir = _unit(_sub(_charge_corner(p, lbox), p), (1.0, -0.5)) if (lbox is not None and frag.charges[i]) else None
    toward_w, apart = (0.8, 0.6 * L) if role == "tail" else (0.3, 1.0 * L)
    bond_dirs = [_unit(_sub(frag.positions[j], p)) for j in frag.neighbors[i] if _dist(frag.positions[j], p) > 1e-6]
    h_dir = _hydrogen_direction(frag, i, bond_dirs) if lbox is None else None
    charge_pt = (0.0, 0.0)
    if charge_dir is not None:
        cc = _charge_corner(p, lbox)
        charge_pt = (cc[0] - (2.5 if cc[0] > p[0] else -2.5), cc[1] + 2.0)
    hs_dir = None                         # side of the label that carries the hydrogens of "H3C" / "OH"
    if lbox is not None and frag.hydrogens[i] > 0 and stage.label_side(key) in ("left", "right"):
        hs_dir = (-1.0, 0.0) if stage.label_side(key) == "left" else (1.0, 0.0)
    if not bond_dirs and other is not None and lbox is None:
        h_dir = None                      # isolated atom: nothing but the partner decides the side
    others = stage.atoms_except([key])
    other_boxes = stage.label_boxes_except([key])
    best_pt: Point = _add(p, (base, 0.0))
    best = -1e9
    for k in range(72):
        th = 2.0 * math.pi * k / 72.0
        d = (math.cos(th), math.sin(th))
        offset = _clamp(max(base, _exit_distance(p, d, lbox) + 0.06 * L), base, MAX_LABEL_OFFSET * L)
        pt = _add(p, _mul(d, offset))
        score = -0.6 * (offset - base) / L
        for b in bond_dirs:
            c = _dot(d, b)
            if c > 0:
                score -= 1.5 * c * c
        if h_dir is not None:
            c = _dot(d, h_dir)
            if c > 0:
                score -= 1.0 * c * c
        if charge_dir is not None:
            c = _dot(d, charge_dir)
            if c > 0:
                score -= 1.2 * c * c                    # do not cover the charge superscript
            cd = _dist(pt, charge_pt)
            if cd < 0.65 * L:
                score -= 4.0 * (0.65 * L - cd) / L      # ... nor touch it
        if hs_dir is not None:
            c = _dot(d, hs_dir)
            if c > 0:
                score -= 2.5 * c * c                    # ends point at the atom's symbol, not at its hydrogens ("H3C")
        if other is not None:
            v = _sub(other, p)
            nv = _norm(v)
            if nv > 1e-6:
                score += toward_w * _dot(d, v) / nv
            dd = _dist(pt, other)
            if dd < apart:
                score -= 2.0 * (apart - dd) / L         # keep tail and head apart (visible arc)
        v = _sub(p, centroid)
        nv = _norm(v)
        if nv > 0.05 * L:
            score += 0.25 * _dot(d, v) / nv
        score += 0.05 * (-d[1])                          # tie-break: up
        if others:
            dmin = min(_dist(pt, o) for o in others)
            if dmin < 0.45 * L:
                score -= 2.0 * (0.45 * L - dmin) / L
        for ob in other_boxes:
            c = _box_dist(pt, ob)
            if c < 0.1 * L:
                score -= 1.5 * (0.1 * L - c) / L + 0.5
        if score > best:
            best, best_pt = score, pt
    return best_pt


def _bond_endpoint(stage: _Stage, a: AtomKey, b: AtomKey, other: Optional[Point], centroid: Point, offset: Optional[float] = None) -> Point:
    """Bond midpoint offset ~0.25 L perpendicular to the bond, on the side away from the centroid / less crowded,
    moved along the bond if a label of one of its atoms covers it."""
    L = stage.L
    offset = BOND_OFFSET * L if offset is None else offset
    pa, pb = stage.point(a), stage.point(b)
    m = _mid(pa, pb)
    u = _unit(_sub(pb, pa))
    others = stage.atoms_except([a, b])
    other_boxes = stage.label_boxes_except([a, b])
    stage_centroid = stage.centroid_of(range(len(stage.frags)))
    best_pt = m
    best = -1e9
    for n in (_perp(u), _mul(_perp(u), -1.0)):
        pt = _add(m, _mul(n, offset))
        s = 0.0
        v = _sub(m, centroid)
        nv = _norm(v)
        if nv > 0.05 * L:
            s += 1.0 * _dot(n, v) / nv
        w = _sub(m, stage_centroid)
        nw = _norm(w)
        if nw > 0.05 * L:
            s += 0.3 * _dot(n, w) / nw
        if other is not None:
            v = _sub(other, m)
            nv = _norm(v)
            if nv > 1e-6:
                s += 0.5 * _dot(n, v) / nv
        if others:
            dmin = min(_dist(pt, o) for o in others)
            if dmin < 0.6 * L:
                s -= 2.5 * (0.6 * L - dmin) / L
        for ob in other_boxes:
            c = _box_dist(pt, ob)
            if c < 0.1 * L:
                s -= 2.0
        s += 0.05 * (-n[1])
        if s > best:
            best, best_pt = s, pt
    # slide along the bond out of the labels of its own atoms ("Cl", "OH")
    for key, direction in ((a, u), (b, _mul(u, -1.0))):
        box = stage.label_box(key)
        if box is not None and _box_dist(best_pt, box) < 1e-9:
            t = _exit_distance(best_pt, direction, box) + 0.05 * L
            best_pt = _add(best_pt, _mul(direction, min(t, 0.3 * L)))
    return best_pt


def _new_bond_endpoint(stage: _Stage, a: AtomKey, b: AtomKey, source_atoms: Sequence[AtomKey], tail: Point, centroid: Point) -> Point:
    """Head for a bond that is being formed: close to the atom being attacked, on the line towards the
    nucleophilic atom (outside its label), nudged slightly sideways away from the centroid."""
    L = stage.L
    pa, pb = stage.point(a), stage.point(b)
    src = set(source_atoms)
    if a in src and b not in src:
        s, e, ekey = pa, pb, b
    elif b in src and a not in src:
        s, e, ekey = pb, pa, a
    else:
        m = _mid(pa, pb)
        u = _unit(_sub(pb, pa))
        n = _perp(u)
        v = _sub(m, centroid)
        if _norm(v) > 0.05 * L and _dot(n, v) < 0:
            n = _mul(n, -1.0)
        elif _norm(v) <= 0.05 * L and n[1] > 0:
            n = _mul(n, -1.0)
        return _add(m, _mul(n, 0.2 * L))
    u = _unit(_sub(s, e))
    lbox = stage.label_box(ekey)
    frag_e = stage.frag(ekey[0])
    charge_pt = None
    if lbox is not None and frag_e.charges[ekey[1]]:
        cc = _charge_corner(e, lbox)
        charge_pt = (cc[0] - (2.5 if cc[0] > e[0] else -2.5), cc[1] + 2.0)
    hs_dir = None
    if lbox is not None and frag_e.hydrogens[ekey[1]] > 0 and stage.label_side(ekey) in ("left", "right"):
        hs_dir = (-1.0, 0.0) if stage.label_side(ekey) == "left" else (1.0, 0.0)
    base_off = _clamp(0.5 * _dist(s, e), 0.4 * L, MAX_ATOM_OFFSET * L)
    if charge_pt is not None or hs_dir is not None:
        # swing the approach direction (up to 75 degrees either way) away from the charge sign and the hydrogens
        best_u, best_sc = u, -1e9
        for k in range(-15, 16):
            th = math.radians(5.0 * k)
            cu = (u[0] * math.cos(th) - u[1] * math.sin(th), u[0] * math.sin(th) + u[1] * math.cos(th))
            o = max(base_off, min(_exit_distance(e, cu, lbox) + 0.06 * L, MAX_LABEL_OFFSET * L))
            cand = _add(e, _mul(cu, o))
            sc = -abs(k) * 0.02
            if charge_pt is not None:
                cd = _dist(cand, charge_pt)
                if cd < 0.7 * L:
                    sc -= 4.0 * (0.7 * L - cd) / L
            if hs_dir is not None:
                c = _dot(cu, hs_dir)
                if c > 0:
                    sc -= 2.5 * c * c
            if sc > best_sc:
                best_sc, best_u = sc, cu
        u = best_u
    off = max(base_off, min(_exit_distance(e, u, lbox) + 0.06 * L, MAX_LABEL_OFFSET * L))
    pt = _add(e, _mul(u, off))
    n = _perp(u)
    v = _sub(pt, centroid)
    if _norm(v) > 0.05 * L:
        if _dot(n, v) < 0:
            n = _mul(n, -1.0)
    elif n[1] > 0:
        n = _mul(n, -1.0)
    side = min(0.12 * L, math.sqrt(max((MAX_ATOM_OFFSET * L) ** 2 - off * off, 0.0)))
    return _add(pt, _mul(n, side))


def _push_clear(pt: Point, other: Point, obstacles: Sequence[Obstacle], min_d: float, max_shift: float,
                anchor: Optional[Point] = None, max_anchor_dist: float = float("inf")) -> Point:
    """Move ``pt`` towards ``other`` until it is at least ``min_d`` from every obstacle (at most ``max_shift``, and
    never further than ``max_anchor_dist`` from ``anchor``)."""
    if not obstacles:
        return pt
    u = _unit(_sub(other, pt), (0.0, 0.0))
    if u == (0.0, 0.0):
        return pt
    step = max(min_d / 10.0, 1e-3)
    moved = 0.0
    while moved < max_shift:
        if all(_clearance(pt, ob) >= min_d for ob in obstacles):
            break
        nxt = _add(pt, _mul(u, step))
        if anchor is not None and _dist(nxt, anchor) > max_anchor_dist:
            break
        pt = nxt
        moved += step
    return pt


def _end_obstacles(stage: _Stage, own: Sequence[AtomKey]) -> List[Obstacle]:
    """Other atoms (as discs) and their label boxes - what an arrow end must stay clear of."""
    L = stage.L
    obs: List[Obstacle] = [(p, (MIN_CLEARANCE - 0.05) * L) for p in stage.atoms_except(own)]
    obs += [b.padded(0.05 * L - MIN_CLEARANCE * L + 0.05 * L) for b in stage.label_boxes_except(own)]
    return obs


def _choose_sweep(stage: _Stage, tail: Point, head: Point, centroid: Point, bulge: str, frag_ids: Sequence[int],
                  exclude: Sequence[AtomKey]) -> float:
    """Signed sweep (90-150 degrees): bulging away from the centroid (or as forced by ``bulge``), keeping clear of
    atoms, labels, '+' signs, reaction arrows, earlier arrows and the fragment box centres."""
    L = stage.L
    chord = _dist(tail, head)
    base = _clamp(135.0 - 30.0 * (chord / L - 1.0), MIN_SWEEP, MAX_SWEEP)
    u = _unit(_sub(head, tail))
    mid = _mid(tail, head)
    centres = [(stage.frag(f).box.cx, stage.frag(f).box.cy) for f in frag_ids]
    obstacles: List[Obstacle] = [(p, 0.3 * L) for p in stage.atoms_except(exclude)]
    obstacles += list(stage.label_boxes_except([]))          # the arc must not cut through any label, own ones included
    obstacles += list(stage.obstacles)
    best_angle = base
    best = -1e9
    for sweep in sorted({base, min(MAX_SWEEP, base + 20.0), max(MIN_SWEEP, base - 20.0)}):
        for sign in (1.0, -1.0):
            a = sign * sweep
            pts = arc_points(tail, head, a, 24)
            inner = [p for p in pts if _dist(p, tail) > 0.35 * L and _dist(p, head) > 0.35 * L] or pts[1:-1]
            bd = _perp(u) if sign > 0 else _mul(_perp(u), -1.0)
            score = -0.3 * abs(sweep - base) / 20.0
            v = _sub(mid, centroid)
            nv = _norm(v)
            if nv > 0.05 * L:
                score += 1.0 * _dot(bd, v) / nv * min(1.0, nv / (0.75 * L))
            score += 0.3 * (-bd[1])                       # tie-break: up
            clear = min([_clearance(p, ob) for ob in obstacles for p in inner] or [1e9])
            if clear < 0.3 * L:
                score -= 3.0 * (0.3 * L - clear) / L
            if any(_dist(p, c) < 0.45 * L for c in centres for p in inner):
                score -= 3.0
            if bulge == "up":
                score += 2.0 if bd[1] < -0.2 else (-2.0 if bd[1] > 0.2 else 0.0)
            elif bulge == "down":
                score += 2.0 if bd[1] > 0.2 else (-2.0 if bd[1] < -0.2 else 0.0)
            if score > best:
                best, best_angle = score, a
    return best_angle


def _plan_arrow(stage: _Stage, arrow: ElectronArrow) -> _Plan:
    L = stage.L
    arrow.validate(stage.n_atoms(), stage.has_bond)
    arrow = arrow.resolved(stage.has_bond)
    src_kind, src_atoms = arrow.source_kind, arrow.source_atoms
    tgt_kind, tgt_atoms = arrow.target_kind, arrow.target_atoms
    frag_ids = arrow.fragments
    centroid = stage.centroid_of(frag_ids)
    src_anchor = stage.point(src_atoms[0]) if len(src_atoms) == 1 else _mid(stage.point(src_atoms[0]), stage.point(src_atoms[1]))
    tgt_anchor = stage.point(tgt_atoms[0]) if len(tgt_atoms) == 1 else _mid(stage.point(tgt_atoms[0]), stage.point(tgt_atoms[1]))
    # tail
    if src_kind in ("lp", "atom"):
        tail = _atom_endpoint(stage, src_atoms[0], tgt_anchor, centroid, role="tail")
    elif src_kind == "h":
        tail = _atom_endpoint(stage, src_atoms[0], tgt_anchor, centroid, base=0.45 * L, role="tail")
    else:
        tail = _bond_endpoint(stage, src_atoms[0], src_atoms[1], tgt_anchor, centroid)
    # head
    if tgt_kind == "atom":
        head = _atom_endpoint(stage, tgt_atoms[0], tail, centroid, role="head")
    elif tgt_kind == "bond":
        head = _bond_endpoint(stage, tgt_atoms[0], tgt_atoms[1], tail, centroid)
    else:
        head = _new_bond_endpoint(stage, tgt_atoms[0], tgt_atoms[1], src_atoms, tail, centroid)
    # keep the ends clear of other atoms and their labels (without drifting away from the referenced atom)
    tail = _push_clear(tail, head, _end_obstacles(stage, src_atoms), 0.05 * L, 0.15 * L,
                       src_anchor if len(src_atoms) == 1 else None, max(MAX_ATOM_OFFSET * L, _dist(tail, src_anchor)))
    head = _push_clear(head, tail, _end_obstacles(stage, tgt_atoms), 0.05 * L, 0.15 * L,
                       tgt_anchor if len(tgt_atoms) == 1 else None, max(MAX_ATOM_OFFSET * L, _dist(head, tgt_anchor)))
    if _dist(tail, head) < 0.3 * L:
        # degenerate (e.g. lone pair pushed onto its own bond): open the arrow up a little
        d = _unit(_sub(head, tail), _unit(_sub(head, centroid), (0.0, -1.0)))
        head = _add(tail, _mul(d, 0.5 * L))
    if arrow.angle is not None:
        angle = arrow.angle
    else:
        angle = _choose_sweep(stage, tail, head, centroid, arrow.bulge, frag_ids, list(src_atoms) + list(tgt_atoms))
    return _Plan(tail, head, angle, arc_apex(tail, head, angle), centroid, arc_points(tail, head, angle, 16), src_kind, tgt_kind)


def add_electron_arrows(doc: CDXMLDocument, ref: Union[ObjRef, Sequence[ObjRef]], arrows: Sequence[Any], mol: Any = None,
                        obstacles: Optional[Sequence[Any]] = None, color: Any = None, label_size: Optional[float] = None
                        ) -> List[ObjRef]:
    """Draw electron-pushing arrows on a placed fragment (or on several fragments of one stage: pass a list of
    ObjRefs; arrows then address atoms as ``fragment:atom``).

    ``arrows``: :class:`ElectronArrow` objects or specs accepted by :func:`parse_arrow_spec`.  ``mol``: the RDKit
    molecule(s) behind the fragment(s) (optional - connectivity is read back from the generated XML otherwise).
    ``obstacles``: boxes / ObjRefs / ``(x, y, radius)`` the arcs should keep away from (plus signs, captions).
    Returns the arrow ObjRefs; ``extra`` carries tail, head, apex, angle, electrons, source/target and spec.
    """
    refs = [ref] if isinstance(ref, ObjRef) else list(ref)
    if not refs:
        raise ValueError("add_electron_arrows needs at least one fragment")
    mols = _mol_list(mol)
    mols = mols + [None] * (len(refs) - len(mols))
    L = doc.bond_length
    geoms = [fragment_geometry(doc, r, i, mols[i]) for i, r in enumerate(refs)]
    stage = _Stage(geoms, L, obstacles)
    out: List[ObjRef] = []
    for spec in arrows:
        arrow = parse_arrow_spec(spec)
        plan = _plan_arrow(stage, arrow)
        aref = doc.add_curved_arrow(plan.tail[0], plan.tail[1], plan.head[0], plan.head[1], plan.angle, kind=arrow.arrow_kind, color=color)
        tight = union_boxes([Box(p[0], p[1], p[0], p[1]) for p in plan.samples])
        if tight is not None:
            aref.box = tight.padded(0.1 * L)
        aref.extra.update({
            "mechanism": True, "arrow_kind": arrow.arrow_kind, "electrons": arrow.electrons, "spec": arrow.to_spec(),
            "source": {"kind": plan.source_kind, "atoms": [list(k) for k in arrow.source_atoms]},
            "target": {"kind": plan.target_kind, "atoms": [list(k) for k in arrow.target_atoms]},
            "tail": [round(plan.tail[0], 2), round(plan.tail[1], 2)], "head": [round(plan.head[0], 2), round(plan.head[1], 2)],
            "apex": [round(plan.apex[0], 2), round(plan.apex[1], 2)], "angle": round(plan.angle, 2), "bulge": arrow.bulge,
            "centroid": [round(plan.centroid[0], 2), round(plan.centroid[1], 2)],
            "fragment_ids": [refs[f].id for f in arrow.fragments], "label": arrow.label,
        })
        stage_idx = refs[arrow.fragments[0]].extra.get("stage")
        if stage_idx is not None:
            aref.extra["stage"] = stage_idx
        if arrow.label:
            size = label_size or max(6.0, doc.caption_size - 2)
            bd = _unit(_sub(plan.apex, _mid(plan.tail, plan.head)), (0.0, -1.0))
            lp = _add(plan.apex, _mul(bd, 0.3 * L + 0.4 * size))
            tref = doc.add_text(arrow.label, lp[0], lp[1], size=size, anchor="center", color=color)
            aref.extra["label_id"] = tref.id
            stage.add_obstacle(tref)
        for p in plan.samples:
            stage.obstacles.append((p, 0.2 * L))
        out.append(aref)
    return out


# ---------------------------------------------------------------------------
# mechanism layout
# ---------------------------------------------------------------------------

@dataclass
class Stage:
    """One snapshot of a mechanism: one or more structures (drawn side by side with '+') with their arrows."""
    mols: List[Any]
    arrows: List[ElectronArrow] = field(default_factory=list)
    label: Optional[str] = None                      # caption below the stage
    arrow: str = "full"                              # reaction arrow leading to the next stage ("resonance", "equilibrium"...)
    above: List[str] = field(default_factory=list)   # text above that arrow
    below: List[str] = field(default_factory=list)   # text below that arrow
    labels: List[Optional[str]] = field(default_factory=list)   # per-structure captions
    colors: List[Any] = field(default_factory=list)
    label_atoms: List[Optional[bool]] = field(default_factory=list)   # per structure; None = automatic (tiny molecules)

    @classmethod
    def coerce(cls, obj: Any) -> "Stage":
        if isinstance(obj, Stage):
            return obj
        if isinstance(obj, dict):
            structures = None
            for key in ("structures", "structure", "smiles", "mols", "mol", "molecules", "molecule"):
                if obj.get(key) is not None:
                    structures = obj[key]
                    break
            if structures is None:
                raise ValueError("mechanism stage needs a 'structure' (SMILES / mol or a list of them): %r" % (obj,))
            arrows = obj.get("arrows", obj.get("electron_arrows")) or []
            label = obj.get("label", obj.get("caption"))
            arrow = obj.get("arrow", obj.get("arrow_to_next")) or "full"
            above = _lines(obj.get("above", obj.get("reagents")))
            below = _lines(obj.get("below", obj.get("conditions")))
            label_atoms_default = obj.get("label_atoms")
        else:
            structures, arrows, label, arrow, above, below, label_atoms_default = obj, [], None, "full", [], [], None
        if isinstance(structures, (str, dict)) or hasattr(structures, "GetNumAtoms"):
            structures = [structures]
        mols: List[Any] = []
        labels: List[Optional[str]] = []
        colors: List[Any] = []
        label_atoms: List[Optional[bool]] = []
        for item in structures:
            if isinstance(item, dict):
                src = item.get("structure", item.get("smiles", item.get("mol")))
                mols.append(_coerce_mol(src))
                labels.append(item.get("label"))
                colors.append(item.get("color"))
                la = item.get("label_atoms", label_atoms_default)
            else:
                mols.append(_coerce_mol(item))
                labels.append(None)
                colors.append(None)
                la = label_atoms_default
            label_atoms.append(None if la is None else bool(la))
        if not mols:
            raise ValueError("mechanism stage has no structures")
        parsed = [parse_arrow_spec(a, mols) for a in (arrows if isinstance(arrows, (list, tuple)) else [arrows])]
        return cls(mols, parsed, label, str(arrow), above, below, labels, colors, label_atoms)


def _lines(v: Any) -> List[str]:
    if v is None:
        return []
    if isinstance(v, str):
        return [v] if v else []
    return [str(x) for x in v if x is not None and str(x)]


def _coerce_mol(obj: Any):
    """SMILES (explicit ``[H]`` atoms kept) / any input chem.mol_from_input understands / RDKit Mol -> depiction-ready Mol."""
    from . import chem
    if obj is None:
        raise ValueError("structure is missing")
    if hasattr(obj, "GetNumAtoms"):
        if any(a.GetAtomicNum() == 1 for a in obj.GetAtoms()):
            return _prepare_keep_h(obj)
        return _ensure_finite_coords(chem.prepare_for_depiction(obj))
    s = str(obj).strip()
    if "[H]" in s and "\n" not in s:
        from rdkit import Chem
        params = Chem.SmilesParserParams()
        params.removeHs = False
        m = Chem.MolFromSmiles(s, params)
        if m is not None:
            return _prepare_keep_h(m)
    return _ensure_finite_coords(chem.mol_from_input(s))


def _ensure_finite_coords(mol):
    """Guarantee usable 2D coordinates: single atoms (``[H-]``, ``[H+]``) can come out of the depictor as NaN."""
    from rdkit.Chem import rdDepictor
    from rdkit.Geometry import Point3D

    def bad(m) -> List[int]:
        if m.GetNumConformers() == 0:
            return list(range(m.GetNumAtoms()))
        conf = m.GetConformer()
        return [i for i in range(m.GetNumAtoms())
                if not all(math.isfinite(v) for v in (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y))]

    n = mol.GetNumAtoms()
    if n == 0:
        raise ValueError("structure has no atoms")
    if not bad(mol):
        return mol
    if n == 1:
        if mol.GetNumConformers() == 0:
            from rdkit import Chem
            conf = Chem.Conformer(1)
            conf.Set3D(False)
            mol.AddConformer(conf, assignId=True)
        mol.GetConformer().SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
        return mol
    try:
        rdDepictor.SetPreferCoordGen(False)
        rdDepictor.Compute2DCoords(mol, clearConfs=True)
    except Exception:
        pass
    finally:
        try:
            rdDepictor.SetPreferCoordGen(True)
        except Exception:
            pass
    if bad(mol):
        raise ValueError("could not compute 2D coordinates for the structure (atoms %s)" % bad(mol))
    mol.GetConformer().Set3D(False)
    return mol


def _prepare_keep_h(mol):
    """Like chem.prepare_for_depiction but explicit hydrogens stay (so H atoms can take part in arrows)."""
    from rdkit import Chem
    m = Chem.Mol(mol)
    Chem.SanitizeMol(m)
    if m.GetNumConformers() == 0 or m.GetConformer().Is3D():
        try:
            from .chem import _compute_2d
            _compute_2d(m)
        except Exception:
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    Chem.Kekulize(m, clearAromaticFlags=True)
    try:
        Chem.WedgeMolBonds(m, m.GetConformer())
    except Exception:
        pass
    m.UpdatePropertyCache(strict=False)
    return _ensure_finite_coords(m)


def _orient_stage(stage: Stage) -> None:
    """Rotate structures (by multiples of 90 degrees - wedges stay valid) so that atoms forming bonds to another
    structure of the stage face that structure (structures are drawn left to right in list order)."""
    import numpy as np
    from rdkit.Geometry import Point3D
    n = len(stage.mols)
    if n < 2:
        return
    want: List[List[Tuple[int, float]]] = [[] for _ in range(n)]     # per fragment: (atom, +1 partner to the right / -1 left)
    for arrow in stage.arrows:
        if arrow.target_kind != "new_bond":
            continue
        (fa, ia), (fb, ib) = arrow.target_atoms
        if fa == fb:
            continue
        want[fa].append((ia, 1.0 if fb > fa else -1.0))
        want[fb].append((ib, 1.0 if fa > fb else -1.0))
    for f, items in enumerate(want):
        mol = stage.mols[f]
        if not items or mol.GetNumAtoms() < 2 or mol.GetNumConformers() == 0:
            continue
        conf = mol.GetConformer()
        pos = np.asarray(conf.GetPositions())[:, :2]
        centre = pos.mean(axis=0)
        rel = pos - centre
        extent = max(float(np.abs(rel).max()), 1e-6)
        best_theta, best_score, base_score = 0.0, -1e9, 0.0
        for theta in (0.0, 90.0, 180.0, 270.0):
            c, s_ = math.cos(math.radians(theta)), math.sin(math.radians(theta))
            score = sum((rel[i, 0] * c - rel[i, 1] * s_) / extent * side for i, side in items) / len(items)
            if theta == 0.0:
                base_score = score
            if score > best_score + 1e-9:
                best_theta, best_score = theta, score
        if best_theta == 0.0 or best_score - base_score < 0.25:
            continue
        c, s_ = math.cos(math.radians(best_theta)), math.sin(math.radians(best_theta))
        for i in range(mol.GetNumAtoms()):
            x, y = rel[i]
            conf.SetAtomPosition(i, Point3D(float(centre[0] + x * c - y * s_), float(centre[1] + x * s_ + y * c), 0.0))


def _tiny(mol) -> bool:
    """Molecules that are unreadable as skeletal drawings (CH3OH -> '-OH', CO2 -> 'O==O'): label every atom."""
    try:
        if mol.GetNumHeavyAtoms() <= 2:
            return True
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 6 and a.GetDegree() >= 2 and all(b.GetBondTypeAsDouble() >= 2 for b in a.GetBonds()):
                return True   # cumulated carbon (CO2, ketene, allene)
    except Exception:
        pass
    return False


def _label_atoms(stage: Stage, i: int) -> bool:
    la = stage.label_atoms[i] if i < len(stage.label_atoms) else None
    return _tiny(stage.mols[i]) if la is None else bool(la)


def _stage_size(doc: CDXMLDocument, stage: Stage, gap: float) -> Tuple[float, float]:
    from .layout import _measure_fragment
    L = doc.bond_length
    sizes = [_measure_fragment(doc, m, _label_atoms(stage, i)) for i, m in enumerate(stage.mols)]
    w = sum(s[0] for s in sizes) + (len(sizes) - 1) * (gap + 0.2 * L) + 0.6 * L
    h = max(s[1] for s in sizes) + (1.0 * L if stage.arrows else 0.0)
    return w, h


def _stage_extents(doc: CDXMLDocument, stage: Stage, text_size: float, has_next: bool = True) -> Tuple[float, float]:
    """(extent above the row baseline, extent below it) of a stage as drawn: structures (with room for electron
    arrows), per-structure captions, the stage caption and - when a reaction arrow follows (``has_next``) - the
    reagent lines stacked above that arrow and the condition lines below it (cf. ``layout._step_extents``)."""
    from .cdxml import text_box
    from .layout import CAPTION_GAP, _arrow_clearance, _caption_font, _measure_fragment, _text_stack
    L = doc.bond_length
    font = _caption_font(doc)
    pad = 0.5 * L if stage.arrows else 0.0
    up = down = 0.0
    for i, m in enumerate(stage.mols):
        half = 0.5 * _measure_fragment(doc, m, _label_atoms(stage, i))[1] + pad
        up = max(up, half)
        lab = stage.labels[i] if i < len(stage.labels) else None
        if lab:
            _, asc, desc = text_box(str(lab), text_size, font)
            half += CAPTION_GAP * L + asc + desc
        down = max(down, half)
    if stage.label:
        _, asc, desc = text_box(str(stage.label), text_size, font)
        down += 0.35 * L + asc + desc
    if has_next:
        clear = _arrow_clearance(doc, stage.arrow or "full")
        if stage.above:
            up = max(up, clear + _text_stack(doc, stage.above, text_size)[1])
        if stage.below:
            down = max(down, clear + _text_stack(doc, stage.below, text_size)[1])
    return up, down


def _transition_width(doc: CDXMLDocument, stage: Stage, gap: float, text_size: float) -> float:
    from .cdxml import text_box
    L = doc.bond_length
    widest = max([text_box(ln, text_size)[0] for ln in list(stage.above) + list(stage.below)] or [0.0])
    return gap + max(3.0 * L, widest + 0.8 * L) + gap


def build_mechanism(stages: Sequence[Any], style: Optional[str] = None, bond_length: Optional[float] = None,
                    name: str = "mechanism", title: Optional[str] = None, max_width: Optional[float] = None,
                    text_size: Optional[float] = None, row_gap: Optional[float] = None,
                    settings: Optional[Dict[str, Any]] = None, origin: Tuple[float, float] = (36.0, 72.0),
                    doc: Optional[CDXMLDocument] = None, arrow_color: Any = None, orient: bool = True) -> CDXMLDocument:
    """Lay out a mechanism: stages left to right joined by reaction arrows, wrapping into rows.

    Each stage is a dict ``{"structure": SMILES | mol | [SMILES, ...], "arrows": [specs], "label": caption,
    "arrow": kind of the arrow to the next stage ("full" default, "resonance", "equilibrium", ...),
    "above"/"below": text on that arrow, "label_atoms": force/suppress atom labels (default: automatic for
    tiny molecules such as CH3Cl, CH3OH, CO2)}`` (or a :class:`Stage`).  Charges and radicals come from the SMILES.
    With several structures in a stage, arrow specs address atoms as ``fragment:atom``; with ``orient=True`` a
    structure is turned by 180 degrees when the atom it bonds to in another structure faces away from it.
    ``origin`` = (left x, centre line y of the first row).
    """
    from . import layout
    doc = doc or CDXMLDocument(style, name=name, bond_length=bond_length, settings=settings)
    L = doc.bond_length
    gap = 1.0 * L
    text_size = text_size or doc.caption_size
    sts = [Stage.coerce(s) for s in stages]
    if not sts:
        raise ValueError("build_mechanism needs at least one stage")
    if orient:
        for st in sts:
            _orient_stage(st)
    page_w = doc.style.get("page", [0, 0, 540, 720])[2]
    max_width = max_width or max(300.0, page_w - 72.0)
    row_gap = row_gap if row_gap is not None else 2.5 * L
    sizes = [_stage_size(doc, s, gap) for s in sts]
    ups = [_stage_extents(doc, s, text_size, k < len(sts) - 1)[0] for k, s in enumerate(sts)]
    trans = [_transition_width(doc, s, gap, text_size) for s in sts]
    # greedy row assignment (the arrow leading to a wrapped stage stays at the end of the previous row)
    rows: List[List[int]] = [[0]]
    width = sizes[0][0]
    for k in range(1, len(sts)):
        need = trans[k - 1] + sizes[k][0]
        if width + need > max_width and rows[-1]:
            rows.append([k])
            width = sizes[k][0]
        else:
            rows[-1].append(k)
            width += need
    x0, y = origin
    up0 = max(ups[k] for k in rows[0])        # tallest structure / reagent stack of the first row
    if title:
        ty = max(18.0, y - up0 - 1.5 * L)
        t = doc.add_text(title, x0, ty, size=text_size + 1, face="bold", justification="Left", anchor="top")
        y = t.box.y1 + 1.5 * L + up0
        doc.mechanism_title = title            # type: ignore[attr-defined]
    else:
        y = max(y, 18.0 + up0)                 # keep the tallest reagent stack below the top margin
    stage_refs: List[List[ObjRef]] = [[] for _ in sts]
    stage_plus: List[List[ObjRef]] = [[] for _ in sts]
    meta: List[Dict[str, Any]] = []
    pending_step: Optional[Tuple[int, ObjRef, List[ObjRef], List[ObjRef]]] = None
    prev_arrow: Optional[ObjRef] = None
    row_bottom = y
    for ri, row in enumerate(rows):
        if ri > 0:
            y = row_bottom + row_gap + max(ups[k] for k in row)
        x = x0
        prev_arrow = None
        for k in row:
            st = sts[k]
            comps = [layout.Component(m, st.labels[i] if i < len(st.labels) else None, st.colors[i] if i < len(st.colors) else None,
                                      None, _label_atoms(st, i)) for i, m in enumerate(st.mols)]
            placed = layout._place_side(doc, comps, x, y, gap, text_size)
            refs, plusses, labels, x = placed[0], placed[1], placed[2], placed[-1]
            for fi, r in enumerate(refs):
                r.extra["stage"] = k
                r.extra["fragment_index"] = fi
            obstacles: List[Any] = list(plusses) + list(labels)
            if prev_arrow is not None:
                obstacles.append(prev_arrow)
            if k < len(sts) - 1:          # the reaction arrow that will follow this stage (drawn afterwards)
                obstacles.append(Box(x + gap, y - 0.05 * L, x + gap + trans[k] - 2 * gap, y + 0.05 * L))
            earrows = add_electron_arrows(doc, refs, st.arrows, st.mols, obstacles=obstacles, color=arrow_color)
            for a in earrows:
                a.extra["stage"] = k
            stage_refs[k], stage_plus[k] = refs, list(plusses)
            boxes = [r.box for r in refs] + [a.box for a in earrows] + [t.box for t in labels] + [p.box for p in plusses]
            for a in earrows:
                lid = a.extra.get("label_id")
                if lid is not None:
                    boxes += [t.box for t in doc.refs if t.id == lid]
            ub = union_boxes(boxes)
            cap = None
            if st.label and ub is not None:
                cap = doc.add_text(st.label, ub.cx, ub.y1 + 0.35 * L, size=text_size, anchor="top")
                boxes.append(cap.box)
            x = max(x, max(b.x1 for b in boxes))
            entry = {"index": k, "label": st.label, "fragment_ids": [r.id for r in refs],
                     "smiles": [r.extra.get("smiles") for r in refs], "electron_arrow_ids": [a.id for a in earrows],
                     "electron_arrows": [a.extra.get("spec") for a in earrows], "arrow_to_next": st.arrow if k < len(sts) - 1 else None,
                     "above": list(st.above), "below": list(st.below), "reaction_arrow_id": None, "row": ri,
                     "caption_id": cap.id if cap else None}
            meta.append(entry)
            if pending_step is not None:
                prev_k, arrow_ref, above_refs, below_refs = pending_step
                doc.add_step(stage_refs[prev_k], refs, [arrow_ref], above_refs, below_refs, stage_plus[prev_k] + list(plusses))
                pending_step = None
            ub = union_boxes(boxes)
            row_bottom = max(row_bottom, ub.y1 if ub else y)
            if k < len(sts) - 1:
                kind = st.arrow or "full"
                arrow_ref, above_refs, below_refs, x = layout._arrow_with_text(doc, x + gap, y, st.above, st.below, kind, 3.0 * L, text_size)
                arrow_ref.extra.update({"mechanism_step": True, "stage_from": k, "stage_to": k + 1, "arrow_kind": kind})
                entry["reaction_arrow_id"] = arrow_ref.id
                prev_arrow = arrow_ref
                if "reson" not in kind.lower():
                    pending_step = (k, arrow_ref, above_refs, below_refs)
                bottom = union_boxes([t.box for t in below_refs])
                row_bottom = max(row_bottom, bottom.y1 if bottom else y)
                x += gap
    doc.mechanism_stages = meta        # type: ignore[attr-defined]
    finish = getattr(layout, "_finish", None)
    return finish(doc) if callable(finish) else doc


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

def _arrow_summary(ref: ObjRef) -> dict:
    e = ref.extra
    return {"id": ref.id, "spec": e.get("spec"), "electrons": e.get("electrons"), "kind": e.get("arrow_kind"),
            "source": e.get("source"), "target": e.get("target"), "tail": e.get("tail"), "head": e.get("head"),
            "apex": e.get("apex"), "angle": e.get("angle"), "bulge": e.get("bulge"), "label": e.get("label")}


def describe_mechanism(doc: CDXMLDocument) -> dict:
    """Summary of a mechanism document: stages (structures, captions, electron arrows), reaction arrows, counts."""
    by_id = {r.id: r for r in doc.refs}
    e_arrows = [r for r in doc.refs if r.kind == "arrow" and r.extra.get("mechanism")]
    r_arrows = [r for r in doc.refs if r.kind == "arrow" and r.extra.get("mechanism_step")]
    content = doc.content_box()
    out: Dict[str, Any] = {
        "name": doc.name, "style": doc.style.get("name"), "bond_length": doc.bond_length,
        "title": getattr(doc, "mechanism_title", None), "stages": [],
        "electron_arrows": len(e_arrows), "reaction_arrows": [], "bounds": content.to_dict() if content else None,
    }
    stages_meta = getattr(doc, "mechanism_stages", None)
    if stages_meta:
        for sm in stages_meta:
            entry = {k: v for k, v in sm.items() if k != "electron_arrows"}
            entry["arrows"] = [_arrow_summary(by_id[i]) for i in sm.get("electron_arrow_ids", []) if i in by_id]
            out["stages"].append(entry)
    else:
        groups: Dict[Any, Dict[str, Any]] = {}
        for r in doc.refs:
            if r.kind != "fragment":
                continue
            key = r.extra.get("stage", r.id)
            g = groups.setdefault(key, {"index": key, "fragment_ids": [], "smiles": [], "arrows": [], "label": None})
            g["fragment_ids"].append(r.id)
            g["smiles"].append(r.extra.get("smiles"))
        for a in e_arrows:
            key = a.extra.get("stage")
            if key is None:
                fids = a.extra.get("fragment_ids") or []
                key = next((k for k, g in groups.items() if fids and fids[0] in g["fragment_ids"]), None)
            if key in groups:
                groups[key]["arrows"].append(_arrow_summary(a))
            else:
                groups.setdefault("unassigned", {"index": None, "fragment_ids": [], "smiles": [], "arrows": [], "label": None})["arrows"].append(_arrow_summary(a))
        out["stages"] = list(groups.values())
    for r in r_arrows:
        out["reaction_arrows"].append({"id": r.id, "kind": r.extra.get("arrow_kind"), "from_stage": r.extra.get("stage_from"),
                                       "to_stage": r.extra.get("stage_to"), "tail": r.extra.get("tail"), "head": r.extra.get("head")})
    warnings = getattr(doc, "warnings", None)
    if warnings:
        out["warnings"] = list(warnings)
    return out
