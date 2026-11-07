"""Text metrics for the CDXML writer and the layout engine.

ChemDraw re-measures every caption and atom label with the real font when it opens a file, so the only job
of this module is to predict those sizes well enough that nothing overlaps.  Two sources are used:

* **real fonts** - when Arial / Helvetica / Liberation Sans / DejaVu Sans (or the requested font itself, e.g.
  Times New Roman) can be found on the machine, glyph advances are measured once with Pillow
  (``ImageFont.truetype(...).getlength()`` at a 1000 px reference size) and cached per character;
* **heuristic** - otherwise an Arial/Helvetica advance-width table for printable ASCII is used, so results on a
  machine without fonts are within rounding of the measured ones for ordinary chemistry text.

Formula-face text (ChemDraw face 96: automatic sub/superscripts) is handled by :func:`formula_runs`, which
splits "NaBH4" into runs of normal / subscript / superscript characters; small runs are measured at
:data:`SMALL` times the size.

Atom labels are drawn by ChemDraw itself from ``Element`` + ``NumHydrogens`` (+ charge, isotope, radical).
:func:`fragment_labels` predicts their text, the side they extend to (away from the neighbours, horizontal)
and their bounding boxes so that fragments get label-aware bounds.

Set ``CHEMDRAWMCP_NO_FONTS=1`` to force the heuristic (reproducible output); ``CHEMDRAWMCP_FONT_DIR`` adds a
directory to the font search path.
"""
from __future__ import annotations

import math
import os
import sys
from functools import lru_cache
from typing import Any, Dict, List, Optional, Sequence, Tuple

__all__ = ["text_width", "text_box", "formula_runs", "glyph_center", "font_path", "fragment_labels",
           "atom_label", "reset_cache", "ARIAL_WIDTHS", "ASCENT", "DESCENT", "LINE_HEIGHT", "SMALL"]

# vertical metrics as fractions of the font size (kept identical to the historical heuristic so that boxes of
# existing documents do not move; Arial's real hhea values are 0.906 / 0.212)
ASCENT = 0.9
DESCENT = 0.25
LINE_HEIGHT = 1.2
# ChemDraw formula face: sub/superscripts are drawn at roughly this fraction of the size
SMALL = 0.6
SUB_SHIFT = 0.25      # subscript baseline shift down (fraction of size)
SUP_SHIFT = 0.35      # superscript baseline shift up

_REF = 1000           # reference pixel size for measuring

# Arial / Helvetica advance widths (em units) for printable ASCII - measured from Arial.ttf (identical to the
# Helvetica AFM metrics).  This is the no-font fallback.
ARIAL_WIDTHS: Dict[str, float] = {
    " ": 0.278, "!": 0.278, '"': 0.355, "#": 0.556, "$": 0.556, "%": 0.889, "&": 0.667, "'": 0.191, "(": 0.333,
    ")": 0.333, "*": 0.389, "+": 0.584, ",": 0.278, "-": 0.333, ".": 0.278, "/": 0.278, "0": 0.556, "1": 0.556,
    "2": 0.556, "3": 0.556, "4": 0.556, "5": 0.556, "6": 0.556, "7": 0.556, "8": 0.556, "9": 0.556, ":": 0.278,
    ";": 0.278, "<": 0.584, "=": 0.584, ">": 0.584, "?": 0.556, "@": 1.015, "A": 0.667, "B": 0.667, "C": 0.722,
    "D": 0.722, "E": 0.667, "F": 0.611, "G": 0.778, "H": 0.722, "I": 0.278, "J": 0.5, "K": 0.667, "L": 0.556,
    "M": 0.833, "N": 0.722, "O": 0.778, "P": 0.667, "Q": 0.778, "R": 0.722, "S": 0.667, "T": 0.611, "U": 0.722,
    "V": 0.667, "W": 0.944, "X": 0.667, "Y": 0.667, "Z": 0.611, "[": 0.278, "\\": 0.278, "]": 0.278, "^": 0.469,
    "_": 0.556, "`": 0.333, "a": 0.556, "b": 0.556, "c": 0.5, "d": 0.556, "e": 0.556, "f": 0.278, "g": 0.556,
    "h": 0.556, "i": 0.222, "j": 0.222, "k": 0.5, "l": 0.222, "m": 0.833, "n": 0.556, "o": 0.556, "p": 0.556,
    "q": 0.556, "r": 0.333, "s": 0.5, "t": 0.278, "u": 0.556, "v": 0.5, "w": 0.722, "x": 0.5, "y": 0.5, "z": 0.5,
    "{": 0.334, "|": 0.26, "}": 0.334, "~": 0.584,
    # a few non-ASCII characters common in reaction conditions
    "\u00b0": 0.4, "\u00b7": 0.333, "\u2022": 0.35, "\u2212": 0.584, "\u2013": 0.556, "\u2014": 1.0, "\u2192": 1.0,
    "\u00b5": 0.556, "\u0394": 0.667, "\u00e9": 0.556, "\u00c5": 0.667, "\u2019": 0.222, "\u00d7": 0.584,
}
# vertical glyph bounds (top, bottom) above the baseline in em for the glyphs whose centring matters
_ARIAL_VBOUNDS: Dict[str, Tuple[float, float]] = {"+": (0.589, 0.0), "-": (0.303, 0.0), "=": (0.503, 0.0),
                                                   "\u2212": (0.503, 0.0), "H": (0.716, 0.0), "x": (0.519, 0.0)}

# font files to try, per requested family; the generic sans fallbacks are appended to every list
_SANS_FALLBACK = ["arial.ttf", "arialmt.ttf", "helvetica.ttc", "helvetica.ttf", "liberationsans-regular.ttf",
                  "dejavusans.ttf", "nimbussans-regular.otf", "freesans.ttf"]
_FAMILY_FILES: Dict[str, List[str]] = {
    "arial": ["arial.ttf", "arialmt.ttf"],
    "helvetica": ["helvetica.ttc", "helvetica.ttf", "arial.ttf"],
    "liberation sans": ["liberationsans-regular.ttf"],
    "dejavu sans": ["dejavusans.ttf"],
    "times new roman": ["times new roman.ttf", "times.ttf", "timesnewroman.ttf", "times.ttc",
                        "liberationserif-regular.ttf", "dejavuserif.ttf"],
    "times": ["times.ttc", "times.ttf", "times new roman.ttf", "liberationserif-regular.ttf", "dejavuserif.ttf"],
}

_FACES: Dict[str, Any] = {}          # family key -> _Face | _HeuristicFace
_INDEX: Optional[Dict[str, str]] = None


# ---------------------------------------------------------------------------
# font discovery
# ---------------------------------------------------------------------------

def _font_dirs() -> List[str]:
    dirs: List[str] = []
    extra = os.environ.get("CHEMDRAWMCP_FONT_DIR")
    if extra:
        dirs.append(extra)
    home = os.path.expanduser("~")
    if sys.platform == "darwin":
        dirs += ["/System/Library/Fonts/Supplemental", "/Library/Fonts", "/System/Library/Fonts",
                 os.path.join(home, "Library", "Fonts")]
    elif sys.platform.startswith("win"):
        windir = os.environ.get("WINDIR") or os.environ.get("SystemRoot") or "C:\\Windows"
        dirs += [os.path.join(windir, "Fonts"),
                 os.path.join(os.environ.get("LOCALAPPDATA", ""), "Microsoft", "Windows", "Fonts")]
    else:
        dirs += ["/usr/share/fonts", "/usr/local/share/fonts", os.path.join(home, ".fonts"),
                 os.path.join(home, ".local", "share", "fonts"), "/usr/share/fonts/truetype"]
    return [d for d in dirs if d and os.path.isdir(d)]


def _font_index() -> Dict[str, str]:
    """lower-case file name -> path for every .ttf/.ttc/.otf under the font directories (depth <= 3)."""
    global _INDEX
    if _INDEX is not None:
        return _INDEX
    index: Dict[str, str] = {}
    if os.environ.get("CHEMDRAWMCP_NO_FONTS"):
        _INDEX = index
        return index
    for base in _font_dirs():
        base_depth = base.rstrip(os.sep).count(os.sep)
        for root, dirnames, files in os.walk(base):
            if root.count(os.sep) - base_depth >= 3:
                dirnames[:] = []
            for fn in files:
                low = fn.lower()
                if low.endswith((".ttf", ".ttc", ".otf")) and low not in index:
                    index[low] = os.path.join(root, fn)
    _INDEX = index
    return index


def _family_key(font: Optional[str]) -> str:
    return (font or "Arial").strip().lower()


def font_path(font: Optional[str] = None) -> Optional[str]:
    """Path of the font file that will be used to measure ``font`` (None when the heuristic is used)."""
    face = _face(font)
    return getattr(face, "path", None)


# ---------------------------------------------------------------------------
# faces
# ---------------------------------------------------------------------------

def _heuristic_advance(ch: str) -> float:
    w = ARIAL_WIDTHS.get(ch)
    if w is not None:
        return w
    if ch.islower():
        return 0.556
    if ch.isupper():
        return 0.667
    if ch.isdigit():
        return 0.556
    if ch.isspace():
        return 0.278
    return 0.584


class _HeuristicFace:
    path = None
    ascent = ASCENT
    descent = DESCENT

    def advance(self, ch: str) -> float:
        return _heuristic_advance(ch)

    def vbounds(self, ch: str) -> Tuple[float, float]:
        return _ARIAL_VBOUNDS.get(ch, (0.716 if ch.isupper() or ch.isdigit() else 0.52, 0.0))


class _Face:
    """A Pillow font measured at the reference size with a per-character advance cache."""

    def __init__(self, path: str, index: int = 0):
        from PIL import ImageFont
        self.path = path
        self.font = ImageFont.truetype(path, _REF, index=index)
        asc, desc = self.font.getmetrics()
        self.ascent = asc / float(_REF)
        self.descent = desc / float(_REF)
        self._adv: Dict[str, float] = {}
        self._vb: Dict[str, Tuple[float, float]] = {}

    def advance(self, ch: str) -> float:
        w = self._adv.get(ch)
        if w is None:
            try:
                w = float(self.font.getlength(ch)) / _REF
            except Exception:
                w = _heuristic_advance(ch)
            if w <= 0 and not ch.isspace():
                w = _heuristic_advance(ch)
            self._adv[ch] = w
        return w

    def vbounds(self, ch: str) -> Tuple[float, float]:
        vb = self._vb.get(ch)
        if vb is None:
            try:
                bb = self.font.getbbox(ch)
                asc_px = self.font.getmetrics()[0]
                vb = ((asc_px - bb[1]) / float(_REF), (asc_px - bb[3]) / float(_REF))
            except Exception:
                vb = _HeuristicFace().vbounds(ch)
            self._vb[ch] = vb
        return vb


def _face(font: Optional[str]):
    key = _family_key(font)
    face = _FACES.get(key)
    if face is not None:
        return face
    face = _HeuristicFace()
    if not os.environ.get("CHEMDRAWMCP_NO_FONTS"):
        index = _font_index()
        candidates = list(_FAMILY_FILES.get(key, [key + ".ttf", key.replace(" ", "") + ".ttf"])) + _SANS_FALLBACK
        for fn in candidates:
            path = index.get(fn)
            if path:
                try:
                    face = _Face(path, 0)
                    break
                except Exception:
                    continue
    _FACES[key] = face
    return face


def reset_cache() -> None:
    """Forget discovered fonts (tests / after changing the environment)."""
    global _INDEX
    _INDEX = None
    _FACES.clear()
    _em_width.cache_clear()
    formula_runs.cache_clear()


# ---------------------------------------------------------------------------
# formula face
# ---------------------------------------------------------------------------

_CHARGE_CHARS = "+-\u2212"


@lru_cache(maxsize=4096)
def formula_runs(text: str) -> Tuple[Tuple[str, str], ...]:
    """Split formula-face text into (chunk, kind) runs, kind in {"normal", "sub", "sup"}.

    ChemDraw's formula face subscripts digits that follow a letter or a closing bracket ("H2SO4", "Pd(PPh3)4",
    "C12H26") and superscripts a trailing charge sign after an element symbol or subscript ("H3O+", "Cl-").
    Digits at the start of a token stay normal ("1. Mg", "2 equiv", "80 °C").  Cached; returns a tuple.
    """
    runs: List[Tuple[str, str]] = []
    n = len(text)
    i = 0
    prev = "other"       # kind of the previous character: letter | close | sub | sup | digit | other
    while i < n:
        c = text[i]
        if c.isdigit():
            j = i
            while j < n and text[j].isdigit():
                j += 1
            if prev in ("letter", "close", "sub"):
                runs.append((text[i:j], "sub"))
                prev = "sub"
            else:
                runs.append((text[i:j], "normal"))
                prev = "digit"
            i = j
            continue
        if c in _CHARGE_CHARS and prev in ("letter", "sub") and (i + 1 >= n or not text[i + 1].isalnum()):
            runs.append((c, "sup"))
            prev = "sup"
            i += 1
            continue
        if c.isalpha():
            prev = "letter"
        elif c in ")]":
            prev = "close"
        else:
            prev = "other"
        runs.append((c, "normal"))
        i += 1
    # merge neighbouring runs of the same kind
    merged: List[Tuple[str, str]] = []
    for chunk, kind in runs:
        if merged and merged[-1][1] == kind:
            merged[-1] = (merged[-1][0] + chunk, kind)
        else:
            merged.append((chunk, kind))
    return tuple(merged)


# ---------------------------------------------------------------------------
# public measurements
# ---------------------------------------------------------------------------

def _runs_width(runs: Sequence[Tuple[str, str]], size: float, face) -> float:
    w = 0.0
    for chunk, kind in runs:
        s = size * (SMALL if kind != "normal" else 1.0)
        w += sum(face.advance(ch) for ch in chunk) * s
    return w


def text_width(text: str, size: float, font: Optional[str] = None, formula: bool = False) -> float:
    """Width in points of one line of ``text`` at ``size`` points (fractional widths)."""
    return _em_width(text, _family_key(font), bool(formula)) * size


@lru_cache(maxsize=8192)
def _em_width(text: str, family: str, formula: bool) -> float:
    face = _face(family)
    if formula:
        return _runs_width(formula_runs(text), 1.0, face)
    return sum(face.advance(ch) for ch in text)


def text_box(text: str, size: float, font: Optional[str] = None, formula: bool = False) -> Tuple[float, float, float]:
    """(width, ascent, descent) in points of a (possibly multi-line) text.

    ascent is the height above the first baseline, descent everything below it including further lines.
    Formula text with subscripts gets a slightly deeper descent (the subscript digits hang below the baseline).
    """
    lines = text.split("\n") or [""]
    w = 0.0
    sub = False
    for ln in lines:
        w = max(w, text_width(ln, size, font, formula))
        if formula and any(kind == "sub" for _, kind in formula_runs(ln)):
            sub = True
    line_h = size * LINE_HEIGHT
    ascent = size * ASCENT
    descent = size * DESCENT + line_h * (len(lines) - 1)
    if sub:
        descent += size * 0.1
    return w, ascent, descent


def glyph_center(ch: str, size: float, font: Optional[str] = None) -> float:
    """Height of the visual centre of glyph ``ch`` above the baseline, in points (e.g. 0.29 * size for "+")."""
    top, bottom = _face(font).vbounds(ch[:1] or "+")
    return (top + bottom) / 2.0 * size


# ---------------------------------------------------------------------------
# atom labels
# ---------------------------------------------------------------------------

def _charge_text(charge: int) -> str:
    if not charge:
        return ""
    sign = "+" if charge > 0 else "-"
    return sign if abs(charge) == 1 else "%d%s" % (abs(charge), sign)


def atom_label(symbol: str, n_h: int = 0, charge: int = 0, isotope: int = 0, side: str = "right",
               show_h: bool = True) -> Tuple[str, List[Tuple[str, str]]]:
    """Predict the label ChemDraw draws for an atom: (text, runs).

    side: "left" -> the label extends to the left of the atom, hydrogens first ("H2N", "HO");
    "right" / "center" -> symbol first ("NH2", "OH").  Isotopes are superscript prefixes ("13C"), charges
    superscript suffixes ("NH3+", "O-", "Fe2+").
    """
    runs: List[Tuple[str, str]] = []
    sym: List[Tuple[str, str]] = ([(str(isotope), "sup")] if isotope else []) + [(symbol, "normal")]
    hs: List[Tuple[str, str]] = []
    if show_h and n_h > 0:
        hs.append(("H", "normal"))
        if n_h > 1:
            hs.append((str(n_h), "sub"))
    if side == "left" and hs:
        runs += hs + sym
    else:
        runs += sym + hs
    ch = _charge_text(charge)
    if ch:
        runs.append((ch, "sup"))
    return "".join(chunk for chunk, _ in runs), runs


def _label_side(deg: int, mean_dx: float) -> str:
    if deg == 0:
        return "center"
    # neighbours mostly to the right -> label extends left (hydrogens first); otherwise it extends right.
    return "left" if mean_dx > 0.15 else "right"


def fragment_labels(mol, pts: Sequence[Tuple[float, float]], size: float, label_atoms: bool = False,
                    font: Optional[str] = None, margin: float = 1.6, show_h: bool = True,
                    show_terminal_carbons: bool = False) -> List[Dict[str, Any]]:
    """Predict ChemDraw's atom labels for ``mol`` whose atoms sit at ``pts`` (page coordinates, y down).

    Returns one dict per labelled atom: {"atom": idx, "text": "H3N+", "side": "left"|"right"|"center",
    "box": [x0, y0, x1, y1]} (box already padded by ``margin``, ChemDraw's label margin).  Atoms with a label:
    non-carbon, charged, isotope, radical, degree 0, or all atoms when ``label_atoms``.  Charged / radical carbons
    that ChemDraw draws without the "C" only get their charge sign.
    """
    out: List[Dict[str, Any]] = []
    n = mol.GetNumAtoms()
    if n == 0 or n != len(pts):
        return out
    face = _face(font)
    for i in range(n):            # GetAtomWithIdx avoids RDKit's slow Python sequence wrapper
        atom = mol.GetAtomWithIdx(i)
        z = atom.GetAtomicNum()
        deg = atom.GetDegree()
        if z == 6 and deg > 0 and not label_atoms and not (show_terminal_carbons and deg == 1):
            # plain carbon: only a charge / isotope / radical gives it a label
            if not (atom.GetFormalCharge() or atom.GetIsotope() or atom.GetNumRadicalElectrons()):
                continue
        charge = atom.GetFormalCharge()
        iso = atom.GetIsotope()
        rad = atom.GetNumRadicalElectrons()
        explicit = bool(label_atoms or z == 0) and z != 1
        full = bool(z != 6 or deg == 0 or iso or explicit or (show_terminal_carbons and deg == 1))
        if not (full or charge or rad):
            continue
        x, y = pts[i]
        # mean direction towards the neighbours (unit vectors)
        sx = sy = 0.0
        for nb in atom.GetNeighbors():
            nx, ny = pts[nb.GetIdx()]
            d = math.hypot(nx - x, ny - y)
            if d > 1e-9:
                sx += (nx - x) / d
                sy += (ny - y) / d
        mean_dx = sx / deg if deg else 0.0
        side = "center" if explicit else _label_side(deg, mean_dx)
        if full:
            if z == 0:
                sym = atom.GetProp("dummyLabel") if atom.HasProp("dummyLabel") else (
                    atom.GetProp("_MolFileRLabel") if atom.HasProp("_MolFileRLabel") else "R")
                n_h = 0
            else:
                sym = atom.GetSymbol()
                n_h = atom.GetTotalNumHs()
            text, runs = atom_label(sym, n_h, charge, iso, side, show_h=show_h or deg == 0)
        else:
            # charged / radical carbon vertex: ChemDraw draws only the charge sign (or a radical dot) next to the vertex
            text = _charge_text(charge) or ("•" if rad else "")
            runs = [(text, "sup")] if text else []
            sym = ""
        if not text:
            continue
        total = _runs_width(runs, size, face)
        sym_w = sum(face.advance(c) for c in sym) * size if sym else 0.0
        charge_w = sum(face.advance(c) for c in _charge_text(charge)) * size * SMALL if charge else 0.0
        iso_w = sum(face.advance(c) for c in str(iso)) * size * SMALL if iso else 0.0
        if side == "center" or not sym:
            left = right = total / 2.0
        elif side == "left":
            # hydrogens first, symbol centred on the atom; be conservative and keep the charge on both sides
            left = total - sym_w / 2.0
            right = sym_w / 2.0 + charge_w
        else:
            left = sym_w / 2.0 + iso_w
            right = total - iso_w - sym_w / 2.0
        top = 0.4 * size + (0.05 * size if iso else 0.0)
        bottom = 0.4 * size + (0.2 * size if any(k == "sub" for _, k in runs) else 0.0)
        if not sym:   # bare charge sign: a small superscript next to the vertex
            top, bottom = 0.5 * size, 0.1 * size
        out.append({"atom": i, "text": text, "side": side,
                    "box": [x - left - margin, y - top - margin, x + right + margin, y + bottom + margin]})
    return out
