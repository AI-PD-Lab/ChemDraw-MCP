"""Nearest-named-colour lookup used by the render analysis.

The table mixes CSS basics with the colour names molecular viewers use by
default (PyMOL element/chain palettes, ChimeraX model colours) so that the
agent reads "marine", "salmon" or "cornflowerblue" rather than raw RGB.
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

NAMED_COLORS: Dict[str, Tuple[int, int, int]] = {
    # neutrals
    "black": (0, 0, 0),
    "white": (255, 255, 255),
    "grey": (128, 128, 128),
    "lightgrey": (200, 200, 200),
    "darkgrey": (64, 64, 64),
    "silver": (192, 192, 192),
    # PyMOL classics
    "red": (255, 0, 0),
    "green": (0, 255, 0),
    "blue": (0, 0, 255),
    "yellow": (255, 255, 0),
    "magenta": (255, 0, 255),
    "cyan": (0, 255, 255),
    "orange": (255, 140, 0),
    "salmon": (255, 153, 153),
    "lime": (128, 255, 0),
    "slate": (128, 128, 255),
    "marine": (0, 128, 255),
    "forest": (51, 153, 51),
    "deepteal": (26, 153, 153),
    "teal": (0, 191, 191),
    "olive": (196, 196, 51),
    "purple": (191, 0, 191),
    "hotpink": (255, 0, 128),
    "wheat": (255, 230, 204),
    "violet": (255, 128, 255),
    "brown": (165, 102, 51),
    "tv_red": (255, 51, 51),
    "tv_green": (51, 255, 51),
    "tv_blue": (77, 77, 255),
    "tv_yellow": (255, 255, 51),
    "tv_orange": (255, 140, 26),
    "lightblue": (191, 204, 255),
    "lightpink": (255, 204, 204),
    "palegreen": (166, 230, 166),
    "paleyellow": (255, 255, 128),
    "palecyan": (204, 255, 255),
    "lightorange": (255, 204, 128),
    "lightmagenta": (255, 128, 255),
    "smudge": (140, 179, 102),
    "skyblue": (84, 171, 236),
    "density": (26, 26, 153),
    "firebrick": (178, 34, 34),
    "chocolate": (210, 105, 30),
    "gold": (255, 215, 0),
    "limon": (191, 255, 128),
    "splitpea": (133, 179, 69),
    "sand": (240, 190, 100),
    "ruby": (153, 26, 51),
    "raspberry": (179, 51, 102),
    "dirtyviolet": (179, 128, 179),
    "greencyan": (64, 255, 191),
    "bluewhite": (217, 217, 255),
    "yelloworange": (255, 222, 0),
    "deepblue": (0, 0, 153),
    "deeppurple": (153, 0, 153),
    "darkred": (139, 0, 0),
    "darkgreen": (0, 100, 0),
    "navy": (0, 0, 128),
    # ChimeraX default model colours and other CSS names
    "tan": (210, 180, 140),
    "skyblue_css": (135, 206, 235),
    "plum": (221, 160, 221),
    "lightgreen": (144, 238, 144),
    "cornflowerblue": (100, 149, 237),
    "goldenrod": (218, 165, 32),
    "coral": (255, 127, 80),
    "orchid": (218, 112, 214),
    "khaki": (240, 230, 140),
    "turquoise": (64, 224, 208),
    "indigo": (75, 0, 130),
    "crimson": (220, 20, 60),
    "steelblue": (70, 130, 180),
    "seagreen": (46, 139, 87),
    "sienna": (160, 82, 45),
    "pink": (255, 192, 203),
    "beige": (245, 245, 220),
    "ivory": (255, 255, 240),
    "lavender": (230, 230, 250),
    "mediumpurple": (147, 112, 219),
    "dodgerblue": (30, 144, 255),
    "tomato": (255, 99, 71),
    "limegreen": (50, 205, 50),
    "chartreuse": (127, 255, 0),
    "springgreen": (0, 255, 127),
    "aquamarine": (127, 255, 212),
    "royalblue": (65, 105, 225),
    "darkorange": (255, 140, 0),
    "peru": (205, 133, 63),
}

# Coarse hue families used for a quick "what colours are in the frame" read.
HUE_FAMILIES = [
    ("red", 345, 15),
    ("orange", 15, 45),
    ("yellow", 45, 70),
    ("green", 70, 170),
    ("cyan", 170, 200),
    ("blue", 200, 260),
    ("purple", 260, 300),
    ("magenta", 300, 345),
]


def _srgb_to_linear(c: float) -> float:
    c = c / 255.0
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def rgb_to_lab(rgb) -> Tuple[float, float, float]:
    r, g, b = (_srgb_to_linear(float(x)) for x in rgb)
    # D65 reference white
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
    y = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 1.0
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883

    def f(t: float) -> float:
        return t ** (1.0 / 3.0) if t > 0.008856 else 7.787 * t + 16.0 / 116.0

    fx, fy, fz = f(x), f(y), f(z)
    return (116.0 * fy - 16.0, 500.0 * (fx - fy), 200.0 * (fy - fz))


_LAB_TABLE: List[Tuple[str, Tuple[float, float, float]]] = [
    (name, rgb_to_lab(rgb)) for name, rgb in NAMED_COLORS.items()
]


def nearest_color_name(rgb) -> Tuple[str, float]:
    """Return (name, distance) of the closest named colour in CIE Lab."""
    lab = rgb_to_lab(rgb)
    best, best_d = "unknown", float("inf")
    for name, ref in _LAB_TABLE:
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(lab, ref)))
        if d < best_d:
            best, best_d = name, d
    if best.endswith("_css"):
        best = best[:-4]
    return best, best_d


def rgb_to_hsv(rgb) -> Tuple[float, float, float]:
    r, g, b = (float(x) / 255.0 for x in rgb)
    mx, mn = max(r, g, b), min(r, g, b)
    d = mx - mn
    if d == 0:
        h = 0.0
    elif mx == r:
        h = (60.0 * ((g - b) / d) + 360.0) % 360.0
    elif mx == g:
        h = 60.0 * ((b - r) / d) + 120.0
    else:
        h = 60.0 * ((r - g) / d) + 240.0
    s = 0.0 if mx == 0 else d / mx
    return h, s, mx


def hue_family(rgb) -> str:
    h, s, v = rgb_to_hsv(rgb)
    if v < 0.12:
        return "black"
    if s < 0.15:
        return "white" if v > 0.85 else "grey"
    for name, lo, hi in HUE_FAMILIES:
        if lo > hi:  # wraps around 360
            if h >= lo or h < hi:
                return name
        elif lo <= h < hi:
            return name
    return "unknown"


def hex_color(rgb) -> str:
    return "#%02x%02x%02x" % tuple(int(x) for x in rgb[:3])
