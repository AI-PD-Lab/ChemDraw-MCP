"""Numerical analysis of rendered images (no vision model required).

Everything here works on RGBA uint8 numpy arrays.  The public entry points
return plain JSON-serialisable dictionaries that an agent can read directly:

* :func:`analyze_render`        - composition, coverage, clipping, palette,
                                  contrast and an ASCII density map.
* :func:`analyze_segmentation`  - per-object / per-chain visibility from an
                                  "ID pass" render where every item was drawn
                                  in a unique flat colour.
* :func:`diff_images`           - what changed between two renders.
"""
from __future__ import annotations

import base64
import io
import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
from PIL import Image

from .colornames import hex_color, hue_family, nearest_color_name

DENSITY_CHARS = " .:-=+*#%@"
LABEL_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

ImageSource = Union[str, bytes, bytearray, Image.Image, np.ndarray]


# --------------------------------------------------------------------------
# loading / helpers
# --------------------------------------------------------------------------

def load_rgba(src: ImageSource) -> np.ndarray:
    """Return an (H, W, 4) uint8 RGBA array from a path, PNG bytes, PIL image
    or numpy array."""
    if isinstance(src, np.ndarray):
        arr = src
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        if arr.shape[2] == 3:
            alpha = np.full(arr.shape[:2] + (1,), 255, dtype=np.uint8)
            arr = np.concatenate([arr.astype(np.uint8), alpha], axis=2)
        return np.ascontiguousarray(arr[:, :, :4].astype(np.uint8))
    if isinstance(src, Image.Image):
        img = src
    elif isinstance(src, (bytes, bytearray)):
        img = Image.open(io.BytesIO(bytes(src)))
    else:
        img = Image.open(str(src))
    return np.array(img.convert("RGBA"), dtype=np.uint8)


def png_bytes_from_b64(b64: str) -> bytes:
    return base64.b64decode(b64)


def thumbnail_png(src: ImageSource, max_dim: int = 512) -> bytes:
    """Downscale an image to at most ``max_dim`` on its long side; return PNG bytes."""
    rgba = load_rgba(src)
    img = Image.fromarray(rgba)
    img.thumbnail((max_dim, max_dim))
    buf = io.BytesIO()
    img.save(buf, format="PNG", optimize=True)
    return buf.getvalue()


def _r(x: float, nd: int = 3) -> float:
    return float(round(float(x), nd))


def _integral(mask: np.ndarray) -> np.ndarray:
    ii = np.zeros((mask.shape[0] + 1, mask.shape[1] + 1), dtype=np.int64)
    ii[1:, 1:] = np.cumsum(np.cumsum(mask.astype(np.int64), axis=0), axis=1)
    return ii


def _cell_bounds(n: int, cells: int) -> List[Tuple[int, int]]:
    edges = np.linspace(0, n, cells + 1).round().astype(int)
    return [(int(edges[i]), int(max(edges[i + 1], edges[i] + 1))) for i in range(cells)]


def grid_occupancy(mask: np.ndarray, cols: int, rows: int) -> np.ndarray:
    """Fraction of True pixels per cell of a rows x cols grid."""
    h, w = mask.shape
    cols = max(1, min(cols, w))
    rows = max(1, min(rows, h))
    ii = _integral(mask)
    out = np.zeros((rows, cols), dtype=np.float64)
    ys = _cell_bounds(h, rows)
    xs = _cell_bounds(w, cols)
    for r, (y0, y1) in enumerate(ys):
        for c, (x0, x1) in enumerate(xs):
            s = ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]
            out[r, c] = s / float((y1 - y0) * (x1 - x0))
    return out


def ascii_map(cells: np.ndarray, chars: str = DENSITY_CHARS) -> List[str]:
    n = len(chars) - 1
    lines = []
    for row in cells:
        lines.append("".join(chars[int(round(min(max(v, 0.0), 1.0) * n))] for v in row))
    return lines


def bbox_of(mask: np.ndarray) -> Optional[Tuple[int, int, int, int]]:
    """(x0, y0, x1, y1) with exclusive upper bounds, or None if empty."""
    rows = np.flatnonzero(mask.any(axis=1))
    cols = np.flatnonzero(mask.any(axis=0))
    if rows.size == 0 or cols.size == 0:
        return None
    return int(cols[0]), int(rows[0]), int(cols[-1]) + 1, int(rows[-1]) + 1


def _bbox_dict(bb: Optional[Tuple[int, int, int, int]], w: int, h: int) -> Optional[dict]:
    if bb is None:
        return None
    x0, y0, x1, y1 = bb
    return {
        "px": [x0, y0, x1, y1],
        "frac": [_r(x0 / w), _r(y0 / h), _r(x1 / w), _r(y1 / h)],
        "width_frac": _r((x1 - x0) / w),
        "height_frac": _r((y1 - y0) / h),
    }


def _centroid(mask: np.ndarray) -> Optional[Tuple[float, float]]:
    ys, xs = np.nonzero(mask)
    if xs.size == 0:
        return None
    return float(xs.mean()), float(ys.mean())


def _edge_coverage(mask: np.ndarray, band: int = 2) -> Dict[str, float]:
    h, w = mask.shape
    band = max(1, min(band, h, w))
    return {
        "top": _r(mask[:band, :].mean()),
        "bottom": _r(mask[-band:, :].mean()),
        "left": _r(mask[:, :band].mean()),
        "right": _r(mask[:, -band:].mean()),
    }


def _luma(rgb: np.ndarray) -> np.ndarray:
    rgb = rgb.astype(np.float64)
    return 0.2126 * rgb[..., 0] + 0.7152 * rgb[..., 1] + 0.0722 * rgb[..., 2]


def _rel_luminance(rgb: Sequence[float]) -> float:
    def lin(c: float) -> float:
        c = c / 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    r, g, b = (lin(float(x)) for x in rgb[:3])
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast_ratio(rgb_a: Sequence[float], rgb_b: Sequence[float]) -> float:
    la, lb = _rel_luminance(rgb_a), _rel_luminance(rgb_b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def estimate_background(rgba: np.ndarray) -> Tuple[Tuple[int, int, int], bool, float]:
    """Return (bg_rgb, transparent, border_agreement).

    If the alpha channel is used, the background is transparent and bg_rgb is
    the colour of fully transparent pixels (usually black).  Otherwise the
    most common colour along the 1-pixel border is taken as background and
    border_agreement says what fraction of the border matched it.
    """
    alpha = rgba[:, :, 3]
    if (alpha < 255).any():
        bg_px = rgba[alpha == 0][:, :3]
        if bg_px.size:
            vals, counts = np.unique(bg_px, axis=0, return_counts=True)
            rgb = tuple(int(v) for v in vals[counts.argmax()])
        else:
            rgb = (0, 0, 0)
        return rgb, True, 1.0
    border = np.concatenate([rgba[0, :, :3], rgba[-1, :, :3], rgba[:, 0, :3], rgba[:, -1, :3]])
    vals, counts = np.unique(border, axis=0, return_counts=True)
    i = counts.argmax()
    rgb = tuple(int(v) for v in vals[i])
    return rgb, False, float(counts[i] / border.shape[0])


def foreground_mask(rgba: np.ndarray, bg_rgb: Sequence[int], transparent: bool, tol: int = 18) -> np.ndarray:
    if transparent:
        return rgba[:, :, 3] > 8
    diff = np.abs(rgba[:, :, :3].astype(np.int16) - np.array(bg_rgb, dtype=np.int16)).max(axis=2)
    return diff > tol


def _sample(arr: np.ndarray, max_n: int = 400_000) -> np.ndarray:
    if arr.shape[0] <= max_n:
        return arr
    step = int(math.ceil(arr.shape[0] / max_n))
    return arr[::step]


def palette(rgb_pixels: np.ndarray, k: int = 8, quant: int = 24) -> List[dict]:
    """Top-k quantised colours of the given (N, 3) pixel array."""
    if rgb_pixels.shape[0] == 0:
        return []
    px = _sample(rgb_pixels)
    q = (px.astype(np.int32) // quant) * quant + quant // 2
    q = np.clip(q, 0, 255)
    keys = (q[:, 0] << 16) | (q[:, 1] << 8) | q[:, 2]
    vals, counts = np.unique(keys, return_counts=True)
    order = counts.argsort()[::-1][:k]
    total = float(px.shape[0])
    # merge quantisation bins that resolve to the same colour name (shading of one colour)
    merged: Dict[str, dict] = {}
    for i in counts.argsort()[::-1][:k * 4]:
        key = int(vals[i])
        sel = px[keys == key]
        mean = sel.mean(axis=0)
        name, dist = nearest_color_name(tuple(int(round(v)) for v in mean))
        n = int(counts[i])
        m = merged.get(name)
        if m is None:
            merged[name] = {"name": name, "sum": mean * n, "n": n, "dist": dist}
        else:
            m["sum"] = m["sum"] + mean * n
            m["n"] += n
            m["dist"] = min(m["dist"], dist)
    out = []
    for m in sorted(merged.values(), key=lambda d: -d["n"])[:k]:
        mean = tuple(int(round(v)) for v in (m["sum"] / m["n"]))
        out.append({
            "hex": hex_color(mean),
            "rgb": list(mean),
            "name": m["name"],
            "name_distance": _r(m["dist"], 1),
            "family": hue_family(mean),
            "share": _r(m["n"] / total),
        })
    return out


def hue_shares(rgb_pixels: np.ndarray) -> Dict[str, float]:
    if rgb_pixels.shape[0] == 0:
        return {}
    px = _sample(rgb_pixels, 200_000).astype(np.float64) / 255.0
    r, g, b = px[:, 0], px[:, 1], px[:, 2]
    mx = px.max(axis=1)
    mn = px.min(axis=1)
    d = mx - mn
    s = np.where(mx > 0, d / np.maximum(mx, 1e-9), 0.0)
    v = mx
    h = np.zeros_like(mx)
    nz = d > 0
    rm = nz & (mx == r)
    gm = nz & (mx == g) & ~rm
    bm = nz & ~rm & ~gm
    h[rm] = (60.0 * ((g[rm] - b[rm]) / d[rm]) + 360.0) % 360.0
    h[gm] = 60.0 * ((b[gm] - r[gm]) / d[gm]) + 120.0
    h[bm] = 60.0 * ((r[bm] - g[bm]) / d[bm]) + 240.0
    fam = np.full(px.shape[0], "unknown", dtype=object)
    fam[(v < 0.12)] = "black"
    grey = (v >= 0.12) & (s < 0.15)
    fam[grey & (v > 0.85)] = "white"
    fam[grey & (v <= 0.85)] = "grey"
    chroma = (v >= 0.12) & (s >= 0.15)
    bands = [("red", 345, 15), ("orange", 15, 45), ("yellow", 45, 70), ("green", 70, 170),
             ("cyan", 170, 200), ("blue", 200, 260), ("purple", 260, 300), ("magenta", 300, 345)]
    for name, lo, hi in bands:
        if lo > hi:
            m = chroma & ((h >= lo) | (h < hi))
        else:
            m = chroma & (h >= lo) & (h < hi)
        fam[m] = name
    names, counts = np.unique(fam.astype(str), return_counts=True)
    total = float(px.shape[0])
    shares = {str(n): _r(c / total) for n, c in zip(names, counts)}
    return dict(sorted(shares.items(), key=lambda kv: -kv[1]))


def _thirds(mask: np.ndarray) -> Dict[str, float]:
    cells = grid_occupancy(mask, 3, 3)
    total = cells.sum()
    if total <= 0:
        return {}
    names = [["upper-left", "upper-centre", "upper-right"],
             ["middle-left", "centre", "middle-right"],
             ["lower-left", "lower-centre", "lower-right"]]
    out = {}
    for r in range(3):
        for c in range(3):
            out[names[r][c]] = _r(cells[r, c] / total)
    return out


def _describe_offset(dx: float, dy: float) -> str:
    parts = []
    if abs(dx) >= 0.04:
        parts.append(("right" if dx > 0 else "left") + " by %d%%" % round(abs(dx) * 100))
    if abs(dy) >= 0.04:
        parts.append(("down" if dy > 0 else "up") + " by %d%%" % round(abs(dy) * 100))
    return "centred" if not parts else "shifted " + " and ".join(parts)


# --------------------------------------------------------------------------
# main analysis
# --------------------------------------------------------------------------

def analyze_render(src: ImageSource, grid_cols: int = 48, grid_rows: int = 20,
                   palette_k: int = 8, bg_rgb: Optional[Sequence[int]] = None,
                   tol: int = 18) -> dict:
    """Describe a rendered frame numerically.

    ``bg_rgb`` can be supplied by the caller (the viewer knows its background
    colour exactly), otherwise it is estimated from the image border.
    """
    rgba = load_rgba(src)
    h, w = rgba.shape[:2]
    est_bg, transparent, agreement = estimate_background(rgba)
    bg = tuple(int(x) for x in bg_rgb[:3]) if (bg_rgb is not None and not transparent) else est_bg
    fg = foreground_mask(rgba, bg, transparent, tol=tol)
    n_fg = int(fg.sum())
    coverage = n_fg / float(h * w)
    bb = bbox_of(fg)
    cen = _centroid(fg)
    edges = _edge_coverage(fg)
    clipped = [side for side, frac in edges.items() if frac >= 0.02]
    touching = [side for side, frac in edges.items() if 0 < frac < 0.02]

    fg_rgb = rgba[fg][:, :3]
    pal = palette(fg_rgb, k=palette_k)
    hues = hue_shares(fg_rgb)

    luma = _luma(rgba[:, :, :3])
    if n_fg:
        fg_mean_rgb = fg_rgb.mean(axis=0)
        fg_luma_mean = float(luma[fg].mean())
        fg_luma_std = float(luma[fg].std())
        contrast = contrast_ratio(fg_mean_rgb, bg)
        # edge density: proportion of foreground pixels on a strong gradient
        gy, gx = np.gradient(luma)
        grad = np.hypot(gx, gy)
        edge_density = float((grad[fg] > 24).mean())
    else:
        fg_mean_rgb = np.array(bg, dtype=float)
        fg_luma_mean = fg_luma_std = 0.0
        contrast = 1.0
        edge_density = 0.0

    cells = grid_occupancy(fg, grid_cols, grid_rows)
    amap = ascii_map(cells)
    thirds = _thirds(fg)

    bg_name, _ = nearest_color_name(bg)
    result = {
        "width": w,
        "height": h,
        "background": {
            "rgb": list(bg),
            "hex": hex_color(bg),
            "name": "transparent" if transparent else bg_name,
            "transparent": transparent,
            "uniform_border": _r(agreement),
        },
        "coverage": _r(coverage, 4),
        "empty": coverage < 0.0005,
        "bbox": _bbox_dict(bb, w, h),
        "centroid_frac": None if cen is None else [_r(cen[0] / w), _r(cen[1] / h)],
        "offset_from_centre_frac": None if cen is None else [_r(cen[0] / w - 0.5), _r(cen[1] / h - 0.5)],
        "margins_frac": None if bb is None else {
            "left": _r(bb[0] / w), "right": _r(1 - bb[2] / w),
            "top": _r(bb[1] / h), "bottom": _r(1 - bb[3] / h),
        },
        "edge_coverage": edges,
        "clipped_edges": clipped,
        "touching_edges": touching,
        "thirds": thirds,
        "palette": pal,
        "hue_shares": hues,
        "luminance": {
            "fg_mean": _r(fg_luma_mean, 1),
            "fg_std": _r(fg_luma_std, 1),
            "bg": _r(float(_luma(np.array(bg, dtype=np.uint8)[None, :])[0]), 1),
            "contrast_ratio_fg_bg": _r(contrast, 2),
        },
        "edge_density": _r(edge_density),
        "ascii": amap,
        "ascii_legend": "each char = one cell of a %dx%d grid; ' ' empty .. '@' full" % (grid_cols, grid_rows),
    }
    result["summary"] = summarize_render(result)
    return result


def summarize_render(a: dict) -> str:
    if a.get("empty"):
        return "Frame is empty: nothing but the %s background is visible." % a["background"]["name"]
    parts = []
    bb = a["bbox"]
    parts.append("Content covers %.0f%% of the %dx%d frame (bounding box %.0f%% x %.0f%%)." % (
        a["coverage"] * 100, a["width"], a["height"], bb["width_frac"] * 100, bb["height_frac"] * 100))
    dx, dy = a["offset_from_centre_frac"]
    parts.append("Centroid is " + _describe_offset(dx, dy) + ".")
    if a["clipped_edges"]:
        parts.append("Content is cut off at the %s edge%s (%s)." % (
            ", ".join(a["clipped_edges"]), "s" if len(a["clipped_edges"]) > 1 else "",
            ", ".join("%s %.0f%%" % (e, a["edge_coverage"][e] * 100) for e in a["clipped_edges"])))
    elif a["touching_edges"]:
        parts.append("Content just touches the %s edge%s." % (
            ", ".join(a["touching_edges"]), "s" if len(a["touching_edges"]) > 1 else ""))
    else:
        m = a["margins_frac"]
        parts.append("Free margins: left %.0f%%, right %.0f%%, top %.0f%%, bottom %.0f%%." % (
            m["left"] * 100, m["right"] * 100, m["top"] * 100, m["bottom"] * 100))
    if a["palette"]:
        top = a["palette"][:4]
        parts.append("Dominant colours: " + ", ".join("%s %.0f%%" % (p["name"], p["share"] * 100) for p in top) + ".")
    lum = a["luminance"]
    parts.append("Background %s; foreground/background contrast %.1f:1." % (a["background"]["name"], lum["contrast_ratio_fg_bg"]))
    if lum["contrast_ratio_fg_bg"] < 1.6:
        parts.append("Warning: low contrast, the content may be hard to distinguish from the background.")
    if a["coverage"] < 0.05:
        parts.append("Content is small in the frame; consider zooming in.")
    return " ".join(parts)


# --------------------------------------------------------------------------
# segmentation (ID pass)
# --------------------------------------------------------------------------

def _parse_hex(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def id_colors(n: int) -> List[Tuple[int, int, int]]:
    """``n`` well-separated flat colours for an ID pass (never pure black)."""
    out = []
    # spread over a 3D lattice of 6 levels per channel, skipping (0,0,0)
    levels = [40, 80, 120, 160, 200, 240]
    import itertools
    combos = list(itertools.product(levels, repeat=3))
    # order combos so successive colours are far apart
    combos.sort(key=lambda c: (sum(c) % 3, c))
    stride = max(1, len(combos) // max(n, 1))
    i = 0
    while len(out) < n:
        c = combos[(i * stride) % len(combos)]
        if c not in out:
            out.append(c)
        i += 1
        if i > len(combos) * 4:
            break
    while len(out) < n:  # fall back to any remaining colour
        for c in combos:
            if c not in out:
                out.append(c)
                if len(out) >= n:
                    break
    return out[:n]


def analyze_segmentation(src: ImageSource, legend: Dict[str, str], bg_rgb: Optional[Sequence[int]] = None,
                         grid_cols: int = 48, grid_rows: int = 20, tol: int = 40) -> dict:
    """Per-item visibility statistics from an ID-pass render.

    ``legend`` maps hex colour -> item label.  Pixels are matched exactly
    first; anti-aliased border pixels are assigned to the nearest legend
    colour within ``tol`` (max channel difference), the rest count as
    "unassigned".
    """
    rgba = load_rgba(src)
    h, w = rgba.shape[:2]
    est_bg, transparent, _ = estimate_background(rgba)
    bg = tuple(int(x) for x in bg_rgb[:3]) if (bg_rgb is not None and not transparent) else est_bg
    labels = list(legend.values())
    colors = np.array([_parse_hex(k) for k in legend.keys()], dtype=np.int16)
    rgb = rgba[:, :, :3].astype(np.int16)
    flat = rgb.reshape(-1, 3)

    assign = np.full(flat.shape[0], -1, dtype=np.int32)  # -1 = background, -2 = unassigned
    if transparent:
        is_bg = rgba[:, :, 3].reshape(-1) <= 8
    else:
        is_bg = (np.abs(flat - np.array(bg, dtype=np.int16)).max(axis=1) <= 12)
    fg_idx = np.flatnonzero(~is_bg)
    if fg_idx.size and colors.shape[0]:
        fgpx = flat[fg_idx]
        # nearest legend colour, chunked to bound memory
        best = np.empty(fgpx.shape[0], dtype=np.int32)
        bestd = np.empty(fgpx.shape[0], dtype=np.int16)
        chunk = max(1, 2_000_000 // max(colors.shape[0], 1))
        for s in range(0, fgpx.shape[0], chunk):
            blk = fgpx[s:s + chunk]
            d = np.abs(blk[:, None, :] - colors[None, :, :]).max(axis=2)
            bi = d.argmin(axis=1)
            best[s:s + chunk] = bi
            bestd[s:s + chunk] = d[np.arange(blk.shape[0]), bi]
        ok = bestd <= tol
        assign[fg_idx[ok]] = best[ok]
        assign[fg_idx[~ok]] = -2
    assign = assign.reshape(h, w)

    n_fg = int((assign != -1).sum())
    items = []
    dominant = np.full((grid_rows, grid_cols), -1, dtype=np.int32)
    dom_score = np.zeros((grid_rows, grid_cols), dtype=np.float64)
    for i, label in enumerate(labels):
        m = assign == i
        cnt = int(m.sum())
        if cnt == 0:
            items.append({"label": label, "visible": False, "pixels": 0, "share_of_frame": 0.0,
                          "share_of_content": 0.0, "bbox": None, "centroid_frac": None,
                          "touches_edges": []})
            continue
        bb = bbox_of(m)
        cen = _centroid(m)
        edges = _edge_coverage(m)
        cells = grid_occupancy(m, grid_cols, grid_rows)
        better = cells > dom_score
        dominant[better] = i
        dom_score[better] = cells[better]
        items.append({
            "label": label,
            "visible": True,
            "pixels": cnt,
            "share_of_frame": _r(cnt / float(h * w), 4),
            "share_of_content": _r(cnt / float(max(n_fg, 1)), 4),
            "bbox": _bbox_dict(bb, w, h),
            "centroid_frac": [_r(cen[0] / w), _r(cen[1] / h)],
            "touches_edges": [s for s, f in edges.items() if f > 0],
            "edge_coverage": edges,
        })

    # 2D contacts between items (4-neighbourhood, ignoring background)
    contacts: Dict[Tuple[int, int], int] = {}
    for a, b in ((assign[:, :-1], assign[:, 1:]), (assign[:-1, :], assign[1:, :])):
        diff = (a != b) & (a >= 0) & (b >= 0)
        if diff.any():
            pa = a[diff]
            pb = b[diff]
            lo = np.minimum(pa, pb)
            hi = np.maximum(pa, pb)
            keys = lo.astype(np.int64) * 100000 + hi
            vals, counts = np.unique(keys, return_counts=True)
            for k, c in zip(vals, counts):
                pair = (int(k // 100000), int(k % 100000))
                contacts[pair] = contacts.get(pair, 0) + int(c)
    contact_list = sorted(
        ({"a": labels[p[0]], "b": labels[p[1]], "boundary_px": c} for p, c in contacts.items()),
        key=lambda d: -d["boundary_px"])[:40]

    # ASCII layout map: letter of dominant item per cell
    letters = {i: LABEL_CHARS[i % len(LABEL_CHARS)] for i in range(len(labels))}
    lines = []
    for r in range(grid_rows):
        line = []
        for c in range(grid_cols):
            i = int(dominant[r, c])
            line.append("." if i < 0 or dom_score[r, c] < 0.02 else letters[i])
        lines.append("".join(line))
    legend_lines = {letters[i]: labels[i] for i in range(len(labels)) if items[i]["visible"]}

    unassigned = int((assign == -2).sum())
    visible_items = [it for it in items if it["visible"]]
    visible_items.sort(key=lambda it: -it["pixels"])
    summary = "%d of %d items visible." % (len(visible_items), len(items))
    if visible_items:
        summary += " Largest on screen: " + ", ".join(
            "%s (%.0f%% of content)" % (it["label"], it["share_of_content"] * 100) for it in visible_items[:5]) + "."
    hidden = [it["label"] for it in items if not it["visible"]]
    if hidden:
        summary += " Not visible at all: " + ", ".join(hidden[:10]) + ("..." if len(hidden) > 10 else "") + "."
    cut = [it["label"] for it in visible_items if it["touches_edges"]]
    if cut:
        summary += " Touching the frame edge: " + ", ".join(cut[:10]) + "."

    return {
        "width": w,
        "height": h,
        "content_pixels": n_fg,
        "coverage": _r(n_fg / float(h * w), 4),
        "unassigned_pixels": unassigned,
        "items": items,
        "contacts_2d": contact_list,
        "ascii": lines,
        "ascii_legend": legend_lines,
        "summary": summary,
    }


def occlusion_report(composite: dict, solo_counts: Dict[str, int]) -> List[dict]:
    """Combine composite visibility with per-item solo renders.

    ``solo_counts`` maps label -> pixels the item covers when drawn alone.
    """
    out = []
    by_label = {it["label"]: it for it in composite["items"]}
    for label, solo in solo_counts.items():
        it = by_label.get(label)
        vis = it["pixels"] if it else 0
        if solo <= 0:
            out.append({"label": label, "solo_pixels": 0, "visible_pixels": vis,
                        "visible_fraction": None, "note": "item renders nothing on its own"})
            continue
        frac = min(1.0, vis / float(solo))
        out.append({"label": label, "solo_pixels": int(solo), "visible_pixels": int(vis),
                    "visible_fraction": _r(frac), "occluded_fraction": _r(1 - frac)})
    out.sort(key=lambda d: (d["visible_fraction"] if d["visible_fraction"] is not None else 2))
    return out


# --------------------------------------------------------------------------
# diff
# --------------------------------------------------------------------------

def diff_images(a_src: ImageSource, b_src: ImageSource, tol: int = 16,
                grid_cols: int = 48, grid_rows: int = 20) -> dict:
    a = load_rgba(a_src)
    b = load_rgba(b_src)
    resized = False
    if a.shape != b.shape:
        b = np.array(Image.fromarray(b).resize((a.shape[1], a.shape[0]), Image.BILINEAR), dtype=np.uint8)
        resized = True
    d = np.abs(a.astype(np.int16) - b.astype(np.int16)).max(axis=2)
    changed = d > tol
    n = int(changed.sum())
    h, w = changed.shape
    bb = bbox_of(changed)
    cells = grid_occupancy(changed, grid_cols, grid_rows)
    frac = n / float(h * w)
    if n == 0:
        summary = "Images are identical within tolerance."
    else:
        summary = "%.1f%% of pixels changed" % (frac * 100)
        if bb:
            summary += ", inside a region spanning %.0f%% x %.0f%% of the frame" % (
                (bb[2] - bb[0]) / w * 100, (bb[3] - bb[1]) / h * 100)
        summary += "."
    return {
        "width": w, "height": h, "resized_b": resized,
        "changed_fraction": _r(frac, 4),
        "changed_pixels": n,
        "mean_abs_diff": _r(float(d.mean()), 2),
        "bbox": _bbox_dict(bb, w, h),
        "ascii": ascii_map(cells),
        "summary": summary,
    }
