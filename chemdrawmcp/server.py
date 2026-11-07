"""chemdraw-mcp: MCP server that gives agents full access to ChemDraw.

Two engines work together:

* the live ChemDraw application (macOS via AppleScript/JXA, Windows via COM) for everything ChemDraw is good at:
  publication-quality rendering in journal styles, Struct=Name / Name=Struct, chemical properties, format
  conversion, NMR prediction, clean-up, its ~1600 named commands and menus;
* a headless CDXML engine (RDKit + our own CDXML writer/reader) that lays out structures, reactions and multi-step
  schemes with reagents, conditions, labels and arrows, and that keeps working when ChemDraw is not installed.
"""
from __future__ import annotations

import contextlib
import datetime as _dt
import json
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import __version__
from . import cdxml as cdx
from . import chem
from . import discovery
from . import imageanalysis as ia
from . import layout
from . import styles as styles_mod
from .backends import ChemDrawError, ChemDrawNotAvailable, NullBackend, get_backend
from .backends.base import EXPORT_FORMATS, INSERT_KINDS, ChemDrawBackend, normalize_format
from .mcpcompat import FastMCP, Image

try:  # keep RDKit's warnings off the server's stderr
    from rdkit import RDLogger as _RDLogger
    _RDLogger.DisableLog("rdApp.*")
except Exception:  # pragma: no cover
    pass

INSTRUCTIONS = """ChemDraw MCP server (chemdrawmcp). You control the ChemDraw desktop application and a headless CDXML engine.

Typical flow:
1. `status` - is ChemDraw installed/running, which documents are open, where files go.  `launch` starts it.
2. Draw: `draw_structures` (grid of molecules with labels), `draw_reaction` (reactants -> products with reagents above
   and conditions below the arrow), `draw_scheme` (multi-step synthesis, wraps into rows).  Inputs are SMILES, InChI,
   molblocks, CDXML, or chemical names (resolved by an offline dictionary of ~3,100 names, then OPSIN if installed,
   then ChemDraw's Name=Struct, then PubChem when CHEMDRAWMCP_ONLINE=1).  Each call
   opens a new ChemDraw document (or pastes into `document`), renders it and returns a numerical description plus a
   thumbnail you can look at.  Journal styles: acs (default), nature, rsc, wiley, science_of_synthesis, ... (`styles`).
3. Look: `render` (PNG/SVG/PDF/EPS/TIFF... export + analysis), `screenshot`, `read_document` (every molecule as SMILES
   with formula/MW, captions, arrows, reaction steps as reaction SMILES).
4. Chemistry: `name_to_structure`, `structure_to_name` (IUPAC name from ChemDraw), `properties` (ChemDraw + RDKit:
   formula, MW, exact mass, elemental analysis, logP, TPSA, InChIKey...), `analyze_reaction` (balance, atom economy),
   `convert` (any format <-> any format: smiles, inchi, mol, cdxml, cdx, sdf, rxn, cml, png, svg, pdf...),
   `predict_nmr`, `clean_structure`, `check_structure`.
5. Raw control: `insert` (paste SMILES/name/CDXML/MOL/InChI/HELM/text into a document), `selection`, `select`,
   `commands` + `run_command` (any of ChemDraw's ~1600 named commands, e.g. cleanStructure, mapReactionAtoms,
   autoNumberReaction, toggleShowingAtomNumber, flipObjectsHorizontal), `menus` + `do_menu`, `document_settings`,
   `save_document`, `new_document`, `open_document`, `close_document`, `dismiss_dialogs`.

Documents are addressed by name (as in the window title) or by window id (number); None = front document.
ChemDraw renames untitled documents after their first edit, so prefer the window id returned by each tool.
Files are written to the renders directory reported by `status` (CHEMDRAWMCP_OUTPUT_DIR).
Without ChemDraw (Linux, or not installed) the draw_* tools still produce CDXML files and RDKit renders.
"""

mcp = FastMCP("chemdraw", instructions=INSTRUCTIONS)

MAX_TEXT = 8000
_BACKEND: Optional[ChemDrawBackend] = None
_NAME_CACHE: Dict[str, str] = {}
_NAME_SOURCE: Dict[str, str] = {}
DEFAULT_STYLE = os.environ.get("CHEMDRAWMCP_STYLE", "acs")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def dumps(obj: Any) -> str:
    return json.dumps(obj, indent=1, ensure_ascii=False, default=str)


def _backend() -> ChemDrawBackend:
    global _BACKEND
    if _BACKEND is None:
        _BACKEND = get_backend()
    return _BACKEND


def _live(autolaunch: bool = True) -> ChemDrawBackend:
    """The backend, with ChemDraw running (launched on demand unless CHEMDRAWMCP_AUTOLAUNCH=0)."""
    b = _backend()
    if isinstance(b, NullBackend):
        raise ChemDrawNotAvailable(b.reason)
    if not b.is_running():
        if autolaunch and os.environ.get("CHEMDRAWMCP_AUTOLAUNCH", "1") != "0":
            b.launch()
        else:
            raise ChemDrawError("ChemDraw is not running; call `launch` first")
    return b


def _has_live() -> bool:
    b = _backend()
    return not isinstance(b, NullBackend) and b.available()


def _doc(document: Any) -> Any:
    """Normalise a document argument: None (front), a name, or a window id (int or '#123')."""
    if document is None or document == "" or document == "front":
        return None
    if isinstance(document, int):
        return document
    s = str(document).strip()
    if re.match(r"^#?\d+$", s):
        return int(s.lstrip("#"))
    return s


_SEQ = [0]


def _timestamp() -> str:
    _SEQ[0] += 1
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3] + "_%02d" % (_SEQ[0] % 100)


def _slug(text: str, n: int = 40) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return (s[:n] or "item").lower()


def _out_path(tag: str, ext: str, path: Optional[str] = None) -> str:
    if path:
        p = os.path.abspath(os.path.expanduser(path))
        os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
        return p
    return os.path.join(discovery.output_dir(), "%s_%s.%s" % (_slug(tag), _timestamp(), ext.lstrip(".")))


def _in_path(path: str) -> str:
    """An existing input file: `~` expanded and made absolute (as every other path-taking tool does)."""
    p = os.path.abspath(os.path.expanduser(str(path).strip()))
    if not os.path.isfile(p):
        raise chem.ChemError("file not found: %s" % p)
    return p


def _trim(text: str, limit: int = MAX_TEXT) -> str:
    if len(text) <= limit:
        return text
    return text[: limit // 2] + "\n... [%d chars omitted] ...\n" % (len(text) - limit) + text[-limit // 2:]


def _err(e: Exception, **extra) -> str:
    out = {"ok": False, "error": str(e), "error_type": type(e).__name__}
    if isinstance(e, ChemDrawNotAvailable):
        out["hint"] = "ChemDraw is not available here; headless tools (draw_* with open_in_chemdraw=false, properties, convert, analyze_reaction) still work."
    out.update(extra)
    return dumps(out)


def _doc_info(b: ChemDrawBackend, name: Any) -> dict:
    """Name + window id + counts for a document (robust to ChemDraw renaming it)."""
    try:
        for d in b.documents():
            if d.name == name or d.window_id == name or (isinstance(name, str) and name.isdigit() and d.window_id == int(name)):
                return d.to_dict()
    except ChemDrawError:
        pass
    return {"name": name}


@contextlib.contextmanager
def _temp_document(b: ChemDrawBackend, style: Optional[str] = None):
    """A scratch document that is closed afterwards; yields its window id (stable across renames)."""
    name = b.new_document(style or DEFAULT_STYLE)
    wid = None
    try:
        for d in b.documents():
            if d.name == name:
                wid = d.window_id
    except ChemDrawError:
        pass
    handle = wid if wid is not None else name
    try:
        yield handle
    finally:
        try:
            b.close_document(handle, save=False)
        except ChemDrawError:
            pass


# -- chemistry input -------------------------------------------------------------

def _dictionary_lookup(name: str) -> Optional[str]:
    """chem's built-in dictionary of common names alone (instant; no OPSIN, no network)."""
    raw = name.strip()
    smi = chem._CASE_SENSITIVE.get(raw)
    if smi is None:
        key = chem._normalize_name(raw)
        smi = chem._resolve_offline(key) if key else None
    return smi


def _resolve_name(name: str) -> Optional[str]:
    """Chemical name -> SMILES.  Order: 1. the curated offline dictionary (instant), 2. OPSIN when installed
    (offline IUPAC parser), 3. ChemDraw's Name=Struct when ChemDraw is available (~5 s per lookup), 4. PubChem when
    CHEMDRAWMCP_ONLINE=1.  _NAME_SOURCE records which one answered ("offline dictionary" | "ChemDraw Name=Struct" |
    "OPSIN" | "PubChem")."""
    key = name.strip().lower()
    if key in _NAME_CACHE:
        return _NAME_CACHE[key]
    raw = name.strip()
    smi = _dictionary_lookup(raw)
    source = "offline dictionary" if smi else None
    if not smi:
        try:
            smi = chem._opsin(raw)       # offline IUPAC parser when installed: fast, no ChemDraw round trip
        except Exception:
            smi = None
        source = "OPSIN" if smi else None
    if not smi and _has_live():
        try:
            b = _live()
            with _temp_document(b) as h:
                r = b.insert(h, raw, "name")
                if r.get("atoms_added", 0) > 0:
                    smi = b.selection_as(h, "smiles").strip() or None
                    source = "ChemDraw Name=Struct" if smi else None
        except ChemDrawError:
            smi = None
    if not smi and os.environ.get("CHEMDRAWMCP_ONLINE", "").strip() == "1":
        try:
            smi = chem._pubchem(raw)
        except Exception:
            smi = None
        source = "PubChem" if smi else None
    if smi:
        _NAME_CACHE[key] = smi
        _NAME_SOURCE[key] = source or "unknown"
    return smi


def _cdx_text(data: str) -> Optional[str]:
    """If ``data`` is a path to a binary .cdx file, return its CDXML text (via the CDX reader), else None."""
    p = os.path.expanduser(data.strip())
    if len(p) < 4096 and os.path.isfile(p):
        with open(p, "rb") as fh:
            raw = fh.read()
        if raw[:8] == b"VjCD0100":
            from . import cdx as cdxbin
            return cdxbin.cdx_to_cdxml(raw)
    return None


def _mol(data: str, fmt: Optional[str] = None):
    """Parse one structure (SMILES/InChI/molblock/CDXML/CDX file/path/name) into a depiction-ready RDKit mol."""
    data = data.strip()
    if fmt in (None, "cdx"):
        cdx_text = _cdx_text(data)
        if cdx_text is not None:
            data, fmt = cdx_text, "cdxml"
    f = fmt or chem.detect_format(data)
    if f == "name":
        smi = _resolve_name(data)
        if not smi:
            raise chem.ChemError("could not interpret %r as SMILES, InChI, molfile, CDXML or a chemical name" % data[:80])
        m = chem.mol_from_input(smi, "smiles")
        m.SetProp("_Name", data)
        return m
    return chem.mol_from_input(data, f)


def _looks_like_reaction(text: str, fmt: Optional[str] = None) -> bool:
    """Reaction SMILES, or arrow text in any form chem.parse_reaction accepts ('->', '-->', '>>', '=>', ' = ', ' > ',
    unicode arrows, '--[conditions]-->').  Molblocks, SDF, CDXML and InChI are never reactions here, and the dative
    '->' inside SMILES does not count."""
    f = fmt or chem.detect_format(text)
    if f == "reaction_smiles":
        return True
    return f in ("name", "smiles") and chem._ARROW_RE.search(chem._strip_dative(text)) is not None


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


def _component(spec: Any) -> layout.Component:
    if isinstance(spec, dict):
        mol = _mol(str(spec.get("structure") or spec.get("smiles") or spec.get("mol")))
        la = spec.get("label_atoms")
        return layout.Component(mol, spec.get("label"), spec.get("color"), None, _tiny(mol) if la is None else bool(la))
    mol = _mol(str(spec))
    return layout.Component(mol, label_atoms=_tiny(mol))


def _reaction_step(reaction: Optional[str], reactants: Optional[Sequence[Any]], products: Optional[Sequence[Any]],
                   above: Optional[Sequence[str]], below: Optional[Sequence[str]], arrow: str,
                   reactant_labels: Optional[Sequence[str]] = None, product_labels: Optional[Sequence[str]] = None,
                   atom_map: bool = False) -> layout.Step:
    above = list(above or [])
    below = list(below or [])
    if reaction:
        # names inside "A + B -> C" syntax are resolved through our resolver first
        text = reaction
        rxn = chem.parse_reaction(_resolve_names_in_reaction(text))
        r_comps = [layout.Component(m, label_atoms=_tiny(m)) for m in rxn.reactants]
        p_comps = [layout.Component(m, label_atoms=_tiny(m)) for m in rxn.products]
        if rxn.agents and not above:
            above = [_agent_text(m) for m in rxn.agents]
        cond = getattr(rxn, "conditions", "") or ""
        if cond and not above:
            above = [cond]
        elif cond:
            below = below + [cond]
    else:
        if not reactants or not products:
            raise chem.ChemError("give either `reaction` or both `reactants` and `products`")
        r_comps = [_component(s) for s in reactants]
        p_comps = [_component(s) for s in products]
    for comps, labels in ((r_comps, reactant_labels), (p_comps, product_labels)):
        if labels:
            for c, lab in zip(comps, labels):
                if lab:
                    c.label = lab
    return layout.Step(r_comps, p_comps, above, below, arrow, atom_map=atom_map)


def _agent_text(mol) -> str:
    try:
        from rdkit.Chem import rdMolDescriptors
        if mol.HasProp("_Name") and mol.GetProp("_Name"):
            return mol.GetProp("_Name")
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return "reagent"


def _resolve_names_in_reaction(text: str) -> str:
    """Replace chemical names in arrow-syntax reactions with SMILES so parse_reaction never needs the network.
    Uses the parser's own arrow set (chem._ARROW_RE) and keeps stoichiometric coefficients ("2 water" -> "2 O") so
    that parse_reaction still expands them."""
    s = text.strip()
    if re.match(r"^[^\s>]+>[^\s>]*>[^\s>]+$", s):
        return text   # proper reaction SMILES
    m = chem._ARROW_RE.search(s)
    if m is None or chem._ARROW_RE.search(s[m.end():]):
        return text   # no arrow, or a multi-step reaction: parse_reaction reports it

    def fix(side: str) -> str:
        parts = [p.strip() for p in re.split(r"\s\+\s|^\+\s|\s\+$", side) if p.strip()]
        out = []
        for p in parts:
            cm = chem._COEFF_RE.match(p)
            coeff, body = (cm.group(1) + " ", cm.group(2)) if cm else ("", p)
            if chem.detect_format(body) == "name":
                smi = _resolve_name(body)
                if smi:
                    out.append(coeff + smi)
                    continue
            out.append(p)
        return " + ".join(out)
    return fix(s[: m.start()]) + m.group(0) + fix(s[m.end():])


# -- materialising CDXML in ChemDraw and rendering ----------------------------------

def _materialize(doc: cdx.CDXMLDocument, document: Any, open_in_chemdraw: bool, tag: str) -> Tuple[dict, str]:
    xml = doc.to_xml()
    path = _out_path(tag, "cdxml")
    with open(path, "w", encoding="utf-8") as f:
        f.write(xml)
    info: Dict[str, Any] = {"cdxml_path": path, "layout": doc.summary()}
    warnings = list(getattr(doc, "warnings", []) or [])
    if warnings:
        info["layout_warnings"] = warnings
    if open_in_chemdraw and _has_live():
        b = _live()
        if document is None:
            name = b.open_document(path)
            info.update({"document": name, "opened": True})
            info.update({k: v for k, v in _doc_info(b, name).items() if k in ("window_id", "atoms", "molecules", "captions", "arrows")})
        else:
            res = b.insert(_doc(document), xml, "cdxml")
            info.update({"document": res.get("document"), "inserted": True, "atoms_added": res.get("atoms_added"),
                         "note": "ChemDraw chooses the paste position; the pasted objects are selected so the user can move them"})
            info.update({k: v for k, v in _doc_info(b, res.get("document")).items() if k in ("window_id",)})
    else:
        info["opened"] = False
        info["note"] = ("ChemDraw not available: CDXML written to disk" if not _has_live() else "open_in_chemdraw=false: CDXML written to disk")
    return info, xml


def _analysis(png_bytes: bytes) -> dict:
    try:
        a = ia.analyze_render(png_bytes)
        a["summary"] = ia.summarize_render(a)
        for k in ("ascii_map", "grid"):
            a.pop(k, None)
        return a
    except Exception as e:  # pragma: no cover
        return {"error": "image analysis failed: %s" % e}


def _render_document(b: ChemDrawBackend, handle: Any, fmt: str, path: Optional[str], tag: str,
                     return_image: bool, width: Optional[int] = None) -> Tuple[dict, Optional[bytes]]:
    fmt = normalize_format(fmt)
    out = _out_path(tag, EXPORT_FORMATS[fmt]["ext"].lstrip("."), path)
    written = b.save_document(handle, out, fmt)
    info: Dict[str, Any] = {"path": written, "format": fmt, "rendered_by": "ChemDraw", "bytes": os.path.getsize(written)}
    png = None
    if fmt in ("png", "jpeg", "gif", "bmp", "tiff"):
        try:
            from PIL import Image as PILImage
            import io
            im = PILImage.open(written)
            info["size_px"] = list(im.size)
            info["dpi"] = [round(v) for v in im.info.get("dpi", (300, 300))]
            if width and im.size[0] != width:
                ratio = width / float(im.size[0])
                mode = "RGB" if fmt in ("jpeg", "bmp") else "RGBA"
                im = im.convert(mode).resize((width, max(1, int(round(im.size[1] * ratio)))), PILImage.LANCZOS)
                im.save(written)
                info["size_px"] = list(im.size)
                info["resized_to_width"] = width
            buf = io.BytesIO()
            rgba = im.convert("RGBA")
            if rgba.getextrema()[3][0] < 255:
                # ChemDraw exports a transparent background: analyse and show the picture on white, like a viewer does
                flat = PILImage.new("RGBA", rgba.size, (255, 255, 255, 255))
                flat.paste(rgba, (0, 0), rgba)
                rgba = flat
            rgba.save(buf, format="PNG")
            png = buf.getvalue()
            info["analysis"] = _analysis(png)
        except Exception as e:  # pragma: no cover
            info["analysis_error"] = str(e)
    return info, (png if return_image else None)


def _headless_png(doc_or_obj, tag: str, path: Optional[str]) -> Tuple[dict, bytes]:
    png = chem.depict_png(doc_or_obj, width=900, height=450)
    out = _out_path(tag, "png", path)
    with open(out, "wb") as f:
        f.write(png)
    return {"path": out, "format": "png", "rendered_by": "RDKit" + ("" if _has_live() else " (ChemDraw not available)"),
            "analysis": _analysis(png)}, png


def _with_image(text: str, png: Optional[bytes]):
    if png:
        return [text, Image(data=ia.thumbnail_png(png, 1024), format="png")]
    return text


def _grid_png(mols, labels) -> bytes:
    from rdkit.Chem import Draw
    import io
    legends = [((labels[i] if labels and i < len(labels) else "") or "") for i in range(len(mols))]
    img = Draw.MolsToGridImage(mols, molsPerRow=min(4, max(1, len(mols))), subImgSize=(300, 250), legends=legends)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# status / lifecycle
# ---------------------------------------------------------------------------

@mcp.tool()
def status() -> str:
    """ChemDraw availability and state: app path/version, running?, open documents (names, window ids, counts),
    pending dialogs, renders directory, RDKit version, available styles."""
    b = _backend()
    out: Dict[str, Any] = {"server_version": __version__, "platform": sys.platform, "backend": b.platform,
                           "chemdraw_app": b.executable(), "chemdraw_available": _has_live(),
                           "renders_dir": discovery.output_dir(), "home": discovery.home_dir(),
                           "rdkit": chem.rdkit_available(), "styles": styles_mod.style_names(), "default_style": DEFAULT_STYLE}
    try:
        out["chemdraw_version"] = discovery.chemdraw_version(b.executable()) if b.executable() else None
    except Exception:
        pass
    if isinstance(b, NullBackend):
        out["headless_reason"] = b.reason
        return dumps(out)
    out["running"] = b.is_running()
    if out["running"]:
        try:
            if hasattr(b, "snapshot"):
                snap = b.snapshot()
                out["version"] = snap.get("version")
                out["documents"] = snap.get("documents", [])
                out["front_document"] = snap.get("front_document")
            else:
                out["version"] = b.version()
                out["documents"] = [d.to_dict() for d in b.documents()]
                out["front_document"] = next((d["name"] for d in out["documents"] if d.get("front")), None)
            dl = b.pending_dialogs() if hasattr(b, "pending_dialogs") else []
            if dl:
                out["pending_dialogs"] = dl
                out["hint"] = "ChemDraw is showing a dialog; call dismiss_dialogs"
        except ChemDrawError as e:
            out["error"] = str(e)
    else:
        out["hint"] = "ChemDraw is installed but not running; tools launch it automatically (or call `launch`)"
    return dumps(out)


@mcp.tool()
def launch(wait_seconds: float = 90) -> str:
    """Start ChemDraw (if not running) and wait until it answers."""
    try:
        b = _backend()
        if isinstance(b, NullBackend):
            raise ChemDrawNotAvailable(b.reason)
        return dumps(b.launch(wait_seconds))
    except Exception as e:
        return _err(e)


@mcp.tool()
def quit_chemdraw(save: bool = False) -> str:
    """Quit ChemDraw (discarding unsaved changes unless save=true)."""
    try:
        _live(autolaunch=False).quit(save)
        return dumps({"ok": True})
    except Exception as e:
        return _err(e)


@mcp.tool()
def activate_chemdraw(document: Optional[str] = None) -> str:
    """Bring ChemDraw (and optionally a document) to the foreground so the user sees it.  Normally tools give
    keyboard focus back to the user's previous application; this one keeps ChemDraw in front."""
    try:
        b = _live()
        b.activate()
        if document is not None:
            b.js("const d = docByName(ARGS.name); d.window().index = 1; return d.name();", name=_doc(document), __restore_focus=False) \
                if hasattr(b, "js") else b.activate_document(_doc(document))
        return dumps({"ok": True})
    except Exception as e:
        return _err(e)


# ---------------------------------------------------------------------------
# documents
# ---------------------------------------------------------------------------

@mcp.tool()
def documents() -> str:
    """Open ChemDraw documents: name, window id, path, modified flag, counts of atoms/bonds/molecules/captions/arrows."""
    try:
        return dumps({"documents": [d.to_dict() for d in _live().documents()]})
    except Exception as e:
        return _err(e)


@mcp.tool()
def new_document(style: str = "acs") -> str:
    """Create a new ChemDraw document from a stationery/style: acs, nature, rsc, wiley, science_of_synthesis,
    j_het_chem, helv_chim_acta, adv_synth_catal, phytomedicine, a4, slide, dark, default... (see `styles`), or a .cds path."""
    try:
        b = _live()
        name = b.new_document(style)
        return dumps({"ok": True, **_doc_info(b, name)})
    except Exception as e:
        return _err(e)


@mcp.tool()
def open_document(path: str) -> str:
    """Open a file in ChemDraw (cdxml, cdx, mol, sdf, rxn, cml, skc, cif, jdx spectra, pdf, png ...)."""
    try:
        b = _live()
        name = b.open_document(path)
        return dumps({"ok": True, **_doc_info(b, name)})
    except Exception as e:
        return _err(e)


@mcp.tool()
def close_document(document: Optional[str] = None, save: bool = False) -> str:
    """Close a document (front document by default), discarding changes unless save=true."""
    try:
        b = _live(autolaunch=False)
        b.close_document(_doc(document), save)
        return dumps({"ok": True, "remaining": [d.name for d in b.documents()]})
    except Exception as e:
        return _err(e)


@mcp.tool()
def save_document(path: str, document: Optional[str] = None, format: Optional[str] = None) -> str:
    """Save/export a document to any format ChemDraw writes: cdxml, cdx, mol (V2000), molv3000, sdf, sdfv3000, rxn (V2000,
    needs a reaction arrow), rxnv3000, rdf, cml, ct, png, svg, pdf, eps, tiff, gif, jpeg, bmp, 3mf.  Format defaults to the
    file extension."""
    try:
        b = _live(autolaunch=False)
        fmt = format or os.path.splitext(path)[1].lstrip(".") or "cdxml"
        written = b.save_document(_doc(document), path, fmt)
        return dumps({"ok": True, "path": written, "format": normalize_format(fmt), "bytes": os.path.getsize(written)})
    except Exception as e:
        return _err(e)


@mcp.tool()
def read_document(document: Optional[str] = None, include_cdxml: bool = False) -> str:
    """Structured content of a document: every molecule as SMILES with formula/MW and position, captions, arrows,
    graphics, reaction steps (as reaction SMILES with reagents), drawing settings.  Use it before editing."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        xml = b.document_cdxml(handle)
        out = cdx.read_cdxml(xml)
        out["document"] = _doc_info(b, handle if handle is not None else b.front_document())
        try:
            out["chemdraw_settings"] = b.document_settings(handle)
        except ChemDrawError:
            pass
        if include_cdxml:
            out["cdxml"] = _trim(xml)
        return dumps(out)
    except Exception as e:
        return _err(e)


@mcp.tool()
def read_file(path: str, include_cdxml: bool = False) -> str:
    """Parse a ChemDraw file on disk WITHOUT opening it: .cdxml or binary .cdx (own reader, works on Linux too).
    Same structured output as read_document: molecules as SMILES with formula/MW, captions, arrows, reaction steps."""
    try:
        p = os.path.abspath(os.path.expanduser(path))
        if not os.path.exists(p):
            raise chem.ChemError("file not found: %s" % p)
        with open(p, "rb") as fh:
            raw = fh.read()
        text = None
        if raw[:8] == b"VjCD0100" or p.lower().endswith(".cdx"):
            try:
                from . import cdx as cdxbin
                text = cdxbin.cdx_to_cdxml(raw)
            except ImportError:
                raise chem.ChemError("binary CDX reading is not available in this build; open the file in ChemDraw with open_document and use read_document")
        else:
            text = raw.decode("utf-8", "replace")
        out = cdx.read_cdxml(text)
        out["path"] = p
        if include_cdxml:
            out["cdxml"] = _trim(text)
        return dumps(out)
    except Exception as e:
        return _err(e)


@mcp.tool()
def document_settings(document: Optional[str] = None, settings: Optional[Dict[str, Any]] = None) -> str:
    """Read (and optionally change) a document's drawing settings in points: bond_length, line_width, bold_width,
    margin_width, hash_spacing, bond_spacing (%), chain_angle, label_size, caption_size, label_font, caption_font,
    magnification (%), fractional_widths, rulers_visible, crosshair_visible.  Example: {"bond_length": 14.4}."""
    try:
        b = _live(autolaunch=False)
        if settings:
            return dumps({"ok": True, "settings": b.set_document_settings(_doc(document), settings)})
        return dumps({"ok": True, "settings": b.document_settings(_doc(document))})
    except Exception as e:
        return _err(e)


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

@mcp.tool()
def draw_structures(structures: List[str], labels: Optional[List[Optional[str]]] = None, style: str = "acs",
                    columns: Optional[int] = None, colors: Optional[List[str]] = None, document: Optional[str] = None,
                    open_in_chemdraw: bool = True, render: bool = True, name: str = "structures",
                    bond_length: Optional[float] = None, return_image: bool = True):
    """Draw one or more molecules in a grid (optionally labelled) in a journal style and open the result in ChemDraw.

    structures: SMILES / InChI / molblock / CDXML / chemical names (e.g. ["aspirin", "CC(=O)O", "InChI=..."]).
    labels: caption under each structure (name, number, yield...).  colors: per structure ("red", "#0055ff"...).
    document: paste into this open document instead of creating a new one.  open_in_chemdraw=false only writes CDXML.
    Returns the document name/window id, the CDXML path, the layout, and (render=true) a PNG path + numerical
    description + thumbnail image.
    """
    try:
        mols = [_mol(s) for s in structures]
        if labels is None:
            labels = [m.GetProp("_Name") if m.HasProp("_Name") and chem.detect_format(s) == "name" else None for m, s in zip(mols, structures)]
        doc = layout.build_structures(mols, labels=labels, style=style, columns=columns, colors=colors, bond_length=bond_length, name=name,
                                      label_atoms=[_tiny(m) for m in mols])
        info, _xml = _materialize(doc, document, open_in_chemdraw, name)
        info["molecules"] = [chem.describe(m) for m in mols]
        png = None
        if render:
            if (info.get("opened") or info.get("inserted")) and _has_live():
                rinfo, png = _render_document(_live(), info.get("window_id") or info["document"], "png", None, name, return_image)
            else:
                png = _grid_png(mols, labels)
                out = _out_path(name, "png")
                with open(out, "wb") as f:
                    f.write(png)
                rinfo = {"path": out, "rendered_by": "RDKit", "analysis": _analysis(png)}
            info["render"] = rinfo
        return _with_image(dumps(info), png if return_image else None)
    except Exception as e:
        return _err(e)


@mcp.tool()
def draw_reaction(reaction: Optional[str] = None, reactants: Optional[List[str]] = None, products: Optional[List[str]] = None,
                  reagents_above: Optional[List[str]] = None, conditions_below: Optional[List[str]] = None,
                  arrow: str = "full", reactant_labels: Optional[List[Optional[str]]] = None, product_labels: Optional[List[Optional[str]]] = None,
                  style: str = "acs", document: Optional[str] = None, open_in_chemdraw: bool = True, render: bool = True,
                  atom_map: bool = False, name: str = "reaction", bond_length: Optional[float] = None, return_image: bool = True):
    """Draw a reaction: reactants + reactants --[reagents above]--> products, conditions below the arrow.

    reaction: reaction SMILES "CC(=O)O.OCC>[H+]>CC(=O)OCC.O", or arrow syntax "acetic acid + ethanol -> ethyl acetate",
    or "A + B --[H2SO4, reflux]--> C" (bracket text goes above the arrow).  Alternatively give `reactants` and
    `products` lists (SMILES or names) plus `reagents_above` / `conditions_below` text lines (formula-style text such
    as "NaBH4" gets automatic subscripts).  arrow: full | equilibrium | retrosynthetic | double | dashed | bold |
    nogo | hollow | resonance.  product_labels e.g. ["85%"].  atom_map=true records atom mapping from map numbers.
    Opens a new ChemDraw document (or pastes into `document`), renders and returns a thumbnail + description.
    """
    try:
        step = _reaction_step(reaction, reactants, products, reagents_above, conditions_below, arrow, reactant_labels,
                              product_labels, atom_map)
        doc = layout.build_reaction(step, style=style, bond_length=bond_length, name=name)
        info, _xml = _materialize(doc, document, open_in_chemdraw, name)
        info["reaction_smiles"] = (cdx.read_cdxml(_xml).get("reaction_smiles") or [None])[0]
        png = None
        if render:
            if (info.get("opened") or info.get("inserted")) and _has_live():
                rinfo, png = _render_document(_live(), info.get("window_id") or info["document"], "png", None, name, return_image)
            else:
                rxn = chem.Reaction([c.mol for c in step.reactants], [], [c.mol for c in step.products])
                rinfo, png = _headless_png(rxn, name, None)
            info["render"] = rinfo
        return _with_image(dumps(info), png if return_image else None)
    except Exception as e:
        return _err(e)


@mcp.tool()
def draw_scheme(steps: List[Dict[str, Any]], title: Optional[str] = None, style: str = "acs", max_width: Optional[float] = None,
                document: Optional[str] = None, open_in_chemdraw: bool = True, render: bool = True, name: str = "scheme",
                bond_length: Optional[float] = None, return_image: bool = True):
    """Draw a multi-step synthesis scheme laid out left to right, wrapping into rows.

    steps: list of {"reactants": [...], "products": [...], "above": ["Br2"], "below": ["FeBr3, 0 °C"], "arrow": "full",
    "reactant_labels": [...], "product_labels": ["1", "72%"]} or {"reaction": "A + B --[reagents]--> C"}.
    Consecutive steps that start with the previous product are chained (the intermediate is drawn once).
    max_width: row width in points (default: page width); title: bold caption above the scheme.
    """
    try:
        if not steps:
            raise chem.ChemError("steps is empty")
        step_objs = []
        for s in steps:
            step_objs.append(_reaction_step(s.get("reaction"), s.get("reactants"), s.get("products"), s.get("above") or s.get("reagents_above"),
                                            s.get("below") or s.get("conditions_below"), s.get("arrow") or "full",
                                            s.get("reactant_labels"), s.get("product_labels"), bool(s.get("atom_map"))))
        doc = layout.build_scheme(step_objs, style=style, bond_length=bond_length, max_width=max_width, title=title, name=name)
        info, _xml = _materialize(doc, document, open_in_chemdraw, name)
        info["reaction_smiles"] = cdx.read_cdxml(_xml).get("reaction_smiles", [])
        png = None
        if render:
            if (info.get("opened") or info.get("inserted")) and _has_live():
                rinfo, png = _render_document(_live(), info.get("window_id") or info["document"], "png", None, name, return_image)
            else:
                mols = [c.mol for st in step_objs for c in st.reactants + st.products]
                png = _grid_png(mols, None)
                out = _out_path(name, "png")
                with open(out, "wb") as f:
                    f.write(png)
                rinfo = {"path": out, "rendered_by": "RDKit (components only; open the CDXML in ChemDraw for the scheme)"}
            info["render"] = rinfo
        return _with_image(dumps(info), png if return_image else None)
    except Exception as e:
        return _err(e)


def _resolve_stage_structures(stages: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Replace chemical names in mechanism stages with SMILES (same resolver chain as name_to_structure)."""
    out = []
    for st in stages:
        st = dict(st)
        for key in ("structure", "structures", "smiles"):
            val = st.get(key)
            if val is None:
                continue
            items = val if isinstance(val, list) else [val]
            fixed = []
            for it in items:
                if isinstance(it, str) and chem.detect_format(it) == "name":
                    smi = _resolve_name(it)
                    if not smi:
                        raise chem.ChemError("could not interpret %r as a structure or chemical name" % it)
                    fixed.append(smi)
                else:
                    fixed.append(it)
            st[key] = fixed if isinstance(val, list) else fixed[0]
        out.append(st)
    return out


@mcp.tool()
def draw_mechanism(stages: List[Dict[str, Any]], title: Optional[str] = None, style: str = "acs", max_width: Optional[float] = None,
                   document: Optional[str] = None, open_in_chemdraw: bool = True, render: bool = True, name: str = "mechanism",
                   bond_length: Optional[float] = None, return_image: bool = True):
    """Draw a reaction mechanism with curved electron-pushing arrows between stages (intermediates, transition states,
    resonance structures).

    stages: list of {"structures": ["[OH-]", "CCl"] (or "structure": one SMILES/name), "arrows": [...], "label": "caption",
    "arrow": kind of the arrow to the next stage (full | equilibrium | resonance ...), "above": [...], "below": [...]}.
    Electron arrows use 0-based RDKit atom indices (SMILES order) with "frag:atom" when a stage has several structures:
    "lp:0:0->0:0-1:0" = lone pair of atom 0 of structure 0 forms a bond to atom 0 of structure 1;
    "1:0-1:1->1:1" = bond 0-1 of structure 1 collapses onto its atom 1; "3=4->4-7" = pi bond to a new bond;
    "lp:5->5-8:1e" = single-electron (fishhook) arrow; "h:3->3-4" = a C-H bond as source.  Add ":1e" or "~>" for radicals.
    Returns the document, layout, the mechanism description (every arrow's tail/head) and a render.
    """
    try:
        from . import mechanism
        if not stages:
            raise chem.ChemError("stages is empty")
        doc = mechanism.build_mechanism(_resolve_stage_structures(stages), style=style, title=title, name=name,
                                        max_width=max_width, bond_length=bond_length)
        info, _xml = _materialize(doc, document, open_in_chemdraw, name)
        info["mechanism"] = mechanism.describe_mechanism(doc)
        png = None
        if render:
            if (info.get("opened") or info.get("inserted")) and _has_live():
                rinfo, png = _render_document(_live(), info.get("window_id") or info["document"], "png", None, name, return_image)
            else:
                mols = cdx.cdxml_mols(_xml)
                png = _grid_png(mols, None)
                out = _out_path(name, "png")
                with open(out, "wb") as f:
                    f.write(png)
                rinfo = {"path": out, "rendered_by": "RDKit (structures only; open the CDXML in ChemDraw for the arrows)"}
            info["render"] = rinfo
        return _with_image(dumps(info), png if return_image else None)
    except Exception as e:
        return _err(e)


@mcp.tool()
def named_reactions(query: Optional[str] = None, reaction_class: Optional[str] = None, name: Optional[str] = None) -> str:
    """Library of 100 named reaction templates (Suzuki, Grignard, Wittig, amide coupling, Boc protection, CuAAC,
    Diels-Alder, oxidations/reductions, ...).  query/reaction_class filter the list; name returns one template's full
    description (roles, typical reagents/conditions, SMARTS, example)."""
    try:
        from . import reactions
        if name:
            t = reactions.get_reaction(name)
            return dumps({"template": t.as_dict(), "description": reactions.describe(name)})
        rows = reactions.list_reactions(reaction_class, query)
        return dumps({"count": len(rows), "classes": reactions.classes(), "reactions": rows})
    except Exception as e:
        return _err(e)


@mcp.tool()
def apply_named_reaction(reaction: str, reactants: List[str], draw: bool = True, product_index: int = 0, style: str = "acs",
                         document: Optional[str] = None, open_in_chemdraw: bool = True, render: bool = True,
                         atom_map: bool = True, max_products: int = 8, return_image: bool = True):
    """Apply a named reaction template to reactants (SMILES or names, any order) to predict the product(s), and
    (draw=true) draw the atom-mapped reaction with the template's typical reagents above and conditions below the arrow.
    Example: apply_named_reaction("suzuki", ["bromobenzene", "phenylboronic acid"])."""
    try:
        from . import reactions
        resolved = []
        for r in reactants:
            if chem.detect_format(r) == "name":
                smi = _resolve_name(r)
                if not smi:
                    raise chem.ChemError("could not interpret %r as a structure or chemical name" % r)
                resolved.append(smi)
            else:
                resolved.append(r)
        outcomes = reactions.apply_reaction(reaction, resolved, max_products=max_products)
        out: Dict[str, Any] = {"template": outcomes[0]["template"], "name": outcomes[0]["name"], "outcomes": outcomes}
        if draw:
            step = reactions.to_step(reaction, resolved, product_index=product_index, atom_map=atom_map)
            doc = layout.build_reaction(step, style=style, name=_slug(outcomes[0]["template"]))
            info, _xml = _materialize(doc, document, open_in_chemdraw, _slug(outcomes[0]["template"]))
            png = None
            if render:
                if (info.get("opened") or info.get("inserted")) and _has_live():
                    rinfo, png = _render_document(_live(), info.get("window_id") or info["document"], "png", None, "reaction", return_image)
                else:
                    rxn = chem.Reaction([c.mol for c in step.reactants], [], [c.mol for c in step.products])
                    rinfo, png = _headless_png(rxn, "reaction", None)
                info["render"] = rinfo
            out["drawn"] = info
            return _with_image(dumps(out), png if return_image else None)
        return dumps(out)
    except Exception as e:
        return _err(e)


@mcp.tool()
def suggest_reactions(reactants: List[str]) -> str:
    """Which named reactions the given reactants (SMILES or names) can undergo together, with predicted products."""
    try:
        from . import reactions
        resolved = [(_resolve_name(r) or r) if chem.detect_format(r) == "name" else r for r in reactants]
        return dumps({"reactants": resolved, "suggestions": reactions.suggest_reactions(resolved)})
    except Exception as e:
        return _err(e)


@mcp.tool()
def insert(data: str, kind: str = "auto", document: Optional[str] = None) -> str:
    """Paste chemistry into a document through ChemDraw's own importers.  kind: auto | smiles | name | cdxml | mol |
    inchi | helm | sequence | peptide | dna | rna | text (plain caption).  The inserted objects stay selected."""
    try:
        b = _live()
        k = kind.lower()
        if k == "auto":
            f = chem.detect_format(data)
            if _looks_like_reaction(data, f):
                raise chem.ChemError("use draw_reaction for reactions")
            k = {"smiles": "smiles", "inchi": "inchi", "molblock": "mol", "cdxml": "cdxml", "name": "name", "sdf": "mol"}.get(f, "text")
        if k not in INSERT_KINDS:
            raise chem.ChemError("kind must be one of auto, %s" % ", ".join(INSERT_KINDS))
        res = b.insert(_doc(document), data, k)
        res.update({k2: v for k2, v in _doc_info(b, res.get("document")).items() if k2 == "window_id"})
        return dumps(res)
    except Exception as e:
        return _err(e)


def _paste_result(b: ChemDrawBackend, res: dict, count_key: str) -> dict:
    """add_caption / add_arrow result: the document's `<count_key>_after` taken from documents() (both backends count
    captions and arrows there; the insert result itself never reports arrows), plus the window id."""
    info = _doc_info(b, res.get("document"))
    n = info.get(count_key)
    if n is None:  # documents() failed: fall back to whatever the insert reported (captions only)
        n = res.get("after", {}).get(count_key, res.get(count_key))
    return {"ok": True, "document": res.get("document"), "window_id": info.get("window_id"), count_key + "_after": n}


@mcp.tool()
def add_caption(text: str, document: Optional[str] = None, size: Optional[float] = None, bold: bool = False,
                italic: bool = False, formula: bool = False, color: Optional[str] = None) -> str:
    """Add a text caption to a document (ChemDraw places it; it stays selected so the user can move it).
    formula=true enables automatic sub/superscripts for formulas like H2SO4."""
    try:
        b = _live()
        doc = cdx.CDXMLDocument(DEFAULT_STYLE)
        faces = [f for f, on in (("bold", bold), ("italic", italic)) if on] or "plain"
        doc.add_text(text, 100, 100, size=size, face=faces, formula=formula, color=color)
        res = b.insert(_doc(document), doc.to_xml(), "cdxml")
        return dumps(_paste_result(b, res, "captions"))
    except Exception as e:
        return _err(e)


@mcp.tool()
def add_arrow(kind: str = "full", length: Optional[float] = None, document: Optional[str] = None,
              text_above: Optional[str] = None, text_below: Optional[str] = None) -> str:
    """Add a reaction arrow (full, equilibrium, retrosynthetic, double, dashed, bold, nogo, hollow) with optional
    text above/below, pasted into a document (ChemDraw places it; it stays selected)."""
    try:
        b = _live()
        doc = cdx.CDXMLDocument(DEFAULT_STYLE)
        L = doc.bond_length
        length = length or 3 * L
        ar = doc.add_arrow(100, 100, 100 + length, 100, kind)
        if text_above:
            doc.add_text(text_above, ar.box.cx, 100 - 0.45 * L - 2, formula=True)
        if text_below:
            doc.add_text(text_below, ar.box.cx, 100 + 0.45 * L, anchor="top")
        res = b.insert(_doc(document), doc.to_xml(), "cdxml")
        return dumps(_paste_result(b, res, "arrows"))
    except Exception as e:
        return _err(e)


# ---------------------------------------------------------------------------
# rendering & looking
# ---------------------------------------------------------------------------

@mcp.tool()
def render(document: Optional[str] = None, format: str = "png", path: Optional[str] = None, width: Optional[int] = None,
           return_image: bool = True):
    """Export a document as an image (png at 300 dpi, svg, pdf, eps, tiff, gif, jpeg, bmp) and describe it numerically:
    size, coverage, bounding box, margins, dominant colours.  Returns a thumbnail you can look at (png formats)."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        info, png = _render_document(b, handle, format, path, "render", return_image, width)
        info["document"] = _doc_info(b, handle if handle is not None else b.front_document()).get("name")
        return _with_image(dumps(info), png)
    except Exception as e:
        return _err(e)


@mcp.tool()
def quick_render(input: str, style: str = "acs", format: str = "png", path: Optional[str] = None, label: Optional[str] = None,
                 return_image: bool = True):
    """Render a structure or reaction (SMILES, name, InChI, molblock, reaction SMILES / "A + B -> C") to an image file
    without leaving a document open in ChemDraw (uses a temporary document; falls back to RDKit without ChemDraw)."""
    try:
        if _looks_like_reaction(input):
            step = _reaction_step(input, None, None, None, None, "full")
            doc = layout.build_reaction(step, style=style, name="quick")
            headless_obj: Any = chem.Reaction([c.mol for c in step.reactants], [], [c.mol for c in step.products])
        else:
            mol = _mol(input)
            doc = layout.build_structures([mol], labels=[label] if label else None, style=style, name="quick")
            headless_obj = mol
        if _has_live():
            b = _live()
            xml_path = _out_path("quick", "cdxml")
            with open(xml_path, "w", encoding="utf-8") as fh:
                fh.write(doc.to_xml())
            name = b.open_document(xml_path)
            wid = _doc_info(b, name).get("window_id") or name
            try:
                info, png = _render_document(b, wid, format, path, "quick", return_image)
            finally:
                try:
                    b.close_document(wid)
                except ChemDrawError:
                    pass
            info["cdxml_path"] = xml_path
        else:
            fmt_h = normalize_format(format)
            if fmt_h not in ("png", "svg"):
                raise chem.ChemError("without ChemDraw only png and svg can be rendered (asked for %s)" % fmt_h)
            if fmt_h == "svg":
                svg = chem.depict_svg(headless_obj)
                out = _out_path("quick", "svg", path)
                with open(out, "w") as fh:
                    fh.write(svg)
                return dumps({"path": out, "format": "svg", "rendered_by": "RDKit"})
            info, png = _headless_png(headless_obj, "quick", path)
        return _with_image(dumps(info), png if return_image else None)
    except Exception as e:
        return _err(e)


@mcp.tool()
def screenshot(document: Optional[str] = None, return_image: bool = True):
    """Screenshot of a document window exactly as the user sees it (needs Screen Recording permission on macOS)."""
    try:
        b = _live(autolaunch=False)
        png = b.screenshot(_doc(document))
        out = _out_path("screenshot", "png")
        with open(out, "wb") as f:
            f.write(png)
        return _with_image(dumps({"path": out, "analysis": _analysis(png)}), png if return_image else None)
    except Exception as e:
        return _err(e)


@mcp.tool()
def analyze_image(path: str) -> str:
    """Numerical description of any PNG/JPEG on disk (coverage, bounding box, margins, colours, contrast, ASCII map)."""
    try:
        p = _in_path(path)
        a = ia.analyze_render(p)
        a["path"] = p
        a["summary"] = ia.summarize_render(a)
        return dumps(a)
    except Exception as e:
        return _err(e)


@mcp.tool()
def compare_renders(path_a: str, path_b: str) -> str:
    """Diff two images (e.g. before/after a change): changed-pixel fraction, changed regions, size differences."""
    try:
        return dumps(ia.diff_images(_in_path(path_a), _in_path(path_b)))
    except Exception as e:
        return _err(e)


# ---------------------------------------------------------------------------
# chemistry
# ---------------------------------------------------------------------------

@mcp.tool()
def name_to_structure(name: str, document: Optional[str] = None, draw: bool = False) -> str:
    """Chemical name -> structure: offline dictionary (~3,100 names) -> OPSIN (if installed) -> ChemDraw's Name=Struct ->
    PubChem (CHEMDRAWMCP_ONLINE=1).
    Returns SMILES, InChIKey, formula, MW...  draw=true also draws it into a (new or given) document."""
    try:
        smi = _resolve_name(name)
        if not smi:
            raise chem.ChemError("could not convert %r to a structure" % name)
        mol = chem.mol_from_input(smi, "smiles")
        out = {"name": name, "smiles": smi, "properties": chem.properties(mol, extended=False),
               "source": _NAME_SOURCE.get(name.strip().lower(), "unknown")}
        if draw:
            doc = layout.build_structures([mol], labels=[name], style=DEFAULT_STYLE, name=_slug(name))
            info, _ = _materialize(doc, document, True, _slug(name))
            out["drawn"] = info
        return dumps(out)
    except Exception as e:
        return _err(e)


@mcp.tool()
def structure_to_name(structure: Optional[str] = None, document: Optional[str] = None) -> str:
    """IUPAC-style name from ChemDraw's Struct=Name.  Give a structure (SMILES/InChI/molblock) or use the current
    selection of a document (select the structure first; whole document if nothing is selected)."""
    try:
        b = _live()
        if structure:
            mol = _mol(structure)
            with _temp_document(b) as h:
                b.insert(h, chem.mol_to_molblock(mol), "mol")
                b.select_all(h)
                before = set(cdx.document_text(b.document_cdxml(h)))
                res = b.run_command("convertStructureToName", h)
                after = cdx.document_text(b.document_cdxml(h))
                names = [t for t in after if t not in before]
                return dumps({"ok": bool(names), "name": names[-1] if names else None, "smiles": chem.properties(mol, False)["canonical_smiles"],
                              "all_new_text": names, "alerts": res.get("alerts")})
        handle = _doc(document)
        info = b.selection_info(handle)
        if info.get("empty"):
            b.select_all(handle)
        before = set(cdx.document_text(b.document_cdxml(handle)))
        res = b.run_command("convertStructureToName", handle)
        after = cdx.document_text(b.document_cdxml(handle))
        names = [t for t in after if t not in before]
        return dumps({"ok": bool(names), "name": names[-1] if names else None, "smiles": info.get("smiles"), "all_new_text": names,
                      "alerts": res.get("alerts"), "note": "the name was added as a caption below the structure"})
    except Exception as e:
        return _err(e)


@mcp.tool()
def properties(structure: Optional[str] = None, document: Optional[str] = None, extended: bool = True,
               use_chemdraw: bool = False) -> str:
    """Chemical properties.  structure: SMILES/InChI/molblock/name -> RDKit properties (formula, MW, exact mass,
    elemental analysis, InChI/InChIKey, logP, TPSA, HBD/HBA, rotatable bonds, stereocentres, Lipinski, QED).
    document: ChemDraw's own values for the current selection (whole document if nothing is selected) plus RDKit.
    use_chemdraw=true also computes ChemDraw's values for `structure` (slower: uses a temporary document)."""
    try:
        out: Dict[str, Any] = {}
        if structure:
            mol = _mol(structure)
            out["rdkit"] = chem.properties(mol, extended)
            if use_chemdraw and _has_live():
                b = _live()
                with _temp_document(b) as h:
                    b.insert(h, chem.mol_to_molblock(mol), "mol")
                    out["chemdraw"] = b.selection_info(h)
        elif _has_live():
            b = _live(autolaunch=False)
            handle = _doc(document)
            info = b.selection_info(handle)
            if info.get("empty"):
                b.select_all(handle)
                info = b.selection_info(handle)
            out["chemdraw"] = info
            if info.get("smiles"):
                try:
                    out["rdkit"] = chem.properties(chem.mol_from_input(info["smiles"], "smiles"), extended)
                except chem.ChemError as e:
                    out["rdkit_error"] = str(e)
        else:
            raise chem.ChemError("give `structure` (ChemDraw is not available for document selections)")
        return dumps(out)
    except Exception as e:
        return _err(e)


@mcp.tool()
def analyze_reaction(reaction: str) -> str:
    """Balance check and bookkeeping for a reaction (reaction SMILES or "A + B -> C" with SMILES/names): per-component
    formula/MW, element and charge balance, atom economy."""
    try:
        rxn = chem.parse_reaction(_resolve_names_in_reaction(reaction))
        return dumps(chem.reaction_properties(rxn))
    except Exception as e:
        return _err(e)


@mcp.tool()
def convert(data: str, to: str, from_format: Optional[str] = None, path: Optional[str] = None, use_chemdraw: bool = True) -> str:
    """Convert chemistry between formats.  Input: SMILES, InChI, molblock, CDXML, a file path, or a chemical name.
    to: smiles | inchi | inchikey | mol | molv3000 | sdf | sdfv3000 | cdxml | cdx | rxn | rxnv3000 | cml | helm | sln | name |
    png | svg | pdf | eps | tiff.  Text results are returned inline; file formats are written to `path` (or the renders dir).
    ChemDraw does the conversion when available (use_chemdraw=false forces RDKit)."""
    try:
        to_l = to.lower().strip().lstrip(".")
        text_formats = {"smiles", "inchi", "inchikey", "mol", "molv3000", "cdxml", "helm", "sln", "name"}
        cdx_text = _cdx_text(data) if not from_format or from_format.lower() == "cdx" else None
        if cdx_text is not None:
            data, from_format = cdx_text, "cdxml"
        norm = getattr(chem, "_norm_fmt", None)
        f = (norm(from_format) if (from_format and norm) else from_format) or chem.detect_format(data)
        if to_l == "name":
            return structure_to_name(data)
        if use_chemdraw and _has_live() and f != "reaction_smiles":
            b = _live()
            kind = {"smiles": "smiles", "inchi": "inchi", "molblock": "mol", "cdxml": "cdxml", "name": "name", "sdf": "mol"}.get(f)
            payload = data
            if f in ("molblock", "sdf", "cdxml") and os.path.exists(os.path.expanduser(data)):
                with open(os.path.expanduser(data), encoding="utf-8", errors="replace") as fh:
                    payload = fh.read()
            if kind is None:
                raise chem.ChemError("cannot feed %s input to ChemDraw" % f)
            with _temp_document(b) as h:
                r = b.insert(h, payload, kind)
                if r.get("atoms_added", 0) <= 0:
                    raise chem.ChemError("ChemDraw did not accept the input as %s" % kind)
                if to_l in text_formats:
                    text = b.selection_as(h, to_l)
                    return dumps({"ok": True, "to": to_l, "result": _trim(text), "converted_by": "ChemDraw"})
                fmt = normalize_format(to_l)
                out = b.save_document(h, _out_path("convert", EXPORT_FORMATS[fmt]["ext"].lstrip("."), path), fmt)
                return dumps({"ok": True, "to": fmt, "path": out, "bytes": os.path.getsize(out), "converted_by": "ChemDraw"})
        # RDKit path
        from rdkit import Chem
        mol = _mol(data, None if from_format is None else from_format)
        if to_l == "smiles":
            res = Chem.MolToSmiles(mol)
        elif to_l == "inchi":
            res = Chem.MolToInchi(mol)
        elif to_l == "inchikey":
            res = Chem.MolToInchiKey(mol)
        elif to_l == "mol":
            res = chem.mol_to_molblock(mol)
        elif to_l == "molv3000":
            res = chem.mol_to_molblock(mol, v3000=True)
        elif to_l == "sdf":
            res = chem.mol_to_sdf([mol])
        elif to_l == "cdxml":
            doc = layout.build_structures([mol], style=DEFAULT_STYLE, name="convert")
            res = doc.to_xml()
        elif to_l in ("png", "svg"):
            if to_l == "svg":
                svg = chem.depict_svg(mol)
                out = _out_path("convert", "svg", path)
                with open(out, "w") as fh:
                    fh.write(svg)
            else:
                out = _out_path("convert", "png", path)
                with open(out, "wb") as fh:
                    fh.write(chem.depict_png(mol))
            return dumps({"ok": True, "to": to_l, "path": out, "converted_by": "RDKit"})
        else:
            raise chem.ChemError("RDKit cannot write %s; ChemDraw is needed for that format" % to_l)
        out = None
        if path and to_l in text_formats:
            out = _out_path("convert", to_l, path)   # expands `~`, creates the directory, absolute
            with open(out, "w", encoding="utf-8") as fh:
                fh.write(res)
        return dumps({"ok": True, "to": to_l, "result": _trim(res), "converted_by": "RDKit", "path": out})
    except Exception as e:
        return _err(e)


@mcp.tool()
def clean_structure(document: Optional[str] = None, reaction: bool = False, select_all: bool = True) -> str:
    """Run ChemDraw's Clean Up Structure (or Clean Up Reaction) on the selection (whole document by default)."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        if select_all:
            b.select_all(handle)
        res = b.run_command("cleanReaction" if reaction else "cleanStructure", handle)
        return dumps(res)
    except Exception as e:
        return _err(e)


@mcp.tool()
def check_structure(document: Optional[str] = None) -> str:
    """Chemical warnings for a document: ChemDraw's per-object warnings (valence, unusual bonds...) from its CDXML
    plus the result of the Check Structure command."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        xml = b.document_cdxml(handle)
        warnings, caption_notes = [], []
        for m in re.finditer(r'<(n|b|t|fragment)\b[^>]*\bWarning="([^"]*)"', xml):
            text = m.group(2).replace("&apos;", "'").replace("&quot;", '"')
            # ChemDraw also flags every ordinary caption it cannot read as chemistry; keep those apart
            if m.group(1) == "t" and ("interpret" in text.lower() or "Interpretation" in text):
                caption_notes.append(text)
            else:
                warnings.append({"object": m.group(1), "warning": text})
        b.select_all(handle)
        res = b.run_command("checkStructure", handle)
        return dumps({"warnings": warnings, "captions_not_chemical": len(caption_notes),
                      "check_structure": res.get("alerts") or res.get("error") or "no dialog shown", "ok": not warnings})
    except Exception as e:
        return _err(e)


@mcp.tool()
def predict_nmr(structure: Optional[str] = None, document: Optional[str] = None, nucleus: str = "1H", keep_open: bool = True,
                return_image: bool = True):
    """Predict a 1H or 13C NMR spectrum with ChemDraw (structure given as SMILES/name, or the selection of a document).
    ChemDraw opens a spectrum document; it is rendered to PNG and its text (shift assignments) is returned."""
    try:
        b = _live()
        nuc = "predict13C_NMR_Shifts" if "13" in nucleus else "predict1H_NMR_Shifts"
        before = {d.name for d in b.documents()}

        def run(handle):
            b.select_all(handle)
            return b.run_command(nuc, handle)

        if structure:
            mol = _mol(structure)
            with _temp_document(b) as h:
                b.insert(h, chem.mol_to_molblock(mol), "mol")
                res = run(h)
        else:
            res = run(_doc(document))
        new = [d for d in b.documents() if d.name not in before]
        if not new:
            return dumps({"ok": False, "error": "no spectrum document appeared", "command": res})
        spec = new[0]
        info, png = _render_document(b, spec.window_id or spec.name, "png", None, "nmr_%s" % nucleus, return_image)
        texts = cdx.document_text(b.document_cdxml(spec.window_id or spec.name))
        if not keep_open:
            b.close_document(spec.window_id or spec.name)
        return _with_image(dumps({"ok": True, "spectrum_document": spec.to_dict(), "render": info, "text": texts[:200]}), png)
    except Exception as e:
        return _err(e)


# ---------------------------------------------------------------------------
# selection, commands, menus
# ---------------------------------------------------------------------------

@mcp.tool()
def selection(document: Optional[str] = None, as_format: str = "info") -> str:
    """The current selection of a document.  as_format: info (SMILES, formula, MW, exact mass, elemental analysis,
    caption text, bounds, counts) | smiles | inchi | inchikey | mol | molv3000 | cdxml | helm | sln."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        if as_format == "info":
            return dumps(b.selection_info(handle))
        return dumps({"format": as_format, "result": _trim(b.selection_as(handle, as_format))})
    except Exception as e:
        return _err(e)


@mcp.tool()
def select(mode: str = "all", document: Optional[str] = None) -> str:
    """Change the selection: all | none | next (next structure) | previous | invert."""
    try:
        b = _live(autolaunch=False)
        handle = _doc(document)
        m = mode.lower()
        if m == "all":
            b.select_all(handle)
        elif m == "none":
            b.deselect(handle)
        else:
            cmdname = {"next": "selectNextStructure", "previous": "selectPreviousStructure", "invert": "invertSelection"}.get(m)
            if not cmdname:
                raise chem.ChemError("mode must be all, none, next, previous or invert")
            b.run_command(cmdname, handle)
        return dumps(b.selection_info(handle))
    except Exception as e:
        return _err(e)


@mcp.tool()
def commands(pattern: Optional[str] = None, enabled_only: bool = False, limit: int = 150) -> str:
    """List ChemDraw's named commands (regex filter), with enabled state and descriptions.  Examples: pattern="clean",
    "arrow", "copyAs", "paste", "NMR", "reaction", "flip|rotate", "toggle".  Run them with run_command."""
    try:
        b = _live(autolaunch=False)
        rows = b.list_commands(pattern, enabled_only, with_descriptions=True)
        out = [r.to_dict() for r in rows[:limit]]
        return dumps({"count": len(rows), "shown": len(out), "commands": out})
    except Exception as e:
        return _err(e)


@mcp.tool()
def run_command(command: str, document: Optional[str] = None) -> str:
    """Execute a ChemDraw named command on a document (e.g. cleanStructure, cleanReaction, mapReactionAtoms,
    autoNumberReaction, convertStructureToName, toggleShowingAtomNumber, flipObjectsHorizontal, centerObjectOnPage,
    zoomToFit, analyzeStoichiometry, expandLabel, toggleKekulization, copyAsSMILES, selectAll, undo).  Reports object
    counts before/after and any alert ChemDraw showed."""
    try:
        return dumps(_live(autolaunch=False).run_command(command, _doc(document)))
    except Exception as e:
        return _err(e)


@mcp.tool()
def menus(menu: Optional[str] = None) -> str:
    """ChemDraw's menus and their items (optionally one menu).  Invoke items with do_menu."""
    try:
        m = _live(autolaunch=False).menus()
        if menu:
            return dumps({menu: m.get(menu) or [k for k in m if k.lower().startswith(menu.lower())]})
        return dumps(m)
    except Exception as e:
        return _err(e)


@mcp.tool()
def do_menu(menu: str, item: str, document: Optional[str] = None) -> str:
    """Invoke a menu item by titles, e.g. do_menu("Structure", "Clean Up Structure"), do_menu("View", "Show Chemical Warnings")."""
    try:
        return dumps(_live(autolaunch=False).do_menu(menu, item, _doc(document)))
    except Exception as e:
        return _err(e)


@mcp.tool()
def dismiss_dialogs() -> str:
    """Close any alert/sheet ChemDraw is showing and return its text (tools do this automatically on errors)."""
    try:
        b = _live(autolaunch=False)
        return dumps({"dismissed": b.dismiss_dialogs()})
    except Exception as e:
        return _err(e)


@mcp.tool()
def styles() -> str:
    """Available document styles (from ChemDraw's stationery): bond length, fonts, sizes, line widths, page size."""
    out = {}
    for k in styles_mod.style_names():
        s = styles_mod.get_style(k)
        out[k] = {kk: s[kk] for kk in ("name", "bond_length", "line_width", "bold_width", "label_font", "label_size", "caption_font",
                                      "caption_size", "bond_spacing", "page") if kk in s}
    return dumps({"default": DEFAULT_STYLE, "styles": out})


def main():
    mcp.run()


if __name__ == "__main__":
    main()
