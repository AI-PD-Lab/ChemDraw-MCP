"""RDKit-based chemistry utilities for chemdrawmcp.

Parsing (SMILES, InChI, molblock, SDF, CDXML, reactions, chemical names),
preparation for depiction (2D coordinates, Kekulé bonds, wedge bonds),
properties, reactions, PNG/SVG rendering, similarity and name resolution.

RDKit is imported lazily so that this module (and the package) imports even
when RDKit is missing; every function that needs it raises ``ChemError``
with an installation hint instead.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import types
import urllib.parse
import urllib.request
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union

if TYPE_CHECKING:  # pragma: no cover - typing only
    from rdkit import Chem

__all__ = [
    "ChemError",
    "Reaction",
    "rdkit_available",
    "detect_format",
    "mol_from_input",
    "mols_from_input",
    "prepare_for_depiction",
    "properties",
    "elemental_analysis",
    "parse_reaction",
    "reaction_smiles",
    "reaction_properties",
    "depict_png",
    "depict_svg",
    "mol_to_molblock",
    "mol_to_sdf",
    "substructure_matches",
    "tanimoto",
    "resolve_name",
    "describe",
]

FORMATS = (
    "smiles", "inchi", "molblock", "sdf", "cdxml", "rxn", "reaction_smiles", "name",
)


class ChemError(ValueError):
    """Raised for invalid chemical input or a failed chemistry operation."""


# --------------------------------------------------------------------------- #
# RDKit access
# --------------------------------------------------------------------------- #
_RD: Any = None


def rdkit_available() -> bool:
    """Return True when RDKit can be imported."""
    try:
        from rdkit import Chem  # noqa: F401
    except Exception:  # pragma: no cover - depends on environment
        return False
    return True


def _rd() -> Any:
    """Import RDKit lazily and cache the modules in a namespace."""
    global _RD
    if _RD is None:
        try:
            from rdkit import Chem, DataStructs, rdBase
            from rdkit.Chem import (
                AllChem, Crippen, Descriptors, Draw, Lipinski, rdChemReactions,
                rdDepictor, rdMolDescriptors,
            )
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise ChemError(
                "RDKit is required for this operation but is not installed "
                "(pip install rdkit)."
            ) from exc
        _RD = types.SimpleNamespace(
            Chem=Chem, DataStructs=DataStructs, rdBase=rdBase, AllChem=AllChem,
            Crippen=Crippen, Descriptors=Descriptors, Draw=Draw, Lipinski=Lipinski,
            rdChemReactions=rdChemReactions, rdDepictor=rdDepictor,
            rdMolDescriptors=rdMolDescriptors,
        )
    return _RD


@contextmanager
def _quiet() -> Iterator[None]:
    """Silence RDKit's C++ log output inside the block."""
    rd = _rd()
    block = getattr(rd.rdBase, "BlockLogs", None)
    if block is not None:
        with block():
            yield
        return
    from rdkit import RDLogger  # pragma: no cover - old RDKit
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")


def _exc_msg(exc: BaseException) -> str:
    """First meaningful line of an exception message."""
    text = str(exc).strip()
    for line in text.splitlines():
        line = line.strip()
        if line and not line.startswith("*"):
            return line
    return text or exc.__class__.__name__


# --------------------------------------------------------------------------- #
# Input handling / format detection
# --------------------------------------------------------------------------- #
_PATH_FORMATS: Dict[str, str] = {
    ".mol": "molblock", ".mdl": "molblock", ".sdf": "sdf", ".sd": "sdf",
    ".cdxml": "cdxml", ".rxn": "rxn", ".smi": "smiles", ".smiles": "smiles",
    ".inchi": "inchi",
}
_FORMAT_ALIASES: Dict[str, str] = {
    "smi": "smiles", "smarts": "smiles", "sma": "smiles", "cxsmiles": "smiles",
    "mol": "molblock", "molfile": "molblock", "mdl": "molblock", "ctab": "molblock",
    "sd": "sdf", "mdl_sdf": "sdf",
    "rxnsmiles": "reaction_smiles", "reaction": "reaction_smiles", "rsmi": "reaction_smiles",
    "rxn_smiles": "reaction_smiles",
    "iupac": "name", "names": "name", "text": "name",
}


def _as_path(data: str) -> Optional[Path]:
    """Return a Path when ``data`` names an existing file, else None."""
    s = data.strip()
    if not s or "\n" in s or len(s) > 1024:
        return None
    try:
        p = Path(os.path.expanduser(s))
        if p.is_file():
            return p
    except (OSError, ValueError):
        pass
    return None


def _read_text(path: Path) -> str:
    """File contents as text; binary ChemDraw .cdx files are converted to CDXML on the fly."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        head = b""
    if head == b"VjCD0100":
        from . import cdx as _cdx
        with open(path, "rb") as fh:
            return _cdx.cdx_to_cdxml(fh.read())
    return _read_text_plain(path)


def _read_text_plain(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ChemError(f"Cannot read file {path}: {exc}") from exc


def _norm_fmt(fmt: str) -> str:
    f = fmt.strip().lower().lstrip(".")
    f = _FORMAT_ALIASES.get(f, f)
    if f == "cdx":
        raise ChemError("Binary CDX is not supported; export CDXML from ChemDraw instead.")
    if f not in FORMATS:
        raise ChemError(f"Unknown format {fmt!r}; expected one of {', '.join(FORMATS)}.")
    return f


def _looks_like_molblock(s: str) -> bool:
    if "M  END" in s or "M  V30 BEGIN CTAB" in s:
        return True
    lines = s.splitlines()
    if len(lines) >= 4 and ("V2000" in lines[3] or "V3000" in lines[3]):
        return True
    return bool(re.search(r"^\s*\d{1,3}\s+\d{1,3}\s.*V2000\s*$", s, re.M))


def _looks_like_cdxml(s: str) -> bool:
    low = s.lstrip().lower()
    return "<cdxml" in low or "<fragment" in low or (low.startswith("<?xml") and "cdx" in low)


_RXN_SMILES_RE = re.compile(r"[^\s>]*>[^\s>]*>[^\s>]*")
# SMILES dative-bond tokens ("[NH3]->[Pt]", "[Fe]<-N") sit directly between atom
# characters; a reaction arrow written in text ("A -> B") is surrounded by spaces.
_DATIVE_RE = re.compile(r"(?<=\S)(?:->|<-)(?=\S)")
# ">" that is not the tail of a dative "->" token
_RXN_SEP_RE = re.compile(r"(?<!-)>")


def _strip_dative(s: str) -> str:
    """Remove dative-bond tokens so that ``>`` counting sees only reaction arrows."""
    return _DATIVE_RE.sub("", s)


def _has_reaction_arrow(s: str) -> bool:
    """True when ``s`` holds a ``>`` (or unicode arrow) outside dative-bond tokens."""
    t = _strip_dative(s)
    return ">" in t or any(a in t for a in _UNICODE_ARROWS)


def _split_reaction_smiles(core: str) -> List[str]:
    """Split ``reactants>agents>products`` on ``>`` while keeping ``->`` dative bonds."""
    return _RXN_SEP_RE.split(core)


def _is_reaction_smiles(s: str) -> bool:
    core = s.split(" |", 1)[0].strip()
    probe = _strip_dative(core)
    return bool(_RXN_SMILES_RE.fullmatch(probe)) and probe != ">>"


def _smiles_ok(token: str) -> bool:
    """True when ``token`` parses as SMILES (quietly)."""
    Chem = _rd().Chem
    with _quiet():
        return Chem.MolFromSmiles(token) is not None


def _smarts_ok(token: str) -> bool:
    Chem = _rd().Chem
    with _quiet():
        return Chem.MolFromSmarts(token) is not None


def _detect_text(text: str) -> str:
    s = text.strip()
    if not s:
        raise ChemError("Empty input.")
    if _looks_like_cdxml(s):
        return "cdxml"
    if s.startswith("$RXN"):
        return "rxn"
    if _looks_like_molblock(s):
        return "sdf" if "$$$$" in s else "molblock"
    if s.startswith("InChI="):
        return "inchi"
    first = s.splitlines()[0].strip()
    if "\n" not in s and _has_reaction_arrow(s):
        return "reaction_smiles"
    token = first.split()[0]
    if _smiles_ok(token):
        return "smiles"
    if " " not in first and _smarts_ok(token):
        return "smiles"
    return "name"


def detect_format(data: str) -> str:
    """Guess the format of ``data``.

    Returns one of ``smiles``, ``inchi``, ``molblock``, ``sdf``, ``cdxml``,
    ``rxn``, ``reaction_smiles`` or ``name``.  ``data`` may also be a path to
    a ``.mol``/``.sdf``/``.cdxml``/``.rxn``/``.smi`` file.  SMARTS strings are
    reported as ``smiles``; anything unparseable is assumed to be a chemical
    name.
    """
    if not isinstance(data, str) or not data.strip():
        raise ChemError(
            "Empty input: expected a SMILES, InChI, molblock, SDF, CDXML, "
            "reaction or chemical name."
        )
    path = _as_path(data)
    if path is not None:
        fmt = _PATH_FORMATS.get(path.suffix.lower())
        return fmt or _detect_text(_read_text(path))
    return _detect_text(data)


def _load(data: Any, fmt: Optional[str]) -> Tuple[str, str]:
    """Resolve file paths and the format; return ``(text, fmt)``."""
    if not isinstance(data, str):
        raise ChemError(f"Input must be a string, got {type(data).__name__}.")
    if not data.strip():
        raise ChemError("Empty input.")
    path = _as_path(data)
    if path is not None:
        text = _read_text(path)
        if fmt is None:
            f = _PATH_FORMATS.get(path.suffix.lower()) or _detect_text(text)
        else:
            f = _norm_fmt(fmt)
        return text, f
    return data, (_norm_fmt(fmt) if fmt else _detect_text(data))


# --------------------------------------------------------------------------- #
# Record parsers
# --------------------------------------------------------------------------- #
_ORGANIC_BRACKET_RE = re.compile(r"\[(Cl|Br|B|C|N|O|P|S|F|I|c|n|o|p|s)\]")


def _parse_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES (falling back to SMARTS); raise ChemError on failure."""
    Chem = _rd().Chem
    smiles = smiles.strip()
    with _quiet():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            query = Chem.MolFromSmarts(smiles)
            if query is not None:
                # SMARTS -> SMILES; unbracket plain organic-subset atoms so they
                # regain implicit hydrogens ([#6][#8] -> CO rather than [C][O]).
                plain = _ORGANIC_BRACKET_RE.sub(r"\1", Chem.MolToSmiles(query))
                mol = Chem.MolFromSmiles(plain)
    if mol is None:
        raise ChemError(
            f"Could not parse SMILES {smiles!r}. Check the syntax, or pass fmt='name' "
            "if this is a chemical name."
        )
    return mol


_CX_RE = re.compile(r"^\S+\s+\|[^|]*\|")


def _smiles_records(text: str, first_only: bool) -> List[Chem.Mol]:
    """Parse newline/whitespace separated SMILES, optionally 'SMILES name'."""
    Chem = _rd().Chem
    mols: List[Chem.Mol] = []
    errors: List[str] = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        cx = _CX_RE.match(line)
        if cx and _smiles_ok(cx.group(0).strip()):
            # CXSMILES ("CC |c:1|"): RDKit parses the extension block itself; the rest of the line is the name
            smiles_list, name = [cx.group(0).strip()], line[cx.end():].strip()
        else:
            tokens = line.split()
            extra = tokens[1:]
            if extra and all(_smiles_ok(t) for t in extra):
                smiles_list, name = tokens, ""
            else:
                smiles_list, name = tokens[:1], " ".join(extra)
        for smi in smiles_list:
            try:
                mol = _parse_smiles(smi)
            except ChemError as exc:
                errors.append(f"line {lineno}: {_exc_msg(exc)}")
                continue
            if name:
                mol.SetProp("_Name", name)
            mols.append(mol)
            if first_only:
                return mols
        if first_only and errors:
            break
    if errors and (first_only or not mols):
        raise ChemError("Invalid SMILES input: " + "; ".join(errors))
    if errors:
        raise ChemError(
            f"{len(errors)} SMILES record(s) could not be parsed: " + "; ".join(errors)
        )
    return mols


def _inchi_records(text: str, first_only: bool) -> List[Chem.Mol]:
    Chem = _rd().Chem
    mols: List[Chem.Mol] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        token = line.split()[0]
        with _quiet():
            mol = Chem.MolFromInchi(token)
        if mol is None:
            raise ChemError(f"Could not parse InChI {token!r}.")
        rest = line[len(token):].strip()
        if rest:
            mol.SetProp("_Name", rest)
        mols.append(mol)
        if first_only:
            break
    return mols


def _molblock_records(text: str, first_only: bool) -> List[Chem.Mol]:
    """Parse molblock/SDF text (all records, or just the first)."""
    Chem = _rd().Chem
    if "$$$$" not in text:
        try:
            with _quiet():
                mol = Chem.MolFromMolBlock(text, removeHs=True)
                if mol is None:
                    mol = Chem.MolFromMolBlock(text, removeHs=True, strictParsing=False)
        except Exception as exc:
            raise ChemError(f"Could not parse molblock: {_exc_msg(exc)}") from exc
        if mol is None:
            raise ChemError("Could not parse molblock (check the counts line and 'M  END').")
        return [mol]
    mols: List[Chem.Mol] = []
    failed = 0
    with _quiet():
        supplier = Chem.SDMolSupplier()
        supplier.SetData(text, removeHs=True)
        for mol in supplier:
            if mol is None:
                failed += 1
                continue
            mols.append(mol)
            if first_only:
                break
    if not mols:
        raise ChemError(
            f"No valid molecule records in SDF input ({failed} record(s) failed to parse)."
        )
    return mols


def _cdxml_records(text: str) -> List[Chem.Mol]:
    Chem = _rd().Chem
    try:
        with _quiet():
            mols = list(Chem.MolsFromCDXML(text))
    except Exception as exc:
        raise ChemError(f"Could not parse CDXML: {_exc_msg(exc)}") from exc
    if not mols:
        raise ChemError("No structures found in the CDXML document.")
    return mols


def _name_records(text: str, first_only: bool) -> List[Chem.Mol]:
    names = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if first_only:
        names = names[:1]
    mols: List[Chem.Mol] = []
    failed: List[str] = []
    for name in names:
        smiles = resolve_name(name)
        if smiles is None:
            failed.append(name)
            continue
        mol = _parse_smiles(smiles)
        mol.SetProp("_Name", name)
        mols.append(mol)
    if failed:
        raise ChemError(
            "Could not resolve chemical name(s): " + ", ".join(repr(n) for n in failed)
            + ". Provide a SMILES/InChI/molblock instead, install OPSIN (py2opsin), "
            "or set CHEMDRAWMCP_ONLINE=1 to allow PubChem lookups."
        )
    return mols


def _records(text: str, fmt: str, first_only: bool, prepare: bool) -> List[Chem.Mol]:
    if fmt == "smiles":
        mols = _smiles_records(text, first_only)
    elif fmt == "inchi":
        mols = _inchi_records(text, first_only)
    elif fmt in ("molblock", "sdf"):
        mols = _molblock_records(text, first_only)
    elif fmt == "cdxml":
        mols = _cdxml_records(text)
        if first_only:
            mols = mols[:1]
    elif fmt == "name":
        mols = _name_records(text, first_only)
    elif fmt in ("rxn", "reaction_smiles"):
        rxn = _parse_reaction(text, fmt, prepare)
        return rxn.components[:1] if first_only else rxn.components
    else:  # pragma: no cover - guarded by _norm_fmt
        raise ChemError(f"Unsupported format {fmt!r}.")
    if prepare:
        mols = [prepare_for_depiction(m) for m in mols]
    return mols


def mol_from_input(data: str, fmt: Optional[str] = None, prepare: bool = True) -> Chem.Mol:
    """Parse one molecule from SMILES, InChI, molblock/SDF, CDXML, a file path or a name.

    For multi-record input (SDF, CDXML, SMILES lists) the first molecule is
    returned.  Reactions raise ``ChemError`` (use ``parse_reaction``).  With
    ``prepare=True`` the result is passed through ``prepare_for_depiction``.
    """
    text, f = _load(data, fmt)
    if f in ("rxn", "reaction_smiles"):
        raise ChemError(
            "Input is a reaction, not a single structure; use parse_reaction() "
            "or mols_from_input() to get its components."
        )
    mols = _records(text, f, first_only=True, prepare=prepare)
    if not mols:
        raise ChemError(f"No molecule found in {f} input.")
    return mols[0]


def mols_from_input(data: str, fmt: Optional[str] = None, prepare: bool = True) -> List[Chem.Mol]:
    """Parse all molecules from multi-record input.

    Handles SDF (all records), CDXML (all fragments), newline- or
    whitespace-separated SMILES lists (``SMILES name`` per line sets ``_Name``),
    InChI lists, lists of names and reactions (all components in the order
    reactants, agents, products).
    """
    text, f = _load(data, fmt)
    mols = _records(text, f, first_only=False, prepare=prepare)
    if not mols:
        raise ChemError(f"No molecules found in {f} input.")
    return mols


# --------------------------------------------------------------------------- #
# Depiction preparation
# --------------------------------------------------------------------------- #
def _is_flat(conf: Any) -> bool:
    """True when the conformer holds usable 2D coordinates (z == 0, not all zero)."""
    import numpy as np

    pos = np.asarray(conf.GetPositions())
    if pos.size == 0 or not np.all(np.isfinite(pos)):
        return False
    if pos.shape[0] > 1 and not np.any(np.abs(pos[:, :2]) > 1e-6):
        return False
    return bool(np.all(np.abs(pos[:, 2]) < 1e-4))


def _separate_fragments(mol: Chem.Mol, gap: float = 2.0) -> None:
    """Lay disconnected fragments out left-to-right so they do not overlap."""
    import numpy as np
    from rdkit.Geometry import Point3D

    Chem = _rd().Chem
    frags = Chem.GetMolFrags(mol)
    if len(frags) < 2:
        return
    conf = mol.GetConformer()
    pos = np.asarray(conf.GetPositions())
    cursor = 0.0
    for frag in frags:
        idx = list(frag)
        sub = pos[idx]
        xmin, ymin = sub[:, 0].min(), sub[:, 1].min()
        xmax, ymax = sub[:, 0].max(), sub[:, 1].max()
        dx, dy = cursor - xmin, -(ymin + ymax) / 2.0
        for i in idx:
            conf.SetAtomPosition(i, Point3D(float(pos[i, 0] + dx), float(pos[i, 1] + dy), 0.0))
        cursor += (xmax - xmin) + gap


def _coords_usable(mol: Chem.Mol) -> bool:
    """True when the molecule has one conformer with a finite position per atom."""
    import numpy as np

    if mol.GetNumConformers() == 0:
        return mol.GetNumAtoms() == 0
    pos = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    return pos.shape[0] == mol.GetNumAtoms() and bool(np.all(np.isfinite(pos)))


def _fallback_coords(mol: Chem.Mol) -> None:
    """Trivial 2D layout: one atom at the origin, further atoms 1.5 apart on x.

    Used when RDKit's depiction code yields NaN or no coordinates, which happens
    for all-hydrogen species such as ``[H][H]``, ``[H+]``, ``[H-]`` or ``[2H][2H]``.
    """
    from rdkit.Geometry import Point3D

    Chem = _rd().Chem
    mol.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(1.5 * i, 0.0, 0.0))
    conf.Set3D(False)
    mol.AddConformer(conf, assignId=True)


def _compute_2d(mol: Chem.Mol) -> None:
    """Compute 2D coordinates (CoordGen preferred) and spread fragments.

    Falls back to :func:`_fallback_coords` when the depiction code produces
    NaN or missing positions (all-hydrogen species).
    """
    rd = _rd()
    dep = rd.rdDepictor
    try:
        dep.SetPreferCoordGen(True)
    except Exception:  # pragma: no cover - old RDKit
        pass
    with _quiet():
        try:
            dep.Compute2DCoords(mol, clearConfs=True)
        except Exception:
            try:
                dep.SetPreferCoordGen(False)
                dep.Compute2DCoords(mol, clearConfs=True)
            finally:
                try:
                    dep.SetPreferCoordGen(True)
                except Exception:  # pragma: no cover
                    pass
    if not _coords_usable(mol):
        _fallback_coords(mol)
    if mol.GetNumConformers():
        mol.GetConformer().Set3D(False)
    _separate_fragments(mol)


def prepare_for_depiction(
    mol: Chem.Mol, kekulize: bool = True, keep_coords: bool = True, wedge: bool = True
) -> Chem.Mol:
    """Return a copy ready for 2D drawing / CDXML export.

    Steps: sanitize (errors become ``ChemError``), perceive stereo from 3D
    coordinates if present, remove non-essential explicit hydrogens
    (``Chem.RemoveHs`` defaults keep isotopic, stereo-defining and isolated
    ones), keep existing 2D coordinates when ``keep_coords`` (3D or missing
    coordinates are recomputed with CoordGen when available), assign
    stereochemistry, kekulize and clear aromatic flags, and set
    ``BEGINWEDGE``/``BEGINDASH`` directions on one bond per stereocentre,
    starting at the chiral atom.  Implicit hydrogen counts are up to date on
    the returned molecule.
    """
    rd = _rd()
    Chem = rd.Chem
    if mol is None:
        raise ChemError("No molecule given.")
    m = Chem.Mol(mol)
    try:
        m.UpdatePropertyCache(strict=False)
        with _quiet():
            Chem.SanitizeMol(m)
    except Exception as exc:
        raise ChemError(f"Sanitization failed: {_exc_msg(exc)}") from exc

    has_conf = m.GetNumConformers() > 0
    flat = has_conf and _is_flat(m.GetConformer())
    if has_conf and not flat and m.GetConformer().Is3D():
        with _quiet():
            Chem.AssignStereochemistryFrom3D(m, replaceExistingTags=True)

    try:
        with _quiet():
            m = Chem.RemoveHs(m)
    except Exception as exc:
        raise ChemError(f"Hydrogen removal failed: {_exc_msg(exc)}") from exc

    has_conf = m.GetNumConformers() > 0
    flat = has_conf and _is_flat(m.GetConformer())
    if keep_coords and flat:
        m.GetConformer().Set3D(False)
    else:
        m.RemoveAllConformers()
        _compute_2d(m)

    with _quiet():
        Chem.AssignStereochemistry(m, cleanIt=True, force=False, flagPossibleStereoCenters=True)

    if kekulize:
        try:
            Chem.Kekulize(m, clearAromaticFlags=True)
        except Exception as exc:
            raise ChemError(f"Kekulization failed: {_exc_msg(exc)}") from exc

    if wedge and m.GetNumConformers():
        wedge_dirs = (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH, Chem.BondDir.UNKNOWN)
        for bond in m.GetBonds():
            if bond.GetBondDir() in wedge_dirs:
                bond.SetBondDir(Chem.BondDir.NONE)
        with _quiet():
            Chem.WedgeMolBonds(m, m.GetConformer())

    m.UpdatePropertyCache(strict=False)
    return m


def _aromatic_copy(mol: Chem.Mol) -> Chem.Mol:
    """Sanitized (aromatic) copy used for SMILES, InChI and descriptors."""
    Chem = _rd().Chem
    m = Chem.Mol(mol)
    try:
        with _quiet():
            Chem.SanitizeMol(m)
    except Exception:
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:
            pass
    return m


def _as_mol(obj: Any, prepare: bool = True) -> Chem.Mol:
    """Accept a Mol or an input string."""
    if isinstance(obj, str):
        return mol_from_input(obj, prepare=prepare)
    if obj is None or not hasattr(obj, "GetNumAtoms"):
        raise ChemError(f"Expected an RDKit Mol or an input string, got {type(obj).__name__}.")
    return obj


# --------------------------------------------------------------------------- #
# Properties
# --------------------------------------------------------------------------- #
def _hill_order(symbols: Sequence[str]) -> List[str]:
    syms = set(symbols)
    if "C" in syms:
        rest = sorted(syms - {"C", "H"})
        return ["C"] + (["H"] if "H" in syms else []) + rest
    return sorted(syms)


def _element_masses(mol: Chem.Mol) -> Dict[str, float]:
    """Mass contribution per element symbol (implicit H included)."""
    Chem = _rd().Chem
    pt = Chem.GetPeriodicTable()
    masses: Dict[str, float] = {}
    h_weight = pt.GetAtomicWeight(1)
    for atom in mol.GetAtoms():
        num = atom.GetAtomicNum()
        if num == 0:
            continue
        iso = atom.GetIsotope()
        mass = pt.GetMassForIsotope(num, iso) if iso else pt.GetAtomicWeight(num)
        sym = atom.GetSymbol()
        masses[sym] = masses.get(sym, 0.0) + mass
        n_h = atom.GetTotalNumHs()
        if n_h:
            masses["H"] = masses.get("H", 0.0) + n_h * h_weight
    return masses


def _element_counts(mol: Chem.Mol) -> Counter:
    counts: Counter = Counter()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            continue
        counts[atom.GetSymbol()] += 1
        n_h = atom.GetTotalNumHs()
        if n_h:
            counts["H"] += n_h
    return counts


def elemental_analysis(mol: Chem.Mol) -> str:
    """ChemDraw-style elemental analysis, e.g. ``"C, 60.00; H, 4.48; O, 35.52"``."""
    mol = _as_mol(mol)
    masses = _element_masses(mol)
    total = sum(masses.values())
    if total <= 0:
        return ""
    return "; ".join(f"{el}, {100.0 * masses[el] / total:.2f}" for el in _hill_order(list(masses)))


def _stereo_summary(mol: Chem.Mol) -> Dict[str, Any]:
    Chem = _rd().Chem
    with _quiet():
        centers = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, includeCIP=True, useLegacyImplementation=False
        )
    labels = {int(i): str(lab) for i, lab in centers}
    assigned = sum(1 for lab in labels.values() if lab != "?")
    return {
        "assigned": assigned,
        "unassigned": len(labels) - assigned,
        "atoms": sorted(labels),
        "labels": labels,
    }


def _double_bond_stereo(mol: Chem.Mol) -> Dict[str, int]:
    Chem = _rd().Chem
    assigned = unassigned = 0
    try:
        with _quiet():
            infos = Chem.FindPotentialStereo(mol)
        for info in infos:
            if info.type == Chem.StereoType.Bond_Double:
                if info.specified == Chem.StereoSpecified.Specified:
                    assigned += 1
                else:
                    unassigned += 1
    except Exception:  # pragma: no cover - very old RDKit
        pass
    return {"assigned": assigned, "unassigned": unassigned}


def properties(mol: Chem.Mol, extended: bool = True) -> Dict[str, Any]:
    """Compute identifiers and properties of a molecule.

    Always: formula (Hill), molecular_weight, exact_mass, canonical_smiles,
    isomeric_smiles, inchi, inchikey, heavy_atoms, atoms (with H), bonds,
    charge, elemental_analysis, rings, aromatic_rings, stereocenters,
    double_bond_stereo.  With ``extended``: logp, tpsa, hbd, hba,
    rotatable_bonds, fraction_csp3, heavy_atom_count, lipinski_violations, qed.
    """
    rd = _rd()
    Chem = rd.Chem
    mol = _as_mol(mol)
    ms = _aromatic_copy(mol)
    D = rd.rdMolDescriptors
    with _quiet():
        try:
            inchi = Chem.MolToInchi(ms) or None
            inchikey = Chem.MolToInchiKey(ms) if inchi else None
        except Exception:
            inchi = inchikey = None
        try:
            mw = float(rd.Descriptors.MolWt(ms))
            exact = float(rd.Descriptors.ExactMolWt(ms))
        except Exception as exc:
            raise ChemError(f"Could not compute molecular weight: {_exc_msg(exc)}") from exc
    out: Dict[str, Any] = {
        "formula": D.CalcMolFormula(ms),
        "molecular_weight": round(mw, 3),
        "exact_mass": round(exact, 5),
        "canonical_smiles": Chem.MolToSmiles(ms, isomericSmiles=False),
        "isomeric_smiles": Chem.MolToSmiles(ms),
        "inchi": inchi,
        "inchikey": inchikey,
        "heavy_atoms": ms.GetNumHeavyAtoms(),
        "atoms": sum(1 + a.GetTotalNumHs() for a in ms.GetAtoms()),
        "bonds": ms.GetNumBonds(),
        "charge": Chem.GetFormalCharge(ms),
        "elemental_analysis": elemental_analysis(ms),
        "rings": D.CalcNumRings(ms),
        "aromatic_rings": D.CalcNumAromaticRings(ms),
        "stereocenters": _stereo_summary(ms),
        "double_bond_stereo": _double_bond_stereo(ms),
    }
    if ms.HasProp("_Name") and ms.GetProp("_Name").strip():
        out["name"] = ms.GetProp("_Name").strip()
    if extended:
        with _quiet():
            logp = float(rd.Crippen.MolLogP(ms))
            hbd = int(rd.Lipinski.NumHDonors(ms))
            hba = int(rd.Lipinski.NumHAcceptors(ms))
            out.update({
                "logp": round(logp, 3),
                "tpsa": round(float(D.CalcTPSA(ms)), 2),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": int(D.CalcNumRotatableBonds(ms)),
                "fraction_csp3": round(float(D.CalcFractionCSP3(ms)), 3),
                "heavy_atom_count": ms.GetNumHeavyAtoms(),
                "lipinski_violations": int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10),
            })
            try:
                from rdkit.Chem import QED

                out["qed"] = round(float(QED.qed(ms)), 3)
            except Exception:
                out["qed"] = None
    return out


def describe(mol: Chem.Mol) -> str:
    """One-line summary: name, formula, MW, SMILES, charge and stereocentres."""
    p = properties(mol, extended=False)
    st = p["stereocenters"]
    parts = [p["formula"], f"MW {p['molecular_weight']:.3f}", f"SMILES {p['isomeric_smiles']}"]
    parts.append(f"charge {p['charge']:+d}" if p["charge"] else "neutral")
    n = st["assigned"] + st["unassigned"]
    if n:
        parts.append(f"{n} stereocentre(s): {st['assigned']} assigned, {st['unassigned']} unassigned")
    else:
        parts.append("no stereocentres")
    prefix = f"{p['name']}: " if p.get("name") else ""
    return prefix + ", ".join(parts)


# --------------------------------------------------------------------------- #
# Reactions
# --------------------------------------------------------------------------- #
@dataclass
class Reaction:
    """A single-step reaction: reactants, agents (above the arrow) and products."""

    reactants: List[Chem.Mol] = field(default_factory=list)
    agents: List[Chem.Mol] = field(default_factory=list)
    products: List[Chem.Mol] = field(default_factory=list)
    name: str = ""
    conditions: str = ""

    @property
    def components(self) -> List[Chem.Mol]:
        """All molecules in the order reactants, agents, products."""
        return [*self.reactants, *self.agents, *self.products]


_UNICODE_ARROWS = ("→", "⟶", "⇒", "⟹", "⇌", "⇄", "⟷", "↔")
_ARROW_RE = re.compile(
    r"\s*-{1,3}\[(?P<cond>[^\]]*)\]-{0,3}>\s*"          # --[conditions]-->
    r"|\s*-{1,3}>\s*"                                     # -> or -->
    r"|\s*>>\s*"                                          # >>
    r"|\s*(?:→|⟶|⇒|⟹|⇌|⇄|⟷|↔|=>)\s*"                     # unicode arrows, =>
    r"|\s+=\s+"                                           # A + B = C
    r"|\s+>\s+"                                           # A + B > C
)
_COEFF_RE = re.compile(r"^(\d+)\s+(\S.*)$")


def _split_plus(side: str) -> List[str]:
    """Split on '+' outside square brackets (so [NH4+] survives)."""
    parts: List[str] = []
    depth = 0
    cur: List[str] = []
    for ch in side:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        if ch == "+" and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


def _split_dot(part: str) -> List[str]:
    """Split reaction-SMILES components on top-level '.', honouring '(A.B)' groups."""
    parts: List[str] = []
    depth = 0
    cur: List[str] = []
    for ch in part:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        if ch == "." and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    out: List[str] = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        if p.startswith("(") and p.endswith(")"):
            inner = p[1:-1]
            if _balanced(inner):
                p = inner
        out.append(p)
    return out


def _balanced(s: str) -> bool:
    depth = 0
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def _component_mols(spec: str, failures: List[str]) -> List[Chem.Mol]:
    """Turn 'SMILES', 'name' or '2 SMILES' into one or more raw Mol objects."""
    Chem = _rd().Chem
    spec = spec.strip()
    count = 1
    m = _COEFF_RE.match(spec)
    if m:
        count, spec = int(m.group(1)), m.group(2).strip()
    mol = None
    with _quiet():
        mol = Chem.MolFromSmiles(spec)
    if mol is None:
        smiles = resolve_name(spec)
        if smiles is not None:
            mol = _parse_smiles(smiles)
            mol.SetProp("_Name", spec)
    if mol is None:
        failures.append(spec)
        return []
    return [Chem.Mol(mol) for _ in range(max(1, count))]


def _components(specs: Sequence[str], failures: List[str]) -> List[Chem.Mol]:
    out: List[Chem.Mol] = []
    for spec in specs:
        out.extend(_component_mols(spec, failures))
    return out


def _finish_reaction(
    reactants: List[Chem.Mol], agents: List[Chem.Mol], products: List[Chem.Mol],
    name: str, conditions: str, prepare: bool,
) -> Reaction:
    Chem = _rd().Chem

    def conv(mols: Sequence[Chem.Mol]) -> List[Chem.Mol]:
        out: List[Chem.Mol] = []
        for m in mols:
            m = Chem.Mol(m)
            if prepare:
                out.append(prepare_for_depiction(m))
            else:
                try:
                    with _quiet():
                        Chem.SanitizeMol(m)
                except Exception:
                    m.UpdatePropertyCache(strict=False)
                out.append(m)
        return out

    return Reaction(conv(reactants), conv(agents), conv(products), name=name, conditions=conditions)


def _from_rdkit_reaction(rxn: Any, name: str, prepare: bool) -> Reaction:
    if rxn is None:
        raise ChemError("Could not parse the reaction.")
    return _finish_reaction(
        list(rxn.GetReactants()), list(rxn.GetAgents()), list(rxn.GetProducts()),
        name, "", prepare,
    )


def _rxn_name(text: str) -> str:
    lines = text.splitlines()
    return lines[1].strip() if len(lines) > 1 and lines[0].startswith("$RXN") else ""


def _reaction_from_cdxml_file(path: Path, prepare: bool) -> Reaction:
    RX = _rd().rdChemReactions
    reader = getattr(RX, "ReactionsFromCDXMLFile", None)
    if reader is None:  # pragma: no cover - old RDKit
        raise ChemError("This RDKit build cannot read reactions from CDXML.")
    try:
        with _quiet():
            rxns = list(reader(str(path)))
    except Exception as exc:
        raise ChemError(f"Could not read reaction from CDXML: {_exc_msg(exc)}") from exc
    if not rxns:
        raise ChemError("No reaction scheme found in the CDXML (it needs a <scheme>/<step>).")
    return _from_rdkit_reaction(rxns[0], "", prepare)


def _parse_reaction(text: str, fmt: Optional[str], prepare: bool) -> Reaction:
    rd = _rd()
    RX = rd.rdChemReactions
    if not isinstance(text, str) or not text.strip():
        raise ChemError("Empty reaction input.")
    f = (fmt or "").strip().lower() or None
    path = _as_path(text)
    if path is not None:
        suffix = path.suffix.lower()
        if suffix == ".rxn" or f == "rxn":
            content = _read_text(path)
            try:
                with _quiet():
                    rxn = RX.ReactionFromRxnFile(str(path))
            except Exception as exc:
                raise ChemError(f"Could not parse RXN file {path}: {_exc_msg(exc)}") from exc
            return _from_rdkit_reaction(rxn, _rxn_name(content), prepare)
        if suffix == ".cdxml" or f == "cdxml":
            return _reaction_from_cdxml_file(path, prepare)
        text = _read_text(path)
    s = text.strip()
    if f is None:
        if s.startswith("$RXN"):
            f = "rxn"
        elif _looks_like_cdxml(s):
            f = "cdxml"
        elif _is_reaction_smiles(s):
            f = "reaction_smiles"
        else:
            f = "text"
    if f == "rxn":
        try:
            with _quiet():
                rxn = RX.ReactionFromRxnBlock(s)
        except Exception as exc:
            raise ChemError(f"Could not parse RXN block: {_exc_msg(exc)}") from exc
        return _from_rdkit_reaction(rxn, _rxn_name(s), prepare)
    if f == "cdxml":
        fd, tmp = tempfile.mkstemp(suffix=".cdxml")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(s)
            return _reaction_from_cdxml_file(Path(tmp), prepare)
        finally:
            try:
                os.unlink(tmp)
            except OSError:  # pragma: no cover
                pass
    failures: List[str] = []
    if f in ("reaction_smiles", "smiles") and _is_reaction_smiles(s):
        core = s.split(" |", 1)[0].strip()
        parts = _split_reaction_smiles(core)
        if len(parts) != 3:
            raise ChemError("Reaction SMILES must have the form reactants>agents>products (got %d '>'-separated parts in %r)"
                            % (len(parts), core[:80]))
        r_part, a_part, p_part = parts
        reactants = _components(_split_dot(r_part), failures)
        agents = _components(_split_dot(a_part), failures)
        products = _components(_split_dot(p_part), failures)
        conditions = ""
    else:
        m = _ARROW_RE.search(s)
        if m is None:
            raise ChemError(
                "No reaction arrow found. Write e.g. 'A + B -> C', "
                "'A + B --[H2SO4, reflux]--> C' or reaction SMILES 'A.B>C>D'."
            )
        lhs, rhs = s[: m.start()], s[m.end():]
        if _ARROW_RE.search(rhs):
            raise ChemError("Multi-step reactions are not supported; give a single step.")
        conditions = (m.group("cond") or "").strip()
        reactants = _components(_split_plus(lhs), failures)
        agents = []
        products = _components(_split_plus(rhs), failures)
    if failures:
        raise ChemError(
            "Could not interpret reaction component(s): "
            + ", ".join(repr(x) for x in failures)
            + ". Use SMILES, or chemical names known to resolve_name()."
        )
    if not reactants or not products:
        raise ChemError("A reaction needs at least one reactant and one product.")
    return _finish_reaction(reactants, agents, products, "", conditions, prepare)


def parse_reaction(text: str, fmt: Optional[str] = None) -> Reaction:
    """Parse a reaction from reaction SMILES, an RXN block/file, CDXML or arrow text.

    Accepted: ``A.B>C>D`` (agents optional, ``A>>B``), ``$RXN`` blocks or a
    ``.rxn`` path, CDXML with a reaction scheme, and ``A + B -> C`` style text
    (arrows ``->``, ``-->``, ``>>``, ``=``, ``→`` ...; ``A + B --[reagents]--> C``
    stores the bracket text in ``Reaction.conditions``; a leading integer is a
    stoichiometric coefficient).  Components may be SMILES or names resolvable
    by ``resolve_name``; failures are reported in the ``ChemError``.  All
    components are passed through ``prepare_for_depiction``.
    """
    return _parse_reaction(text, fmt, prepare=True)


def _smiles_of(mol: Chem.Mol) -> str:
    return _rd().Chem.MolToSmiles(_aromatic_copy(mol))


def reaction_smiles(rxn: Reaction) -> str:
    """Reaction SMILES ``reactants>agents>products`` (canonical, aromatic)."""
    return ">".join(
        ".".join(_smiles_of(m) for m in group)
        for group in (rxn.reactants, rxn.agents, rxn.products)
    )


def reaction_properties(rxn: Reaction) -> Dict[str, Any]:
    """Per-component formula/MW/SMILES, element balance, masses and atom economy.

    ``atom_economy`` is the classic value: MW of the first (desired) product
    divided by the total reactant MW, in percent; it is ``None`` when the
    reaction is not balanced.
    """
    def info(m: Chem.Mol) -> Dict[str, Any]:
        ms = _aromatic_copy(m)
        rd = _rd()
        return {
            "formula": rd.rdMolDescriptors.CalcMolFormula(ms),
            "molecular_weight": round(float(rd.Descriptors.MolWt(ms)), 3),
            "smiles": rd.Chem.MolToSmiles(ms),
            "name": ms.GetProp("_Name") if ms.HasProp("_Name") else "",
        }

    Chem = _rd().Chem
    reactants = [info(m) for m in rxn.reactants]
    agents = [info(m) for m in rxn.agents]
    products = [info(m) for m in rxn.products]
    r_counts: Counter = Counter()
    p_counts: Counter = Counter()
    for m in rxn.reactants:
        r_counts.update(_element_counts(m))
    for m in rxn.products:
        p_counts.update(_element_counts(m))
    balance = {
        el: (int(r_counts.get(el, 0)), int(p_counts.get(el, 0)))
        for el in _hill_order(list(set(r_counts) | set(p_counts)))
    }
    q_r = sum(Chem.GetFormalCharge(m) for m in rxn.reactants)
    q_p = sum(Chem.GetFormalCharge(m) for m in rxn.products)
    unbalanced = [el for el, (a, b) in balance.items() if a != b]
    balanced = not unbalanced and q_r == q_p
    mass_r = round(sum(x["molecular_weight"] for x in reactants), 3)
    mass_p = round(sum(x["molecular_weight"] for x in products), 3)
    atom_economy = None
    if balanced and products and mass_r > 0:
        atom_economy = round(100.0 * products[0]["molecular_weight"] / mass_r, 1)
    return {
        "reactants": reactants,
        "agents": agents,
        "products": products,
        "reaction_smiles": reaction_smiles(rxn),
        "conditions": rxn.conditions,
        "name": rxn.name,
        "element_balance": balance,
        "unbalanced_elements": unbalanced,
        "charge_balance": (q_r, q_p),
        "balanced": balanced,
        "mass_reactants": mass_r,
        "mass_products": mass_p,
        "atom_economy": atom_economy,
    }


def _to_rdkit_reaction(rxn: Reaction) -> Any:
    rd = _rd()
    out = rd.rdChemReactions.ChemicalReaction()
    for m in rxn.reactants:
        out.AddReactantTemplate(rd.Chem.Mol(m))
    for m in rxn.agents:
        out.AddAgentTemplate(rd.Chem.Mol(m))
    for m in rxn.products:
        out.AddProductTemplate(rd.Chem.Mol(m))
    return out


# --------------------------------------------------------------------------- #
# Rendering and export
# --------------------------------------------------------------------------- #
def _drawable(obj: Any) -> Union[Chem.Mol, Reaction]:
    if isinstance(obj, Reaction):
        return obj
    if isinstance(obj, str):
        if detect_format(obj) in ("rxn", "reaction_smiles"):
            return parse_reaction(obj)
        return mol_from_input(obj)
    return _as_mol(obj)


def _draw(drawer: Any, obj: Any, legend: str, bond_length_px: Optional[float]) -> None:
    target = _drawable(obj)
    opts = drawer.drawOptions()
    if bond_length_px:
        opts.fixedBondLength = float(bond_length_px)
    with _quiet():
        if isinstance(target, Reaction):
            drawer.DrawReaction(_to_rdkit_reaction(target))
        else:
            drawer.DrawMolecule(_rd().Chem.Mol(target), legend=legend or "")
        drawer.FinishDrawing()


def depict_png(
    obj: Union[Chem.Mol, Reaction, str], width: int = 800, height: int = 400,
    legend: str = "", bond_length_px: Optional[float] = None,
) -> bytes:
    """Render a molecule or Reaction to PNG bytes with RDKit's Cairo drawer.

    ``bond_length_px`` fixes the drawn bond length; the legend is only drawn
    for molecules.
    """
    Draw = _rd().Draw
    if not hasattr(Draw, "MolDraw2DCairo"):  # pragma: no cover - build dependent
        raise ChemError("This RDKit build has no Cairo support; use depict_svg().")
    try:
        drawer = Draw.MolDraw2DCairo(int(width), int(height))
        _draw(drawer, obj, legend, bond_length_px)
        return bytes(drawer.GetDrawingText())
    except ChemError:
        raise
    except Exception as exc:
        raise ChemError(f"PNG rendering failed: {_exc_msg(exc)}") from exc


def depict_svg(
    obj: Union[Chem.Mol, Reaction, str], width: int = 800, height: int = 400,
    legend: str = "", bond_length_px: Optional[float] = None,
) -> str:
    """Render a molecule or Reaction to an SVG string."""
    Draw = _rd().Draw
    try:
        drawer = Draw.MolDraw2DSVG(int(width), int(height))
        _draw(drawer, obj, legend, bond_length_px)
        return str(drawer.GetDrawingText())
    except ChemError:
        raise
    except Exception as exc:
        raise ChemError(f"SVG rendering failed: {_exc_msg(exc)}") from exc


def _with_coords(mol: Chem.Mol) -> Chem.Mol:
    return mol if mol.GetNumConformers() else prepare_for_depiction(mol)


def mol_to_molblock(mol: Chem.Mol, v3000: bool = False) -> str:
    """MDL molblock (V2000, or V3000 when requested) with 2D coordinates and wedges."""
    Chem = _rd().Chem
    m = _with_coords(_as_mol(mol))
    try:
        with _quiet():
            return Chem.MolToMolBlock(m, forceV3000=bool(v3000))
    except Exception as exc:
        raise ChemError(f"Could not write molblock: {_exc_msg(exc)}") from exc


def mol_to_sdf(mols: Union[Chem.Mol, Sequence[Chem.Mol]]) -> str:
    """SDF text for one or more molecules (properties and ``_Name`` included)."""
    Chem = _rd().Chem
    if hasattr(mols, "GetNumAtoms"):
        mols = [mols]  # type: ignore[list-item]
    buf = StringIO()
    writer = Chem.SDWriter(buf)
    try:
        with _quiet():
            for m in mols:
                writer.write(_with_coords(_as_mol(m)))
    except ChemError:
        raise
    except Exception as exc:
        raise ChemError(f"Could not write SDF: {_exc_msg(exc)}") from exc
    finally:
        writer.close()
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# Substructure / similarity
# --------------------------------------------------------------------------- #
def substructure_matches(mol: Chem.Mol, smarts: str) -> List[List[int]]:
    """Atom-index lists of every match of a SMARTS (or SMILES) pattern."""
    Chem = _rd().Chem
    mol = _as_mol(mol)
    with _quiet():
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            query = Chem.MolFromSmiles(smarts)
    if query is None:
        raise ChemError(f"Invalid SMARTS pattern {smarts!r}.")
    ms = _aromatic_copy(mol)
    matches = ms.GetSubstructMatches(query, uniquify=True, maxMatches=100000)
    return [list(map(int, t)) for t in matches]


def tanimoto(a: Union[Chem.Mol, str], b: Union[Chem.Mol, str], radius: int = 2) -> float:
    """Tanimoto similarity of Morgan (ECFP-like) fingerprints, 2048 bits."""
    rd = _rd()
    ma, mb = _aromatic_copy(_as_mol(a)), _aromatic_copy(_as_mol(b))
    try:
        from rdkit.Chem import rdFingerprintGenerator

        gen = rdFingerprintGenerator.GetMorganGenerator(radius=int(radius), fpSize=2048)
        fa, fb = gen.GetFingerprint(ma), gen.GetFingerprint(mb)
    except Exception:  # pragma: no cover - old RDKit
        fa = rd.AllChem.GetMorganFingerprintAsBitVect(ma, int(radius), nBits=2048)
        fb = rd.AllChem.GetMorganFingerprintAsBitVect(mb, int(radius), nBits=2048)
    return float(rd.DataStructs.TanimotoSimilarity(fa, fb))


# --------------------------------------------------------------------------- #
# Name resolution
# --------------------------------------------------------------------------- #
# Offline dictionary of common names: ((alias, ...), SMILES), grouped in sections.
#
# Conventions
# * Aliases are lower-case and already normalised (see ``_normalize_name``);
#   lookups also try a "loose" key without spaces, hyphens, dots and quotes, so
#   "n-BuLi", "nBuLi", "EDC·HCl", "EDC HCl" and "Hünig's base" all resolve.
# * Amino acids are L (CIP S, cysteine/selenocysteine R) unless prefixed "d-";
#   "dl-"/"rac-" names return the structure without tetrahedral stereo.
# * Sugars are the usual ring forms: D-glucose = beta-D-glucopyranose, galactose /
#   mannose / allose / talose = beta-D-pyranoses, ribose = beta-D-ribofuranose,
#   fructose = beta-D-fructofuranose; the pentoses arabinose, xylose and lyxose
#   are open-chain aldehydes.  "d-glucose (open chain)" is available explicitly.
# * Hydrate water is not part of the structure: "monohydrate", "·5H2O" etc. are
#   stripped on lookup (a few entries such as CuSO4·5H2O are listed with water).
# * Drugs and natural products carry stereocentres only where they are certain;
#   racemates and uncertain cases are flat.  Metal complexes are written as
#   simple salts / fragment lists (e.g. Pd(PPh3)4 as [Pd] + 4 PPh3); solid
#   catalysts such as Pd/C resolve to the metal.
_COMMON: Sequence[Tuple[Tuple[str, ...], str]] = (
    # ------------------------------------------------------------------ #
    # Solvents (protic, aprotic, chlorinated, ethers, hydrocarbons)
    # ------------------------------------------------------------------ #
    (("water", "h2o", "wasser"), "O"),
    (("heavy water", "d2o", "deuterium oxide", "water-d2"), "[2H]O[2H]"),
    (("methanol", "meoh", "methyl alcohol", "wood alcohol"), "CO"),
    (("ethanol", "etoh", "ethyl alcohol", "alcohol"), "CCO"),
    (("isopropanol", "isopropyl alcohol", "2-propanol", "ipa", "iproh", "propan-2-ol"), "CC(C)O"),
    (("n-propanol", "1-propanol", "propanol", "propan-1-ol", "nproh"), "CCCO"),
    (("n-butanol", "1-butanol", "butanol", "butan-1-ol", "nbuoh"), "CCCCO"),
    (("2-butanol", "sec-butanol", "butan-2-ol"), "CCC(C)O"),
    (("isobutanol", "2-methyl-1-propanol", "isobutyl alcohol"), "CC(C)CO"),
    (("tert-butanol", "t-butanol", "tbuoh", "t-buoh", "2-methyl-2-propanol", "tert-butyl alcohol"),
     "CC(C)(C)O"),
    (("tert-amyl alcohol", "t-amyl alcohol", "2-methyl-2-butanol", "tamyloh"), "CCC(C)(C)O"),
    (("1-pentanol", "n-pentanol", "pentanol", "amyl alcohol"), "CCCCCO"),
    (("1-hexanol", "n-hexanol", "hexanol"), "CCCCCCO"),
    (("1-octanol", "n-octanol", "octanol"), "CCCCCCCCO"),
    (("cyclohexanol",), "OC1CCCCC1"),
    (("allyl alcohol", "2-propen-1-ol"), "C=CCO"),
    (("propargyl alcohol", "2-propyn-1-ol"), "C#CCO"),
    (("trifluoroethanol", "tfe", "2,2,2-trifluoroethanol"), "OCC(F)(F)F"),
    (("hexafluoroisopropanol", "hfip", "1,1,1,3,3,3-hexafluoro-2-propanol", "hexafluoro-2-propanol"),
     "OC(C(F)(F)F)C(F)(F)F"),
    (("ethylene glycol", "1,2-ethanediol", "glycol"), "OCCO"),
    (("propylene glycol", "1,2-propanediol"), "CC(O)CO"),
    (("1,3-propanediol",), "OCCCO"),
    (("1,4-butanediol",), "OCCCCO"),
    (("diethylene glycol", "deg"), "OCCOCCO"),
    (("triethylene glycol",), "OCCOCCOCCO"),
    (("2-methoxyethanol", "methyl cellosolve", "egme", "ethylene glycol monomethyl ether"), "COCCO"),
    (("2-ethoxyethanol", "cellosolve"), "CCOCCO"),
    (("glycerol", "glycerin", "glycerine", "propane-1,2,3-triol"), "OCC(O)CO"),
    (("acetone", "propanone", "propan-2-one", "dimethyl ketone"), "CC(C)=O"),
    (("2-butanone", "methyl ethyl ketone", "mek", "butanone"), "CCC(C)=O"),
    (("methyl isobutyl ketone", "mibk", "4-methyl-2-pentanone"), "CC(C)CC(C)=O"),
    (("acetonitrile", "mecn", "acn", "methyl cyanide", "cyanomethane", "ch3cn"), "CC#N"),
    (("propionitrile", "propanenitrile"), "CCC#N"),
    (("benzonitrile", "phcn"), "N#Cc1ccccc1"),
    (("ethyl acetate", "etoac", "ethyl ethanoate", "ea"), "CCOC(C)=O"),
    (("methyl acetate", "meoac"), "COC(C)=O"),
    (("isopropyl acetate", "ipac", "iproac"), "CC(C)OC(C)=O"),
    (("butyl acetate", "n-butyl acetate", "buoac"), "CCCCOC(C)=O"),
    (("ethyl formate",), "CCOC=O"),
    (("methyl formate",), "COC=O"),
    (("dmf", "n,n-dimethylformamide", "dimethylformamide"), "CN(C)C=O"),
    (("dma", "dmac", "n,n-dimethylacetamide", "dimethylacetamide"), "CC(=O)N(C)C"),
    (("nmp", "n-methyl-2-pyrrolidone", "n-methylpyrrolidone", "1-methyl-2-pyrrolidinone", "n-methylpyrrolidinone"),
     "CN1CCCC1=O"),
    (("dmso", "dimethyl sulfoxide", "dimethylsulfoxide", "methylsulfinylmethane", "dimethylsulphoxide", "dimethyl sulphoxide", "me2so"),
     "CS(C)=O"),
    (("sulfolane", "tetramethylene sulfone", "tetrahydrothiophene 1,1-dioxide"), "O=S1(=O)CCCC1"),
    (("hmpa", "hexamethylphosphoramide", "hexamethylphosphoric triamide"), "CN(C)P(=O)(N(C)C)N(C)C"),
    (("dmpu", "n,n'-dimethylpropyleneurea", "1,3-dimethyl-3,4,5,6-tetrahydro-2(1h)-pyrimidinone"),
     "CN1CCCN(C)C1=O"),
    (("dmi", "1,3-dimethyl-2-imidazolidinone", "n,n'-dimethylethyleneurea"), "CN1CCN(C)C1=O"),
    (("formamide",), "NC=O"),
    (("n-methylformamide", "nmf"), "CNC=O"),
    (("acetamide",), "CC(N)=O"),
    (("propylene carbonate", "4-methyl-1,3-dioxolan-2-one"), "CC1COC(=O)O1"),
    (("ethylene carbonate", "1,3-dioxolan-2-one"), "O=C1OCCO1"),
    (("dimethyl carbonate", "dmc", "(meo)2co"), "COC(=O)OC"),
    (("diethyl carbonate", "dec", "(eto)2co"), "CCOC(=O)OCC"),
    (("thf", "tetrahydrofuran", "oxolane"), "C1CCOC1"),
    (("2-methyltetrahydrofuran", "2-methf", "2-me-thf", "2-methyl-thf", "2-methyloxolane", "methf"),
     "CC1CCCO1"),
    (("1,4-dioxane", "dioxane", "p-dioxane"), "C1COCCO1"),
    (("1,3-dioxolane", "dioxolane"), "C1COCO1"),
    (("tetrahydropyran", "thp", "oxane"), "C1CCOCC1"),
    (("3,4-dihydro-2h-pyran", "dihydropyran", "dhp"), "C1CC=COC1"),
    (("diethyl ether", "ether", "et2o", "ethyl ether", "ethoxyethane"), "CCOCC"),
    (("diisopropyl ether", "ipr2o", "isopropyl ether"), "CC(C)OC(C)C"),
    (("dibutyl ether", "di-n-butyl ether", "bu2o", "butyl ether"), "CCCCOCCCC"),
    (("mtbe", "methyl tert-butyl ether", "tert-butyl methyl ether", "tbme", "methyl t-butyl ether"),
     "COC(C)(C)C"),
    (("cpme", "cyclopentyl methyl ether"), "COC1CCCC1"),
    (("dme", "1,2-dimethoxyethane", "glyme", "monoglyme", "ethylene glycol dimethyl ether"), "COCCOC"),
    (("diglyme", "diethylene glycol dimethyl ether", "bis(2-methoxyethyl) ether"), "COCCOCCOC"),
    (("triglyme", "triethylene glycol dimethyl ether"), "COCCOCCOCCOC"),
    (("dimethyl ether", "methoxymethane", "me2o"), "COC"),
    (("dichloromethane", "dcm", "methylene chloride", "ch2cl2", "methylene dichloride"), "ClCCl"),
    (("chloroform", "trichloromethane", "chcl3"), "ClC(Cl)Cl"),
    (("carbon tetrachloride", "tetrachloromethane", "ccl4", "tetrachloride"), "ClC(Cl)(Cl)Cl"),
    (("1,2-dichloroethane", "dce", "ethylene dichloride"), "ClCCCl"),
    (("1,1,2,2-tetrachloroethane",), "ClC(Cl)C(Cl)Cl"),
    (("chloroethane", "ethyl chloride"), "CCCl"),
    (("1,2-dibromoethane", "ethylene dibromide"), "BrCCBr"),
    (("dibromomethane", "methylene bromide"), "BrCBr"),
    (("bromoform", "tribromomethane"), "BrC(Br)Br"),
    (("iodoform", "triiodomethane"), "IC(I)I"),
    (("trichloroethylene", "tce", "trichloroethene"), "ClC=C(Cl)Cl"),
    (("perchloroethylene", "tetrachloroethylene", "perc"), "ClC(Cl)=C(Cl)Cl"),
    (("1,2-dichlorobenzene", "o-dichlorobenzene", "odcb"), "Clc1ccccc1Cl"),
    (("chlorobenzene", "phcl"), "Clc1ccccc1"),
    (("fluorobenzene", "phf"), "Fc1ccccc1"),
    (("bromobenzene", "phbr"), "Brc1ccccc1"),
    (("iodobenzene", "phi"), "Ic1ccccc1"),
    (("benzotrifluoride", "trifluorotoluene", "(trifluoromethyl)benzene", "alpha,alpha,alpha-trifluorotoluene", "α,α,α-trifluorotoluene", "phcf3"),
     "FC(F)(F)c1ccccc1"),
    (("pentane", "n-pentane"), "CCCCC"),
    (("isopentane", "2-methylbutane"), "CCC(C)C"),
    (("neopentane", "2,2-dimethylpropane"), "CC(C)(C)C"),
    (("hexane", "n-hexane", "hexanes"), "CCCCCC"),
    (("heptane", "n-heptane"), "CCCCCCC"),
    (("octane", "n-octane"), "CCCCCCCC"),
    (("isooctane", "2,2,4-trimethylpentane"), "CC(C)CC(C)(C)C"),
    (("nonane", "n-nonane"), "CCCCCCCCC"),
    (("decane", "n-decane"), "CCCCCCCCCC"),
    (("dodecane", "n-dodecane"), "CCCCCCCCCCCC"),
    (("hexadecane", "cetane", "n-hexadecane"), "CCCCCCCCCCCCCCCC"),
    (("cyclopentane",), "C1CCCC1"),
    (("cyclohexane",), "C1CCCCC1"),
    (("methylcyclohexane",), "CC1CCCCC1"),
    (("cycloheptane",), "C1CCCCCC1"),
    (("cyclooctane",), "C1CCCCCCC1"),
    (("decalin", "decahydronaphthalene"), "C1CCC2CCCCC2C1"),
    (("tetralin", "1,2,3,4-tetrahydronaphthalene"), "c1ccc2CCCCc2c1"),
    (("benzene", "phh", "c6h6"), "c1ccccc1"),
    (("toluene", "methylbenzene", "phme", "phch3", "toluol"), "Cc1ccccc1"),
    (("o-xylene", "1,2-dimethylbenzene"), "Cc1ccccc1C"),
    (("m-xylene", "1,3-dimethylbenzene"), "Cc1cccc(C)c1"),
    (("p-xylene", "xylene", "1,4-dimethylbenzene", "xylenes"), "Cc1ccc(C)cc1"),
    (("mesitylene", "1,3,5-trimethylbenzene"), "Cc1cc(C)cc(C)c1"),
    (("ethylbenzene",), "CCc1ccccc1"),
    (("cumene", "isopropylbenzene"), "CC(C)c1ccccc1"),
    (("tert-butylbenzene", "t-butylbenzene"), "CC(C)(C)c1ccccc1"),
    (("anisole", "methoxybenzene", "methyl phenyl ether"), "COc1ccccc1"),
    (("phenetole", "ethoxybenzene"), "CCOc1ccccc1"),
    (("diphenyl ether", "diphenyl oxide", "ph2o"), "O(c1ccccc1)c1ccccc1"),
    (("veratrole", "1,2-dimethoxybenzene"), "COc1ccccc1OC"),
    (("pyridine", "py"), "c1ccncc1"),
    (("2-picoline", "2-methylpyridine", "alpha-picoline", "α-picoline"), "Cc1ccccn1"),
    (("3-picoline", "3-methylpyridine"), "Cc1cccnc1"),
    (("4-picoline", "4-methylpyridine"), "Cc1ccncc1"),
    (("2,6-lutidine", "2,6-dimethylpyridine", "lutidine"), "Cc1cccc(C)n1"),
    (("2,4,6-collidine", "collidine", "2,4,6-trimethylpyridine", "sym-collidine"), "Cc1cc(C)nc(C)c1"),
    (("nitromethane", "meno2"), "C[N+](=O)[O-]"),
    (("nitroethane",), "CC[N+](=O)[O-]"),
    (("carbon disulfide", "carbon disulphide", "cs2"), "S=C=S"),
    (("dimethyl sulfide", "dms", "methyl sulfide", "dimethylsulfide", "me2s"), "CSC"),
    (("thioanisole", "methyl phenyl sulfide"), "CSc1ccccc1"),
    (("ethyl lactate",), "CCOC(=O)C(C)O"),
    (("triacetin", "glyceryl triacetate"), "CC(=O)OCC(COC(C)=O)OC(C)=O"),
    (("cyrene", "dihydrolevoglucosenone"), "O=C1C[C@@H]2OC[C@H]1O2"),
    (("gamma-valerolactone", "gvl", "γ-valerolactone"), "CC1CCC(=O)O1"),
    (("gamma-butyrolactone", "gbl", "γ-butyrolactone", "butyrolactone"), "O=C1CCCO1"),
    (("delta-valerolactone", "δ-valerolactone"), "O=C1CCCCO1"),
    (("epsilon-caprolactone", "ε-caprolactone", "caprolactone"), "O=C1CCCCCO1"),
    (("epsilon-caprolactam", "ε-caprolactam", "caprolactam"), "O=C1CCCCCN1"),
    (("2-pyrrolidone", "2-pyrrolidinone", "pyrrolidone"), "O=C1CCCN1"),
    # ------------------------------------------------------------------ #
    # Deuterated NMR solvents
    # ------------------------------------------------------------------ #
    (("cdcl3", "chloroform-d", "deuterochloroform", "deuterated chloroform", "chloroform-d1"),
     "[2H]C(Cl)(Cl)Cl"),
    (("dmso-d6", "d6-dmso", "deuterated dmso", "dimethyl sulfoxide-d6", "(cd3)2so"),
     "[2H]C([2H])([2H])S(=O)C([2H])([2H])[2H]"),
    (("methanol-d4", "cd3od", "d4-methanol", "deuterated methanol", "meod", "meoh-d4"),
     "[2H]OC([2H])([2H])[2H]"),
    (("methanol-d3", "cd3oh"), "OC([2H])([2H])[2H]"),
    (("acetone-d6", "d6-acetone", "(cd3)2co", "deuterated acetone"),
     "[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H]"),
    (("acetonitrile-d3", "cd3cn", "d3-acetonitrile", "mecn-d3"), "[2H]C([2H])([2H])C#N"),
    (("benzene-d6", "c6d6", "d6-benzene", "deuterated benzene"), "[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H]"),
    (("toluene-d8", "d8-toluene"), "[2H]C([2H])([2H])c1c([2H])c([2H])c([2H])c([2H])c1[2H]"),
    (("thf-d8", "d8-thf", "tetrahydrofuran-d8"), "[2H]C1([2H])OC([2H])([2H])C([2H])([2H])C1([2H])[2H]"),
    (("dichloromethane-d2", "cd2cl2", "dcm-d2", "methylene chloride-d2"), "[2H]C([2H])(Cl)Cl"),
    (("pyridine-d5", "d5-pyridine", "c5d5n"), "[2H]c1nc([2H])c([2H])c([2H])c1[2H]"),
    (("dmf-d7", "d7-dmf", "n,n-dimethylformamide-d7"), "[2H]C(=O)N(C([2H])([2H])[2H])C([2H])([2H])[2H]"),
    (("trifluoroacetic acid-d", "tfa-d", "cf3cood", "tfa-d1"), "[2H]OC(=O)C(F)(F)F"),
    (("acetic acid-d4", "cd3cood", "d4-acetic acid"), "[2H]OC(=O)C([2H])([2H])[2H]"),
    (("deuterium", "d2"), "[2H][2H]"),
    (("methane", "ch4", "natural gas"), "C"),
    (("ethane", "c2h6"), "CC"),
    (("propane", "c3h8"), "CCC"),
    (("butane", "n-butane", "c4h10"), "CCCC"),
    (("isobutane", "2-methylpropane"), "CC(C)C"),
    (("ethylene", "ethene", "c2h4"), "C=C"),
    (("acetylene", "ethyne", "c2h2"), "C#C"),
    (("styrene", "vinylbenzene", "phenylethylene"), "C=Cc1ccccc1"),
    (("phenol", "hydroxybenzene", "phoh", "carbolic acid"), "Oc1ccccc1"),
    (("imidazole", "1h-imidazole", "1,3-diazole"), "c1cnc[nH]1"),
    (("pyrrole", "1h-pyrrole"), "c1cc[nH]c1"),
    (("furan",), "c1ccoc1"),
    (("thiophene",), "c1ccsc1"),
    (("indole", "1h-indole"), "c1ccc2[nH]ccc2c1"),
    (("indoline", "2,3-dihydroindole"), "C1Cc2ccccc2N1"),
    (("thiazolidine",), "C1CSCN1"),
    (("oxazolidine",), "C1COCN1"),
    (("2-methylfuran", "sylvan"), "Cc1ccco1"),
    (("2,5-dimethylfuran",), "Cc1ccc(C)o1"),
    (("tetrahydrothiophene", "thiolane"), "C1CCSC1"),
    (("benzofuran",), "c1ccc2occc2c1"),
    (("benzothiophene", "benzo[b]thiophene"), "c1ccc2sccc2c1"),
    (("isoindoline",), "C1NCc2ccccc21"),
    (("1,2,3,4-tetrahydroisoquinoline", "tetrahydroisoquinoline", "thiq"), "C1Cc2ccccc2CN1"),
    (("1,2,3,4-tetrahydroquinoline", "tetrahydroquinoline"), "C1CNc2ccccc2C1"),
    (("chroman", "chromane"), "C1COc2ccccc2C1"),
    (("coumaran", "2,3-dihydrobenzofuran"), "C1Oc2ccccc2C1"),
    (("2h-chromene", "chromene"), "C1C=Cc2ccccc2O1"),
    (("1,3-benzodioxole", "methylenedioxybenzene"), "c1ccc2OCOc2c1"),
    (("quinuclidinone", "3-quinuclidinone"), "O=C1CN2CCC1CC2"),
    (("tropinone",), "CN1[C@H]2CC[C@@H]1CC(=O)C2"),
    (("nortropinone",), "O=C1C[C@H]2CC[C@@H](C1)N2"),
    (("homopiperazine", "1,4-diazepane"), "C1CNCCNC1"),
    (("azetidin-3-ol hydrochloride",), "Cl.OC1CNC1"),
    (("n-methylpiperazine dihydrochloride",), "Cl.Cl.CN1CCNCC1"),
    # ------------------------------------------------------------------ #
    # Gases, elements and small inorganic molecules
    # ------------------------------------------------------------------ #
    (("hydrogen", "h2", "dihydrogen"), "[H][H]"),
    (("oxygen", "o2", "dioxygen"), "O=O"),
    (("nitrogen", "n2", "dinitrogen"), "N#N"),
    (("ozone", "o3"), "[O-][O+]=O"),
    (("fluorine", "f2"), "FF"),
    (("chlorine", "cl2"), "ClCl"),
    (("bromine", "br2"), "BrBr"),
    (("iodine", "i2"), "II"),
    (("sulfur", "sulphur", "s8", "elemental sulfur"), "S1SSSSSSS1"),
    (("argon", "ar"), "[Ar]"),
    (("helium", "he"), "[He]"),
    (("neon",), "[Ne]"),
    (("xenon", "xe"), "[Xe]"),
    (("carbon dioxide", "co2", "dry ice"), "O=C=O"),
    (("carbon monoxide", "co"), "[C-]#[O+]"),
    (("nitric oxide", "no", "nitrogen monoxide"), "[N]=O"),
    (("nitrous oxide", "n2o", "laughing gas", "dinitrogen monoxide"), "[N-]=[N+]=O"),
    (("nitrogen dioxide", "no2"), "[O-][N+]=O"),
    (("sulfur dioxide", "sulphur dioxide", "so2"), "O=S=O"),
    (("sulfur trioxide", "sulphur trioxide", "so3"), "O=S(=O)=O"),
    (("hydrogen sulfide", "hydrogen sulphide", "h2s", "sulfane"), "S"),
    (("hydrogen cyanide", "hcn", "prussic acid", "hydrocyanic acid"), "C#N"),
    (("cyanogen bromide", "cnbr", "brcn"), "BrC#N"),
    (("phosgene", "cocl2", "carbonyl chloride", "carbonyl dichloride"), "O=C(Cl)Cl"),
    (("triphosgene", "bis(trichloromethyl) carbonate", "btc"), "ClC(Cl)(Cl)OC(=O)OC(Cl)(Cl)Cl"),
    (("diphosgene", "trichloromethyl chloroformate"), "ClC(=O)OC(Cl)(Cl)Cl"),
    (("ammonia", "nh3", "ammonia gas"), "N"),
    (("ammonium hydroxide", "nh4oh", "aqueous ammonia", "ammonia water"), "[NH4+].[OH-]"),
    (("hydrazine", "n2h4"), "NN"),
    (("hydrazine hydrate", "hydrazine monohydrate", "n2h4.h2o"), "NN.O"),
    (("hydroxylamine", "nh2oh"), "NO"),
    (("hydroxylamine hydrochloride", "nh2oh.hcl", "hydroxylammonium chloride"), "Cl.NO"),
    (("hydrogen peroxide", "h2o2", "peroxide"), "OO"),
    (("hydrochloric acid", "hcl", "hydrogen chloride", "muriatic acid"), "Cl"),
    (("hydrobromic acid", "hbr", "hydrogen bromide"), "Br"),
    (("hydroiodic acid", "hi", "hydrogen iodide", "hydriodic acid"), "I"),
    (("hydrofluoric acid", "hf", "hydrogen fluoride"), "F"),
    (("sulfuric acid", "sulphuric acid", "h2so4", "oil of vitriol"), "OS(=O)(=O)O"),
    (("sulfurous acid", "sulphurous acid", "h2so3"), "OS(=O)O"),
    (("chlorosulfonic acid", "chlorosulphonic acid", "clso3h"), "OS(=O)(=O)Cl"),
    (("nitric acid", "hno3", "aqua fortis"), "O[N+](=O)[O-]"),
    (("nitrous acid", "hno2"), "ON=O"),
    (("phosphoric acid", "h3po4", "orthophosphoric acid"), "OP(=O)(O)O"),
    (("phosphorous acid", "h3po3"), "OP(O)=O"),
    (("hypophosphorous acid", "h3po2", "phosphinic acid"), "OP=O"),
    (("pyrophosphoric acid", "diphosphoric acid"), "OP(=O)(O)OP(=O)(O)O"),
    (("polyphosphoric acid", "ppa"), "OP(=O)(O)OP(=O)(O)OP(=O)(O)O"),
    (("perchloric acid", "hclo4"), "OCl(=O)(=O)=O"),
    (("periodic acid", "hio4", "h5io6"), "OI(=O)(=O)=O"),
    (("boric acid", "h3bo3"), "OB(O)O"),
    (("carbonic acid", "h2co3"), "OC(O)=O"),
    (("hypochlorous acid", "hocl"), "OCl"),
    (("hexafluorophosphoric acid", "hpf6"), "F[P-](F)(F)(F)(F)F.[H+]"),
    (("tetrafluoroboric acid", "hbf4", "fluoroboric acid"), "F[B-](F)(F)F.[H+]"),
    (("silica", "silicon dioxide", "sio2", "silica gel"), "O=[Si]=O"),
    (("alumina", "aluminium oxide", "aluminum oxide", "al2o3"), "[O-2].[O-2].[O-2].[Al+3].[Al+3]"),
    (("titanium dioxide", "titania", "tio2"), "O=[Ti]=O"),
    (("manganese dioxide", "mno2", "manganese(iv) oxide"), "O=[Mn]=O"),
    (("zinc oxide", "zno"), "[Zn]=O"),
    (("magnesium oxide", "mgo", "magnesia"), "[Mg]=O"),
    (("calcium oxide", "cao", "quicklime", "lime"), "[Ca]=O"),
    (("silver oxide", "ag2o", "silver(i) oxide"), "[Ag]O[Ag]"),
    (("copper(ii) oxide", "cuo", "cupric oxide"), "[Cu]=O"),
    (("iron(iii) oxide", "fe2o3", "ferric oxide", "rust"), "O=[Fe]O[Fe]=O"),
    (("phosphorus pentoxide", "p2o5", "p4o10", "phosphorus(v) oxide"), "O=P12OP3(=O)OP(=O)(O1)OP(=O)(O2)O3"),
    (("phosphorus pentachloride", "pcl5"), "ClP(Cl)(Cl)(Cl)Cl"),
    (("phosphorus trichloride", "pcl3"), "ClP(Cl)Cl"),
    (("phosphorus tribromide", "pbr3"), "BrP(Br)Br"),
    (("phosphorus oxychloride", "pocl3", "phosphoryl chloride"), "O=P(Cl)(Cl)Cl"),
    (("thionyl chloride", "socl2"), "O=S(Cl)Cl"),
    (("sulfuryl chloride", "sulphuryl chloride", "so2cl2"), "O=S(=O)(Cl)Cl"),
    (("sulfur tetrafluoride", "sf4"), "FS(F)(F)F"),
    (("sulfur hexafluoride", "sf6"), "FS(F)(F)(F)(F)F"),
    (("boron trifluoride", "bf3"), "FB(F)F"),
    (("boron tribromide", "bbr3"), "BrB(Br)Br"),
    (("boron trichloride", "bcl3"), "ClB(Cl)Cl"),
    (("borane", "bh3"), "B"),
    (("silane", "sih4"), "[SiH4]"),
    # ------------------------------------------------------------------ #
    # Inorganic bases and salts
    # ------------------------------------------------------------------ #
    (("sodium hydroxide", "naoh", "caustic soda", "lye"), "[Na+].[OH-]"),
    (("potassium hydroxide", "koh", "caustic potash"), "[K+].[OH-]"),
    (("lithium hydroxide", "lioh"), "[Li+].[OH-]"),
    (("cesium hydroxide", "caesium hydroxide", "csoh"), "[Cs+].[OH-]"),
    (("calcium hydroxide", "ca(oh)2", "slaked lime"), "[Ca+2].[OH-].[OH-]"),
    (("barium hydroxide", "ba(oh)2"), "[Ba+2].[OH-].[OH-]"),
    (("magnesium hydroxide", "mg(oh)2", "milk of magnesia"), "[Mg+2].[OH-].[OH-]"),
    (("tetramethylammonium hydroxide", "tmah", "tmaoh"), "C[N+](C)(C)C.[OH-]"),
    (("tetrabutylammonium hydroxide", "tbaoh"), "CCCC[N+](CCCC)(CCCC)CCCC.[OH-]"),
    (("benzyltrimethylammonium hydroxide", "triton b"), "C[N+](C)(C)Cc1ccccc1.[OH-]"),
    (("sodium chloride", "nacl", "table salt", "salt", "kochsalz"), "[Na+].[Cl-]"),
    (("potassium chloride", "kcl"), "[K+].[Cl-]"),
    (("lithium chloride", "licl"), "[Li+].[Cl-]"),
    (("calcium chloride", "cacl2"), "[Ca+2].[Cl-].[Cl-]"),
    (("magnesium chloride", "mgcl2"), "[Mg+2].[Cl-].[Cl-]"),
    (("ammonium chloride", "nh4cl", "sal ammoniac"), "[NH4+].[Cl-]"),
    (("sodium bromide", "nabr"), "[Na+].[Br-]"),
    (("potassium bromide", "kbr"), "[K+].[Br-]"),
    (("lithium bromide", "libr"), "[Li+].[Br-]"),
    (("sodium iodide", "nai"), "[Na+].[I-]"),
    (("potassium iodide", "ki"), "[K+].[I-]"),
    (("lithium iodide", "lii"), "[Li+].[I-]"),
    (("sodium fluoride", "naf"), "[Na+].[F-]"),
    (("potassium fluoride", "kf"), "[K+].[F-]"),
    (("cesium fluoride", "caesium fluoride", "csf"), "[Cs+].[F-]"),
    (("sodium bicarbonate", "sodium hydrogen carbonate", "nahco3", "baking soda", "bicarb"),
     "[Na+].OC(=O)[O-]"),
    (("potassium bicarbonate", "khco3", "potassium hydrogen carbonate"), "[K+].OC(=O)[O-]"),
    (("sodium carbonate", "na2co3", "soda ash", "washing soda"), "[Na+].[Na+].[O-]C(=O)[O-]"),
    (("potassium carbonate", "k2co3", "potash"), "[K+].[K+].[O-]C(=O)[O-]"),
    (("cesium carbonate", "caesium carbonate", "cs2co3"), "[Cs+].[Cs+].[O-]C(=O)[O-]"),
    (("lithium carbonate", "li2co3"), "[Li+].[Li+].[O-]C(=O)[O-]"),
    (("calcium carbonate", "caco3", "chalk", "limestone"), "[Ca+2].[O-]C(=O)[O-]"),
    (("sodium acetate", "naoac", "sodium ethanoate"), "[Na+].CC(=O)[O-]"),
    (("potassium acetate", "koac"), "[K+].CC(=O)[O-]"),
    (("ammonium acetate", "nh4oac"), "[NH4+].CC(=O)[O-]"),
    (("sodium formate", "hcoona", "nahco2"), "[Na+].[O-]C=O"),
    (("ammonium formate", "hcoonh4"), "[NH4+].[O-]C=O"),
    (("sodium benzoate",), "[Na+].[O-]C(=O)c1ccccc1"),
    (("sodium sulfate", "sodium sulphate", "na2so4", "glauber's salt"), "[Na+].[Na+].[O-]S(=O)(=O)[O-]"),
    (("magnesium sulfate", "magnesium sulphate", "mgso4", "epsom salt", "epsom salts"),
     "[Mg+2].[O-]S(=O)(=O)[O-]"),
    (("sodium bisulfate", "sodium hydrogen sulfate", "nahso4"), "[Na+].OS(=O)(=O)[O-]"),
    (("potassium bisulfate", "khso4", "potassium hydrogen sulfate"), "[K+].OS(=O)(=O)[O-]"),
    (("ammonium sulfate", "ammonium sulphate", "(nh4)2so4"), "[NH4+].[NH4+].[O-]S(=O)(=O)[O-]"),
    (("sodium sulfite", "sodium sulphite", "na2so3"), "[Na+].[Na+].[O-]S(=O)[O-]"),
    (("sodium bisulfite", "sodium hydrogen sulfite", "nahso3"), "[Na+].OS(=O)[O-]"),
    (("sodium metabisulfite", "na2s2o5", "sodium disulfite"), "[Na+].[Na+].[O-]S(=O)S(=O)(=O)[O-]"),
    (("sodium thiosulfate", "sodium thiosulphate", "na2s2o3", "hypo"), "[Na+].[Na+].[O-]S(=S)(=O)[O-]"),
    (("sodium dithionite", "na2s2o4", "sodium hydrosulfite"), "[Na+].[Na+].[O-]S(=O)S(=O)[O-]"),
    (("sodium sulfide", "na2s"), "[Na+].[Na+].[S-2]"),
    (("sodium hydrosulfide", "nash", "sodium bisulfide"), "[Na+].[SH-]"),
    (("sodium nitrite", "nano2"), "[Na+].[O-]N=O"),
    (("sodium nitrate", "nano3"), "[Na+].[O-][N+](=O)[O-]"),
    (("potassium nitrate", "kno3", "saltpetre", "saltpeter"), "[K+].[O-][N+](=O)[O-]"),
    (("silver nitrate", "agno3"), "[Ag+].[O-][N+](=O)[O-]"),
    (("ammonium nitrate", "nh4no3"), "[NH4+].[O-][N+](=O)[O-]"),
    (("sodium phosphate", "trisodium phosphate", "na3po4"), "[Na+].[Na+].[Na+].[O-]P(=O)([O-])[O-]"),
    (("disodium hydrogen phosphate", "na2hpo4", "sodium phosphate dibasic", "disodium phosphate"),
     "[Na+].[Na+].OP(=O)([O-])[O-]"),
    (("sodium dihydrogen phosphate", "nah2po4", "sodium phosphate monobasic", "monosodium phosphate"),
     "[Na+].OP(=O)(O)[O-]"),
    (("dipotassium hydrogen phosphate", "k2hpo4", "potassium phosphate dibasic"),
     "[K+].[K+].OP(=O)([O-])[O-]"),
    (("potassium dihydrogen phosphate", "kh2po4", "potassium phosphate monobasic"), "[K+].OP(=O)(O)[O-]"),
    (("potassium phosphate", "tripotassium phosphate", "k3po4"), "[K+].[K+].[K+].[O-]P(=O)([O-])[O-]"),
    (("sodium azide", "nan3"), "[Na+].[N-]=[N+]=[N-]"),
    (("sodium cyanide", "nacn"), "[Na+].[C-]#N"),
    (("potassium cyanide", "kcn"), "[K+].[C-]#N"),
    (("sodium thiocyanate", "nascn"), "[Na+].[S-]C#N"),
    (("potassium thiocyanate", "kscn"), "[K+].[S-]C#N"),
    (("sodium perchlorate", "naclo4"), "[Na+].[O-]Cl(=O)(=O)=O"),
    (("lithium perchlorate", "liclo4"), "[Li+].[O-]Cl(=O)(=O)=O"),
    (("sodium hypochlorite", "naocl", "bleach"), "[Na+].[O-]Cl"),
    (("sodium chlorite", "naclo2"), "[Na+].[O-]Cl=O"),
    (("sodium chlorate", "naclo3"), "[Na+].[O-]Cl(=O)=O"),
    (("sodium periodate", "sodium metaperiodate", "naio4"), "[Na+].[O-]I(=O)(=O)=O"),
    (("potassium permanganate", "kmno4"), "[K+].[O-][Mn](=O)(=O)=O"),
    (("potassium dichromate", "k2cr2o7"), "[K+].[K+].[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-]"),
    (("sodium dichromate", "na2cr2o7"), "[Na+].[Na+].[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-]"),
    (("chromium trioxide", "cro3", "chromic anhydride"), "O=[Cr](=O)=O"),
    (("potassium ferricyanide", "k3fe(cn)6", "potassium hexacyanoferrate(iii)"),
     "[K+].[K+].[K+].N#C[Fe-3](C#N)(C#N)(C#N)(C#N)C#N"),
    (("potassium ferrocyanide", "k4fe(cn)6", "potassium hexacyanoferrate(ii)"),
     "[K+].[K+].[K+].[K+].N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N"),
    (("sodium tetrafluoroborate", "nabf4"), "[Na+].F[B-](F)(F)F"),
    (("sodium hexafluorophosphate", "napf6"), "[Na+].F[P-](F)(F)(F)(F)F"),
    (("ammonium hexafluorophosphate", "nh4pf6"), "[NH4+].F[P-](F)(F)(F)(F)F"),
    (("potassium hexafluorophosphate", "kpf6"), "[K+].F[P-](F)(F)(F)(F)F"),
    (("lithium hexafluorophosphate", "lipf6"), "[Li+].F[P-](F)(F)(F)(F)F"),
    (("sodium tetraphenylborate", "nabph4"), "[Na+].c1ccc(cc1)[B-](c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("lithium triflate", "lithium trifluoromethanesulfonate", "liotf"), "[Li+].[O-]S(=O)(=O)C(F)(F)F"),
    (("sodium triflate", "naotf"), "[Na+].[O-]S(=O)(=O)C(F)(F)F"),
    (("potassium triflate", "kotf"), "[K+].[O-]S(=O)(=O)C(F)(F)F"),
    (("silver triflate", "agotf", "silver trifluoromethanesulfonate"), "[Ag+].[O-]S(=O)(=O)C(F)(F)F"),
    (("silver tetrafluoroborate", "agbf4"), "[Ag+].F[B-](F)(F)F"),
    (("silver hexafluoroantimonate", "agsbf6"), "[Ag+].F[Sb-](F)(F)(F)(F)F"),
    (("silver acetate", "agoac"), "[Ag+].CC(=O)[O-]"),
    (("silver carbonate", "ag2co3"), "[Ag+].[Ag+].[O-]C(=O)[O-]"),
    (("silver(i) fluoride", "agf", "silver fluoride"), "[Ag+].[F-]"),
    (("silver chloride", "agcl"), "[Ag+].[Cl-]"),
    (("calcium sulfate", "caso4", "gypsum", "drierite", "plaster of paris"), "[Ca+2].[O-]S(=O)(=O)[O-]"),
    (("barium sulfate", "baso4"), "[Ba+2].[O-]S(=O)(=O)[O-]"),
    (("copper(ii) sulfate", "copper sulfate", "cuso4", "cupric sulfate", "copper sulphate", "copper(ii) sulphate"),
     "[Cu+2].[O-]S(=O)(=O)[O-]"),
    (("iron(ii) sulfate", "ferrous sulfate", "feso4"), "[Fe+2].[O-]S(=O)(=O)[O-]"),
    (("zinc sulfate", "znso4"), "[Zn+2].[O-]S(=O)(=O)[O-]"),
    (("nickel(ii) chloride", "nicl2", "nickel chloride"), "[Ni+2].[Cl-].[Cl-]"),
    (("cobalt(ii) chloride", "cobalt chloride", "cobaltous chloride", "cobalt dichloride"),
     "[Co+2].[Cl-].[Cl-]"),
    (("tin(ii) chloride", "stannous chloride", "sncl2"), "Cl[Sn]Cl"),
    (("tin(iv) chloride", "stannic chloride", "sncl4", "tin tetrachloride"), "Cl[Sn](Cl)(Cl)Cl"),
    (("cerium(iii) chloride", "cecl3", "cerium chloride", "cerous chloride"), "[Ce+3].[Cl-].[Cl-].[Cl-]"),
    (("cerium ammonium nitrate", "can", "ceric ammonium nitrate", "(nh4)2ce(no3)6"),
     "[NH4+].[NH4+].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[Ce+4]"),
    (("lanthanum(iii) chloride", "lacl3"), "[La+3].[Cl-].[Cl-].[Cl-]"),
    (("sodium amide", "nanh2", "sodamide"), "[Na+].[NH2-]"),
    (("potassium amide", "knh2"), "[K+].[NH2-]"),
    (("sodium hydride", "nah"), "[Na+].[H-]"),
    (("potassium hydride", "kh"), "[K+].[H-]"),
    (("lithium hydride", "lih"), "[Li+].[H-]"),
    (("calcium hydride", "cah2"), "[Ca+2].[H-].[H-]"),
    (("lithium", "li metal", "lithium metal"), "[Li]"),
    (("sodium", "na metal", "sodium metal"), "[Na]"),
    (("potassium", "k metal", "potassium metal"), "[K]"),
    (("magnesium", "mg", "magnesium metal", "magnesium turnings"), "[Mg]"),
    (("zinc", "zn", "zinc dust", "zinc metal"), "[Zn]"),
    (("iron", "fe", "iron powder"), "[Fe]"),
    (("copper", "cu", "copper metal"), "[Cu]"),
    (("tin", "sn", "tin metal"), "[Sn]"),
    (("palladium", "pd", "pd/c", "pd-c", "palladium on carbon", "palladium on charcoal", "pd on c", "10% pd/c", "5% pd/c", "palladium black"),
     "[Pd]"),
    (("platinum", "pt", "pt/c", "platinum on carbon"), "[Pt]"),
    (("rhodium", "rh", "rh/c"), "[Rh]"),
    (("ruthenium", "ru"), "[Ru]"),
    (("nickel", "ni", "raney nickel", "raney ni", "ra-ni"), "[Ni]"),
    (("platinum oxide", "pto2", "adams' catalyst", "adams catalyst"), "O=[Pt]=O"),
    (("palladium hydroxide", "pd(oh)2", "pearlman's catalyst", "pearlmans catalyst", "pd(oh)2/c"),
     "[Pd+2].[OH-].[OH-]"),
    (("samarium", "sm"), "[Sm]"),
    (("indium", "indium metal"), "[In]"),
    (("selenium", "se"), "[Se]"),
    (("mercury", "hg", "quicksilver"), "[Hg]"),
    # ------------------------------------------------------------------ #
    # Organic acids
    # ------------------------------------------------------------------ #
    (("formic acid", "methanoic acid", "hcooh", "hco2h"), "OC=O"),
    (("acetic acid", "ethanoic acid", "acoh", "glacial acetic acid", "hoac", "ch3cooh", "essigsäure", "essigsaure"),
     "CC(=O)O"),
    (("propionic acid", "propanoic acid"), "CCC(=O)O"),
    (("butyric acid", "butanoic acid"), "CCCC(=O)O"),
    (("isobutyric acid", "2-methylpropanoic acid"), "CC(C)C(=O)O"),
    (("valeric acid", "pentanoic acid"), "CCCCC(=O)O"),
    (("pivalic acid", "trimethylacetic acid", "pivoh", "2,2-dimethylpropanoic acid"), "CC(C)(C)C(=O)O"),
    (("hexanoic acid", "caproic acid"), "CCCCCC(=O)O"),
    (("octanoic acid", "caprylic acid"), "CCCCCCCC(=O)O"),
    (("decanoic acid", "capric acid"), "CCCCCCCCCC(=O)O"),
    (("acrylic acid", "propenoic acid"), "C=CC(=O)O"),
    (("methacrylic acid", "2-methylpropenoic acid"), "CC(=C)C(=O)O"),
    (("crotonic acid", "(e)-2-butenoic acid"), "C/C=C/C(=O)O"),
    (("propiolic acid", "propynoic acid"), "C#CC(=O)O"),
    (("trifluoroacetic acid", "tfa", "cf3cooh", "cf3co2h"), "OC(=O)C(F)(F)F"),
    (("trichloroacetic acid", "tca", "ccl3cooh"), "OC(=O)C(Cl)(Cl)Cl"),
    (("dichloroacetic acid", "dca"), "OC(=O)C(Cl)Cl"),
    (("chloroacetic acid", "monochloroacetic acid"), "OC(=O)CCl"),
    (("bromoacetic acid",), "OC(=O)CBr"),
    (("iodoacetic acid", "iodoacetate"), "OC(=O)CI"),
    (("glycolic acid", "hydroxyacetic acid"), "OCC(=O)O"),
    (("lactic acid", "2-hydroxypropanoic acid", "dl-lactic acid"), "CC(O)C(=O)O"),
    (("l-lactic acid", "(s)-lactic acid"), "C[C@H](O)C(=O)O"),
    (("d-lactic acid", "(r)-lactic acid"), "C[C@@H](O)C(=O)O"),
    (("pyruvic acid", "2-oxopropanoic acid"), "CC(=O)C(=O)O"),
    (("glyoxylic acid",), "O=CC(=O)O"),
    (("levulinic acid", "4-oxopentanoic acid"), "CC(=O)CCC(=O)O"),
    (("oxalic acid", "ethanedioic acid"), "OC(=O)C(=O)O"),
    (("malonic acid", "propanedioic acid"), "OC(=O)CC(=O)O"),
    (("succinic acid", "butanedioic acid"), "OC(=O)CCC(=O)O"),
    (("glutaric acid", "pentanedioic acid"), "OC(=O)CCCC(=O)O"),
    (("adipic acid", "hexanedioic acid"), "OC(=O)CCCCC(=O)O"),
    (("pimelic acid", "heptanedioic acid"), "OC(=O)CCCCCC(=O)O"),
    (("suberic acid", "octanedioic acid"), "OC(=O)CCCCCCC(=O)O"),
    (("azelaic acid", "nonanedioic acid"), "OC(=O)CCCCCCCC(=O)O"),
    (("sebacic acid", "decanedioic acid"), "OC(=O)CCCCCCCCC(=O)O"),
    (("fumaric acid", "(e)-butenedioic acid", "trans-butenedioic acid"), "OC(=O)/C=C/C(=O)O"),
    (("maleic acid", "(z)-butenedioic acid", "cis-butenedioic acid"), "OC(=O)/C=C\\C(=O)O"),
    (("malic acid", "hydroxybutanedioic acid", "dl-malic acid"), "OC(CC(=O)O)C(=O)O"),
    (("l-malic acid", "(s)-malic acid"), "O[C@@H](CC(=O)O)C(=O)O"),
    (("tartaric acid", "l-tartaric acid", "l-(+)-tartaric acid", "(2r,3r)-tartaric acid", "(+)-tartaric acid"),
     "O=C(O)[C@H](O)[C@@H](O)C(=O)O"),
    (("d-tartaric acid", "d-(-)-tartaric acid", "(2s,3s)-tartaric acid", "(-)-tartaric acid"),
     "O=C(O)[C@@H](O)[C@H](O)C(=O)O"),
    (("meso-tartaric acid", "(2r,3s)-tartaric acid"), "O=C(O)[C@@H](O)[C@@H](O)C(=O)O"),
    (("citric acid", "2-hydroxypropane-1,2,3-tricarboxylic acid"), "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
    (("isocitric acid",), "OC(=O)C(O)C(CC(=O)O)C(=O)O"),
    (("alpha-ketoglutaric acid", "α-ketoglutaric acid", "2-oxoglutaric acid", "akg", "2-oxopentanedioic acid"),
     "OC(=O)CCC(=O)C(=O)O"),
    (("oxaloacetic acid", "oxaloacetate"), "OC(=O)CC(=O)C(=O)O"),
    (("aconitic acid", "cis-aconitic acid"), "OC(=O)C/C(=C\\C(=O)O)C(=O)O"),
    (("mevalonic acid",), "CC(O)(CCO)CC(=O)O"),
    (("shikimic acid",), "O[C@@H]1CC(=C[C@@H]([C@H]1O)O)C(=O)O"),
    (("quinic acid",), "O[C@@H]1C[C@](O)(C[C@@H](O)[C@H]1O)C(=O)O"),
    (("gluconic acid", "d-gluconic acid"), "OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)O"),
    (("glucuronic acid", "d-glucuronic acid"), "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C(=O)O"),
    (("benzoic acid", "phcooh", "phco2h", "bzoh"), "OC(=O)c1ccccc1"),
    (("salicylic acid", "2-hydroxybenzoic acid"), "OC(=O)c1ccccc1O"),
    (("4-hydroxybenzoic acid", "p-hydroxybenzoic acid"), "OC(=O)c1ccc(O)cc1"),
    (("anthranilic acid", "2-aminobenzoic acid"), "Nc1ccccc1C(=O)O"),
    (("4-aminobenzoic acid", "p-aminobenzoic acid", "paba"), "Nc1ccc(cc1)C(=O)O"),
    (("4-nitrobenzoic acid", "p-nitrobenzoic acid"), "OC(=O)c1ccc(cc1)[N+](=O)[O-]"),
    (("phthalic acid", "benzene-1,2-dicarboxylic acid"), "OC(=O)c1ccccc1C(=O)O"),
    (("isophthalic acid", "benzene-1,3-dicarboxylic acid"), "OC(=O)c1cccc(c1)C(=O)O"),
    (("terephthalic acid", "benzene-1,4-dicarboxylic acid", "tpa"), "OC(=O)c1ccc(cc1)C(=O)O"),
    (("phthalic anhydride",), "O=C1OC(=O)c2ccccc12"),
    (("maleic anhydride",), "O=C1OC(=O)C=C1"),
    (("succinic anhydride",), "O=C1CCC(=O)O1"),
    (("glutaric anhydride",), "O=C1CCCC(=O)O1"),
    (("trifluoroacetic anhydride", "tfaa"), "FC(F)(F)C(=O)OC(=O)C(F)(F)F"),
    (("acetic anhydride", "ac2o", "ethanoic anhydride"), "CC(=O)OC(C)=O"),
    (("propionic anhydride",), "CCC(=O)OC(=O)CC"),
    (("benzoic anhydride", "bz2o"), "O=C(OC(=O)c1ccccc1)c1ccccc1"),
    (("pivalic anhydride", "piv2o"), "CC(C)(C)C(=O)OC(=O)C(C)(C)C"),
    (("mandelic acid", "2-hydroxy-2-phenylacetic acid", "dl-mandelic acid"), "OC(C(=O)O)c1ccccc1"),
    (("(r)-mandelic acid", "d-mandelic acid", "(r)-(-)-mandelic acid"), "O[C@@H](C(=O)O)c1ccccc1"),
    (("(s)-mandelic acid", "l-mandelic acid", "(s)-(+)-mandelic acid"), "O[C@H](C(=O)O)c1ccccc1"),
    (("phenylacetic acid",), "OC(=O)Cc1ccccc1"),
    (("cinnamic acid", "trans-cinnamic acid", "(e)-cinnamic acid"), "OC(=O)/C=C/c1ccccc1"),
    (("gallic acid", "3,4,5-trihydroxybenzoic acid"), "OC(=O)c1cc(O)c(O)c(O)c1"),
    (("caffeic acid",), "OC(=O)/C=C/c1ccc(O)c(O)c1"),
    (("ferulic acid",), "COc1cc(/C=C/C(=O)O)ccc1O"),
    (("nicotinic acid", "niacin", "vitamin b3", "pyridine-3-carboxylic acid"), "OC(=O)c1cccnc1"),
    (("isonicotinic acid", "pyridine-4-carboxylic acid"), "OC(=O)c1ccncc1"),
    (("picolinic acid", "pyridine-2-carboxylic acid"), "OC(=O)c1ccccn1"),
    (("2-furoic acid", "furan-2-carboxylic acid"), "OC(=O)c1ccco1"),
    (("methanesulfonic acid", "msoh", "mesylic acid", "mesic acid", "ch3so3h"), "CS(=O)(=O)O"),
    (("p-toluenesulfonic acid", "tsoh", "ptsa", "tosic acid", "p-tsoh", "para-toluenesulfonic acid", "toluenesulfonic acid", "pts"),
     "Cc1ccc(cc1)S(=O)(=O)O"),
    (("camphorsulfonic acid", "csa", "(1s)-(+)-10-camphorsulfonic acid", "(+)-csa", "10-camphorsulfonic acid"),
     "CC1(C)[C@@H]2CC[C@@]1(CS(=O)(=O)O)C(=O)C2"),
    (("benzenesulfonic acid", "phso3h", "besylic acid"), "OS(=O)(=O)c1ccccc1"),
    (("trifluoromethanesulfonic acid", "triflic acid", "tfoh", "cf3so3h", "hotf"), "OS(=O)(=O)C(F)(F)F"),
    (("sulfamic acid",), "NS(=O)(=O)O"),
    (("taurine", "2-aminoethanesulfonic acid"), "NCCS(=O)(=O)O"),
    (("pyridinium p-toluenesulfonate", "ppts"), "Cc1ccc(cc1)S(=O)(=O)[O-].c1cc[nH+]cc1"),
    (("pyridinium hydrochloride", "pyridine hydrochloride", "py.hcl"), "[Cl-].c1cc[nH+]cc1"),
    (("diphenyl phosphate", "dpp", "diphenyl hydrogen phosphate"), "O=P(O)(Oc1ccccc1)Oc1ccccc1"),
    (("phenylphosphonic acid",), "O=P(O)(O)c1ccccc1"),
    (("phenylboronic acid", "phb(oh)2", "benzeneboronic acid"), "OB(O)c1ccccc1"),
    (("saccharin",), "O=C1NS(=O)(=O)c2ccccc12"),
    (("squaric acid", "3,4-dihydroxycyclobut-3-ene-1,2-dione"), "O=C1C(=O)C(O)=C1O"),
    (("meldrum's acid", "meldrums acid", "2,2-dimethyl-1,3-dioxane-4,6-dione"), "CC1(C)OC(=O)CC(=O)O1"),
    (("barbituric acid",), "O=C1CC(=O)NC(=O)N1"),
    (("picric acid", "2,4,6-trinitrophenol"), "Oc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]"),
    (("uric acid",), "O=c1[nH]c(=O)c2[nH]c(=O)[nH]c2[nH]1"),
    (("hippuric acid", "n-benzoylglycine"), "OC(=O)CNC(=O)c1ccccc1"),
    (("thioacetic acid", "acsh"), "CC(=O)S"),
    (("thioglycolic acid", "mercaptoacetic acid"), "OC(=O)CS"),
    (("lipoic acid", "alpha-lipoic acid", "α-lipoic acid", "thioctic acid"), "OC(=O)CCCCC1CCSS1"),
    # ------------------------------------------------------------------ #
    # Organic bases and amines
    # ------------------------------------------------------------------ #
    (("methylamine", "menh2", "ch3nh2"), "CN"),
    (("dimethylamine", "me2nh"), "CNC"),
    (("trimethylamine", "me3n", "nme3", "tma"), "CN(C)C"),
    (("ethylamine", "etnh2"), "CCN"),
    (("diethylamine", "et2nh", "hnet2"), "CCNCC"),
    (("triethylamine", "et3n", "tea", "net3", "nete3", "triethyl amine"), "CCN(CC)CC"),
    (("n-propylamine", "propylamine"), "CCCN"),
    (("isopropylamine", "iprnh2"), "CC(C)N"),
    (("diisopropylamine", "ipr2nh", "dipa"), "CC(C)NC(C)C"),
    (("n-butylamine", "butylamine", "bunh2"), "CCCCN"),
    (("tert-butylamine", "t-butylamine", "tbunh2"), "CC(C)(C)N"),
    (("dibutylamine", "di-n-butylamine", "bu2nh"), "CCCCNCCCC"),
    (("tributylamine", "bu3n", "tri-n-butylamine"), "CCCCN(CCCC)CCCC"),
    (("cyclohexylamine", "cynh2"), "NC1CCCCC1"),
    (("dicyclohexylamine", "cy2nh"), "C1CCC(CC1)NC1CCCCC1"),
    (("allylamine",), "C=CCN"),
    (("propargylamine",), "C#CCN"),
    (("benzylamine", "bnnh2"), "NCc1ccccc1"),
    (("dibenzylamine", "bn2nh"), "C(NCc1ccccc1)c1ccccc1"),
    (("n-methylbenzylamine", "bnnhme"), "CNCc1ccccc1"),
    (("phenethylamine", "2-phenylethylamine", "beta-phenylethylamine", "β-phenylethylamine"), "NCCc1ccccc1"),
    (("ethylenediamine", "en", "1,2-diaminoethane", "ethane-1,2-diamine"), "NCCN"),
    (("1,3-diaminopropane", "propane-1,3-diamine"), "NCCCN"),
    (("putrescine", "1,4-diaminobutane", "butane-1,4-diamine"), "NCCCCN"),
    (("cadaverine", "1,5-diaminopentane"), "NCCCCCN"),
    (("hexamethylenediamine", "1,6-diaminohexane", "hmda"), "NCCCCCCN"),
    (("spermidine",), "NCCCCNCCCN"),
    (("spermine",), "NCCCNCCCCNCCCN"),
    (("tmeda", "n,n,n',n'-tetramethylethylenediamine", "tetramethylethylenediamine"), "CN(C)CCN(C)C"),
    (("ethanolamine", "2-aminoethanol", "monoethanolamine", "mea"), "NCCO"),
    (("diethanolamine", "dea"), "OCCNCCO"),
    (("triethanolamine", "tea-oh", "trolamine"), "OCCN(CCO)CCO"),
    (("2-amino-2-methyl-1-propanol", "2-amino-2-methylpropan-1-ol"), "CC(C)(N)CO"),
    (("aniline", "phnh2", "aminobenzene", "phenylamine"), "Nc1ccccc1"),
    (("n-methylaniline", "phnhme"), "CNc1ccccc1"),
    (("n,n-dimethylaniline", "phnme2"), "CN(C)c1ccccc1"),
    (("n,n-diethylaniline", "phnet2"), "CCN(CC)c1ccccc1"),
    (("diphenylamine", "ph2nh"), "N(c1ccccc1)c1ccccc1"),
    (("triphenylamine", "ph3n"), "N(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("p-toluidine", "4-methylaniline", "4-toluidine"), "Cc1ccc(N)cc1"),
    (("o-toluidine", "2-methylaniline"), "Cc1ccccc1N"),
    (("p-anisidine", "4-methoxyaniline", "4-anisidine"), "COc1ccc(N)cc1"),
    (("o-anisidine", "2-methoxyaniline"), "COc1ccccc1N"),
    (("4-chloroaniline", "p-chloroaniline"), "Nc1ccc(Cl)cc1"),
    (("4-bromoaniline", "p-bromoaniline"), "Nc1ccc(Br)cc1"),
    (("4-fluoroaniline", "p-fluoroaniline"), "Nc1ccc(F)cc1"),
    (("4-nitroaniline", "p-nitroaniline", "pna"), "Nc1ccc(cc1)[N+](=O)[O-]"),
    (("2-nitroaniline", "o-nitroaniline"), "Nc1ccccc1[N+](=O)[O-]"),
    (("p-phenylenediamine", "1,4-diaminobenzene", "benzene-1,4-diamine", "ppd"), "Nc1ccc(N)cc1"),
    (("o-phenylenediamine", "1,2-diaminobenzene", "benzene-1,2-diamine", "opd"), "Nc1ccccc1N"),
    (("m-phenylenediamine", "1,3-diaminobenzene"), "Nc1cccc(N)c1"),
    (("4-aminophenol", "p-aminophenol"), "Nc1ccc(O)cc1"),
    (("2-aminopyridine",), "Nc1ccccn1"),
    (("3-aminopyridine",), "Nc1cccnc1"),
    (("4-aminopyridine", "4-ap", "fampridine"), "Nc1ccncc1"),
    (("2-aminothiazole",), "Nc1nccs1"),
    (("pyrrolidine",), "C1CCNC1"),
    (("n-methylpyrrolidine",), "CN1CCCC1"),
    (("piperidine",), "C1CCNCC1"),
    (("n-methylpiperidine", "1-methylpiperidine"), "CN1CCCCC1"),
    (("4-hydroxypiperidine", "piperidin-4-ol"), "OC1CCNCC1"),
    (("4-aminopiperidine",), "NC1CCNCC1"),
    (("piperazine",), "C1CNCCN1"),
    (("n-methylpiperazine", "1-methylpiperazine"), "CN1CCNCC1"),
    (("n-boc-piperazine", "1-boc-piperazine", "tert-butyl piperazine-1-carboxylate", "boc-piperazine"),
     "CC(C)(C)OC(=O)N1CCNCC1"),
    (("morpholine",), "C1COCCN1"),
    (("n-methylmorpholine", "nmm", "4-methylmorpholine"), "CN1CCOCC1"),
    (("thiomorpholine",), "C1CSCCN1"),
    (("azetidine",), "C1CNC1"),
    (("aziridine", "ethyleneimine"), "C1CN1"),
    (("azepane", "hexamethyleneimine"), "C1CCCNCC1"),
    (("dipea", "diea", "diisopropylethylamine", "n,n-diisopropylethylamine", "hunig's base", "hunigs base", "hünig's base", "hünigs base", "ethyldiisopropylamine", "ipr2net"),
     "CCN(C(C)C)C(C)C"),
    (("dbu", "1,8-diazabicyclo[5.4.0]undec-7-ene", "1,8-diazabicycloundec-7-ene"), "C1CCC2=NCCCN2CC1"),
    (("dbn", "1,5-diazabicyclo[4.3.0]non-5-ene"), "C1CN2CCCC2=NC1"),
    (("dabco", "1,4-diazabicyclo[2.2.2]octane", "triethylenediamine", "teda"), "C1CN2CCN1CC2"),
    (("dmap", "4-dimethylaminopyridine", "4-(dimethylamino)pyridine", "n,n-dimethylaminopyridine"),
     "CN(C)c1ccncc1"),
    (("ppy", "4-pyrrolidinopyridine", "4-(pyrrolidin-1-yl)pyridine"), "c1cc(ccn1)N1CCCC1"),
    (("proton sponge", "proton-sponge", "1,8-bis(dimethylamino)naphthalene", "dman"),
     "CN(C)c1cccc2cccc(N(C)C)c12"),
    (("tmg", "1,1,3,3-tetramethylguanidine", "tetramethylguanidine"), "CN(C)C(=N)N(C)C"),
    (("mtbd", "7-methyl-1,5,7-triazabicyclo[4.4.0]dec-5-ene"), "CN1CCCN2CCCN=C12"),
    (("tbd", "1,5,7-triazabicyclo[4.4.0]dec-5-ene"), "C1CN=C2NCCCN2C1"),
    (("bemp", "2-tert-butylimino-2-diethylamino-1,3-dimethylperhydro-1,3,2-diazaphosphorine"),
     "CCN(CC)P1(=NC(C)(C)C)N(C)CCCN1C"),
    (("bttp", "tert-butylimino-tri(pyrrolidino)phosphorane"), "CC(C)(C)N=P(N1CCCC1)(N1CCCC1)N1CCCC1"),
    (("quinuclidine", "1-azabicyclo[2.2.2]octane"), "C1CN2CCC1CC2"),
    (("3-quinuclidinol", "quinuclidin-3-ol"), "OC1CN2CCC1CC2"),
    (("2,6-di-tert-butylpyridine", "2,6-di-t-butylpyridine", "dtbp"), "CC(C)(C)c1cccc(n1)C(C)(C)C"),
    (("2,6-di-tert-butyl-4-methylpyridine", "dtbmp"), "Cc1cc(nc(c1)C(C)(C)C)C(C)(C)C"),
    (("n-methylimidazole", "1-methylimidazole", "nmi", "methylimidazole"), "Cn1ccnc1"),
    (("2-methylimidazole",), "Cc1ncc[nH]1"),
    (("benzimidazole",), "c1ccc2[nH]cnc2c1"),
    (("1,2,4-triazole", "1h-1,2,4-triazole"), "c1nc[nH]n1"),
    (("1,2,3-triazole", "1h-1,2,3-triazole"), "c1c[nH]nn1"),
    (("tetrazole", "1h-tetrazole"), "c1nnn[nH]1"),
    (("pyrazole",), "c1cn[nH]c1"),
    (("3,5-dimethylpyrazole",), "Cc1cc(C)[nH]n1"),
    (("oxazole",), "c1cocn1"),
    (("isoxazole",), "c1conc1"),
    (("thiazole",), "c1cscn1"),
    (("benzothiazole",), "c1ccc2scnc2c1"),
    (("benzoxazole",), "c1ccc2ocnc2c1"),
    (("pyrimidine",), "c1cncnc1"),
    (("pyrazine",), "c1cnccn1"),
    (("pyridazine",), "c1ccnnc1"),
    (("1,3,5-triazine", "s-triazine"), "c1ncncn1"),
    (("quinoline",), "c1ccc2ncccc2c1"),
    (("isoquinoline",), "c1ccc2cnccc2c1"),
    (("quinoxaline",), "c1ccc2nccnc2c1"),
    (("quinazoline",), "c1ccc2ncncc2c1"),
    (("purine",), "c1ncc2[nH]cnc2n1"),
    (("acridine",), "c1ccc2nc3ccccc3cc2c1"),
    (("phenanthroline", "1,10-phenanthroline", "phen", "o-phenanthroline", "ortho-phenanthroline"),
     "c1cnc2c(c1)ccc1cccnc12"),
    (("2,2'-bipyridine", "bipyridine", "bipy", "bpy", "2,2'-bipyridyl"), "c1ccc(nc1)-c1ccccn1"),
    (("4,4'-bipyridine", "4,4'-bipyridyl"), "c1cc(ccn1)-c1ccncc1"),
    (("terpyridine", "2,2':6',2''-terpyridine", "terpy", "tpy"), "c1ccc(nc1)-c1cccc(n1)-c1ccccn1"),
    (("guanidine",), "NC(N)=N"),
    (("guanidine hydrochloride", "guanidinium chloride", "guanidine hcl", "gdmcl", "gdncl", "gdnhcl", "guhcl", "gnhcl", "gdm.hcl", "gdn.hcl", "gu.hcl"),
     "NC(N)=[NH2+].[Cl-]"),
    (("guanidine thiocyanate", "guanidinium thiocyanate", "gitc", "gscn"), "NC(N)=[NH2+].[S-]C#N"),
    (("arginine methyl ester",), "COC(=O)[C@@H](N)CCCNC(N)=N"),
    (("urea", "carbamide", "harnstoff"), "NC(N)=O"),
    (("thiourea", "thiocarbamide"), "NC(N)=S"),
    (("n,n'-dimethylurea", "1,3-dimethylurea"), "CNC(=O)NC"),
    (("tetramethylurea", "tmu"), "CN(C)C(=O)N(C)C"),
    (("biuret",), "NC(=O)NC(N)=O"),
    (("semicarbazide",), "NNC(N)=O"),
    (("semicarbazide hydrochloride",), "Cl.NNC(N)=O"),
    (("thiosemicarbazide",), "NNC(N)=S"),
    (("cyanamide",), "NC#N"),
    (("dicyandiamide", "cyanoguanidine"), "NC(=N)NC#N"),
    (("melamine",), "Nc1nc(N)nc(N)n1"),
    (("acetamidine hydrochloride",), "Cl.CC(N)=N"),
    (("formamidine acetate",), "CC(=O)O.NC=N"),
    (("hydrazine hydrochloride",), "Cl.NN"),
    (("methylhydrazine",), "CNN"),
    (("phenylhydrazine",), "NNc1ccccc1"),
    (("phenylhydrazine hydrochloride",), "Cl.NNc1ccccc1"),
    (("2,4-dinitrophenylhydrazine", "dnph", "2,4-dnp"), "NNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]"),
    (("hydrazine monohydrobromide",), "Br.NN"),
    (("o-methylhydroxylamine hydrochloride", "methoxyamine hydrochloride", "meonh2.hcl", "nh2ome.hcl"),
     "Cl.CON"),
    (("n,o-dimethylhydroxylamine hydrochloride", "weinreb amine", "menh(ome).hcl", "n,o-dimethylhydroxylamine"),
     "Cl.CNOC"),
    (("o-benzylhydroxylamine", "bnonh2"), "NOCc1ccccc1"),
    (("o-benzylhydroxylamine hydrochloride", "bnonh2.hcl"), "Cl.NOCc1ccccc1"),
    # ------------------------------------------------------------------ #
    # Organometallic bases and nucleophiles
    # ------------------------------------------------------------------ #
    (("n-butyllithium", "butyllithium", "n-buli", "buli", "nbuli"), "CCCC[Li]"),
    (("sec-butyllithium", "sec-buli", "s-buli", "secbuli"), "CCC(C)[Li]"),
    (("tert-butyllithium", "t-buli", "tert-buli", "tbuli"), "CC(C)(C)[Li]"),
    (("methyllithium", "meli"), "C[Li]"),
    (("phenyllithium", "phli"), "[Li]c1ccccc1"),
    (("hexyllithium", "n-hexyllithium"), "CCCCCC[Li]"),
    (("lda", "lithium diisopropylamide"), "[Li+].CC(C)[N-]C(C)C"),
    (("lihmds", "lithium bis(trimethylsilyl)amide", "lithium hexamethyldisilazide", "lin(tms)2", "lithium hexamethyldisilazane"),
     "[Li+].C[Si](C)(C)[N-][Si](C)(C)C"),
    (("nahmds", "sodium bis(trimethylsilyl)amide", "sodium hexamethyldisilazide", "nan(tms)2"),
     "[Na+].C[Si](C)(C)[N-][Si](C)(C)C"),
    (("khmds", "potassium bis(trimethylsilyl)amide", "potassium hexamethyldisilazide", "kn(tms)2"),
     "[K+].C[Si](C)(C)[N-][Si](C)(C)C"),
    (("hmds", "hexamethyldisilazane", "bis(trimethylsilyl)amine", "1,1,1,3,3,3-hexamethyldisilazane"),
     "C[Si](C)(C)N[Si](C)(C)C"),
    (("litmp", "lithium tetramethylpiperidide", "lithium 2,2,6,6-tetramethylpiperidide"),
     "[Li+].CC1(C)CCCC(C)(C)[N-]1"),
    (("tmp", "2,2,6,6-tetramethylpiperidine", "tetramethylpiperidine"), "CC1(C)CCCC(C)(C)N1"),
    (("sodium methoxide", "naome", "sodium methylate", "meona"), "[Na+].C[O-]"),
    (("sodium ethoxide", "naoet", "sodium ethylate", "etona"), "[Na+].CC[O-]"),
    (("potassium tert-butoxide", "kotbu", "potassium t-butoxide", "tbuok", "kobut", "potassium tert-butylate"),
     "[K+].CC(C)(C)[O-]"),
    (("sodium tert-butoxide", "naotbu", "sodium t-butoxide", "tbuona"), "[Na+].CC(C)(C)[O-]"),
    (("lithium tert-butoxide", "liotbu", "tbuoli"), "[Li+].CC(C)(C)[O-]"),
    (("potassium tert-pentoxide", "potassium tert-amylate", "kotam", "potassium t-amyloxide"),
     "[K+].CCC(C)(C)[O-]"),
    (("sodium phenoxide", "naoph", "sodium phenolate"), "[Na+].[O-]c1ccccc1"),
    (("lithium acetylide ethylenediamine complex",), "[Li+].[C-]#C.NCCN"),
    (("methylmagnesium bromide", "memgbr", "methyl grignard"), "C[Mg]Br"),
    (("methylmagnesium chloride", "memgcl"), "C[Mg]Cl"),
    (("methylmagnesium iodide", "memgi"), "C[Mg]I"),
    (("ethylmagnesium bromide", "etmgbr"), "CC[Mg]Br"),
    (("isopropylmagnesium chloride", "iprmgcl"), "CC(C)[Mg]Cl"),
    (("isopropylmagnesium chloride lithium chloride complex", "turbo grignard", "iprmgcl.licl"),
     "CC(C)[Mg]Cl.[Li+].[Cl-]"),
    (("phenylmagnesium bromide", "phmgbr", "phenyl grignard"), "Br[Mg]c1ccccc1"),
    (("vinylmagnesium bromide", "vinyl grignard"), "C=C[Mg]Br"),
    (("allylmagnesium bromide", "allyl grignard"), "C=CC[Mg]Br"),
    (("benzylmagnesium chloride",), "Cl[Mg]Cc1ccccc1"),
    (("tert-butylmagnesium chloride", "tbumgcl"), "CC(C)(C)[Mg]Cl"),
    (("dimethylzinc", "me2zn", "znme2"), "C[Zn]C"),
    (("diethylzinc", "et2zn", "znet2"), "CC[Zn]CC"),
    (("trimethylaluminium", "trimethylaluminum", "alme3", "me3al"), "C[Al](C)C"),
    (("triethylaluminium", "triethylaluminum", "alet3", "et3al"), "CC[Al](CC)CC"),
    (("diethylaluminium chloride", "diethylaluminum chloride", "et2alcl"), "CC[Al](Cl)CC"),
    (("lithium dimethylcuprate", "me2culi", "gilman reagent"), "[Li+].C[Cu-]C"),
    (("trimethylsilylacetylene", "tms-acetylene", "ethynyltrimethylsilane", "tmsa"), "C#C[Si](C)(C)C"),
    (("trimethylsilyl cyanide", "tmscn", "cyanotrimethylsilane"), "C[Si](C)(C)C#N"),
    (("trimethylsilyl azide", "tmsn3", "azidotrimethylsilane"), "C[Si](C)(C)N=[N+]=[N-]"),
    (("trimethylsilyldiazomethane", "tmschn2", "tmsd"), "C[Si](C)(C)C=[N+]=[N-]"),
    (("tetramethylsilane", "tms"), "C[Si](C)(C)C"),
    (("triethylsilane", "et3sih", "hsiet3"), "CC[SiH](CC)CC"),
    (("triisopropylsilane", "tis", "tips-h", "ipr3sih"), "CC(C)[SiH](C(C)C)C(C)C"),
    (("phenylsilane", "phsih3"), "[SiH3]c1ccccc1"),
    (("polymethylhydrosiloxane", "pmhs"), "C[SiH](C)O[SiH](C)O[SiH](C)O[SiH](C)C"),
    (("tributyltin hydride", "bu3snh", "tbth"), "CCCC[SnH](CCCC)CCCC"),
    (("tributyltin chloride", "bu3sncl"), "CCCC[Sn](Cl)(CCCC)CCCC"),
    (("hexamethylditin", "me3snsnme3"), "C[Sn](C)(C)[Sn](C)(C)C"),
    (("hexabutylditin", "bu3snsnbu3"), "CCCC[Sn](CCCC)(CCCC)[Sn](CCCC)(CCCC)CCCC"),
    (("tributyl(vinyl)tin", "vinyltributyltin", "tributylvinylstannane"), "CCCC[Sn](CCCC)(CCCC)C=C"),
    (("allyltributyltin", "allyltributylstannane"), "CCCC[Sn](CCCC)(CCCC)CC=C"),
    (("tetrakis(trimethylsilyl)silane", "(tms)3sih", "tris(trimethylsilyl)silane", "ttmss"),
     "C[Si](C)(C)[SiH]([Si](C)(C)C)[Si](C)(C)C"),
    # ------------------------------------------------------------------ #
    # Oxidants
    # ------------------------------------------------------------------ #
    (("dess-martin periodinane", "dmp", "dess martin periodinane", "dess-martin reagent", "dess martin"),
     "CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc12"),
    (("ibx", "2-iodoxybenzoic acid", "o-iodoxybenzoic acid"), "O=C1OI(=O)(O)c2ccccc12"),
    (("pcc", "pyridinium chlorochromate"), "[O-][Cr](=O)(=O)Cl.c1cc[nH+]cc1"),
    (("pdc", "pyridinium dichromate"), "[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-].c1cc[nH+]cc1.c1cc[nH+]cc1"),
    (("jones reagent", "chromic acid", "h2cro4"), "O[Cr](=O)(=O)O"),
    (("tempo", "2,2,6,6-tetramethylpiperidin-1-oxyl", "2,2,6,6-tetramethylpiperidine-1-oxyl", "tempo radical"),
     "CC1(C)CCCC(C)(C)N1[O]"),
    (("4-acetamido-tempo", "acnh-tempo"), "CC(=O)NC1CC(C)(C)N([O])C(C)(C)C1"),
    (("4-hydroxy-tempo", "tempol", "4-oh-tempo"), "CC1(C)CC(O)CC(C)(C)N1[O]"),
    (("abno", "9-azabicyclo[3.3.1]nonane n-oxyl"), "[O]N1C2CCCC1CCC2"),
    (("mcpba", "m-cpba", "meta-chloroperoxybenzoic acid", "3-chloroperoxybenzoic acid", "3-chloroperbenzoic acid", "m-chloroperbenzoic acid", "meta-chloroperbenzoic acid"),
     "OOC(=O)c1cccc(Cl)c1"),
    (("peracetic acid", "peroxyacetic acid", "acooh", "paa"), "CC(=O)OO"),
    (("oxone", "potassium peroxymonosulfate", "potassium monopersulfate", "khso5"),
     "[K+].[O-]S(=O)(=O)OO.[K+].[K+].[O-]S(=O)(=O)[O-].[K+].[O-]S(=O)(=O)O"),
    (("potassium persulfate", "potassium peroxydisulfate", "k2s2o8"),
     "[K+].[K+].[O-]S(=O)(=O)OOS(=O)(=O)[O-]"),
    (("ammonium persulfate", "aps", "ammonium peroxydisulfate", "(nh4)2s2o8"),
     "[NH4+].[NH4+].[O-]S(=O)(=O)OOS(=O)(=O)[O-]"),
    (("sodium persulfate", "na2s2o8"), "[Na+].[Na+].[O-]S(=O)(=O)OOS(=O)(=O)[O-]"),
    (("sodium percarbonate",), "[Na+].[Na+].[O-]C(=O)[O-].OO"),
    (("sodium perborate",), "[Na+].[O-]B1OOB([O-])OO1.[Na+]"),
    (("urea hydrogen peroxide", "uhp", "urea peroxide", "carbamide peroxide"), "NC(N)=O.OO"),
    (("tbhp", "tert-butyl hydroperoxide", "t-butyl hydroperoxide", "tbuooh"), "CC(C)(C)OO"),
    (("cumene hydroperoxide", "chp"), "CC(C)(OO)c1ccccc1"),
    (("di-tert-butyl peroxide", "dtbp peroxide", "tbuootbu"), "CC(C)(C)OOC(C)(C)C"),
    (("benzoyl peroxide", "bpo", "dibenzoyl peroxide", "(bzo)2"), "O=C(OOC(=O)c1ccccc1)c1ccccc1"),
    (("osmium tetroxide", "oso4", "osmium(viii) oxide"), "O=[Os](=O)(=O)=O"),
    (("potassium osmate", "k2oso4", "potassium osmate dihydrate", "k2oso2(oh)4"),
     "[K+].[K+].[O-][Os](=O)(=O)[O-]"),
    (("ruthenium tetroxide", "ruo4"), "O=[Ru](=O)(=O)=O"),
    (("ruthenium(iii) chloride", "rucl3"), "[Ru+3].[Cl-].[Cl-].[Cl-]"),
    (("selenium dioxide", "seo2"), "O=[Se]=O"),
    (("lead tetraacetate", "pb(oac)4", "lta"), "CC(=O)O[Pb](OC(C)=O)(OC(C)=O)OC(C)=O"),
    (("iodobenzene diacetate", "phi(oac)2", "(diacetoxyiodo)benzene", "diacetoxyiodobenzene", "pida"),
     "CC(=O)OI(OC(C)=O)c1ccccc1"),
    (("pifa", "[bis(trifluoroacetoxy)iodo]benzene", "phi(ocf3)2", "phi(tfa)2", "bis(trifluoroacetoxy)iodobenzene"),
     "FC(F)(F)C(=O)OI(OC(=O)C(F)(F)F)c1ccccc1"),
    (("koser's reagent", "kosers reagent", "hydroxy(tosyloxy)iodobenzene", "hti"),
     "Cc1ccc(cc1)S(=O)(=O)OI(O)c1ccccc1"),
    (("iodosobenzene", "iodosylbenzene", "phio"), "O=Ic1ccccc1"),
    (("togni reagent", "togni reagent ii", "togni's reagent", "1-trifluoromethyl-1,2-benziodoxol-3(1h)-one"),
     "O=C1OI(C(F)(F)F)c2ccccc12"),
    (("selectfluor", "f-tedac", "1-chloromethyl-4-fluoro-1,4-diazoniabicyclo[2.2.2]octane bis(tetrafluoroborate)"),
     "F[N+]12CC[N+](CCl)(CC1)CC2.F[B-](F)(F)F.F[B-](F)(F)F"),
    (("nfsi", "n-fluorobenzenesulfonimide", "n-fluorobis(benzenesulfonyl)imide"),
     "FN(S(=O)(=O)c1ccccc1)S(=O)(=O)c1ccccc1"),
    (("n-bromosuccinimide", "nbs"), "O=C1CCC(=O)N1Br"),
    (("n-chlorosuccinimide", "ncs"), "O=C1CCC(=O)N1Cl"),
    (("n-iodosuccinimide", "nis"), "O=C1CCC(=O)N1I"),
    (("dbdmh", "1,3-dibromo-5,5-dimethylhydantoin"), "CC1(C)N(Br)C(=O)N(Br)C1=O"),
    (("dcdmh", "1,3-dichloro-5,5-dimethylhydantoin"), "CC1(C)N(Cl)C(=O)N(Cl)C1=O"),
    (("trichloroisocyanuric acid", "tcca", "tcica"), "O=c1n(Cl)c(=O)n(Cl)c(=O)n1Cl"),
    (("chloramine-t", "chloramine t"), "Cc1ccc(cc1)S(=O)(=O)[N-]Cl.[Na+]"),
    (("t-butyl hypochlorite", "tert-butyl hypochlorite", "tbuocl"), "CC(C)(C)OCl"),
    (("ddq", "2,3-dichloro-5,6-dicyano-1,4-benzoquinone"), "N#CC1=C(C#N)C(=O)C(Cl)=C(Cl)C1=O"),
    (("chloranil", "tetrachloro-1,4-benzoquinone", "p-chloranil"), "O=C1C(Cl)=C(Cl)C(=O)C(Cl)=C1Cl"),
    (("benzoquinone", "1,4-benzoquinone", "p-benzoquinone", "quinone", "bq"), "O=C1C=CC(=O)C=C1"),
    (("hydroquinone", "1,4-dihydroxybenzene", "hq"), "Oc1ccc(O)cc1"),
    (("catechol", "1,2-dihydroxybenzene", "pyrocatechol"), "Oc1ccccc1O"),
    (("resorcinol", "1,3-dihydroxybenzene"), "Oc1cccc(O)c1"),
    (("naphthoquinone", "1,4-naphthoquinone"), "O=C1C=CC(=O)c2ccccc12"),
    (("menadione", "vitamin k3", "2-methyl-1,4-naphthoquinone"), "CC1=CC(=O)c2ccccc2C1=O"),
    (("anthraquinone", "9,10-anthraquinone"), "O=C1c2ccccc2C(=O)c2ccccc12"),
    (("nitrosobenzene", "phno"), "O=Nc1ccccc1"),
    (("isoamyl nitrite", "isopentyl nitrite", "amyl nitrite"), "CC(C)CCON=O"),
    (("tert-butyl nitrite", "t-butyl nitrite", "tbuono"), "CC(C)(C)ON=O"),
    (("n-methylmorpholine n-oxide", "nmo", "nmmo", "4-methylmorpholine n-oxide"), "C[N+]1([O-])CCOCC1"),
    (("trimethylamine n-oxide", "tmao"), "C[N+](C)(C)[O-]"),
    (("pyridine n-oxide",), "[O-][n+]1ccccc1"),
    (("dimethyldioxirane", "dmdo"), "CC1(C)OO1"),
    (("nitrosonium tetrafluoroborate", "nobf4", "nitrosyl tetrafluoroborate"), "[N+]#[O+].F[B-](F)(F)F"),
    (("copper(ii) acetate", "cu(oac)2", "cupric acetate", "copper acetate"), "[Cu+2].CC(=O)[O-].CC(=O)[O-]"),
    (("ferric chloride", "iron(iii) chloride", "fecl3", "iron trichloride"), "Cl[Fe](Cl)Cl"),
    (("ferrous chloride", "iron(ii) chloride", "fecl2"), "[Fe+2].[Cl-].[Cl-]"),
    (("iron(iii) nitrate", "ferric nitrate", "fe(no3)3"),
     "[Fe+3].[O-][N+](=O)[O-].[O-][N+](=O)[O-].[O-][N+](=O)[O-]"),
    (("iron(iii) acetylacetonate", "fe(acac)3"), "CC(=O)/C=C(\\C)O[Fe](O/C(C)=C\\C(C)=O)O/C(C)=C\\C(C)=O"),
    (("manganese(iii) acetate", "mn(oac)3"), "[Mn+3].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-]"),
    (("mercury(ii) acetate", "hg(oac)2", "mercuric acetate"), "CC(=O)O[Hg]OC(C)=O"),
    (("thallium(iii) acetate",), "CC(=O)O[Tl](OC(C)=O)OC(C)=O"),
    (("sulfur trioxide pyridine complex", "so3.py", "so3-pyridine", "pyridine sulfur trioxide", "py.so3", "parikh-doering reagent"),
     "O=S(=O)([O-])[n+]1ccccc1"),
    # ------------------------------------------------------------------ #
    # Reductants
    # ------------------------------------------------------------------ #
    (("sodium borohydride", "nabh4", "sodium tetrahydroborate"), "[Na+].[BH4-]"),
    (("lithium borohydride", "libh4"), "[Li+].[BH4-]"),
    (("potassium borohydride", "kbh4"), "[K+].[BH4-]"),
    (("sodium cyanoborohydride", "nabh3cn", "sodium cyanotrihydroborate"), "[Na+].[BH3-]C#N"),
    (("sodium triacetoxyborohydride", "stab", "nabh(oac)3", "na(oac)3bh"),
     "[Na+].CC(=O)O[BH-](OC(C)=O)OC(C)=O"),
    (("lithium aluminium hydride", "lithium aluminum hydride", "lialh4", "lah"), "[Li+].[AlH4-]"),
    (("dibal-h", "dibal", "diisobutylaluminium hydride", "diisobutylaluminum hydride", "ibu2alh"),
     "CC(C)C[AlH]CC(C)C"),
    (("red-al", "sodium bis(2-methoxyethoxy)aluminium hydride", "sodium bis(2-methoxyethoxy)aluminum hydride", "vitride"),
     "[Na+].COCCO[AlH-](OCCOC)"),
    (("l-selectride", "lithium tri-sec-butylborohydride", "l selectride"), "[Li+].CCC(C)[BH-](C(C)CC)C(C)CC"),
    (("k-selectride", "potassium tri-sec-butylborohydride"), "[K+].CCC(C)[BH-](C(C)CC)C(C)CC"),
    (("super-hydride", "superhydride", "lithium triethylborohydride", "liet3bh"), "[Li+].CC[BH-](CC)CC"),
    (("borane thf complex", "bh3.thf", "borane-thf", "bh3-thf", "borane tetrahydrofuran complex"),
     "B.C1CCOC1"),
    (("borane dimethyl sulfide complex", "bh3.sme2", "bh3.dms", "borane-dms", "bms", "borane dimethylsulfide"),
     "B.CSC"),
    (("borane ammonia complex", "ammonia borane", "nh3.bh3", "bh3.nh3"), "[NH3+][BH3-]"),
    (("catecholborane", "hbcat"), "B1Oc2ccccc2O1"),
    (("pinacolborane", "hbpin"), "CC1(C)OB(OC1(C)C)"),
    (("9-bbn", "9-borabicyclo[3.3.1]nonane", "9-bbn dimer"), "B1C2CCCC1CCC2"),
    (("disiamylborane", "sia2bh"), "CC(C)C(C)BC(C)C(C)C"),
    (("dicyclohexylborane", "cy2bh"), "C1CCC(CC1)BC1CCCCC1"),
    (("thexylborane",), "CC(C)C(C)(C)B"),
    (("bis(pinacolato)diboron", "b2pin2", "bis(pinacolato)diborane"), "CC1(C)OB(OC1(C)C)B1OC(C)(C)C(C)(C)O1"),
    (("bis(catecholato)diboron", "b2cat2"), "B1(Oc2ccccc2O1)B1Oc2ccccc2O1"),
    (("bis(neopentyl glycolato)diboron", "b2neop2", "b2(neop)2"), "B1(OCC(C)(C)CO1)B1OCC(C)(C)CO1"),
    (("tetrahydroxydiboron", "b2(oh)4", "bba"), "OB(O)B(O)O"),
    (("samarium(ii) iodide", "samarium diiodide", "smi2", "kagan's reagent", "kagans reagent"), "I[Sm]I"),
    (("titanium(iii) chloride", "ticl3"), "Cl[Ti](Cl)Cl"),
    (("chromium(ii) chloride", "crcl2", "chromous chloride"), "[Cr+2].[Cl-].[Cl-]"),
    (("tin(ii) chloride dihydrate", "sncl2.2h2o"), "Cl[Sn]Cl.O.O"),
    (("diimide", "diazene", "n2h2"), "N=N"),
    (("tosylhydrazide", "p-toluenesulfonyl hydrazide", "tsnhnh2"), "Cc1ccc(cc1)S(=O)(=O)NN"),
    (("hantzsch ester", "diethyl 1,4-dihydro-2,6-dimethyl-3,5-pyridinedicarboxylate", "hantzsch dihydropyridine"),
     "CCOC(=O)C1=C(C)NC(C)=C(C1)C(=O)OCC"),
    (("lithium naphthalenide", "linaph"), "[Li+].c1ccc2cc[c-]cc2c1"),
    (("sodium naphthalenide",), "[Na+].c1ccc2cc[c-]cc2c1"),
    (("thiourea dioxide", "formamidinesulfinic acid"), "NC(=N)S(=O)O"),
    (("zinc amalgam", "zn(hg)", "zn/hg"), "[Zn].[Hg]"),
    (("sodium amalgam", "na(hg)", "na/hg"), "[Na].[Hg]"),
    (("wilkinson's catalyst", "wilkinsons catalyst", "wilkinson catalyst", "rhcl(pph3)3", "chlorotris(triphenylphosphine)rhodium(i)"),
     "Cl[Rh](<-P(c1ccccc1)(c1ccccc1)c1ccccc1)(<-P(c1ccccc1)(c1ccccc1)c1ccccc1)<-P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("crabtree's catalyst", "crabtrees catalyst", "crabtree catalyst", "[ir(cod)(pcy3)(py)]pf6"),
     "F[P-](F)(F)(F)(F)F.[Ir+].C1CC=CCCC=C1.C1CCC(CC1)P(C1CCCCC1)C1CCCCC1.c1ccncc1"),
    (("lindlar catalyst", "lindlar's catalyst", "pd/caco3/pb", "lindlar"), "[Pd].[Pb].[Ca+2].[O-]C(=O)[O-]"),
    (("zinc borohydride", "zn(bh4)2"), "[Zn+2].[BH4-].[BH4-]"),
    (("sodium hydrogen telluride", "nahte"), "[Na+].[TeH-]"),
    (("diethyl phosphite", "dep", "(eto)2p(o)h"), "CCOP(=O)OCC"),
    (("dimethyl phosphite", "(meo)2p(o)h"), "COP(=O)OC"),
    # ------------------------------------------------------------------ #
    # Amide/peptide coupling reagents and additives
    # ------------------------------------------------------------------ #
    (("dcc", "n,n'-dicyclohexylcarbodiimide", "dicyclohexylcarbodiimide"), "C1CCC(CC1)N=C=NC1CCCCC1"),
    (("dic", "n,n'-diisopropylcarbodiimide", "diisopropylcarbodiimide", "dipc"), "CC(C)N=C=NC(C)C"),
    (("edc", "edci", "edac", "1-ethyl-3-(3-dimethylaminopropyl)carbodiimide", "n-(3-dimethylaminopropyl)-n'-ethylcarbodiimide"),
     "CCN=C=NCCCN(C)C"),
    (("edc.hcl", "edci.hcl", "edc hydrochloride", "edci hydrochloride", "edc-hcl", "1-ethyl-3-(3-dimethylaminopropyl)carbodiimide hydrochloride"),
     "Cl.CCN=C=NCCCN(C)C"),
    (("dcu", "n,n'-dicyclohexylurea", "dicyclohexylurea"), "O=C(NC1CCCCC1)NC1CCCCC1"),
    (("hobt", "1-hydroxybenzotriazole", "hydroxybenzotriazole", "n-hydroxybenzotriazole"), "On1nnc2ccccc12"),
    (("hoat", "1-hydroxy-7-azabenzotriazole", "7-aza-1-hydroxybenzotriazole"), "On1nnc2cccnc12"),
    (("6-cl-hobt", "6-chloro-1-hydroxybenzotriazole"), "On1nnc2ccc(Cl)cc12"),
    (("oxyma", "oxyma pure", "ethyl cyanohydroxyiminoacetate", "ethyl 2-cyano-2-(hydroxyimino)acetate"),
     "CCOC(=O)C(=NO)C#N"),
    (("hosu", "n-hydroxysuccinimide", "nhs", "1-hydroxypyrrolidine-2,5-dione"), "O=C1CCC(=O)N1O"),
    (("sulfo-nhs", "n-hydroxysulfosuccinimide", "sulfo-nhs sodium salt", "n-hydroxysulfosuccinimide sodium salt"),
     "[Na+].O=C1CC(C(=O)N1O)S(=O)(=O)[O-]"),
    (("hopo", "2-hydroxypyridine n-oxide", "1-hydroxy-2-pyridone", "1-hydroxypyridin-2(1h)-one"),
     "On1ccccc1=O"),
    (("pentafluorophenol", "pfp-oh", "pfpoh", "c6f5oh"), "Oc1c(F)c(F)c(F)c(F)c1F"),
    (("hatu", "1-[bis(dimethylamino)methylene]-1h-1,2,3-triazolo[4,5-b]pyridinium 3-oxide hexafluorophosphate", "o-(7-azabenzotriazol-1-yl)-n,n,n',n'-tetramethyluronium hexafluorophosphate"),
     "CN(C)C(=[N+](C)C)On1nnc2cccnc12.F[P-](F)(F)(F)(F)F"),
    (("hbtu", "o-(benzotriazol-1-yl)-n,n,n',n'-tetramethyluronium hexafluorophosphate", "n,n,n',n'-tetramethyl-o-(1h-benzotriazol-1-yl)uronium hexafluorophosphate"),
     "CN(C)C(=[N+](C)C)On1nnc2ccccc12.F[P-](F)(F)(F)(F)F"),
    (("tbtu", "o-(benzotriazol-1-yl)-n,n,n',n'-tetramethyluronium tetrafluoroborate"),
     "CN(C)C(=[N+](C)C)On1nnc2ccccc12.F[B-](F)(F)F"),
    (("hctu", "o-(6-chlorobenzotriazol-1-yl)-n,n,n',n'-tetramethyluronium hexafluorophosphate"),
     "CN(C)C(=[N+](C)C)On1nnc2ccc(Cl)cc12.F[P-](F)(F)(F)(F)F"),
    (("tatu", "o-(7-azabenzotriazol-1-yl)-n,n,n',n'-tetramethyluronium tetrafluoroborate"),
     "CN(C)C(=[N+](C)C)On1nnc2cccnc12.F[B-](F)(F)F"),
    (("comu", "(1-cyano-2-ethoxy-2-oxoethylidenaminooxy)dimethylamino-morpholino-carbenium hexafluorophosphate"),
     "CCOC(=O)C(=NOC(=[N+](C)C)N1CCOCC1)C#N.F[P-](F)(F)(F)(F)F"),
    (("tstu", "n,n,n',n'-tetramethyl-o-(n-succinimidyl)uronium tetrafluoroborate"),
     "CN(C)C(=[N+](C)C)ON1C(=O)CCC1=O.F[B-](F)(F)F"),
    (("tptu", "o-(2-oxo-1(2h)pyridyl)-n,n,n',n'-tetramethyluronium tetrafluoroborate"),
     "CN(C)C(=[N+](C)C)On1ccccc1=O.F[B-](F)(F)F"),
    (("tffh", "fluoro-n,n,n',n'-tetramethylformamidinium hexafluorophosphate", "tetramethylfluoroformamidinium hexafluorophosphate"),
     "CN(C)C(F)=[N+](C)C.F[P-](F)(F)(F)(F)F"),
    (("pybop", "benzotriazol-1-yloxytripyrrolidinophosphonium hexafluorophosphate"),
     "C1CCN(C1)[P+](N1CCCC1)(N1CCCC1)On1nnc2ccccc12.F[P-](F)(F)(F)(F)F"),
    (("pyaop", "(7-azabenzotriazol-1-yloxy)tripyrrolidinophosphonium hexafluorophosphate"),
     "C1CCN(C1)[P+](N1CCCC1)(N1CCCC1)On1nnc2cccnc12.F[P-](F)(F)(F)(F)F"),
    (("pybrop", "bromotripyrrolidinophosphonium hexafluorophosphate"),
     "Br[P+](N1CCCC1)(N1CCCC1)N1CCCC1.F[P-](F)(F)(F)(F)F"),
    (("bop", "(benzotriazol-1-yloxy)tris(dimethylamino)phosphonium hexafluorophosphate", "castro's reagent"),
     "CN(C)[P+](N(C)C)(N(C)C)On1nnc2ccccc12.F[P-](F)(F)(F)(F)F"),
    (("bop-cl", "bis(2-oxo-3-oxazolidinyl)phosphinic chloride", "bopcl"), "O=C1OCCN1P(=O)(Cl)N1CCOC1=O"),
    (("t3p", "propylphosphonic anhydride", "propanephosphonic acid anhydride", "1-propanephosphonic acid cyclic anhydride", "2,4,6-tripropyl-1,3,5,2,4,6-trioxatriphosphorinane-2,4,6-trioxide"),
     "CCCP1(=O)OP(=O)(CCC)OP(=O)(CCC)O1"),
    (("cdi", "1,1'-carbonyldiimidazole", "carbonyldiimidazole", "n,n'-carbonyldiimidazole"),
     "O=C(n1ccnc1)n1ccnc1"),
    (("tcdi", "1,1'-thiocarbonyldiimidazole", "thiocarbonyldiimidazole"), "S=C(n1ccnc1)n1ccnc1"),
    (("dppa", "diphenylphosphoryl azide", "diphenyl phosphorazidate"),
     "[N-]=[N+]=NP(=O)(Oc1ccccc1)Oc1ccccc1"),
    (("diethyl cyanophosphonate", "diethyl phosphorocyanidate", "(eto)2p(o)cn"), "CCOP(=O)(OCC)C#N"),
    (("dmtmm", "4-(4,6-dimethoxy-1,3,5-triazin-2-yl)-4-methylmorpholinium chloride"),
     "COc1nc(OC)nc(n1)[N+]1(C)CCOCC1.[Cl-]"),
    (("cdmt", "2-chloro-4,6-dimethoxy-1,3,5-triazine"), "COc1nc(Cl)nc(OC)n1"),
    (("ibcf", "isobutyl chloroformate", "ibucocl"), "CC(C)COC(=O)Cl"),
    (("ethyl chloroformate", "etococl", "ethyl chlorocarbonate"), "CCOC(=O)Cl"),
    (("methyl chloroformate", "meococl"), "COC(=O)Cl"),
    (("pivaloyl chloride", "pivcl", "trimethylacetyl chloride"), "CC(C)(C)C(=O)Cl"),
    (("2,4,6-trichlorobenzoyl chloride", "yamaguchi reagent"), "O=C(Cl)c1c(Cl)cc(Cl)cc1Cl"),
    (("mnba", "2-methyl-6-nitrobenzoic anhydride", "shiina reagent"),
     "Cc1cccc(c1C(=O)OC(=O)c1c(C)cccc1[N+](=O)[O-])[N+](=O)[O-]"),
    (("ghosez's reagent", "ghosez reagent", "1-chloro-n,n,2-trimethyl-1-propenylamine", "1-chloro-n,n,2-trimethylpropenylamine"),
     "CC(C)=C(Cl)N(C)C"),
    (("mukaiyama reagent", "mukaiyama's reagent", "2-chloro-1-methylpyridinium iodide"),
     "C[n+]1ccccc1Cl.[I-]"),
    (("2-mercaptobenzothiazole", "mbt"), "Sc1nc2ccccc2s1"),
    (("dpds", "2,2'-dipyridyl disulfide", "2,2'-dithiodipyridine", "aldrithiol-2", "2-pyridyl disulfide"),
     "c1ccc(nc1)SSc1ccccn1"),
    # ------------------------------------------------------------------ #
    # Protecting-group reagents
    # ------------------------------------------------------------------ #
    (("di-tert-butyl dicarbonate", "boc anhydride", "boc2o", "(boc)2o", "boc-anhydride", "di-t-butyl dicarbonate"),
     "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"),
    (("boc-on", "2-(tert-butoxycarbonyloxyimino)-2-phenylacetonitrile"), "CC(C)(C)OC(=O)ON=C(C#N)c1ccccc1"),
    (("fmoc-cl", "fmoc chloride", "9-fluorenylmethyl chloroformate", "9-fluorenylmethoxycarbonyl chloride"),
     "O=C(Cl)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-osu", "fmoc-nhs", "n-(9-fluorenylmethoxycarbonyloxy)succinimide", "fmoc n-hydroxysuccinimide ester"),
     "O=C1CCC(=O)N1OC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("cbz-cl", "cbzcl", "benzyl chloroformate", "z-cl", "carbobenzoxy chloride", "benzyloxycarbonyl chloride"),
     "O=C(Cl)OCc1ccccc1"),
    (("cbz-osu", "z-osu", "n-(benzyloxycarbonyloxy)succinimide"), "O=C1CCC(=O)N1OC(=O)OCc1ccccc1"),
    (("alloc-cl", "allyl chloroformate"), "C=CCOC(=O)Cl"),
    (("troc-cl", "2,2,2-trichloroethyl chloroformate"), "ClC(Cl)(Cl)COC(=O)Cl"),
    (("teoc-osu", "2-(trimethylsilyl)ethyl succinimidyl carbonate"), "C[Si](C)(C)CCOC(=O)ON1C(=O)CCC1=O"),
    (("tbscl", "tbdmscl", "tert-butyldimethylsilyl chloride", "tert-butyldimethylchlorosilane", "tbs-cl", "tbdms-cl"),
     "CC(C)(C)[Si](C)(C)Cl"),
    (("tbsotf", "tbdmsotf", "tert-butyldimethylsilyl triflate", "tert-butyldimethylsilyl trifluoromethanesulfonate", "tbs-otf"),
     "CC(C)(C)[Si](C)(C)OS(=O)(=O)C(F)(F)F"),
    (("tipscl", "triisopropylsilyl chloride", "chlorotriisopropylsilane", "tips-cl"),
     "CC(C)[Si](Cl)(C(C)C)C(C)C"),
    (("tipsotf", "triisopropylsilyl triflate", "tips-otf"), "CC(C)[Si](OS(=O)(=O)C(F)(F)F)(C(C)C)C(C)C"),
    (("tescl", "triethylsilyl chloride", "chlorotriethylsilane", "tes-cl"), "CC[Si](Cl)(CC)CC"),
    (("tesotf", "triethylsilyl triflate", "tes-otf"), "CC[Si](CC)(CC)OS(=O)(=O)C(F)(F)F"),
    (("tbdpscl", "tert-butyldiphenylsilyl chloride", "tert-butylchlorodiphenylsilane", "tbdps-cl"),
     "CC(C)(C)[Si](Cl)(c1ccccc1)c1ccccc1"),
    (("tmscl", "chlorotrimethylsilane", "trimethylsilyl chloride", "trimethylchlorosilane", "tms-cl"),
     "C[Si](C)(C)Cl"),
    (("tmsotf", "trimethylsilyl triflate", "trimethylsilyl trifluoromethanesulfonate", "tms-otf"),
     "C[Si](C)(C)OS(=O)(=O)C(F)(F)F"),
    (("tmsi", "iodotrimethylsilane", "trimethylsilyl iodide", "tms-i"), "C[Si](C)(C)I"),
    (("tmsbr", "bromotrimethylsilane", "trimethylsilyl bromide"), "C[Si](C)(C)Br"),
    (("bstfa", "n,o-bis(trimethylsilyl)trifluoroacetamide"), "C[Si](C)(C)OC(=N[Si](C)(C)C)C(F)(F)F"),
    (("bsa", "n,o-bis(trimethylsilyl)acetamide"), "C[Si](C)(C)OC(C)=N[Si](C)(C)C"),
    (("mstfa", "n-methyl-n-(trimethylsilyl)trifluoroacetamide"), "CN([Si](C)(C)C)C(=O)C(F)(F)F"),
    (("pmbcl", "pmb-cl", "4-methoxybenzyl chloride", "p-methoxybenzyl chloride"), "COc1ccc(CCl)cc1"),
    (("pmbbr", "pmb-br", "4-methoxybenzyl bromide"), "COc1ccc(CBr)cc1"),
    (("pmb alcohol", "4-methoxybenzyl alcohol", "anisyl alcohol", "pmboh"), "COc1ccc(CO)cc1"),
    (("pmb trichloroacetimidate", "4-methoxybenzyl trichloroacetimidate"), "COc1ccc(COC(=N)C(Cl)(Cl)Cl)cc1"),
    (("benzyl trichloroacetimidate", "bn trichloroacetimidate"), "N=C(OCc1ccccc1)C(Cl)(Cl)Cl"),
    (("benzyl bromide", "bnbr", "benzylbromide", "(bromomethyl)benzene"), "BrCc1ccccc1"),
    (("benzyl chloride", "bncl", "(chloromethyl)benzene"), "ClCc1ccccc1"),
    (("benzyl alcohol", "bnoh", "phenylmethanol"), "OCc1ccccc1"),
    (("mom-cl", "momcl", "chloromethyl methyl ether", "methoxymethyl chloride"), "COCCl"),
    (("mem-cl", "memcl", "2-methoxyethoxymethyl chloride"), "COCCOCCl"),
    (("sem-cl", "semcl", "2-(trimethylsilyl)ethoxymethyl chloride"), "C[Si](C)(C)CCOCCl"),
    (("bom-cl", "bomcl", "benzyl chloromethyl ether"), "ClCOCc1ccccc1"),
    (("dmb-cl", "3,4-dimethoxybenzyl chloride"), "COc1ccc(CCl)cc1OC"),
    (("trityl chloride", "trcl", "tr-cl", "chlorotriphenylmethane", "triphenylmethyl chloride", "triphenylchloromethane"),
     "ClC(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("dmt-cl", "dmtcl", "4,4'-dimethoxytrityl chloride"), "COc1ccc(cc1)C(Cl)(c1ccccc1)c1ccc(OC)cc1"),
    (("mmt-cl", "4-methoxytrityl chloride"), "COc1ccc(cc1)C(Cl)(c1ccccc1)c1ccccc1"),
    (("acetyl chloride", "accl", "ethanoyl chloride"), "CC(=O)Cl"),
    (("benzoyl chloride", "bzcl", "phcocl"), "O=C(Cl)c1ccccc1"),
    (("tosyl chloride", "p-toluenesulfonyl chloride", "tscl", "ts-cl", "4-toluenesulfonyl chloride", "4-methylbenzenesulfonyl chloride"),
     "Cc1ccc(cc1)S(=O)(=O)Cl"),
    (("mesyl chloride", "methanesulfonyl chloride", "mscl", "ms-cl", "ch3so2cl"), "CS(=O)(=O)Cl"),
    (("mesyl anhydride", "methanesulfonic anhydride", "ms2o"), "CS(=O)(=O)OS(C)(=O)=O"),
    (("nosyl chloride", "nscl", "2-nitrobenzenesulfonyl chloride", "o-nosyl chloride", "ns-cl"),
     "O=S(=O)(Cl)c1ccccc1[N+](=O)[O-]"),
    (("4-nitrobenzenesulfonyl chloride", "p-nosyl chloride", "p-nscl"), "O=S(=O)(Cl)c1ccc(cc1)[N+](=O)[O-]"),
    (("benzenesulfonyl chloride", "phso2cl", "bscl"), "O=S(=O)(Cl)c1ccccc1"),
    (("triflic anhydride", "tf2o", "trifluoromethanesulfonic anhydride"),
     "FC(F)(F)S(=O)(=O)OS(=O)(=O)C(F)(F)F"),
    (("n-phenyl-bis(trifluoromethanesulfonimide)", "phntf2", "n-phenyltriflimide", "phn(tf)2", "tf2nph"),
     "FC(F)(F)S(=O)(=O)N(S(=O)(=O)C(F)(F)F)c1ccccc1"),
    (("comins' reagent", "comins reagent", "n-(5-chloro-2-pyridyl)triflimide"),
     "FC(F)(F)S(=O)(=O)N(S(=O)(=O)C(F)(F)F)c1ccc(Cl)cn1"),
    (("trifluoroacetyl chloride", "cf3cocl"), "FC(F)(F)C(=O)Cl"),
    (("ethyl trifluoroacetate", "cf3co2et", "tfa-oet"), "CCOC(=O)C(F)(F)F"),
    (("phthalimide",), "O=C1NC(=O)c2ccccc12"),
    (("potassium phthalimide", "phthnk"), "[K+].O=C1[N-]C(=O)c2ccccc12"),
    (("n-boc-glycine", "boc-gly-oh", "boc-glycine", "boc-gly"), "CC(C)(C)OC(=O)NCC(=O)O"),
    (("n-fmoc-glycine", "fmoc-gly-oh", "fmoc-glycine"), "OC(=O)CNC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-ala-oh", "fmoc-l-alanine", "fmoc-alanine"), "C[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("boc-ala-oh", "boc-l-alanine", "boc-alanine"), "C[C@H](NC(=O)OC(C)(C)C)C(=O)O"),
    (("fmoc-lys(boc)-oh", "fmoc-l-lys(boc)-oh"),
     "CC(C)(C)OC(=O)NCCCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("boc-lys(fmoc)-oh",), "CC(C)(C)OC(=O)N[C@@H](CCCCNC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-cys(trt)-oh",), "OC(=O)[C@H](CSC(c1ccccc1)(c1ccccc1)c1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-ser(tbu)-oh",), "CC(C)(C)OC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-asp(otbu)-oh",), "CC(C)(C)OC(=O)C[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-glu(otbu)-oh",), "CC(C)(C)OC(=O)CC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-arg(pbf)-oh",),
     "Cc1c(C)c2OCC(C)(C)c2c(C)c1S(=O)(=O)NC(=N)NCCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-his(trt)-oh",),
     "OC(=O)[C@H](Cc1cn(cn1)C(c1ccccc1)(c1ccccc1)c1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-asn(trt)-oh",), "OC(=O)[C@H](CC(=O)NC(c1ccccc1)(c1ccccc1)c1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-gln(trt)-oh",), "OC(=O)[C@H](CCC(=O)NC(c1ccccc1)(c1ccccc1)c1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-thr(tbu)-oh",), "C[C@@H](OC(C)(C)C)[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-tyr(tbu)-oh",), "CC(C)(C)Oc1ccc(C[C@H](NC(=O)OCC2c3ccccc3-c3ccccc23)C(=O)O)cc1"),
    (("fmoc-trp(boc)-oh",), "CC(C)(C)OC(=O)n1cc(C[C@H](NC(=O)OCC2c3ccccc3-c3ccccc23)C(=O)O)c2ccccc12"),
    (("fmoc-pro-oh", "fmoc-l-proline"), "OC(=O)[C@@H]1CCCN1C(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-phe-oh", "fmoc-l-phenylalanine"), "OC(=O)[C@H](Cc1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-leu-oh", "fmoc-l-leucine"), "CC(C)C[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-ile-oh", "fmoc-l-isoleucine"), "CC[C@H](C)[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-val-oh", "fmoc-l-valine"), "CC(C)[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-met-oh", "fmoc-l-methionine"), "CSCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-aib-oh",), "CC(C)(NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("dibenzofulvene", "dbf", "9-methylenefluorene"), "C=C1c2ccccc2-c2ccccc12"),
    (("fluorene", "9h-fluorene"), "C1c2ccccc2-c2ccccc12"),
    (("fluorenone", "9-fluorenone"), "O=C1c2ccccc2-c2ccccc12"),
    (("fmoc-oh", "9-fluorenylmethanol", "fluorenylmethanol"), "OCC1c2ccccc2-c2ccccc12"),
    (("2-chlorotrityl chloride", "2-cl-trt-cl", "2-chlorotriphenylmethyl chloride", "ctc"),
     "ClC(c1ccccc1Cl)(c1ccccc1)c1ccccc1"),
    (("rink amide linker", "rink linker", "4-[(2,4-dimethoxyphenyl)(aminomethyl)]phenoxyacetic acid"),
     "COc1ccc(C(N)c2ccc(OCC(=O)O)cc2)c(OC)c1"),
    (("4-hydroxymethylphenoxyacetic acid", "hmpa linker", "wang linker"), "OCc1ccc(OCC(=O)O)cc1"),
    (("4-methylpiperidine", "4-mp"), "CC1CCNCC1"),
    (("acetic acid n-hydroxysuccinimide ester", "nhs acetate", "acosu"), "CC(=O)ON1C(=O)CCC1=O"),
    (("2,2-dimethoxypropane", "acetone dimethyl acetal", "2,2-dmp"), "COC(C)(C)OC"),
    (("trimethyl orthoformate", "tmof", "hc(ome)3"), "COC(OC)OC"),
    (("triethyl orthoformate", "teof", "hc(oet)3"), "CCOC(OCC)OCC"),
    (("trimethyl orthoacetate", "mec(ome)3"), "COC(C)(OC)OC"),
    (("triethyl orthoacetate", "mec(oet)3"), "CCOC(C)(OCC)OCC"),
    (("1,3-propanedithiol",), "SCCCS"),
    (("1,2-ethanedithiol", "edt"), "SCCS"),
    (("1,3-dithiane",), "C1CSCSC1"),
    (("2-methyl-2-butene", "amylene"), "CC=C(C)C"),
    (("pinacol", "2,3-dimethyl-2,3-butanediol"), "CC(C)(O)C(C)(C)O"),
    (("neopentyl glycol", "2,2-dimethyl-1,3-propanediol"), "CC(C)(CO)CO"),
    (("mida", "n-methyliminodiacetic acid"), "CN(CC(=O)O)CC(=O)O"),
    (("phenylboronic acid mida ester", "ph-b(mida)"), "CN1CC(=O)OB(OC(=O)C1)c1ccccc1"),
    (("potassium phenyltrifluoroborate", "phbf3k"), "[K+].F[B-](F)(F)c1ccccc1"),
    (("allyl bromide", "3-bromopropene"), "C=CCBr"),
    (("allyl chloride", "3-chloropropene"), "C=CCCl"),
    (("propargyl bromide", "3-bromopropyne"), "C#CCBr"),
    (("1,3-dimethylbarbituric acid", "nndmba", "n,n'-dimethylbarbituric acid"), "CN1C(=O)CC(=O)N(C)C1=O"),
    (("thiophenol", "phsh", "benzenethiol"), "Sc1ccccc1"),
    (("4-methoxythiophenol", "4-methoxybenzenethiol"), "COc1ccc(S)cc1"),
    (("2-mercaptoethanol", "beta-mercaptoethanol", "β-mercaptoethanol", "bme", "b-me", "2-me", "βme", "mercaptoethanol", "2-hydroxyethanethiol"),
     "OCCS"),
    (("dodecanethiol", "1-dodecanethiol", "lauryl mercaptan"), "CCCCCCCCCCCCS"),
    # ------------------------------------------------------------------ #
    # Catalysts, ligands, Lewis acids and metal salts
    # ------------------------------------------------------------------ #
    (("pd(pph3)4", "tetrakis(triphenylphosphine)palladium(0)", "tetrakis", "tetrakis(triphenylphosphine)palladium", "pd(ph3p)4", "pd tetrakis"),
     "[Pd].c1ccc(cc1)P(c1ccccc1)c1ccccc1.c1ccc(cc1)P(c1ccccc1)c1ccccc1.c1ccc(cc1)P(c1ccccc1)c1ccccc1.c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    (("pd(oac)2", "palladium(ii) acetate", "palladium acetate", "palladium diacetate"), "CC(=O)O[Pd]OC(C)=O"),
    (("pdcl2", "palladium(ii) chloride", "palladium chloride", "palladium dichloride"), "Cl[Pd]Cl"),
    (("pd(tfa)2", "palladium(ii) trifluoroacetate"), "FC(F)(F)C(=O)O[Pd]OC(=O)C(F)(F)F"),
    (("pd2(dba)3", "tris(dibenzylideneacetone)dipalladium(0)", "tris(dibenzylideneacetone)dipalladium"),
     "[Pd].[Pd].O=C(/C=C/c1ccccc1)/C=C/c1ccccc1.O=C(/C=C/c1ccccc1)/C=C/c1ccccc1.O=C(/C=C/c1ccccc1)/C=C/c1ccccc1"),
    (("pd(dba)2", "bis(dibenzylideneacetone)palladium(0)"),
     "[Pd].O=C(/C=C/c1ccccc1)/C=C/c1ccccc1.O=C(/C=C/c1ccccc1)/C=C/c1ccccc1"),
    (("dba", "dibenzylideneacetone", "1,5-diphenyl-1,4-pentadien-3-one"), "O=C(/C=C/c1ccccc1)/C=C/c1ccccc1"),
    (("pdcl2(pph3)2", "bis(triphenylphosphine)palladium(ii) dichloride", "pd(pph3)2cl2", "bis(triphenylphosphine)palladium dichloride"),
     "Cl[Pd]Cl.c1ccc(cc1)P(c1ccccc1)c1ccccc1.c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    (("pdcl2(dppf)", "pd(dppf)cl2", "[1,1'-bis(diphenylphosphino)ferrocene]dichloropalladium(ii)", "pd(dppf)cl2.ch2cl2", "pd(dppf)cl2.dcm", "pdcl2(dppf).ch2cl2"),
     "Cl[Pd]Cl.[Fe+2].c1ccc(cc1)P(c1ccccc1)[c-]1cccc1.c1ccc(cc1)P(c1ccccc1)[c-]1cccc1"),
    (("pdcl2(mecn)2", "bis(acetonitrile)palladium(ii) dichloride", "pdcl2(ch3cn)2"), "Cl[Pd]Cl.CC#N.CC#N"),
    (("pdcl2(dppe)",), "Cl[Pd]Cl.c1ccc(cc1)P(CCP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("pd(amphos)cl2", "pd(amphos)2cl2", "bis(di-tert-butyl(4-dimethylaminophenyl)phosphine)dichloropalladium(ii)"),
     "Cl[Pd]Cl.CN(C)c1ccc(cc1)P(C(C)(C)C)C(C)(C)C.CN(C)c1ccc(cc1)P(C(C)(C)C)C(C)(C)C"),
    (("pd(acac)2", "palladium(ii) acetylacetonate"), "CC(=O)/C=C(\\C)O[Pd]O/C(C)=C\\C(C)=O"),
    (("pd(dtbpf)cl2", "[1,1'-bis(di-tert-butylphosphino)ferrocene]dichloropalladium(ii)"),
     "Cl[Pd]Cl.[Fe+2].CC(C)(C)P(C(C)(C)C)[c-]1cccc1.CC(C)(C)P(C(C)(C)C)[c-]1cccc1"),
    (("xphos", "2-dicyclohexylphosphino-2',4',6'-triisopropylbiphenyl"),
     "CC(C)c1cc(C(C)C)c(c(c1)C(C)C)-c1ccccc1P(C1CCCCC1)C1CCCCC1"),
    (("sphos", "2-dicyclohexylphosphino-2',6'-dimethoxybiphenyl"),
     "COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1"),
    (("ruphos", "2-dicyclohexylphosphino-2',6'-diisopropoxybiphenyl"),
     "CC(C)Oc1cccc(OC(C)C)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1"),
    (("brettphos", "2-(dicyclohexylphosphino)3,6-dimethoxy-2',4',6'-triisopropyl-1,1'-biphenyl"),
     "COc1ccc(OC)c(c1-c1c(cc(cc1C(C)C)C(C)C)C(C)C)P(C1CCCCC1)C1CCCCC1"),
    (("tbuxphos", "t-buxphos", "2-di-tert-butylphosphino-2',4',6'-triisopropylbiphenyl"),
     "CC(C)c1cc(C(C)C)c(c(c1)C(C)C)-c1ccccc1P(C(C)(C)C)C(C)(C)C"),
    (("davephos", "2-dicyclohexylphosphino-2'-(n,n-dimethylamino)biphenyl"),
     "CN(C)c1ccccc1-c1ccccc1P(C1CCCCC1)C1CCCCC1"),
    (("johnphos", "2-(di-tert-butylphosphino)biphenyl"), "CC(C)(C)P(C(C)(C)C)c1ccccc1-c1ccccc1"),
    (("cyjohnphos", "2-(dicyclohexylphosphino)biphenyl"), "C1CCC(CC1)P(C1CCCCC1)c1ccccc1-c1ccccc1"),
    (("xantphos", "4,5-bis(diphenylphosphino)-9,9-dimethylxanthene"),
     "CC1(C)c2cccc(P(c3ccccc3)c3ccccc3)c2Oc2c1cccc2P(c1ccccc1)c1ccccc1"),
    (("dppf", "1,1'-bis(diphenylphosphino)ferrocene"),
     "[Fe+2].c1ccc(cc1)P(c1ccccc1)[c-]1cccc1.c1ccc(cc1)P(c1ccccc1)[c-]1cccc1"),
    (("dppe", "1,2-bis(diphenylphosphino)ethane", "diphos"), "c1ccc(cc1)P(CCP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("dppp", "1,3-bis(diphenylphosphino)propane"), "c1ccc(cc1)P(CCCP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("dppb", "1,4-bis(diphenylphosphino)butane"), "c1ccc(cc1)P(CCCCP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("dppm", "bis(diphenylphosphino)methane"), "c1ccc(cc1)P(CP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("dpephos", "bis[(2-diphenylphosphino)phenyl] ether"),
     "c1ccc(cc1)P(c1ccccc1)c1ccccc1Oc1ccccc1P(c1ccccc1)c1ccccc1"),
    (("binap", "2,2'-bis(diphenylphosphino)-1,1'-binaphthyl", "rac-binap", "(r)-binap", "(s)-binap", "r-binap", "s-binap"),
     "c1ccc(cc1)P(c1ccccc1)c1ccc2ccccc2c1-c1c(ccc2ccccc12)P(c1ccccc1)c1ccccc1"),
    (("binol", "1,1'-bi-2-naphthol", "2,2'-dihydroxy-1,1'-binaphthyl", "rac-binol"),
     "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12"),
    (("triphenylphosphine", "pph3", "ph3p"), "c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    (("triphenylphosphine oxide", "pph3o", "ph3po", "tppo", "o=pph3"), "O=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("tricyclohexylphosphine", "pcy3", "cy3p"), "C1CCC(CC1)P(C1CCCCC1)C1CCCCC1"),
    (("tri-tert-butylphosphine", "p(tbu)3", "ptbu3", "tbu3p"), "CC(C)(C)P(C(C)(C)C)C(C)(C)C"),
    (("tri-tert-butylphosphonium tetrafluoroborate", "[htbu3p]bf4", "p(tbu)3.hbf4"),
     "CC(C)(C)[PH+](C(C)(C)C)C(C)(C)C.F[B-](F)(F)F"),
    (("tri-o-tolylphosphine", "p(o-tol)3", "tri(o-tolyl)phosphine"), "Cc1ccccc1P(c1ccccc1C)c1ccccc1C"),
    (("tri(2-furyl)phosphine", "tfp", "trifurylphosphine"), "c1coc(c1)P(c1ccco1)c1ccco1"),
    (("tributylphosphine", "pbu3", "bu3p", "tri-n-butylphosphine"), "CCCCP(CCCC)CCCC"),
    (("trimethylphosphine", "pme3", "me3p"), "CP(C)C"),
    (("triethylphosphine", "pet3", "et3p"), "CCP(CC)CC"),
    (("triphenylphosphite", "p(oph)3", "triphenyl phosphite"), "c1ccc(cc1)OP(Oc1ccccc1)Oc1ccccc1"),
    (("trimethyl phosphite", "p(ome)3"), "COP(OC)OC"),
    (("triethyl phosphite", "p(oet)3"), "CCOP(OCC)OCC"),
    (("triphenylarsine", "asph3"), "c1ccc(cc1)[As](c1ccccc1)c1ccccc1"),
    (("ipr", "1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene", "ipr carbene", "ipr nhc"),
     "CC(C)c1cccc(C(C)C)c1N1C=CN([C]1)c1c(cccc1C(C)C)C(C)C"),
    (("imes", "1,3-dimesitylimidazol-2-ylidene", "imes carbene", "imes nhc"),
     "Cc1cc(C)c(c(C)c1)N1C=CN([C]1)c1c(C)cc(C)cc1C"),
    (("ipr.hcl", "1,3-bis(2,6-diisopropylphenyl)imidazolium chloride", "ipr hcl"),
     "CC(C)c1cccc(C(C)C)c1-n1cc[n+](c1)-c1c(cccc1C(C)C)C(C)C.[Cl-]"),
    (("imes.hcl", "1,3-dimesitylimidazolium chloride", "imes hcl"),
     "Cc1cc(C)c(c(C)c1)-n1cc[n+](c1)-c1c(C)cc(C)cc1C.[Cl-]"),
    (("peppsi-ipr", "peppsi", "pd-peppsi-ipr", "[1,3-bis(2,6-diisopropylphenyl)imidazol-2-ylidene](3-chloropyridyl)palladium(ii) dichloride"),
     "CC(C)c1cccc(C(C)C)c1N1C=CN(c2c(cccc2C(C)C)C(C)C)[C]1->[Pd](Cl)(Cl)<-n1ccc(Cl)cc1"),
    (("grubbs i", "grubbs 1", "grubbs first generation", "grubbs' first generation catalyst", "grubbs first-generation catalyst", "grubbs catalyst 1st generation", "grubbs catalyst first generation", "g1", "grubbs-i"),
     "C1CCC(CC1)P(C1CCCCC1)(C1CCCCC1)->[Ru](Cl)(Cl)(=Cc1ccccc1)<-P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"),
    (("grubbs ii", "grubbs 2", "grubbs second generation", "grubbs' second generation catalyst", "grubbs second-generation catalyst", "grubbs catalyst 2nd generation", "grubbs catalyst second generation", "g2", "grubbs-ii"),
     "Cc1cc(C)c(c(C)c1)N1CCN(c2c(C)cc(C)cc2C)[C]1[Ru](Cl)(Cl)(=Cc1ccccc1)<-P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"),
    (("hoveyda-grubbs ii", "hoveyda-grubbs 2", "hoveyda grubbs ii", "hoveyda-grubbs second generation", "hoveyda-grubbs second-generation catalyst", "hg-ii", "hgii", "hg2", "hoveyda grubbs 2nd generation"),
     "Cc1cc(C)c(c(C)c1)N1CCN(c2c(C)cc(C)cc2C)[C]1[Ru](Cl)(Cl)=Cc1ccccc1OC(C)C"),
    (("hoveyda-grubbs i", "hoveyda-grubbs 1", "hoveyda grubbs i", "hoveyda-grubbs first generation", "hg-i", "hgi"),
     "CC(C)Oc1ccccc1C=[Ru](Cl)(Cl)<-P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"),
    (("schrock catalyst", "schrock's catalyst", "schrock molybdenum catalyst"),
     "CC(C)c1cccc(C(C)C)c1N=[Mo](=CC(C)(C)c1ccccc1)(OC(C)(C(F)(F)F)C(F)(F)F)OC(C)(C(F)(F)F)C(F)(F)F"),
    (("cui", "copper(i) iodide", "cuprous iodide", "copper iodide"), "[Cu]I"),
    (("cubr", "copper(i) bromide", "cuprous bromide", "copper bromide"), "[Cu]Br"),
    (("cucl", "copper(i) chloride", "cuprous chloride"), "[Cu]Cl"),
    (("cubr.sme2", "copper(i) bromide dimethyl sulfide complex", "cubr.dms", "cubr-sme2"), "[Cu]Br.CSC"),
    (("cucn", "copper(i) cyanide", "cuprous cyanide"), "[Cu]C#N"),
    (("copper(i) thiophene-2-carboxylate", "cutc"), "[Cu+].[O-]C(=O)c1cccs1"),
    (("cu(otf)2", "copper(ii) triflate", "copper(ii) trifluoromethanesulfonate"),
     "[Cu+2].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("cucl2", "copper(ii) chloride", "cupric chloride"), "Cl[Cu]Cl"),
    (("cubr2", "copper(ii) bromide", "cupric bromide"), "Br[Cu]Br"),
    (("cu(acac)2", "copper(ii) acetylacetonate"), "CC(=O)/C=C(\\C)O[Cu]O/C(C)=C\\C(C)=O"),
    (("tbaf", "tetrabutylammonium fluoride", "tetra-n-butylammonium fluoride", "bu4nf", "nbu4nf"),
     "CCCC[N+](CCCC)(CCCC)CCCC.[F-]"),
    (("tbab", "tetrabutylammonium bromide", "tetra-n-butylammonium bromide", "bu4nbr", "nbu4nbr"),
     "CCCC[N+](CCCC)(CCCC)CCCC.[Br-]"),
    (("tbai", "tetrabutylammonium iodide", "tetra-n-butylammonium iodide", "bu4ni", "nbu4ni"),
     "CCCC[N+](CCCC)(CCCC)CCCC.[I-]"),
    (("tbac", "tetrabutylammonium chloride", "bu4ncl", "nbu4ncl"), "CCCC[N+](CCCC)(CCCC)CCCC.[Cl-]"),
    (("tbahs", "tetrabutylammonium hydrogen sulfate", "tetrabutylammonium bisulfate", "bu4nhso4"),
     "CCCC[N+](CCCC)(CCCC)CCCC.OS(=O)(=O)[O-]"),
    (("tbapf6", "tetrabutylammonium hexafluorophosphate", "bu4npf6"),
     "CCCC[N+](CCCC)(CCCC)CCCC.F[P-](F)(F)(F)(F)F"),
    (("tbabf4", "tetrabutylammonium tetrafluoroborate", "bu4nbf4"), "CCCC[N+](CCCC)(CCCC)CCCC.F[B-](F)(F)F"),
    (("tbaoac", "tetrabutylammonium acetate"), "CCCC[N+](CCCC)(CCCC)CCCC.CC(=O)[O-]"),
    (("tbat", "tetrabutylammonium difluorotriphenylsilicate"),
     "CCCC[N+](CCCC)(CCCC)CCCC.F[Si-](F)(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("teba", "tebac", "benzyltriethylammonium chloride", "triethylbenzylammonium chloride"),
     "CC[N+](CC)(CC)Cc1ccccc1.[Cl-]"),
    (("aliquat 336", "methyltrioctylammonium chloride", "aliquat"),
     "CCCCCCCC[N+](C)(CCCCCCCC)CCCCCCCC.[Cl-]"),
    (("18-crown-6", "18-c-6"), "C1COCCOCCOCCOCCOCCO1"),
    (("15-crown-5", "15-c-5"), "C1COCCOCCOCCOCCO1"),
    (("12-crown-4", "12-c-4"), "C1COCCOCCOCCO1"),
    (("dibenzo-18-crown-6",), "c1ccc2OCCOCCOc3ccccc3OCCOCCOc2c1"),
    (("ti(oipr)4", "titanium(iv) isopropoxide", "titanium tetraisopropoxide", "titanium isopropoxide", "ti(o-ipr)4"),
     "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C"),
    (("ticl4", "titanium(iv) chloride", "titanium tetrachloride"), "Cl[Ti](Cl)(Cl)Cl"),
    (("bf3.oet2", "boron trifluoride etherate", "boron trifluoride diethyl etherate", "bf3-oet2", "bf3 etherate", "bf3.et2o"),
     "CC[O+](CC)[B-](F)(F)F"),
    (("alcl3", "aluminium chloride", "aluminum chloride", "aluminium trichloride", "aluminum trichloride"),
     "Cl[Al](Cl)Cl"),
    (("albr3", "aluminium bromide", "aluminum bromide"), "Br[Al](Br)Br"),
    (("zncl2", "zinc chloride", "zinc(ii) chloride"), "Cl[Zn]Cl"),
    (("znbr2", "zinc bromide"), "Br[Zn]Br"),
    (("zni2", "zinc iodide"), "I[Zn]I"),
    (("zn(otf)2", "zinc triflate", "zinc trifluoromethanesulfonate"),
     "[Zn+2].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("zn(oac)2", "zinc acetate"), "[Zn+2].CC(=O)[O-].CC(=O)[O-]"),
    (("mgbr2", "magnesium bromide"), "Br[Mg]Br"),
    (("mgbr2.oet2", "magnesium bromide etherate", "magnesium bromide diethyl etherate"), "Br[Mg]Br.CCOCC"),
    (("mg(otf)2", "magnesium triflate"), "[Mg+2].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("mgi2", "magnesium iodide"), "I[Mg]I"),
    (("sc(otf)3", "scandium(iii) triflate", "scandium triflate", "scandium trifluoromethanesulfonate"),
     "[Sc+3].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("yb(otf)3", "ytterbium(iii) triflate", "ytterbium triflate"),
     "[Yb+3].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("la(otf)3", "lanthanum(iii) triflate"),
     "[La+3].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("in(otf)3", "indium(iii) triflate"),
     "[In+3].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("bi(otf)3", "bismuth(iii) triflate"),
     "[Bi+3].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("incl3", "indium(iii) chloride"), "Cl[In](Cl)Cl"),
    (("bicl3", "bismuth(iii) chloride"), "Cl[Bi](Cl)Cl"),
    (("sbcl5", "antimony pentachloride"), "Cl[Sb](Cl)(Cl)(Cl)Cl"),
    (("gacl3", "gallium(iii) chloride"), "Cl[Ga](Cl)Cl"),
    (("et2alcn", "diethylaluminium cyanide", "diethylaluminum cyanide", "nagata's reagent"), "CC[Al](CC)C#N"),
    (("rh2(oac)4", "rhodium(ii) acetate dimer", "dirhodium tetraacetate", "rhodium acetate"),
     "[Rh+2].[Rh+2].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-]"),
    (("rh2(esp)2", "bis[rhodium(α,α,α',α'-tetramethyl-1,3-benzenedipropionic acid)]"),
     "[Rh+2].[Rh+2].CC(C)(CC(=O)[O-])c1cccc(c1)C(C)(C)CC(=O)[O-].CC(C)(CC(=O)[O-])c1cccc(c1)C(C)(C)CC(=O)[O-]"),
    (("[rh(cod)cl]2", "chloro(1,5-cyclooctadiene)rhodium(i) dimer", "rh(cod)cl dimer"),
     "Cl[Rh]Cl[Rh].C1CC=CCCC=C1.C1CC=CCCC=C1"),
    (("[rh(cod)2]bf4", "bis(1,5-cyclooctadiene)rhodium(i) tetrafluoroborate"),
     "[Rh+].C1CC=CCCC=C1.C1CC=CCCC=C1.F[B-](F)(F)F"),
    (("[ir(cod)cl]2", "chloro(1,5-cyclooctadiene)iridium(i) dimer"),
     "Cl[Ir]Cl[Ir].C1CC=CCCC=C1.C1CC=CCCC=C1"),
    (("[ir(cod)(ome)]2", "(1,5-cyclooctadiene)(methoxy)iridium(i) dimer"),
     "CO[Ir]O(C)[Ir].C1CC=CCCC=C1.C1CC=CCCC=C1"),
    (("ir(ppy)3", "fac-ir(ppy)3", "tris(2-phenylpyridine)iridium(iii)"),
     "[Ir+3].[c-]1ccccc1-c1ccccn1.[c-]1ccccc1-c1ccccn1.[c-]1ccccc1-c1ccccn1"),
    (("ir[df(cf3)ppy]2(dtbbpy)pf6", "ir(dfcf3ppy)2(dtbbpy)pf6", "ir[df(cf3)ppy]2(dtbpy)pf6", "ir photocatalyst"),
     "[Ir+3].F[P-](F)(F)(F)(F)F.CC(C)(C)c1ccnc(c1)-c1cc(ccn1)C(C)(C)C.[c-]1c(F)cc(F)cc1-c1ccc(cn1)C(F)(F)F.[c-]1c(F)cc(F)cc1-c1ccc(cn1)C(F)(F)F"),
    (("ru(bpy)3cl2", "tris(2,2'-bipyridyl)ruthenium(ii) chloride", "ru(bpy)3", "[ru(bpy)3]cl2"),
     "[Ru+2].[Cl-].[Cl-].c1ccc(nc1)-c1ccccn1.c1ccc(nc1)-c1ccccn1.c1ccc(nc1)-c1ccccn1"),
    (("4czipn", "2,4,5,6-tetra(9h-carbazol-9-yl)isophthalonitrile"),
     "N#Cc1c(-n2c3ccccc3c3ccccc23)c(-n2c3ccccc3c3ccccc23)c(-n2c3ccccc3c3ccccc23)c(C#N)c1-n1c2ccccc2c2ccccc12"),
    (("eosin y", "eosin", "eosin y disodium", "acid red 87", "tetrabromofluorescein disodium"),
     "[Na+].[Na+].[O-]C(=O)c1ccccc1-c1c2cc(Br)c(=O)c(Br)c-2oc2c(Br)c([O-])c(Br)cc12"),
    (("rucl2(pph3)3", "dichlorotris(triphenylphosphine)ruthenium(ii)"),
     "Cl[Ru](Cl)(<-P(c1ccccc1)(c1ccccc1)c1ccccc1)(<-P(c1ccccc1)(c1ccccc1)c1ccccc1)<-P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("[rucl2(p-cymene)]2", "dichloro(p-cymene)ruthenium(ii) dimer", "ru(p-cymene)cl2 dimer"),
     "Cl[Ru](Cl)Cl[Ru]Cl.CC(C)c1ccc(C)cc1.CC(C)c1ccc(C)cc1"),
    (("noyori catalyst", "rucl[(s,s)-tsdpen](p-cymene)", "(s,s)-ru-tsdpen"),
     "Cc1ccc(cc1)S(=O)(=O)N[C@@H](c1ccccc1)[C@@H](N)c1ccccc1.Cl[Ru].CC(C)c1ccc(C)cc1"),
    (("ni(cod)2", "bis(1,5-cyclooctadiene)nickel(0)"), "[Ni].C1CC=CCCC=C1.C1CC=CCCC=C1"),
    (("nicl2.glyme", "nicl2(dme)", "nickel(ii) chloride ethylene glycol dimethyl ether complex"),
     "Cl[Ni]Cl.COCCOC"),
    (("nicl2(dppp)", "[1,3-bis(diphenylphosphino)propane]nickel(ii) dichloride"),
     "Cl[Ni]Cl.c1ccc(cc1)P(CCCP(c1ccccc1)c1ccccc1)c1ccccc1"),
    (("nicl2(pph3)2",), "Cl[Ni]Cl.c1ccc(cc1)P(c1ccccc1)c1ccccc1.c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    (("ni(acac)2", "nickel(ii) acetylacetonate"), "CC(=O)/C=C(\\C)O[Ni]O/C(C)=C\\C(C)=O"),
    (("dtbbpy", "4,4'-di-tert-butyl-2,2'-bipyridine", "dtbpy"), "CC(C)(C)c1ccnc(c1)-c1cc(ccn1)C(C)(C)C"),
    (("dmbpy", "4,4'-dimethoxy-2,2'-bipyridine", "4,4'-dimethoxybipyridine"), "COc1ccnc(c1)-c1cc(OC)ccn1"),
    (("neocuproine", "2,9-dimethyl-1,10-phenanthroline"), "Cc1ccc2ccc3ccc(C)nc3c2n1"),
    (("bathophenanthroline", "4,7-diphenyl-1,10-phenanthroline"),
     "c1ccc(cc1)-c1ccnc2c1ccc1c(ccnc12)-c1ccccc1"),
    (("bathocuproine", "2,9-dimethyl-4,7-diphenyl-1,10-phenanthroline"),
     "Cc1cc(-c2ccccc2)c2ccc3c(-c4ccccc4)cc(C)nc3c2n1"),
    (("salen", "n,n'-bis(salicylidene)ethylenediamine"), "Oc1ccccc1/C=N/CC/N=C/c1ccccc1O"),
    (("jacobsen's catalyst", "jacobsen catalyst", "(r,r)-jacobsen catalyst", "(r,r)-(-)-n,n'-bis(3,5-di-tert-butylsalicylidene)-1,2-cyclohexanediaminomanganese(iii) chloride"),
     "CC(C)(C)c1cc(C=N[C@@H]2CCCC[C@H]2N=Cc2cc(cc(c2O)C(C)(C)C)C(C)(C)C)c(O)c(c1)C(C)(C)C.Cl[Mn]"),
    (("acetylacetone", "acac", "hacac", "2,4-pentanedione", "pentane-2,4-dione"), "CC(=O)CC(C)=O"),
    (("dibenzoylmethane", "dbm", "1,3-diphenyl-1,3-propanedione"), "O=C(CC(=O)c1ccccc1)c1ccccc1"),
    (("macmillan catalyst", "macmillan's catalyst", "(s)-5-benzyl-2,2,3-trimethylimidazolidin-4-one", "macmillan imidazolidinone"),
     "CN1C(=O)[C@H](Cc2ccccc2)NC1(C)C"),
    (("jørgensen-hayashi catalyst", "jorgensen-hayashi catalyst", "hayashi-jorgensen catalyst", "(s)-diphenylprolinol tms ether", "(s)-2-(diphenyl((trimethylsilyl)oxy)methyl)pyrrolidine"),
     "C[Si](C)(C)OC(c1ccccc1)(c1ccccc1)[C@@H]1CCCN1"),
    (("(s)-diphenylprolinol", "(s)-α,α-diphenyl-2-pyrrolidinemethanol", "(s)-alpha,alpha-diphenylprolinol"),
     "OC(c1ccccc1)(c1ccccc1)[C@@H]1CCCN1"),
    (("cbs catalyst", "(s)-cbs", "(s)-2-methyl-cbs-oxazaborolidine", "corey-bakshi-shibata catalyst"),
     "CB1OC(c2ccccc2)(c2ccccc2)[C@@H]2CCCN12"),
    (("(r)-cbs", "(r)-2-methyl-cbs-oxazaborolidine"), "CB1OC(c2ccccc2)(c2ccccc2)[C@H]2CCCN12"),
    (("(-)-sparteine", "sparteine"), "C1CCN2C[C@H]3C[C@H](C2C1)[C@@H]1CCCCN1C3"),
    (("ferrocene", "fecp2", "cp2fe", "bis(cyclopentadienyl)iron"), "[Fe+2].[cH-]1cccc1.[cH-]1cccc1"),
    (("ferrocenium hexafluorophosphate", "fcpf6"), "[Fe+3].[cH-]1cccc1.[cH-]1cccc1.F[P-](F)(F)(F)(F)F"),
    (("cyclopentadiene", "cp-h", "1,3-cyclopentadiene"), "C1C=CC=C1"),
    (("dicyclopentadiene", "dcpd"), "C1C=CC2C1C1C=CC2C1"),
    (("cp2ticl2", "titanocene dichloride", "bis(cyclopentadienyl)titanium dichloride"),
     "Cl[Ti+2]Cl.[cH-]1cccc1.[cH-]1cccc1"),
    (("cp2zrcl2", "zirconocene dichloride"), "Cl[Zr+2]Cl.[cH-]1cccc1.[cH-]1cccc1"),
    (("schwartz's reagent", "schwartz reagent", "cp2zr(h)cl", "zirconocene chloride hydride"),
     "Cl[ZrH+2].[cH-]1cccc1.[cH-]1cccc1"),
    (("tebbe reagent", "tebbe's reagent", "tebbes reagent"), "C[Al](C)Cl.C=[Ti+2].[cH-]1cccc1.[cH-]1cccc1"),
    (("petasis reagent", "dimethyltitanocene", "cp2time2"), "C[Ti+2]C.[cH-]1cccc1.[cH-]1cccc1"),
    (("mo(co)6", "molybdenum hexacarbonyl"),
     "[Mo].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+]"),
    (("co2(co)8", "dicobalt octacarbonyl"),
     "[Co].[Co].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+]"),
    (("fe(co)5", "iron pentacarbonyl"), "[Fe].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+]"),
    (("cr(co)6", "chromium hexacarbonyl"),
     "[Cr].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+]"),
    (("au(pph3)cl", "chloro(triphenylphosphine)gold(i)", "ph3paucl", "(ph3p)aucl"),
     "Cl[Au]<-P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("aucl3", "gold(iii) chloride"), "Cl[Au](Cl)Cl"),
    (("haucl4", "chloroauric acid", "tetrachloroauric acid", "gold(iii) chloride hydrate"),
     "Cl[Au-](Cl)(Cl)Cl.[H+]"),
    (("ptcl2", "platinum(ii) chloride"), "Cl[Pt]Cl"),
    (("k2ptcl4", "potassium tetrachloroplatinate(ii)"), "[K+].[K+].Cl[Pt-2](Cl)(Cl)Cl"),
    (("karstedt's catalyst", "karstedt catalyst"), "[Pt].C=C[Si](C)(C)O[Si](C)(C)C=C"),
    (("hgcl2", "mercury(ii) chloride", "mercuric chloride"), "Cl[Hg]Cl"),
    (("hg(otf)2", "mercury(ii) triflate"), "[Hg+2].[O-]S(=O)(=O)C(F)(F)F.[O-]S(=O)(=O)C(F)(F)F"),
    (("cocl2.6h2o", "cobalt(ii) chloride hexahydrate"), "[Co+2].[Cl-].[Cl-].O.O.O.O.O.O"),
    (("co(acac)2", "cobalt(ii) acetylacetonate"), "CC(=O)/C=C(\\C)O[Co]O/C(C)=C\\C(C)=O"),
    (("vo(acac)2", "vanadyl acetylacetonate", "vanadium(iv) oxide acetylacetonate"),
     "CC(=O)/C=C(\\C)O[V](=O)O/C(C)=C\\C(C)=O"),
    (("mn(dpm)3", "tris(dipivaloylmethanato)manganese(iii)", "manganese(iii) dipivaloylmethanate"),
     "CC(C)(C)C(=O)/C=C(\\O[Mn](O/C(=C\\C(=O)C(C)(C)C)C(C)(C)C)O/C(=C\\C(=O)C(C)(C)C)C(C)(C)C)C(C)(C)C"),
    (("fe(acac)2", "iron(ii) acetylacetonate"), "CC(=O)/C=C(\\C)O[Fe]O/C(C)=C\\C(C)=O"),
    (("zeolite", "molecular sieves", "4a molecular sieves", "3a molecular sieves", "ms 4a", "4å ms", "4a ms", "mol sieves"),
     "O=[Si]=O.[O-2].[O-2].[O-2].[Al+3].[Al+3]"),
    (("amberlyst-15", "amberlyst 15", "amberlyst", "dowex 50", "dowex 50w", "dowex"),
     "CC(c1ccc(cc1)S(=O)(=O)O)c1ccc(cc1)S(=O)(=O)O"),
    (("polystyrene", "ps"), "CC(c1ccccc1)Cc1ccccc1"),
    (("merrifield resin", "chloromethyl polystyrene"), "ClCc1ccc(cc1)C(C)Cc1ccccc1"),
    (("wang resin",), "OCc1ccc(OCc2ccc(cc2)C(C)Cc2ccccc2)cc1"),
    # ------------------------------------------------------------------ #
    # Phosphorus, sulfur and fluorination reagents
    # ------------------------------------------------------------------ #
    (("oxalyl chloride", "(cocl)2", "ethanedioyl dichloride", "oxalyl dichloride"), "O=C(Cl)C(=O)Cl"),
    (("oxalyl bromide", "(cobr)2"), "O=C(Br)C(=O)Br"),
    (("lawesson's reagent", "lawessons reagent", "lawesson reagent", "lr", "2,4-bis(4-methoxyphenyl)-1,3,2,4-dithiadiphosphetane 2,4-disulfide"),
     "COc1ccc(cc1)P1(=S)SP(=S)(S1)c1ccc(OC)cc1"),
    (("burgess reagent", "burgess' reagent", "methyl n-(triethylammoniumsulfonyl)carbamate", "(methoxycarbonylsulfamoyl)triethylammonium hydroxide inner salt"),
     "CC[N+](CC)(CC)S(=O)(=O)[N-]C(=O)OC"),
    (("martin sulfurane", "martin's sulfurane", "bis[α,α-bis(trifluoromethyl)benzenemethanolato]diphenylsulfur", "martin's reagent"),
     "FC(F)(F)C(OS(c1ccccc1)(c1ccccc1)OC(c1ccccc1)(C(F)(F)F)C(F)(F)F)(c1ccccc1)C(F)(F)F"),
    (("dead", "diethyl azodicarboxylate"), "CCOC(=O)/N=N/C(=O)OCC"),
    (("diad", "diisopropyl azodicarboxylate"), "CC(C)OC(=O)/N=N/C(=O)OC(C)C"),
    (("dtad", "di-tert-butyl azodicarboxylate", "dbad"), "CC(C)(C)OC(=O)/N=N/C(=O)OC(C)(C)C"),
    (("addp", "1,1'-(azodicarbonyl)dipiperidine"), "O=C(/N=N/C(=O)N1CCCCC1)N1CCCCC1"),
    (("tmad", "n,n,n',n'-tetramethylazodicarboxamide"), "CN(C)C(=O)/N=N/C(=O)N(C)C"),
    (("ptad", "4-phenyl-1,2,4-triazoline-3,5-dione"), "O=C1N=NC(=O)N1c1ccccc1"),
    (("tsunoda reagent", "cmbp", "(cyanomethylene)tributylphosphorane"), "CCCCP(=CC#N)(CCCC)CCCC"),
    (("deoxo-fluor", "deoxofluor", "bis(2-methoxyethyl)aminosulfur trifluoride"), "COCCN(CCOC)S(F)(F)F"),
    (("dast", "diethylaminosulfur trifluoride", "(diethylamino)sulfur trifluoride"), "CCN(CC)S(F)(F)F"),
    (("xtalfluor-e", "xtalfluor e"), "CC[N+](CC)=S(F)F.F[B-](F)(F)F"),
    (("fluolead", "4-tert-butyl-2,6-dimethylphenylsulfur trifluoride"), "Cc1cc(cc(C)c1S(F)(F)F)C(C)(C)C"),
    (("pybrox", "pyfluor", "2-pyridinesulfonyl fluoride"), "O=S(=O)(F)c1ccccn1"),
    (("tmscf3", "ruppert-prakash reagent", "ruppert's reagent", "trimethyl(trifluoromethyl)silane", "(trifluoromethyl)trimethylsilane", "cf3tms", "cf3sime3"),
     "C[Si](C)(C)C(F)(F)F"),
    (("sodium trifluoromethanesulfinate", "langlois reagent", "langlois' reagent", "cf3so2na", "naso2cf3"),
     "[Na+].[O-]S(=O)C(F)(F)F"),
    (("umemoto reagent", "umemoto's reagent", "s-(trifluoromethyl)dibenzothiophenium tetrafluoroborate"),
     "FC(F)(F)[s+]1c2ccccc2c2ccccc12.F[B-](F)(F)F"),
    (("trifluoromethyl iodide", "cf3i", "iodotrifluoromethane"), "FC(F)(F)I"),
    (("bromotrifluoromethane", "cf3br"), "FC(F)(F)Br"),
    (("diazomethane", "ch2n2"), "C=[N+]=[N-]"),
    (("ethyl diazoacetate", "eda", "n2chco2et"), "CCOC(=O)C=[N+]=[N-]"),
    (("diazald", "n-methyl-n-nitroso-p-toluenesulfonamide"), "Cc1ccc(cc1)S(=O)(=O)N(C)N=O"),
    (("p-absa", "p-acetamidobenzenesulfonyl azide", "4-acetamidobenzenesulfonyl azide"),
     "CC(=O)Nc1ccc(cc1)S(=O)(=O)N=[N+]=[N-]"),
    (("tosyl azide", "tsn3", "p-toluenesulfonyl azide"), "Cc1ccc(cc1)S(=O)(=O)N=[N+]=[N-]"),
    (("mesyl azide", "msn3", "methanesulfonyl azide"), "CS(=O)(=O)N=[N+]=[N-]"),
    (("imidazole-1-sulfonyl azide hydrochloride", "diazo transfer reagent", "stick's reagent"),
     "Cl.[N-]=[N+]=NS(=O)(=O)n1ccnc1"),
    (("trifluoromethanesulfonyl azide", "triflyl azide", "tfn3"), "FC(F)(F)S(=O)(=O)N=[N+]=[N-]"),
    # ------------------------------------------------------------------ #
    # Alkylating agents, alkyl halides and epoxides
    # ------------------------------------------------------------------ #
    (("methyl triflate", "meotf", "methyl trifluoromethanesulfonate"), "COS(=O)(=O)C(F)(F)F"),
    (("meerwein's salt", "meerwein salt", "trimethyloxonium tetrafluoroborate", "me3obf4", "me3o+bf4-"),
     "C[O+](C)C.F[B-](F)(F)F"),
    (("triethyloxonium tetrafluoroborate", "et3obf4", "meerwein's reagent"), "CC[O+](CC)CC.F[B-](F)(F)F"),
    (("dimethyl sulfate", "me2so4", "(meo)2so2"), "COS(=O)(=O)OC"),
    (("diethyl sulfate", "et2so4"), "CCOS(=O)(=O)OCC"),
    (("iodomethane", "methyl iodide", "mei", "ch3i"), "CI"),
    (("bromomethane", "methyl bromide", "mebr", "ch3br"), "CBr"),
    (("chloromethane", "methyl chloride", "mecl"), "CCl"),
    (("iodoethane", "ethyl iodide", "eti"), "CCI"),
    (("bromoethane", "ethyl bromide"), "CCBr"),
    (("1-bromopropane", "n-propyl bromide"), "CCCBr"),
    (("2-bromopropane", "isopropyl bromide", "iprbr"), "CC(C)Br"),
    (("1-bromobutane", "n-butyl bromide", "bubr"), "CCCCBr"),
    (("2-bromobutane", "sec-butyl bromide"), "CCC(C)Br"),
    (("tert-butyl bromide", "2-bromo-2-methylpropane", "tbubr"), "CC(C)(C)Br"),
    (("tert-butyl chloride", "2-chloro-2-methylpropane", "tbucl"), "CC(C)(C)Cl"),
    (("1-iodobutane", "n-butyl iodide", "bui"), "CCCCI"),
    (("1-bromohexane", "hexyl bromide"), "CCCCCCBr"),
    (("1-bromooctane", "octyl bromide"), "CCCCCCCCBr"),
    (("1-bromododecane", "dodecyl bromide", "lauryl bromide"), "CCCCCCCCCCCCBr"),
    (("cyclohexyl bromide", "bromocyclohexane"), "BrC1CCCCC1"),
    (("cyclopropyl bromide", "bromocyclopropane"), "BrC1CC1"),
    (("1,3-dibromopropane",), "BrCCCBr"),
    (("1,4-dibromobutane",), "BrCCCCBr"),
    (("1,5-dibromopentane",), "BrCCCCCBr"),
    (("2-bromoethanol",), "OCCBr"),
    (("2-chloroethanol",), "OCCCl"),
    (("3-bromopropanol", "3-bromo-1-propanol"), "OCCCBr"),
    (("epichlorohydrin", "2-(chloromethyl)oxirane"), "ClCC1CO1"),
    (("ethylene oxide", "oxirane"), "C1CO1"),
    (("propylene oxide", "2-methyloxirane"), "CC1CO1"),
    (("styrene oxide", "2-phenyloxirane"), "c1ccc(cc1)C1CO1"),
    (("cyclohexene oxide", "1,2-epoxycyclohexane"), "C1CCC2OC2C1"),
    (("glycidol", "oxiran-2-ylmethanol"), "OCC1CO1"),
    (("oxetane",), "C1COC1"),
    (("3-oxetanone", "oxetan-3-one"), "O=C1COC1"),
    (("3-oxetanol", "oxetan-3-ol"), "OC1COC1"),
    (("chloroacetyl chloride",), "ClCC(=O)Cl"),
    (("bromoacetyl bromide",), "BrCC(=O)Br"),
    (("chloroacetonitrile",), "ClCC#N"),
    (("bromoacetonitrile",), "BrCC#N"),
    (("ethyl bromoacetate", "brch2co2et"), "CCOC(=O)CBr"),
    (("methyl bromoacetate",), "COC(=O)CBr"),
    (("tert-butyl bromoacetate",), "CC(C)(C)OC(=O)CBr"),
    (("ethyl chloroacetate",), "CCOC(=O)CCl"),
    (("ethyl 2-bromopropionate", "ethyl 2-bromopropanoate"), "CCOC(=O)C(C)Br"),
    (("ethyl 2-bromoisobutyrate", "ethyl α-bromoisobutyrate", "ebib"), "CCOC(=O)C(C)(C)Br"),
    (("bromoacetaldehyde diethyl acetal",), "CCOC(CBr)OCC"),
    (("chloroacetaldehyde dimethyl acetal",), "COC(CCl)OC"),
    (("2-bromoacetophenone", "phenacyl bromide"), "BrCC(=O)c1ccccc1"),
    (("chloromethyl methyl sulfide",), "CSCCl"),
    (("pph3.br2", "dibromotriphenylphosphorane", "triphenylphosphine dibromide"),
     "Br[P](Br)(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("appel reagent", "pph3/ccl4"), "ClC(Cl)(Cl)Cl.c1ccc(cc1)P(c1ccccc1)c1ccccc1"),
    (("vilsmeier reagent", "chloromethylenedimethyliminium chloride", "(chloromethylene)dimethylammonium chloride"),
     "C[N+](C)=CCl.[Cl-]"),
    # ------------------------------------------------------------------ #
    # Alkenes, alkynes, monomers and radical initiators
    # ------------------------------------------------------------------ #
    (("ethyl vinyl ether", "eve", "vinyl ethyl ether"), "C=COCC"),
    (("butyl vinyl ether",), "CCCCOC=C"),
    (("methyl vinyl ketone", "mvk", "3-buten-2-one"), "C=CC(C)=O"),
    (("ethyl vinyl ketone", "evk", "1-penten-3-one"), "CCC(=O)C=C"),
    (("acrolein", "propenal", "2-propenal"), "C=CC=O"),
    (("methacrolein",), "C=C(C)C=O"),
    (("crotonaldehyde", "(e)-2-butenal"), "C/C=C/C=O"),
    (("acrylonitrile", "an"), "C=CC#N"),
    (("methyl acrylate",), "COC(=O)C=C"),
    (("ethyl acrylate",), "CCOC(=O)C=C"),
    (("tert-butyl acrylate", "t-butyl acrylate"), "CC(C)(C)OC(=O)C=C"),
    (("butyl acrylate", "n-butyl acrylate"), "CCCCOC(=O)C=C"),
    (("methyl methacrylate", "mma"), "COC(=O)C(C)=C"),
    (("2-hydroxyethyl methacrylate", "hema"), "CC(=C)C(=O)OCCO"),
    (("ethylene glycol dimethacrylate", "egdma"), "CC(=C)C(=O)OCCOC(=O)C(C)=C"),
    (("acrylamide",), "NC(=O)C=C"),
    (("n,n'-methylenebisacrylamide", "bisacrylamide", "bis-acrylamide", "mba"), "C=CC(=O)NCNC(=O)C=C"),
    (("n-isopropylacrylamide", "nipam", "nipaam"), "CC(C)NC(=O)C=C"),
    (("n,n-dimethylacrylamide", "dmaa"), "CN(C)C(=O)C=C"),
    (("vinyl acetate", "vam"), "CC(=O)OC=C"),
    (("vinyl chloride", "chloroethene"), "C=CCl"),
    (("vinylidene chloride", "1,1-dichloroethene"), "C=C(Cl)Cl"),
    (("vinylidene fluoride", "1,1-difluoroethene", "vdf"), "C=C(F)F"),
    (("tetrafluoroethylene", "tetrafluoroethene", "c2f4"), "FC(F)=C(F)F"),
    (("n-vinylpyrrolidone", "nvp", "1-vinyl-2-pyrrolidone"), "C=CN1CCCC1=O"),
    (("4-vinylpyridine",), "C=Cc1ccncc1"),
    (("divinylbenzene", "dvb", "1,4-divinylbenzene"), "C=Cc1ccc(C=C)cc1"),
    (("alpha-methylstyrene", "α-methylstyrene", "2-phenylpropene"), "CC(=C)c1ccccc1"),
    (("4-methylstyrene", "p-methylstyrene"), "C=Cc1ccc(C)cc1"),
    (("4-chlorostyrene",), "C=Cc1ccc(Cl)cc1"),
    (("4-methoxystyrene",), "C=Cc1ccc(OC)cc1"),
    (("isoprene", "2-methyl-1,3-butadiene"), "C=CC(C)=C"),
    (("1,3-butadiene", "butadiene"), "C=CC=C"),
    (("chloroprene", "2-chloro-1,3-butadiene"), "C=CC(Cl)=C"),
    (("propylene", "propene"), "CC=C"),
    (("1-butene", "but-1-ene"), "CCC=C"),
    (("isobutylene", "isobutene", "2-methylpropene"), "CC(C)=C"),
    (("1-hexene", "hex-1-ene"), "CCCCC=C"),
    (("1-octene", "oct-1-ene"), "CCCCCCC=C"),
    (("cyclohexene",), "C1CCC=CC1"),
    (("cyclopentene",), "C1CC=CC1"),
    (("cyclooctene", "cis-cyclooctene"), "C1CCC/C=C\\CC1"),
    (("1,5-cyclooctadiene", "cod", "cycloocta-1,5-diene"), "C1C/C=C\\CC/C=C\\1"),
    (("norbornene", "2-norbornene", "bicyclo[2.2.1]hept-2-ene"), "C1CC2C=CC1C2"),
    (("norbornadiene", "nbd"), "C1=CC2C=CC1C2"),
    (("1-hexyne", "hex-1-yne"), "CCCCC#C"),
    (("1-octyne", "oct-1-yne"), "CCCCCCC#C"),
    (("phenylacetylene", "ethynylbenzene"), "C#Cc1ccccc1"),
    (("propyne", "methylacetylene"), "CC#C"),
    (("2-butyne", "dimethylacetylene"), "CC#CC"),
    (("diphenylacetylene", "tolan"), "c1ccc(cc1)C#Cc1ccccc1"),
    (("dimethyl acetylenedicarboxylate", "dmad"), "COC(=O)C#CC(=O)OC"),
    (("ethyl propiolate",), "CCOC(=O)C#C"),
    (("methyl propiolate",), "COC(=O)C#C"),
    (("3-butyn-1-ol", "homopropargyl alcohol"), "C#CCCO"),
    (("4-pentyn-1-ol",), "C#CCCCO"),
    (("5-hexynoic acid",), "C#CCCCC(=O)O"),
    # ------------------------------------------------------------------ #
    # 1,3-Dicarbonyls, esters and nitriles
    # ------------------------------------------------------------------ #
    (("ethyl acetoacetate", "eaa", "ethyl 3-oxobutanoate"), "CCOC(=O)CC(C)=O"),
    (("methyl acetoacetate",), "COC(=O)CC(C)=O"),
    (("tert-butyl acetoacetate",), "CC(=O)CC(=O)OC(C)(C)C"),
    (("diethyl malonate", "dem", "malonic ester"), "CCOC(=O)CC(=O)OCC"),
    (("dimethyl malonate",), "COC(=O)CC(=O)OC"),
    (("di-tert-butyl malonate",), "CC(C)(C)OC(=O)CC(=O)OC(C)(C)C"),
    (("diethyl oxalate",), "CCOC(=O)C(=O)OCC"),
    (("dimethyl oxalate",), "COC(=O)C(=O)OC"),
    (("diethyl succinate",), "CCOC(=O)CCC(=O)OCC"),
    (("dimethyl succinate",), "COC(=O)CCC(=O)OC"),
    (("diethyl phthalate",), "CCOC(=O)c1ccccc1C(=O)OCC"),
    (("dimethyl terephthalate",), "COC(=O)c1ccc(cc1)C(=O)OC"),
    (("ethyl cyanoacetate",), "CCOC(=O)CC#N"),
    (("methyl cyanoacetate",), "COC(=O)CC#N"),
    (("malononitrile",), "N#CCC#N"),
    (("ethyl benzoate",), "CCOC(=O)c1ccccc1"),
    (("methyl benzoate",), "COC(=O)c1ccccc1"),
    (("methyl salicylate", "oil of wintergreen"), "COC(=O)c1ccccc1O"),
    (("ethyl 4-hydroxybenzoate", "ethylparaben"), "CCOC(=O)c1ccc(O)cc1"),
    (("methylparaben", "methyl 4-hydroxybenzoate"), "COC(=O)c1ccc(O)cc1"),
    (("propylparaben", "propyl 4-hydroxybenzoate"), "CCCOC(=O)c1ccc(O)cc1"),
    (("benzyl benzoate",), "O=C(OCc1ccccc1)c1ccccc1"),
    (("methyl 4-bromobenzoate",), "COC(=O)c1ccc(Br)cc1"),
    (("methyl 4-iodobenzoate",), "COC(=O)c1ccc(I)cc1"),
    (("4-bromobenzoic acid",), "OC(=O)c1ccc(Br)cc1"),
    (("4-iodobenzoic acid",), "OC(=O)c1ccc(I)cc1"),
    (("4-fluorobenzoic acid",), "OC(=O)c1ccc(F)cc1"),
    (("4-chlorobenzoic acid",), "OC(=O)c1ccc(Cl)cc1"),
    (("4-methoxybenzoic acid", "p-anisic acid", "anisic acid"), "COc1ccc(cc1)C(=O)O"),
    (("4-methylbenzoic acid", "p-toluic acid"), "Cc1ccc(cc1)C(=O)O"),
    (("4-(trifluoromethyl)benzoic acid",), "OC(=O)c1ccc(cc1)C(F)(F)F"),
    (("4-cyanobenzoic acid",), "OC(=O)c1ccc(cc1)C#N"),
    (("3-nitrobenzoic acid", "m-nitrobenzoic acid"), "OC(=O)c1cccc(c1)[N+](=O)[O-]"),
    (("2-nitrobenzoic acid", "o-nitrobenzoic acid"), "OC(=O)c1ccccc1[N+](=O)[O-]"),
    # ------------------------------------------------------------------ #
    # Aldehydes and ketones
    # ------------------------------------------------------------------ #
    (("formaldehyde", "methanal", "formalin", "hcho", "paraformaldehyde"), "C=O"),
    (("acetaldehyde", "ethanal", "ch3cho"), "CC=O"),
    (("propionaldehyde", "propanal"), "CCC=O"),
    (("butyraldehyde", "butanal"), "CCCC=O"),
    (("isobutyraldehyde", "2-methylpropanal"), "CC(C)C=O"),
    (("valeraldehyde", "pentanal"), "CCCCC=O"),
    (("hexanal",), "CCCCCC=O"),
    (("heptanal",), "CCCCCCC=O"),
    (("octanal",), "CCCCCCCC=O"),
    (("decanal",), "CCCCCCCCCC=O"),
    (("pivaldehyde", "pivalaldehyde", "trimethylacetaldehyde", "2,2-dimethylpropanal"), "CC(C)(C)C=O"),
    (("cyclohexanecarboxaldehyde", "cyclohexanecarbaldehyde"), "O=CC1CCCCC1"),
    (("glyoxal", "ethanedial"), "O=CC=O"),
    (("glutaraldehyde", "pentanedial", "1,5-pentanedial"), "O=CCCCC=O"),
    (("chloral", "trichloroacetaldehyde"), "ClC(Cl)(Cl)C=O"),
    (("chloral hydrate",), "OC(O)C(Cl)(Cl)Cl"),
    (("benzaldehyde", "phcho"), "O=Cc1ccccc1"),
    (("4-methoxybenzaldehyde", "anisaldehyde", "p-anisaldehyde", "4-anisaldehyde"), "COc1ccc(C=O)cc1"),
    (("4-nitrobenzaldehyde", "p-nitrobenzaldehyde"), "O=Cc1ccc(cc1)[N+](=O)[O-]"),
    (("4-chlorobenzaldehyde", "p-chlorobenzaldehyde"), "O=Cc1ccc(Cl)cc1"),
    (("4-bromobenzaldehyde", "p-bromobenzaldehyde"), "O=Cc1ccc(Br)cc1"),
    (("4-fluorobenzaldehyde", "p-fluorobenzaldehyde"), "O=Cc1ccc(F)cc1"),
    (("4-hydroxybenzaldehyde", "p-hydroxybenzaldehyde"), "O=Cc1ccc(O)cc1"),
    (("4-methylbenzaldehyde", "p-tolualdehyde", "tolualdehyde"), "Cc1ccc(C=O)cc1"),
    (("4-(dimethylamino)benzaldehyde", "p-dimethylaminobenzaldehyde", "ehrlich's reagent", "dmab"),
     "CN(C)c1ccc(C=O)cc1"),
    (("4-(trifluoromethyl)benzaldehyde",), "O=Cc1ccc(cc1)C(F)(F)F"),
    (("2-nitrobenzaldehyde", "o-nitrobenzaldehyde"), "O=Cc1ccccc1[N+](=O)[O-]"),
    (("2-hydroxybenzaldehyde", "salicylaldehyde"), "O=Cc1ccccc1O"),
    (("2-bromobenzaldehyde",), "O=Cc1ccccc1Br"),
    (("terephthalaldehyde", "1,4-benzenedicarboxaldehyde"), "O=Cc1ccc(C=O)cc1"),
    (("isophthalaldehyde", "1,3-benzenedicarboxaldehyde"), "O=Cc1cccc(C=O)c1"),
    (("piperonal", "heliotropin", "3,4-methylenedioxybenzaldehyde"), "O=Cc1ccc2OCOc2c1"),
    (("veratraldehyde", "3,4-dimethoxybenzaldehyde"), "COc1ccc(C=O)cc1OC"),
    (("furfural", "2-furaldehyde", "furan-2-carbaldehyde"), "O=Cc1ccco1"),
    (("5-hydroxymethylfurfural", "hmf", "5-hmf"), "OCc1ccc(C=O)o1"),
    (("2-thiophenecarboxaldehyde", "thiophene-2-carbaldehyde"), "O=Cc1cccs1"),
    (("2-pyridinecarboxaldehyde", "picolinaldehyde", "pyridine-2-carbaldehyde", "2pca"), "O=Cc1ccccn1"),
    (("3-pyridinecarboxaldehyde", "nicotinaldehyde"), "O=Cc1cccnc1"),
    (("4-pyridinecarboxaldehyde", "isonicotinaldehyde"), "O=Cc1ccncc1"),
    (("indole-3-carboxaldehyde", "indole-3-carbaldehyde"), "O=Cc1c[nH]c2ccccc12"),
    (("cinnamaldehyde", "trans-cinnamaldehyde", "(e)-cinnamaldehyde"), "O=C/C=C/c1ccccc1"),
    (("phenylacetaldehyde",), "O=CCc1ccccc1"),
    (("hydrocinnamaldehyde", "3-phenylpropanal"), "O=CCCc1ccccc1"),
    (("acetophenone", "phcome", "methyl phenyl ketone"), "CC(=O)c1ccccc1"),
    (("4-methoxyacetophenone",), "COc1ccc(cc1)C(C)=O"),
    (("4-nitroacetophenone",), "CC(=O)c1ccc(cc1)[N+](=O)[O-]"),
    (("4-bromoacetophenone",), "CC(=O)c1ccc(Br)cc1"),
    (("4-chloroacetophenone",), "CC(=O)c1ccc(Cl)cc1"),
    (("4-fluoroacetophenone",), "CC(=O)c1ccc(F)cc1"),
    (("4-hydroxyacetophenone",), "CC(=O)c1ccc(O)cc1"),
    (("2-acetylpyridine",), "CC(=O)c1ccccn1"),
    (("3-acetylpyridine",), "CC(=O)c1cccnc1"),
    (("2-acetylfuran",), "CC(=O)c1ccco1"),
    (("2-acetylthiophene",), "CC(=O)c1cccs1"),
    (("propiophenone",), "CCC(=O)c1ccccc1"),
    (("benzophenone", "ph2co", "diphenyl ketone"), "O=C(c1ccccc1)c1ccccc1"),
    (("4,4'-dimethoxybenzophenone",), "COc1ccc(cc1)C(=O)c1ccc(OC)cc1"),
    (("michler's ketone", "michlers ketone", "4,4'-bis(dimethylamino)benzophenone"),
     "CN(C)c1ccc(cc1)C(=O)c1ccc(cc1)N(C)C"),
    (("benzil",), "O=C(C(=O)c1ccccc1)c1ccccc1"),
    (("benzoin",), "OC(C(=O)c1ccccc1)c1ccccc1"),
    (("chalcone", "trans-chalcone", "(e)-chalcone"), "O=C(/C=C/c1ccccc1)c1ccccc1"),
    (("cyclopropanone ethyl hemiketal",), "CCOC1(O)CC1"),
    (("cyclobutanone",), "O=C1CCC1"),
    (("cyclopentanone",), "O=C1CCCC1"),
    (("cyclohexanone",), "O=C1CCCCC1"),
    (("2-methylcyclohexanone",), "CC1CCCCC1=O"),
    (("4-methylcyclohexanone",), "CC1CCC(=O)CC1"),
    (("cycloheptanone",), "O=C1CCCCCC1"),
    (("cyclooctanone",), "O=C1CCCCCCC1"),
    (("2-cyclohexen-1-one", "cyclohexenone", "cyclohex-2-en-1-one"), "O=C1C=CCCC1"),
    (("2-cyclopenten-1-one", "cyclopentenone"), "O=C1C=CCC1"),
    (("isophorone", "3,5,5-trimethyl-2-cyclohexen-1-one"), "CC1=CC(=O)CC(C)(C)C1"),
    (("1,3-cyclohexanedione",), "O=C1CCCC(=O)C1"),
    (("dimedone", "5,5-dimethyl-1,3-cyclohexanedione"), "CC1(C)CC(=O)CC(=O)C1"),
    (("1,4-cyclohexanedione",), "O=C1CCC(=O)CC1"),
    (("1,4-cyclohexanedione monoethylene ketal", "1,4-dioxaspiro[4.5]decan-8-one"), "O=C1CCC2(CC1)OCCO2"),
    # ------------------------------------------------------------------ #
    # Saturated N-heterocycle and amine building blocks
    # ------------------------------------------------------------------ #
    (("n-boc-4-piperidone", "1-boc-4-piperidone", "tert-butyl 4-oxopiperidine-1-carboxylate"),
     "CC(C)(C)OC(=O)N1CCC(=O)CC1"),
    (("n-boc-3-pyrrolidinone", "1-boc-3-pyrrolidinone"), "CC(C)(C)OC(=O)N1CCC(=O)C1"),
    (("n-boc-4-piperidinecarboxylic acid", "boc-isonipecotic acid", "1-boc-piperidine-4-carboxylic acid"),
     "CC(C)(C)OC(=O)N1CCC(CC1)C(=O)O"),
    (("n-boc-4-aminopiperidine", "4-amino-1-boc-piperidine"), "CC(C)(C)OC(=O)N1CCC(N)CC1"),
    (("n-boc-4-hydroxypiperidine", "1-boc-4-hydroxypiperidine"), "CC(C)(C)OC(=O)N1CCC(O)CC1"),
    (("n-boc-ethylenediamine", "tert-butyl n-(2-aminoethyl)carbamate", "boc-en"), "CC(C)(C)OC(=O)NCCN"),
    (("n-boc-glycinol", "tert-butyl n-(2-hydroxyethyl)carbamate"), "CC(C)(C)OC(=O)NCCO"),
    (("n-boc-proline", "boc-pro-oh", "boc-l-proline"), "CC(C)(C)OC(=O)N1CCC[C@H]1C(=O)O"),
    (("n-boc-azetidin-3-one", "1-boc-3-azetidinone", "tert-butyl 3-oxoazetidine-1-carboxylate"),
     "CC(C)(C)OC(=O)N1CC(=O)C1"),
    (("n-cbz-piperazine", "1-cbz-piperazine", "benzyl piperazine-1-carboxylate"), "O=C(OCc1ccccc1)N1CCNCC1"),
    (("2-acetylcyclohexanone",), "CC(=O)C1CCCCC1=O"),
    (("diethyl ethylmalonate",), "CCOC(=O)C(CC)C(=O)OCC"),
    # ------------------------------------------------------------------ #
    # Isocyanides, isocyanates and polymer chemistry
    # ------------------------------------------------------------------ #
    (("ethyl isocyanoacetate",), "CCOC(=O)C[N+]#[C-]"),
    (("tosylmethyl isocyanide", "tosmic"), "Cc1ccc(cc1)S(=O)(=O)C[N+]#[C-]"),
    (("tert-butyl isocyanide", "tbunc"), "CC(C)(C)[N+]#[C-]"),
    (("cyclohexyl isocyanide", "cync"), "[C-]#[N+]C1CCCCC1"),
    (("benzyl isocyanide", "bnnc"), "[C-]#[N+]Cc1ccccc1"),
    (("phenyl isocyanate", "phnco"), "O=C=Nc1ccccc1"),
    (("phenyl isothiocyanate", "phncs", "pitc"), "S=C=Nc1ccccc1"),
    (("methyl isocyanate", "mic"), "CN=C=O"),
    (("tosyl isocyanate", "tsnco", "p-toluenesulfonyl isocyanate"), "Cc1ccc(cc1)S(=O)(=O)N=C=O"),
    (("chlorosulfonyl isocyanate", "csi"), "O=C=NS(=O)(=O)Cl"),
    (("trichloroacetyl isocyanate",), "ClC(Cl)(Cl)C(=O)N=C=O"),
    (("hexamethylene diisocyanate", "hdi"), "O=C=NCCCCCCN=C=O"),
    (("toluene diisocyanate", "tdi", "2,4-toluene diisocyanate"), "Cc1ccc(cc1N=C=O)N=C=O"),
    (("methylene diphenyl diisocyanate", "mdi", "4,4'-mdi"), "O=C=Nc1ccc(Cc2ccc(cc2)N=C=O)cc1"),
    (("bisphenol a", "bpa", "2,2-bis(4-hydroxyphenyl)propane"), "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1"),
    (("bisphenol a diglycidyl ether", "badge", "dgeba"), "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"),
    (("lauryl methacrylate",), "CCCCCCCCCCCCOC(=O)C(C)=C"),
    (("adipoyl chloride", "adipoyl dichloride"), "O=C(Cl)CCCCC(=O)Cl"),
    (("sebacoyl chloride", "sebacoyl dichloride"), "O=C(Cl)CCCCCCCCC(=O)Cl"),
    (("terephthaloyl chloride", "terephthaloyl dichloride"), "O=C(Cl)c1ccc(cc1)C(=O)Cl"),
    (("lactide", "l-lactide", "(3s,6s)-3,6-dimethyl-1,4-dioxane-2,5-dione"), "C[C@H]1OC(=O)[C@H](C)OC1=O"),
    (("glycolide", "1,4-dioxane-2,5-dione"), "O=C1COC(=O)CO1"),
    (("aibn", "azobisisobutyronitrile", "2,2'-azobis(2-methylpropionitrile)", "2,2'-azobisisobutyronitrile"),
     "CC(C)(C#N)/N=N/C(C)(C)C#N"),
    (("v-70", "2,2'-azobis(4-methoxy-2,4-dimethylvaleronitrile)"),
     "COCC(C)(C)CC(C)(C#N)/N=N/C(C)(CC(C)(C)COC)C#N"),
    (("acva", "4,4'-azobis(4-cyanovaleric acid)", "v-501"), "OC(=O)CCC(C)(C#N)/N=N/C(C)(CCC(=O)O)C#N"),
    (("irgacure 2959", "2-hydroxy-4'-(2-hydroxyethoxy)-2-methylpropiophenone"),
     "CC(C)(O)C(=O)c1ccc(OCCO)cc1"),
    (("irgacure 651", "dmpa", "2,2-dimethoxy-2-phenylacetophenone"), "COC(OC)(C(=O)c1ccccc1)c1ccccc1"),
    (("lap", "lithium phenyl-2,4,6-trimethylbenzoylphosphinate"),
     "[Li+].Cc1cc(C)c(c(C)c1)C(=O)P(=O)([O-])c1ccccc1"),
    (("bht", "butylated hydroxytoluene", "2,6-di-tert-butyl-4-methylphenol"),
     "Cc1cc(c(O)c(c1)C(C)(C)C)C(C)(C)C"),
    (("bha", "butylated hydroxyanisole", "2-tert-butyl-4-methoxyphenol"), "COc1ccc(O)c(c1)C(C)(C)C"),
    (("galvinoxyl", "galvinoxyl free radical"),
     "CC(C)(C)C1=CC(=CC2=CC(=C([O])C(=C2)C(C)(C)C)C(C)(C)C)C=C(C1=O)C(C)(C)C"),
    (("dpph", "2,2-diphenyl-1-picrylhydrazyl"),
     "[O-][N+](=O)c1cc(cc([N+](=O)[O-])c1[N]N(c1ccccc1)c1ccccc1)[N+](=O)[O-]"),
    (("tbaf trihydrate", "tbaf.3h2o"), "CCCC[N+](CCCC)(CCCC)CCCC.[F-].O.O.O"),
    # ------------------------------------------------------------------ #
    # Aromatic building blocks
    # ------------------------------------------------------------------ #
    (("hexamethylbenzene",), "Cc1c(C)c(C)c(C)c(C)c1C"),
    (("1,3,5-trimethoxybenzene",), "COc1cc(OC)cc(OC)c1"),
    (("1,4-dimethoxybenzene",), "COc1ccc(OC)cc1"),
    (("4-bromoanisole",), "COc1ccc(Br)cc1"),
    (("4-iodoanisole",), "COc1ccc(I)cc1"),
    (("4-bromotoluene", "p-bromotoluene"), "Cc1ccc(Br)cc1"),
    (("4-iodotoluene", "p-iodotoluene"), "Cc1ccc(I)cc1"),
    (("4-chlorotoluene", "p-chlorotoluene"), "Cc1ccc(Cl)cc1"),
    (("1-bromo-4-fluorobenzene", "4-bromofluorobenzene", "4-fluorobromobenzene"), "Fc1ccc(Br)cc1"),
    (("1-bromo-4-chlorobenzene", "4-bromochlorobenzene"), "Clc1ccc(Br)cc1"),
    (("1,4-dibromobenzene", "p-dibromobenzene"), "Brc1ccc(Br)cc1"),
    (("1,4-diiodobenzene",), "Ic1ccc(I)cc1"),
    (("1,3-dibromobenzene",), "Brc1cccc(Br)c1"),
    (("1,2-dibromobenzene",), "Brc1ccccc1Br"),
    (("4-bromobenzonitrile",), "N#Cc1ccc(Br)cc1"),
    (("4-fluorobenzonitrile",), "N#Cc1ccc(F)cc1"),
    (("4-bromonitrobenzene", "1-bromo-4-nitrobenzene"), "[O-][N+](=O)c1ccc(Br)cc1"),
    (("4-fluoronitrobenzene", "1-fluoro-4-nitrobenzene"), "[O-][N+](=O)c1ccc(F)cc1"),
    (("4-chloronitrobenzene", "1-chloro-4-nitrobenzene"), "[O-][N+](=O)c1ccc(Cl)cc1"),
    (("2,4-dinitrofluorobenzene", "sanger's reagent", "sangers reagent", "1-fluoro-2,4-dinitrobenzene", "dnfb", "fdnb"),
     "[O-][N+](=O)c1ccc(F)c(c1)[N+](=O)[O-]"),
    (("2,4-dinitrochlorobenzene", "1-chloro-2,4-dinitrobenzene", "cdnb", "dncb"),
     "[O-][N+](=O)c1ccc(Cl)c(c1)[N+](=O)[O-]"),
    (("nitrobenzene", "phno2"), "O=[N+]([O-])c1ccccc1"),
    (("1,3-dinitrobenzene", "m-dinitrobenzene"), "[O-][N+](=O)c1cccc(c1)[N+](=O)[O-]"),
    (("4-nitrophenol", "p-nitrophenol", "pnp"), "Oc1ccc(cc1)[N+](=O)[O-]"),
    (("2-nitrophenol", "o-nitrophenol"), "Oc1ccccc1[N+](=O)[O-]"),
    (("4-nitrophenyl chloroformate", "p-nitrophenyl chloroformate", "pnp-cl"),
     "O=C(Cl)Oc1ccc(cc1)[N+](=O)[O-]"),
    (("bis(4-nitrophenyl) carbonate",), "O=C(Oc1ccc(cc1)[N+](=O)[O-])Oc1ccc(cc1)[N+](=O)[O-]"),
    (("disuccinimidyl carbonate", "dsc", "n,n'-disuccinimidyl carbonate"),
     "O=C(ON1C(=O)CCC1=O)ON1C(=O)CCC1=O"),
    (("4-nitrophenyl acetate", "p-nitrophenyl acetate", "pnpa"), "CC(=O)Oc1ccc(cc1)[N+](=O)[O-]"),
    (("4-cresol", "p-cresol", "4-methylphenol"), "Cc1ccc(O)cc1"),
    (("2-cresol", "o-cresol", "2-methylphenol"), "Cc1ccccc1O"),
    (("3-cresol", "m-cresol", "3-methylphenol"), "Cc1cccc(O)c1"),
    (("4-chlorophenol", "p-chlorophenol"), "Oc1ccc(Cl)cc1"),
    (("4-bromophenol", "p-bromophenol"), "Oc1ccc(Br)cc1"),
    (("4-fluorophenol", "p-fluorophenol"), "Oc1ccc(F)cc1"),
    (("4-methoxyphenol", "mequinol", "p-methoxyphenol", "mehq", "hydroquinone monomethyl ether"),
     "COc1ccc(O)cc1"),
    (("2-naphthol", "beta-naphthol", "β-naphthol", "naphthalen-2-ol"), "Oc1ccc2ccccc2c1"),
    (("1-naphthol", "alpha-naphthol", "α-naphthol", "naphthalen-1-ol"), "Oc1cccc2ccccc12"),
    (("4-tert-butylphenol",), "CC(C)(C)c1ccc(O)cc1"),
    (("2,6-dimethylphenol", "2,6-xylenol"), "Cc1cccc(C)c1O"),
    (("guaiacol", "2-methoxyphenol"), "COc1ccccc1O"),
    (("eugenol",), "C=CCc1ccc(O)c(OC)c1"),
    (("thymol",), "Cc1ccc(C(C)C)c(O)c1"),
    (("carvacrol",), "Cc1ccc(C(C)C)cc1O"),
    (("pyrogallol", "1,2,3-trihydroxybenzene"), "Oc1cccc(O)c1O"),
    (("phloroglucinol", "1,3,5-trihydroxybenzene"), "Oc1cc(O)cc(O)c1"),
    (("4-aminobenzonitrile",), "Nc1ccc(cc1)C#N"),
    (("4-bromobenzyl bromide",), "BrCc1ccc(Br)cc1"),
    (("4-nitrobenzyl bromide",), "BrCc1ccc(cc1)[N+](=O)[O-]"),
    (("4-methoxybenzylamine", "pmb-nh2", "pmbnh2"), "COc1ccc(CN)cc1"),
    (("2,4-dimethoxybenzylamine", "dmb-nh2"), "COc1ccc(CN)c(OC)c1"),
    (("tryptamine",), "NCCc1c[nH]c2ccccc12"),
    (("2-phenylethanol", "phenethyl alcohol"), "OCCc1ccccc1"),
    (("1-phenylethanol",), "CC(O)c1ccccc1"),
    (("(r)-1-phenylethanol", "(r)-(+)-1-phenylethanol"), "C[C@@H](O)c1ccccc1"),
    (("(s)-1-phenylethanol", "(s)-(-)-1-phenylethanol"), "C[C@H](O)c1ccccc1"),
    (("1-phenylethylamine", "alpha-methylbenzylamine", "α-methylbenzylamine", "dl-1-phenylethylamine"),
     "CC(N)c1ccccc1"),
    (("(r)-1-phenylethylamine", "(r)-(+)-α-methylbenzylamine", "(r)-alpha-methylbenzylamine", "(r)-(+)-1-phenylethylamine"),
     "C[C@@H](N)c1ccccc1"),
    (("(s)-1-phenylethylamine", "(s)-(-)-α-methylbenzylamine", "(s)-alpha-methylbenzylamine", "(s)-(-)-1-phenylethylamine"),
     "C[C@H](N)c1ccccc1"),
    (("benzhydrol", "diphenylmethanol"), "OC(c1ccccc1)c1ccccc1"),
    (("benzhydrylamine", "aminodiphenylmethane"), "NC(c1ccccc1)c1ccccc1"),
    (("triphenylmethanol", "trityl alcohol", "troh"), "OC(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("diphenylmethane",), "C(c1ccccc1)c1ccccc1"),
    (("triphenylmethane",), "C(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("4-methylbenzonitrile", "p-tolunitrile"), "Cc1ccc(cc1)C#N"),
    # ------------------------------------------------------------------ #
    # Heteroaromatic building blocks
    # ------------------------------------------------------------------ #
    (("4-bromopyridine hydrochloride",), "Cl.Brc1ccncc1"),
    (("2-bromopyridine",), "Brc1ccccn1"),
    (("3-bromopyridine",), "Brc1cccnc1"),
    (("2-chloropyridine",), "Clc1ccccn1"),
    (("4-chloropyridine hydrochloride",), "Cl.Clc1ccncc1"),
    (("2-chloropyrimidine",), "Clc1ncccn1"),
    (("2,4-dichloropyrimidine",), "Clc1ccnc(Cl)n1"),
    (("4,6-dichloropyrimidine",), "Clc1cc(Cl)ncn1"),
    (("2-chloropyrazine",), "Clc1cnccn1"),
    (("cyanuric chloride", "2,4,6-trichloro-1,3,5-triazine", "tct"), "Clc1nc(Cl)nc(Cl)n1"),
    (("2-aminopyrimidine",), "Nc1ncccn1"),
    (("2-aminobenzimidazole",), "Nc1nc2ccccc2[nH]1"),
    (("6-chloropurine",), "Clc1ncnc2[nH]cnc12"),
    (("2-chlorobenzoxazole",), "Clc1nc2ccccc2o1"),
    (("2-chlorobenzothiazole",), "Clc1nc2ccccc2s1"),
    (("2-bromothiazole",), "Brc1nccs1"),
    (("2-bromothiophene",), "Brc1cccs1"),
    (("3-bromothiophene",), "Brc1ccsc1"),
    (("2-bromofuran",), "Brc1ccco1"),
    (("4-bromoindole",), "Brc1cccc2[nH]ccc12"),
    (("5-bromoindole",), "Brc1ccc2[nH]ccc2c1"),
    (("5-bromoindazole", "5-bromo-1h-indazole"), "Brc1ccc2[nH]ncc2c1"),
    (("indazole", "1h-indazole"), "c1ccc2[nH]ncc2c1"),
    (("7-azaindole", "1h-pyrrolo[2,3-b]pyridine"), "c1cnc2[nH]ccc2c1"),
    (("5-bromo-2-methoxypyridine",), "COc1ccc(Br)cn1"),
    (("2-bromo-5-methylpyridine",), "Cc1ccc(Br)nc1"),
    (("3-bromoquinoline",), "Brc1cnc2ccccc2c1"),
    (("4-bromo-1h-pyrazole", "4-bromopyrazole"), "Brc1c[nH]nc1"),
    (("4-iodo-1h-pyrazole", "4-iodopyrazole"), "Ic1c[nH]nc1"),
    # ------------------------------------------------------------------ #
    # Boron reagents
    # ------------------------------------------------------------------ #
    (("1-methylpyrazole-4-boronic acid pinacol ester", "1-methyl-4-(4,4,5,5-tetramethyl-1,3,2-dioxaborolan-2-yl)pyrazole"),
     "Cn1cc(cn1)B1OC(C)(C)C(C)(C)O1"),
    (("4-pyridylboronic acid", "pyridine-4-boronic acid"), "OB(O)c1ccncc1"),
    (("3-pyridylboronic acid", "pyridine-3-boronic acid"), "OB(O)c1cccnc1"),
    (("2-thienylboronic acid", "thiophene-2-boronic acid"), "OB(O)c1cccs1"),
    (("4-methoxyphenylboronic acid", "4-anisylboronic acid"), "COc1ccc(cc1)B(O)O"),
    (("4-fluorophenylboronic acid",), "OB(O)c1ccc(F)cc1"),
    (("4-chlorophenylboronic acid",), "OB(O)c1ccc(Cl)cc1"),
    (("4-tolylboronic acid", "p-tolylboronic acid", "4-methylphenylboronic acid"), "Cc1ccc(cc1)B(O)O"),
    (("4-(trifluoromethyl)phenylboronic acid",), "OB(O)c1ccc(cc1)C(F)(F)F"),
    (("4-cyanophenylboronic acid",), "OB(O)c1ccc(cc1)C#N"),
    (("3-nitrophenylboronic acid",), "OB(O)c1cccc(c1)[N+](=O)[O-]"),
    (("vinylboronic acid pinacol ester", "vinyl-bpin"), "C=CB1OC(C)(C)C(C)(C)O1"),
    (("phenylboronic acid pinacol ester", "ph-bpin", "phbpin", "2-phenyl-4,4,5,5-tetramethyl-1,3,2-dioxaborolane"),
     "CC1(C)OB(OC1(C)C)c1ccccc1"),
    (("cyclopropylboronic acid",), "OB(O)C1CC1"),
    (("methylboronic acid", "meb(oh)2"), "CB(O)O"),
    (("potassium vinyltrifluoroborate",), "[K+].C=C[B-](F)(F)F"),
    (("potassium cyclopropyltrifluoroborate",), "[K+].F[B-](F)(F)C1CC1"),
    (("allylboronic acid pinacol ester", "allyl-bpin", "allylbpin"), "C=CCB1OC(C)(C)C(C)(C)O1"),
    (("isopropoxyboronic acid pinacol ester", "2-isopropoxy-4,4,5,5-tetramethyl-1,3,2-dioxaborolane", "ipro-bpin", "iprobpin"),
     "CC(C)OB1OC(C)(C)C(C)(C)O1"),
    (("trimethyl borate", "b(ome)3"), "COB(OC)OC"),
    (("triisopropyl borate", "b(oipr)3"), "CC(C)OB(OC(C)C)OC(C)C"),
    (("triethylborane", "et3b", "bet3"), "CCB(CC)CC"),
    (("tributyltin azide",), "CCCC[Sn](CCCC)(CCCC)N=[N+]=[N-]"),
    # ------------------------------------------------------------------ #
    # Ylides, phosphonates and sulfur ylides
    # ------------------------------------------------------------------ #
    (("wittig reagent", "methylenetriphenylphosphorane", "ph3p=ch2"), "C=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("methyltriphenylphosphonium bromide", "ph3pmebr", "mepph3br", "ch3pph3br"),
     "C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[Br-]"),
    (("methyltriphenylphosphonium iodide", "ph3pmei", "mepph3i"), "C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[I-]"),
    (("ethyltriphenylphosphonium bromide", "etpph3br"), "CC[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[Br-]"),
    (("benzyltriphenylphosphonium chloride", "bnpph3cl"),
     "[Cl-].c1ccc(cc1)C[P+](c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("(methoxymethyl)triphenylphosphonium chloride", "ph3pch2omecl", "meoch2pph3cl"),
     "COC[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[Cl-]"),
    (("(carbethoxymethylene)triphenylphosphorane", "ph3p=chco2et", "ethyl (triphenylphosphoranylidene)acetate"),
     "CCOC(=O)C=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("(carbomethoxymethylene)triphenylphosphorane", "ph3p=chco2me", "methyl (triphenylphosphoranylidene)acetate"),
     "COC(=O)C=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("(triphenylphosphoranylidene)acetaldehyde", "ph3p=chcho"), "O=CC=P(c1ccccc1)(c1ccccc1)c1ccccc1"),
    (("triethyl phosphonoacetate", "(eto)2p(o)ch2co2et", "hwe reagent"), "CCOC(=O)CP(=O)(OCC)OCC"),
    (("trimethyl phosphonoacetate", "(meo)2p(o)ch2co2me", "methyl 2-(dimethoxyphosphoryl)acetate"),
     "COC(=O)CP(=O)(OC)OC"),
    (("diethyl (cyanomethyl)phosphonate",), "CCOP(=O)(CC#N)OCC"),
    (("bestmann-ohira reagent", "ohira-bestmann reagent", "dimethyl (1-diazo-2-oxopropyl)phosphonate"),
     "CC(=O)C(=[N+]=[N-])P(=O)(OC)OC"),
    (("diethyl benzylphosphonate",), "CCOP(=O)(Cc1ccccc1)OCC"),
    (("tetraethyl methylenediphosphonate",), "CCOP(=O)(CP(=O)(OCC)OCC)OCC"),
    (("ethyl 2-(diethoxyphosphoryl)propanoate", "triethyl 2-phosphonopropionate"),
     "CCOC(=O)C(C)P(=O)(OCC)OCC"),
    (("still-gennari reagent", "bis(2,2,2-trifluoroethyl) (methoxycarbonylmethyl)phosphonate"),
     "COC(=O)CP(=O)(OCC(F)(F)F)OCC(F)(F)F"),
    (("dimethyl sulfoxonium methylide", "corey's ylide", "corey ylide", "corey-chaykovsky reagent"),
     "C[S+](C)([CH2-])=O"),
    (("trimethylsulfoxonium iodide", "me3soi", "me3s(o)i"), "C[S+](C)(C)=O.[I-]"),
    (("trimethylsulfonium iodide",), "C[S+](C)C.[I-]"),
    # ------------------------------------------------------------------ #
    # Nitrosyl, diazonium, iodonium and amination reagents
    # ------------------------------------------------------------------ #
    (("nitrosyl chloride", "nocl"), "O=NCl"),
    (("benzenediazonium tetrafluoroborate", "phn2bf4", "phn2+bf4-"), "[N+](#N)c1ccccc1.F[B-](F)(F)F"),
    (("4-methoxybenzenediazonium tetrafluoroborate",), "COc1ccc(cc1)[N+]#N.F[B-](F)(F)F"),
    (("diphenyliodonium triflate", "ph2iotf"), "[I+](c1ccccc1)c1ccccc1.[O-]S(=O)(=O)C(F)(F)F"),
    (("diphenyliodonium tetrafluoroborate",), "[I+](c1ccccc1)c1ccccc1.F[B-](F)(F)F"),
    (("triphenylsulfonium triflate", "ph3sotf"), "[S+](c1ccccc1)(c1ccccc1)c1ccccc1.[O-]S(=O)(=O)C(F)(F)F"),
    (("n-tosylimidazole", "1-tosylimidazole", "tsim"), "Cc1ccc(cc1)S(=O)(=O)n1ccnc1"),
    (("n-acetylimidazole", "1-acetylimidazole"), "CC(=O)n1ccnc1"),
    (("n,n'-carbonyldi(1,2,4-triazole)", "cdt", "1,1'-carbonyldi(1,2,4-triazole)"), "O=C(n1cncn1)n1cncn1"),
    (("di-tert-butyl carbonate",), "CC(C)(C)OC(=O)OC(C)(C)C"),
    (("dimethyl pyrocarbonate", "dimethyl dicarbonate"), "COC(=O)OC(=O)OC"),
    (("diethyl pyrocarbonate", "depc", "diethyl dicarbonate"), "CCOC(=O)OC(=O)OCC"),
    (("acetone cyanohydrin", "2-hydroxy-2-methylpropanenitrile"), "CC(C)(O)C#N"),
    (("acetone oxime", "propan-2-one oxime"), "CC(C)=NO"),
    (("acetaldoxime",), "CC=NO"),
    (("methoxyamine", "o-methylhydroxylamine", "meonh2"), "CON"),
    (("hydroxylamine-o-sulfonic acid", "hosa", "hydroxylamine o-sulfonic acid"), "NOS(=O)(=O)O"),
    (("o-(mesitylsulfonyl)hydroxylamine", "msh"), "Cc1cc(C)c(c(C)c1)S(=O)(=O)ON"),
    (("o-(diphenylphosphinyl)hydroxylamine", "diphenylphosphinyl hydroxylamine"),
     "NOP(=O)(c1ccccc1)c1ccccc1"),
    (("ethyl nitroacetate",), "CCOC(=O)C[N+](=O)[O-]"),
    (("nitroethylene",), "C=C[N+](=O)[O-]"),
    (("beta-nitrostyrene", "β-nitrostyrene", "trans-beta-nitrostyrene", "(e)-nitrostyrene", "nitrostyrene"),
     "[O-][N+](=O)/C=C/c1ccccc1"),
    # ------------------------------------------------------------------ #
    # Sulfur and thiol reagents
    # ------------------------------------------------------------------ #
    (("phenyl vinyl sulfone",), "C=CS(=O)(=O)c1ccccc1"),
    (("methyl vinyl sulfone",), "C=CS(C)(=O)=O"),
    (("phenyl vinyl sulfoxide",), "C=CS(=O)c1ccccc1"),
    (("divinyl sulfone",), "C=CS(=O)(=O)C=C"),
    (("methyl phenyl sulfone", "phso2me"), "CS(=O)(=O)c1ccccc1"),
    (("dimethyl sulfone", "msm", "methylsulfonylmethane"), "CS(C)(=O)=O"),
    (("diphenyl sulfide", "ph2s"), "S(c1ccccc1)c1ccccc1"),
    (("diphenyl disulfide", "ph2s2", "phssph"), "c1ccc(cc1)SSc1ccccc1"),
    (("dimethyl disulfide", "me2s2", "dmds"), "CSSC"),
    (("benzyl mercaptan", "benzyl thiol", "phenylmethanethiol", "bnsh"), "SCc1ccccc1"),
    (("ethanethiol", "ethyl mercaptan", "etsh"), "CCS"),
    (("methanethiol", "methyl mercaptan", "mesh"), "CS"),
    (("tert-butyl mercaptan", "2-methyl-2-propanethiol", "tbush"), "CC(C)(C)S"),
    (("thioacetamide",), "CC(N)=S"),
    (("sodium thiomethoxide", "nasme", "sodium methanethiolate"), "[Na+].C[S-]"),
    (("sodium thiophenoxide", "nasph", "sodium benzenethiolate"), "[Na+].[S-]c1ccccc1"),
    (("potassium thioacetate", "ksac", "potassium thiolacetate", "acsk"), "[K+].CC([S-])=O"),
    (("potassium ethyl xanthate", "potassium o-ethyl xanthate"), "[K+].CCOC([S-])=S"),
    (("s-methyl methanethiosulfonate", "mmts", "methyl methanethiosulfonate"), "CSS(C)(=O)=O"),
    (("n-methylmaleimide", "1-methylmaleimide"), "CN1C(=O)C=CC1=O"),
    # ------------------------------------------------------------------ #
    # Imides, amidines, acetals and orthoesters
    # ------------------------------------------------------------------ #
    (("maleimide",), "O=C1NC(=O)C=C1"),
    (("n-phenylmaleimide",), "O=C1C=CC(=O)N1c1ccccc1"),
    (("n-ethylmaleimide", "nem"), "CCN1C(=O)C=CC1=O"),
    (("succinimide",), "O=C1CCC(=O)N1"),
    (("glutarimide",), "O=C1CCCC(=O)N1"),
    (("dimethylformamide dimethyl acetal", "dmf-dma", "dmfdma", "n,n-dimethylformamide dimethyl acetal"),
     "COC(OC)N(C)C"),
    (("bredereck's reagent", "bredereck reagent", "tert-butoxybis(dimethylamino)methane"),
     "CC(C)(C)OC(N(C)C)N(C)C"),
    (("dmf-dea", "n,n-dimethylformamide diethyl acetal"), "CCOC(OCC)N(C)C"),
    (("gold's reagent", "golds reagent"), "CN(C)C=[N+](C)C.[Cl-]"),
    (("eschenmoser's salt", "eschenmoser salt", "n,n-dimethylmethyleneiminium iodide", "dimethylmethyleneammonium iodide"),
     "C=[N+](C)C.[I-]"),
    (("acetal", "1,1-diethoxyethane", "acetaldehyde diethyl acetal"), "CCOC(C)OCC"),
    (("dimethoxymethane", "methylal", "formaldehyde dimethyl acetal"), "COCOC"),
    (("paraldehyde",), "CC1OC(C)OC(C)O1"),
    (("trioxane", "1,3,5-trioxane"), "C1OCOCO1"),
    (("tetramethylethylene", "2,3-dimethyl-2-butene"), "CC(C)=C(C)C"),
    (("methylenecyclohexane",), "C=C1CCCCC1"),
    # ------------------------------------------------------------------ #
    # Carbocycles, cages and sp3-rich scaffolds
    # ------------------------------------------------------------------ #
    (("methyl cyclopropyl ketone", "cyclopropyl methyl ketone"), "CC(=O)C1CC1"),
    (("cyclopropylamine",), "NC1CC1"),
    (("cyclopropanecarboxylic acid",), "OC(=O)C1CC1"),
    (("cyclobutylamine",), "NC1CCC1"),
    (("cyclopentylamine",), "NC1CCCC1"),
    (("cyclohexanecarboxylic acid",), "OC(=O)C1CCCCC1"),
    (("adamantane",), "C1C2CC3CC1CC(C2)C3"),
    (("1-adamantanol",), "OC12CC3CC(CC(C3)C1)C2"),
    (("1-adamantylamine", "amantadine", "adamantan-1-amine"), "NC12CC3CC(CC(C3)C1)C2"),
    (("1-adamantanecarboxylic acid",), "OC(=O)C12CC3CC(CC(C3)C1)C2"),
    (("memantine",), "CC12CC3CC(C)(C1)CC(N)(C3)C2"),
    (("cubane",), "C12C3C4C1C1C2C3C41"),
    (("bicyclo[1.1.1]pentane", "bcp"), "C1C2CC1C2"),
    (("1,4-bcp-dicarboxylic acid", "bicyclo[1.1.1]pentane-1,3-dicarboxylic acid"),
     "OC(=O)C12CC(C1)(C2)C(=O)O"),
    (("spiro[3.3]heptane",), "C1CC2(C1)CCC2"),
    (("3-azabicyclo[3.1.0]hexane",), "C1NCC2CC12"),
    (("2-oxa-6-azaspiro[3.3]heptane",), "C1OCC12CNC2"),
    # ------------------------------------------------------------------ #
    # Amino acids (L unless stated; L = S except cysteine/selenocysteine R), D-enantiomers,
    # non-canonical residues, protected amino acids and short peptides
    # ------------------------------------------------------------------ #
    (("glycine", "gly", "g", "aminoacetic acid", "2-aminoacetic acid", "l-glycine", "d-glycine"), "NCC(=O)O"),
    (("alanine", "l-alanine", "ala", "l-ala", "(s)-alanine", "(s)-2-aminopropanoic acid"), "C[C@H](N)C(=O)O"),
    (("d-alanine", "d-ala", "(r)-alanine"), "C[C@@H](N)C(=O)O"),
    (("dl-alanine", "rac-alanine", "racemic alanine"), "CC(N)C(=O)O"),
    (("valine", "l-valine", "val", "l-val", "(s)-valine"), "CC(C)[C@H](N)C(=O)O"),
    (("d-valine", "d-val"), "CC(C)[C@@H](N)C(=O)O"),
    (("leucine", "l-leucine", "leu", "l-leu", "(s)-leucine"), "CC(C)C[C@H](N)C(=O)O"),
    (("d-leucine", "d-leu"), "CC(C)C[C@@H](N)C(=O)O"),
    (("isoleucine", "l-isoleucine", "ile", "l-ile", "(2s,3s)-isoleucine"), "CC[C@H](C)[C@H](N)C(=O)O"),
    (("d-isoleucine", "d-ile", "(2r,3r)-isoleucine"), "CC[C@@H](C)[C@@H](N)C(=O)O"),
    (("allo-isoleucine", "l-allo-isoleucine", "l-alloisoleucine", "(2s,3r)-isoleucine"),
     "CC[C@@H](C)[C@H](N)C(=O)O"),
    (("proline", "l-proline", "pro", "l-pro", "(s)-proline"), "OC(=O)[C@@H]1CCCN1"),
    (("d-proline", "d-pro", "(r)-proline"), "OC(=O)[C@H]1CCCN1"),
    (("phenylalanine", "l-phenylalanine", "phe", "l-phe", "(s)-phenylalanine"), "N[C@@H](Cc1ccccc1)C(=O)O"),
    (("d-phenylalanine", "d-phe"), "N[C@H](Cc1ccccc1)C(=O)O"),
    (("tyrosine", "l-tyrosine", "tyr", "l-tyr", "(s)-tyrosine"), "N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
    (("d-tyrosine", "d-tyr"), "N[C@H](Cc1ccc(O)cc1)C(=O)O"),
    (("tryptophan", "l-tryptophan", "trp", "l-trp", "(s)-tryptophan"), "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    (("d-tryptophan", "d-trp"), "N[C@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    (("serine", "l-serine", "ser", "l-ser", "(s)-serine"), "N[C@@H](CO)C(=O)O"),
    (("d-serine", "d-ser"), "N[C@H](CO)C(=O)O"),
    (("threonine", "l-threonine", "thr", "l-thr", "(2s,3r)-threonine"), "C[C@@H](O)[C@H](N)C(=O)O"),
    (("d-threonine", "d-thr", "(2r,3s)-threonine"), "C[C@H](O)[C@@H](N)C(=O)O"),
    (("allo-threonine", "l-allo-threonine", "l-allothreonine", "(2s,3s)-threonine"),
     "C[C@H](O)[C@H](N)C(=O)O"),
    (("cysteine", "l-cysteine", "cys", "l-cys", "(r)-cysteine"), "N[C@@H](CS)C(=O)O"),
    (("d-cysteine", "d-cys"), "N[C@H](CS)C(=O)O"),
    (("cysteine hydrochloride", "l-cysteine hydrochloride", "l-cysteine hcl", "cys.hcl"),
     "Cl.N[C@@H](CS)C(=O)O"),
    (("cystine", "l-cystine", "cys-cys", "(cys)2"), "N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O"),
    (("methionine", "l-methionine", "met", "l-met", "(s)-methionine"), "CSCC[C@H](N)C(=O)O"),
    (("d-methionine", "d-met"), "CSCC[C@@H](N)C(=O)O"),
    (("asparagine", "l-asparagine", "asn", "l-asn", "(s)-asparagine"), "NC(=O)C[C@H](N)C(=O)O"),
    (("d-asparagine", "d-asn"), "NC(=O)C[C@@H](N)C(=O)O"),
    (("glutamine", "l-glutamine", "gln", "l-gln", "(s)-glutamine"), "NC(=O)CC[C@H](N)C(=O)O"),
    (("d-glutamine", "d-gln"), "NC(=O)CC[C@@H](N)C(=O)O"),
    (("aspartic acid", "l-aspartic acid", "aspartate", "asp", "l-asp", "l-aspartate", "(s)-aspartic acid"),
     "N[C@@H](CC(=O)O)C(=O)O"),
    (("d-aspartic acid", "d-asp", "d-aspartate"), "N[C@H](CC(=O)O)C(=O)O"),
    (("glutamic acid", "l-glutamic acid", "glutamate", "glu", "l-glu", "l-glutamate", "(s)-glutamic acid"),
     "N[C@@H](CCC(=O)O)C(=O)O"),
    (("d-glutamic acid", "d-glu", "d-glutamate"), "N[C@H](CCC(=O)O)C(=O)O"),
    (("monosodium glutamate", "msg", "sodium glutamate", "sodium l-glutamate"),
     "[Na+].N[C@@H](CCC(=O)O)C(=O)[O-]"),
    (("lysine", "l-lysine", "lys", "l-lys", "(s)-lysine"), "NCCCC[C@H](N)C(=O)O"),
    (("d-lysine", "d-lys"), "NCCCC[C@@H](N)C(=O)O"),
    (("lysine hydrochloride", "l-lysine hydrochloride", "l-lysine hcl", "lys.hcl"), "Cl.NCCCC[C@H](N)C(=O)O"),
    (("arginine", "l-arginine", "arg", "l-arg", "(s)-arginine"), "N=C(N)NCCC[C@H](N)C(=O)O"),
    (("d-arginine", "d-arg"), "N=C(N)NCCC[C@@H](N)C(=O)O"),
    (("arginine hydrochloride", "l-arginine hydrochloride", "l-arginine hcl", "arg.hcl"),
     "Cl.N=C(N)NCCC[C@H](N)C(=O)O"),
    (("histidine", "l-histidine", "his", "l-his", "(s)-histidine"), "N[C@@H](Cc1c[nH]cn1)C(=O)O"),
    (("d-histidine", "d-his"), "N[C@H](Cc1c[nH]cn1)C(=O)O"),
    # ------------------------------------------------------------------ #
    # Non-canonical and modified amino acids
    # ------------------------------------------------------------------ #
    (("selenocysteine", "l-selenocysteine", "sec", "u", "(r)-selenocysteine"), "N[C@@H](C[SeH])C(=O)O"),
    (("pyrrolysine", "l-pyrrolysine", "pyl", "o"), "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N)C(=O)O"),
    (("ornithine", "l-ornithine", "orn", "l-orn", "(s)-ornithine", "2,5-diaminopentanoic acid"),
     "NCCC[C@H](N)C(=O)O"),
    (("d-ornithine", "d-orn"), "NCCC[C@@H](N)C(=O)O"),
    (("ornithine hydrochloride", "l-ornithine hydrochloride", "l-ornithine hcl"), "Cl.NCCC[C@H](N)C(=O)O"),
    (("citrulline", "l-citrulline", "cit", "(s)-citrulline"), "NC(=O)NCCC[C@H](N)C(=O)O"),
    (("homoserine", "l-homoserine", "hse", "(s)-homoserine"), "N[C@@H](CCO)C(=O)O"),
    (("homocysteine", "l-homocysteine", "hcy", "(s)-homocysteine"), "N[C@@H](CCS)C(=O)O"),
    (("homoarginine", "l-homoarginine", "har"), "N=C(N)NCCCC[C@H](N)C(=O)O"),
    (("homophenylalanine", "l-homophenylalanine", "hphe"), "N[C@@H](CCc1ccccc1)C(=O)O"),
    (("norleucine", "l-norleucine", "nle", "(s)-2-aminohexanoic acid"), "CCCC[C@H](N)C(=O)O"),
    (("norvaline", "l-norvaline", "nva", "(s)-2-aminopentanoic acid"), "CCC[C@H](N)C(=O)O"),
    (("tert-leucine", "l-tert-leucine", "tle", "tert-butylglycine", "(s)-tert-leucine", "l-tle"),
     "CC(C)(C)[C@H](N)C(=O)O"),
    (("aib", "2-aminoisobutyric acid", "alpha-aminoisobutyric acid", "α-aminoisobutyric acid", "2-methylalanine", "alpha-methylalanine", "α-methylalanine"),
     "CC(C)(N)C(=O)O"),
    (("dap", "2,3-diaminopropionic acid", "l-dap", "(s)-2,3-diaminopropanoic acid", "2,3-diaminopropanoic acid", "l-2,3-diaminopropionic acid"),
     "NC[C@H](N)C(=O)O"),
    (("dab", "2,4-diaminobutyric acid", "l-dab", "(s)-2,4-diaminobutanoic acid", "2,4-diaminobutanoic acid", "l-2,4-diaminobutyric acid"),
     "NCC[C@H](N)C(=O)O"),
    (("beta-alanine", "β-alanine", "b-alanine", "beta alanine", "3-aminopropanoic acid", "3-aminopropionic acid", "bala", "β-ala"),
     "NCCC(=O)O"),
    (("gaba", "gamma-aminobutyric acid", "γ-aminobutyric acid", "4-aminobutanoic acid", "4-aminobutyric acid", "gamma aminobutyric acid"),
     "NCCCC(=O)O"),
    (("5-aminovaleric acid", "5-aminopentanoic acid", "ava"), "NCCCCC(=O)O"),
    (("6-aminohexanoic acid", "6-aminocaproic acid", "aminocaproic acid", "ahx", "epsilon-aminocaproic acid", "ε-aminocaproic acid", "eaca"),
     "NCCCCCC(=O)O"),
    (("8-aminooctanoic acid",), "NCCCCCCCC(=O)O"),
    (("sarcosine", "n-methylglycine", "sar"), "CNCC(=O)O"),
    (("n-methylalanine", "l-n-methylalanine", "meala", "n-me-ala"), "CN[C@@H](C)C(=O)O"),
    (("n,n-dimethylglycine", "dimethylglycine", "dmg"), "CN(C)CC(=O)O"),
    (("betaine", "trimethylglycine", "glycine betaine", "tmg (betaine)"), "C[N+](C)(C)CC(=O)[O-]"),
    (("creatine",), "CN(CC(=O)O)C(N)=N"),
    (("creatinine",), "CN1CC(=O)N=C1N"),
    (("carnitine", "l-carnitine", "(r)-carnitine"), "C[N+](C)(C)C[C@H](O)CC(=O)[O-]"),
    (("acetylcarnitine", "acetyl-l-carnitine", "alcar"), "CC(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C"),
    (("phenylglycine", "l-phenylglycine", "(s)-phenylglycine", "phg", "2-phenylglycine"),
     "N[C@H](C(=O)O)c1ccccc1"),
    (("d-phenylglycine", "(r)-phenylglycine"), "N[C@@H](C(=O)O)c1ccccc1"),
    (("4-hydroxyphenylglycine", "d-4-hydroxyphenylglycine", "(r)-4-hydroxyphenylglycine"),
     "N[C@@H](C(=O)O)c1ccc(O)cc1"),
    (("cyclohexylglycine", "l-cyclohexylglycine", "chg"), "N[C@H](C(=O)O)C1CCCCC1"),
    (("cyclohexylalanine", "l-cyclohexylalanine", "cha", "3-cyclohexyl-l-alanine"),
     "N[C@@H](CC1CCCCC1)C(=O)O"),
    (("1-naphthylalanine", "l-1-naphthylalanine", "1-nal"), "N[C@@H](Cc1cccc2ccccc12)C(=O)O"),
    (("2-naphthylalanine", "l-2-naphthylalanine", "2-nal"), "N[C@@H](Cc1ccc2ccccc2c1)C(=O)O"),
    (("4-fluorophenylalanine", "l-4-fluorophenylalanine", "4-f-phe", "phe(4-f)", "p-fluorophenylalanine"),
     "N[C@@H](Cc1ccc(F)cc1)C(=O)O"),
    (("4-chlorophenylalanine", "l-4-chlorophenylalanine", "4-cl-phe", "phe(4-cl)"),
     "N[C@@H](Cc1ccc(Cl)cc1)C(=O)O"),
    (("4-bromophenylalanine", "l-4-bromophenylalanine", "4-br-phe"), "N[C@@H](Cc1ccc(Br)cc1)C(=O)O"),
    (("4-iodophenylalanine", "l-4-iodophenylalanine", "4-i-phe"), "N[C@@H](Cc1ccc(I)cc1)C(=O)O"),
    (("4-nitrophenylalanine", "l-4-nitrophenylalanine", "4-no2-phe"),
     "N[C@@H](Cc1ccc(cc1)[N+](=O)[O-])C(=O)O"),
    (("4-azidophenylalanine", "l-4-azidophenylalanine", "azf", "p-azido-l-phenylalanine", "azpha"),
     "N[C@@H](Cc1ccc(cc1)N=[N+]=[N-])C(=O)O"),
    (("4-benzoylphenylalanine", "p-benzoyl-l-phenylalanine", "l-bpa", "benzoylphenylalanine"),
     "N[C@@H](Cc1ccc(cc1)C(=O)c1ccccc1)C(=O)O"),
    (("4-acetylphenylalanine", "p-acetyl-l-phenylalanine", "acf"), "CC(=O)c1ccc(C[C@H](N)C(=O)O)cc1"),
    (("pentafluorophenylalanine", "l-pentafluorophenylalanine"), "N[C@@H](Cc1c(F)c(F)c(F)c(F)c1F)C(=O)O"),
    (("3,4-dihydroxyphenylalanine", "l-dopa", "levodopa", "dopa", "l-3,4-dihydroxyphenylalanine"),
     "N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O"),
    (("3-iodotyrosine", "3-iodo-l-tyrosine"), "N[C@@H](Cc1ccc(O)c(I)c1)C(=O)O"),
    (("3-nitrotyrosine", "3-nitro-l-tyrosine"), "N[C@@H](Cc1ccc(O)c(c1)[N+](=O)[O-])C(=O)O"),
    (("o-methyltyrosine", "4-methoxyphenylalanine", "tyr(me)", "o-methyl-l-tyrosine"),
     "COc1ccc(C[C@H](N)C(=O)O)cc1"),
    (("phosphotyrosine", "o-phospho-l-tyrosine", "ptyr", "p-tyr", "tyr(po3h2)"),
     "N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O"),
    (("phosphoserine", "o-phospho-l-serine", "pser", "p-ser", "ser(po3h2)"), "N[C@@H](COP(=O)(O)O)C(=O)O"),
    (("phosphothreonine", "o-phospho-l-threonine", "pthr", "p-thr"), "C[C@@H](OP(=O)(O)O)[C@H](N)C(=O)O"),
    (("5-hydroxytryptophan", "5-htp", "l-5-hydroxytryptophan", "oxitriptan"),
     "N[C@@H](Cc1c[nH]c2ccc(O)cc12)C(=O)O"),
    (("hydroxyproline", "l-hydroxyproline", "trans-4-hydroxy-l-proline", "4-hydroxyproline", "hyp", "(2s,4r)-4-hydroxyproline"),
     "O[C@H]1CN[C@@H](C1)C(=O)O"),
    (("cis-4-hydroxy-l-proline", "(2s,4s)-4-hydroxyproline"), "O[C@@H]1CN[C@@H](C1)C(=O)O"),
    (("4-fluoroproline", "(2s,4r)-4-fluoroproline", "trans-4-fluoro-l-proline"), "F[C@H]1CN[C@@H](C1)C(=O)O"),
    (("pipecolic acid", "l-pipecolic acid", "pip", "(s)-piperidine-2-carboxylic acid", "homoproline"),
     "OC(=O)[C@@H]1CCCCN1"),
    (("azetidine-2-carboxylic acid", "l-azetidine-2-carboxylic acid", "aze"), "OC(=O)[C@@H]1CCN1"),
    (("thiazolidine-4-carboxylic acid", "l-thioproline", "thioproline", "thz"), "OC(=O)[C@@H]1CSCN1"),
    (("tic", "tetrahydroisoquinoline-3-carboxylic acid", "(s)-tic", "l-tic"), "OC(=O)[C@@H]1Cc2ccccc2CN1"),
    (("allylglycine", "l-allylglycine", "(s)-2-aminopent-4-enoic acid"), "C=CC[C@H](N)C(=O)O"),
    (("propargylglycine", "l-propargylglycine", "(s)-2-aminopent-4-ynoic acid", "pra"), "C#CC[C@H](N)C(=O)O"),
    (("homopropargylglycine", "l-homopropargylglycine", "hpg", "(s)-2-aminohex-5-ynoic acid"),
     "C#CCC[C@H](N)C(=O)O"),
    (("azidohomoalanine", "l-azidohomoalanine", "aha", "(s)-2-amino-4-azidobutanoic acid"),
     "N[C@@H](CCN=[N+]=[N-])C(=O)O"),
    (("azidolysine", "l-azidolysine", "6-azido-l-lysine", "lys(n3)", "azk"),
     "N[C@@H](CCCCN=[N+]=[N-])C(=O)O"),
    (("azidoalanine", "l-azidoalanine", "3-azido-l-alanine"), "N[C@@H](CN=[N+]=[N-])C(=O)O"),
    (("acetyllysine", "n-epsilon-acetyl-l-lysine", "nε-acetyllysine", "n-ε-acetyl-l-lysine", "ack", "lys(ac)", "acetyl-lysine"),
     "CC(=O)NCCCC[C@H](N)C(=O)O"),
    (("methyllysine", "n-epsilon-methyl-l-lysine", "lys(me)", "monomethyllysine"), "CNCCCC[C@H](N)C(=O)O"),
    (("dimethyllysine", "n-epsilon-dimethyl-l-lysine", "lys(me2)"), "CN(C)CCCC[C@H](N)C(=O)O"),
    (("trimethyllysine", "n-epsilon-trimethyl-l-lysine", "lys(me3)"), "C[N+](C)(C)CCCC[C@H](N)C(=O)[O-]"),
    (("asymmetric dimethylarginine", "adma", "n,n-dimethylarginine"), "CN(C)C(=N)NCCC[C@H](N)C(=O)O"),
    (("symmetric dimethylarginine", "sdma"), "CN=C(NC)NCCC[C@H](N)C(=O)O"),
    (("monomethylarginine", "n-methylarginine", "l-nmma", "nmma"), "CNC(=N)NCCC[C@H](N)C(=O)O"),
    (("nitroarginine", "n-omega-nitro-l-arginine", "nω-nitro-l-arginine", "l-nna"),
     "N[C@@H](CCCNC(=N)N[N+](=O)[O-])C(=O)O"),
    (("l-name", "n-omega-nitro-l-arginine methyl ester", "nω-nitro-l-arginine methyl ester"),
     "COC(=O)[C@@H](N)CCCNC(=N)N[N+](=O)[O-]"),
    (("selenomethionine", "l-selenomethionine", "sem", "semet"), "C[Se]CC[C@H](N)C(=O)O"),
    (("methionine sulfoxide", "l-methionine sulfoxide", "met(o)"), "CS(=O)CC[C@H](N)C(=O)O"),
    (("methionine sulfone", "met(o2)"), "CS(=O)(=O)CC[C@H](N)C(=O)O"),
    (("s-methylcysteine", "s-methyl-l-cysteine"), "CSC[C@H](N)C(=O)O"),
    (("s-benzylcysteine", "s-benzyl-l-cysteine", "cys(bzl)"), "N[C@@H](CSCc1ccccc1)C(=O)O"),
    (("penicillamine", "d-penicillamine", "(s)-penicillamine", "3-mercapto-d-valine", "cuprimine"),
     "CC(C)(S)[C@@H](N)C(=O)O"),
    (("n-acetylcysteine", "nac", "n-acetyl-l-cysteine", "acetylcysteine"), "CC(=O)N[C@@H](CS)C(=O)O"),
    (("n-acetylglycine", "aceturic acid", "ac-gly-oh"), "CC(=O)NCC(=O)O"),
    (("n-acetylalanine", "ac-ala-oh", "n-acetyl-l-alanine"), "CC(=O)N[C@@H](C)C(=O)O"),
    (("acetylphenylalanine", "ac-phe-oh", "n-acetyl-l-phenylalanine"), "CC(=O)N[C@@H](Cc1ccccc1)C(=O)O"),
    (("n-acetyltryptophan", "ac-trp-oh"), "CC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    # ------------------------------------------------------------------ #
    # Protected amino acids, esters and short peptides
    # ------------------------------------------------------------------ #
    (("glycinamide", "gly-nh2", "h-gly-nh2"), "NCC(N)=O"),
    (("glycinamide hydrochloride", "gly-nh2.hcl"), "Cl.NCC(N)=O"),
    (("glycine methyl ester hydrochloride", "h-gly-ome.hcl", "gly-ome.hcl"), "Cl.COC(=O)CN"),
    (("glycine ethyl ester hydrochloride", "h-gly-oet.hcl", "gly-oet.hcl"), "Cl.CCOC(=O)CN"),
    (("glycine tert-butyl ester hydrochloride", "h-gly-otbu.hcl"), "Cl.CC(C)(C)OC(=O)CN"),
    (("glycine benzyl ester p-toluenesulfonate", "h-gly-obzl.tsoh"),
     "Cc1ccc(cc1)S(=O)(=O)O.NCC(=O)OCc1ccccc1"),
    (("alanine methyl ester hydrochloride", "h-ala-ome.hcl", "l-alanine methyl ester hydrochloride"),
     "Cl.COC(=O)[C@H](C)N"),
    (("phenylalanine methyl ester hydrochloride", "h-phe-ome.hcl"), "Cl.COC(=O)[C@@H](N)Cc1ccccc1"),
    (("leucine methyl ester hydrochloride", "h-leu-ome.hcl"), "Cl.COC(=O)[C@@H](N)CC(C)C"),
    (("proline methyl ester hydrochloride", "h-pro-ome.hcl"), "Cl.COC(=O)[C@@H]1CCCN1"),
    (("valine methyl ester hydrochloride", "h-val-ome.hcl"), "Cl.COC(=O)[C@@H](N)C(C)C"),
    (("serine methyl ester hydrochloride", "h-ser-ome.hcl"), "Cl.COC(=O)[C@@H](N)CO"),
    (("glutamic acid diethyl ester hydrochloride",), "Cl.CCOC(=O)CC[C@H](N)C(=O)OCC"),
    (("boc-phe-oh", "boc-l-phenylalanine", "n-boc-phenylalanine"), "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O"),
    (("boc-leu-oh", "boc-l-leucine"), "CC(C)(C)OC(=O)N[C@@H](CC(C)C)C(=O)O"),
    (("boc-val-oh", "boc-l-valine"), "CC(C)(C)OC(=O)N[C@@H](C(C)C)C(=O)O"),
    (("boc-ser-oh", "boc-l-serine"), "CC(C)(C)OC(=O)N[C@@H](CO)C(=O)O"),
    (("boc-lys(z)-oh", "boc-lys(cbz)-oh"), "CC(C)(C)OC(=O)N[C@@H](CCCCNC(=O)OCc1ccccc1)C(=O)O"),
    (("boc-cys(trt)-oh",), "CC(C)(C)OC(=O)N[C@@H](CSC(c1ccccc1)(c1ccccc1)c1ccccc1)C(=O)O"),
    (("boc-trp-oh", "boc-l-tryptophan"), "CC(C)(C)OC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O"),
    (("boc-tyr-oh", "boc-l-tyrosine"), "CC(C)(C)OC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)O"),
    (("boc-asp(obzl)-oh",), "CC(C)(C)OC(=O)N[C@@H](CC(=O)OCc1ccccc1)C(=O)O"),
    (("boc-glu(obzl)-oh",), "CC(C)(C)OC(=O)N[C@@H](CCC(=O)OCc1ccccc1)C(=O)O"),
    (("boc-beta-alanine", "boc-β-ala-oh", "boc-beta-ala-oh", "3-(boc-amino)propionic acid"),
     "CC(C)(C)OC(=O)NCCC(=O)O"),
    (("boc-gaba-oh", "4-(boc-amino)butyric acid"), "CC(C)(C)OC(=O)NCCCC(=O)O"),
    (("boc-ahx-oh", "6-(boc-amino)hexanoic acid", "boc-6-aminohexanoic acid"), "CC(C)(C)OC(=O)NCCCCCC(=O)O"),
    (("fmoc-beta-alanine", "fmoc-β-ala-oh", "fmoc-beta-ala-oh"), "OC(=O)CCNC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-ahx-oh", "fmoc-6-aminohexanoic acid"), "OC(=O)CCCCCNC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-lys(mtt)-oh",), "Cc1ccc(cc1)C(NCCCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O)(c1ccccc1)c1ccccc1"),
    (("fmoc-lys(alloc)-oh",), "C=CCOC(=O)NCCCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-lys(dde)-oh",), "CC(=NCCCC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O)C1=C(O)CC(C)(C)CC1=O"),
    (("fmoc-glu(oall)-oh",), "C=CCOC(=O)CC[C@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-d-ala-oh", "fmoc-d-alanine"), "C[C@@H](NC(=O)OCC1c2ccccc2-c2ccccc12)C(=O)O"),
    (("fmoc-d-phe-oh", "fmoc-d-phenylalanine"), "OC(=O)[C@@H](Cc1ccccc1)NC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-sar-oh", "fmoc-sarcosine", "fmoc-n-methylglycine"), "CN(CC(=O)O)C(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-hyp(tbu)-oh",), "CC(C)(C)O[C@H]1CN(C(=O)OCC2c3ccccc3-c3ccccc23)[C@@H](C1)C(=O)O"),
    (("fmoc-peg2-oh", "fmoc-aeea", "fmoc-8-amino-3,6-dioxaoctanoic acid", "fmoc-o2oc-oh"),
     "OC(=O)COCCOCCNC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("fmoc-pna-t", "fmoc-pna-t-oh"), "Cc1cn(CC(=O)N(CCNC(=O)OCC2c3ccccc3-c3ccccc23)CC(=O)O)c(=O)[nH]c1=O"),
    (("z-gly-oh", "cbz-glycine", "n-cbz-glycine", "cbz-gly-oh"), "OC(=O)CNC(=O)OCc1ccccc1"),
    (("z-phe-oh", "cbz-l-phenylalanine", "cbz-phe-oh"), "OC(=O)[C@H](Cc1ccccc1)NC(=O)OCc1ccccc1"),
    (("z-ala-oh", "cbz-l-alanine", "cbz-ala-oh"), "C[C@H](NC(=O)OCc1ccccc1)C(=O)O"),
    (("z-pro-oh", "cbz-l-proline"), "OC(=O)[C@@H]1CCCN1C(=O)OCc1ccccc1"),
    (("z-lys(boc)-oh",), "CC(C)(C)OC(=O)NCCCC[C@H](NC(=O)OCc1ccccc1)C(=O)O"),
    (("ac-gly-nh2", "n-acetylglycinamide"), "CC(=O)NCC(N)=O"),
    (("glycylglycine", "gly-gly", "glygly", "diglycine"), "NCC(=O)NCC(=O)O"),
    (("triglycine", "gly-gly-gly", "glycylglycylglycine", "ggg"), "NCC(=O)NCC(=O)NCC(=O)O"),
    (("alanylalanine", "ala-ala", "l-alanyl-l-alanine"), "C[C@H](N)C(=O)N[C@@H](C)C(=O)O"),
    (("aspartame", "l-aspartyl-l-phenylalanine methyl ester", "asp-phe-ome"),
     "COC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CC(=O)O"),
    (("carnosine", "beta-alanyl-l-histidine", "β-alanyl-l-histidine"), "NCCC(=O)N[C@@H](Cc1c[nH]cn1)C(=O)O"),
    (("glutathione", "gsh", "reduced glutathione", "γ-glu-cys-gly", "gamma-glu-cys-gly", "l-glutathione"),
     "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O"),
    (("glutathione disulfide", "gssg", "oxidized glutathione", "oxidised glutathione"),
     "N[C@@H](CCC(=O)N[C@@H](CSSC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O)C(=O)NCC(=O)O)C(=O)O"),
    (("rgd", "arg-gly-asp", "arginylglycylaspartic acid"),
     "N=C(N)NCCC[C@H](N)C(=O)NCC(=O)N[C@@H](CC(=O)O)C(=O)O"),
    (("leu-enkephalin", "leucine enkephalin", "[leu5]-enkephalin", "tyr-gly-gly-phe-leu"),
     "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)CNC(=O)CNC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)O"),
    (("met-enkephalin", "methionine enkephalin", "[met5]-enkephalin", "tyr-gly-gly-phe-met"),
     "CSCC[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)CNC(=O)CNC(=O)[C@@H](N)Cc1ccc(O)cc1)C(=O)O"),
    (("oxytocin",),
     "CC[C@H](C)[C@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@@H](C(=O)N2CCC[C@H]2C(=O)N[C@@H](CC(C)C)C(=O)NCC(N)=O)NC(=O)[C@H](CC(N)=O)NC(=O)[C@H](CCC(N)=O)NC1=O"),
    (("thyrotropin-releasing hormone", "trh", "pglu-his-pro-nh2", "protirelin"),
     "NC(=O)[C@@H]1CCCN1C(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@@H]1CCC(=O)N1"),
    (("pyroglutamic acid", "l-pyroglutamic acid", "pglu", "5-oxoproline", "pca"), "OC(=O)[C@@H]1CCC(=O)N1"),
    (("kainic acid", "kainate"), "C=C(C)[C@@H]1CN[C@@H]([C@H]1CC(=O)O)C(=O)O"),
    (("theanine", "l-theanine", "n-ethyl-l-glutamine"), "CCNC(=O)CC[C@H](N)C(=O)O"),
    (("ectoine",), "CC1=N[C@@H](CCN1)C(=O)O"),
    (("hypotaurine",), "NCCS(=O)O"),
    (("canavanine", "l-canavanine"), "N=C(N)NOCC[C@H](N)C(=O)O"),
    (("mimosine", "l-mimosine"), "N[C@@H](Cn1ccc(=O)c(O)c1)C(=O)O"),
    (("statine", "(3s,4s)-statine", "sta"), "CC(C)C[C@H](N)[C@@H](O)CC(=O)O"),
    (("4-aminomethylbenzoic acid", "pamba", "4-(aminomethyl)benzoic acid"), "NCc1ccc(cc1)C(=O)O"),
    (("tranexamic acid", "trans-4-(aminomethyl)cyclohexanecarboxylic acid", "txa"),
     "NC[C@H]1CC[C@@H](CC1)C(=O)O"),
    (("iminodiacetic acid", "ida"), "OC(=O)CNCC(=O)O"),
    (("nitrilotriacetic acid", "nta"), "OC(=O)CN(CC(=O)O)CC(=O)O"),
    (("aminomalonic acid",), "NC(C(=O)O)C(=O)O"),
    # ------------------------------------------------------------------ #
    # Sugars (pyranose forms where common; D-glucose = beta-D-glucopyranose)
    # ------------------------------------------------------------------ #
    (("glucose", "d-glucose", "dextrose", "beta-d-glucose", "β-d-glucose", "beta-d-glucopyranose", "β-d-glucopyranose", "grape sugar", "blood sugar", "glc"),
     "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("alpha-d-glucose", "α-d-glucose", "alpha-d-glucopyranose", "α-d-glucopyranose"),
     "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("d-glucose (open chain)", "open-chain d-glucose", "aldehydo-d-glucose"),
     "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"),
    (("l-glucose",), "OC[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),
    (("fructose", "d-fructose", "levulose", "fruit sugar", "beta-d-fructofuranose", "β-d-fructofuranose", "fru"),
     "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O"),
    (("galactose", "d-galactose", "beta-d-galactose", "β-d-galactose", "beta-d-galactopyranose", "gal"),
     "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"),
    (("mannose", "d-mannose", "beta-d-mannose", "β-d-mannose", "beta-d-mannopyranose", "man"),
     "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"),
    (("alpha-d-mannose", "α-d-mannose", "alpha-d-mannopyranose"),
     "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O"),
    (("allose", "d-allose", "beta-d-allose"), "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O"),
    (("talose", "d-talose", "beta-d-talose"), "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"),
    (("ribose", "d-ribose", "beta-d-ribofuranose", "β-d-ribofuranose", "rib"),
     "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O"),
    (("2-deoxyribose", "deoxyribose", "2-deoxy-d-ribose", "d-2-deoxyribose", "beta-d-2-deoxyribofuranose"),
     "OC[C@H]1O[C@@H](O)C[C@@H]1O"),
    (("arabinose", "l-arabinose", "ara", "l-arabinose (open chain)"), "O=C[C@H](O)[C@@H](O)[C@@H](O)CO"),
    (("d-arabinose",), "O=C[C@@H](O)[C@H](O)[C@H](O)CO"),
    (("xylose", "d-xylose", "wood sugar", "xyl", "d-xylose (open chain)"), "O=C[C@H](O)[C@@H](O)[C@H](O)CO"),
    (("xylitol",), "OCC(O)C(O)C(O)CO"),
    (("lyxose", "d-lyxose"), "O=C[C@@H](O)[C@@H](O)[C@H](O)CO"),
    (("rhamnose", "l-rhamnose", "6-deoxy-l-mannose", "alpha-l-rhamnopyranose", "rha"),
     "C[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O"),
    (("fucose", "l-fucose", "6-deoxy-l-galactose", "alpha-l-fucopyranose", "fuc"),
     "C[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),
    (("erythrose", "d-erythrose"), "O=C[C@H](O)[C@H](O)CO"),
    (("threose", "d-threose"), "O=C[C@@H](O)[C@H](O)CO"),
    (("glyceraldehyde", "d-glyceraldehyde", "(r)-glyceraldehyde"), "O=C[C@H](O)CO"),
    (("l-glyceraldehyde", "(s)-glyceraldehyde"), "O=C[C@@H](O)CO"),
    (("dihydroxyacetone", "1,3-dihydroxyacetone"), "OCC(=O)CO"),
    (("sucrose", "table sugar", "saccharose", "cane sugar", "beet sugar"),
     "OC[C@H]1O[C@@](CO)(O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O"),
    (("lactose", "milk sugar", "d-lactose", "beta-lactose", "β-lactose"),
     "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@H]1O"),
    (("maltose", "malt sugar", "d-maltose"),
     "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("cellobiose", "d-cellobiose"),
     "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("trehalose", "alpha,alpha-trehalose", "α,α-trehalose", "d-trehalose", "mycose"),
     "OC[C@H]1O[C@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("isomaltose",), "OC[C@H]1O[C@H](OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("gentiobiose",),
     "OC[C@H]1O[C@@H](OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("melibiose",), "OC[C@H]1O[C@H](OC[C@H]2O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
    (("raffinose",),
     "OC[C@H]1O[C@H](OC[C@H]2O[C@H](O[C@]3(CO)O[C@H](CO)[C@@H](O)[C@@H]3O)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@H]1O"),
    (("n-acetylglucosamine", "glcnac", "n-acetyl-d-glucosamine", "nag", "2-acetamido-2-deoxy-d-glucose", "beta-d-glcnac"),
     "CC(=O)N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"),
    (("n-acetylgalactosamine", "galnac", "n-acetyl-d-galactosamine"),
     "CC(=O)N[C@@H]1[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]1O"),
    (("glucosamine", "d-glucosamine", "2-amino-2-deoxy-d-glucose", "glcn"),
     "N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"),
    (("glucosamine hydrochloride", "d-glucosamine hydrochloride"),
     "Cl.N[C@@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O"),
    (("galactosamine", "d-galactosamine"), "N[C@@H]1[C@@H](O)[C@@H](O)[C@@H](CO)O[C@H]1O"),
    (("sialic acid", "n-acetylneuraminic acid", "neu5ac", "neuac", "nana", "neuraminic acid n-acetyl"),
     "CC(=O)N[C@@H]1[C@@H](O)C[C@@](O)(O[C@H]1[C@H](O)[C@H](O)CO)C(=O)O"),
    (("kdo", "3-deoxy-d-manno-oct-2-ulosonic acid"), "OC[C@@H](O)[C@@H]1O[C@@](O)(C[C@@H](O)[C@@H]1O)C(=O)O"),
    (("muramic acid", "n-acetylmuramic acid", "murnac", "nam"),
     "CC(=O)N[C@@H]1[C@@H](O[C@@H](C)C(=O)O)[C@H](O)[C@@H](CO)O[C@H]1O"),
    (("beta-d-glucuronic acid", "β-d-glucuronic acid", "glca"),
     "O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(=O)O"),
    (("galacturonic acid", "d-galacturonic acid", "gala"), "O[C@@H]1O[C@@H]([C@H](O)[C@H](O)[C@H]1O)C(=O)O"),
    (("sorbitol", "d-sorbitol", "glucitol", "d-glucitol"), "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"),
    (("mannitol", "d-mannitol"), "OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"),
    (("erythritol",), "OC[C@@H](O)[C@H](O)CO"),
    (("inositol", "myo-inositol", "cyclohexane-1,2,3,4,5,6-hexol"),
     "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("ascorbic acid", "vitamin c", "l-ascorbic acid", "ascorbate"), "OC[C@H](O)[C@H]1OC(=O)C(O)=C1O"),
    (("dehydroascorbic acid",), "OC[C@H](O)[C@H]1OC(=O)C(=O)C1=O"),
    (("sodium ascorbate", "naasc"), "[Na+].OC[C@H](O)[C@H]1OC(=O)C([O-])=C1O"),
    (("gluconolactone", "d-glucono-1,5-lactone", "glucono delta-lactone", "gdl"),
     "OC[C@H]1OC(=O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("methyl alpha-d-glucopyranoside", "methyl glucoside", "methyl α-d-glucopyranoside"),
     "CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),
    (("methyl beta-d-glucopyranoside", "methyl β-d-glucopyranoside"),
     "CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),
    (("glucose pentaacetate", "beta-d-glucose pentaacetate", "1,2,3,4,6-penta-o-acetyl-beta-d-glucopyranose"),
     "CC(=O)OC[C@H]1O[C@@H](OC(C)=O)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@H]1OC(C)=O"),
    (("acetobromoglucose", "2,3,4,6-tetra-o-acetyl-alpha-d-glucopyranosyl bromide"),
     "CC(=O)OC[C@H]1O[C@H](Br)[C@H](OC(C)=O)[C@@H](OC(C)=O)[C@@H]1OC(C)=O"),
    (("diacetone glucose", "dag", "1,2:5,6-di-o-isopropylidene-alpha-d-glucofuranose", "1,2:5,6-di-o-isopropylidene-α-d-glucofuranose"),
     "CC1(C)OC[C@H](O1)[C@H]1O[C@@H]2OC(C)(C)O[C@@H]2[C@H]1O"),
    (("levoglucosan", "1,6-anhydro-beta-d-glucopyranose"), "O[C@@H]1[C@@H](O)[C@H]2CO[C@@H](O2)[C@@H]1O"),
    (("glucose 6-phosphate", "g6p", "d-glucose 6-phosphate", "glucose-6-phosphate"),
     "O[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H](O)[C@H]1O"),
    (("glucose 1-phosphate", "g1p", "alpha-d-glucose 1-phosphate"),
     "OC[C@H]1O[C@H](OP(=O)(O)O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("fructose 6-phosphate", "f6p"), "OC[C@@]1(O)O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O"),
    (("fructose 1,6-bisphosphate", "fbp", "f1,6bp", "fructose-1,6-bisphosphate"),
     "O[C@@]1(COP(=O)(O)O)O[C@H](COP(=O)(O)O)[C@@H](O)[C@@H]1O"),
    (("ribose 5-phosphate", "r5p"), "O[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("prpp", "phosphoribosyl pyrophosphate", "5-phospho-alpha-d-ribose 1-diphosphate"),
     "O[C@H]1[C@@H](O)[C@H](OP(=O)(O)OP(=O)(O)O)O[C@@H]1COP(=O)(O)O"),
    (("glyceraldehyde 3-phosphate", "g3p", "gap", "d-glyceraldehyde 3-phosphate"), "O=C[C@H](O)COP(=O)(O)O"),
    (("dihydroxyacetone phosphate", "dhap"), "OCC(=O)COP(=O)(O)O"),
    (("phosphoenolpyruvate", "pep", "phosphoenolpyruvic acid"), "C=C(OP(=O)(O)O)C(=O)O"),
    (("3-phosphoglycerate", "3-pg", "3-phosphoglyceric acid"), "O[C@H](COP(=O)(O)O)C(=O)O"),
    (("2-phosphoglycerate", "2-pg"), "OC[C@@H](OP(=O)(O)O)C(=O)O"),
    (("glycerol 3-phosphate", "sn-glycerol 3-phosphate", "g3p (glycerol)", "glycerophosphate"),
     "OC[C@@H](O)COP(=O)(O)O"),
    (("udp-glucose", "uridine diphosphate glucose", "udp-glc"),
     "OC[C@H]1O[C@H](OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3ccc(=O)[nH]c3=O)[C@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("gdp-mannose",),
     "OC[C@H]1O[C@H](OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c3nc(N)[nH]c4=O)[C@H](O)[C@@H]2O)[C@@H](O)[C@@H](O)[C@@H]1O"),
    (("cyclodextrin", "beta-cyclodextrin", "β-cyclodextrin", "bcd"),
     "OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@H]4[C@H](O)[C@@H](O)[C@@H](O[C@H]5[C@H](O)[C@@H](O)[C@@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O[C@H]7[C@H](O)[C@@H](O)[C@@H](O[C@H]8[C@H](O)[C@@H](O)[C@@H](O[C@H]1[C@H](O)[C@H]2O)O[C@@H]8CO)O[C@@H]7CO)O[C@@H]6CO)O[C@@H]5CO)O[C@@H]4CO)O[C@@H]3CO"),
    # ------------------------------------------------------------------ #
    # Nucleobases, nucleosides, nucleotides and cofactors
    # ------------------------------------------------------------------ #
    (("adenine", "ade", "6-aminopurine", "vitamin b4"), "Nc1ncnc2[nH]cnc12"),
    (("guanine", "gua", "2-amino-6-oxopurine"), "Nc1nc2[nH]cnc2c(=O)[nH]1"),
    (("cytosine", "cyt"), "Nc1cc[nH]c(=O)n1"),
    (("thymine", "thy", "5-methyluracil"), "Cc1c[nH]c(=O)[nH]c1=O"),
    (("uracil", "ura"), "O=c1cc[nH]c(=O)[nH]1"),
    (("hypoxanthine", "6-oxopurine"), "O=c1[nH]cnc2[nH]cnc12"),
    (("xanthine",), "O=c1[nH]c(=O)c2[nH]cnc2[nH]1"),
    (("5-methylcytosine", "5mc"), "Cc1c[nH]c(=O)nc1N"),
    (("5-fluorouracil", "5-fu", "fluorouracil"), "O=c1[nH]cc(F)c(=O)[nH]1"),
    (("6-mercaptopurine", "6-mp", "mercaptopurine"), "S=c1[nH]cnc2[nH]cnc12"),
    (("6-thioguanine", "thioguanine", "6-tg"), "Nc1nc(=S)c2[nH]cnc2[nH]1"),
    (("5-azacytidine", "azacitidine", "5-aza-c"), "Nc1ncn(c(=O)n1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("adenosine", "ado", "a (nucleoside)"), "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("guanosine", "guo"), "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("cytidine", "cyd"), "Nc1ccn(c(=O)n1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("uridine", "urd"), "OC[C@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H]1O"),
    (("thymidine", "deoxythymidine", "2'-deoxythymidine", "dthd", "dt"),
     "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"),
    (("inosine", "ino"), "O=c1[nH]cnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("deoxyadenosine", "2'-deoxyadenosine", "dado"), "Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](CO)O1"),
    (("deoxyguanosine", "2'-deoxyguanosine", "dguo", "dg"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@H]1C[C@H](O)[C@@H](CO)O1"),
    (("deoxycytidine", "2'-deoxycytidine", "dcyd", "dc"), "Nc1ccn(c(=O)n1)[C@H]1C[C@H](O)[C@@H](CO)O1"),
    (("deoxyuridine", "2'-deoxyuridine", "durd", "du"), "OC[C@H]1O[C@@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"),
    (("bromodeoxyuridine", "brdu", "5-bromo-2'-deoxyuridine"),
     "OC[C@H]1O[C@@H](C[C@@H]1O)n1cc(Br)c(=O)[nH]c1=O"),
    (("edu", "5-ethynyl-2'-deoxyuridine"), "C#Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O"),
    (("amp", "adenosine monophosphate", "adenosine 5'-monophosphate", "adenylic acid", "5'-amp"),
     "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("adp", "adenosine diphosphate", "adenosine 5'-diphosphate"),
     "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("atp", "adenosine triphosphate", "adenosine 5'-triphosphate"),
     "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("atp disodium salt", "na2atp", "atp disodium", "adenosine 5'-triphosphate disodium salt"),
     "[Na+].[Na+].Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OP(=O)([O-])O)[C@@H](O)[C@H]1O"),
    (("gmp", "guanosine monophosphate", "guanosine 5'-monophosphate", "guanylic acid"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("gdp", "guanosine diphosphate", "guanosine 5'-diphosphate"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("gtp", "guanosine triphosphate", "guanosine 5'-triphosphate"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("cmp", "cytidine monophosphate", "cytidine 5'-monophosphate"),
     "Nc1ccn(c(=O)n1)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("ctp", "cytidine triphosphate", "cytidine 5'-triphosphate"),
     "Nc1ccn(c(=O)n1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"),
    (("ump", "uridine monophosphate", "uridine 5'-monophosphate", "uridylic acid"),
     "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1"),
    (("udp", "uridine diphosphate", "uridine 5'-diphosphate"),
     "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1"),
    (("utp", "uridine triphosphate", "uridine 5'-triphosphate"),
     "O=c1ccn([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]2O)c(=O)[nH]1"),
    (("ttp", "dttp", "thymidine triphosphate", "thymidine 5'-triphosphate", "deoxythymidine triphosphate"),
     "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O2)c(=O)[nH]c1=O"),
    (("datp", "deoxyadenosine triphosphate", "2'-deoxyadenosine 5'-triphosphate"),
     "Nc1ncnc2n(cnc12)[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O1"),
    (("dgtp", "deoxyguanosine triphosphate"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O1"),
    (("dctp", "deoxycytidine triphosphate"),
     "Nc1ccn(c(=O)n1)[C@H]1C[C@H](O)[C@@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O1"),
    (("camp", "cyclic amp", "3',5'-cyclic amp", "adenosine 3',5'-cyclic monophosphate", "cyclic adenosine monophosphate"),
     "Nc1ncnc2n(cnc12)[C@@H]1O[C@@H]2COP(=O)(O)O[C@H]2[C@H]1O"),
    (("cgmp", "cyclic gmp", "guanosine 3',5'-cyclic monophosphate"),
     "Nc1nc2n(cnc2c(=O)[nH]1)[C@@H]1O[C@@H]2COP(=O)(O)O[C@H]2[C@H]1O"),
    (("nad+", "nad", "nicotinamide adenine dinucleotide", "beta-nad", "β-nad", "oxidized nad", "nad(+)"),
     "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O"),
    (("nadh", "reduced nad", "1,4-dihydronicotinamide adenine dinucleotide"),
     "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O"),
    (("nadp+", "nadp", "nicotinamide adenine dinucleotide phosphate", "nadp(+)"),
     "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@H](OP(=O)(O)O)[C@@H]2O)[C@@H](O)[C@H]1O"),
    (("nadph", "reduced nadp"),
     "NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc34)[C@H](OP(=O)(O)O)[C@@H]2O)[C@@H](O)[C@H]1O"),
    (("fad", "flavin adenine dinucleotide"),
     "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc45)[C@H](O)[C@@H]3O)c2cc1C"),
    (("fmn", "flavin mononucleotide", "riboflavin 5'-phosphate"),
     "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)O)c2cc1C"),
    (("riboflavin", "vitamin b2"), "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C"),
    (("coenzyme a", "coa", "coash", "hs-coa"),
     "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)(O)O)[C@@H](O)C(=O)NCCC(=O)NCCS"),
    (("acetyl-coa", "acetyl coenzyme a", "acetyl coa", "ac-coa", "acetyl-coenzyme a"),
     "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)(O)O"),
    (("malonyl-coa",),
     "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)(O)O"),
    (("s-adenosylmethionine", "sam", "adomet", "s-adenosyl-l-methionine", "same"),
     "C[S+](CC[C@H](N)C(=O)[O-])C[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1O"),
    (("s-adenosylhomocysteine", "sah", "adohcy", "s-adenosyl-l-homocysteine"),
     "N[C@@H](CCSC[C@H]1O[C@@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1O)C(=O)O"),
    (("pantothenic acid", "vitamin b5", "d-pantothenic acid", "pantothenate"),
     "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)O"),
    (("calcium pantothenate",),
     "[Ca+2].CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)[O-].CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)[O-]"),
    (("biotin", "vitamin b7", "vitamin h", "d-biotin", "(+)-biotin"),
     "O=C1N[C@@H]2[C@H](CCCCC(=O)O)SC[C@@H]2N1"),
    (("biotin-nhs", "nhs-biotin", "biotin n-hydroxysuccinimide ester", "biotin nhs ester", "biotin-osu"),
     "O=C1CCC(=O)N1OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("sulfo-nhs-biotin", "biotin sulfo-nhs ester"),
     "[Na+].O=C1CC(S(=O)(=O)[O-])C(=O)N1OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("nhs-lc-biotin", "biotin-lc-nhs", "biotinamidohexanoic acid n-hydroxysuccinimide ester", "biotin-x-nhs"),
     "O=C1CCC(=O)N1OC(=O)CCCCCNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("biotin-peg3-azide", "biotin-peg3-n3"),
     "[N-]=[N+]=NCCOCCOCCOCCNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("biotin-peg4-alkyne", "biotin-peg4-propargyl"),
     "C#CCNC(=O)CCOCCOCCOCCOCCNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("biotin hydrazide",), "NNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("biocytin", "epsilon-n-biotinyl-l-lysine", "ε-n-biotinyl-l-lysine"),
     "N[C@@H](CCCCNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12)C(=O)O"),
    (("desthiobiotin", "dethiobiotin"), "C[C@@H]1NC(=O)N[C@@H]1CCCCCC(=O)O"),
    (("iminobiotin", "2-iminobiotin"), "N=C1N[C@@H]2[C@H](CCCCC(=O)O)SC[C@@H]2N1"),
    (("folic acid", "folate", "vitamin b9", "pteroylglutamic acid", "pteroyl-l-glutamic acid"),
     "Nc1nc2ncc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)nc2c(=O)[nH]1"),
    (("tetrahydrofolate", "tetrahydrofolic acid", "5,6,7,8-tetrahydrofolate"),
     "Nc1nc2NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)Nc2c(=O)[nH]1"),
    (("5-methyltetrahydrofolate", "5-mthf", "levomefolic acid", "methylfolate", "l-methylfolate"),
     "CN1[C@@H](CNc2ccc(cc2)C(=O)N[C@@H](CCC(=O)O)C(=O)O)CNc2nc(N)[nH]c(=O)c12"),
    (("leucovorin", "folinic acid", "5-formyltetrahydrofolate", "calcium folinate"),
     "Nc1nc2NC[C@H](CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)N(C=O)c2c(=O)[nH]1"),
    (("methotrexate", "mtx", "amethopterin"),
     "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(cc1)C(=O)N[C@@H](CCC(=O)O)C(=O)O"),
    (("pyridoxal", "vitamin b6 aldehyde"), "Cc1ncc(CO)c(C=O)c1O"),
    (("pyridoxal 5'-phosphate", "plp", "pyridoxal phosphate", "pyridoxal-5-phosphate", "p5p"),
     "Cc1ncc(COP(=O)(O)O)c(C=O)c1O"),
    (("pyridoxine", "vitamin b6", "pyridoxol"), "Cc1ncc(CO)c(CO)c1O"),
    (("pyridoxine hydrochloride", "vitamin b6 hydrochloride"), "Cl.Cc1ncc(CO)c(CO)c1O"),
    (("pyridoxamine",), "Cc1ncc(CO)c(CN)c1O"),
    (("thiamine", "vitamin b1", "aneurine"), "Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"),
    (("thiamine hydrochloride", "vitamin b1 hydrochloride", "thiamine hcl"),
     "Cl.[Cl-].Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1"),
    (("thiamine pyrophosphate", "tpp", "thiamine diphosphate", "thdp", "cocarboxylase"),
     "Cc1ncc(C[n+]2csc(CCOP(=O)([O-])OP(=O)(O)O)c2C)c(N)n1"),
    (("nicotinamide", "niacinamide", "vitamin b3 amide", "nicotinic acid amide"), "NC(=O)c1cccnc1"),
    (("nicotinamide riboside", "nr"), "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("nicotinamide mononucleotide", "nmn", "β-nmn", "beta-nmn"),
     "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])O)[C@@H](O)[C@H]1O"),
    (("lipoamide", "alpha-lipoamide", "thioctic acid amide"), "NC(=O)CCCCC1CCSS1"),
    (("coenzyme q10", "coq10", "ubiquinone", "ubidecarenone", "q10"),
     "COC1=C(OC)C(=O)C(C/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)=C(C)C1=O"),
    (("coenzyme q0", "coq0", "2,3-dimethoxy-5-methyl-1,4-benzoquinone"), "COC1=C(OC)C(=O)C=C(C)C1=O"),
    (("tetrahydrobiopterin", "bh4", "sapropterin", "(6r)-tetrahydrobiopterin"),
     "C[C@@H](O)[C@@H](O)[C@H]1CNc2nc(N)[nH]c(=O)c2N1"),
    (("biopterin",), "C[C@@H](O)[C@@H](O)c1cnc2nc(N)[nH]c(=O)c2n1"),
    (("pqq", "pyrroloquinoline quinone", "methoxatin"),
     "OC(=O)c1cc2C(=O)C(=O)c3[nH]c(cc3-c2nc1C(=O)O)C(=O)O"),
    (("menaquinone-4", "vitamin k2", "mk-4"),
     "CC1=C(C/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)C(=O)c2ccccc2C1=O"),
    (("phylloquinone", "vitamin k1", "phytonadione", "phytomenadione"),
     "CC1=C(C/C=C(\\C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C(=O)c2ccccc2C1=O"),
    (("retinol", "vitamin a", "all-trans-retinol"), "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/CO)C(C)(C)CCC1"),
    (("retinal", "retinaldehyde", "all-trans-retinal", "vitamin a aldehyde"),
     "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=O)C(C)(C)CCC1"),
    (("retinoic acid", "all-trans-retinoic acid", "tretinoin", "atra"),
     "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C(=O)O)C(C)(C)CCC1"),
    (("isotretinoin", "13-cis-retinoic acid"), "CC1=C(/C=C/C(C)=C/C=C/C(C)=C\\C(=O)O)C(C)(C)CCC1"),
    (("beta-carotene", "β-carotene", "b-carotene", "carotene", "provitamin a"),
     "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C)C(C)(C)CCC1"),
    (("lycopene",), "CC(C)=CCC/C(C)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(C)/CCC=C(C)C"),
    (("alpha-tocopherol", "α-tocopherol", "vitamin e", "tocopherol", "d-alpha-tocopherol", "(+)-α-tocopherol"),
     "Cc1c(C)c2O[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)CCc2c(C)c1O"),
    (("tocopherol acetate", "vitamin e acetate", "alpha-tocopheryl acetate"),
     "CC(=O)Oc1c(C)c(C)c2O[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)CCc2c1C"),
    (("trolox", "6-hydroxy-2,5,7,8-tetramethylchroman-2-carboxylic acid"),
     "Cc1c(C)c2OC(C)(CCc2c(C)c1O)C(=O)O"),
    (("cholecalciferol", "vitamin d3", "calciol"),
     "C=C1CC[C@@H](O)C/C1=C/C=C1\\CCC[C@]2(C)[C@@H]1CC[C@@H]2[C@H](C)CCCC(C)C"),
    (("ergocalciferol", "vitamin d2"),
     "C=C1CC[C@@H](O)C/C1=C/C=C1\\CCC[C@]2(C)[C@@H]1CC[C@@H]2[C@H](C)/C=C/[C@H](C)C(C)C"),
    (("calcitriol", "1,25-dihydroxyvitamin d3", "1α,25-dihydroxycholecalciferol", "1,25-(oh)2d3"),
     "C=C1[C@@H](O)C[C@@H](O)C/C1=C/C=C1\\CCC[C@]2(C)[C@@H]1CC[C@@H]2[C@H](C)CCCC(C)(C)O"),
    (("heme b", "heme", "haem", "protoheme ix", "iron protoporphyrin ix"),
     "CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C)CCC(=O)O)CCC(=O)O.[Fe+2]"),
    (("protoporphyrin ix", "ppix"),
     "CC1=C(C2=CC3=C(C(=C(N3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1N2)C)C=C)C)C=C)C)CCC(=O)O)CCC(=O)O"),
    (("porphin", "porphine"), "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3"),
    (("tetraphenylporphyrin", "h2tpp", "meso-tetraphenylporphyrin"),
     "c1ccc(cc1)-c1c2ccc(n2)c(-c2ccccc2)c2ccc([nH]2)c(-c2ccccc2)c2ccc(n2)c(-c2ccccc2)c2ccc1[nH]2"),
    (("hemin", "chlorohemin", "ferriprotoporphyrin ix chloride"),
     "CC1=C(C2=CC3=C(C(=C([N-]3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C)CCC(=O)O)CCC(=O)O.[Fe+3].[Cl-]"),
    # ------------------------------------------------------------------ #
    # Lipids and fatty acids
    # ------------------------------------------------------------------ #
    (("lauric acid", "dodecanoic acid", "c12:0"), "CCCCCCCCCCCC(=O)O"),
    (("myristic acid", "tetradecanoic acid", "c14:0"), "CCCCCCCCCCCCCC(=O)O"),
    (("palmitic acid", "hexadecanoic acid", "c16:0", "palmitate"), "CCCCCCCCCCCCCCCC(=O)O"),
    (("palmitoleic acid", "(z)-hexadec-9-enoic acid", "c16:1"), "CCCCCC/C=C\\CCCCCCCC(=O)O"),
    (("stearic acid", "octadecanoic acid", "c18:0", "stearate"), "CCCCCCCCCCCCCCCCCC(=O)O"),
    (("oleic acid", "(z)-octadec-9-enoic acid", "c18:1", "oleate", "cis-9-octadecenoic acid"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)O"),
    (("elaidic acid", "(e)-octadec-9-enoic acid", "trans-9-octadecenoic acid"), "CCCCCCCC/C=C/CCCCCCCC(=O)O"),
    (("linoleic acid", "c18:2", "(9z,12z)-octadecadienoic acid", "linoleate"),
     "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O"),
    (("alpha-linolenic acid", "α-linolenic acid", "linolenic acid", "c18:3 n-3", "(9z,12z,15z)-octadecatrienoic acid"),
     "CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O"),
    (("gamma-linolenic acid", "γ-linolenic acid", "gla", "(6z,9z,12z)-octadecatrienoic acid"),
     "CCCCC/C=C\\C/C=C\\C/C=C\\CCCCC(=O)O"),
    (("arachidic acid", "eicosanoic acid", "c20:0"), "CCCCCCCCCCCCCCCCCCCC(=O)O"),
    (("arachidonic acid", "c20:4", "(5z,8z,11z,14z)-eicosatetraenoic acid", "arachidonate"),
     "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O"),
    (("eicosapentaenoic acid", "epa", "c20:5", "timnodonic acid"),
     "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O"),
    (("docosahexaenoic acid", "dha", "c22:6", "cervonic acid"),
     "CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC(=O)O"),
    (("behenic acid", "docosanoic acid", "c22:0"), "CCCCCCCCCCCCCCCCCCCCCC(=O)O"),
    (("erucic acid", "(z)-docos-13-enoic acid"), "CCCCCCCC/C=C\\CCCCCCCCCCCC(=O)O"),
    (("lignoceric acid", "tetracosanoic acid", "c24:0"), "CCCCCCCCCCCCCCCCCCCCCCCC(=O)O"),
    (("sodium stearate", "sodium octadecanoate"), "[Na+].CCCCCCCCCCCCCCCCCC(=O)[O-]"),
    (("sodium oleate",), "[Na+].CCCCCCCC/C=C\\CCCCCCCC(=O)[O-]"),
    (("methyl oleate",), "CCCCCCCC/C=C\\CCCCCCCC(=O)OC"),
    (("methyl palmitate",), "CCCCCCCCCCCCCCCC(=O)OC"),
    (("methyl stearate",), "CCCCCCCCCCCCCCCCCC(=O)OC"),
    (("ethyl oleate",), "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC"),
    (("palmitoyl chloride",), "CCCCCCCCCCCCCCCC(=O)Cl"),
    (("stearoyl chloride",), "CCCCCCCCCCCCCCCCCC(=O)Cl"),
    (("oleoyl chloride",), "CCCCCCCC/C=C\\CCCCCCCC(=O)Cl"),
    (("lauroyl chloride",), "CCCCCCCCCCCC(=O)Cl"),
    (("cetyl alcohol", "1-hexadecanol", "hexadecanol", "palmityl alcohol"), "CCCCCCCCCCCCCCCCO"),
    (("stearyl alcohol", "1-octadecanol", "octadecanol"), "CCCCCCCCCCCCCCCCCCO"),
    (("oleyl alcohol",), "CCCCCCCC/C=C\\CCCCCCCCO"),
    (("lauryl alcohol", "1-dodecanol", "dodecanol"), "CCCCCCCCCCCCO"),
    (("tripalmitin", "glyceryl tripalmitate"),
     "CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"),
    (("tristearin", "glyceryl tristearate"),
     "CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCCCC"),
    (("triolein", "glyceryl trioleate"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(COC(=O)CCCCCCC/C=C\\CCCCCCCC)OC(=O)CCCCCCC/C=C\\CCCCCCCC"),
    (("glyceryl monostearate", "monostearin", "1-monostearin", "gms"), "CCCCCCCCCCCCCCCCCC(=O)OCC(O)CO"),
    (("monoolein", "glyceryl monooleate", "1-monoolein"), "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(O)CO"),
    (("dppc", "dipalmitoylphosphatidylcholine", "1,2-dipalmitoyl-sn-glycero-3-phosphocholine"),
     "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC"),
    (("popc", "1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine"),
     "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"),
    (("dopc", "1,2-dioleoyl-sn-glycero-3-phosphocholine"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC"),
    (("dope", "1,2-dioleoyl-sn-glycero-3-phosphoethanolamine"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)(O)OCCN)OC(=O)CCCCCCC/C=C\\CCCCCCCC"),
    (("dmpc", "1,2-dimyristoyl-sn-glycero-3-phosphocholine"),
     "CCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC"),
    (("dspc", "1,2-distearoyl-sn-glycero-3-phosphocholine"),
     "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCCCC"),
    (("pops", "1-palmitoyl-2-oleoyl-sn-glycero-3-phospho-l-serine"),
     "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)OC[C@H](N)C(=O)O)OC(=O)CCCCCCC/C=C\\CCCCCCCC"),
    (("dotap", "1,2-dioleoyl-3-trimethylammonium-propane"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(C[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC.[Cl-]"),
    (("1,2-dimyristoyl-rac-glycerol",), "CCCCCCCCCCCCCC(=O)OCC(CO)OC(=O)CCCCCCCCCCCCC"),
    (("sphingosine", "d-erythro-sphingosine", "(2s,3r,4e)-2-aminooctadec-4-ene-1,3-diol"),
     "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](N)CO"),
    (("sphingosine-1-phosphate", "s1p"), "CCCCCCCCCCCCC/C=C/[C@@H](O)[C@@H](N)COP(=O)(O)O"),
    (("sphinganine", "dihydrosphingosine"), "CCCCCCCCCCCCCCC[C@@H](O)[C@@H](N)CO"),
    (("ceramide", "c16 ceramide", "n-palmitoyl-d-erythro-sphingosine", "c16-ceramide"),
     "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)/C=C/CCCCCCCCCCCCC"),
    (("c2 ceramide", "n-acetylsphingosine", "c2-ceramide"), "CC(=O)N[C@@H](CO)[C@H](O)/C=C/CCCCCCCCCCCCC"),
    (("sphingomyelin", "c16 sphingomyelin", "n-palmitoyl-sphingomyelin"),
     "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)[C@H](O)/C=C/CCCCCCCCCCCCC"),
    (("phosphocholine", "phosphorylcholine", "choline phosphate"), "C[N+](C)(C)CCOP(=O)([O-])O"),
    (("glycerophosphocholine", "sn-glycero-3-phosphocholine", "gpc", "alpha-gpc", "choline alfoscerate"),
     "C[N+](C)(C)CCOP(=O)([O-])OC[C@H](O)CO"),
    (("citicoline", "cdp-choline", "cytidine diphosphate choline"),
     "C[N+](C)(C)CCOP(=O)([O-])OP(=O)(O)OC[C@H]1O[C@@H](n2ccc(N)nc2=O)[C@H](O)[C@@H]1O"),
    (("choline",), "C[N+](C)(C)CCO"),
    (("choline chloride",), "C[N+](C)(C)CCO.[Cl-]"),
    (("acetylcholine", "ach"), "CC(=O)OCC[N+](C)(C)C"),
    (("acetylcholine chloride", "ach chloride"), "CC(=O)OCC[N+](C)(C)C.[Cl-]"),
    (("phytol",), "CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC/C(C)=C/CO"),
    (("squalene",), "CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(C)/CC/C=C(C)/CCC=C(C)C"),
    (("geraniol", "(e)-3,7-dimethylocta-2,6-dien-1-ol"), "CC(C)=CCC/C(C)=C/CO"),
    (("nerol", "(z)-3,7-dimethylocta-2,6-dien-1-ol"), "CC(C)=CCC/C(C)=C\\CO"),
    (("farnesol", "(2e,6e)-farnesol"), "CC(C)=CCC/C(C)=C/CC/C(C)=C/CO"),
    (("citral", "geranial", "(e)-citral"), "CC(C)=CCC/C(C)=C/C=O"),
    (("citronellol",), "CC(CCO)CCC=C(C)C"),
    (("linalool",), "CC(C)=CCCC(C)(O)C=C"),
    (("isopentenyl pyrophosphate", "ipp", "isopentenyl diphosphate"), "C=C(C)CCOP(=O)(O)OP(=O)(O)O"),
    (("dimethylallyl pyrophosphate", "dmapp", "dimethylallyl diphosphate"), "CC(C)=CCOP(=O)(O)OP(=O)(O)O"),
    (("geranyl pyrophosphate", "gpp", "geranyl diphosphate"), "CC(C)=CCC/C(C)=C/COP(=O)(O)OP(=O)(O)O"),
    (("farnesyl pyrophosphate", "fpp", "farnesyl diphosphate"),
     "CC(C)=CCC/C(C)=C/CC/C(C)=C/COP(=O)(O)OP(=O)(O)O"),
    # ------------------------------------------------------------------ #
    # Steroids, bile acids and eicosanoids
    # ------------------------------------------------------------------ #
    (("cholesterol", "cholest-5-en-3β-ol", "cholest-5-en-3beta-ol"),
     "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"),
    (("7-dehydrocholesterol", "provitamin d3"),
     "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C"),
    (("cholesteryl hemisuccinate", "chems"),
     "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)OC(=O)CCC(=O)O)C)C"),
    (("cholesteryl acetate",),
     "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)OC(C)=O)C)C"),
    (("ergosterol", "provitamin d2"),
     "C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C"),
    (("stigmasterol",),
     "CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C"),
    (("beta-sitosterol", "β-sitosterol", "sitosterol"),
     "CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C"),
    (("lanosterol",),
     "CC(C)=CCC[C@@H](C)[C@H]1CC[C@@]2(C)C1=CC[C@H]1[C@@]2(C)CC[C@@H]2C(C)(C)[C@@H](O)CC[C@@]12C"),
    (("testosterone", "17β-hydroxyandrost-4-en-3-one", "17beta-hydroxyandrost-4-en-3-one"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C"),
    (("testosterone propionate",), "CCC(=O)O[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C"),
    (("testosterone enanthate",), "CCCCCCC(=O)O[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C"),
    (("dihydrotestosterone", "dht", "androstanolone", "5α-dihydrotestosterone", "5alpha-dihydrotestosterone", "stanolone"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CC[C@@H]4[C@@]3(CCC(=O)C4)C"),
    (("androstenedione", "androst-4-ene-3,17-dione", "4-androstenedione"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CCC4=CC(=O)CC[C@]34C"),
    (("dehydroepiandrosterone", "dhea", "prasterone", "androstenolone"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CC=C4[C@@]3(CC[C@@H](C4)O)C"),
    (("nandrolone", "19-nortestosterone"), "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@H]34"),
    (("estradiol", "17β-estradiol", "17beta-estradiol", "oestradiol", "e2", "estra-1,3,5(10)-triene-3,17β-diol"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O"),
    (("estrone", "oestrone", "e1"), "C[C@]12CC[C@H]3[C@H]([C@@H]1CCC2=O)CCC4=C3C=CC(=C4)O"),
    (("estriol", "oestriol", "e3"), "C[C@]12CC[C@H]3[C@H]([C@@H]1C[C@@H]([C@@H]2O)O)CCC4=C3C=CC(=C4)O"),
    (("ethinylestradiol", "ethinyl estradiol", "ethynylestradiol", "ee2", "17α-ethinylestradiol"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)O"),
    (("progesterone", "pregn-4-ene-3,20-dione"),
     "CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C"),
    (("pregnenolone",), "CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"),
    (("17α-hydroxyprogesterone", "17-hydroxyprogesterone", "17-ohp", "17alpha-hydroxyprogesterone"),
     "CC(=O)[C@@]1(O)CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C"),
    (("cortisol", "hydrocortisone", "11β,17α,21-trihydroxypregn-4-ene-3,20-dione", "11beta,17alpha,21-trihydroxypregn-4-ene-3,20-dione"),
     "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O"),
    (("cortisone",), "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2C(=O)C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C"),
    (("corticosterone",), "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C)O"),
    (("aldosterone",), "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C=O)O"),
    (("11-deoxycortisol", "cortexolone", "cortodoxone"),
     "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)CC[C@]34C"),
    (("prednisolone",), "C[C@]12C=CC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O"),
    (("prednisone",), "C[C@]12C=CC(=O)C=C1CC[C@@H]3[C@@H]2C(=O)C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C"),
    (("methylprednisolone",),
     "C[C@H]1C[C@H]2[C@@H]3CC[C@@]([C@]3(C[C@@H]([C@@H]2[C@@]4(C)C=CC(=O)C=C14)O)C)(C(=O)CO)O"),
    (("dexamethasone",),
     "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C"),
    (("betamethasone",),
     "C[C@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)CO)O)C)O)F)C"),
    (("triamcinolone",),
     "C[C@]12C[C@@H]([C@]3([C@H]([C@@H]1C[C@@H]([C@]2(C(=O)CO)O)O)CCC4=CC(=O)C=C[C@@]43C)F)O"),
    (("fludrocortisone",), "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@]2(F)[C@H](C[C@]4([C@H]3CC[C@@]4(C(=O)CO)O)C)O"),
    (("budesonide",),
     "CCCC1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@@]5([C@H]4[C@H](C[C@@]3([C@@]2(O1)C(=O)CO)C)O)C"),
    (("fluticasone propionate", "flovent", "flonase"),
     "C[C@@H]1C[C@H]2[C@@H]3C[C@@H](C4=CC(=O)C=C[C@@]4([C@]3([C@H](C[C@@]2([C@]1(C(=O)SCF)OC(=O)CC)C)O)F)C)F"),
    (("cholic acid",),
     "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    (("deoxycholic acid",),
     "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"),
    (("sodium deoxycholate",),
     "[Na+].C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"),
    (("chenodeoxycholic acid", "cdca", "chenodiol"),
     "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"),
    (("ursodeoxycholic acid", "udca", "ursodiol"),
     "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2[C@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)C"),
    (("lithocholic acid", "lca"),
     "C[C@H](CCC(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@H]4[C@@]3(CC[C@H](C4)O)C)C"),
    (("taurocholic acid", "taurocholate"),
     "C[C@H](CCC(=O)NCCS(=O)(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    (("glycocholic acid", "glycocholate"),
     "C[C@H](CCC(=O)NCC(=O)O)[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    (("prostaglandin e2", "pge2", "dinoprostone"),
     "CCCCC[C@H](O)/C=C/[C@H]1[C@H](O)CC(=O)[C@@H]1C/C=C\\CCCC(=O)O"),
    (("prostaglandin f2alpha", "pgf2α", "pgf2a", "pgf2alpha", "dinoprost"),
     "CCCCC[C@H](O)/C=C/[C@H]1[C@H](O)C[C@H](O)[C@@H]1C/C=C\\CCCC(=O)O"),
    (("prostaglandin e1", "pge1", "alprostadil"),
     "CCCCC[C@H](O)/C=C/[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(=O)O"),
    (("prostaglandin d2", "pgd2"), "CCCCC[C@H](O)/C=C/[C@H]1[C@@H](C/C=C\\CCCC(=O)O)C(=O)C[C@H]1O"),
    (("thromboxane b2", "txb2"), "CCCCC[C@H](O)/C=C/[C@H]1OC(O)C[C@@H](O)[C@@H]1C/C=C\\CCCC(=O)O"),
    (("leukotriene b4", "ltb4"), "CCCCC/C=C\\CC(O)/C=C/C=C/C=C\\C(O)CCCC(=O)O"),
    (("anandamide", "aea", "n-arachidonoylethanolamine"), "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)NCCO"),
    (("2-arachidonoylglycerol", "2-ag"), "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)OC(CO)CO"),
    (("platelet-activating factor", "paf", "paf-c16"),
     "CCCCCCCCCCCCCCCCOC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(C)=O"),
    (("lysophosphatidic acid", "lpa", "18:1 lpa", "oleoyl-lpa"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](O)COP(=O)(O)O"),
    (("dipalmitoylphosphatidic acid", "1,2-dipalmitoyl-sn-glycero-3-phosphate"),
     "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)(O)O)OC(=O)CCCCCCCCCCCCCCC"),
    (("ethanolamine phosphate", "phosphoethanolamine", "o-phosphoethanolamine"), "NCCOP(=O)(O)O"),
    # ------------------------------------------------------------------ #
    # Neurotransmitters, metabolites and signalling molecules
    # ------------------------------------------------------------------ #
    (("melatonin", "n-acetyl-5-methoxytryptamine"), "COc1ccc2[nH]cc(CCNC(C)=O)c2c1"),
    (("n-acetylserotonin",), "CC(=O)NCCc1c[nH]c2ccc(O)cc12"),
    (("serotonin", "5-hydroxytryptamine", "5-ht", "enteramine"), "NCCc1c[nH]c2ccc(O)cc12"),
    (("serotonin hydrochloride", "5-ht hcl"), "Cl.NCCc1c[nH]c2ccc(O)cc12"),
    (("dopamine", "3,4-dihydroxyphenethylamine", "3-hydroxytyramine"), "NCCc1ccc(O)c(O)c1"),
    (("dopamine hydrochloride",), "Cl.NCCc1ccc(O)c(O)c1"),
    (("norepinephrine", "noradrenaline", "(r)-norepinephrine", "l-norepinephrine"),
     "NC[C@H](O)c1ccc(O)c(O)c1"),
    (("epinephrine", "adrenaline", "(r)-epinephrine", "l-epinephrine", "l-adrenaline", "(-)-adrenaline"),
     "CNC[C@H](O)c1ccc(O)c(O)c1"),
    (("tyramine", "4-hydroxyphenethylamine"), "NCCc1ccc(O)cc1"),
    (("octopamine", "p-octopamine", "norsynephrine"), "NCC(O)c1ccc(O)cc1"),
    (("synephrine", "p-synephrine"), "CNCC(O)c1ccc(O)cc1"),
    (("hordenine", "n,n-dimethyltyramine"), "CN(C)CCc1ccc(O)cc1"),
    (("histamine", "2-(1h-imidazol-4-yl)ethanamine"), "NCCc1c[nH]cn1"),
    (("histamine dihydrochloride",), "Cl.Cl.NCCc1c[nH]cn1"),
    (("n,n-dimethyltryptamine", "dmt"), "CN(C)CCc1c[nH]c2ccccc12"),
    (("psilocin", "4-hydroxy-dmt", "4-ho-dmt"), "CN(C)CCc1c[nH]c2cccc(O)c12"),
    (("psilocybin",), "CN(C)CCc1c[nH]c2cccc(OP(=O)(O)O)c12"),
    (("bufotenine", "5-hydroxy-dmt", "5-ho-dmt"), "CN(C)CCc1c[nH]c2ccc(O)cc12"),
    (("mescaline", "3,4,5-trimethoxyphenethylamine"), "COc1cc(CCN)cc(OC)c1OC"),
    (("amphetamine", "alpha-methylphenethylamine", "α-methylphenethylamine"), "CC(N)Cc1ccccc1"),
    (("dextroamphetamine", "d-amphetamine", "(s)-amphetamine", "dexamfetamine"), "C[C@H](N)Cc1ccccc1"),
    (("methamphetamine", "metamfetamine", "n-methylamphetamine"), "CNC(C)Cc1ccccc1"),
    (("mdma", "3,4-methylenedioxymethamphetamine", "ecstasy", "midomafetamine"), "CNC(C)Cc1ccc2OCOc2c1"),
    (("cathinone",), "CC(N)C(=O)c1ccccc1"),
    (("ephedrine", "(-)-ephedrine", "(1r,2s)-ephedrine", "l-ephedrine"), "CN[C@@H](C)[C@H](O)c1ccccc1"),
    (("pseudoephedrine", "(+)-pseudoephedrine", "(1s,2s)-pseudoephedrine"), "CN[C@@H](C)[C@@H](O)c1ccccc1"),
    (("phenylephrine", "(r)-phenylephrine"), "CNC[C@H](O)c1cccc(O)c1"),
    (("isoproterenol", "isoprenaline"), "CC(C)NCC(O)c1ccc(O)c(O)c1"),
    (("homovanillic acid", "hva"), "COc1cc(CC(=O)O)ccc1O"),
    (("vanillylmandelic acid", "vma"), "COc1cc(C(O)C(=O)O)ccc1O"),
    (("3,4-dihydroxyphenylacetic acid", "dopac"), "OC(=O)Cc1ccc(O)c(O)c1"),
    (("5-hydroxyindoleacetic acid", "5-hiaa"), "OC(=O)Cc1c[nH]c2ccc(O)cc12"),
    (("indole-3-acetic acid", "iaa", "auxin", "heteroauxin"), "OC(=O)Cc1c[nH]c2ccccc12"),
    (("indole-3-butyric acid", "iba"), "OC(=O)CCCc1c[nH]c2ccccc12"),
    (("kynurenine", "l-kynurenine"), "N[C@@H](CC(=O)c1ccccc1N)C(=O)O"),
    (("kynurenic acid",), "OC(=O)c1cc(O)c2ccccc2n1"),
    (("quinolinic acid", "pyridine-2,3-dicarboxylic acid"), "OC(=O)c1cccnc1C(=O)O"),
    (("xanthurenic acid",), "OC(=O)c1cc(O)c2cccc(O)c2n1"),
    (("indoxyl sulfate", "3-indoxyl sulfate", "indican"), "OS(=O)(=O)Oc1c[nH]c2ccccc12"),
    (("agmatine",), "N=C(N)NCCCCN"),
    (("gamma-hydroxybutyric acid", "ghb", "γ-hydroxybutyric acid", "4-hydroxybutanoic acid"), "OCCCC(=O)O"),
    (("sodium oxybate", "ghb sodium salt"), "[Na+].OCCCC(=O)[O-]"),
    (("succinic semialdehyde", "4-oxobutanoic acid"), "O=CCCC(=O)O"),
    (("acetoacetic acid", "acetoacetate", "3-oxobutanoic acid"), "CC(=O)CC(=O)O"),
    (("beta-hydroxybutyric acid", "β-hydroxybutyric acid", "3-hydroxybutyric acid", "bhb", "(r)-3-hydroxybutyrate", "d-3-hydroxybutyrate"),
     "C[C@@H](O)CC(=O)O"),
    (("n-acetylaspartate", "naa", "n-acetyl-l-aspartic acid"), "CC(=O)N[C@@H](CC(=O)O)C(=O)O"),
    (("allantoin",), "NC(=O)NC1NC(=O)NC1=O"),
    (("orotic acid", "vitamin b13"), "OC(=O)c1cc(=O)[nH]c(=O)[nH]1"),
    (("carbamoyl phosphate",), "NC(=O)OP(=O)(O)O"),
    (("cystathionine", "l-cystathionine"), "N[C@@H](CCSC[C@H](N)C(=O)O)C(=O)O"),
    (("nmda", "n-methyl-d-aspartic acid", "n-methyl-d-aspartate"), "CN[C@H](CC(=O)O)C(=O)O"),
    (("ampa", "α-amino-3-hydroxy-5-methyl-4-isoxazolepropionic acid", "alpha-amino-3-hydroxy-5-methyl-4-isoxazolepropionic acid"),
     "Cc1onc(O)c1CC(N)C(=O)O"),
    (("muscimol",), "NCc1cc(=O)[nH]o1"),
    (("ibotenic acid",), "NC(C(=O)O)c1cc(=O)[nH]o1"),
    (("muscarine", "(+)-muscarine", "l-(+)-muscarine"), "C[C@@H]1O[C@H](C[N+](C)(C)C)C[C@H]1O"),
    (("carbachol", "carbamylcholine"), "C[N+](C)(C)CCOC(N)=O.[Cl-]"),
    (("hexamethonium", "hexamethonium bromide"), "C[N+](C)(C)CCCCCC[N+](C)(C)C.[Br-].[Br-]"),
    (("6-hydroxydopamine", "6-ohda", "oxidopamine"), "NCCc1cc(O)c(O)cc1O"),
    (("mptp", "1-methyl-4-phenyl-1,2,3,6-tetrahydropyridine"), "CN1CCC(=CC1)c1ccccc1"),
    (("mpp+", "1-methyl-4-phenylpyridinium", "cyperquat"), "C[n+]1ccc(cc1)-c1ccccc1"),
    (("paraquat", "methyl viologen", "1,1'-dimethyl-4,4'-bipyridinium dichloride"),
     "C[n+]1ccc(cc1)-c1cc[n+](C)cc1.[Cl-].[Cl-]"),
    (("diquat",), "c1cc[n+]2CC[n+]3ccccc3c2c1.[Br-].[Br-]"),
    (("cotinine", "(s)-cotinine", "(-)-cotinine"), "CN1CC[C@H](C1=O)c1cccnc1"),
    (("theophylline", "1,3-dimethylxanthine"), "Cn1c2nc[nH]c2c(=O)n(C)c1=O"),
    (("theobromine", "3,7-dimethylxanthine"), "Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
    (("paraxanthine", "1,7-dimethylxanthine"), "Cn1cnc2c1c(=O)n(C)c(=O)[nH]2"),
    (("1-methylxanthine",), "Cn1c(=O)c2[nH]cnc2[nH]c1=O"),
    (("pentoxifylline", "trental"), "CC(=O)CCCCn1c(=O)c2c(ncn2C)n(C)c1=O"),
    (("thyroxine", "t4", "levothyroxine", "l-thyroxine", "3,5,3',5'-tetraiodothyronine"),
     "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O"),
    (("triiodothyronine", "t3", "liothyronine", "3,5,3'-triiodothyronine"),
     "N[C@@H](Cc1cc(I)c(Oc2ccc(O)c(I)c2)c(I)c1)C(=O)O"),
    (("reverse t3", "rt3", "3,3',5'-triiodothyronine"), "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)cc1)C(=O)O"),
    (("thyronine",), "N[C@@H](Cc1ccc(Oc2ccc(O)cc2)cc1)C(=O)O"),
    (("3-iodothyronamine", "t1am"), "NCCc1cc(I)c(Oc2ccc(O)cc2)cc1"),
    (("isosorbide dinitrate", "isdn"), "[O-][N+](=O)O[C@H]1CO[C@@H]2[C@H]1OC[C@H]2O[N+](=O)[O-]"),
    (("isosorbide mononitrate", "ismn", "isosorbide 5-mononitrate"),
     "O[C@H]1CO[C@@H]2[C@H]1OC[C@H]2O[N+](=O)[O-]"),
    (("sodium nitroprusside", "snp", "nitroprusside"), "[Na+].[Na+].N#C[Fe-2](C#N)(C#N)(C#N)(C#N)N=O"),
    (("s-nitrosoglutathione", "gsno"), "N[C@@H](CCC(=O)N[C@@H](CSN=O)C(=O)NCC(=O)O)C(=O)O"),
    (("s-nitroso-n-acetylpenicillamine", "snap"), "CC(=O)N[C@@H](C(C)(C)SN=O)C(=O)O"),
    (("dea nonoate", "dea-nonoate", "deanonoate"), "[Na+].CCN(CC)N([O-])N=O"),
    (("l-nio", "n-iminoethyl-l-ornithine"), "CC(=N)NCCC[C@H](N)C(=O)O"),
    # ------------------------------------------------------------------ #
    # Biochemistry buffers and chelators
    # ------------------------------------------------------------------ #
    (("tris", "tris base", "tris(hydroxymethyl)aminomethane", "trometamol", "tromethamine", "thame", "2-amino-2-(hydroxymethyl)propane-1,3-diol"),
     "NC(CO)(CO)CO"),
    (("tris hydrochloride", "tris-hcl", "tris hcl", "tris.hcl", "trizma hydrochloride"), "Cl.NC(CO)(CO)CO"),
    (("hepes", "4-(2-hydroxyethyl)-1-piperazineethanesulfonic acid", "4-(2-hydroxyethyl)piperazine-1-ethanesulfonic acid"),
     "OCCN1CCN(CCS(=O)(=O)O)CC1"),
    (("hepes sodium salt", "hepes-na", "sodium hepes"), "[Na+].OCCN1CCN(CCS(=O)(=O)[O-])CC1"),
    (("heppso", "4-(2-hydroxyethyl)piperazine-1-(2-hydroxypropanesulfonic acid)"),
     "OCCN1CCN(CC(O)CS(=O)(=O)O)CC1"),
    (("epps", "hepps", "4-(2-hydroxyethyl)-1-piperazinepropanesulfonic acid (epps)"),
     "OCCN1CCN(CCCS(=O)(=O)O)CC1"),
    (("mes", "2-(n-morpholino)ethanesulfonic acid", "2-morpholinoethanesulfonic acid", "4-morpholineethanesulfonic acid"),
     "OS(=O)(=O)CCN1CCOCC1"),
    (("mops", "3-(n-morpholino)propanesulfonic acid", "3-morpholinopropane-1-sulfonic acid"),
     "OS(=O)(=O)CCCN1CCOCC1"),
    (("mopso", "3-morpholino-2-hydroxypropanesulfonic acid"), "OC(CS(=O)(=O)O)CN1CCOCC1"),
    (("pipes", "piperazine-n,n'-bis(2-ethanesulfonic acid)", "1,4-piperazinediethanesulfonic acid"),
     "OS(=O)(=O)CCN1CCN(CCS(=O)(=O)O)CC1"),
    (("bis-tris", "bistris", "bis-tris methane", "2-[bis(2-hydroxyethyl)amino]-2-(hydroxymethyl)propane-1,3-diol"),
     "OCCN(CCO)C(CO)(CO)CO"),
    (("bis-tris propane", "1,3-bis[tris(hydroxymethyl)methylamino]propane", "btp"),
     "OCC(CO)(CO)NCCCNC(CO)(CO)CO"),
    (("caps", "3-(cyclohexylamino)-1-propanesulfonic acid", "n-cyclohexyl-3-aminopropanesulfonic acid"),
     "OS(=O)(=O)CCCNC1CCCCC1"),
    (("capso", "3-(cyclohexylamino)-2-hydroxy-1-propanesulfonic acid"), "OC(CS(=O)(=O)O)CNC1CCCCC1"),
    (("ches", "2-(cyclohexylamino)ethanesulfonic acid", "n-cyclohexyl-2-aminoethanesulfonic acid"),
     "OS(=O)(=O)CCNC1CCCCC1"),
    (("taps", "n-tris(hydroxymethyl)methyl-3-aminopropanesulfonic acid", "3-[[tris(hydroxymethyl)methyl]amino]propanesulfonic acid"),
     "OCC(CO)(CO)NCCCS(=O)(=O)O"),
    (("tapso",), "OCC(CO)(CO)NCC(O)CS(=O)(=O)O"),
    (("tes", "n-tris(hydroxymethyl)methyl-2-aminoethanesulfonic acid", "2-[[tris(hydroxymethyl)methyl]amino]ethanesulfonic acid"),
     "OCC(CO)(CO)NCCS(=O)(=O)O"),
    (("bicine", "n,n-bis(2-hydroxyethyl)glycine"), "OCCN(CCO)CC(=O)O"),
    (("tricine", "n-[tris(hydroxymethyl)methyl]glycine"), "OCC(CO)(CO)NCC(=O)O"),
    (("aces", "n-(2-acetamido)-2-aminoethanesulfonic acid"), "NC(=O)CNCCS(=O)(=O)O"),
    (("ada", "n-(2-acetamido)iminodiacetic acid"), "NC(=O)CN(CC(=O)O)CC(=O)O"),
    (("bes", "n,n-bis(2-hydroxyethyl)-2-aminoethanesulfonic acid"), "OCCN(CCO)CCS(=O)(=O)O"),
    (("ampso", "3-[(1,1-dimethyl-2-hydroxyethyl)amino]-2-hydroxypropanesulfonic acid"),
     "CC(C)(CO)NCC(O)CS(=O)(=O)O"),
    (("cabs", "4-(cyclohexylamino)-1-butanesulfonic acid"), "OS(=O)(=O)CCCCNC1CCCCC1"),
    (("imidazole hydrochloride",), "Cl.c1cnc[nH]1"),
    (("sodium cacodylate", "cacodylate", "sodium dimethylarsinate"), "[Na+].C[As](C)(=O)[O-]"),
    (("cacodylic acid", "dimethylarsinic acid"), "C[As](C)(=O)O"),
    (("sodium citrate", "trisodium citrate", "citrate buffer salt"),
     "[Na+].[Na+].[Na+].[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]"),
    (("potassium citrate", "tripotassium citrate"), "[K+].[K+].[K+].[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]"),
    (("ammonium bicarbonate", "nh4hco3", "ammonium hydrogen carbonate"), "[NH4+].OC(=O)[O-]"),
    (("ammonium carbonate", "(nh4)2co3"), "[NH4+].[NH4+].[O-]C(=O)[O-]"),
    (("triethylammonium acetate", "teaa"), "CC(=O)[O-].CC[NH+](CC)CC"),
    (("triethylammonium bicarbonate", "teab"), "OC(=O)[O-].CC[NH+](CC)CC"),
    (("triethylamine hydrochloride", "et3n.hcl"), "[Cl-].CC[NH+](CC)CC"),
    (("edta", "ethylenediaminetetraacetic acid", "edetic acid", "(ethylenedinitrilo)tetraacetic acid", "versene"),
     "OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O"),
    (("edta disodium", "disodium edta", "na2edta", "edta disodium salt", "edta-na2", "edetate disodium"),
     "[Na+].[Na+].OC(=O)CN(CCN(CC(=O)[O-])CC(=O)O)CC(=O)[O-]"),
    (("edta tetrasodium", "na4edta", "edta tetrasodium salt"),
     "[Na+].[Na+].[Na+].[Na+].[O-]C(=O)CN(CCN(CC(=O)[O-])CC(=O)[O-])CC(=O)[O-]"),
    (("egta", "ethylene glycol-bis(β-aminoethyl ether)-n,n,n',n'-tetraacetic acid", "ethylene glycol tetraacetic acid", "egtazic acid"),
     "OC(=O)CN(CCOCCOCCN(CC(=O)O)CC(=O)O)CC(=O)O"),
    (("bapta", "1,2-bis(o-aminophenoxy)ethane-n,n,n',n'-tetraacetic acid"),
     "OC(=O)CN(CC(=O)O)c1ccccc1OCCOc1ccccc1N(CC(=O)O)CC(=O)O"),
    (("dtpa", "diethylenetriaminepentaacetic acid", "pentetic acid"),
     "OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CCN(CC(=O)O)CC(=O)O"),
    (("dota", "1,4,7,10-tetraazacyclododecane-1,4,7,10-tetraacetic acid"),
     "OC(=O)CN1CCN(CC(=O)O)CCN(CC(=O)O)CCN(CC(=O)O)CC1"),
    (("nota", "1,4,7-triazacyclononane-1,4,7-triacetic acid"), "OC(=O)CN1CCN(CC(=O)O)CCN(CC(=O)O)CC1"),
    (("cyclen", "1,4,7,10-tetraazacyclododecane"), "C1CNCCNCCNCCN1"),
    (("cyclam", "1,4,8,11-tetraazacyclotetradecane"), "C1CNCCNCCCNCCNC1"),
    (("tacn", "1,4,7-triazacyclononane"), "C1CNCCNCCN1"),
    (("hedta", "n-(2-hydroxyethyl)ethylenediaminetriacetic acid"), "OCCN(CCN(CC(=O)O)CC(=O)O)CC(=O)O"),
    (("deferoxamine", "desferrioxamine", "desferal", "dfo", "deferoxamine b"),
     "CC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN"),
    (("deferoxamine mesylate", "desferrioxamine mesylate", "dfo mesylate"),
     "CS(=O)(=O)O.CC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN"),
    (("deferiprone", "1,2-dimethyl-3-hydroxypyridin-4-one"), "Cc1n(C)ccc(=O)c1O"),
    (("deferasirox",), "OC(=O)c1ccc(cc1)-n1nc(-c2ccccc2O)nc1-c1ccccc1O"),
    (("8-hydroxyquinoline", "oxine", "8-quinolinol"), "Oc1cccc2cccnc12"),
    (("dipicolinic acid", "pyridine-2,6-dicarboxylic acid"), "OC(=O)c1cccc(n1)C(=O)O"),
    (("bathocuproine disulfonate", "bcs"),
     "Cc1cc(-c2ccc(cc2)S(=O)(=O)[O-])c2ccc3c(-c4ccc(cc4)S(=O)(=O)[O-])cc(C)nc3c2n1.[Na+].[Na+]"),
    (("ferrozine",), "[Na+].[O-]S(=O)(=O)c1ccc(cc1)-c1cc(nnc1-c1ccccn1)-c1ccc(cc1)S(=O)(=O)O"),
    (("dithizone", "diphenylthiocarbazone"), "S=C(NNc1ccccc1)N=Nc1ccccc1"),
    (("dimercaprol", "bal", "british anti-lewisite", "2,3-dimercaptopropanol"), "OCC(S)CS"),
    (("dmsa", "succimer", "meso-2,3-dimercaptosuccinic acid"), "OC(=O)[C@@H](S)[C@H](S)C(=O)O"),
    # ------------------------------------------------------------------ #
    # Reducing agents, alkylating agents, denaturants and detergents
    # ------------------------------------------------------------------ #
    (("dtt", "dithiothreitol", "cleland's reagent", "clelands reagent", "(2s,3s)-1,4-dimercaptobutane-2,3-diol", "threo-1,4-dimercapto-2,3-butanediol", "d-dtt"),
     "O[C@H](CS)[C@H](O)CS"),
    (("dte", "dithioerythritol", "erythro-1,4-dimercapto-2,3-butanediol"), "OC(CS)C(O)CS"),
    (("tcep", "tris(2-carboxyethyl)phosphine", "tris(2-carboxyethyl)phosphine hydrochloride", "tcep hydrochloride", "tcep-hcl", "tcep.hcl"),
     "OC(=O)CCP(CCC(=O)O)CCC(=O)O"),
    (("2-mercaptoethylamine", "cysteamine", "2-aminoethanethiol", "mercaptamine"), "NCCS"),
    (("cysteamine hydrochloride",), "Cl.NCCS"),
    (("cystamine", "2,2'-dithiobis(ethylamine)"), "NCCSSCCN"),
    (("cystamine dihydrochloride",), "Cl.Cl.NCCSSCCN"),
    (("dtnb", "ellman's reagent", "ellmans reagent", "5,5'-dithiobis(2-nitrobenzoic acid)"),
     "OC(=O)c1cc(SSc2ccc(c(c2)C(=O)O)[N+](=O)[O-])ccc1[N+](=O)[O-]"),
    (("tnb", "5-thio-2-nitrobenzoic acid", "2-nitro-5-thiobenzoic acid"), "OC(=O)c1cc(S)ccc1[N+](=O)[O-]"),
    (("iodoacetamide", "iam", "2-iodoacetamide"), "NC(=O)CI"),
    (("chloroacetamide", "2-chloroacetamide", "caa"), "NC(=O)CCl"),
    (("sds", "sodium dodecyl sulfate", "sodium dodecyl sulphate", "sodium lauryl sulfate", "sls", "sodium lauryl sulphate"),
     "[Na+].CCCCCCCCCCCCOS(=O)(=O)[O-]"),
    (("lithium dodecyl sulfate", "lds"), "[Li+].CCCCCCCCCCCCOS(=O)(=O)[O-]"),
    (("chaps", "3-[(3-cholamidopropyl)dimethylammonio]-1-propanesulfonate"),
     "C[C@H](CCC(=O)NCCC[N+](C)(C)CCCS(=O)(=O)[O-])[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    (("chapso",),
     "C[C@H](CCC(=O)NCCC[N+](C)(C)CC(O)CS(=O)(=O)[O-])[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    (("triton x-100", "triton x100", "octylphenol ethoxylate", "t-octylphenoxypolyethoxyethanol", "polyethylene glycol tert-octylphenyl ether"),
     "CC(C)(C)CC(C)(C)c1ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCO)cc1"),
    (("tween-20", "tween 20", "polysorbate 20"),
     "CCCCCCCCCCCC(=O)OCCOCCOC[C@H]1O[C@@H](OCCOCCO)[C@@H](OCCO)[C@@H]1OCCO"),
    (("tween-80", "tween 80", "polysorbate 80", "polyoxyethylene (80) sorbitan monooleate"),
     "CCCCCCCC/C=C\\CCCCCCCC(=O)OCCOCCOC[C@H]1O[C@@H](OCCOCCO)[C@@H](OCCO)[C@@H]1OCCO"),
    (("span 80", "sorbitan monooleate"), "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](O)[C@H]1OC[C@H](O)[C@H]1O"),
    (("np-40", "nonidet p-40", "igepal ca-630", "nonylphenol ethoxylate"),
     "CCCCCCCCCc1ccc(OCCOCCOCCOCCOCCOCCOCCOCCOCCO)cc1"),
    (("brij-35", "brij 35", "polyoxyethylene (23) lauryl ether"), "CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
    (("brij-58", "brij 58", "polyoxyethylene (20) cetyl ether"),
     "CCCCCCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
    (("octyl glucoside", "n-octyl-β-d-glucopyranoside", "og", "n-octyl-beta-d-glucopyranoside", "octyl-β-d-glucoside", "β-og"),
     "CCCCCCCCO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),
    (("dodecyl maltoside", "n-dodecyl-β-d-maltoside", "ddm", "n-dodecyl-beta-d-maltoside", "lauryl maltoside"),
     "CCCCCCCCCCCCO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O"),
    (("decyl maltoside", "n-decyl-β-d-maltoside", "dm", "n-decyl-beta-d-maltoside"),
     "CCCCCCCCCCO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O"),
    (("lmng", "lauryl maltose neopentyl glycol"),
     "CCCCCCCCCCCCC(CO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O)(CCCCCCCCCCCC)CO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@H](O)[C@H]1O"),
    (("cetyltrimethylammonium bromide", "ctab", "cetrimonium bromide", "hexadecyltrimethylammonium bromide"),
     "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]"),
    (("cetyltrimethylammonium chloride", "ctac", "cetrimonium chloride"),
     "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Cl-]"),
    (("benzalkonium chloride", "bac", "benzyldodecyldimethylammonium chloride", "bzk"),
     "CCCCCCCCCCCC[N+](C)(C)Cc1ccccc1.[Cl-]"),
    (("dodecyltrimethylammonium bromide", "dtab"), "CCCCCCCCCCCC[N+](C)(C)C.[Br-]"),
    (("didodecyldimethylammonium bromide", "ddab"), "CCCCCCCCCCCC[N+](C)(C)CCCCCCCCCCCC.[Br-]"),
    (("sarkosyl", "sodium lauroyl sarcosinate", "n-lauroylsarcosine sodium salt", "sodium n-lauroylsarcosinate"),
     "[Na+].CCCCCCCCCCCC(=O)N(C)CC(=O)[O-]"),
    (("fos-choline-12", "n-dodecylphosphocholine", "dpc"), "CCCCCCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C"),
    (("zwittergent 3-12", "n-dodecyl-n,n-dimethyl-3-ammonio-1-propanesulfonate", "sb-12", "sb3-12"),
     "CCCCCCCCCCCC[N+](C)(C)CCCS(=O)(=O)[O-]"),
    (("lauryldimethylamine n-oxide", "ldao", "dodecyldimethylamine oxide", "dda0"),
     "CCCCCCCCCCCC[N+](C)(C)[O-]"),
    (("sodium cholate", "cholate"),
     "[Na+].C[C@H](CCC(=O)[O-])[C@H]1CC[C@@H]2[C@@]1([C@H](C[C@H]3[C@H]2[C@@H](C[C@H]4[C@@]3(CC[C@H](C4)O)C)O)O)C"),
    # ------------------------------------------------------------------ #
    # PEGs, linkers and protease inhibitors
    # ------------------------------------------------------------------ #
    (("peg-400", "peg 400", "polyethylene glycol 400", "polyethylene glycol", "peg", "poly(ethylene glycol)"),
     "OCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
    (("peg-200", "peg 200", "polyethylene glycol 200", "tetraethylene glycol", "peg-4", "peg4", "tetraglycol"),
     "OCCOCCOCCOCCO"),
    (("pentaethylene glycol", "peg-5"), "OCCOCCOCCOCCOCCO"),
    (("hexaethylene glycol", "peg-6", "peg6"), "OCCOCCOCCOCCOCCOCCO"),
    (("peg-3350", "peg 3350", "polyethylene glycol 3350", "macrogol 3350", "macrogol"),
     "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO"),
    (("mpeg", "methoxypolyethylene glycol", "peg monomethyl ether", "mpeg-oh"), "COCCOCCOCCOCCOCCOCCOCCO"),
    (("peg dimethyl ether", "polyethylene glycol dimethyl ether", "pegdme"), "COCCOCCOCCOCCOCCOCCOCCOC"),
    (("peg-diacrylate", "pegda", "poly(ethylene glycol) diacrylate"), "C=CC(=O)OCCOCCOCCOCCOCCOCCOC(=O)C=C"),
    (("peg-dimethacrylate", "pegdma"), "CC(=C)C(=O)OCCOCCOCCOCCOCCOCCOC(=O)C(C)=C"),
    (("azido-peg4-amine", "n3-peg4-nh2", "azido-peg4-nh2"), "NCCOCCOCCOCCOCCN=[N+]=[N-]"),
    (("azido-peg4-acid", "n3-peg4-cooh", "azido-peg4-cooh"), "OC(=O)CCOCCOCCOCCOCCN=[N+]=[N-]"),
    (("azido-peg4-nhs ester", "n3-peg4-nhs"), "[N-]=[N+]=NCCOCCOCCOCCOCCC(=O)ON1C(=O)CCC1=O"),
    (("amino-peg4-acid", "nh2-peg4-cooh", "amino-peg4-cooh"), "NCCOCCOCCOCCOCCC(=O)O"),
    (("boc-nh-peg4-cooh", "boc-amino-peg4-acid", "boc-peg4-acid"), "CC(C)(C)OC(=O)NCCOCCOCCOCCOCCC(=O)O"),
    (("fmoc-nh-peg4-cooh", "fmoc-amino-peg4-acid"), "OC(=O)CCOCCOCCOCCOCCNC(=O)OCC1c2ccccc2-c2ccccc12"),
    (("mal-peg4-nhs", "maleimide-peg4-nhs ester", "nhs-peg4-maleimide", "mal-peg4-nhs ester", "smpeg4"),
     "O=C1C=CC(=O)N1CCC(=O)NCCOCCOCCOCCOCCC(=O)ON1C(=O)CCC1=O"),
    (("biotin-peg4-acid", "biotin-peg4-cooh"),
     "OC(=O)CCOCCOCCOCCOCCNC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@@H]12"),
    (("pmsf", "phenylmethylsulfonyl fluoride", "phenylmethanesulfonyl fluoride"), "O=S(=O)(F)Cc1ccccc1"),
    (("aebsf", "4-(2-aminoethyl)benzenesulfonyl fluoride", "pefabloc"), "NCCc1ccc(cc1)S(=O)(=O)F"),
    (("aebsf hydrochloride", "pefabloc sc"), "Cl.NCCc1ccc(cc1)S(=O)(=O)F"),
    (("benzamidine",), "NC(=N)c1ccccc1"),
    (("benzamidine hydrochloride",), "Cl.NC(=N)c1ccccc1"),
    (("leupeptin",), "CC(C)C[C@H](NC(C)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCNC(N)=N)C=O"),
    (("pepstatin", "pepstatin a"),
     "CC(C)C[C@H](NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)CC(C)C)C(C)C)C(C)C)[C@@H](O)CC(=O)N[C@@H](C)C(=O)N[C@@H](CC(C)C)[C@@H](O)CC(=O)O"),
    (("e-64",), "CC(C)C[C@H](NC(=O)[C@@H]1O[C@H]1C(=O)O)C(=O)NCCCCNC(N)=N"),
    (("bestatin", "ubenimex"), "CC(C)C[C@H](NC(=O)[C@H](O)[C@H](N)Cc1ccccc1)C(=O)O"),
    (("sodium orthovanadate", "na3vo4", "sodium vanadate"), "[Na+].[Na+].[Na+].[O-][V](=O)([O-])[O-]"),
    (("beta-glycerophosphate", "β-glycerophosphate", "glycerol 2-phosphate", "sodium β-glycerophosphate", "sodium beta-glycerophosphate"),
     "[Na+].[Na+].OCC(CO)OP(=O)([O-])[O-]"),
    (("spectinomycin",),
     "CN[C@@H]1[C@@H](O)[C@@H](NC)[C@@H]2O[C@@]3(O)[C@@H](O[C@@H]2[C@@H]1O)CC(=O)[C@H](C)O3"),
    (("neomycin", "neomycin b"),
     "NC[C@H]1O[C@H](O[C@@H]2[C@@H](N)C[C@@H](N)[C@H](O[C@H]3O[C@H](CO)[C@@H](O[C@@H]4O[C@H](CN)[C@@H](O)[C@H](O)[C@H]4N)[C@H]3O)[C@@H]2O)[C@H](N)[C@@H](O)[C@@H]1O"),
    # ------------------------------------------------------------------ #
    # Enzyme / reporter substrates and assay reagents
    # ------------------------------------------------------------------ #
    (("iptg", "isopropyl β-d-1-thiogalactopyranoside", "isopropyl beta-d-1-thiogalactopyranoside", "isopropylthiogalactoside"),
     "CC(C)S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"),
    (("x-gal", "5-bromo-4-chloro-3-indolyl β-d-galactopyranoside", "5-bromo-4-chloro-3-indolyl beta-d-galactopyranoside", "bcig"),
     "OC[C@H]1O[C@@H](Oc2c[nH]c3ccc(Br)c(Cl)c23)[C@H](O)[C@@H](O)[C@H]1O"),
    (("onpg", "o-nitrophenyl β-d-galactopyranoside", "2-nitrophenyl beta-d-galactopyranoside"),
     "OC[C@H]1O[C@@H](Oc2ccccc2[N+](=O)[O-])[C@H](O)[C@@H](O)[C@H]1O"),
    (("pnpg", "p-nitrophenyl β-d-glucopyranoside", "4-nitrophenyl beta-d-glucopyranoside"),
     "OC[C@H]1O[C@@H](Oc2ccc(cc2)[N+](=O)[O-])[C@H](O)[C@@H](O)[C@@H]1O"),
    (("pnpp", "p-nitrophenyl phosphate", "4-nitrophenyl phosphate", "pnp-phosphate"),
     "O=P(O)(O)Oc1ccc(cc1)[N+](=O)[O-]"),
    (("4-methylumbelliferone", "4-mu", "hymecromone", "7-hydroxy-4-methylcoumarin"),
     "Cc1cc(=O)oc2cc(O)ccc12"),
    (("4-methylumbelliferyl phosphate", "mup", "4-mup"), "Cc1cc(=O)oc2cc(OP(=O)(O)O)ccc12"),
    (("mug", "4-methylumbelliferyl β-d-glucuronide", "4-methylumbelliferyl beta-d-glucuronide"),
     "Cc1cc(=O)oc2cc(O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)C(=O)O)ccc12"),
    (("luminol", "5-amino-2,3-dihydrophthalazine-1,4-dione"), "Nc1cccc2C(=O)NNC(=O)c12"),
    (("d-luciferin", "firefly luciferin", "luciferin"), "OC(=O)[C@@H]1CSC(=N1)c1nc2ccc(O)cc2s1"),
    (("mtt", "3-(4,5-dimethylthiazol-2-yl)-2,5-diphenyltetrazolium bromide", "thiazolyl blue"),
     "CC1=C(N=C(S1)[N+]2=NC(=NN2C3=CC=CC=C3)C4=CC=CC=C4)C.[Br-]"),
    (("resazurin", "resazurin sodium", "alamar blue", "alamarblue"),
     "[Na+].[O-]c1ccc2[n+]([O-])c3ccc(=O)cc3oc2c1"),
    (("resorufin", "7-hydroxy-3h-phenoxazin-3-one"), "Oc1ccc2nc3ccc(=O)cc3oc2c1"),
    (("amplex red", "n-acetyl-3,7-dihydroxyphenoxazine", "10-acetyl-3,7-dihydroxyphenoxazine"),
     "CC(=O)N1c2ccc(O)cc2Oc2cc(O)ccc12"),
    (("nbt", "nitro blue tetrazolium", "nitrotetrazolium blue chloride"),
     "COc1cc(ccc1-c1ccc(c(OC)c1)-[n+]1nc(nn1-c1ccccc1)-c1ccc(cc1)[N+](=O)[O-])-[n+]1nc(nn1-c1ccccc1)-c1ccc(cc1)[N+](=O)[O-].[Cl-].[Cl-]"),
    (("bcip", "5-bromo-4-chloro-3-indolyl phosphate", "x-phosphate"), "O=P(O)(O)Oc1c[nH]c2ccc(Br)c(Cl)c12"),
    (("tmb", "3,3',5,5'-tetramethylbenzidine"), "Cc1cc(cc(C)c1N)-c1cc(C)c(N)c(C)c1"),
    (("abts", "2,2'-azino-bis(3-ethylbenzothiazoline-6-sulfonic acid)"),
     "CCN1/C(=N/N=C2\\Sc3cc(ccc3N2CC)S(=O)(=O)O)Sc2cc(ccc12)S(=O)(=O)O"),
    (("o-phenylenediamine dihydrochloride", "opd dihydrochloride"), "Cl.Cl.Nc1ccccc1N"),
    (("3,3'-diaminobenzidine", "diaminobenzidine"), "Nc1ccc(cc1N)-c1ccc(N)c(N)c1"),
    (("bromophenol blue", "bpb", "3',3'',5',5''-tetrabromophenolsulfonphthalein"),
     "Oc1c(Br)cc(cc1Br)C1(OS(=O)(=O)c2ccccc12)c1cc(Br)c(O)c(Br)c1"),
    (("orange g", "acid orange 10"), "[Na+].[Na+].Oc1ccc2cc(cc(c2c1N=Nc1ccccc1)S(=O)(=O)[O-])S(=O)(=O)[O-]"),
    # ------------------------------------------------------------------ #
    # Dyes, indicators and stains
    # ------------------------------------------------------------------ #
    (("ethidium bromide", "etbr", "ethidium"), "CC[n+]1c(-c2ccccc2)c2cc(N)ccc2c2ccc(N)cc12.[Br-]"),
    (("propidium iodide", "pi"), "CC[N+](C)(CC)CCC[n+]1c(-c2ccccc2)c2cc(N)ccc2c2ccc(N)cc12.[I-].[I-]"),
    (("dapi", "4',6-diamidino-2-phenylindole"), "NC(=N)c1ccc(cc1)-c1cc2cc(ccc2[nH]1)C(N)=N"),
    (("dapi dihydrochloride",), "Cl.Cl.NC(=N)c1ccc(cc1)-c1cc2cc(ccc2[nH]1)C(N)=N"),
    (("hoechst 33258", "bisbenzimide", "pibenzimol"),
     "CN1CCN(CC1)c1ccc2[nH]c(nc2c1)-c1ccc2[nH]c(nc2c1)-c1ccc(O)cc1"),
    (("hoechst 33342",), "CCOc1ccc(cc1)-c1nc2cc(ccc2[nH]1)-c1nc2cc(ccc2[nH]1)N1CCN(C)CC1"),
    (("acridine orange",), "CN(C)c1ccc2cc3ccc(cc3nc2c1)N(C)C"),
    (("methylene blue", "methylthioninium chloride", "mb"), "CN(C)c1ccc2nc3ccc(cc3[s+]c2c1)N(C)C.[Cl-]"),
    (("toluidine blue", "toluidine blue o", "tolonium chloride"), "Cc1cc2nc3ccc(cc3[s+]c2cc1N)N(C)C.[Cl-]"),
    (("crystal violet", "gentian violet", "methyl violet 10b", "basic violet 3"),
     "CN(C)c1ccc(cc1)C(=C1C=CC(=[N+](C)C)C=C1)c1ccc(cc1)N(C)C.[Cl-]"),
    (("malachite green", "malachite green chloride", "basic green 4"),
     "CN(C)c1ccc(cc1)C(=C1C=CC(=[N+](C)C)C=C1)c1ccccc1.[Cl-]"),
    (("malachite green oxalate",), "CN(C)c1ccc(cc1)C(=C1C=CC(=[N+](C)C)C=C1)c1ccccc1.OC(=O)C(=O)[O-]"),
    (("safranin", "safranin o", "basic red 2"), "Cc1cc2nc3cc(C)c(N)cc3[n+](-c3ccccc3)c2cc1N.[Cl-]"),
    (("congo red",),
     "[Na+].[Na+].Nc1c(N=Nc2ccc(cc2)-c2ccc(cc2)N=Nc2c(N)c3ccccc3c(c2)S(=O)(=O)[O-])cc(c2ccccc12)S(=O)(=O)[O-]"),
    (("phenol red", "phenolsulfonphthalein", "psp"), "Oc1ccc(cc1)C1(OS(=O)(=O)c2ccccc12)c1ccc(O)cc1"),
    (("phenolphthalein",), "Oc1ccc(cc1)C1(OC(=O)c2ccccc12)c1ccc(O)cc1"),
    (("bromothymol blue", "btb"), "CC(C)c1cc(C2(OS(=O)(=O)c3ccccc23)c2cc(C(C)C)c(O)c(Br)c2C)c(C)c(Br)c1O"),
    (("bromocresol green", "bcg"), "Cc1c(Br)c(O)c(Br)cc1C1(OS(=O)(=O)c2ccccc12)c1cc(Br)c(O)c(Br)c1C"),
    (("bromocresol purple",), "Cc1cc(C2(OS(=O)(=O)c3ccccc23)c2cc(C)c(O)c(Br)c2)cc(Br)c1O"),
    (("methyl orange", "helianthine", "acid orange 52"), "[Na+].CN(C)c1ccc(cc1)N=Nc1ccc(cc1)S(=O)(=O)[O-]"),
    (("methyl red",), "CN(C)c1ccc(cc1)N=Nc1ccccc1C(=O)O"),
    (("thymol blue", "thymolsulfonphthalein"),
     "CC(C)c1cc(C2(OS(=O)(=O)c3ccccc23)c2cc(C(C)C)c(O)cc2C)c(C)cc1O"),
    (("alizarin", "1,2-dihydroxyanthraquinone"), "Oc1ccc2C(=O)c3ccccc3C(=O)c2c1O"),
    (("alizarin red s", "sodium alizarinsulfonate"), "[Na+].Oc1c(O)c2C(=O)c3ccccc3C(=O)c2cc1S(=O)(=O)[O-]"),
    (("murexide", "ammonium purpurate"), "[NH4+].O=C1NC(=O)C(=NC2=C([O-])NC(=O)NC2=O)C(=O)N1"),
    (("ninhydrin", "2,2-dihydroxyindane-1,3-dione", "triketohydrindene hydrate"), "OC1(O)C(=O)c2ccccc2C1=O"),
    (("opa", "o-phthalaldehyde", "phthalaldehyde", "ortho-phthalaldehyde"), "O=Cc1ccccc1C=O"),
    (("bicinchoninic acid", "bca", "2,2'-biquinoline-4,4'-dicarboxylic acid"),
     "OC(=O)c1cc(nc2ccccc12)-c1cc(C(=O)O)c2ccccc2n1"),
    # ------------------------------------------------------------------ #
    # Fluorophores, dyes for biology, bioconjugation and click reagents
    # ------------------------------------------------------------------ #
    (("fluorescein", "resorcinolphthalein", "solvent yellow 94"),
     "O=C1OC2(c3ccccc13)c1ccc(O)cc1Oc1cc(O)ccc12"),
    (("fluorescein sodium", "uranine", "sodium fluorescein", "fluorescein disodium", "acid yellow 73"),
     "[Na+].[Na+].[O-]C(=O)c1ccccc1-c1c2ccc(=O)cc2oc2cc([O-])ccc12"),
    (("fitc", "fluorescein isothiocyanate", "fluorescein 5-isothiocyanate", "5-fitc", "fluorescein-5-isothiocyanate"),
     "O=C1OC2(c3cc(N=C=S)ccc13)c1ccc(O)cc1Oc1cc(O)ccc12"),
    (("5-fam", "5-carboxyfluorescein", "carboxyfluorescein", "fam", "5(6)-carboxyfluorescein", "5(6)-fam", "fam-5(6)"),
     "O=C1OC2(c3cc(C(=O)O)ccc13)c1ccc(O)cc1Oc1cc(O)ccc12"),
    (("6-fam", "6-carboxyfluorescein"), "O=C1OC2(c3ccc(C(=O)O)cc13)c1ccc(O)cc1Oc1cc(O)ccc12"),
    (("fluorescein diacetate", "fda"), "CC(=O)Oc1ccc2c(c1)Oc1cc(OC(C)=O)ccc1C21OC(=O)c2ccccc21"),
    (("dcfh-da", "2',7'-dichlorodihydrofluorescein diacetate", "h2dcfda", "dcfda", "dcf-da"),
     "CC(=O)Oc1ccc2c(c1Cl)Oc1c(Cl)c(OC(C)=O)ccc1C2c1ccccc1C(=O)O"),
    (("2',7'-dichlorofluorescein", "dcf"), "O=C1OC2(c3ccccc13)c1ccc(O)c(Cl)c1Oc1c(Cl)c(O)ccc12"),
    (("calcein",), "OC(=O)CN(CC(=O)O)Cc1c(O)ccc2c1Oc1c(CN(CC(=O)O)CC(=O)O)c(O)ccc1C21OC(=O)c2ccccc21"),
    (("erythrosin b", "erythrosine", "tetraiodofluorescein disodium", "fd&c red 3"),
     "[Na+].[Na+].[O-]C(=O)c1ccccc1-c1c2cc(I)c(=O)c(I)c-2oc2c(I)c([O-])c(I)cc12"),
    (("rose bengal",),
     "[Na+].[Na+].[O-]C(=O)c1c(Cl)c(Cl)c(Cl)c(Cl)c1-c1c2cc(I)c(=O)c(I)c-2oc2c(I)c([O-])c(I)cc12"),
    (("rhodamine b", "rhodamine b chloride", "basic violet 10", "rhb"),
     "CCN(CC)c1ccc2c(-c3ccccc3C(=O)O)c3ccc(=[N+](CC)CC)cc3oc2c1.[Cl-]"),
    (("rhodamine 6g", "r6g", "rhodamine 6g chloride", "basic red 1"),
     "CCNc1cc2oc3cc(=[NH+]CC)c(C)cc3c(-c3ccccc3C(=O)OCC)c2cc1C.[Cl-]"),
    (("rhodamine 110", "rh110", "rhodamine 110 chloride"),
     "Nc1ccc2c(-c3ccccc3C(=O)O)c3ccc(=[NH2+])cc3oc2c1.[Cl-]"),
    (("rhodamine 123", "rh123"), "COC(=O)c1ccccc1-c1c2ccc(N)cc2oc2cc(=[NH2+])ccc12.[Cl-]"),
    (("tamra", "5-tamra", "5-carboxytetramethylrhodamine", "carboxytetramethylrhodamine", "tetramethylrhodamine carboxylic acid"),
     "CN(C)c1ccc2c(-c3ccc(cc3C(=O)[O-])C(=O)O)c3ccc(=[N+](C)C)cc3oc2c1"),
    (("tmrm", "tetramethylrhodamine methyl ester"),
     "CN(C)c1ccc2c(-c3ccccc3C(=O)OC)c3ccc(=[N+](C)C)cc3oc2c1.[Cl-]"),
    (("tmre", "tetramethylrhodamine ethyl ester"),
     "CCOC(=O)c1ccccc1-c1c2ccc(cc2oc2cc(=[N+](C)C)ccc12)N(C)C.[Cl-]"),
    (("tritc", "tetramethylrhodamine isothiocyanate", "tetramethylrhodamine-5-isothiocyanate"),
     "CN(C)c1ccc2c(-c3ccc(N=C=S)cc3C(=O)[O-])c3ccc(=[N+](C)C)cc3oc2c1"),
    (("sulforhodamine b", "kiton red", "acid red 52"),
     "CCN(CC)c1ccc2c(-c3cc(ccc3S(=O)(=O)[O-])S(=O)(=O)O)c3ccc(=[N+](CC)CC)cc3oc2c1"),
    (("cy3", "cy3 core", "cyanine 3", "1,1',3,3,3',3'-hexamethylindocarbocyanine iodide"),
     "C[N+]1=C(C=CC=C2N(C)c3ccccc3C2(C)C)C(C)(C)c2ccccc21.[I-]"),
    (("cy5", "cy5 core", "cyanine 5", "1,1',3,3,3',3'-hexamethylindodicarbocyanine iodide"),
     "C[N+]1=C(C=CC=CC=C2N(C)c3ccccc3C2(C)C)C(C)(C)c2ccccc21.[I-]"),
    (("cy7", "cy7 core", "cyanine 7"), "C[N+]1=C(C=CC=CC=CC=C2N(C)c3ccccc3C2(C)C)C(C)(C)c2ccccc21.[I-]"),
    (("dii", "diic18(3)", "1,1'-dioctadecyl-3,3,3',3'-tetramethylindocarbocyanine perchlorate"),
     "CCCCCCCCCCCCCCCCCC[N+]1=C(C=CC=C2N(CCCCCCCCCCCCCCCCCC)c3ccccc3C2(C)C)C(C)(C)c2ccccc21.[O-]Cl(=O)(=O)=O"),
    (("dio", "dioc18(3)", "3,3'-dioctadecyloxacarbocyanine perchlorate"),
     "CCCCCCCCCCCCCCCCCC[n+]1c(C=CC=C2N(CCCCCCCCCCCCCCCCCC)c3ccccc3O2)oc2ccccc21.[O-]Cl(=O)(=O)=O"),
    (("indocyanine green", "icg", "cardiogreen"),
     "[Na+].CC1(C)C(=CC=CC=CC=CC2=[N+](CCCCS(=O)(=O)[O-])c3ccc4ccccc4c3C2(C)C)N(CCCCS(=O)(=O)[O-])c2ccc3ccccc3c21"),
    (("bodipy", "bodipy core", "4,4-difluoro-4-bora-3a,4a-diaza-s-indacene"),
     "F[B-]1(F)n2cccc2C=C2C=CC=[N+]21"),
    (("bodipy fl", "bodipy fl acid", "bodipy fl carboxylic acid", "4,4-difluoro-5,7-dimethyl-4-bora-3a,4a-diaza-s-indacene-3-propionic acid"),
     "Cc1cc(C)n2c1C=C1C=CC(CCC(=O)O)=[N+]1[B-]2(F)F"),
    (("bodipy 493/503",), "Cc1cc(C)c2C=C3C(C)=C(C)C(C)=[N+]3[B-](F)(F)n12"),
    (("coumarin", "2h-chromen-2-one", "1,2-benzopyrone"), "O=c1ccc2ccccc2o1"),
    (("umbelliferone", "7-hydroxycoumarin"), "Oc1ccc2ccc(=O)oc2c1"),
    (("7-amino-4-methylcoumarin", "amc", "coumarin 120"), "Cc1cc(=O)oc2cc(N)ccc12"),
    (("coumarin 1", "7-diethylamino-4-methylcoumarin", "coumarin 460"), "CCN(CC)c1ccc2c(C)cc(=O)oc2c1"),
    (("coumarin 6", "3-(2-benzothiazolyl)-7-(diethylamino)coumarin"),
     "CCN(CC)c1ccc2cc(-c3nc4ccccc4s3)c(=O)oc2c1"),
    (("coumarin 343",), "OC(=O)c1cc2cc3CCCN4CCCc(c34)c2oc1=O"),
    (("7-diethylaminocoumarin-3-carboxylic acid", "deac"), "CCN(CC)c1ccc2cc(C(=O)O)c(=O)oc2c1"),
    (("dansyl chloride", "5-(dimethylamino)naphthalene-1-sulfonyl chloride", "dns-cl"),
     "CN(C)c1cccc2c(cccc12)S(=O)(=O)Cl"),
    (("dansyl amide", "dansylamide"), "CN(C)c1cccc2c(cccc12)S(N)(=O)=O"),
    (("dansyl glycine",), "CN(C)c1cccc2c(cccc12)S(=O)(=O)NCC(=O)O"),
    (("nbd-cl", "4-chloro-7-nitrobenzofurazan", "4-chloro-7-nitro-2,1,3-benzoxadiazole"),
     "[O-][N+](=O)c1ccc(Cl)c2nonc12"),
    (("nbd-f", "4-fluoro-7-nitrobenzofurazan"), "[O-][N+](=O)c1ccc(F)c2nonc12"),
    (("pyrene",), "c1cc2ccc3cccc4ccc(c1)c2c34"),
    (("1-pyrenebutyric acid", "pyrenebutyric acid", "4-(1-pyrenyl)butyric acid"),
     "OC(=O)CCCc1ccc2ccc3cccc4ccc1c2c34"),
    (("1-pyrenebutyric acid nhs ester", "pyrenebutyric acid n-hydroxysuccinimide ester", "pba-nhs"),
     "O=C1CCC(=O)N1OC(=O)CCCc1ccc2ccc3cccc4ccc1c2c34"),
    (("1-aminopyrene",), "Nc1ccc2ccc3cccc4ccc1c2c34"),
    (("perylene",), "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54"),
    (("9,10-diphenylanthracene",), "c1ccc(cc1)-c1c2ccccc2c(-c2ccccc2)c2ccccc12"),
    (("prodan", "6-propionyl-2-(dimethylamino)naphthalene"), "CCC(=O)c1ccc2cc(ccc2c1)N(C)C"),
    (("laurdan", "6-dodecanoyl-2-dimethylaminonaphthalene"), "CCCCCCCCCCCC(=O)c1ccc2cc(ccc2c1)N(C)C"),
    (("dph", "1,6-diphenyl-1,3,5-hexatriene", "diphenylhexatriene"), "c1ccc(cc1)/C=C/C=C/C=C/c1ccccc1"),
    (("ans", "8-anilino-1-naphthalenesulfonic acid", "1-anilinonaphthalene-8-sulfonic acid", "1,8-ans"),
     "OS(=O)(=O)c1cccc2c(Nc3ccccc3)cccc12"),
    (("thioflavin t", "tht", "thioflavine t", "basic yellow 1"),
     "Cc1ccc2[n+](C)c(sc2c1)-c1ccc(cc1)N(C)C.[Cl-]"),
    (("nile red", "9-diethylamino-5h-benzo[a]phenoxazin-5-one"), "CCN(CC)c1ccc2nc3c(cc(=O)c4ccccc34)oc2c1"),
    (("nile blue", "nile blue a", "nile blue chloride", "basic blue 12"),
     "CCN(CC)c1ccc2nc3c(cc(N)c4ccccc34)[o+]c2c1.[Cl-]"),
    (("dihydroethidium", "dhe", "hydroethidine"), "CCN1c2ccc(N)cc2-c2cc(N)ccc2C1c1ccccc1"),
    (("lucigenin", "bis-n-methylacridinium nitrate"),
     "C[n+]1c2ccccc2c(-c2c3ccccc3[n+](C)c3ccccc23)c2ccccc12.[O-][N+](=O)[O-].[O-][N+](=O)[O-]"),
    (("dabcyl", "4-((4-(dimethylamino)phenyl)azo)benzoic acid", "dabcyl acid"),
     "CN(C)c1ccc(cc1)N=Nc1ccc(cc1)C(=O)O"),
    (("dabcyl nhs ester", "dabcyl-se"), "CN(C)c1ccc(cc1)N=Nc1ccc(cc1)C(=O)ON1C(=O)CCC1=O"),
    # bioconjugation crosslinkers
    (("smcc", "succinimidyl 4-(n-maleimidomethyl)cyclohexane-1-carboxylate"),
     "O=C1CCC(=O)N1OC(=O)C1CCC(CC1)CN1C(=O)C=CC1=O"),
    (("sulfo-smcc",), "[Na+].O=C1CC(S(=O)(=O)[O-])C(=O)N1OC(=O)C1CCC(CC1)CN1C(=O)C=CC1=O"),
    (("spdp", "n-succinimidyl 3-(2-pyridyldithio)propionate"), "O=C1CCC(=O)N1OC(=O)CCSSc1ccccn1"),
    (("dss", "disuccinimidyl suberate"), "O=C1CCC(=O)N1OC(=O)CCCCCCC(=O)ON1C(=O)CCC1=O"),
    (("bs3", "bis(sulfosuccinimidyl) suberate", "bs(3)", "sulfo-dss"),
     "[Na+].[Na+].O=C1CC(S(=O)(=O)[O-])C(=O)N1OC(=O)CCCCCCC(=O)ON1C(=O)CC(S(=O)(=O)[O-])C1=O"),
    (("dsg", "disuccinimidyl glutarate"), "O=C1CCC(=O)N1OC(=O)CCCC(=O)ON1C(=O)CCC1=O"),
    (("dsp", "dithiobis(succinimidyl propionate)", "lomant's reagent", "lomants reagent"),
     "O=C1CCC(=O)N1OC(=O)CCSSCCC(=O)ON1C(=O)CCC1=O"),
    (("bmoe", "bismaleimidoethane", "1,2-bismaleimidoethane"), "O=C1C=CC(=O)N1CCN1C(=O)C=CC1=O"),
    (("bmh", "bismaleimidohexane", "1,6-bismaleimidohexane"), "O=C1C=CC(=O)N1CCCCCCN1C(=O)C=CC1=O"),
    (("sata", "n-succinimidyl s-acetylthioacetate"), "CC(=O)SCC(=O)ON1C(=O)CCC1=O"),
    (("traut's reagent", "trauts reagent", "2-iminothiolane", "2-iminothiolane hydrochloride"),
     "Cl.N=C1CCCS1"),
    (("maleimidopropionic acid", "3-maleimidopropionic acid"), "OC(=O)CCN1C(=O)C=CC1=O"),
    (("maleimidohexanoic acid", "6-maleimidohexanoic acid", "mhc acid", "emca"), "OC(=O)CCCCCN1C(=O)C=CC1=O"),
    (("mal-nhs", "n-succinimidyl 6-maleimidohexanoate", "emcs"), "O=C1CCC(=O)N1OC(=O)CCCCCN1C(=O)C=CC1=O"),
    (("n-(2-aminoethyl)maleimide", "aminoethylmaleimide", "aem"), "NCCN1C(=O)C=CC1=O"),
    (("n-(2-aminoethyl)maleimide trifluoroacetate", "aem tfa"), "OC(=O)C(F)(F)F.NCCN1C(=O)C=CC1=O"),
    (("n-succinimidyl iodoacetate", "sia"), "O=C1CCC(=O)N1OC(=O)CI"),
    (("n-succinimidyl bromoacetate", "sba"), "O=C1CCC(=O)N1OC(=O)CBr"),
    (("ethyl acetimidate hydrochloride", "ethyl acetimidate"), "Cl.CCOC(C)=N"),
    (("dimethyl suberimidate", "dimethyl suberimidate dihydrochloride"), "Cl.Cl.COC(=N)CCCCCCC(=N)OC"),
    (("dimethyl adipimidate",), "Cl.Cl.COC(=N)CCCCC(=N)OC"),
    (("4-formylbenzoic acid", "p-formylbenzoic acid"), "OC(=O)c1ccc(C=O)cc1"),
    (("4-formylbenzoate nhs ester", "4-formylbenzoic acid nhs ester"), "O=Cc1ccc(cc1)C(=O)ON1C(=O)CCC1=O"),
    (("aminooxyacetic acid", "aoa", "carboxymethoxylamine"), "NOCC(=O)O"),
    (("aminooxyacetic acid hemihydrochloride", "carboxymethoxylamine hemihydrochloride"),
     "Cl.NOCC(=O)O.NOCC(=O)O"),
    (("aminoguanidine",), "NNC(N)=N"),
    (("aminoguanidine hydrochloride", "aminoguanidine hcl"), "Cl.NNC(N)=N"),
    (("acetohydrazide",), "CC(=O)NN"),
    (("benzhydrazide",), "NNC(=O)c1ccccc1"),
    (("adipic dihydrazide", "adh"), "NNC(=O)CCCCC(=O)NN"),
    (("girard's reagent t", "girard reagent t"), "C[N+](C)(C)CC(=O)NN.[Cl-]"),
    (("2-hydrazinopyridine",), "NNc1ccccn1"),
    (("6-hydrazinonicotinic acid", "hynic"), "NNc1ccc(cn1)C(=O)O"),
    # click chemistry
    (("dbco-acid", "dibenzocyclooctyne-acid", "dbco acid", "4-(11,12-didehydrodibenzo[b,f]azocin-5(6h)-yl)-4-oxobutanoic acid"),
     "OC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-amine", "dibenzocyclooctyne-amine", "dbco amine"), "NCCNC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-nhs ester", "dibenzocyclooctyne-n-hydroxysuccinimidyl ester", "dbco-nhs", "dbco nhs ester"),
     "O=C1CCC(=O)N1OC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-peg4-nhs ester", "dbco-peg4-nhs"),
     "O=C1CCC(=O)N1OC(=O)CCOCCOCCOCCOCCNC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-peg4-acid",), "OC(=O)CCOCCOCCOCCOCCNC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-peg4-amine",), "NCCOCCOCCOCCOCCNC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco-maleimide",), "O=C1C=CC(=O)N1CCNC(=O)CCC(=O)N1Cc2ccccc2C#Cc2ccccc21"),
    (("dbco", "dibenzocyclooctyne", "dibenzoazacyclooctyne", "aza-dibenzocyclooctyne", "adibo", "dibac"),
     "C1Nc2ccccc2C#Cc2ccccc21"),
    (("bcn", "bicyclo[6.1.0]non-4-yne", "bcn-oh", "(1r,8s,9s)-bicyclo[6.1.0]non-4-yn-9-ylmethanol", "bcn alcohol"),
     "OCC1C2CCC#CCCC12"),
    (("bcn-nhs", "bcn-nhs carbonate", "bcn n-hydroxysuccinimidyl carbonate"),
     "O=C1CCC(=O)N1OC(=O)OCC1C2CCC#CCCC12"),
    (("cyclooctyne",), "C1CCC#CCCC1"),
    (("tco", "trans-cyclooctene", "(e)-cyclooctene"), "C1CC/C=C/CCC1"),
    (("tco-oh", "trans-cyclooct-4-enol", "(e)-cyclooct-4-en-1-ol", "tco alcohol"), "OC1CC/C=C/CCC1"),
    (("tco-nhs ester", "tco-nhs", "trans-cyclooctene nhs carbonate"), "O=C1CCC(=O)N1OC(=O)OC1CC/C=C/CCC1"),
    (("tco-peg4-nhs ester", "tco-peg4-nhs"), "O=C1CCC(=O)N1OC(=O)CCOCCOCCOCCOCCNC(=O)OC1CC/C=C/CCC1"),
    (("tco-amine", "tco-nh2"), "NCCNC(=O)OC1CC/C=C/CCC1"),
    (("methyltetrazine-amine", "methyltetrazine amine", "(4-(6-methyl-1,2,4,5-tetrazin-3-yl)phenyl)methanamine"),
     "Cc1nnc(nn1)-c1ccc(CN)cc1"),
    (("methyltetrazine-nhs ester", "methyltetrazine-nhs"), "Cc1nnc(nn1)-c1ccc(cc1)C(=O)ON1C(=O)CCC1=O"),
    (("methyltetrazine-peg4-nhs ester",), "Cc1nnc(nn1)-c1ccc(cc1)C(=O)NCCOCCOCCOCCOCCC(=O)ON1C(=O)CCC1=O"),
    (("methyltetrazine-acid", "methyltetrazine-cooh", "4-(6-methyl-1,2,4,5-tetrazin-3-yl)benzoic acid"),
     "Cc1nnc(nn1)-c1ccc(cc1)C(=O)O"),
    (("tetrazine", "1,2,4,5-tetrazine", "s-tetrazine"), "c1nncnn1"),
    (("3,6-dimethyl-1,2,4,5-tetrazine",), "Cc1nnc(C)nn1"),
    (("3,6-di(2-pyridyl)-1,2,4,5-tetrazine", "dipyridyl tetrazine", "dptz"),
     "c1ccc(nc1)-c1nnc(nn1)-c1ccccn1"),
    (("3,6-diphenyl-1,2,4,5-tetrazine",), "c1ccc(cc1)-c1nnc(nn1)-c1ccccc1"),
    (("propargyl-peg4-nhs ester", "alkyne-peg4-nhs"), "C#CCOCCOCCOCCOCCC(=O)ON1C(=O)CCC1=O"),
    (("propargyl-peg4-acid", "alkyne-peg4-acid"), "C#CCOCCOCCOCCOCCC(=O)O"),
    (("propargyl-peg4-amine", "alkyne-peg4-amine"), "C#CCOCCOCCOCCOCCN"),
    (("4-pentynoic acid", "pentynoic acid"), "C#CCCC(=O)O"),
    (("4-pentynoic acid nhs ester",), "C#CCCC(=O)ON1C(=O)CCC1=O"),
    (("azidoacetic acid", "2-azidoacetic acid"), "OC(=O)CN=[N+]=[N-]"),
    (("azidoacetic acid nhs ester", "azidoacetic nhs"), "[N-]=[N+]=NCC(=O)ON1C(=O)CCC1=O"),
    (("3-azidopropylamine", "3-azido-1-propanamine"), "NCCCN=[N+]=[N-]"),
    (("3-azidopropanoic acid", "3-azidopropionic acid"), "OC(=O)CCN=[N+]=[N-]"),
    (("6-azidohexanoic acid",), "OC(=O)CCCCCN=[N+]=[N-]"),
    (("4-azidobenzoic acid",), "OC(=O)c1ccc(cc1)N=[N+]=[N-]"),
    (("4-azidoaniline",), "Nc1ccc(cc1)N=[N+]=[N-]"),
    (("benzyl azide", "(azidomethyl)benzene"), "[N-]=[N+]=NCc1ccccc1"),
    (("phenyl azide", "azidobenzene"), "[N-]=[N+]=Nc1ccccc1"),
    (("1-azidoadamantane", "adamantyl azide"), "[N-]=[N+]=NC12CC3CC(CC(C3)C1)C2"),
    (("picolyl azide", "2-(azidomethyl)pyridine"), "[N-]=[N+]=NCc1ccccn1"),
    (("azido-peg3-amine", "n3-peg3-nh2"), "NCCOCCOCCOCCN=[N+]=[N-]"),
    (("thpta", "tris(3-hydroxypropyltriazolylmethyl)amine", "tris((1-hydroxypropyl-1h-1,2,3-triazol-4-yl)methyl)amine"),
     "OCCCn1cc(CN(Cc2cn(CCCO)nn2)Cc2cn(CCCO)nn2)nn1"),
    (("tbta", "tris(benzyltriazolylmethyl)amine", "tris((1-benzyl-1h-1,2,3-triazol-4-yl)methyl)amine"),
     "C(n1cc(CN(Cc2cn(Cc3ccccc3)nn2)Cc2cn(Cc3ccccc3)nn2)nn1)c1ccccc1"),
    (("bttaa",), "OC(=O)Cn1cc(CN(Cc2cn(CC(C)(C)C)nn2)Cc2cn(CC(C)(C)C)nn2)nn1"),
    (("copper(ii) sulfate pentahydrate", "cuso4.5h2o", "copper sulfate pentahydrate", "blue vitriol"),
     "[Cu+2].[O-]S(=O)(=O)[O-].O.O.O.O.O"),
    (("bathophenanthroline disulfonate", "bphen disulfonate", "bathophenanthrolinedisulfonic acid disodium"),
     "[Na+].[Na+].[O-]S(=O)(=O)c1ccc(cc1)-c1ccnc2c1ccc1c(ccnc12)-c1ccc(cc1)S(=O)(=O)[O-]"),
    (("2-(diphenylphosphino)benzoic acid", "dppba"), "OC(=O)c1ccccc1P(c1ccccc1)c1ccccc1"),
    (("methyl 2-(diphenylphosphino)benzoate", "staudinger-bertozzi ester"),
     "COC(=O)c1ccccc1P(c1ccccc1)c1ccccc1"),
    (("norbornene carboxylic acid", "5-norbornene-2-carboxylic acid"), "OC(=O)C1CC2C=CC1C2"),
    # ------------------------------------------------------------------ #
    # Alkaloids and other natural products
    # ------------------------------------------------------------------ #
    (("caffeine", "1,3,7-trimethylxanthine", "guaranine", "theine", "methyltheobromine"),
     "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    (("nicotine", "(s)-nicotine", "l-nicotine", "(-)-nicotine"), "CN1CCC[C@H]1c1cccnc1"),
    (("morphine", "(-)-morphine"), "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"),
    (("codeine", "methylmorphine"), "CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"),
    (("heroin", "diacetylmorphine", "diamorphine"),
     "CN1CC[C@]23c4c5ccc(OC(C)=O)c4O[C@H]2[C@@H](OC(C)=O)C=C[C@H]3[C@H]1C5"),
    (("oxycodone",), "COc1ccc2CC3N(C)CCC45C(Oc1c24)C(=O)CCC35O"),
    (("hydrocodone",), "COc1ccc2CC3N(C)CCC45C(Oc1c24)C(=O)CCC35"),
    (("hydromorphone",), "Oc1ccc2CC3N(C)CCC45C(Oc1c24)C(=O)CCC35"),
    (("naloxone",), "C=CCN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=O)CC[C@@]3(O)[C@H]1C5"),
    (("naltrexone",), "O=C1CC[C@@]2(O)[C@H]3Cc4ccc(O)c5O[C@@H]1[C@]2(CCN3CC1CC1)c45"),
    (("methadone",), "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1"),
    (("fentanyl",), "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1"),
    (("tramadol",), "COc1cccc(c1)C1(O)CCCCC1CN(C)C"),
    (("cocaine", "(-)-cocaine", "benzoylmethylecgonine"),
     "CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c1ccccc1)C(=O)OC"),
    (("atropine", "dl-hyoscyamine", "tropine tropate"), "CN1[C@H]2CC[C@@H]1CC(C2)OC(=O)C(CO)c1ccccc1"),
    (("scopolamine", "hyoscine", "(-)-scopolamine"),
     "CN1[C@@H]2C[C@@H](C[C@H]1[C@@H]1O[C@H]21)OC(=O)[C@H](CO)c1ccccc1"),
    (("tropine", "3-tropanol"), "CN1[C@H]2CC[C@@H]1CC(O)C2"),
    (("pilocarpine", "(+)-pilocarpine"), "CC[C@H]1[C@H](COC1=O)Cc1cn(C)cn1"),
    (("physostigmine", "eserine"), "CNC(=O)Oc1ccc2N(C)[C@@H]3N(C)CC[C@@]3(C)c2c1"),
    (("neostigmine", "neostigmine bromide"), "CN(C)C(=O)Oc1cccc(c1)[N+](C)(C)C.[Br-]"),
    (("pyridostigmine", "pyridostigmine bromide"), "CN(C)C(=O)Oc1ccc[n+](C)c1.[Br-]"),
    (("strychnine",), "C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C1=CC=CC=C15"),
    (("brucine",), "COC1=C(C=C2C(=C1)[C@]34CCN5[C@H]3C[C@@H]6[C@@H]7[C@@H]4N2C(=O)C[C@@H]7OCC=C6C5)OC"),
    (("quinine", "(-)-quinine"), "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1"),
    (("quinidine", "(+)-quinidine"), "COc1ccc2nccc([C@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1"),
    (("cinchonidine", "(-)-cinchonidine"), "C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccccc12"),
    (("cinchonine", "(+)-cinchonine"), "C=C[C@H]1CN2CC[C@H]1C[C@@H]2[C@@H](O)c1ccnc2ccccc12"),
    (("lsd", "lysergic acid diethylamide", "lysergide", "d-lsd"),
     "CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34"),
    (("lysergic acid",), "CN1C[C@H](C(=O)O)C=C2c3cccc4[nH]cc(c34)C[C@H]21"),
    (("tetrahydrocannabinol", "thc", "delta-9-thc", "δ9-thc", "dronabinol", "(-)-trans-δ9-tetrahydrocannabinol"),
     "CCCCCc1cc(O)c2[C@@H]3C=C(C)CC[C@H]3C(C)(C)Oc2c1"),
    (("cannabidiol", "cbd"), "CCCCCc1cc(O)c(c(O)c1)[C@@H]1C=C(C)CC[C@H]1C(C)=C"),
    (("cannabinol", "cbn"), "CCCCCc1cc(O)c2c(c1)OC(C)(C)c1ccc(C)cc12"),
    (("colchicine",), "CC(=O)N[C@H]1CCc2cc(OC)c(OC)c(OC)c2-c2ccc(OC)c(=O)cc12"),
    (("camptothecin", "(s)-camptothecin", "cpt"), "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccccc4cc3Cn1c2=O"),
    (("topotecan",), "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccc(O)c(CN(C)C)c4cc3Cn1c2=O"),
    (("irinotecan",), "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccc(OC(=O)N5CCC(CC5)N5CCCCC5)cc4c(CC)c3Cn1c2=O"),
    (("paclitaxel", "taxol"),
     "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)C5=CC=CC=C5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)C6=CC=CC=C6)C7=CC=CC=C7)O)(C)C"),
    (("docetaxel",),
     "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)O)OC(=O)C5=CC=CC=C5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)OC(C)(C)C)C6=CC=CC=C6)O)(C)C"),
    (("artemisinin", "qinghaosu"), "C[C@@H]1CC[C@H]2[C@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4"),
    (("artesunate",), "C[C@@H]1CC[C@H]2[C@H](C)[C@H](OC(=O)CCC(=O)O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4"),
    (("artemether",), "C[C@@H]1CC[C@H]2[C@H](C)[C@@H](OC)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4"),
    (("dihydroartemisinin", "artenimol"),
     "C[C@@H]1CC[C@H]2[C@H](C)[C@H](O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@@]23OO4"),
    # ------------------------------------------------------------------ #
    # Antibiotics, macrolides and cytotoxic natural products
    # ------------------------------------------------------------------ #
    (("penicillin g", "benzylpenicillin", "penicillin", "pen g"),
     "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"),
    (("penicillin g sodium", "benzylpenicillin sodium"),
     "[Na+].CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)[O-]"),
    (("penicillin v", "phenoxymethylpenicillin"), "CC1(C)S[C@@H]2[C@H](NC(=O)COc3ccccc3)C(=O)N2[C@H]1C(=O)O"),
    (("amoxicillin",), "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O"),
    (("ampicillin",), "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O"),
    (("carbenicillin",), "CC1(C)S[C@@H]2[C@H](NC(=O)C(C(=O)O)c3ccccc3)C(=O)N2[C@H]1C(=O)O"),
    (("6-apa", "6-aminopenicillanic acid"), "CC1(C)S[C@@H]2[C@H](N)C(=O)N2[C@H]1C(=O)O"),
    (("clavulanic acid",), "OC(=O)[C@@H]1N2C(=O)C[C@H]2O/C1=C\\CO"),
    (("cephalexin", "cefalexin"), "CC1=C(N2[C@@H]([C@@H](C2=O)NC(=O)[C@@H](c3ccccc3)N)SC1)C(=O)O"),
    (("cefazolin",), "Cc1nnc(SCC2=C(N3[C@@H]([C@@H](C3=O)NC(=O)Cn3cnnn3)SC2)C(=O)O)s1"),
    (("ceftriaxone",),
     "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(CSc3nc(=O)c(=O)[nH]n3C)CS[C@H]12)c1csc(N)n1"),
    (("cefuroxime",), "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(COC(N)=O)CS[C@H]12)c1ccco1"),
    (("meropenem",), "C[C@@H]1[C@@H]2[C@H](C(=O)N2C(=C1S[C@@H]1CN[C@@H](C1)C(=O)N(C)C)C(=O)O)[C@@H](C)O"),
    (("imipenem",), "C[C@@H](O)[C@@H]1C(=O)N2C(C(=O)O)=C(SCCNC=N)C[C@@H]12"),
    (("aztreonam",), "C[C@@H]1[C@H](NC(=O)/C(=N\\OC(C)(C)C(=O)O)c2csc(N)n2)C(=O)N1S(=O)(=O)O"),
    (("erythromycin", "erythromycin a"),
     "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"),
    (("azithromycin",),
     "CC[C@@H]1[C@@]([C@@H]([C@H](N(C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O"),
    (("clarithromycin",),
     "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)OC)C)C)O)(C)O"),
    (("tetracycline",),
     "CN(C)[C@H]1[C@@H]2C[C@H]3C(=C(O)c4c(O)cccc4[C@@]3(C)O)C(=O)[C@]2(O)C(O)=C(C(N)=O)C1=O"),
    (("doxycycline",),
     "CN(C)[C@H]1[C@@H]2[C@@H](O)[C@H]3C(=C(O)c4c(O)cccc4[C@H]3C)C(=O)[C@]2(O)C(O)=C(C(N)=O)C1=O"),
    (("minocycline",), "CN(C)c1ccc(O)c2c1C[C@H]1C[C@H]3[C@H](N(C)C)C(O)=C(C(N)=O)C(=O)[C@@]3(O)C(O)=C1C2=O"),
    (("chloramphenicol",), "O=C(N[C@H](CO)[C@@H](O)c1ccc(cc1)[N+](=O)[O-])C(Cl)Cl"),
    (("kanamycin", "kanamycin a"),
     "C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)N)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)O)N"),
    (("kanamycin sulfate",),
     "OS(=O)(=O)O.C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)N)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CN)O)O)O)N"),
    (("rifampicin", "rifampin"),
     "CO[C@H]1/C=C/O[C@@]2(C)Oc3c(C)c(O)c4c(O)c(/C=N/N5CCN(C)CC5)c(NC(=O)/C(C)=C\\C=C\\[C@H](C)[C@H](O)[C@@H](C)[C@@H](O)[C@@H](C)[C@H](OC(C)=O)[C@@H]1C)c(O)c4c3C2=O"),
    (("puromycin",), "CN(C)c1ncnc2n(cnc12)[C@H]1O[C@H](CO)[C@@H](NC(=O)[C@H](N)Cc2ccc(OC)cc2)[C@H]1O"),
    (("puromycin dihydrochloride",),
     "Cl.Cl.CN(C)c1ncnc2n(cnc12)[C@H]1O[C@H](CO)[C@@H](NC(=O)[C@H](N)Cc2ccc(OC)cc2)[C@H]1O"),
    (("cycloheximide", "chx"), "C[C@H]1C[C@@H](C)C(=O)[C@@H](C1=O)[C@H](O)CC1CC(=O)NC(=O)C1"),
    (("staurosporine",), "C[C@@]12[C@@H]([C@@H](C[C@@H](O1)N1c3ccccc3-c3c4C(=O)NCc4c4c5ccccc5n2c4c31)NC)OC"),
    (("doxorubicin", "adriamycin"),
     "C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O"),
    (("daunorubicin",),
     "C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)C)O)N)O"),
    (("doxorubicin hydrochloride",),
     "Cl.C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C=CC=C5OC)O)(C(=O)CO)O)N)O"),
    (("mitomycin c",), "CO[C@]12[C@H](COC(N)=O)C3=C(N1C[C@@H]1N[C@@H]21)C(=O)C(C)=C(N)C3=O"),
    # ------------------------------------------------------------------ #
    # Statins, flavonoids and polyphenols
    # ------------------------------------------------------------------ #
    (("lovastatin", "mevinolin"),
     "CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C"),
    (("simvastatin",),
     "CCC(C)(C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@@H]3C[C@H](CC(=O)O3)O)C"),
    (("pravastatin",),
     "CC[C@H](C)C(=O)O[C@H]1C[C@H](C=C2[C@H]1[C@H]([C@H](C=C2)C)CC[C@H](C[C@H](CC(=O)O)O)O)O"),
    (("atorvastatin",), "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@H](O)C[C@H](O)CC(=O)O"),
    (("rosuvastatin",), "CC(C)c1nc(nc(c1/C=C/[C@H](O)C[C@H](O)CC(=O)O)-c1ccc(F)cc1)N(C)S(C)(=O)=O"),
    (("fluvastatin",), "CC(C)n1c(/C=C/[C@H](O)C[C@H](O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc12"),
    (("quercetin",), "Oc1cc(O)c2c(c1)oc(-c1ccc(O)c(O)c1)c(O)c2=O"),
    (("kaempferol",), "Oc1ccc(cc1)-c1oc2cc(O)cc(O)c2c(=O)c1O"),
    (("myricetin",), "Oc1cc(O)c2c(c1)oc(-c1cc(O)c(O)c(O)c1)c(O)c2=O"),
    (("luteolin",), "Oc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(O)c(O)c1"),
    (("apigenin",), "Oc1ccc(cc1)-c1cc(=O)c2c(O)cc(O)cc2o1"),
    (("chrysin",), "Oc1cc(O)c2c(c1)oc(cc2=O)-c1ccccc1"),
    (("genistein",), "Oc1ccc(cc1)-c1coc2cc(O)cc(O)c2c1=O"),
    (("daidzein",), "Oc1ccc(cc1)-c1coc2cc(O)ccc2c1=O"),
    (("naringenin", "(s)-naringenin"), "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"),
    (("hesperetin",), "COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"),
    (("catechin", "(+)-catechin"), "Oc1cc(O)c2C[C@H](O)[C@H](Oc2c1)c1ccc(O)c(O)c1"),
    (("epicatechin", "(-)-epicatechin"), "Oc1cc(O)c2C[C@H](O)[C@@H](Oc2c1)c1ccc(O)c(O)c1"),
    (("epigallocatechin gallate", "egcg", "epigallocatechin-3-gallate"),
     "C1[C@@H]([C@H](Oc2cc(O)cc(O)c21)c1cc(O)c(O)c(O)c1)OC(=O)c1cc(O)c(O)c(O)c1"),
    (("rutin", "rutoside"),
     "C[C@H]1O[C@@H](OC[C@H]2O[C@@H](Oc3c(oc4cc(O)cc(O)c4c3=O)-c3ccc(O)c(O)c3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O"),
    (("resveratrol", "trans-resveratrol", "(e)-resveratrol", "3,5,4'-trihydroxystilbene"),
     "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1"),
    (("curcumin", "diferuloylmethane"), "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"),
    (("capsaicin",), "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O"),
    (("dihydrocapsaicin",), "COc1cc(CNC(=O)CCCCCCC(C)C)ccc1O"),
    (("piperine",), "O=C(/C=C/C=C/c1ccc2OCOc2c1)N1CCCCC1"),
    (("gingerol", "[6]-gingerol", "6-gingerol"), "CCCCC[C@H](O)CC(=O)CCc1ccc(O)c(OC)c1"),
    (("allicin",), "C=CCS(=O)SCC=C"),
    (("alliin",), "C=CCS(=O)C[C@H](N)C(=O)O"),
    (("diallyl disulfide",), "C=CCSSCC=C"),
    (("sulforaphane",), "CS(=O)CCCCN=C=S"),
    (("allyl isothiocyanate", "mustard oil", "aitc"), "C=CCN=C=S"),
    # ------------------------------------------------------------------ #
    # Terpenes, flavour and fragrance compounds
    # ------------------------------------------------------------------ #
    (("menthol", "(-)-menthol", "l-menthol", "(1r,2s,5r)-menthol"), "C[C@@H]1CC[C@H]([C@@H](C1)O)C(C)C"),
    (("menthone", "(-)-menthone"), "C[C@@H]1CC[C@H](C(C)C)C(=O)C1"),
    (("limonene", "(r)-(+)-limonene", "d-limonene", "(r)-limonene"), "CC1=CC[C@@H](CC1)C(=C)C"),
    (("(s)-limonene", "l-limonene", "(s)-(-)-limonene"), "CC1=CC[C@H](CC1)C(=C)C"),
    (("camphor", "(1r)-(+)-camphor", "d-camphor", "(+)-camphor"), "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2"),
    (("(-)-camphor", "l-camphor", "(1s)-(-)-camphor"), "CC1(C)[C@H]2CC[C@]1(C)C(=O)C2"),
    (("borneol", "(-)-borneol"), "CC1(C)[C@@H]2CC[C@@]1(C)[C@@H](O)C2"),
    (("alpha-pinene", "α-pinene", "pinene", "(-)-alpha-pinene"), "CC1=CC[C@@H]2C[C@H]1C2(C)C"),
    (("beta-pinene", "β-pinene"), "C=C1CC[C@@H]2C[C@H]1C2(C)C"),
    (("carvone", "(r)-(-)-carvone", "l-carvone", "(r)-carvone"), "CC1=CC[C@H](CC1=O)C(C)=C"),
    (("(s)-carvone", "(s)-(+)-carvone", "d-carvone"), "CC1=CC[C@@H](CC1=O)C(C)=C"),
    (("eucalyptol", "1,8-cineole", "cineole"), "CC1(C)OC2(C)CCC1CC2"),
    (("thujone", "(-)-alpha-thujone", "alpha-thujone", "α-thujone"), "C[C@@H]1C(=O)C[C@]2(C[C@@H]12)C(C)C"),
    (("terpinen-4-ol",), "CC1=CCC(O)(CC1)C(C)C"),
    (("alpha-terpineol", "α-terpineol", "terpineol"), "CC1=CCC(CC1)C(C)(C)O"),
    (("nerolidol",), "CC(C)=CCC/C(C)=C/CCC(C)(O)C=C"),
    (("beta-caryophyllene", "β-caryophyllene", "caryophyllene"), "C=C1CC/C=C(\\C)CC[C@@H]2[C@@H]1CC2(C)C"),
    (("humulene", "α-humulene", "alpha-humulene"), "C/C1=C\\CC(C)(C)/C=C/C/C(C)=C/CC1"),
    (("ocimene", "beta-ocimene"), "C=C/C(C)=C/CC=C(C)C"),
    (("myrcene", "beta-myrcene", "β-myrcene"), "C=CC(=C)CCC=C(C)C"),
    (("vanillin", "4-hydroxy-3-methoxybenzaldehyde"), "COc1cc(C=O)ccc1O"),
    (("ethyl vanillin", "bourbonal"), "CCOc1cc(C=O)ccc1O"),
    (("vanillic acid",), "COc1cc(ccc1O)C(=O)O"),
    (("syringaldehyde",), "COc1cc(C=O)cc(OC)c1O"),
    (("cinnamyl alcohol",), "OC/C=C/c1ccccc1"),
    (("anethole", "trans-anethole"), "C/C=C/c1ccc(OC)cc1"),
    (("estragole", "methyl chavicol"), "C=CCc1ccc(OC)cc1"),
    (("safrole",), "C=CCc1ccc2OCOc2c1"),
    (("isoeugenol",), "C/C=C/c1ccc(O)c(OC)c1"),
    (("benzyl acetate",), "CC(=O)OCc1ccccc1"),
    (("isoamyl acetate", "isopentyl acetate", "banana oil"), "CC(C)CCOC(C)=O"),
    (("ethyl butyrate", "ethyl butanoate"), "CCCC(=O)OCC"),
    (("methyl anthranilate",), "COC(=O)c1ccccc1N"),
    (("methyl cinnamate",), "COC(=O)/C=C/c1ccccc1"),
    (("ethyl cinnamate",), "CCOC(=O)/C=C/c1ccccc1"),
    (("maltol",), "Cc1occc(=O)c1O"),
    (("ethyl maltol",), "CCc1occc(=O)c1O"),
    (("diacetyl", "2,3-butanedione", "butane-2,3-dione"), "CC(=O)C(C)=O"),
    (("acetoin", "3-hydroxy-2-butanone"), "CC(O)C(C)=O"),
    (("raspberry ketone", "4-(4-hydroxyphenyl)-2-butanone", "frambinone"), "CC(=O)CCc1ccc(O)cc1"),
    (("skatole", "3-methylindole"), "Cc1c[nH]c2ccccc12"),
    (("2-acetyl-1-pyrroline",), "CC(=O)C1=NCCC1"),
    (("furaneol", "4-hydroxy-2,5-dimethyl-3(2h)-furanone"), "CC1OC(C)=C(O)C1=O"),
    (("sotolon",), "CC1OC(=O)C(O)=C1C"),
    (("geosmin",), "C[C@@]12CCCC[C@]1(C)[C@@H](O)CCC2"),
    (("2-methylisoborneol", "2-mib"), "CC1(C)[C@@H]2CC[C@@]1(C)[C@](C)(O)C2"),
    (("psoralen",), "O=c1ccc2cc3ccoc3cc2o1"),
    (("8-methoxypsoralen", "methoxsalen", "xanthotoxin", "8-mop"), "COc1c2occc2cc2ccc(=O)oc12"),
    (("salicin",), "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("arbutin", "hydroquinone glucoside"), "OC[C@H]1O[C@@H](Oc2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("esculin",), "OC[C@H]1O[C@@H](Oc2cc3oc(=O)ccc3cc2O)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("aflatoxin b1",), "COC1=C2C3=C(C(=O)CC3)C(=O)OC2=C4[C@@H]5C=CO[C@@H]5OC4=C1"),
    (("ochratoxin a",), "C[C@@H]1Cc2c(Cl)cc(C(=O)N[C@@H](Cc3ccccc3)C(=O)O)c(O)c2C(=O)O1"),
    (("patulin",), "O=C1C=C2OC(O)CC2=CCO1"),
    # ------------------------------------------------------------------ #
    # Phytohormones, agrochemicals and environmental chemicals
    # ------------------------------------------------------------------ #
    (("kinetin", "6-furfurylaminopurine"), "C(Nc1ncnc2[nH]cnc12)c1ccco1"),
    (("zeatin", "trans-zeatin"), "C/C(CO)=C/CNc1ncnc2[nH]cnc12"),
    (("6-benzylaminopurine", "bap", "6-bap", "benzyladenine"), "C(Nc1ncnc2[nH]cnc12)c1ccccc1"),
    (("abscisic acid", "aba", "(s)-abscisic acid", "(+)-abscisic acid"),
     "CC1=CC(=O)CC(C)(C)[C@]1(O)/C=C/C(C)=C\\C(=O)O"),
    (("gibberellic acid", "ga3", "gibberellin a3"),
     "C[C@@]12[C@H](C=C[C@@]3([C@@H]1[C@@H]([C@]45[C@H]3CC[C@](C4)(C(=C)C5)O)C(=O)O)OC2=O)O"),
    (("jasmonic acid",), "CC/C=C\\C[C@H]1[C@@H](CC(=O)O)CCC1=O"),
    (("methyl jasmonate",), "CC/C=C\\C[C@H]1[C@@H](CC(=O)OC)CCC1=O"),
    (("glyphosate", "n-(phosphonomethyl)glycine", "roundup"), "OC(=O)CNCP(=O)(O)O"),
    (("glufosinate", "phosphinothricin"), "CP(=O)(O)CCC(N)C(=O)O"),
    (("atrazine",), "CCNc1nc(Cl)nc(NC(C)C)n1"),
    (("simazine",), "CCNc1nc(Cl)nc(NCC)n1"),
    (("2,4-d", "2,4-dichlorophenoxyacetic acid"), "OC(=O)COc1ccc(Cl)cc1Cl"),
    (("2,4,5-t", "2,4,5-trichlorophenoxyacetic acid"), "OC(=O)COc1cc(Cl)c(Cl)cc1Cl"),
    (("dicamba",), "COc1c(Cl)ccc(Cl)c1C(=O)O"),
    (("ddt", "dichlorodiphenyltrichloroethane", "p,p'-ddt"), "Clc1ccc(cc1)C(c1ccc(Cl)cc1)C(Cl)(Cl)Cl"),
    (("dde", "p,p'-dde"), "Clc1ccc(cc1)C(=C(Cl)Cl)c1ccc(Cl)cc1"),
    (("lindane", "gamma-hexachlorocyclohexane", "γ-hch"), "ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl"),
    (("imidacloprid",), "[O-][N+](=O)/N=C1\\NCCN1Cc1ccc(Cl)nc1"),
    (("thiamethoxam",), "CN1COCN(Cc2cnc(Cl)s2)/C1=N/[N+](=O)[O-]"),
    (("fipronil",), "N#Cc1c(N)n(nc1S(=O)C(F)(F)F)-c1c(Cl)cc(cc1Cl)C(F)(F)F"),
    (("chlorpyrifos",), "CCOP(=S)(OCC)Oc1nc(Cl)c(Cl)cc1Cl"),
    (("malathion",), "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC"),
    (("parathion",), "CCOP(=S)(OCC)Oc1ccc(cc1)[N+](=O)[O-]"),
    (("carbaryl", "sevin"), "CNC(=O)Oc1cccc2ccccc12"),
    (("permethrin",), "CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1"),
    (("deltamethrin",), "CC1(C)[C@@H](C=C(Br)Br)[C@H]1C(=O)O[C@H](C#N)c1cccc(Oc2ccccc2)c1"),
    (("triclosan",), "Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl"),
    (("deet", "n,n-diethyl-m-toluamide", "diethyltoluamide"), "CCN(CC)C(=O)c1cccc(C)c1"),
    (("picaridin", "icaridin"), "CCC(C)OC(=O)N1CCCCC1CCO"),
    (("dehp", "bis(2-ethylhexyl) phthalate", "di(2-ethylhexyl) phthalate", "dioctyl phthalate"),
     "CCCCC(CC)COC(=O)c1ccccc1C(=O)OCC(CC)CCCC"),
    (("dibutyl phthalate", "dbp"), "CCCCOC(=O)c1ccccc1C(=O)OCCCC"),
    (("pfoa", "perfluorooctanoic acid"), "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    (("pfos", "perfluorooctanesulfonic acid"),
     "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    (("tnt", "2,4,6-trinitrotoluene", "trinitrotoluene"),
     "Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]"),
    (("rdx", "cyclotrimethylenetrinitramine", "hexogen"), "[O-][N+](=O)N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1"),
    (("petn", "pentaerythritol tetranitrate"),
     "[O-][N+](=O)OCC(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]"),
    (("nitroglycerin", "glyceryl trinitrate", "gtn", "trinitroglycerin"),
     "[O-][N+](=O)OCC(CO[N+](=O)[O-])O[N+](=O)[O-]"),
    (("pentaerythritol",), "OCC(CO)(CO)CO"),
    (("acesulfame potassium", "acesulfame k", "ace-k"), "CC1=CC(=O)[N-]S(=O)(=O)O1.[K+]"),
    (("sucralose",), "OC[C@H]1O[C@@](CCl)(O[C@H]2O[C@H](CCl)[C@@H](O)[C@H](O)[C@H]2Cl)[C@@H](O)[C@@H]1O"),
    (("sodium cyclamate", "cyclamate"), "[Na+].[O-]S(=O)(=O)NC1CCCCC1"),
    (("sorbic acid", "(2e,4e)-hexa-2,4-dienoic acid"), "C/C=C/C=C/C(=O)O"),
    (("potassium sorbate",), "[K+].C/C=C/C=C/C(=O)[O-]"),
    (("bisphenol s", "bps"), "Oc1ccc(cc1)S(=O)(=O)c1ccc(O)cc1"),
    (("nonylphenol", "4-nonylphenol"), "CCCCCCCCCc1ccc(O)cc1"),
    (("octylphenol", "4-tert-octylphenol"), "CC(C)(C)CC(C)(C)c1ccc(O)cc1"),
    # ------------------------------------------------------------------ #
    # Polyaromatics, materials chemistry and ionic liquids
    # ------------------------------------------------------------------ #
    (("benzo[a]pyrene", "3,4-benzpyrene"), "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34"),
    (("naphthalene", "mothballs"), "c1ccc2ccccc2c1"),
    (("anthracene",), "c1ccc2cc3ccccc3cc2c1"),
    (("phenanthrene",), "c1ccc2c(c1)ccc1ccccc12"),
    (("fluoranthene",), "c1ccc2c(c1)-c1cccc3cccc2c13"),
    (("chrysene",), "c1ccc2c(c1)ccc1c3ccccc3ccc21"),
    (("tetracene", "naphthacene"), "c1ccc2cc3cc4ccccc4cc3cc2c1"),
    (("pentacene",), "c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1"),
    (("coronene",), "c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61"),
    (("fullerene c60", "c60", "buckminsterfullerene"),
     "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"),
    (("biphenyl",), "c1ccc(cc1)-c1ccccc1"),
    (("terphenyl", "p-terphenyl"), "c1ccc(cc1)-c1ccc(cc1)-c1ccccc1"),
    (("triphenylene",), "c1ccc2c(c1)c1ccccc1c1ccccc21"),
    (("azulene",), "c1ccc2cccc2cc1"),
    (("indene",), "C1C=Cc2ccccc12"),
    (("indane", "indan"), "C1Cc2ccccc2C1"),
    (("cyclooctatetraene", "cot"), "C1=CC=CC=CC=C1"),
    (("tetrathiafulvalene", "ttf"), "C1=CSC(=C2SC=CS2)S1"),
    (("tcnq", "7,7,8,8-tetracyanoquinodimethane"), "N#CC(C#N)=C1C=CC(=C(C#N)C#N)C=C1"),
    (("tcne", "tetracyanoethylene"), "N#CC(C#N)=C(C#N)C#N"),
    (("carbazole", "9h-carbazole"), "c1ccc2c(c1)[nH]c1ccccc12"),
    (("dibenzofuran",), "c1ccc2c(c1)oc1ccccc12"),
    (("dibenzothiophene",), "c1ccc2c(c1)sc1ccccc12"),
    (("xanthene",), "C1c2ccccc2Oc2ccccc12"),
    (("thioxanthone",), "O=C1c2ccccc2Sc2ccccc12"),
    (("xanthone",), "O=C1c2ccccc2Oc2ccccc12"),
    (("phenothiazine",), "c1ccc2c(c1)Nc1ccccc1S2"),
    (("phenoxazine",), "c1ccc2c(c1)Nc1ccccc1O2"),
    (("phenazine",), "c1ccc2nc3ccccc3nc2c1"),
    (("acridone", "9(10h)-acridone"), "O=c1c2ccccc2[nH]c2ccccc12"),
    (("2,2'-bithiophene", "bithiophene"), "c1csc(c1)-c1cccs1"),
    (("edot", "3,4-ethylenedioxythiophene"), "c1cc2OCCOc2s1"),
    (("thieno[3,2-b]thiophene",), "c1cc2sccc2s1"),
    (("benzothiadiazole", "2,1,3-benzothiadiazole", "btd"), "c1ccc2nsnc2c1"),
    (("4,7-dibromo-2,1,3-benzothiadiazole",), "Brc1ccc(Br)c2nsnc12"),
    (("ferrocenecarboxylic acid", "ferrocene carboxylic acid"), "[Fe+2].OC(=O)[c-]1cccc1.[cH-]1cccc1"),
    (("ferrocenemethanol",), "[Fe+2].OC[c-]1cccc1.[cH-]1cccc1"),
    (("alq3", "tris(8-hydroxyquinolinato)aluminium", "tris(8-hydroxyquinolinato)aluminum"),
     "[Al+3].[O-]c1cccc2cccnc12.[O-]c1cccc2cccnc12.[O-]c1cccc2cccnc12"),
    (("ionic liquid bmim bf4", "[bmim]bf4", "1-butyl-3-methylimidazolium tetrafluoroborate", "bmim bf4"),
     "CCCCn1cc[n+](C)c1.F[B-](F)(F)F"),
    (("[bmim]pf6", "1-butyl-3-methylimidazolium hexafluorophosphate", "bmim pf6"),
     "CCCCn1cc[n+](C)c1.F[P-](F)(F)(F)(F)F"),
    (("[bmim]cl", "1-butyl-3-methylimidazolium chloride", "bmim cl"), "CCCCn1cc[n+](C)c1.[Cl-]"),
    (("[emim]oac", "1-ethyl-3-methylimidazolium acetate", "emim acetate"), "CCn1cc[n+](C)c1.CC(=O)[O-]"),
    (("[emim]ntf2", "1-ethyl-3-methylimidazolium bis(trifluoromethylsulfonyl)imide", "emim tfsi"),
     "CCn1cc[n+](C)c1.FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F"),
    (("litfsi", "lithium bis(trifluoromethanesulfonyl)imide", "lithium bistriflimide", "lintf2"),
     "[Li+].FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F"),
    (("tf2nh", "bis(trifluoromethanesulfonyl)imide", "triflimide", "hntf2", "bistriflimide"),
     "FC(F)(F)S(=O)(=O)NS(=O)(=O)C(F)(F)F"),
    (("choline chloride/urea des", "reline", "deep eutectic solvent"),
     "C[N+](C)(C)CCO.[Cl-].NC(N)=O.NC(N)=O"),
    # ------------------------------------------------------------------ #
    # Drugs (international nonproprietary names; common trade names as aliases)
    # Stereochemistry only where known; racemates / unknown left flat.
    # ------------------------------------------------------------------ #
    (("abacavir", "ziagen"), "Nc1nc(NC2CC2)c2ncn([C@@H]3C=C[C@H](CO)C3)c2n1"),
    (("acetazolamide", "diamox"), "CC(=O)Nc1nnc(s1)S(N)(=O)=O"),
    (("acyclovir", "aciclovir", "zovirax"), "Nc1nc2n(COCCO)cnc2c(=O)[nH]1"),
    (("albuterol", "salbutamol", "ventolin"), "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1"),
    (("alendronic acid", "alendronate", "fosamax"), "NCCCC(O)(P(=O)(O)O)P(=O)(O)O"),
    (("allopurinol", "zyloprim"), "O=c1[nH]cnc2[nH]ncc12"),
    (("alprazolam", "xanax"), "Cc1nnc2CN=C(c3ccccc3)c3cc(Cl)ccc3-n12"),
    (("amiloride",), "NC(=N)NC(=O)c1nc(Cl)c(N)nc1N"),
    (("amiodarone", "cordarone"), "CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1"),
    (("amitriptyline", "elavil"), "CN(C)CCC=C1c2ccccc2CCc2ccccc12"),
    (("amlodipine", "norvasc"), "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC"),
    (("amlodipine besylate", "amlodipine besilate"),
     "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC.OS(=O)(=O)c1ccccc1"),
    (("anastrozole", "arimidex"), "CC(C)(C#N)c1cc(Cn2cncn2)cc(c1)C(C)(C)C#N"),
    (("apixaban", "eliquis"), "COc1ccc(cc1)-n1nc(C(N)=O)c2CCN(C(=O)c12)c1ccc(cc1)N1CCCCC1=O"),
    (("aripiprazole", "abilify"), "Clc1cccc(N2CCN(CCCCOc3ccc4CCC(=O)Nc4c3)CC2)c1Cl"),
    (("aspirin", "acetylsalicylic acid", "asa", "2-acetoxybenzoic acid"), "CC(=O)Oc1ccccc1C(=O)O"),
    (("atenolol", "tenormin"), "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
    (("atomoxetine", "strattera", "(r)-atomoxetine"), "CNCC[C@@H](Oc1ccccc1C)c1ccccc1"),
    (("azathioprine", "imuran"), "Cn1cnc(c1Sc1ncnc2[nH]cnc12)[N+](=O)[O-]"),
    (("baclofen", "lioresal"), "NCC(CC(=O)O)c1ccc(Cl)cc1"),
    (("benzocaine", "ethyl 4-aminobenzoate"), "CCOC(=O)c1ccc(N)cc1"),
    (("benztropine", "cogentin"), "CN1[C@H]2CC[C@@H]1CC(C2)OC(c1ccccc1)c1ccccc1"),
    (("bisacodyl", "dulcolax"), "CC(=O)Oc1ccc(cc1)C(c1ccccn1)c1ccc(OC(C)=O)cc1"),
    (("bisoprolol", "concor", "zebeta"), "CC(C)NCC(O)COc1ccc(COCCOC(C)C)cc1"),
    (("bortezomib", "velcade"), "CC(C)C[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)c1cnccn1)B(O)O"),
    (("brimonidine", "alphagan"), "Brc1c(NC2=NCCN2)ccc2nccnc12"),
    (("bumetanide", "bumex"), "CCCCNc1cc(cc(c1Oc1ccccc1)S(N)(=O)=O)C(=O)O"),
    (("bupivacaine", "marcaine"), "CCCCN1CCCCC1C(=O)Nc1c(C)cccc1C"),
    (("bupropion", "wellbutrin", "zyban"), "CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1"),
    (("buspirone", "buspar"), "O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(CC1)c1ncccn1"),
    (("carbamazepine", "tegretol"), "NC(=O)N1c2ccccc2C=Cc2ccccc12"),
    (("carbidopa", "lodosyn"), "C[C@](Cc1ccc(O)c(O)c1)(NN)C(=O)O"),
    (("carboplatin", "paraplatin"), "[NH3]->[Pt]1(<-[NH3])OC(=O)C2(CCC2)C(=O)O1"),
    (("carvedilol", "coreg"), "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12"),
    (("celecoxib", "celebrex"), "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F"),
    (("cetirizine", "zyrtec"), "OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1"),
    (("chloroquine",), "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
    (("chlorpheniramine", "chlorphenamine"), "CN(C)CCC(c1ccc(Cl)cc1)c1ccccn1"),
    (("chlorpromazine", "thorazine", "largactil"), "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc12"),
    (("chlorthalidone", "chlortalidone"), "NS(=O)(=O)c1cc(ccc1Cl)C1(O)NC(=O)c2ccccc12"),
    (("cimetidine", "tagamet"), "CNC(=NC#N)NCCSCc1[nH]cnc1C"),
    (("ciprofloxacin", "cipro"), "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O"),
    (("cisplatin", "cis-diamminedichloroplatinum(ii)", "cddp"), "[NH3]->[Pt](Cl)(Cl)<-[NH3]"),
    (("citalopram", "celexa"), "CN(C)CCCC1(OCc2cc(ccc12)C#N)c1ccc(F)cc1"),
    (("clonazepam", "klonopin", "rivotril"), "O=C1CN=C(c2ccccc2Cl)c2cc(ccc2N1)[N+](=O)[O-]"),
    (("clonidine", "catapres"), "Clc1cccc(Cl)c1NC1=NCCN1"),
    (("clopidogrel", "plavix", "(s)-clopidogrel"), "COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1"),
    (("clotrimazole", "canesten"), "Clc1ccccc1C(n1ccnc1)(c1ccccc1)c1ccccc1"),
    (("clozapine", "clozaril"), "CN1CCN(CC1)C1=Nc2cc(Cl)ccc2Nc2ccccc12"),
    (("cyclobenzaprine", "flexeril"), "CN(C)CCC=C1c2ccccc2C=Cc2ccccc12"),
    (("cyclophosphamide", "cytoxan"), "ClCCN(CCCl)P1(=O)NCCCO1"),
    (("cytarabine", "ara-c", "cytosine arabinoside"), "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@@H]2O)c(=O)n1"),
    (("dapagliflozin", "farxiga", "forxiga"),
     "CCOC1=CC=C(C=C1)CC2=C(C=CC(=C2)[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)O)O)Cl"),
    (("dapsone",), "Nc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1"),
    (("desloratadine", "clarinex"), "Clc1ccc2C(=C3CCNCC3)c3ncccc3CCc2c1"),
    (("dexmedetomidine", "precedex"), "Cc1cccc(C)c1[C@H](C)c1c[nH]cn1"),
    (("dextromethorphan", "dxm"), "CN1CC[C@]23CCCC[C@H]2[C@H]1CC4=C3C=C(C=C4)OC"),
    (("diazepam", "valium"), "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12"),
    (("diclofenac", "voltaren"), "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    (("diclofenac sodium",), "[Na+].[O-]C(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    (("diltiazem", "cardizem"), "CC(=O)O[C@@H]1[C@@H](SC2=CC=CC=C2N(C1=O)CCN(C)C)C3=CC=C(C=C3)OC"),
    (("diphenhydramine", "benadryl"), "CN(C)CCOC(c1ccccc1)c1ccccc1"),
    (("disulfiram", "antabuse"), "CCN(CC)C(=S)SSC(=S)N(CC)CC"),
    (("dobutamine",), "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
    (("donepezil", "aricept"), "COc1cc2CC(CC3CCN(Cc4ccccc4)CC3)C(=O)c2cc1OC"),
    (("doxazosin", "cardura"), "COc1cc2nc(nc(N)c2cc1OC)N1CCN(CC1)C(=O)C1COc2ccccc2O1"),
    (("duloxetine", "cymbalta", "(s)-duloxetine"), "CNCC[C@H](Oc1cccc2ccccc12)c1cccs1"),
    (("efavirenz", "sustiva"), "C1CC1C#C[C@]2(C3=C(C=CC(=C3)Cl)NC(=O)O2)C(F)(F)F"),
    (("empagliflozin", "jardiance"),
     "OC[C@H]1O[C@@H](c2ccc(Cl)c(Cc3ccc(OC4CCOC4)cc3)c2)[C@H](O)[C@@H](O)[C@@H]1O"),
    (("emtricitabine", "ftc"), "Nc1nc(=O)n(cc1F)[C@@H]1CS[C@H](CO)O1"),
    (("enalapril", "vasotec"), "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)O"),
    (("erlotinib", "tarceva"), "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC"),
    (("escitalopram", "lexapro", "(s)-citalopram"), "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2cc(C#N)ccc21"),
    (("esomeprazole", "nexium", "(s)-omeprazole"), "COc1ccc2[nH]c(nc2c1)[S@](=O)Cc1ncc(C)c(OC)c1C"),
    (("ethambutol",), "CC[C@@H](CO)NCCN[C@@H](CC)CO"),
    (("ezetimibe", "zetia"), "C1=CC(=CC=C1[C@@H]2[C@H](C(=O)N2C3=CC=C(C=C3)F)CC[C@@H](C4=CC=C(C=C4)F)O)O"),
    (("famotidine", "pepcid"), "NC(N)=Nc1nc(CSCCC(N)=NS(N)(=O)=O)cs1"),
    (("favipiravir", "avigan"), "NC(=O)c1nc(F)cnc1O"),
    (("febuxostat", "uloric"), "Cc1nc(sc1C(=O)O)-c1ccc(OCC(C)C)c(c1)C#N"),
    (("fenofibrate", "tricor"), "CC(C)OC(=O)C(C)(C)Oc1ccc(cc1)C(=O)c1ccc(Cl)cc1"),
    (("fexofenadine", "allegra"), "CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN1CCC(CC1)C(O)(c1ccccc1)c1ccccc1"),
    (("finasteride", "propecia", "proscar"),
     "CC(C)(C)NC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(C=CC(=O)N4)C)C"),
    (("fluconazole", "diflucan"), "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"),
    (("fluoxetine", "prozac"), "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1"),
    (("fluoxetine hydrochloride", "prozac hcl"), "Cl.CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1"),
    (("furosemide", "frusemide", "lasix"), "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl"),
    (("gabapentin", "neurontin"), "NCC1(CC(=O)O)CCCCC1"),
    (("gefitinib", "iressa"), "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"),
    (("gemcitabine", "gemzar"), "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)C2(F)F)c(=O)n1"),
    (("gemfibrozil", "lopid"), "Cc1ccc(C)c(OCCCC(C)(C)C(=O)O)c1"),
    (("glimepiride", "amaryl"), "CCC1=C(C)CN(C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC2CCC(C)CC2)C1=O"),
    (("glipizide", "glucotrol"), "Cc1cnc(cn1)C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"),
    (("glyburide", "glibenclamide", "diabeta"), "COc1ccc(Cl)cc1C(=O)NCCc1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1"),
    (("guaifenesin", "mucinex"), "COc1ccccc1OCC(O)CO"),
    (("guanfacine", "intuniv"), "NC(=N)NC(=O)Cc1c(Cl)cccc1Cl"),
    (("haloperidol", "haldol"), "OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1"),
    (("hydralazine",), "NNc1nncc2ccccc12"),
    (("hydrochlorothiazide", "hctz", "microzide"), "NS(=O)(=O)c1cc2NCNS(=O)(=O)c2cc1Cl"),
    (("hydroxychloroquine", "plaquenil"), "CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
    (("hydroxychloroquine sulfate",), "OS(=O)(=O)O.CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
    (("hydroxyzine", "atarax", "vistaril"), "OCCOCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1"),
    (("ibrutinib", "imbruvica"), "C=CC(=O)N1CCC[C@H](C1)n1nc(-c2ccc(Oc3ccccc3)cc2)c2c(N)ncnc12"),
    (("ibuprofen", "advil", "motrin", "brufen", "2-(4-isobutylphenyl)propanoic acid"),
     "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    (("(s)-ibuprofen", "s-ibuprofen", "dexibuprofen", "(s)-(+)-ibuprofen"), "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
    (("(r)-ibuprofen", "r-ibuprofen", "(r)-(-)-ibuprofen"), "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"),
    (("imatinib", "gleevec", "glivec"), "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1"),
    (("imatinib mesylate", "imatinib mesilate"),
     "CS(=O)(=O)O.Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1"),
    (("imipramine", "tofranil"), "CN(C)CCCN1c2ccccc2CCc2ccccc12"),
    (("indomethacin", "indometacin", "indocin"), "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1"),
    (("irbesartan", "avapro"), "CCCCC1=NC2(CCCC2)C(=O)N1Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1"),
    (("isoniazid", "isonicotinic acid hydrazide", "inh"), "NNC(=O)c1ccncc1"),
    (("ketamine", "ketalar"), "CNC1(CCCCC1=O)c1ccccc1Cl"),
    (("esketamine", "(s)-ketamine", "spravato"), "CN[C@]1(c2ccccc2Cl)CCCCC1=O"),
    (("ketoprofen",), "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1"),
    (("ketorolac", "toradol"), "OC(=O)C1CCn2c1ccc2C(=O)c1ccccc1"),
    (("labetalol", "trandate"), "CC(CCc1ccccc1)NCC(O)c1ccc(O)c(C(N)=O)c1"),
    (("lamivudine", "3tc", "epivir"), "Nc1ccn([C@@H]2CS[C@H](CO)O2)c(=O)n1"),
    (("lamotrigine", "lamictal"), "Nc1nnc(c(N)n1)-c1cccc(Cl)c1Cl"),
    (("lansoprazole", "prevacid"), "Cc1c(OCC(F)(F)F)ccnc1CS(=O)c1nc2ccccc2[nH]1"),
    (("latanoprost", "xalatan"),
     "CC(C)OC(=O)CCC/C=C\\C[C@H]1[C@@H](O)C[C@@H](O)[C@@H]1CC[C@@H](O)CCc1ccccc1"),
    (("leflunomide", "arava"), "Cc1oncc1C(=O)Nc1ccc(cc1)C(F)(F)F"),
    (("lenalidomide", "revlimid"), "Nc1cccc2C(=O)N(Cc12)C1CCC(=O)NC1=O"),
    (("letrozole", "femara"), "N#Cc1ccc(cc1)C(c1ccc(cc1)C#N)n1cncn1"),
    (("levetiracetam", "keppra"), "CC[C@H](N1CCCC1=O)C(N)=O"),
    (("levofloxacin", "levaquin"), "C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(cn1c23)C(=O)O"),
    (("levonorgestrel", "plan b"), "CC[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC1=CC(=O)CC[C@H]31"),
    (("lidocaine", "lignocaine", "xylocaine"), "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    (("lidocaine hydrochloride",), "Cl.CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    (("linezolid", "zyvox"), "CC(=O)NC[C@H]1CN(C(=O)O1)c1ccc(N2CCOCC2)c(F)c1"),
    (("lisinopril", "zestril", "prinivil"), "NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O"),
    (("loperamide", "imodium"), "CN(C)C(=O)C(CCN1CCC(O)(CC1)c1ccc(Cl)cc1)(c1ccccc1)c1ccccc1"),
    (("loratadine", "claritin"), "CCOC(=O)N1CCC(CC1)=C1c2ccc(Cl)cc2CCc2cccnc12"),
    (("lorazepam", "ativan"), "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O"),
    (("losartan", "cozaar"), "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1"),
    (("losartan potassium",), "[K+].CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nn[n-]n1"),
    (("meclizine", "meclozine", "antivert"), "Cc1cccc(CN2CCN(CC2)C(c2ccccc2)c2ccc(Cl)cc2)c1"),
    (("meloxicam", "mobic"), "Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1"),
    (("mesalamine", "mesalazine", "5-aminosalicylic acid", "5-asa"), "Nc1ccc(O)c(c1)C(=O)O"),
    (("metformin", "glucophage"), "CN(C)C(=N)NC(N)=N"),
    (("metformin hydrochloride", "metformin hcl"), "Cl.CN(C)C(=N)NC(N)=N"),
    (("methimazole", "thiamazole", "tapazole"), "Cn1cc[nH]c1=S"),
    (("methocarbamol", "robaxin"), "COc1ccccc1OCC(O)COC(N)=O"),
    (("methylphenidate", "ritalin"), "COC(=O)C(c1ccccc1)C1CCCCN1"),
    (("metoclopramide", "reglan"), "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC"),
    (("metoprolol", "lopressor", "toprol"), "COCCc1ccc(OCC(O)CNC(C)C)cc1"),
    (("metronidazole", "flagyl"), "Cc1ncc(n1CCO)[N+](=O)[O-]"),
    (("midazolam", "versed"), "Cc1ncc2CN=C(c3ccccc3F)c3cc(Cl)ccc3-n12"),
    (("mirtazapine", "remeron"), "CN1CCN2c3ncccc3Cc3ccccc3C2C1"),
    (("misoprostol", "cytotec"), "CCCCC(C)(O)C/C=C/[C@H]1[C@H](O)CC(=O)[C@@H]1CCCCCCC(=O)OC"),
    (("modafinil", "provigil"), "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1"),
    (("montelukast", "singulair"),
     "CC(C)(C1=CC=CC=C1CC[C@H](C2=CC=CC(=C2)/C=C/C3=NC4=C(C=CC(=C4)Cl)C=C3)SCC5(CC5)CC(=O)O)O"),
    (("mycophenolic acid",), "COc1c(C)c2COC(=O)c2c(O)c1C/C=C(\\C)CCC(=O)O"),
    (("naproxen", "aleve", "naprosyn", "(s)-naproxen"), "COc1ccc2cc(ccc2c1)[C@H](C)C(=O)O"),
    (("naproxen sodium",), "[Na+].COc1ccc2cc(ccc2c1)[C@H](C)C(=O)[O-]"),
    (("nifedipine", "adalat", "procardia"), "COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC"),
    (("nitrofurantoin", "macrobid"), "O=C1CN(N=Cc2ccc(o2)[N+](=O)[O-])C(=O)N1"),
    (("norethindrone", "norethisterone"), "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC1=CC(=O)CC[C@H]31"),
    (("nortriptyline", "pamelor"), "CNCCC=C1c2ccccc2CCc2ccccc12"),
    (("olanzapine", "zyprexa"), "Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1"),
    (("olaparib", "lynparza"), "C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C45)F"),
    (("omeprazole", "prilosec", "losec"), "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C"),
    (("ondansetron", "zofran"), "Cc1nccn1CC1CCc2c(C1=O)c1ccccc1n2C"),
    (("oseltamivir", "tamiflu"), "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(C)=O)N)C(=O)OCC"),
    (("oseltamivir phosphate",), "OP(=O)(O)O.CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(C)=O)N)C(=O)OCC"),
    (("oxcarbazepine", "trileptal"), "NC(=O)N1c2ccccc2CC(=O)c2ccccc12"),
    (("oxybutynin", "ditropan"), "CCN(CC)CC#CCOC(=O)C(O)(C1CCCCC1)c1ccccc1"),
    (("pantoprazole", "protonix"), "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC"),
    (("paracetamol", "acetaminophen", "tylenol", "panadol", "apap", "n-acetyl-p-aminophenol", "4-acetamidophenol"),
     "CC(=O)Nc1ccc(O)cc1"),
    (("paroxetine", "paxil", "seroxat"), "C1CNC[C@H]([C@@H]1C2=CC=C(C=C2)F)COC3=CC4=C(C=C3)OCO4"),
    (("phenobarbital", "phenobarbitone", "luminal"), "CCC1(C(=O)NC(=O)NC1=O)c1ccccc1"),
    (("phentermine", "adipex"), "CC(C)(N)Cc1ccccc1"),
    (("phenytoin", "dilantin"), "O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1"),
    (("pioglitazone", "actos"), "CCc1ccc(OCCc2ccc(CC3SC(=O)NC3=O)cc2)nc1"),
    (("pramipexole", "mirapex"), "CCCN[C@H]1CCc2nc(N)sc2C1"),
    (("praziquantel", "biltricide"), "O=C(C1CCCCC1)N1CC(=O)N2CCc3ccccc3C2C1"),
    (("prazosin", "minipress"), "COc1cc2nc(nc(N)c2cc1OC)N1CCN(CC1)C(=O)c1ccco1"),
    (("pregabalin", "lyrica", "(s)-pregabalin"), "CC(C)C[C@H](CN)CC(=O)O"),
    (("primaquine",), "COc1cc(NC(C)CCCN)c2ncccc2c1"),
    (("procaine", "novocaine"), "CCN(CC)CCOC(=O)c1ccc(N)cc1"),
    (("promethazine", "phenergan"), "CC(CN1c2ccccc2Sc2ccccc12)N(C)C"),
    (("propofol", "diprivan", "2,6-diisopropylphenol"), "CC(C)c1cccc(C(C)C)c1O"),
    (("propranolol", "inderal"), "CC(C)NCC(O)COc1cccc2ccccc12"),
    (("propylthiouracil", "ptu"), "CCCc1cc(=O)[nH]c(=S)[nH]1"),
    (("pyrazinamide",), "NC(=O)c1cnccn1"),
    (("quetiapine", "seroquel"), "OCCOCCN1CCN(CC1)C1=Nc2ccccc2Sc2ccccc12"),
    (("ramipril", "altace"), "CCOC(=O)[C@H](CCC1=CC=CC=C1)N[C@@H](C)C(=O)N2[C@@H]3CCC[C@@H]3C[C@H]2C(=O)O"),
    (("ranitidine", "zantac"), "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1"),
    (("rasagiline", "azilect"), "C#CCN[C@@H]1CCc2ccccc12"),
    (("ribavirin", "virazole"), "NC(=O)c1ncn(n1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"),
    (("risedronic acid", "risedronate", "actonel"), "OC(Cc1cccnc1)(P(=O)(O)O)P(=O)(O)O"),
    (("risperidone", "risperdal"), "Cc1nc2CCCCn2c(=O)c1CCN1CCC(CC1)c1noc2cc(F)ccc12"),
    (("rivaroxaban", "xarelto"), "O=C1COCCN1c1ccc(cc1)N1C[C@H](CNC(=O)c2ccc(Cl)s2)OC1=O"),
    (("rivastigmine", "exelon"), "CCN(C)C(=O)Oc1cccc(c1)[C@H](C)N(C)C"),
    (("rizatriptan", "maxalt"), "CN(C)CCc1c[nH]c2ccc(Cn3cncn3)cc12"),
    (("ropinirole", "requip"), "CCCN(CCC)CCc1cccc2NC(=O)Cc12"),
    (("salmeterol", "serevent"), "CCCCCCOCCCCNCC(O)c1ccc(O)c(CO)c1"),
    (("selegiline", "deprenyl", "eldepryl", "l-deprenyl"), "C#CCN(C)[C@H](C)Cc1ccccc1"),
    (("sertraline", "zoloft"), "CN[C@H]1CC[C@H](C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl"),
    (("sertraline hydrochloride",), "Cl.CN[C@H]1CC[C@H](C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl"),
    (("sildenafil", "viagra"), "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"),
    (("sildenafil citrate",),
     "OC(=O)CC(O)(CC(=O)O)C(=O)O.CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"),
    (("sitagliptin", "januvia"), "N[C@@H](CC(=O)N1CCn2c(C1)nnc2C(F)(F)F)Cc1cc(F)c(F)cc1F"),
    (("sorafenib", "nexavar"), "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(c3)C(F)(F)F)cc2)ccn1"),
    (("spironolactone", "aldactone"),
     "CC(=O)S[C@H]1CC2=CC(=O)CC[C@@]2([C@H]3[C@@H]1[C@@H]4CC[C@]5([C@@]4(CC3)C)CCC(=O)O5)C"),
    (("succinylcholine", "suxamethonium", "succinylcholine chloride"),
     "C[N+](C)(C)CCOC(=O)CCC(=O)OCC[N+](C)(C)C.[Cl-].[Cl-]"),
    (("sulfamethoxazole", "smx"), "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"),
    (("sulfasalazine", "azulfidine"), "OC(=O)c1cc(ccc1O)N=Nc1ccc(cc1)S(=O)(=O)Nc1ccccn1"),
    (("sumatriptan", "imitrex"), "CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1"),
    (("sunitinib", "sutent"), "CCN(CC)CCNC(=O)c1c(C)[nH]c(/C=C2\\C(=O)Nc3ccc(F)cc23)c1C"),
    (("tadalafil", "cialis"), "CN1CC(=O)N2[C@H](C1=O)CC3=C([C@H]2C4=CC5=C(C=C4)OCO5)NC6=CC=CC=C36"),
    (("tamoxifen", "nolvadex", "(z)-tamoxifen"), "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
    (("tamsulosin", "flomax"), "CCOc1ccccc1OCCN[C@@H](C)Cc1ccc(OC)c(c1)S(N)(=O)=O"),
    (("telmisartan", "micardis"), "CCCc1nc2c(C)cc(cc2n1Cc1ccc(cc1)-c1ccccc1C(=O)O)-c1nc2ccccc2n1C"),
    (("temazepam", "restoril"), "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc12"),
    (("temozolomide", "temodar"), "CN1N=NC2=C(N=CN2C1=O)C(N)=O"),
    (("tenofovir", "pmpa"), "C[C@H](Cn1cnc2c(N)ncnc12)OCP(=O)(O)O"),
    (("terbinafine", "lamisil"), "CN(C/C=C/C#CC(C)(C)C)Cc1cccc2ccccc12"),
    (("thalidomide",), "O=C1CCC(N2C(=O)c3ccccc3C2=O)C1=O"),
    (("timolol",), "CC(C)(C)NCC(O)COc1nsnc1N1CCOCC1"),
    (("tizanidine", "zanaflex"), "Clc1ccc2nsnc2c1NC1=NCCN1"),
    (("tofacitinib", "xeljanz"), "C[C@@H]1CCN(C[C@@H]1N(C)c1ncnc2[nH]ccc12)C(=O)CC#N"),
    (("topiramate", "topamax"), "CC1(C)O[C@@H]2CO[C@@]3(COS(N)(=O)=O)OC(C)(C)O[C@H]3[C@@H]2O1"),
    (("torsemide", "torasemide", "demadex"), "Cc1cccc(Nc2ccncc2S(=O)(=O)NC(=O)NC(C)C)c1"),
    (("trazodone", "desyrel"), "Clc1cccc(c1)N1CCN(CCCn2nc3ccccn3c2=O)CC1"),
    (("triamterene", "dyrenium"), "Nc1nc(N)c2nc(c(N)nc2n1)-c1ccccc1"),
    (("trimethoprim",), "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC"),
    (("valacyclovir", "valaciclovir", "valtrex"), "CC(C)[C@H](N)C(=O)OCCOCn1cnc2c(=O)[nH]c(N)nc12"),
    (("valproic acid", "2-propylpentanoic acid", "depakene"), "CCCC(CCC)C(=O)O"),
    (("sodium valproate", "valproate sodium", "depakote"), "[Na+].CCCC(CCC)C(=O)[O-]"),
    (("valsartan", "diovan"), "CCCCC(=O)N(Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1)[C@@H](C(C)C)C(=O)O"),
    (("varenicline", "chantix", "champix"), "C1C2CNCC1C3=CC4=NC=CN=C4C=C23"),
    (("venlafaxine", "effexor"), "COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1"),
    (("verapamil", "calan", "isoptin"), "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC"),
    (("voriconazole", "vfend"), "C[C@@H](c1ncncc1F)[C@](O)(Cn1cncn1)c1ccc(F)cc1F"),
    (("warfarin", "coumadin"), "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    (("warfarin sodium",), "[Na+].CC(=O)CC(c1ccccc1)c1c([O-])c2ccccc2oc1=O"),
    (("zidovudine", "azt", "azidothymidine", "retrovir"),
     "Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O"),
    (("zoledronic acid", "zoledronate", "zometa", "reclast"), "OC(Cn1ccnc1)(P(=O)(O)O)P(=O)(O)O"),
    (("zolpidem", "ambien"), "Cc1ccc(cc1)-c1nc2ccc(C)cn2c1CC(=O)N(C)C"),
    (("albendazole",), "CCCSc1ccc2[nH]c(NC(=O)OC)nc2c1"),
    (("mebendazole",), "COC(=O)Nc1nc2cc(ccc2[nH]1)C(=O)c1ccccc1"),
    (("acetylsalicylic acid lysine salt", "lysine acetylsalicylate"),
     "CC(=O)Oc1ccccc1C(=O)[O-].NCCCC[C@H]([NH3+])C(=O)O"),
    (("ferrous gluconate", "iron(ii) gluconate"),
     "[Fe+2].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]"),
    (("calcium gluconate",),
     "[Ca+2].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]"),
    (("magnesium citrate",),
     "[Mg+2].[Mg+2].[Mg+2].[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-].[O-]C(=O)CC(O)(CC(=O)[O-])C(=O)[O-]"),
    (("potassium gluconate",), "[K+].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)[O-]"),
    (("sodium lactate",), "[Na+].CC(O)C(=O)[O-]"),
    (("dimethyl fumarate", "tecfidera"), "COC(=O)/C=C/C(=O)OC"),
    (("hydroxyurea", "hydroxycarbamide"), "NC(=O)NO"),
    (("busulfan",), "CS(=O)(=O)OCCCCOS(C)(=O)=O"),
    (("melphalan",), "N[C@@H](Cc1ccc(cc1)N(CCCl)CCCl)C(=O)O"),
    (("chlorambucil",), "OC(=O)CCCc1ccc(cc1)N(CCCl)CCCl"),
    (("mechlorethamine", "nitrogen mustard", "hn2"), "CN(CCCl)CCCl"),
    (("procarbazine",), "CNNCc1ccc(cc1)C(=O)NC(C)C"),
    (("dacarbazine",), "CN(C)N=Nc1[nH]cnc1C(N)=O"),
    (("fludarabine",), "Nc1nc(F)nc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@@H]1O"),
    (("cladribine",), "Nc1nc(Cl)nc2n(cnc12)[C@H]1C[C@H](O)[C@@H](CO)O1"),
    (("decitabine", "5-aza-2'-deoxycytidine"), "Nc1ncn(c(=O)n1)[C@H]1C[C@H](O)[C@@H](CO)O1"),
    (("pemetrexed",), "Nc1nc(=O)c2c([nH]1)[nH]cc2CCc1ccc(cc1)C(=O)N[C@@H](CCC(=O)O)C(=O)O"),
    (("vorinostat", "saha", "suberoylanilide hydroxamic acid"), "ONC(=O)CCCCCCC(=O)Nc1ccccc1"),
    (("ruxolitinib",), "N#CC[C@H](C1CCCC1)n1cc(cn1)-c1ncnc2[nH]ccc12"),
    (("baricitinib",), "CCS(=O)(=O)N1CC(CC#N)(C1)n1cc(cn1)-c1ncnc2[nH]ccc12"),
    (("osimertinib",), "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(n1)-c1cn(C)c2ccccc12"),
    (("dasatinib",), "Cc1nc(Nc2ncc(s2)C(=O)Nc2c(C)cccc2Cl)cc(n1)N1CCN(CCO)CC1"),
    (("nilotinib",), "Cc1cn(cn1)-c1cc(NC(=O)c2ccc(C)c(Nc3nccc(n3)-c3cccnc3)c2)cc(c1)C(F)(F)F"),
    (("crizotinib",), "C[C@@H](Oc1cc(cnc1N)-c1cnn(C2CCNCC2)c1)c1c(Cl)ccc(F)c1Cl"),
    (("vemurafenib",), "CCCS(=O)(=O)Nc1ccc(F)c(c1F)C(=O)c1c[nH]c2ncc(cc12)-c1ccc(Cl)cc1"),
    (("palbociclib",), "CC(=O)c1c(C)c2cnc(Nc3ccc(cn3)N3CCNCC3)nc2n(C2CCCC2)c1=O"),
    (("enzalutamide",), "CNC(=O)c1ccc(cc1F)N1C(=S)N(C(=O)C1(C)C)c1ccc(C#N)c(c1)C(F)(F)F"),
    (("abiraterone",), "C[C@]12CC[C@H]3[C@H]([C@@H]1CC=C2c1cccnc1)CC=C1C[C@@H](O)CC[C@]31C"),
    (("bicalutamide",), "CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(c1)C(F)(F)F"),
    (("flutamide",), "CC(C)C(=O)Nc1ccc(cc1C(F)(F)F)[N+](=O)[O-]"),
    (("raloxifene",), "Oc1ccc(cc1)-c1sc2cc(O)ccc2c1C(=O)c1ccc(OCCN2CCCCC2)cc1"),
    (("exemestane",), "C=C1C[C@H]2[C@@H](CC[C@]3(C)C(=O)CC[C@@H]23)[C@@]2(C)C=CC(=O)C=C12"),
    (("megestrol acetate",), "CC(=O)O[C@@]1(CC[C@H]2[C@@H]3CCC4=CC(=O)C=C(C)[C@@H]4[C@@H]3CC[C@@]21C)C(C)=O"),
    (("medroxyprogesterone acetate", "depo-provera"),
     "C[C@H]1C[C@H]2[C@@H]3CC[C@]([C@]3(CC[C@@H]2[C@@]4(C1=CC(=O)CC4)C)C)(C(=O)C)OC(=O)C"),
    (("mifepristone", "ru-486"),
     "CC#C[C@]1(O)CC[C@H]2[C@@H]3CCC4=CC(=O)CCC4=C3[C@@H](C[C@@]21C)c1ccc(cc1)N(C)C"),
    (("clomiphene", "clomifene"), "CCN(CC)CCOc1ccc(cc1)C(=C(Cl)c1ccccc1)c1ccccc1"),
    (("tolterodine",), "CC(C)N(CC[C@H](c1ccccc1)c1cc(C)ccc1O)C(C)C"),
    (("mirabegron",), "Nc1nc(CC(=O)Nc2ccc(CCNC[C@H](O)c3ccccc3)cc2)cs1"),
    (("glycopyrrolate", "glycopyrronium bromide"), "C[N+]1(C)CCC(C1)OC(=O)C(O)(C1CCCC1)c1ccccc1.[Br-]"),
    (("ipratropium", "ipratropium bromide", "atrovent"),
     "CC(C)[N+]1(C)[C@H]2CC[C@@H]1CC(C2)OC(=O)C(CO)c1ccccc1.[Br-]"),
    (("formoterol",), "COc1ccc(cc1)CC(C)NCC(O)c1ccc(O)c(NC=O)c1"),
    (("cromolyn", "cromoglicic acid", "sodium cromoglycate acid"),
     "OC(=O)c1cc(=O)c2cccc(OCC(O)COc3cccc4oc(cc(=O)c34)C(=O)O)c2o1"),
    (("miconazole",), "Clc1ccc(COC(Cn2ccnc2)c2ccc(Cl)cc2Cl)c(Cl)c1"),
    (("griseofulvin",), "C[C@@H]1CC(=O)C=C([C@]12C(=O)C3=C(O2)C(=C(C=C3OC)OC)Cl)OC"),
    (("moxifloxacin",), "COc1c(N2C[C@@H]3CCCN[C@@H]3C2)c(F)cc2c(=O)c(cn(C3CC3)c12)C(=O)O"),
    (("norfloxacin",), "CCn1cc(C(=O)O)c(=O)c2cc(F)c(N3CCNCC3)cc12"),
    (("ofloxacin",), "CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(cn1c23)C(=O)O"),
    (("nalidixic acid",), "CCn1cc(C(=O)O)c(=O)c2ccc(C)nc12"),
    (("fosfomycin",), "C[C@@H]1O[C@@H]1P(=O)(O)O"),
    (("nitrofurazone", "nitrofural"), "NC(=O)N/N=C/c1ccc(o1)[N+](=O)[O-]"),
    (("sulfadiazine",), "Nc1ccc(cc1)S(=O)(=O)Nc1ncccn1"),
    (("silver sulfadiazine",), "[Ag+].Nc1ccc(cc1)S(=O)(=O)[N-]c1ncccn1"),
    (("sulfanilamide",), "Nc1ccc(cc1)S(N)(=O)=O"),
    (("sulfapyridine",), "Nc1ccc(cc1)S(=O)(=O)Nc1ccccn1"),
    (("prontosil",), "Nc1ccc(N=Nc2ccc(cc2)S(N)(=O)=O)c(N)c1"),
    (("ethionamide",), "CCc1cc(C(N)=S)ccn1"),
    (("rimantadine",), "CC(N)C12CC3CC(CC(C3)C1)C2"),
    (("ganciclovir", "cytovene"), "Nc1nc2n(COC(CO)CO)cnc2c(=O)[nH]1"),
    (("valganciclovir",), "CC(C)[C@H](N)C(=O)OCC(CO)OCn1cnc2c(=O)[nH]c(N)nc12"),
    (("famciclovir",), "CC(=O)OCC(CCn1cnc2c(N)ncnc12)COC(C)=O"),
    (("penciclovir",), "Nc1nc2n(CCC(CO)CO)cnc2c(=O)[nH]1"),
    (("cidofovir",), "Nc1ccn(C[C@@H](CO)OCP(=O)(O)O)c(=O)n1"),
    (("foscarnet", "phosphonoformic acid"), "OC(=O)P(=O)(O)O"),
    (("didanosine", "ddi"), "O=c1[nH]cnc2n(cnc12)[C@H]1CC[C@@H](CO)O1"),
    (("stavudine", "d4t"), "Cc1cn([C@H]2C=C[C@@H](CO)O2)c(=O)[nH]c1=O"),
    (("abacavir sulfate",),
     "OS(=O)(=O)O.Nc1nc(NC2CC2)c2ncn([C@@H]3C=C[C@H](CO)C3)c2n1.Nc1nc(NC2CC2)c2ncn([C@@H]3C=C[C@H](CO)C3)c2n1"),
    (("nevirapine",), "Cc1ccnc2N(C3CC3)c3ncccc3C(=O)Nc12"),
    (("raltegravir",), "Cc1nnc(o1)C(=O)NC(C)(C)c1nc(C(=O)NCc2ccc(F)cc2)c(O)c(=O)n1C"),
    (("dolutegravir",), "C[C@@H]1CCO[C@@H]2Cn3cc(C(=O)NCc4ccc(F)cc4F)c(=O)c(O)c3C(=O)N12"),
    (("atovaquone",), "O=C1C(O)=C(C2CCC(CC2)c2ccc(Cl)cc2)C(=O)c2ccccc12"),
    (("proguanil",), "CC(C)NC(=N)NC(=N)Nc1ccc(Cl)cc1"),
    (("pyrimethamine",), "CCc1nc(N)nc(N)c1-c1ccc(Cl)cc1"),
    (("miltefosine",), "CCCCCCCCCCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C"),
    (("benznidazole",), "O=C(Cn1ccnc1[N+](=O)[O-])NCc1ccccc1"),
    (("pentamidine",), "NC(=N)c1ccc(OCCCCCOc2ccc(cc2)C(N)=N)cc1"),
)
COMMON_NAMES: Dict[str, str] = {alias: smi for aliases, smi in _COMMON for alias in aliases}

_DASHES = str.maketrans({
    "‐": "-", "‑": "-", "–": "-", "—": "-", "−": "-", "’": "'", "‘": "'", "`": "'", "ʼ": "'", "´": "'",
    "ü": "u", "ö": "o", "ä": "a", "é": "e", "è": "e", "ë": "e", "å": "a", "ß": "ss", " ": " ",
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon", "κ": "kappa", "λ": "lambda",
    "μ": "mu", "µ": "mu", "ω": "omega", "σ": "sigma", "τ": "tau", "±": "+/-", "×": ".", "·": ".",
    "•": ".", "⋅": ".", "∙": ".", "∗": ".", "*": ".",
})
_SALT_DOT_RE = re.compile(r"\s*\.\s*")
_TIMES_RE = re.compile(r"(?<=\S)\s+x\s+(?=\S)")   # "NiCl2 x 6H2O" -> "nicl2.6h2o"
_LOOSE_RE = re.compile(r"[\s\-.'’,]")
_SPELLINGS = (("sulph", "sulf"), ("aluminium", "aluminum"), ("caesium", "cesium"), ("oestr", "estr"),
              ("haem", "hem"), ("amphetamine", "amfetamine"), ("ketal", "acetal"))
# Hydrate / solvate suffixes are dropped on lookup: the water or solvent of
# crystallisation is not part of the structure ("MgSO4·7H2O" -> MgSO4,
# "PdCl2(dppf)·CH2Cl2" -> PdCl2(dppf)).  Salts (".HCl", ".TFA") are kept.
_SOLVENTS = r"(?:h2o|d2o|hydrate|ch2cl2|dcm|thf|et2o|ether|meoh|etoh|dmf|dmso|dioxane|toluene|acoh|hoac|mecn|nh3)"
_HYDRATE_RE = re.compile(
    r"(?:[\s.]*\(?(?:mono|hemi|sesqui|di|tri|tetra|penta|hexa|hepta|octa|nona|deca)?hydrate\)?"
    r"|\.\s*(?:\d+(?:\.\d+)?|x|n)?\s*" + _SOLVENTS + r"|\s+(?:\d+(?:\.\d+)?|x|n)\s*" + _SOLVENTS
    + r"|\s*\(anhydrous\)|\s+anhydrous)$"
)
_LEADING_QUALIFIER_RE = re.compile(r"^(?:anhydrous|dry|aqueous|aq\.?|conc\.?|concentrated|\d+(?:\.\d+)?\s*%|\d+(?:\.\d+)?\s*m)\s+")
_TRAILING_QUALIFIER_RE = re.compile(
    r"\s*\((?=[^)]*(?:\d|\baq\b|aqueous|anhydrous|dry|solid|solution|\bin\s))[^)]*\)$|\s+\d+(?:\.\d+)?\s*%$"
)
_STEREO_PREFIX_RE = re.compile(
    r"^(?:(?P<rac>dl-|rac-|racemic\s+|\(\+/-\)-|\(\+-\)-|\(rs\)-|\(r/s\)-|\(r,s\)-)"
    r"|(?P<ld>[ld]-)"
    r"|(?P<rs>\((?P<rslabel>[rs])\)-)"
    r"|(?P<multi>\((?P<multilabels>\d*[rs](?:,\s*\d*[rs])+)\)-)"
    r"|(?P<sign>\([+-]\)-|[+-]-)"
    r"|(?P<geo>cis-|trans-|\([ez]\)-))"
)
# Element-symbol collisions that only differ by case (resolved before lower-casing).
_CASE_SENSITIVE: Dict[str, str] = {
    "CoCl2": "[Co+2].[Cl-].[Cl-]",
    "COCl2": "O=C(Cl)Cl",
    "Co": "[Co]",
    "CO": "[C-]#[O+]",
    "In": "[In]",
    "NO": "[N]=O",
    "Nb": "[Nb]",
    "NB": "[Nb]",
}
_LOOSE_INDEX: Optional[Dict[str, str]] = None


def _normalize_name(name: str) -> str:
    """Lower-case, unify dashes/quotes/accents/Greek letters, collapse whitespace, ``·`` -> ``.``."""
    s = name.strip().strip("\"'").lower().translate(_DASHES)
    s = _TIMES_RE.sub(".", s)
    s = _SALT_DOT_RE.sub(".", s)
    s = re.sub(r"\s+", " ", s).strip(" .;,")
    return s


def _loose_key(key: str) -> str:
    """Key without separators and with British/German spellings unified ("n-buli" -> "nbuli")."""
    s = _LOOSE_RE.sub("", key)
    for a, b in _SPELLINGS:
        s = s.replace(a, b)
    return s


def _loose_index() -> Dict[str, str]:
    """Loose-key -> SMILES table, built on first use so that importing stays cheap.

    Keys shared by different structures (e.g. "pd-c" vs "pdc") are ambiguous and
    are left out, so only an exact alias resolves them.
    """
    global _LOOSE_INDEX
    if _LOOSE_INDEX is None:
        table: Dict[str, str] = {}
        ambiguous = set()
        for aliases, smi in _COMMON:
            for alias in aliases:
                key = _loose_key(_normalize_name(alias))
                if table.setdefault(key, smi) != smi:
                    ambiguous.add(key)
        for key in ambiguous:
            del table[key]
        _LOOSE_INDEX = table
    return _LOOSE_INDEX


def _lookup_exact(key: str) -> Optional[str]:
    smi = COMMON_NAMES.get(key)
    if smi is None and key:
        smi = _loose_index().get(_loose_key(key))
    return smi


def _without_tetrahedral(smiles: str) -> str:
    """SMILES with tetrahedral stereo removed (double-bond geometry is kept)."""
    if "@" not in smiles:
        return smiles
    try:
        Chem = _rd().Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles
        for atom in mol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def _cip_labels(mol: Any) -> Dict[int, str]:
    """Atom index -> CIP label for every assigned tetrahedral centre."""
    try:
        from rdkit.Chem import rdCIPLabeler

        rdCIPLabeler.AssignCIPLabels(mol)
        return {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}
    except Exception:  # pragma: no cover - old RDKit
        Chem = _rd().Chem
        return dict(Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False))


def _with_single_center(smiles: str, label: str) -> Optional[str]:
    """Set the one stereocentre of ``smiles`` to CIP ``label`` (R/S); None when ambiguous."""
    try:
        Chem = _rd().Chem
        from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        if not centres:
            return smiles
        if len(centres) != 1:
            return None
        if _cip_labels(mol).get(centres[0][0]) == label:
            return smiles
        for atom in mol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        opts = StereoEnumerationOptions(onlyUnassigned=True, unique=True)
        with _quiet():
            for iso in EnumerateStereoisomers(mol, options=opts):
                if _cip_labels(iso).get(centres[0][0]) == label:
                    return Chem.MolToSmiles(iso)
    except Exception:
        return None
    return None


def _resolve_offline(key: str, depth: int = 0) -> Optional[str]:
    """Dictionary lookup with hydrate suffixes, qualifiers and stereo prefixes handled."""
    if not key or depth > 6:
        return None
    smi = _lookup_exact(key)
    if smi is not None:
        return smi
    stripped = _HYDRATE_RE.sub("", key)
    if stripped != key:
        return _resolve_offline(stripped.strip(" ."), depth + 1)
    stripped = _LEADING_QUALIFIER_RE.sub("", key)
    if stripped == key:
        stripped = _TRAILING_QUALIFIER_RE.sub("", key)
    if stripped != key:
        return _resolve_offline(stripped.strip(), depth + 1)
    m = _STEREO_PREFIX_RE.match(key)
    if not m:
        return None
    base = _resolve_offline(key[m.end():].strip(), depth + 1)
    if base is None:
        return None
    if m.group("rac"):
        return _without_tetrahedral(base)
    if m.group("rs"):
        return _with_single_center(base, m.group("rslabel").upper())
    if m.group("multi") and "@" in base:
        # "(2S,3R)-threonine": accept when the descriptors are exactly the base's CIP labels
        wanted = sorted(re.findall(r"[rs]", m.group("multilabels").lower()))
        try:
            got = sorted(c.lower() for c in _cip_labels(_rd().Chem.MolFromSmiles(base)).values())
        except Exception:
            got = []
        return base if wanted == got else None
    # L-/D-/(+)/(-)/cis/trans/multi-descriptor prefixes on an unknown name: only safe
    # when the base structure carries no stereo information of its own.
    if "@" in base or "/" in base or "\\" in base:
        return None
    return base


def _validated(smiles: Optional[str]) -> Optional[str]:
    if not smiles:
        return None
    smiles = smiles.strip().split()[0] if smiles.strip() else ""
    return smiles if smiles and _smiles_ok(smiles) else None


def _opsin(name: str) -> Optional[str]:
    """Try py2opsin, an ``opsin`` executable or ``$OPSIN_JAR``; None on failure."""
    try:
        from py2opsin import py2opsin  # type: ignore

        try:
            out = py2opsin(name, output_format="SMILES")
            smi = _validated(out if isinstance(out, str) else None)
            if smi:
                return smi
        except Exception:
            pass
    except Exception:
        pass
    commands: List[List[str]] = []
    exe = shutil.which("opsin")
    if exe:
        commands.append([exe, "-osmi"])
    jar = os.environ.get("OPSIN_JAR")
    java = shutil.which("java")
    if jar and java and os.path.isfile(jar):
        commands.append([java, "-jar", jar, "-osmi"])
    for cmd in commands:
        try:
            proc = subprocess.run(
                cmd, input=name + "\n", capture_output=True, text=True, timeout=60, check=False,
            )
        except Exception:
            continue
        for line in proc.stdout.splitlines():
            smi = _validated(line)
            if smi:
                return smi
    return None


_PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/{prop}/TXT"


def _pubchem(name: str) -> Optional[str]:
    """Look a name up in PubChem PUG REST; never raises."""
    quoted = urllib.parse.quote(name.strip(), safe="")
    for prop in ("IsomericSMILES", "SMILES"):
        url = _PUBCHEM_URL.format(name=quoted, prop=prop)
        try:
            with urllib.request.urlopen(url, timeout=10) as resp:
                text = resp.read().decode("utf-8", "replace")
        except Exception:
            continue
        for line in text.splitlines():
            smi = _validated(line)
            if smi:
                return smi
    return None


def resolve_name(name: str, online: Optional[bool] = None) -> Optional[str]:
    """Resolve a chemical name to SMILES, or return None.

    Order: built-in dictionary of common names (exact, then without separators,
    then with hydrate suffixes / purity qualifiers / ``L-``, ``D-``, ``DL-``,
    ``(S)-``, ``(R)-`` prefixes handled), OPSIN (py2opsin, ``opsin`` on PATH or
    ``$OPSIN_JAR``), then PubChem PUG REST only when ``online`` is True or the
    environment variable ``CHEMDRAWMCP_ONLINE`` is ``"1"``.  Network or tool
    failures never raise.
    """
    if not isinstance(name, str):
        return None
    raw = name.strip()
    smi = _CASE_SENSITIVE.get(raw)
    if smi is not None:
        return smi
    key = _normalize_name(name)
    if not key:
        return None
    smi = _resolve_offline(key)
    if smi is not None:
        return smi
    try:
        smi = _opsin(raw)
    except Exception:
        smi = None
    if smi:
        return smi
    if online is None:
        online = os.environ.get("CHEMDRAWMCP_ONLINE", "").strip() == "1"
    if not online:
        return None
    try:
        return _pubchem(name)
    except Exception:
        return None
