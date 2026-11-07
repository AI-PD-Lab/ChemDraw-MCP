"""Backend interface: how the MCP server drives a ChemDraw application instance.

A backend wraps one platform's automation channel to the ChemDraw desktop app
(AppleScript/JXA on macOS, COM automation on Windows).  Everything that does not
need the live application (CDXML generation, RDKit chemistry, image analysis)
lives outside the backends so the server also works "headless" when ChemDraw
is not installed.

All methods raise :class:`ChemDrawError` on failure with a message that is
meaningful to an agent (including the text of any modal alert ChemDraw showed).
Document names are the names ChemDraw shows in the window title
(e.g. ``"Untitled ACS Document 1996-1"`` or ``"aspirin.mol"``); ``None`` means
"the front document".
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional


class ChemDrawError(RuntimeError):
    """Any failure talking to ChemDraw (including modal alerts, timeouts)."""


class ChemDrawNotAvailable(ChemDrawError):
    """ChemDraw is not installed / not reachable on this platform."""


# Export formats supported by every backend via ``save_document``.  Keys are
# what agents pass; values: (file extension, human name used by the app).
EXPORT_FORMATS: Dict[str, Dict[str, str]] = {
    "cdxml": {"ext": ".cdxml", "name": "ChemDraw XML", "mime": "chemical/x-cdxml"},
    "cdx": {"ext": ".cdx", "name": "ChemDraw", "mime": "chemical/x-cdx"},
    # ChemDraw 23's plain "MDL Molfile"/"SDfile"/"RXN" types write V3000; the V2000 variants are the compatible ones
    "mol": {"ext": ".mol", "name": "MDL Molfile V2000", "mime": "chemical/x-mdl-molfile"},
    "molv3000": {"ext": ".mol", "name": "MDL Molfile", "mime": "chemical/x-mdl-molfile"},
    "sdf": {"ext": ".sdf", "name": "MDL SDfile V2000", "mime": "chemical/x-mdl-sdfile"},
    "sdfv3000": {"ext": ".sdf", "name": "MDL SDfile", "mime": "chemical/x-mdl-sdfile"},
    "rxn": {"ext": ".rxn", "name": "Reaction Molfile V2000 (RXN)", "mime": "chemical/x-mdl-rxnfile"},
    "rxnv3000": {"ext": ".rxn", "name": "Reaction Molfile (RXN)", "mime": "chemical/x-mdl-rxnfile"},
    "rdf": {"ext": ".rdf", "name": "MDL RDfile", "mime": "chemical/x-mdl-rdfile"},
    "cml": {"ext": ".cml", "name": "Chemical Markup Language (CML)", "mime": "chemical/x-cml"},
    "ct": {"ext": ".ct", "name": "Connection Table", "mime": "chemical/x-ct"},
    "png": {"ext": ".png", "name": "PNG", "mime": "image/png"},
    "svg": {"ext": ".svg", "name": "Scalable Vector Graphics (SVG)", "mime": "image/svg+xml"},
    "pdf": {"ext": ".pdf", "name": "PDF", "mime": "application/pdf"},
    "eps": {"ext": ".eps", "name": "Encapsulated PostScript", "mime": "application/postscript"},
    "tiff": {"ext": ".tiff", "name": "TIFF", "mime": "image/tiff"},
    "gif": {"ext": ".gif", "name": "GIF", "mime": "image/gif"},
    "jpeg": {"ext": ".jpeg", "name": "JPEG", "mime": "image/jpeg"},
    "bmp": {"ext": ".bmp", "name": "BMP", "mime": "image/bmp"},
    "3mf": {"ext": ".3mf", "name": "3D Manufacturing Format for Display (3MF)",
            "mime": "application/vnd.ms-package.3dmanufacturing-3dmodel+xml"},
}
FORMAT_ALIASES = {"jpg": "jpeg", "tif": "tiff", "molfile": "mol", "molv2000": "mol", "sd": "sdf", "sdfv2000": "sdf",
                  "rxnv2000": "rxn", "mol3000": "molv3000", "v3000": "molv3000", "chemdraw": "cdx", "xml": "cdxml"}


def normalize_format(fmt: str) -> str:
    f = (fmt or "").strip().lower().lstrip(".")
    f = FORMAT_ALIASES.get(f, f)
    if f not in EXPORT_FORMATS:
        raise ChemDrawError("unknown export format %r; use one of %s" % (fmt, ", ".join(sorted(EXPORT_FORMATS))))
    return f


# Kinds of data that ``insert`` can paste into a document.
INSERT_KINDS = ("cdxml", "smiles", "inchi", "mol", "name", "helm", "sequence", "peptide", "dna", "rna", "text")


@dataclass
class DocumentInfo:
    name: str
    path: Optional[str] = None
    modified: bool = False
    front: bool = False
    atoms: int = 0
    bonds: int = 0
    molecules: int = 0
    captions: int = 0
    arrows: int = 0
    window_id: Optional[int] = None
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        if not d["extra"]:
            d.pop("extra")
        return d


@dataclass
class CommandInfo:
    name: str
    enabled: Optional[bool] = None
    short_description: str = ""
    long_description: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


class ChemDrawBackend(ABC):
    """Abstract automation channel to one ChemDraw application."""

    platform: str = "none"           # "macos" | "windows" | "none"

    # -- application lifecycle ------------------------------------------------
    @abstractmethod
    def executable(self) -> Optional[str]:
        """Path of the ChemDraw application/executable, or None if not found."""

    def available(self) -> bool:
        return bool(self.executable())

    @abstractmethod
    def is_running(self) -> bool:
        """True if a ChemDraw process is running that this backend can talk to."""

    @abstractmethod
    def launch(self, wait_seconds: float = 90.0) -> dict:
        """Start ChemDraw if it is not running and wait until it answers.  Returns {pid, version, ...}."""

    @abstractmethod
    def quit(self, save: bool = False) -> None:
        """Quit ChemDraw (discarding unsaved changes unless save=True)."""

    @abstractmethod
    def version(self) -> str:
        """ChemDraw version string (e.g. '23.1.2.7')."""

    @abstractmethod
    def activate(self) -> None:
        """Bring ChemDraw to the foreground."""

    # -- documents ------------------------------------------------------------
    @abstractmethod
    def documents(self) -> List[DocumentInfo]:
        """All open documents (hidden auto-save documents excluded), front one first."""

    @abstractmethod
    def front_document(self) -> Optional[str]:
        """Name of the front document or None if no document is open."""

    @abstractmethod
    def activate_document(self, name: Optional[str]) -> str:
        """Bring the document to the front; returns its (resolved) name."""

    @abstractmethod
    def new_document(self, stationery: Optional[str] = None) -> str:
        """Create a document (optionally from a stationery/style file or style name); returns its name."""

    @abstractmethod
    def open_document(self, path: str) -> str:
        """Open any file ChemDraw can read (cdxml, cdx, mol, sdf, rxn, cml, ...); returns the document name."""

    @abstractmethod
    def close_document(self, name: Optional[str] = None, save: bool = False) -> None:
        """Close a document, discarding changes unless save=True."""

    @abstractmethod
    def save_document(self, name: Optional[str], path: str, fmt: str) -> str:
        """Export/save a document to ``path`` in format ``fmt`` (see EXPORT_FORMATS).  Returns the actual path written."""

    @abstractmethod
    def document_cdxml(self, name: Optional[str] = None) -> str:
        """The whole document as CDXML text."""

    @abstractmethod
    def document_settings(self, name: Optional[str] = None) -> Dict[str, Any]:
        """Drawing settings of a document: bond length, line widths, fonts, sizes... (points where applicable)."""

    @abstractmethod
    def set_document_settings(self, name: Optional[str], settings: Dict[str, Any]) -> Dict[str, Any]:
        """Change document settings (same keys as ``document_settings``); returns the new settings."""

    # -- content --------------------------------------------------------------
    @abstractmethod
    def insert(self, name: Optional[str], data: str, kind: str = "cdxml") -> dict:
        """Insert chemistry into a document from text: kind in INSERT_KINDS.

        Returns {"atoms_before", "atoms_after", "molecules", ...}.  The inserted objects are left selected.
        """

    @abstractmethod
    def select_all(self, name: Optional[str] = None) -> None:
        """Select everything in the document."""

    @abstractmethod
    def deselect(self, name: Optional[str] = None) -> None:
        """Clear the selection."""

    @abstractmethod
    def selection_info(self, name: Optional[str] = None) -> Dict[str, Any]:
        """Chemistry of the current selection: smiles, formula, molecular_weight, exact_mass,
        elemental_analysis, text (captions), bounds, atoms, bonds, molecules."""

    @abstractmethod
    def selection_cdxml(self, name: Optional[str] = None) -> str:
        """The current selection as CDXML text (via Copy As > CDXML)."""

    @abstractmethod
    def selection_as(self, name: Optional[str], what: str) -> str:
        """Selection converted by ChemDraw: what in {smiles, inchi, inchikey, mol, molv3000, cdxml, helm, sln}."""

    # -- commands & menus -----------------------------------------------------
    @abstractmethod
    def list_commands(self, pattern: Optional[str] = None, enabled_only: bool = False,
                      with_descriptions: bool = True) -> List[CommandInfo]:
        """ChemDraw's named commands (e.g. cleanStructure, convertNameToStructure, copyAsSMILES)."""

    @abstractmethod
    def run_command(self, command: str, name: Optional[str] = None) -> dict:
        """Execute a named command on the (front) document.  Returns {ok, command, enabled_before, alerts}."""

    @abstractmethod
    def menus(self) -> Dict[str, List[str]]:
        """Top-level menus and their item titles."""

    @abstractmethod
    def do_menu(self, menu: str, item: str, name: Optional[str] = None) -> dict:
        """Invoke a menu item by menu title and item title."""

    # -- pictures -------------------------------------------------------------
    @abstractmethod
    def screenshot(self, name: Optional[str] = None) -> bytes:
        """PNG bytes of the document window as the user sees it."""

    # -- dialogs --------------------------------------------------------------
    @abstractmethod
    def dismiss_dialogs(self) -> List[str]:
        """Close any modal alert/sheet ChemDraw is showing; returns their texts."""

    # -- optional conveniences ------------------------------------------------
    def ping(self) -> bool:
        try:
            self.version()
            return True
        except Exception:
            return False

    def capabilities(self) -> Dict[str, Any]:
        return {"platform": self.platform, "executable": self.executable(), "available": self.available()}


class NullBackend(ChemDrawBackend):
    """Used where no ChemDraw automation exists (Linux, or ChemDraw missing): every call raises."""

    platform = "none"

    def __init__(self, reason: str = "ChemDraw automation is not available on this platform"):
        self.reason = reason

    def executable(self) -> Optional[str]:
        return None

    def _raise(self, *a, **k):
        raise ChemDrawNotAvailable(self.reason)

    is_running = lambda self: False  # noqa: E731
    launch = quit = version = activate = _raise
    documents = front_document = activate_document = new_document = open_document = _raise
    close_document = save_document = document_cdxml = document_settings = set_document_settings = _raise
    insert = select_all = deselect = selection_info = selection_cdxml = selection_as = _raise
    list_commands = run_command = menus = do_menu = screenshot = _raise

    def dismiss_dialogs(self) -> List[str]:
        return []
