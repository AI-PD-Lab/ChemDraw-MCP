"""Reader for ChemDraw's binary CDX format (``.cdx``).

CDX is a tagged little-endian binary stream (CambridgeSoft, "CDX file format").  This module parses it into a
tree of :class:`CDXObject` and converts that tree into CDXML text that :func:`chemdrawmcp.cdxml.read_cdxml`
understands, so ``.cdx`` files can be read without ChemDraw::

    from chemdrawmcp.cdx import read_cdx, cdx_to_cdxml, is_cdx
    info = read_cdx(open("scheme.cdx", "rb").read())      # same dict as cdxml.read_cdxml()
    xml = cdx_to_cdxml(data)                               # CDXML text

Everything below marked (v) was verified empirically on six documents saved by ChemDraw 23.1.2 (macOS) in both
formats (object ids are identical in the ``.cdx`` and the ``.cdxml`` of the same document, so every property
could be aligned with its CDXML attribute value by id); entries marked (s) come from the published spec /
``CDXConstants.h`` and were not exercised by the samples.

File layout (v)
---------------
* 8-byte magic ``VjCD0100``, 4 bytes ``04 03 02 01``, then 10 reserved zero bytes: the Document object tag
  (0x8000) starts at offset 22.  (The spec text describes 16 reserved bytes = offset 28; both are accepted.)
* A stream of tagged items.  A UINT16 tag with the high bit set opens an *object*: a UINT32 object id follows,
  then the object's own properties and child objects until a 0x0000 tag closes it.  Any other tag is a
  *property*: UINT16 length (0xFFFF means a UINT32 length follows), then the value bytes.
* CDXCoordinate = INT32 in 1/65536 point (y grows downwards).  2D positions are stored (y, x); 3D points
  (x, y, z); rectangles (top, left, bottom, right) -> CDXML ``BoundingBox="left top right bottom"``.
* CDXString = UINT16 style-run count, runs of {UINT16 start, UINT16 font, UINT16 face, UINT16 size (1/20 pt),
  UINT16 color}, then the text bytes in the font's charset (10000 = Mac Roman, 1252 = Windows Latin);
  line breaks are CR.  A few strings (objecttag Name) are written without the run-count header.
* CDXFontStyle (LabelStyle/CaptionStyle) = {UINT16 font, UINT16 face, UINT16 size (1/20 pt), UINT16 color}.
* FontTable = UINT16 platform, UINT16 count, then {UINT16 id, UINT16 charset, UINT16 length, name}.
  ColorTable = UINT16 count, then {UINT16 r, UINT16 g, UINT16 b} (0..65535).

Object tags
-----------
(v) 0x8000 Document, 0x8001 page, 0x8002 group, 0x8003 fragment, 0x8004 n (node/atom), 0x8005 b (bond),
    0x8006 t (text), 0x8007 graphic, 0x8008 curve, 0x8009 embeddedobject, 0x800D scheme, 0x800E step,
    0x8011 objecttag, 0x8017 bracketedgroup, 0x8018 bracketattachment, 0x8020 chemicalproperty,
    0x8021 arrow, 0x802B annotation.
(s) 0x800A altgroup, 0x800B templategrid, 0x800C regnum, 0x800F objectdefinition, 0x8010 spectrum,
    0x8012 oleclientitem, 0x8013 sequence, 0x8014 crossreference, 0x8015 splitter, 0x8016 table,
    0x8019 crossingbond, 0x801A border, 0x801B geometry, 0x801C constraint, 0x801D..0x801F tlc*,
    0x8022..0x8024 stoichiometry grid, 0x8025 bioshape, 0x8026..0x8028 plasmid*, 0x8029/0x802A rlogic*,
    0x802C..0x802F gel electrophoresis, 0x8030 documentproperties, 0x8031 property, 0x8032 coloredmoleculararea.

Nesting seen (v): Document > page > {group, fragment, t, graphic, arrow, curve, scheme, chemicalproperty,
annotation, bracketedgroup, embeddedobject}; group > {fragment, t, graphic, arrow}; fragment > {n, b, graphic};
n > {fragment (nickname/abbreviation contents), t (atom label, id 0), objecttag, annotation}; objecttag > t;
scheme > step; bracketedgroup > bracketattachment.

Property tags (see ``PROPERTY_TAGS`` for the complete table used by the converter)
----------------------------------------------------------------------------------
General (v): 0x0003 CreationProgram, 0x0008 Name, 0x000A Z, 0x000E RepresentsProperty (-> <represent>),
  0x000F IgnoreWarnings, 0x0010 Warning, 0x0013 SupersededBy (old-style graphic arrow -> new arrow object).
Position (v): 0x0200 p (y,x), 0x0204 BoundingBox, 0x0205 RotationAngle (raw INT32), 0x0207 Head3D,
  0x0208 Tail3D, 0x020D Center3D, 0x020E MajorAxisEnd3D, 0x020F MinorAxisEnd3D.
Colors/fonts (v): 0x0100 FontTable, 0x0300 ColorTable, 0x0301 color, 0x0302 bgcolor, 0x0303 FadePercent.
Atoms (v): 0x0400 NodeType (0 Unspecified, 4 Nickname, 5 Fragment, 7 GenericNickname,
  8 AnonymousAlternativeGroup, 12 ExternalConnectionPoint), 0x0402 Element, 0x0420 Isotope, 0x0421 Charge
  (INT8), 0x0422 Radical (1 Singlet), 0x042B NumHydrogens, 0x0430 Geometry (7 Tetrahedral),
  0x0431 BondOrdering, 0x0433 GenericNickname, 0x0437 AS (1 N, 2 R, 3 S, 6 u), 0x043B ShowAtomStereo,
  0x0440 ExternalConnectionType, 0x0448 NeedsClean, 0x0449 ResidueID, 0x044B ExternalConnectionNum,
  0x044D AtomID.  (s): 0x0401 LabelDisplay, 0x0403 ElementList, 0x0404 Formula, 0x0423..0x0429 query
  restrictions, 0x0432 Attachments, 0x0434 AltGroupID, 0x0439 AtomNumber, 0x0446/0x0447 enhanced stereo.
Bonds (v): 0x0600 Order (bitmask: 1,2,4 = single/double/triple, 128 = 1.5, 256 = 2.5), 0x0601 Display
  (2 Hash, 3 WedgedHashBegin, 4 WedgedHashEnd, 5 Bold, 6 WedgeBegin, 7 WedgeEnd, 8 Wavy), 0x0603
  DoublePosition, 0x0604 B, 0x0605 E, 0x0608 BeginAttach, 0x0609 EndAttach, 0x060A BS (1 N, 2 E, 3 Z),
  0x060B BondCircularOrdering, 0x060D ShowBondStereo, 0x060E CrossingBonds, 0x0610 Connectivity.
Text (v): 0x0700 Text (CDXString -> <s> runs), 0x0701 Justification, 0x0702 LineHeight (1 = auto, else
  1/20 pt), 0x0703 WordWrapWidth, 0x0704 LineStarts, 0x0705 LabelAlignment, 0x0706 LabelLineHeight,
  0x0707 CaptionLineHeight, 0x0708 InterpretChemically, 0x0823 LabelJustification (-1 Right, 0 Left,
  1 Center, 5 Auto), 0x080C CaptionJustification.  0x0709 duplicates 0x0700 and is ignored.
Document (v): 0x0800 MacPrintInfo, 0x0802 PrintMargins, 0x0803 ChainAngle, 0x0804 BondSpacing (1/10 %),
  0x0805 BondLength, 0x0806 BoldWidth, 0x0807 LineWidth, 0x0808 MarginWidth, 0x0809 HashSpacing,
  0x080A LabelStyle (-> LabelFont/LabelSize/LabelFace), 0x080B CaptionStyle, 0x080D FractionalWidths,
  0x080E Magnification, 0x043A/0x043C/0x0442..0x0445/0x044A/0x060C/0x060D/0x060F Show* flags,
  0x0828 AminoAcidTermini, 0x0829 ShowSequenceTermini, 0x082A ShowSequenceBonds, 0x082B ResidueWrapCount,
  0x082C ResidueBlockCount, 0x0832 ShowSequenceUnlinkedBranches, 0x0BB4..0x0BCA ChemProp* captions,
  0x0C09..0x0C0C RxnAutonumber*.  (s) 0x0901 WindowPosition, 0x0902 WindowSize (ChemDraw writes "0 0" in its
  own CDXML regardless of the stored values).  Page (v): 0x080F WidthPages, 0x0810 HeightPages, 0x0811 DrawingSpace,
  0x0812 Width, 0x0813 Height, 0x0814 PageOverlap, 0x0816 HeaderPosition, 0x0818 FooterPosition,
  0x0819 PrintTrimMarks.  Fragment (v): 0x0505 ConnectionOrder, 0x0506 SequenceType.
Graphics/arrows (v): 0x0A00 GraphicType (1 Line, 2 Arc, 3 Rectangle, 4 Oval, 5 Orbital, 6 Bracket,
  7 Symbol), 0x0A02 ArrowType (2 FullHead, 16 Hollow), 0x0A03 RectangleType, 0x0A04 OvalType,
  0x0A05 OrbitalType, 0x0A06 BracketType (3 Square), 0x0A07 SymbolType (0 LonePair, 5 CircleMinus),
  0x0A20 HeadSize, 0x0A21 AngularSize (1/10 degree), 0x0A22 LipSize, 0x0A23 CurvePoints, 0x0A27
  BracketedObjectIDs, 0x0A2B GraphicID, 0x0A2F ArrowheadType (1 Solid, 2 Hollow), 0x0A30
  ArrowheadCenterSize, 0x0A31 ArrowheadWidth, 0x0A32 ShadowSize, 0x0A33 ArrowShaftSpacing,
  0x0A35 ArrowheadHead (2 Full), 0x0A36 ArrowheadTail, 0x0A37 FillType (1 None), 0x0A3C CornerRadius.
  (s): 0x0A01 LineType, 0x0A24..0x0A2A bracket/polymer attributes, 0x0A34 ArrowEquilibriumRatio,
  0x0A3A Dipole, 0x0A3B NoGo, 0x0A3D FrameType.
Reaction steps (v): 0x0C01 ReactionStepReactants, 0x0C02 ReactionStepProducts, 0x0C04 ReactionStepArrows,
  0x0C05 ReactionStepObjectsAboveArrow, 0x0C06 ReactionStepObjectsBelowArrow (UINT32 id lists).
  (s): 0x0C00 ReactionStepAtomMap, 0x0C03 ReactionStepPlusses.
Other (v): 0x0B82 BasisObjects, 0x0BB0 ChemicalPropertyType, 0x0BB1 ChemicalPropertyDisplayID,
  0x0BB2 ChemicalPropertyIsActive, 0x0D00 TagType, 0x0D06 PositioningType, 0x0D08 PositioningOffset,
  0x1500 Keyword, 0x1501 Content (annotation).

Limitations: unknown property tags are skipped (they are listed by :func:`dump_cdx`), unknown object tags are
emitted as ``<obj_XXXX>`` elements with their children; embedded pictures are emitted only when their format
can be sniffed (PNG/JPEG/GIF/BMP/PDF); spectra, tables, TLC plates and bio shapes are passed through
structurally without interpreting their specialised properties.
"""
from __future__ import annotations

import re
import struct
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Optional, Tuple
from xml.sax.saxutils import escape, quoteattr

__all__ = ["CDXError", "CDXObject", "CDXProperty", "CDXStyleRun", "OBJECT_TAGS", "PROPERTY_TAGS", "is_cdx",
           "parse_cdx", "cdx_to_cdxml", "read_cdx", "dump_cdx", "decode_cdx_string"]

CDX_MAGIC = b"VjCD0100"
# Deeper object nesting than this is treated as a corrupt stream (real ChemDraw documents nest < 10 levels); it keeps
# the recursive consumers (emit, iter, dump) well inside the interpreter's recursion limit.
MAX_OBJECT_DEPTH = 256
# Characters that XML 1.0 forbids anywhere in a document (C0 controls except TAB/LF/CR, surrogates, U+FFFE/U+FFFF).
_XML_ILLEGAL = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f\ud800-\udfff\ufffe\uffff]")


class CDXError(ValueError):
    """The data is not a CDX stream (or is too damaged to read)."""


# ---------------------------------------------------------------------------
# tag tables
# ---------------------------------------------------------------------------

OBJECT_TAGS: Dict[int, str] = {
    0x8000: "CDXML", 0x8001: "page", 0x8002: "group", 0x8003: "fragment", 0x8004: "n", 0x8005: "b", 0x8006: "t",
    0x8007: "graphic", 0x8008: "curve", 0x8009: "embeddedobject", 0x800A: "altgroup", 0x800B: "templategrid",
    0x800C: "regnum", 0x800D: "scheme", 0x800E: "step", 0x800F: "objectdefinition", 0x8010: "spectrum",
    0x8011: "objecttag", 0x8012: "oleclientitem", 0x8013: "sequence", 0x8014: "crossreference", 0x8015: "splitter",
    0x8016: "table", 0x8017: "bracketedgroup", 0x8018: "bracketattachment", 0x8019: "crossingbond", 0x801A: "border",
    0x801B: "geometry", 0x801C: "constraint", 0x801D: "tlcplate", 0x801E: "tlclane", 0x801F: "tlcspot",
    0x8020: "chemicalproperty", 0x8021: "arrow", 0x8022: "stoichiometrygrid", 0x8023: "sgcomponent",
    0x8024: "sgdatum", 0x8025: "bioshape", 0x8026: "plasmidmap", 0x8027: "plasmidmarker", 0x8028: "plasmidregion",
    0x8029: "rlogic", 0x802A: "rlogicitem", 0x802B: "annotation", 0x802C: "gepplate", 0x802D: "gepband",
    0x802E: "geplane", 0x802F: "marker", 0x8030: "documentproperties", 0x8031: "property",
    0x8032: "coloredmoleculararea",
}

# value kinds used by PROPERTY_TAGS:
#   string    CDXString (styled runs dropped)        int/uint   1/2/4/8-byte integer
#   coord     CDXCoordinate -> points                 point2d    (y, x) -> "x y"
#   point3d   (x, y, z)                               rect       (t, l, b, r) -> "l t r b"
#   bool      empty or non-zero byte = yes            ids        UINT32 list -> "1 2 3"
#   hex       raw bytes as hex                        lineheight 0 variable, 1 auto, else 1/20 pt
#   tenths    INT16 / 10                              u16list    UINT16 count + UINT16 values
#   curvepts  UINT16 count + (y, x) coordinate pairs  order      bond-order bitmask
#   enum:X    integer -> ENUMS[X]                     flags:X    bitmask -> space-joined names
#   text/fonttable/colortable/fontstyle/represent/blob/skip  handled specially
PROPERTY_TAGS: Dict[int, Tuple[str, str]] = {
    # general
    0x0001: ("CreationUserName", "string"), 0x0002: ("CreationDate", "hex"), 0x0003: ("CreationProgram", "string"),
    0x0004: ("ModificationUserName", "string"), 0x0005: ("ModificationDate", "hex"),
    0x0006: ("ModificationProgram", "string"), 0x0008: ("Name", "string"), 0x0009: ("Comment", "string"),
    0x000A: ("Z", "uint"), 0x000B: ("RegistryNumber", "string"), 0x000C: ("RegistryAuthority", "string"),
    0x000E: ("RepresentsProperty", "represent"), 0x000F: ("IgnoreWarnings", "bool"), 0x0010: ("Warning", "string"),
    0x0011: ("Visible", "bool"), 0x0012: ("Transparent", "bool"), 0x0013: ("SupersededBy", "uint"),
    0x0014: ("StructurePerspective", "bool"),
    # fonts and colors
    0x0100: ("fonttable", "fonttable"), 0x0300: ("colortable", "colortable"), 0x0301: ("color", "uint"),
    0x0302: ("bgcolor", "uint"), 0x0303: ("FadePercent", "uint"), 0x0304: ("ForegroundAlpha", "uint"),
    0x0305: ("BackgroundAlpha", "uint"), 0x0306: ("HighlightColor", "uint"),
    # coordinates
    0x0200: ("p", "point2d"), 0x0201: ("xyz", "point3d"), 0x0202: ("extent", "point2d"),
    0x0203: ("extent3D", "point3d"), 0x0204: ("BoundingBox", "rect"), 0x0205: ("RotationAngle", "int"),
    0x0206: ("BoundsInParent", "rect"), 0x0207: ("Head3D", "point3d"), 0x0208: ("Tail3D", "point3d"),
    0x0209: ("TopLeft", "point2d"), 0x020A: ("TopRight", "point2d"), 0x020B: ("BottomRight", "point2d"),
    0x020C: ("BottomLeft", "point2d"), 0x020D: ("Center3D", "point3d"), 0x020E: ("MajorAxisEnd3D", "point3d"),
    0x020F: ("MinorAxisEnd3D", "point3d"),
    # atoms
    0x0400: ("NodeType", "enum:NodeType"), 0x0401: ("LabelDisplay", "enum:LabelAlignment"), 0x0402: ("Element", "uint"),
    0x0403: ("ElementList", "elementlist"), 0x0404: ("Formula", "string"), 0x0420: ("Isotope", "uint"),
    0x0421: ("Charge", "int"), 0x0422: ("Radical", "enum:Radical"), 0x0423: ("FreeSites", "uint"),
    0x0424: ("ImplicitHydrogens", "bool"), 0x0425: ("RingBondCount", "enum:RingBondCount"),
    0x0426: ("UnsaturatedBonds", "enum:Unsaturation"), 0x0427: ("RxnChange", "bool"),
    0x0428: ("RxnStereo", "enum:RxnStereo"), 0x0429: ("AbnormalValence", "bool"), 0x042B: ("NumHydrogens", "uint"),
    0x042E: ("HDot", "bool"), 0x042F: ("HDash", "bool"), 0x0430: ("Geometry", "enum:Geometry"),
    0x0431: ("BondOrdering", "ids"), 0x0432: ("Attachments", "ids"), 0x0433: ("GenericNickname", "string"),
    0x0434: ("AltGroupID", "uint"), 0x0435: ("SubstituentsUpTo", "uint"), 0x0436: ("SubstituentsExactly", "uint"),
    0x0437: ("AS", "enum:AtomCIP"), 0x0438: ("Translation", "enum:Translation"), 0x0439: ("AtomNumber", "string"),
    0x043A: ("ShowAtomQuery", "bool"), 0x043B: ("ShowAtomStereo", "bool"), 0x043C: ("ShowAtomNumber", "bool"),
    0x043D: ("LinkCountLow", "uint"), 0x043E: ("LinkCountHigh", "uint"),
    0x043F: ("IsotopicAbundance", "enum:IsotopicAbundance"),
    0x0440: ("ExternalConnectionType", "enum:ExternalConnectionType"), 0x0441: ("GenericList", "string"),
    0x0442: ("ShowTerminalCarbonLabels", "bool"), 0x0443: ("ShowNonTerminalCarbonLabels", "bool"),
    0x0444: ("HideImplicitHydrogens", "bool"), 0x0445: ("ShowAtomEnhancedStereo", "bool"),
    0x0446: ("EnhancedStereoType", "enum:EnhancedStereoType"), 0x0447: ("EnhancedStereoGroupNum", "uint"),
    0x0448: ("NeedsClean", "bool"), 0x0449: ("ResidueID", "string"), 0x044A: ("ShowResidueID", "bool"),
    0x044B: ("ExternalConnectionNum", "int"), 0x044C: ("ShowAtomID", "bool"), 0x044D: ("AtomID", "uint"),
    # fragments
    0x0505: ("ConnectionOrder", "ids"), 0x0506: ("SequenceType", "enum:SequenceType"),
    # bonds
    0x0600: ("Order", "order"), 0x0601: ("Display", "enum:BondDisplay"), 0x0602: ("Display2", "enum:BondDisplay"),
    0x0603: ("DoublePosition", "enum:DoublePosition"), 0x0604: ("B", "uint"), 0x0605: ("E", "uint"),
    0x0606: ("Topology", "enum:Topology"), 0x0607: ("RxnParticipation", "enum:RxnParticipation"),
    0x0608: ("BeginAttach", "uint"), 0x0609: ("EndAttach", "uint"), 0x060A: ("BS", "enum:BondCIP"),
    0x060B: ("BondCircularOrdering", "ids"), 0x060C: ("ShowBondQuery", "bool"), 0x060D: ("ShowBondStereo", "bool"),
    0x060E: ("CrossingBonds", "ids"), 0x060F: ("ShowBondRxn", "bool"), 0x0610: ("Connectivity", "enum:Connectivity"),
    0x0611: ("BeginExternalNum", "uint"), 0x0612: ("EndExternalNum", "uint"),
    # text
    0x0700: ("Text", "text"), 0x0701: ("Justification", "enum:Justification"), 0x0702: ("LineHeight", "lineheight"),
    0x0703: ("WordWrapWidth", "int"), 0x0704: ("LineStarts", "u16list"), 0x0705: ("LabelAlignment", "enum:LabelAlignment"),
    0x0706: ("LabelLineHeight", "lineheight"), 0x0707: ("CaptionLineHeight", "lineheight"),
    0x0708: ("InterpretChemically", "bool"), 0x0709: ("TextInterpreted", "skip"),
    # document and page
    0x0800: ("MacPrintInfo", "hex"), 0x0801: ("WinPrintInfo", "hex"), 0x0802: ("PrintMargins", "rect"),
    0x0803: ("ChainAngle", "coord"), 0x0804: ("BondSpacing", "tenths"), 0x0805: ("BondLength", "coord"),
    0x0806: ("BoldWidth", "coord"), 0x0807: ("LineWidth", "coord"), 0x0808: ("MarginWidth", "coord"),
    0x0809: ("HashSpacing", "coord"), 0x080A: ("Label", "fontstyle"), 0x080B: ("Caption", "fontstyle"),
    0x080C: ("CaptionJustification", "enum:Justification"), 0x080D: ("FractionalWidths", "bool"),
    0x080E: ("Magnification", "uint"), 0x080F: ("WidthPages", "uint"), 0x0810: ("HeightPages", "uint"),
    0x0811: ("DrawingSpace", "enum:DrawingSpace"), 0x0812: ("Width", "coord"), 0x0813: ("Height", "coord"),
    0x0814: ("PageOverlap", "coord"), 0x0815: ("Header", "string"), 0x0816: ("HeaderPosition", "coord"),
    0x0817: ("Footer", "string"), 0x0818: ("FooterPosition", "coord"), 0x0819: ("PrintTrimMarks", "bool"),
    0x0822: ("BondSpacingAbs", "coord"), 0x0823: ("LabelJustification", "enum:Justification"),
    0x0828: ("AminoAcidTermini", "enum:AminoAcidTermini"), 0x0829: ("ShowSequenceTermini", "bool"),
    0x082A: ("ShowSequenceBonds", "bool"), 0x082B: ("ResidueWrapCount", "uint"), 0x082C: ("ResidueBlockCount", "uint"),
    0x0832: ("ShowSequenceUnlinkedBranches", "bool"), 0x0901: ("WindowPosition", "point2d"),
    0x0902: ("WindowSize", "point2d"),
    # graphics, arrows, brackets, curves
    0x0A00: ("GraphicType", "enum:GraphicType"), 0x0A01: ("LineType", "flags:LineType"),
    0x0A02: ("ArrowType", "enum:ArrowType"), 0x0A03: ("RectangleType", "flags:RectangleType"),
    0x0A04: ("OvalType", "flags:OvalType"), 0x0A05: ("OrbitalType", "enum:OrbitalType"),
    0x0A06: ("BracketType", "enum:BracketType"), 0x0A07: ("SymbolType", "enum:SymbolType"),
    0x0A20: ("HeadSize", "uint"), 0x0A21: ("AngularSize", "tenths"), 0x0A22: ("LipSize", "uint"),
    0x0A23: ("CurvePoints", "curvepts"), 0x0A24: ("BracketUsage", "enum:BracketUsage"),
    0x0A25: ("PolymerRepeatPattern", "enum:PolymerRepeat"), 0x0A26: ("PolymerFlipType", "enum:PolymerFlip"),
    0x0A27: ("BracketedObjectIDs", "ids"), 0x0A28: ("RepeatCount", "int"), 0x0A29: ("ComponentOrder", "int"),
    0x0A2A: ("SRULabel", "string"), 0x0A2B: ("GraphicID", "uint"), 0x0A2C: ("BondID", "uint"),
    0x0A2D: ("InnerAtomID", "uint"), 0x0A2F: ("ArrowheadType", "enum:ArrowheadType"),
    0x0A30: ("ArrowheadCenterSize", "uint"), 0x0A31: ("ArrowheadWidth", "uint"), 0x0A32: ("ShadowSize", "uint"),
    0x0A33: ("ArrowShaftSpacing", "uint"), 0x0A34: ("ArrowEquilibriumRatio", "uint"),
    0x0A35: ("ArrowheadHead", "enum:Arrowhead"), 0x0A36: ("ArrowheadTail", "enum:Arrowhead"),
    0x0A37: ("FillType", "enum:FillType"), 0x0A38: ("CurveSpacing", "uint"), 0x0A39: ("Closed", "bool"),
    0x0A3A: ("Dipole", "bool"), 0x0A3B: ("NoGo", "enum:NoGo"), 0x0A3C: ("CornerRadius", "uint"),
    0x0A3D: ("FrameType", "enum:FrameType"), 0x0A70: ("Picture", "blob"), 0x0A71: ("Picture", "blob"),
    # chemical properties and annotations
    0x0B82: ("BasisObjects", "ids"), 0x0BB0: ("ChemicalPropertyType", "uint"),
    0x0BB1: ("ChemicalPropertyDisplayID", "uint"), 0x0BB2: ("ChemicalPropertyIsActive", "bool"),
    0x0BB4: ("ChemPropName", "string"), 0x0BB5: ("ChemPropFormula", "string"), 0x0BB6: ("ChemPropExactMass", "string"),
    0x0BB7: ("ChemPropMolWt", "string"), 0x0BB8: ("ChemPropMOverZ", "string"), 0x0BB9: ("ChemPropAnalysis", "string"),
    0x0BBA: ("ChemPropBoilingPt", "string"), 0x0BBB: ("ChemPropMeltingPt", "string"),
    0x0BBC: ("ChemPropCritTemp", "string"), 0x0BBD: ("ChemPropCritPres", "string"),
    0x0BBE: ("ChemPropCritVol", "string"), 0x0BBF: ("ChemPropGibbs", "string"), 0x0BC0: ("ChemPropLogP", "string"),
    0x0BC1: ("ChemPropMR", "string"), 0x0BC2: ("ChemPropHenry", "string"), 0x0BC3: ("ChemPropEForm", "string"),
    0x0BC4: ("ChemProptPSA", "string"), 0x0BC9: ("ChemPropID", "string"), 0x0BCA: ("ChemPropFragmentLabel", "string"),
    # reaction steps and autonumbering
    0x0C00: ("ReactionStepAtomMap", "ids"), 0x0C01: ("ReactionStepReactants", "ids"),
    0x0C02: ("ReactionStepProducts", "ids"), 0x0C03: ("ReactionStepPlusses", "ids"),
    0x0C04: ("ReactionStepArrows", "ids"), 0x0C05: ("ReactionStepObjectsAboveArrow", "ids"),
    0x0C06: ("ReactionStepObjectsBelowArrow", "ids"), 0x0C07: ("ReactionStepAtomMapManual", "ids"),
    0x0C08: ("ReactionStepAtomMapAuto", "ids"), 0x0C09: ("RxnAutonumberStyle", "enum:RxnAutonumberStyle"),
    0x0C0A: ("RxnAutonumberConditions", "bool"), 0x0C0B: ("RxnAutonumberStart", "uint"),
    0x0C0C: ("RxnAutonumberFormat", "string"),
    # object tags and positioning
    0x0D00: ("TagType", "enum:TagType"), 0x0D02: ("Tracking", "bool"), 0x0D03: ("Persistent", "bool"),
    0x0D06: ("PositioningType", "enum:PositioningType"), 0x0D07: ("PositioningAngle", "int"),
    0x0D08: ("PositioningOffset", "point2d"),
    # annotations
    0x1500: ("Keyword", "string"), 0x1501: ("Content", "string"),
}

ENUMS: Dict[str, Dict[int, Optional[str]]] = {
    "NodeType": {0: "Unspecified", 1: "Element", 2: "ElementList", 3: "ElementListNickname", 4: "Nickname",
                 5: "Fragment", 6: "Formula", 7: "GenericNickname", 8: "AnonymousAlternativeGroup",
                 9: "NamedAlternativeGroup", 10: "MultiAttachment", 11: "VariableAttachment",
                 12: "ExternalConnectionPoint", 13: "LinkNode", 14: "Monomer"},
    "LabelAlignment": {0: "Auto", 1: "Left", 2: "Center", 3: "Right", 4: "Above", 5: "Below", 6: "Best"},
    "Radical": {0: "None", 1: "Singlet", 2: "Doublet", 3: "Triplet"},
    "RingBondCount": {0: "Unspecified", 1: "NoRingBonds", 2: "AsDrawn", 3: "SimpleRing", 4: "Fusion", 5: "SpiroOrHigher"},
    "Unsaturation": {0: "Unspecified", 1: "MustBeAbsent", 2: "MustBePresent"},
    "RxnStereo": {0: "Unspecified", 1: "Inversion", 2: "Retention"},
    "Geometry": {0: "Unknown", 1: "1", 2: "Linear", 3: "Bent", 4: "TrigonalPlanar", 5: "TrigonalPyramidal",
                 6: "SquarePlanar", 7: "Tetrahedral", 8: "TrigonalBipyramidal", 9: "SquarePyramidal", 10: "5",
                 11: "Octahedral", 12: "6", 13: "7", 14: "8", 15: "9", 16: "10"},
    "AtomCIP": {0: "U", 1: "N", 2: "R", 3: "S", 4: "r", 5: "s", 6: "u"},
    "BondCIP": {0: "U", 1: "N", 2: "E", 3: "Z"},
    "Translation": {0: "Equal", 1: "Broad", 2: "Narrow", 3: "Any"},
    "IsotopicAbundance": {0: "Unspecified", 1: "Any", 2: "Natural", 3: "Enriched", 4: "Deficient", 5: "Nonnatural"},
    "ExternalConnectionType": {0: "Unspecified", 1: "Diamond", 2: "Star", 3: "PolymerBead", 4: "Wavy", 5: "Residue",
                               6: "Peptide", 7: "DNA", 8: "RNA", 9: "Terminal", 10: "Sulfur", 11: "Nucleophile"},
    "EnhancedStereoType": {0: "Unspecified", 1: "None", 2: "Absolute", 3: "Or", 4: "And"},
    "SequenceType": {0: "Unknown", 1: "Peptide", 2: "Peptide1", 3: "Peptide3", 4: "DNA", 5: "RNA", 6: "Biopolymer"},
    "BondDisplay": {0: "Solid", 1: "Dash", 2: "Hash", 3: "WedgedHashBegin", 4: "WedgedHashEnd", 5: "Bold",
                    6: "WedgeBegin", 7: "WedgeEnd", 8: "Wavy", 9: "HollowWedgeBegin", 10: "HollowWedgeEnd",
                    11: "WavyWedgeBegin", 12: "WavyWedgeEnd", 13: "Dot", 14: "DashDot", 15: "DativeBegin",
                    16: "DativeEnd"},
    # automatic double-bond placement (0..2) is not written to CDXML by ChemDraw; user choices are
    "DoublePosition": {0: None, 1: None, 2: None, 256: "Center", 257: "Right", 258: "Left"},
    "Topology": {0: "Unspecified", 1: "Ring", 2: "Chain", 3: "RingOrChain"},
    "RxnParticipation": {0: "Unspecified", 1: "ReactionCenter", 2: "MakeOrBreak", 3: "ChangeType", 4: "MakeAndChange",
                         5: "NotReactionCenter", 6: "NoChange", 7: "Unmapped"},
    "Connectivity": {0: "Unspecified", 1: "Linear", 2: "Bridged", 3: "Staggered", 4: "Cyclic"},
    "Justification": {-1: "Right", 0: "Left", 1: "Center", 2: "Full", 3: "Above", 4: "Below", 5: "Auto",
                      6: "BestInitial"},
    "DrawingSpace": {0: "pages", 1: "poster"},
    "AminoAcidTermini": {1: "HOH", 2: "NH2/COOH"},
    "GraphicType": {0: "Undefined", 1: "Line", 2: "Arc", 3: "Rectangle", 4: "Oval", 5: "Orbital", 6: "Bracket",
                    7: "Symbol", 8: "Spectrum"},
    "ArrowType": {0: "NoHead", 1: "HalfHead", 2: "FullHead", 4: "Resonance", 8: "Equilibrium", 16: "Hollow",
                  32: "RetroSynthetic", 64: "NoGo", 128: "Dipole"},
    "OrbitalType": {0: "s", 1: "oval", 2: "lobe", 3: "p", 4: "hybridPlus", 5: "hybridMinus", 6: "dz2Plus",
                    7: "dz2Minus", 8: "dxy", 256: "sShaded", 257: "ovalShaded", 258: "lobeShaded", 259: "pShaded",
                    512: "sFilled", 513: "ovalFilled", 514: "lobeFilled", 515: "pFilled", 516: "hybridPlusFilled",
                    517: "hybridMinusFilled", 518: "dz2PlusFilled", 519: "dz2MinusFilled", 520: "dxyFilled"},
    "BracketType": {0: "RoundPair", 1: "SquarePair", 2: "CurlyPair", 3: "Square", 4: "Curly", 5: "Round"},
    "SymbolType": {0: "LonePair", 1: "Electron", 2: "RadicalCation", 3: "RadicalAnion", 4: "CirclePlus",
                   5: "CircleMinus", 6: "Dagger", 7: "DoubleDagger", 8: "Plus", 9: "Minus", 10: "Racemic",
                   11: "Absolute", 12: "Relative", 13: "LonePairBar"},
    "BracketUsage": {0: "Unspecified", 1: "Unused1", 2: "Unused2", 3: "SRU", 4: "Monomer", 5: "Mer", 6: "Copolymer",
                     7: "CopolymerAlternating", 8: "CopolymerRandom", 9: "CopolymerBlock", 10: "Crosslink",
                     11: "Graft", 12: "Modification", 13: "Component", 14: "MixtureUnordered", 15: "MixtureOrdered",
                     16: "MultipleGroup", 17: "Generic", 18: "Anypolymer"},
    "PolymerRepeat": {0: "HeadToTail", 1: "HeadToHead", 2: "EitherUnknown"},
    "PolymerFlip": {0: "Unspecified", 1: "NoFlip", 2: "Flip"},
    "ArrowheadType": {0: "Unspecified", 1: "Solid", 2: "Hollow", 3: "Angle"},
    "Arrowhead": {0: "Unspecified", 1: "None", 2: "Full", 3: "HalfLeft", 4: "HalfRight"},
    "FillType": {0: "Unspecified", 1: "None", 2: "Solid", 4: "Shaded", 8: "Gradient", 16: "Pattern"},
    "NoGo": {0: "Unspecified", 1: "None", 2: "Cross", 3: "Hash"},
    "FrameType": {0: "Unspecified", 1: "None", 2: "Plain", 4: "RoundEdge", 8: "Shadow", 16: "Shaded", 32: "Filled",
                  64: "Dashed", 128: "Bold"},
    "RxnAutonumberStyle": {0: "Roman", 1: "Arabic", 2: "Alpha"},
    "TagType": {0: "Unknown", 1: "Double", 2: "Long", 3: "String"},
    "PositioningType": {0: "auto", 1: "angle", 2: "offset", 3: "absolute"},
}

# bit flags: the 0 entry names the value with no bits set
FLAGS: Dict[str, Dict[int, str]] = {
    "LineType": {0: "Solid", 1: "Dashed", 2: "Bold", 4: "Wavy"},
    "RectangleType": {0: "Plain", 1: "RoundEdge", 2: "Shadow", 4: "Shaded", 8: "Filled", 16: "Dashed", 32: "Bold"},
    "OvalType": {0: "Plain", 1: "Circle", 2: "Shaded", 4: "Filled", 8: "Dashed", 16: "Bold"},
}

BOND_ORDER_BITS: Dict[int, str] = {1: "1", 2: "2", 4: "3", 8: "4", 16: "5", 32: "6", 64: "0.5", 128: "1.5", 256: "2.5",
                                   512: "3.5", 1024: "4.5", 2048: "5.5", 4096: "dative", 8192: "ionic",
                                   16384: "hydrogen", 32768: "threecenter"}

# CDX charset number -> (CDXML charset name, Python codec)
CHARSETS: Dict[int, Tuple[str, str]] = {
    0: ("Unknown", "latin-1"), 1252: ("iso-8859-1", "cp1252"), 1250: ("iso-8859-2", "cp1250"),
    1251: ("iso-8859-5", "cp1251"), 1253: ("iso-8859-7", "cp1253"), 1254: ("iso-8859-9", "cp1254"),
    1255: ("iso-8859-8", "cp1255"), 1256: ("iso-8859-6", "cp1256"), 1257: ("iso-8859-4", "cp1257"),
    932: ("Shift_JIS", "cp932"), 936: ("GB2312", "gbk"), 949: ("EUC-KR", "cp949"), 950: ("big5", "big5"),
    10000: ("x-mac-roman", "mac_roman"), 10001: ("x-mac-japanese", "shift_jis"), 10006: ("x-mac-greek", "mac_greek"),
    10007: ("x-mac-cyrillic", "mac_cyrillic"), 10029: ("x-mac-ce", "mac_latin2"), 10079: ("x-mac-icelandic", "mac_iceland"),
    10081: ("x-mac-turkish", "mac_turkish"), 65001: ("utf-8", "utf-8"),
}

PICTURE_MAGIC: Tuple[Tuple[bytes, str], ...] = ((b"\x89PNG", "PNG"), (b"\xff\xd8\xff", "JPEG"), (b"GIF8", "GIF"),
                                                (b"BM", "BMP"), (b"%PDF", "PDF"), (b"\xd7\xcd\xc6\x9a", "WindowsMetafile"),
                                                (b"\x01\x00\x00\x00", "EnhancedMetafile"))


# ---------------------------------------------------------------------------
# object tree
# ---------------------------------------------------------------------------

@dataclass
class CDXProperty:
    tag: int
    data: bytes

    @property
    def name(self) -> str:
        spec = PROPERTY_TAGS.get(self.tag)
        return spec[0] if spec else "0x%04X" % self.tag


@dataclass
class CDXObject:
    tag: int
    id: int
    properties: List[CDXProperty] = field(default_factory=list)
    children: List["CDXObject"] = field(default_factory=list)

    @property
    def name(self) -> str:
        return OBJECT_TAGS.get(self.tag, "obj_%04x" % self.tag)

    def prop(self, tag: int) -> Optional[bytes]:
        """Value bytes of the first property with this tag (None if absent)."""
        for p in self.properties:
            if p.tag == tag:
                return p.data
        return None

    def iter(self) -> Iterator["CDXObject"]:
        """Depth-first iteration over this object and all descendants."""
        yield self
        for c in self.children:
            yield from c.iter()

    def find_all(self, tag: int) -> List["CDXObject"]:
        return [o for o in self.iter() if o.tag == tag]


@dataclass
class CDXStyleRun:
    start: int
    font: int
    face: int
    size: int       # in 1/20 point
    color: int
    text: str = ""


# ---------------------------------------------------------------------------
# low-level decoding
# ---------------------------------------------------------------------------

def _u16(b: bytes, pos: int = 0) -> int:
    return struct.unpack_from("<H", b, pos)[0]


def _u32(b: bytes, pos: int = 0) -> int:
    return struct.unpack_from("<I", b, pos)[0]


def _int(b: bytes, signed: bool = True) -> Optional[int]:
    fmt = {1: "b", 2: "h", 4: "i", 8: "q"}.get(len(b))
    if fmt is None:
        return None
    return struct.unpack("<" + (fmt if signed else fmt.upper()), b)[0]


def _coord(b: bytes, pos: int = 0) -> float:
    return struct.unpack_from("<i", b, pos)[0] / 65536.0


def _coords(b: bytes) -> List[float]:
    n = len(b) // 4
    return [v / 65536.0 for v in struct.unpack("<%di" % n, b[:n * 4])]


def _f(v: float) -> str:
    """Compact number formatting (2 decimals, like ChemDraw's CDXML writer)."""
    if abs(v - round(v)) < 0.005:
        return str(int(round(v)))
    return ("%.2f" % v).rstrip("0").rstrip(".")


def _parse_fonttable(b: bytes) -> Dict[int, Tuple[str, int]]:
    """{font id: (name, charset)} from the 0x0100 FontTable value."""
    fonts: Dict[int, Tuple[str, int]] = {}
    if len(b) < 4:
        return fonts
    count = _u16(b, 2)
    pos = 4
    for _ in range(count):
        if pos + 6 > len(b):
            break
        fid, charset, length = struct.unpack_from("<HHH", b, pos)
        pos += 6
        name = b[pos:pos + length].decode(CHARSETS.get(charset, ("", "latin-1"))[1], "replace")
        pos += length
        fonts[fid] = (name, charset)
    return fonts


def decode_cdx_string(data: bytes, fonts: Optional[Dict[int, Tuple[str, int]]] = None,
                      default_codec: str = "cp1252") -> Tuple[str, List[CDXStyleRun]]:
    """Decode a CDXString value into (text, style runs); CR line breaks become LF.

    ``fonts`` ({id: (name, charset)}) selects the codec per run; strings written without the run-count header
    (ChemDraw does that for objecttag names) are decoded as plain text.
    """
    fonts = fonts or {}
    if len(data) < 2:
        return data.decode(default_codec, "replace"), []
    count = _u16(data, 0)
    if 2 + 10 * count > len(data):
        return data.decode(default_codec, "replace").replace("\r\n", "\n").replace("\r", "\n"), []
    runs = [CDXStyleRun(*struct.unpack_from("<HHHHH", data, 2 + 10 * i)) for i in range(count)]
    raw = data[2 + 10 * count:]
    if not runs:
        return raw.decode(default_codec, "replace").replace("\r\n", "\n").replace("\r", "\n"), []
    runs.sort(key=lambda r: r.start)
    pieces = []
    for i, run in enumerate(runs):
        end = runs[i + 1].start if i + 1 < len(runs) else len(raw)
        chunk = raw[min(run.start, len(raw)):max(min(end, len(raw)), min(run.start, len(raw)))]
        charset = fonts.get(run.font, ("", -1))[1]
        codec = CHARSETS.get(charset, ("", default_codec))[1]
        run.text = chunk.decode(codec, "replace").replace("\r\n", "\n").replace("\r", "\n")
        pieces.append(run.text)
    return "".join(pieces), runs


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def is_cdx(data: bytes) -> bool:
    """True if ``data`` starts with the CDX magic ``VjCD0100``."""
    try:
        return bytes(data[:8]) == CDX_MAGIC
    except TypeError:
        return False


def _document_offset(data: bytes) -> int:
    for off in (22, 28):
        if data[off:off + 2] == b"\x00\x80":
            return off
    off = data.find(b"\x00\x80\x00\x00\x00\x00", 12, 96)   # Document tag followed by object id 0
    return off


def parse_cdx(data: bytes) -> CDXObject:
    """Parse CDX bytes into an object tree; the result is the Document object (tag 0x8000).

    Tolerant of truncation and unknown tags: unknown properties are kept as raw bytes, unknown objects are
    parsed recursively like any other object.
    """
    data = bytes(data)
    if not is_cdx(data):
        raise CDXError("not a CDX stream (expected it to start with %r)" % CDX_MAGIC)
    start = _document_offset(data)
    wrapper = CDXObject(-1, 0)
    stack = [wrapper]
    pos = start if start >= 0 else 28
    n = len(data)
    while pos + 2 <= n:
        tag = _u16(data, pos)
        pos += 2
        if tag == 0:
            if len(stack) > 1:
                stack.pop()
            if len(stack) == 1 and wrapper.children and wrapper.children[-1].tag == 0x8000:
                break   # end of the document object; ignore trailing bytes
            continue
        if tag & 0x8000:
            if pos + 4 > n:
                break
            if len(stack) > MAX_OBJECT_DEPTH:
                raise CDXError("object nesting deeper than %d levels (corrupt CDX stream?)" % MAX_OBJECT_DEPTH)
            obj = CDXObject(tag, _u32(data, pos))
            pos += 4
            stack[-1].children.append(obj)
            stack.append(obj)
            continue
        if pos + 2 > n:
            break
        length = _u16(data, pos)
        pos += 2
        if length == 0xFFFF:
            if pos + 4 > n:
                break
            length = _u32(data, pos)
            pos += 4
        value = data[pos:pos + length]
        pos += length
        stack[-1].properties.append(CDXProperty(tag, value))
        if len(value) < length:
            break   # truncated file
    if len(wrapper.children) == 1 and wrapper.children[0].tag == 0x8000 and not wrapper.properties:
        return wrapper.children[0]
    return CDXObject(0x8000, 0, wrapper.properties, wrapper.children)


# ---------------------------------------------------------------------------
# CDXML conversion
# ---------------------------------------------------------------------------

def _xml_text(s: str) -> str:
    """Drop characters that may not appear in XML 1.0 at all (not even escaped)."""
    return _XML_ILLEGAL.sub("", s)


def _attr_str(attrs: List[Tuple[str, str]]) -> str:
    return "".join(" %s=%s" % (k, quoteattr(_xml_text(v))) for k, v in attrs)


class _Converter:
    def __init__(self, doc: CDXObject):
        self.doc = doc
        raw = doc.prop(0x0100)
        self.fonts = _parse_fonttable(raw) if raw is not None else {}
        self.default_codec = "mac_roman" if any(cs == 10000 for _, cs in self.fonts.values()) else "cp1252"
        self.unknown: Dict[int, int] = {}

    # -- values -------------------------------------------------------------
    def string(self, data: bytes) -> str:
        return decode_cdx_string(data, self.fonts, self.default_codec)[0]

    def value(self, kind: str, data: bytes) -> Optional[str]:
        if kind == "string":
            return self.string(data)
        if kind in ("int", "uint"):
            v = _int(data, signed=(kind == "int"))
            return None if v is None else str(v)
        if kind == "coord":
            return _f(_coord(data)) if len(data) >= 4 else None
        if kind == "point2d":
            c = _coords(data)
            return "%s %s" % (_f(c[1]), _f(c[0])) if len(c) >= 2 else None
        if kind == "point3d":
            c = _coords(data)
            return " ".join(_f(v) for v in c[:3]) if len(c) >= 3 else None
        if kind == "rect":
            c = _coords(data)
            return "%s %s %s %s" % (_f(c[1]), _f(c[0]), _f(c[3]), _f(c[2])) if len(c) >= 4 else None
        if kind == "bool":
            return "yes" if len(data) == 0 or data[0] != 0 else "no"
        if kind == "ids":
            n = len(data) // 4
            return " ".join(str(v) for v in struct.unpack("<%dI" % n, data[:n * 4]))
        if kind == "hex":
            return data.hex().upper()
        if kind == "lineheight":
            v = _int(data)
            if v is None:
                return None
            if v == 0:          # "variable": ChemDraw omits the attribute in CDXML
                return None
            return "auto" if v == 1 else _f(v / 20.0)
        if kind == "tenths":
            v = _int(data)
            return None if v is None else _f(v / 10.0)
        if kind == "u16list":
            if len(data) < 2:
                return None
            n = min(_u16(data, 0), (len(data) - 2) // 2)
            return " ".join(str(v) for v in struct.unpack_from("<%dH" % n, data, 2))
        if kind == "curvepts":
            if len(data) < 2:
                return None
            n = min(_u16(data, 0), (len(data) - 2) // 8)
            c = _coords(data[2:2 + 8 * n])
            return " ".join("%s %s" % (_f(c[2 * i + 1]), _f(c[2 * i])) for i in range(n))
        if kind == "elementlist":
            if len(data) < 2:
                return None
            n = struct.unpack_from("<h", data, 0)[0]
            m = min(abs(n), (len(data) - 2) // 2)
            elems = " ".join(str(v) for v in struct.unpack_from("<%dH" % m, data, 2))
            return ("NOT " + elems) if n < 0 else elems
        if kind == "order":
            v = _int(data, signed=False)
            if v is None:
                return None
            names = [name for bit, name in BOND_ORDER_BITS.items() if v & bit]
            return " ".join(names) if names else str(v)
        if kind.startswith("enum:"):
            v = _int(data)
            if v is None:
                return None
            table = ENUMS[kind[5:]]
            return table[v] if v in table else str(v)
        if kind.startswith("flags:"):
            v = _int(data, signed=False)
            if v is None:
                return None
            table = FLAGS[kind[6:]]
            names = [name for bit, name in table.items() if bit and v & bit]
            return " ".join(names) if names else table.get(0, str(v))
        return None

    # -- objects ------------------------------------------------------------
    def attributes(self, obj: CDXObject) -> Tuple[List[Tuple[str, str]], List[str]]:
        """(attribute list, extra child XML snippets) for one object."""
        attrs: List[Tuple[str, str]] = []
        extra: List[str] = []
        for p in obj.properties:
            spec = PROPERTY_TAGS.get(p.tag)
            if spec is None:
                self.unknown[p.tag] = self.unknown.get(p.tag, 0) + 1
                continue
            name, kind = spec
            try:
                self._convert(p, name, kind, attrs, extra)
            except (struct.error, IndexError, ValueError, KeyError):
                continue
        # A (corrupt) object may carry the same property tag twice; duplicate XML attributes are not well-formed, so
        # keep one value per name (the last one wins, as it would when ChemDraw reads the stream).
        return list(dict(attrs).items()), extra

    def _convert(self, p: CDXProperty, name: str, kind: str, attrs: List[Tuple[str, str]], extra: List[str]) -> None:
        if kind == "skip":
            return
        if kind == "text":
            extra.extend(self.style_runs_xml(p.data))
        elif kind == "fonttable":
            extra.append(self.fonttable_xml())
        elif kind == "colortable":
            extra.append(self.colortable_xml(p.data))
        elif kind == "fontstyle":
            if len(p.data) >= 8:
                font, face, size, _color = struct.unpack_from("<HHHH", p.data, 0)
                attrs.append((name + "Font", str(font)))
                attrs.append((name + "Size", _f(size / 20.0)))
                if face:
                    attrs.append((name + "Face", str(face)))
        elif kind == "represent":
            if len(p.data) >= 6:
                target, ptag = struct.unpack_from("<IH", p.data, 0)
                pname = PROPERTY_TAGS.get(ptag, ("0x%04X" % ptag, ""))[0]
                extra.append('<represent attribute=%s object="%d"/>' % (quoteattr(pname), target))
        elif kind == "blob":
            for magic, fmt in PICTURE_MAGIC:
                if p.data.startswith(magic):
                    attrs.append((fmt, p.data.hex().upper()))
                    break
        else:
            v = self.value(kind, p.data)
            if v is not None:
                attrs.append((name, v))

    def style_runs_xml(self, data: bytes) -> List[str]:
        text, runs = decode_cdx_string(data, self.fonts, self.default_codec)
        if not runs:
            return ["<s>%s</s>" % escape(_xml_text(text))]
        out = []
        for r in runs:
            a = [("font", str(r.font)), ("size", _f(r.size / 20.0))]
            if r.color != 3:          # 3 is ChemDraw's default text colour (black); it is omitted in CDXML
                a.append(("color", str(r.color)))
            if r.face:
                a.append(("face", str(r.face)))
            out.append("<s%s>%s</s>" % (_attr_str(a), escape(_xml_text(r.text))))
        return out

    def fonttable_xml(self) -> str:
        lines = ["<fonttable>"]
        for fid, (fname, charset) in self.fonts.items():
            cname = CHARSETS.get(charset, ("iso-8859-1", ""))[0]
            lines.append('<font id="%d" charset=%s name=%s/>' % (fid, quoteattr(cname), quoteattr(_xml_text(fname))))
        lines.append("</fonttable>")
        return "\n".join(lines)

    @staticmethod
    def colortable_xml(data: bytes) -> str:
        lines = ["<colortable>"]
        if len(data) >= 2:
            n = min(_u16(data, 0), (len(data) - 2) // 6)
            for i in range(n):
                rgb = struct.unpack_from("<HHH", data, 2 + 6 * i)
                # ChemDraw writes 4 decimals for non-integral components (0x8080 -> "0.5020")
                vals = [str(v // 65535) if v % 65535 == 0 else "%.4f" % (v / 65535.0) for v in rgb]
                lines.append('<color r="%s" g="%s" b="%s"/>' % tuple(vals))
        lines.append("</colortable>")
        return "\n".join(lines)

    def emit(self, obj: CDXObject, out: List[str]) -> None:
        name = obj.name
        attrs, extra = self.attributes(obj)
        if obj.id:
            attrs.insert(0, ("id", str(obj.id)))
        if not extra and not obj.children:
            out.append("<%s%s/>" % (name, _attr_str(attrs)))
            return
        if obj.tag == 0x8006 and not obj.children:
            out.append("<%s%s>%s</%s>" % (name, _attr_str(attrs), "".join(extra), name))
            return
        out.append("<%s%s>" % (name, _attr_str(attrs)))
        out.extend(extra)
        for c in obj.children:
            self.emit(c, out)
        out.append("</%s>" % name)

    def render(self) -> str:
        attrs, extra = self.attributes(self.doc)
        out = ['<?xml version="1.0" encoding="UTF-8" ?>',
               '<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >',
               "<CDXML%s>" % _attr_str(attrs)]
        out.extend(extra)
        for c in self.doc.children:
            self.emit(c, out)
        out.append("</CDXML>")
        return "\n".join(out) + "\n"


def cdx_to_cdxml(data: bytes) -> str:
    """Convert CDX bytes into CDXML text (parseable by :func:`chemdrawmcp.cdxml.read_cdxml` and by ChemDraw)."""
    return _Converter(parse_cdx(data)).render()


def read_cdx(data: bytes, want_mols: bool = True) -> Dict[str, Any]:
    """Structured description of a CDX document (same layout as :func:`chemdrawmcp.cdxml.read_cdxml`)."""
    from .cdxml import read_cdxml
    xml = cdx_to_cdxml(data)
    try:
        return read_cdxml(xml, want_mols=want_mols)
    except ET.ParseError as e:
        raise CDXError("converted CDXML is not well-formed: %s" % e) from e


# ---------------------------------------------------------------------------
# debugging
# ---------------------------------------------------------------------------

def dump_cdx(data: bytes, max_depth: Optional[int] = None) -> str:
    """Human-readable dump of the object tree with decoded properties (unknown tags shown as hex)."""
    doc = parse_cdx(data)
    conv = _Converter(doc)
    lines: List[str] = []

    def rec(obj: CDXObject, depth: int) -> None:
        attrs, _extra = conv.attributes(obj)
        parts = ["%s id=%d" % (obj.name, obj.id)]
        parts += ["%s=%s" % (k, v if len(v) <= 60 else v[:57] + "...") for k, v in attrs]
        txt = obj.prop(0x0700)
        if txt is not None:
            parts.append("Text=%r" % conv.string(txt))
        for p in obj.properties:
            if p.tag not in PROPERTY_TAGS:
                parts.append("0x%04X=<%d bytes%s>" % (p.tag, len(p.data), (" " + p.data[:8].hex()) if p.data else ""))
        lines.append("  " * depth + " ".join(parts))
        if max_depth is not None and depth >= max_depth:
            if obj.children:
                lines.append("  " * (depth + 1) + "... %d children" % len(obj.children))
            return
        for c in obj.children:
            rec(c, depth + 1)

    rec(doc, 0)
    return "\n".join(lines)
