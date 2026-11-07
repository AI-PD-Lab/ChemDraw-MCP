"""Tests for the binary CDX reader (chemdrawmcp.cdx).

Two layers:

* hand-assembled CDX byte strings (always run, also on CI);
* round-trip comparisons against documents that ChemDraw 23 saved in both formats.  The pair directory is
  taken from ``CHEMDRAWMCP_CDX_PAIRS`` (default ``~/.chemdrawmcp/testdata/cdxpairs``); the tests are skipped
  when no ``X.cdx`` + ``X.cdxml`` pairs are found there.
"""
import collections
import glob
import os
import re
import struct
import xml.etree.ElementTree as ET

import pytest

from chemdrawmcp.cdx import (CDXError, CDXObject, cdx_to_cdxml, decode_cdx_string, dump_cdx, is_cdx, parse_cdx,
                             read_cdx)
from chemdrawmcp.cdxml import read_cdxml

_HOME_PAIRS = os.path.join(os.path.expanduser("~"), ".chemdrawmcp", "testdata", "cdxpairs")
PAIR_DIR = os.environ.get("CHEMDRAWMCP_CDX_PAIRS") or _HOME_PAIRS
PAIRS = sorted(p for p in glob.glob(os.path.join(PAIR_DIR, "*.cdx")) if os.path.exists(p[:-4] + ".cdxml"))


# ---------------------------------------------------------------------------
# tiny hand-assembled document
# ---------------------------------------------------------------------------

def _prop(tag: int, data: bytes) -> bytes:
    return struct.pack("<HH", tag, len(data)) + data


def _long_prop(tag: int, data: bytes) -> bytes:
    """Property written with the 0xFFFF + UINT32 length form."""
    return struct.pack("<HHI", tag, 0xFFFF, len(data)) + data


def _obj(tag: int, oid: int, body: bytes = b"") -> bytes:
    return struct.pack("<HI", tag, oid) + body + b"\x00\x00"


def _coord(v: float) -> bytes:
    return struct.pack("<i", int(round(v * 65536)))


def _pos(x: float, y: float) -> bytes:          # 2D positions are stored y first
    return _coord(y) + _coord(x)


def _rect(left: float, top: float, right: float, bottom: float) -> bytes:
    return _coord(top) + _coord(left) + _coord(bottom) + _coord(right)


def _cdxstring(text: str, runs=()) -> bytes:
    out = struct.pack("<H", len(runs))
    for start, font, face, size, color in runs:
        out += struct.pack("<HHHHH", start, font, face, size, color)
    return out + text.encode("cp1252")


def _u16(v: int) -> bytes:
    return struct.pack("<H", v)


def _u32(v: int) -> bytes:
    return struct.pack("<I", v)


def tiny_cdx(bond_extra: bytes = b"", reserved: int = 10, extra_page_objects: bytes = b"",
             extra_doc_props: bytes = b"", num_h: int = 1) -> bytes:
    """Methanol (O-C, one fragment, two nodes, one bond), a caption, an arrow and a reaction step."""
    header = b"VjCD0100" + b"\x04\x03\x02\x01" + b"\x00" * reserved
    fonttable = struct.pack("<HH", 1, 1) + struct.pack("<HHH", 3, 1252, 5) + b"Arial"
    colortable = _u16(3) + struct.pack("<HHH", 0xFFFF, 0xFFFF, 0xFFFF) + struct.pack("<HHH", 0, 0, 0) \
        + struct.pack("<HHH", 0x8080, 0x8080, 0x8080)
    doc_props = (_prop(0x0003, _cdxstring("tiny test")) + _prop(0x0008, _cdxstring("tiny.cdx"))
                 + _prop(0x0204, _rect(0, 0, 200, 100)) + _prop(0x0100, fonttable) + _prop(0x0300, colortable)
                 + _prop(0x0805, _coord(30)) + _prop(0x0807, _coord(1)) + _prop(0x0804, _u16(120))
                 + _prop(0x080A, struct.pack("<HHHH", 3, 96, 200, 3)) + _prop(0x080B, struct.pack("<HHHH", 3, 0, 240, 3))
                 + _prop(0x0708, b"\x01") + extra_doc_props)
    label = _obj(0x8006, 0, _prop(0x0200, _pos(52, 54)) + _prop(0x0701, b"\xff")
                 + _prop(0x0700, _cdxstring("OH", [(0, 3, 96, 200, 3)])))
    node_o = _obj(0x8004, 11, _prop(0x0200, _pos(50, 100)) + _prop(0x0402, _u16(8)) + _prop(0x042B, _u16(num_h))
                  + _prop(0x0437, b"\x01") + _prop(0x044D, _u32(1)) + label)
    node_c = _obj(0x8004, 12, _prop(0x0200, _pos(80, 100)) + _prop(0x000A, _u16(5)) + _prop(0x000F, b""))
    bond = _obj(0x8005, 13, _prop(0x0604, _u32(11)) + _prop(0x0605, _u32(12)) + _prop(0x060A, b"\x01") + bond_extra)
    frag = _obj(0x8003, 10, _prop(0x0204, _rect(50, 95, 80, 105)) + node_o + node_c + bond)
    caption = _obj(0x8006, 20, _prop(0x0200, _pos(50, 120)) + _prop(0x0204, _rect(50, 110, 90, 122))
                   + _prop(0x0701, b"\x01") + _prop(0x0700, _cdxstring("methanol\rline two", [(0, 3, 0, 240, 3)])))
    arrow = _obj(0x8021, 30, _prop(0x0204, _rect(100, 95, 140, 105)) + _prop(0x0207, _coord(140) + _coord(100) + _coord(0))
                 + _prop(0x0208, _coord(100) + _coord(100) + _coord(0)) + _prop(0x0A35, _u16(2)) + _prop(0x0A2F, _u16(1))
                 + _prop(0x0A20, _u16(1500)))
    step = _obj(0x800E, 41, _prop(0x0C01, _u32(10)) + _prop(0x0C04, _u32(30)) + _prop(0x0C05, _u32(20)))
    scheme = _obj(0x800D, 40, step)
    page = _obj(0x8001, 1, _prop(0x0204, _rect(0, 0, 540, 720)) + _prop(0x0819, b"") + frag + caption + arrow + scheme
                + extra_page_objects)
    return header + _obj(0x8000, 0, doc_props + page)


def _root(xml: str) -> ET.Element:
    return ET.fromstring(re.sub(r"<!DOCTYPE[^>]*>", "", xml))


# ---------------------------------------------------------------------------
# unit tests (run everywhere)
# ---------------------------------------------------------------------------

def test_is_cdx():
    assert is_cdx(tiny_cdx())
    assert is_cdx(b"VjCD0100")
    assert not is_cdx(b"<?xml version='1.0'?><CDXML/>")
    assert not is_cdx(b"")
    assert not is_cdx(None)


def test_parse_rejects_non_cdx():
    with pytest.raises(CDXError):
        parse_cdx(b"this is not a cdx file")
    assert issubclass(CDXError, ValueError)


def test_parse_tree_structure():
    doc = parse_cdx(tiny_cdx())
    assert isinstance(doc, CDXObject)
    assert doc.tag == 0x8000 and doc.name == "CDXML" and doc.id == 0
    assert [c.name for c in doc.children] == ["page"]
    page = doc.children[0]
    assert page.id == 1
    assert [c.name for c in page.children] == ["fragment", "t", "arrow", "scheme"]
    frag = page.children[0]
    assert frag.id == 10
    assert [(c.name, c.id) for c in frag.children] == [("n", 11), ("n", 12), ("b", 13)]
    node_o = frag.children[0]
    assert node_o.prop(0x0402) == b"\x08\x00"
    assert node_o.prop(0x9999) is None
    assert [c.name for c in node_o.children] == ["t"]        # the atom label lives inside the node
    assert [p.name for p in node_o.properties][:3] == ["p", "Element", "NumHydrogens"]
    assert len(doc.find_all(0x8004)) == 2
    assert page.children[3].children[0].name == "step"


def test_header_with_16_reserved_bytes_is_accepted():
    data = tiny_cdx(reserved=16)
    assert parse_cdx(data).children[0].name == "page"
    assert read_cdx(data)["molecules"][0]["smiles"] == "CO"


def test_cdxml_output_content():
    xml = cdx_to_cdxml(tiny_cdx())
    assert xml.startswith('<?xml version="1.0" encoding="UTF-8" ?>')
    assert "<!DOCTYPE CDXML" in xml
    root = _root(xml)
    assert root.tag == "CDXML"
    assert root.get("CreationProgram") == "tiny test"
    assert root.get("BoundingBox") == "0 0 200 100"
    assert root.get("BondLength") == "30"
    assert root.get("BondSpacing") == "12"
    assert root.get("LabelFont") == "3" and root.get("LabelSize") == "10" and root.get("LabelFace") == "96"
    assert root.get("CaptionFont") == "3" and root.get("CaptionSize") == "12" and root.get("CaptionFace") is None
    assert root.get("InterpretChemically") == "yes"
    fonts = root.findall("fonttable/font")
    assert [(f.get("id"), f.get("charset"), f.get("name")) for f in fonts] == [("3", "iso-8859-1", "Arial")]
    colors = root.findall("colortable/color")
    assert [(c.get("r"), c.get("g"), c.get("b")) for c in colors] == [("1", "1", "1"), ("0", "0", "0"),
                                                                       ("0.5020", "0.5020", "0.5020")]
    page = root.find("page")
    assert page.get("id") == "1" and page.get("BoundingBox") == "0 0 540 720" and page.get("PrintTrimMarks") == "yes"
    frag = page.find("fragment")
    assert frag.get("id") == "10" and frag.get("BoundingBox") == "50 95 80 105"
    nodes = frag.findall("n")
    assert nodes[0].attrib == {"id": "11", "p": "50 100", "Element": "8", "NumHydrogens": "1", "AS": "N", "AtomID": "1"}
    assert nodes[1].attrib == {"id": "12", "p": "80 100", "Z": "5", "IgnoreWarnings": "yes"}
    label = nodes[0].find("t")
    assert label.get("id") is None                       # id 0 is not written, like ChemDraw
    assert label.get("p") == "52 54" and label.get("Justification") == "Right"
    assert [(s.attrib, s.text) for s in label.findall("s")] == [({"font": "3", "size": "10", "face": "96"}, "OH")]
    bond = frag.find("b")
    assert bond.attrib == {"id": "13", "B": "11", "E": "12", "BS": "N"}
    caption = page.find("t")
    assert caption.get("id") == "20" and caption.get("p") == "50 120" and caption.get("Justification") == "Center"
    assert caption.get("BoundingBox") == "50 110 90 122"
    assert [(s.attrib, s.text) for s in caption.findall("s")] == [({"font": "3", "size": "12"}, "methanol\nline two")]
    arrow = page.find("arrow")
    assert arrow.get("id") == "30"
    assert arrow.get("Head3D") == "140 100 0" and arrow.get("Tail3D") == "100 100 0"
    assert arrow.get("ArrowheadHead") == "Full" and arrow.get("ArrowheadType") == "Solid" and arrow.get("HeadSize") == "1500"
    step = page.find("scheme/step")
    assert step.attrib == {"id": "41", "ReactionStepReactants": "10", "ReactionStepArrows": "30",
                           "ReactionStepObjectsAboveArrow": "20"}


def test_read_cdx_round_trip():
    info = read_cdx(tiny_cdx())
    assert info["pages"] == 1
    assert len(info["molecules"]) == 1
    assert info["molecules"][0]["smiles"] == "CO"
    assert info["molecules"][0]["formula"] == "CH4O"
    assert [t["text"] for t in info["texts"]] == ["methanol\nline two"]
    assert len(info["arrows"]) == 1 and info["arrows"][0]["head3D"] == "140 100 0"
    assert len(info["steps"]) == 1
    assert info["steps"][0]["reactants"][0]["smiles"] == "CO"
    assert info["steps"][0]["above_arrow"][0]["text"] == "methanol\nline two"
    assert info["reaction_smiles"] == ["CO>>"]
    assert info["settings"]["BondLength"] == 30.0
    assert info["settings"]["LabelFont"] == "Arial"
    assert info["settings"]["CreationProgram"] == "tiny test"
    assert read_cdx(tiny_cdx(), want_mols=False)["molecules"][0].get("smiles") is None


@pytest.mark.parametrize("order, display, expected", [
    (_u16(2), None, {"Order": "2"}),
    (_u16(4), None, {"Order": "3"}),
    (_u16(128), None, {"Order": "1.5"}),
    (_u16(4096), None, {"Order": "dative"}),
    (_u16(3), None, {"Order": "1 2"}),                       # query bond: single or double
    (None, _u16(6), {"Display": "WedgeBegin"}),
    (None, _u16(3), {"Display": "WedgedHashBegin"}),
    (None, _u16(8), {"Display": "Wavy"}),
    (None, _u16(5), {"Display": "Bold"}),
])
def test_bond_order_and_display(order, display, expected):
    extra = b""
    if order is not None:
        extra += _prop(0x0600, order)
    if display is not None:
        extra += _prop(0x0601, display)
    bond = _root(cdx_to_cdxml(tiny_cdx(bond_extra=extra))).find("page/fragment/b")
    for k, v in expected.items():
        assert bond.get(k) == v


def test_bond_order_feeds_rdkit():
    info = read_cdx(tiny_cdx(bond_extra=_prop(0x0600, _u16(2)), num_h=0))
    assert info["molecules"][0]["smiles"] == "C=O"
    info = read_cdx(tiny_cdx(bond_extra=_prop(0x0600, _u16(2))))       # explicit NumHydrogens=1 is kept
    assert info["molecules"][0]["smiles"] == "C=[OH]"


def test_node_attributes_and_nested_fragment():
    # an abbreviation node (NodeType Fragment) with a nested fragment (ethyl: CH3-CH2-*) and a charged, isotopic atom
    inner_c1 = _obj(0x8004, 101, _prop(0x0200, _pos(0, 0)))
    inner_c2 = _obj(0x8004, 102, _prop(0x0200, _pos(10, 0)))
    inner_ext = _obj(0x8004, 103, _prop(0x0200, _pos(20, 0)) + _prop(0x0400, _u16(12)) + _prop(0x044B, b"\x01"))
    inner_b1 = _obj(0x8005, 104, _prop(0x0604, _u32(101)) + _prop(0x0605, _u32(102)))
    inner_b2 = _obj(0x8005, 105, _prop(0x0604, _u32(102)) + _prop(0x0605, _u32(103)))
    nested = _obj(0x8003, 100, inner_c1 + inner_c2 + inner_ext + inner_b1 + inner_b2)
    et_label = _obj(0x8006, 0, _prop(0x0700, _cdxstring("Et", [(0, 3, 96, 200, 3)])))
    abbrev = _obj(0x8004, 51, _prop(0x0200, _pos(30, 0)) + _prop(0x0400, _u16(5)) + nested + et_label)
    nplus = _obj(0x8004, 52, _prop(0x0200, _pos(0, 0)) + _prop(0x0402, _u16(7)) + _prop(0x0421, b"\x01")
                 + _prop(0x0420, _u16(15)) + _prop(0x042B, _u16(3)) + _prop(0x0422, b"\x01") + _prop(0x0430, b"\x07"))
    rnode = _obj(0x8004, 53, _prop(0x0200, _pos(0, 30)) + _prop(0x0400, _u16(7)) + _prop(0x0433, _cdxstring("R")))
    b1 = _obj(0x8005, 54, _prop(0x0604, _u32(52)) + _prop(0x0605, _u32(51)))
    b2 = _obj(0x8005, 55, _prop(0x0604, _u32(52)) + _prop(0x0605, _u32(53)))
    frag = _obj(0x8003, 50, abbrev + nplus + rnode + b1 + b2)
    xml = cdx_to_cdxml(tiny_cdx(extra_page_objects=frag))
    root = _root(xml)
    frags = root.findall("page/fragment")
    assert [f.get("id") for f in frags] == ["10", "50"]
    n51 = frags[1].find("n[@id='51']")
    assert n51.get("NodeType") == "Fragment"
    assert n51.find("fragment").get("id") == "100"
    assert [c.tag for c in n51] == ["fragment", "t"]
    assert n51.find("fragment").find("n[@id='103']").get("NodeType") == "ExternalConnectionPoint"
    assert n51.find("fragment").find("n[@id='103']").get("ExternalConnectionNum") == "1"
    n52 = frags[1].find("n[@id='52']")
    assert n52.attrib == {"id": "52", "p": "0 0", "Element": "7", "Charge": "1", "Isotope": "15", "NumHydrogens": "3",
                          "Radical": "Singlet", "Geometry": "Tetrahedral"}
    n53 = frags[1].find("n[@id='53']")
    assert n53.get("NodeType") == "GenericNickname" and n53.get("GenericNickname") == "R"
    info = read_cdxml(xml)
    assert len(info["molecules"]) == 2
    mol2 = info["molecules"][1]
    assert mol2["atoms"] == 3
    assert "[15NH3+]" in mol2["smiles"] and "CC" in mol2["smiles"]     # abbreviation expanded to an ethyl group


def test_graphic_arrow_and_old_style_line_arrow():
    graphic = _obj(0x8007, 60, _prop(0x0204, _rect(100, 50, 140, 50)) + _prop(0x0A00, _u16(1)) + _prop(0x0A02, _u16(2))
                   + _prop(0x0A20, _u16(1500)) + _prop(0x0013, _u32(30)))
    symbol = _obj(0x8007, 61, _prop(0x0204, _rect(10, 10, 14, 14)) + _prop(0x0A00, _u16(7)) + _prop(0x0A07, _u16(5))
                  + _prop(0x000E, _u32(11) + _u16(0x0421)))
    rect = _obj(0x8007, 62, _prop(0x0204, _rect(0, 0, 10, 10)) + _prop(0x0A00, _u16(3)) + _prop(0x0A03, _u16(3))
                + _prop(0x0A3C, _u16(600)))
    arc = _obj(0x8021, 63, _prop(0x0A21, struct.pack("<h", -1200)) + _prop(0x0A36, _u16(2)) + _prop(0x0A37, _u16(1)))
    root = _root(cdx_to_cdxml(tiny_cdx(extra_page_objects=graphic + symbol + rect + arc)))
    g = root.find("page/graphic[@id='60']")
    assert g.get("GraphicType") == "Line" and g.get("ArrowType") == "FullHead" and g.get("HeadSize") == "1500"
    assert g.get("SupersededBy") == "30" and g.get("BoundingBox") == "100 50 140 50"
    s = root.find("page/graphic[@id='61']")
    assert s.get("GraphicType") == "Symbol" and s.get("SymbolType") == "CircleMinus"
    assert s.find("represent").attrib == {"attribute": "Charge", "object": "11"}
    r = root.find("page/graphic[@id='62']")
    assert r.get("RectangleType") == "RoundEdge Shadow" and r.get("CornerRadius") == "600"
    a = root.find("page/arrow[@id='63']")
    assert a.get("AngularSize") == "-120" and a.get("ArrowheadTail") == "Full" and a.get("FillType") == "None"
    info = read_cdxml(cdx_to_cdxml(tiny_cdx(extra_page_objects=graphic)))
    assert len(info["arrows"]) == 1  # the SupersededBy legacy line is hidden by read_cdxml          # the <arrow> and the superseded line graphic, as with ChemDraw's CDXML


def test_unknown_tags_are_tolerated():
    unknown_prop = _prop(0x7F77, b"\x01\x02\x03")
    unknown_obj = _obj(0x8077, 90, _prop(0x7F78, b"xyz") + _obj(0x8006, 91, _prop(0x0700, _cdxstring("inside"))))
    data = tiny_cdx(extra_page_objects=unknown_obj, extra_doc_props=unknown_prop)
    doc = parse_cdx(data)
    assert doc.prop(0x7F77) == b"\x01\x02\x03"
    assert doc.children[0].children[-1].name == "obj_8077"
    xml = cdx_to_cdxml(data)
    root = _root(xml)
    unk = root.find("page/obj_8077")
    assert unk.get("id") == "90" and unk.find("t").find("s").text == "inside"
    info = read_cdxml(xml)
    assert info["molecules"][0]["smiles"] == "CO"
    dump = dump_cdx(data)
    assert "0x7F77=<3 bytes 010203>" in dump and "obj_8077 id=90" in dump


def test_long_length_form_and_truncation():
    warning = "x" * 70000
    node_extra = _long_prop(0x0010, _cdxstring(warning))
    page_obj = _obj(0x8003, 70, _obj(0x8004, 71, _prop(0x0200, _pos(1, 1)) + node_extra))
    data = tiny_cdx(extra_page_objects=page_obj)
    root = _root(cdx_to_cdxml(data))
    assert root.find("page/fragment[@id='70']/n").get("Warning") == warning
    truncated = data[:len(data) - 30000]
    doc = parse_cdx(truncated)                      # no exception, partial tree
    assert doc.children[0].children[0].id == 10
    xml = cdx_to_cdxml(truncated)
    assert read_cdxml(xml)["molecules"][0]["smiles"] == "CO"


def test_decode_cdx_string():
    fonts = {5: ("Arial", 10000), 6: ("Symbol", 10000), 7: ("Arial", 1252)}
    data = struct.pack("<H", 1) + struct.pack("<HHHHH", 0, 5, 0, 200, 3) + b"1. BH3\xe1THF"
    text, runs = decode_cdx_string(data, fonts)
    assert text == "1. BH3·THF" and len(runs) == 1 and runs[0].text == text
    assert (runs[0].font, runs[0].face, runs[0].size, runs[0].color) == (5, 0, 200, 3)
    two = (struct.pack("<H", 2) + struct.pack("<HHHHH", 0, 7, 1, 200, 0) + struct.pack("<HHHHH", 3, 5, 64, 160, 0)
           + b"abc\xa9 d\rx")
    text, runs = decode_cdx_string(two, fonts)
    assert text == "abc© d\nx" and [r.text for r in runs] == ["abc", "© d\nx"]
    assert decode_cdx_string(b"stereo") == ("stereo", [])        # no run-count header (objecttag names)
    assert decode_cdx_string(b"\x00\x00") == ("", [])
    assert decode_cdx_string(b"") == ("", [])
    assert decode_cdx_string(b"\x00\x00caf\xe9", {}, "cp1252") == ("café", [])


def test_dump_cdx():
    dump = dump_cdx(tiny_cdx())
    lines = dump.splitlines()
    assert lines[0].startswith("CDXML id=0 CreationProgram=tiny test")
    assert any(ln.strip().startswith("fragment id=10") for ln in lines)
    assert any("n id=11 p=50 100 Element=8" in ln for ln in lines)
    assert any("Text='methanol\\nline two'" in ln for ln in lines)
    shallow = dump_cdx(tiny_cdx(), max_depth=1)
    assert "... 4 children" in shallow and "n id=11" not in shallow


# ---------------------------------------------------------------------------
# ChemDraw-written pairs (skipped when the sample directory is absent)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pair_results():
    results = {}
    for cdx_path in PAIRS:
        with open(cdx_path, "rb") as fh:
            data = fh.read()
        with open(cdx_path[:-4] + ".cdxml", encoding="utf-8") as fh:
            xml_ref = fh.read()
        results[os.path.basename(cdx_path)] = (read_cdxml(cdx_to_cdxml(data)), read_cdxml(xml_ref))
    return results


# One parameter per pair; with no pairs a single placeholder parameter carries a skip mark with a clear reason.
_NO_PAIRS = pytest.mark.skip(reason="no CDX/CDXML sample pairs in %s (set CHEMDRAWMCP_CDX_PAIRS)" % PAIR_DIR)
pairs_param = pytest.mark.parametrize("name", [os.path.basename(p) for p in PAIRS]
                                      or [pytest.param(None, marks=_NO_PAIRS)])


@pairs_param
def test_pair_molecule_count(pair_results, name):
    ours, ref = pair_results[name]
    assert len(ours["molecules"]) == len(ref["molecules"])
    assert ours["pages"] == ref["pages"]


@pairs_param
def test_pair_smiles(pair_results, name):
    ours, ref = pair_results[name]
    sa = collections.Counter(m.get("smiles") for m in ours["molecules"])
    sb = collections.Counter(m.get("smiles") for m in ref["molecules"])
    common = sum((sa & sb).values())
    rate = 1.0 - common / max(1, len(ref["molecules"]))
    print("%s: SMILES mismatch rate %.1f%% (%d/%d)" % (name, 100 * rate, len(ref["molecules"]) - common, len(ref["molecules"])))
    assert rate <= 0.10, "SMILES mismatch %.1f%%: only in cdxml %s / only in cdx %s" % (100 * rate, sb - sa, sa - sb)


@pairs_param
def test_pair_captions(pair_results, name):
    ours, ref = pair_results[name]
    ta = collections.Counter(t["text"] for t in ours["texts"])
    tb = collections.Counter(t["text"] for t in ref["texts"])
    if not tb:
        assert not ta
        return
    overlap = sum((ta & tb).values()) / sum(tb.values())
    assert overlap >= 0.90, "caption overlap %.1f%%: %s" % (100 * overlap, tb - ta)


@pairs_param
def test_pair_arrows_steps_bondlength(pair_results, name):
    ours, ref = pair_results[name]
    assert abs(len(ours["arrows"]) - len(ref["arrows"])) <= 1
    assert abs(len(ours["steps"]) - len(ref["steps"])) <= 1
    assert ours["settings"].get("BondLength") == ref["settings"].get("BondLength")
    for key in ("LabelSize", "CaptionSize", "LineWidth", "LabelFont", "CaptionFont"):
        assert ours["settings"].get(key) == ref["settings"].get(key)


@pairs_param
def test_pair_cdxml_is_well_formed_and_dumpable(name):
    with open(os.path.join(PAIR_DIR, name), "rb") as fh:
        data = fh.read()
    assert is_cdx(data)
    root = _root(cdx_to_cdxml(data))
    assert root.tag == "CDXML" and root.find("page") is not None
    assert root.get("CreationProgram", "").startswith("ChemDraw")
    assert "page id=" in dump_cdx(data, max_depth=1)


# ---------------------------------------------------------------------------
# corrupt-input robustness (every reader fails with CDXError or degrades gracefully, never with ParseError /
# RecursionError)
# ---------------------------------------------------------------------------

def test_xml_illegal_control_characters_are_stripped():
    fonttable = struct.pack("<HH", 1, 1) + struct.pack("<HHH", 3, 1252, 6) + b"Ar\x01ial"
    data = tiny_cdx(extra_page_objects=_obj(0x8006, 99, _prop(0x0700, _cdxstring("bad\x01\x0cchar\ttab")))
                    + _obj(0x8004, 98, _prop(0x0008, b"Me\x01thyl")),
                    extra_doc_props=_prop(0x0009, _cdxstring("com\x02ment")) + _prop(0x0100, fonttable))
    xml = cdx_to_cdxml(data)
    root = _root(xml)                               # well-formed: no ParseError
    assert root.find("page/t[@id='99']/s").text == "badchar\ttab"
    assert root.find("page/n[@id='98']").get("Name") == "Methyl"
    assert root.get("Comment") == "comment"
    assert root.find("fonttable/font").get("name") == "Arial"
    assert read_cdx(data)["molecules"][0]["smiles"] == "CO"


def test_duplicate_property_tag_yields_one_attribute():
    obj = _obj(0x8006, 98, _prop(0x0200, _pos(1, 1)) + _prop(0x0200, _pos(2, 2)) + _prop(0x0700, _cdxstring("x")))
    data = tiny_cdx(extra_page_objects=obj)
    root = _root(cdx_to_cdxml(data))                # no "duplicate attribute" ParseError
    assert root.find("page/t[@id='98']").get("p") == "2 2"      # last occurrence wins
    assert read_cdx(data)["molecules"][0]["smiles"] == "CO"
    assert "p=2 2" in dump_cdx(data)


def test_read_cdx_wraps_parse_errors_as_cdxerror(monkeypatch):
    import chemdrawmcp.cdx as cdx_module
    monkeypatch.setattr(cdx_module, "cdx_to_cdxml", lambda data: "<CDXML><page></CDXML>")
    with pytest.raises(CDXError, match="not well-formed"):
        read_cdx(tiny_cdx())


def test_deep_nesting_raises_cdxerror_instead_of_recursionerror():
    def nested(n, closed=True):
        return b"".join(struct.pack("<HI", 0x8002, 1000 + i) for i in range(n)) + (b"\x00\x00" * n if closed else b"")
    ok = tiny_cdx(extra_page_objects=nested(200))
    assert len(parse_cdx(ok).find_all(0x8002)) == 200
    assert read_cdx(ok)["molecules"][0]["smiles"] == "CO"
    assert dump_cdx(ok).count("group id=") == 200
    for data in (tiny_cdx(extra_page_objects=nested(1200)), tiny_cdx(extra_page_objects=nested(1200, closed=False))):
        for fn in (parse_cdx, cdx_to_cdxml, read_cdx, dump_cdx):
            with pytest.raises(CDXError, match="nesting deeper than"):
                fn(data)
