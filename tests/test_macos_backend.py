"""Headless tests for the macOS (JXA) backend.  Nothing here launches or talks to ChemDraw.

* The JXA scripts the backend generates (``JS_PRELUDE`` + body) are executed with **node** against a fake
  ChemDraw object model (documents, windows, commands, an NSPasteboard-like ``$`` bridge) defined in ``FAKE_JS``.
  Those tests are skipped when node is not installed.
* ``save_document`` is exercised with ``js`` stubbed to write (or not write) files.
* The pasteboard helpers (``JS_PASTEBOARD``) are run for real with ``osascript`` on macOS, on a private *named*
  pasteboard; the general pasteboard (the user's clipboard) is never touched.
"""
from __future__ import annotations

import base64
import json
import os
import shutil
import subprocess
import sys
import time

import pytest

from chemdrawmcp.backends.base import ChemDrawError
from chemdrawmcp.backends.macos import JS_PASTEBOARD, MacBackend

NODE = shutil.which("node")
needs_node = pytest.mark.skipif(not NODE, reason="node is needed to run the backend's JXA scripts against a fake ChemDraw")


def b64(text: str) -> str:
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


# ---------------------------------------------------------------------------
# a fake ChemDraw + JXA environment for node
# ---------------------------------------------------------------------------

FAKE_JS = r"""
const vm = require("vm"), fs = require("fs"), pathmod = require("path");
const SCRIPT = fs.readFileSync(process.argv[2], "utf8");
const WORLD = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
// WORLD = {docs: [{name, path, wid, atoms}], front: <wid>, pasteboard: [{type: base64}], commands: {name: enabled}, copy_result}
const CALLS = [];
globalThis.Path = function (p) { return {toString() { return String(p); }}; };
globalThis.delay = function () {};
function rec(wid) { return WORLD.docs.find(r => r.wid === wid); }
// windows / documents are ordered front to back, like ChemDraw's
function ordered() { return WORLD.docs.slice().sort((a, b) => (a.wid === WORLD.front ? -1 : 0) - (b.wid === WORLD.front ? -1 : 0)); }
function mkDoc(r) {
  return {
    name: () => r.name, modified: () => false, file: () => (r.path ? Path(r.path) : null), window: () => mkWin(r),
    atoms: {get length() { return r.atoms || 0; }}, bonds: {length: 0}, molecules: {get length() { return r.atoms ? 1 : 0; }},
    captions: {length: 0}, arrows: {length: 0},
    selection: {smiles: () => { throw new Error("no selection"); }, bounds: () => { throw new Error("no selection"); }},
    close: (o) => { CALLS.push(["close", r.name, o]); },
  };
}
function mkWin(r) {
  const w = {id: () => r.wid, document: () => mkDoc(r)};
  Object.defineProperty(w, "index", {set(v) { if (v === 1) WORLD.front = r.wid; }});
  return w;
}
function unb64(s) { return Buffer.from(s, "base64").toString("utf8"); }
function b64(s) { return Buffer.from(String(s), "utf8").toString("base64"); }
const NIL = {isNil: () => true};
function nsdata(b) { return {b64: b, isNil: () => false, base64EncodedStringWithOptions: () => b}; }
function nsstr(s) { const o = new String(s); o.isNil = () => false; return o; }
function mkItem(entry) {
  return {
    get types() { return Object.keys(entry); },
    dataForType: (t) => (t in entry ? nsdata(entry[t]) : NIL),
    setDataForType: (d, t) => { entry[t] = d.b64; return true; },
    _entry: entry,
  };
}
const PB = {   // the general pasteboard: a list of items, each {type: base64 data}
  items: (WORLD.pasteboard || []).map(e => Object.assign({}, e)),
  get pasteboardItems() { const its = this.items.map(mkItem); return {isNil: () => false, count: String(its.length), objectAtIndex: (i) => its[i]}; },
  get clearContents() { this.items = []; return 0; },
  setStringForType(s, t) { this.items.push({[t]: b64(s)}); return true; },
  stringForType(t) { for (const e of this.items) if (t in e) return nsstr(unb64(e[t])); return NIL; },
  writeObjects(arr) { for (const it of arr) this.items.push(it._entry); return true; },
  text() { const s = this.stringForType("public.utf8-plain-text"); return s.isNil() ? null : s.valueOf(); },
};
globalThis.ObjC = {import: () => {}, unwrap: (x) => (x instanceof String ? x.valueOf() : x), deepUnwrap: (x) => x};
globalThis.$ = Object.assign((x) => x, {
  NSPasteboard: {generalPasteboard: PB},
  NSPasteboardItem: {alloc: {get init() { return mkItem({}); }}},
  NSData: {alloc: {initWithBase64EncodedStringOptions: (s) => nsdata(s)}},
});
const documents = () => ordered().map(mkDoc);                     // app.documents() and app.documents.name() / .byName()
Object.defineProperty(documents, "name", {value: () => ordered().map(r => r.name)});
documents.byName = (n) => { const r = ordered().find(r => r.name === n); if (!r) throw new Error("fake: no document " + n); return mkDoc(r); };
const app = {
  includeStandardAdditions: false, frontmost: () => true, activate: () => {}, version: () => "fake 0.0",
  windows: () => ordered().map(mkWin),
  documents: documents,
  open: (p) => {
    const path = p.toString(); CALLS.push(["open", path]);
    const open = WORLD.docs.find(r => r.path === path);
    if (open) { WORLD.front = open.wid; return; }          // already open: ChemDraw brings its window forward
    const wid = Math.max(0, ...WORLD.docs.map(r => r.wid)) + 1;
    WORLD.docs.push({name: pathmod.basename(path), path: path, wid: wid, atoms: 0}); WORLD.front = wid;
  },
  commands: {whose: (q) => (q.name in (WORLD.commands || {})) ? [{name: () => q.name, enabled: () => !!WORLD.commands[q.name]}] : []},
  do: (c) => {
    const n = c.name(); CALLS.push(["do", n, PB.text()]);
    const r = rec(WORLD.front);
    if (n === "paste" || n.startsWith("pasteAs")) r.atoms += (PB.text() || "").length;
    if (n.startsWith("copyAs")) PB.items = [{"public.utf8-plain-text": b64(WORLD.copy_result || "")}];
  },
};
globalThis.Application = function (id) { if (id === "System Events") return {processes: {whose: () => []}}; return app; };
const out = vm.runInThisContext(SCRIPT, {filename: "jxa.js"});
process.stdout.write(String(out) + "\n" + JSON.stringify({world: WORLD, calls: CALLS, pasteboard: PB.items}) + "\n");
"""


class FakeChemDraw:
    """A ``MacBackend`` whose scripts run in node against the fake object model above instead of osascript."""

    def __init__(self, tmp_path, world: dict):
        self.dir = tmp_path
        self.world = {"docs": [], "front": None, **world}
        self.pasteboard = [dict(e) for e in world.get("pasteboard", [])]
        self.calls: list = []
        self.backend = MacBackend(app_path=str(tmp_path / "NoChemDraw.app"))   # never discovers or launches anything
        self.backend._run_js = self._run_js
        self.backend.dismiss_dialogs = lambda: []

    def _run_js(self, script: str, timeout=None):
        stub, src, world = self.dir / "fake.js", self.dir / "script.js", self.dir / "world.json"
        stub.write_text(FAKE_JS)
        src.write_text(script)
        world.write_text(json.dumps(dict(self.world, pasteboard=self.pasteboard)))
        r = subprocess.run([NODE, str(stub), str(src), str(world)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, r.stderr
        lines = r.stdout.strip().splitlines()
        payload, state = json.loads(lines[0]), json.loads(lines[1])
        self.world, self.pasteboard = state["world"], state["pasteboard"]
        self.calls.extend(state["calls"])
        if not payload.get("ok"):
            raise ChemDrawError(payload.get("error", "unknown error"))
        return payload.get("result")


def docs_world(names, front=1, **extra):
    return {"docs": [{"name": n, "path": None, "wid": i + 1, "atoms": 0} for i, n in enumerate(names)], "front": front, **extra}


# ---------------------------------------------------------------------------
# docByName: the fuzzy rules are tiers, not one OR-ed filter
# ---------------------------------------------------------------------------

@needs_node
def test_doc_by_name_cdxml_less_match_wins_over_a_prefix_sibling(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world(["scheme.cdxml", "scheme.cdx"]))
    assert fake.backend.activate_document("scheme") == "scheme.cdxml"
    assert fake.world["front"] == 1
    assert fake.backend.activate_document("SCHEME.CDX") == "scheme.cdx"      # case-insensitive
    assert fake.world["front"] == 2
    assert fake.backend.activate_document("scheme.cdx") == "scheme.cdx"      # exact


@needs_node
def test_doc_by_name_case_insensitive_match_wins_over_a_prefix_sibling(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world(["Untitled ACS Document 1996-1", "Untitled ACS Document 1996-10"]))
    assert fake.backend.activate_document("untitled acs document 1996-1") == "Untitled ACS Document 1996-1"
    assert fake.backend.activate_document("untitled acs document 1996-10") == "Untitled ACS Document 1996-10"
    fake = FakeChemDraw(tmp_path, docs_world(["scheme.cdxml", "scheme2.cdxml"]))
    assert fake.backend.activate_document("scheme") == "scheme.cdxml"
    assert fake.backend.activate_document("Scheme2") == "scheme2.cdxml"


@needs_node
def test_doc_by_name_unique_prefix_ambiguity_and_window_id(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world(["alpha.cdxml", "beta.cdxml"]))
    assert fake.backend.activate_document("alp") == "alpha.cdxml"
    assert fake.backend.activate_document(2) == "beta.cdxml"                # window id
    assert fake.backend.activate_document("#2") == "beta.cdxml"
    fake = FakeChemDraw(tmp_path, docs_world(["alpha1.cdxml", "alpha2.cdxml"]))
    with pytest.raises(ChemDrawError, match="ambiguous document name"):
        fake.backend.activate_document("alpha")
    with pytest.raises(ChemDrawError, match="no document named"):
        fake.backend.activate_document("gamma")


# ---------------------------------------------------------------------------
# open_document / documents / snapshot: documents are identified by window and file, not by name
# ---------------------------------------------------------------------------

@needs_node
def test_open_second_file_with_same_basename_targets_the_new_document(tmp_path):
    a, b = tmp_path / "a" / "scheme.cdxml", tmp_path / "b" / "scheme.cdxml"
    for p in (a, b):
        p.parent.mkdir()
        p.write_text("<CDXML/>")
    fake = FakeChemDraw(tmp_path, {"docs": [{"name": "scheme.cdxml", "path": str(a), "wid": 1, "atoms": 3}], "front": 1})
    assert fake.backend.open_document(str(b)) == "scheme.cdxml"
    assert fake.calls == [["open", str(b)]]
    assert fake.world["front"] == 2 and [d["path"] for d in fake.world["docs"]] == [str(a), str(b)]
    # only the opened document is the front one (both share the name)
    docs = fake.backend.documents()
    assert [(d.name, d.window_id, d.front, d.path) for d in docs] == [("scheme.cdxml", 2, True, str(b)), ("scheme.cdxml", 1, False, str(a))]
    snap = fake.backend.snapshot()
    assert [(d["window_id"], d["front"]) for d in snap["documents"]] == [(2, True), (1, False)]
    assert snap["front_document"] == "scheme.cdxml"


@needs_node
def test_open_already_open_file_brings_its_window_forward(tmp_path):
    a = tmp_path / "a" / "scheme.cdxml"
    a.parent.mkdir()
    a.write_text("<CDXML/>")
    fake = FakeChemDraw(tmp_path, {"docs": [{"name": "scheme.cdxml", "path": str(a), "wid": 1, "atoms": 0},
                                            {"name": "scheme.cdxml", "path": str(tmp_path / "elsewhere" / "scheme.cdxml"), "wid": 2, "atoms": 0},
                                            {"name": "other.cdxml", "path": None, "wid": 3, "atoms": 0}], "front": 3})
    assert fake.backend.open_document(str(a)) == "scheme.cdxml"
    assert len(fake.world["docs"]) == 3
    assert fake.world["front"] == 1                                         # the one with that file, not the other scheme.cdxml
    assert [d.front for d in fake.backend.documents()] == [True, False, False]
    assert fake.backend.documents()[0].window_id == 1


@needs_node
def test_open_document_missing_file(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world([]))
    with pytest.raises(ChemDrawError, match="file not found"):
        fake.backend.open_document(str(tmp_path / "nope.cdxml"))
    assert fake.calls == []


# ---------------------------------------------------------------------------
# insert / selection_as: the whole clipboard survives, not only its text
# ---------------------------------------------------------------------------

USER_CLIPBOARD = [
    {"public.png": b64("\x89PNG not really"), "com.cambridgesoft.cdx": b64("VjCD binary-ish")},
    {"public.utf8-plain-text": b64("hello"), "public.rtf": b64("{\\rtf1 hello}")},
]


@needs_node
@pytest.mark.parametrize("clipboard", [USER_CLIPBOARD, []], ids=["image+cdx+rtf", "empty"])
def test_insert_pastes_the_payload_and_restores_every_flavour(tmp_path, clipboard):
    fake = FakeChemDraw(tmp_path, docs_world(["Untitled ACS Document 1996-1"], pasteboard=clipboard, commands={"pasteAsSMILES": True}))
    res = fake.backend.insert(None, "CCO", "smiles")
    assert res["document"] == "Untitled ACS Document 1996-1" and res["kind"] == "smiles" and res["atoms_added"] == 3
    assert [c for c in fake.calls if c[0] == "do"] == [["do", "pasteAsSMILES", "CCO"]]   # the paste saw the payload
    assert fake.pasteboard == clipboard                                                 # and the user's clipboard is back


@needs_node
def test_insert_restores_the_clipboard_when_the_paste_command_is_disabled(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world(["doc.cdxml"], pasteboard=USER_CLIPBOARD, commands={"pasteAsCdxml": False}))
    with pytest.raises(ChemDrawError, match="pasteAsCdxml is not enabled"):
        fake.backend.insert("doc.cdxml", "<CDXML/>", "cdxml")
    assert fake.pasteboard == USER_CLIPBOARD
    with pytest.raises(ChemDrawError, match="kind must be one of"):
        fake.backend.insert("doc.cdxml", "x", "bogus")


@needs_node
def test_selection_as_reads_the_copy_and_restores_the_clipboard(tmp_path):
    fake = FakeChemDraw(tmp_path, docs_world(["doc.cdxml"], pasteboard=USER_CLIPBOARD, commands={"copyAsSMILES": True, "copyAsMol": False},
                                             copy_result="c1ccccc1\r\n"))
    assert fake.backend.selection_as(None, "smiles") == "c1ccccc1\n"
    assert [c for c in fake.calls if c[0] == "do"] == [["do", "copyAsSMILES", None]]     # pasteboard cleared before the copy
    assert fake.pasteboard == USER_CLIPBOARD
    with pytest.raises(ChemDrawError, match="copyAsMol is not enabled"):
        fake.backend.selection_as(None, "mol")
    assert fake.pasteboard == USER_CLIPBOARD


# ---------------------------------------------------------------------------
# save_document: the "ChemDraw appended its own extension" fallback
# ---------------------------------------------------------------------------

def headless_saver(tmp_path, write=None):
    """A backend whose ``js`` writes ``write(path)`` (or nothing) instead of asking ChemDraw to save."""
    b = MacBackend(app_path=str(tmp_path / "NoChemDraw.app"))
    b.js = lambda body, timeout=None, **args: (write(args["path"]) if write else None) or True
    b.dismiss_dialogs = lambda: []
    return b


def test_save_document_exact_path(tmp_path):
    b = headless_saver(tmp_path, lambda p: open(p, "w").write("<svg/>"))
    assert b.save_document(None, str(tmp_path / "scheme.svg"), "svg") == str(tmp_path / "scheme.svg")
    assert b.save_document(None, str(tmp_path / "scheme"), "png") == str(tmp_path / "scheme.png")    # extension added


@pytest.mark.parametrize("fmt, ext", [("svg", "svg"), ("rxn", "rxn"), ("png", "png")])
def test_save_document_fallback_never_returns_a_sibling_of_another_format(tmp_path, fmt, ext):
    (tmp_path / "scheme.cdxml").write_text("<CDXML/>")          # written a moment ago, like a previous export
    b = headless_saver(tmp_path)                                  # ChemDraw "saves" without writing anything
    with pytest.raises(ChemDrawError, match="ChemDraw did not write"):
        b.save_document(None, str(tmp_path / ("scheme." + ext)), fmt)


@pytest.mark.parametrize("folder", ["run1", "run[1]", "a b (c)"])
def test_save_document_fallback_accepts_an_appended_extension(tmp_path, folder):
    d = tmp_path / folder
    d.mkdir()
    b = headless_saver(tmp_path, lambda p: open(p + ".png", "w").write("png"))     # ChemDraw writes x.png.png
    assert b.save_document(None, str(d / "x.png"), "png") == str(d / "x.png.png")


def test_save_document_fallback_accepts_the_formats_other_spelling(tmp_path):
    b = headless_saver(tmp_path, lambda p: open(os.path.splitext(p)[0] + ".jpeg", "w").write("jpg"))   # x.jpg -> x.jpeg
    assert b.save_document(None, str(tmp_path / "x.jpg"), "jpeg") == str(tmp_path / "x.jpeg")
    b = headless_saver(tmp_path, lambda p: open(os.path.splitext(p)[0] + ".tif", "w").write("tif"))
    assert b.save_document(None, str(tmp_path / "y.tiff"), "tiff") == str(tmp_path / "y.tif")


def test_save_document_fallback_ignores_a_stale_file_of_the_right_format(tmp_path):
    stale = tmp_path / "x.png.png"
    stale.write_text("old")
    old = time.time() - 120
    os.utime(stale, (old, old))
    b = headless_saver(tmp_path)
    with pytest.raises(ChemDrawError, match="ChemDraw did not write"):
        b.save_document(None, str(tmp_path / "x.png"), "png")
    assert stale.exists()


# ---------------------------------------------------------------------------
# the real pasteboard helpers, on a private named pasteboard (macOS only)
# ---------------------------------------------------------------------------

PASTEBOARD_PROBE = JS_PASTEBOARD + r"""
(function () {
  const pb = $.NSPasteboard.pasteboardWithName("org.chemdrawmcp.test." + $.NSProcessInfo.processInfo.processIdentifier);
  const out = {};
  try {
    const a = $.NSPasteboardItem.alloc.init;
    a.setDataForType($.NSData.alloc.initWithBase64EncodedStringOptions("iVBORw0KGgo=", 0), "public.png");
    a.setDataForType($.NSData.alloc.initWithBase64EncodedStringOptions("AAECAwQF/w==", 0), "com.cambridgesoft.cdx");
    const b = $.NSPasteboardItem.alloc.init;
    b.setStringForType("hello é world", PB_TEXT);
    b.setStringForType("{\\rtf1 hi}", "public.rtf");
    pb.clearContents; pb.writeObjects($([a, b]));
    const saved = pbSnapshot(pb);
    out.saved = saved;
    pbSetText(pb, "C1=CC=CC=C1 é");
    out.items_while_pasting = Number(pb.pasteboardItems.count);
    out.read_back = pbGetText(pb);
    pbRestore(pb, saved);
    out.restored = pbSnapshot(pb);
    out.restored_types = ObjC.deepUnwrap(pb.types);
    pb.clearContents;
    pbRestore(pb, pbSnapshot(pb));
    out.empty_items = Number(pb.pasteboardItems.count);
    out.empty_text = pbGetText(pb);
  } finally { pb.releaseGlobally; }
  return JSON.stringify(out);
})();
"""


@pytest.mark.skipif(sys.platform != "darwin" or not shutil.which("osascript"), reason="needs macOS (JXA + NSPasteboard)")
def test_pasteboard_helpers_round_trip_every_flavour_jxa():
    r = subprocess.run(["osascript", "-l", "JavaScript", "-"], input=PASTEBOARD_PROBE, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    out = json.loads(r.stdout.strip().splitlines()[-1])
    assert len(out["saved"]) == 2
    assert out["saved"][0]["public.png"] == "iVBORw0KGgo=" and out["saved"][0]["com.cambridgesoft.cdx"] == "AAECAwQF/w=="
    assert base64.b64decode(out["saved"][1]["public.utf8-plain-text"]).decode("utf-8") == "hello é world"
    assert "public.rtf" in out["saved"][1]
    assert out["items_while_pasting"] == 1 and out["read_back"] == "C1=CC=CC=C1 é"
    assert out["restored"] == out["saved"]
    assert {"public.png", "com.cambridgesoft.cdx", "public.utf8-plain-text", "public.rtf"} <= set(out["restored_types"])
    assert out["empty_items"] == 0 and out["empty_text"] == ""
