"""macOS backend: drives ChemDraw through its AppleScript dictionary using JXA (JavaScript for Automation).

Verified against ChemDraw 23.1.2 (bundle id com.revvity.ChemDraw):

* ``Application(id)`` objects: ``documents`` (name, modified, file, window, atoms/bonds/molecules/captions/arrows,
  selection with smiles / molecularFormula / molecularWeight / exactMass / elementalAnalysis / characters / bounds,
  drawing settings), ``windows`` (ordered front to back), ``menus`` / ``menuItems`` and ``commands`` (about 1600
  named commands such as selectAll, cleanStructure, convertNameToStructure, copyAsSMILES, pasteAsCdxml ...).
* ``app.open(Path)``, ``app.save(doc, {in: Path, as: "<file type name>"})`` for every format in EXPORT_FORMATS,
  ``doc.close({saving: "no"})``, ``app.do(command | menuItem)``.
* ChemDraw is never activated: the target document is brought to the front inside ChemDraw (window index / AXRaise),
  which is enough for commands and pasteboard operations, so the user keeps keyboard focus (hotkeys would otherwise
  relabel selected atoms).  ``CHEMDRAWMCP_ALLOW_ACTIVATE=1`` re-enables activation as a last resort.
* Modal alerts (e.g. "No reaction center was found" when exporting RXN) block the Apple event; they are detected
  and dismissed through System Events accessibility, and their text is returned as the error message.
* Data goes in through the general pasteboard: the text is written with NSPasteboard (ObjC bridge, same flavour as
  ``pbcopy``) and read by ``pasteAsSMILES`` / ``pasteAsCdxml`` / ``pasteAsMol`` / ``pasteAsInChi`` /
  ``pasteAsStructureToName`` (name -> structure) / ``paste`` (caption); and out through ``copyAsSMILES`` etc. +
  NSPasteboard.  The user's clipboard is snapshotted with every item and flavour (images, files, ChemDraw's own
  CDX/CDXML ...) and put back afterwards, inside the same script (see JS_PASTEBOARD).
"""
from __future__ import annotations

import glob
import json
import os
import plistlib
import re
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Optional

from .base import (EXPORT_FORMATS, INSERT_KINDS, ChemDrawBackend, ChemDrawError, ChemDrawNotAvailable, CommandInfo,
                   DocumentInfo, normalize_format)

DEFAULT_BUNDLE_ID = "com.revvity.ChemDraw"
KNOWN_BUNDLE_IDS = ("com.revvity.ChemDraw", "com.perkinelmer.ChemDraw", "com.cambridgesoft.ChemDraw")
AUTOSAVE_PREFIX = "ChemDraw Auto-Save"

# JXA property names for document settings and the conversion to points / plain values.
# (value in the app, factor to points)  1440ths of an inch -> pt: /20 ; 20ths of a point -> pt: /20
SETTINGS = {
    "bond_length": ("fixedLength", 20.0),
    "line_width": ("lineWidth", 20.0),
    "bold_width": ("boldWidth", 20.0),
    "margin_width": ("marginWidth", 20.0),
    "hash_spacing": ("hashSpacing", 20.0),
    "bond_spacing": ("bondSpacing", 1.0),        # percent
    "chain_angle": ("chainAngle", 1.0),          # degrees
    "label_size": ("labelSize", 20.0),
    "caption_size": ("captionSize", 20.0),
    "label_font": ("labelFont", None),
    "caption_font": ("captionFont", None),
    "magnification": ("magnification", 1.0),
    "fractional_widths": ("fractionalWidths", None),
    "rulers_visible": ("rulersVisible", None),
    "crosshair_visible": ("crosshairVisible", None),
}

COPY_COMMANDS = {
    "smiles": "copyAsSMILES", "inchi": "copyAsInChi", "inchikey": "copyAsInChIKey", "mol": "copyAsMol",
    "molv3000": "copyAsMolV3000", "cdxml": "copyAsCdxml", "helm": "copyAsHELM", "sln": "copyAsSLN",
    "helm_natural": "copyAsHELMUsingNaturalAnalog",
}
PASTE_COMMANDS = {
    "cdxml": "pasteAsCdxml", "smiles": "pasteAsSMILES", "inchi": "pasteAsInChi", "mol": "pasteAsMol",
    "name": "pasteAsStructureToName", "helm": "pasteAsHELM", "sequence": "pasteAsSequence", "peptide": "pasteAsPeptide",
    "dna": "pasteAsDNA", "rna": "pasteAsRNA", "text": "paste",
}

JS_PRELUDE = """
const BUNDLE = %(bundle)s;
const ARGS = %(args)s;
const app = Application(BUNDLE);
app.includeStandardAdditions = true;
function cmd(name) { const c = app.commands.whose({name: name}); if (!c.length) throw new Error("no ChemDraw command named " + name); return c[0]; }
function docByName(name) {
  if (name === null || name === undefined || name === "") {
    const ws = app.windows();
    for (const w of ws) { try { const d = w.document(); const n = d.name(); if (n && !n.startsWith(%(autosave)s)) return d; } catch (e) {} }
    throw new Error("no front document");
  }
  if (typeof name === "number" || /^#?\\d+$/.test(String(name))) {
    // a window id (stable even when ChemDraw renames an untitled document after its first edit)
    const wid = parseInt(String(name).replace("#", ""), 10);
    for (const w of app.windows()) { try { if (w.id() === wid) return w.document(); } catch (e) {} }
    throw new Error("no document window with id " + wid);
  }
  const names = app.documents.name();
  if (names.indexOf(name) < 0) {
    // tiers, most specific first: case-insensitive, then without ".cdxml", then a unique prefix -- so that "scheme" still
    // resolves to scheme.cdxml while scheme.cdx is open, and "untitled acs document 1996-1" is not spoiled by "...1996-10"
    const low = name.toLowerCase();
    let m = names.filter(n => n.toLowerCase() === low);
    if (!m.length) m = names.filter(n => n.toLowerCase() === low + ".cdxml");
    if (!m.length) m = names.filter(n => n.toLowerCase().startsWith(low));
    if (m.length === 1) name = m[0];
    else if (m.length) throw new Error("ambiguous document name " + JSON.stringify(name) + "; it matches: " + JSON.stringify(m));
    else throw new Error("no document named " + JSON.stringify(name) + "; open documents: " + JSON.stringify(names));
  }
  return app.documents.byName(name);
}
// Commands act on the front document window *inside* ChemDraw.  Reordering windows (index = 1, or AXRaise through
// System Events) works while ChemDraw stays in the background, so the app is never activated: activating would steal
// keyboard focus from the user, and stray keystrokes become ChemDraw hotkeys that relabel selected atoms.
let PREV = null;
function frontDocName() {
  try { for (const w of app.windows()) { try { const n = w.document().name(); if (n && !n.startsWith(%(autosave)s)) return n; } catch (e) {} } } catch (e) {}
  return null;
}
function front(d) {
  let name = null; try { name = d.name(); } catch (e) { return; }
  if (frontDocName() === name) return;
  try { d.window().index = 1; } catch (e) {}
  if (frontDocName() === name) return;
  try {
    const se = Application("System Events");
    const procs = se.processes.whose({bundleIdentifier: BUNDLE});
    if (procs.length) { procs[0].windows.byName(name).actions.byName("AXRaise").perform(); }
  } catch (e) {}
  if (frontDocName() === name) return;
  if (ARGS.__allow_activate) {   // last resort, only when the caller permits it
    try { const se = Application("System Events"); const fp = se.processes.whose({frontmost: true}); if (fp.length) PREV = fp[0].bundleIdentifier(); } catch (e) {}
    app.activate();
    for (let i = 0; i < 40 && !app.frontmost(); i++) delay(0.05);
    try { d.window().index = 1; } catch (e) {}
  }
}
function restoreFocus() {
  if (PREV && PREV !== BUNDLE && ARGS.__restore_focus !== false) { try { Application(PREV).activate(); } catch (e) {} }
}
function doCmd(name) { const c = cmd(name); if (!c.enabled()) return {ok: false, enabled: false}; app.do(c); return {ok: true, enabled: true}; }
function main() {
%(body)s
}
(function () { let out; try { out = JSON.stringify({ok: true, result: main()}); } catch (e) { out = JSON.stringify({ok: false, error: String(e)}); } restoreFocus(); return out; })();
"""

JS_UI = """
const ARGS = %(args)s;
const se = Application("System Events");
function main() {
  const procs = se.processes.whose({bundleIdentifier: %(bundle)s});
  const out = [];
  if (!procs.length) return out;
  const p = procs[0];
  const dismiss = ARGS.dismiss;
  const pref = ["OK", "Yes", "Don't Save", "Don’t Save", "Discard", "Cancel", "Close", "No"];
  function handle(el, where) {
    let text = [], buttons = [];
    try { text = el.staticTexts.value().filter(t => t); } catch (e) {}
    try { buttons = el.buttons.name(); } catch (e) {}
    const entry = {where: where, text: text, buttons: buttons};
    if (dismiss && buttons.length) {
      let chosen = null;
      for (const b of pref) { if (buttons.indexOf(b) >= 0) { chosen = b; break; } }
      if (!chosen) chosen = buttons[0];
      try { el.buttons.byName(chosen).click(); entry.clicked = chosen; } catch (e) { entry.clickError = String(e); }
    }
    out.push(entry);
  }
  for (const w of p.windows()) {
    let sub = ""; try { sub = w.subrole(); } catch (e) {}
    let name = ""; try { name = w.name(); } catch (e) {}
    try { for (const sh of w.sheets()) handle(sh, "sheet of " + name); } catch (e) {}
    if (sub === "AXDialog" || sub === "AXSystemDialog" || sub === "AXFloatingWindow" && !name) handle(w, "dialog " + name);
  }
  return out;
}
(function () { try { return JSON.stringify({ok: true, result: main()}); } catch (e) { return JSON.stringify({ok: false, error: String(e)}); } })();
"""

# Pasteboard helpers (ObjC bridge), prepended to the body of the scripts that paste into / copy out of ChemDraw.
# ``pbSnapshot`` records every item of a pasteboard with every flavour it carries (base64 data per type) and
# ``pbRestore`` writes them back, so an image, Finder files, rich text or a structure copied inside ChemDraw (binary
# CDX/CDXML flavours) on the user's clipboard survives an insert / selection round trip -- a pbpaste/pbcopy round trip
# would keep only the plain text.  Everything happens inside one script, so nothing can slip in between the snapshot,
# the paste and the restore.
JS_PASTEBOARD = """
ObjC.import("AppKit");
const PB_TEXT = "public.utf8-plain-text";
function pbSnapshot(pb) {
  const items = []; const arr = pb.pasteboardItems; const n = (arr && !arr.isNil()) ? Number(arr.count) : 0;
  for (let i = 0; i < n; i++) {
    const it = arr.objectAtIndex(i); const entry = {};
    for (const t of (ObjC.deepUnwrap(it.types) || [])) {
      try { const d = it.dataForType(t); if (d && !d.isNil()) entry[t] = ObjC.unwrap(d.base64EncodedStringWithOptions(0)); } catch (e) {}
    }
    items.push(entry);
  }
  return items;
}
function pbRestore(pb, items) {
  const objs = [];
  for (const entry of items) {
    const it = $.NSPasteboardItem.alloc.init; let any = false;
    for (const t in entry) { const d = $.NSData.alloc.initWithBase64EncodedStringOptions(entry[t], 0); if (d && !d.isNil() && it.setDataForType(d, t)) any = true; }
    if (any) objs.push(it);
  }
  pb.clearContents;
  if (objs.length) pb.writeObjects($(objs));
}
function pbSetText(pb, text) { pb.clearContents; if (!pb.setStringForType(text, PB_TEXT)) throw new Error("could not write to the pasteboard"); }
function pbGetText(pb) {
  // like pbpaste: the plain-text flavour (in any of its spellings)
  for (const t of [PB_TEXT, "public.utf16-external-plain-text", "public.plain-text", "com.apple.traditional-mac-plain-text"]) {
    try { const s = pb.stringForType(t); if (s && !s.isNil()) return ObjC.unwrap(s); } catch (e) {}
  }
  return "";
}
"""


def _plist(app_path: str) -> dict:
    try:
        with open(os.path.join(app_path, "Contents", "Info.plist"), "rb") as f:
            return plistlib.load(f)
    except Exception:
        return {}


def find_app() -> Optional[str]:
    try:
        from .. import discovery
        p = discovery.find_chemdraw()
        if p:
            return p
    except Exception:
        pass
    for env in ("CHEMDRAW_APP", "CHEMDRAW_EXE"):
        p = os.environ.get(env)
        if p and os.path.exists(p):
            return p
    cands = [p for p in glob.glob("/Applications/ChemDraw*.app") + glob.glob(os.path.expanduser("~/Applications/ChemDraw*.app"))
             if "Collections" not in p]
    return sorted(cands)[-1] if cands else None


class MacBackend(ChemDrawBackend):
    platform = "macos"

    def __init__(self, app_path: Optional[str] = None, timeout: float = 60.0):
        self.app_path = app_path or find_app()
        self.bundle_id = (_plist(self.app_path).get("CFBundleIdentifier") if self.app_path else None) or DEFAULT_BUNDLE_ID
        self.timeout = timeout
        self.last_alerts: List[str] = []

    # -- plumbing ---------------------------------------------------------------
    def executable(self) -> Optional[str]:
        return self.app_path

    def _run_js(self, script: str, timeout: Optional[float] = None) -> Any:
        """Run a JXA script.  While it runs, ChemDraw is polled for modal alerts (which block Apple events):
        an alert is dismissed within a couple of seconds and its text ends up in the error / ``last_alerts``."""
        limit = timeout or self.timeout
        proc = subprocess.Popen(["osascript", "-l", "JavaScript", "-"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        alerts: List[str] = []
        t0 = time.time()
        next_check = t0 + 1.5
        import threading
        result: Dict[str, Any] = {}

        def pump():
            try:
                result["out"], result["err"] = proc.communicate(input=script)
            except Exception as e:  # pragma: no cover
                result["exc"] = e

        th = threading.Thread(target=pump, daemon=True)
        th.start()
        while th.is_alive():
            th.join(0.2)
            now = time.time()
            if th.is_alive() and now >= next_check:
                next_check = now + 2.0
                found = self.dismiss_dialogs()
                if found:
                    alerts.extend(found)
            if th.is_alive() and now - t0 > limit:
                proc.kill()
                th.join(2)
                found = self.dismiss_dialogs()
                alerts.extend(found)
                msg = "ChemDraw did not answer within %.0fs" % limit
                if alerts:
                    msg += "; it showed an alert which was dismissed: " + " | ".join(alerts)
                self.last_alerts = alerts
                raise ChemDrawError(msg)
        self.last_alerts = alerts
        if os.environ.get("CHEMDRAWMCP_DEBUG"):
            sys.stderr.write("[chemdrawmcp] jxa %.2fs polls=%d alerts=%d pid=%d tid=%s\n" % (
                time.time() - t0, int((time.time() - t0) // 2), len(alerts), os.getpid(), threading.get_ident()))
            sys.stderr.flush()
        out = (result.get("out") or "").strip()
        stderr = result.get("err") or ""
        if proc.returncode != 0 and not out:
            raise ChemDrawError("osascript failed: %s%s" % (stderr.strip(), (" (alert: %s)" % " | ".join(alerts)) if alerts else ""))
        try:
            payload = json.loads(out.splitlines()[-1])
        except Exception:
            raise ChemDrawError("unexpected osascript output: %s %s" % (out[:500], stderr[:300]))
        if not payload.get("ok"):
            err = payload.get("error", "unknown error")
            if not alerts and ("timed out" in err or "not available" in err or "handler failed" in err):
                alerts = self.dismiss_dialogs()
            if alerts:
                err += " (ChemDraw alert: %s)" % " | ".join(alerts)
            raise ChemDrawError(_clean_error(err))
        return payload.get("result")

    def js(self, body: str, timeout: Optional[float] = None, **args) -> Any:
        """Run JavaScript in the ChemDraw automation context.  ``body`` is the body of ``main()``; use ``return``.
        Available: app, ARGS (your kwargs), cmd(name), docByName(name), front(doc), doCmd(name)."""
        args.setdefault("__restore_focus", os.environ.get("CHEMDRAWMCP_RESTORE_FOCUS", "1") != "0")
        args.setdefault("__allow_activate", os.environ.get("CHEMDRAWMCP_ALLOW_ACTIVATE", "0") == "1")
        script = JS_PRELUDE % {"bundle": json.dumps(self.bundle_id), "args": json.dumps(args), "body": body,
                               "autosave": json.dumps(AUTOSAVE_PREFIX)}
        return self._run_js(script, timeout)

    def _ui(self, dismiss: bool) -> List[dict]:
        script = JS_UI % {"bundle": json.dumps(self.bundle_id), "args": json.dumps({"dismiss": dismiss})}
        t0 = time.time()
        try:
            r = subprocess.run(["osascript", "-l", "JavaScript", "-"], input=script, capture_output=True, text=True, timeout=20)
            payload = json.loads((r.stdout or "").strip().splitlines()[-1])
            if os.environ.get("CHEMDRAWMCP_DEBUG"):
                sys.stderr.write("[chemdrawmcp] ui(dismiss=%s) %.2fs rc=%s ok=%s\n" % (dismiss, time.time() - t0, r.returncode, payload.get("ok")))
            if payload.get("ok"):
                return payload["result"] or []
        except Exception as e:
            if os.environ.get("CHEMDRAWMCP_DEBUG"):
                sys.stderr.write("[chemdrawmcp] ui(dismiss=%s) failed after %.2fs: %s\n" % (dismiss, time.time() - t0, e))
        return []

    # -- lifecycle --------------------------------------------------------------
    def is_running(self) -> bool:
        try:
            r = subprocess.run(["osascript", "-l", "JavaScript", "-e", "Application(%s).running()" % json.dumps(self.bundle_id)],
                               capture_output=True, text=True, timeout=15)
            return (r.stdout or "").strip() == "true"
        except Exception:
            return False

    def launch(self, wait_seconds: float = 90.0) -> dict:
        if not self.app_path:
            raise ChemDrawNotAvailable("ChemDraw.app not found in /Applications (set CHEMDRAW_APP)")
        was_running = self.is_running()
        if not was_running:
            # -g: do not bring ChemDraw to the foreground, -j: launch hidden (its windows still open behind other apps)
            subprocess.run(["open", "-g", "-j", "-a", self.app_path], check=False)
        t0 = time.time()
        last = None
        while time.time() - t0 < wait_seconds:
            try:
                v = self.js("return app.version();", timeout=20)
                return {"running": True, "was_running": was_running, "version": v, "app": self.app_path,
                        "seconds": round(time.time() - t0, 1)}
            except ChemDrawError as e:
                last = e
                time.sleep(1.0)
        raise ChemDrawError("ChemDraw did not start within %.0fs (%s)" % (wait_seconds, last))

    def quit(self, save: bool = False) -> None:
        self.js("app.quit({saving: ARGS.save ? 'yes' : 'no'}); return true;", save=save)

    def version(self) -> str:
        return str(self.js("return app.version();", timeout=20))

    def activate(self) -> None:
        self.js("app.activate(); return true;")

    # -- documents ---------------------------------------------------------------
    def documents(self) -> List[DocumentInfo]:
        rows = self.js("""
            const out = [];
            let frontName = null, frontWid = null;
            try { for (const w of app.windows()) { try { const n = w.document().name(); if (n && !n.startsWith(ARGS.autosave)) { frontName = n; try { frontWid = w.id(); } catch (e) {} break; } } catch (e) {} } } catch (e) {}
            for (const d of app.documents()) {
              const n = d.name();
              if (n.startsWith(ARGS.autosave)) continue;
              const e = {name: n};
              try { e.modified = d.modified(); } catch (err) {}
              try { const f = d.file(); e.path = f ? f.toString() : null; } catch (err) { e.path = null; }
              try { e.atoms = d.atoms.length; e.bonds = d.bonds.length; e.molecules = d.molecules.length; e.captions = d.captions.length; e.arrows = d.arrows.length; } catch (err) {}
              try { e.window_id = d.window().id(); } catch (err) {}
              // by window id, not by name: two open files with the same basename share a document name
              e.front = (frontWid !== null && e.window_id !== undefined) ? e.window_id === frontWid : (n === frontName);
              out.push(e);
            }
            out.sort((a, b) => (b.front ? 1 : 0) - (a.front ? 1 : 0));
            return out;""", autosave=AUTOSAVE_PREFIX)
        return [DocumentInfo(name=r["name"], path=r.get("path"), modified=bool(r.get("modified")), front=bool(r.get("front")),
                             atoms=int(r.get("atoms", 0) or 0), bonds=int(r.get("bonds", 0) or 0),
                             molecules=int(r.get("molecules", 0) or 0), captions=int(r.get("captions", 0) or 0),
                             arrows=int(r.get("arrows", 0) or 0), window_id=r.get("window_id")) for r in rows]

    def snapshot(self) -> Dict[str, Any]:
        """version + documents + front document in one round trip (used by the status tool)."""
        rows = self.js("""
            const out = {version: app.version(), frontmost: app.frontmost(), documents: []};
            let frontName = null, frontWid = null;
            try { for (const w of app.windows()) { try { const n = w.document().name(); if (n && !n.startsWith(ARGS.autosave)) { frontName = n; try { frontWid = w.id(); } catch (e) {} break; } } catch (e) {} } } catch (e) {}
            for (const d of app.documents()) {
              const n = d.name();
              if (n.startsWith(ARGS.autosave)) continue;
              const e = {name: n};
              try { e.modified = d.modified(); } catch (err) {}
              try { const f = d.file(); e.path = f ? f.toString() : null; } catch (err) { e.path = null; }
              try { e.atoms = d.atoms.length; e.molecules = d.molecules.length; e.captions = d.captions.length; e.arrows = d.arrows.length; } catch (err) {}
              try { e.window_id = d.window().id(); } catch (err) {}
              e.front = (frontWid !== null && e.window_id !== undefined) ? e.window_id === frontWid : (n === frontName);
              out.documents.push(e);
            }
            out.documents.sort((a, b) => (b.front ? 1 : 0) - (a.front ? 1 : 0));
            out.front_document = frontName;
            return out;""", autosave=AUTOSAVE_PREFIX, timeout=45)
        return rows

    def front_document(self) -> Optional[str]:
        try:
            return self.js("return docByName(null).name();")
        except ChemDrawError:
            return None

    def activate_document(self, name: Optional[str]) -> str:
        return self.js("const d = docByName(ARGS.name); front(d); return d.name();", name=name)

    def stationery_path(self, style: str) -> Optional[str]:
        """Path of a stationery file: an absolute path, or a style/stationery name inside the app bundle."""
        if os.path.isabs(style) and os.path.exists(style):
            return style
        if not self.app_path:
            return None
        folder = os.path.join(self.app_path, "Contents", "Resources", "SpecialPurpose", "Stationery")
        try:
            from ..styles import get_style
            title = get_style(style).get("stationery")
            p = os.path.join(folder, title)
            if os.path.exists(p):
                return p
        except KeyError:
            pass
        low = style.lower().replace(".cds", "")
        for p in sorted(glob.glob(os.path.join(folder, "*.cds"))):
            base = os.path.basename(p)[:-4].lower()
            if base == low or base.startswith(low):
                return p
        return None

    def new_document(self, stationery: Optional[str] = None) -> str:
        if stationery:
            path = self.stationery_path(stationery)
            if not path:
                raise ChemDrawError("stationery %r not found (use a style name like 'acs', 'nature', 'rsc', 'wiley' or a .cds path)" % stationery)
            return self.js("""
                const before = app.documents.name();
                app.open(Path(ARGS.path));
                for (let i = 0; i < 40; i++) {
                  const after = app.documents.name().filter(n => before.indexOf(n) < 0 && !n.startsWith(ARGS.autosave));
                  if (after.length) { const d = app.documents.byName(after[0]); front(d); return after[0]; }
                  delay(0.1);
                }
                throw new Error("no new document appeared after opening stationery");""",
                           path=path, autosave=AUTOSAVE_PREFIX, timeout=60)
        return self.js("""
            const before = app.documents.name();
            try { const d = app.Document(); app.documents.push(d); } catch (e) { app.do(cmd("newDocument")); }
            for (let i = 0; i < 40; i++) {
              const after = app.documents.name().filter(n => before.indexOf(n) < 0 && !n.startsWith(ARGS.autosave));
              if (after.length) { const d = app.documents.byName(after[0]); front(d); return after[0]; }
              delay(0.1);
            }
            throw new Error("no new document appeared");""", autosave=AUTOSAVE_PREFIX)

    def open_document(self, path: str) -> str:
        path = os.path.abspath(os.path.expanduser(path))
        if not os.path.exists(path):
            raise ChemDrawError("file not found: %s" % path)
        base = os.path.basename(path)
        real = os.path.realpath(path)
        # spellings under which ChemDraw may report the file (symlinks resolved, with or without /private)
        paths = list(dict.fromkeys([path, real] + [p[len("/private"):] for p in (path, real) if p.startswith("/private/")]))
        return self.js("""
            // The opened document is identified by its window (a window id that did not exist before) or by its file,
            // never by its name alone: a second file with the basename of an open document gets the same document name.
            const widsBefore = [];
            try { for (const w of app.windows()) { try { widsBefore.push(w.id()); } catch (e) {} } } catch (e) {}
            const already = app.documents.name().indexOf(ARGS.base) >= 0;
            function byPath() {
              for (const d of app.documents()) { try { const f = d.file(); if (f && ARGS.paths.indexOf(f.toString()) >= 0) return d; } catch (e) {} }
              return null;
            }
            app.open(Path(ARGS.path));
            for (let i = 0; i < 50; i++) {
              for (const w of app.windows()) {
                try {
                  const wid = w.id();
                  if (widsBefore.indexOf(wid) >= 0) continue;
                  const d = w.document(); const n = d.name();
                  if (n && !n.startsWith(ARGS.autosave)) { front(d); return n; }
                } catch (e) {}
              }
              const d = byPath();   // the same file was already open: ChemDraw just brings its window forward
              if (d) { front(d); return d.name(); }
              if (already && i >= 10) { const d2 = app.documents.byName(ARGS.base); front(d2); return ARGS.base; }
              delay(0.1);
            }
            throw new Error("ChemDraw did not open " + ARGS.path);""", path=path, paths=paths, base=base, autosave=AUTOSAVE_PREFIX, timeout=90)

    def close_document(self, name: Optional[str] = None, save: bool = False) -> None:
        self.js("const d = docByName(ARGS.name); d.close({saving: ARGS.save ? 'yes' : 'no'}); return true;", name=name, save=save)

    def save_document(self, name: Optional[str], path: str, fmt: str) -> str:
        fmt = normalize_format(fmt)
        spec = EXPORT_FORMATS[fmt]
        exts = {fmt, spec["ext"].lstrip("."), {"jpeg": "jpg", "tiff": "tif"}.get(fmt, fmt)}   # extensions of this format
        path = os.path.abspath(os.path.expanduser(path))
        root, ext = os.path.splitext(path)
        if ext.lower().lstrip(".") not in exts:
            path = path + spec["ext"]
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        if os.path.exists(path):
            os.remove(path)
        t_save = time.time()
        self.js("const d = docByName(ARGS.name); app.save(d, {in: Path(ARGS.path), as: ARGS.type}); return true;",
                name=name, path=path, type=spec["name"], timeout=120)
        if not os.path.exists(path):
            # ChemDraw may append its own extension (x.png -> x.png.png) or spell it its way (x.jpg -> x.jpeg).  Accept
            # only a file of the requested format written by this save: an older sibling of another format (scheme.cdxml
            # next to the requested scheme.svg) must never be reported as the exported file.  glob.escape keeps a
            # directory such as "run[1]" from being read as a character class.
            cands = glob.glob(glob.escape(path) + ".*") + glob.glob(glob.escape(root) + ".*")
            cands = [c for c in dict.fromkeys(cands)
                     if os.path.isfile(c) and os.path.splitext(c)[1].lower().lstrip(".") in exts
                     and os.path.getmtime(c) >= int(t_save)]
            if cands:
                return sorted(cands, key=os.path.getmtime)[-1]
            alerts = self.dismiss_dialogs()
            raise ChemDrawError("ChemDraw did not write %s%s" % (path, (" (alert: %s)" % " | ".join(alerts)) if alerts else ""))
        return path

    def document_cdxml(self, name: Optional[str] = None) -> str:
        tmpdir = tempfile.mkdtemp(prefix="chemdrawmcp_")
        try:
            p = self.save_document(name, os.path.join(tmpdir, "document.cdxml"), "cdxml")
            with open(p, encoding="utf-8", errors="replace") as f:
                return f.read().replace("\r\n", "\n")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def document_settings(self, name: Optional[str] = None) -> Dict[str, Any]:
        raw = self.js("""
            const d = docByName(ARGS.name); const out = {name: d.name()};
            for (const k of ARGS.keys) { try { out[k] = d[k](); } catch (e) { out[k] = null; } }
            try { out.selectionBounds = d.selection.bounds(); } catch (e) {}
            return out;""", name=name, keys=[v[0] for v in SETTINGS.values()])
        out: Dict[str, Any] = {"name": raw.get("name")}
        for key, (prop, factor) in SETTINGS.items():
            v = raw.get(prop)
            if v is None:
                continue
            out[key] = (v / factor) if factor and isinstance(v, (int, float)) else v
        return out

    def set_document_settings(self, name: Optional[str], settings: Dict[str, Any]) -> Dict[str, Any]:
        assignments = {}
        for key, val in settings.items():
            if key not in SETTINGS:
                raise ChemDrawError("unknown setting %r; known: %s" % (key, ", ".join(SETTINGS)))
            prop, factor = SETTINGS[key]
            assignments[prop] = int(round(val * factor)) if factor and isinstance(val, (int, float)) else val
        self.js("""
            const d = docByName(ARGS.name);
            for (const k in ARGS.set) { d[k] = ARGS.set[k]; }
            return true;""", name=name, set=assignments)
        return self.document_settings(name)

    # -- content -------------------------------------------------------------------
    def insert(self, name: Optional[str], data: str, kind: str = "cdxml") -> dict:
        """Paste ``data`` into the document through the general pasteboard.  The user's clipboard (every item and
        flavour) is snapshotted and put back inside the same script, whether or not the paste succeeds."""
        kind = kind.lower()
        if kind not in PASTE_COMMANDS:
            raise ChemDrawError("kind must be one of %s" % ", ".join(INSERT_KINDS))
        res = self.js(JS_PASTEBOARD + """
            const d = docByName(ARGS.name); front(d);
            const before = {atoms: d.atoms.length, molecules: d.molecules.length, captions: d.captions.length};
            const pb = $.NSPasteboard.generalPasteboard; const saved = pbSnapshot(pb);
            let r;
            let after = before;
            try {
              pbSetText(pb, ARGS.data); r = doCmd(ARGS.command);
              // under load ChemDraw registers the pasted objects a moment after the command returns: wait for a change
              for (let i = 0; i < 40 && r.enabled; i++) {
                after = {atoms: d.atoms.length, molecules: d.molecules.length, captions: d.captions.length};
                if (after.atoms !== before.atoms || after.captions !== before.captions || after.molecules !== before.molecules) break;
                delay(0.05);
              }
            } finally { pbRestore(pb, saved); }
            if (!r.enabled) throw new Error("ChemDraw command " + ARGS.command + " is not enabled (is the clipboard content valid " + ARGS.kind + "?)");
            let sel = null; try { sel = {smiles: d.selection.smiles(), bounds: d.selection.bounds()}; } catch (e) {}
            return {document: d.name(), before: before, after: after, selection: sel};""",
                      name=name, data=data, command=PASTE_COMMANDS[kind], kind=kind, timeout=90)
        res["atoms_added"] = res["after"]["atoms"] - res["before"]["atoms"]
        res["kind"] = kind
        if res["atoms_added"] <= 0 and kind != "text":
            alerts = self.dismiss_dialogs()
            res["warning"] = "no atoms were added" + ((": " + " | ".join(alerts)) if alerts else "")
        return res

    def select_all(self, name: Optional[str] = None) -> None:
        self.js("const d = docByName(ARGS.name); front(d); doCmd('selectAll'); return true;", name=name)

    def deselect(self, name: Optional[str] = None) -> None:
        self.js("""const d = docByName(ARGS.name); front(d);
                 const r = doCmd('selectAll'); if (cmd('invertSelection').enabled()) app.do(cmd('invertSelection')); return true;""", name=name)

    def selection_info(self, name: Optional[str] = None) -> Dict[str, Any]:
        return self.js("""
            const d = docByName(ARGS.name); const s = d.selection; const out = {document: d.name()};
            const get = (k, f) => { try { out[k] = f(); } catch (e) { out[k] = null; } };
            get('smiles', () => s.smiles()); get('formula', () => s.molecularFormula()); get('molecular_weight', () => s.molecularWeight());
            get('exact_mass', () => s.exactMass()); get('elemental_analysis', () => s.elementalAnalysis()); get('text', () => s.characters());
            get('bounds', () => s.bounds()); get('atoms', () => s.atoms.length); get('bonds', () => s.bonds.length);
            get('molecules', () => s.molecules.length); get('captions', () => s.captions.length); get('arrows', () => s.arrows.length);
            out.empty = !(out.atoms || out.captions || out.arrows);
            return out;""", name=name)

    def selection_cdxml(self, name: Optional[str] = None) -> str:
        return self.selection_as(name, "cdxml")

    def selection_as(self, name: Optional[str], what: str) -> str:
        what = what.lower()
        if what not in COPY_COMMANDS:
            raise ChemDrawError("what must be one of %s" % ", ".join(COPY_COMMANDS))
        out = self.js(JS_PASTEBOARD + """
            const d = docByName(ARGS.name); front(d);
            const pb = $.NSPasteboard.generalPasteboard; const saved = pbSnapshot(pb);
            let out = "";
            try {
              pb.clearContents;
              const count0 = Number(pb.changeCount);
              const r = doCmd(ARGS.command);
              if (!r.enabled) throw new Error(ARGS.command + " is not enabled: select one or more structures first");
              // ChemDraw fills the pasteboard asynchronously: wait until its change count advances (up to ~3 s)
              for (let i = 0; i < 60 && Number(pb.changeCount) <= count0; i++) delay(0.05);
              delay(0.05);
              out = pbGetText(pb);
            } finally { pbRestore(pb, saved); }
            return out;""", name=name, command=COPY_COMMANDS[what])
        return str(out or "").replace("\r\n", "\n")

    # -- commands / menus ---------------------------------------------------------
    def list_commands(self, pattern: Optional[str] = None, enabled_only: bool = False,
                      with_descriptions: bool = True) -> List[CommandInfo]:
        rows = self.js("""
            const names = app.commands.name();
            const pat = ARGS.pattern ? new RegExp(ARGS.pattern, 'i') : null;
            let idx = []; names.forEach((n, i) => { if (!pat || pat.test(n)) idx.push(i); });
            if (idx.length > 400 && (ARGS.desc || ARGS.enabled_only)) idx = idx.slice(0, 400);
            const out = [];
            if (ARGS.desc || ARGS.enabled_only) {
              const en = app.commands.enabled(); const sd = ARGS.desc ? app.commands.shortDescription() : null; const ld = ARGS.desc ? app.commands.longDescription() : null;
              for (const i of idx) { if (ARGS.enabled_only && !en[i]) continue; out.push({name: names[i], enabled: en[i], sd: sd ? sd[i] : '', ld: ld ? ld[i] : ''}); }
            } else { for (const i of idx) out.push({name: names[i]}); }
            return out;""", pattern=pattern, desc=with_descriptions, enabled_only=enabled_only, timeout=120)
        return [CommandInfo(r["name"], r.get("enabled"), r.get("sd", "") or "", r.get("ld", "") or "") for r in rows]

    def run_command(self, command: str, name: Optional[str] = None) -> dict:
        res = self.js("""
            const d = docByName(ARGS.name); front(d);
            const c = cmd(ARGS.command); const en = c.enabled();
            if (!en) return {ok: false, command: ARGS.command, enabled_before: false, document: d.name(),
                             error: "command is not enabled right now (" + c.shortDescription() + "; " + c.longDescription() + ")"};
            const before = {atoms: d.atoms.length, molecules: d.molecules.length, captions: d.captions.length, arrows: d.arrows.length};
            app.do(c);
            const after = {atoms: d.atoms.length, molecules: d.molecules.length, captions: d.captions.length, arrows: d.arrows.length};
            return {ok: true, command: ARGS.command, enabled_before: true, document: d.name(), before: before, after: after,
                    description: c.shortDescription()};""", name=name, command=command, timeout=120)
        res["alerts"] = self.dismiss_dialogs()
        return res

    def menus(self) -> Dict[str, List[str]]:
        # fetch names first: evaluating the menu objects themselves ("app.menus()") fails with
        # "AppleEvent handler failed" while byName()/name() work reliably
        return self.js("""
            const out = {};
            const names = app.menus.name();
            for (const n of names) {
              if (!n) continue;
              try { out[n] = app.menus.byName(n).menuItems.name().filter(x => x); } catch (e) { out[n] = []; }
            }
            return out;""")

    def do_menu(self, menu: str, item: str, name: Optional[str] = None) -> dict:
        res = self.js("""
            if (ARGS.name) { const d = docByName(ARGS.name); front(d); }
            const m = app.menus.byName(ARGS.menu); const it = m.menuItems.byName(ARGS.item);
            const en = it.enabled(); if (!en) return {ok: false, error: "menu item is disabled"};
            app.do(it); return {ok: true, menu: ARGS.menu, item: ARGS.item};""", name=name, menu=menu, item=item, timeout=120)
        res["alerts"] = self.dismiss_dialogs()
        return res

    # -- pictures ---------------------------------------------------------------------
    def screenshot(self, name: Optional[str] = None) -> bytes:
        wid = self.js("const d = docByName(ARGS.name); front(d); return d.window().id();", name=name)
        fd, path = tempfile.mkstemp(suffix=".png", prefix="chemdrawmcp_shot_")
        os.close(fd)
        try:
            r = subprocess.run(["screencapture", "-l", str(wid), "-x", "-o", path], capture_output=True, text=True, timeout=30)
            if r.returncode != 0 or not os.path.getsize(path):
                raise ChemDrawError("screencapture failed (grant Screen Recording permission to your terminal): %s" % (r.stderr or "").strip())
            with open(path, "rb") as f:
                return f.read()
        finally:
            try:
                os.remove(path)
            except OSError:
                pass

    # -- dialogs -----------------------------------------------------------------------
    def dismiss_dialogs(self) -> List[str]:
        entries = self._ui(dismiss=True)
        texts = []
        for e in entries:
            t = " ".join(str(x) for x in e.get("text", []) if x).strip()
            if e.get("clicked"):
                t += " [dismissed with '%s']" % e["clicked"]
            texts.append(t or "(dialog without text)")
        self.last_alerts = texts
        return texts

    def pending_dialogs(self) -> List[dict]:
        return self._ui(dismiss=False)

    def capabilities(self) -> Dict[str, Any]:
        info = super().capabilities()
        info.update({"bundle_id": self.bundle_id, "version": _plist(self.app_path).get("CFBundleShortVersionString") if self.app_path else None,
                     "running": self.is_running()})
        return info


def _clean_error(err: str) -> str:
    err = re.sub(r"^Error:\s*", "", err.strip())
    return err
