"""Tests for the "Agent Chat" ChemDraw add-in package (chemdrawmcp/addin) and chemdrawmcp.addin_builder.

Nothing here talks to ChemDraw or the bridge: the add-in files are checked statically (metadata, balanced
braces/tags, no external URLs, `node --check` when node is available) and the builder runs in temp dirs.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from html.parser import HTMLParser
from urllib.parse import urlsplit

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from chemdrawmcp import addin_builder as ab  # noqa: E402

ADDIN = ab.ADDIN_SRC
WEB_FILES = ["main.html", "app.js", "style.css", "config.example.js"]
VOID_TAGS = {"area", "base", "br", "col", "embed", "hr", "img", "input", "link", "meta", "param", "source", "track", "wbr"}


def read(name, folder=ADDIN):
    with open(os.path.join(folder, name), encoding="utf-8") as f:
        return f.read()


# ---------------------------------------------------------------------------
# cheap syntax checkers
# ---------------------------------------------------------------------------

_OPEN = {"(": ")", "[": "]", "{": "}"}
_CLOSE = {")": "(", "]": "[", "}": "{"}
_REGEX_PREV = set("(,=:[!&|?{};+-*%<>~^")
_REGEX_WORDS = {"return", "typeof", "case", "in", "of", "new", "delete", "void", "throw", "instanceof", "else", "do"}


def js_balance(src):
    """None when (), [] and {} balance outside strings/comments/regex literals, else a message."""
    stack = []
    i, n, line = 0, len(src), 1
    last_sig = ""      # last significant character
    last_word = ""     # identifier/keyword ending at last_sig
    while i < n:
        c = src[i]
        if c == "\n":
            line += 1
        if src.startswith("//", i):
            j = src.find("\n", i)
            i = n if j < 0 else j
            continue
        if src.startswith("/*", i):
            j = src.find("*/", i + 2)
            if j < 0:
                return "unterminated block comment at line %d" % line
            line += src.count("\n", i, j)
            i = j + 2
            continue
        if c in "'\"`":
            j = i + 1
            while j < n and src[j] != c:
                if src[j] == "\\":
                    j += 1
                elif src[j] == "\n" and c != "`":
                    return "unterminated string at line %d" % line
                j += 1
            if j >= n:
                return "unterminated string at line %d" % line
            line += src.count("\n", i, j)
            i = j + 1
            last_sig, last_word = c, ""
            continue
        if c == "/" and (last_sig == "" or last_sig in _REGEX_PREV or last_word in _REGEX_WORDS):
            j, in_class = i + 1, False
            while j < n:
                ch = src[j]
                if ch == "\\":
                    j += 2
                    continue
                if ch == "\n":
                    return "unterminated regex at line %d" % line
                if in_class:
                    if ch == "]":
                        in_class = False
                elif ch == "[":
                    in_class = True
                elif ch == "/":
                    break
                j += 1
            if j >= n:
                return "unterminated regex at line %d" % line
            j += 1
            while j < n and src[j].isalpha():
                j += 1
            i = j
            last_sig, last_word = "/", ""
            continue
        if c in _OPEN:
            stack.append((c, line))
        elif c in _CLOSE:
            if not stack or stack[-1][0] != _CLOSE[c]:
                return "unexpected %r at line %d" % (c, line)
            stack.pop()
        if c.isalnum() or c in "_$":
            last_word = (last_word + c) if (last_sig and (last_sig.isalnum() or last_sig in "_$")) else c
        elif not c.isspace():
            last_word = ""
        if not c.isspace():
            last_sig = c
        i += 1
    if stack:
        return "unclosed %r from line %d" % stack[-1]
    return None


class TagChecker(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.errors = []
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        for k, v in attrs:
            if k == "id":
                self.ids.add(v)
        if tag not in VOID_TAGS:
            self.stack.append((tag, self.getpos()[0]))

    def handle_startendtag(self, tag, attrs):
        for k, v in attrs:
            if k == "id":
                self.ids.add(v)

    def handle_endtag(self, tag):
        if tag in VOID_TAGS:
            return
        if not self.stack or self.stack[-1][0] != tag:
            self.errors.append("unexpected </%s> at line %d (open: %s)" % (tag, self.getpos()[0], [t for t, _ in self.stack[-3:]]))
            return
        self.stack.pop()


# ---------------------------------------------------------------------------
# the package
# ---------------------------------------------------------------------------

def test_metadata_json_valid_and_complete():
    meta = json.loads(read(ab.METADATA_FILE))
    for key in ("name", "description", "version", "menuItemText", "minimumAPIVersion", "url", "isModalDialog", "canBeUninstalled"):
        assert key in meta, key
    assert meta["url"] == "main.html"
    assert meta["minimumAPIVersion"] == "1.0"
    assert meta["isModalDialog"] is False
    assert meta["canBeUninstalled"] is True
    assert meta["name"] == ab.DEFAULT_NAME == meta["menuItemText"]
    assert re.match(r"^\d+\.\d+", meta["version"])


def test_package_files_present():
    for fname in ab.PACKAGE_FILES:
        assert os.path.isfile(os.path.join(ADDIN, fname)), fname
    assert not os.path.exists(os.path.join(ADDIN, "config.js")), "config.js is generated by the installer, never shipped"


def test_js_balance_checker_itself():
    assert js_balance("var a = {b: [1, (2)]}; var r = /[)]/g; s = 'x)'; t = `y}`; // }\n/* ( */") is None
    assert js_balance("function f() { return [1, 2; }") is not None
    assert js_balance("if (x) { y = (1 + 2; }") is not None
    assert js_balance("var s = 'open") is not None


def test_app_js_and_config_example_balanced():
    for fname in ("app.js", "config.example.js"):
        assert js_balance(read(fname)) is None, fname


def test_app_js_node_check():
    node = shutil.which("node")
    if not node:
        pytest.skip("node is not on PATH")
    for fname in ("app.js", "config.example.js"):
        r = subprocess.run([node, "--check", os.path.join(ADDIN, fname)], capture_output=True, text=True, timeout=60)
        assert r.returncode == 0, r.stderr


def test_main_html_tags_balanced_and_wired():
    html = read("main.html")
    assert html.lstrip().lower().startswith("<!doctype html>")
    assert '<meta charset="utf-8">' in html
    p = TagChecker()
    p.feed(html)
    p.close()
    assert not p.errors, p.errors
    assert not p.stack, "unclosed tags: %s" % p.stack
    for needed in ("prompt", "transcript", "status", "model-btn", "bridge-dot", "bridge-text", "suggest", "chooser", "sel-chip"):
        assert needed in p.ids, needed
    assert re.search(r'<script src="config\.js"></script>', html)
    assert re.search(r'<script src="app\.js"></script>', html)
    assert re.search(r'<link rel="stylesheet" href="style\.css">', html)
    assert html.index('src="config.js"') < html.index('src="app.js"'), "config.js must load before app.js"


def test_no_external_urls_in_web_files():
    for fname in WEB_FILES:
        text = read(fname)
        for m in re.finditer(r"https?://[^\s'\"<>)]+", text):
            host = urlsplit(m.group(0)).hostname or ""
            assert host in ("127.0.0.1", "localhost"), "%s references %s" % (fname, m.group(0))
        assert "<script src=\"http" not in text and "@import" not in text


def test_app_js_contract_strings():
    js = read("app.js")
    assert "47325" in js                           # default port when config.js is missing
    assert "CHEMDRAWMCP_BRIDGE" in js
    assert "X-Bridge-Token" in js
    assert "install.py --addin" in js              # the "bridge not configured" instruction
    assert "chemdraw-agent-bridge" in js           # how to start the bridge
    for route in ("/health", "/run", "/events/", "/cancel/", "/settings", "/clear", "/sessions"):
        assert route in js, route
    for kind in ("'session'", "'delta'", "'text'", "'tool'", "'tool_result'", "'stats'", "'done'", "'end'"):
        assert kind in js, kind
    for cmd in ("/model", "/clear", "/claude", "/codex", "/stop", "/help", "/skills"):
        assert "'%s'" % cmd in js, cmd
    for model in ("claude-opus-5-5", "claude-sonnet-5-5", "claude-haiku-5-5", "claude-fable-5-1", "gpt-6-sol", "gpt-6"):
        assert model in js, model


def test_app_js_feature_detects_chemdraw_api():
    js = read("app.js")
    assert "typeof ChemDrawAPI !== 'undefined'" in js
    assert "setDefaultSize" in js and "resizeTo" in js
    # every ChemDrawAPI.window call sits behind a typeof check
    for method in ("setDefaultSize", "resizeTo", "close", "onClose"):
        assert re.search(r"typeof w\.%s === 'function'" % method, js), method
    assert "window.localStorage" in js
    assert "EventSource" in js and "new window.EventSource" in js
    # no frameworks / bundler leftovers
    for bad in ("require(", "import ", "export ", "React", "Vue", "jQuery", "$("):
        assert bad not in js.replace("function $(sel)", "").replace("$('#", "").replace("\"_$\"", ""), bad


def test_config_example_defines_bridge_object():
    cfg = read("config.example.js")
    assert "window.CHEMDRAWMCP_BRIDGE" in cfg
    assert "port: 47325" in cfg
    assert "token:" in cfg


def test_readme_documents_install_and_bridge():
    readme = read("README.md")
    assert "install.py --addin" in readme
    assert "config.js" in readme
    assert "chemdraw-agent-bridge" in readme
    assert ".chemdrawaddin" in readme


# ---------------------------------------------------------------------------
# the builder
# ---------------------------------------------------------------------------

def test_build_writes_config_js(tmp_path):
    folder = ab.build_addin(str(tmp_path), 47999, "s3cret-token", name="Agent Chat")
    assert folder == str(tmp_path / "Agent Chat")
    for fname in ab.PACKAGE_FILES + ["config.js"]:
        assert os.path.isfile(os.path.join(folder, fname)), fname
    cfg = read("config.js", folder)
    assert cfg.startswith("//")
    m = re.search(r"window\.CHEMDRAWMCP_BRIDGE = (\{.*?\});", cfg, re.S)
    assert m, cfg
    assert json.loads(m.group(1)) == {"port": 47999, "token": "s3cret-token"}
    assert js_balance(cfg) is None
    meta = json.loads(read(ab.METADATA_FILE, folder))
    assert meta["name"] == "Agent Chat" and meta["menuItemText"] == "Agent Chat" and meta["url"] == "main.html"
    assert read("app.js", folder) == read("app.js")
    assert read("main.html", folder) == read("main.html")


def test_build_is_idempotent_and_renames(tmp_path):
    ab.build_addin(str(tmp_path), 47325, "one")
    folder = ab.build_addin(str(tmp_path), 47326, "two", name="Agent Chat Dev")
    assert os.path.basename(folder) == "Agent Chat Dev"
    meta = json.loads(read(ab.METADATA_FILE, folder))
    assert meta["name"] == meta["menuItemText"] == "Agent Chat Dev"
    assert '"port": 47326' in read("config.js", folder) and '"token": "two"' in read("config.js", folder)
    # rebuilding the default one overwrites config.js
    ab.build_addin(str(tmp_path), 47327, "three")
    assert '"token": "three"' in read("config.js", str(tmp_path / "Agent Chat"))


def test_build_rejects_bad_name_and_port(tmp_path):
    with pytest.raises(ValueError):
        ab.build_addin(str(tmp_path), 47325, "t", name="bad/name")
    with pytest.raises(ValueError):
        ab.build_addin(str(tmp_path), 0, "t")
    with pytest.raises(ValueError):
        ab.config_js(70000, "t")


def test_package_zip_has_main_html_at_top_level(tmp_path):
    ab.build_addin(str(tmp_path), 47325, "tok")
    out = ab.package_addin(str(tmp_path))
    assert out == str(tmp_path / ("Agent Chat" + ab.PACKAGE_SUFFIX))
    with zipfile.ZipFile(out) as zf:
        assert zf.testzip() is None
        names = zf.namelist()
        assert "main.html" in names
        assert "chemdraw-addin-metadata.json" in names
        assert "config.js" in names and "app.js" in names and "style.css" in names
        assert not any(n.startswith("Agent Chat/") for n in names), names
        assert json.loads(zf.read("chemdraw-addin-metadata.json"))["url"] == "main.html"
        assert b"47325" in zf.read("config.js")


def test_package_without_build_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ab.package_addin(str(tmp_path), "Nope")


def test_remove_addin(tmp_path):
    ab.build_addin(str(tmp_path), 47325, "tok")
    assert ab.remove_addin(str(tmp_path)) is True
    assert not os.path.exists(tmp_path / "Agent Chat")
    assert ab.remove_addin(str(tmp_path)) is False
    os.makedirs(tmp_path / "Other")
    with pytest.raises(RuntimeError):
        ab.remove_addin(str(tmp_path), "Other")   # not an add-in folder


def test_addin_install_dir_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("CHEMDRAWMCP_ADDIN_DIR", str(tmp_path))
    assert ab.addin_install_dir() == str(tmp_path)
    assert ab.installed_addin_dir("X") == os.path.join(str(tmp_path), "X")
    assert ab.is_installed() is False
    ab.build_addin(ab.addin_install_dir(), 47325, "tok")
    assert ab.is_installed() is True


def test_addin_install_dir_platform_default(monkeypatch):
    monkeypatch.delenv("CHEMDRAWMCP_ADDIN_DIR", raising=False)
    d = ab.addin_install_dir()
    assert os.path.isabs(d)
    if sys.platform == "darwin":
        assert d.endswith(os.path.join("com.revvity.ChemDraw", "Add-ins"))
        assert "Library/Application Support" in d
    elif sys.platform.startswith("win"):
        assert "RevvitySignalsSoftware" in d and d.endswith("Add-ins")
    else:
        assert d.endswith("addins")


def test_bridge_connection_resolution(monkeypatch, tmp_path):
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path))      # no bridge.json here
    monkeypatch.setenv("CHEMDRAWMCP_BRIDGE_PORT", "48000")
    monkeypatch.setenv("CHEMDRAWMCP_BRIDGE_TOKEN", "env-token")
    conn = ab.bridge_connection()
    assert conn["port"] == 48000 and conn["token"] == "env-token"
    conn = ab.bridge_connection(port=48001, token="explicit")
    assert conn["port"] == 48001 and conn["token"] == "explicit"
    monkeypatch.delenv("CHEMDRAWMCP_BRIDGE_PORT")
    monkeypatch.delenv("CHEMDRAWMCP_BRIDGE_TOKEN")
    conn = ab.bridge_connection()
    assert conn["port"] == ab.DEFAULT_PORT and conn["token"] == ""
    # a bridge.json written by the bridge wins over the defaults
    with open(tmp_path / "bridge.json", "w", encoding="utf-8") as f:
        json.dump({"port": 48123, "token": "from-bridge-json-0123456789", "pid": 1}, f)
    conn = ab.bridge_connection()
    assert conn["port"] == 48123 and conn["token"] == "from-bridge-json-0123456789" and conn["source"] == "bridge.json"


def test_cli_install_check_package_uninstall(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("CHEMDRAWMCP_ADDIN_DIR", str(tmp_path / "Add-ins"))
    monkeypatch.setenv("CHEMDRAWMCP_HOME", str(tmp_path / "home"))
    assert ab.main(["--check"]) == 0
    out = capsys.readouterr().out
    assert "not installed" in out and "MISSING" in out
    assert ab.main(["--port", "47555", "--token", "cli-token", "--package", str(tmp_path / "pkg")]) == 0
    out = capsys.readouterr().out
    assert "Installed add-in" in out and "Package:" in out
    assert os.path.isfile(tmp_path / "Add-ins" / "Agent Chat" / "config.js")
    assert '"token": "cli-token"' in read("config.js", str(tmp_path / "Add-ins" / "Agent Chat"))
    assert os.path.isfile(tmp_path / "pkg" / "Agent Chat.chemdrawaddin")
    assert ab.main(["--check"]) == 0
    assert "(installed)" in capsys.readouterr().out
    assert ab.main(["--uninstall"]) == 0
    assert not os.path.exists(tmp_path / "Add-ins" / "Agent Chat")


def test_version_is_two_part():
    """ChemDraw 23 silently drops add-ins whose version has more than two components (verified: "0.1.0" was not
    listed in the Add-ins menu, "0.1" was)."""
    from chemdrawmcp import addin_builder as _ab
    assert _ab.metadata()["version"].count(".") <= 1
