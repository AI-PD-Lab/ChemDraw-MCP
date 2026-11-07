/* Agent Chat for ChemDraw - add-in UI.
 *
 * Plain ES2017 for the WebKit web view inside ChemDraw (also runs in any browser for development).
 * No frameworks, no bundler, no network except the local agent bridge on 127.0.0.1.
 * Look & behaviour mirror molmcp/inapp/shared/agentchat_ui.py (the Claude Code terminal UI):
 * `>` user prompts, `⏺` assistant / tool lines, `⎿` result previews, animated status line, bordered `>` input
 * with a slash-command popup, a provider · model · effort chooser and a shortcuts footer.
 *
 * Bridge contract (chemdrawmcp/bridge/server.py):
 *   GET  /health                 -> {ok, version, backends:{claude,codex}, default_backend, settings, runs_active}
 *   POST /run {backend, prompt, context?} -> {run_id}
 *   GET  /events/<run_id>?token= -> SSE: session|delta|text|tool|tool_result|stats|done|error|end (data: JSON)
 *   POST /cancel/<run_id>, GET/POST /settings, POST /clear {backend}, GET /sessions
 * Auth: header X-Bridge-Token on fetch, ?token= on EventSource.  Port/token come from config.js
 * (window.CHEMDRAWMCP_BRIDGE = {port, token}, written by `python3 install.py --addin`).
 *
 * ChemDrawAPI (feature-detected, never assumed): activeDocument.selection.{isEmpty,getSMILES,getCDXML,
 * numberOfObjects,numberOfStructures,containsPartialStructure,onChange}, activeDocument.add{SMILES,CDXML,MolV2000,...},
 * window.{setDefaultSize,resizeTo,close,onClose}.
 */
(function () {
'use strict';

// ------------------------------------------------------------------------------------------------
// constants
// ------------------------------------------------------------------------------------------------
var DEFAULT_PORT = 47325;
var SPINNER = ['·', '✢', '✳', '✶', '✻', '✽', '✻', '✶', '✳', '✢'];
var VERBS = ['Thinking', 'Pondering', 'Working', 'Inspecting', 'Rendering', 'Reasoning', 'Composing'];
var BACKENDS = ['claude', 'codex'];
var LABELS = { claude: 'Claude Code', codex: 'Codex' };
var MODELS = {
  claude: [['', 'default'], ['claude-opus-5-5', 'opus-5-5'], ['claude-sonnet-5-5', 'sonnet-5-5'],
           ['claude-haiku-5-5', 'haiku-5-5'], ['claude-fable-5-1', 'fable-5-1']],
  codex: [['', 'default'], ['gpt-6-sol', 'gpt-6-sol'], ['gpt-6', 'gpt-6']]
};
var EFFORTS = ['', 'low', 'medium', 'high', 'xhigh', 'max'];
var SLASH = [
  ['/model', 'choose provider, model and effort'],
  ['/clear', 'new conversation'],
  ['/claude', 'switch to Claude Code'],
  ['/codex', 'switch to Codex'],
  ['/skills', 'ask the agent to list its skills'],
  ['/stop', 'interrupt the current run'],
  ['/help', 'shortcuts and commands']
];
var HELP = '/model   provider · model · effort (or the button bottom right)\n' +
  '/claude  /codex   switch provider        /clear   new conversation\n' +
  '/skills  list the agent\'s skill files     /stop    interrupt\n' +
  'Enter sends · Shift+Enter newline · Esc interrupts · ↑/↓ history · Tab accepts a / suggestion\n' +
  'A non-empty ChemDraw selection is attached to every prompt (SMILES, or CDXML for partial structures).';
var MAX_ENTRIES = 200;
var MAX_HISTORY = 100;
var HEALTH_EVERY_MS = 5000;
var CONTEXT_LIMIT = 200000;
var STORE_PREFIX = 'chemdrawmcp.addin.';
var START_BRIDGE = 'Start it with `chemdraw-agent-bridge` in a terminal, or use ChemDraw ▸ Scripts ▸ Start Agent Bridge.';
var CONFIGURE = 'Run `python3 install.py --addin` in the chemdraw-mcp folder (it writes config.js with the bridge port and token next to main.html), then reopen Agent Chat.';
var INSERT_METHODS = {
  smiles: 'addSMILES', smi: 'addSMILES',
  cdxml: 'addCDXML',
  mol: 'addMolV2000', molfile: 'addMolV2000', sdf: 'addMolV2000', molv2000: 'addMolV2000',
  molv3000: 'addMolV3000', mol3000: 'addMolV3000',
  inchi: 'addInChI',
  rxn: 'addRXNV2000', rxnv2000: 'addRXNV2000', rxnv3000: 'addRXNV3000'
};

// ------------------------------------------------------------------------------------------------
// configuration (config.js is optional; without it: port 47325, no token -> the bridge answers 401)
// ------------------------------------------------------------------------------------------------
var rawCfg = (typeof window.CHEMDRAWMCP_BRIDGE === 'object' && window.CHEMDRAWMCP_BRIDGE) || null;
var cfg = {
  port: (rawCfg && parseInt(rawCfg.port, 10)) || DEFAULT_PORT,
  token: (rawCfg && typeof rawCfg.token === 'string') ? rawCfg.token : ''
};
var hasConfig = !!(rawCfg && cfg.token);
var BASE = 'http://127.0.0.1:' + cfg.port;

// ------------------------------------------------------------------------------------------------
// DOM
// ------------------------------------------------------------------------------------------------
function $(sel) { return document.querySelector(sel); }
var el = {
  bridge: $('#bridge'), dot: $('#bridge-dot'), bridgeText: $('#bridge-text'), help: $('#bridge-help'),
  transcript: $('#transcript'), status: $('#status'), suggest: $('#suggest'), prompt: $('#prompt'),
  modelBtn: $('#model-btn'), selChip: $('#sel-chip'), version: $('#version'), chooser: $('#chooser'),
  close: $('#btn-close')
};

// ------------------------------------------------------------------------------------------------
// state
// ------------------------------------------------------------------------------------------------
var backend = 'claude';
var entries = [];            // transcript of the current backend: {k, t, name, args, c, ctx}
var welcomeEl = null;
var run = null;              // active run: {id, es, t0, verb, spin, stats, draft, finished, cancelled, timer}
var history = [];
var histPos = 0;
var suggestItems = [];
var suggestIdx = 0;
var chooserBackend = null;
var stickToBottom = true;
var helpDismissedFor = '';
var everConnected = false;
var saveTimer = null;
var draftTimer = null;
var codeSeq = 0;
var codeStore = {};          // id -> {method, code} for the "insert into document" buttons
var bridge = {
  state: 'unknown',          // unknown | ok | down | unauthorized
  version: '', backends: {}, default_backend: '', runs_active: 0, sessions: {},
  settings: { claude: { model: '', effort: '' }, codex: { model: '', effort: '' } }
};

// ------------------------------------------------------------------------------------------------
// small helpers
// ------------------------------------------------------------------------------------------------
var store = {
  get: function (key, fallback) {
    try {
      var v = window.localStorage.getItem(STORE_PREFIX + key);
      return v === null ? fallback : JSON.parse(v);
    } catch (e) { return fallback; }
  },
  set: function (key, value) {
    try { window.localStorage.setItem(STORE_PREFIX + key, JSON.stringify(value)); } catch (e) { /* no storage */ }
  },
  del: function (key) {
    try { window.localStorage.removeItem(STORE_PREFIX + key); } catch (e) { /* no storage */ }
  }
};

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}
function escAttr(s) {
  return esc(s).replace(/"/g, '&quot;');
}
function isHidden(node) { return node.classList.contains('hidden'); }
function show(node) { node.classList.remove('hidden'); }
function hide(node) { node.classList.add('hidden'); }
function shortModel(m) { return (m || 'default').replace(/^claude-/, ''); }
function settingsFor(b) { return bridge.settings[b] || { model: '', effort: '' }; }
function modelLabel(b) {
  b = b || backend;
  var st = settingsFor(b);
  return b + ' · ' + shortModel(st.model) + (st.effort ? ' · ' + st.effort : '');
}
function parseData(s) {
  if (typeof s !== 'string') return s;
  try { return JSON.parse(s); } catch (e) { return s; }
}

// ------------------------------------------------------------------------------------------------
// ChemDrawAPI (guarded: the page also runs in a plain browser)
// ------------------------------------------------------------------------------------------------
function cdApi() {
  try { return (typeof ChemDrawAPI !== 'undefined' && ChemDrawAPI) ? ChemDrawAPI : null; } catch (e) { return null; }
}
function cdDocument() {
  var a = cdApi();
  try { return (a && a.activeDocument) || null; } catch (e) { return null; }
}
function cdSelection() {
  var d = cdDocument();
  try { return (d && d.selection) || null; } catch (e) { return null; }
}
function cdCall(obj, name) {
  // call obj[name]() if it is a function (or return the property), swallowing ChemDraw's exceptions
  try {
    if (!obj) return undefined;
    var v = obj[name];
    if (typeof v === 'function') return v.apply(obj, Array.prototype.slice.call(arguments, 2));
    return v;
  } catch (e) { return undefined; }
}

var windowReady = false;
function setupWindow() {
  var a = cdApi();
  if (!a || !a.window || windowReady) return !!windowReady;
  windowReady = true;
  var w = a.window;
  try {
    if (typeof w.setDefaultSize === 'function') w.setDefaultSize(560, 720);
    var tooSmall = (window.innerWidth && window.innerWidth < 420) || (window.innerHeight && window.innerHeight < 480);
    if (typeof w.resizeTo === 'function' && (tooSmall || typeof w.setDefaultSize !== 'function')) w.resizeTo(560, 720);
  } catch (e) { /* ignore */ }
  if (typeof w.close === 'function') show(el.close);
  try { if (typeof w.onClose === 'function') w.onClose(function () { saveTranscript(true); }); } catch (e) { /* ignore */ }
  return true;
}

var selectionHooked = false;
function hookSelection() {
  var sel = cdSelection();
  if (!sel || selectionHooked) return;
  selectionHooked = true;
  try { if (typeof sel.onChange === 'function') sel.onChange(function () { updateSelectionChip(); }); } catch (e) { /* ignore */ }
  updateSelectionChip();
}

function updateSelectionChip() {
  var sel = cdSelection();
  if (!sel) { hide(el.selChip); return; }
  var text = 'no selection';
  var title = 'Nothing is selected in ChemDraw; prompts are sent without a selection';
  try {
    var empty = cdCall(sel, 'isEmpty');
    if (empty === false) {
      var s = cdCall(sel, 'getSMILES');
      s = (typeof s === 'string') ? s.trim() : '';
      var n = cdCall(sel, 'numberOfObjects');
      var what = s ? (s.length > 36 ? s.slice(0, 35) + '…' : s)
                   : (typeof n === 'number' ? n + ' object' + (n === 1 ? '' : 's') : 'objects');
      text = 'selection: ' + what;
      title = 'Attached to the next prompt' + (s ? ' as SMILES: ' + s : ' (CDXML)');
    }
  } catch (e) { text = 'selection: ?'; }
  el.selChip.textContent = text;
  el.selChip.title = title;
  show(el.selChip);
}

function collectContext() {
  // the current selection (SMILES, or CDXML when SMILES is not available) for the agent
  var sel = cdSelection();
  if (!sel) return null;
  if (cdCall(sel, 'isEmpty') !== false) return null;
  var smiles = cdCall(sel, 'getSMILES');
  smiles = (typeof smiles === 'string') ? smiles.trim() : '';
  var cdxml = '';
  if (!smiles) {
    cdxml = cdCall(sel, 'getCDXML');
    cdxml = (typeof cdxml === 'string') ? cdxml.trim() : '';
  }
  var payload = smiles || cdxml;
  if (!payload) return null;
  if (payload.length > CONTEXT_LIMIT) {
    pushEntry({ k: 'note', t: '(selection too large to attach — ' + payload.length + ' characters)' });
    return null;
  }
  var n = cdCall(sel, 'numberOfStructures');
  var partial = cdCall(sel, 'containsPartialStructure');
  var fmt = smiles ? 'SMILES' : 'CDXML';
  var text = 'The user has selected ' + (typeof n === 'number' ? n + ' structure(s)' : 'objects') +
    ' in the active ChemDraw document' + (partial === true ? ' (the selection contains a partial structure)' : '') +
    '. Selection as ' + fmt + ':\n' + payload;
  return { text: text, label: 'selection attached (' + fmt + ')' };
}

function insertMethodFor(lang) {
  var method = INSERT_METHODS[String(lang || '').toLowerCase()];
  if (!method) return null;
  var doc = cdDocument();
  try { return (doc && typeof doc[method] === 'function') ? method : null; } catch (e) { return null; }
}

function insertCode(id) {
  var item = codeStore[id];
  var doc = cdDocument();
  if (!item || !doc) return;
  try {
    doc[item.method](item.code);
    pushEntry({ k: 'note', t: '⤓ inserted into the ChemDraw document (' + item.method + ')' });
  } catch (e) {
    pushEntry({ k: 'error', t: 'ChemDraw refused the insert: ' + ((e && e.message) || e) });
  }
  updateSelectionChip();
}

// ------------------------------------------------------------------------------------------------
// bridge API
// ------------------------------------------------------------------------------------------------
function api(method, path, body) {
  var headers = {};
  if (cfg.token) headers['X-Bridge-Token'] = cfg.token;
  var init = { method: method, headers: headers, cache: 'no-store' };
  if (body !== undefined) {
    headers['Content-Type'] = 'application/json';
    init.body = JSON.stringify(body);
  }
  return window.fetch(BASE + path, init).then(function (r) {
    return r.text().then(function (t) {
      var data = null;
      try { data = t ? JSON.parse(t) : null; } catch (e) { data = { raw: t }; }
      if (!r.ok) {
        var msg = (data && (data.error || data.message || data.detail)) || (method + ' ' + path + ' → HTTP ' + r.status);
        var err = new Error(String(msg));
        err.status = r.status;
        err.data = data;
        throw err;
      }
      return data;
    });
  });
}

function describeError(e) {
  if (e && (e.status === 401 || e.status === 403)) {
    return 'the bridge rejected the token (HTTP ' + e.status + ') — ' + CONFIGURE;
  }
  if (e && e.status) return e.message || ('HTTP ' + e.status);
  return 'cannot reach the agent bridge at 127.0.0.1:' + cfg.port + ' — ' + START_BRIDGE;
}

function applySettings(s) {
  if (!s || typeof s !== 'object') return;
  var src = (s.settings && typeof s.settings === 'object') ? s.settings : s;
  BACKENDS.forEach(function (b) {
    var st = src[b];
    if (st && typeof st === 'object') {
      bridge.settings[b] = {
        model: typeof st.model === 'string' ? st.model : (bridge.settings[b].model || ''),
        effort: typeof st.effort === 'string' ? st.effort : (bridge.settings[b].effort || '')
      };
    }
  });
  if (typeof s.default_backend === 'string') bridge.default_backend = s.default_backend;
}

function sessionIdOf(data) {
  if (!data) return '';
  if (typeof data === 'string') return data;
  if (typeof data !== 'object') return String(data);
  return String(data.session_id || data.id || data.session || '');
}

function fetchSessions() {
  return api('GET', '/sessions').then(function (s) {
    var src = (s && s.sessions && typeof s.sessions === 'object') ? s.sessions : s;
    var out = {};
    if (src && typeof src === 'object') {
      BACKENDS.forEach(function (b) { out[b] = sessionIdOf(src[b]); });
    }
    bridge.sessions = out;
    renderWelcome();
  }).catch(function () { /* optional endpoint */ });
}

function checkHealth() {
  return api('GET', '/health').then(function (h) {
    var was = bridge.state;
    bridge.state = 'ok';
    bridge.version = (h && h.version) ? String(h.version) : '';
    bridge.backends = (h && h.backends && typeof h.backends === 'object') ? h.backends : {};
    bridge.runs_active = (h && h.runs_active) || 0;
    if (h && typeof h.default_backend === 'string') bridge.default_backend = h.default_backend;
    applySettings(h && h.settings);
    if (!everConnected) { everConnected = true; firstConnect(); }
    if (was !== 'ok') { renderBridge(); renderWelcome(); }
    updateFooter();
  }).catch(function (e) {
    var was = bridge.state;
    bridge.state = (e && (e.status === 401 || e.status === 403)) ? 'unauthorized' : 'down';
    if (was !== bridge.state) { renderBridge(); renderWelcome(); }
  });
}

function isInstalled(b) {
  if (bridge.state !== 'ok') return true;            // unknown until the bridge answers
  return !!bridge.backends[b];
}

function firstConnect() {
  var stored = store.get('backend', '');
  if (!stored && bridge.default_backend && BACKENDS.indexOf(bridge.default_backend) >= 0) {
    setBackend(bridge.default_backend);
  } else if (!isInstalled(backend)) {
    var other = BACKENDS.filter(isInstalled)[0];
    if (other) {
      setBackend(other);
      pushEntry({ k: 'note', t: 'agent: ' + modelLabel() + ' (the previous provider is not installed on this machine)' });
    }
  }
  fetchSessions();
}

// ------------------------------------------------------------------------------------------------
// bridge status in the header
// ------------------------------------------------------------------------------------------------
function problemText(state) {
  var where = '127.0.0.1:' + cfg.port;
  if (state === 'unauthorized') {
    return 'The agent bridge at ' + where + ' rejected this add-in\'s token (HTTP 401). ' +
      (hasConfig ? 'config.js next to main.html is stale. ' : 'config.js is missing next to main.html. ') + CONFIGURE;
  }
  if (state === 'down') {
    return 'No agent bridge is listening on ' + where + '. ' + START_BRIDGE + ' Retrying every 5 s…' +
      (hasConfig ? '' : '\n\nconfig.js is also missing — ' + CONFIGURE);
  }
  return hasConfig ? 'Connecting to the agent bridge at ' + where + '…' : 'config.js is missing next to main.html. ' + CONFIGURE;
}

function availableLabel() {
  var names = BACKENDS.filter(function (b) { return !!bridge.backends[b]; });
  return names.length ? names.join(' · ') : 'no agent CLI found';
}

function renderBridge() {
  var s = bridge.state;
  el.dot.className = 'dot ' + (s === 'ok' ? 'ok' : s === 'unauthorized' ? 'warn' : s === 'down' ? 'down' : '');
  el.bridgeText.textContent = s === 'ok' ? 'connected · ' + availableLabel()
    : s === 'unauthorized' ? 'bridge not configured'
    : s === 'down' ? 'bridge unreachable' : 'connecting…';
  el.version.textContent = bridge.version ? 'bridge v' + bridge.version : '';
  el.help.className = (s === 'unauthorized') ? 'warn' : '';
  if (s === 'ok') { hide(el.help); return; }
  el.help.innerHTML = mdInline(problemText(s)).replace(/\n/g, '<br>');
  if (helpDismissedFor === s) hide(el.help); else show(el.help);
}

function toggleHelp() {
  if (isHidden(el.help)) {
    el.help.className = bridge.state === 'unauthorized' ? 'warn' : (bridge.state === 'ok' ? 'ok' : '');
    el.help.innerHTML = mdInline(bridge.state === 'ok'
      ? 'Connected to the agent bridge at 127.0.0.1:' + cfg.port + (bridge.version ? ' (v' + bridge.version + ')' : '') +
        ' · agents: ' + availableLabel() + ' · active runs: ' + bridge.runs_active + '. ' +
        'Without the bridge: ' + START_BRIDGE
      : problemText(bridge.state)).replace(/\n/g, '<br>');
    helpDismissedFor = '';
    show(el.help);
  } else {
    helpDismissedFor = bridge.state;
    hide(el.help);
  }
}

// ------------------------------------------------------------------------------------------------
// markdown-ish rendering (inline code, bold, italics, bullet lists, headings, fenced code windows)
// ------------------------------------------------------------------------------------------------
function mdInline(s) {
  s = esc(s);
  s = s.replace(/`([^`\n]+)`/g, '<code>$1</code>');
  s = s.replace(/\*\*([^*\n]+)\*\*/g, '<b>$1</b>');
  s = s.replace(/(^|[^\w*])\*([^*\n]+)\*(?![\w*])/g, '$1<i>$2</i>');
  return s;
}

function codeWindow(code, lang, title) {
  var lines = String(code).replace(/\n+$/, '').split('\n');
  var rows = '';
  for (var i = 0; i < lines.length; i++) {
    rows += '<tr><td class="n">' + (i + 1) + '</td><td>' + esc(lines[i]) + '</td></tr>';
  }
  var method = insertMethodFor(lang);
  var head = '';
  if (title || lang || method) {
    head = '<div class="code-title">' + esc(title || '') + (title && lang ? ' · ' : '') + esc(lang || '');
    if (method) {
      var id = ++codeSeq;
      codeStore[id] = { method: method, code: lines.join('\n') };
      head += '<button type="button" class="insert" data-code="' + id + '" title="' + method + '">⤓ insert into document</button>';
    }
    head += '</div>';
  }
  return head + '<div class="code"><table>' + rows + '</table></div>';
}

function mdText(chunk) {
  var lines = chunk.replace(/^\n+|\n+$/g, '').split('\n');
  if (lines.length === 1 && !lines[0].trim()) return '';
  var html = '';
  var list = null;
  var para = [];
  function flushPara() { if (para.length) { html += '<p>' + para.join('<br>') + '</p>'; para = []; } }
  function flushList() { if (list) { html += '</' + list + '>'; list = null; } }
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i];
    var m = /^\s*(?:[-*•]|\d+[.)])\s+(.*)$/.exec(line);
    if (m) {
      flushPara();
      var kind = /^\s*\d/.test(line) ? 'ol' : 'ul';
      if (list !== kind) { flushList(); list = kind; html += '<' + kind + '>'; }
      html += '<li>' + mdInline(m[1]) + '</li>';
      continue;
    }
    m = /^\s*#{1,6}\s+(.*)$/.exec(line);
    if (m) { flushPara(); flushList(); html += '<strong class="h">' + mdInline(m[1]) + '</strong>'; continue; }
    if (!line.trim()) { flushPara(); flushList(); continue; }
    flushList();
    para.push(mdInline(line));
  }
  flushPara();
  flushList();
  return html;
}

function mdBlock(text) {
  // split on ``` fences: [text, lang, code, text, lang, code, ...]
  var parts = String(text).split(/```([a-zA-Z0-9_+-]*)[ \t]*\n?/);
  var out = '';
  for (var i = 0; i < parts.length; i++) {
    if (i % 3 === 0) {
      out += mdText(parts[i]);
    } else {
      out += codeWindow(parts[i + 1] || '', parts[i], '');
      i += 1;
    }
  }
  return out;
}

function preview(text, maxLines, width) {
  maxLines = maxLines || 2;
  width = width || 110;
  text = String(text == null ? '' : text);
  var stripped = text.trim();
  var first = stripped.charAt(0);
  if (first === '{' || first === '[') {
    try {
      var compact = JSON.stringify(JSON.parse(stripped));
      var total = text.split('\n').length;
      var cut = compact.length > width;
      return { shown: [cut ? compact.slice(0, width - 1) + '…' : compact], more: Math.max(0, total - 1), truncated: cut || total > 1 };
    } catch (e) { /* not JSON */ }
  }
  var lines = text.split('\n').filter(function (l) { return l.trim(); });
  var truncated = false;
  var shown = lines.slice(0, maxLines).map(function (l) {
    if (l.length <= width) return l;
    truncated = true;
    return l.slice(0, width - 1) + '…';
  });
  return { shown: shown, more: lines.length - shown.length, truncated: truncated || lines.length > shown.length };
}

function toolCallHtml(name, args) {
  var obj = null;
  var argText = '';
  if (args && typeof args === 'object') {
    obj = args;
    argText = JSON.stringify(args);
  } else {
    argText = String(args == null ? '' : args).trim();
    try { obj = argText ? JSON.parse(argText) : {}; } catch (e) { obj = null; }
  }
  var windows = '';
  var compact = argText;
  if (obj && typeof obj === 'object' && !Array.isArray(obj)) {
    var rest = {};
    Object.keys(obj).forEach(function (k) {
      var v = obj[k];
      if (typeof v === 'string' && (v.length > 160 || v.indexOf('\n') >= 0)) {
        windows += codeWindow(v, k === 'code' ? 'python' : '', k);
      } else if (Array.isArray(v) && k === 'commands') {
        windows += codeWindow(v.map(String).join('\n'), '', k);
      } else {
        rest[k] = v;
      }
    });
    compact = Object.keys(rest).length ? JSON.stringify(rest) : '';
  }
  var html = '<span class="bullet">⏺</span><div class="body"><span class="name">' + esc(name) + '</span>';
  if (compact || !windows) {
    var shortText = compact.length > 400 ? compact.slice(0, 399) + '…' : compact;
    html += '<span class="args" data-full="' + escAttr(compact) + '" title="click to expand">(' + esc(shortText) + ')</span>';
  }
  return html + windows + '</div>';
}

// ------------------------------------------------------------------------------------------------
// transcript
// ------------------------------------------------------------------------------------------------
function renderEntry(e) {
  var div = document.createElement('div');
  div.className = 'entry';
  switch (e.k) {
    case 'user':
      div.className += ' e-user';
      div.innerHTML = '<span class="gt">&gt;</span> ' + mdInline(e.t) +
        (e.ctx ? ' <span class="ctx">⌁ ' + esc(e.ctx) + '</span>' : '');
      break;
    case 'assistant':
      div.className += ' e-assistant';
      div.innerHTML = '<span class="bullet">⏺</span><div class="md">' + mdBlock(e.t) +
        (e.live ? '<span class="cursor">▌</span>' : '') + '</div>';
      break;
    case 'tool':
      div.className += ' e-tool';
      div.innerHTML = toolCallHtml(e.name, e.args);
      break;
    case 'tool_result': {
      div.className += ' e-result';
      var p = preview(e.t);
      var body = p.shown.map(esc).join('\n   ') || '(no output)';
      if (p.more > 0) body += '\n   <span class="more">… +' + p.more + ' lines</span>';
      div.innerHTML = '⎿  ' + body;
      div.title = p.truncated ? 'click to expand' : '';
      break;
    }
    case 'note':
      div.className += ' e-note' + (e.c ? ' ' + e.c : '');
      if (e.c === 'ok') div.innerHTML = '<span class="tick">✓</span>' + esc(String(e.t).replace(/^✓/, ''));
      else div.textContent = e.t;
      break;
    case 'error':
      div.className += ' e-error';
      div.textContent = '! ' + e.t;
      break;
    default:
      div.textContent = String(e.t || '');
  }
  div._entry = e;
  return div;
}

function plainEntry(e) {
  var out = {};
  ['k', 't', 'name', 'args', 'c', 'ctx'].forEach(function (key) { if (e[key] !== undefined) out[key] = e[key]; });
  return out;
}

function validEntry(e) {
  return e && typeof e === 'object' && typeof e.k === 'string' &&
    ['user', 'assistant', 'tool', 'tool_result', 'note', 'error'].indexOf(e.k) >= 0;
}

function scrollToBottom(force) {
  if (!force && !stickToBottom) return;
  el.transcript.scrollTop = el.transcript.scrollHeight;
}

function pushEntry(e) {
  entries.push(e);
  var node = renderEntry(e);
  e.el = node;
  el.transcript.appendChild(node);
  while (entries.length > MAX_ENTRIES) {
    var old = entries.shift();
    if (old.el && old.el.parentNode) old.el.parentNode.removeChild(old.el);
  }
  scrollToBottom();
  saveTranscript();
  return e;
}

function rerenderEntry(e) {
  var node = renderEntry(e);
  if (e.el && e.el.parentNode) e.el.parentNode.replaceChild(node, e.el);
  else el.transcript.appendChild(node);
  e.el = node;
  scrollToBottom();
}

function saveTranscript(now) {
  if (now) {
    if (saveTimer) { clearTimeout(saveTimer); saveTimer = null; }
    store.set('transcript.' + backend, entries.map(plainEntry));
    return;
  }
  if (saveTimer) return;
  saveTimer = setTimeout(function () { saveTimer = null; store.set('transcript.' + backend, entries.map(plainEntry)); }, 300);
}

function loadTranscript() {
  var saved = store.get('transcript.' + backend, []);
  entries = Array.isArray(saved) ? saved.filter(validEntry).slice(-MAX_ENTRIES) : [];
}

function welcomeHtml() {
  var who;
  if (bridge.state === 'ok') {
    who = isInstalled(backend) ? LABELS[backend] + ' (' + modelLabel() + ')' : LABELS[backend] + ' is not installed on this machine';
  } else if (bridge.state === 'unauthorized') {
    who = 'bridge not configured — run `python3 install.py --addin`';
  } else if (bridge.state === 'down') {
    who = 'agent bridge not running — `chemdraw-agent-bridge`';
  } else {
    who = LABELS[backend] + ' · connecting to the bridge';
  }
  var session = bridge.sessions && bridge.sessions[backend];
  var sessionText = session ? '<br>session ' + esc(String(session).slice(0, 8)) + '… continues (' + entries.length + ' entries restored)'
    : (entries.length ? '<br>' + entries.length + ' entries restored — the next prompt starts a new conversation on the bridge' : '');
  return '<span class="accent">✻</span> <b>Welcome to Agent Chat for ChemDraw</b><br><br>' +
    mdInline(who) + ' · bridge at 127.0.0.1:' + cfg.port + sessionText +
    '<br>Type a request below · /model to switch provider or model · /help for shortcuts';
}

function renderWelcome() {
  if (!welcomeEl) {
    welcomeEl = document.createElement('div');
    welcomeEl.className = 'welcome';
    el.transcript.insertBefore(welcomeEl, el.transcript.firstChild);
  }
  welcomeEl.innerHTML = welcomeHtml();
}

function renderAll() {
  el.transcript.innerHTML = '';
  welcomeEl = null;
  renderWelcome();
  entries.forEach(function (e) {
    e.el = renderEntry(e);
    el.transcript.appendChild(e.el);
  });
  stickToBottom = true;
  scrollToBottom(true);
}

// ------------------------------------------------------------------------------------------------
// backend / settings
// ------------------------------------------------------------------------------------------------
function updateFooter() {
  el.modelBtn.textContent = modelLabel() + ' ▾';
}

function setBackend(b) {
  if (BACKENDS.indexOf(b) < 0) return;
  if (b !== backend) {
    saveTranscript(true);
    backend = b;
    loadTranscript();
    renderAll();
  }
  store.set('backend', backend);
  updateFooter();
}

function choose(b, opts) {
  opts = opts || {};
  if (BACKENDS.indexOf(b) < 0) return Promise.resolve();
  var body = { backend: b, default_backend: b };
  var st = settingsFor(b);
  if (opts.model !== undefined) { body.model = opts.model; st.model = opts.model; }
  if (opts.effort !== undefined) { body.effort = opts.effort; st.effort = opts.effort; }
  bridge.settings[b] = st;
  bridge.default_backend = b;
  setBackend(b);
  return api('POST', '/settings', body).then(function (r) {
    applySettings(r);
  }).catch(function (e) {
    pushEntry({ k: 'error', t: 'settings were not saved on the bridge: ' + describeError(e) });
  }).then(function () {
    updateFooter();
    pushEntry({ k: 'note', t: 'agent: ' + modelLabel() });
    renderWelcome();
    renderChooser();
  });
}

function clearConversation() {
  api('POST', '/clear', { backend: backend }).catch(function () { /* bridge may be down */ });
  entries = [];
  store.del('transcript.' + backend);
  if (bridge.sessions) bridge.sessions[backend] = '';
  renderAll();
  pushEntry({ k: 'note', t: 'conversation cleared' });
}

// ------------------------------------------------------------------------------------------------
// model chooser (provider · model · effort)
// ------------------------------------------------------------------------------------------------
function openChooser() {
  chooserBackend = backend;
  show(el.chooser);
  renderChooser();
}
function closeChooser() { hide(el.chooser); }

function renderChooser() {
  if (isHidden(el.chooser)) return;
  var b = (chooserBackend && BACKENDS.indexOf(chooserBackend) >= 0) ? chooserBackend : backend;
  var st = settingsFor(b);
  var html = '<div class="col"><div class="col-title">Provider</div>';
  BACKENDS.forEach(function (k) {
    var installed = isInstalled(k);
    html += '<button type="button" class="row' + (k === b ? ' on' : '') + '" data-provider="' + k + '"' +
      (installed ? '' : ' disabled') + '>' + (k === backend ? '● ' : '  ') + LABELS[k] +
      (installed ? '' : '  (not installed)') + '</button>';
  });
  html += '<div class="hint">' + (b === backend ? 'current provider' : 'click to switch') + '</div></div>';
  html += '<div class="col"><div class="col-title">Model</div>';
  var known = false;
  MODELS[b].forEach(function (m) {
    var on = m[0] === st.model;
    known = known || on;
    html += '<button type="button" class="row' + (on ? ' on' : '') + '" data-model="' + escAttr(m[0]) + '">' +
      (on ? '● ' : '  ') + esc(m[1]) + '</button>';
  });
  if (!known && st.model) {
    html += '<button type="button" class="row on" data-model="' + escAttr(st.model) + '">● ' + esc(st.model) + '</button>';
  }
  html += '<div class="hint">other: /model &lt;name&gt;</div></div>';
  html += '<div class="col"><div class="col-title">Effort</div>';
  EFFORTS.forEach(function (e) {
    var on = e === (st.effort || '');
    html += '<button type="button" class="row' + (on ? ' on' : '') + '" data-effort="' + escAttr(e) + '">' +
      (on ? '● ' : '  ') + (e || 'default') + '</button>';
  });
  html += '</div>';
  el.chooser.innerHTML = html;
}

el.chooser.addEventListener('click', function (ev) {
  var btn = ev.target.closest('button');
  if (!btn || btn.disabled) return;
  var b = (chooserBackend && BACKENDS.indexOf(chooserBackend) >= 0) ? chooserBackend : backend;
  if (btn.dataset.provider !== undefined) {
    chooserBackend = btn.dataset.provider;
    if (chooserBackend !== backend) choose(chooserBackend); else renderChooser();
  } else if (btn.dataset.model !== undefined) {
    choose(b, { model: btn.dataset.model });
  } else if (btn.dataset.effort !== undefined) {
    choose(b, { effort: btn.dataset.effort });
  }
});

el.modelBtn.addEventListener('click', function (ev) {
  ev.stopPropagation();
  if (isHidden(el.chooser)) openChooser(); else closeChooser();
});

document.addEventListener('click', function (ev) {
  if (!isHidden(el.chooser) && !el.chooser.contains(ev.target) && ev.target !== el.modelBtn) closeChooser();
});

// ------------------------------------------------------------------------------------------------
// slash-command suggestions
// ------------------------------------------------------------------------------------------------
function hideSuggest() { suggestItems = []; hide(el.suggest); }

function renderSuggest() {
  el.suggest.innerHTML = suggestItems.map(function (it, i) {
    return '<div class="item' + (i === suggestIdx ? ' sel' : '') + '" role="option" data-i="' + i + '">' +
      '<span class="cmd">' + esc(it[0]) + '</span> <span class="desc">' + esc(it[1]) + '</span></div>';
  }).join('');
}

function updateSuggestions() {
  var t = el.prompt.value;
  if (t.charAt(0) === '/' && t.indexOf('\n') < 0 && t.indexOf(' ') < 0) {
    suggestItems = SLASH.filter(function (it) { return it[0].indexOf(t) === 0; });
    if (suggestItems.length) {
      suggestIdx = Math.min(suggestIdx, suggestItems.length - 1);
      renderSuggest();
      show(el.suggest);
      return;
    }
  }
  hideSuggest();
}

function acceptSuggestion() {
  if (!suggestItems.length) return false;
  el.prompt.value = suggestItems[suggestIdx][0] + ' ';
  hideSuggest();
  autogrow();
  el.prompt.focus();
  return true;
}

el.suggest.addEventListener('mousedown', function (ev) {
  var item = ev.target.closest('.item');
  if (!item) return;
  ev.preventDefault();
  suggestIdx = parseInt(item.dataset.i, 10) || 0;
  acceptSuggestion();
});

// ------------------------------------------------------------------------------------------------
// prompt input
// ------------------------------------------------------------------------------------------------
function autogrow() {
  el.prompt.style.height = 'auto';
  el.prompt.style.height = Math.min(160, Math.max(54, el.prompt.scrollHeight + 2)) + 'px';
}

function caretOnFirstLine() {
  return el.prompt.value.lastIndexOf('\n', Math.max(0, el.prompt.selectionStart - 1)) < 0;
}
function caretOnLastLine() {
  return el.prompt.value.indexOf('\n', el.prompt.selectionEnd) < 0;
}

function historyStep(step) {
  if (!history.length) return;
  histPos = Math.max(0, Math.min(history.length, histPos + step));
  el.prompt.value = histPos < history.length ? history[histPos] : '';
  el.prompt.selectionStart = el.prompt.selectionEnd = el.prompt.value.length;
  autogrow();
  updateSuggestions();
}

function onSubmit() {
  var text = el.prompt.value;
  el.prompt.value = '';
  hideSuggest();
  autogrow();
  submit(text);
}

el.prompt.addEventListener('keydown', function (ev) {
  var sVis = !isHidden(el.suggest) && suggestItems.length > 0;
  if (ev.key === 'Enter' && !ev.shiftKey) {
    ev.preventDefault();
    if (sVis && el.prompt.value.trim() !== suggestItems[suggestIdx][0]) { acceptSuggestion(); return; }
    onSubmit();
    return;
  }
  if (ev.key === 'Tab' && sVis) { ev.preventDefault(); acceptSuggestion(); return; }
  if (ev.key === 'Escape') {
    ev.preventDefault();
    if (!isHidden(el.chooser)) closeChooser();
    else if (sVis) hideSuggest();
    else stop();
    return;
  }
  if (ev.key === 'ArrowUp') {
    if (sVis) { ev.preventDefault(); suggestIdx = (suggestIdx - 1 + suggestItems.length) % suggestItems.length; renderSuggest(); return; }
    if (caretOnFirstLine()) { ev.preventDefault(); historyStep(-1); }
    return;
  }
  if (ev.key === 'ArrowDown') {
    if (sVis) { ev.preventDefault(); suggestIdx = (suggestIdx + 1) % suggestItems.length; renderSuggest(); return; }
    if (caretOnLastLine()) { ev.preventDefault(); historyStep(+1); }
  }
});
el.prompt.addEventListener('input', function () { updateSuggestions(); autogrow(); });
el.prompt.addEventListener('focus', function () { updateSelectionChip(); });

document.addEventListener('keydown', function (ev) {
  if (ev.key === 'Escape' && document.activeElement !== el.prompt) {
    if (!isHidden(el.chooser)) closeChooser(); else stop();
  }
});

// ------------------------------------------------------------------------------------------------
// transcript interactions: expand tool args / results, insert buttons
// ------------------------------------------------------------------------------------------------
el.transcript.addEventListener('click', function (ev) {
  var insert = ev.target.closest('button.insert');
  if (insert) { ev.preventDefault(); insertCode(insert.dataset.code); return; }
  var args = ev.target.closest('.args');
  if (args) {
    var full = args.classList.toggle('full');
    var text = args.dataset.full || '';
    args.textContent = '(' + (full || text.length <= 400 ? text : text.slice(0, 399) + '…') + ')';
    return;
  }
  var result = ev.target.closest('.e-result');
  if (result && result._entry) {
    if (window.getSelection && String(window.getSelection()).length) return;   // let the user copy text
    var open = result.classList.toggle('open');
    if (open) {
      result.textContent = '⎿  ' + String(result._entry.t || '(no output)');
    } else {
      var fresh = renderEntry(result._entry);
      result.parentNode.replaceChild(fresh, result);
      result._entry.el = fresh;
    }
  }
});

el.transcript.addEventListener('scroll', function () {
  stickToBottom = el.transcript.scrollHeight - el.transcript.scrollTop - el.transcript.clientHeight < 40;
});

el.bridge.addEventListener('click', toggleHelp);
el.close.addEventListener('click', function () {
  saveTranscript(true);
  var a = cdApi();
  try { if (a && a.window && typeof a.window.close === 'function') a.window.close(); } catch (e) { /* ignore */ }
});

// ------------------------------------------------------------------------------------------------
// running a prompt
// ------------------------------------------------------------------------------------------------
function submit(text) {
  text = (text || '').trim();
  if (!text) return;
  if (text.charAt(0) === '/') { slash(text); return; }
  if (run) { pushEntry({ k: 'note', t: '(still working — Esc to interrupt)' }); return; }
  history.push(text);
  if (history.length > MAX_HISTORY) history = history.slice(-MAX_HISTORY);
  histPos = history.length;
  store.set('history', history);
  var ctx = collectContext();
  pushEntry({ k: 'user', t: text, ctx: ctx ? ctx.label : '' });
  beginRun(text, ctx ? ctx.text : null);
}

function tick() {
  if (!run) return;
  run.spin = (run.spin + 1) % SPINNER.length;
  var secs = Math.floor((Date.now() - run.t0) / 1000);
  el.status.innerHTML = '<span class="accent">' + SPINNER[run.spin] + ' ' + run.verb + '…</span> ' +
    '<span class="dim">(esc to interrupt · ' + secs + 's)</span>';
}

function beginRun(prompt, context) {
  var me = { id: null, es: null, t0: Date.now(), verb: VERBS[history.length % VERBS.length], spin: 0,
             stats: null, draft: null, finished: false, cancelled: false, timer: null };
  run = me;
  show(el.status);
  tick();
  me.timer = setInterval(tick, 120);
  var body = { backend: backend, prompt: prompt };
  if (context) body.context = context;
  api('POST', '/run', body).then(function (r) {
    var id = r && (r.run_id || r.id);
    if (!id) throw new Error('the bridge did not return a run_id');
    me.id = String(id);
    if (me.cancelled || me.finished) {
      api('POST', '/cancel/' + encodeURIComponent(me.id)).catch(function () { /* ignore */ });
      return;
    }
    openEvents(me);
  }).catch(function (e) {
    if (me.finished) return;
    pushEntry({ k: 'error', t: describeError(e) });
    finishRun(true, me);
    if (bridge.state !== 'ok' || (e && e.status === 401)) checkHealth();
  });
}

function openEvents(me) {
  var url = BASE + '/events/' + encodeURIComponent(me.id) + (cfg.token ? '?token=' + encodeURIComponent(cfg.token) : '');
  var es;
  try { es = new window.EventSource(url); } catch (e) {
    pushEntry({ k: 'error', t: 'cannot open the event stream: ' + ((e && e.message) || e) });
    finishRun(true, me);
    return;
  }
  me.es = es;
  ['session', 'delta', 'text', 'tool', 'tool_result', 'stats', 'done', 'end'].forEach(function (kind) {
    es.addEventListener(kind, function (ev) { onEvent(me, kind, parseData(ev.data)); });
  });
  // "error" is both the SSE kind from the bridge (has .data) and the EventSource connection error
  es.addEventListener('error', function (ev) {
    if (ev && typeof ev.data === 'string') { onEvent(me, 'error', parseData(ev.data)); return; }
    if (me.finished) { es.close(); return; }
    if (es.readyState === window.EventSource.CLOSED) {
      pushEntry({ k: 'error', t: 'lost the event stream from the bridge' });
      finishRun(true, me);
    }
    // CONNECTING: the browser reconnects by itself
  });
  es.onmessage = function (ev) {
    // unnamed events: {event|kind|type, data|payload} envelopes
    var d = parseData(ev.data);
    if (d && typeof d === 'object') {
      var kind = d.event || d.kind || d.type;
      if (kind) onEvent(me, String(kind), d.data !== undefined ? d.data : (d.payload !== undefined ? d.payload : d));
    }
  };
}

var TEXT_KEYS = ['text', 'delta', 'content', 'message', 'error', 'output', 'result', 'value', 'data', 'payload'];
function textOf(data, fallbackJson) {
  if (data == null) return '';
  if (typeof data === 'string') return data;
  if (typeof data !== 'object') return String(data);
  if (Array.isArray(data)) return data.map(function (p) { return textOf(p, true); }).join('\n');
  if (data.type === 'image' || data.type === 'image_url' || data.type === 'input_image') return '[image]';
  for (var i = 0; i < TEXT_KEYS.length; i++) {
    var v = data[TEXT_KEYS[i]];
    if (typeof v === 'string') return v;
    if (Array.isArray(v)) return textOf(v, true);
  }
  return fallbackJson ? JSON.stringify(data, null, 1) : '';
}

function resultTextOf(data) {
  var t = textOf(data, true);
  if (data && typeof data === 'object' && !Array.isArray(data) &&
      (data.image || data.images || data.png || data.png_base64 || data.image_base64)) {
    t = (t ? t + '\n' : '') + '[image]';
  }
  t = t.replace(/data:image\/[a-zA-Z0-9.+-]+;base64,[A-Za-z0-9+/=]+/g, '[image]');
  if (/^[A-Za-z0-9+/=\s]{1500,}$/.test(t)) t = '[image]';
  return t;
}

function toolOf(data) {
  if (typeof data === 'string') {
    var sp = data.indexOf(' ');
    return sp < 0 ? [data, ''] : [data.slice(0, sp), data.slice(sp + 1).trim()];
  }
  if (data && typeof data === 'object') {
    var name = data.name || data.tool || data.tool_name || 'tool';
    var args = data.args;
    if (args === undefined) args = data.input;
    if (args === undefined) args = data.arguments;
    if (args === undefined) args = data.params;
    if (args === undefined) args = data.text;
    return [String(name), args === undefined ? '' : args];
  }
  return ['tool', ''];
}

function onEvent(me, kind, data) {
  if (me !== run || me.finished) return;
  switch (kind) {
    case 'session': {
      var id = sessionIdOf(data);
      if (id) { bridge.sessions[backend] = id; store.set('session.' + backend, id); }
      break;
    }
    case 'delta':
      addDelta(me, textOf(data));
      break;
    case 'text':
      finalizeDraft(me, textOf(data));
      break;
    case 'tool': {
      finalizeDraft(me, '');
      var t = toolOf(data);
      pushEntry({ k: 'tool', name: t[0], args: t[1] });
      break;
    }
    case 'tool_result':
      pushEntry({ k: 'tool_result', t: resultTextOf(data) });
      break;
    case 'stats':
      if (data && typeof data === 'object') me.stats = data;
      break;
    case 'error':
      pushEntry({ k: 'error', t: textOf(data, true) || 'unknown error' });
      break;
    case 'done': {
      var text = textOf(data);
      if (data && typeof data === 'object' && !me.stats && (data.stats || data.usage || data.total_cost_usd !== undefined)) {
        me.stats = data.stats || data;
      }
      finalizeDraft(me, text);
      finishRun(false, me);
      break;
    }
    case 'end':
      finishRun(false, me);
      break;
    default:
      break;
  }
}

function addDelta(me, text) {
  if (!text) return;
  if (!me.draft) me.draft = pushEntry({ k: 'assistant', t: '', live: true });
  me.draft.t += text;
  if (!draftTimer) {
    draftTimer = setTimeout(function () {
      draftTimer = null;
      if (me.draft) { rerenderEntry(me.draft); saveTranscript(); }
    }, 30);
  }
}

function finalizeDraft(me, text) {
  if (me.draft) {
    if (text) me.draft.t = text;
    me.draft.live = false;
    rerenderEntry(me.draft);
    me.draft = null;
    saveTranscript();
  } else if (text) {
    pushEntry({ k: 'assistant', t: text });
  }
}

function finishRun(failed, me) {
  me = me || run;
  if (!me || me.finished) return;
  me.finished = true;
  if (me.timer) clearInterval(me.timer);
  if (me.es) { try { me.es.close(); } catch (e) { /* ignore */ } }
  if (me.draft) { me.draft.live = false; rerenderEntry(me.draft); me.draft = null; }
  hide(el.status);
  el.status.textContent = '';
  var secs = ((Date.now() - me.t0) / 1000).toFixed(1);
  if (!failed) {
    var st = me.stats || {};
    var usage = st.usage || {};
    var toks = usage.output_tokens || st.output_tokens;
    var cost = st.total_cost_usd;
    if (cost == null) cost = st.cost_usd;
    if (cost == null) cost = st.cost;
    var extra = '';
    if (toks) extra += ' · ' + toks + ' output tokens';
    if (typeof cost === 'number') extra += ' · $' + cost.toFixed(2);
    pushEntry({ k: 'note', c: 'ok', t: '✓ ' + modelLabel() + ' · ' + secs + 's' + extra });
  }
  if (run === me) run = null;
  saveTranscript(true);
  updateSelectionChip();
  el.prompt.focus();
}

function stop() {
  if (!run) return;
  var me = run;
  me.cancelled = true;
  if (me.id) api('POST', '/cancel/' + encodeURIComponent(me.id)).catch(function () { /* ignore */ });
  pushEntry({ k: 'note', c: 'red', t: 'interrupted' });
  finishRun(true, me);
}

// ------------------------------------------------------------------------------------------------
// slash commands
// ------------------------------------------------------------------------------------------------
function slash(text) {
  var sp = text.indexOf(' ');
  var cmd = (sp < 0 ? text : text.slice(0, sp)).toLowerCase();
  var rest = sp < 0 ? '' : text.slice(sp + 1).trim();
  switch (cmd) {
    case '/clear':
      clearConversation();
      break;
    case '/claude':
    case '/codex':
      choose(cmd.slice(1)).then(function () { if (rest) submit(rest); });
      break;
    case '/model':
      if (!rest) { openChooser(); break; }
      var b = backend;
      var opts = {};
      rest.split(/\s+/).forEach(function (tok) {
        if (BACKENDS.indexOf(tok) >= 0) b = tok;
        else if (tok && EFFORTS.indexOf(tok) > 0) opts.effort = tok;
        else if (tok === 'default' || tok === 'none') opts.model = '';
        else if (tok) opts.model = tok;
      });
      choose(b, opts);
      break;
    case '/skills':
      submit('List your available skills (call the skills tool) and summarise them in one line each.');
      break;
    case '/stop':
      stop();
      break;
    case '/help':
    case '/?':
      pushEntry({ k: 'note', t: HELP });
      break;
    default:
      pushEntry({ k: 'error', t: 'unknown command ' + cmd + ' — /help' });
  }
}

// ------------------------------------------------------------------------------------------------
// start
// ------------------------------------------------------------------------------------------------
function init() {
  var storedBackend = store.get('backend', '');
  backend = BACKENDS.indexOf(storedBackend) >= 0 ? storedBackend : 'claude';
  var savedHistory = store.get('history', []);
  history = Array.isArray(savedHistory) ? savedHistory.filter(function (h) { return typeof h === 'string'; }).slice(-MAX_HISTORY) : [];
  histPos = history.length;
  bridge.sessions = {};
  BACKENDS.forEach(function (b) { var s = store.get('session.' + b, ''); if (typeof s === 'string') bridge.sessions[b] = ''; });
  loadTranscript();
  renderAll();
  updateFooter();
  renderBridge();
  if (!hasConfig) { el.help.innerHTML = mdInline(problemText('unknown')); show(el.help); }
  checkHealth();
  setInterval(checkHealth, HEALTH_EVERY_MS);
  // ChemDrawAPI may be injected a moment after the page loads
  [0, 250, 1000, 3000].forEach(function (ms) {
    setTimeout(function () { setupWindow(); hookSelection(); }, ms);
  });
  window.addEventListener('pagehide', function () { saveTranscript(true); });
  window.addEventListener('beforeunload', function () { saveTranscript(true); });
  autogrow();
  el.prompt.focus();
}

if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
else init();

// a few internals for development / debugging in a plain browser
window.AgentChat = { submit: submit, stop: stop, choose: choose, checkHealth: checkHealth, bridge: bridge,
                     entries: function () { return entries; }, backend: function () { return backend; } };
})();
