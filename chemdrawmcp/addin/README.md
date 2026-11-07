# Agent Chat — ChemDraw add-in

A ChemDraw add-in (macOS and Windows) that puts a Claude-Code-style chat panel inside ChemDraw.
It talks to the local **agent bridge** (`chemdraw-agent-bridge`, `chemdrawmcp/bridge/server.py`), which
runs Claude Code or Codex with the `chemdraw` MCP server, so the agent can draw, read, name and convert
in the open document while you chat.

```
✻ Agent Chat for ChemDraw                      ● connected · claude · codex
┌──────────────────────────────────────────────────────────────────────────┐
│ > draw the Suzuki coupling of 4-bromotoluene and phenylboronic acid      │
│ ⏺ chemdraw/draw_reaction({"reactants": ["Cc1ccc(Br)cc1", ...], ...})    │
│   ⎿  {"document": "Untitled-1", "window_id": 12, ...}                    │
│ ⏺ Drawn in a new document in ACS style: ...                              │
│ ✓ claude · sonnet-5-5 · high · 22.9s · 328 output tokens · $0.01         │
└──────────────────────────────────────────────────────────────────────────┘
> ▌                                                     claude · sonnet-5-5 · high ▾
```

## Files

| file | purpose |
| --- | --- |
| `chemdraw-addin-metadata.json` | what ChemDraw reads: name, menu item, `url: main.html`, `minimumAPIVersion: 1.0` |
| `main.html` | the page shown in ChemDraw's add-in web view (loads `config.js`, then `app.js`) |
| `app.js` | the whole UI: transcript, streaming, slash commands, model chooser, bridge client (plain ES2017, no dependencies) |
| `style.css` | the dark terminal look (`#1F1F1F` background, orange `#D97757` accent) |
| `config.example.js` | template for `config.js` |
| `config.js` | **written by the installer, not shipped** — `window.CHEMDRAWMCP_BRIDGE = {port, token}` |

## Install

```
python3 install.py --addin          # from the chemdraw-mcp checkout (creates the venv if needed)
python -m chemdrawmcp.addin_builder # inside the venv: same thing
```

Both copy this folder to the user's ChemDraw add-in directory and write `config.js` with the bridge port
and token (taken from the running bridge's `~/.chemdrawmcp/bridge.json`, else `$CHEMDRAWMCP_BRIDGE_PORT` /
`$CHEMDRAWMCP_BRIDGE_TOKEN`; pass `--port` / `--token` explicitly otherwise):

* macOS: `~/Library/Application Support/com.revvity.ChemDraw/Add-ins/Agent Chat/`
* Windows: `C:\ProgramData\RevvitySignalsSoftware\ChemDrawApplications\ChemDraw\ChemDraw Items\Add-ins\Agent Chat\`
* override with `CHEMDRAWMCP_ADDIN_DIR=/path/to/Add-ins`

Then start the bridge (`chemdraw-agent-bridge`, or ChemDraw ▸ Scripts ▸ *Start Agent Bridge*) and pick
**Add-ins ▸ Agent Chat** in ChemDraw (restart ChemDraw once if the menu does not list it yet).

`python -m chemdrawmcp.addin_builder --package DIR` additionally zips the folder as `Agent Chat.chemdrawaddin`
(files at the top level of the archive) for ChemDraw's Add-in Manager; `--uninstall` removes the folder;
`--check` prints the locations.

If the bridge's token changes (`chemdraw-agent-bridge --new-token`) rerun the installer: the panel shows
**bridge not configured** with that instruction whenever the bridge answers 401.

## How `config.js` is consumed

`main.html` loads `config.js` before `app.js`.  `app.js` reads `window.CHEMDRAWMCP_BRIDGE`:

```js
window.CHEMDRAWMCP_BRIDGE = { port: 47325, token: "…" };
```

* base URL `http://127.0.0.1:<port>` (the host is fixed — the bridge only binds loopback)
* `X-Bridge-Token: <token>` on every `fetch`, `?token=<token>` on the `EventSource` (which cannot set headers)
* when `config.js` is missing (ChemDraw's web view ignores the 404) the page falls back to port 47325 and
  no token; the bridge answers 401 and the header shows a red dot with *bridge not configured* and
  `python3 install.py --addin`

## Bridge protocol used by the page

| call | use |
| --- | --- |
| `GET /health` | every 5 s: status dot, available backends, settings snapshot (models/efforts), session ids |
| `POST /run {backend, prompt, context?}` | one turn; `context` is the ChemDraw selection (SMILES, or CDXML) |
| `GET /events/<run_id>?token=` | SSE `session · delta · text · tool · tool_result · stats · done · error · end` |
| `POST /cancel/<run_id>` | Esc / `/stop` |
| `POST /settings {backend, model, effort, default_backend}` | `/model`, `/claude`, `/codex`, the chooser |
| `POST /clear {backend}` | `/clear` (forgets the session on the bridge; the transcript is wiped locally) |
| `GET /sessions` | only when `/health` carries no `sessions` |

Everything else is local: the transcript (last 200 entries per provider), the prompt history and the
chosen provider live in `localStorage`; the bridge keeps the session ids so follow-up turns continue the
conversation after the panel was closed and reopened.

## Keyboard

* **Enter** send · **Shift+Enter** newline · **Esc** interrupt (or close the chooser/popup)
* **↑ / ↓** prompt history (or move in the slash popup) · **Tab** accept the slash suggestion
* `/model` (chooser: provider · model · effort), `/claude`, `/codex`, `/clear`, `/skills`, `/stop`, `/help`
* click a `⎿` result or a tool's `(args)` to expand it

## ChemDrawAPI usage (all feature-detected)

The page works in a plain browser too (open `main.html` from disk for development); every
`ChemDrawAPI` access sits behind `typeof` checks and `try/catch`.

* `ChemDrawAPI.window.setDefaultSize(560, 720)` on load (`resizeTo` as fallback / when the window is tiny),
  `window.close()` for the × button, `window.onClose` to flush the transcript
* `ChemDrawAPI.activeDocument.selection`: `isEmpty()`, `getSMILES()` (falls back to `getCDXML()` when SMILES
  throws, e.g. partial structures), `numberOfObjects()`, `numberOfStructures()`, `containsPartialStructure()`,
  `onChange(cb)` → footer chip *selection: C9H8O4…* / *no selection*; a non-empty selection is sent as
  `context` with each prompt and the `>` line gets a blue *selection attached* tag
* fenced blocks in an answer tagged `smiles`, `cdxml`, `mol`/`molfile`, `molv3000`, `inchi`, `rxn` get an
  **⤓ insert into document** button (`addSMILES`, `addCDXML`, `addMolV2000`, `addMolV3000`, `addInChI`,
  `addRXNV2000`) when the document exposes that method

## Things to verify inside ChemDraw

* the add-in appears in the Add-ins menu and opens non-modal at 560×720
* CORS: the page runs from `file://` (origin `null`); the bridge allows that and answers the `OPTIONS`
  preflight caused by the `X-Bridge-Token` header — the status dot must turn green within 5 s
* `EventSource` streaming through the web view (deltas appear live; keep-alive comments are ignored)
* `localStorage` persists across reopenings of the panel (file:// origins are sometimes sandboxed; the
  page falls back to in-memory state without errors)
* the selection chip updates when the selection changes and `selection.getSMILES()` throws harmlessly on
  partial structures (CDXML is sent instead)
