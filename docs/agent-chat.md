# Agent Chat inside ChemDraw

`chemdrawmcp` ships an in-app chat panel for ChemDraw modelled on the Agent Chat panels of
[mol-mcp](https://github.com/Jetzu/Structural-Biology-MCPs) for PyMOL and ChimeraX.

## Pieces

| piece | where | role |
|---|---|---|
| **Agent Chat add-in** | `chemdrawmcp/addin/` → installed to `~/Library/Application Support/com.revvity.ChemDraw/Add-ins/Agent Chat/` (macOS) or `…\ChemDraw Items\Add-ins\Agent Chat\` (Windows) | HTML/JS panel shown by ChemDraw (**Add-ins ▸ Agent Chat**). Claude-Code look: dark transcript, `>` prompts, `⏺` assistant text and tool calls, `⎿` result previews, status spinner, slash commands, model switcher, shortcuts footer. Uses the ChemDraw add-in JavaScript API (`ChemDrawAPI` 1.6) to attach the current selection as SMILES and to resize its window. |
| **Agent bridge** | `chemdrawmcp/bridge/`, console script `chemdraw-agent-bridge` | Local HTTP + SSE server on `127.0.0.1:47325` that runs `claude -p … --output-format stream-json` or `codex exec --json …` and streams the events to the panel. Keeps per-provider session ids and history in `~/.chemdrawmcp/agent_sessions.json`, model/effort settings in `agent_settings.json`. |
| **Scripts ▸ Start Agent Bridge** | `~/Library/Application Support/com.revvity.ChemDraw/Scripts/Start Agent Bridge.scpt` | AppleScript that starts the bridge from inside ChemDraw (macOS). |
| **launchd agent** (optional) | `~/Library/LaunchAgents/com.chemdrawmcp.bridge.plist` (`python3 install.py --launch-agent`) | keeps the bridge running at login. |

The CLIs already have the `chemdraw` MCP server registered by the installer, so the agent reads and edits the
document you are looking at through MCP tools (`read_document`, `draw_reaction`, `render`, `run_command`, …).

## Install / uninstall

```bash
python3 install.py                 # registers the MCP server and, when ChemDraw is found, installs add-in + Scripts item
python3 install.py --addin         # only the add-in and the Scripts item
python3 install.py --launch-agent  # plus the launchd agent
python3 install.py --uninstall     # removes everything (registrations, add-in, Scripts item, launchd agent)
```

Restart ChemDraw after installing: add-ins and scripts are read at launch.

## Using the panel

* Enter sends, Shift+Enter inserts a newline, Esc interrupts the running turn, ↑/↓ recall earlier prompts, Tab accepts a
  slash suggestion.
* `/model` – provider · model · effort chooser (also the button bottom right); `/claude`, `/codex` – switch provider;
  `/clear` – new conversation; `/stop` – interrupt; `/skills` – ask the agent to list its skills; `/help`.
* The header dot shows the bridge state. Red means the bridge is not running: use **Scripts ▸ Start Agent Bridge** or run
  `chemdraw-agent-bridge` (`--check` lists the CLIs it found).
* "selection attached" means the SMILES of your current selection was sent with the prompt.
* Transcripts persist per provider across reopenings; the bridge resumes the CLI session, so follow-up questions keep
  their context.

## Security

* The bridge binds to `127.0.0.1` only.  Every request needs the per-install token (`~/.chemdrawmcp/bridge_token`,
  mode 0600) as `X-Bridge-Token` header or `?token=` query; the installer bakes it into the add-in's `config.js`.
* Requests carrying an `http(s)` `Origin` other than localhost are refused, so a web page cannot drive the bridge even
  if it guessed the port.
* Claude runs with `--allowedTools mcp__chemdraw` (plus `Bash,Edit,Write` only when a run is started with
  `allow_shell`); Codex runs with `--approve-for-me`.

## Environment variables

| variable | meaning |
|---|---|
| `CHEMDRAWMCP_BRIDGE_PORT` | bridge port (default 47325) |
| `CHEMDRAWMCP_CLAUDE_EXE`, `CHEMDRAWMCP_CODEX_EXE` | explicit CLI executables |
| `CHEMDRAWMCP_BRIDGE_CWD` | working directory for the CLIs (default `~/.chemdrawmcp/workspace`) |
| `CHEMDRAWMCP_AGENT_SETTINGS` | settings file (default `~/.chemdrawmcp/agent_settings.json`) |
| `CHEMDRAWMCP_ADDIN_DIR` | where to install the add-in instead of ChemDraw's user add-in folder |

## Bridge API (for other front ends)

`GET /health`, `GET/POST /settings`, `POST /run {backend, prompt, context?, new?}` → `{run_id}`,
`GET /events/<run_id>?token=` (SSE: `session`, `delta`, `text`, `tool`, `tool_result`, `stats`, `error`, `done`, `end`),
`POST /cancel/<run_id>`, `GET /runs`, `GET /history`, `POST /clear {backend}`, `GET /sessions`.
