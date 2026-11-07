// Agent Chat for ChemDraw - bridge connection settings.
//
// The installer (`python3 install.py --addin`, or `python -m chemdrawmcp.addin_builder`)
// writes the real `config.js` next to main.html with the port and the secret token of the
// local agent bridge (`chemdraw-agent-bridge`).  Copy this file to `config.js` and fill in
// the token only when installing by hand; without config.js the add-in assumes port 47325
// and no token, and the bridge answers 401 ("bridge not configured").
window.CHEMDRAWMCP_BRIDGE = {
  port: 47325,
  token: "paste-the-bridge-token-here"
};
