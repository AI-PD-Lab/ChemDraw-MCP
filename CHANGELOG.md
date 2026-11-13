# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Mechanism arrows: label boxes of explicitly labelled atoms now follow ChemDraw's hydrogens-first layout ("H3C"), so heads no longer land on the hydrogens; heads keep clear of charge superscripts (also for bond-forming arrows); bond-sourced tails start closer to the bond.
- README rewritten around a gallery of ChemDraw renders (mechanisms, reactions, schemes, NMR); reference material moved to `docs/detailed.md`; `examples/make_gallery.py` regenerates the images.

## 0.1.0 - 2026-11-07

First release of `chemdrawmcp`, an MCP server that gives Claude Code, Codex and any other MCP client access to
ChemDraw (Revvity Signals, formerly PerkinElmer/CambridgeSoft).

### Added
- **Agent Chat inside ChemDraw**: a ChemDraw add-in (`Add-ins ▸ Agent Chat`, Claude-Code-style transcript with
  `⏺` tool calls, `⎿` result previews, status spinner, `/model`, `/claude`, `/codex`, `/clear`, `/stop`, `/help`) backed by
  `chemdraw-agent-bridge`, a localhost HTTP/SSE server that runs the installed Claude Code / Codex CLI with the
  `chemdraw` MCP server; the panel attaches the current selection via the ChemDraw add-in JavaScript API.  Installed by
  `python3 install.py` (`--addin`, optional `--launch-agent`) together with a `Scripts ▸ Start Agent Bridge` item.
- macOS backend no longer activates ChemDraw (windows are raised inside the app); ChemDraw is launched hidden.

- **Live ChemDraw backend** for macOS, driving the application's AppleScript dictionary through JXA: documents,
  atoms/bonds/molecules, captions, arrows, the selection (SMILES, formula, MW, exact mass, elemental analysis),
  drawing settings, `open`/`save as` for every export format, menus and the ~1600 named commands. Data goes in and
  out through the general pasteboard via NSPasteboard (`pasteAsSMILES`, `pasteAsCdxml`, `copyAsCdxml`, ...); the
  user's clipboard is snapshotted with every item and flavour and restored inside the same script. While a script
  runs, ChemDraw is polled every 2 s for modal alerts, which are detected via Accessibility, dismissed and returned
  as the error text. `status` reads version, documents and the front document in one `snapshot()` round trip.
  ChemDraw is driven in the background: documents are raised inside ChemDraw without activating the app, so the
  user keeps keyboard focus (`CHEMDRAWMCP_ALLOW_ACTIVATE=1` permits activation as a last resort).
- **Windows backend** (best effort) driving `ChemDraw.Application` through COM (pywin32 or comtypes) on one STA
  worker thread with dialog polling, the candidate member names in editable tables at the top of
  `chemdrawmcp/backends/windows.py`, and an optional ChemScript bridge.
- **Headless CDXML engine** (`cdxml.py`, `layout.py`, `textmetrics.py`, `styles.py`): RDKit plus an in-house CDXML
  writer/reader that lays out molecules, reactions (reagents above the arrow, conditions below, plus signs, yields,
  labels, atom mapping) and wrapped multi-step schemes using the settings of ChemDraw's stationery styles (ACS 1996,
  Nature, RSC, Wiley, Science of Synthesis, ...), with label-aware fragment boxes measured from real fonts, overlap
  detection (`layout_warnings`), automatic page tiling, and the reader's nickname expansion (`abbreviations.py`).
  Works without ChemDraw (Linux, CI) via `CHEMDRAWMCP_BACKEND=none`.
- **Binary CDX reader** (`cdx.py`): parses `.cdx` files (verified against ChemDraw 23 documents) and converts them to
  CDXML, so `read_file` and `convert` handle `.cdx` without ChemDraw.
- **Mechanisms** (`mechanism.py`, `draw_mechanism`): curved electron-pushing arrows (`lp:5->5-8`, `3=4->4-7`,
  `h:3->3-4`, `:1e` fishhooks, `fragment:atom` addressing) planned around atom labels, stages joined by
  reaction/resonance arrows and wrapped into rows.
- **Named reactions** (`reactions.py`): 100 templates (Suzuki, Grignard, Wittig, amide coupling, Boc, CuAAC,
  Diels–Alder, aldol, oxidations, reductions, ...) with RDKit reaction SMARTS, typical reagents/conditions, roles,
  aliases, by-products and a worked example each; product prediction with atom-mapped reaction SMILES
  (`named_reactions`, `apply_named_reaction`, `suggest_reactions`).
- **Name resolution**: built-in offline dictionary (`chem.py`, ~3,100 structures under ~7,800 names and aliases:
  reagents, solvents, drugs, amino acids, sugars, cofactors, natural products) -> OPSIN if installed -> ChemDraw's
  Name=Struct -> PubChem (only with `CHEMDRAWMCP_ONLINE=1`); results report their `source`.
- **43 tools**: `status`, `launch`, `quit_chemdraw`, `activate_chemdraw`; `documents`, `new_document`, `open_document`,
  `close_document`, `save_document` (cdxml, cdx, mol / sdf / rxn as V2000, molv3000 / sdfv3000 / rxnv3000, rdf, cml,
  ct, png, svg, pdf, eps, tiff, gif, jpeg, bmp, 3mf), `read_document`, `read_file`, `document_settings`;
  `draw_structures`, `draw_reaction`, `draw_scheme`, `draw_mechanism`, `named_reactions`, `apply_named_reaction`,
  `suggest_reactions`, `insert`, `add_caption`, `add_arrow`; `render`, `quick_render`, `screenshot`, `analyze_image`,
  `compare_renders`; `name_to_structure`, `structure_to_name`, `properties`, `analyze_reaction`, `convert`,
  `clean_structure`, `check_structure`, `predict_nmr`; `selection`, `select`, `commands`, `run_command`, `menus`,
  `do_menu`; `dismiss_dialogs`, `styles`.
- Every render returns a numerical description (size, coverage, bounding box, margins, dominant colours) and a
  thumbnail; `read_document` / `read_file` turn a document back into SMILES / formulas / captions / reaction SMILES.
- **Styles**: `acs` (default), `nature`, `rsc`, `rsc_1col`, `wiley`, `science_of_synthesis`, `j_het_chem`,
  `j_chin_chem_soc`, `j_mol_mod_1col`, `j_mol_mod_2col`, `helv_chim_acta`, `adv_synth_catal`, `phytomedicine`,
  `default`, `a4`, `slide`, `dark`, `idraw`, `terminal_carbons`.
- **Installer** (`python install.py` / `chemdrawmcp-install`): creates `.venv`, installs the package (+ rdkit) and
  registers the `chemdraw` server with Claude Code (`~/.claude.json`), Codex CLI (`~/.codex/config.toml`) and
  Claude Desktop; `--check`, `--dry-run`, `--uninstall`. Config files are backed up as `*.chemdrawmcp-backup`.
- **Examples** (`examples/`): seven runnable scripts (library API, CDXML reader, properties, an MCP stdio client,
  a hand-built CDXML document) and a prompt collection for agents.
- Environment variables: `CHEMDRAW_APP` / `CHEMDRAW_EXE`, `CHEMDRAWMCP_HOME`, `CHEMDRAWMCP_OUTPUT_DIR`,
  `CHEMDRAWMCP_AUTOLAUNCH`, `CHEMDRAWMCP_BACKEND`, `CHEMDRAWMCP_STYLE`, `CHEMDRAWMCP_ALLOW_ACTIVATE`,
  `CHEMDRAWMCP_ONLINE`, `CHEMDRAWMCP_DEBUG`, `CHEMDRAWMCP_NO_FONTS`, `CHEMDRAWMCP_FONT_DIR`, `CHEMDRAWMCP_LIVE`,
  `CHEMDRAWMCP_CDX_PAIRS`, `CHEMDRAWMCP_SAMPLES`, `OPSIN_JAR`, `CHEMSCRIPT_PATH` / `CHEMSCRIPT_MODULE`.
- About 2,070 unit tests in twelve modules (server tools headless, CDXML writer/reader and robustness, binary CDX,
  layout quality, mechanisms, named reactions, chemistry, the name dictionary, installer, macOS backend against a
  fake object model, Windows backend against a fake COM application) that run without ChemDraw or network; 31
  `live` tests against a running ChemDraw are opt-in via `CHEMDRAWMCP_LIVE=1 pytest -m live`.
- Repository hygiene: PolyForm Noncommercial licence, GitHub Actions CI (pytest on Linux/macOS/Windows x Python 3.10/3.12, ruff lint),
  `.editorconfig`, this changelog.

### Known limits

- ChemDraw's PNG export is always 300 dpi; `render(width=...)` resamples. SVG/PDF/EPS are vector.
- Struct=Name, NMR prediction and some import formats depend on the ChemDraw edition/licence.
- The Windows backend was written from the published object model without a Windows machine at hand; it offers a
  curated table of 19 commands, and `menus`/`do_menu`, `check_structure` and `predict_nmr` are macOS-only.

[0.1.0]: https://github.com/Jetzu/ChemDraw-MCP/releases/tag/v0.1.0
