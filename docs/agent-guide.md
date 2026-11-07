# Agent guide

How an AI agent (Claude Code, Codex, Claude Desktop or any MCP client) should work with the `chemdraw`
server.  The tool reference is in [tools.md](tools.md); this page is about *sequencing* the tools.

## 1. Recommended workflow

```
status ──► draw_* / insert / new_document ──► read_document ──► render / screenshot ──► fix ──► save_document
  │                                                                                   ▲
  └── headless? (chemdraw_available=false) ──► draw_* still write CDXML + RDKit PNG ───┘
```

1. **`status` first.**  It tells you whether ChemDraw is live (`chemdraw_available`, `running`), which
   documents are open (`documents` with `window_id`s), whether a dialog is blocking the app
   (`pending_dialogs`), where files go (`renders_dir`) and which styles exist.  Do not call `launch`
   explicitly unless the user asked; every live tool starts ChemDraw on demand.
2. **Draw with the high-level tools** (`draw_structures`, `draw_reaction`, `draw_scheme`, `draw_mechanism`,
   `apply_named_reaction`).  They lay the page out for you, open a new document and come back with a thumbnail and a numerical analysis.  Use
   `insert` / `add_caption` / `add_arrow` only to add to an existing document, and `run_command` for
   ChemDraw operations the high-level tools do not cover (atom mapping, autonumbering, flipping, zoom).
3. **Keep the `window_id`** from the response and pass it as `document` from then on.  Names of untitled
   documents change after the first edit; window ids do not.
4. **Verify** (section 4) before telling the user it is done.
5. **Export** with `save_document` (SVG/PDF/EPS for manuscripts, CDXML to keep editing, PNG for chat) and
   give the user the path.  If they want to see it in ChemDraw, `activate_chemdraw`.

Behavioural rules:

* The user's ChemDraw is a live application with their own documents in it.  Do not close documents you did
  not open (`close_document` discards changes by default); do not `quit_chemdraw` unless asked.
* Prefer **new documents** for your drawings.  Pasting into the user's document (`document=...`) lets
  ChemDraw choose the position; the pasted objects stay selected so the user can move them, and you must
  say so.
* Commands act on the selection of the key window.  `insert` leaves its objects selected; otherwise call
  `select` first.
* A `.cdxml` or `.cdx` file the user points you to can be inspected with `read_file(path)` without opening
  it (no ChemDraw needed); `open_document` it only when you need to render or edit it.
* ChemDraw stays in the background while tools run (the backend never activates it); only `activate_chemdraw`
  brings it forward.  Do not activate it yourself during long operations: ChemDraw hotkeys relabel selected atoms.
* Set `return_image=false` when you only need the data (for example a batch of structures you will not
  inspect) — thumbnails cost context.
* Tool results are JSON; a failure is `{"ok": false, "error": ..., "error_type": ...}`.  Check it.

## 2. Picking a style

`style` on the drawing tools reproduces the settings of a ChemDraw stationery (bond length, fonts, sizes,
line widths, page size); `new_document(style)` opens the real stationery.  Call `styles` for the table.

| the user says | use |
|---------------|-----|
| nothing, "a quick drawing", "for the group meeting" | `acs` (default; ACS Document 1996, the most common journal style) |
| JACS, Org. Lett., J. Org. Chem., any ACS journal, "ACS style" | `acs` |
| Nature, Nature Chemistry, Nat. Commun. | `nature` (small: 10.8 pt bonds, 6 pt Arial) |
| Chem. Commun., Chem. Sci., Org. Biomol. Chem., any RSC journal | `rsc` (two-column) or `rsc_1col` |
| Angew. Chem., Chem. Eur. J., Adv. Synth. Catal., any Wiley journal | `wiley` (large: 17 pt bonds, 12 pt Arial); `adv_synth_catal` for that journal |
| Science of Synthesis / Thieme | `science_of_synthesis` |
| Helv. Chim. Acta | `helv_chim_acta` |
| a slide, a poster, a presentation | `slide` (white on dark, 30 pt bonds) or `default` (ChemDraw's New Document, 30 pt bonds) |
| dark mode | `dark` |
| "show the methyl groups", "show terminal carbons as CH3" | `terminal_carbons` |
| a specific bond length ("1.5 cm bonds") | any style plus `bond_length` in points (1 cm = 28.35 pt) |

Everything inside one document should use one style.  If the user already has a document open in a
journal style and wants additions, read `document_settings` and pass the same `style`/`bond_length` so the
pasted drawing matches.

## 3. Describing reactions

Use the text syntax for anything simple; use component lists when you need labels, colours or control over
the text lines.

```
A + B --[reagents]--> C              bracket text goes ABOVE the arrow
A + B -> C + D                       also  -->  >>  =>  →  ⟶  ⇌  ⇄  and  " = "
CC(=O)O.OCC>[H+]>CC(=O)OCC.O         reaction SMILES; agents are drawn above the arrow
2 H2O -> 2 H2 + O2                   leading integers are stoichiometric coefficients
```

* Components are SMILES or names; names must resolve (the built-in dictionary first, then OPSIN if installed,
  then ChemDraw's Name=Struct when live, then PubChem only with `CHEMDRAWMCP_ONLINE=1`).  When a name fails, give
  the SMILES yourself rather than guessing a synonym.
* Conditions that are not reagents (temperature, time, solvent, yield) go in `conditions_below`; reagents in
  `reagents_above` or the bracket.  Formula-looking text (`NaBH4`, `H3O+`, `Pd(PPh3)4`) automatically gets
  subscripts and superscripts.
* Yields and compound numbers are captions under the components: `product_labels: ["85%", null]`,
  `reactant_labels: ["1", "2"]`.
* Arrow types: `equilibrium` for reversible steps, `retrosynthetic` for disconnections, `nogo` for a reaction
  that does not work, `dashed`/`bold`/`hollow`/`double`/`resonance` as needed.
* One arrow per `draw_reaction`.  A sequence is a `draw_scheme` with one step per arrow; write each step as
  its own `{"reaction": ...}` or `{"reactants": ..., "products": ...}` object.  Consecutive steps that start
  from the previous product are chained automatically (the intermediate is drawn once); the scheme wraps into
  rows at the page width (`max_width` to change it) and `title` puts a bold caption above.
* Atom mapping: give map numbers in the SMILES (`[CH3:1]C(=O)[OH:2]`) and set `atom_map=true`, or draw the
  reaction and then `run_command("mapReactionAtoms")` to let ChemDraw map it.
* To hand ChemDraw a reaction it recognises as such (for RXN export, *Autonumber Reaction*, stoichiometry
  analysis) you need the `<scheme>/<step>` record that `draw_reaction`/`draw_scheme` write.  A reaction the
  user drew by hand may lack it; `read_document` shows `steps: []` in that case.

**Named reactions.**  When the user names a reaction and gives the starting materials but not the product
("Suzuki coupling of bromobenzene and phenylboronic acid"), let the template library do the chemistry:
`apply_named_reaction(reaction, reactants)` predicts the product(s), draws the atom-mapped reaction with the
template's standard reagents/conditions and returns the mapped reaction SMILES; `named_reactions(name=...)`
shows what a template expects (roles, reagents) and `suggest_reactions(reactants)` lists what a set of
starting materials can do.  Mention the template's reagents to the user and replace them with
`draw_reaction(..., reagents_above=[...])` when the user specified different conditions.

**Mechanisms** (curved electron-pushing arrows) are a separate tool, `draw_mechanism`: one stage per
snapshot of the mechanism, each with its structures and arrow specs.  Write the SMILES yourself and number
the atoms in your head in SMILES order, 0-based (`[OH-]` is atom 0; in `CCl` the carbon is 0 and chlorine 1);
with several structures in a stage prefix the structure index (`1:0`).  `lp:` for a lone pair, `h:` for a
C–H bond, `a-b` / `a=b` for bonds, `->` to an atom, a new bond or an existing bond, `:1e` for fishhooks.
Put protons you need to move as explicit `[H]` in the SMILES.  Use `"arrow": "resonance"` between resonance
structures and `"equilibrium"` for reversible steps.  Verify with the thumbnail and with
`mechanism.stages[*].arrows[*]` (`source`/`target` say what the server understood).

## 4. Verifying output

Never report success from the drawing call alone.  Use at least one of:

| check | tool | what to look at |
|-------|------|-----------------|
| **Chemistry is right** | `read_document(window_id)` | `molecules[*].smiles` / `formula` match what was asked; `steps[*]` reactants/products/above/below are the intended ones; `reaction_smiles`; captions in `texts` |
| **Picture is right** | the `render.analysis` of the drawing call, or `render(document)` | `summary`: content not clipped (`clipped_edges` empty), reasonable coverage, contrast; the thumbnail itself for overlaps, wrong wedges, missing charges |
| **Layout is clean** | `layout_warnings` in the drawing result | absent when nothing overlaps; otherwise it lists the colliding objects — redraw with fewer `columns`, a larger `max_width` or shorter text |
| **What the user sees** | `screenshot(document)` | magnification, selection handles, a dialog in the way, where ChemDraw pasted something |
| **A change took effect** | `compare_renders(before.png, after.png)` | `changed_fraction` > 0 and the `bbox` of the change where you expected it |
| **ChemDraw is happy** | `check_structure(document)` | `warnings` empty, no alert text |
| **A command did something** | the `before`/`after` counts and `alerts` in `run_command` results | `alerts` non-empty means ChemDraw complained |

Workflow for an edit: `render` -> change -> `render` -> `compare_renders`; or `read_document` before and after
and compare `molecules`/`texts`.

Compare the SMILES semantically (canonicalise in your head or with `convert(to="smiles")`) rather than as
strings: ChemDraw and RDKit write different but equivalent SMILES.

## 5. Error handling

| symptom | meaning | what to do |
|---------|---------|------------|
| `error_type: "ChemDrawNotAvailable"` (+ `hint`) | no live ChemDraw on this machine (Linux, not installed, or `CHEMDRAWMCP_BACKEND=none`) | stay headless: `draw_*` still write CDXML and RDKit PNGs; `properties`, `convert` (RDKit formats), `analyze_reaction`, `name_to_structure` (offline) work; tell the user the CDXML path to open in ChemDraw later |
| `"ChemDraw is not running; call `launch` first"` | a no-autolaunch tool (`render`, `read_document`, `run_command`, ...) was called while the app is closed | call `launch` (or any drawing tool) first |
| `"ChemDraw did not answer within 60s; it showed an alert which was dismissed: ..."` | a modal dialog blocked the Apple event; it has already been dismissed | read the alert text, fix the cause, retry once |
| `status.pending_dialogs` non-empty | a dialog is open (maybe the user's) | `dismiss_dialogs`, then continue; mention it to the user |
| `"no document named ..."; open documents: [...]` | the name changed or the document was closed | use the listed names / `documents` and switch to `window_id` |
| `"command is not enabled right now (...)"` | the command needs a selection or a different document state | `select(mode="all")` or select the structure, then retry; `commands(pattern, enabled_only=true)` shows what is possible now |
| `"copyAsSMILES is not enabled: select one or more structures first"` | `selection(as_format=...)` with nothing chemical selected | `select` first |
| `"No reaction center was found"` on `save_document(format="rxn")` | the document has no reaction ChemDraw recognises | draw it with `draw_reaction` (writes the step record) or ask ChemDraw to map it (`run_command("mapReactionAtoms")`) |
| `"could not interpret 'X' as SMILES, InChI, molfile, CDXML or a chemical name"` | name resolution failed | provide SMILES; for IUPAC names ChemDraw's Name=Struct usually works when live |
| `"Multi-step reactions are not supported; give a single step"` | several arrows in one `reaction` string | `draw_scheme` |
| `insert` returns `warning: "no atoms were added"` | ChemDraw rejected the pasted text (wrong `kind`, invalid SMILES, unlicensed importer) | check `kind`, try `convert(..., to="mol")` and `insert(kind="mol")` |
| `structure_to_name` returns `ok: false` with `alerts` | Struct=Name refused or is not licensed | report the alert; fall back to the SMILES/InChIKey from `properties` |
| `predict_nmr` returns `"no spectrum document appeared"` | NMR prediction unavailable in this edition, or the structure was rejected | report; do not retry blindly |
| `"ChemDraw COM call ... failed"` / `"none of the candidate members ... worked"` (Windows) | the COM object model differs from the table in `backends/windows.py` | tell the user which member failed; see `docs/windows.md` |

General: errors are already cleaned up (dialogs dismissed, temporary documents closed).  Retry once after
fixing the cause; if the same error repeats, report the exact `error` text to the user instead of trying
variations.

## 6. Example prompts and the tool sequences that answer them

1. **"Draw aspirin, ibuprofen and paracetamol side by side with their names."**
   `draw_structures(structures=["aspirin", "ibuprofen", "paracetamol"], columns=3)` (labels default to the
   names) -> inspect thumbnail/`analysis` -> report `window_id` and the PNG path.

2. **"Draw the Fischer esterification of acetic acid with ethanol, H2SO4 cat., reflux, ACS style."**
   `draw_reaction(reaction="acetic acid + ethanol --[H2SO4 (cat.)]--> ethyl acetate + water", conditions_below=["reflux"], style="acs")`
   -> `read_document(window_id)` to confirm `reaction_smiles` -> done.

3. **"Make a three-step scheme: benzene -> bromobenzene (Br2/FeBr3) -> phenylmagnesium bromide (Mg, THF) -> benzoic acid (CO2 then H3O+), number the compounds."**
   `draw_scheme(steps=[{"reaction": "benzene --[Br2, FeBr3]--> bromobenzene", "product_labels": ["2"]}, {"reaction": "bromobenzene --[Mg, THF]--> Brc1ccccc1 ... ` — Grignard reagents are easier as SMILES:
   `{"reactants": ["Brc1ccccc1"], "products": ["Br[Mg]c1ccccc1"], "above": ["Mg, THF"], "product_labels": ["3"]}, {"reactants": ["Br[Mg]c1ccccc1"], "products": ["OC(=O)c1ccccc1"], "above": ["1. CO2", "2. H3O+"], "product_labels": ["4"]}]`, `title="Scheme 1"`)
   -> check `reaction_smiles` list has 3 entries and the render is not clipped.

4. **"What is the IUPAC name of this structure?" (user pastes a SMILES)**
   `structure_to_name(structure="<smiles>")` -> if `ok` false, `properties(structure=...)` for formula/InChIKey and explain.

5. **"Give me logP, TPSA and the Lipinski status for these five SMILES."**
   `properties(structure=s)` for each (headless-safe; consider `extended=true` default) -> tabulate `rdkit.logp`, `rdkit.tpsa`, `rdkit.lipinski_violations`.  No drawing needed.

6. **"Is this reaction balanced? What is the atom economy?"**
   `analyze_reaction(reaction=...)` -> report `balanced`, `unbalanced_elements`, `atom_economy`.

7. **"Convert this molfile to a ChemDraw file I can open."**
   `convert(data="<path or molblock>", to="cdx", path="~/Desktop/x.cdx")` (ChemDraw writes real CDX; headless use `to="cdxml"`).

8. **"Name this: 2-chloro-N-(4-methoxyphenyl)acetamide, draw it and export an SVG for the SI."**
   `name_to_structure(name=..., draw=true)` -> `render(document=<drawn.window_id>, format="svg", path="~/Desktop/compound.svg")`.

9. **"Clean up the structure I just drew and check it for errors."**
   `status` (find the front document / `window_id`) -> `clean_structure(document=wid)` -> `check_structure(document=wid)` -> `render(document=wid)` and show the thumbnail.

10. **"Add a yield of 92% under the product in my open reaction."**
    `read_document()` to find the product's `bounds` -> `add_caption(text="92%", document=wid)` -> tell the user the caption is pasted and selected so they can drop it under the product (ChemDraw places pasted objects), or `screenshot` to show where it landed.

11. **"Predict the 1H NMR of ethyl acetate."**
    `predict_nmr(structure="CC(=O)OCC", nucleus="1H")` -> relay the `text` (shift assignments) and the spectrum PNG; mention the spectrum document is open in ChemDraw (`keep_open=true`).

12. **"Map the atoms of the reaction in my current document and number the compounds."**
    `run_command("mapReactionAtoms", document=wid)` -> `run_command("autoNumberReaction", document=wid)` -> check `alerts` empty -> `read_document` shows `atom_map_pairs` > 0 -> `render`.

13. **"Put this peptide in: H-Ala-Gly-Ser-Phe-OH."**
    `insert(data="Ala-Gly-Ser-Phe", kind="peptide", document=wid)` -> `selection()` to read back SMILES/formula -> `render`.

14. **"Turn the ChemDraw file on my desktop into a PNG and tell me what is in it."**
    `read_file(path="~/Desktop/file.cdx")` for molecules/captions/steps (works headless, also for binary CDX) -> `open_document(path)` -> `render(document=window_id, format="png")` -> summarise; `close_document(window_id)` only if the user does not want it left open.

15. **"Show me the ACS vs Nature rendering of the same reaction."**
    `quick_render(input="CC(=O)O + OCC -> CC(=O)OCC + O", style="acs")` and `quick_render(..., style="nature")` (no documents left open) -> show both thumbnails and the bond length / font differences from `styles`.

16. **"Flip the molecule horizontally and zoom to fit."**
    `select(mode="all", document=wid)` -> `run_command("flipObjectsHorizontal", document=wid)` -> `run_command("zoomToFit", document=wid)` -> `screenshot(document=wid)` (zoom is a view property; `render` would not show it).

17. **"I'm on Linux, just give me a CDXML of these structures."**
    `draw_structures(structures=[...], open_in_chemdraw=false)` -> return `cdxml_path` (and the RDKit PNG).  Nothing else is needed; `status` will have said `chemdraw_available: false`.

18. **"Draw the SN2 mechanism of hydroxide attacking chloromethane with curved arrows."**
    `draw_mechanism(stages=[{"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"], "label": "SN2"}, {"structures": ["CO", "[Cl-]"]}], title="SN2")`
    -> check the thumbnail (arrow from the oxygen lone pair to the carbon, C–Cl bond onto Cl) and `mechanism.stages[0].arrows` -> `save_document(format="svg")` if asked.

19. **"Show the resonance structures of the acetate anion."**
    `draw_mechanism(stages=[{"structure": "CC(=O)[O-]", "arrows": ["lp:3->3-1", "1=2->2"], "arrow": "resonance"}, {"structure": "CC([O-])=O"}])`
    (atom 1 is the carbonyl carbon, 2 the carbonyl oxygen, 3 the anionic oxygen) -> thumbnail -> done.

20. **"Draw the Suzuki coupling of 4-bromoanisole with phenylboronic acid."**
    `apply_named_reaction(reaction="suzuki", reactants=["COc1ccc(Br)cc1", "OB(O)c1ccccc1"])` -> read `outcomes[0].products` and tell the user the predicted product and the template's reagents (`Pd(PPh3)4, K2CO3; dioxane/H2O, 80 °C`); the drawing is already open (`drawn.window_id`).

21. **"What could I do with benzaldehyde and nitromethane?"**
    `suggest_reactions(reactants=["O=Cc1ccccc1", "C[N+](=O)[O-]"])` -> list the matching templates with their products -> offer to draw one with `apply_named_reaction`.
