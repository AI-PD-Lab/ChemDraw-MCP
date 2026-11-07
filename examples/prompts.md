# Prompts for Claude Code / Codex with the `chemdraw` MCP server

Twenty-four prompts that exercise the server the way a chemist would ask.  Paste them as written; the tools an agent
typically reaches for are noted after each one.  With ChemDraw installed the drawings open as documents and are
rendered by ChemDraw; without it (`CHEMDRAWMCP_BACKEND=none`, Linux, CI) the same prompts produce CDXML files
and RDKit renders, except for the ChemDraw-only tasks (Struct=Name, NMR, clean-up, cdx/eps/pdf export).

## Drawing structures

1. "Draw aspirin, ibuprofen and naproxen side by side in ACS style, label each with its name and molecular
   weight, and show me the render."
   `draw_structures(structures, labels, style="acs")`, `properties`, `render`

2. "Make a 4-column grid of these SMILES in Nature style with the compound numbers 1-8 underneath, then save it
   as `compounds.cdxml` and as a 300 dpi PNG on my Desktop."
   `draw_structures(columns=4, style="nature")`, `save_document(format="cdxml")`, `render(format="png")`

3. "Draw (S)-naproxen with the stereocentre shown as a wedge and colour the carboxylic acid red.  Then tell me
   what is actually in the document."
   `draw_structures(colors=...)` (or `insert` + `select`), `read_document`

4. "Insert this molfile into the open ChemDraw document next to what is already there, and don't start a new
   document."  (paste the molblock)
   `documents`, `insert(data, kind="mol", document=...)`, `screenshot`

## Reactions and schemes

5. "Draw the Fischer esterification of acetic acid with ethanol, H2SO4 (cat.) above the arrow, 'reflux, 4 h'
   below, equilibrium arrow, 85% yield under the ester."
   `draw_reaction(reaction="acetic acid + ethanol --[H2SO4 (cat.)]--> ethyl acetate + water", conditions_below, product_labels, arrow="equilibrium")`

6. "Draw this reaction SMILES with atom mapping and number the reaction with ChemDraw's auto-numbering:
   `[CH3:1][C:2](=[O:3])[OH:4].[OH:5][CH2:6][CH3:7]>>[CH3:1][C:2](=[O:3])[O:5][CH2:6][CH3:7].[OH2:4]`"
   `draw_reaction(atom_map=true)`, `run_command("autoNumberReaction")`

7. "Lay out a five-step synthesis of ibuprofen from isobutylbenzene as a scheme titled 'Scheme 2', wrapped to
   a two-column RSC width, with reagents above and conditions below every arrow and yields under the products."
   `draw_scheme(steps=[...], title="Scheme 2", style="rsc", max_width=...)`, `read_document`

8. "Show a retrosynthetic analysis of paracetamol: paracetamol => 4-aminophenol + acetic anhydride, with
   retrosynthetic arrows, then export it as SVG."
   `draw_reaction(arrow="retrosynthetic")` (or `draw_scheme`), `render(format="svg")`

## Naming

9. "What is the IUPAC name of `CC(C)Cc1ccc(cc1)C(C)C(=O)O`?  Also give me the name ChemDraw generates for the
   molecule currently selected in my document."
   `structure_to_name(structure=...)`, `structure_to_name(document=None)` (ChemDraw Struct=Name)

10. "Turn these names into structures and draw them: 2-acetoxybenzoic acid, (2S)-2-amino-3-phenylpropanoic
    acid, 1,3,7-trimethylxanthine.  Flag any name you could not interpret."
    `name_to_structure` (offline dictionary, then OPSIN, then ChemDraw's Name=Struct, then PubChem with `CHEMDRAWMCP_ONLINE=1`), `draw_structures`

## Properties and analysis

11. "Give me formula, exact mass, elemental analysis, logP, TPSA, H-bond donors/acceptors and Lipinski
    violations for sildenafil, and compare ChemDraw's values with RDKit's."
    `properties(structure="sildenafil", use_chemdraw=true)`

12. "Is `C6H6 + Br2 -> C6H5Br` balanced?  Fix it, draw the balanced equation and report the atom economy."
    `analyze_reaction`, `draw_reaction`

13. "For every molecule in the front document, list SMILES, formula, MW and the number of stereocentres, and
    tell me which ones have unassigned stereocentres."
    `read_document`, `properties`, `check_structure`

## Conversions and batch export

14. "Convert this InChI to a V3000 molfile and to canonical SMILES: `InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3`"
    `convert(data, to="molv3000")`, `convert(data, to="smiles")`

15. "Open every .cdxml in ~/Projects/paper/figures, export each one as PDF and as a 2000 px wide PNG into
    ~/Projects/paper/export, and give me a table of file names and image sizes."
    `open_document`, `render(format="pdf")`, `render(format="png", width=2000)`, `close_document` in a loop

16. "Take the SDF at ~/data/hits.sdf, draw the first twelve compounds as a labelled grid in the slide style and
    save it as `hits.cdx` for the group meeting."
    `convert`/`draw_structures` from the SDF, `save_document(format="cdx")`

## NMR

17. "Predict the 1H and 13C NMR spectra of ethyl acetate, keep the spectrum windows open and show me the PNGs;
    then summarise the predicted shifts in a table."
    `predict_nmr(structure="ethyl acetate", nucleus="1H")`, `predict_nmr(nucleus="13C")`

## Clean-up and checking

18. "The structure I just sketched in ChemDraw looks crooked.  Clean it up, check it for chemical warnings and
    tell me what ChemDraw complains about."
    `clean_structure`, `check_structure`, `screenshot`

19. "Change the front document to ACS 1996 settings (bond length 14.4 pt, Arial 10), then clean the reaction
    so the arrows and components line up."
    `document_settings(settings={...})`, `clean_structure(reaction=true)`

## Document control and raw commands

20. "Which ChemDraw commands have to do with atom numbering?  Turn atom numbers on for the current document,
    flip the molecule horizontally, take a screenshot, and then save the document as `numbered.cdxml`."
    `commands(pattern="number")`, `run_command("toggleShowingAtomNumber")`, `run_command("flipObjectsHorizontal")`, `screenshot`, `save_document`

## Mechanisms, named reactions, files on disk

21. "Draw the SN2 mechanism of hydroxide attacking chloromethane with curved electron-pushing arrows, then the
    products."
    `draw_mechanism(stages=[{"structures": ["[OH-]", "CCl"], "arrows": ["lp:0:0->0:0-1:0", "1:0-1:1->1:1"]}, {"structures": ["CO", "[Cl-]"]}], title="SN2")`

22. "Suzuki coupling of 4-bromoanisole with phenylboronic acid: predict the product, draw it with the usual
    conditions and give me the atom-mapped reaction SMILES."
    `apply_named_reaction(reaction="suzuki", reactants=["COc1ccc(Br)cc1", "OB(O)c1ccccc1"])`, `named_reactions(name="suzuki")`

23. "What named reactions could benzaldehyde and nitromethane undergo together?"
    `suggest_reactions(reactants=["O=Cc1ccccc1", "C[N+](=O)[O-]"])`, then `apply_named_reaction` for the one the user picks

24. "Read `~/Desktop/scheme.cdx` without opening it and list every molecule, caption and reaction step."
    `read_file(path="~/Desktop/scheme.cdx")` (binary CDX and CDXML, no ChemDraw needed)

### Tips for agents

* Start with `status` (is ChemDraw running? which documents are open? where do files go?) and prefer the window
  id it returns over document names, which ChemDraw changes after the first edit.
* `read_document` before editing an existing document, `render` or `screenshot` after, `dismiss_dialogs` if a
  tool reports a pending dialog.
* Reaction text accepts names and SMILES mixed: `"A + B --[reagents]--> C + D"`; the bracket goes above the arrow.
