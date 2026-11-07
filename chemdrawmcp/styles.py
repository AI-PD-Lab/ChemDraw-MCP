"""Document styles extracted from the stationery files shipped with ChemDraw 23.1 (settings only).

All lengths are in points (1/72 inch).  ``bond_length`` is the fixed bond length; ``label_*`` apply to atom
labels and ``caption_*`` to captions / free text.  ``colortable`` lists RGB triplets (0-1); ``color`` and
``bgcolor`` index into it (foreground, background).  ``page`` is the page bounding box (x0 y0 x1 y1).
"""
from __future__ import annotations

from typing import Any, Dict

STYLES: Dict[str, Dict[str, Any]] = {
 "acs": {
  "name": "ACS Document 1996",
  "stationery": "ACS Document 1996.cds",
  "bond_length": 14.4,
  "line_width": 0.6,
  "bold_width": 2,
  "margin_width": 1.6,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 10,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "adv_synth_catal": {
  "name": "Adv. Synth. Catal. Document",
  "stationery": "Adv. Synth. Catal. Document.cds",
  "bond_length": 17,
  "line_width": 1,
  "bold_width": 2,
  "margin_width": 1.6,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   121.69,
   36.0,
   121.69,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   368.5,
   720.91
  ]
 },
 "dark": {
  "name": "Dark Mode",
  "stationery": "Dark Mode.cds",
  "bond_length": 14.4,
  "line_width": 1.61,
  "bold_width": 2.72,
  "margin_width": 1.6,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 10,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    0.1216,
    0.1059,
    0.1412
   ],
   [
    0.9961,
    1.0,
    1.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 3,
  "bgcolor": 2,
  "page": [
   0.0,
   0.0,
   1440.0,
   1080.0
  ]
 },
 "idraw": {
  "name": "I Draw Styles",
  "stationery": "I Draw Styles.cds",
  "bond_length": 30,
  "line_width": 1,
  "bold_width": 4,
  "margin_width": 2,
  "hash_spacing": 2.7,
  "bond_spacing": 12,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": True,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "j_chin_chem_soc": {
  "name": "J. Chin. Chem. Soc. Doc",
  "stationery": "J. Chin. Chem. Soc. Doc.cds",
  "bond_length": 18,
  "line_width": 1,
  "bold_width": 2.5,
  "margin_width": 2,
  "hash_spacing": 2.5,
  "bond_spacing": 20,
  "chain_angle": 120,
  "label_size": 12,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   675.0,
   1170.0
  ]
 },
 "j_het_chem": {
  "name": "J. Het. Chem",
  "stationery": "J. Het. Chem.cds",
  "bond_length": 10.8,
  "line_width": 0.6,
  "bold_width": 2.06,
  "margin_width": 1.08,
  "hash_spacing": 1.3,
  "bond_spacing": 16,
  "chain_angle": 120,
  "label_size": 7,
  "caption_size": 7,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Times New Roman",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "j_mol_mod_1col": {
  "name": "J. Mol. Mod. (1 Column)",
  "stationery": "J. Mol. Mod. (1 Column).cds",
  "bond_length": 14.4,
  "line_width": 0.6,
  "bold_width": 2,
  "margin_width": 1.6,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 10,
  "label_face": 96,
  "print_margins": [
   185.55,
   36.0,
   185.55,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Times New Roman",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   240.9,
   720.0
  ]
 },
 "j_mol_mod_2col": {
  "name": "J. Mol. Mod. (2 Column)",
  "stationery": "J. Mol. Mod. (2 Column).cds",
  "bond_length": 14.4,
  "line_width": 0.6,
  "bold_width": 2,
  "margin_width": 1.6,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 10,
  "label_face": 96,
  "print_margins": [
   53.7,
   36.0,
   53.7,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Times New Roman",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   504.6,
   720.0
  ]
 },
 "nature": {
  "name": "Nature Research Document",
  "stationery": "Nature Research Document.cds",
  "bond_length": 10.8,
  "line_width": 0.6,
  "bold_width": 1.56,
  "margin_width": 1.19,
  "hash_spacing": 1.7,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 6,
  "caption_size": 6,
  "label_face": 96,
  "caption_face": 96,
  "print_margins": [
   50.88,
   55.84,
   50.88,
   55.84
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    0.0745,
    0.0745,
    0.0745
   ],
   [
    0.7294,
    0.051,
    0.1255
   ],
   [
    0.9608,
    0.749,
    0.3686
   ],
   [
    0.0,
    0.4431,
    0.2039
   ],
   [
    0.0,
    0.4549,
    0.5373
   ],
   [
    0.0039,
    0.1647,
    0.3451
   ],
   [
    0.3451,
    0.298,
    0.5098
   ],
   [
    0.4,
    0.4,
    1.0
   ],
   [
    0.7961,
    0.902,
    1.0
   ],
   [
    0.9333,
    0.8667,
    0.9333
   ],
   [
    1.0,
    0.0,
    1.0
   ],
   [
    0.2431,
    0.5255,
    0.4039
   ],
   [
    0.7529,
    0.7529,
    0.7529
   ],
   [
    0.9216,
    0.9216,
    0.9216
   ],
   [
    0.502,
    0.0,
    0.502
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   510.24,
   680.31
  ]
 },
 "a4": {
  "name": "New A4 Document",
  "stationery": "New A4 Document.cds",
  "bond_length": 30,
  "line_width": 1,
  "bold_width": 4,
  "margin_width": 2,
  "hash_spacing": 2.7,
  "bond_spacing": 12,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   523.28,
   769.89
  ]
 },
 "default": {
  "name": "New Document",
  "stationery": "New Document.cds",
  "bond_length": 30,
  "line_width": 1,
  "bold_width": 4,
  "margin_width": 2,
  "hash_spacing": 2.7,
  "bond_spacing": 12,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Helvetica",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "slide": {
  "name": "New Slide",
  "stationery": "New Slide.cds",
  "bond_length": 30,
  "line_width": 1.6,
  "bold_width": 4,
  "margin_width": 2,
  "hash_spacing": 2.7,
  "bond_spacing": 15,
  "chain_angle": 120,
  "label_size": 12,
  "caption_size": 16,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Helvetica",
  "colortable": [
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    1.0,
    1.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 1,
  "bgcolor": 2,
  "page": [
   0.0,
   0.0,
   1101.18,
   635.29
  ]
 },
 "phytomedicine": {
  "name": "Phytomedicine Document",
  "stationery": "Phytomedicine Document.cds",
  "bond_length": 20,
  "line_width": 1,
  "bold_width": 1.33,
  "margin_width": 1.25,
  "hash_spacing": 3,
  "bond_spacing": 8,
  "chain_angle": 120,
  "label_size": 12,
  "caption_size": 12,
  "label_face": 97,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "rsc_1col": {
  "name": "RSC (1 Column) Document",
  "stationery": "RSC (1 Column) Document.cds",
  "bond_length": 12.2,
  "line_width": 0.45,
  "bold_width": 1.6,
  "margin_width": 1.25,
  "hash_spacing": 1.75,
  "bond_spacing": 20,
  "chain_angle": 120,
  "label_size": 7,
  "caption_size": 7,
  "label_face": 96,
  "caption_face": 96,
  "print_margins": [
   162.0,
   36.0,
   162.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Helvetica",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   235.28,
   646.3
  ]
 },
 "rsc": {
  "name": "RSC (2 Column) Document",
  "stationery": "RSC (2 Column) Document.cds",
  "bond_length": 12.2,
  "line_width": 0.45,
  "bold_width": 1.6,
  "margin_width": 1.25,
  "hash_spacing": 1.75,
  "bond_spacing": 20,
  "chain_angle": 120,
  "label_size": 7,
  "caption_size": 7,
  "label_face": 96,
  "print_margins": [
   28.35,
   28.35,
   28.35,
   28.35
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Helvetica",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   538.58,
   785.2
  ]
 },
 "science_of_synthesis": {
  "name": "Science of Synthesis",
  "stationery": "Science of Synthesis.cds",
  "bond_length": 17,
  "line_width": 0.8,
  "bold_width": 2,
  "margin_width": 2.25,
  "hash_spacing": 2.5,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 10,
  "label_face": 96,
  "print_margins": [
   24.95,
   36.0,
   9.35,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0275,
    0.0078
   ],
   [
    1.0,
    0.4706,
    0.0
   ],
   [
    0.9882,
    1.0,
    0.0118
   ],
   [
    0.0118,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.9843
   ],
   [
    0.0078,
    0.0235,
    1.0
   ],
   [
    0.6314,
    0.0078,
    1.0
   ],
   [
    1.0,
    0.0078,
    0.651
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   560.97,
   769.89
  ]
 },
 "terminal_carbons": {
  "name": "Show Terminal Carbon Labels",
  "stationery": "Show Terminal Carbon Labels.cds",
  "bond_length": 30,
  "line_width": 1,
  "bold_width": 4,
  "margin_width": 2,
  "hash_spacing": 2.7,
  "bond_spacing": 12,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": True,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Times New Roman",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   540.0,
   720.0
  ]
 },
 "helv_chim_acta": {
  "name": "Verlag Helv. Chim. Acta Doc",
  "stationery": "Verlag Helv. Chim. Acta Doc.cds",
  "bond_length": 17,
  "line_width": 0.6,
  "bold_width": 2,
  "margin_width": 2,
  "hash_spacing": 2,
  "bond_spacing": 15,
  "chain_angle": 120,
  "label_size": 10,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   59.39,
   60.55,
   59.39,
   60.55
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Helvetica",
  "caption_font": "Helvetica",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   510.24,
   720.91
  ]
 },
 "wiley": {
  "name": "Wiley Document",
  "stationery": "Wiley Document.cds",
  "bond_length": 17,
  "line_width": 0.75,
  "bold_width": 2.6,
  "margin_width": 2,
  "hash_spacing": 2.6,
  "bond_spacing": 18,
  "chain_angle": 120,
  "label_size": 12,
  "caption_size": 12,
  "label_face": 96,
  "print_margins": [
   36.0,
   36.0,
   36.0,
   36.0
  ],
  "show_terminal_carbon_labels": False,
  "hide_implicit_hydrogens": False,
  "label_justification": "Auto",
  "caption_justification": "Left",
  "fractional_widths": True,
  "label_font": "Arial",
  "caption_font": "Arial",
  "colortable": [
   [
    1.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    0.0
   ],
   [
    1.0,
    0.0,
    0.0
   ],
   [
    1.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    0.0
   ],
   [
    0.0,
    1.0,
    1.0
   ],
   [
    0.0,
    0.0,
    1.0
   ],
   [
    1.0,
    0.0,
    1.0
   ]
  ],
  "color": 0,
  "bgcolor": 1,
  "page": [
   0.0,
   0.0,
   769.89,
   523.28
  ]
 }
}

ALIASES = {
    "acs1996": "acs", "acs_1996": "acs", "jacs": "acs", "new": "default", "chemdraw": "default",
    "rsc_2col": "rsc", "nature_research": "nature", "science": "science_of_synthesis",
    "helvetica_chimica_acta": "helv_chim_acta", "hca": "helv_chim_acta", "advanced_synthesis": "adv_synth_catal",
    "dark_mode": "dark", "a4_document": "a4", "presentation": "slide",
}

DEFAULT_STYLE = "acs"


def style_names():
    return sorted(STYLES)


def get_style(name=None) -> Dict[str, Any]:
    """Return a copy of a style by key, alias or (case-insensitive) stationery title.  Unknown -> KeyError."""
    if not name:
        name = DEFAULT_STYLE
    key = str(name).strip().lower().replace(" ", "_").replace("-", "_")
    key = ALIASES.get(key, key)
    if key in STYLES:
        return dict(STYLES[key])
    low = str(name).strip().lower()
    for k, s in STYLES.items():
        if s["name"].lower() == low or s["stationery"].lower() == low or s["name"].lower().startswith(low):
            return dict(s)
    raise KeyError("unknown style %r; known: %s" % (name, ", ".join(style_names())))
