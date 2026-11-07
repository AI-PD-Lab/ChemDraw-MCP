"""Named-reaction templates for chemdrawmcp.

A registry of ~100 textbook / medicinal-chemistry reactions written as RDKit
reaction SMARTS, each with display name, class, typical reagents (text above
the arrow), conditions (text below the arrow), reactant role names, aliases
and a worked example.  The module can

* look templates up fuzzily (:func:`get_reaction`, :func:`list_reactions`),
* apply one to concrete reactants and return deduplicated products together
  with an atom-mapped reaction SMILES (:func:`apply_reaction`),
* tell which templates fit a set of reactants (:func:`suggest_reactions`),
* build a drawable :class:`chemdrawmcp.layout.Step` with the template's
  reagent / condition text (:func:`to_step`), and
* describe a template as text (:func:`describe`).

Conventions used in the SMARTS
------------------------------
* Only the drawn products are written in the product side; small
  by-products (H2O, HX, Ph3P=O, C2H4, N2 ...) are listed in
  ``ReactionTemplate.byproducts`` as text instead.
* Reactions that may open a ring (ester hydrolysis, ozonolysis ...) write
  their products as one grouped template ``(A.B)`` and the runner splits the
  result into fragments; this avoids RDKit duplicating ring atoms.
* Stoichiometry > 1 for a role (e.g. 2 equiv RMgX on an ester) is encoded by
  repeating the reactant template and declaring ``stoichiometry``; callers
  still pass each reactant once.
* Stereochemistry is only written where the reaction defines it (Lindlar
  → Z, Heck / HWE → E, SN2 / Mitsunobu → inversion).  Everything else is
  left unspecified.
* Regiochemistry: electrophilic aromatic substitutions return every
  positional isomer (deduplicated); Zaitsev / Markovnikov preferences are
  encoded where a simple SMARTS can express them.

RDKit is imported lazily (through :mod:`chemdrawmcp.chem`), so the registry
itself can be listed without RDKit; applying a template needs it.
"""
from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from difflib import get_close_matches
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Sequence, Tuple, Union

from .chem import ChemError, _quiet, _rd, mol_from_input, prepare_for_depiction

if TYPE_CHECKING:  # pragma: no cover - typing only

    from .layout import Step

__all__ = [
    "ReactionError",
    "ReactionTemplate",
    "REACTIONS",
    "classes",
    "list_reactions",
    "get_reaction",
    "apply_reaction",
    "suggest_reactions",
    "to_step",
    "describe",
]


class ReactionError(ChemError):
    """Raised when a template cannot be found or applied (a ``ValueError``)."""


# --------------------------------------------------------------------------- #
# Template dataclass
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ReactionTemplate:
    """One named reaction as an RDKit reaction SMARTS plus presentation data."""

    key: str                                 # snake_case identifier
    name: str                                # display name
    cls: str                                 # reaction class ("C–C coupling", "oxidation" ...)
    smarts: str                              # RDKit reaction SMARTS
    reagents: Tuple[str, ...]                # text above the arrow
    conditions: Tuple[str, ...]              # text below the arrow
    description: str                         # one line
    roles: Tuple[str, ...]                   # one name per reactant the caller passes
    aliases: Tuple[str, ...] = ()
    byproducts: Tuple[str, ...] = ()         # not drawn, text only
    stoichiometry: Tuple[int, ...] = ()      # equivalents per role (default 1 each)
    example_reactants: Tuple[str, ...] = ()  # SMILES, in role order
    example_products: Tuple[str, ...] = ()   # SMILES of the expected drawn products

    # -- derived ----------------------------------------------------------- #
    @property
    def n_reactants(self) -> int:
        """Number of reactants a caller passes (one per role)."""
        return len(self.roles)

    @property
    def equivalents(self) -> Tuple[int, ...]:
        return self.stoichiometry or tuple(1 for _ in self.roles)

    @property
    def rxn(self) -> Any:
        """The compiled ``rdChemReactions.ChemicalReaction`` (cached)."""
        return _compiled(self)

    @property
    def n_products(self) -> int:
        """Number of drawn product templates (after splitting grouped products)."""
        return max(1, self.smarts.split(">>", 1)[1].count(".") + 1)

    def role_template(self, role: int) -> Any:
        """Query molecule of the reactant template that fills ``role``."""
        offset = sum(self.equivalents[:role])
        return self.rxn.GetReactantTemplate(offset)

    def matches_role(self, mol: Any, role: int) -> bool:
        """True when ``mol`` contains the substructure required for ``role``."""
        with _quiet():
            return bool(mol.HasSubstructMatch(self.role_template(role)))

    def as_dict(self) -> Dict[str, Any]:
        return {
            "key": self.key,
            "name": self.name,
            "class": self.cls,
            "roles": list(self.roles),
            "n_reactants": self.n_reactants,
            "equivalents": list(self.equivalents),
            "reagents": list(self.reagents),
            "conditions": list(self.conditions),
            "byproducts": list(self.byproducts),
            "description": self.description,
            "aliases": list(self.aliases),
            "smarts": self.smarts,
            "example": {
                "reactants": list(self.example_reactants),
                "products": list(self.example_products),
            },
        }


_COMPILED: Dict[str, Any] = {}


def _compiled(tmpl: ReactionTemplate) -> Any:
    rxn = _COMPILED.get(tmpl.key)
    if rxn is None:
        rd = _rd()
        try:
            with _quiet():
                rxn = rd.AllChem.ReactionFromSmarts(tmpl.smarts)
                rxn.Initialize()
        except Exception as exc:  # pragma: no cover - templates are tested
            raise ReactionError(f"Template {tmpl.key!r} has invalid SMARTS: {exc}") from exc
        expected = sum(tmpl.equivalents)
        if rxn.GetNumReactantTemplates() != expected:  # pragma: no cover - tested
            raise ReactionError(
                f"Template {tmpl.key!r}: SMARTS has {rxn.GetNumReactantTemplates()} reactant "
                f"templates but roles/stoichiometry imply {expected}."
            )
        _COMPILED[tmpl.key] = rxn
    return rxn


# --------------------------------------------------------------------------- #
# SMARTS building blocks
# --------------------------------------------------------------------------- #
def _m(fragment: str, n: int) -> str:
    """Add atom map ``n`` to a single bracket-atom SMARTS fragment."""
    assert fragment.startswith("[") and fragment.endswith("]"), fragment
    return f"{fragment[:-1]}:{n}]"


# ammonia + 1°/2° amines (incl. anilines); no amides, sulfonamides, N–N, N–O
_AMINE = "[NX3;H3,H2,H1;+0;!$(N[C,S,P]=[O,S,N]);!$(N[#7,#8,#16])]"
_AMINE12 = "[NX3;H2,H1;+0;!$(N[C,S,P]=[O,S,N]);!$(N[#7,#8,#16])]"
_PRIM_AMINE = "[NX3;H2;+0;$(N[#6]);!$(N[C,S,P]=[O,S,N]);!$(N[#7,#8,#16])]"
_ACID = "[C;$(C(=O)[OH]);!$(C(=O)(O)[!#6])]"                 # carboxylic acid carbon
_KETALD = "[CX3;$(C=O);!$(C(=O)[!#6])]"                     # aldehyde or ketone carbon
_ALD = "[CX3;H1,H2;$(C=O);!$(C(=O)[!#6])]"                  # aldehyde carbon (incl. HCHO)
_KET = "[CX3;H0;$(C(=O)([#6])[#6])]"                        # ketone carbon
_ENOL_KET = "[C;$(C(=O)([#6])[#6]);!$(C(=O)[!#6])]"          # ketone carbon next to [CH2,CH3]
_ARYL = "[c,$(C=C)]"                                         # aryl or vinyl carbon
_ARLG = "[Cl,Br,I,$(OS(=O)(=O)C(F)(F)F)]"                   # halide / triflate on sp2 C
_LG = "[Cl,Br,I,$(OS(=O)(=O)[#6])]"                         # halide / sulfonate on sp3 C
_ALKYL = "[CX4;H3,H2,H1]"                                    # 1°/2° sp3 carbon (not tertiary)
_ALKYL_INV = "[C@;X4;H3,H2,H1]"                              # same, inverted in the product ([C@@:n])
_HYDROXYL = "[OX2H1;$(O[C,c]);!$(OC=O)]"                     # alcohol / phenol OH
_ALCOHOL_C = "[CX4;$(C[OH])]"
_AZIDE = "[#6:4][N;$(N=[N+]=[N-]),$([N-][N+]#N):1]~[N+:2]~[NX1:3]"   # R-N3 (both resonance forms)
_AZIDE_ION = "[N-:2]=[N+:3]=[N-:4]"
_ISOCYANIDE = "[C-:7]#[N+:6][#6:9]"
_MICHAEL_DONOR = "[C;H1,H2,H3;$(C(C=O)C=O),$(C(C=O)C#N),$(C(C#N)C#N),$(C[N+](=O)[O-])]"
_MICHAEL_ACCEPTOR = "[CX3:2]=[CX3:3][C,N,S;$(C=O),$(C#N),$([N+](=O)[O-]),$(S(=O)=O):4]"
_ACTIVE_CH2 = "[CH2;$(C(C=O)C=O),$(C(C=O)C#N),$(C(C#N)C#N)]"
_ANTI_MARK = "[C;$([CH2]=C),$([CH1]=[CH1,CH0]),$([CH0]=[CH0])]"   # less substituted alkene C
_EPOX_C1 = "[C;$([CH2]1OC1),$([CH1]1O[CH1,CH0]1),$([CH0]1O[CH0]1)]"  # less substituted epoxide C
_MITSUNOBU_NU = "[O,N,S;H1;$([OH]C=O),$([OH]c),$([SH]),$([NH](C=O)C=O),$([NH]S(=O)=O)]"


def _t(key: str, name: str, cls: str, smarts: str, reagents: Sequence[str], conditions: Sequence[str],
       description: str, roles: Sequence[str], example: Tuple[Sequence[str], Sequence[str]],
       aliases: Sequence[str] = (), byproducts: Sequence[str] = (),
       stoichiometry: Sequence[int] = ()) -> ReactionTemplate:
    return ReactionTemplate(
        key=key, name=name, cls=cls, smarts=smarts, reagents=tuple(reagents),
        conditions=tuple(conditions), description=description, roles=tuple(roles),
        aliases=tuple(aliases), byproducts=tuple(byproducts), stoichiometry=tuple(stoichiometry),
        example_reactants=tuple(example[0]), example_products=tuple(example[1]),
    )


# --------------------------------------------------------------------------- #
# The registry
# --------------------------------------------------------------------------- #
_TEMPLATES: List[ReactionTemplate] = [
    # ------------------------------------------------------------- acylation
    _t("fischer_esterification", "Fischer esterification", "acylation",
       f"{_m(_ACID, 1)}(=[O:2])[OH].[OH:3][CX4:4]>>[C:1](=[O:2])[O:3][C:4]",
       ["H2SO4 (cat.)"], ["ROH, reflux"],
       "Acid-catalysed condensation of a carboxylic acid with an alcohol to the ester.",
       ["carboxylic acid", "alcohol"], (["CC(=O)O", "CCO"], ["CCOC(C)=O"]),
       aliases=["fischer", "esterification", "fischer speier", "acid catalysed esterification"],
       byproducts=["H2O"]),
    _t("steglich_esterification", "Steglich esterification", "acylation",
       f"{_m(_ACID, 1)}(=[O:2])[OH].[OH:3][CX4:4]>>[C:1](=[O:2])[O:3][C:4]",
       ["DCC, DMAP (cat.)"], ["DCM, rt"],
       "Carbodiimide-mediated esterification; mild, works for hindered or acid-sensitive partners.",
       ["carboxylic acid", "alcohol"], (["OC(=O)c1ccccc1", "CC(C)(C)O"], ["CC(C)(C)OC(=O)c1ccccc1"]),
       aliases=["steglich", "dcc esterification", "edc esterification"], byproducts=["DCU"]),
    _t("acid_chloride_alcohol", "Acylation of an alcohol with an acid chloride", "acylation",
       f"[C:1](=[O:2])[Cl].{_m(_HYDROXYL, 3)}>>[C:1](=[O:2])[O:3]",
       ["pyridine or Et3N"], ["DCM, 0 °C → rt"],
       "Ester formation from an acyl chloride and an alcohol or phenol.",
       ["acid chloride", "alcohol/phenol"], (["CC(=O)Cl", "OCc1ccccc1"], ["CC(=O)OCc1ccccc1"]),
       aliases=["acyl chloride alcohol", "o-acylation", "ester from acid chloride"], byproducts=["HCl"]),
    _t("amide_coupling", "Amide coupling (acid + amine)", "acylation",
       f"{_m(_ACID, 1)}(=[O:2])[OH].{_m(_AMINE, 3)}>>[C:1](=[O:2])[N:3]",
       ["EDC·HCl, HOBt", "DIPEA"], ["DMF, rt"],
       "Coupling-reagent mediated amide formation from a carboxylic acid and an amine (EDC/HOBt, HATU, T3P ...).",
       ["carboxylic acid", "amine"], (["OC(=O)c1ccccc1", "CN"], ["CNC(=O)c1ccccc1"]),
       aliases=["amidation", "peptide coupling", "hatu coupling", "edc coupling", "amide formation",
                "amide bond formation"],
       byproducts=["H2O"]),
    _t("acid_chloride_amine", "Acid chloride + amine (Schotten–Baumann)", "acylation",
       f"[C:1](=[O:2])[Cl].{_m(_AMINE, 3)}>>[C:1](=[O:2])[N:3]",
       ["Et3N or aq. NaOH"], ["DCM, 0 °C → rt"],
       "Amide formation from an acyl chloride and an amine.",
       ["acid chloride", "amine"], (["CC(=O)Cl", "Nc1ccccc1"], ["CC(=O)Nc1ccccc1"]),
       aliases=["schotten baumann", "acyl chloride amine", "n-acylation", "acylation of amine"],
       byproducts=["HCl"]),
    _t("acid_chloride_formation", "Acid chloride formation", "acylation",
       f"{_m(_ACID, 1)}(=[O:2])[OH]>>[C:1](=[O:2])Cl",
       ["SOCl2 or (COCl)2, DMF (cat.)"], ["DCM, 0 °C → rt"],
       "Conversion of a carboxylic acid into the acyl chloride.",
       ["carboxylic acid"], (["OC(=O)c1ccccc1"], ["O=C(Cl)c1ccccc1"]),
       aliases=["thionyl chloride", "acyl chloride formation", "oxalyl chloride"], byproducts=["SO2", "HCl"]),
    _t("sulfonamide_formation", "Sulfonamide formation", "acylation",
       f"[SX4:1](=[O:2])(=[O:3])[Cl].{_m(_AMINE, 4)}>>[S:1](=[O:2])(=[O:3])[N:4]",
       ["pyridine or Et3N"], ["DCM, 0 °C → rt"],
       "Sulfonyl chloride + amine → sulfonamide.",
       ["sulfonyl chloride", "amine"], (["Cc1ccc(S(=O)(=O)Cl)cc1", "Nc1ccccc1"],
                                        ["Cc1ccc(S(=O)(=O)Nc2ccccc2)cc1"]),
       aliases=["sulfonylation", "sulfonyl chloride amine", "tosylation of amine"], byproducts=["HCl"]),
    _t("urea_formation", "Urea formation (isocyanate + amine)", "acylation",
       f"[N:1]=[C:2]=[O:3].{_m(_AMINE, 4)}>>[N:1]-[C:2](=[O:3])[N:4]",
       [], ["DCM, rt"],
       "Addition of an amine to an isocyanate giving a urea.",
       ["isocyanate", "amine"], (["O=C=Nc1ccccc1", "CN"], ["CNC(=O)Nc1ccccc1"]),
       aliases=["isocyanate amine", "urea synthesis"]),
    # ----------------------------------------------------- C–C cross-coupling
    _t("suzuki_coupling", "Suzuki–Miyaura coupling", "C–C coupling",
       f"{_m(_ARYL, 1)}{_ARLG}.{_m(_ARYL, 2)}[#5;$(B(O)O),$([B-](F)(F)F)]>>[*:1]-[*:2]",
       ["Pd(PPh3)4 (5 mol%)", "K2CO3"], ["dioxane/H2O, 80 °C"],
       "Pd-catalysed cross-coupling of an aryl/vinyl halide or triflate with a boronic acid, boronate ester or trifluoroborate.",
       ["aryl/vinyl halide", "boronic acid/ester"], (["Brc1ccccc1", "OB(O)c1ccccc1"], ["c1ccc(-c2ccccc2)cc1"]),
       aliases=["suzuki", "suzuki miyaura", "boronic acid coupling", "sm coupling"],
       byproducts=["B(OH)2X"]),
    _t("heck_reaction", "Heck reaction", "C–C coupling",
       f"{_m(_ARYL, 1)}{_ARLG}.[CH2:2]=[CH:3][*:4]>>[*:1]/[CH:2]=[CH:3]/[*:4]",
       ["Pd(OAc)2, PPh3", "Et3N"], ["DMF, 100 °C"],
       "Pd-catalysed arylation of a terminal alkene; the aryl group ends up on the terminal carbon, (E)-selective.",
       ["aryl/vinyl halide", "terminal alkene"], (["Brc1ccccc1", "C=CC(=O)OC"], ["COC(=O)/C=C/c1ccccc1"]),
       aliases=["heck", "mizoroki heck", "heck coupling"], byproducts=["HX"]),
    _t("sonogashira_coupling", "Sonogashira coupling", "C–C coupling",
       f"{_m(_ARYL, 1)}{_ARLG}.[CH:2]#[C:3]>>[*:1][C:2]#[C:3]",
       ["Pd(PPh3)2Cl2, CuI", "Et3N"], ["THF, rt"],
       "Pd/Cu-catalysed coupling of an aryl/vinyl halide with a terminal alkyne.",
       ["aryl/vinyl halide", "terminal alkyne"], (["Ic1ccccc1", "C#Cc1ccccc1"], ["c1ccc(C#Cc2ccccc2)cc1"]),
       aliases=["sonogashira", "alkyne coupling"], byproducts=["HX"]),
    _t("negishi_coupling", "Negishi coupling", "C–C coupling",
       f"{_m(_ARYL, 1)}{_ARLG}.[#6:2][Zn][Cl,Br,I,#6]>>[*:1][#6:2]",
       ["Pd(PPh3)4 (5 mol%)"], ["THF, rt → 60 °C"],
       "Pd-catalysed coupling of an aryl/vinyl halide with an organozinc reagent (RZnX or R2Zn).",
       ["aryl/vinyl halide", "organozinc reagent"], (["Brc1ccccc1", "CC[Zn]Cl"], ["CCc1ccccc1"]),
       aliases=["negishi", "organozinc coupling"], byproducts=["ZnX2"]),
    _t("stille_coupling", "Stille coupling", "C–C coupling",
       f"{_m(_ARYL, 1)}{_ARLG}.[#6;$([#6]=*),$([#6]#*),c:2][Sn]([#6])([#6])[#6]>>[*:1][#6:2]",
       ["Pd(PPh3)4 (5 mol%)"], ["toluene, 110 °C"],
       "Pd-catalysed coupling with an organostannane; the sp2/sp carbon on tin is transferred.",
       ["aryl/vinyl halide", "organostannane"], (["Ic1ccccc1", "C=C[Sn](CCCC)(CCCC)CCCC"], ["C=Cc1ccccc1"]),
       aliases=["stille", "organotin coupling"], byproducts=["Bu3SnX"]),
    # ------------------------------------------------------ organometallics
    _t("grignard_aldehyde", "Grignard addition to an aldehyde", "C–C bond formation",
       f"{_m(_ALD, 2)}=[O:3].[#6:1][Mg][Cl,Br,I]>>[#6:1][C:2]-[OH:3]",
       ["Et2O or THF"], ["0 °C → rt; then aq. NH4Cl"],
       "RMgX adds to an aldehyde; aqueous work-up gives the secondary alcohol (primary from formaldehyde).",
       ["aldehyde", "Grignard reagent"], (["O=Cc1ccccc1", "C[Mg]Br"], ["CC(O)c1ccccc1"]),
       aliases=["grignard aldehyde", "grignard reaction aldehyde", "rmgx aldehyde"]),
    _t("grignard_ketone", "Grignard addition to a ketone", "C–C bond formation",
       f"{_m(_KET, 2)}=[O:3].[#6:1][Mg][Cl,Br,I]>>[#6:1][C:2]-[OH:3]",
       ["Et2O or THF"], ["0 °C → rt; then aq. NH4Cl"],
       "RMgX adds to a ketone; work-up gives the tertiary alcohol.",
       ["ketone", "Grignard reagent"], (["CC(=O)c1ccccc1", "C[Mg]Br"], ["CC(C)(O)c1ccccc1"]),
       aliases=["grignard ketone", "grignard reaction ketone", "rmgx ketone"]),
    _t("grignard_ester", "Grignard addition to an ester (2 equiv)", "C–C bond formation",
       "[C;$(C(=O)([#6])O[#6]):2](=O)O[#6].[#6:1][Mg][Cl,Br,I].[#6:5][Mg][Cl,Br,I]>>[C:2]([#6:1])([#6:5])O",
       ["2 equiv RMgX, THF"], ["0 °C → rt; then aq. NH4Cl"],
       "Two equivalents of RMgX add to an ester via the ketone, giving a tertiary alcohol with two identical R groups.",
       ["ester", "Grignard reagent (2 equiv)"], (["CCOC(=O)c1ccccc1", "C[Mg]Br"], ["CC(C)(O)c1ccccc1"]),
       aliases=["grignard ester", "double grignard", "rmgx ester"], byproducts=["ROMgX"], stoichiometry=[1, 2]),
    # ---------------------------------------------------------- olefination
    _t("wittig_reaction", "Wittig reaction", "olefination",
       f"[C;$(C=P),$([C-][P+]):1]~[P;$(P(c)(c)c)].{_m(_KETALD, 2)}=[O]>>[C+0:1]=[C:2]",
       ["THF"], ["0 °C → rt"],
       "A phosphonium ylide (Ph3P=CR2, drawn as ylene or ylide) olefinates an aldehyde or ketone; E/Z not specified (non-stabilised ylides → Z, stabilised → E).",
       ["phosphonium ylide", "aldehyde/ketone"],
       (["C=P(c1ccccc1)(c1ccccc1)c1ccccc1", "O=Cc1ccccc1"], ["C=Cc1ccccc1"]),
       aliases=["wittig", "wittig olefination", "ylide olefination"], byproducts=["Ph3P=O"]),
    _t("horner_wadsworth_emmons", "Horner–Wadsworth–Emmons olefination", "olefination",
       "[#6;$(C=O),$(C#N):4][CH2,CH1:1]P(=O)(O[#6])O[#6].[CH:2](=O)[#6:5]>>[#6:4]/[C:1]=[C:2]/[#6:5]",
       ["NaH"], ["THF, 0 °C → rt"],
       "A stabilised phosphonate carbanion olefinates an aldehyde, (E)-selectively.",
       ["phosphonate (α-EWG)", "aldehyde"],
       (["CCOC(=O)CP(=O)(OCC)OCC", "O=Cc1ccccc1"], ["CCOC(=O)/C=C/c1ccccc1"]),
       aliases=["hwe", "horner wadsworth emmons", "horner emmons", "wadsworth emmons", "hwe olefination"],
       byproducts=["(RO)2PO2–"]),
    _t("cross_metathesis", "Olefin cross-metathesis", "olefination",
       "[CH2]=[C:1].[CH2]=[C:2]>>[C:1]=[C:2]",
       ["Grubbs II (5 mol%)"], ["DCM, 40 °C"],
       "Ru-catalysed cross-metathesis of two terminal alkenes releasing ethylene; E/Z not specified (usually E).",
       ["terminal alkene A", "terminal alkene B"],
       (["C=CCc1ccccc1", "C=CC(=O)OC"], ["COC(=O)C=CCc1ccccc1"]),
       aliases=["cross metathesis", "cm", "olefin metathesis", "grubbs"], byproducts=["C2H4"]),
    _t("ring_closing_metathesis", "Ring-closing metathesis", "olefination",
       "([CH2]=[C:1].[CH2]=[C:2])>>[C:1]=[C:2]",
       ["Grubbs II (5 mol%)"], ["DCM, 40 °C, –C2H4"],
       "Intramolecular metathesis of an α,ω-diene to a cycloalkene.",
       ["α,ω-diene"], (["C=CCC(CC=C)(C(=O)OCC)C(=O)OCC"], ["CCOC(=O)C1(C(=O)OCC)CC=CC1"]),
       aliases=["rcm", "ring closing metathesis", "grubbs rcm"], byproducts=["C2H4"]),
    # ------------------------------------------------- carbonyl C–C chemistry
    _t("aldol_addition", "Aldol addition", "C–C bond formation",
       f"{_m(_ENOL_KET, 1)}(=[O:2])[CH2,CH3:3].[CH:4](=[O:5])[#6:6]>>[C:1](=[O:2])[C:3][C:4](-[OH:5])[#6:6]",
       ["LDA"], ["THF, –78 °C; then RCHO"],
       "The enolate of a ketone adds to an aldehyde giving the β-hydroxy ketone (crossed/directed aldol).",
       ["enolisable ketone", "aldehyde"], (["CC(C)=O", "O=Cc1ccccc1"], ["CC(=O)CC(O)c1ccccc1"]),
       aliases=["aldol", "aldol reaction", "directed aldol", "crossed aldol"]),
    _t("aldol_condensation", "Aldol condensation (Claisen–Schmidt)", "C–C bond formation",
       f"{_m(_ENOL_KET, 1)}(=[O:2])[CH2,CH3:3].[CH:4](=O)[#6:6]>>[C:1](=[O:2])[C:3]=[C:4][#6:6]",
       ["NaOH"], ["EtOH/H2O, rt"],
       "Base-mediated aldol followed by dehydration to the α,β-unsaturated ketone (enone); E/Z not specified (usually E).",
       ["enolisable ketone", "aldehyde (non-enolisable)"],
       (["CC(C)=O", "O=Cc1ccccc1"], ["CC(=O)C=Cc1ccccc1"]),
       aliases=["claisen schmidt", "aldol condensation", "enone synthesis"], byproducts=["H2O"]),
    _t("claisen_condensation", "Claisen condensation", "C–C bond formation",
       "[#6:8][O:7][C:1](=[O:2])[CH2,CH3:3].[#6:6][C:4](=[O:5])O[#6]>>[#6:8][O:7][C:1](=[O:2])[C:3][C:4](=[O:5])[#6:6]",
       ["NaOEt"], ["EtOH, reflux"],
       "Ester enolate acylation by a second ester giving a β-keto ester.",
       ["ester (enolisable)", "ester (acceptor)"], (["CCOC(C)=O", "CCOC(C)=O"], ["CCOC(=O)CC(C)=O"]),
       aliases=["claisen", "claisen ester condensation", "beta keto ester synthesis"], byproducts=["ROH"]),
    _t("knoevenagel_condensation", "Knoevenagel condensation", "C–C bond formation",
       f"[CH:1](=O)[#6:2].{_m(_ACTIVE_CH2, 3)}>>[#6:2][CH:1]=[C:3]",
       ["piperidine (cat.)"], ["EtOH, reflux"],
       "Condensation of an aldehyde with an active-methylene compound (malonate, malononitrile, cyanoacetate ...) to the alkylidene.",
       ["aldehyde", "active methylene compound"], (["O=Cc1ccccc1", "N#CCC#N"], ["N#CC(C#N)=Cc1ccccc1"]),
       aliases=["knoevenagel", "doebner knoevenagel"], byproducts=["H2O"]),
    _t("michael_addition", "Michael addition", "C–C bond formation",
       f"{_m(_MICHAEL_DONOR, 1)}.{_MICHAEL_ACCEPTOR}>>[C:1][C:2]-[C:3][*:4]",
       ["NaOEt (cat.)"], ["EtOH, rt"],
       "Conjugate (1,4-) addition of a stabilised carbanion (1,3-dicarbonyl, malononitrile, nitroalkane) to an enone / acrylate / nitroalkene.",
       ["Michael donor", "Michael acceptor"],
       (["CCOC(=O)CC(=O)OCC", "C=CC(C)=O"], ["CCOC(=O)C(CCC(C)=O)C(=O)OCC"]),
       aliases=["michael", "conjugate addition", "1,4-addition", "michael reaction"]),
    _t("mannich_reaction", "Mannich reaction", "C–C bond formation",
       f"{_m(_ENOL_KET, 1)}(=[O:2])[CH2,CH3:3].{_m(_ALD, 4)}=O.{_m(_AMINE12, 5)}>>[C:1](=[O:2])[C:3][C:4][N:5]",
       ["HCl (cat.)"], ["EtOH, reflux"],
       "Three-component aminoalkylation of an enolisable ketone with an aldehyde (usually HCHO) and an amine to the β-amino ketone.",
       ["enolisable ketone", "aldehyde (HCHO)", "amine (1°/2°)"],
       (["CC(=O)c1ccccc1", "C=O", "CNC"], ["CN(C)CCC(=O)c1ccccc1"]),
       aliases=["mannich", "mannich base", "aminomethylation"], byproducts=["H2O"]),
    _t("robinson_annulation", "Robinson annulation", "C–C bond formation",
       f"{_m(_ENOL_KET, 1)}(=O)[CH2,CH1:3].[CH2:4]=[CH:5][C:6](=[O:7])[CH3:8]>>[C:1]1=[C:8][C:6](=[O:7])[C:5]-[C:4][C:3]1",
       ["KOH or pyrrolidine"], ["MeOH, reflux; –H2O"],
       "Michael addition of a ketone enolate to methyl vinyl ketone followed by intramolecular aldol condensation to the cyclohexenone.",
       ["enolisable ketone", "methyl vinyl ketone"], (["O=C1CCCCC1", "CC(=O)C=C"], ["O=C1CCC2CCCCC2=C1"]),
       aliases=["robinson", "annulation", "robinson ring annulation"], byproducts=["H2O"]),
    _t("decarboxylation", "Decarboxylation of a β-keto / malonic acid", "C–C bond formation",
       "[C;$(C=O):1][CX4:2]C(=O)[OH]>>[C:1][C:2]",
       [], ["Δ, –CO2"],
       "Thermal loss of CO2 from a β-keto acid or malonic acid.",
       ["β-keto acid / malonic acid"], (["CC(=O)CC(=O)O"], ["CC(C)=O"]),
       aliases=["krapcho", "thermal decarboxylation"], byproducts=["CO2"]),
    # --------------------------------------------------------- cycloadditions
    _t("diels_alder", "Diels–Alder cycloaddition", "cycloaddition",
       "[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1-[C:2]=[C:3]-[C:4][C:5]-[C:6]1",
       [], ["toluene, Δ"],
       "[4+2] cycloaddition of a 1,3-diene with a dienophile giving a cyclohexene; stereochemistry not specified.",
       ["diene", "dienophile"], (["C1=CCC=C1", "O=C1OC(=O)C=C1"], ["O=C1OC(=O)C2C3C=CC(C3)C12"]),
       aliases=["diels alder", "da reaction", "4+2 cycloaddition", "[4+2]"]),
    _t("cuaac_click", "CuAAC click reaction (azide–alkyne)", "cycloaddition",
       f"{_AZIDE}.[CH:5]#[C:6]>>[#6:4]-[N+0:1]1-[N+0:2]=[N+0:3]-[C:6]=[C:5]-1",
       ["CuSO4·5H2O, sodium ascorbate"], ["tBuOH/H2O, rt"],
       "Cu(I)-catalysed azide–alkyne cycloaddition giving the 1,4-disubstituted 1,2,3-triazole.",
       ["azide", "terminal alkyne"], (["[N-]=[N+]=NCc1ccccc1", "C#Cc1ccccc1"], ["c1ccc(Cn2cc(-c3ccccc3)nn2)cc1"]),
       aliases=["cuaac", "click", "click reaction", "huisgen", "azide alkyne cycloaddition", "triazole formation"]),
    _t("spaac_click", "Strain-promoted azide–alkyne cycloaddition (SPAAC)", "cycloaddition",
       f"{_AZIDE}.[C;R:5]#[C;R:6]>>[#6:4]-[N+0:1]1-[N+0:2]=[N+0:3]-[C:6]=[C:5]-1",
       [], ["rt, no catalyst"],
       "Copper-free click of an azide with a cyclooctyne; one of the two triazole regioisomers is drawn.",
       ["azide", "cyclooctyne"], (["[N-]=[N+]=NCc1ccccc1", "C1CCCC#CCC1"], ["c1ccc(Cn2nnc3c2CCCCCC3)cc1"]),
       aliases=["spaac", "copper free click", "strain promoted click", "cyclooctyne click"]),
    # ------------------------------------------------- nucleophilic substitution
    _t("sn2_amine", "SN2 alkylation of an amine", "substitution",
       f"{_m(_ALKYL_INV, 1)}{_LG}.{_m(_AMINE, 2)}>>[C@@:1][N:2]",
       ["K2CO3"], ["MeCN, 60 °C"],
       "N-alkylation of ammonia / a 1° or 2° amine with a primary or secondary alkyl halide or sulfonate, with inversion (over-alkylation not modelled).",
       ["alkyl halide", "amine"], (["BrCc1ccccc1", "CN"], ["CNCc1ccccc1"]),
       aliases=["n-alkylation", "amine alkylation", "sn2 with amine"], byproducts=["HX"]),
    _t("sn2_alkoxide", "SN2 with an alkoxide", "substitution",
       f"{_m(_ALKYL_INV, 1)}{_LG}.[O-;$(O[#6]):2]>>[C@@:1][O+0:2]",
       [], ["DMF, rt"],
       "O-alkylation of a preformed alkoxide / phenoxide (e.g. NaOEt) with an alkyl halide to the ether.",
       ["alkyl halide", "alkoxide"], (["CCBr", "[Na+].CC[O-]"], ["CCOCC"]),
       aliases=["alkoxide alkylation", "sn2 with alkoxide", "sn2 ether"], byproducts=["NaX"]),
    _t("sn2_azide", "SN2 with azide", "substitution",
       f"{_m(_ALKYL_INV, 1)}{_LG}.{_AZIDE_ION}>>[C@@:1][N+0:2]=[N+:3]=[N-:4]",
       [], ["DMF, 60 °C"],
       "Displacement of a halide / sulfonate by azide (NaN3) giving the alkyl azide, with inversion.",
       ["alkyl halide", "azide (NaN3)"], (["BrCc1ccccc1", "[Na+].[N-]=[N+]=[N-]"], ["[N-]=[N+]=NCc1ccccc1"]),
       aliases=["azide substitution", "sn2 with azide", "nan3", "azidation"], byproducts=["NaX"]),
    _t("sn2_cyanide", "SN2 with cyanide", "substitution",
       f"{_m(_ALKYL_INV, 1)}{_LG}.[C-:2]#[N:3]>>[C@@:1][C+0:2]#[N:3]",
       [], ["DMSO, 60 °C"],
       "Displacement of a halide / sulfonate by cyanide (NaCN/KCN) giving the nitrile (one-carbon homologation).",
       ["alkyl halide", "cyanide (NaCN)"], (["BrCc1ccccc1", "[Na+].[C-]#N"], ["N#CCc1ccccc1"]),
       aliases=["cyanation", "sn2 with cyanide", "nitrile from halide", "kolbe nitrile"], byproducts=["NaX"]),
    _t("sn2_thiol", "SN2 with a thiol / thiolate", "substitution",
       f"{_m(_ALKYL_INV, 1)}{_LG}.[S;$([SH][#6]),$([S-][#6]):2]>>[C@@:1][S+0:2]",
       ["K2CO3"], ["DMF, rt"],
       "S-alkylation of a thiol (or thiolate) with an alkyl halide / sulfonate giving the thioether.",
       ["alkyl halide", "thiol"], (["BrCc1ccccc1", "Sc1ccccc1"], ["c1ccc(CSc2ccccc2)cc1"]),
       aliases=["s-alkylation", "thioether synthesis", "sn2 with thiol", "thiol alkylation"], byproducts=["HX"]),
    _t("williamson_ether_synthesis", "Williamson ether synthesis", "substitution",
       f"{_m(_HYDROXYL, 1)}.{_m(_ALKYL_INV, 2)}{_LG}>>[O:1][C@@:2]",
       ["NaH"], ["DMF, 0 °C → rt"],
       "Deprotonated alcohol / phenol displaces a primary or secondary alkyl halide or sulfonate giving the ether.",
       ["alcohol/phenol", "alkyl halide"], (["Oc1ccccc1", "CI"], ["COc1ccccc1"]),
       aliases=["williamson", "o-alkylation", "ether synthesis", "etherification"], byproducts=["NaX"]),
    _t("mitsunobu_reaction", "Mitsunobu reaction", "substitution",
       f"[C@;X4;H1,H2:1][OH].{_m(_MITSUNOBU_NU, 3)}>>[C@@:1][*+0:3]",
       ["PPh3, DIAD"], ["THF, 0 °C → rt"],
       "Redox condensation of a 1°/2° alcohol with an acidic nucleophile (RCO2H, ArOH, imide, RSH, sulfonamide) with inversion at carbon.",
       ["alcohol (1°/2°)", "acidic nucleophile"],
       (["C[C@H](O)CC", "OC(=O)c1ccccc1"], ["CC[C@@H](C)OC(=O)c1ccccc1"]),
       aliases=["mitsunobu", "mitsunobu inversion", "diad pph3"], byproducts=["Ph3P=O", "DIAD-H2"]),
    _t("appel_reaction", "Appel reaction", "substitution",
       f"{_m(_ALKYL, 1)}[OH]>>[C:1]Br",
       ["CBr4, PPh3"], ["DCM, 0 °C → rt"],
       "Conversion of an alcohol into the alkyl bromide (CCl4 gives the chloride).",
       ["alcohol (1°/2°)"], (["OCc1ccccc1"], ["BrCc1ccccc1"]),
       aliases=["appel", "alcohol to bromide", "cbr4 pph3"], byproducts=["Ph3P=O", "CHBr3"]),
    _t("tosylation", "Tosylation of an alcohol", "substitution",
       "[CX4:1][OH:2]>>[C:1][O:2]S(=O)(=O)c1ccc(C)cc1",
       ["TsCl, pyridine"], ["DCM, 0 °C → rt"],
       "Alcohol → tosylate (activation of the hydroxyl as a leaving group).",
       ["alcohol"], (["CCO"], ["CCOS(=O)(=O)c1ccc(C)cc1"]),
       aliases=["tscl", "tosyl", "tosylate formation", "alcohol tosylation"], byproducts=["HCl"]),
    _t("mesylation", "Mesylation of an alcohol", "substitution",
       "[CX4:1][OH:2]>>[C:1][O:2]S(C)(=O)=O",
       ["MsCl, Et3N"], ["DCM, 0 °C"],
       "Alcohol → mesylate (activation of the hydroxyl as a leaving group).",
       ["alcohol"], (["CCO"], ["CCOS(C)(=O)=O"]),
       aliases=["mscl", "mesyl", "mesylate formation", "alcohol mesylation"], byproducts=["HCl"]),
    _t("e2_elimination", "E2 elimination", "elimination",
       f"[CX4:1]({_LG})[CX4;H1,H2,H3:2]>>[C:1]=[C:2]",
       ["KOtBu"], ["tBuOH, Δ"],
       "Base-induced β-elimination of HX from an alkyl halide / sulfonate; every regioisomeric alkene is generated (Zaitsev product usually dominates).",
       ["alkyl halide (with β-H)"], (["CC(C)(C)Br"], ["C=C(C)C"]),
       aliases=["e2", "elimination", "dehydrohalogenation", "beta elimination"], byproducts=["HX"]),
    # ------------------------------------------------- aromatic substitution
    _t("friedel_crafts_acylation", "Friedel–Crafts acylation", "aromatic substitution",
       "[cH:1].[C:2](=[O:3])[Cl,Br,$(OC(=O)[#6])]>>[c:1][C:2]=[O:3]",
       ["AlCl3"], ["DCM, 0 °C → rt"],
       "Lewis-acid catalysed acylation of an arene with an acyl chloride or anhydride; all positional isomers are generated.",
       ["arene", "acyl chloride/anhydride"], (["c1ccccc1", "CC(=O)Cl"], ["CC(=O)c1ccccc1"]),
       aliases=["friedel crafts acylation", "fc acylation", "aryl ketone synthesis"], byproducts=["HCl"]),
    _t("friedel_crafts_alkylation", "Friedel–Crafts alkylation", "aromatic substitution",
       "[cH:1].[CX4:2][Cl,Br,I]>>[c:1][C:2]",
       ["AlCl3"], ["DCM, rt"],
       "Lewis-acid catalysed alkylation of an arene with an alkyl halide (carbocation rearrangements not modelled); all positional isomers are generated.",
       ["arene", "alkyl halide"], (["c1ccccc1", "CC(C)Cl"], ["CC(C)c1ccccc1"]),
       aliases=["friedel crafts alkylation", "fc alkylation"], byproducts=["HX"]),
    _t("aromatic_nitration", "Electrophilic aromatic nitration", "aromatic substitution",
       "[cH:1]>>[c:1][N+](=O)[O-]",
       ["HNO3, H2SO4"], ["0 °C → rt"],
       "Nitration of an arene with mixed acid; all positional isomers are generated.",
       ["arene"], (["c1ccccc1"], ["O=[N+]([O-])c1ccccc1"]),
       aliases=["nitration", "eas nitration", "mixed acid nitration"], byproducts=["H2O"]),
    _t("aromatic_bromination", "Electrophilic aromatic bromination", "aromatic substitution",
       "[cH:1]>>[c:1]Br",
       ["Br2, FeBr3 (or NBS)"], ["DCM, rt"],
       "Bromination of an arene; all positional isomers are generated.",
       ["arene"], (["c1ccccc1"], ["Brc1ccccc1"]),
       aliases=["bromination", "eas bromination", "nbs bromination", "aryl bromide synthesis"], byproducts=["HBr"]),
    # -------------------------------------------------------------- reduction
    _t("nitro_reduction", "Nitro reduction to amine", "reduction",
       "[#6:1][N+:2](=O)[O-]>>[#6:1][NH2+0:2]",
       ["H2, Pd/C (or Fe/NH4Cl, SnCl2)"], ["MeOH, rt"],
       "Reduction of an aromatic or aliphatic nitro group to the primary amine.",
       ["nitro compound"], (["O=[N+]([O-])c1ccccc1"], ["Nc1ccccc1"]),
       aliases=["nitro to amine", "nitro group reduction", "bechamp", "sncl2 reduction"], byproducts=["H2O"]),
    _t("alkene_hydrogenation", "Hydrogenation of an alkene", "reduction",
       "[C:1]=[C:2]>>[C:1]-[C:2]",
       ["H2 (1 atm), Pd/C"], ["EtOAc, rt"],
       "Catalytic hydrogenation of a C=C double bond.",
       ["alkene"], (["C=Cc1ccccc1"], ["CCc1ccccc1"]),
       aliases=["hydrogenation", "alkene reduction", "h2 pd/c", "olefin hydrogenation"]),
    _t("alkyne_hydrogenation", "Hydrogenation of an alkyne to the alkane", "reduction",
       "[C:1]#[C:2]>>[C:1]-[C:2]",
       ["H2 (excess), Pd/C"], ["EtOAc, rt"],
       "Full catalytic hydrogenation of a C≡C triple bond.",
       ["alkyne"], (["c1ccc(C#Cc2ccccc2)cc1"], ["c1ccc(CCc2ccccc2)cc1"]),
       aliases=["alkyne reduction", "full hydrogenation", "alkyne to alkane"]),
    _t("lindlar_hydrogenation", "Lindlar hydrogenation to the (Z)-alkene", "reduction",
       "[#6:3][C:1]#[C:2][#6:4]>>[#6:3]/[C:1]=[C:2]\\[#6:4]",
       ["H2, Lindlar Pd (Pd/CaCO3, Pb)", "quinoline"], ["hexane, rt"],
       "Partial hydrogenation of an internal alkyne to the cis (Z) alkene.",
       ["internal alkyne"], (["CC#Cc1ccccc1"], ["C/C=C\\c1ccccc1"]),
       aliases=["lindlar", "alkyne to cis alkene", "alkyne to z alkene", "partial hydrogenation"]),
    _t("carbonyl_reduction", "Carbonyl reduction (NaBH4)", "reduction",
       f"{_m(_KETALD, 1)}=[O:2]>>[C:1]-[OH:2]",
       ["NaBH4"], ["MeOH, 0 °C"],
       "Hydride reduction of an aldehyde or ketone to the alcohol.",
       ["aldehyde/ketone"], (["CC(=O)c1ccccc1"], ["CC(O)c1ccccc1"]),
       aliases=["nabh4 reduction", "ketone reduction", "aldehyde reduction", "luche reduction", "borohydride reduction"]),
    _t("ester_reduction", "Ester reduction (LiAlH4)", "reduction",
       "[C;!$(C(=O)(O)[!#6]):1](=[O:4])[O:2][#6:3]>>([C:1]-[OH:4].[OH:2][#6:3])",
       ["LiAlH4"], ["THF, 0 °C → rt; then H2O"],
       "Reduction of an ester to the primary alcohol (plus the released alcohol).",
       ["ester"], (["COC(=O)c1ccccc1"], ["OCc1ccccc1", "CO"]),
       aliases=["lah ester reduction", "ester to alcohol", "lialh4 ester", "dibal"]),
    _t("amide_reduction", "Amide reduction (LiAlH4)", "reduction",
       "[C;$(C(=O)(N)[#6]):1](=O)[N:2]>>[C:1][N:2]",
       ["LiAlH4"], ["THF, reflux"],
       "Reduction of an amide to the amine (C=O → CH2).",
       ["amide"], (["CNC(=O)c1ccccc1"], ["CNCc1ccccc1"]),
       aliases=["lah amide reduction", "amide to amine", "lialh4 amide"]),
    _t("nitrile_reduction", "Nitrile reduction to the amine", "reduction",
       "[#6:2][C:1]#[N:3]>>[#6:2][CH2:1]-[NH2:3]",
       ["LiAlH4 (or H2, Raney Ni)"], ["THF, 0 °C → reflux"],
       "Reduction of a nitrile to the primary amine.",
       ["nitrile"], (["N#Cc1ccccc1"], ["NCc1ccccc1"]),
       aliases=["nitrile to amine", "lah nitrile", "raney nickel nitrile"]),
    _t("staudinger_reduction", "Staudinger reduction", "reduction",
       f"{_AZIDE}>>[#6:4][NH2+0:1]",
       ["PPh3; then H2O"], ["THF, rt"],
       "Azide → primary amine via the iminophosphorane.",
       ["azide"], (["[N-]=[N+]=NCc1ccccc1"], ["NCc1ccccc1"]),
       aliases=["staudinger", "azide reduction", "azide to amine"], byproducts=["N2", "Ph3P=O"]),
    _t("wolff_kishner_reduction", "Wolff–Kishner reduction", "reduction",
       f"{_m(_KETALD, 1)}=O>>[C:1]",
       ["N2H4, KOH"], ["ethylene glycol, 200 °C"],
       "Deoxygenation of an aldehyde or ketone to the methylene via the hydrazone (basic conditions).",
       ["aldehyde/ketone"], (["CC(=O)c1ccccc1"], ["CCc1ccccc1"]),
       aliases=["wolff kishner", "huang minlon", "carbonyl deoxygenation"], byproducts=["N2", "H2O"]),
    _t("clemmensen_reduction", "Clemmensen reduction", "reduction",
       f"{_m(_KETALD, 1)}=O>>[C:1]",
       ["Zn(Hg), HCl (conc.)"], ["reflux"],
       "Deoxygenation of an aldehyde or ketone to the methylene (acidic conditions).",
       ["aldehyde/ketone"], (["CC(=O)c1ccccc1"], ["CCc1ccccc1"]),
       aliases=["clemmensen", "zinc amalgam reduction"], byproducts=["H2O"]),
    _t("birch_reduction", "Birch reduction", "reduction",
       "[cH:1]1[cH:2][c:3][cH:4][cH:5][c:6]1>>[CH2:1]1[CH:2]=[C:3][CH2:4][CH:5]=[C:6]1",
       ["Na, NH3 (l)", "EtOH"], ["–78 °C"],
       "Dissolving-metal reduction of a benzene ring to the 1,4-cyclohexadiene; substituents are placed on the sp2 carbons (correct for donors such as OMe; acceptors would end up on sp3 carbons).",
       ["arene"], (["COc1ccccc1"], ["COC1=CCC=CC1"]),
       aliases=["birch", "dissolving metal reduction", "sodium ammonia reduction"]),
    # -------------------------------------------------------------- oxidation
    _t("swern_oxidation", "Swern oxidation", "oxidation",
       "[CX4;H1,H2;$(C[#6]):1][OH:2]>>[C:1]=[O:2]",
       ["(COCl)2, DMSO; Et3N"], ["DCM, –78 °C → rt"],
       "Oxidation of a primary alcohol to the aldehyde (or a secondary alcohol to the ketone) under activated-DMSO conditions.",
       ["alcohol (1°/2°)"], (["OCc1ccccc1"], ["O=Cc1ccccc1"]),
       aliases=["swern", "alcohol to aldehyde", "dmso oxidation"], byproducts=["Me2S", "CO", "CO2"]),
    _t("dess_martin_oxidation", "Dess–Martin oxidation", "oxidation",
       "[CX4;H1,H2;$(C[#6]):1][OH:2]>>[C:1]=[O:2]",
       ["Dess–Martin periodinane"], ["DCM, rt"],
       "Hypervalent-iodine oxidation of a primary alcohol to the aldehyde or a secondary alcohol to the ketone.",
       ["alcohol (1°/2°)"], (["CC(O)c1ccccc1"], ["CC(=O)c1ccccc1"]),
       aliases=["dess martin", "dmp", "dmp oxidation", "periodinane oxidation"]),
    _t("secondary_alcohol_oxidation", "Secondary alcohol → ketone (PCC)", "oxidation",
       "[CX4;H1;$(C([#6])[#6]):1][OH:2]>>[C:1]=[O:2]",
       ["PCC"], ["DCM, rt"],
       "Oxidation of a secondary alcohol to the ketone (PCC, TPAP/NMO, Jones ...).",
       ["secondary alcohol"], (["OC1CCCCC1"], ["O=C1CCCCC1"]),
       aliases=["pcc oxidation", "alcohol to ketone", "ketone from alcohol", "tpap"]),
    _t("jones_oxidation", "Jones oxidation (primary alcohol → acid)", "oxidation",
       "[CH2;$(C[#6]):1][OH:2]>>[C:1](=O)[OH:2]",
       ["CrO3, H2SO4 (Jones)"], ["acetone, 0 °C"],
       "Oxidation of a primary alcohol all the way to the carboxylic acid.",
       ["primary alcohol"], (["OCc1ccccc1"], ["O=C(O)c1ccccc1"]),
       aliases=["jones", "alcohol to acid", "chromic acid oxidation"]),
    _t("pinnick_oxidation", "Pinnick oxidation (aldehyde → acid)", "oxidation",
       "[CH:1](=[O:2])[#6:3]>>[#6:3][C:1](=[O:2])O",
       ["NaClO2, NaH2PO4", "2-methyl-2-butene"], ["tBuOH/H2O, rt"],
       "Chlorite oxidation of an aldehyde to the carboxylic acid.",
       ["aldehyde"], (["O=Cc1ccccc1"], ["O=C(O)c1ccccc1"]),
       aliases=["pinnick", "lindgren oxidation", "aldehyde to acid", "chlorite oxidation"]),
    _t("baeyer_villiger_oxidation", "Baeyer–Villiger oxidation", "oxidation",
       "[#6:3][C:1](=[O:2])[#6;!$([CH3]):4]>>[#6:3][C:1](=[O:2])O[#6:4]",
       ["mCPBA"], ["DCM, rt"],
       "Peracid oxidation of a ketone to the ester / lactone; oxygen is inserted next to the non-methyl (more substituted) group.",
       ["ketone"], (["O=C1CCCCC1"], ["O=C1CCCCCO1"]),
       aliases=["baeyer villiger", "bv oxidation", "lactone from ketone", "ketone to ester"]),
    _t("ozonolysis_reductive", "Ozonolysis (reductive work-up)", "oxidation",
       "[C:1]=[C:2]>>([C:1]=O.[C:2]=O)",
       ["1. O3, DCM, –78 °C", "2. Me2S"], [],
       "Oxidative cleavage of an alkene to two carbonyl compounds (aldehydes/ketones).",
       ["alkene"], (["CC1=CCCCC1"], ["CC(=O)CCCCC=O"]),
       aliases=["ozonolysis", "alkene cleavage", "o3 me2s"], byproducts=["Me2SO"]),
    _t("epoxidation", "Epoxidation (mCPBA)", "oxidation",
       "[C:1]=[C:2]>>[C:1]1-[C:2]O1",
       ["mCPBA"], ["DCM, 0 °C → rt"],
       "Peracid epoxidation of an alkene (stereospecific, cis-alkene → cis-epoxide; not drawn).",
       ["alkene"], (["C1CCC=CC1"], ["C1CCC2OC2C1"]),
       aliases=["mcpba epoxidation", "prilezhaev", "alkene epoxidation"], byproducts=["mCBA"]),
    _t("sharpless_dihydroxylation", "Sharpless / Upjohn dihydroxylation", "oxidation",
       "[C:1]=[C:2]>>O[C:1]-[C:2]O",
       ["AD-mix-β (or OsO4, NMO)"], ["tBuOH/H2O, 0 °C"],
       "syn-Dihydroxylation of an alkene to the 1,2-diol; drawn without stereochemistry.",
       ["alkene"], (["C=Cc1ccccc1"], ["OCC(O)c1ccccc1"]),
       aliases=["sharpless", "dihydroxylation", "upjohn", "oso4", "ad mix", "asymmetric dihydroxylation",
                "syn dihydroxylation"]),
    _t("alkene_dibromination", "Alkene dibromination", "addition",
       "[C:1]=[C:2]>>Br[C:1]-[C:2]Br",
       ["Br2"], ["DCM, 0 °C"],
       "Anti addition of Br2 across an alkene giving the vicinal dibromide (stereochemistry not drawn).",
       ["alkene"], (["C1CCC=CC1"], ["BrC1CCCCC1Br"]),
       aliases=["bromine addition", "dibromination", "vicinal dibromide"]),
    _t("hydroboration_oxidation", "Hydroboration–oxidation", "addition",
       f"{_m(_ANTI_MARK, 1)}=[C:2]>>[OH][C:1]-[C:2]",
       ["1. BH3·THF", "2. H2O2, NaOH"], [],
       "Anti-Markovnikov hydration of an alkene: OH on the less substituted carbon.",
       ["alkene"], (["C=Cc1ccccc1"], ["OCCc1ccccc1"]),
       aliases=["hydroboration", "anti markovnikov hydration", "bh3 h2o2"]),
    _t("markovnikov_hydration", "Markovnikov hydration (oxymercuration)", "addition",
       f"{_m(_ANTI_MARK, 1)}=[C:2]>>[C:1]-[C:2][OH]",
       ["1. Hg(OAc)2, H2O", "2. NaBH4"], [],
       "Markovnikov hydration of an alkene: OH on the more substituted carbon (oxymercuration–demercuration or H3O+).",
       ["alkene"], (["C=Cc1ccccc1"], ["CC(O)c1ccccc1"]),
       aliases=["oxymercuration", "markovnikov", "alkene hydration", "acid catalysed hydration"]),
    _t("epoxide_opening_amine", "Epoxide opening by an amine", "addition",
       f"{_m(_EPOX_C1, 1)}1[O:3][C:2]1.{_m(_AMINE12, 4)}>>[N:4][C:1][C:2][OH:3]",
       [], ["EtOH, 60 °C"],
       "Nucleophilic ring opening of an epoxide by a 1°/2° amine at the less hindered carbon giving the β-amino alcohol.",
       ["epoxide", "amine (1°/2°)"], (["C1OC1c1ccccc1", "CN"], ["CNCC(O)c1ccccc1"]),
       aliases=["epoxide aminolysis", "amino alcohol from epoxide", "epoxide ring opening"]),
    # --------------------------------------------------- C–N bond formation
    _t("buchwald_hartwig_amination", "Buchwald–Hartwig amination", "C–N coupling",
       f"[c:1]{_ARLG}.{_m(_AMINE12, 2)}>>[c:1][N:2]",
       ["Pd2(dba)3, XPhos", "NaOtBu"], ["toluene, 100 °C"],
       "Pd-catalysed C–N coupling of an aryl halide / triflate with a 1° or 2° amine.",
       ["aryl halide", "amine (1°/2°)"], (["Brc1ccccc1", "C1COCCN1"], ["c1ccc(N2CCOCC2)cc1"]),
       aliases=["buchwald", "buchwald hartwig", "c-n coupling", "pd amination", "aryl amination"],
       byproducts=["HX"]),
    _t("reductive_amination", "Reductive amination", "C–N coupling",
       f"{_m(_KETALD, 1)}=O.{_m(_AMINE12, 2)}>>[C:1][N:2]",
       ["NaBH(OAc)3"], ["DCE, rt"],
       "Condensation of an aldehyde / ketone with a 1° or 2° amine and in situ reduction of the imine / iminium.",
       ["aldehyde/ketone", "amine (1°/2°)"], (["O=Cc1ccccc1", "CN"], ["CNCc1ccccc1"]),
       aliases=["reductive amination", "stab", "nabh(oac)3", "borch reduction", "nabh3cn amination"],
       byproducts=["H2O"]),
    _t("imine_formation", "Imine (Schiff base) formation", "condensation",
       f"{_m(_KETALD, 1)}=O.{_m(_PRIM_AMINE, 2)}>>[C:1]=[N:2]",
       ["MgSO4 or 4 Å MS"], ["DCM, rt"],
       "Condensation of an aldehyde / ketone with a primary amine to the imine.",
       ["aldehyde/ketone", "primary amine"], (["O=Cc1ccccc1", "Nc1ccccc1"], ["C(=Nc1ccccc1)c1ccccc1"]),
       aliases=["imine", "schiff base", "imine condensation"], byproducts=["H2O"]),
    _t("oxime_formation", "Oxime formation", "condensation",
       f"{_m(_KETALD, 1)}=O.[NH2:2][OH:3]>>[C:1]=[N:2][OH:3]",
       ["NH2OH·HCl, NaOAc"], ["EtOH, rt"],
       "Condensation of an aldehyde / ketone with hydroxylamine to the oxime.",
       ["aldehyde/ketone", "hydroxylamine"], (["CC(C)=O", "NO"], ["CC(C)=NO"]),
       aliases=["oxime", "hydroxylamine condensation"], byproducts=["H2O"]),
    # ------------------------------------------------------------ protection
    _t("boc_protection", "Boc protection of an amine", "protection",
       f"{_m(_AMINE12, 1)}>>[N:1]C(=O)OC(C)(C)C",
       ["Boc2O", "Et3N"], ["DCM, rt"],
       "N-tert-butoxycarbonylation of a 1° or 2° amine.",
       ["amine (1°/2°)"], (["Nc1ccccc1"], ["CC(C)(C)OC(=O)Nc1ccccc1"]),
       aliases=["boc", "boc2o", "boc amine", "n-boc protection"], byproducts=["tBuOH", "CO2"]),
    _t("boc_deprotection", "Boc deprotection", "deprotection",
       "[N:1]C(=O)OC(C)(C)C>>[N:1]",
       ["TFA (or 4 M HCl/dioxane)"], ["DCM, rt"],
       "Acidolytic removal of a Boc group.",
       ["Boc-amine"], (["CC(C)(C)OC(=O)Nc1ccccc1"], ["Nc1ccccc1"]),
       aliases=["boc removal", "tfa deprotection", "boc cleavage", "n-boc deprotection"],
       byproducts=["isobutene", "CO2"]),
    _t("fmoc_protection", "Fmoc protection of an amine", "protection",
       f"{_m(_AMINE12, 1)}>>[N:1]C(=O)OCC1c2ccccc2-c2ccccc12",
       ["Fmoc-OSu (or Fmoc-Cl)", "NaHCO3"], ["dioxane/H2O, rt"],
       "N-(9-fluorenylmethoxycarbonyl)ation of a 1° or 2° amine.",
       ["amine (1°/2°)"], (["NCC(=O)O"], ["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc12"]),
       aliases=["fmoc", "fmoc amine", "n-fmoc protection", "fmoc osu"]),
    _t("fmoc_deprotection", "Fmoc deprotection", "deprotection",
       "[N:1]C(=O)OCC1c2ccccc2-c2ccccc12>>[N:1]",
       ["20% piperidine"], ["DMF, rt"],
       "Base-mediated (β-elimination) removal of an Fmoc group.",
       ["Fmoc-amine"], (["O=C(O)CNC(=O)OCC1c2ccccc2-c2ccccc12"], ["NCC(=O)O"]),
       aliases=["fmoc removal", "piperidine deprotection", "fmoc cleavage"], byproducts=["dibenzofulvene", "CO2"]),
    _t("tbs_protection", "TBS protection of an alcohol", "protection",
       f"{_m(_HYDROXYL, 1)}>>[O:1][Si](C)(C)C(C)(C)C",
       ["TBSCl, imidazole"], ["DMF, rt"],
       "Silylation of an alcohol / phenol as the tert-butyldimethylsilyl ether.",
       ["alcohol/phenol"], (["OCc1ccccc1"], ["C[Si](C)(OCc1ccccc1)C(C)(C)C"]),
       aliases=["tbs", "tbdms", "silyl protection", "tbscl", "silylation"], byproducts=["HCl"]),
    _t("tbs_deprotection", "Silyl ether deprotection (TBAF)", "deprotection",
       "[O:1][Si]([#6])([#6])[#6]>>[O:1]",
       ["TBAF"], ["THF, rt"],
       "Fluoride-mediated cleavage of a silyl ether (TBS, TES, TIPS, TMS) back to the alcohol.",
       ["silyl ether"], (["C[Si](C)(OCc1ccccc1)C(C)(C)C"], ["OCc1ccccc1"]),
       aliases=["tbaf", "tbs removal", "desilylation", "silyl deprotection", "tbs cleavage"],
       byproducts=["R3SiF"]),
    _t("benzyl_protection", "Benzyl protection of an alcohol", "protection",
       f"{_m(_HYDROXYL, 1)}>>[O:1]Cc1ccccc1",
       ["BnBr, NaH"], ["DMF, 0 °C → rt"],
       "O-benzylation of an alcohol / phenol.",
       ["alcohol/phenol"], (["OCCO"], ["OCCOCc1ccccc1"]),
       aliases=["benzylation", "bn protection", "benzyl ether formation", "bnbr"], byproducts=["NaBr"]),
    _t("benzyl_deprotection", "Benzyl deprotection (hydrogenolysis)", "deprotection",
       "[OX2;H0:1]Cc1ccccc1>>[O:1]",
       ["H2, Pd/C"], ["MeOH, rt"],
       "Hydrogenolysis of a benzyl ether / ester to the free alcohol / acid.",
       ["benzyl ether/ester"], (["COc1ccc(OCc2ccccc2)cc1"], ["COc1ccc(O)cc1"]),
       aliases=["hydrogenolysis", "bn removal", "debenzylation", "benzyl ether cleavage"], byproducts=["toluene"]),
    _t("acetylation_alcohol", "O-Acetylation (Ac2O)", "protection",
       f"{_m(_HYDROXYL, 1)}>>[O:1]C(C)=O",
       ["Ac2O, pyridine", "DMAP (cat.)"], ["DCM, rt"],
       "Acetylation of an alcohol / phenol with acetic anhydride.",
       ["alcohol/phenol"], (["OCc1ccccc1"], ["CC(=O)OCc1ccccc1"]),
       aliases=["acetylation", "ac2o", "acetate protection", "o-acetylation"], byproducts=["AcOH"]),
    _t("acetal_formation", "Cyclic acetal (dioxolane) formation", "protection",
       f"{_m(_KETALD, 1)}=O.[OH:2][C:3][C:4][OH:5]>>[C:1]1[O:2][C:3][C:4][O:5]1",
       ["p-TsOH (cat.)"], ["toluene, Dean–Stark"],
       "Protection of an aldehyde / ketone as the 1,3-dioxolane with a 1,2-diol.",
       ["aldehyde/ketone", "1,2-diol"], (["O=C1CCCCC1", "OCCO"], ["C1CCC2(CC1)OCCO2"]),
       aliases=["acetal", "ketal", "dioxolane", "ketal formation", "acetalisation", "ethylene glycol protection"],
       byproducts=["H2O"]),
    _t("acetal_formation_dialkyl", "Dialkyl acetal formation (2 ROH)", "protection",
       f"{_m(_KETALD, 1)}=O.[OH:2][CX4:3].[OH:4][CX4:5]>>[C:1]([O:2][C:3])[O:4][C:5]",
       ["HC(OMe)3, p-TsOH (cat.)"], ["MeOH, rt"],
       "Acetalisation of an aldehyde / ketone with two equivalents of a simple alcohol.",
       ["aldehyde/ketone", "alcohol (2 equiv)"], (["O=Cc1ccccc1", "CO"], ["COC(OC)c1ccccc1"]),
       aliases=["dimethyl acetal", "diethyl acetal", "dialkyl acetal", "acetal from alcohol"],
       byproducts=["H2O"], stoichiometry=[1, 2]),
    _t("acetal_hydrolysis", "Acetal hydrolysis", "deprotection",
       "[CX4;!$(C(O)(O)O):1]([OX2;$(O[#6])])[OX2;$(O[#6])]>>[C:1]=O",
       ["aq. HCl (or PPTS)"], ["acetone/H2O, rt"],
       "Acid-catalysed cleavage of an acetal / ketal back to the carbonyl compound (the alcohol / diol is not drawn).",
       ["acetal/ketal"], (["COC(OC)c1ccccc1"], ["O=Cc1ccccc1"]),
       aliases=["acetal deprotection", "ketal hydrolysis", "ketal deprotection", "acetal cleavage"],
       byproducts=["ROH / diol"]),
    # ------------------------------------------------------------ hydrolysis
    _t("ester_hydrolysis", "Ester hydrolysis (saponification)", "hydrolysis",
       "[C;!$(C(=O)(O)[!#6]):1](=[O:2])[O:3][#6:4]>>([C:1](=[O:2])O.[OH:3][#6:4])",
       ["LiOH (or NaOH)"], ["THF/H2O, rt; then HCl"],
       "Base-mediated ester cleavage to the carboxylic acid and the alcohol (lactones give the hydroxy acid).",
       ["ester"], (["CCOC(C)=O"], ["CC(=O)O", "CCO"]),
       aliases=["saponification", "ester cleavage", "ester saponification", "lioh hydrolysis"]),
    _t("amide_hydrolysis", "Amide hydrolysis", "hydrolysis",
       "[C;$(C(=O)(N)[#6]):1](=[O:2])[N:3]>>([C:1](=[O:2])O.[N:3])",
       ["6 M HCl (or NaOH)"], ["reflux"],
       "Hydrolysis of an amide to the carboxylic acid and the amine (lactams give the amino acid).",
       ["amide"], (["CC(=O)Nc1ccccc1"], ["CC(=O)O", "Nc1ccccc1"]),
       aliases=["amide cleavage", "amide to acid"]),
    _t("nitrile_hydrolysis", "Nitrile hydrolysis", "hydrolysis",
       "[#6:2][C:1]#[N]>>[#6:2][C:1](=O)[OH]",
       ["H2SO4 (aq.)"], ["reflux"],
       "Hydrolysis of a nitrile to the carboxylic acid (via the primary amide).",
       ["nitrile"], (["N#Cc1ccccc1"], ["O=C(O)c1ccccc1"]),
       aliases=["nitrile to acid", "cyano hydrolysis"], byproducts=["NH3"]),
    # ---------------------------------------------------------- rearrangement
    _t("curtius_rearrangement", "Curtius rearrangement (acyl azide → isocyanate)", "rearrangement",
       "[#6:4][C:1](=[O:2])[N;$(N=[N+]=[N-]),$([N-][N+]#N):3]~[N+]~[NX1]>>[#6:4][N+0:3]=[C:1]=[O:2]",
       [], ["toluene, 80 °C, –N2"],
       "Thermal rearrangement of an acyl azide to the isocyanate with loss of N2.",
       ["acyl azide"], (["[N-]=[N+]=NC(=O)c1ccccc1"], ["O=C=Nc1ccccc1"]),
       aliases=["curtius", "acyl azide rearrangement"], byproducts=["N2"]),
    _t("isocyanate_hydrolysis", "Isocyanate hydrolysis to the amine", "hydrolysis",
       "[#6:4][N:3]=C=O>>[#6:4][NH2:3]",
       ["H2O"], ["Δ, –CO2"],
       "Hydrolysis of an isocyanate via the carbamic acid to the primary amine (second step of Curtius / Hofmann degradations).",
       ["isocyanate"], (["O=C=Nc1ccccc1"], ["Nc1ccccc1"]),
       aliases=["isocyanate to amine", "carbamic acid decarboxylation"], byproducts=["CO2"]),
    _t("curtius_to_amine", "Curtius degradation (acid → amine, one pot)", "rearrangement",
       "[#6:4]C(=O)[OH]>>[#6:4][NH2]",
       ["DPPA, Et3N", "then H2O, Δ"], ["toluene, 80 °C"],
       "Carboxylic acid → acyl azide → isocyanate → amine with loss of one carbon (DPPA one-pot variant).",
       ["carboxylic acid"], (["OC(=O)c1ccccc1"], ["Nc1ccccc1"]),
       aliases=["curtius degradation", "dppa amination", "acid to amine"], byproducts=["N2", "CO2"]),
    _t("hofmann_rearrangement", "Hofmann rearrangement", "rearrangement",
       "[#6:4]C(=O)[NH2:3]>>[#6:4][NH2:3]",
       ["Br2, NaOH"], ["H2O, 0 °C → 70 °C"],
       "Primary amide → amine with loss of the carbonyl carbon via the isocyanate.",
       ["primary amide"], (["NC(=O)c1ccccc1"], ["Nc1ccccc1"]),
       aliases=["hofmann", "hofmann degradation", "amide to amine (hofmann)"], byproducts=["CO2", "NaBr"]),
    # --------------------------------------------------- heterocycle synthesis
    _t("fischer_indole_synthesis", "Fischer indole synthesis", "heterocycle synthesis",
       "[NH2][NH:1][c:2]1[cH:3][c:4][c:5][c:6][c:7]1.[CX3;$(C=O);!$(C(=O)[!#6]):8](=O)[CH2,CH3:9]"
       ">>[NH:1]1[C:8]=[C:9][c:3]2[c:4][c:5][c:6][c:7][c:2]21",
       ["H2SO4 or ZnCl2"], ["AcOH, reflux"],
       "Aryl hydrazine + enolisable ketone / aldehyde → indole ([3,3]-sigmatropic rearrangement of the ene-hydrazine).",
       ["aryl hydrazine", "enolisable ketone/aldehyde"], (["NNc1ccccc1", "O=C1CCCCC1"], ["c1ccc2c(c1)[nH]c1c2CCCC1"]),
       aliases=["fischer indole", "indole synthesis", "fischer indolisation"], byproducts=["NH3", "H2O"]),
    _t("paal_knorr_pyrrole", "Paal–Knorr pyrrole synthesis", "heterocycle synthesis",
       f"[#6:7][C:1](=O)[CH2:2][CH2:3][C:4](=O)[#6:8].{_m(_PRIM_AMINE, 5)}>>[#6:7][C:1]1=[C:2][C:3]=[C:4]([#6:8])[N:5]1",
       ["AcOH (cat.)"], ["toluene, reflux"],
       "1,4-Diketone + primary amine → N-substituted pyrrole.",
       ["1,4-diketone", "primary amine"], (["CC(=O)CCC(C)=O", "NCc1ccccc1"], ["Cc1ccc(C)n1Cc1ccccc1"]),
       aliases=["paal knorr", "pyrrole synthesis"], byproducts=["2 H2O"]),
    _t("paal_knorr_furan", "Paal–Knorr furan synthesis", "heterocycle synthesis",
       "[#6:7][C:1](=[O:9])[CH2:2][CH2:3][C:4](=O)[#6:8]>>[#6:7][C:1]1=[C:2][C:3]=[C:4]([#6:8])[O:9]-1",
       ["p-TsOH or P2O5"], ["toluene, reflux"],
       "Acid-catalysed cyclodehydration of a 1,4-diketone to the furan.",
       ["1,4-diketone"], (["CC(=O)CCC(C)=O"], ["Cc1ccc(C)o1"]),
       aliases=["furan synthesis", "paal knorr furan"], byproducts=["H2O"]),
    # ---------------------------------------------------------- multicomponent
    _t("biginelli_reaction", "Biginelli reaction", "multicomponent",
       "[CH:1](=O)[#6:2].[CH3:4][C:3](=O)[CH2:5][C:6](=[O:12])[O:7][#6:13].[NH2:8][C:9](=[O,S:11])[NH2:10]"
       ">>[C:1]1([#6:2])[NH:8][C:9](=[O,S:11])[NH:10][C:3]([CH3:4])=[C:5]1[C:6](=[O:12])[O:7][#6:13]",
       ["HCl (cat.)"], ["EtOH, reflux"],
       "Aldehyde + β-keto ester + (thio)urea → 3,4-dihydropyrimidin-2(1H)-one.",
       ["aldehyde", "β-keto ester", "urea/thiourea"],
       (["O=Cc1ccccc1", "CCOC(=O)CC(C)=O", "NC(N)=O"], ["CCOC(=O)C1=C(C)NC(=O)NC1c1ccccc1"]),
       aliases=["biginelli", "dihydropyrimidinone synthesis", "dhpm"], byproducts=["2 H2O"]),
    _t("ugi_reaction", "Ugi four-component reaction", "multicomponent",
       f"{_m(_ALD, 1)}=O.{_m(_PRIM_AMINE, 3)}.{_m(_ACID, 5)}(=[O:8])[OH].{_ISOCYANIDE}"
       ">>[C:1]([N:3][C:5]=[O:8])[C+0:7](=O)[NH+0:6][#6:9]",
       [], ["MeOH, rt"],
       "Aldehyde + primary amine + carboxylic acid + isocyanide → α-acylamino amide (bis-amide).",
       ["aldehyde", "primary amine", "carboxylic acid", "isocyanide"],
       (["O=Cc1ccccc1", "CN", "CC(=O)O", "CC(C)(C)[N+]#[C-]"], ["CC(=O)N(C)C(C(=O)NC(C)(C)C)c1ccccc1"]),
       aliases=["ugi", "ugi 4cr", "u-4cr", "four component"], byproducts=["H2O"]),
    _t("passerini_reaction", "Passerini three-component reaction", "multicomponent",
       f"{_m(_ALD, 1)}=[O:10].{_m(_ACID, 5)}(=[O:8])[OH].{_ISOCYANIDE}"
       ">>[C:1]([O:10][C:5]=[O:8])[C+0:7](=O)[NH+0:6][#6:9]",
       [], ["DCM, rt"],
       "Aldehyde + carboxylic acid + isocyanide → α-acyloxy amide.",
       ["aldehyde", "carboxylic acid", "isocyanide"],
       (["O=Cc1ccccc1", "CC(=O)O", "CC(C)(C)[N+]#[C-]"], ["CC(=O)OC(C(=O)NC(C)(C)C)c1ccccc1"]),
       aliases=["passerini", "passerini 3cr", "p-3cr", "three component"]),
]

REACTIONS: Dict[str, ReactionTemplate] = {}
for _tmpl in _TEMPLATES:
    if _tmpl.key in REACTIONS:  # pragma: no cover - guarded by tests
        raise RuntimeError(f"duplicate reaction key {_tmpl.key!r}")
    REACTIONS[_tmpl.key] = _tmpl
del _tmpl


def classes() -> List[str]:
    """Sorted list of the reaction classes used in the registry."""
    return sorted({t.cls for t in REACTIONS.values()})


# --------------------------------------------------------------------------- #
# Lookup
# --------------------------------------------------------------------------- #
_DASHES = str.maketrans({"‐": "-", "‑": "-", "–": "-", "—": "-", "−": "-", "’": "'", "`": "'"})


def _norm(text: str) -> str:
    """Lower-case, ASCII dashes, strip everything but letters and digits."""
    s = str(text).translate(_DASHES).lower()
    return re.sub(r"[^a-z0-9]+", "", s)


def _lookup_index() -> List[Tuple[str, str]]:
    """``(normalised label, key)`` for every key, display name and alias."""
    out: List[Tuple[str, str]] = []
    for key, t in REACTIONS.items():
        labels = [key, t.name, *t.aliases]
        for label in labels:
            n = _norm(label)
            if n:
                out.append((n, key))
    return out


_INDEX: List[Tuple[str, str]] = _lookup_index()
_EXACT: Dict[str, str] = {}
for _n, _k in _INDEX:
    _EXACT.setdefault(_n, _k)
del _n, _k


def get_reaction(name_or_alias: Union[str, ReactionTemplate]) -> ReactionTemplate:
    """Find a template by key, display name or alias.

    Matching ignores case, whitespace, punctuation and dash variants
    (``"Suzuki-Miyaura"``, ``"suzuki miyaura coupling"``, ``"SUZUKI_COUPLING"``
    all work).  After exact matching, a unique substring match (the query
    contained in a key, name or alias) and then a close fuzzy match are
    tried; failures raise :class:`ReactionError` listing the nearest
    candidates.
    """
    if isinstance(name_or_alias, ReactionTemplate):
        return name_or_alias
    if not isinstance(name_or_alias, str) or not name_or_alias.strip():
        raise ReactionError("Give a reaction name, key or alias (see list_reactions()).")
    q = _norm(name_or_alias)
    if not q:
        raise ReactionError(f"Reaction name {name_or_alias!r} contains no letters or digits.")
    key = _EXACT.get(q)
    if key is not None:
        return REACTIONS[key]
    # drop trivial suffixes ("reaction", "coupling", "synthesis") and retry exact
    for suffix in ("reaction", "coupling", "synthesis", "oxidation", "reduction", "rearrangement"):
        if q.endswith(suffix) and len(q) > len(suffix):
            key = _EXACT.get(q[: -len(suffix)])
            if key is not None:
                return REACTIONS[key]
    # unique substring match (query contained in a key / name / alias).  The
    # reverse direction (label contained in the query) is deliberately not
    # used: short generic aliases such as "bromination" or "acetylation" would
    # otherwise hijack unrelated queries like "alkene bromination" and silently
    # select the wrong template.
    hits = {k for n, k in _INDEX if len(q) >= 4 and q in n}
    if len(hits) == 1:
        return REACTIONS[hits.pop()]
    # fuzzy
    close = get_close_matches(q, [n for n, _ in _INDEX], n=5, cutoff=0.8)
    fuzzy = []
    for n in close:
        k = _EXACT[n]
        if k not in fuzzy:
            fuzzy.append(k)
    if len(fuzzy) == 1 and not hits:
        return REACTIONS[fuzzy[0]]
    candidates = sorted(hits) or fuzzy
    if not candidates:
        candidates = [_EXACT[n] for n in get_close_matches(q, [n for n, _ in _INDEX], n=5, cutoff=0.5)]
        candidates = list(dict.fromkeys(candidates))
    hint = ", ".join(candidates[:6]) if candidates else "use list_reactions() to browse"
    raise ReactionError(
        f"Unknown reaction {name_or_alias!r}"
        + (f"; did you mean: {hint}?" if candidates else f" ({hint}).")
    )


def list_reactions(cls: Optional[str] = None, query: Optional[str] = None) -> List[Dict[str, Any]]:
    """Registry entries as dicts, optionally filtered by class and/or free-text query.

    ``cls`` matches the reaction class (case/punctuation-insensitive substring,
    e.g. ``"oxidation"``, ``"C-C coupling"``); ``query`` is matched against
    key, name, class, roles, aliases and description.
    """
    out: List[Dict[str, Any]] = []
    c = _norm(cls) if cls else ""
    q = _norm(query) if query else ""
    for t in REACTIONS.values():
        if c and c not in _norm(t.cls):
            continue
        if q:
            hay = _norm(" ".join([t.key, t.name, t.cls, t.description, *t.roles, *t.aliases]))
            if q not in hay:
                continue
        out.append(t.as_dict())
    return out


# --------------------------------------------------------------------------- #
# Applying templates
# --------------------------------------------------------------------------- #
def _reactant_mol(obj: Any, position: int) -> Any:
    """Parse/copy one reactant into a sanitised, map-number-free RDKit Mol."""
    Chem = _rd().Chem
    if isinstance(obj, str):
        try:
            mol = mol_from_input(obj, prepare=False)
        except ChemError as exc:
            raise ReactionError(f"Reactant {position + 1} ({obj!r}): {exc}") from exc
    elif obj is not None and hasattr(obj, "GetNumAtoms"):
        mol = obj
    else:
        raise ReactionError(
            f"Reactant {position + 1}: expected a SMILES/name/molblock string or an RDKit Mol, "
            f"got {type(obj).__name__}."
        )
    mol = Chem.Mol(mol)
    try:
        with _quiet():
            Chem.SanitizeMol(mol)
    except Exception as exc:
        raise ReactionError(f"Reactant {position + 1} could not be sanitised: {exc}") from exc
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol


def _canon(mol: Any) -> str:
    """Canonical isomeric SMILES without atom-map numbers."""
    Chem = _rd().Chem
    m = Chem.Mol(mol)
    for atom in m.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


@dataclass
class _Outcome:
    """One distinct product set of a template run."""

    order: Tuple[int, ...]              # input reactant index per role
    reactants: List[Any]                # mapped copies, one per role (first equivalent)
    expanded: List[Any]                 # mapped copies, one per reactant template
    products: List[Any]                 # sanitised, mapped product mols
    smiles: Tuple[str, ...]             # canonical product SMILES (map-free)

    @property
    def key(self) -> Tuple[str, ...]:
        return tuple(sorted(self.smiles))

    @property
    def mapped(self) -> str:
        Chem = _rd().Chem
        return ">>".join((
            ".".join(Chem.MolToSmiles(m) for m in self.expanded),
            ".".join(Chem.MolToSmiles(p) for p in self.products),
        ))


_PRODUCT_PROPS = ("react_idx", "react_atom_idx", "old_mapno")


def _drop_spectator_fragments(mol: Any) -> Any:
    """Remove fragments (counter-ions such as Na+) that contribute no mapped atom."""
    Chem = _rd().Chem
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) < 2:
        return mol
    keep = [f for f in frags if any(a.GetAtomMapNum() for a in f.GetAtoms())]
    if not keep or len(keep) == len(frags):
        return mol
    out = keep[0]
    for f in keep[1:]:
        out = Chem.CombineMols(out, f)
    return Chem.Mol(out)


def _run_template(tmpl: ReactionTemplate, role_mols: Sequence[Any], order: Tuple[int, ...],
                  limit: int) -> List[_Outcome]:
    """Run ``tmpl`` on reactants given in role order; returns distinct outcomes."""
    Chem = _rd().Chem
    rxn = tmpl.rxn
    expanded: List[Any] = []
    first_index: List[int] = []
    for role, mol in enumerate(role_mols):
        first_index.append(len(expanded))
        for _ in range(tmpl.equivalents[role]):
            expanded.append(Chem.Mol(mol))
    with _quiet():
        try:
            raw = rxn.RunReactants(tuple(expanded), maxProducts=max(50, 10 * limit))
        except Exception as exc:
            raise ReactionError(f"{tmpl.name}: RDKit could not run the template: {exc}") from exc
    outcomes: List[_Outcome] = []
    seen = set()
    for product_set in raw:
        mols: List[Any] = []
        ok = True
        for p in product_set:
            p = Chem.Mol(p)
            try:
                with _quiet():
                    Chem.SanitizeMol(p)
            except Exception:
                ok = False
                break
            with _quiet():
                frags = Chem.GetMolFrags(p, asMols=True, sanitizeFrags=False)
            mols.extend(Chem.Mol(f) for f in frags)
        if not ok or not mols:
            continue
        # drop exact duplicates inside one outcome (ring-opening gives two copies)
        uniq: List[Any] = []
        uniq_smiles: List[str] = []
        for m in mols:
            s = _canon(m)
            if s not in uniq_smiles:
                uniq.append(m)
                uniq_smiles.append(s)
        key = tuple(sorted(uniq_smiles))
        if key in seen:
            continue
        # atom mapping: number every product atom that came from a reactant
        r_copies = [Chem.Mol(m) for m in expanded]
        assigned: Dict[Tuple[int, int], int] = {}
        counter = 1
        valid = True
        for pm in uniq:
            for atom in pm.GetAtoms():
                if atom.HasProp("react_idx") and atom.HasProp("react_atom_idx"):
                    src = (atom.GetIntProp("react_idx"), atom.GetIntProp("react_atom_idx"))
                    if src in assigned:
                        valid = False      # RDKit duplicated a reactant atom: bogus outcome
                        break
                    assigned[src] = counter
                    atom.SetAtomMapNum(counter)
                    r_copies[src[0]].GetAtomWithIdx(src[1]).SetAtomMapNum(counter)
                    counter += 1
                for prop in _PRODUCT_PROPS:
                    if atom.HasProp(prop):
                        atom.ClearProp(prop)
            if not valid:
                break
        if not valid:
            continue
        seen.add(key)
        r_copies = [_drop_spectator_fragments(m) for m in r_copies]
        outcomes.append(_Outcome(
            order=tuple(order), reactants=[r_copies[i] for i in first_index], expanded=r_copies,
            products=uniq, smiles=tuple(uniq_smiles),
        ))
        if len(outcomes) >= limit:
            break
    return outcomes


def _role_matrix(tmpl: ReactionTemplate, mols: Sequence[Any]) -> List[List[bool]]:
    return [[tmpl.matches_role(m, r) for r in range(tmpl.n_reactants)] for m in mols]


def _orders(tmpl: ReactionTemplate, mols: Sequence[Any], matrix: Sequence[Sequence[bool]]
            ) -> Iterator[Tuple[int, ...]]:
    """Reactant orderings (input index per role) in which every role is satisfied."""
    n = tmpl.n_reactants
    smiles = [_canon(m) for m in mols]
    tried = set()
    for perm in itertools.permutations(range(n)):
        if not all(matrix[perm[r]][r] for r in range(n)):
            continue
        sig = tuple(smiles[i] for i in perm)
        if sig in tried:
            continue
        tried.add(sig)
        yield perm


def _apply(tmpl: ReactionTemplate, mols: Sequence[Any], max_products: int) -> List[_Outcome]:
    if len(mols) != tmpl.n_reactants:
        roles = ", ".join(f"{i + 1}. {r}" for i, r in enumerate(tmpl.roles))
        raise ReactionError(
            f"{tmpl.name} takes {tmpl.n_reactants} reactant(s) ({roles}); got {len(mols)}."
        )
    limit = max(1, int(max_products))
    matrix = _role_matrix(tmpl, mols)
    outcomes: List[_Outcome] = []
    seen = set()
    tried = False
    for order in _orders(tmpl, mols, matrix):
        tried = True
        for oc in _run_template(tmpl, [mols[i] for i in order], order, limit):
            if oc.key in seen:
                continue
            seen.add(oc.key)
            outcomes.append(oc)
            if len(outcomes) >= limit:
                return outcomes
    if outcomes:
        return outcomes
    raise ReactionError(_no_match_message(tmpl, mols, matrix, tried))


def _no_match_message(tmpl: ReactionTemplate, mols: Sequence[Any], matrix: Sequence[Sequence[bool]],
                      tried: bool) -> str:
    smiles = [_canon(m) for m in mols]
    if tried:
        return (
            f"{tmpl.name}: the reactants fit the template roles but RDKit produced no valid "
            f"product (sanitisation failed or the transformation is not applicable to "
            f"{', '.join(repr(s) for s in smiles)})."
        )
    lines = [f"{tmpl.name}: the reactants do not fit the template roles."]
    for r, role in enumerate(tmpl.roles):
        q = _rd().Chem.MolToSmarts(tmpl.role_template(r))
        who = [f"reactant {i + 1} {smiles[i]!r}" for i in range(len(mols)) if matrix[i][r]]
        lines.append(
            f"  role {r + 1} '{role}' ({q}): " + (", ".join(who) if who else "matched by none of the reactants")
        )
    for i, s in enumerate(smiles):
        if not any(matrix[i]):
            lines.append(f"  reactant {i + 1} {s!r} matches no role.")
    return "\n".join(lines)


def _roles_dict(tmpl: ReactionTemplate, mols: Sequence[Any]) -> Dict[str, str]:
    """role name -> reactant SMILES (duplicate role names get a numeric suffix)."""
    out: Dict[str, str] = {}
    for role, m in zip(tmpl.roles, mols):
        label, n = role, 2
        while label in out:
            label, n = f"{role} ({n})", n + 1
        out[label] = _canon(m)
    return out


def _result_dict(tmpl: ReactionTemplate, oc: _Outcome) -> Dict[str, Any]:
    return {
        "template": tmpl.key,
        "name": tmpl.name,
        "class": tmpl.cls,
        "reactant_order": list(oc.order),
        "reactants": [_canon(m) for m in oc.reactants],
        "roles": _roles_dict(tmpl, oc.reactants),
        "products": list(oc.smiles),
        "mapped": oc.mapped,
        "reagents": list(tmpl.reagents),
        "conditions": list(tmpl.conditions),
        "byproducts": list(tmpl.byproducts),
    }


def apply_reaction(name: Union[str, ReactionTemplate], reactants: Sequence[Union[str, Any]],
                   max_products: int = 8) -> List[Dict[str, Any]]:
    """Apply a named template to concrete reactants.

    ``reactants`` are SMILES / names / molblocks (anything ``mol_from_input``
    accepts) or RDKit Mols, one per template role in any order: every
    ordering in which each reactant fits its role pattern is tried.  Products
    are sanitised and deduplicated by canonical SMILES; at most
    ``max_products`` distinct outcomes are returned, each a dict with

    * ``products`` - canonical SMILES of the drawn products,
    * ``mapped`` - atom-mapped reaction SMILES ``reactants>>products`` (the
      map numbers link reactant atoms to product atoms for drawing),
    * ``template``, ``name``, ``class``, ``reagents``, ``conditions``,
      ``byproducts``, ``reactants`` (canonical, in role order),
      ``reactant_order`` (input index per role) and ``roles``.

    Raises :class:`ReactionError` for an unknown template, a wrong reactant
    count, or when no reactant arrangement matches (the message says which
    role each reactant fails).
    """
    tmpl = get_reaction(name)
    if isinstance(reactants, (str, bytes)) or not isinstance(reactants, Sequence):
        raise ReactionError("reactants must be a list of SMILES/names or RDKit Mols.")
    mols = [_reactant_mol(r, i) for i, r in enumerate(reactants)]
    return [_result_dict(tmpl, oc) for oc in _apply(tmpl, mols, max_products)]


def suggest_reactions(reactants: Sequence[Union[str, Any]], max_products: int = 4) -> List[Dict[str, Any]]:
    """Templates that can be applied to exactly these reactants.

    Every template with the same number of roles is tried; those where some
    ordering of the reactants satisfies all role patterns *and* yields a
    sanitisable product are returned (registry order) as dicts with
    ``template``, ``name``, ``class``, ``roles`` (role → reactant SMILES),
    ``products`` (first outcome), ``n_outcomes`` (up to ``max_products``),
    ``reagents`` and ``conditions``.
    """
    if isinstance(reactants, (str, bytes)) or not isinstance(reactants, Sequence):
        raise ReactionError("reactants must be a list of SMILES/names or RDKit Mols.")
    mols = [_reactant_mol(r, i) for i, r in enumerate(reactants)]
    out: List[Dict[str, Any]] = []
    for tmpl in REACTIONS.values():
        if tmpl.n_reactants != len(mols):
            continue
        try:
            matrix = _role_matrix(tmpl, mols)
        except ReactionError:  # pragma: no cover - bad template
            continue
        outcomes: List[_Outcome] = []
        seen = set()
        for order in _orders(tmpl, mols, matrix):
            try:
                found = _run_template(tmpl, [mols[i] for i in order], order, max(1, int(max_products)))
            except ReactionError:
                continue
            for oc in found:
                if oc.key not in seen:
                    seen.add(oc.key)
                    outcomes.append(oc)
            if len(outcomes) >= max_products:
                break
        if not outcomes:
            continue
        first = outcomes[0]
        out.append({
            "template": tmpl.key,
            "name": tmpl.name,
            "class": tmpl.cls,
            "roles": _roles_dict(tmpl, first.reactants),
            "products": list(first.smiles),
            "mapped": first.mapped,
            "n_outcomes": len(outcomes),
            "reagents": list(tmpl.reagents),
            "conditions": list(tmpl.conditions),
        })
    return out


def to_step(name: Union[str, ReactionTemplate], reactants: Sequence[Union[str, Any]],
            product_index: int = 0, include_reagents: bool = True, atom_map: bool = True,
            labels: bool = False) -> "Step":
    """Build a drawable :class:`chemdrawmcp.layout.Step` for one outcome.

    ``product_index`` selects among the distinct outcomes of
    :func:`apply_reaction` (0 = first).  Reactant and product molecules are
    passed through ``prepare_for_depiction`` and keep matching atom-map
    numbers so ``Step.atom_map`` can emit a ReactionStepAtomMap; set
    ``atom_map=False`` to strip them.  With ``include_reagents`` the
    template's reagents go above and its conditions below the arrow; with
    ``labels`` each reactant is captioned with its role name.
    """
    from . import layout  # lazy: pulls in the CDXML engine

    tmpl = get_reaction(name)
    if isinstance(reactants, (str, bytes)) or not isinstance(reactants, Sequence):
        raise ReactionError("reactants must be a list of SMILES/names or RDKit Mols.")
    mols = [_reactant_mol(r, i) for i, r in enumerate(reactants)]
    idx = int(product_index)
    if idx < 0:
        raise ReactionError("product_index must be >= 0.")
    outcomes = _apply(tmpl, mols, max_products=idx + 1)
    if idx >= len(outcomes):
        raise ReactionError(
            f"{tmpl.name}: only {len(outcomes)} distinct outcome(s) for these reactants; "
            f"product_index {idx} is out of range."
        )
    oc = outcomes[idx]

    def prep(m: Any) -> Any:
        m = _rd().Chem.Mol(m)
        if not atom_map:
            for atom in m.GetAtoms():
                atom.SetAtomMapNum(0)
        return prepare_for_depiction(m)

    r_comps = [
        layout.Component(prep(m), label=(tmpl.roles[i] if labels else None))
        for i, m in enumerate(oc.reactants)
    ]
    p_comps = [layout.Component(prep(p)) for p in oc.products]
    above = list(tmpl.reagents) if include_reagents else []
    below = list(tmpl.conditions) if include_reagents else []
    return layout.Step(r_comps, p_comps, above, below, "full", atom_map=atom_map)


def describe(name: Union[str, ReactionTemplate]) -> str:
    """Human-readable multi-line summary of a template."""
    t = get_reaction(name)
    lines = [f"{t.name} [{t.key}] — class: {t.cls}", f"  {t.description}"]
    eq = t.equivalents
    roles = []
    for i, role in enumerate(t.roles):
        roles.append(f"{i + 1}. {role}" + (f" ({eq[i]} equiv)" if eq[i] != 1 else ""))
    lines.append("  reactants: " + "; ".join(roles))
    if t.reagents:
        lines.append("  reagents (above arrow): " + ", ".join(t.reagents))
    if t.conditions:
        lines.append("  conditions (below arrow): " + ", ".join(t.conditions))
    if t.byproducts:
        lines.append("  by-products (not drawn): " + ", ".join(t.byproducts))
    lines.append("  SMARTS: " + t.smarts)
    if t.example_reactants:
        lines.append(
            "  example: " + " + ".join(t.example_reactants) + " -> " + " + ".join(t.example_products)
        )
    if t.aliases:
        lines.append("  aliases: " + ", ".join(t.aliases))
    return "\n".join(lines)
