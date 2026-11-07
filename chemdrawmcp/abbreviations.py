"""Chemical abbreviations used in ChemDraw labels, and a small atom-label parser.

ChemDraw stores nicknames ("OMe", "Boc", "Ph" ...) as ``<n NodeType="Nickname">`` nodes that normally carry an
expanded ``<fragment>`` child.  Files written by other tools (or labels ChemDraw could not interpret) have only
the text, so the reader falls back to the tables here.

* :data:`SUBSTITUENTS` map a label to a SMILES whose **first atom is the attachment point**.
* :data:`MOLECULES` map reagent / solvent abbreviations to complete molecules (no attachment point).
* :func:`parse_atom_label` derives element, hydrogen count, charge and isotope from a plain atom label such as
  ``"NH3+"``, ``"O-"``, ``"13C"`` or ``"D"``.
"""
from __future__ import annotations

import re
from typing import Dict, Optional, Tuple

# ---------------------------------------------------------------------------
# substituents: first atom of the SMILES is the atom bonded to the rest of the molecule
# ---------------------------------------------------------------------------

SUBSTITUENTS: Dict[str, str] = {
    # alkyl / aryl
    "Me": "C", "CH3": "C", "Et": "CC", "CH2CH3": "CC", "Pr": "CCC", "nPr": "CCC", "n-Pr": "CCC", "iPr": "C(C)C",
    "i-Pr": "C(C)C", "Bu": "CCCC", "nBu": "CCCC", "n-Bu": "CCCC", "iBu": "CC(C)C", "i-Bu": "CC(C)C", "sBu": "C(C)CC",
    "s-Bu": "C(C)CC", "tBu": "C(C)(C)C", "t-Bu": "C(C)(C)C", "Pent": "CCCCC", "Am": "CCCCC", "Hex": "CCCCCC",
    "Oct": "CCCCCCCC", "Cy": "C1CCCCC1", "cHex": "C1CCCCC1", "cPr": "C1CC1", "c-Pr": "C1CC1", "cBu": "C1CCC1",
    "cPent": "C1CCCC1", "Allyl": "CC=C", "Vinyl": "C=C", "Propargyl": "CC#C", "Bn": "Cc1ccccc1", "CH2Ph": "Cc1ccccc1",
    "Ph": "c1ccccc1", "C6H5": "c1ccccc1", "Tol": "c1ccc(C)cc1", "pTol": "c1ccc(C)cc1", "p-Tol": "c1ccc(C)cc1",
    "Mes": "c1c(C)cc(C)cc1C", "Xyl": "c1c(C)cccc1C", "Np": "c1cccc2ccccc12", "An": "c1ccc(OC)cc1",
    "Fc": "[Fe]", "Cp": "C1C=CC=C1",
    # acyl / carbamate protecting groups
    "Ac": "C(C)=O", "COCH3": "C(C)=O", "Bz": "C(=O)c1ccccc1", "Piv": "C(=O)C(C)(C)C", "Boc": "C(=O)OC(C)(C)C",
    "Fmoc": "C(=O)OCC1c2ccccc2-c2ccccc21", "Cbz": "C(=O)OCc1ccccc1",
    "Alloc": "C(=O)OCC=C", "Troc": "C(=O)OCC(Cl)(Cl)Cl", "Teoc": "C(=O)OCC[Si](C)(C)C", "Moz": "C(=O)OCc1ccc(OC)cc1",
    "Tfa": "C(=O)C(F)(F)F", "COCF3": "C(=O)C(F)(F)F", "Fmoc-": "C(=O)OCC1c2ccccc2-c2ccccc21",
    # sulfonyl
    "Ts": "S(=O)(=O)c1ccc(C)cc1", "Tos": "S(=O)(=O)c1ccc(C)cc1", "Ms": "S(=O)(=O)C", "Tf": "S(=O)(=O)C(F)(F)F",
    "Ns": "S(=O)(=O)c1ccc([N+](=O)[O-])cc1", "Bs": "S(=O)(=O)c1ccccc1", "SES": "S(=O)(=O)CC[Si](C)(C)C",
    "SO2Me": "S(=O)(=O)C", "SO2Ph": "S(=O)(=O)c1ccccc1", "SO2NH2": "S(=O)(=O)N", "SO2Cl": "S(=O)(=O)Cl",
    "SO3H": "S(=O)(=O)O", "SO3Na": "S(=O)(=O)[O-].[Na+]", "SOMe": "S(C)=O", "SO2CF3": "S(=O)(=O)C(F)(F)F",
    # silyl
    "TMS": "[Si](C)(C)C", "SiMe3": "[Si](C)(C)C", "TES": "[Si](CC)(CC)CC", "TBS": "[Si](C)(C)C(C)(C)C",
    "TBDMS": "[Si](C)(C)C(C)(C)C", "TIPS": "[Si](C(C)C)(C(C)C)C(C)C", "TBDPS": "[Si](c1ccccc1)(c1ccccc1)C(C)(C)C",
    "TPS": "[Si](c1ccccc1)(c1ccccc1)c1ccccc1",
    # ethers / acetals
    "PMB": "Cc1ccc(OC)cc1", "MPM": "Cc1ccc(OC)cc1", "MOM": "COC", "MEM": "COCCOC", "BOM": "COCc1ccccc1",
    "SEM": "COCC[Si](C)(C)C", "THP": "C1CCCCO1", "Tr": "C(c1ccccc1)(c1ccccc1)c1ccccc1", "Trt": "C(c1ccccc1)(c1ccccc1)c1ccccc1",
    "DMT": "C(c1ccccc1)(c1ccc(OC)cc1)c1ccc(OC)cc1", "Dmb": "Cc1ccc(OC)cc1OC",
    # O-linked
    "OH": "O", "OMe": "OC", "OCH3": "OC", "OEt": "OCC", "OPr": "OCCC", "OiPr": "OC(C)C", "OnBu": "OCCCC", "OBu": "OCCCC",
    "OtBu": "OC(C)(C)C", "OBn": "OCc1ccccc1", "OPh": "Oc1ccccc1", "OAc": "OC(C)=O", "OBz": "OC(=O)c1ccccc1",
    "OPiv": "OC(=O)C(C)(C)C", "OTf": "OS(=O)(=O)C(F)(F)F", "OTs": "OS(=O)(=O)c1ccc(C)cc1", "OMs": "OS(=O)(=O)C",
    "ONs": "OS(=O)(=O)c1ccc([N+](=O)[O-])cc1", "OTMS": "O[Si](C)(C)C", "OTES": "O[Si](CC)(CC)CC",
    "OTBS": "O[Si](C)(C)C(C)(C)C", "OTBDMS": "O[Si](C)(C)C(C)(C)C", "OTIPS": "O[Si](C(C)C)(C(C)C)C(C)C",
    "OTBDPS": "O[Si](c1ccccc1)(c1ccccc1)C(C)(C)C", "OTHP": "OC1CCCCO1", "OMOM": "OCOC", "OMEM": "OCOCCOC",
    "OPMB": "OCc1ccc(OC)cc1", "OBoc": "OC(=O)OC(C)(C)C", "OCbz": "OC(=O)OCc1ccccc1", "OCF3": "OC(F)(F)F",
    "OOH": "OO", "OTr": "OC(c1ccccc1)(c1ccccc1)c1ccccc1", "OAllyl": "OCC=C", "OSO3H": "OS(=O)(=O)O",
    "OPO3H2": "OP(=O)(O)O", "ONO2": "O[N+](=O)[O-]", "OCOPh": "OC(=O)c1ccccc1", "OCOMe": "OC(C)=O",
    "ONa": "[O-].[Na+]", "OK": "[O-].[K+]", "OLi": "[O-].[Li+]",
    # S-linked
    "SH": "S", "SMe": "SC", "SEt": "SCC", "SPh": "Sc1ccccc1", "SBn": "SCc1ccccc1", "SAc": "SC(C)=O", "StBu": "SC(C)(C)C",
    "SCF3": "SC(F)(F)F", "SCN": "SC#N", "STrt": "SC(c1ccccc1)(c1ccccc1)c1ccccc1",
    # N-linked
    "NH2": "N", "NHMe": "NC", "NMe2": "N(C)C", "NEt2": "N(CC)CC", "NHEt": "NCC", "NiPr2": "N(C(C)C)C(C)C",
    "NHAc": "NC(C)=O", "NHBoc": "NC(=O)OC(C)(C)C", "NHFmoc": "NC(=O)OCC1c2ccccc2-c2ccccc21", "NHCbz": "NC(=O)OCc1ccccc1",
    "NHZ": "NC(=O)OCc1ccccc1", "NHTs": "NS(=O)(=O)c1ccc(C)cc1", "NHMs": "NS(=O)(=O)C", "NHBn": "NCc1ccccc1",
    "NHPh": "Nc1ccccc1", "NHBz": "NC(=O)c1ccccc1", "NHTf": "NS(=O)(=O)C(F)(F)F", "NBoc2": "N(C(=O)OC(C)(C)C)C(=O)OC(C)(C)C",
    "NPhth": "N1C(=O)c2ccccc2C1=O", "NHOH": "NO", "NHNH2": "NN", "NHOMe": "NOC", "NMe3+": "[N+](C)(C)C",
    "N3": "N=[N+]=[N-]", "NO2": "[N+](=O)[O-]", "NO": "N=O", "NC": "[N+]#[C-]", "NCO": "N=C=O", "NCS": "N=C=S",
    "NH3+": "[NH3+]", "NHNHBoc": "NNC(=O)OC(C)(C)C", "NHCOMe": "NC(C)=O", "NHCO2Me": "NC(=O)OC",
    "N(Me)Boc": "N(C)C(=O)OC(C)(C)C", "NMeBoc": "N(C)C(=O)OC(C)(C)C",
    # C-linked functional groups
    "CN": "C#N", "CHO": "C=O", "COOH": "C(=O)O", "CO2H": "C(=O)O", "COOMe": "C(=O)OC", "CO2Me": "C(=O)OC",
    "COOEt": "C(=O)OCC", "CO2Et": "C(=O)OCC", "CO2tBu": "C(=O)OC(C)(C)C", "COOtBu": "C(=O)OC(C)(C)C",
    "CO2Bn": "C(=O)OCc1ccccc1", "CO2iPr": "C(=O)OC(C)C", "CO2Ph": "C(=O)Oc1ccccc1", "CO2-": "C(=O)[O-]",
    "COO-": "C(=O)[O-]", "CO2Na": "C(=O)[O-].[Na+]", "COMe": "C(C)=O", "COEt": "C(=O)CC", "COPh": "C(=O)c1ccccc1",
    "CONH2": "C(N)=O", "CONHMe": "C(=O)NC", "CONMe2": "C(=O)N(C)C", "CONHPh": "C(=O)Nc1ccccc1", "COCl": "C(=O)Cl",
    "COBr": "C(=O)Br", "COSMe": "C(=O)SC", "CF3": "C(F)(F)F", "CCl3": "C(Cl)(Cl)Cl", "CBr3": "C(Br)(Br)Br",
    "CHF2": "C(F)F", "CH2F": "CF", "CH2Cl": "CCl", "CH2Br": "CBr", "CH2I": "CI", "CH2OH": "CO", "CH2OMe": "COC",
    "CH2OAc": "COC(C)=O", "CH2OBn": "COCc1ccccc1", "CH2NH2": "CN", "CH2NMe2": "CN(C)C", "CH2CO2H": "CC(=O)O",
    "CH2CO2Me": "CC(=O)OC", "CH2CO2Et": "CC(=O)OCC", "CH2CN": "CC#N", "CH2SH": "CS", "CH2OTBS": "CO[Si](C)(C)C(C)(C)C",
    "CH2CH2OH": "CCO", "CH2CH2NH2": "CCN", "C2H5": "CC", "C3H7": "CCC", "C4H9": "CCCC", "CH=CH2": "C=C", "C≡CH": "C#C",
    "CCH": "C#C", "CH(OH)": "CO", "C(O)": "C=O", "CHOH": "CO", "CH2": "C", "CH": "C", "C(CH3)3": "C(C)(C)C", "CH(CH3)2": "C(C)C",
    # P / B / Sn / other
    "PPh2": "P(c1ccccc1)c1ccccc1", "PPh3": "P(c1ccccc1)(c1ccccc1)c1ccccc1", "PPh3+": "[P+](c1ccccc1)(c1ccccc1)c1ccccc1",
    "PMe2": "P(C)C", "PCy2": "P(C1CCCCC1)C1CCCCC1", "PtBu2": "P(C(C)(C)C)C(C)(C)C", "P(O)(OEt)2": "P(=O)(OCC)OCC",
    "PO(OEt)2": "P(=O)(OCC)OCC", "PO3H2": "P(=O)(O)O", "P(OEt)2": "P(OCC)OCC", "B(OH)2": "B(O)O", "Bpin": "B1OC(C)(C)C(C)(C)O1",
    "BPin": "B1OC(C)(C)C(C)(C)O1", "B(pin)": "B1OC(C)(C)C(C)(C)O1", "BF3K": "[B-](F)(F)F.[K+]", "Bcat": "B1Oc2ccccc2O1",
    "SnBu3": "[Sn](CCCC)(CCCC)CCCC", "SnMe3": "[Sn](C)(C)C", "SiMe2Ph": "[Si](C)(C)c1ccccc1", "SiPh3": "[Si](c1ccccc1)(c1ccccc1)c1ccccc1",
    "SeMe": "[Se]C", "SePh": "[Se]c1ccccc1", "MgBr": "[Mg]Br", "MgCl": "[Mg]Cl", "MgI": "[Mg]I", "ZnBr": "[Zn]Br",
    "ZnCl": "[Zn]Cl", "Li": "[Li]", "Na": "[Na]", "K": "[K]",
    # ring-ish
    "Ind": "c1cc2ccccc2[nH]1", "Pyr": "c1ccncc1", "Py": "c1ccncc1", "Fur": "c1ccco1", "Thi": "c1cccs1",
    "Im": "N1C=CN=C1", "Pip": "N1CCCCC1", "Mor": "N1CCOCC1", "Pyrr": "N1CCCC1",
    # terminators
    "D": "[2H]", "T": "[3H]", "H": "[H]", "OD": "O[2H]", "CD3": "C([2H])([2H])[2H]", "ND2": "N([2H])[2H]",
}

# ---------------------------------------------------------------------------
# whole molecules (reagents, solvents, catalysts): no attachment point
# ---------------------------------------------------------------------------

MOLECULES: Dict[str, str] = {
    # small inorganics
    "H2O": "O", "D2O": "[2H]O[2H]", "H2": "[H][H]", "D2": "[2H][2H]", "O2": "O=O", "N2": "N#N", "O3": "[O-][O+]=O",
    "CO2": "O=C=O", "CO": "[C-]#[O+]", "H2O2": "OO", "H2S": "S", "NH3": "N", "HCl": "Cl", "HBr": "Br", "HI": "I", "HF": "F",
    "H2SO4": "OS(=O)(=O)O", "HNO3": "O[N+](=O)[O-]", "H3PO4": "OP(=O)(O)O", "HClO4": "OCl(=O)(=O)=O", "HCO2H": "OC=O",
    "HCOOH": "OC=O", "HCHO": "C=O", "CH2O": "C=O", "SO2": "O=S=O", "SO3": "O=S(=O)=O", "NO": "[N]=O", "N2O": "[N-]=[N+]=O",
    "Cl2": "ClCl", "Br2": "BrBr", "I2": "II", "F2": "FF", "NH4Cl": "[NH4+].[Cl-]", "NH4OH": "[NH4+].[OH-]",
    "NaOH": "[Na+].[OH-]", "KOH": "[K+].[OH-]", "LiOH": "[Li+].[OH-]", "CsOH": "[Cs+].[OH-]", "NaH": "[Na+].[H-]",
    "KH": "[K+].[H-]", "LiH": "[Li+].[H-]", "CaH2": "[Ca+2].[H-].[H-]", "NaCl": "[Na+].[Cl-]", "KCl": "[K+].[Cl-]",
    "LiCl": "[Li+].[Cl-]", "NaBr": "[Na+].[Br-]", "KBr": "[K+].[Br-]", "LiBr": "[Li+].[Br-]", "NaI": "[Na+].[I-]",
    "KI": "[K+].[I-]", "NaF": "[Na+].[F-]", "KF": "[K+].[F-]", "CsF": "[Cs+].[F-]", "NaN3": "[Na+].[N-]=[N+]=[N-]",
    "NaCN": "[Na+].[C-]#N", "KCN": "[K+].[C-]#N", "NaNO2": "[Na+].[O-]N=O", "NaOMe": "[Na+].[O-]C", "NaOEt": "[Na+].[O-]CC",
    "KOtBu": "[K+].[O-]C(C)(C)C", "NaOtBu": "[Na+].[O-]C(C)(C)C", "NaOAc": "[Na+].CC(=O)[O-]", "KOAc": "[K+].CC(=O)[O-]",
    "NaHCO3": "[Na+].OC(=O)[O-]", "Na2CO3": "[Na+].[Na+].[O-]C(=O)[O-]", "K2CO3": "[K+].[K+].[O-]C(=O)[O-]",
    "Cs2CO3": "[Cs+].[Cs+].[O-]C(=O)[O-]", "K3PO4": "[K+].[K+].[K+].[O-]P(=O)([O-])[O-]", "NaNH2": "[Na+].[NH2-]",
    "KHSO4": "[K+].OS(=O)(=O)[O-]", "Na2SO4": "[Na+].[Na+].[O-]S(=O)(=O)[O-]", "MgSO4": "[Mg+2].[O-]S(=O)(=O)[O-]",
    "CaCl2": "[Ca+2].[Cl-].[Cl-]", "NaIO4": "[Na+].[O-]I(=O)(=O)=O", "NaOCl": "[Na+].[O-]Cl", "KMnO4": "[K+].[O-][Mn](=O)(=O)=O",
    "OsO4": "O=[Os](=O)(=O)=O", "SeO2": "O=[Se]=O", "CrO3": "O=[Cr](=O)=O", "MnO2": "O=[Mn]=O", "PtO2": "O=[Pt]=O",
    "Na2S2O3": "[Na+].[Na+].[O-]S(=S)(=O)[O-]", "Na2SO3": "[Na+].[Na+].[O-]S(=O)[O-]", "NaHSO3": "[Na+].OS(=O)[O-]",
    "LiClO4": "[Li+].[O-]Cl(=O)(=O)=O", "AgNO3": "[Ag+].[O-][N+](=O)[O-]", "AgOTf": "[Ag+].[O-]S(=O)(=O)C(F)(F)F",
    "Ag2O": "[Ag+].[Ag+].[O-2]", "Ag2CO3": "[Ag+].[Ag+].[O-]C(=O)[O-]", "CuSO4": "[Cu+2].[O-]S(=O)(=O)[O-]",
    "ZnCl2": "Cl[Zn]Cl", "ZnBr2": "Br[Zn]Br", "ZnI2": "I[Zn]I", "MgBr2": "Br[Mg]Br", "MgCl2": "Cl[Mg]Cl",
    "AlCl3": "Cl[Al](Cl)Cl", "AlBr3": "Br[Al](Br)Br", "BF3": "FB(F)F", "BCl3": "ClB(Cl)Cl", "BBr3": "BrB(Br)Br",
    "BF3·OEt2": "FB(F)F.CCOCC", "BF3.OEt2": "FB(F)F.CCOCC", "BF3·Et2O": "FB(F)F.CCOCC", "BF3-OEt2": "FB(F)F.CCOCC",
    "TiCl4": "Cl[Ti](Cl)(Cl)Cl", "Ti(OiPr)4": "CC(C)O[Ti](OC(C)C)(OC(C)C)OC(C)C", "SnCl2": "Cl[Sn]Cl", "SnCl4": "Cl[Sn](Cl)(Cl)Cl",
    "FeCl3": "Cl[Fe](Cl)Cl", "FeCl2": "Cl[Fe]Cl", "FeBr3": "Br[Fe](Br)Br", "CuI": "[Cu]I", "CuBr": "[Cu]Br", "CuCl": "[Cu]Cl",
    "CuCN": "[Cu]C#N", "Cu(OAc)2": "CC(=O)[O-].CC(=O)[O-].[Cu+2]", "CuCl2": "Cl[Cu]Cl", "CuBr2": "Br[Cu]Br",
    "Pd(OAc)2": "CC(=O)[O-].CC(=O)[O-].[Pd+2]", "PdCl2": "Cl[Pd]Cl", "Pd/C": "[Pd]", "Pd(OH)2": "[OH-].[OH-].[Pd+2]",
    "Pd(PPh3)4": "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "PdCl2(PPh3)2": "Cl[Pd]Cl.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "Pd(dppf)Cl2": "Cl[Pd]Cl.[Fe].c1cc[c-](P(c2ccccc2)c2ccccc2)c1.c1cc[c-](P(c2ccccc2)c2ccccc2)c1",
    "Ni(COD)2": "[Ni].C1=CCCC=CCC1.C1=CCCC=CCC1", "NiCl2": "Cl[Ni]Cl", "RuCl3": "Cl[Ru](Cl)Cl", "RhCl3": "Cl[Rh](Cl)Cl",
    "PtCl2": "Cl[Pt]Cl", "HgCl2": "Cl[Hg]Cl", "Hg(OAc)2": "CC(=O)[O-].CC(=O)[O-].[Hg+2]", "Zn": "[Zn]", "Mg": "[Mg]",
    "Fe": "[Fe]", "Na": "[Na]", "K": "[K]", "Li": "[Li]", "Cu": "[Cu]", "Pt": "[Pt]", "Pd": "[Pd]", "Ni": "[Ni]", "Rh": "[Rh]",
    "Ru": "[Ru]", "Sn": "[Sn]", "Zn/Cu": "[Zn].[Cu]", "Zn(Cu)": "[Zn].[Cu]", "Zn/Hg": "[Zn].[Hg]", "Na/Hg": "[Na].[Hg]",
    "Na/NH3": "[Na].N", "Li/NH3": "[Li].N", "Mg/MeOH": "[Mg].CO", "H2/Pd": "[H][H].[Pd]", "H2/Pd-C": "[H][H].[Pd]",
    "H2/Pd/C": "[H][H].[Pd]", "H2, Pd/C": "[H][H].[Pd]", "H2/PtO2": "[H][H].O=[Pt]=O", "H2/Raney Ni": "[H][H].[Ni]",
    "Raney Ni": "[Ni]", "Lindlar": "[Pd].[Ca+2].[O-]C(=O)[O-].[Pb]",
    # hydrides and organometallics
    "NaBH4": "[Na+].[BH4-]", "LiBH4": "[Li+].[BH4-]", "KBH4": "[K+].[BH4-]", "LiAlH4": "[Li+].[AlH4-]", "LAH": "[Li+].[AlH4-]",
    "NaBH3CN": "[Na+].[BH3-]C#N", "NaBH(OAc)3": "[Na+].[BH-](OC(C)=O)(OC(C)=O)OC(C)=O", "Red-Al": "[Na+].COCCO[AlH2-]OCCOC",
    "DIBAL": "CC(C)C[AlH]CC(C)C", "DIBAL-H": "CC(C)C[AlH]CC(C)C", "DIBALH": "CC(C)C[AlH]CC(C)C", "BH3": "B",
    "BH3·THF": "B.C1CCOC1", "BH3·SMe2": "B.CSC", "BH3.THF": "B.C1CCOC1", "9-BBN": "B1C2CCCC1CCC2", "Et3SiH": "CC[SiH](CC)CC",
    "Bu3SnH": "CCCC[SnH](CCCC)CCCC", "Ph3SnH": "c1ccc([SnH](c2ccccc2)c2ccccc2)cc1", "HSiEt3": "CC[SiH](CC)CC",
    "LDA": "[Li+].CC(C)[N-]C(C)C", "LiHMDS": "[Li+].C[Si](C)(C)[N-][Si](C)(C)C", "NaHMDS": "[Na+].C[Si](C)(C)[N-][Si](C)(C)C",
    "KHMDS": "[K+].C[Si](C)(C)[N-][Si](C)(C)C", "LiTMP": "[Li+].CC1(C)CCCC(C)(C)[N-]1", "nBuLi": "[Li]CCCC", "n-BuLi": "[Li]CCCC",
    "BuLi": "[Li]CCCC", "sBuLi": "[Li]C(C)CC", "s-BuLi": "[Li]C(C)CC", "tBuLi": "[Li]C(C)(C)C", "t-BuLi": "[Li]C(C)(C)C",
    "MeLi": "[Li]C", "PhLi": "[Li]c1ccccc1", "MeMgBr": "C[Mg]Br", "CH3MgBr": "C[Mg]Br", "MeMgCl": "C[Mg]Cl", "MeMgI": "C[Mg]I",
    "EtMgBr": "CC[Mg]Br", "PhMgBr": "Br[Mg]c1ccccc1", "PhMgCl": "Cl[Mg]c1ccccc1", "iPrMgCl": "CC(C)[Mg]Cl",
    "Me2Zn": "C[Zn]C", "Et2Zn": "CC[Zn]CC", "Me3Al": "C[Al](C)C", "AlMe3": "C[Al](C)C", "Et3Al": "CC[Al](CC)CC",
    "Me2CuLi": "[Li+].C[Cu-]C", "Me3SiCl": "C[Si](C)(C)Cl", "TMSCl": "C[Si](C)(C)Cl", "TMSCN": "C[Si](C)(C)C#N",
    "TMSOTf": "C[Si](C)(C)OS(=O)(=O)C(F)(F)F", "TMSN3": "C[Si](C)(C)N=[N+]=[N-]", "TMSCHN2": "C[Si](C)(C)[CH-][N+]#N",
    "TMSI": "C[Si](C)(C)I", "TMSBr": "C[Si](C)(C)Br", "TESCl": "CC[Si](CC)(CC)Cl", "TBSCl": "C[Si](C)(Cl)C(C)(C)C",
    "TBDMSCl": "C[Si](C)(Cl)C(C)(C)C", "TBSOTf": "C[Si](C)(OS(=O)(=O)C(F)(F)F)C(C)(C)C", "TIPSCl": "CC(C)[Si](Cl)(C(C)C)C(C)C",
    "TIPSOTf": "CC(C)[Si](OS(=O)(=O)C(F)(F)F)(C(C)C)C(C)C", "TBDPSCl": "CC(C)(C)[Si](Cl)(c1ccccc1)c1ccccc1",
    "HMDS": "C[Si](C)(C)N[Si](C)(C)C", "TBAF": "CCCC[N+](CCCC)(CCCC)CCCC.[F-]", "TBAB": "CCCC[N+](CCCC)(CCCC)CCCC.[Br-]",
    "TBAI": "CCCC[N+](CCCC)(CCCC)CCCC.[I-]", "Bu4NF": "CCCC[N+](CCCC)(CCCC)CCCC.[F-]",
    # bases / amines
    "Et3N": "CCN(CC)CC", "NEt3": "CCN(CC)CC", "TEA": "CCN(CC)CC", "DIPEA": "CCN(C(C)C)C(C)C", "DIEA": "CCN(C(C)C)C(C)C",
    "iPr2NEt": "CCN(C(C)C)C(C)C", "Hünig's base": "CCN(C(C)C)C(C)C", "DBU": "C1CCC2=NCCCN2CC1", "DBN": "C1CN=C2CCCN2C1",
    "DABCO": "C1CN2CCN1CC2", "DMAP": "CN(C)c1ccncc1", "pyridine": "c1ccncc1", "Pyridine": "c1ccncc1", "Py": "c1ccncc1",
    "pyr": "c1ccncc1", "2,6-lutidine": "Cc1cccc(C)n1", "lutidine": "Cc1cccc(C)n1", "collidine": "Cc1cc(C)nc(C)c1",
    "imidazole": "c1c[nH]cn1", "Imidazole": "c1c[nH]cn1", "NMM": "CN1CCOCC1", "TMEDA": "CN(C)CCN(C)C", "morpholine": "C1COCCN1",
    "piperidine": "C1CCNCC1", "pyrrolidine": "C1CCNC1", "proton sponge": "CN(C)c1cccc2cccc(N(C)C)c12", "Et2NH": "CCNCC",
    "iPr2NH": "CC(C)NC(C)C", "NH2NH2": "NN", "N2H4": "NN", "NH2OH": "NO", "NH2OH·HCl": "NO.Cl", "MeNH2": "CN", "Me2NH": "CNC",
    "NaOH aq": "[Na+].[OH-].O", "aq. NaOH": "[Na+].[OH-].O",
    # solvents
    "DCM": "ClCCl", "CH2Cl2": "ClCCl", "CHCl3": "ClC(Cl)Cl", "CCl4": "ClC(Cl)(Cl)Cl", "DCE": "ClCCCl", "THF": "C1CCOC1",
    "2-MeTHF": "CC1CCCO1", "DMF": "CN(C)C=O", "DMSO": "CS(C)=O", "DMA": "CC(=O)N(C)C", "DMAc": "CC(=O)N(C)C", "NMP": "CN1CCCC1=O",
    "MeCN": "CC#N", "ACN": "CC#N", "CH3CN": "CC#N", "MeOH": "CO", "CH3OH": "CO", "EtOH": "CCO", "iPrOH": "CC(C)O", "i-PrOH": "CC(C)O",
    "IPA": "CC(C)O", "nBuOH": "CCCCO", "tBuOH": "CC(C)(C)O", "t-BuOH": "CC(C)(C)O", "AcOH": "CC(=O)O", "HOAc": "CC(=O)O",
    "CH3COOH": "CC(=O)O", "Ac2O": "CC(=O)OC(C)=O", "TFA": "OC(=O)C(F)(F)F", "TFAA": "O=C(OC(=O)C(F)(F)F)C(F)(F)F",
    "TFE": "OCC(F)(F)F", "HFIP": "OC(C(F)(F)F)C(F)(F)F", "EtOAc": "CCOC(C)=O", "AcOEt": "CCOC(C)=O", "Et2O": "CCOCC",
    "ether": "CCOCC", "MTBE": "COC(C)(C)C", "TBME": "COC(C)(C)C", "dioxane": "C1COCCO1", "1,4-dioxane": "C1COCCO1",
    "toluene": "Cc1ccccc1", "Toluene": "Cc1ccccc1", "PhMe": "Cc1ccccc1", "PhH": "c1ccccc1", "benzene": "c1ccccc1",
    "xylene": "Cc1ccccc1C", "PhCl": "Clc1ccccc1", "PhCF3": "FC(F)(F)c1ccccc1", "hexane": "CCCCCC", "hexanes": "CCCCCC",
    "heptane": "CCCCCCC", "pentane": "CCCCC", "cyclohexane": "C1CCCCC1", "DME": "COCCOC", "glyme": "COCCOC",
    "diglyme": "COCCOCCOC", "HMPA": "CN(C)P(=O)(N(C)C)N(C)C", "DMPU": "CN1CCCN(C)C1=O", "acetone": "CC(C)=O",
    "Acetone": "CC(C)=O", "MEK": "CCC(C)=O", "ethylene glycol": "OCCO", "PEG": "OCCO", "H2O/THF": "O.C1CCOC1",
    "THF/H2O": "C1CCOC1.O", "MeOH/H2O": "CO.O", "DMSO-H2O": "CS(C)=O.O", "DMSO/H2O": "CS(C)=O.O", "EtOH/H2O": "CCO.O",
    # electrophiles / reagents
    "MsCl": "CS(=O)(=O)Cl", "TsCl": "Cc1ccc(S(=O)(=O)Cl)cc1", "Tf2O": "O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F",
    "Ms2O": "CS(=O)(=O)OS(C)(=O)=O", "Tf2NPh": "O=S(=O)(N(S(=O)(=O)C(F)(F)F)c1ccccc1)C(F)(F)F", "PhNTf2": "O=S(=O)(N(S(=O)(=O)C(F)(F)F)c1ccccc1)C(F)(F)F",
    "SOCl2": "O=S(Cl)Cl", "SO2Cl2": "O=S(=O)(Cl)Cl", "POCl3": "O=P(Cl)(Cl)Cl", "PCl3": "ClP(Cl)Cl", "PCl5": "ClP(Cl)(Cl)(Cl)Cl",
    "PBr3": "BrP(Br)Br", "(COCl)2": "O=C(Cl)C(Cl)=O", "oxalyl chloride": "O=C(Cl)C(Cl)=O",
    "AcCl": "CC(=O)Cl", "BzCl": "O=C(Cl)c1ccccc1", "PivCl": "CC(C)(C)C(=O)Cl", "Boc2O": "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C",
    "(Boc)2O": "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C", "FmocCl": "O=C(Cl)OCC1c2ccccc2-c2ccccc21", "Fmoc-Cl": "O=C(Cl)OCC1c2ccccc2-c2ccccc21",
    "FmocOSu": "O=C(ON1C(=O)CCC1=O)OCC1c2ccccc2-c2ccccc21", "Fmoc-OSu": "O=C(ON1C(=O)CCC1=O)OCC1c2ccccc2-c2ccccc21",
    "CbzCl": "O=C(Cl)OCc1ccccc1", "Cbz-Cl": "O=C(Cl)OCc1ccccc1", "ZCl": "O=C(Cl)OCc1ccccc1", "ClCO2Me": "COC(=O)Cl",
    "ClCO2Et": "CCOC(=O)Cl", "MOMCl": "COCCl", "MEMCl": "COCCOCCl", "SEMCl": "C[Si](C)(C)CCOCCl", "BnBr": "BrCc1ccccc1",
    "BnCl": "ClCc1ccccc1", "PMBCl": "COc1ccc(CCl)cc1", "PMBBr": "COc1ccc(CBr)cc1", "MeI": "CI", "CH3I": "CI", "EtI": "CCI",
    "MeBr": "CBr", "EtBr": "CCBr", "allyl bromide": "C=CCBr", "AllylBr": "C=CCBr", "propargyl bromide": "C#CCBr",
    "MeOTf": "COS(=O)(=O)C(F)(F)F", "Me2SO4": "COS(=O)(=O)OC", "Et2SO4": "CCOS(=O)(=O)OCC", "Me3OBF4": "C[O+](C)C.F[B-](F)(F)F",
    "Meerwein's salt": "CC[O+](CC)CC.F[B-](F)(F)F", "CH2N2": "C=[N+]=[N-]", "CH2I2": "ICI", "CH2Br2": "BrCBr", "CHI3": "IC(I)I",
    "CBr4": "BrC(Br)(Br)Br", "NBS": "O=C1CCC(=O)N1Br", "NCS": "O=C1CCC(=O)N1Cl", "NIS": "O=C1CCC(=O)N1I",
    "DBDMH": "CC1(C)N(Br)C(=O)N(Br)C1=O", "Selectfluor": "[B-](F)(F)(F)F.[B-](F)(F)(F)F.F[N+]12CC[N+](CCl)(CC1)CC2",
    "NFSI": "O=S(=O)(N(F)S(=O)(=O)c1ccccc1)c1ccccc1", "DAST": "CCN(CC)S(F)(F)F", "Deoxo-Fluor": "COCCN(CCOC)S(F)(F)F",
    "AIBN": "CC(C)(C#N)N=NC(C)(C)C#N", "BPO": "O=C(OOC(=O)c1ccccc1)c1ccccc1", "PPh3": "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "Ph3P": "c1ccc(P(c2ccccc2)c2ccccc2)cc1", "PBu3": "CCCCP(CCCC)CCCC", "Ph3P=O": "O=P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "Ph3PO": "O=P(c1ccccc1)(c1ccccc1)c1ccccc1", "DEAD": "CCOC(=O)N=NC(=O)OCC", "DIAD": "CC(C)OC(=O)N=NC(=O)OC(C)C",
    "DCC": "C1CCC(CC1)N=C=NC1CCCCC1", "DIC": "CC(C)N=C=NC(C)C", "EDC": "CCN=C=NCCCN(C)C", "EDCI": "CCN=C=NCCCN(C)C",
    "EDC·HCl": "CCN=C=NCCCN(C)C.Cl", "HOBt": "On1nnc2ccccc21", "HOAt": "On1nnc2cccnc21", "HOSu": "ON1C(=O)CCC1=O",
    "NHS": "ON1C(=O)CCC1=O", "CDI": "O=C(n1ccnc1)n1ccnc1", "HATU": "CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F",
    "HBTU": "CN(C)C(On1nnc2ccccc21)=[N+](C)C.F[P-](F)(F)(F)(F)F", "PyBOP": "C1CCN(CC1)[P+](N1CCCC1)(N1CCCC1)On1nnc2ccccc21.F[P-](F)(F)(F)(F)F",
    "T3P": "CCCP1(=O)OP(=O)(CCC)OP(=O)(CCC)O1", "mCPBA": "O=C(OO)c1cccc(Cl)c1", "m-CPBA": "O=C(OO)c1cccc(Cl)c1", "MCPBA": "O=C(OO)c1cccc(Cl)c1",
    "TBHP": "CC(C)(C)OO", "t-BuOOH": "CC(C)(C)OO", "tBuOOH": "CC(C)(C)OO", "DMDO": "CC1(C)OO1", "Oxone": "[K+].[K+].[K+].[K+].[K+].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)[O-].[O-]S(=O)(=O)OO.[O-]S(=O)(=O)OO.OS(=O)(=O)[O-]",
    "PCC": "[O-][Cr](=O)(=O)Cl.c1cc[nH+]cc1", "PDC": "[O-][Cr](=O)(=O)O[Cr](=O)(=O)[O-].c1cc[nH+]cc1.c1cc[nH+]cc1",
    "DMP": "CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21", "Dess-Martin": "CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21",
    "IBX": "O=C1OI(=O)(O)c2ccccc21", "TEMPO": "CC1(C)CCCC(C)(C)N1[O]", "TPAP": "CCCC[N+](CCCC)(CCCC)CCCC.[O-][Ru](=O)(=O)=O",
    "NMO": "C[N+]1([O-])CCOCC1", "Jones reagent": "O=[Cr](=O)=O.OS(=O)(=O)O", "KMnO4/H2O": "[K+].[O-][Mn](=O)(=O)=O.O",
    "LiAlH4/THF": "[Li+].[AlH4-].C1CCOC1", "Ph3P=CH2": "C=P(c1ccccc1)(c1ccccc1)c1ccccc1", "Ph3PCH3Br": "C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[Br-]",
    "MePPh3Br": "C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[Br-]",
    "Grubbs II": "CC1=CC(C)=C(C(C)=C1)N1CCN(C1=[Ru](Cl)(Cl)(=Cc1ccccc1)P(C1CCCCC1)(C1CCCCC1)C1CCCCC1)c1c(C)cc(C)cc1C",
    "NaH/DMF": "[Na+].[H-].CN(C)C=O", "KOH/EtOH": "[K+].[OH-].CCO", "HCl/MeOH": "Cl.CO", "HCl/dioxane": "Cl.C1COCCO1",
    "TFA/DCM": "OC(=O)C(F)(F)F.ClCCl", "H2SO4 (cat.)": "OS(=O)(=O)O", "H+": "[H+]", "H3O+": "[OH3+]", "OH-": "[OH-]",
    "HO-": "[OH-]", "CN-": "[C-]#N", "N3-": "[N-]=[N+]=[N-]", "MeO-": "C[O-]", "EtO-": "CC[O-]", "tBuO-": "CC(C)(C)[O-]",
    "AcO-": "CC(=O)[O-]", "OAc-": "CC(=O)[O-]", "HS-": "[SH-]", "hν": "", "hv": "", "Δ": "", "heat": "",
    "CS2": "S=C=S", "COCl2": "O=C(Cl)Cl", "phosgene": "O=C(Cl)Cl", "triphosgene": "O=C(OC(Cl)(Cl)Cl)OC(Cl)(Cl)Cl",
    "HC(OEt)3": "CCOC(OCC)OCC", "HC(OMe)3": "COC(OC)OC", "DMF-DMA": "CN(C)C(OC)OC", "CH(OEt)3": "CCOC(OCC)OCC",
    "Me3SiCN": "C[Si](C)(C)C#N", "KCN/18-crown-6": "[K+].[C-]#N.C1COCCOCCOCCOCCOCCO1", "18-crown-6": "C1COCCOCCOCCOCCOCCO1",
    "glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "urea": "NC(N)=O", "thiourea": "NC(N)=S", "guanidine": "NC(N)=N",
    "ethylene": "C=C", "acetylene": "C#C", "HCN": "C#N", "MeCHO": "CC=O", "CH3CHO": "CC=O", "PhCHO": "O=Cc1ccccc1",
    "EtOCOCl": "CCOC(=O)Cl", "MeOH/HCl": "CO.Cl", "NaOMe/MeOH": "[Na+].[O-]C.CO", "EtONa": "[Na+].[O-]CC",
    "NaOEt/EtOH": "[Na+].[O-]CC.CCO", "KOtBu/THF": "[K+].[O-]C(C)(C)C.C1CCOC1",
}

# Labels that mean "nothing chemical" (conditions) map to an empty SMILES above; filter them out on lookup.

_SUBSCRIPT = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_SUPERSCRIPT = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")


def normalize_label(text: str) -> str:
    """Collapse whitespace and unicode sub/superscripts, minus signs and dots so table lookups are forgiving."""
    s = (text or "").strip()
    s = s.translate(_SUBSCRIPT).translate(_SUPERSCRIPT)
    s = s.replace("−", "-").replace("–", "-").replace("—", "-").replace("·", "·").replace("•", "·")
    s = re.sub(r"\s+", " ", s)
    return s


def lookup(label: str) -> Optional[Tuple[str, bool]]:
    """Return (smiles, is_substituent) for a known abbreviation, else None.

    Substituent SMILES have their attachment atom first.  Lookup tries the exact label, then the label with
    whitespace/dots removed, then a few conventional variants ("n-Bu" -> "nBu", trailing dashes stripped).
    """
    if not label:
        return None
    cands = []
    s = normalize_label(label)
    cands.append(s)
    cands.append(s.replace(" ", ""))
    cands.append(s.replace("-", "").replace(" ", ""))
    cands.append(s.rstrip("-").strip())
    for c in cands:
        if c in SUBSTITUENTS:
            return SUBSTITUENTS[c], True
        if c in MOLECULES:
            smi = MOLECULES[c]
            if smi:
                return smi, False
            return None
    for c in cands:   # case-insensitive fall back for reagents/solvents only (substituent labels are case sensitive)
        lc = c.lower()
        for k, v in MOLECULES.items():
            if k.lower() == lc and v:
                return v, False
    return None


_ELEMENTS = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
             "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
             "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
             "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
             "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
             "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
             "Lv", "Ts", "Og"}

_LABEL_RE = re.compile(
    r"^(?P<iso>\d{1,3})?(?P<sym>[A-Z][a-z]?)(?:H(?P<hn>\d*))?(?P<rad>[·•]*)"
    r"(?:(?P<cn1>\d*)(?P<cs1>[+\-])|(?P<cs2>[+\-])(?P<cn2>\d*))?(?P<rad2>[·•]*)$")


def parse_atom_label(text: str) -> Optional[dict]:
    """Interpret a single-atom label.

    ``"NH3+"`` -> {symbol: "N", hcount: 3, charge: 1}; ``"O-"`` -> {symbol: "O", charge: -1};
    ``"13C"`` -> {symbol: "C", isotope: 13}; ``"D"`` -> {symbol: "H", isotope: 2}; ``"Fe3+"`` -> charge 3.
    Keys are present only when the label states them.  Returns None for anything that is not a plain atom label.
    """
    s = normalize_label(text).replace(" ", "")
    if not s:
        return None
    m = _LABEL_RE.match(s)
    if not m:
        return None
    sym = m.group("sym")
    out: dict = {}
    if sym in ("D", "T") and not m.group("iso"):
        out["symbol"] = "H"
        out["isotope"] = 2 if sym == "D" else 3
    elif sym in _ELEMENTS:
        out["symbol"] = sym
        if m.group("iso"):
            out["isotope"] = int(m.group("iso"))
    else:
        return None
    if m.group("hn") is not None:
        out["hcount"] = int(m.group("hn") or "1")
    sign = m.group("cs1") or m.group("cs2")
    if sign:
        n = m.group("cn1") or m.group("cn2") or "1"
        out["charge"] = int(n) * (1 if sign == "+" else -1)
    if m.group("rad") or m.group("rad2"):
        out["radical"] = len(m.group("rad") or "") + len(m.group("rad2") or "")
    return out


def all_labels():
    """Every known label (for tests)."""
    return list(SUBSTITUENTS) + list(MOLECULES)
