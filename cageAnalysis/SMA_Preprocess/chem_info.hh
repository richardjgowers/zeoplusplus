/*
 * chem_info.hh
 *
 *  Created on: Apr 18, 2016
 *      Author: ismael.gomez
 */

#ifndef CHEM_INFO_HH_
#define CHEM_INFO_HH_

#include <string.h>

#include <iostream>

#define	H_RAD		1.09
#define H_SYMB		"H"

#define D_RAD		1.09
#define D_SYMB		"D"

#define HE_RAD		1.4
#define HE_SYMB		"He"

#define LI_RAD		1.82
#define LI_SYMB		"Li"

#define BE_RAD		2.
#define BE_SYMB		"Be"

#define B_RAD		2.
#define B_SYMB		"B"

#define C_RAD		1.7
#define C_SYMB		"C"

#define N_RAD		1.55
#define N_SYMB		"N"

#define O_RAD		1.52
#define O_SYMB		"O"

#define F_RAD		1.47
#define F_SYMB		"F"

#define NE_RAD		1.54
#define NE_SYMB		"Ne"

#define NA_RAD		2.27
#define NA_SYMB		"Na"

#define MG_RAD		1.73
#define MG_SYMB		"Mg"

#define AL_RAD		2.
#define AL_SYMB		"Al"

#define SI_RAD		2.1
#define SI_SYMB		"Si"

#define P_RAD		1.8
#define P_SYMB		"P"

#define S_RAD		1.8
#define S_SYMB		"S"

#define CL_RAD		1.75
#define CL_SYMB		"Cl"

#define AR_RAD		1.88
#define AR_SYMB		"Ar"

#define K_RAD		2.75
#define K_SYMB		"K"

#define CA_RAD		2.
#define CA_SYMB		"Ca"

#define SC_RAD		2.
#define SC_SYMB		"Sc"

#define TI_RAD		2.
#define TI_SYMB		"Ti"

#define V_RAD		2.
#define V_SYMB		"V"

#define CR_RAD		2.
#define CR_SYMB		"Cr"

#define MN_RAD		2.
#define MN_SYMB		"Mn"

#define FE_RAD		2.
#define FE_SYMB		"Fe"

#define CO_RAD		2.
#define CO_SYMB		"Co"

#define NI_RAD		1.63
#define NI_SYMB		"Ni"

#define CU_RAD		1.4
#define CU_SYMB		"Cu"

#define ZN_RAD		1.39
#define ZN_SYMB		"Zn"

#define GA_RAD		1.87
#define GA_SYMB		"Ga"

#define GE_RAD		2.
#define GE_SYMB		"Ge"

#define AS_RAD		1.85
#define AS_SYMB		"As"

#define SE_RAD		1.9
#define SE_SYMB		"Se"

#define BR_RAD		1.85
#define BR_SYMB		"Br"

#define KR_RAD		2.02
#define KR_SYMB		"Kr"

#define RB_RAD		2.
#define RB_SYMB		"Rb"

#define SR_RAD		2.
#define SR_SYMB		"Sr"

#define Y_RAD		2.
#define Y_SYMB		"Y"

#define ZR_RAD		2.
#define ZR_SYMB		"Zr"

#define NB_RAD		2.
#define NB_SYMB		"Nb"

#define MO_RAD		2.
#define MO_SYMB		"Mo"

#define TC_RAD		2.
#define TC_SYMB		"Tc"

#define RU_RAD		2.
#define RU_SYMB		"Ru"

#define RH_RAD		2.
#define RH_SYMB		"Rh"

#define PD_RAD		1.63
#define PD_SYMB		"Pd"

#define AG_RAD		1.72
#define AG_SYMB		"Ag"

#define CD_RAD		1.58
#define CD_SYMB		"Cd"

#define IN_RAD		1.93
#define IN_SYMB		"In"

#define SN_RAD		2.17
#define SN_SYMB		"Sn"

#define SB_RAD		2.
#define SB_SYMB		"Sb"

#define TE_RAD		2.06
#define TE_SYMB		"Te"

#define I_RAD		1.98
#define I_SYMB		"I"

#define XE_RAD		2.16
#define XE_SYMB		"Xe"

#define CS_RAD		2.
#define CS_SYMB		"Cs"

#define BA_RAD		2.
#define BA_SYMB		"Ba"

#define LA_RAD		2.
#define LA_SYMB		"La"

#define CE_RAD		2.
#define CE_SYMB		"Ce"

#define PR_RAD		2.
#define PR_SYMB		"Pr"

#define ND_RAD		2.
#define ND_SYMB		"Nd"

#define PM_RAD		2.
#define PM_SYMB		"Pm"

#define SM_RAD		2.
#define SM_SYMB		"Sm"

#define EU_RAD		2.
#define EU_SYMB		"Eu"

#define GD_RAD		2.
#define GD_SYMB		"Gd"

#define TB_RAD		2.
#define TB_SYMB		"Tb"

#define DY_RAD		2.
#define DY_SYMB		"Dy"

#define HO_RAD		2.
#define HO_SYMB		"Ho"

#define ER_RAD		2.
#define ER_SYMB		"Er"

#define TM_RAD		2.
#define TM_SYMB		"Tm"

#define YB_RAD		2.
#define YB_SYMB		"Yb"

#define LU_RAD		2.
#define LU_SYMB		"Lu"

#define HF_RAD		2.
#define HF_SYMB		"Hf"

#define TA_RAD		2.
#define TA_SYMB		"Ta"

#define W_RAD		2.
#define W_SYMB		"W"

#define RE_RAD		2.
#define RE_SYMB		"Re"

#define OS_RAD		2.
#define OS_SYMB		"Os"

#define IR_RAD		2.
#define IR_SYMB		"Ir"

#define PT_RAD		1.72
#define PT_SYMB		"Pt"

#define AU_RAD		1.66
#define AU_SYMB		"Au"

#define HG_RAD		1.55
#define HG_SYMB		"Hg"

#define TL_RAD		1.96
#define TL_SYMB		"Tl"

#define PB_RAD		2.02
#define PB_SYMB		"Pb"

#define BI_RAD		2.
#define BI_SYMB		"Bi"

#define PO_RAD		2.
#define PO_SYMB		"Po"

#define AT_RAD		2.
#define AT_SYMB		"At"

#define RN_RAD		2.
#define RN_SYMB		"Rn"

#define FR_RAD		2.
#define FR_SYMB		"Fr"

#define RA_RAD		2.
#define RA_SYMB		"Ra"

#define AC_RAD		2.
#define AC_SYMB		"Ac"

#define TH_RAD		2.
#define TH_SYMB		"Th"

#define PA_RAD		2.
#define PA_SYMB		"Pa"

#define U_RAD		1.86
#define U_SYMB		"U"

#define NP_RAD		2.
#define NP_SYMB		"Np"

#define PU_RAD		2.
#define PU_SYMB		"Pu"

#define AM_RAD		2.
#define AM_SYMB		"Am"

#define CM_RAD		2.
#define CM_SYMB		"Cm"

#define BK_RAD		2.
#define BK_SYMB		"Bk"

#define CF_RAD		2.
#define CF_SYMB		"Cf"

#define ES_RAD		2.
#define ES_SYMB		"Es"

#define FM_RAD		2.
#define FM_SYMB		"Fm"

#define MD_RAD		2.
#define MD_SYMB		"Md"

#define NO_RAD		2.
#define NO_SYMB		"No"

#define LR_RAD		2.
#define LR_SYMB		"Lr"

#define RF_RAD		2.
#define RF_SYMB		"Rf"

#define DB_RAD		2.
#define DB_SYMB		"Db"

#define SG_RAD		2.
#define SG_SYMB		"Sg"

#define BH_RAD		2.
#define BH_SYMB		"Bh"

#define HS_RAD		2.
#define HS_SYMB		"Hs"

#define MT_RAD		2.
#define MT_SYMB		"Mt"

#define DS_RAD		2.
#define DS_SYMB		"DS"

/*static std::map <std::string,double> radTable;

radTable.insert(pair <string,double> ("H",  1.09));
radTable.insert(pair <string,double> ("D",  1.09));
radTable.insert(pair <string,double> ("He",  1.4));
radTable.insert(pair <string,double> ("Li",  1.82));
radTable.insert(pair <string,double> ("Be",  2));
radTable.insert(pair <string,double> ("B",  2));
radTable.insert(pair <string,double> ("C",  1.7));
radTable.insert(pair <string,double> ("N",  1.55));
radTable.insert(pair <string,double> ("O",  1.52));
radTable.insert(pair <string,double> ("F",  1.47));
radTable.insert(pair <string,double> ("Ne",  1.54));
radTable.insert(pair <string,double> ("Na",  2.27));
radTable.insert(pair <string,double> ("Mg",  1.73));
radTable.insert(pair <string,double> ("Al",  2));
radTable.insert(pair <string,double> ("Si",  2.1));
radTable.insert(pair <string,double> ("P",  1.8));
radTable.insert(pair <string,double> ("S",  1.8));
radTable.insert(pair <string,double> ("Cl",  1.75));
radTable.insert(pair <string,double> ("Ar",  1.88));
radTable.insert(pair <string,double> ("K",  2.75));
radTable.insert(pair <string,double> ("Ca",  2));
radTable.insert(pair <string,double> ("Sc",  2));
radTable.insert(pair <string,double> ("Ti",  2));
radTable.insert(pair <string,double> ("V",  2));
radTable.insert(pair <string,double> ("Cr",  2));
radTable.insert(pair <string,double> ("Mn",  2));
radTable.insert(pair <string,double> ("Fe",  2));
radTable.insert(pair <string,double> ("Co",  2));
radTable.insert(pair <string,double> ("Ni",  1.63));
radTable.insert(pair <string,double> ("Cu",  1.4));
radTable.insert(pair <string,double> ("Zn",  1.39));
radTable.insert(pair <string,double> ("Ga",  1.87));
radTable.insert(pair <string,double> ("Ge",  2));
radTable.insert(pair <string,double> ("As",  1.85));
radTable.insert(pair <string,double> ("Se",  1.9));
radTable.insert(pair <string,double> ("Br",  1.85));
radTable.insert(pair <string,double> ("Kr",  2.02));
radTable.insert(pair <string,double> ("Rb",  2));
radTable.insert(pair <string,double> ("Sr",  2));
radTable.insert(pair <string,double> ("Y",  2));
radTable.insert(pair <string,double> ("Zr",  2));
radTable.insert(pair <string,double> ("Nb",  2));
radTable.insert(pair <string,double> ("Mo",  2));
radTable.insert(pair <string,double> ("Tc",  2));
radTable.insert(pair <string,double> ("Ru",  2));
radTable.insert(pair <string,double> ("Rh",  2));
radTable.insert(pair <string,double> ("Pd",  1.63));
radTable.insert(pair <string,double> ("Ag",  1.72));
radTable.insert(pair <string,double> ("Cd",  1.58));
radTable.insert(pair <string,double> ("In",  1.93));
radTable.insert(pair <string,double> ("Sn",  2.17));
radTable.insert(pair <string,double> ("Sb",  2));
radTable.insert(pair <string,double> ("Te",  2.06));
radTable.insert(pair <string,double> ("I",  1.98));
radTable.insert(pair <string,double> ("Xe",  2.16));
radTable.insert(pair <string,double> ("Cs",  2));
radTable.insert(pair <string,double> ("Ba",  2));
radTable.insert(pair <string,double> ("La",  2));
radTable.insert(pair <string,double> ("Ce",  2));
radTable.insert(pair <string,double> ("Pr",  2));
radTable.insert(pair <string,double> ("Nd",  2));
radTable.insert(pair <string,double> ("Pm",  2));
radTable.insert(pair <string,double> ("Sm",  2));
radTable.insert(pair <string,double> ("Eu",  2));
radTable.insert(pair <string,double> ("Gd",  2));
radTable.insert(pair <string,double> ("Tb",  2));
radTable.insert(pair <string,double> ("Dy",  2));
radTable.insert(pair <string,double> ("Ho",  2));
radTable.insert(pair <string,double> ("Er",  2));
radTable.insert(pair <string,double> ("Tm",  2));
radTable.insert(pair <string,double> ("Yb",  2));
radTable.insert(pair <string,double> ("Lu",  2));
radTable.insert(pair <string,double> ("Hf",  2));
radTable.insert(pair <string,double> ("Ta",  2));
radTable.insert(pair <string,double> ("W",  2));
radTable.insert(pair <string,double> ("Re",  2));
radTable.insert(pair <string,double> ("Os",  2));
radTable.insert(pair <string,double> ("Ir",  2));
radTable.insert(pair <string,double> ("Pt",  1.72));
radTable.insert(pair <string,double> ("Au",  1.66));
radTable.insert(pair <string,double> ("Hg",  1.55));
radTable.insert(pair <string,double> ("Tl",  1.96));
radTable.insert(pair <string,double> ("Pb",  2.02));
radTable.insert(pair <string,double> ("Bi",  2));
radTable.insert(pair <string,double> ("Po",  2));
radTable.insert(pair <string,double> ("At",  2));
radTable.insert(pair <string,double> ("Rn",  2));
radTable.insert(pair <string,double> ("Fr",  2));
radTable.insert(pair <string,double> ("Ra",  2));
radTable.insert(pair <string,double> ("Ac",  2));
radTable.insert(pair <string,double> ("Th",  2));
radTable.insert(pair <string,double> ("Pa",  2));
radTable.insert(pair <string,double> ("U",  1.86));
radTable.insert(pair <string,double> ("Np",  2));
radTable.insert(pair <string,double> ("Pu",  2));
radTable.insert(pair <string,double> ("Am",  2));
radTable.insert(pair <string,double> ("Cm",  2));
radTable.insert(pair <string,double> ("Bk",  2));
radTable.insert(pair <string,double> ("Cf",  2));
radTable.insert(pair <string,double> ("Es",  2));
radTable.insert(pair <string,double> ("Fm",  2));
radTable.insert(pair <string,double> ("Md",  2));
radTable.insert(pair <string,double> ("No",  2));
radTable.insert(pair <string,double> ("Lr",  2));
radTable.insert(pair <string,double> ("Rf",  2));
radTable.insert(pair <string,double> ("Db",  2));
radTable.insert(pair <string,double> ("Sg",  2));
radTable.insert(pair <string,double> ("Bh",  2));
radTable.insert(pair <string,double> ("Hs",  2));
radTable.insert(pair <string,double> ("Mt",  2));
radTable.insert(pair <string,double> ("Ds",  2));*/


using namespace std;

double AtomRadius(char *Element);

#endif /* CHEM_INFO_HH_ */
