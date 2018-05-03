/*
 * chem_info.cc
 *
 *  Created on: Apr 18, 2016
 *      Author: ismael.gomez
 */

#include "chem_info.hh"

//////////// TAKEN FROM ZEO++ (TODO: LINK instead)



double AtomRadius(char *Element)
{


	if (strncmp(Element, H_SYMB, 2) == 0) { return H_RAD;}

	if (strncmp(Element, D_SYMB, 2) == 0) { return D_RAD; }

	if (strncmp(Element, HE_SYMB, 2) == 0) { return HE_RAD; }

	if (strncmp(Element, LI_SYMB, 2) == 0) { return LI_RAD; }

	if (strncmp(Element, BE_SYMB, 2) == 0) { return BE_RAD; }

	if (strncmp(Element, B_SYMB, 2) == 0) { return B_RAD; }

	if (strncmp(Element, C_SYMB, 2) == 0) { return C_RAD; }

	if (strncmp(Element, N_SYMB, 2) == 0) { return N_RAD; }

	if (strncmp(Element, O_SYMB, 2) == 0) { return O_RAD; }

	if (strncmp(Element, F_SYMB, 2) == 0) { return F_RAD; }

	if (strncmp(Element, NE_SYMB, 2) == 0) { return NE_RAD; }

	if (strncmp(Element, NA_SYMB, 2) == 0) { return NA_RAD; }

	if (strncmp(Element, MG_SYMB, 2) == 0) { return MG_RAD; }

	if (strncmp(Element, AL_SYMB, 2) == 0) { return AL_RAD; }

	if (strncmp(Element, SI_SYMB, 2) == 0) { return SI_RAD; }

	if (strncmp(Element, P_SYMB, 2) == 0) { return P_RAD; }

	if (strncmp(Element, S_SYMB, 2) == 0) { return S_RAD; }

	if (strncmp(Element, CL_SYMB, 2) == 0) { return CL_RAD; }

	if (strncmp(Element, AR_SYMB, 2) == 0) { return AR_RAD; }

	if (strncmp(Element, K_SYMB, 2) == 0) { return K_RAD; }

	if (strncmp(Element, CA_SYMB, 2) == 0) { return CA_RAD; }

	if (strncmp(Element, SC_SYMB, 2) == 0) { return SC_RAD; }

	if (strncmp(Element, TI_SYMB, 2) == 0) { return TI_RAD; }

	if (strncmp(Element, V_SYMB, 2) == 0) { return V_RAD; }

	if (strncmp(Element, CR_SYMB, 2) == 0) { return CR_RAD; }

	if (strncmp(Element, MN_SYMB, 2) == 0) { return MN_RAD; }

	if (strncmp(Element, FE_SYMB, 2) == 0) { return FE_RAD; }

	if (strncmp(Element, CO_SYMB, 2) == 0) { return CO_RAD; }

	if (strncmp(Element, NI_SYMB, 2) == 0) { return NI_RAD; }

	if (strncmp(Element, CU_SYMB, 2) == 0) { return CU_RAD; }

	if (strncmp(Element, ZN_SYMB, 2) == 0) { return ZN_RAD; }

	if (strncmp(Element, GA_SYMB, 2) == 0) { return GA_RAD; }

	if (strncmp(Element, GE_SYMB, 2) == 0) { return GE_RAD; }

	if (strncmp(Element, AS_SYMB, 2) == 0) { return AS_RAD; }

	if (strncmp(Element, SE_SYMB, 2) == 0) { return SE_RAD; }

	if (strncmp(Element, BR_SYMB, 2) == 0) { return BR_RAD; }

	if (strncmp(Element, KR_SYMB, 2) == 0) { return KR_RAD; }

	if (strncmp(Element, PD_SYMB, 2) == 0) { return PD_RAD; }

	if (strncmp(Element, AG_SYMB, 2) == 0) { return AG_RAD; }

	if (strncmp(Element, CD_SYMB, 2) == 0) { return CD_RAD; }

	if (strncmp(Element, IN_SYMB, 2) == 0) { return IN_RAD; }

	if (strncmp(Element, SN_SYMB, 2) == 0) { return SN_RAD; }

	if (strncmp(Element, SB_SYMB, 2) == 0) { return SB_RAD; }

	if (strncmp(Element, TE_SYMB, 2) == 0) { return TE_RAD; }

	if (strncmp(Element, I_SYMB, 2) == 0) { return I_RAD; }

	if (strncmp(Element, XE_SYMB, 2) == 0) { return XE_RAD; }

	if (strncmp(Element, PT_SYMB, 2) == 0) { return PT_RAD; }

	if (strncmp(Element, AU_SYMB, 2) == 0) { return AU_RAD; }

	if (strncmp(Element, HG_SYMB, 2) == 0) { return HG_RAD; }

	if (strncmp(Element, TL_SYMB, 2) == 0) { return TL_RAD; }

	if (strncmp(Element, PB_SYMB, 2) == 0) { return PB_RAD; }

	if (strncmp(Element, BI_SYMB, 2) == 0) { return BI_RAD; }

	if (strncmp(Element, PO_SYMB, 2) == 0) { return PO_RAD; }

	if (strncmp(Element, AT_SYMB, 2) == 0) { return AT_RAD; }

	if (strncmp(Element, RN_SYMB, 2) == 0) { return RN_RAD; }

	if (strncmp(Element, FR_SYMB, 2) == 0) { return FR_RAD; }

	if (strncmp(Element, RA_SYMB, 2) == 0) { return RA_RAD; }

	if (strncmp(Element, AC_SYMB, 2) == 0) { return AC_RAD; }

	if (strncmp(Element, TH_SYMB, 2) == 0) { return TH_RAD; }

	if (strncmp(Element, PA_SYMB, 2) == 0) { return PA_RAD; }

	if (strncmp(Element, U_SYMB, 2) == 0) { return U_RAD; }

	if (strncmp(Element, NP_SYMB, 2) == 0) { return NP_RAD; }

	if (strncmp(Element, PU_SYMB, 2) == 0) { return PU_RAD; }

	if (strncmp(Element, AM_SYMB, 2) == 0) { return AM_RAD; }

	if (strncmp(Element, CM_SYMB, 2) == 0) { return CM_RAD; }

	if (strncmp(Element, BK_SYMB, 2) == 0) { return BK_RAD; }

	if (strncmp(Element, CF_SYMB, 2) == 0) { return CF_RAD; }

	if (strncmp(Element, ES_SYMB, 2) == 0) { return ES_RAD; }

	if (strncmp(Element, FM_SYMB, 2) == 0) { return FM_RAD; }

	if (strncmp(Element, MD_SYMB, 2) == 0) { return MD_RAD; }

	if (strncmp(Element, NO_SYMB, 2) == 0) { return NO_RAD; }

	if (strncmp(Element, LR_SYMB, 2) == 0) { return LR_RAD; }

	if (strncmp(Element, RF_SYMB, 2) == 0) { return RF_RAD; }

	if (strncmp(Element, DB_SYMB, 2) == 0) { return DB_RAD; }

	if (strncmp(Element, SG_SYMB, 2) == 0) { return SG_RAD; }

	if (strncmp(Element, BH_SYMB, 2) == 0) { return BH_RAD; }

	if (strncmp(Element, HS_SYMB, 2) == 0) { return HS_RAD; }

	if (strncmp(Element, MT_SYMB, 2) == 0) { return MT_RAD; }

	if (strncmp(Element, DS_SYMB, 2) == 0) { return DS_RAD; }

	if (strncmp(Element, RB_SYMB, 2) == 0) { return RB_RAD; }

	if (strncmp(Element, SR_SYMB, 2) == 0) { return SR_RAD; }

	if (strncmp(Element, Y_SYMB, 2) == 0) { return Y_RAD; }

	if (strncmp(Element, ZR_SYMB, 2) == 0) { return ZR_RAD; }

	if (strncmp(Element, NB_SYMB, 2) == 0) { return NB_RAD; }

	if (strncmp(Element, MO_SYMB, 2) == 0) { return MO_RAD; }

	if (strncmp(Element, TC_SYMB, 2) == 0) { return TC_RAD; }

	if (strncmp(Element, RU_SYMB, 2) == 0) { return RU_RAD; }

	if (strncmp(Element, RH_SYMB, 2) == 0) { return RH_RAD; }

	if (strncmp(Element, CS_SYMB, 2) == 0) { return CS_RAD; }

	if (strncmp(Element, BA_SYMB, 2) == 0) { return BA_RAD; }

	if (strncmp(Element, LA_SYMB, 2) == 0) { return LA_RAD; }

	if (strncmp(Element, CE_SYMB, 2) == 0) { return CE_RAD; }

	if (strncmp(Element, PR_SYMB, 2) == 0) { return PR_RAD; }

	if (strncmp(Element, ND_SYMB, 2) == 0) { return ND_RAD; }

	if (strncmp(Element, PM_SYMB, 2) == 0) { return PM_RAD; }

	if (strncmp(Element, SM_SYMB, 2) == 0) { return SM_RAD; }

	if (strncmp(Element, EU_SYMB, 2) == 0) { return EU_RAD; }

	if (strncmp(Element, GD_SYMB, 2) == 0) { return GD_RAD; }

	if (strncmp(Element, TB_SYMB, 2) == 0) { return TB_RAD; }

	if (strncmp(Element, DY_SYMB, 2) == 0) { return DY_RAD; }

	if (strncmp(Element, HO_SYMB, 2) == 0) { return HO_RAD; }

	if (strncmp(Element, ER_SYMB, 2) == 0) { return ER_RAD; }

	if (strncmp(Element, TM_SYMB, 2) == 0) { return TM_RAD; }

	if (strncmp(Element, YB_SYMB, 2) == 0) { return YB_RAD; }

	if (strncmp(Element, LU_SYMB, 2) == 0) { return LU_RAD; }

	if (strncmp(Element, HF_SYMB, 2) == 0) { return HF_RAD; }

	if (strncmp(Element, TA_SYMB, 2) == 0) { return TA_RAD; }

	if (strncmp(Element, W_SYMB, 2) == 0) { return W_RAD; }

	if (strncmp(Element, RE_SYMB, 2) == 0) { return RE_RAD; }

	if (strncmp(Element, OS_SYMB, 2) == 0) { return OS_RAD; }

	if (strncmp(Element, IR_SYMB, 2) == 0) { return IR_RAD; }

	/*if (strncmp(Element, H_SYMB, 2) == 0) { return H_RAD/2.;}

	if (strncmp(Element, D_SYMB, 2) == 0) { return D_RAD/2.; }

	if (strncmp(Element, HE_SYMB, 2) == 0) { return HE_RAD/2.; }

	if (strncmp(Element, LI_SYMB, 2) == 0) { return LI_RAD/2.; }

	if (strncmp(Element, BE_SYMB, 2) == 0) { return BE_RAD/2.; }

	if (strncmp(Element, B_SYMB, 2) == 0) { return B_RAD/2.; }

	if (strncmp(Element, C_SYMB, 2) == 0) { return C_RAD; }

	if (strncmp(Element, N_SYMB, 2) == 0) { return N_RAD/2.; }

	if (strncmp(Element, O_SYMB, 2) == 0) { return O_RAD/2.; }

	if (strncmp(Element, F_SYMB, 2) == 0) { return F_RAD/2.; }

	if (strncmp(Element, NE_SYMB, 2) == 0) { return NE_RAD/2.; }

	if (strncmp(Element, NA_SYMB, 2) == 0) { return NA_RAD/2.; }

	if (strncmp(Element, MG_SYMB, 2) == 0) { return MG_RAD/2.; }

	if (strncmp(Element, AL_SYMB, 2) == 0) { return AL_RAD/2.; }

	if (strncmp(Element, SI_SYMB, 2) == 0) { return SI_RAD/2.; }

	if (strncmp(Element, P_SYMB, 2) == 0) { return P_RAD/2.; }

	if (strncmp(Element, S_SYMB, 2) == 0) { return S_RAD/2.; }

	if (strncmp(Element, CL_SYMB, 2) == 0) { return CL_RAD/2.; }

	if (strncmp(Element, AR_SYMB, 2) == 0) { return AR_RAD/2.; }

	if (strncmp(Element, K_SYMB, 2) == 0) { return K_RAD/2.; }

	if (strncmp(Element, CA_SYMB, 2) == 0) { return CA_RAD/2.; }

	if (strncmp(Element, SC_SYMB, 2) == 0) { return SC_RAD/2.; }

	if (strncmp(Element, TI_SYMB, 2) == 0) { return TI_RAD/2.; }

	if (strncmp(Element, V_SYMB, 2) == 0) { return V_RAD/2.; }

	if (strncmp(Element, CR_SYMB, 2) == 0) { return CR_RAD/2.; }

	if (strncmp(Element, MN_SYMB, 2) == 0) { return MN_RAD/2.; }

	if (strncmp(Element, FE_SYMB, 2) == 0) { return FE_RAD/2.; }

	if (strncmp(Element, CO_SYMB, 2) == 0) { return CO_RAD/2.; }

	if (strncmp(Element, NI_SYMB, 2) == 0) { return NI_RAD/2.; }

	if (strncmp(Element, CU_SYMB, 2) == 0) { return CU_RAD/2.; }

	if (strncmp(Element, ZN_SYMB, 2) == 0) { return ZN_RAD/2.; }

	if (strncmp(Element, GA_SYMB, 2) == 0) { return GA_RAD/2.; }

	if (strncmp(Element, GE_SYMB, 2) == 0) { return GE_RAD/2.; }

	if (strncmp(Element, AS_SYMB, 2) == 0) { return AS_RAD/2.; }

	if (strncmp(Element, SE_SYMB, 2) == 0) { return SE_RAD/2.; }

	if (strncmp(Element, BR_SYMB, 2) == 0) { return BR_RAD/2.; }

	if (strncmp(Element, KR_SYMB, 2) == 0) { return KR_RAD/2.; }

	if (strncmp(Element, PD_SYMB, 2) == 0) { return PD_RAD/2.; }

	if (strncmp(Element, AG_SYMB, 2) == 0) { return AG_RAD/2.; }

	if (strncmp(Element, CD_SYMB, 2) == 0) { return CD_RAD/2.; }

	if (strncmp(Element, IN_SYMB, 2) == 0) { return IN_RAD/2.; }

	if (strncmp(Element, SN_SYMB, 2) == 0) { return SN_RAD/2.; }

	if (strncmp(Element, SB_SYMB, 2) == 0) { return SB_RAD/2.; }

	if (strncmp(Element, TE_SYMB, 2) == 0) { return TE_RAD/2.; }

	if (strncmp(Element, I_SYMB, 2) == 0) { return I_RAD/2.; }

	if (strncmp(Element, XE_SYMB, 2) == 0) { return XE_RAD/2.; }

	if (strncmp(Element, PT_SYMB, 2) == 0) { return PT_RAD/2.; }

	if (strncmp(Element, AU_SYMB, 2) == 0) { return AU_RAD/2.; }

	if (strncmp(Element, HG_SYMB, 2) == 0) { return HG_RAD/2.; }

	if (strncmp(Element, TL_SYMB, 2) == 0) { return TL_RAD/2.; }

	if (strncmp(Element, PB_SYMB, 2) == 0) { return PB_RAD/2.; }

	if (strncmp(Element, BI_SYMB, 2) == 0) { return BI_RAD/2.; }

	if (strncmp(Element, PO_SYMB, 2) == 0) { return PO_RAD/2.; }

	if (strncmp(Element, AT_SYMB, 2) == 0) { return AT_RAD/2.; }

	if (strncmp(Element, RN_SYMB, 2) == 0) { return RN_RAD/2.; }

	if (strncmp(Element, FR_SYMB, 2) == 0) { return FR_RAD/2.; }

	if (strncmp(Element, RA_SYMB, 2) == 0) { return RA_RAD/2.; }

	if (strncmp(Element, AC_SYMB, 2) == 0) { return AC_RAD/2.; }

	if (strncmp(Element, TH_SYMB, 2) == 0) { return TH_RAD/2.; }

	if (strncmp(Element, PA_SYMB, 2) == 0) { return PA_RAD/2.; }

	if (strncmp(Element, U_SYMB, 2) == 0) { return U_RAD/2.; }

	if (strncmp(Element, NP_SYMB, 2) == 0) { return NP_RAD/2.; }

	if (strncmp(Element, PU_SYMB, 2) == 0) { return PU_RAD/2.; }

	if (strncmp(Element, AM_SYMB, 2) == 0) { return AM_RAD/2.; }

	if (strncmp(Element, CM_SYMB, 2) == 0) { return CM_RAD/2.; }

	if (strncmp(Element, BK_SYMB, 2) == 0) { return BK_RAD/2.; }

	if (strncmp(Element, CF_SYMB, 2) == 0) { return CF_RAD/2.; }

	if (strncmp(Element, ES_SYMB, 2) == 0) { return ES_RAD/2.; }

	if (strncmp(Element, FM_SYMB, 2) == 0) { return FM_RAD/2.; }

	if (strncmp(Element, MD_SYMB, 2) == 0) { return MD_RAD/2.; }

	if (strncmp(Element, NO_SYMB, 2) == 0) { return NO_RAD/2.; }

	if (strncmp(Element, LR_SYMB, 2) == 0) { return LR_RAD/2.; }

	if (strncmp(Element, RF_SYMB, 2) == 0) { return RF_RAD/2.; }

	if (strncmp(Element, DB_SYMB, 2) == 0) { return DB_RAD/2.; }

	if (strncmp(Element, SG_SYMB, 2) == 0) { return SG_RAD/2.; }

	if (strncmp(Element, BH_SYMB, 2) == 0) { return BH_RAD/2.; }

	if (strncmp(Element, HS_SYMB, 2) == 0) { return HS_RAD/2.; }

	if (strncmp(Element, MT_SYMB, 2) == 0) { return MT_RAD/2.; }

	if (strncmp(Element, DS_SYMB, 2) == 0) { return DS_RAD/2.; }

	if (strncmp(Element, RB_SYMB, 2) == 0) { return RB_RAD/2.; }

	if (strncmp(Element, SR_SYMB, 2) == 0) { return SR_RAD/2.; }

	if (strncmp(Element, Y_SYMB, 2) == 0) { return Y_RAD/2.; }

	if (strncmp(Element, ZR_SYMB, 2) == 0) { return ZR_RAD/2.; }

	if (strncmp(Element, NB_SYMB, 2) == 0) { return NB_RAD/2.; }

	if (strncmp(Element, MO_SYMB, 2) == 0) { return MO_RAD/2.; }

	if (strncmp(Element, TC_SYMB, 2) == 0) { return TC_RAD/2.; }

	if (strncmp(Element, RU_SYMB, 2) == 0) { return RU_RAD/2.; }

	if (strncmp(Element, RH_SYMB, 2) == 0) { return RH_RAD/2.; }

	if (strncmp(Element, CS_SYMB, 2) == 0) { return CS_RAD/2.; }

	if (strncmp(Element, BA_SYMB, 2) == 0) { return BA_RAD/2.; }

	if (strncmp(Element, LA_SYMB, 2) == 0) { return LA_RAD/2.; }

	if (strncmp(Element, CE_SYMB, 2) == 0) { return CE_RAD/2.; }

	if (strncmp(Element, PR_SYMB, 2) == 0) { return PR_RAD/2.; }

	if (strncmp(Element, ND_SYMB, 2) == 0) { return ND_RAD/2.; }

	if (strncmp(Element, PM_SYMB, 2) == 0) { return PM_RAD/2.; }

	if (strncmp(Element, SM_SYMB, 2) == 0) { return SM_RAD/2.; }

	if (strncmp(Element, EU_SYMB, 2) == 0) { return EU_RAD/2.; }

	if (strncmp(Element, GD_SYMB, 2) == 0) { return GD_RAD/2.; }

	if (strncmp(Element, TB_SYMB, 2) == 0) { return TB_RAD/2.; }

	if (strncmp(Element, DY_SYMB, 2) == 0) { return DY_RAD/2.; }

	if (strncmp(Element, HO_SYMB, 2) == 0) { return HO_RAD/2.; }

	if (strncmp(Element, ER_SYMB, 2) == 0) { return ER_RAD/2.; }

	if (strncmp(Element, TM_SYMB, 2) == 0) { return TM_RAD/2.; }

	if (strncmp(Element, YB_SYMB, 2) == 0) { return YB_RAD/2.; }

	if (strncmp(Element, LU_SYMB, 2) == 0) { return LU_RAD/2.; }

	if (strncmp(Element, HF_SYMB, 2) == 0) { return HF_RAD/2.; }

	if (strncmp(Element, TA_SYMB, 2) == 0) { return TA_RAD/2.; }

	if (strncmp(Element, W_SYMB, 2) == 0) { return W_RAD/2.; }

	if (strncmp(Element, RE_SYMB, 2) == 0) { return RE_RAD/2.; }

	if (strncmp(Element, OS_SYMB, 2) == 0) { return OS_RAD/2.; }

	if (strncmp(Element, IR_SYMB, 2) == 0) { return IR_RAD/2.; }*/

	//cout << "ERROR: Unknown atom radius: " << Element << endl;

	return 0.0;





}
