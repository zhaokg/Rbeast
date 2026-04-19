#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include "beastv2_header.h"
 

static void DD_CalcBasisKsKeK_prec0123(BEAST2_BASIS_PTR  basis) {

	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//  KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	

	I32		   period   = basis->bConst.dummy.period; // period is a field of the dummy basis only!
	TKNOT_PTR  KNOT     = basis->KNOT;
	I16PTR     KS       = basis->ks;
	I16PTR     KE       = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; // Number of seasonal segment		
	int kCounter   = 1L;                  // A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		//Loop through all the terms of each segment: order is the number of the temrs
		// order must be less than or equal to period
		I32 Nseg   = (KNOT[(i)-1] - 1) - KNOT[(i-1)-1];
		I32 Kterms = Nseg >= period ? period : Nseg;
		kCounter   += Kterms;
		// the end index of the seasonal segment		
		*KE++ = kCounter - 1L;
	}
	// Total number of terms for the season component
	basis->K = kCounter - 1;

}

static void VV_CalcBasisKsKeK_prec0123(BEAST2_BASIS_PTR  basis)
{
	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//  KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	

	//U08PTR     TERM_TYPE = basis->termType;
	TORDER_PTR ORDER = basis->ORDER;
	I16PTR     KS    = basis->ks;
	I16PTR     KE    = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; // Number of seasonal segment		
	int kCounter   = 1L;                  // A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		//Loop through all the terms of each segment: Number of terms is 2*Harmonic Order		
		I32 order = ORDER[i - 1];
		for (I32 j = 1; j <= order; j++) {
			//The first term of each order is SIN
			///////////*TERM_TYPE++ = j; // 0 stands for a seasonal term
			kCounter++;
			//The second term of each order is COS
			///////////*TERM_TYPE++ = j; // 0 stands for a seasonal term
			//kCounter++;
		}
		// the end index of the seasonal segment		
		*KE++ = kCounter - 1L;
	}
	// Total number of terms for the season component
	basis->K = kCounter - 1;

}
static void SS_CalcBasisKsKeK_prec012(BEAST2_BASIS_PTR  basis)
{
	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//  KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	

	//U08PTR     TERM_TYPE = basis->termType;
	TORDER_PTR ORDER     = basis->ORDER;
	I16PTR     KS        = basis->ks;
	I16PTR     KE        = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; // Number of seasonal segment		
	int kCounter   = 1L;                  // A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		//Loop through all the terms of each segment: Number of terms is 2*Harmonic Order		
		I32 order = ORDER[i - 1];
		for (I32 j = 1; j <= order; j++) {
			//The first term of each order is SIN
			///////////*TERM_TYPE++ = j; // 0 stands for a seasonal term
			kCounter++;
			//The second term of each order is COS
			///////////*TERM_TYPE++ = j; // 0 stands for a seasonal term
			kCounter++;
		}
		// the end index of the seasonal segment		
		*KE++ = kCounter - 1L;
	}
	// Total number of terms for the season component
	basis->K = kCounter - 1;

}
static void TT_CalcBasisKsKeK_prec012(BEAST2_BASIS_PTR  basis)
{

	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS	  
    // Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	
	//U08PTR TERM_TYPE = basis->termType;
	TORDER_PTR ORDER = basis->ORDER;
	I16PTR KS = basis->ks;
	I16PTR KE = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; //number of seasonal segment		
	int kCounter = 1L;                            //A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		I32 order = ORDER[i - 1];		
		for (I32 j = 1; j <= order+1; j++) {
			//*TERM_TYPE++ = j; //j=1 corresponds to the zero-th order, j=2 to the linear term
			kCounter++;
		}
		*KE++ = kCounter - 1L;
	}
	//total number of terms for the season component
	basis->K = kCounter - 1;

}
static void OO_CalcBasisKsKeK_prec012(BEAST2_BASIS_PTR  basis)
{
	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS	  
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	
	U08PTR TERM_TYPE = basis->termType;
	TORDER_PTR ORDER = basis->ORDER;
	I16PTR KS = basis->ks;
	I16PTR KE = basis->ke;

	int NUM_OF_SEG = basis->nKnot; //number of seasonal segment		
	int kCounter   = 1L;                            //A counter for basis terms
	for (rI32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		*KE++ = kCounter;
		//*TERM_TYPE++ = 1;
		kCounter++;
	}
	//total number of terms for the season component
	basis->K = kCounter - 1;

}


static void SS_CalcBasisKsKeK_prec3(BEAST2_BASIS_PTR  basis)
{
	  
	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//  KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	

	U08PTR     TERM_TYPE = basis->termType;
	TORDER_PTR ORDER = basis->ORDER;
	I16PTR     KS = basis->ks;
	I16PTR     KE = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; // Number of seasonal segment		
	int kCounter   = 1L;                  // A counter for basis terms
	for (rI32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		//Loop through all the terms of each segment: Number of terms is 2*Harmonic Order		
		I32 order = ORDER[i - 1];
		for (rI32 j = 1; j <= order; j++) {
			//The first term of each order is SIN
			*TERM_TYPE++ = j; // 0 stands for a seasonal term
			kCounter++;
			//The second term of each order is COS
			*TERM_TYPE++ = j; // 0 stands for a seasonal term
			kCounter++;
		}
		// the end index of the seasonal segment		
		*KE++ = kCounter - 1L;
	}
	//total number of terms for the season component
	basis->K = kCounter - 1;

}
static void TT_CalcBasisKsKeK_prec3(BEAST2_BASIS_PTR  basis)
{

	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS	  
    // Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	
	U08PTR     TERM_TYPE = basis->termType;
	TORDER_PTR ORDER     = basis->ORDER;
	I16PTR KS = basis->ks;
	I16PTR KE = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L; //number of seasonal segment		
	int kCounter = 1L;                            //A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		I32 order = ORDER[i - 1];		
		for (I32 j = 1; j <= order+1; j++) {
			*TERM_TYPE++ = j; //j=1 corresponds to the zero-th order, j=2 to the linear term
			kCounter++;
		}
		*KE++ = kCounter - 1L;
	}
	//total number of terms for the season component
	basis->K = kCounter - 1;

}
static void OO_CalcBasisKsKeK_prec3(BEAST2_BASIS_PTR  basis)
{
	/***************************************************************
	// Determine the start and end indices for each seasonal segment
	***************************************************************/
	//KS, KE, TERMTYPE, & K_CONST_TERMS are temporary varialbes to access to the pointers of TREND_BASIS	  
	// Use of them to reduce the memory access costs assocaited with basis->xxx: basis->sY.ks[i] is more expensive than pt++	
	U08PTR TERM_TYPE = basis->termType;
	TORDER_PTR ORDER = basis->ORDER;
	I16PTR KS = basis->ks;
	I16PTR KE = basis->ke;

	int NUM_OF_SEG = basis->nKnot; //number of seasonal segment		
	int kCounter   = 1L;                            //A counter for basis terms
	for (rI32 i = 1; i <= NUM_OF_SEG; i++) {
		//The start index is k
		*KS++ = kCounter;
		*KE++ = kCounter;
		*TERM_TYPE++ = 1;
		kCounter++;
	}
	//total number of terms for the season component
	basis->K = kCounter - 1;

}

void* Get_CalcBasisKsKeK(I08 id, I08 precPriorType) {
	if (precPriorType == 0 || precPriorType == 1 || precPriorType == 2)
	{
		switch (id) {
		case DUMMYID:   return DD_CalcBasisKsKeK_prec0123;
		case SVDID:	    return VV_CalcBasisKsKeK_prec0123;
		case SEASONID:  return SS_CalcBasisKsKeK_prec012;
		case TRENDID:   return TT_CalcBasisKsKeK_prec012;
		case OUTLIERID: return OO_CalcBasisKsKeK_prec012;
		}
	}
	else if (precPriorType == 3) {
		switch (id) {
		case SEASONID:  return SS_CalcBasisKsKeK_prec3;
		case TRENDID:   return TT_CalcBasisKsKeK_prec3;
		case OUTLIERID: return OO_CalcBasisKsKeK_prec3;
		}
	}

	return NULL;
}

#include "abc_000_warning.h"