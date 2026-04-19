#include "math.h"
#include "abc_000_warning.h"
#include "abc_mcmc.h"
#include "abc_vec.h"
#include "beastv2_header.h"

static void DSVT(BEAST2_BASIS_PTR basis, I32 N, BEAST2_RNDSTREAM* PRND)
{
	/*
	I32 minKnotNum = basis->prior.minKnotNum;
	I32 rndOrder   = RANDINT(basis->prior.minOrder, basis->prior.maxOrder, *(PRAND->rnd08)++);
	//basis->sY.ORDER[0] = (unsigned char)(minSeasonOrder - 1) + (unsigned char)ceilfunc((*rnd++)* (maxSeasonOrder - minSeasonOrder + 1));
	if (minKnotNum == 0){
		basis->nKnot				= 0;     
		basis->ORDER[0]				= rndOrder;
		basis->KNOT[basis->nKnot] = N + 1L;
	} else {
		basis->nKnot	= minKnotNum;     		
		I32 SEP		 = N / (minKnotNum + 1);
		I32 initKnot = 1 + SEP;
		for (I32 i = 1; i <= minKnotNum; i++) {
			basis->ORDER[i - 1]	= rndOrder;
			basis->KNOT[i-   1]    = initKnot;
			initKnot += SEP;
		}
		basis->ORDER[basis->nKnot] = rndOrder;
		basis->KNOT[basis->nKnot]  = N + 1L;
	}
	*/

	//rndOrder is not used at all for DUMMY basis but kepted here to ensure the code re-usability
	//KMAX must be appropriately chosen to accomodate the initial basis

	I32 rndOrder = ceil((basis->prior.maxOrder + basis->prior.minOrder) / 2.0); // Not needed for DUMMY

	I32 maxKnotNum = basis->prior.maxKnotNum;
	I32 minKnotNum = basis->prior.minKnotNum;
	//basis->nKnot = ceil((minKnotNum + maxKnotNum) / 2.0);   // nKnot is zero if maxKnot is zero	
	// minKnotNum is 0 by default, so nKnot is also set to 0: this is chosen to be consistent with early versions
	// But if minKnotN>0, a min num of knots must be allocated; otherwise, the program crashed.
	basis->nKnot   = minKnotNum;  
	
	I32 leftMargin = basis->prior.leftMargin;
	I32 rightargin = basis->prior.rightMargin;
	I32 minSepDist = basis->prior.minSepDist;
	
	I32 Nvalid    = (N - rightargin) - (1 + leftMargin + 1) + 1;
	I32 SEP       = Nvalid / max(basis->nKnot, 1);
	I32 initKnot  = 1 + leftMargin + 1;
	for (I32 i = 1; i <= basis->nKnot; ++i) {
		basis->ORDER[i - 1] = rndOrder;       //ORDER is allocated for DUMMY but not used at all: a waste of memoery
		basis->KNOT[i - 1]  = initKnot;
		initKnot            += SEP;
	}
	basis->ORDER[basis->nKnot] = rndOrder;

	// Added to accomodate left and aright margsins
	I32 fakeStart =  1L     + (leftMargin - minSepDist);
	I32 fakeEnd   = (N + 1) -  rightargin + minSepDist;
	//The first changepont is fixed at the ts start, and the pt has been moved forward by 1
	//  so the fixed brk has a index of -1. there is one extra chngpt fixed at N+1, so the skip lenghth here is MAX_KNOTNUM+1
		
	basis->KNOT[INDEX_FakeStart]  = fakeStart;
	basis->KNOT[INDEX_FakeEnd]    = fakeEnd;	 
 
	basis->KNOT[-1]              = 1L;	    
	basis->KNOT[basis->nKnot]    = N + 1L;

	//basis->Kbase = 0; // It must be re-calcuated in UpdateKbase
	// Get Ks andKe for individula segments of each components AND the total number of terms K
	basis->CalcBasisKsKeK_TermType(basis); 
}


 static  void OO(BEAST2_BASIS_PTR basis, I32 N, BEAST2_RNDSTREAM* PRND,  BEAST2_YINFO_PTR yInfo)
{
	basis->nKnot = basis->prior.minKnotNum;    
	//basis->sY.ORDER[0] = (unsigned char)(minSeasonOrder - 1) + (unsigned char)ceilfunc((*rnd++)* (maxSeasonOrder - minSeasonOrder + 1));
	//basis->Kbase = 0; // It must be re-calcuated in UpdateKbase

	if (basis->nKnot > 0) {	
		I32    nMissing     =  yInfo->nMissing;
		I32PTR rowsMissing  =  yInfo->rowsMissing;
	 
		/////////////////////////////////////////
		memset(basis->goodvec, 1, N);
		for (int i = 0; i < nMissing; i++) {
			basis->goodvec[rowsMissing[i]] = 0;  //rowMissing is zero-based
		}
 
		//I32  Npad16  = ((N + 15) / 16) * 16;
		//I32  goodNum = i08_sum_binvec(basis->goodvec, Npad16);
		I32  goodNum = yInfo->n;
	
		for (I32 i = 1; i <= basis->nKnot; ++i) {		
			int   randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++);
			int   newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, N, randLoc);
			basis->KNOT[i - 1] = newKnot;
			basis->goodvec[newKnot - 1] = 0;
			goodNum--;
		}	
	}


	// These are not used at all and made to provide consistenst API 
	// (e.g., around UpdateGoodVec_KnotList in beast_core)
	int fakeStart = 1;
	int fakeEnd   = N + 1L;
	basis->KNOT[INDEX_FakeStart] = fakeStart;
	basis->KNOT[INDEX_FakeEnd]   = fakeEnd;
	basis->KNOT[-1]              = 1L;
	basis->KNOT[basis->nKnot]    = N + 1L;

	// Get Ks andKe for individula segments of each components AND the total number of terms K
	basis->CalcBasisKsKeK_TermType(basis);
}


void GenarateRandomBasis(BEAST2_BASIS_PTR basis, I32 NUMBASIS, I32 N, BEAST2_RNDSTREAM* PRAND, BEAST2_YINFO_PTR yInfo)
{
	//yInfo.missingRows needed only for the outllier component

	for (I32 i = 0; i < NUMBASIS; i++) {

		I32 id = basis[i].type;
		switch (id) {
		case DUMMYID:    
		case SEASONID:    
		case SVDID:
		case TRENDID:   
			DSVT(basis + i, N, PRAND); break;
		case OUTLIERID: 
			OO(basis + i, N, PRAND, yInfo); break;
		}
	}
}
	
#include "abc_000_warning.h"