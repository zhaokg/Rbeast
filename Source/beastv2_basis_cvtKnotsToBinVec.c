#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include <string.h>
#include "abc_vec.h"   // for i08_sum_binvec only
#include "abc_blas_lapack_lib.h" //r_ippsSet_8u


#include "beastv2_header.h"

static int DSVT(U08PTR goodVec, I32 N, TKNOT_PTR KNOT, I64 nKnot, KNOT2BINVEC * info)
{
	I32 minSepDist = info->minSepDist;

	memset(goodVec, 1, N);
	for (int i = 1; i <= nKnot; i++) {
		r_ippsSet_8u(0L, goodVec + (KNOT[i-1] - minSepDist) - 1, 2 * minSepDist + 1);
	}

	I32 leftMargin  = info->leftMargin;
	I32 rightMargin = info->rightMargin;
	r_ippsSet_8u(0, goodVec,                            leftMargin + 1);
	r_ippsSet_8u(0, goodVec + (N- rightMargin +1L) - 1, rightMargin);

	I32  Npad16  = ((N + 15) / 16) *16;
	I32  goodNum = i08_sum_binvec(goodVec, Npad16);
	return goodNum;
}

static int OO(U08PTR goodVec, I32 N,TKNOT_PTR KNOT, I64 nKnot, KNOT2BINVEC* info)
{
	I32    nMissing    = info->yInfo->nMissing;
	I32PTR rowsMissing = info->yInfo->rowsMissing;

	/////////////////////////////////////////
	memset(goodVec, 1, N);
	for (int i = 0; i < nMissing; i++) { 
		goodVec[rowsMissing[i]] = 0;  //rowMissing is zero-based
	}
	for (int i = 0; i < nKnot; i++) {
		I32 idx       = KNOT[i];
		goodVec[idx - 1] = 0;
	}

	I32  Npad16  =  ((N + 15) / 16) *16;
	I32  goodNum = i08_sum_binvec(goodVec, Npad16);
	return goodNum;
}
 

void CvtKnotsToBinVec(BEAST2_BASIS_PTR b, I32 NUMBASIS, I32 N, BEAST2_YINFO_PTR yInfo ) {
	// nMissing & rowsMissing used for the outlier function
	KNOT2BINVEC info = {.yInfo=yInfo};

	for (int i = 0; i < NUMBASIS; i++) {

		info.minSepDist = b[i].prior.minSepDist;			
		info.leftMargin  = b[i].prior.leftMargin;  // not needed for outlier compnt
		info.rightMargin = b[i].prior.rightMargin; // not needed for outlier compnt

		switch (b[i].type) {
			case DUMMYID:  
			case SEASONID: 
			case SVDID:
			case TRENDID: 
				//CHANGE:(MODEL.b[i].goodvec, (MODEL.b[i].goodNum	
				b[i].goodNum = DSVT(b[i].goodvec, N, b[i].KNOT, b[i].nKnot, &info);
				break;
			case OUTLIERID: 
				b[i].goodNum = OO(b[i].goodvec, N, b[i].KNOT, b[i].nKnot, &info);
				break;
		}
		
	}

}

#include "abc_000_warning.h"