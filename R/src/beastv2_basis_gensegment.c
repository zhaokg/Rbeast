#include <math.h> //sqrtf
#include "abc_000_macro.h" 
#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_blas_lapack_lib.h" //r_ippsSet_32 r_cblas_scopy
#include "abc_common.h"          // normalize_x_factor normalize
#include "beastv2_header.h"

 
static int TT_03(F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST * ptr) {
	#define TREND  (*((TREND_CONST*)ptr))	
	
	I32    Npad   = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32    Nseg   = seg->R2 - seg->R1 + 1L;		
	I32    Kterms = (seg->ORDER2 - seg->ORDER1) + 1; //the trend order starts from 0: 0, 1, 2, ..	
	F32PTR TERMS  = TREND.TERMS + seg->ORDER1 * N + seg->R1 - 1L;
	F32    scale  = TREND.INV_SQR[(Nseg)-1];

	r_ippsSet_32f(0, X, Kterms * Npad);
		
	//...Calcuate new terms........		
	I32  k = 0;
	for (I32 order = seg->ORDER1; order <= seg->ORDER2; order++) {

		if (order == 0) {
			//If there is only a global trend (.tNum=1) & it'sY a consant term (j=1), do not normalize it.
			//if (numOfSeg != 1) f32_normalize_inplace(Xt, N);
			//if (numOfSeg != 1) r_ippsMulC_32f_I(scale, Xt + r1 - 1, segLength);
			r_ippsSet_32f(scale, X + (seg->R1)-1, Nseg);
		} else if (order == 1) {			
			F32      b = TREND.COEFF_B[Nseg - 1];
			F32      a = TREND.COEFF_A[Nseg - 1];
			f32_seq(X+ seg->R1 -1,a,b,Nseg);
		} else { // higer oder than 1th
			//f32_normalize_inplace(Xt + r1 - 1, segLength );r_cblas_sscal(segLength , scale, Xt + r1 - 1, 1);
			r_cblas_scopy(Nseg, TERMS, 1, X + seg->R1 - 1, 1);
			f32_normalize_x_factor_inplace(X+ seg->R1 -1, Nseg, scale);
		}

		k++;
		X     += Npad;
		TERMS += N;
	} //for (rI32 j = 1; j <= ORDER; j++)

	return k;
	/*

	//...Calcuate new terms........			
#define TREND  (*((TREND_CONST*)ptr))
	F32PTR TREND_TERMS = TREND.TERMS ;
	I32  Npad = ((N + 7L) / 8L) * 8L;
	I32  Nseg = r2 - r1 + 1L;
	//the trend order starts from 0: 0, 1, 2, ..
	I32 Kterms = (order2 - order1) + 1;
	r_ippsSet_32f(0, X, Kterms * Npad);
 
		// 1..(1)....(2)...(s-1)...N
 
	I32 k = 0;
		int TORDER = order2 + 1;
		int segLength = r2 - r1 + 1;
		for (int j = order1+1; j <= TORDER; j++)
		{
 
			F32 scale = sqrtf((F32)N / (F32)segLength);
			// if there is only a global trend (.tNum=1) and it is the consant term (j=1), do not normalize it.
			if (j == 1)			{
				f32_fill_val(1., X + k * Npad + r1 - 1, segLength);
				f32_mul_val_inplace(scale, X + k * Npad + r1 - 1, segLength);
			}
			else
			{				//if (numOfSeg != 1) f32_normalize_inplace(Xt_mars + (k - 1)*Npad, N);
				f32_seq(X + k * Npad + r1 - 1, 1,1,segLength);
				f32_normalize_inplace(X + k * Npad + r1 - 1, segLength);
				f32_mul_val_inplace(scale, X + k * Npad + r1 - 1, segLength);
			}
			
			k++;
 
		}
 

		return k;
		*/
	#undef TREND
}
static int TT_1( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {
	#define TREND  (*((TREND_CONST*)ptr))	
	
	I32    Npad   = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32    Nseg   = seg->R2 - seg->R1 + 1L;		
	I32    Kterms = (seg->ORDER2 - seg->ORDER1) + 1; //the trend order starts from 0: 0, 1, 2, ..	
	F32PTR TERMS  = TREND.TERMS + seg->ORDER1 * N + seg->R1 - 1L;
	F32    scale  = TREND.INV_SQR[(Nseg)-1];

	r_ippsSet_32f(0, X, Kterms * Npad);
		
	//...Calcuate new terms........		
	I32  k = 0;
	for (I32 order = seg->ORDER1; order <= seg->ORDER2; order++) {
		r_cblas_scopy(Nseg, TERMS, 1, X + seg->R1 - 1, 1);
		// if there is only a global trend (.tNum=1) and it is the consant term (j=1), do not normalize it.
		if (order != 0) {
			/*
			F32 sum;
			r_ippsSum_32f(X + seg->R1 - 1, Nseg, &sum, ippAlgHintAccurate);
			r_ippsSubC_32f_I(sum / Nseg, X + seg->R1 - 1, Nseg); //centering the data vector
			*/
			r_ippsSubC_32f_I( *(X + seg->R1 - 1), X + seg->R1 - 1, Nseg); //put the start to zero
		}
		k++;
		X     += Npad;
		TERMS += N;
	} //for (rI32 j = 1; j <= ORDER; j++)

	return k;
#undef TREND
}
static int TT_2( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {
	#define TREND  (*((TREND_CONST*)ptr))	
	
	I32    Npad   = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32    Nseg   = seg->R2 - seg->R1 + 1L;		
	I32    Kterms = (seg->ORDER2 - seg->ORDER1) + 1; //the trend order starts from 0: 0, 1, 2, ..	
	F32PTR TERMS  = TREND.TERMS + seg->ORDER1 * N + seg->R1 - 1L;
	F32    scale  = TREND.INV_SQR[(Nseg)-1];

	r_ippsSet_32f(0, X, Kterms * Npad);
		
	//...Calcuate new terms........		
	I32  k = 0;
	for (I32 order = seg->ORDER1; order <= seg->ORDER2; order++) {
 		r_cblas_scopy(Nseg, TERMS, 1, X + seg->R1 - 1, 1);
		// if there is only a global trend (.tNum=1) and it is the consant term (j=1), do not normalize it.
		if ( Nseg != N || order != 0) f32_normalize_inplace(X, N);

 		k++;
		X     += Npad;
		TERMS += N;
	} //for (rI32 j = 1; j <= ORDER; j++)

	return k;
#undef TREND
}

static int DD_0( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {

	// seg->ORDER1 and ORDER2 are not used
	#define dummy (*((DUMMY_CONST*)ptr))

	I32  Npad   = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32  Nseg   = seg->R2 - seg->R1 + 1L;
	I32  period = dummy.period;

    //the dummy order must be less than or equal to period
	I32 Kterms = Nseg >= period ? period : Nseg;

	//...Calcuate new terms........			
	r_ippsSet_32f(0, X, Kterms * Npad);

	F32PTR SCALE        = dummy.SQRT_N_div_n;
	I32    nFullPeriod  = Nseg / period;
	I32    nRemainder   = Nseg - nFullPeriod*period;
	for (I32 j = 1; j <= Kterms; j++) {

        // 1,2,..,nRemainder|--PERIOD1--|--PERIOD2--|--PERIOD3-|..|--PERIOD-nFullperiod-|
		// a, 0,0,0,..      |0000000a000|0000000a000|0000000a000|0000000a000|
		// 0, b,0,0,..      |00000000b00|00000000b00|00000000b00|00000000b00|
		I32 n    =  (j <= nRemainder) ? (nFullPeriod + 1) : nFullPeriod;
		F32 value = SCALE[n];

		for (I32 r = seg->R1+(j-1) - 1; r<seg->R2-1; r += period) {		
			X[r] = value;		
		}
		X += Npad;
	} //for (rI32 j = 1; j <= ORDER; j++)


	return Kterms;
#undef dummy
}
static int VV_0( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {

#define SVD  (*((SVD_CONST*)ptr))

	//...Calcuate new terms........	
	I32     Npad = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32     Nseg = seg->R2 - seg->R1 + 1L;
	//the season order starts from 1: 
	I32 kTerms    = (seg->ORDER2 - seg->ORDER1) + 1 ;
	r_ippsSet_32f(0, X, kTerms * Npad);

	F32PTR TERM       = SVD.TERMS     + N * (seg->ORDER1 - 1)  + seg->R1 - 1;
	F32PTR svd_csum   = SVD.SQR_CSUM  + 1L + (N + 1) * (seg->ORDER1 - 1);
	// !L is needed bcz each col is  padded with an exta zero at the beginning so 
	// that season_csum[r2]-season_csum[r1-1] will be working

	I32    k = 0;
	for (I32 j = seg->ORDER1; j <= seg->ORDER2; j++){
		// there is only one term per order
		r_cblas_scopy(Nseg, TERM, 1L, X + seg->R1 - 1, 1L);
		F32 csum_diff     = (svd_csum[seg->R2 - 1] - svd_csum[(seg->R1 - 1) - 1]);
		F32 scalingFactor = csum_diff==0?0.0: sqrtf(N / csum_diff);  // Some hihger SVD terms may be zeros
		//r_cblas_sscal(Nseg, scalingFactor, X + seg->R1 - 1, 1L);
		
		k    += 1;
		TERM += N ;
		X    += Npad;
		svd_csum += (N + 1L) ;	
	}

	
	return k;
#undef SEASON
}

static int SS_0( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {

	#define SEASON  (*((SEASON_CONST*)ptr))

	//...Calcuate new terms........	
	I32     Npad = ((N + 7L) / 8L) * 8L;  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32     Nseg = seg->R2 - seg->R1 + 1L;
	//the season order starts from 1: sin/cos
	I32 kTerms = ((seg->ORDER2 - seg->ORDER1) + 1) * 2;
	r_ippsSet_32f(0, X, kTerms * Npad);

	F32PTR TERM        = SEASON.TERMS    + N * (seg->ORDER1-1) * 2 +  seg->R1 - 1;
    F32PTR season_csum = SEASON.SQR_CSUM + 1L + (N + 1) * (seg->ORDER1-1)*2;
	// !L is needed bcz each col is  padded with an exta zero at the beginning so 
	// that season_csum[r2]-season_csum[r1-1] will be working

	I32    k = 0;
	for (I32 j = seg->ORDER1; j <= seg->ORDER2; j++) 	{
		//memset(Xt, 0, (r1 - 1)*sizeof(F32));
		//f32_normalize_inplace(Xt + r1 - 1, segLength );r_cblas_sscal(segLength , scale, Xt + r1 - 1, 1);
		//memset(Xt + (r2 + 1) - 1, 0, (N - r2)*sizeof(F32));

		//r_cblas_scopy(segLength, tmpTERM, 1, Xt + r1 - 1, 1);						
		//f32_normalize_x_factor_inplace(Xt + r1 - 1, segLength, scale);

		/*
		r_cblas_scopy(segLength,  tmpTERM, 1L, Xt + r1 - 1, 1); //r_ippsMulC_32f_I(0.1*scale, Xt + r1 - 1,segLength);
		//rF32 dotval = DOT(segLength, Xt + r1 - 1, Xt + r1 - 1);
		//r_ippsMulC_32f_I(sqrt(N/dotval), Xt + r1 - 1, segLength);
		r_ippsMulC_32f_I(scale, Xt + r1 - 1, segLength);
		*/
		F32   scalingFactor;
		//sin
		r_cblas_scopy(Nseg, TERM, 1L, X + seg->R1 - 1, 1L);
		scalingFactor = sqrtf(N / (season_csum[seg->R2 - 1] - season_csum[(seg->R1 - 1) - 1]));
		r_cblas_sscal(Nseg, scalingFactor, X + seg->R1 - 1, 1L);
		//cos
		r_cblas_scopy(Nseg, TERM + N, 1L, (X + Npad) + seg->R1 - 1, 1L);
		scalingFactor = sqrtf(N / (season_csum[(N + 1) + seg->R2 - 1] - season_csum[(N + 1) + (seg->R1 - 1) - 1]));
		r_cblas_sscal(Nseg, scalingFactor, (X + Npad) + seg->R1 - 1, 1L);

		k += 2;
		TERM += N * 2;
		X += Npad * 2;
		season_csum += (N + 1L) * 2;
	}

	return k;
#undef SEASON
}
static int SS_1( F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {

	#define SEASON  (*((SEASON_CONST*)ptr))

	//...Calcuate new terms........	
	I32     Npad = ((N + 7L) / 8L) * 8L; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32     Nseg = seg->R2 - seg->R1 + 1L;
	//the season order starts from 1: sin/cos
	I32 kTerms = ((seg->ORDER2 - seg->ORDER1) + 1) * 2;
	r_ippsSet_32f(0, X, kTerms * Npad);

	F32PTR TERM        = SEASON.TERMS    + N * (seg->ORDER1-1) * 2 +  seg->R1 - 1;
    F32PTR season_csum = SEASON.SQR_CSUM + 1L + (N + 1) * (seg->ORDER1-1)*2;
	// !L is needed bcz each col is  padded with an exta zero at the beginning so 
	// that season_csum[r2]-season_csum[r1-1] will be working

	I32    k = 0;
	for (I32 order = seg->ORDER1; order <= seg->ORDER2; order++)
	{
		r_cblas_scopy(Nseg, TERM,   1, X+seg->R1 - 1,       1);
		r_cblas_scopy(Nseg, TERM+N, 1, X+Npad+seg->R1 - 1, 1);

		k    += 2;
		TERM += N*2;
		X    += Npad*2;
		season_csum += (N+1L) * 2;
	}

	return k;
#undef SEASON
}

/*
// Get rid of this basis type because it is almost equivalent to SS_0 and also because
// it makes less sense to introduce an offset shift to the sin/cos terms
static int SS_2(F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* ptr) {

	#define SEASON  (*((SEASON_CONST*)ptr))

	//...Calcuate new terms........	
	I32     Npad = ((N + 7L) / 8L) * 8L; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32     Nseg = seg->R2 - seg->R1 + 1L;
	//the season order starts from 1: sin/cos
	I32 kTerms = ((seg->ORDER2 - seg->ORDER1) + 1) * 2;
	r_ippsSet_32f(0, X, kTerms * Npad);

	F32PTR TERM        = SEASON.TERMS    + N * (seg->ORDER1-1) * 2 +  seg->R1 - 1;
    F32PTR season_csum = SEASON.SQR_CSUM + 1L + (N + 1) * (seg->ORDER1-1)*2;
	// !L is needed bcz each col is  padded with an exta zero at the beginning so 
	// that season_csum[r2]-season_csum[r1-1] will be working

	I32    k = 0;
	for (I32 order = seg->ORDER1; order <= seg->ORDER2; order++) {  
		r_cblas_scopy(Nseg, TERM,     1, X + seg->R1 - 1,        1);  		f32_normalize_inplace(X, N);
		r_cblas_scopy(Nseg, TERM + N, 1, X + Npad + seg->R1 - 1, 1);		f32_normalize_inplace(X+Npad, N);
		k    += 2;
		TERM += N * 2;
		X    += Npad * 2;
		season_csum += (N + 1L) * 2;
	}

	return k;
#undef SEASON
}
 */

static int OO_0(F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* bConst) {
	I32 Npad        = ((N + 7L) / 8L) * 8L; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 Kterms      = 1L; 
	I32 knotOutlier = seg->outlierKnot;

	r_ippsSet_32f(0, X, Kterms * N);	
	X[knotOutlier - 1] =bConst->outlier.SQRTN;;
	return 1;
}
static int OO_1(F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* bConst) {
	I32 Npad        = ((N + 7L) / 8L) * 8L; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 Kterms      = 1L;
	I32 knotOutlier = seg->outlierKnot;

	r_ippsSet_32f(0, X, Kterms * N);
	X[knotOutlier - 1] = 1.0;
	return 1;
}
/*
//Get rid of this one becaust it introduces a global shift to the  outlier term
static int OO_2(F32PTR X, I32 N, BEAST2_BASESEG_PTR seg, BASIS_CONST* bConst) {

	I32 Npad   = ((N + 7L) / 8L) * 8L; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 Kterms = 1;

	F32 sqrt_n1     = bConst->outlier.SQRTN_1;
	I32 knotOutlier = seg->outlierKnot;

	r_ippsSet_32f(-1/ sqrt_n1, X, Kterms * N);
	X[knotOutlier-1] = sqrt_n1; 
 
	return 1;
}
*/
pfnGenTerms Get_GenTerms(I08 id, BEAST2_PRIOR_PTR prior) {
	switch (id) {
	case DUMMYID:
		 return DD_0; 
	case SVDID:
		 return VV_0;
	case SEASONID:
		if      (prior->seasonBasisFuncType==0)		    return SS_0;
		else if (prior->seasonBasisFuncType == 1) 		return SS_1;
		//else if (prior->seasonBasisFuncType == 2)		return SS_2; 
	case TRENDID:  
		if		(prior->trendBasisFuncType == 0)		return TT_03;
		else if (prior->trendBasisFuncType == 1)		return TT_1;
		else if (prior->trendBasisFuncType == 2)		return TT_2;
		else if (prior->trendBasisFuncType == 3)		return TT_03;
		
	case OUTLIERID: 
		if		(prior->outlierBasisFuncType == 0)		return OO_0;
		else if (prior->outlierBasisFuncType == 1)		return OO_1;			
	}
	return NULL;
}
#include "abc_000_warning.h"