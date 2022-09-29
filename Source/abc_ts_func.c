#include "abc_000_warning.h"

#include "abc_001_config.h"
 
#include "abc_ts_func.h"
#include "abc_vec.h"   // for f32_seq only
#include "abc_math.h"  // for fastsqrt only
#include "abc_date.h"  // for fastsqrt only
#include "abc_ide_util.h"  //printf
#include "abc_common.h" //normalize quicksort
#include "abc_blas_lapack_lib.h"

#if defined(PI)
	#undef PI	
#endif
#define PI   (3.141592653589793)

void preCalc_terms_season(F32PTR SEASON_TERMS, F32PTR SEASON_SQR_CSUM, F32PTR SCALE_FACTOR, int N, F32 PERIOD, int maxSeasonOrder)
{
	// At the first run of the function, pre-calcalte terms for SEASON AND TRENDS Bases	

	// SIN or COS: vmdCos(n, a, y, mode); 
	// https: //software.intel.com/en-us/node/521751 All the VM mathematical functions can perform 
	// in - place operations, where the input and output arrays are at the same memory locations.
	//Multiply by scalar: cblas_dscal(const MKL_INT n, const double a, double *x, const MKL_INT incx);
	//r_cblas_sscal(N, freq_factor * (F32)order, ptr, 1);
	//r_cblas_scopy(N, ptr, 1, ptr + N, 1);
	// vmsSin(N, ptr, ptr, VML_EP); VML_HA
	// //vmsCos(N, ptr + N, ptr + N, VML_EP); VML_HA
	// r_ippsMulC_32f_I(1 / sqrtf(dotProduct / N), ptr + N, N);

	// Precomputing SEASONAL TERMS
	if (SEASON_TERMS == NULL) 	return;	

	F32   freq_factor = 2.0f *  3.141592653589793f / PERIOD;

	F32PTR ptr  = SEASON_TERMS;
	F32PTR ptr1 = SEASON_SQR_CSUM;
	F32PTR ptr2 = SEASON_SQR_CSUM + (N + 1);

	for (I32 order = 1; order <= maxSeasonOrder; order++)	{
		f32_seq(ptr,    1.f, 1.f, N);	
		f32_mul_val_inplace(freq_factor * (F32)order, ptr, N);
		
		f32_copy(ptr, ptr + N, N);
		                          
		f32_sincos_vec_inplace(ptr+N, ptr, N); //(sin,cos): we want the cos term first and then the sin term.
		//because for the highest order period/2, sins are all zeros.

	 
		F32 dotProduct, scale1, scale2;
		dotProduct = DOT(N, ptr, ptr);   	   scale1 = 1/sqrtf(dotProduct / N);	f32_mul_val_inplace(scale1, ptr,  N);
		dotProduct = DOT(N, ptr + N, ptr + N); scale2 = 1/sqrtf(dotProduct / N);	f32_mul_val_inplace(scale2, ptr+N, N);
		// Commented out here : needed for the computation of samp because we need sin(x) and cos(x) has 
		// an amplitude of 1.0. The scaling by dotProduct here will screw up the amplitude
		if (SCALE_FACTOR) {
			SCALE_FACTOR[(order-1)*2   ] = scale1;
			SCALE_FACTOR[(order-1)*2+1L] = scale2;
		}


		if (SEASON_SQR_CSUM) {
			*ptr1 = 0.f;
			*ptr2 = 0.f;
			f32_copy(ptr,     ptr1 + 1L, N);      f32_cumsumsqr_inplace(ptr1 + 1L, N);
			f32_copy(ptr + N, ptr2 + 1L, N);      f32_cumsumsqr_inplace(ptr2 + 1L, N);

			ptr  += 2 * N;
			ptr1 += 2 * (N + 1);
			ptr2 += 2 * (N + 1);
		}	else {
			ptr += 2 * N;
		}
		
	} // or (rI32 order = 1; order <= maxSeasonOrder; order++)

}

void preCalc_terms_trend(F32PTR TREND_TERMS, F32PTR INV_SQR, int N, int maxTrendOrder)
{	
	if (TREND_TERMS != NULL) {
		//PRECOMPUTING TREND TERMS
		F32PTR	ptr = TREND_TERMS;
		for (I32 i = 1; i <= (maxTrendOrder + 1); i++) {
			if (i == 1) // the constant term
				r_ippsSet_32f(1.0f, ptr, N); //IppStatus r_ippsSet_32f(Ipp64f val, Ipp64f* pDst, int len);	
			else {
				f32_seq(ptr, 1.0f, 1.0f, N);
				r_vsPowx(N, ptr, (F32)(i - 1), ptr);
				f32_normalize_inplace(ptr, N);
			}
			ptr += N;
		}
	}

	//Pre - compute sqrt(N / i), i = 1, ..., N, which will be used later
	if (INV_SQR != NULL) {
		F32 sqrt_N = fastsqrt((F32)N);
		for (I32 i = 0; i < N; i++)
			INV_SQR[i] = sqrt_N / fastsqrt((F32)(i + 1));
	}
}

void preCalc_XmarsTerms_extra(F32PTR COEFF_A, F32PTR COEFF_B, I32 N)
{
	//PRECOMPUTING TREND TERMS

	if (COEFF_A != NULL &&  COEFF_B != NULL) {
		COEFF_B[1 - 1] = 0;
		COEFF_A[1 - 1] = fastsqrt(N);

		for (I32 n = 2; n <= N; n++) {
			F32 sum   = (1L + n) / 2.f;
			F32 b     = 1.f / ((n + 1L) * (2L * n + 1) / 6.f - sum * sum);
			F32 bsqrt = fastsqrt(b * N / n);
			COEFF_B[n - 1] = bsqrt;
			COEFF_A[n - 1] = -bsqrt * sum;
		}
	}

}

void preCalc_XmarsTerms_extra_fmt3(F32PTR COEFF_A, F32PTR COEFF_B, I32 N)
{
	//PRECOMPUTING TREND TERMS

	if (COEFF_A != NULL &&  COEFF_B != NULL) {
		COEFF_B[1 - 1] = 0;
		COEFF_A[1 - 1] = fastsqrt(N);

		for (I32 n = 2; n <= N; n++) {
			F32 b     = 6.f / (n * (F32)(n + 1) * (F32)(2*n + 1));
			F32 bsqrt = fastsqrt(b);
			COEFF_B[n - 1] = bsqrt;
			COEFF_A[n - 1] = bsqrt;
		}
	}


}

void preCalc_scale_factor(F32PTR sclFactor, I32 N, I32 maxKnotNum, I32 minSepDist, F32PTR mem1, F32PTR mem2)
{
	if (sclFactor == NULL) {
		return;
	}

	F32 N_tmp, tmp1, tmp2;

	for (int k = 0; k <= maxKnotNum; k++) {

		N_tmp = N - (k+1)*(minSepDist - 1) - 1.f;

		if (k == 0)	{
			*mem1 = 1.0f;
			tmp1  = logf(1.0f);
		} else {			
			f32_seq(mem1, (F32)1, (F32)1, k);
			r_ippsSubC_32f_I(1.f, mem1, k);
			r_ippsSubCRev_32f_I(N_tmp, mem1, k);
			r_ippsLn_32f_I(mem1, k);
			r_ippsSum_32f(mem1, k, &tmp1, ippAlgHintAccurate);
		}

		N_tmp    = N - (k + 2)*(minSepDist - 1) - 1.f;		
		f32_seq(mem2, 1.f, 1.f, k+1);
		r_ippsSubC_32f_I(1.f, mem2, k + 1);
		r_ippsSubCRev_32f_I(N_tmp, mem2, k + 1);
		r_ippsLn_32f_I(mem2, k + 1);	
		r_ippsSum_32f(mem2, k + 1, &tmp2, ippAlgHintAccurate);

		sclFactor[k] = (N - (k + 2)*minSepDist + 1) *expf(tmp1 - tmp2);
	}	

}


void cvt_knotList_to_bincode(U08PTR  good, I32 N, U16 minSepDist, U16PTR knotList, I64 knotNum)
{
	r_ippsSet_8u(1, good, N);	
	for (int i = 1; i <= knotNum; i++)	{		
		r_ippsSet_8u(0L, good + (knotList[i-1] - minSepDist) - 1, 2*minSepDist+1);			
	}
	r_ippsSet_8u(0, good,                      (minSepDist+1) );
	r_ippsSet_8u(0, good+ (N-minSepDist+1) - 1, minSepDist);	 
}

I32 tsAggegrationPrepare(F32PTR oldTime, I32 Nold, F32 dT, I32PTR *SortedTimeIdx, I32PTR *NumPtsPerInterval,
					   I32 *startIdxOfFirsInterval, F32 *startTime)
{	 
	//Return the number of newly aggregrated ts
	/************************************************/
	I32PTR  SORTED_IDX = malloc(sizeof(I32)*Nold);	
	for (I32 i = 0; i < Nold; i++) SORTED_IDX[i] = i;
	QuickSortA(oldTime, SORTED_IDX, 0, Nold - 1);	
	*SortedTimeIdx = SORTED_IDX;
	/************************************************/

	// Now oldTime is the sorted times
	/************************************/
	F32PTR	SortedTimes = oldTime;
	F32 T0   = SortedTimes[0],     T1 = SortedTimes[Nold-1];
	I32 i0   = round(T0 / dT),	   i1 = round(T1 / dT);
	//F32 i0 = floor(T0/period);
	//T0     = i0*period + floor( (T0-i0)/dT )*dT;		
	I32 Nnew = ((i1 - i0) + 1);	
	*startTime = i0*dT; // the midpoint of the first internval
	/************************************/

	/************************************/
	I32PTR NUM_PER_INTERVAL = malloc(sizeof(I32)*Nnew);
	memset(NUM_PER_INTERVAL, 0L, sizeof(F32)*Nnew);
	*NumPtsPerInterval = NUM_PER_INTERVAL;
	/************************************/

	/************************************/
	I32 idxTime      = 0;	
	F32 UpperEndInterval = i0*dT+0.5*dT;
	while (SortedTimes[idxTime] < (UpperEndInterval - dT) && idxTime < Nold) {
		idxTime++;
	}	
	*startIdxOfFirsInterval = idxTime;
	/************************************/
	
	for (I32 i = 0; i < Nnew; i++) {		 
		I32 nptsPerInterval  = 0;
		F32 time = SortedTimes[idxTime];
		while (time <= UpperEndInterval && idxTime < Nold) {
			++nptsPerInterval;
			time = SortedTimes[++idxTime];
		}
		NUM_PER_INTERVAL[i] = nptsPerInterval;
		UpperEndInterval   += dT;
	}

	return Nnew;
}

void tsAggegrationPerform(F32PTR RegularTS, I32 Nnew, F32PTR IrregularTS, I32 Nold, I32PTR NumPerSeg,I32PTR SorteTimeIdx){
	
	F32 nan = getNaN();	
	I32 idx = 0;	
	for (I32 i = 0; i < Nnew; i++) {				
		F32 sum = 0;
		I32 num = 0;
		I32 nPts = NumPerSeg[i];
		for (I32 j = 0; j < nPts; j++) {
			I32 id = SorteTimeIdx[idx++];
			F32 Y  = IrregularTS[id];
			if (Y == Y)	{sum += Y; num++;}
		}
		RegularTS[i] = num == 0 ? nan : sum / num;
	}
}


void tsRemoveNaNs(F32PTR x, int N) {

	I32 preGoodIdx  = -1;
	I32 postGoodIdx = -1;
	for (int i = 0; i < N; ++i) {
		if (x[i]==x[i]) {
			preGoodIdx = i;
			continue;
		}

		if (postGoodIdx<=i) {
			for (int j = i + 1; j < N; j++) {
				if (x[j] == x[j]) {
					postGoodIdx = j;
					break;
				}
			}
		}

		if (preGoodIdx < 0 ) {
			if (postGoodIdx <= i) {
				break;
			}	else {
				x[i] = x[postGoodIdx];				
			}			 
		} else {
			if (postGoodIdx <= i) {
				x[i]      = x[preGoodIdx];				
			}
			else {
				x[i] = (x[preGoodIdx]*(postGoodIdx-i)+x[postGoodIdx]*(i- preGoodIdx))/(postGoodIdx-preGoodIdx);
			}
		}
		
		preGoodIdx = i;		 
	}

}

#include "abc_000_warning.h"