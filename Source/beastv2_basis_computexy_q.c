#include "abc_000_macro.h" 
#include "abc_000_warning.h"

#include "abc_blas_lapack_lib.h" //r_cblas_sgemv
#include "beastv2_header.h"


static void ANY(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{
	X     +=basis->Kbase * Npad;
	beta  +=basis->Kbase;
	I32 K = basis->K;
	r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad,beta, 1L, 0.f,	Y, 1L);
}

static void ST(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{
	//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K_SN, 1.f, Xt_mars, Npad, MODEL.curr.beta_mean, 1L, 0.f, MEMBUF1, 1L);
	X    += basis->Kbase * Npad;
	beta += basis->Kbase;
	TKNOT_PTR KNOT = basis->KNOT;
 
	for (I32 i = 0; i < basis->nKnot + 1; i++) {
		I32 order = basis->ORDER[i];
		I32 Kseg  = (basis->type == TRENDID) ? order + 1 : order * 2;
		I32 r1	  = KNOT[i - 1];
		I32 r2	  = KNOT[i] - 1;
		I32 Nseg  = r2 - r1 + 1L;

		r_cblas_sgemv(CblasColMajor, CblasNoTrans,Nseg, Kseg, 1.0f, X+r1 - 1, Npad, beta, 1L, 0.f, Y + r1 - 1L, 1L);

		X    += Kseg * Npad;
		beta += Kseg;
	}

}

static void DD_0(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{

	memset(Y, 0, sizeof(F32) * Npad);

	X	  = X    + basis->Kbase * Npad;
	beta  = beta + basis->Kbase;
	I32 K = basis->K;

	I32		   period   = basis->bConst.dummy.period; // period is a field of the dummy basis only!
	TKNOT_PTR  KNOT     = basis->KNOT;
	I16PTR     KS       = basis->ks;
	I16PTR     KE       = basis->ke;

	int NUM_OF_SEG = basis->nKnot + 1L;  // Number of seasonal segment		
	int kCounter   = 1L;                  // A counter for basis terms
	for (I32 i = 1; i <= NUM_OF_SEG; i++) {
	 
		I32 k1 = KS[i-1]-1, k2 = KE[i-1]-1;
		I32 r1 = KNOT[(i - 1) - 1];
		I32 r2 = KNOT[(i)-1] - 1;
		for (I32 j = 0; j < k2-k1+1; j++) {
			F32 coeff = *beta++;
			F32 x0     = X[r1+j-1];
			F32 y      =x0*coeff;

			for (I32 r = r1 + j; r <= r2; r += period) {
				Y[r - 1] = y;
			}
			
			X += Npad;
		} // Loop through indiviual terms of each segment

	} //Loop through individual segments
 
}

static void OO_0(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{
	memset(Y, 0, sizeof(F32) * Npad);
	beta += basis->Kbase;
	
	TKNOT_PTR	knotList	= basis->KNOT;
	F32			sqrtN		= basis->bConst.outlier.SQRTN;
	for (I32 i = 0; i < basis->nKnot; i++) {
		Y[knotList[i] - 1] = beta[i] * sqrtN;
	}
}
static void OO_1(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{
	memset(Y, 0, sizeof(F32) * Npad);
	beta += basis->Kbase;

	TKNOT_PTR	knotList = basis->KNOT;
	for (I32 i = 0; i < basis->nKnot; i++) {
		Y[knotList[i] - 1] = beta[i] ;
	}
}
/*
static void OO_2(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad)
{

	X    += basis->Kbase * Npad;
	beta += basis->Kbase;
	I32 K = basis->K;
	
	if (K == 0)
		memset(Y, 0, sizeof(F32) * Npad);
	else
		r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad,beta, 1L, 0.f, Y, 1L);
	
}
*/
void* Get_ComputeY(I08 id, BEAST2_OPTIONS_PTR opt) {
	switch (id) {
		case DUMMYID:   return DD_0;
		case SVDID:     return ANY; //ST//TODO11
		case SEASONID:  return ANY; //ST//TODO11
		case TRENDID:   return ANY; //ST//TODO11
		case OUTLIERID: {
			if      (opt->prior.outlierBasisFuncType==0)    return OO_0;
			else if (opt->prior.outlierBasisFuncType == 1)	return OO_1;
			//else if (opt->prior.outlierBasisFuncType == 2)	return OO_2;

		}
	}
	return NULL;
}
#include "abc_000_warning.h"