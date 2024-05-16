#pragma once

#include <inttypes.h>  //#include <stdint.h>

#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_mem.h"
#include "abc_blas_lapack_lib.h"
#include "beastv2_header.h"

extern F32  GetPercentileNcp(F32PTR prob, I32 N, F32 pctile);

extern void SetupPointersForCoreResults(CORESULT* coreResults, BEAST2_BASIS_PTR b, I32 NumBasis, BEAST2_RESULT* resultChain);

extern void CvtKnotsToBinVec(BEAST2_BASIS_PTR b, I32 NUMBASIS, I32 N, BEAST2_YINFO_PTR yInfo);
extern void GenarateRandomBasis(BEAST2_BASIS_PTR basis, I32 NUMBASIS, I32 N, BEAST2_RNDSTREAM* PRAND);
 
void BEAST2_EvaluateModel(
    BEAST2_MODELDATA* curmodel, BEAST2_BASIS_PTR b, F32PTR Xt_mars, I32 N, I32 NUMBASIS,
    BEAST2_YINFO_PTR  pyInfo, BEAST2_HyperPar* hyperPar, F32PTR precVec, VOID_PTR stream);

void MatxVec(BEAST2_BASESEG* SEG, I32 numSeg, F32PTR X, F32PTR Y, F32PTR XtY, I32 N);
I32  GetInfoBandList(BEAST2_BASESEG* info,      BEAST2_MODEL_PTR model, I32 Klastcol);
I32  GetInfoBandList_post(BEAST2_BASESEG* info, BEAST2_MODEL_PTR model, I32 Kstartcol);

void MatxMat(BEAST2_BASESEG* infoX, I32 numBandsX, F32PTR X, BEAST2_BASESEG* infoY,
             I32 numBandsY, F32PTR Y, F32PTR XtY, I32 N, I32 Knew);
void XtX_ByGroup(BEAST2_BASESEG* SEG, I32 numSeg, F32PTR X, F32PTR XtX, I32 N, I32 Knew);

static INLINE int GetTotalColNumOfXmars(BEAST2_BASIS_PTR b, I32 NUM_BASIS) {
    I32 K;
    if      (NUM_BASIS == 1)     K = b[0].K;
    else if (NUM_BASIS == 2)     K = b[0].K+ b[1].K;
    else if (NUM_BASIS == 3)     K = b[0].K+ b[1].K + b[2].K;
    else   {
        K=0;
        for (int i=0; i < NUM_BASIS; i++)    K += b[i].K;
    }
    return K;
}
static INLINE void UpdateBasisKbase(BEAST2_BASIS_PTR b, I32 NUMBASIS, int startBasisIdx0) {

    if (NUMBASIS == 1) 
        b[0].Kbase = 0;    
    else if (NUMBASIS == 2) {
        if (startBasisIdx0 == 0) {
           // b[0].Kbase = 0;
             b[1].Kbase = b[0].K;
        }  
    } 
	else if (NUMBASIS == 3) {
        if (startBasisIdx0 == 0) {
            b[1].Kbase = b[0].K;
            b[2].Kbase = b[1].Kbase + b[1].K;
        } else if (startBasisIdx0 == 1) {
            b[2].Kbase = b[1].Kbase + b[1].K;
        } else if (startBasisIdx0 == 2) {
        }          
    }
	/*
	* commented out bcz currently, MAX_NUMBAIS is set to 3 (e.g., S+T+O)
	else if (NUMBASIS > 3)  {           
        for (rI32 i = startBasisIdx0 +1; i < MAX_NUM_BASIS; i++) {
            b[i].Kbase  = b[i-1].Kbase + b[i-1].K;
        }   
    }
	*/
}


 
 static void MatxVec_FULL(BEAST2_BASESEG* SEG, I32 numSeg, F32PTR X, F32PTR Y, F32PTR XtY, I32 N)
{
	// New components
	 I32 Kterms = 0;
	 for (int i = 0; i < numSeg; i++) {			 Kterms += SEG[i].K; 	 }
	 I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
		// cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, 1, N, 1, X_mars_prop + (k1_new - 1)*N, N, Y, N, 0, buff1, K_newTerm);								  

	//?????????????????????????????????????????????????????? Npad->N
	r_cblas_sgemv(CblasColMajor, CblasTrans,
		N, Kterms, 1.f,
		X , Npad,
		Y, 1L, 0.f,
		XtY, 1L);	//MEMBUF1, 1L);
	//r_cblas_scopy(K_newTerm, MEMBUF1, 1, XtY_prop + k1_new - 1, 1);

 


}
 static  void MatxMat_FULL(BEAST2_BASESEG*infoX,I32 numBandsX, F32PTR X,BEAST2_BASESEG*infoY,I32 numBandsY, F32PTR Y, F32PTR XtY, I32 N, I32 Knew)
{

	I32 Kterms1 = 0;
	I32 Kterms2 = 0;
	for (int i = 0; i < numBandsX; i++) { Kterms1 += infoX[i].K; }
	for (int i = 0; i < numBandsY; i++) { Kterms2 += infoY[i].K; }

	I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
 

	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		Kterms1, Kterms2, N, 1.0f, \
		X, Npad,
		Y, Npad, 0.f, \
		XtY, Knew); //0.f, MEMBUF1, k1_new - 1);	
 
}
static void XtX_ByGroup_FULL(BEAST2_BASESEG* SEG, I32 numSeg,F32PTR X,F32PTR XtX,I32 N, I32 Knew) {
	/*
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	NEW.Knewterm, NEW.Knewterm, Nseg, 1.0,
	Xnewterm + r1 - 1, Npad,
	Xnewterm + r1 - 1, Npad, 0.f,
	MODEL.prop.XtX + (NEW.k1 - 1) * Knew + NEW.k1 - 1, Knew);// MEMBUF2, K_newTerm);
	*/

	I32 Kterms = 0;
	for (int i = 0; i < numSeg; i++) { Kterms += SEG[i].K; }
	I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv

	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
					Kterms, Kterms, N, 1.0,
					X , Npad,
					X , Npad, 0.,
					XtX, Knew);

}


void MR_EvaluateModel(
	BEAST2_MODELDATA* curmodel, BEAST2_BASIS_PTR b, F32PTR Xt_mars, I32 N, I32 NUMBASIS,
	BEAST2_YINFO_PTR pyInfo, BEAST2_HyperPar* hyperPar, F32PTR precVec, VOID_PTR stream);