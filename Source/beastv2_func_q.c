#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "abc_001_config.h"
#include "abc_rand.h"
#include "abc_mat.h"
#include "abc_math.h"  //sum_log_diag_v2
#include "abc_vec.h"  //sum_log_diag_v2
#include "beastv2_func.h"

/*
F32  GetPercentileNcp_old(F32PTR prob, I32 N, F32 pctile) {
// THis versio is abandoned bcz it can give riduous results
// For example, p=[0, 1.0] will give a median of 0.5 if doing
// interpoltation. The actual median will be 1.0
	F32 preNcp = 0;
	F32 ncp    = 0;
	F32 cumProb=0.f;
	for (int i = 0; i < N; i++) {
		cumProb += prob[i];
		if (cumProb > pctile) {
			// Linear interpolation betweren prevNcp and i
			F32 d2 = cumProb - pctile;
			F32 d1 = pctile - (cumProb-prob[i]);
			ncp = (preNcp*d2  + i*d1)/prob[i];
			break;
		}
		preNcp = i;
	}
	return ncp;
}
*/

F32  GetPercentileNcp(F32PTR prob, I32 N, F32 pctile) {
 
	F32 cumProb=0.f;
	for (int i = 0; i < N; i++) {
		cumProb += prob[i];
		if (cumProb > pctile) {
			return i;		
		}	 
	}
	
	return N - 1; //this branch should never be visited (or only if pctile >1.0)
}

void SetupPointersForCoreResults(CORESULT* coreResults, BEAST2_BASIS_PTR b, I32 NumBasis, BEAST2_RESULT* resultChain) {
	for (I32 i = 0; i < NumBasis; i++) {
		if (b[i].type == SEASONID || b[i].type == DUMMYID || b[i].type == SVDID) {
			coreResults[i].xNProb = resultChain->sncpPr,
			coreResults[i].xProb  = resultChain->scpOccPr,
			coreResults[i].xorder = resultChain->sorder, //sorder=NULL if extra.computerSeasonOrder=0
			coreResults[i].x      = resultChain->sY,
			coreResults[i].xSD    = resultChain->sSD;
		}
		else if (b[i].type == TRENDID) {
			coreResults[i].xNProb = resultChain->tncpPr,
			coreResults[i].xProb  = resultChain->tcpOccPr,
			coreResults[i].xorder = resultChain->torder, //torder=NULL if extra.computerTrendOrder=0
			coreResults[i].x      = resultChain->tY,
			coreResults[i].xSD    = resultChain->tSD;
		}
		else if (b[i].type == OUTLIERID) {
			coreResults[i].xNProb = resultChain->oncpPr,
			coreResults[i].xProb  = resultChain->ocpOccPr,
			coreResults[i].xorder = NULL,				//torder=NULL if extra.computerSeasonOrder=0
			coreResults[i].x      = resultChain->oY,
			coreResults[i].xSD    = resultChain->oSD;
		}
	}
}

void BEAST2_EvaluateModel(
	BEAST2_MODELDATA *curmodel, BEAST2_BASIS_PTR b, F32PTR Xt_mars, I32 N, I32 NUMBASIS,
	BEAST2_YINFO_PTR  yInfo,    BEAST2_HyperPar *hyperPar, F32PTR precVec, VOID_PTR stream )
{
	//TO RUN THIS FUNCTION, MUST FIRST RUN CONVERTBAIS so to GET NUMBERS OF BASIS TERMS	
	
	I32 Npad	= (I32)ceil((F32)N / 8.0f) * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 K		= 0;	 
	for (I32 basisID = 0; basisID < NUMBASIS; basisID++) {
		BEAST2_BASIS_PTR basis = b + basisID;
		if (basis->type != OUTLIERID) {
			int         NUM_SEG = basis->nKnot + 1;
			TKNOT_PTR   KNOT	= basis->KNOT;
			TORDER_PTR  ORDER   = basis->ORDER;

			BEAST2_BASESEG seg;
			seg.ORDER1 = basis->type == TRENDID ? 0 : 1;
			//For season terms: 1..(1)....(2)...(|sSegNum-1)...N (|N+1)
			for (int i = 1; i <= NUM_SEG; i++) {
				///The first segment starts from 1; other segments start one interval after the previous segment
				//The last segment ends at N'; others end exactly at the current breakpoint minus one;
				//the last bk is (N+1)
				seg.R1     = KNOT[(i - 1) - 1L];
				seg.R2     = KNOT[i - 1L] - 1L;
				seg.ORDER2 = basis->type==DUMMYID? 0 :  ORDER[i - 1L];
				I32 k  = basis->GenTerms(Xt_mars + Npad * K, N, &seg, &(basis->bConst));
				K += k;
			}
		} 	else	{
			int         numOfSeg = basis->nKnot;
			TKNOT_PTR   knotList = basis->KNOT;

			BEAST2_BASESEG seg;
			seg.ORDER1 = seg.ORDER2 = 0; // orders are not used at all
			//For season terms: 1..(1)....(2)...(|sSegNum-1)...N (|N+1)
			for (int i = 1; i <= numOfSeg; i++) {
				///The first segment starts from 1; other segments start one interval after the previous segment
				//The last segment ends at N'; others end exactly at the current breakpoint minus one;
				//the last bk is (N+1)
				seg.R1  = knotList[(i)-1L];
				seg.R2  = knotList[(i)-1L];
				I32 k = basis->GenTerms(Xt_mars + Npad * K, N, &seg, &(basis->bConst));
				K += k;
			}
		}


	}
	curmodel->K = K;

	//--------------------------------------------------------------------------------------------------
	// Set those rows of X_mars specfied by rowsMissing  to zeros 
	//	TODO: this could be buggy: Xt_mars may be not large enough to backup the X values at the missing rows
	// for the full initial model
	F32PTR	GlobalMEMBuf  = Xt_mars + K * Npad;
	F32PTR	Xt_zeroBackup = GlobalMEMBuf;
	if (yInfo->nMissing > 0) {
		F32 fillvalue = 0.f;
		f32_mat_multirows_extract_set_by_scalar(Xt_mars, Npad, K, Xt_zeroBackup, yInfo->rowsMissing, yInfo->nMissing, fillvalue);
	}

	//Calcuate X'*X and Calcuate X'*Y
	//DGEMM('T', 'N', &m, &n, &K, &alpha, X_mars, &K, X_mars, &K, &beta, XtX, &m);	
	//cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);

	
	F32PTR XtX = curmodel->XtX;
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.f, Xt_mars, Npad, Xt_mars, Npad, 0.f, XtX, K);
	//cblas_dgemmt (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
	//cblas_dgemmt(CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, k, N, 1, X_mars, N, X_mars, N, 0, XtY, k);				

	F32PTR XtY = curmodel->XtY;
	//cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);
	r_cblas_sgemv(CblasColMajor, CblasTrans, Npad, K, 1, Xt_mars, Npad, yInfo->Y, 1, 0, XtY, 1);

	// Restore the backuped zeros at the missing rows
	if (yInfo->nMissing > 0) {
		f32_mat_multirows_set_by_submat(Xt_mars, Npad, K, Xt_zeroBackup, yInfo->rowsMissing, yInfo->nMissing);
	}
	// Now GlobalMEMBuf_1st is free to use;

	//X_mars has been constructed. Now use it to calcuate its margin likelihood
	//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec)	
	F32PTR cholXtX = curmodel->cholXtX;
	//f32_copy( XtX,  cholXtX, K * K);
	//f32_add_val_matrixdiag(cholXtX, precVec[0], K);

	//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
	F32PTR beta_mean = curmodel->beta_mean;	
	///////////////////////////////////////////////////////////////
	/*	
		//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
		r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, cholXtX, K); // Choleskey decomposition; only the upper triagnle elements are used
		//LAPACKE_spotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
		r_cblas_scopy(K, XtY, 1, beta_mean, 1);
		r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1, cholXtX, K, beta_mean, K);	
	*/
	
	chol_addCol_skipleadingzeros_prec_invdiag(XtX, cholXtX, precVec, K, 1, K);
	solve_U_as_LU_invdiag_sqrmat(cholXtX, XtY, beta_mean, K);
	//////////////////////////////////////////////////////////////////


	//Compute beta = beta_mean + Rsig2 * randn(p, 1);
	//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
	//status = vdRngGaussian( method, stream, n, r, a, sigma );

	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/
	/*
	F32PTR beta = model->beta;
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, (*(VSLStreamStatePtr *)stream), K, beta, 0, 1);

	// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
	{
		//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, post_P_U, K, beta, K);
		solve_U_as_U_invdiag(cholXtX, beta, K, K);
	}
	r_ippsMulC_32f_I(sqrtf(model->sig2), beta, K);
	r_ippsAdd_32f_I(beta_mean, beta, K);
	*/
	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/

	//Compute alpha2_star	 
	/*
	r_cblas_scopy(K, beta_mean, 1, GlobalMEMBuf_1st, 1);
	//cblas_dtrmv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE trans, const CBLAS_DIAG diag, const MKL_INT n, const double *a, const MKL_INT lda, double *x, const MKL_INT incx);
	r_cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, K, cholXtX, K, GlobalMEMBuf_1st, 1);
	F32 alpha2_star = pyInfo->YtY - DOT(K, GlobalMEMBuf_1st, GlobalMEMBuf_1st);
	 */
	F32 alpha2_star = (yInfo->YtY_plus_alpha2Q[0] - DOT(K,  XtY,  beta_mean)) * 0.5;

	//half_log_det_post; = sum(log(1. / diag(U)))
	F32 half_log_det_post = sum_log_diagv2(cholXtX, K);
	
	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));		
	F32 half_log_det_prior	= -.5f * K*logf(precVec[0]);

	//log_ML    = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik = half_log_det_post - half_log_det_prior - yInfo->alpha1_star * logf(alpha2_star);

	curmodel->alpha2Q_star[0] = alpha2_star;
	curmodel->marg_lik        = marg_lik;

	return;
}

void MatxVec(BEAST2_BASESEG* SEG, I32 numSeg, F32PTR X, F32PTR Y, F32PTR XtY, I32 N)
{
	I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	// New components
	for (int i=1; i<=numSeg; i++) {  		
		// cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, 1, N, 1, X_mars_prop + (k1_new - 1)*N, N, Y, N, 0, buff1, K_newTerm);								  
		I32 r1		= SEG[i -1].R1;
		I32 r2		= SEG[i -1].R2;
		I32 Nseg	= r2 - r1 + 1;
		I32 Kseg	= SEG[i - 1].K;
		//?????????????????????????????????????????????????????? Npad->N
		r_cblas_sgemv(CblasColMajor, CblasTrans,
			Nseg,      Kseg,   1.f,
			X+r1 - 1,  Npad,
			Y+r1 - 1,  1L, 0.f,
			XtY, 1L);	//MEMBUF1, 1L);
		//r_cblas_scopy(K_newTerm, MEMBUF1, 1, XtY_prop + k1_new - 1, 1);

		X   += Kseg*Npad;
		XtY += Kseg;
	}

 
}

I32  GetInfoBandList(BEAST2_BASESEG* info, BEAST2_MODEL_PTR model, I32 Klastcol) {
	I32 numBands = 0;
	I32 QUITFLAG = 0;
	for (int basisID = 0; basisID < model->NUMBASIS; basisID++) {

		BEAST2_BASIS_PTR b = model->b + basisID;

		if (b->type != OUTLIERID)
		{
			I32 numSeg = b->nKnot + 1;
			for (int j = 0; j < numSeg; j++) {
				I32 Kbase = b->Kbase;
				if (Klastcol >= (Kbase + b->ks[j])) {
					info->R1 = b->KNOT[j - 1];
					info->R2 = b->KNOT[j] - 1L;
					info->K = min(Kbase + b->ke[j], Klastcol) - (Kbase + b->ks[j]) + 1;
					info++;
					numBands++;
				}
				else {
					QUITFLAG = 1;	break;
				}
			}
		}
		else
		{
			I32 numSeg = b->nKnot;
			for (int j = 0; j < numSeg; j++) {
				I32 Kbase = b->Kbase;
				if (Klastcol >= (Kbase + b->ks[j])) {
					info->R1 = b->KNOT[j];
					info->R2 = b->KNOT[j] ;
					info->K  = min(Kbase + b->ke[j], Klastcol) - (Kbase + b->ks[j]) + 1;
					info++;
					numBands++;
				}
				else {
					QUITFLAG = 1;	break;
				}
			}
		}


		if (QUITFLAG == 1) { break; }
	}
	return numBands;
}
I32  GetInfoBandList_post(BEAST2_BASESEG* info, BEAST2_MODEL_PTR model, I32 Kstartcol) {
	I32 numBands = 0;
	for (int basisID = 0; basisID < model->NUMBASIS; basisID++) {

		BEAST2_BASIS_PTR b = model->b + basisID;
		if (b->type != OUTLIERID)
		{
			for (int j = 0; j < b->nKnot + 1; j++) {
				I32  Kbase = b->Kbase;
				if (Kstartcol <= (Kbase + b->ke[j])) {
					info->R1 = b->KNOT[j - 1];
					info->R2 = b->KNOT[j] - 1L;
					info->K = (Kbase + b->ke[j]) - max(Kbase + b->ks[j], Kstartcol) + 1;
					info++;
					numBands++;
				}
			}
		}
		else
		{
			for (int j = 0; j < b->nKnot; j++) {
				I32  Kbase = b->Kbase;
				if (Kstartcol <= (Kbase + b->ke[j])) {
					info->R1 = b->KNOT[j];
					info->R2 = b->KNOT[j];
					info->K = (Kbase + b->ke[j]) - max(Kbase + b->ks[j], Kstartcol) + 1;
					info++;
					numBands++;
				}
			}
		}


	}
	return numBands;
}

void MatxMat(BEAST2_BASESEG*infoX,I32 numBandsX, F32PTR X,BEAST2_BASESEG*infoY,I32 numBandsY, F32PTR Y, F32PTR XtY, I32 N, I32 Knew)
{
	I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv

	I32 Ksegcsum = 0;
	for (int i = 0; i < numBandsY; i++) {

		I32 new_r1 = infoY[i].R1;
		I32 new_r2 = infoY[i].R2;
		I32 Kseg   = infoY[i].K;

		I32 Kbandcsum = 0;	
		for (int j = 0; j < numBandsX; j++) {
			I32 old_r1 = infoX[j].R1;
			I32 old_r2 = infoX[j].R2;
              
			I32 r1 = max(new_r1, old_r1);
			I32 r2 = min(new_r2, old_r2);

			I32 Kband = infoX[j].K;
			if (r2>= r1) {
				I32 Nseg = r2 - r1 + 1;				
				
				r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
					Kband, Kseg, Nseg, 1.0f, \
					X + r1 - 1L, Npad,
					Y + r1 - 1L, Npad, 0.f, \
					XtY + Kbandcsum, Knew); //0.f, MEMBUF1, k1_new - 1);	
			}
			Kbandcsum += Kband;
			X         += Kband * Npad;
		}
		Ksegcsum += Kseg;
		Y        += Kseg * Npad;
		XtY      += Kseg *Knew;

		X -= Kbandcsum * Npad;
	}
 
}
void XtX_ByGroup(BEAST2_BASESEG* SEG, I32 numSeg,F32PTR X,F32PTR XtX,I32 N, I32 Knew) {
	/*
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
	NEW.Knewterm, NEW.Knewterm, Nseg, 1.0,
	Xnewterm + r1 - 1, Npad,
	Xnewterm + r1 - 1, Npad, 0.f,
	MODEL.prop.XtX + (NEW.k1 - 1) * Knew + NEW.k1 - 1, Knew);// MEMBUF2, K_newTerm);
	*/
	I32 Npad = (N + 7) / 8 * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 Ksegcolcsum = 0;
	for (int i = 1; i <= numSeg; i++) {
		I32 Ksegcol =SEG[i- 1].K;
		I32 col_r1  =SEG[i - 1].R1;
		I32 col_r2  =SEG[i - 1].R2;
		I32 Ksegrowcsum = 0;
		for (int j = 1; j <= i; j++) {
			I32 Ksegrow = SEG[j - 1].K;
			I32 row_r1  = SEG[j - 1].R1;
			I32 row_r2  = SEG[j - 1].R2;

			I32 r1 = max(col_r1, row_r1);
			I32 r2 = max(col_r2, row_r2);
			I32 Nseg = r2 - r1 + 1;
			if (r2 >= r1) {
				r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
					Ksegrow, Ksegcol, Nseg, 1.0,
					X + Ksegrowcsum * Npad + r1 - 1, Npad,
					X + Ksegcolcsum * Npad + r1 - 1, Npad, 0.,
					XtX + Ksegrowcsum, Knew);
			}
			Ksegrowcsum += Ksegrow;			
		}

		Ksegcolcsum += Ksegcol;
		XtX += Ksegcol*Knew;
	}
	
}
void MoveCOLsWithinMatrix(F32PTR X, I32 N, I32 Kstart, I32 Kend, I32 Knewstart) {

	rI32 j = Knewstart - Kstart;
	if (j < 0 || Knewstart > Kend)
		// dst(k2_new):-----123455----
		// src(k2_old):----------123455----

		// dst(k2_new):----------------123455----
		// src(k2_old):-----123455----
		r_cblas_scopy((Kend-Kstart+1)*N, X+(Kstart-1)*N, 1, X+(Knewstart-1)*N, 1);
	else
	{
		// dst(k2_new):-------------123456----
		// src(k2_old):---------123456----

		rI32 segStartIdx = Kend+1;
		while (true) {
			segStartIdx = segStartIdx - j;
			if (segStartIdx > Kstart) {
				r_cblas_scopy(j * N, X + (segStartIdx-1) * N, 1L, X + ((segStartIdx+j)- 1) * N, 1);
			}
			else {
				j = (segStartIdx+j) - Kstart;
				r_cblas_scopy(j *N, X + (Kstart - 1)*N, 1L, X+(Knewstart-1) * N, 1);
				break;
			}
		}//while (true)

	}//if (j < 0 || k2_new + 1 > Kold)
}
 
//https: //stackoverflow.com/questions/3174850/what-is-the-correct-type-for-array-indexes-in-c#
//https: //stackoverflow.com/questions/797318/how-to-split-a-string-literal-across-multiple-lines-in-c-objective-c
 

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MR_EvaluateModel(
	BEAST2_MODELDATA *curmodel,  BEAST2_BASIS_PTR b, F32PTR Xt_mars, I32 N, I32 NUMBASIS,
	BEAST2_YINFO_PTR yInfo,	 BEAST2_HyperPar *hyperPar, F32PTR precVec, VOID_PTR stream )
{
	//TO RUN THIS FUNCTION, MUST FIRST RUN CONVERTBAIS so to GET NUMBERS OF BASIS TERMS	

	
	I32 Npad	= (I32)ceil((F32)N / 8.0f) * 8; Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 K		= 0;	 
	for (I32 basisID = 0; basisID < NUMBASIS; basisID++) {
		BEAST2_BASIS_PTR basis = b + basisID;
		if (basis->type != OUTLIERID) {
			int         NUM_SEG = basis->nKnot + 1;
			TKNOT_PTR   KNOT	= basis->KNOT;
			TORDER_PTR  ORDER   = basis->ORDER;

			BEAST2_BASESEG seg;
			seg.ORDER1 = basis->type == TRENDID ? 0 : 1;
			//For season terms: 1..(1)....(2)...(|sSegNum-1)...N (|N+1)
			for (int i = 1; i <= NUM_SEG; i++) {
				///The first segment starts from 1; other segments start one interval after the previous segment
				//The last segment ends at N'; others end exactly at the current breakpoint minus one;
				//the last bk is (N+1)
				seg.R1     = KNOT[(i - 1) - 1L];
				seg.R2     = KNOT[i - 1L] - 1L;
				seg.ORDER2 = basis->type==DUMMYID? 0 :  ORDER[i - 1L];
				I32 k  = basis->GenTerms(Xt_mars + Npad * K, N, &seg, &(basis->bConst));
				K += k;
			}
		} 	else	{
			int         numOfSeg = basis->nKnot;
			TKNOT_PTR   knotList = basis->KNOT;

			BEAST2_BASESEG seg;
			seg.ORDER1 = seg.ORDER2 = 0; // orders are not used at all
			//For season terms: 1..(1)....(2)...(|sSegNum-1)...N (|N+1)
			for (int i = 1; i <= numOfSeg; i++) {
				///The first segment starts from 1; other segments start one interval after the previous segment
				//The last segment ends at N'; others end exactly at the current breakpoint minus one;
				//the last bk is (N+1)
				seg.R1  = knotList[(i)-1L];
				seg.R2  = knotList[(i)-1L];
				I32 k = basis->GenTerms(Xt_mars + Npad * K, N, &seg, &(basis->bConst));
				K += k;
			}
		}


	}
	curmodel->K = K;

	//--------------------------------------------------------------------------------------------------
	// Set those rows of X_mars specfied by rowsMissing  to zeros 
	//	TODO: this could be buggy: Xt_mars may be not large enough to backup the X values at the missing rows
	// for the full initial model
	F32PTR	GlobalMEMBuf  = Xt_mars + K * Npad;
	F32PTR	Xt_zeroBackup = GlobalMEMBuf;
	if (yInfo->nMissing > 0) {
		F32 fillvalue = 0.f;
		f32_mat_multirows_extract_set_by_scalar(Xt_mars, Npad, K, Xt_zeroBackup, yInfo->rowsMissing, yInfo->nMissing, fillvalue);
	}

	//Calcuate X'*X and Calcuate X'*Y
	//DGEMM('T', 'N', &m, &n, &K, &alpha, X_mars, &K, X_mars, &K, &beta, XtX, &m);	
	//cblas_dgemm(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);

	
	F32PTR XtX = curmodel->XtX;
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.f, Xt_mars, Npad, Xt_mars, Npad, 0.f, XtX, K);
	//cblas_dgemmt (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
	//cblas_dgemmt(CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, k, N, 1, X_mars, N, X_mars, N, 0, XtY, k);				

    //MRBEAST
    I32 q = yInfo->q;

	F32PTR XtY = curmodel->XtY;
	//cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);
	//BEAST: r_cblas_sgemv(CblasColMajor, CblasTrans, Npad, K, 1, Xt_mars, Npad, pyInfo->Y, 1, 0, XtY, 1);	
	
	//MRBEAST--------start
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, q, N, 1.f, Xt_mars, Npad, yInfo->Y, N, 0.f, XtY, K);
	//MRBEAST---------end

	// Restore the backuped zeros at the missing rows
	if (yInfo->nMissing > 0) {
		f32_mat_multirows_set_by_submat(Xt_mars, Npad, K, Xt_zeroBackup, yInfo->rowsMissing, yInfo->nMissing);
	}
	// Now GlobalMEMBuf_1st is free to use;

	//X_mars has been constructed. Now use it to calcuate its margin likelihood
	//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec)	
	F32PTR cholXtX = curmodel->cholXtX;
	//f32_copy( XtX,  cholXtX, K * K);
	//f32_add_val_matrixdiag(cholXtX, precVal, K);

	//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
	F32PTR beta_mean = curmodel->beta_mean;	
	///////////////////////////////////////////////////////////////
	/*	
		//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
		r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, cholXtX, K); // Choleskey decomposition; only the upper triagnle elements are used
		//LAPACKE_spotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
		r_cblas_scopy(K, XtY, 1, beta_mean, 1);
		r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1, cholXtX, K, beta_mean, K);	
	*/
	
	chol_addCol_skipleadingzeros_prec_invdiag(XtX, cholXtX, precVec, K, 1, K);
	//BEAST: solve_U_as_LU_invdiag_sqrmat(cholXtX, XtY, beta_mean, K);
	//MRBEAST--------start
	solve_U_as_LU_invdiag_sqrmat_multicols(cholXtX, XtY, beta_mean, K, q);
   	//MRBEAST--------end


	//////////////////////////////////////////////////////////////////


	//Compute beta = beta_mean + Rsig2 * randn(p, 1);
	//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
	//status = vdRngGaussian( method, stream, n, r, a, sigma );

	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/
	/*
	F32PTR beta = model->beta;
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, (*(VSLStreamStatePtr *)stream), K, beta, 0, 1);

	// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
	{
		//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, post_P_U, K, beta, K);
		solve_U_as_U_invdiag(cholXtX, beta, K, K);
	}
	r_ippsMulC_32f_I(sqrtf(model->sig2), beta, K);
	r_ippsAdd_32f_I(beta_mean, beta, K);
	*/
	/**********************************************************/
	// Sample beta from beta_mean and cholXtX
	/********************************************************/


    //BEAST: F32 alpha2_star = (pyInfo->YtY_plus_alpha2Q[0] - DOT(K,  XtY,  beta_mean)) * 0.5;
	
	//MRBEAST--------start
    r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, K, 1.f, beta_mean, K, XtY, K, 0.f, curmodel->alpha2Q_star, q);
	r_ippsSub_32f(curmodel->alpha2Q_star, yInfo->YtY_plus_alpha2Q, curmodel->alpha2Q_star, q * q);
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', q, curmodel->alpha2Q_star, q); // Choleskey decomposition; only the upper triagnle elements are used

	F32 log_det_alphaQ = sum_log_diagv2(curmodel->alpha2Q_star, q);
	//MRBEAST--------end

	//half_log_det_post; = sum(log(1. / diag(U)))
	F32 half_log_det_post = sum_log_diagv2(cholXtX, K);
	
	//half_log_det_prior = -0.5 * (k_SN * log(prec(1)) + length(k_const_terms) * log(prec(2)) + length(linear_terms) * log(prec(3)));		
	F32 half_log_det_prior	= -.5f * K*logf(precVec[0]);

 	//log_ML = half_log_det_post - half_log_det_prior - (n / 2 + b) * log(a + a_star / 2);
	F32 marg_lik = q*(half_log_det_post - half_log_det_prior) - yInfo->alpha1_star * log_det_alphaQ*2.0f;

    curmodel->marg_lik    = marg_lik;
	//r_printf("Eval: det_post %f\n", f32_abs_sum(beta_mean,K*q));
	//r_printf("Lik: det_post %f\n", f32_sum_matrixdiag(cholXtX, K));
 	return;

	//MRBEAST--
 
}
#include "abc_000_warning.h"
