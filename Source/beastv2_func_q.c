#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include "abc_001_config.h"
#include "abc_rand.h"
#include "abc_mat.h"
#include "abc_math.h"  //sum_log_diag_v2
#include "abc_vec.h"   //sum_log_diag_v2
#include "beastv2_func.h"
#include "beastv2_prior_precfunc.h"


#include <math.h>
#include <string.h>
#include <stdio.h>	 

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

int BEAST2_Basis_To_XmarsXtX_XtY(BEAST2_BASIS_PTR b, I32 NUMBASIS, F32PTR Xt_mars, I32 N, F32PTR XtX, F32PTR XtY, BEAST2_YINFO_PTR  yInfo) {
	
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
	//curmodel->K = K;

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

	
	//F32PTR XtX = curmodel->XtX;
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.f, Xt_mars, Npad, Xt_mars, Npad, 0.f, XtX, K);
	//cblas_dgemmt (const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const CBLAS_TRANSPOSE transa, const CBLAS_TRANSPOSE transb, const MKL_INT n, const MKL_INT k, const double alpha, const double *a, const MKL_INT lda, const double *b, const MKL_INT ldb, const double beta, double *c, const MKL_INT ldc);
	//cblas_dgemmt(CblasColMajor, CblasLower, CblasTrans, CblasNoTrans, k, N, 1, X_mars, N, X_mars, N, 0, XtY, k);	


	//F32PTR XtY = curmodel->XtY;
	if (yInfo->q == 1) {	
		//cblas_dgemv(const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE trans, const MKL_INT m, const MKL_INT n, const double alpha, const double *a, const MKL_INT lda, const double *x, const MKL_INT incx, const double beta, double *y, const MKL_INT incy);
		r_cblas_sgemv(CblasColMajor, CblasTrans, Npad, K, 1, Xt_mars, Npad, yInfo->Y, 1, 0, XtY, 1);
	}	else {
		//MRBEAST--------start
		r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, yInfo->q, N, 1.f, Xt_mars, Npad, yInfo->Y, N, 0.f, XtY, K);
		//MRBEAST---------end
	}

	// Restore the backuped zeros at the missing rows
	if (yInfo->nMissing > 0) {
		f32_mat_multirows_set_by_submat(Xt_mars, Npad, K, Xt_zeroBackup, yInfo->rowsMissing, yInfo->nMissing);
	}
	// Now GlobalMEMBuf_1st is free to use;

	return K;
}

void BEAST2_EvaluateModel(	BEAST2_MODELDATA *curmodel, BEAST2_BASIS_PTR b, F32PTR Xt_mars, I32 N, I32 NUMBASIS,
	      BEAST2_YINFO_PTR  yInfo,    BEAST2_HyperPar *hyperPar, PRECSTATE_PTR precState, PREC_FUNCS * precFunc )
{
	//TO RUN THIS FUNCTION, MUST FIRST RUN CONVERTBAIS so to GET NUMBERS OF BASIS TERMS	
	
	// urmodel->XtX is filled in the function below
	curmodel->K = BEAST2_Basis_To_XmarsXtX_XtY(b, NUMBASIS, Xt_mars, N, curmodel->XtX, curmodel->XtY, yInfo);

	I32 Npad    = N;
	I32 K       = curmodel->K;
	F32PTR XtX  = curmodel->XtX;
	F32PTR XtY  = curmodel->XtY;

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
	
   
	/////////////////////////////////////////////////////////
	//                Get PrecXttDiag
	///////////////////////////////////////////////////////// 
	// [0:ConstPrec, 1:UniformPrec, 2:ComponentWise], 3:OrderWise		 
	precFunc->SetPrecXtXDiag(curmodel->precXtXDiag, b, NUMBASIS, precState);    //CHANGE: (1)nothing  or (2) MODEL.curr.precXtXDiag
	precFunc->chol_addCol(XtX, cholXtX, curmodel->precXtXDiag, K, 1, K);


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
	 

	precFunc->SetNtermsPerPrecGrp(curmodel->nTermsPerPrecGrp, b, NUMBASIS, precState);
	precFunc->ComputeMargLik(curmodel, precState, yInfo, hyperPar);

	return;
}

  
/////////////////////////////////////////////////////////////////////////////////////////////////


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

		if (b->type != OUTLIERID) {
 
			I32 numSeg = b->nKnot + 1;
			for (int j = 0; j < numSeg; j++) {
				I32 Kbase = b->Kbase;
				if (Klastcol >= (Kbase + b->ks[j])) {
					info->R1 = b->KNOT[j - 1];
					info->R2 = b->KNOT[j] - 1L;
					info->K  = min(Kbase + b->ke[j], Klastcol) - (Kbase + b->ks[j]) + 1;
					info++;
					numBands++;
				} else {
					QUITFLAG = 1;	
					break;
				}
			}

		} 
		else {
			I32 numSeg = b->nKnot;
			for (int j = 0; j < numSeg; j++) {
				I32 Kbase = b->Kbase;
				if (Klastcol >= (Kbase + b->ks[j])) {
					info->R1 = b->KNOT[j];
					info->R2 = b->KNOT[j] ;
					info->K  = min(Kbase + b->ke[j], Klastcol) - (Kbase + b->ks[j]) + 1;
					info++;
					numBands++;
				}	else {
					QUITFLAG = 1;	
					break;
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
		if (b->type != OUTLIERID) {

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
		else {
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
 

void Update_XtX_from_Xnewterm_ByGroup(F32PTR X, F32PTR Xnewterm, F32PTR XtX, F32PTR XtXnew, NEWTERM* NEW , BEAST2_MODEL_PTR model) {

    const I32 k1       = NEW->newcols.k1;
	const I32 k2_old   = NEW->newcols.k2_old;
	const I32 k2_new   = NEW->newcols.k2_new;
	const I32 N        = NEW->newcols.N;
	const I32 Nlda     = NEW->newcols.Nlda;
	const I32 Knewterm = NEW->newcols.Knewterm;
	const I32 KOLD     = NEW->newcols.KOLD;
	const I32 KNEW     = NEW->newcols.KNEW;

	const I32  Npad    = N;//Correct for the inconsitency of X and Y in gemm and gemv
	/*************************************************************************/
	//               The FIRST component:	
	/*************************************************************************/
	// There'sY no first component if k1_old/k1_new=1 for SEASON	
	for (I32 i = 1; i < k1; i++) SCPY(i, XtX + (i - 1L) * KOLD, XtXnew + (i - 1L) * KNEW);

	/*************************************************************************/
	//              The SECOND component
	/*************************************************************************/
	// No new cols/terms if flag=ChORDER && isInsert=0:the resampled basis has a higher order than the old basis
	if (Knewterm != 0) {

		FILL0(XtXnew + (k1 - 1) * KNEW, (KNEW - k1 + 1) * KNEW); // zero out the cols from k1-th to the end
		if (k1 > 1) {
			/*
			I32 r1  = NEW.r1[0];
			I32 r2  = NEW.r2[NEW.numSeg - 1];
			I32 Nseg = r2 - r1 + 1;
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
				NEW.k1 - 1L, Knewterm, Nseg, 1.0f, \
				Xt_mars  + r1-1L, Npad,
				Xnewterm + r1-1L, Npad, 0.f, \
				MODEL.prop.XtX + (NEW.k1 - 1L) * KNEW, KNEW); //0.f, MEMBUF1, k1_new - 1);
			*/

			// Xnewterm is pre-allocated with sufficent mem to ensure segInfo won't overflow in __GetMAXNumElemXnewTerm
			BEAST2_BASESEG*  OLD_SEG    = (BEAST2_BASESEG*)(Xnewterm + Knewterm * Npad);
			I32              OLD_numSeg = GetInfoBandList(OLD_SEG, model, k1 - 1);
			MatxMat( OLD_SEG,  OLD_numSeg,  X,
				     NEW->SEG, NEW->numSeg, Xnewterm,
				     XtXnew + (k1 - 1L) * KNEW, N, KNEW);
		}

		/*
		  r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
		   Knewterm, Knewterm, Nseg, 1.0,
		   Xnewterm + r1 - 1, Npad,
		   Xnewterm + r1 - 1, Npad, 0.f,
		   MODEL.prop.XtX + (NEW.k1-1) * KNEW + NEW.k1 - 1, KNEW);// MEMBUF2, K_newTerm);
		*/
		//XnewtermTXnewterm(&NEW, Xnewterm, MODEL.prop.XtX + (NEW.k1 - 1) * KNEW + NEW.k1 - 1, Npad, KNEW);
		XtX_ByGroup(NEW->SEG, NEW->numSeg, Xnewterm, XtXnew + (k1 - 1) * KNEW + k1 - 1, N, KNEW);

		/* //After obtaining Xnewterm'*Xnewterm, insert it into XtX_prop at appropriate locations
			for (rI32 i = k1_new, j = 1; i <= k2_new; i++, j++) {
				if (k1_new != 1) r_cblas_scopy(k1_new - 1, MEMBUF1 + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
				r_cblas_scopy(j, MEMBUF2 + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);}
		*/
	}
	/*{
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1_new - 1, K_newTerm, N, 1, X_mars, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, k1_new - 1);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
		 r_cblas_scopy(k1_new - 1, GlobalMEMBuf_1st + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, K_newTerm, N, 1, X_mars_prop + (k1_new - 1)*N, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, K_newTerm);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
		 r_cblas_scopy(j, GlobalMEMBuf_1st + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);
	}*/

	/*************************************************************************/
	//                  The THRID component: 
	/*************************************************************************/
	//There is no third componet if k2_old=KOLD 
	if (k2_old != KOLD) {
		/*for (rI32  j = 1; i <= KOLD; i++, j++) {r_cblas_scopy(K_newTerm,  MEMBUF1 + (j - 1)*K_newTerm, 1, XtX_prop + (k - 1)*KNEW + k1_new - 1, 1),					*/
		for (I32 kold = k2_old + 1, knew = k2_new + 1; kold <= KOLD; kold++, knew++) {
			F32PTR ColStart_old = XtX    + (kold - 1) * KOLD;
			F32PTR ColStart_new = XtXnew + (knew - 1) * KNEW;
			SCPY(k1 - 1,        ColStart_old,                    ColStart_new                   ); //the upper part of the third componet
			SCPY(kold - k2_old, ColStart_old + (k2_old + 1) - 1, ColStart_new + (k2_new + 1) - 1); // the bottom part of the 3rd cmpnt
		}

		// If there is a MIDDLE part of the componet (i.e, Knewterm>0); this part 
		// will be missing if flag is resmaplingOder and isInsert = 0.
		if (Knewterm != 0) {
			/*
			  rI32 r1 = NEW.r1[0];
			  rI32 r2 = NEW.r2[NEW.numSeg - 1];
			  rI32 Nseg = r2 - r1 + 1;
			  r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
				  Knewterm, (KOLD - NEW.k2_old), Nseg, 1,
				  Xnewterm + r1 - 1, Npad,
				  Xt_mars + (NEW.k2_old + 1 - 1) * Npad + r1 - 1, Npad, 0,
				  MODEL.prop.XtX + (NEW.k2_new+1L-1L) * KNEW + NEW.k1 - 1, KNEW );//MEMBUF1, K_newTerm);
			 */
			BEAST2_BASESEG* OLD_SEG    = (BEAST2_BASESEG*)(Xnewterm + Knewterm * Npad);
			I32             OLD_numSeg = GetInfoBandList_post(OLD_SEG, model, k2_old + 1);
			MatxMat(NEW->SEG, NEW->numSeg, Xnewterm,
				    OLD_SEG, OLD_numSeg, X + k2_old * Npad,
				    XtXnew + (k2_new + 1 - 1) * KNEW + k1 - 1, N, KNEW);
		}

	}

}

void Update_XtY_from_Xnewterm_ByGroup(F32PTR Y, F32PTR Xnewterm, F32PTR XtY, F32PTR XtYnew, NEWTERM* NEW , int q) {

    const I32 k1       = NEW->newcols.k1;
	const I32 k2_old   = NEW->newcols.k2_old;
	const I32 k2_new   = NEW->newcols.k2_new;
	const I32 N        = NEW->newcols.N;
	const I32 Nlda     = NEW->newcols.Nlda;
	const I32 Knewterm = NEW->newcols.Knewterm;
	const I32 KOLD     = NEW->newcols.KOLD;
	const I32 KNEW     = NEW->newcols.KNEW;

	const I32  Npad    = N;//Correct for the inconsitency of X and Y in gemm and gemv
	 
	/*********************************************************************************/
	//                 Compute XtY_prop from XtY
	/********************************************************************************/
	if (q == 1) {

		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1)          SCPY(k1 - 1,XtY, XtYnew);
		// New components : XnewTemrm*Y
		if (Knewterm > 0)  	 MatxVec(NEW->SEG, NEW->numSeg, Xnewterm, Y,  XtYnew + k1 - 1, N);
		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD)  SCPY(KNEW - k2_new, XtY + (k2_old + 1L) - 1L, XtYnew + (k2_new + 1) - 1);

	} else {
		// FOR MrBEAST

		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(k1 - 1,XtY + KOLD * c, XtYnew + KNEW * c);
			}
		}
		// New components : XnewTemrm*Y
		if (Knewterm > 0) {
			//MatxVec(NEW.SEG, NEW.numSeg, Xnewterm, yInfo.Y, MODEL.prop.XtY + NEW.k1 - 1, N);
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, Knewterm, q, N, 1.f, Xnewterm, Npad, Y, N, 0.f, XtYnew + k1 - 1, KNEW);
		}
		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(KNEW -k2_new, XtY + (k2_old + 1L) - 1L + KOLD * c, XtYnew + (k2_new + 1) - 1 + KNEW * c);
			}
		}

	} //if (q == 1) 

}

void Update_XtX_from_Xnewterm_NoGroup(F32PTR X, F32PTR Xnewterm, F32PTR XtX, F32PTR XtXnew, NEWTERM* NEW, BEAST2_MODEL* MODEL_not_used) {

    const I32 k1       = NEW->newcols.k1;
	const I32 k2_old   = NEW->newcols.k2_old;
	const I32 k2_new   = NEW->newcols.k2_new;
	const I32 N        = NEW->newcols.N;
	const I32 Nlda     = NEW->newcols.Nlda;
	const I32 Knewterm = NEW->newcols.Knewterm; // // k2_new - k1 + 1L;
	const I32 KOLD     = NEW->newcols.KOLD;
	const I32 KNEW     = NEW->newcols.KNEW;
 
	/*************************************************************************/
	//               The FIRST component:	
	/*************************************************************************/
	// There'sY no first component if k1_old/k1_new=1 for SEASON	
	for (I32 i = 1; i < k1; i++) SCPY(i, XtX + (i - 1L) * KOLD, XtXnew + (i - 1L) * KNEW);

	/*************************************************************************/
	//              The SECOND component
	/*************************************************************************/
	// No new cols/terms if flag=ChORDER && isInsert=0:the resampled basis has a higher order than the old basis
	if (Knewterm != 0) {

		FILL0(XtXnew + (k1 - 1) * KNEW, (KNEW - k1 + 1) * KNEW); // zero out the cols from k1-th to the end
		if (k1 > 1) {
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1 - 1, Knewterm, N, 1.0f,
				X, Nlda,
				Xnewterm, Nlda, 0.f,
				XtXnew + (k1 - 1L) * KNEW, KNEW);
		}


		// Three alternative ways to compute Xnewterm'*XnewTerm. Note that the resulting matrix is symmetric

		// THE FIRST WAY:
		r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
			Knewterm, Knewterm, N, 1.0,
			Xnewterm, Nlda,
			Xnewterm, Nlda, 0.f,
			XtXnew + (k1 - 1) * KNEW + k1 - 1, KNEW);


		//THE SECOND WAY: 
		//sgemmt only updates the upper triangular part of the resulting matrix, which is supposed to be faster than sgemm, but it is NOT
		//cblas_sgemmt(CblasColMajor, CblasUpper, CblasTrans, CblasNoTrans, K_newTerm, Npad, 1.0f, Xnewterm, Npad, Xnewterm, Npad, 0.f, GlobalMEMBuf_2nd, K_newTerm);

		//THE THIRD WAY: 
		//This is the fastest way when using Intel'sY MKL
		/*
		{  for (int i = 1; i <= K_newTerm; i++)
		   for (int j = 1; j <= i; j++)
		   GlobalMEMBuf_2nd[K_newTerm*(i - 1) + j - 1] = DOT(N, Xnewterm + (j - 1)*Npad, Xnewterm + (i - 1)*Npad);
		} */

		/*
		//After obtaining Xnewterm'*Xnewterm, insert it into XtX_prop at appropriate locations
		for (rI32 i = k1_new, j = 1; i <= k2_new; i++, j++) {
			if (k1_new != 1) r_cblas_scopy(k1_new - 1, MEMBUF1 + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
			r_cblas_scopy(j, MEMBUF2 + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);
		}*/
	}
	/*{
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1_new - 1, K_newTerm, N, 1, X_mars, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, k1_new - 1);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
	r_cblas_scopy(k1_new - 1, GlobalMEMBuf_1st + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, K_newTerm, N, 1, X_mars_prop + (k1_new - 1)*N, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, K_newTerm);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
	r_cblas_scopy(j, GlobalMEMBuf_1st + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);
	}*/

	/*************************************************************************/
	//                  The THRID component: 
	/*************************************************************************/
	//There is no third componet if k2_old=KOLD 
	if (k2_old != KOLD) {
		/*for (rI32  j = 1; i <= KOLD; i++, j++) {r_cblas_scopy(K_newTerm,  MEMBUF1 + (j - 1)*K_newTerm, 1, XtX_prop + (k - 1)*KNEW + k1_new - 1, 1),					*/
		for (I32 kold = k2_old + 1, knew = k2_new + 1; kold <= KOLD; kold++, knew++) {
			F32PTR ColStart_old = XtX + (kold - 1) * KOLD;
			F32PTR ColStart_new = XtXnew + (knew - 1) * KNEW;
			SCPY(k1 - 1,       ColStart_old, ColStart_new); //the upper part of the third componet
			SCPY(kold - k2_old, ColStart_old + (k2_old + 1) - 1, ColStart_new + (k2_new + 1) - 1); // the bottom part of the 3rd cmpnt
		}

		// If there is a MIDDLE part of the componet (i.e, Knewterm>0); this part 
		// will be missing if flag is resmaplingOder and isInsert = 0.
		if (Knewterm != 0) {
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
				Knewterm, (KOLD - k2_old), N, 1.0,
				Xnewterm, Nlda,
				X + (k2_old + 1 - 1) * Nlda, Nlda, 0.0,
				XtXnew + (k2_new + 1 - 1) * KNEW + k1 - 1, KNEW);
		}

	}

			 
}

void Update_XtY_from_Xnewterm_NoGroup(F32PTR Y, F32PTR Xnewterm, F32PTR XtY, F32PTR XtYnew, NEWTERM* NEW, I32 q) {

	// X and Xnewterm has a leading dimesnion of new.Nlada
	// Y has a leading dimension of new.N

    const I32 k1       = NEW->newcols.k1;
	const I32 k2_old   = NEW->newcols.k2_old;
	const I32 k2_new   = NEW->newcols.k2_new;
	const I32 N        = NEW->newcols.N;
	const I32 Nlda     = NEW->newcols.Nlda;
	const I32 Knewterm = NEW->newcols.Knewterm; // // k2_new - k1 + 1L;
	const I32 KOLD     = NEW->newcols.KOLD;
	const I32 KNEW     = NEW->newcols.KNEW;

/*********************************************************************************/
//                 Compute XtY_prop from XtY
/********************************************************************************/
	if (q == 1) {
		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1)       SCPY(k1 - 1, XtY,XtYnew);
		// New components : XnewTemrm*Y
		if (Knewterm > 0) { 
				r_cblas_sgemv(CblasColMajor, CblasTrans, N, Knewterm, 1.f,
						Xnewterm,   Nlda,
						Y,         1L, 0.f,
					    XtYnew + k1 - 1, 1L);
		}
		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD) SCPY(KNEW - k2_new, XtY + (k2_old + 1L) - 1L, XtYnew + (k2_new + 1) - 1);

	}
	else {
		// FOR MrBEAST

		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(k1 - 1, XtY + KOLD * c, XtYnew + KNEW * c);
			}
		}
		// New components : XnewTemrm*Y
		if (Knewterm > 0) {
			//MatxVec(NEW.SEG, NEW.numSeg, Xnewterm, yInfo.Y, MODEL.prop.XtY + NEW.k1 - 1, N);
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
				Knewterm, q, N, 1.f, 
				Xnewterm, Nlda,
				Y,        N, 0.f,
				XtYnew + k1 -1, KNEW);
		}

		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(KNEW - k2_new, XtY + (k2_old + 1L) - 1L + KOLD * c, XtYnew + (k2_new + 1) - 1 + KNEW * c);
			}
		}


	}
}
//https: //stackoverflow.com/questions/3174850/what-is-the-correct-type-for-array-indexes-in-c#
//https: //stackoverflow.com/questions/797318/how-to-split-a-string-literal-across-multiple-lines-in-c-objective-c
#include "abc_000_warning.h"
