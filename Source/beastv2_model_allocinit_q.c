#include <string.h>  //memset
#include <math.h>    //sqrt

#include "abc_000_warning.h"

#include "abc_vec.h"
#include "abc_ts_func.h"
#include "abc_mem.h"
#include "beastv2_header.h"
#include "beastv2_model_allocinit.h"
#include "beastv2_prior_precfunc.h" 

extern void* Get_Propose(I08, BEAST2_OPTIONS_PTR);
extern pfnGenTerms * Get_GenTerms(I08 id, BEAST2_PRIOR_PTR prior);
extern void* Get_CalcBasisKsKeK(I08 id, I08 precPriorType);
extern void* Get_UpdateGoodVec(I08);
extern void* Get_ComputeY(I08, BEAST2_OPTIONS_PTR);
extern void* Get_ModelPrior(I08);
extern void* Get_GenRandomBasis(I08);
extern void* Get_AllocInitBasis(I08);
extern void* Get_PickBasisID(I08, I08, I32PTR);
extern void* Get_CvtKnotsToBinVec(I08 id);

extern void PreCaclModelNumber(I32 minOrder, I32 maxOrder, I32 maxNumseg, I32 N, I32 minSep, F64PTR TNUM, F64PTR totalNum);


#define MODEL (*model)

void AllocInitModelMEM(BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt, MemPointers* MEM)
{	
	// not needed bcz it is initialized to 0 when model is declared
	// memset(model, 0, sizeof(BEAST2_MODEL));
	I32 N      = opt->io.N;
	I32 K_MAX  = opt->prior.K_MAX;
	I32 q      = opt->io.q;

	// sig2 is intialized only once here. This default value is used only for the the first chain
	// of the first pixel. All other chains or pixels will be using the existing value left from 
	// previous runs. This may have some unintended consequences
	if (q == 1) {
		// BEASTV2
		MODEL.sig2 = opt->prior.sig2; // sig2 is needed only for generating beta,  
	} else { 
		//for MRBEAST		
		// Allocate 3 arrays altogether:  SIG2(q2)+curr.alpah2Q+prop.alpha2Q
		I32   qq = q * q;
		MODEL.SIG2               = MyALLOC0(*MEM, 2*qq + qq + qq, F32, 64);
		MODEL.curr.alphaQ_star   = MODEL.SIG2+ 2 * qq;
		MODEL.prop.alphaQ_star   = MODEL.curr.alphaQ_star + qq;

		// SIG2 has a lenght of 2 * (q*q); the first q*q elements are used to save SIG2_Upper;
		// The second q*q elements are for inv(SIG2_Upper), which is needed to compute  Beta*inv(SIG2)*Beta' when 
		// resampling Precision
		
		// VERY IMPORTANT--alphaQ_star needs to have lower triangles zeroed out when it used to resample SIG2:  WLower*SIG2_upper=alphaQ_star_chol

		
		// FIll the diagonal of SIG2 with an inital value
		f32_fill_val_matrixdiag(MODEL.SIG2, opt->prior.sig2, q);
	}
	

	MODEL.beta			  = MyALLOC(*MEM,  K_MAX * q,      F32, 64);

	MODEL.curr.XtX        = MyALLOC(*MEM, K_MAX * K_MAX,   F32, 64);
	MODEL.curr.XtY        = MyALLOC(*MEM, K_MAX*q,         F32, 64);
	MODEL.curr.cholXtX    = MyALLOC(*MEM, K_MAX * K_MAX,   F32, 64);
	MODEL.curr.beta_mean  = MyALLOC(*MEM, K_MAX * q,       F32, 64);
	
	MODEL.prop.XtX        = MyALLOC(*MEM, K_MAX * K_MAX,   F32, 64);
	MODEL.prop.XtY        = MyALLOC(*MEM, K_MAX * q,       F32, 64);
	MODEL.prop.cholXtX    = MyALLOC(*MEM, K_MAX * K_MAX,   F32, 64);
	MODEL.prop.beta_mean  = MyALLOC(*MEM, K_MAX * q,       F32, 64);

	// Allocate 3 arrays altogether: Deviation (N*q) + avgDeivation (q) + extremePosVec (Npad16)
	MODEL.deviation       = MyALLOC(*MEM, (N*q) + (q),   F32, 64);	  // changed for MRBEAST
	MODEL.avgDeviation    = MODEL.deviation    + N*q;
	{
	//	MODEL.extremePosVec needs to be 8-byte aligned; otherwise,there is an 
	// run-time error in Debian Linux, R-devel, GCC ASAN/UBSAN.
	I32 Npad16 = (N + 15) / 16 * 16;
	MODEL.extremePosVec   = MyALLOC(*MEM, Npad16,  I08,   8);
	// MODEL.extremPosNum    = 0; // No need to initialize bcz it will be reset to yInfo.n in beast_core
	}
	
	

	MODEL.NUMBASIS = opt->prior.numBasis;
	I32   NumBasis = MODEL.NUMBASIS;
	for (int i = 0; i < NumBasis; i++)
		MODEL.b[i].type = opt->prior.basisType[i];

	I32 isComponentFixed[3] = {0,0,0};
	MODEL.did = MODEL.sid = MODEL.tid = MODEL.oid = -1;
	for (int i = 0; i < NumBasis; i++) 	{

		BEAST2_BASIS_PTR	basis	= MODEL.b + i;
		I08					type	= basis->type;

		void (*AllocInitBasisMEM)(BEAST2_BASIS_PTR, I32,I32, MemPointers*)= Get_AllocInitBasis(type);

		if      (type == TRENDID)
		{
			MODEL.tid = i;
			
			basis->prior.minKnotNum = opt->prior.trendMinKnotNum;
			basis->prior.maxKnotNum = opt->prior.trendMaxKnotNum;
			basis->prior.minSepDist = opt->prior.trendMinSepDist;
			basis->prior.minOrder   = opt->prior.trendMinOrder;
			basis->prior.maxOrder   = opt->prior.trendMaxOrder;
			isComponentFixed[i] = basis->prior.minOrder == basis->prior.maxOrder && basis->prior.minKnotNum == 0 && basis->prior.maxKnotNum == 0;

			AllocInitBasisMEM(basis, N, K_MAX, MEM);

			basis->mcmc_Kstopping = K_MAX - (basis->prior.maxOrder + 1);
			basis->mcmc_MoveStep  = opt->mcmc.maxMoveStepSize;

			/****************************************************************************/
			// Prepare and re-adjust the propoal probabilities of indiviudal MOVEs for the
			// trend and seasonal components
			/****************************************************************************/
            #define MAX_RAND_BYTE_VALUE 255 //127 
		   // If resamplingOrderProb is zero,  resampling the order will be NOT conducted; this
		   // happens if minOrder=maxOrder (checked in PostCheckArgs in beastv2_in_args  where 
		   // resamplingOderProb is set to 0 if needed)
			F32 cumProbExcludeResamplingOrder = 1 - opt->mcmc.trendResamplingOrderProb;
			basis->propprob = (PROP_PROB_STRUCT){
				  .birth = (U08)(cumProbExcludeResamplingOrder * 0.33 * MAX_RAND_BYTE_VALUE),
				  .move = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33) * MAX_RAND_BYTE_VALUE) ,
				  .death = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33 + 0.165) * MAX_RAND_BYTE_VALUE),
				  .merge = (U08)(cumProbExcludeResamplingOrder * 1. * MAX_RAND_BYTE_VALUE) };

			/*************************************************************************************************/
			// Pre-calcalte terms for SEASON AND TRENDS Bases:Precaluating SEASON AND TREND terms 
			/*************************************************************************************************/
			I32 MAX_ORDER = basis->prior.maxOrder;

			TREND_CONST *TREND = &basis->bConst;
			TREND->INV_SQR = MyALLOC(*MEM, N, F32, 0),  // sqrt(1/1), sqrt(1/2), ....,sqrt(1/N)
			TREND->COEFF_A = MyALLOC(*MEM, N, F32, 0),  // intercpet and slope for a linear segment of length ith: i=1, ,..N
 			TREND->COEFF_B = MyALLOC(*MEM, N, F32, 0),
			TREND->TERMS   = MyALLOC(*MEM, N*(MAX_ORDER + 1L), F32, 64);

			preCalc_terms_trend(TREND->TERMS, TREND->INV_SQR, N, MAX_ORDER);

			
			if (opt->prior.trendBasisFuncType==4)
				preCalc_XmarsTerms_extra_fmt4(TREND->COEFF_A, TREND->COEFF_B, N);
			else
				preCalc_XmarsTerms_extra(TREND->COEFF_A, TREND->COEFF_B, N);

		
			/****************************************************************************/
			//				Pre-compute the scaling factor
			/****************************************************************************/
			I32 tMAXNUMKNOT = basis->prior.maxKnotNum;
			I32 tMINSEPDIST = basis->prior.minSepDist;
			F32PTR scaleFactor = MyALLOC(*MEM, tMAXNUMKNOT + 1, F32, 0);
			{
				F32PTR MEMBUF1 = MODEL.curr.XtX;  //MODEL.curr.XtX used here a temp MEM buffer
				F32PTR MEMBUF2 = MODEL.curr.XtX + (tMAXNUMKNOT + 1);
				preCalc_scale_factor(scaleFactor, N, tMAXNUMKNOT, tMINSEPDIST, MEMBUF1, MEMBUF2);
			}
			basis->scalingFactor = scaleFactor;

			/////////////////////////////////////////////////////////////
			//TODO: the model prior adjustment factors are needed only for certain modelPriorType
			I32 NUMSEG    = opt->prior.trendMaxKnotNum + 1;
			I32 MINORDER  = opt->prior.trendMinOrder + 1;
			I32 MAXORDER  = opt->prior.trendMaxOrder + 1;

			if (opt->prior.modelPriorType ==4) {
				F64PTR priorVec = MyALLOC(*MEM, NUMSEG * MAXORDER, F64, 64);
				F64PTR priorMat = MyALLOC(*MEM, NUMSEG * MAXORDER * NUMSEG, F64, 64);

				PreCaclModelNumber(MINORDER, MAXORDER, NUMSEG, N, opt->prior.trendMinSepDist,	priorMat, priorVec);

				basis->priorMat = priorMat;
				basis->priorVec = priorVec;	
			}
			
		}
		else if (type == SEASONID)
		{
			MODEL.sid = i;

			basis->prior.minKnotNum = opt->prior.seasonMinKnotNum;
			basis->prior.maxKnotNum = opt->prior.seasonMaxKnotNum;
			basis->prior.minSepDist = opt->prior.seasonMinSepDist;
			basis->prior.minOrder   = opt->prior.seasonMinOrder;
			basis->prior.maxOrder   = opt->prior.seasonMaxOrder;
			isComponentFixed[i]     = basis->prior.minOrder == basis->prior.maxOrder && basis->prior.minKnotNum == 0 && basis->prior.maxKnotNum == 0;

			AllocInitBasisMEM(basis, N, K_MAX, MEM);

			basis->mcmc_Kstopping = K_MAX - (2 * basis->prior.maxOrder);
			basis->mcmc_MoveStep  = opt->mcmc.maxMoveStepSize;

			/****************************************************************************/
			// Prepare and re-adjust the propoal probabilities of indiviudal MOVEs for the
			// trend and seasonal components
			/****************************************************************************/
		   // If resamplingOrderProb is zero,  resampling the order will be NOT conducted; this
		   // happens if minOrder=maxOrder (checked in PostCheckArgs in beastv2_in_args  where 
		   // resamplingOderProb is set to 0 if needed)
			F32 cumProbExcludeResamplingOrder = 1 - opt->mcmc.seasonResamplingOrderProb;
			basis->propprob = (PROP_PROB_STRUCT){
					  .birth = (U08)(cumProbExcludeResamplingOrder * 0.33 * MAX_RAND_BYTE_VALUE),
					  .move = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33) * MAX_RAND_BYTE_VALUE) ,
					  .death = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33 + 0.165) * MAX_RAND_BYTE_VALUE),
					  .merge = (U08)(cumProbExcludeResamplingOrder * 1. * MAX_RAND_BYTE_VALUE) };

			/*************************************************************************************************/
			// Pre-calcalte terms for SEASON AND TRENDS Bases:Precaluating SEASON AND TREND terms 
			/*************************************************************************************************/
			// When deseaonalize=TRUE, the full set of cos/sin are needed to compute seasonality 
			I32 sMAXORDER = opt->io.meta.deseasonalize? opt->io.meta.period/2 : basis->prior.maxOrder;
			//SEASON_CONST  *SEASON; // Must be global;otherwise &SEASON is transient and will be lost.

			SEASON_CONST* SEASON = &basis->bConst;
			SEASON->TERMS        = MyALLOC(*MEM,  N * sMAXORDER * 2L, F32, 64),
			SEASON->SQR_CSUM     = MyALLOC(*MEM,  (N + 1L)*sMAXORDER * 2L, F32, 64) ;  
			//term:       1, 2,  3,  . ..,   N
			//sqr_cum: 0, 1, 12, 123,     , 1-N
			//so the sum of term_i_to_j can be computed as  sqr_csum(j)-sqrt_csum(i-1).

	
			preCalc_terms_season(SEASON->TERMS, SEASON->SQR_CSUM, N, opt->io.meta.period, sMAXORDER);	 

			/****************************************************************************/
			//				Pre-compute the scaling factor
			/****************************************************************************/
			I32 sMAXNUMKNOT = basis->prior.maxKnotNum;
			I32 sMINSEPDIST = basis->prior.minSepDist;
			F32PTR scaleFactor = MyALLOC(*MEM, sMAXNUMKNOT + 1, F32, 0);
			{
				F32PTR MEMBUF1 = MODEL.curr.XtX;  //MODEL.curr.XtX used here a temp MEM buffer
				F32PTR MEMBUF2 = MODEL.curr.XtX + (sMAXNUMKNOT + 1);
				preCalc_scale_factor(scaleFactor, N, sMAXNUMKNOT, sMINSEPDIST, MEMBUF1, MEMBUF2);
			}
			basis->scalingFactor = scaleFactor;

			/////////////////////////////////////////////////////////////
			// TODO: the model prior adjustment factors are needed only for certain modelPriorType
			I32 NUMSEG    = opt->prior.seasonMaxKnotNum + 1;
			I32 MINORDER  = opt->prior.seasonMinOrder;
			I32 MAXORDER  = opt->prior.seasonMaxOrder;

			if (opt->prior.modelPriorType == 4) {
				F64PTR priorVec = MyALLOC(*MEM, NUMSEG * MAXORDER, F64, 64);
				F64PTR priorMat = MyALLOC(*MEM, NUMSEG * MAXORDER * NUMSEG, F64, 64);
				PreCaclModelNumber(MINORDER, MAXORDER, NUMSEG, N, opt->prior.seasonMinSepDist, priorMat, priorVec);
				basis->priorMat = priorMat;
				basis->priorVec = priorVec;
			}

		}
		else if (type == DUMMYID)
		{
			MODEL.did = i;
			MODEL.sid = i;

			basis->prior.minKnotNum = opt->prior.seasonMinKnotNum;
			basis->prior.maxKnotNum = opt->prior.seasonMaxKnotNum;
			basis->prior.minSepDist = opt->prior.seasonMinSepDist;
			basis->prior.minOrder   = -1;//opt->prior.seasonMinOrder;
			basis->prior.maxOrder   = -1;//opt->prior.seasonMaxOrder;
			isComponentFixed[i]     = basis->prior.minOrder == basis->prior.maxOrder && basis->prior.minKnotNum == 0 && basis->prior.maxKnotNum == 0;

			AllocInitBasisMEM(basis, N, K_MAX, MEM);
		 	
			
			basis->mcmc_Kstopping = K_MAX - (2 * opt->io.meta.period);
			basis->mcmc_MoveStep  = opt->mcmc.maxMoveStepSize;

			/****************************************************************************/
			// Prepare and re-adjust the propoal probabilities of indiviudal MOVEs for the
			// trend and seasonal components
			/****************************************************************************/
			// seasonResamplingOrderProb has been set to zero in the PostCheckArgs
			F32 resamplingOrderProb = 1 - opt->mcmc.seasonResamplingOrderProb;
			basis->propprob = (PROP_PROB_STRUCT){
					  .birth = (U08)(resamplingOrderProb * 0.33 * MAX_RAND_BYTE_VALUE),
					  .move = (U08)(resamplingOrderProb * (0.33 + 0.33) * MAX_RAND_BYTE_VALUE) ,
					  .death = (U08)(resamplingOrderProb * (0.33 + 0.33 + 0.165) * MAX_RAND_BYTE_VALUE),
					  .merge = (U08)(resamplingOrderProb * 1. * MAX_RAND_BYTE_VALUE) };

			/*************************************************************************************************/
			// Pre-calcalte terms for SEASON AND TRENDS Bases:Precaluating SEASON AND TREND terms 
			/*************************************************************************************************/
			I32 period = opt->io.meta.period;
			I32 nElem  = (N + period - 1) / period * period + 1;
			//DUMMY_CONST  *DUMMY; // Must be global;otherwise &SEASON is transient and will be lost.

			DUMMY_CONST* DUMMY		= &basis->bConst;
			DUMMY->period           = period;
			DUMMY->SQRT_N_div_n     = MyALLOC(*MEM, nElem, F32, 64) ;
			for (I32 n = 1; n < nElem; n++) {
				DUMMY->SQRT_N_div_n[n] = sqrtf(N) /sqrtf(n);
			}

			// When deseaonalize=TRUE, the full set of cos/sin are needed to compute seasonality in the pre-processing step
			if (opt->io.meta.deseasonalize) 	{
				I32 sMAXORDER = opt->io.meta.period / 2;
				DUMMY->TERMS  = MyALLOC(*MEM, N * sMAXORDER * 2L, F32, 64), 
				//term:       1, 2,  3,  . ..,   N
				//sqr_cum: 0, 1, 12, 123,     , 1-N
				//so the sum of term_i_to_j can be computed as  sqr_csum(j)-sqrt_csum(i-1).
				preCalc_terms_season(DUMMY->TERMS, NULL, N, opt->io.meta.period, sMAXORDER);
			}
			
			
			/****************************************************************************/
			//				Pre-compute the scaling factor
			/****************************************************************************/
			I32 sMAXNUMKNOT = basis->prior.maxKnotNum;
			I32 sMINSEPDIST = basis->prior.minSepDist;
			F32PTR scaleFactor = MyALLOC(*MEM, sMAXNUMKNOT + 1, F32, 0);
			{
				F32PTR MEMBUF1 = MODEL.curr.XtX;  //MODEL.curr.XtX used here a temp MEM buffer
				F32PTR MEMBUF2 = MODEL.curr.XtX + (sMAXNUMKNOT + 1);
				preCalc_scale_factor(scaleFactor, N, sMAXNUMKNOT, sMINSEPDIST, MEMBUF1, MEMBUF2);
			}
			basis->scalingFactor = scaleFactor;

			/////////////////////////////////////////////////////////////
			//TODO:::
			I32		NUMSEG    = opt->prior.seasonMaxKnotNum + 1;
			I32		MINORDER  = opt->prior.seasonMinOrder=1;
			I32		MAXORDER  = opt->prior.seasonMaxOrder=2;
			if (opt->prior.modelPriorType == 4) {
				F64PTR	priorVec = MyALLOC(*MEM, NUMSEG * MAXORDER, F64, 64);
				F64PTR	priorMat = MyALLOC(*MEM, NUMSEG * MAXORDER * NUMSEG, F64, 64);
				PreCaclModelNumber(MINORDER, MAXORDER, NUMSEG, N, opt->prior.seasonMinSepDist, priorMat, priorVec);

				basis->priorMat = priorMat = NULL;
				basis->priorVec = priorVec = NULL;
			}
	
		}
		else if (type == SVDID)
		{
			MODEL.vid = i;	
			MODEL.sid = i;

			basis->prior.minKnotNum = opt->prior.seasonMinKnotNum;
			basis->prior.maxKnotNum = opt->prior.seasonMaxKnotNum;
			basis->prior.minSepDist = opt->prior.seasonMinSepDist;
			basis->prior.minOrder   = opt->prior.seasonMinOrder;
			basis->prior.maxOrder   = opt->prior.seasonMaxOrder;
			isComponentFixed[i]     = basis->prior.minOrder == basis->prior.maxOrder && basis->prior.minKnotNum == 0 && basis->prior.maxKnotNum == 0;

			AllocInitBasisMEM(basis, N, K_MAX, MEM);

	
			basis->mcmc_Kstopping = K_MAX - ( basis->prior.maxOrder); 		//Changed for SVD 
			basis->mcmc_MoveStep  = opt->mcmc.maxMoveStepSize;

			/****************************************************************************/
			// Prepare and re-adjust the propoal probabilities of indiviudal MOVEs for the
			// trend and seasonal components
			/****************************************************************************/
		   // If resamplingOrderProb is zero,  resampling the order will be NOT conducted; this
		   // happens if minOrder=maxOrder (checked in PostCheckArgs in beastv2_in_args  where 
		   // resamplingOderProb is set to 0 if needed)
			F32 cumProbExcludeResamplingOrder = 1 - opt->mcmc.seasonResamplingOrderProb;
			basis->propprob = (PROP_PROB_STRUCT){
					  .birth = (U08)(cumProbExcludeResamplingOrder * 0.33 * MAX_RAND_BYTE_VALUE),
					  .move = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33) * MAX_RAND_BYTE_VALUE) ,
					  .death = (U08)(cumProbExcludeResamplingOrder * (0.33 + 0.33 + 0.165) * MAX_RAND_BYTE_VALUE),
					  .merge = (U08)(cumProbExcludeResamplingOrder * 1. * MAX_RAND_BYTE_VALUE) };

			/*************************************************************************************************/
			// Pre-calcalte terms for SEASON AND TRENDS Bases:Precaluating SEASON AND TREND terms 
			/*************************************************************************************************/
			// When deseaonalize=TRUE, the full set of cos/sin are needed to compute seasonality 		
			SVD_CONST* SVD = &basis->bConst;
			I32 sMAXORDER  = opt->io.meta.deseasonalize? opt->io.meta.period: basis->prior.maxOrder;
			SVD->TERMS     = MyALLOC(*MEM, N*sMAXORDER ,              F32, 64);
			SVD->SQR_CSUM  = MyALLOC(*MEM, (N + 1L) * sMAXORDER ,     F32, 64);
			CopyNumericArrToF32Arr(SVD->TERMS, opt->io.meta.svdTerms, N* sMAXORDER);
			{
				F32PTR ptr  = SVD->TERMS;
				F32PTR ptr1 = SVD->SQR_CSUM; 

				for (I32 order = 1; order <= sMAXORDER; order++) {				 
						*ptr1 = 0.f;						
						f32_copy(ptr, ptr1 + 1, N);         f32_cumsumsqr_inplace(ptr1 + 1, N);
						ptr  +=  N;
						ptr1 += N + 1;
				} // or (rI32 order = 1; order <= maxSeasonOrder; order++)
			}
			/****************************************************************************/
			//				Pre-compute the scaling factor
			/****************************************************************************/
			I32 sMAXNUMKNOT = basis->prior.maxKnotNum;
			I32 sMINSEPDIST = basis->prior.minSepDist;
			F32PTR scaleFactor = MyALLOC(*MEM, sMAXNUMKNOT + 1, F32, 0);
			{
				F32PTR MEMBUF1 = MODEL.curr.XtX;  //MODEL.curr.XtX used here a temp MEM buffer
				F32PTR MEMBUF2 = MODEL.curr.XtX + (sMAXNUMKNOT + 1);
				preCalc_scale_factor(scaleFactor, N, sMAXNUMKNOT, sMINSEPDIST, MEMBUF1, MEMBUF2);
			}
			basis->scalingFactor = scaleFactor;

			/////////////////////////////////////////////////////////////
			// TODO: the model prior adjustment factors are needed only for certain modelPriorType
			I32 NUMSEG    = opt->prior.seasonMaxKnotNum + 1;
			I32 MINORDER  = opt->prior.seasonMinOrder;
			I32 MAXORDER  = opt->prior.seasonMaxOrder;

			if (opt->prior.modelPriorType == 4) {
				F64PTR priorVec = MyALLOC(*MEM, NUMSEG * MAXORDER, F64, 64);
				F64PTR priorMat = MyALLOC(*MEM, NUMSEG * MAXORDER * NUMSEG, F64, 64);
				PreCaclModelNumber(MINORDER, MAXORDER, NUMSEG, N, opt->prior.seasonMinSepDist, priorMat, priorVec);
				basis->priorMat = priorMat;
				basis->priorVec = priorVec;
			}

		}
		else if (type == OUTLIERID) {

			MODEL.oid = i;

			basis->prior.maxKnotNum = opt->prior.outlierMaxKnotNum;

			isComponentFixed[i] = 0;

			AllocInitBasisMEM(basis, N, K_MAX, MEM);
	 
			//MODEL.b[OUTLIERID].TERMS = &TREND;
			basis->mcmc_Kstopping = K_MAX - 1;
			basis->mcmc_MoveStep  = opt->mcmc.maxMoveStepSize;

			/*************************************************************************************************/
			// Pre-calcalte terms for SEASON AND TRENDS Bases:Precaluating SEASON AND TREND terms 
			/*************************************************************************************************/
			OUTLIER_CONST* OUTLIER = &basis->bConst;
			OUTLIER->SQRTN   = sqrtf(N);
			OUTLIER->SQRTN_1 = sqrtf(N-1);
						
			/****************************************************************************/
			// Prepare and re-adjust the propoal probabilities of indiviudal MOVEs for the
			// trend and seasonal components
			/****************************************************************************/
			basis->propprob = (PROP_PROB_STRUCT){
				.birth = (U08)(0.33 * MAX_RAND_BYTE_VALUE),
				.move  = (U08)((0.33 + 0.33) * MAX_RAND_BYTE_VALUE),
				.death = (U08)(1 * MAX_RAND_BYTE_VALUE), };


		}

		basis->Propose			       = Get_Propose(type, opt);		
		basis->UpdateGoodVec	       = Get_UpdateGoodVec(type);
		basis->ComputeY			       = Get_ComputeY(type, opt);
		basis->ModelPrior		       = Get_ModelPrior(opt->prior.modelPriorType);
		basis->GenTerms				   = Get_GenTerms(type, &opt->prior);
		basis->CalcBasisKsKeK_TermType = Get_CalcBasisKsKeK(type, opt->prior.precPriorType);

	}// loop through individual bases


	// Choose the routinue to sample type based on the nubmer of candidate bases
	MODEL.PickBasisID = Get_PickBasisID(NumBasis, MODEL.oid >=0, isComponentFixed);

	// Initialze mem/values for precison prior parameters: Must run after setting up the period
	InitPrecPriorMEM(model,  opt, MEM);


	

#undef MODEL
}


#include "abc_000_warning.h"