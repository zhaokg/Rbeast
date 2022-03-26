#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if defined(MSVC_COMPILER)
#include "intrin.h"                //_rdstc
#endif

#include <stdio.h>	               //fprintf fopen FILE
#include <string.h>	               //memset memcpy
#include <time.h>
#include <math.h>

#include "abc_001_config.h"
#include "abc_mem.h"              // Independent of R/Matlab,  VC/GNU, or MY/MKL.
#include "abc_blas_lapack_lib.h"  // Slight dependence on the choice of VC/GNU. Dependence on MY/MKL. Independacne of R/Matlab.
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_timer.h"
#include "abc_mcmc.h"
#include "abc_mat.h"
#include "abc_rand.h"
#include "abc_vec.h"   // for f32_add_v_v2_vec_in_place, f32_diff_back,i32_increment_bycon_inplace i32_to_f32_scaelby_inplace, f32_sx_sxx_toavstd_inplace 
#include "abc_math.h"  // for fastexp, fastsqrt only

#include "globalvars.h"  

#include "beastv2_header.h"
#include "beastv2_func.h" 
#include "beastv2_model_allocinit.h" 
#include "beastv2_prior_precfunc.h" 
#include "beastv2_xxyy_allocmem.h" 
#include "beastv2_io.h" 

#define LOCAL(...) do{ __VA_ARGS__ } while(0);
//extern MemPointers* mem;
//time_t start, end; 
/*---------WINDOW-------------------------*/
#ifdef WIN64_OS	
#include "abc_win32_demo.h"
#endif
/*---------WINDOW-------------------------*/


/*---------WINDOW-------------------------*/
void beast2_main_corev4_gui()
{
#ifdef WIN64_OS	
	/*---------WINDOW-------------------------*/

	// A struct to track allocated pointers   
	// const MemPointers MEM;
	// do not use 'const' bcz Clang will asssume other fields as zeros (e.g., alloc, and alloc0).
	MemPointers MEM = (MemPointers){.init = mem_init,};
	MEM.init(&MEM);
	//MEM.checkHeader = 1;
	//mem = &MEM;

	VSLStreamStatePtr stream;

	// Get the Option parameters from the global pointer GLOBAL_OPTIONS, which was already 
	// set prior to this point (e.g., time-series length N, and number of pixels).	
	const BEAST2_OPTIONS_PTR	opt   = GLOBAL_OPTIONS;
	const BEAST2_EXTRA          extra = opt->extra;

	/*---------WINDOW-------------------------*/
	if (gData.N <= 0)
	{
		EnterCriticalSection(&gData.cs);
		gData.N			= opt->io.N;		
		gData.optStatus = NoUpdate;
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);
	}

	//TODO:useless
	if (gData.optStatus == NeedUpdate)
	{
		EnterCriticalSection(&gData.cs);
		//F32PTR  tmp = opt.yInputData;
		//memcpy(&opt, &beastOpt, sizeof(opt));
		//opt.yInputData = tmp;
		LeaveCriticalSection(&gData.cs);
	}

	if (gData.plotData[0][0] == NULL)
	{
		EnterCriticalSection(&gData.cs);
		while (gData.plotData[0][0] == NULL)
			SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);
		LeaveCriticalSection(&gData.cs);
	}

	/*---------WINDOW-------------------------*/
	
	typedef int QINT;
	const   QINT  q = 1L;
	//const QINT  q   = opt->io.q;
	

	// Pre-allocate memory to save samples for calculating credibile intervals	
	CI_PARAM     ciParam = {0,};
	CI_RESULT    ci[MAX_NUM_BASIS];
	if (extra.computeCredible) {
		ConstructCIStruct(	opt->mcmc.credIntervalAlphaLevel, opt->mcmc.samples, opt->io.N*opt->io.q,  //for MRBEAST
							opt->prior.numBasis,&MEM, &extra.fastCIComputation, &ciParam, ci );  
	}

	// Allocate MEMORY FOR BASIS VARIABLE: Initialzie two pointers to BASIS
	BEAST2_MODEL  MODEL = {0,};
	AllocInitModelMEM(&MODEL, opt, &MEM);

	//Initializing the random number generaotr	
	LOCAL( 	
		U64 seed = (opt->mcmc.seed == 0) ? TimerGetTickCount() : (opt->mcmc.seed+0x4f352a3dc);
		r_vslNewStream(&stream, VSL_BRNG_MT19937, seed);   	
	)
	

	const U32PTR  RND32        = MyALLOC(MEM, MAX_RAND_NUM,     U32, 64);
	const U16PTR  RND16        = MyALLOC(MEM, MAX_RAND_NUM * 2, U16, 64);
	const U08PTR  RND08        = MyALLOC(MEM, MAX_RAND_NUM * 4, U08, 64);	
	const F32PTR  RNDGAMMA     = MyALLOC(MEM, MAX_RAND_NUM,     F32, 64);	
	const U32PTR  RND32_END    = RND32		+ MAX_RAND_NUM - 7;
	const U16PTR  RND16_END    = RND16		+ MAX_RAND_NUM * 2 - 7;
	const U08PTR  RND08_END    = RND08		+ MAX_RAND_NUM * 4 - 7 -3;     //-3 bcz GenRandomBasis will also consume the stream bits
	const F32PTR  RNDGAMMA_END = RNDGAMMA	+ MAX_RAND_NUM - MODEL.nPrec-1L;
		
	// Allocate mem for current covariate/design matrix (Xt_mars), proposed new terms (Xnewterm),
	// and subset matrix corresponding to rows of missing values.		
	const F32PTR Xt_mars;
	const F32PTR Xnewterm;      //will be re-used as a temp mem for multiple purposes
	const F32PTR Xt_zeroBackup; //mem for saving Xrows of the missing rows	
	AllocateXXXMEM(&Xt_mars, &Xnewterm, &Xt_zeroBackup,&MODEL,opt,&MEM);

	// yInfo used to save the current time series to be processed
	BEAST2_YINFO     yInfo;
	AllocateYinfoMEM(&yInfo, opt, &MEM);
	
	// Allocate the output memory for a single chain (resultChain) and the averaged
	// result of all chains ( result)
	BEAST2_RESULT resultChain = { NULL,}, result={ NULL, };
	BEAST2_Result_AllocMEM(&resultChain, opt, &MEM); 	
	BEAST2_Result_AllocMEM(&result,      opt, &MEM);
	
	if (extra.computeCredible) {
		I32  Npad   = (opt->io.N + 7) / 8 * 8; Npad = opt->io.N;//Correct for the inconsitency of X and Y in gemm and gemv
		I32  XnewOffset = 0;
		for (I32 i = 0; i < MODEL.NUMBASIS; i++) {
			if (MODEL.b[i].type == SEASONID|| MODEL.b[i].type == DUMMYID || MODEL.b[i].type == SVDID)
				ci[i].result     = resultChain.sCI,		//season		
			    ci[i].newDataRow = Xnewterm + XnewOffset;			//season		
			else if (MODEL.b[i].type == TRENDID)
				ci[i].result     = resultChain.tCI,     //trend			
			    ci[i].newDataRow = Xnewterm  + XnewOffset;   //trend		
			else if (MODEL.b[i].type == OUTLIERID)
				ci[i].result     = resultChain.oCI,      //outlier
			    ci[i].newDataRow = Xnewterm + XnewOffset; //outlier       

			XnewOffset += Npad * q;   //FOR MRBEAST
		}
	} //NUMVAR_FOR_CI=3

	const CORESULT coreResults[MAX_NUM_BASIS];
	SetupPointersForCoreResults(coreResults, MODEL.b, MODEL.NUMBASIS, &resultChain);
		 
	const BEAST2_HyperPar  hyperPar = { .alpha_1 = opt->prior.alpha1,  .alpha_2 = opt->prior.alpha2,
										.del_1   = opt->prior.delta1,  .del_2   = opt->prior.delta2};

	/****************************************************************************/
	//		THE OUTERMOST LOOP: Loop through all the time series one by one
	/****************************************************************************/
	// Get conversion factors from counts to seceonds
	InitTimerFunc();
	StartTimer();
	SetBreakPointForStartedTimer();

	const PREC_FUNCS precFunc;
	SetUpPrecFunctions(opt->prior.precPriorType, opt->io.q, &precFunc);

	// Print a blank line to be backspaced by the follow
	if (extra.printProgressBar) {
		F32 frac = 0.0; I32 firstTimeRun = 1;
		printProgress(frac, extra.consoleWidth, Xnewterm, firstTimeRun);
	}
		
		
	
	//#define __DEBUG__
	#undef  __DEBUG__ 

	#ifdef __DEBUG__
		// Allocate a mem block and memset it to zero
		I32    N          = opt->io.N;
		I32    Npad       = (N + 7) / 8 * 8; Npad =  N;//Correct for the inconsitency of X and Y in gemm and gemv
		F32PTR flagSat    = MyALLOC0(MEM, N, I32, 64);
		F32PTR Xdebug    = MyALLOC0(MEM, Npad*(opt->prior.K_MAX+ opt->prior.K_MAX), I32, 64); // one Kmax for Xt_mars, and another for Xtbackup
	#endif


	//#define XtX_ByGroup XtX_ByGroup_FULL
	//#define MatxMat     MatxMat_FULL
	//#define MatxVec     MatxVec_FULL

	// Calculate the total number of time sereies to be processed
	const U32  NUM_PIXELS    = opt->io.numOfPixels;
	const U32  MCMC_SAMPLES  = opt->mcmc.samples;
	const U32  MCMC_THINNING = opt->mcmc.thinningFactor;
	const U32  MCMC_BURNIN   = opt->mcmc.burnin;
	const U32  MCMC_CHAINNUM = opt->mcmc.chainNumber;

	NUM_OF_PROCESSED_GOOD_PIXELS  = 0; //this is a global variable.
	NUM_OF_PROCESSED_PIXELS       =0;  //this is also a global variable.
	for (U32 pixelIndex = 1; pixelIndex <= NUM_PIXELS; pixelIndex++)
	{
		// Fecth a new time-series: set up Y, nMissing,  n, rowsMissing		
		F32PTR MEMBUF           = Xnewterm; // Xnewterm is a temp mem buf.
		BEAST2_fetch_next_timeSeries(&yInfo, pixelIndex,  MEMBUF, &(opt->io));


		F32PTR  Xtmp            = Xt_mars;
		U08     skipCurrentPixel = BEAST2_preprocess_timeSeries(&yInfo, MODEL.b, Xtmp, opt);		
	

		#ifdef __DEBUG__
			I32 accS[5] = { 0, 0, 0, 0, 0 },  accT[5] = { 0, 0, 0, 0, 0 };
			I32 flagS[5] = { 0, 0, 0, 0, 0 }, flagT[5] = { 0, 0, 0, 0, 0 };
			for (int i = 0; i < yInfo.nMissing; i++) { flagSat[yInfo.rowsMissing[i]] = getNaN();}
		#endif

		#define __START_IF_NOT_SKIP_TIMESESIRIES__    
		#define __END_IF_NOT_SKIP_TIMESESIRIES__                        


		__START_IF_NOT_SKIP_TIMESESIRIES__  
		if (!skipCurrentPixel) {

		if (q == 1) {
				// for BEASTV4

				// alpha2_star  = alpha_2 + 0.5(YtY-beta*X'Y) = 0.5* (  [YtY+2*alpha_2] - beta*X'Y  )
				// YtY+2*alpha_2  is pre-cacluated here. THe "2" before alpha_2 is to undo the division
				// later in the calcution of alpha2_star
				yInfo.YtY_plus_alpha2Q[0] = yInfo.YtY_plus_alpha2Q[0] + 2 * hyperPar.alpha_2;
				//Pre-compute alpha1_start, which depends only on n.
				yInfo.alpha1_star = yInfo.n * 0.5 + hyperPar.alpha_1;

		}	else {
				// For MRBEAST

				 //YtY was computed from MV_Fetch; here alpaha_Q is added to its diagonal.			
				f32_add_val_matrixdiag(yInfo.YtY_plus_alpha2Q, hyperPar.alpha_2, q);
				yInfo.alpha1_star = yInfo.n * 0.5 + (hyperPar.alpha_1 + q - 1) * 0.5;
		}
		/*---------WINDOW-------------------------*/
			EnterCriticalSection(&gData.cs);
		{
			int idx;
			int N = opt->io.N;
			I32  Npad = (opt->io.N + 7) / 8 * 8;
			r_ippsMaxIndx_32f(yInfo.Y, N, &gData.yMax, &idx);
			r_ippsMinIndx_32f(yInfo.Y, N, &gData.yMin, &idx);

			gData.yMin = gData.yMin - (gData.yMax - gData.yMin) / 10;
			gData.yMax = gData.yMax + (gData.yMax - gData.yMin) / 10;
			gData.y           = yInfo.Y;
			gData.rowsMissing = yInfo.rowsMissing;
			gData.nMissing    = yInfo.nMissing;

			if (opt->io.meta.hasSeasonCmpnt) {
				gData.s    = coreResults[0].x;
				gData.curs = Xnewterm;
				gData.S    = MODEL.b[0].KNOT;
				gData.sCI   = ci[0].result;
				gData.sProb = resultChain.scpOccPr;

				gData.t   = coreResults[1].x;
				gData.ct  = Xnewterm + Npad;
				gData.T   = MODEL.b[1].KNOT;
				gData.tCI = ci[1].result;
				gData.tProb = resultChain.tcpOccPr;
			} else {
			    gData.s     = NULL;
				gData.curs  = NULL;
				gData.S     = NULL;
				gData.sCI   = NULL;
				gData.sProb = NULL;

				gData.t     = coreResults[0].x;
				gData.ct    = Xnewterm + Npad;
				gData.T     = MODEL.b[0].KNOT;
				gData.tCI   = ci[0].result;
				gData.tProb = resultChain.tcpOccPr;
			}
		}
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);



		if (gData.hwnd == NULL)
		{
			EnterCriticalSection(&gData.cs);
			while (gData.hwnd == NULL)
				SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);
			LeaveCriticalSection(&gData.cs);
		}
		//For seme reason, the stewindowpos in beast2_winmain is not always working. Put hre
		//Again
		SetWindowPos(gData.hwnd, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
		/*---------WINDOW-------------------------*/

		/****************************************************************************/
		// GENERATE A STREAM OF RANDOM NUMBERS FOR FUTURE USE
		/****************************************************************************/
		BEAST2_RNDSTREAM  RND;
		{
			RND.rnd32 = RND32, RND.rnd16 = RND16, RND.rnd08 = RND08, RND.rndgamma = RNDGAMMA;
			//vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE,  stream, MAX_RAND_NUM, rnd, 0, 1.0);		
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND32);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND16);
			r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, MAX_RAND_NUM, (U32PTR)RND08);
			// Depends on n and alpha_1 to generate rndgamma
			// Needed to resmaple Sig2	
			// Cann't be used to sample prec bcz the degree of freedom there depends on the number of terms Kterms
			r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, MAX_RAND_NUM, RNDGAMMA, ( hyperPar.alpha_1 + yInfo.n * 0.5f), 0, 1);
		}

		 
		// Clear up and zero out A(RESULT) for initialization	 
		BEAST2_Result_FillMEM(&result, opt, 0);		

		/****************************************************************************/
		//Iterate all the chains. The individual chain result will be saved into resultChain
		/****************************************************************************/
		for ( U32 chainNumber =0;  chainNumber < MCMC_CHAINNUM; chainNumber++)
		{
			/*---------WINDOW-------------------------*/
			//A timer to compute elapsing time interval
			LARGE_INTEGER tStart, tEnd, tFrequency;
			QueryPerformanceCounter(&tStart);
			QueryPerformanceFrequency(&tFrequency);
			/*---------WINDOW-------------------------*/
			const I32  N      = opt->io.N; 
			//const I32  Npad   = (N + 7) / 8 * 8; 
			const I32  Npad = N;//Correct for the inconsitency of X and Y in gemm and gemv
			const I32  Npad16 = (N + 15) /16 * 16;	
			/****************************************************************************/
			//                 GENERATE AN INITIAL MODEL
			/****************************************************************************/			
			{   
				// Generate random knots (numknot,KNOTs and ORDERS)				
				GenarateRandomBasis(MODEL.b, MODEL.NUMBASIS, N, &RND);//CHANGE: nunKnot, ORDER, KNOT, K, KBase, Ks, Ke, or TERM_TYPE				
				 
				/*
				MODEL.b[0].nKnot = 1;
				MODEL.b[0].ORDER[0] = MODEL.b[0].ORDER[1] = 1;
				MODEL.b[0].KNOT[0] = 500;
				MODEL.b[0].KNOT[1] = N + 1;
				// Get Ks and Ke for individula segments of each components
				MODEL.b[0].CalcBasisKsKeK_TermType(&MODEL.b[0]);
				*/

				/*
				for (int i = 0; i < MODEL.NUMBASIS; ++i) {
					MODEL.b[i].nKnot = 0;
					MODEL.b[i].KNOT[0] = N + 1;
					MODEL.b[i].ORDER[0] = MODEL.b[i].prior.minOrder;
					// Get Ks and Ke for individula segments of each components
					MODEL.b[i].CalcBasisKsKeK_TermType(&MODEL.b[i]);
				}
				*/
 
 		
				/////////////////////////////////
				// Update the Kbase for the bases after the 1st one. The first one is fixed at Kbase=0.					
				// Must be called after GenarateRandomBasis and before BEAST2_EvaluateModel
				MODEL.b[0].Kbase = 0;                           // This is the first time and also the only time b[0].Kbase is specified
				UpdateBasisKbase(MODEL.b, MODEL.NUMBASIS, 0);	//CHANGE: MODEL.b[i].Kbase for i>basisID				
								
				precFunc.GetNumTermsPerPrecGrp(&MODEL); //CHANGE: (1) nothing or (2) MODEL.curr.nTermsPerPrecGrp + MODEL.b[id].offsetPrec								
				precFunc.GetXtXPrecDiag(&MODEL);        //CHANGE: (1)nothing or (2) MODEL.curr.precXtXDiag
				
				// Find candidate positions for SEASON AND TREND				
				// nMissing & rowsMissing used for the outlier function
				CvtKnotsToBinVec(MODEL.b, MODEL.NUMBASIS, N, &yInfo);

				//Evaluate the initial model and compute its marg lik. CHANGE: DERIVE XMARS, K, BETA, BETA_MEAN, MARG_LIK
				// Xtmars is a cotinguous mem block consiting of Xtmars, Xnewterm, and Xt_zerobackup. The first part will
				// be filled with Xtmars, and the rest will be used as a temp block in this function call.
				// We don't use Xnewterm  or Xt_zerobackup as a temp mem buf bcz the zeroOutXmars function may need a much
				// larger  mem due to the many terms of the inital random model	
				// basis->K won't be updated inside the function and the old values from CalcBasisKsKeK is kept
				if (q == 1) {
					BEAST2_EvaluateModel(&MODEL.curr, MODEL.b, Xt_mars, N, MODEL.NUMBASIS, &yInfo, &hyperPar, opt->prior.precValue, &stream); 
				} else 	{
					MR_EvaluateModel(    &MODEL.curr, MODEL.b, Xt_mars, N, MODEL.NUMBASIS, &yInfo, &hyperPar, opt->prior.precValue, &stream);
				}
		
			}
		

			{
				// Clear up and zero out resultChain for initialization
				BEAST2_Result_FillMEM(&resultChain, opt, 0);

				// Reset all nots to ones bcz the real extrem positiosn won'tY be updated till samples > 1
				memset(MODEL.extremePosVec, 1, N);
				for (I32 i = 0; i < yInfo.nMissing; ++i) MODEL.extremePosVec[yInfo.rowsMissing[i]] = 0;
				MODEL.extremPosNum = yInfo.n;
			}

			/**********************************************************************************************/
			// PREPARE FOR THE START OF THE MAIN LOOP
			// The proposed basis needs XtX and XtY of the current basis to compute XtX_prop and XtY_prop.
			// Other varialbes (e.g., beta and cholXtX) are not cross-used by the basis and basis_prop
			/**********************************************************************************************/
			U32 ite            = 0;
			U32 sample         = 0;
			U32 subSampleIndex = 0;

			PROP_DATA PROPINFO = {.N=N,                   .Npad16 = Npad16,   .samples=&sample,
				                  .keyresult=coreResults, .mem    = Xnewterm,  .model  =&MODEL, 
				                  .pRND =&RND,            .yInfo  =&yInfo,    .nSample_ExtremVecNeedUpdate =1L,       
								  .sigFactor = opt->prior.sigFactor,          .outlierSigFactor = opt->prior.outlierSigFactor,
			}; 
			I32 numBadIterations = 0;
			while (sample < MCMC_SAMPLES)
			{
				ite++;
				/**********************************************************************/
				/*     Re-generate a pool of random numbers if almost used up          */
				/***********************************************************************/
				//vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD_ACCURATE, stream, MAX_RAND_NUM, rnd32, 0.f, 1.0f);
				if (RND.rnd32    >= RND32_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, (RND.rnd32 - RND32), (U32PTR)RND32);	                               RND.rnd32    = RND32;   }
				if (RND.rnd16    >= RND16_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd16 - (char*)RND16 + 3) / sizeof(U32), (U32PTR)RND16); RND.rnd16    = RND16;   }
				if (RND.rnd08    >= RND08_END)    {r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD, stream, ((char*)RND.rnd08 - (char*)RND08 + 3) / sizeof(U32), (U32PTR)RND08); RND.rnd08    = RND08;   }
				if (RND.rndgamma >= RNDGAMMA_END) {r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,      stream, MAX_RAND_NUM, RNDGAMMA, (hyperPar.alpha_1+yInfo.n*0.5f), 0.f, 1.f);    RND.rndgamma = RNDGAMMA;}
	
				//  Choose a basis type	
				BEAST2_BASIS_PTR basis = MODEL.b + MODEL.PickBasisID(&PROPINFO); //basisID = *RND.rnd08++ < 128;				
				
				// IMPLEMENT THE NEW PROPOSED BASIS		       
				NEWTERM    NEW;
				//CHANGE: new.newKnot, numSeg, SEG/R1/2, orders2, newIdx, nKnot_new, jumpType,propinfo.pRND.rnd8/32				
				basis->Propose(basis, &NEW, &PROPINFO); //info->mem=Xnewterm is used as a temp membuf here

				#ifdef __DEBUG__
					I32 basisIdx = basis - MODEL.b;
					if (basisIdx == 0 && (NEW.jumpType == BIRTH || NEW.jumpType == MOVE)) {
						flagSat[NEW.newKnot - 1] += 1;
					}
				#endif

				/**********************************************************************/
				// Generate the new terms for the propsed step: To add or move a bk, two segments 
				// are re-generated; to remove a bk or merge two bks, one segment are re-generated.	
				/**********************************************************************/								
				I32 Knewterm = 0;				
				for (I32 i = 0; i < NEW.numSeg; i++) { 
					// NEW.numSeg=0 if removing terms for ChORDER(newOrder<oldeTerm)
					// [1]..(1)....(2)...(sY-1)...N [N+1] 
					I32 kterms   = basis->GenTerms(Xnewterm+Npad*Knewterm, N, &NEW.SEG[i], &(basis->bConst));
					NEW.SEG[i].K = kterms;
					Knewterm    += kterms;
				} // Iterate through all the new segments
				#ifdef SOLARIS_COMPILER
					NEW.k1= NEW.k1_new=NEW.k1_old;
				#endif
				NEW.k2_new = NEW.k1 + Knewterm - 1L;	// if Knewterm=0 (i.e., delete terms), k2_new < k1_new
				
				//Get k1_old, k2_old, k1_new and k2_new
				NEW.k1     += basis->Kbase;			// k1=k1_new=k1_old
				NEW.k2_old += basis->Kbase;		
				NEW.k2_new += basis->Kbase;

				I32 KOLD     = MODEL.curr.K;                   // Total number of basis for the current model				
				I32 KNEW     = KOLD + NEW.k2_new - NEW.k2_old; // Total number of bases in the proposed model	
				MODEL.prop.K = KNEW;

				/*************************************************************************/
				// Get XtX_prop: Copy parts of XtX to XtT_prop and fill new components 
				/*************************************************************************/
				// Set those rows of Xt_mars_newterms specfied by rowsMissing  to zeros
				if (yInfo.nMissing > 0 && Knewterm > 0 /*&& basis->type != OUTLIERID*/)  //needed for basisFunction_OUliter=1
				f32_mat_multirows_extract_set_by_scalar(Xnewterm,Npad,Knewterm,Xt_zeroBackup, yInfo.rowsMissing, yInfo.nMissing, 0.0f);

				/*************************************************************************/
				//               The FIRST component:	
				/*************************************************************************/
				// There'sY no first component if k1_old/k1_new=1 for SEASON	
				for (I32 i = 1; i<NEW.k1; i++) SCPY(i, MODEL.curr.XtX + (i-1L)*KOLD, MODEL.prop.XtX + (i-1L)*KNEW);
				/*************************************************************************/
				//              The Second component
				/*************************************************************************/
				// No new cols/terms if flag=ChORDER && isInsert=0:the resampled basis has a higher order than the old basis
				if (Knewterm != 0) {
					FILL0(MODEL.prop.XtX + (NEW.k1-1)*KNEW, (KNEW-NEW.k1+1)*KNEW); // zero out the cols from k1-th to the end
					if (NEW.k1 > 1) {						
						/*
						rI32 r1 = NEW.r1[0];
					    rI32 r2 = NEW.r2[NEW.numSeg - 1];
					    rI32 Nseg = r2 - r1 + 1;
						r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
							NEW.k1 - 1L, Knewterm, Nseg, 1.0f, \
							Xt_mars  + r1-1L, Npad,
							Xnewterm + r1-1L, Npad, 0.f, \
							MODEL.prop.XtX + (NEW.k1 - 1L) * KNEW, KNEW); //0.f, MEMBUF1, k1_new - 1);
					   */
					  //r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1_new - 1, K_newTerm, N, 1.0f, Xt_mars , Npad,Xnewterm , Npad, 0.f, GlobalMEMBuf_1st, k1_new - 1);
					
						// Xnewterm is pre-allocated with sufficent mem to ensure segInfo won't overflow in __GetMAXNumElemXnewTerm
						BEAST2_BASESEG* _segInfo  = (BEAST2_BASESEG *)(Xnewterm + Knewterm * Npad);
						I32             _numBands = GetInfoBandList(_segInfo, &MODEL, NEW.k1 - 1);
						MatxMat(_segInfo, _numBands, Xt_mars,
							    NEW.SEG, NEW.numSeg, Xnewterm,
							    MODEL.prop.XtX+(NEW.k1-1L)*KNEW, N, KNEW );
					}

					//Three alternative ways to compute Xnewterm'*XnewTerm. Note that the resulting matrix is symmetric
					XtX_ByGroup(NEW.SEG, NEW.numSeg, Xnewterm, MODEL.prop.XtX + (NEW.k1 - 1) * KNEW + NEW.k1 - 1, N, KNEW);
					{
						//THE FIRST WAY:
						//r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, K_newTerm, N, 1.0, Xnewterm, Npad, Xnewterm, Npad, 0.f, GlobalMEMBuf_2nd, K_newTerm);
						/*
						r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
							Knewterm, Knewterm, Nseg, 1.0,
							Xnewterm + r1 - 1, Npad,
							Xnewterm + r1 - 1, Npad, 0.f,
							MODEL.prop.XtX + (NEW.k1-1) * KNEW + NEW.k1 - 1, KNEW);// MEMBUF2, K_newTerm);
						 */
						//XnewtermTXnewterm(&NEW, Xnewterm, MODEL.prop.XtX + (NEW.k1 - 1) * KNEW + NEW.k1 - 1, Npad, KNEW);
						
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
					}
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
				//                  The third component: 
				/*************************************************************************/
				//There is no third componet if k2_old=KOLD 
				if (NEW.k2_old != KOLD) {
					/*for (rI32  j = 1; i <= KOLD; i++, j++) {r_cblas_scopy(K_newTerm,  MEMBUF1 + (j - 1)*K_newTerm, 1, XtX_prop + (k - 1)*KNEW + k1_new - 1, 1),					*/
					for (I32 kold=NEW.k2_old+1, knew=NEW.k2_new+1; kold <= KOLD; kold++, knew++) {
						F32PTR ColStart_old = MODEL.curr.XtX + (kold - 1) * KOLD;
						F32PTR ColStart_new = MODEL.prop.XtX + (knew - 1) * KNEW;
						SCPY(NEW.k1 - 1,        ColStart_old,                        ColStart_new); //the upper part of the third componet
						SCPY(kold - NEW.k2_old, ColStart_old + (NEW.k2_old + 1) - 1, ColStart_new + (NEW.k2_new + 1) - 1); // the bottom part of the 3rd cmpnt
					}

					// If there is a middle part of the componet (i.e, Knewterm>0); this part 
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
						BEAST2_BASESEG* _segInfo  = (BEAST2_BASESEG*) (Xnewterm + Knewterm * Npad);
						I32             _numBands = GetInfoBandList_post(_segInfo, &MODEL, NEW.k2_old + 1);
			 
						MatxMat(NEW.SEG,NEW.numSeg,Xnewterm, _segInfo, _numBands, Xt_mars+NEW.k2_old*Npad,
							    MODEL.prop.XtX + (NEW.k2_new+1 -1)*KNEW + NEW.k1- 1, N, KNEW);				
					} 

				}

				/*********************************************************************************/
				//                 Compute XtY_prop from XtY
				/********************************************************************************/
				if (q == 1) {
					// Skipped if k1_old=1 when dealing with SEASON
					if (NEW.k1 > 1)       SCPY(NEW.k1 - 1, MODEL.curr.XtY, MODEL.prop.XtY);
					// New components : XnewTemrm*Y
					if (Knewterm > 0)       MatxVec(NEW.SEG, NEW.numSeg, Xnewterm, yInfo.Y, MODEL.prop.XtY + NEW.k1 - 1, N);
					//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
					if (NEW.k2_old != KOLD) SCPY(KNEW - NEW.k2_new, MODEL.curr.XtY + (NEW.k2_old + 1L) - 1L, MODEL.prop.XtY + (NEW.k2_new + 1) - 1);

				}	else {
				
					// Skipped if k1_old=1 when dealing with SEASON
					if (NEW.k1 > 1) {
						for (I32 c = 0; c < q; ++c) {
							SCPY( NEW.k1-1, MODEL.curr.XtY+MODEL.curr.K*c, MODEL.prop.XtY +KNEW * c);
						}						
					}
					// New components : XnewTemrm*Y
					if (Knewterm > 0) {
						MatxVec(NEW.SEG, NEW.numSeg, Xnewterm, yInfo.Y, MODEL.prop.XtY + NEW.k1 - 1, N);
						r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, Knewterm, q, N, 1.f, Xnewterm, Npad, yInfo.Y, N, 0.f,
							MODEL.prop.XtY+ NEW.k1 - 1, KNEW);
					}
					//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
					if (NEW.k2_old != KOLD) {
						for (I32 c = 0; c < q; ++c) {
							SCPY(KNEW - NEW.k2_new, MODEL.curr.XtY + (NEW.k2_old + 1L) - 1L + MODEL.curr.K * c, MODEL.prop.XtY + (NEW.k2_new + 1) - 1 + KNEW * c);
						}
					}
				 

				}

				/*************************************************************/
				// XtX_prop has been constructed. Now use it to get the marg lik
				/***************************************************************/			
				if (1L) {		
					//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec) 
					/*
					//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
					//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
					r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', KNEW, MODEL.prop.cholXtX, KNEW); // Choleskey decomposition; only the upper triagnle elements are used
					*/
					for (I32 i=1; i<NEW.k1; i++) 
						SCPY(i, MODEL.curr.cholXtX+(i-1)*KOLD, MODEL.prop.cholXtX+(i-1)*KNEW);

					precFunc.UpdateXtXPrec_nTermsPerGrp(&MODEL, basis, &NEW);
					precFunc.chol_addCol(  MODEL.prop.XtX  + (NEW.k1-1)*KNEW,
							               MODEL.prop.cholXtX,
							               MODEL.prop.precXtXDiag, KNEW, NEW.k1, KNEW);
					//chol_full_v2(MODEL.prop.XtX, MODEL.prop.cholXtX, KNEW, KNEW);
			       /*
					for (rI32 i = 1; i <= (NEW.k1_new - 1L); i++) 	r_cblas_scopy(i, MODEL.curr.cholXtX + (i - 1L) * KOLD, 1L, MODEL.prop.cholXtX + (i - 1L) * KNEW, 1L);
					chol_addCol(MODEL.prop.cholXtX+ (NEW.k1_new - 1L)*KNEW, MODEL.prop.cholXtX, KNEW, NEW.k1_new, KNEW);
					//chol_addCol(MODEL.prop.cholXtX + (1 - 1L) * KNEW, MODEL.prop.cholXtX, KNEW, 1, KNEW);
					*/

					/*{
					//LAPACKE_dpotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );			
					SCPY(KNEW, MODEL.prop.XtY, MODEL.prop.beta_mean);
					r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', KNEW, 1, MODEL.prop.cholXtX, KNEW, MODEL.prop.beta_mean, KNEW);
					}*/

				   precFunc.ComputeMargLik(&MODEL.prop, &MODEL, &yInfo, &hyperPar);
				   //if (MODEL.prop.marg_lik != MODEL.prop.marg_lik || fabs(MODEL.prop.marg_lik )>FLOAT_MAX || MODEL.prop.alpha2_star <0.f) {					   				  

				   if ( IsNaN(MODEL.prop.marg_lik) || IsInf(MODEL.prop.marg_lik ) ) {
					   if (++numBadIterations < 20) {
					    	   
						   precFunc.IncreasePrecValues(&MODEL);
						   precFunc.GetXtXPrecDiag(&MODEL);
						   precFunc.chol_addCol(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, MODEL.curr.K, 1L, MODEL.curr.K);
						   precFunc.ComputeMargLik(&MODEL.curr, &MODEL, &yInfo, &hyperPar);

						   #if !(defined(R_RELEASE) || defined(M_RELEASE))
						   r_printf("prec: %.4f| marg_lik_prop: %.4f | marg_like_curr: %.4f \n", MODEL.precVal, MODEL.prop.marg_lik, MODEL.curr.marg_lik);
						   #endif

						   continue;
					   }  else {
						   skipCurrentPixel = 2;
						   break;
					   }					   
				   }  else {
					   numBadIterations = 0;
				   } //if (marg_lik_prop != marg_lik_prop || alpha2_star_prop <0.f) 

				   if (q == 1) {// added for MRBEAST
					   MODEL.prop.alpha2_star = max(MODEL.prop.alpha2_star,  MIN_ALPHA2_VALUE);
				   }
				}

			   /****************************************************************************************/
			   /*    DETERMINE WHETHER OR NOT TO ACCCEPT THE PROPSOED STEP                              */
			   /****************************************************************************************/
			
				// First, calcuate a factor adjusting the likelihood change
				F32  factor;				
				if   ( NEW.jumpType ==MOVE || basis->type ==OUTLIERID) 	factor = 0.;
				else { factor = basis->ModelPrior(basis, &NEW, N); }

				F32 delta_lik = MODEL.prop.marg_lik - MODEL.curr.marg_lik + factor;
				
				//acceptTheProposal = *(RND.rnd16)++ < fastexp(delta_lik) * 65535.0f;				
				I08     acceptTheProposal;
				if      (delta_lik >   0)   acceptTheProposal = 1;
				else if (delta_lik < -23) 	acceptTheProposal = 0;				
				else {				 
					F32    expValue = fastexp(delta_lik);
					if     (delta_lik > -0.5) 	acceptTheProposal = *(RND.rnd08)++ < expValue * 255.0f;
					else if(delta_lik > -5  )   acceptTheProposal = *(RND.rnd16)++ < expValue * 65535.0f;
					else						acceptTheProposal = *(RND.rnd32)++ < expValue * 4.294967296e+09;
					 
				}
		 
				#ifdef __DEBUG__
					if (basisIdx == 0) ++(flagS[NEW.jumpType]);
					else 		       ++(flagT[NEW.jumpType]);
				#endif

				if(acceptTheProposal)
				{
					#ifdef __DEBUG__
						if (basisIdx == 0) ++(accS[NEW.jumpType]);
						else 		       ++(accT[NEW.jumpType]);
					#endif

					//Recover the orignal vaules for those rows corresponding to missing Y values
					if (yInfo.nMissing > 0 && Knewterm > 0 /*&& basis->type != OUTLIERID*/)  //needed for basisFunction_OUliter=1						
						f32_mat_multirows_set_by_submat(Xnewterm, Npad, Knewterm, Xt_zeroBackup, yInfo.rowsMissing, yInfo.nMissing);

					// Inserting XnewTerms into Xt_mars
					if (NEW.k2_old != KOLD && NEW.k2_new != NEW.k2_old)
						MoveCOLsWithinMatrix(Xt_mars, Npad, NEW.k2_old+1, KOLD, NEW.k2_new+1);
					if (Knewterm != 0)
						SCPY(Knewterm*Npad, Xnewterm, Xt_mars + (NEW.k1-1) * Npad);
					

					/****************************************************/
					//Find the good positions of the proposed MOVE
					//Then update the knotLists and order
					/****************************************************/
					basis->UpdateGoodVec(basis, &NEW, Npad16);
					basis->CalcBasisKsKeK_TermType(basis);
					UpdateBasisKbase(MODEL.b, MODEL.NUMBASIS, basis-MODEL.b);//basisIdx=basis-b Re-compute the K indices of the bases after the basisID 
	
					//Switching between Basis and Basis_prop
					{
						//http: //stackoverflow.com/questions/3647331/how-to-swap-two-numbers-without-using-temp-variables-or-arithmetic-operations
						//basis  = ((I64)basis ^ (I64)basis_prop);basis_prop = ((I64)basis ^ (I64)basis_prop); basis      =  ((I64)basis ^ (I64)basis_prop);						
					
						#define Exchange(x,y)   {void * _restrict tmp; tmp=MODEL.x;  MODEL.x=MODEL.y; MODEL.y=tmp;}
						Exchange(curr.XtX,       prop.XtX);
						Exchange(curr.XtY,       prop.XtY);
						Exchange(curr.beta_mean, prop.beta_mean);
						//Exchange(curr.beta,      prop.beta); //no need to exchange
						Exchange(curr.cholXtX,          prop.cholXtX);
						Exchange(curr.precXtXDiag,      prop.precXtXDiag);
						Exchange(curr.nTermsPerPrecGrp, prop.nTermsPerPrecGrp);

						if (q == 1) MODEL.curr.alpha2_star = MODEL.prop.alpha2_star; //BEASTV4
						else 	    Exchange(curr.alphaQ_star, prop.alphaQ_star);    //MRBEAST 
						
						MODEL.curr.marg_lik    = MODEL.prop.marg_lik;
						MODEL.curr.K		   = MODEL.prop.K;  //GetNumOfXmarCols(&MODEL): this function should also give KNEW; if not, there must be something wrong!
						#undef Exchange											
					}
					#ifdef __DEBUG__	
					if (q == 1) {
						//BEAST2_EvaluateModel(&MODEL.prop, MODEL.b, Xdebug, N, MODEL.NUMBASIS, &yInfo, &hyperPar, MODEL.precVal, &stream);											
						//r_printf("ite:%d K: |%f|%f|diff:%f\n", ite, MODEL.curr.K, MODEL.curr.marg_lik, MODEL.prop.marg_lik, MODEL.prop.marg_lik - MODEL.curr.marg_lik);
					    //r_printf(" %f[%f]-%f %f\n", (MODEL.curr.alpha2_star), (MODEL.prop.alpha2_star),	yInfo.alpha1_star* (log(MODEL.curr.alpha2_star) - log(MODEL.prop.alpha2_star)), yInfo.alpha1_star);
						
						/*
						I32 K = MODEL.prop.K;
						for (int i = 0; i < K; ++i) {

							r_printf("%f %f %f | %f %f %f\n",
								MODEL.prop.cholXtX[i * K + i], MODEL.curr.cholXtX[i * K + i],
								MODEL.prop.cholXtX[i * K + i] - MODEL.curr.cholXtX[i * K + i],
								MODEL.prop.XtX[i * K + i], MODEL.curr.XtX[i * K + i],
								MODEL.prop.XtX[i * K + i] - MODEL.curr.XtX[i * K + i]);

						}
						r_printf("ite----%d\n", ite);
						int a = 1;
						 */
					}
					else {
						MR_EvaluateModel(&MODEL.prop, MODEL.b, Xdebug, N, MODEL.NUMBASIS, &yInfo, &hyperPar, MODEL.precVal, &stream);
						//r_printf("MRite%d |%f|%f|diff:%f -prec %f\n", ite, MODEL.curr.marg_lik, MODEL.prop.marg_lik, MODEL.prop.marg_lik - MODEL.curr.marg_lik, MODEL.precVal);
					 
	                    /*
						I32 K = MODEL.prop.K;
						for (int i = 0; i < MODEL.prop.K; ++i) {
						 
								r_printf("%f %f %f | %f %f %f\n", 
									MODEL.prop.cholXtX[i*K + i], MODEL.curr.cholXtX[i * K + i],
									MODEL.prop.cholXtX[i * K + i] - MODEL.curr.cholXtX[i * K + i],
									MODEL.prop.XtX[i * K + i], MODEL.curr.XtX[i * K + i],
									MODEL.prop.XtX[i * K + i] - MODEL.curr.XtX[i * K + i]);
						 
						}
						
							r_printf("ite----%d\n",ite);
							int a = 1;
							*/ 
					}
					

					#endif


				} //(*rnd32++ < exp(marg_lik_prop - basis->marg_lik))
				
				

				/****************************************************************************************/
				// For buin-in iterations, no need to sample parameters or make posterior inference.
				/****************************************************************************************/
				if (ite <= MCMC_BURNIN) continue;
				
				U08 bResampleParameter  = (ite % 20            == 0);
				U08 bStoreCurrentSample = (ite % MCMC_THINNING == 0);
			
				/**********************/
				//First, Re-SAMPLING SIG2
				/**********************/
				if (bResampleParameter || bStoreCurrentSample)	{		

					if (q == 1) {
						//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_1+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
						MODEL.sig2 = (*RND.rndgamma++) * 1.f / MODEL.curr.alpha2_star;
						MODEL.sig2 = 1.0f / MODEL.sig2;
						MODEL.sig2 = max(MODEL.sig2, MIN_SIG2_VALUE);
						//r_printf("ite-%d SIG %f %f\n", ite, MODEL.sig2,  MODEL.sig2*yInfo.sd*yInfo.sd);
					}	else {
						// For MRBEAST
						F32PTR MEMBUF = Xnewterm;
						local_pcg_invwishart_upper( &stream, MODEL.SIG2, MODEL.SIG2+q*q, MEMBUF, q,
							                        MODEL.curr.alphaQ_star, hyperPar.alpha_1+ yInfo.n + q - 1);					
					}

				}
				/**********************/
				//Re-sample beta to be used for either re-sampling prec (ite%20=0) or predicting Y (ite%thiningFactor=0)
				/**********************/
				if (bResampleParameter || (bStoreCurrentSample && extra.useMeanOrRndBeta)) {

					if (q == 1) {
						//Compute beta = beta_mean + Rsig2 * randn(p, 1);
						//Usig2 = (1 / sqrt(sig2)) * U; 		beta = beta_mean + linsolve(Usig2, randn(p, 1), opts);
						//status = vdRngGaussian( method, stream, n, r, a, sigma );
						I32 K = MODEL.curr.K;
						r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, MODEL.beta, 0, 1);
						//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, MODEL.curr.cholXtX, K, MODEL.curr.beta, K); // LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
						solve_U_as_U_invdiag(MODEL.curr.cholXtX, MODEL.beta, K, K);
						r_ippsMulC_32f_I(fastsqrt(MODEL.sig2), MODEL.beta, K);
						r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K);
					} else {
					    // for MRBEAST
						F32PTR MEMBUF = Xnewterm;
						I32    K      = MODEL.curr.K;

						r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K * q, MEMBUF, 0., 1.);
						r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, K, q, q, 1.0, MEMBUF, K, MODEL.SIG2, q, 0.f, MODEL.beta, K);
						//LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
						//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, basis->post_P_U, K, beta, K);
						//r_ippsMulC_32f_I(sqrtf(modelPar.sig2), beta, K);
						//r_ippsAdd_32f_I(basis->beta_mean, beta, K);
						
						//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, q, post_P_U, K, beta, K);
						//r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K * q);
						solve_U_as_U_invdiag_multicols(MODEL.curr.cholXtX, MODEL.beta, K, K, q);
						r_ippsAdd_32f_I(MODEL.curr.beta_mean, MODEL.beta, K * q);
					}

				}

				/**********************/
				// Re-sample the precison parameters and re-calcuate marg_lik and beta
				/**********************/
				if (bResampleParameter && q==1) 
				{
				   /*
					I32   K     = MODEL.K;
					F32   sumq  = DOT(K, MODEL.curr.beta, MODEL.curr.beta);
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (hyperPar.del_1 + K * 0.5f), 0, 1.f);
					modelPar.prec[2]     = modelPar.prec[1] = modelPar.prec[0] = (*modelPar.prec) / (hyperPar.del_2 + 0.5f * sumq / MODEL.sig2);
					modelPar.LOG_PREC[2] = modelPar.LOG_PREC[1] = modelPar.LOG_PREC[0] = logf(modelPar.prec[0]);	
				
					if (ite % 100 == 0) {
						F32PTR beta = MODEL.curr.beta_mean;
						F32 f[4] = { 0, };
						I32 n[4] = { 0, };
						for (int j = 0; j < MODEL.b[0].K; j++) {
							rF32 tY = *beta++;
							int idx = MODEL.b[0].termType[j];
							f[idx - 1] += tY * tY;
							n[idx - 1] += 1;
						}
						r_printf("ite: %f %f %f %f", f[0]/(n[0]+0.0001),
							f[1] / (n[1] + 0.0001), f[2] / (n[2] + 0.0001), f[3] / (n[3] + 0.0001)			);

						 f[0] = f[1]= 0;
						 n[0] = n[1] = 0;
						for (int j = 0; j < MODEL.b[1].K; j++) {
							rF32 tY = *beta++;
							int idx = MODEL.b[1].termType[j];
							f[idx] += tY * tY;
							n[idx] += 1;
						}
						r_printf("|%f %f", f[0] / (n[0] + 0.0001),			f[1] / (n[1] + 0.0001));


						f[0] = f[1] = 0;
						n[0] = n[1] = 0;
						for (int j = 0; j < MODEL.b[2].K; j++) {
							rF32 tY = *beta++;							
							f[0] += tY * tY;
							n[0] += 1;
						}
						r_printf("|%f\n", f[0] / (n[0] + 0.0001));

					}
					*/
		   
			       #ifdef DEBUG_BEAST2
					if (ite %10==0) {
						 float sum=0;
						 int   n  =0;
							r_printf("ite %d: ",ite);
							for (int i = 1; i <= MODEL.nPrec; i++) {
								if (opt->prior.precPriorType < 2) {
									sum += MODEL.precVal;
									n++;
									r_printf("%7.4f | %7.4f ", MODEL.precVal, sum/n );
								}								
								else
									r_printf("%7.4f ", MODEL.precVec[i - 1]);
							}
							r_printf("\n ");
					}
				   #endif

				    /*
					rF32 f0 = 0, f1 = 0, f2 = 0;
					rU08PTR type = basis->termType;
					for (rI32 i = 0; i < K; i++)
					{
					rF32 f = basis->beta[i];
					f = f*f;
					U08  type_cur = *type++;
					if (type_cur == 0)
					f0 += f;
					else if (type_cur == 1)
					f1 += f;
					else
					f2 += f;
					}


					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (modelPar.del_1 + K_SN *0.5f), 0, 1.f);
					modelPar.prec[0] = modelPar.prec[0] / (modelPar.del_2 + 0.5f*f0 / modelPar.sig2);

					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec + 1, (modelPar.del_1 + basis->k_const *0.5f), 0, 1.f);
					modelPar.prec[1] = modelPar.prec[1] / (modelPar.del_2 + 0.5f*f1 / modelPar.sig2);

					rI32 K_1st = (K - K_SN - basis->k_const);
					if (K_1st > 0)
					{
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec + 2, (modelPar.del_1 + (K - K_SN - basis->k_const) *0.5f), 0, 1.f);
					modelPar.prec[2] = modelPar.prec[2] / (modelPar.del_2 + 0.5f*f2 / modelPar.sig2);

					}


					//r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (modelPar.del_1 + (K-basis->k_const) *0.5f), 0, 1.f);
					//modelPar.prec[2]=modelPar.prec[0] = modelPar.prec[0] / (modelPar.del_2 + 0.5f*(f0+f2) / modelPar.sig2);

					//r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec + 1, (modelPar.del_1 + (basis->k_const) *0.5f), 0, 1.f);
					//modelPar.prec[1] = modelPar.prec[1] / (modelPar.del_2 + 0.5f*(f1) / modelPar.sig2);

					//----------------------------------------------------------------------------------------------
					////r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (modelPar.del_1 + ( basis->k_SN) *0.5f), 0, 1.f);
					////modelPar.prec[0] = modelPar.prec[0] / (modelPar.del_2 + 0.5f*(f0 ) / modelPar.sig2);

					//r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec + 1, (modelPar.del_1 + (K - basis->k_SN) *0.5f), 0, 1.f);
					////modelPar.prec[2]=modelPar.prec[1] = modelPar.prec[1] / (modelPar.del_2 + 0.5f*(f1 + f2) / modelPar.sig2);

					if ( ite%200==0)
					{
					r_printf("%d %f %f %f \n", (K - K_SN - basis->k_const), 100 * f0 / K_SN, 100 * f1 / basis->k_const, 100 * f2 / (K - K_SN - basis->k_const));
					}
					*/

					/*
					//X_mars_prop has been constructed. Now use it to calcuate its marginal likelihood
					//Add precison values to the diagonal of XtX: post_P=XtX +diag(prec)				
					SCPY(K * K, MODEL.curr.XtX, MODEL.curr.cholXtX);

					{//Add precison values to the SEASONAL diagonal compoents		
						//rU08PTR termType = MODEL.termType;			
						for (rI32 i = 1, j = 0; i <= K; i++)						{
							//MODEL.curr.cholXtX[j + (i)-1] += modelPar.prec[*termType++];
							MODEL.curr.cholXtX[j + (i)-1] += modelPar.prec[0];
							j += K;
						}
					}//Add precison values to the SEASONAL diagonal compoents

					//Solve inv(Post_P)*XtY using  Post_P*b=XtY to get beta_mean
					//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, double * a, lapack_int lda);
					r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U', K, MODEL.curr.cholXtX, K); // Choleskey decomposition; only the upper triagnle elements are used
					//chol_addCol(MODEL.curr.cholXtX, MODEL.curr.cholXtX, K, 1, K);

					//LAPACKE_spotrs (int matrix_layout , char uplo , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );			
					SCPY(K, MODEL.curr.XtY, MODEL.curr.beta_mean);
					r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U', K, 1, MODEL.curr.cholXtX, K, MODEL.curr.beta_mean, K);
					*/

					I32 ntries = 0;
					do {
						if (ntries++ == 0)	precFunc.ResamplePrecValues(&MODEL, &hyperPar, &stream);							
						else				precFunc.IncreasePrecValues(&MODEL);											
						precFunc.GetXtXPrecDiag( &MODEL);
						precFunc.chol_addCol(    MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, MODEL.curr.K, 1L, MODEL.curr.K);		
						precFunc.ComputeMargLik( &MODEL.curr, &MODEL, &yInfo, &hyperPar);

					} while (  IsNaN(MODEL.curr.marg_lik) && ntries < 20 );

					if ( IsNaN(MODEL.curr.marg_lik) ) {
						#if !(defined(R_RELEASE) || defined(M_RELEASE)) 
						r_printf("skip3 | prec: %.4f| marg_lik_cur: %.4f \n",  MODEL.precVal, MODEL.curr.marg_lik);
						#endif
						skipCurrentPixel = 3;
						break;
					} 

					/* No need to re-sample beta because it is not really used
					r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, K, beta, 0, 1);
					// LAPACKE_strtrs (int matrix_layout , char uplo , char trans , char diag , lapack_int n , lapack_int nrhs , const double * a , lapack_int lda , double * b , lapack_int ldb );
					r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', K, 1, cholXtX, K, beta, K);
					r_ippsMulC_32f_I(sqrtf(modelPar.sig2), beta, K);
					r_ippsAdd_32f_I(beta_mean, beta, K);
					*/
					/* /FInally, re-sample sig2 based on the lastest alpha2_star
					//vdRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &sig2, (alpha_2+n/2), 0, 1.0/(alpha_1+basis->alpha2_star *0.5));
					modelPar.sig2 = (*rndgamma++)*1.0f / (modelPar.alpha_1 + basis->alpha2_star *0.5f);
					modelPar.sig2 = 1.f / modelPar.sig2;
					*/
				}

				if (bResampleParameter && q>1) 
				{
					F32PTR MEMBUF = Xnewterm;					
					I32    K      = MODEL.curr.K;
					//FLOAT_SHARE.sumq = DOT(K, basis->beta, basis->beta);
					//r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, modelPar.prec, (modelPar.alpha_2 + K *0.5f), 0, 1.f);
					//modelPar.prec[2] = modelPar.prec[1] = (*modelPar.prec) / (modelPar.alpha_1 + 0.5f*FLOAT_SHARE.sumq / modelPar.sig2);

					// Get trace( B*inv(SIG2)*B')						
					//r_LAPACKE_strtrs(LAPACK_COL_MAJOR, 'U', 'N', 'N', q, q, basis->alpha_Q_star, q, W_L, q);	
					r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, K, q, q, 1.0, MODEL.beta, K, MODEL.SIG2+q*q, q, 0.f, MEMBUF, K);
					F32 sumq = DOT(K * q, MEMBUF, MEMBUF);
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, stream, 1, &MODEL.precVal, (hyperPar.del_1 + K * q * 0.5f), 0.f, 1.f);
					MODEL.precVal    = MODEL.precVal / (hyperPar.del_2 + 0.5f * sumq);
					MODEL.logPrecVal = logf(MODEL.precVal);
					 	
					
					I32 ntries = 0;
					do {
						if (ntries++ != 0) {
							precFunc.IncreasePrecValues(&MODEL);
						}
						//precFunc.ResamplePrecValues( &MODEL, &hyperPar,&stream);
						precFunc.GetXtXPrecDiag(&MODEL);
						precFunc.chol_addCol(MODEL.curr.XtX, MODEL.curr.cholXtX, MODEL.curr.precXtXDiag, K, 1, K);
						precFunc.ComputeMargLik(&MODEL.curr, &MODEL, &yInfo, &hyperPar);
					} while (IsNaN(MODEL.curr.marg_lik) && ntries < 20);

					
					if ( IsNaN(MODEL.curr.marg_lik) ) {
						#if !(defined(R_RELEASE) || defined(M_RELEASE)) 
						r_printf("skip4 | prec: %.4f| marg_lik_cur: %.4f \n",  MODEL.precVal, MODEL.curr.marg_lik);
						#endif
						skipCurrentPixel = 4;
						break;
					}  

				}

				if (!bStoreCurrentSample) continue;
				
				sample++;
				if (extra.printProgressBar && NUM_PIXELS == 1 && sample % 1000 == 0) {
					F32 frac = (F32)(chainNumber * MCMC_SAMPLES + sample) / (MCMC_SAMPLES * MCMC_CHAINNUM);
					printProgress(frac, extra.consoleWidth, Xnewterm, 0);
				}

				/**********************************************/
				//
				//      Start to compute final results
				//
				/**********************************************/
				*resultChain.marg_lik += MODEL.curr.marg_lik;
				if (q == 1) {
					*resultChain.sig2 += MODEL.sig2;
				}	else {
					// For MRBEAST
					F32PTR MEMBUF = Xnewterm;
					r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, q, 1.f, MODEL.SIG2, q, MODEL.SIG2, q, 0.f, MEMBUF, q);
					r_ippsAdd_32f_I(MEMBUF, resultChain.sig2,  q*q);
				}

				F32PTR BETA = extra.useMeanOrRndBeta == 0 ? MODEL.curr.beta_mean : MODEL.beta;
				{
					F32PTR MEMBUF1  = Xnewterm;					

					for (I32 i = 0; i < MODEL.NUMBASIS; ++i) 
					{
						BEAST2_BASIS_PTR  basis   = MODEL.b   +i;
						CORESULT        * result  = coreResults+i;

						I32        nKnot  = basis->nKnot;
						TKNOT_PTR  KNOT   = basis->KNOT;

						result->xNProb[nKnot] += 1L;
						//Counting probability of being breakpoints				
						for (I32 i = 0; i < nKnot; i++) result->xProb[ KNOT[i]-1 ] += 1L;

						//Summng up the harmonic orders or trend orders for individual seeasonal segments
						if (result->xorder != NULL) {
							TORDER_PTR  orderList = basis->ORDER;
							for (I32 i = 0; i <= nKnot; ++i) {
								I16 r1 = KNOT[i-1], r2 = KNOT[i]-1;
								r_ippsAddC_32s_ISfs(orderList[i], result->xorder+r1 - 1, r2 - r1 + 1, 0);							
							}
						}

						//Compute the averaged  signals
						if (q == 1) {
							//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K_SN, 1.f, Xt_mars, Npad, MODEL.curr.beta_mean, 1L, 0.f, MEMBUF1, 1L);
							basis->ComputeY(Xt_mars, BETA, MEMBUF1, basis, Npad);
							f32_add_v_v2_vec_inplace(MEMBUF1, result->x, result->xSD, N);
							MEMBUF1 += Npad;
						}	else {
							// for MRBEAST
							F32PTR 	X    = Xt_mars + basis->Kbase * Npad;
							F32PTR  beta = BETA+basis->Kbase;
							I32     K    = basis->K;
							//r_cblas_sgemv(CblasColMajor, CblasNoTrans, Npad, K, 1.f, X, Npad, beta, 1L, 0.f, Y, 1L);
							r_cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, q, K, 1.0f,
									     	      X, Npad,beta, MODEL.curr.K, 0.f, MEMBUF1, N);
							f32_add_v_v2_vec_inplace(MEMBUF1, result->x, result->xSD, N*q);
							MEMBUF1 += Npad*q;
						}
	
						
					}//for (rI32 i = 0; i < MODEL.NUMBASIS; i++)
				
				}
				/*********************************************************************/
				// At this point, MEMBUF1=season, MEMBUF1+Npad=trend, and MEMBUF1+Npad+Npad=outlier;
				/*********************************************************************/

				/********************************************/
				// Compute results for the seasonal cmpnt
				/********************************************/
				if(extra.computeSeasonAmp)
				{
					F32PTR           MEMBUF1 = Xnewterm + 3*Npad;
					F32PTR           MEMBUF2 = MODEL.prop.beta_mean; //re-used here as a temp mem buf.

					BEAST2_BASIS_PTR basis    = &MODEL.b[MODEL.sid];
					I32             knotNum  = basis->nKnot;
					TKNOT_PTR       knotList = basis->KNOT;
					
					//Summng up the per-segment harmonic magnitudes  	
					F32PTR       beta            = BETA;
					TORDER_PTR   orderList       = basis->ORDER;
					F32PTR       SEASON_SQR_CSUM = basis->bConst.season.SQR_CSUM;
					for (I32 i = 0; i <= knotNum; i++) {
						I32 r1 = knotList[i - 1];
						I32 r2 = knotList[i] - 1;
						I32    segLength     = r2-r1 + 1L;
						F32PTR seasonSqrCsum = SEASON_SQR_CSUM + 1L;
						I32    order2        = orderList[i] * 2L;
						F32  amp = 0;
						for (I32 j = 0; j < order2; j++) {
							//TODO: re-check here
							F32 scalingFactor = N / (seasonSqrCsum[r2 - 1] - seasonSqrCsum[(r1 - 1) - 1]);
							amp               += (beta[j] * beta[j]) * scalingFactor;
							seasonSqrCsum += (N + 1LL);
						}						
						//r_ippsSubC_32f_I(-amp, resultChain.samp + r1 - 1, segLength, 0);

						r_ippsSet_32f(amp, MEMBUF1 + r1 - 1, segLength);
						beta += order2;
					}
					r_ippsAdd_32f_I(MEMBUF1, resultChain.samp,   N);
					r_ippsMul_32f_I(MEMBUF1, MEMBUF1,            N);
					r_ippsAdd_32f_I(MEMBUF1, resultChain.sampSD, N); //added to the square of the samp for computering SD

					
					if (extra.tallyPosNegSeasonJump)
					{//NEWLY ADDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD					 
						I32  posKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) { // It must be i<KnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (MEMBUF1[knot - 1] > MEMBUF1[knot - 1 - 1]) {
								resultChain.spos_cpOccPr[knot - 1] += 1;
								posKnotNum++;
							}
						}
						resultChain.spos_ncpPr[posKnotNum]           += 1L;
						resultChain.sneg_ncpPr[knotNum - posKnotNum] += 1L;
					}
				}
				
				/********************************************/
				// Compute results for the trend cmpnt
				/********************************************/
				if(extra.computeTrendSlope)
				{
					BEAST2_BASIS_PTR basis   = &MODEL.b[MODEL.tid];
					I32             knotNum  = basis->nKnot;
					TKNOT_PTR       knotList = basis->KNOT;

					F32PTR TREND = Xnewterm + Npad * MODEL.tid;     //trend signal
					F32PTR SLP   = Xnewterm + 3 * Npad;				//temp mem

																	// Compute the rate of change in trend based on beta. 
					f32_diff_back(TREND, SLP, N);
					f32_add_v_v2_vec_inplace(SLP, resultChain.tslp, resultChain.tslpSD, N); //added to the square of the trend signal for computing SD
					//i32_increment_bycond_inplace(resultChain.tslpSgnPosPr, SLP, N); //increamnent tslpSingPr if SLP is larger than 0				 
					i32_increment_vec2_bycond_inplace(resultChain.tslpSgnPosPr, resultChain.tslpSgnZeroPr, SLP, N); //increamnent tslpSingPr if SLP is larger than 0				 
					if (extra.tallyPosNegTrendJump ) {//NEWLY ADDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD					 
						I32  posKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (SLP[knot - 1] > 0) {
								resultChain.tpos_cpOccPr[knot - 1] += 1;
								posKnotNum++;
							}
						}
						resultChain.tpos_ncpPr[posKnotNum]           += 1L;
						resultChain.tneg_ncpPr[knotNum - posKnotNum] += 1L;
					}

					if (extra.tallyIncDecTrendJump ){//NEWLY ADDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD					 						
						I32  incKnotNum = 0;
						for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
							I64 knot = knotList[i];
							if (knot >= 2 && SLP[(knot + 1) - 1] > SLP[(knot - 1) - 1]) {
								resultChain.tinc_cpOccPr[knot - 1] += 1;
								incKnotNum++;
							}
						}
						resultChain.tinc_ncpPr[incKnotNum]         += 1L;
						resultChain.tdec_ncpPr[knotNum-incKnotNum] += 1L;
					}
			
				}

				/********************************************/
				// Compute results for the outlier cmpnt
				/********************************************/
				if(extra.tallyPosNegOutliers)
				{//NEWLY ADDEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD	
					BEAST2_BASIS_PTR basis    = &MODEL.b[MODEL.oid];
					rI32             knotNum  = basis->nKnot;
					rTKNOT_PTR       knotList = basis->KNOT;

					const F32PTR OUTLIIER  = Xnewterm + Npad* MODEL.oid;
	 
					I32  posKnotNum = 0;
					for (I32 i = 0; i < knotNum; i++) {  // It must be i<sKnotNum bcz of dealing with knots only 
						I64 knot = knotList[i];
						if (OUTLIIER[knot - 1] > 0) {
							resultChain.opos_cpOccPr[knot - 1] += 1;
							posKnotNum++;
						}							
					}
					resultChain.opos_ncpPr[posKnotNum]			 += 1L;
					resultChain.oneg_ncpPr[knotNum - posKnotNum] += 1L;								
				}

				/*************************************************/
				// Compute ci intervals: new row of data have already calculated
				// and saved in Xnewterm, Xnewterm+Npad, and Xnweterm+2*Npad
				/*************************************************/
				if (extra.computeCredible)	{ 	
					if (extra.fastCIComputation) {
						//if (*rnd32++ < subsampleFraction*4.294967296000000e+09)
						if (*RND.rnd16++ < ciParam.subsampleFraction_x_INT16MAX)
							++subSampleIndex;
						else // The current sample not included. No need to insert it
							 // into the ci strips. So, just skip to the next iteration							
							continue;
					} else {
						subSampleIndex = sample;
					}
							
					/*???????At this point, MEMBUF1 stores the rate-of-change signal???*/
					if (subSampleIndex <= ciParam.nSamples)	{	
						for (int i=0; i<MODEL.NUMBASIS;i++)      
							InsertInitialRows(&ciParam, &ci[i], subSampleIndex); //season 				
					} else { //(sample > numCISample) 							
						 // New row of data for slope, seasonal, and trend components: MEMBUF1=slope over time	 					    
						for (int i = 0; i < MODEL.NUMBASIS; i++) 
							InsertNewRowToUpdateCI(&ciParam, &ci[i]);						
					}//smaple <=nunumCISample				

				} // if (extra.computeCredible)


			/*---------WINDOW-------------------------*/
				EnterCriticalSection(&gData.cs);

				gData.ite		= ite;
				gData.sample	= sample;
				gData.tKnotNum	= MODEL.b[1].nKnot;
				gData.sKnotNum  = MODEL.b[0].nKnot;;
				
 				if (gData.yMaxT == gData.yMinT || ite %200==0)
				{
					int idx;
					r_ippsMaxIndx_32f(gData.ct, N, &gData.yMaxT, &idx);
					r_ippsMinIndx_32f(gData.ct, N, &gData.yMinT, &idx);
					gData.yMinT = gData.yMinT - (gData.yMaxT - gData.yMinT) / 5;
					gData.yMaxT = gData.yMaxT + (gData.yMaxT - gData.yMinT) / 5;

					if (gData.curs != NULL) {
						r_ippsMaxIndx_32f(gData.curs, N, &gData.yMaxS, &idx);
						r_ippsMinIndx_32f(gData.curs, N, &gData.yMinS, &idx);
						gData.yMinS = gData.yMinS - (gData.yMaxS - gData.yMinS) / 5;
						gData.yMaxS = gData.yMaxS + (gData.yMaxS - gData.yMinS) / 5;
					}					
				}


				if (gData.sleepInterval > 5)
				{
					Sleep_ms(gData.sleepInterval);
					BEAST2_GeneratePlotData();
				} else	{
					QueryPerformanceCounter(&tEnd);
					if ((tEnd.QuadPart - tStart.QuadPart) * 1000 / tFrequency.QuadPart > gData.timerInterval) 					{
						BEAST2_GeneratePlotData();
						tStart = tEnd;
					}
				}

				while (gData.status == PAUSE && gData.quit == 0)
					SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);

				if (gData.quit)
				{
					LeaveCriticalSection(&gData.cs);
					MEM.free_all(&MEM);
					r_vslDeleteStream(&stream);

					return; 				}

				LeaveCriticalSection(&gData.cs);
				/*---------WINDOW-------------------------*/
			}//WHILE(sample<SAMPLE)

			/*****************************************************************************/
			//
			// One chain is done; then start post-processing the result of the chain     
			//
			/*****************************************************************************/

			if (!skipCurrentPixel)
			{
				int   sum;
				#define GetSum(arr) (r_ippsSum_32s_Sfs(arr, N, &sum, 0), sum) // parenthesis operator

				I32   sMAXNUMKNOT = MODEL.sid >=0? MODEL.b[MODEL.sid].prior.maxKnotNum:-9999999;
				I32   tMAXNUMKNOT = MODEL.tid>=0?  MODEL.b[MODEL.tid].prior.maxKnotNum:-9999999;
				I32   oMAXNUMKNOT = MODEL.oid>=0?  MODEL.b[MODEL.oid].prior.maxKnotNum:- 9999999;
				F32   inv_sample  = 1.f / sample;			
				

				*resultChain.marg_lik = *resultChain.marg_lik * inv_sample;
				// FOR MRBEAST
				for (int col = 0; col < q; col++)	{ 
					for (int i = 0; i < q; i++) {
						resultChain.sig2[col*q+i] = resultChain.sig2[col*q+i] * inv_sample * yInfo.sd[col] * yInfo.sd[i];
					}
				}
				

				if (MODEL.sid >= 0) {

						*resultChain.sncp = GetSum(resultChain.scpOccPr)* inv_sample; 
						i32_to_f32_scaleby_inplace(resultChain.sncpPr,	(sMAXNUMKNOT + 1), inv_sample);
						i32_to_f32_scaleby_inplace(resultChain.scpOccPr, N,		           inv_sample);	
					
						//FOR MRBEAST
						for (int i = 0; i < q; i++) {
							F32 offset = 0.0f;
							f32_sx_sxx_to_avgstd_inplace(resultChain.sY + i * N, resultChain.sSD + i * N, sample, yInfo.sd[i], offset, N);
						}

						if (extra.computeSeasonOrder) i32_to_f32_scaleby_inplace(resultChain.sorder, N, inv_sample);
						if (extra.computeSeasonAmp) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								F32 offset = 0.0f;
								f32_sx_sxx_to_avgstd_inplace(resultChain.samp+i*N, resultChain.sampSD + i*N, sample, yInfo.sd[i], offset, N);
							}							
						}
						if (extra.computeCredible) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {								
								r_ippsMulC_32f_I(yInfo.sd[i], resultChain.sCI+      N*i, N);
								r_ippsMulC_32f_I(yInfo.sd[i], resultChain.sCI+N*q + N*i, N);
							}							
						} 
				}

				///////////////TREND////////////////////////////////////////////////////
				if (MODEL.tid >= 0) {
						*resultChain.tncp = GetSum(resultChain.tcpOccPr)*inv_sample ;
						i32_to_f32_scaleby_inplace(resultChain.tncpPr, (tMAXNUMKNOT + 1),	inv_sample);
						i32_to_f32_scaleby_inplace(resultChain.tcpOccPr, N,					inv_sample);
						//FOR MRBEAST
						for (int i = 0; i < q; i++) {
							F32 offset = 0.0f;
							f32_sx_sxx_to_avgstd_inplace(resultChain.tY + i * N, resultChain.tSD + i * N, sample, yInfo.sd[i], yInfo.mean[i], N);
						}


						if (extra.computeTrendOrder) 	i32_to_f32_scaleby_inplace(resultChain.torder, N, inv_sample);						
						if (extra.computeTrendSlope) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								f32_sx_sxx_to_avgstd_inplace(resultChain.tslp+i*N, resultChain.tslpSD+ i*N, sample, yInfo.sd[i]/opt->io.meta.deltaTime, 0, N);
							}							
							i32_to_f32_scaleby_inplace(resultChain.tslpSgnPosPr, N*q, inv_sample);
							i32_to_f32_scaleby_inplace(resultChain.tslpSgnZeroPr, N*q, inv_sample);
						}
						if (extra.computeCredible) {
							//FOR MRBEAST
							for (int i = 0; i < q; i++) {
								f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.tCI +        N * i, N);
								f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], resultChain.tCI + N*q +  N * i, N);
								//r_ippsMulC_32f_I(,    resultChain.tCI+(2*N)*i, N + N),
								//r_ippsSubC_32f_I(-yInfo.mean[i], ); //ippsAddC_32f_I(yInfo.mean, result.tCI, N + N);
							}							
						}
						
				}

				///////////////OUTLIR////////////////////////////////////////////////////
				if (MODEL.oid >= 0) {					
					 *resultChain.oncp = inv_sample * GetSum(resultChain.ocpOccPr);
					 i32_to_f32_scaleby_inplace(resultChain.oncpPr,  (oMAXNUMKNOT + 1), inv_sample);
					 i32_to_f32_scaleby_inplace(resultChain.ocpOccPr, N,                inv_sample);
					 //FOR MRBEAST
					 for (int i = 0; i < q; i++) {
						 f32_sx_sxx_to_avgstd_inplace(resultChain.oY+i*N, resultChain.oSD+i*N, sample, yInfo.sd[i], 0, N);
					 }					 
					 if (extra.computeCredible) {
						 //FOR MRBEAST
						 for (int i = 0; i < q; i++) {
							 r_ippsMulC_32f_I(yInfo.sd[i], resultChain.oCI +      i*N, N );
							 r_ippsMulC_32f_I(yInfo.sd[i], resultChain.oCI + N*q+ i*N, N);
						 }						 
					 }	
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegSeasonJump && MODEL.sid>=0) {

					GetSum(resultChain.spos_cpOccPr);	
					*resultChain.spos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.sneg_ncp = *resultChain.sncp - *resultChain.spos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.spos_ncpPr, (sMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.sneg_ncpPr, (sMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.spos_cpOccPr, N, inv_sample);
 
					SCPY(N, resultChain.scpOccPr, resultChain.sneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.spos_cpOccPr, (F32PTR)resultChain.sneg_cpOccPr, N);   //NEWLY ADDED
				}

				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegTrendJump) {

					GetSum(resultChain.tpos_cpOccPr);	
					*resultChain.tpos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.tneg_ncp = *resultChain.tncp - *resultChain.tpos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.tpos_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tneg_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tpos_cpOccPr, N,				  inv_sample);
 
					SCPY(N, resultChain.tcpOccPr, resultChain.tneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.tpos_cpOccPr, (F32PTR)resultChain.tneg_cpOccPr, N);   //NEWLY ADDED
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyIncDecTrendJump) {

					GetSum(resultChain.tinc_cpOccPr);	   
					*resultChain.tinc_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.tdec_ncp = *resultChain.tncp - *resultChain.tinc_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.tinc_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tdec_ncpPr, (tMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.tinc_cpOccPr, N, inv_sample);

					SCPY(N, resultChain.tcpOccPr, resultChain.tdec_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.tinc_cpOccPr, (F32PTR)resultChain.tdec_cpOccPr, N);   //NEWLY ADDED
				}


				//..tallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJumtallyPosNegSeasonJum.................
				if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {

					GetSum(resultChain.opos_cpOccPr);	
					*resultChain.opos_ncp = inv_sample * sum; //NEWLY ADDED
					*resultChain.oneg_ncp = *resultChain.oncp - *resultChain.opos_ncp; //NEWLY ADDED	

					i32_to_f32_scaleby_inplace(resultChain.opos_ncpPr, (oMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.oneg_ncpPr, (oMAXNUMKNOT + 1), inv_sample);
					i32_to_f32_scaleby_inplace(resultChain.opos_cpOccPr, N, inv_sample);

					SCPY(N, resultChain.ocpOccPr, resultChain.oneg_cpOccPr); //NEWLY ADDED
					r_ippsSub_32f_I((F32PTR)resultChain.opos_cpOccPr, (F32PTR)resultChain.oneg_cpOccPr, N);   //NEWLY ADDED
				}

			}// Finish computing the result of the single chiain


			/**************************************************/
			//   Add up the individual chain to the Result
			/**************************************************/
			if (!skipCurrentPixel) 
			{
				I32   sMAXNUMKNOT = MODEL.sid >=0? MODEL.b[MODEL.sid].prior.maxKnotNum:-9999999;
				I32   tMAXNUMKNOT = MODEL.tid>=0?  MODEL.b[MODEL.tid].prior.maxKnotNum:-9999999;
				I32   oMAXNUMKNOT = MODEL.oid>=0?  MODEL.b[MODEL.oid].prior.maxKnotNum:- 9999999;

				//https://stackoverflow.com/questions/13216423/error-pasting-and-red-does-not-give-a-valid-preprocessing-token
				//#define _1(x)      *(result.##x) += *(resultChain.##x) // Working with MSVC but not GCC

				#define _1(x)      *(result.x) += *(resultChain.x)
				#define _N(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N)
				#define _Nq(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q)
				#define _q(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q)
				#define _q2(x)      r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, q*q)
				#define _2N(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N+N)
				#define _2Nq(x)     r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, N*q+N*q)
				#define _skn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, sMAXNUMKNOT + 1)
				#define _tkn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, tMAXNUMKNOT + 1)
				#define _okn_1(x)  r_ippsAdd_32f_I((F32PTR)resultChain.x, (F32PTR)result.x, oMAXNUMKNOT + 1)

				_1(marg_lik);
				_q2(sig2);  //Fpr MRBEAST
				
				if (MODEL.sid >= 0) {					
					_1(sncp); _skn_1(sncpPr);	     _N(scpOccPr); _Nq(sY); _Nq(sSD);
					if (extra.computeSeasonOrder)    _N(sorder);
					if (extra.computeSeasonAmp)      _N(samp), _N(sampSD);
					if (extra.computeCredible)       _2Nq(sCI);
				}

				if (MODEL.tid >= 0) {					
					_1(tncp); _tkn_1(tncpPr);	   _N(tcpOccPr); _Nq(tY); _Nq(tSD);
					if (extra.computeTrendOrder)   _N(torder);
					if (extra.computeTrendSlope)   _N(tslp), _N(tslpSD),_N(tslpSgnPosPr), _N(tslpSgnZeroPr);
					if (extra.computeCredible)     _2Nq(tCI);
				}

				if (MODEL.oid >= 0) {
					_1(oncp); _okn_1(oncpPr);	_N(ocpOccPr); _Nq(oY); _Nq(oSD);
					if (extra.computeCredible)   _2Nq(oCI);
				}

				if (extra.tallyPosNegSeasonJump && MODEL.sid >=0) {
					_1(spos_ncp);         _1(sneg_ncp); //NEWLY ADDED					
					_skn_1(spos_ncpPr);   _skn_1(sneg_ncpPr); //NEWLY ADDED
					_N(spos_cpOccPr);     _N(sneg_cpOccPr); //NEWLY ADDED	
				}

				if (extra.tallyPosNegTrendJump ) {
					_1(tpos_ncp);            _1(tneg_ncp); //NEWLY ADDED					
				    _tkn_1(tpos_ncpPr); _tkn_1(tneg_ncpPr); //NEWLY ADDED
					_N(tpos_cpOccPr);         _N(tneg_cpOccPr); //NEWLY ADDED 
				}
				
				if (extra.tallyIncDecTrendJump) {
					_1(tinc_ncp);         _1(tdec_ncp); //NEWLY ADDED					
					_tkn_1(tinc_ncpPr); _tkn_1(tdec_ncpPr); //NEWLY ADDED
					_N(tinc_cpOccPr);     _N(tdec_cpOccPr); //NEWLY ADDED
				}
				
				if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {
					_1(opos_ncp);            _1(oneg_ncp); //NEWLY ADDED					
					_okn_1(opos_ncpPr);      _okn_1(oneg_ncpPr); //NEWLY ADDED
					_N(opos_cpOccPr);        _N(oneg_cpOccPr); //NEWLY ADDED 
				}	
 
				#undef _1
				#undef _N
				#undef _Nq 
				#undef _q 
				#undef _q2
				#undef _2N
				#undef _2Nq
				#undef _skn_1
				#undef _tkn_1
			    #undef _okn_1
			}

			// Jump out of the chainumber loop
			if (skipCurrentPixel) {
				r_warning("\nWARNING(#%d):The max number of bad iterations exceeded. Can't decompose the current time series\n", skipCurrentPixel);
				break;
			}


			/*---------WINDOW-------------------------*/
			EnterCriticalSection(&gData.cs);

			gData.curChainNumber = chainNumber;

			if (opt->io.meta.hasSeasonCmpnt) {
				gData.sN = *resultChain.sncp;
			} else {
				gData.sN = 0;
			}
			gData.tN = *resultChain.tncp;
			PostMessage(gData.hwnd, WM_USER + 2, 0, 0);

			LeaveCriticalSection(&gData.cs);
			/*---------WINDOW-------------------------*/
		}
		/*********************************/
		// WHILE(chainNumber<chainNumber)
		/*********************************/

		__END_IF_NOT_SKIP_TIMESESIRIES__  
	}

	/******************************************************/
	//
	// Finish all the chains and now acverage all of them
	//
	/******************************************************/

		// Average the results from multiple chains
		if (MCMC_CHAINNUM >= 2 && !skipCurrentPixel)
		{
			I32  N = opt->io.N;			
	
			I32   sMAXNUMKNOT = MODEL.sid >= 0 ? MODEL.b[MODEL.sid].prior.maxKnotNum : -9999999;
			I32   tMAXNUMKNOT = MODEL.tid >= 0 ? MODEL.b[MODEL.tid].prior.maxKnotNum : -9999999;
			I32   oMAXNUMKNOT = MODEL.oid >= 0 ? MODEL.b[MODEL.oid].prior.maxKnotNum : -9999999;

			const F32 invChainNumber = 1.f /(F32)MCMC_CHAINNUM;

			#define _1(x)      *((F32PTR)result.x)*=invChainNumber
			#define _N(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N)
			#define _Nq(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,N*q)
			#define _q(x)      r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q)
			#define _q2(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x,q*q)
			#define _2N(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N+N)
			#define _2Nq(x)     r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, N*q+N*q)
			#define _skn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, sMAXNUMKNOT + 1)
			#define _tkn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, tMAXNUMKNOT + 1)
			#define _okn_1(x)  r_ippsMulC_32f_I(invChainNumber, (F32PTR)result.x, oMAXNUMKNOT + 1)

			F32 maxncpProb;	 
			_1(marg_lik);			
			_q2(sig2);
			if (MODEL.sid >= 0) {
				_1(sncp); _skn_1(sncpPr);	     _N(scpOccPr); _Nq(sY); _Nq(sSD); 
				if (extra.computeSeasonOrder)    _N(sorder);
				if (extra.computeSeasonAmp)     {_N(samp), _N(sampSD);}
				if (extra.computeCredible)       _2Nq(sCI);
 
				*result.sncp_mode   = f32_maxidx(result.sncpPr,       sMAXNUMKNOT + 1, &maxncpProb);
				*result.sncp_median = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.5);
				*result.sncp_pct90  = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.9);
				*result.sncp_pct10  = GetPercentileNcp(result.sncpPr, sMAXNUMKNOT + 1, 0.1);
			}

			if (MODEL.tid >= 0) {
				_1(tncp); _tkn_1(tncpPr);	     _N(tcpOccPr); _Nq(tY); _Nq(tSD); 
				if (extra.computeTrendOrder)     _N(torder);
				if (extra.computeTrendSlope)    { _N(tslp), _N(tslpSD), _N(tslpSgnPosPr), _N(tslpSgnZeroPr);}
				if (extra.computeCredible)       _2Nq(tCI);

				*result.tncp_mode   = f32_maxidx(result.tncpPr,       tMAXNUMKNOT + 1, &maxncpProb);
				*result.tncp_median = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.5);
				*result.tncp_pct90  = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.9);
				*result.tncp_pct10  = GetPercentileNcp(result.tncpPr, tMAXNUMKNOT + 1, 0.1);
			}

			if (MODEL.oid >= 0) {
				_1(oncp); _okn_1(oncpPr);	   
				_N(ocpOccPr); _Nq(oY); _Nq(oSD);
				if (extra.computeCredible)      _2Nq(oCI);

				*result.oncp_mode   = f32_maxidx(result.oncpPr,       oMAXNUMKNOT + 1, &maxncpProb);
				*result.oncp_median = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.5);
				*result.oncp_pct90  = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.9);
				*result.oncp_pct10  = GetPercentileNcp(result.oncpPr, oMAXNUMKNOT + 1, 0.1);
			}

			/****************************************************************/

			if (extra.tallyPosNegSeasonJump && MODEL.sid >= 0) {
				_1(spos_ncp);             _1(sneg_ncp); //NEWLY ADDED					
				_skn_1(spos_ncpPr);    _skn_1(sneg_ncpPr); //NEWLY ADDED
				_N(spos_cpOccPr);         _N(sneg_cpOccPr); //NEWLY ADDED	
			}


			if (extra.tallyPosNegTrendJump) {
				_1(tpos_ncp);            _1(tneg_ncp); //NEWLY ADDED					
				_tkn_1(tpos_ncpPr); _tkn_1(tneg_ncpPr); //NEWLY ADDED
				_N(tpos_cpOccPr);         _N(tneg_cpOccPr); //NEWLY ADDED 
			}


			if (extra.tallyIncDecTrendJump) {
				_1(tinc_ncp);         _1(tdec_ncp); //NEWLY ADDED					
				_tkn_1(tinc_ncpPr); _tkn_1(tdec_ncpPr); //NEWLY ADDED
				_N(tinc_cpOccPr);     _N(tdec_cpOccPr); //NEWLY ADDED
			}

			if (extra.tallyPosNegOutliers && MODEL.oid >= 0) {
				_1(opos_ncp);            _1(oneg_ncp); //NEWLY ADDED					
				_okn_1(opos_ncpPr); _okn_1(oneg_ncpPr); //NEWLY ADDED
				_N(opos_cpOccPr);         _N(oneg_cpOccPr); //NEWLY ADDED 
			}

 
			#undef _1
			#undef _N
			#undef _2N
			#undef _skn_1
			#undef _tkn_1
			#undef _okn_1
		}
		if (!skipCurrentPixel) {		
			I32  N = opt->io.N;  
			if (MODEL.sid >= 0) 							tsRemoveNaNs(result.sSD, N);
			if (MODEL.tid >= 0)								tsRemoveNaNs(result.tSD, N);
			if (MODEL.tid >= 0 && extra.computeTrendSlope)  tsRemoveNaNs(result.tslpSD, N);		 
			if (MODEL.oid >= 0) 							tsRemoveNaNs(result.oSD, N);
		}

		// Compute Changepoints and their confidence intervals
		if (!skipCurrentPixel)
		{
			I32     N         = opt->io.N;
			F32     nan       = getNaN();     //A sloppy way to get a NAN

			F32  	threshold = 0.001f;
			F32PTR	mem       = Xnewterm;  //Xnewterm must have a length larger than or equla to 5*N + max(sncp, tncp);
			I32PTR  cptList   = (I32PTR)mem + 5LL * N;
			F32PTR  cptCIList = (F32PTR)mem + 6LL * N;

			I32   cptNumber;
			I32   trueCptNumber;
			F32   maxncpProb;
			const F32 T0 = (F32)opt->io.meta.startTime;
			const F32 dT = (F32)opt->io.meta.deltaTime;

			I32   sMAXNUMKNOT = MODEL.sid >= 0 ? MODEL.b[MODEL.sid].prior.maxKnotNum : -9999999;
			I32   tMAXNUMKNOT = MODEL.tid >= 0 ? MODEL.b[MODEL.tid].prior.maxKnotNum : -9999999;
			I32   oMAXNUMKNOT = MODEL.oid >= 0 ? MODEL.b[MODEL.oid].prior.maxKnotNum : -9999999;

			I32   sMINSEPDIST = MODEL.sid >= 0 ? MODEL.b[MODEL.sid].prior.minSepDist : -9999999;
			I32   tMINSEPDIST = MODEL.tid >= 0 ? MODEL.b[MODEL.tid].prior.minSepDist : -9999999;
			I32   oMINSEPDIST = MODEL.oid >= 0 ? 1 : -9999999;

			//--------------Season----------------------------
			if (extra.computeSeasonChngpt && MODEL.sid >=0)  	{
				cptNumber      = sMAXNUMKNOT;
				trueCptNumber = FindChangepoint((F32PTR)result.scpOccPr, mem, threshold, cptList, cptCIList, N, sMINSEPDIST, cptNumber);
				//In the returned result, cptList is a list of detected changepoints, which are all zero-based indices.
				//That is, a changepoint occurring at t1 has a value of 0.

				for (int i = 0; i < trueCptNumber; i++) {
					*(result.scp + i)	  = (F32)(*(cptList + i)) * dT + T0;
					*(result.scpPr + i)   = (F32)mem[i];
					I32 cptLoc = cptList[i] == 0 ? 1 : cptList[i];
					for (int j = 0; j < q; ++j) {
						*(result.scpAbruptChange +j*sMAXNUMKNOT + i) =	result.sY[j*N+cptLoc] - result.sY[j*N+cptLoc - 1];
					}					
				}
				for (int i = trueCptNumber; i < sMAXNUMKNOT; i++) {
					*(result.scp + i)             = nan,
					*(result.scpPr + i)           = nan;
					for (int j = 0; j < q; ++j) {
						*(result.scpAbruptChange + j * sMAXNUMKNOT+ i) = nan;
					}					
				}
				for (int i = 0; i < trueCptNumber; i++)
					*(result.scpCI + i) = (F32)(*(cptCIList + i)) * dT + T0,
					*(result.scpCI + sMAXNUMKNOT + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;
				for (int i = trueCptNumber; i < sMAXNUMKNOT; i++)
					*(result.scpCI + i) = nan,
					*(result.scpCI + sMAXNUMKNOT + i) = nan;
			}

			//--------------Trend----------------------------
			if (extra.computeTrendChngpt) {
				cptNumber     = tMAXNUMKNOT;
				trueCptNumber = FindChangepoint((F32PTR)result.tcpOccPr, mem, threshold, cptList, cptCIList, N, tMINSEPDIST, cptNumber);
				for (int i = 0; i < trueCptNumber; i++) {
					*(result.tcp + i)          = (F32)(*(cptList + i)) * dT + T0,
					*(result.tcpPr + i)         = (F32)mem[i];
					I32 cptLoc = cptList[i] == 0 ? 1 : cptList[i];
					for (int j = 0; j < q; ++j) {
						*(result.tcpAbruptChange + j*tMAXNUMKNOT + i) = result.tY[j * N + cptLoc] - result.tY[j * N + cptLoc - 1];
					}
				}
				for (int i = trueCptNumber; i < tMAXNUMKNOT; i++) {
					*(result.tcp + i)   = nan,
					*(result.tcpPr + i) = nan;
					for (int j = 0; j < q; ++j) {
						*(result.tcpAbruptChange + j * tMAXNUMKNOT + i) = nan;
					}					
				}					
				for (int i = 0; i < trueCptNumber; i++)
					*(result.tcpCI + i) = (F32)(*(cptCIList + i)) * dT + T0,
					*(result.tcpCI + tMAXNUMKNOT + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;
				for (int i = trueCptNumber; i < tMAXNUMKNOT; i++)
					*(result.tcpCI + i) = nan,
					*(result.tcpCI + tMAXNUMKNOT + i) = nan;
			}
	 

			/**************************************************************************************************/
			#define GET_CHANGPOINTS(NcpProb, KNOTNUM, MINSEP, MAX_KNOTNUM, Y, CpOccPr, CP, CPPROB, CP_CHANGE, CP_CI)    \
			cptNumber     = MAX_KNOTNUM;  \
			trueCptNumber = FindChangepoint((F32PTR)CpOccPr, mem, threshold, cptList, cptCIList, N, MINSEP, cptNumber);\
			for (int i = 0; i < trueCptNumber; i++) {\
				*(CP + i)        = (F32) cptList[i]* dT + T0,\
				*(CPPROB+ i)     = (F32) mem[i];\
		         I32 cptLoc      = cptList[i] == 0 ? 1 : cptList[i];\
				 *(CP_CHANGE + i) = Y[cptLoc] - Y[cptLoc - 1];\
			}\
			for (int i = trueCptNumber; i <MAX_KNOTNUM; i++) {\
				*(CP        + i) = nan;\
				*(CPPROB    + i) = nan;\
				  *(CP_CHANGE + i) = nan;\
			}	\
			for (int i = 0; i < trueCptNumber; i++)\
				*(CP_CI+ i)                = (F32)cptCIList[i] * dT + T0,\
				*(CP_CI + MAX_KNOTNUM + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;\
			for (int i = trueCptNumber; i < MAX_KNOTNUM; i++)\
				*(CP_CI + i)               = nan,\
				*(CP_CI + MAX_KNOTNUM + i) = nan;
			/**************************************************************************************************/

			if (extra.tallyPosNegSeasonJump && MODEL.sid >= 0) {
				GET_CHANGPOINTS(result.spos_ncpPr ,result.spos_ncp, sMINSEPDIST, sMAXNUMKNOT, result.samp,
					result.spos_cpOccPr, result.spos_cp, result.spos_cpPr, result.spos_cpAbruptChange, result.spos_cpCI);

				GET_CHANGPOINTS(result.sneg_ncpPr,result.sneg_ncp, sMINSEPDIST, sMAXNUMKNOT, result.samp,
					result.sneg_cpOccPr, result.sneg_cp, result.sneg_cpPr, result.sneg_cpAbruptChange, result.sneg_cpCI);
			}			
			if (extra.tallyPosNegTrendJump) {
				GET_CHANGPOINTS(result.tpos_ncpPr, result.tpos_ncp, tMINSEPDIST, tMAXNUMKNOT, result.tY,
					result.tpos_cpOccPr, result.tpos_cp, result.tpos_cpPr, result.tpos_cpAbruptChange, result.tpos_cpCI);

				GET_CHANGPOINTS(result.tneg_ncpPr, result.tneg_ncp, tMINSEPDIST, tMAXNUMKNOT, result.tY,
					result.tneg_cpOccPr, result.tneg_cp, result.tneg_cpPr, result.tneg_cpAbruptChange, result.tneg_cpCI);
			}
			if (extra.tallyIncDecTrendJump) {
				GET_CHANGPOINTS(result.tinc_ncpPr,result.tinc_ncp, tMINSEPDIST, tMAXNUMKNOT, result.tslp,
					result.tinc_cpOccPr, result.tinc_cp, result.tinc_cpPr, result.tinc_cpAbruptChange, result.tinc_cpCI);

				GET_CHANGPOINTS(result.tdec_ncpPr,result.tdec_ncp, tMINSEPDIST, tMAXNUMKNOT, result.tslp,
					result.tdec_cpOccPr, result.tdec_cp, result.tdec_cpPr, result.tdec_cpAbruptChange, result.tdec_cpCI);
			}


			/**************************************************************************************************/
			#define OGET_CHANGPOINTS(NcpProb, KNOTNUM,MINSEP, MAX_KNOTNUM, PROBCURVE, CP, CPPROB,CP_CI)    \
			cptNumber     = MAX_KNOTNUM;  \
			trueCptNumber = FindChangepoint((F32PTR)PROBCURVE, mem, threshold, cptList, cptCIList, N, MINSEP, cptNumber);\
			for (int i = 0; i < trueCptNumber; i++) {\
				*(CP + i)        = (F32) cptList[i]* dT + T0;\
				*(CPPROB+ i)     = (F32) mem[i];\
		         I32 cptLoc      = cptList[i] == 0 ? 1 : cptList[i];\
			}\
			for (int i = trueCptNumber; i <MAX_KNOTNUM; i++) {\
				*(CP        + i) = nan;\
				*(CPPROB    + i) = nan;\
		    }\
			for (int i = 0; i < trueCptNumber; i++) \
				*(CP_CI+ i)                = (F32)cptCIList[i] * dT + T0,\
				*(CP_CI + MAX_KNOTNUM + i) = (F32)(*(cptCIList + trueCptNumber + i)) * dT + T0;\
			for (int i = trueCptNumber; i < MAX_KNOTNUM; i++)\
				*(CP_CI + i)               = nan,\
				*(CP_CI + MAX_KNOTNUM + i) = nan;
			/**************************************************************************************************/

			if (extra.computeOutlierChngpt) {
				OGET_CHANGPOINTS(result.oncpPr,result.oncp, oMINSEPDIST, oMAXNUMKNOT, result.ocpOccPr, result.ocp, result.ocpPr, result.ocpCI);
			}

			if (extra.tallyPosNegOutliers && MODEL.oid >=0) {
				OGET_CHANGPOINTS(result.opos_ncpPr,result.opos_ncp, oMINSEPDIST, oMAXNUMKNOT,
					result.opos_cpOccPr, result.opos_cp, result.opos_cpPr, result.opos_cpCI);

				OGET_CHANGPOINTS(result.oneg_ncpPr,result.oneg_ncp, oMINSEPDIST, oMAXNUMKNOT,
					result.oneg_cpOccPr, result.oneg_cp, result.oneg_cpPr, result.oneg_cpCI);
			}


		}
		

		// Recover the orignal Y values if the pixel is not skipped; the value will be needed below 
		// to dump into the output and compute R2 and RMSE. Here Xnewterm is used a buff
		// should not be touched until after the computation of R2 and RMSE
		if ( !skipCurrentPixel) {	
			I32  N  = opt->io.N;	
			I32  Nq = N * q;  // For MRBEAST

			for (int i = 0; i < q; ++i) {
				memcpy(Xnewterm+N*i, yInfo.Y + N*i, sizeof(F32)* N);
				f32_scale_inplace(yInfo.sd[i], yInfo.mean[i], Xnewterm+N*i, N);

				F32PTR NULL_BUF_FOR_VALUES = (Xnewterm + N * i) + N;
				if (yInfo.nMissing > 0) {
					f32_gatherVec_scatterVal_byindex(Xnewterm + N*i, yInfo.rowsMissing, NULL_BUF_FOR_VALUES, getNaN(), yInfo.nMissing);
				}
			}			
		}

		if (opt->extra.dumpInputData && result.data != NULL) {
			I32  N = opt->io.N;
			I32  Nq = N * q;  // For MRBEAST
			memcpy(result.data, Xnewterm, sizeof(F32)* Nq);
		}

		// Compute R2 and RMSE
		// At this point, yInfo.Y is not used any longer, so is re-used
		// to stote the sume of fitted S,T and/or O.
		// Xnewterm still contains the orignal data and shiuldn't not be touched
  		if (!skipCurrentPixel) { 
			
			I32  N  = opt->io.N ;
			I32  Nq = N * q;  // For MRBEAST

			I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
			I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
			I08 hasTrendCmpnt   = 1;
			//Xnewterm is used at this point and should not be touched!
			F32PTR BUF      = yInfo.Y;
			f32_fill_val(0., BUF, Nq);
	
			if (hasTrendCmpnt)   f32_add_vec_inplace(result.tY, BUF, Nq);
			if (hasSeasonCmpnt)  f32_add_vec_inplace(result.sY, BUF, Nq);
			if (hasOutlierCmpnt) f32_add_vec_inplace(result.oY, BUF, Nq);
		
			for (int j = 0; j < q; ++j) {
				//For MRBEAST
				F32 r = f32_corr_rmse_nan(BUF+N*j, Xnewterm + N*j, N, &result.RMSE[j]);
				result.R2[j] = r * r;
			}
			
		}

		// Smooth the changepoint occurance probability curve
		if (!skipCurrentPixel) {
			I32  N = opt->io.N;

			I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
			I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
			I08 hasTrendCmpnt   = 1;
 
			F32PTR BUF = Xnewterm;
			if (hasTrendCmpnt && opt->extra.smoothCpOccPrCurve) {	
				memcpy(BUF, result.tcpOccPr, sizeof(F32)* N);
				f32_sumfilter(BUF, result.tcpOccPr,N, opt->prior.trendMinSepDist);
			}
			if (hasSeasonCmpnt && opt->extra.smoothCpOccPrCurve) {
				memcpy(BUF, result.scpOccPr, sizeof(F32)* N);
				f32_sumfilter(BUF, result.scpOccPr, N, opt->prior.seasonMinSepDist);
			}	
		}


	    if (skipCurrentPixel) BEAST2_Result_FillMEM(&result, opt, getNaN());
		/*********************************************/
		//Write outputs to the mem array or files	
		/*********************************************/	
		if (!skipCurrentPixel) {
			int N = opt->io.N;
			for (int i = 0; i < q; ++i) {
				if (yInfo.Yseason) {
				//If Y has been deseaonalized, add it back
						r_ippsAdd_32f_I(yInfo.Yseason + N*i,   result.sY+N* i, N);
					if (result.sCI) {
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.sCI + 2*N*i, N);
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.sCI + 2*N*i+ N, N);
					}
					if (extra.dumpInputData)
						r_ippsAdd_32f_I(yInfo.Yseason + N * i, result.data+N*i, N);
				}
				if (yInfo.Ytrend) {
				//If Y has been detrended, add it back
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tY + N*i, N);
					if (result.tCI) {
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tCI+ 2*N*i,   N);
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.tCI+ 2*N*i+N, N);
					}
					if (extra.dumpInputData)
						r_ippsAdd_32f_I(yInfo.Ytrend + N * i, result.data + N * i, N);
				}
			} //for (int i = 0; i < q; ++i)
		}

		if(q ==1)
			BEAST2_WriteOutput(opt, &result, pixelIndex);
		else
		    MR_WriteOutput(opt, &result, pixelIndex);

		//if (!skipCurrentPixel)	NUM_OF_PROCESSED_GOOD_PIXELS++; //avoid the branch
		NUM_OF_PROCESSED_GOOD_PIXELS += !skipCurrentPixel;  //this is a global variable.
		NUM_OF_PROCESSED_PIXELS++;							//this is also a global variable.


		F64 elaspedTime = GetElaspedTimeFromBreakPoint();
		if (NUM_OF_PROCESSED_GOOD_PIXELS > 0 && NUM_PIXELS > 1 && (pixelIndex % 50 == 0 || elaspedTime > 1)) 		{
			F64 estTimeForCompletion = GetElapsedSecondsSinceStart()/NUM_OF_PROCESSED_GOOD_PIXELS * (NUM_PIXELS - pixelIndex);
			printProgress2((F32)pixelIndex / NUM_PIXELS, estTimeForCompletion, extra.consoleWidth, Xnewterm, 0);
			if (elaspedTime > 1) SetBreakPointForStartedTimer();
		}

		#ifdef __DEBUG__
		r_printf("TREND: birth%4d/%-5d|death%4d/%-5d|merge%4d/%-5d|move%4d/%-5d|chorder%4d/%-5d\n", 
			      accT[0], flagT[0] , accT[1], flagT[1], accT[2], flagT[2], accT[3], flagT[3], accT[4], flagT[4]);
		r_printf("SEASN: birth%4d/%-5d|death%4d/%-5d|merge%4d/%-5d|move%4d/%-5d|chorder%4d/%-5d\n",
			      accS[0], flagS[0], accS[1], flagS[1], accS[2], flagS[2], accS[3], flagS[3], accS[4], flagS[4]);
		#endif

	} //for (U32 pixelIndex = 1; pixelIndex <= TOTALNUMPIXELS; pixelIndex++)


	/***********************************************************/
	// This is the ending bracekt of the iteration through pixels
	/***********************************************************/

	r_vslDeleteStream(&stream);
	MEM.free_all(&MEM);
 
	/*---------WINDOW-------------------------*/
	PostMessage(gData.hwnd, WM_USER + 1, 0, 0);

	EnterCriticalSection(&gData.cs);
	gData.t = NULL;
	gData.y = NULL;
	gData.s = NULL;
	gData.rowsMissing = NULL;
	while (gData.status != DONE)
		SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);
	LeaveCriticalSection(&gData.cs);
#endif
	return;
	/*---------WINDOW-------------------------*/
	
} /* End of beastST() */


#include "abc_000_warning.h"