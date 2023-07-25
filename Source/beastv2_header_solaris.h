#pragma once

#include "abc_000_macro.h"
#include "abc_datatype.h"   //#include <inttypes.h>  #include <stdint.h>
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_mat.h"   //NEWXCOLINFO


#if R_INTERFACE==1 
#ifdef beta
#undef beta  //beta is refined to Rf_beta in Rmath.h, which causes errors in expanding model->beta
#endif
#endif

#define _IN_ 
#define _OUT_

#define MAX_RAND_NUM	  5000L
#define MAX_NUM_BASIS     3
#define MIN_PREC_VALUE    0.001
#define MIN_SIG2_VALUE    0.001
#define MIN_ALPHA2_VALUE  0.0001

#define PROB_SAMPLE_EXTREME_VECTOR 0.5

#define SEASONID   0
#define TRENDID    1
#define OUTLIERID  2 //oooooooooooooooooooo
#define DUMMYID    3
#define SVDID      4
#define A(xxx)     BEAST2_##xxx

typedef U32 TKNOT, * _restrict TKNOT_PTR;
typedef U08 TORDER, * _restrict TORDER_PTR;
#define rTKNOT_PTR   register  TKNOT_PTR
#define rTORDER_PTR  register  TORDER_PTR
#define rTKNOT       register  TKNOT
#define rTORDER      register  TORDER

/*************************************************/
// Declarations of data types for input options
/*************************************************/

typedef struct BEAST2_METADATA {
	I08      detrend;
	I08      deseasonalize;
	I08		 nrhs;
	I08      seasonForm;
	I08      hasSeasonCmpnt;
	I08      hasOutlierCmpnt;

	I08     whichDimIsTime;
	F32     period;
	F32     missingValue;
	F32     startTime, deltaTime;
	F32     maxMissingRate;
	VOIDPTR svdTerms_Object;
	VOIDPTR svdYseason_Object;
} BEAST2_METADATA, * _restrict BEAST2_METADATA_PTR;


// 0: Use a constant precison value
// 1: Use a uniform precision value for all terms
typedef enum { ConstPrec, UniformPrec, ComponentWise, OrderWise } PRECPRIOR_TYPE;

typedef struct BEAST2_HyperPar {
	F32 alpha_1, alpha_2, del_1, del_2;
} BEAST2_HyperPar, * _restrict BEAST2_HyperPar_PTR;

typedef struct BEAST2_PRIOR {
	I08	  basisType[MAX_NUM_BASIS];			//'ST','STO','TO','T', ...
	I08	  numBasis;

	U08	  seasonBasisFuncType;
	U08	  trendBasisFuncType;
	U08	  outlierBasisFuncType;
	U08	  modelPriorType;
	U08   precPriorType;

	I16   seasonMinOrder, seasonMaxOrder;
	I16   trendMinOrder, trendMaxOrder;
	I32   trendMinSepDist, seasonMinSepDist;

	I16   trendMinKnotNum, seasonMinKnotNum;
	I16   trendMaxKnotNum, seasonMaxKnotNum;
	I16   outlierMaxKnotNum;

	U16   K_MAX;

	F32   sigFactor;
	F32   outlierSigFactor;

	F32   sig2;
	F32	  precValue;
	F32   alpha1, alpha2, delta1, delta2;
} BEAST2_PRIOR, * _restrict BEAST2_PRIOR_PTR;

typedef struct BEAST2_MCMC {
	U64   seed;	                  // Unsigned long long seed;
	F32   credIntervalAlphaLevel;
	F32   ridgeFactor;
	F32   trendResamplingOrderProb;
	F32   seasonResamplingOrderProb;
	U16   maxMoveStepSize;
	U32   burnin, samples, chainNumber;
	U16   thinningFactor;
} BEAST2_MCMC, * _restrict BEAST2_MCMC_PTR;

typedef struct BEAST2_EXTRA {
	I08   smoothCpOccPrCurve;
	I08   useMeanOrRndBeta;
	I08   dumpInputData;
	U08   numThreadsPerCPU;
	U16   numParThreads;
	U16   numCPUCoresToUse;
	U16   consoleWidth;
	I08   whichOutputDimIsTime;
	I08   removeSingletonDims;

	I08   ncpStatMethod;
	Bool  computeCredible;
	Bool  fastCIComputation;
	Bool  computeSeasonOrder;
	Bool  computeTrendOrder;

	Bool  computeSeasonChngpt, computeTrendChngpt, computeOutlierChngpt;
	Bool  computeSeasonAmp;
	Bool  computeTrendSlope;

	Bool tallyPosNegSeasonJump;
	Bool tallyPosNegTrendJump;
	Bool tallyIncDecTrendJump;
	Bool tallyPosNegOutliers;

	Bool  printOptions;
	Bool  printProgressBar;

} BEAST2_EXTRA, * _restrict BEAST2_EXTRA_PTR;


struct BEAST2_RESULT;
typedef struct BEAST2_RESULT BEAST2_RESULT, * _restrict BEAST2_RESULT_PTR;
typedef struct BEAST2_IO {
	BEAST2_METADATA	meta;
	TimeVecInfo  	T;

	VOID_PTR* pdata;
	DATA_TYPE* dtype;
	I08				ndim;
	I08             rowdim, coldim, timedim;
	I32				imgdims[2];
	I32				dims[3];
	I32				numOfPixels;
	// q is the number of time series;q=1 is for univaraite TS in BEASTv4.
	// q is added to ensure a consistent API inteface btw BEAST and MRBEAST
	// For irregular time seires, Nraw is obtained from ndims[whichDimIsTime=1]
	I32				N, q;

	struct {
		BEAST2_RESULT* result;
		DATA_TYPE      dtype;
		U08            whichDimIsTime;
	} out;

} BEAST2_IO, * _restrict BEAST2_IO_PTR;

typedef struct BEAST2_OPTIONS {
	BEAST2_IO		    io;
	BEAST2_MCMC			mcmc;
	BEAST2_EXTRA		extra;
	BEAST2_PRIOR		prior;
} BEAST2_OPTIONS, * _restrict BEAST2_OPTIONS_PTR;

typedef struct BEAST2_YINFO {
	// maen, sd, and YtT are all scalars but pointers are used to maintan 
	// a consitent API with MRBEAST
	F32PTR    Yseason;
	F32PTR    Ytrend;

	F32PTR     mean, sd;
	F32PTR     YtY_plus_alpha2Q;
	F32        alpha1_star;
	TKNOT      n, nMissing;
	//q is added for MRBEAST and used in ComputeLik, PropseNew/__CalcAbsDeviation(compute deviaiton and extrempos)
	//, MR_EvaluateModel,
	I32        q;
	I32PTR     rowsMissing;
	F32PTR     Y;

} BEAST2_YINFO, * _restrict BEAST2_YINFO_PTR;

typedef struct BEAST2_RESULT {
	F32PTR  time;
	F32PTR  data;
	F32PTR  marg_lik, sig2, R2, RMSE;
	F32PTR  sncp, sncp_median, sncp_mode, sncp_pct90, sncp_pct10;
	I32PTR  sncpPr, scpOccPr;
	F32PTR  sY, sSD;
	F32PTR  sCI;            //computeCI
	I32PTR  sorder;         //computeOrder
	F32PTR  samp, sampSD;   //computeAmp
	F32PTR  scp, scpCI, scpPr, scpAbruptChange; //computeChangpt

	F32PTR  tncp, tncp_median, tncp_mode, tncp_pct90, tncp_pct10;
	I32PTR  tncpPr, tcpOccPr;
	F32PTR  tY, tSD;
	F32PTR  tCI;            //computeCI
	I32PTR  torder;         //computeOrder
	F32PTR  tslp, tslpSD;   //computeSlp
	I32PTR  tslpSgnPosPr;   //computeSlp
	I32PTR  tslpSgnZeroPr;   //computeSlp: THere is a sizable probability that slp ==0.

	F32PTR  tcp, tcpCI, tcpPr, tcpAbruptChange; //computeChangpt	

	F32PTR  oncp, oncp_median, oncp_mode, oncp_pct90, oncp_pct10;
	I32PTR  oncpPr, ocpOccPr;
	F32PTR  oY, oSD;
	F32PTR  oCI;                 //computeCI
	F32PTR  ocp, ocpCI, ocpPr;  //computeChangpt	

	//tallyPosNegSeasonJump
	F32PTR  spos_ncp, sneg_ncp;
	I32PTR  spos_ncpPr, sneg_ncpPr;
	I32PTR  spos_cpOccPr, sneg_cpOccPr;
	F32PTR  spos_cp, sneg_cp;
	F32PTR  spos_cpPr, sneg_cpPr;
	F32PTR  spos_cpAbruptChange, sneg_cpAbruptChange;
	F32PTR  spos_cpCI, sneg_cpCI;

	//tallyPosNegTrendJump
	F32PTR  tpos_ncp, tneg_ncp;
	I32PTR  tpos_ncpPr, tneg_ncpPr;
	I32PTR  tpos_cpOccPr, tneg_cpOccPr;
	F32PTR  tpos_cp, tneg_cp;
	F32PTR  tpos_cpPr, tneg_cpPr;
	F32PTR  tpos_cpAbruptChange, tneg_cpAbruptChange;
	F32PTR  tpos_cpCI, tneg_cpCI;

	//tallyIncDecTrendJump
	F32PTR  tinc_ncp, tdec_ncp;
	I32PTR  tinc_ncpPr, tdec_ncpPr;
	I32PTR  tinc_cpOccPr, tdec_cpOccPr;
	F32PTR  tinc_cp, tdec_cp;
	F32PTR  tinc_cpPr, tdec_cpPr;
	F32PTR  tinc_cpAbruptChange, tdec_cpAbruptChange;
	F32PTR  tinc_cpCI, tdec_cpCI;

	//tallyPosNegOutliers
	F32PTR  opos_ncp, oneg_ncp;
	I32PTR  opos_ncpPr, oneg_ncpPr;
	I32PTR  opos_cpOccPr, oneg_cpOccPr;
	F32PTR  opos_cp, oneg_cp;
	F32PTR  opos_cpPr, oneg_cpPr;
	F32PTR  opos_cpCI, oneg_cpCI;


} BEAST2_RESULT, * _restrict BEAST2_RESULT_PTR;

typedef struct CORESULT {
	I32PTR xNProb, xProb, xorder;
	F32PTR x, xSD;
} CORESULT, * _restrict CORESULT_PTR;

typedef struct BEAST2_RNDSTREAM {
	F32PTR  rndgamma;
	U32PTR  rnd32;
	U16PTR  rnd16;
	U08PTR  rnd08;
} BEAST2_RNDSTREAM, * _restrict BEAST2_RANDSEEDPTR;


struct BEAST2_BASIS;
struct BEAST2_MODEL;
typedef struct BEAST2_BASIS* _restrict BEAST2_BASIS_PTR;
typedef struct BEAST2_MODEL BEAST2_MODEL, * _restrict BEAST2_MODEL_PTR;
typedef struct PROPOSE_STRUCT {
	I32PTR             samples;
	CORESULT_PTR       keyresult;
	F32PTR             mem;
	BEAST2_MODEL_PTR   model;
	BEAST2_RANDSEEDPTR pRND;
	BEAST2_YINFO_PTR   yInfo;
	I32                nSample_ExtremVecNeedUpdate;
	I32                N, Npad16;
	F32                sigFactor;
	F32                outlierSigFactor;
} PROP_DATA, * _restrict PROP_DATA_PTR;



typedef struct BEAST2_BASESEG {
	I32 R1, R2, K;
	//solaris: union {
		//solaris: struct {
	I16 ORDER1, ORDER2;
	//solaris: };
	I32     outlierKnot; //used only for the outlier component
	//solaris:};
} BEAST2_BASESEG, * _restrict BEAST2_BASESEG_PTR;

typedef struct _NEWTERM {

	BEAST2_BASESEG  SEG[2];
	//solaris: union {
	TKNOT          newKnot;
	//solaris: struct { 
	TORDER oldOrder, newOrder;
	//solaris: };
//solaris: };

	I16 nKnot_new;
	I16 newIdx;

	NEWCOLINFO xcols;

	U08 numSeg;
	I08 jumpType;
	//I16 Knewterm; //not used
	//U08 basisID;	//not used
} NEWTERM, * _restrict NEWTERM_PTR;



//period is used as "basis->bConst.dummy.period" in the following functons: computeXY,Alloc_Init_PrecPrior,DD_CalcBasisKsKeK_prec0123
//TERMS needed only if meta->io.deseasonlaized=TRUE
typedef struct DUMMY_CONS { F32PTR TERMS; F32PTR SQRT_N_div_n; I32 period; }          DUMMY_CONST;
typedef struct SVD_CONS { F32PTR TERMS; F32PTR SQR_CSUM; }                            SVD_CONST;
typedef struct SEASON_CONST { F32PTR TERMS, SQR_CSUM, SCALE_FACTOR; }				    SEASON_CONST;
typedef struct TREND_CONS { F32PTR TERMS, COEFF_A, COEFF_B, INV_SQR; }	               TREND_CONST;
typedef struct OUTLIER_CONS { F32    SQRTN, SQRTN_1; }		                           OUTLIER_CONST;

typedef union {
	SVD_CONST     svd;
	DUMMY_CONST   dummy;
	SEASON_CONST  season;
	TREND_CONST   trend;
	OUTLIER_CONST outlier;
} BASIS_CONST;

//https://gcc.gnu.org/onlinedocs/gcc/Unnamed-Fields.html

typedef int  (*pfnGenTerms)(F32PTR X, I32 N, BEAST2_BASESEG*, BASIS_CONST* ptr);
typedef struct PROP_PROB_STRUCT { U08 birth, death, merge, move; } PROP_PROB_STRUCT;
typedef struct BEAST2_BASIS {

	BASIS_CONST  bConst;
	struct {
		void        (*Propose)(BEAST2_BASIS_PTR, NEWTERM_PTR, NEWCOLINFO_PTR, PROP_DATA*);
		pfnGenTerms GenTerms;
		int         (*CalcBasisKsKeK_TermType)(BEAST2_BASIS_PTR  basis);
		void        (*UpdateGoodVec)(BEAST2_BASIS_PTR basis, NEWTERM_PTR new, I32 Npad16_not_used);
		void        (*ComputeY)(F32PTR X, F32PTR beta, F32PTR Y, BEAST2_BASIS_PTR basis, I32 Npad);
		F32(*ModelPrior)(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new);
	};

	F32PTR   scalingFactor;
	F64PTR   priorMat;
	F64PTR   priorVec;

	struct {
		TKNOT  minSepDist;
		I16    minKnotNum;
		I16    maxKnotNum;
		TORDER minOrder, maxOrder;
	} prior;

	PROP_PROB_STRUCT	propprob;
	I16					mcmc_Kstopping;
	I16					mcmc_MoveStep;

	TKNOT_PTR   KNOT;
	TORDER_PTR	ORDER;  // not used for the dummpy basis
	/*
	union {
		TORDER_PTR	ORDER;  // not used for the dummpy basis
		I32			period; // only used for the dummy basiss
	};
	*/
	I16PTR     ks, ke;
	U08PTR     termType;
	U08PTR08   goodvec;

	I16      nPrec;
	I16      offsetPrec;

	I32      goodNum;
	I16      nKnot;
	I16      K, Kbase;
	U08      type;

} BEAST2_BASIS, * _restrict BEAST2_BASIS_PTR;


typedef struct {
	F32PTR XtX, XtY, cholXtX, beta_mean;
	F32PTR precXtXDiag;
	I16PTR nTermsPerPrecGrp;

	F32PTR   alpha2Q_star;  // made to be a pointer for MRBEAST

	F32    marg_lik;
	I32    K;

	/*
	union {
		F32      alpha2_star;
		F32PTR   alphaQ_star; // added for MRBEAST
	};
	*/
} BEAST2_MODELDATA, * _restrict  BEAST2_MODELDATA_PTR;

typedef struct BEAST2_MODEL {
	I32(*PickBasisID)(PROP_DATA_PTR);

	F32PTR  beta;
	F32PTR	sig2; // for MRBEAST

	/////////////////////////////	
	I08PTR08 extremePosVec;
	F32PTR   deviation;
	F32PTR   avgDeviation;  // Changed to a pointer for a consistent API with MRBEAST
	I32      extremPosNum;

	I16     nPrec;
	F32PTR  precVec;
	F32PTR  logPrecVec;

	/////////////////
	/*
	F32PTR XtX,  XtY, cholXtX,  beta_mean,  beta;
	F32PTR precXtXDiag;
	I08PTR nTermsPerPrecGrp;
	*/
	//F32PTR XtX_prop, XtY_prop, cholXtX_prop, beta_mean_prop, beta_prop;
	BEAST2_MODELDATA curr, prop;

	I32          NUMBASIS;
	I08          vid, did, sid, tid, oid;
	BEAST2_BASIS* b;

} BEAST2_MODEL, * _restrict BEAST2_MODEL_PTR;



typedef struct {
	I32              minSepDist;
	BEAST2_YINFO_PTR yInfo;
} KNOT2BINVEC, * _restrict KNOT2BINVEC_PTR;



