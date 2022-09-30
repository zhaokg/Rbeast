#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "abc_001_config.h"
#include "abc_rand.h"
#include "abc_mat.h"
#include "beastv2_func.h"

#include "abc_ide_util.h" //printf

#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define max3(a,b,c)         max( max(a,b), c)
#define max4(a,b,c,d)       max( max3(a,b,c),d)
#define min(a,b)            (((a) < (b)) ? (a) : (b))



static I32 __GetNumElem_of_XnewTerm(BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt,   I32 * MAX_COL) {
	#define MODEL (*model)
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// For some proposal moves, two new segments will be generated..
     #define MAX_NUM_NEW_SEG 2
	I32 MXNCOL_PERSEG1 = MODEL.sid >= 0 ? (MODEL.b[MODEL.sid].prior.maxOrder)*2L : -9999;  	// mem to save one season term	
	I32 MXNCOL_PERSEG2 = MODEL.tid >= 0 ? (MODEL.b[MODEL.tid].prior.maxOrder+1L) : -9999;    // mem to save one trend term	
	I32 MXNCOL_PERSEG3 = MODEL.did >= 0 ? opt->io.meta.period : -9999;                       // mem to save one dummyterm   

	I32 MAXNUMCOL_Xnewterm = max3(MXNCOL_PERSEG1, MXNCOL_PERSEG2, MXNCOL_PERSEG3) * MAX_NUM_NEW_SEG;
	/////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Extra mem for storing BASESEG
	I32 MAX_TOTAL_SEGNUM = 0;
	for (int i = 0; i < MODEL.NUMBASIS; i++)	MAX_TOTAL_SEGNUM += (MODEL.b[i].prior.maxKnotNum + 1);
	I32 MAX_NUMELEM_SEGINFO = MAX_TOTAL_SEGNUM * (sizeof(BEAST2_BASESEG)/4);
	/////////////////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Xnewterm must have a length larger than or equla to 5*N + max(sncp, tncp);
	I32 MAX_MEM_FOR_CHANGEPOINTS = 6 * opt->io.N;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	// The length of the raw data of input time series  which differs from N for irregular TS
	// Xnewterm is used as temp mem for fetch_next_timeSeries)(&yInfo, pixelIndex, MEMBUF, opt);
	I32 Nraw = opt->io.dims[ opt->io.meta.whichDimIsTime-1L];
	/////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////
	// For MRBEAST, Xnewterm is used as a buffer to store the predicted Y curves that are fed
	// into the ComputeConfidenceInterval as new rows of Y values

	//Xnewterm is also needed to store the recoverd original Y data values to be dumped into the output
	//varaible o.data and to compute R2 and RMSE.
	I32 MAX_COLS_YPRED = opt->io.q * MODEL.NUMBASIS;
	/////////////////////////////////////////////////////////////////////////////////////////////////


	I32 Npad = ((opt->io.N + 7) / 8) * 8;  Npad = opt->io.N;//Correct for the inconsitency of X and Y in gemm and gemv
	I32 TOTAL_NUM = max4( MAXNUMCOL_Xnewterm* Npad + MAX_NUMELEM_SEGINFO, 
		                  MAX_MEM_FOR_CHANGEPOINTS,
		                  MAX_COLS_YPRED*Npad,
		                  Nraw);
		
		
	*MAX_COL = MAXNUMCOL_Xnewterm;
	return TOTAL_NUM;

    #undef MODEL
}
/*************************************************/
//Arrays that have a leading dimension of Npad:
// Xt_mars, and Xnewterm
// Y has a leading dimension of N, as in fetech_next_time-series
/*************************************************/

void AllocateXXXMEM( F32PTR * Xt_mars, F32PTR*  Xnewterm, F32PTR*  Xt_zeroBackup,
	                 BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt, MemPointers * MEM) 
{	
#define RoundTo64(N)  ((N+63)/64*64)

	I32 N         = opt->io.N;
	I32 Npad      = ((N + 7) / 8) * 8;
	I32 K_MAX     = opt->prior.K_MAX;
	I64 szXtmars  = Npad * K_MAX;

	I32 XNEW_MAX_NUMCOL;
	I32 XNEW_TOTAL_NUM    = __GetNumElem_of_XnewTerm(model, opt, &XNEW_MAX_NUMCOL);
	I32 MAXNUM_MISSINGROW = N * opt->io.meta.maxMissingRate + 1;

	I64 szXnewterm     =  XNEW_TOTAL_NUM;
	I64 szXtzeroBackup =  MAXNUM_MISSINGROW * XNEW_MAX_NUMCOL;


	I32 szTotal = RoundTo64(szXtmars) + RoundTo64(szXnewterm) + RoundTo64(szXtzeroBackup);

	*Xt_mars       = MyALLOC(*MEM, szTotal, F32, 64);
	*Xnewterm      = *Xt_mars  + RoundTo64(szXtmars) ;
	*Xt_zeroBackup = *Xnewterm + RoundTo64(szXnewterm);
	//Xt_mars		= MyALLOC(MEM, Npad * MAX_K, F32, 64); 
	//Xnewterm		= MyALLOC(MEM, TOTAL_NUM, F32, 64);
	//Xt_zeroBackup = MyALLOC(MEM, MAXNUM_MISSINGROW *MAX_NUMCOL, F32, 64);
}


void AllocateYinfoMEM(BEAST2_YINFO_PTR yInfo, BEAST2_OPTIONS_PTR opt, MemPointers* MEM) {
	I32  N    = opt->io.N;	//I32  Npad = ((N+7)/8) * 8;
	I32  q    = opt->io.q;  //q=1 for BEASTV2
	I32  qxq = q*q;
	//Important: In fetch_next_timeSeries (__NormalizeF32ArrayWithNaNOmitted_zerobased), it is possible
	// that all the points are bad indices, so the space for rowMissing should be at least of length N.
	I32  MAXNUM_MISSINGROW = N;  //N * opt->io.meta.maxMissingRate + 1L;

	// Sizes: Y[Npad], sd[q], mean[q], YtY_plus_alpha2Q[q*q]
	yInfo->Y                = MyALLOC(*MEM, N*q+ qxq +q +q,              F32, 64);
	yInfo->YtY_plus_alpha2Q = yInfo->Y                + N*q;
	yInfo->mean             = yInfo->YtY_plus_alpha2Q + qxq;
	yInfo->sd               = yInfo->mean + q;

	// yInfo.rowsMissing: THE INDEX OF ROWS WITH MISSING DATA	
	yInfo->rowsMissing     = MyALLOC(*MEM, MAXNUM_MISSINGROW, U32, 64);	 	

	//added for MRBEAST and used in ComputeLik, PropseNew/__CalcAbsDeviation(compute deviaiton and extrempos)
	//, MR_EvaluateModel,
	yInfo->q               = q;

	if (opt->io.meta.deseasonalize) {
		yInfo->Yseason= MyALLOC(*MEM, N * q , F32, 0);
	} else {
		yInfo->Yseason = NULL;
	}

	if (opt->io.meta.detrend) {
		yInfo->Ytrend = MyALLOC(*MEM, N * q, F32, 0);
	}
	else {
		yInfo->Ytrend = NULL;
	}
}

#include "abc_000_warning.h"