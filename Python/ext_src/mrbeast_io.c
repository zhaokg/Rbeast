#include "abc_001_config.h"

DISABLE_MANY_WARNINGS

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "abc_datatype.h"

#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_vec.h"  // for f32_seq only

#include "mrbeast_header.h"
#include "mrbeast_io.h"
#include "mrbeast_func.h"   


#define BEGIN {
#define END   }
int  MV_Get1stArg_Data(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{
	// Check the number of inputs
	if (nrhs == 1L + 0) { // the first arg is the algorithm name.
		r_error("ERROR: There should be at least one input argument!.\n");
		return 0;
	}

	if (nrhs >= 1L + 1) 
    BEGIN
		// Check if the input 'data' is numeric and lengthy enough
		VOID_PTR DATA   = prhs[1L]; 
	    I32      q   = GetNumberOfElements(DATA);
		if (q < 2) {
			r_error("ERROR: There should be at least two time series variables.\n");
			return 0;
		}
		if (!IsStruct(DATA) && !IsCell(DATA)) {
			r_error("ERROR: The input data should be a LIST variable containing multiple elements!\n");
			return 0;
		}


		MV_IO_PTR io = opt->io;
		VOID_PTR  Y;

		io->pdata = malloc(q*sizeof(VOID_PTR));

		/*************************************************/
		/* If y is not a string but a matrix or vector  */
		/* find the dimesions, and configure io-*/
		/*************************************************/

		Y = GetFieldByIdx(DATA, 0);
		I32  numel = GetNumberOfElements(Y);
		if (!(IsDouble(Y) && numel > 2) && !(IsSingle(Y) && numel > 2) &&
			!(IsInt32(Y) && numel > 2))    { // if (!(TYPEOF(Y) == INTSXP  && XLENGTH(Y) > 2)  &&  !(TYPEOF(Y) == REALSXP && XLENGTH(Y) > 2) && !(TYPEOF(Y) == INTSXP  && isMatrix(Y))     &&  !(TYPEOF(Y) == REALSXP && isMatrix(Y)) && 	!(TYPEOF(Y) == STRSXP  && XLENGTH(Y) == 1) &&  !(TYPEOF(Y) == REALSXP  && isArray(Y)) 		)	
			r_error("ERROR: The input data should be a numeric type and must be long enough.\n");
			return 0;
		}

		if (IsInt32(Y))
			io->dtype = DATA_INT32;
		else if (IsDouble(Y))  // isReal(pY)
			io->dtype = DATA_DOUBLE;
		else if (IsSingle(Y))  // isReal(pY)	
			io->dtype = DATA_FLOAT;
		else { // io->dtype = DATA_UNKNOWN;
			r_error("ERROR: The input data has an uknown numeric type!\n");
			return 0;
		}

		io->pdata[0] = GetData(Y);

		I32 ndims = GetNumOfDim(Y);
		if (ndims == 0) {
			// If the input is a vector: this branch is possible only for R vectors (not  Matlab)
			I32 N = GetNumberOfElements(Y);
			io->numOfPixels = 1L;
			io->ndim = 2L;
			io->dims[0] = N;
			io->dims[1] = 1L;
		}
		else if (ndims == 2){ //ndims is impossible to be 1L
			// Matlab: a vector or matrix; R: a matrix or a vector with a dim attribute.

			int N = GetDim1(Y);
			int M = GetDim2(Y);
			if (min(N, M) == 1L)  //PY is a vector
				N = max(N, M),
				io->numOfPixels = 1;
			else                 //PY is a matrix	
				N = N,
				io->numOfPixels = M;

			io->ndim = 2L;
			io->dims[0] = N;
			io->dims[1] = io->numOfPixels;

		}
		else if (ndims == 3){
			// If the input is a 3D array
			io->ndim = 3L;
			GetDimensions(Y, io->dims, 3L);
			//numOfPixels needs to be assigned later in AllocateOutput
		}
		else {
			r_error("ERROR: The max dim allowed is 3 for a 3D image stack, but the input data cannot have a dimension"
				" of %d .\n", ndims);
			return 0;
		}

		//TODO: here the data type an dimensions of all elments are assumed to be the same
		// A sanity check needs to be done to ensure that is the case.
		for (rI32 i = 1; i < q; i++) {			
			Y = GetFieldByIdx(DATA, i);
			io->pdata[i] = GetData(Y);
		}

		// Added for MRbeast
		opt->q = q;
	END

	return 1L;
}

static int  Get2ndArg_MetaData_Regular_Dim2(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{
	/******************************************/
	/*  If the input is regular time series */
	/******************************************/
	MV_IO_PTR       io   = opt->io;
	MV_METADATA_PTR meta = &(opt->io->meta);

	meta->period = -1;
	if (nrhs < 3) {
		meta->deltaTime = 1.f;
		meta->startTime = 1.f;
		meta->period    =  -1;
		meta->missingValue = FLOAT_MAX;
		meta->maxMissingRate = 0.75;
		meta->whichDimIsTime        = 9999; //Not used
	}
	if (nrhs >= 3) {
		VOID_PTR pmeta = prhs[2L];

		if (IsNumeric(pmeta)) 	{
			meta->period = GetScalar(pmeta);
			meta->deltaTime = 1.f;
			meta->startTime = 1.f;
			meta->missingValue = FLOAT_MAX;
			meta->maxMissingRate = 0.75;
			meta->whichDimIsTime = 255; //Not used
		} else if (IsStruct(pmeta)) {
			VOID_PTR tmp;
			meta->missingValue = (tmp = GetField(pmeta, "missingValue")) ? GetScalar(tmp) : FLOAT_MAX;
			meta->period        = (tmp = GetField(pmeta, "period")) ? GetScalar(tmp) : -1;
			meta->deltaTime     = (tmp = GetField(pmeta, "deltaTime")) ? GetScalar(tmp) : 1;
			meta->startTime     = (tmp = GetField(pmeta, "startTime")) ? GetScalar(tmp) : 1;
			meta->maxMissingRate = (tmp = GetField(pmeta, "maxMissingRate")) ? GetScalar(tmp) : 0.75;
			meta->whichDimIsTime        = (tmp = GetField(pmeta, "whichDimIsTime")) ? GetScalar(tmp) : 255;
		} else{
			r_error("ERROR: The 'metadata' parameter given is of unsupported type.\n");
			return 0;
		}
	} //else if (nrhs >= 3)

	int    N = io->dims[0];

	// If period is missing or negative, estimate its value
	if (meta->period <= 0) {

		F32PTR yData = malloc(sizeof(F32)*N);	

		CopyStrideMEMToF32Arr(yData,io->pdata[0],N, 1L/*stride*/, 0L/*offset*/, io->dtype);
		f32_set_nan_by_value(yData, N, meta->missingValue);
		meta->period = DeterminePeriod(yData, N);

		free(yData);

		r_printf("WARNING: The \"metadata$period\" parameter should be known in advance and supplied by the user but it is missing. A best guess of it is "
			" %5d AND will be used in the decomposition. Please make sure this estimate makes sense; otherwise, the BEAST decomposition result will be "
			"incorrect.\n", (int) meta->period  );
		
	} //if (meta->period <= 0 && io->ndim ==2)

	if (meta->period <= 0 || fabsf(meta->period - ceil(meta->period)) > 1e-4f) {
		r_error("The value of period must be an positive integer!");
		return 0;
	}

	opt->N = N;
	opt->Npad = (I32)ceil((F32)N / 8.0f) * 8;
	opt->Npad16 = (I32)ceil((F32)N / 16.f) * 16;

	return 1;
}
static int  Get2ndArg_MetaData_Regular_Dim3(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{
	/******************************************/
	/*  If the input is regular time series */
	/******************************************/
	MV_IO_PTR       io   = opt->io;
	MV_METADATA_PTR meta = &(opt->io->meta);
	
	if (nrhs < 3) {
		r_error("ERROR: For a 3D stack of time series, \"metadata$whichDimIsTime\" must be given "
			    "to specify which dim of the array refers to time.\n");
		return 0;
	}

	VOID_PTR pmeta = prhs[2L];
	if (!IsStruct(pmeta)) {
		r_error("ERROR: For a 3D stack input, the 'metadata' argument must be a LIST variable describing "
			    " the data strucutre, such as metedata$period and metadata$whichDimIsTime.\n");
		return 0;
	}
	if (GetField(pmeta, "period") == NULL) {
		r_error("ERROR: For a 3D stack input of regular time series, metadata$period must be given.\n");
		return 0;
	}
	if (GetField(pmeta, "whichDimIsTime") == NULL) {
		r_error("ERROR: For a 3D stack input of regular time series, metadata$whichDimIsTime must be given.\n");
		return 0;
	}
 
	VOID_PTR tmp;
	meta->period         = (tmp = GetField(pmeta, "period")) ? GetScalar(tmp) : -1;
	meta->whichDimIsTime = (tmp = GetField(pmeta, "whichDimIsTime")) ? GetScalar(tmp) : 255;
	meta->missingValue  = (tmp = GetField(pmeta, "missingValue")) ? GetScalar(tmp) : FLOAT_MAX;
	
	meta->deltaTime      = (tmp = GetField(pmeta, "deltaTime")) ? GetScalar(tmp) : 1;
	meta->startTime      = (tmp = GetField(pmeta, "startTime")) ? GetScalar(tmp) : 1;
	meta->maxMissingRate = (tmp = GetField(pmeta, "maxMissingRate")) ? GetScalar(tmp) : 0.75;
	
	if (meta->period  <= 0) {
		r_error("ERROR: metadata$period should be postive!\n");
		return 0;
	}
	if (meta->whichDimIsTime != 1 && meta->whichDimIsTime != 2 && meta->whichDimIsTime != 3) {
		r_error("ERROR: metadata$whichDimIsTime should take a value out of 1, 2, and 3 only.");
		return 0;
	}

	int N       = io->dims[meta->whichDimIsTime-1];
	opt->N      = N;
	opt->Npad   = ((N + 8L - 1L) / 8L) * 8L;
	opt->Npad16 = (I32)ceil((F32)N / 16.f) * 16;
	return 1L;	
}

void MV_DeallocateTimeSeriesIO(MV_IO*o) {
	
 
		if (o->T.numPtsPerInterval != NULL) {
			free(o->T.numPtsPerInterval);
			o->T.numPtsPerInterval = NULL;
		}
		if (o->T.sortedTimeIdx != NULL) {
			free(o->T.sortedTimeIdx);
			o->T.sortedTimeIdx = NULL;
		}
		if (o->pdata != NULL) {
			free(o->pdata);
			o->pdata = NULL;
		}
 
}


static I32 __ReadRawTime(F32PTR outtime, VOID_PTR timeField, const I32 N)
{
	if (!IsStruct(timeField)) {
		if (CopyNumericArrToF32Arr(outtime, timeField, N))
			return 1L;
		else {
			r_printf("ERROR: metadata$time has an unsupported data format or type.\n");
			return 0;
		}
	}

	return 0;

}

static Get2ndArg_MetaData_Irregular(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{
	/******************************************/
	/*  Check if the input is IRREGULAR or not   */
	/******************************************/
	/*
	U08 isRegularOrdered=1;

	if (nrhs < 1 + 2L) //algname, data: no metadata supplied 
		isRegularOrdered = 1;
	else { // metadata is supplied 
		VOID_PTR metadata = prhs[2L];
		if (IsNumeric(metadata) && GetScalar(metadata) > 1L)
			isRegularOrdered = 1;
		else if (IsStruct(metadata)) {
			VOID_PTR tmp = GetField(metadata, "isRegularOrdered");
			if (tmp != NULL) {
				if (IsNumeric(tmp))
					isRegularOrdered = GetScalar(tmp);
				else {
					r_error("ERROR: metadata$isRegularOrdered must be numeric.\n");
					return 0;
				}
			} //if (tmp != NULL) 
			else {
				isRegularOrdered = 1;
			}

		}// else if (IsStruct(metadata))
	}//else {   metadata is supplied 
	*/

	/******************************************/
	/* If the input time series data are irregular
	/******************************************/
	MV_IO_PTR       io   = opt->io;
	MV_METADATA_PTR meta = &(opt->io->meta);
	
	VOID_PTR pmeta = prhs[2L];
	if ( nrhs < 3 || !IsStruct(pmeta) ) {
		r_printf("ERROR: For irregular time series, the argument 'metadata' must be supplied and it should be a LIST variable!");
		return 0;
	}
	if ( GetField(pmeta, "timeVec") == NULL ) {
		r_printf("ERROR: For irregular time series, metadata$time must be given to specify the times at which"
			     " data points are collected.");
		return 0;
	}
	
	VOID_PTR tmp;

	meta->whichDimIsTime = 1;
	if (io->ndim == 2)
	{		
		if (GetField(pmeta, "timeDimensionIndex") != NULL) 
			r_printf("WARNING: The input data is either a 1D vector or 2D matrix, for which "
				" metadata$timeDimensionIndex is not used and ignored.\n");			
	} else { //(io->ndim == 3)
		tmp = GetField(pmeta, "timeDimensionIndex");
		if (tmp == NULL) {
			r_printf("ERROR: For a 3D array input, metadata$timeDimensionIndex must be specified to "
				" tell which dim of the array refers to time.\n");
			return 0;
		}		
		meta->whichDimIsTime = GetScalar(tmp);
		if (meta->whichDimIsTime != 1 && meta->whichDimIsTime != 2 && meta->whichDimIsTime != 3) {
			r_error("ERROR: metadata$whichDimIsTime should take a value out of 1, 2, and 3 only.\n");
			return 0;
		}
	}

	int    N    = io->dims[meta->whichDimIsTime - 1];
	void * TIME = GetField(pmeta, "timeVec");

	if ( GetNumberOfElements(TIME)!= N ) {
		r_printf("ERROR: metadata$time must be of the same length as the time dim of the input data.\n");
		return 0;
	}

	tmp = GetField(pmeta, "deltaTime");
	if (tmp == NULL || GetScalar(tmp) <= 0) {
		r_printf("ERROR: metadata$timeInterval must be provided and it should be positive!\n"); 		
		return 0;
	}
	
	F32 dT = GetScalar(tmp);
	meta->deltaTime = dT;

	tmp = GetField(pmeta, "period_in_time");
	if (tmp == NULL || GetScalar(tmp) <= 0) {
		r_printf("ERROR: metadata$period_in_time must be provided and it should be positive!\n"); 	
		return 0;
	}
	F32 period = GetScalar(tmp);

	F32 n = (period / dT);
	if (abs(n - round(n)) > 1e-3) 	{
		r_printf("ERROR: metadata$period_in_time must be divided by metadata$deltaTime for an integer"
				" number of times.\n");		
		return 0;
	}

	tmp = GetField(pmeta, "period");
	if (tmp != NULL && GetScalar(tmp) != n)    {
		r_printf("WARNING:metadata$period_in_time is not equal to metadata$period * metadata$deltaTime, "
				" the specified value of metadata$period is ignored, and a value of %d is used instead.\n"
				, (I32)n);
	}
	meta->period = n;

	//meta->startTime = BuildInfoToAggegrateTS(&(io->T), TIME, dT, period, &N);	
	{
		I32    Nold = GetNumberOfElements(TIME);
		F32PTR time = malloc(sizeof(F32)*Nold);

		__ReadRawTime(time, TIME, Nold);
		N=tsAggegrationPrepare(
			time, Nold, dT,
			&io->T.sortedTimeIdx, &io->T.numPtsPerInterval, &io->T.startIdxOfFirsInterval, &meta->startTime);

		free(time);
	}
	
	//Npad:N_extended; the mutiples of 8 closest to N, defined for 32-byte alginment
	opt->N      = N;
	opt->Npad   = (I32)ceil((F32)N / 8.0f) * 8;
	opt->Npad16 = (I32)ceil((F32)N / 16.f) * 16;

	tmp = GetField(pmeta, "missingValue");
	meta->missingValue = tmp != NULL ? GetScalar(tmp) : FLOAT_MAX;
	meta->maxMissingRate = (tmp = GetField(pmeta, "maxMissingRate")) ? GetScalar(tmp) : 0.75;
	return 1L;	
}
int  MV_Get2ndArg_MetaData(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{
	MV_IO_PTR io = opt->io;
	if (io->isRegularOrdered)
	{
		if (io->ndim==2)
			return Get2ndArg_MetaData_Regular_Dim2(prhs, nrhs, opt);
		else
			return Get2ndArg_MetaData_Regular_Dim3(prhs, nrhs, opt);
	}
	else 
	{
		return Get2ndArg_MetaData_Irregular(prhs, nrhs, opt);	
	}
	

}
int  MV_Get3rdArg_Prior(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{		
	struct PRIOR_MISSING {
		U08   minSeasonOrder, maxSeasonOrder, minTrendOrder, maxTrendOrder;
		U08   trendMinSepDist, seasonMinSepDist;
		U08   trendMaxKnotNum, seasonMaxKnotNum;
	} m;
	memset(&m, 0L, sizeof(struct PRIOR_MISSING));

#define o (*prior)
	MV_PRIOR_PTR prior = &(opt->prior);
	VOID_PTR tmp, S;
	S = prhs[3L];
	if (nrhs < 4)
	{
		memset(&m, 1L, sizeof(struct PRIOR_MISSING));
	}
	if (nrhs >= 4)
	{
		if (!IsStruct(S))
		{
			r_printf("WARNING: The arg 'prior' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct PRIOR_MISSING));
		} else {

			o.minSeasonOrder = (tmp = GetField(S, "minSeasonOrder")) ? GetScalar(tmp) : (m.minSeasonOrder = 1);
			o.maxSeasonOrder = (tmp = GetField(S, "maxSeasonOrder")) ? GetScalar(tmp) : (m.maxSeasonOrder = 1);
			o.minTrendOrder = (tmp = GetField(S, "minTrendOrder")) ? GetScalar(tmp) : (m.minTrendOrder = 1);
			o.maxTrendOrder = (tmp = GetField(S, "maxTrendOrder")) ? GetScalar(tmp) : (m.maxTrendOrder = 1);
			o.trendMinSepDist = (tmp = GetField(S, "trendMinSepDist")) ? GetScalar(tmp) : (m.trendMinSepDist = 1);
			o.seasonMinSepDist = (tmp = GetField(S, "seasonMinSepDist")) ? GetScalar(tmp) : (m.seasonMinSepDist = 1);
			o.trendMaxKnotNum = (tmp = GetField(S, "trendMaxKnotNum")) ? GetScalar(tmp) : (m.trendMaxKnotNum = 1);
			o.seasonMaxKnotNum = (tmp = GetField(S, "seasonMaxKnotNum")) ? GetScalar(tmp) : (m.seasonMaxKnotNum = 1);
		}
		
	} // if (nrhs >= 4)

	F32 period = opt->io->meta.period;
	F32 N      = opt->N;
	if (m.minSeasonOrder)    o.minSeasonOrder = 1L;				    o.minSeasonOrder = min(o.minSeasonOrder, period / 2 - 1);    o.minSeasonOrder = max(o.minSeasonOrder, 1L);
	if (m.maxSeasonOrder)    o.maxSeasonOrder = (period / 2 - 1); o.maxSeasonOrder = min(o.maxSeasonOrder, (period / 2 - 1));  o.maxSeasonOrder = max(o.maxSeasonOrder, o.minSeasonOrder);
	if (m.minTrendOrder)     o.minTrendOrder = 0L;				    o.minTrendOrder = max(o.minTrendOrder, 0L);
	if (m.maxTrendOrder)     o.maxTrendOrder = 1L;				    o.maxTrendOrder = max(o.maxTrendOrder, o.minTrendOrder);
	if (m.trendMinSepDist)  o.trendMinSepDist  = period / 2;	    o.trendMinSepDist = max(o.trendMinSepDist, 0L);   o.trendMinSepDist = min(o.trendMinSepDist, N / 2 - 1);
	if (m.seasonMinSepDist) o.seasonMinSepDist = period / 2;    o.seasonMinSepDist = max(o.seasonMinSepDist, 2L); o.seasonMinSepDist = min(o.seasonMinSepDist, N / 2 - 1);

	if (m.trendMaxKnotNum)  o.trendMaxKnotNum = floor(N / (o.trendMinSepDist + 1) - 1);	o.trendMaxKnotNum = min(o.trendMaxKnotNum, floor(N / (o.trendMinSepDist + 1) - 1));
	if (m.seasonMaxKnotNum) o.seasonMaxKnotNum = floor(N / (o.seasonMinSepDist + 1) - 1);	o.seasonMaxKnotNum = min(o.seasonMaxKnotNum, floor(N / (o.seasonMinSepDist + 1) - 1));

#undef o


	return 1;
	
}
int  MV_Get4thArg_MCMC(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{		
	struct MCMC_MISSING {
		U08   seed;	                  // Unsigned long long seed;
		U08   credIntervalAlphaLevel;
		U08   trendResamplingOrderProb;
		U08   seasonResamplingOrderProb;
		U08   ridgeFactor;		
		U08   burnin, samples, chainNumber;
		U08   maxMoveStepSize;
		U08   thinningFactor;		
	} m;
	memset(&m, 0L, sizeof(struct MCMC_MISSING));

#define o (*mcmc)
	MV_MCMC_PTR mcmc = &(opt->mcmc);
	VOID_PTR tmp, S;
	S = prhs[4L];

	if (nrhs < 5) {
		memset(&m, 1L, sizeof(struct MCMC_MISSING));
	}

	if (nrhs >= 5)	{
		if (!IsStruct(S))
		{
			r_printf("WARNING: The arg 'mcmc' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct MCMC_MISSING));
		} else {

			o.maxMoveStepSize   = (tmp = GetField(S, "maxMoveStepSize"))   ? GetScalar(tmp) : (m.maxMoveStepSize = 1);
			o.samples           = (tmp = GetField(S, "samples"))           ? GetScalar(tmp) : (m.samples = 1);
			o.thinningFactor    = (tmp = GetField(S, "thinningFactor"))    ? GetScalar(tmp) : (m.thinningFactor = 1);
			o.burnin            = (tmp = GetField(S, "burnin"))            ? GetScalar(tmp) : (m.burnin = 1);
			o.chainNumber       = (tmp = GetField(S, "chainNumber"))       ? GetScalar(tmp) : (m.chainNumber = 1);
			o.trendResamplingOrderProb  = (tmp = GetField(S, "trendResamplingOrderProb"))  ? GetScalar(tmp) : (m.trendResamplingOrderProb = 1);
			o.seasonResamplingOrderProb = (tmp = GetField(S, "seasonResamplingOrderProb")) ? GetScalar(tmp) : (m.seasonResamplingOrderProb = 1);
			o.seed                  = (tmp = GetField(S, "seed"))              ? GetScalar(tmp) : (m.seed = 1);
			o.credIntervalAlphaLevel  = (tmp = GetField(S, "credIntervalAlphaLevel")) ? GetScalar(tmp) : (m.credIntervalAlphaLevel = 1);
			o.ridgeFactor = (tmp = GetField(S, "ridgeFactor")) ? GetScalar(tmp) : (m.ridgeFactor = 1);

			}
		
	} // if (nrhs >= 5)

	int trendMinSepDist = opt->prior.trendMinSepDist;

	if (m.maxMoveStepSize) o.maxMoveStepSize = (trendMinSepDist / 3 + 1);
	if (m.samples)         o.samples = 3000;          o.samples = max(o.samples, 800);
	if (m.thinningFactor)  o.thinningFactor = 1L;     o.thinningFactor = max(o.thinningFactor, 1L);
	if (m.burnin)          o.burnin = 200L;			  o.burnin = max(o.burnin, 200L);
	if (m.chainNumber)     o.chainNumber = 3;		  o.chainNumber = max(o.chainNumber, 1L);
	if (m.trendResamplingOrderProb)  o.trendResamplingOrderProb = .1f;
	if (m.seasonResamplingOrderProb) o.seasonResamplingOrderProb = .17f;
	if (m.seed)                  o.seed = 0L;
	if (m.credIntervalAlphaLevel)            o.credIntervalAlphaLevel = .95;
	if (m.ridgeFactor)           o.ridgeFactor = 0.0001f;

#undef o
	return 1;
	
}
int  MV_Get5thArg_FLAGS(VOID_PTR prhs[], int nrhs, MV_OPTIONS_PTR opt)
{		
	struct OUTFLAGS_MISSING {
		I08  whichOutputDimIsTime;
		I08  numCPUCoresToUse;
		I08  consoleWidth;
		I08  computeSlopeSign;
		I08  computeHarmonicOrder;
		I08  computeTrendOrder;
		I08  computeChangepoints;
		I08  computeCredible;
		I08  fastCIComputation;
		I08  printOptions;
		I08  printProgressBar;
	} m;
	memset(&m, 0L, sizeof(struct OUTFLAGS_MISSING));

#define o (*extra)
	MV_EXTRA_PTR extra = &(opt->extra);
	VOID_PTR tmp, S;
	S = prhs[5L];

	if (nrhs < 6) {
		memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));
	}

	if (nrhs >= 6)	{
		if (!IsStruct(S)) 	{
			r_printf("WARNING: The arg 'extra' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));
		} else {

			o.whichOutputDimIsTime = (tmp = GetField(S, "whichOutputDimIsTime")) ? GetScalar(tmp) : (m.numCPUCoresToUse = 1);
			o.numCPUCoresToUse     = (tmp = GetField(S, "numCPUCoresToUse")) ? GetScalar(tmp) : (m.numCPUCoresToUse = 1);
			o.consoleWidth         = (tmp = GetField(S, "consoleWidth")) ? GetScalar(tmp) : (m.consoleWidth = 1);

			o.printProgressBar     = (tmp = GetField(S, "printProgressBar")) ? GetScalar(tmp) : (m.printProgressBar = 1);
			o.printOptions         = (tmp = GetField(S, "printOptions")) ? GetScalar(tmp) : (m.printOptions = 1);

			o.computeCredible  = (tmp = GetField(S, "computeCredible")) ? GetScalar(tmp) : (m.computeCredible = 1);
			o.fastCIComputation = (tmp = GetField(S, "fastCIComputation")) ? GetScalar(tmp) : (m.fastCIComputation = 1);
			o.computeSlopeSign = (tmp = GetField(S, "computeSlopeSign")) ? GetScalar(tmp) : (m.computeSlopeSign = 1);
			o.computeHarmonicOrder = (tmp = GetField(S, "computeHarmonicOrder")) ? GetScalar(tmp) : (m.computeHarmonicOrder = 1);
			o.computeTrendOrder = (tmp = GetField(S, "computeTrendOrder")) ? GetScalar(tmp) : (m.computeTrendOrder = 1);
			o.computeChangepoints = (tmp = GetField(S, "computeChangepoints")) ? GetScalar(tmp) : (m.computeChangepoints = 1);

			}
		
	} // if (nrhs >= 5)

 
	if (m.numCPUCoresToUse) o.numCPUCoresToUse = 0;	 

	if (m.whichOutputDimIsTime)  o.whichOutputDimIsTime = opt->io->meta.whichDimIsTime;

	if (m.consoleWidth)          o.consoleWidth = GetConsoleWidth();
	if (m.printProgressBar)      o.printProgressBar = 0L;
	if (m.computeCredible)       o.computeCredible = 0L;
	if (m.fastCIComputation)     o.fastCIComputation = 1L;

	if (m.computeSlopeSign)      o.computeSlopeSign = 0L;
	if (m.computeHarmonicOrder)  o.computeHarmonicOrder = 0L;
	if (m.computeTrendOrder)     o.computeTrendOrder = 0L;
	if (m.computeChangepoints)   o.computeChangepoints = 1L;
	if (m.printOptions)          o.printOptions = 1; 
#undef o
	return 1;
	
}

void * MV_AllocateOutput(MV_OPTIONS_PTR  opt)
{

#if   R_INTERFACE==1
	  DATA_TYPE dtype = DATA_DOUBLE;
#elif M_INTERFACE==1
	  DATA_TYPE dtype = DATA_FLOAT;
#elif P_INTERFACE==1
	  DATA_TYPE dtype = DATA_FLOAT;
#endif

	const MV_IO_PTR      io  = opt->io;
	const MV_RESULT_PTR  mat = io->out.result;

	/*************************************************************/
	if (io->ndim == 3) {
		io->numOfPixels = (I64)io->dims[0] * io->dims[1] * io->dims[2]/\
			io->dims[io->meta.whichDimIsTime - 1L];
	}		
	io->out.dtype       = dtype;                         // DATA_FLOAT or DATA_DOUBLE 
	io->out.whichDimIsTime = opt->extra.whichOutputDimIsTime;		
	/*************************************************************/

	MV_PRIOR_PTR prior = &(opt->prior);

	const int   N = opt->N;
	const int   M = io->numOfPixels;
	const int   q = opt->q;
	const int   mxKnotNumSeason = prior->seasonMaxKnotNum;
	const int   mxKnotNumTrend  = prior->trendMaxKnotNum;

	FIELD_ITEM  fieldList[35];
	I32         nfields = 0;	
	if (io->ndim < 3)  {  // if the input is not a 3D stack
		FIELD_ITEM fldList[] = {
			{ "time",    dtype, 2, { N, 1, 0L, 0L }, &mat->time },
			{ "sN",      dtype, 2, { 1, M, },        &mat->sN },
			{ "tN",      dtype, 2, { 1, M, },        &mat->tN },
			{ "sNProb",  dtype, 2, { mxKnotNumSeason + 1, M, }, &mat->sNProb },
			{ "tNProb",  dtype, 2, { mxKnotNumTrend + 1, M, },  &mat->tNProb },
			{ "sProb",   dtype, 2, { N, M, },        &mat->sProb },
			{ "tProb",   dtype, 2, { N, M, },        &mat->tProb },

			{ "s",       dtype, 3, { N, q, M, },        &mat->s },    //MrBeast
			{ "sCI",     dtype, 4, { N, 2, q, M, },     &mat->sCI },  //MrBeast
			{ "sSD",     dtype, 3, { N, q, M, },        &mat->sSD },  //MrBeast

			{ "t",       dtype, 3, { N, q, M, },        &mat->t },   //MrBeast
			{ "tCI",     dtype, 4, { N, 2, q, M, },     &mat->tCI }, //MrBeast
			{ "tSD",     dtype, 3, { N, q, M, },        &mat->tSD }, //MrBeast

			{ "b",       dtype, 3, { N, q, M, },        &mat->b },   //MrBeast
			{ "bCI",     dtype, 4, { N, 2, q, M, },     &mat->bCI }, //MrBeast
			{ "bSD",     dtype, 3, { N, q, M, },        &mat->bSD }, //MrBeast

			{ "marg_lik",dtype, 2, { 1, M, },        &mat->marg_lik },
			{ "sig2",    dtype, 3, { q, q, M, },        &mat->sig2 }, //MrBeast

			{ "bsign",   dtype, 3, { N, q, M, },        &mat->bsign }, //MrBeast
			{ "horder",  dtype, 2, { N, M, },        &mat->horder },
			{ "torder",  dtype, 2, { N, M, },        &mat->torder },

			{ "scp",             dtype, 2, { mxKnotNumSeason, M, },    &mat->scp },
			{ "scpProb",         dtype, 2, { mxKnotNumSeason, M, },    &mat->scpProb },
			{ "scpAbruptChange", dtype, 2, { mxKnotNumSeason, M, },    &mat->scpAbruptChange },
			{ "scpCI",           dtype, 3, { mxKnotNumSeason, 2, M, }, &mat->scpCI },
			{ "tcp",             dtype, 2, { mxKnotNumTrend, M, },     &mat->tcp },
			{ "tcpProb",         dtype, 2, { mxKnotNumTrend, M, },     &mat->tcpProb },
			{ "tcpAbruptChange", dtype, 2, { mxKnotNumTrend, M, },     &mat->tcpAbruptChange },
			{ "tcpCI",           dtype, 3, { mxKnotNumTrend, 2, M, },  &mat->tcpCI }
		};
		nfields = sizeof(fldList) / sizeof(FIELD_ITEM);
		memcpy(fieldList, fldList, nfields *sizeof(FIELD_ITEM));
	}


	if (io->ndim == 3) 	
	BEGIN  // if the input is a 3D stack
		int   ROW, COL;
 
		switch (io->meta.whichDimIsTime) {
			case 1: 
				ROW = io->dims[1], COL = io->dims[2]; break;
			case 2:
				ROW = io->dims[0], COL = io->dims[2]; break;
			case 3:
				ROW = io->dims[0], COL = io->dims[1]; break;
		}

		if (io->out.whichDimIsTime == 1) {
			FIELD_ITEM fldList[] =	{
				{ "time",			dtype, 2, {N, 1L,0L,0L },                 &mat->time },
				{ "sN",				dtype, 2, {ROW, COL,       },              &mat->sN },
				{ "tN",				dtype, 2, {ROW, COL,       },              &mat->tN },
				{ "sNProb",			dtype, 3, { mxKnotNumSeason + 1, ROW, COL}, &mat->sNProb },
				{ "tNProb",			dtype, 3, { mxKnotNumTrend  + 1, ROW, COL}, &mat->tNProb },
				{ "sProb",			dtype, 3, { N, ROW, COL,    },              &mat->sProb },
				{ "tProb",			dtype, 3, { N, ROW, COL,    },              &mat->tProb },

				{ "s",				dtype, 4, { N, q, ROW, COL,    },              &mat->s },  //MrBeast
				{ "sCI",			dtype, 5, { N, 2, q, ROW, COL  },              &mat->sCI },
				{ "sSD",			dtype, 4, { N, q, ROW, COL     },              &mat->sSD },

				{ "t",				dtype, 4, { N, q, ROW, COL     },              &mat->t },
				{ "tCI",			dtype, 5, { N, 2, q, ROW, COL   },             &mat->tCI },
				{ "tSD",			dtype, 4, { N, q, ROW, COL   },                &mat->tSD },

				{ "b",				dtype, 4, { N, q, ROW, COL       },            &mat->b },
				{ "bCI",			dtype, 5, { N, 2, q, ROW, COL    },            &mat->bCI },
				{ "bSD",			dtype, 4, { N, q, ROW, COL  },                 &mat->bSD },

				{ "marg_lik",		dtype, 2, { ROW, COL       },               &mat->marg_lik },
				{ "sig2",			dtype, 4, { q,q, ROW, COL       },          &mat->sig2 },
		
				{ "bsign",			dtype, 4, { N, q, ROW, COL       },            &mat->bsign },
				{ "horder",			dtype, 3, { N, ROW, COL      },             &mat->horder },
				{ "torder",			dtype, 3, { N, ROW, COL      },             &mat->torder },

				{ "scp",            dtype, 3, { mxKnotNumSeason, ROW, COL},     &mat->scp },
				{ "scpProb",        dtype, 3, { mxKnotNumSeason, ROW, COL},     &mat->scpProb },
				{ "scpAbruptChange",dtype, 3, { mxKnotNumSeason, ROW, COL},     &mat->scpAbruptChange },
				{ "scpCI",          dtype, 4, { mxKnotNumSeason,  2, ROW, COL}, &mat->scpCI },

				{ "tcp",            dtype, 3, { mxKnotNumTrend, ROW, COL},      &mat->tcp },
				{ "tcpProb",        dtype, 3, { mxKnotNumTrend, ROW, COL},      &mat->tcpProb },
				{ "tcpAbruptChange",dtype, 3, { mxKnotNumTrend, ROW, COL},      &mat->tcpAbruptChange },
				{ "tcpCI",          dtype, 4, { mxKnotNumTrend, 2, ROW, COL },  &mat->tcpCI }
			};
			nfields = sizeof(fldList) / sizeof(FIELD_ITEM);
			memcpy(fieldList, fldList, nfields *sizeof(FIELD_ITEM));
		}
		else if (io->out.whichDimIsTime == 2) {
				FIELD_ITEM fldList[] =	{
				{ "time",			dtype, 2, {N, 1L,0L,0L },                 &mat->time },
				{ "sN",				dtype, 2, {ROW, COL,       },              &mat->sN },
				{ "tN",				dtype, 2, {ROW, COL,       },              &mat->tN },
				{ "sNProb",			dtype, 3, { ROW, mxKnotNumSeason + 1,COL,}, &mat->sNProb },
				{ "tNProb",			dtype, 3, { ROW, mxKnotNumTrend + 1, COL,}, &mat->tNProb },
				{ "sProb",			dtype, 3, {  ROW,N, COL,    },              &mat->sProb },
				{ "tProb",			dtype, 3, {  ROW,N, COL,    },              &mat->tProb },

				{ "s",				dtype, 4, {  ROW, N, q, COL,    },              &mat->s },
				{ "sCI",			dtype, 5, {  ROW, N, 2, q, COL  },              &mat->sCI },
				{ "sSD",			dtype, 4, {  ROW, N, q, COL     },              &mat->sSD },

				{ "t",				dtype, 4, { ROW, N, q, COL     },              &mat->t },
				{ "tCI",			dtype, 5, { ROW, N, 2, q, COL   },             &mat->tCI },
				{ "tSD",			dtype, 4, { ROW, N, q, COL   },                &mat->tSD },

				{ "b",				dtype, 4, { ROW, N, q, COL       },            &mat->b },
				{ "bCI",			dtype, 5, { ROW, N, 2, q,COL    },            &mat->bCI },
				{ "bSD",			dtype, 4, { ROW, N, q, COL  },                 &mat->bSD },

				{ "marg_lik",		dtype, 2, { ROW, COL       },               &mat->marg_lik },
				{ "sig2",			dtype, 4, { ROW, q, q, COL       },               &mat->sig2 },
		
				{ "bsign",			dtype, 4, { ROW,N, q, COL       },            &mat->bsign },
				{ "horder",			dtype, 3, { ROW,N,  COL      },             &mat->horder },
				{ "torder",			dtype, 3, { ROW,N,  COL      },             &mat->torder },

				{ "scp",            dtype, 3, { ROW,mxKnotNumSeason, COL},     &mat->scp },
				{ "scpProb",        dtype, 3, { ROW,mxKnotNumSeason, COL},     &mat->scpProb },
				{ "scpAbruptChange",dtype, 3, { ROW, mxKnotNumSeason,COL},     &mat->scpAbruptChange },
				{ "scpCI",          dtype, 4, { ROW, mxKnotNumSeason,  2,COL}, &mat->scpCI },

				{ "tcp",            dtype, 3, { ROW, mxKnotNumTrend, COL},      &mat->tcp },
				{ "tcpProb",        dtype, 3, { ROW, mxKnotNumTrend,COL},      &mat->tcpProb },
				{ "tcpAbruptChange",dtype, 3, { ROW, mxKnotNumTrend,COL},      &mat->tcpAbruptChange },
				{ "tcpCI",          dtype, 4, { ROW, mxKnotNumTrend, 2, COL },  &mat->tcpCI }
			};
			nfields = sizeof(fldList) / sizeof(FIELD_ITEM);
			memcpy(fieldList, fldList, nfields *sizeof(FIELD_ITEM));
		}
		else if (io->out.whichDimIsTime == 3) {
			FIELD_ITEM fldList[] =	{
				{ "time",			dtype, 2, {N, 1L,0L,0L },                 &mat->time },
				{ "sN",				dtype, 2, {ROW, COL,       },              &mat->sN },
				{ "tN",				dtype, 2, {ROW, COL,       },              &mat->tN },
				{ "sNProb",			dtype, 3, {  ROW, COL,mxKnotNumSeason + 1,}, &mat->sNProb },
				{ "tNProb",			dtype, 3, { ROW, COL ,mxKnotNumTrend  + 1,}, &mat->tNProb },
				{ "sProb",			dtype, 3, {  ROW, COL,N,    },              &mat->sProb },
				{ "tProb",			dtype, 3, {  ROW, COL,N,    },              &mat->tProb },

				{ "s",				dtype, 4, {  ROW, COL,  N, q,  },              &mat->s },
				{ "sCI",			dtype, 5, {  ROW, COL, N, 2,q },              &mat->sCI },
				{ "sSD",			dtype, 4, {  ROW, COL,N, q    },              &mat->sSD },

				{ "t",				dtype, 4, { ROW,  COL, N, q,    },           &mat->t },
				{ "tCI",			dtype, 5, { ROW,  COL, N, 2,q   },             &mat->tCI },
				{ "tSD",			dtype, 4, { ROW,  COL, N, q  },                &mat->tSD },

				{ "b",				dtype, 4, { ROW, COL, N, q,       },            &mat->b },
				{ "bCI",			dtype, 5, { ROW, COL, N, 2, q    },            &mat->bCI },
				{ "bSD",			dtype, 4, { ROW, COL, N, q },                 &mat->bSD },

				{ "marg_lik",		dtype, 2, { ROW, COL       },               &mat->marg_lik },
				{ "sig2",			dtype, 4, { ROW, COL, q,q       },               &mat->sig2 },
		
				{ "bsign",			dtype, 4, { ROW, COL, N, q,     },            &mat->bsign },
				{ "horder",			dtype, 3, { ROW, COL, N,    },             &mat->horder },
				{ "torder",			dtype, 3, { ROW, COL  ,N,  },             &mat->torder },

				{ "scp",            dtype, 3, {  ROW, COL,mxKnotNumSeason,},     &mat->scp },
				{ "scpProb",        dtype, 3, {  ROW, COL,mxKnotNumSeason,},     &mat->scpProb },
				{ "scpAbruptChange",dtype, 3, {  ROW, COL,mxKnotNumSeason,},     &mat->scpAbruptChange },
				{ "scpCI",          dtype, 4, {  ROW, COL, mxKnotNumSeason,  2,}, &mat->scpCI },

				{ "tcp",            dtype, 3, { ROW, COL,mxKnotNumTrend, },      &mat->tcp },
				{ "tcpProb",        dtype, 3, { ROW, COL,mxKnotNumTrend, },      &mat->tcpProb },
				{ "tcpAbruptChange",dtype, 3, { ROW, COL,mxKnotNumTrend, },      &mat->tcpAbruptChange },
				{ "tcpCI",          dtype, 4, { ROW, COL,mxKnotNumTrend, 2, },  &mat->tcpCI }
			};
			nfields = sizeof(fldList) / sizeof(FIELD_ITEM);
			memcpy(fieldList, fldList, nfields *sizeof(FIELD_ITEM));
		}
	
		END

	MV_EXTRA_PTR flag = &(opt->extra);
	if (!flag->computeCredible)
		RemoveField(fieldList, nfields, "sCI"), mat->sCI = NULL,
		RemoveField(fieldList, nfields, "tCI"), mat->tCI = NULL,
		RemoveField(fieldList, nfields, "bCI"), mat->bCI = NULL;

	if (!flag->computeSlopeSign) 	    RemoveField(fieldList, nfields, "bsign"), mat->bsign = NULL;
	if (!flag->computeHarmonicOrder) 	RemoveField(fieldList, nfields, "horder"), mat->horder = NULL;
	if (!flag->computeTrendOrder) 	    RemoveField(fieldList, nfields, "torder"), mat->torder = NULL;

	if (!flag->computeChangepoints) {
		RemoveField(fieldList, nfields, "scp"),    mat->scp = NULL,
		RemoveField(fieldList, nfields, "scpProb"), mat->scpProb = NULL,
		RemoveField(fieldList, nfields, "scpAbruptChange"), mat->scpAbruptChange = NULL,
		RemoveField(fieldList, nfields, "scpCI"), mat->scpCI = NULL,
		RemoveField(fieldList, nfields, "tcp"), mat->tcp = NULL,
		RemoveField(fieldList, nfields, "tcpProb"), mat->tcpProb = NULL,
		RemoveField(fieldList, nfields, "tcpAbruptChange"), mat->tcpAbruptChange = NULL,
		RemoveField(fieldList, nfields, "tcpCI"), mat->tcpCI = NULL;
	}

 	VOID_PTR  out;	
	out = PROTECT(CreateStructVar(fieldList, nfields));

	f32_seq(mat->time, io->meta.startTime, io->meta.deltaTime, N);
	//f32_seq(mat->time, 1.f, 1.f, N);
	//r_ippsMulC_32f_I(io->meta.deltaTime,                         mat->time, N);
	//r_ippsSubC_32f_I(-(io->meta.startTime - io->meta.deltaTime), mat->time, N);
	
	if (dtype == DATA_DOUBLE) {
		double *f64Ptr = (double*)mat->time;
		for (ptrdiff_t i = (ptrdiff_t)(opt->N - 1); i >= 0; i--)
			f64Ptr[i] = (double)mat->time[i];
	}

	UNPROTECT(1L);
	return out;
}


	
void MV_memfill_result_output(MV_RESULT_PTR  result, MV_OPTIONS_PTR  opt, const F32 nan)
{
	*result->sN		  = nan;
	*result->tN		  = nan;
	*result->marg_lik = nan;	

	I64  N  = opt->N;
	I64  q  = opt->q;
	I64  Nq = N*q;

	for (rI64 i = 0; i < N; i++)
		*(result->sProb + i) = nan,
		*(result->tProb + i) = nan;

	for (rI64 i = 0; i < Nq; i++)
		*(result->s     + i) = nan,
		*(result->t     + i) = nan,
		*(result->b     + i) = nan,
		*(result->sSD + i)   = nan,
		*(result->tSD + i)   = nan,
		*(result->bSD + i)   = nan;

	for (rI64 i = 0; i < q*q; i++)
		*(result->sig2 + i) = nan;

	MV_PRIOR_PTR prior = &(opt->prior);
	for (rI64 i = 0; i < (prior->seasonMaxKnotNum + 1); i++)
		*(result->sNProb + i) = nan;

	for (rI64 i = 0; i < (prior->trendMaxKnotNum + 1); i++)
		*(result->tNProb + i) = nan;

	MV_EXTRA_PTR flag = &(opt->extra);
	if (flag->computeCredible) 	{
		for (rI64 i = 0; i < (Nq*2); i++)
			*(result->sCI + i) = nan,
			*(result->tCI + i) = nan,
			*(result->bCI + i) = nan;
	}
	 
	if (flag->computeSlopeSign) {
		for (rI64 i = 0; i < Nq; i++) 	*(result->bsign + i) = nan;
	}
 
	if (flag->computeHarmonicOrder) 	{
		for (rI64 i = 0; i < N; i++) *(result->horder + i) = nan;
	}
 
	if (flag->computeTrendOrder) {
		for (rI64 i = 0; i < N; i++) 	*(result->torder + i) = nan;
	}

	if (flag->computeChangepoints)
	{
		for (rI32 i = 0; i < prior->seasonMaxKnotNum; i++)
			*(result->scp + i)             = nan,
			*(result->scpProb + i)         = nan,
			*(result->scpAbruptChange + i) = nan,
			*(result->scpCI + i)           = nan,
			*(result->scpCI + prior->seasonMaxKnotNum + i) = nan;


		for (rI32 i = 0; i < prior->trendMaxKnotNum; i++)
			*(result->tcp + i)             = nan,
			*(result->tcpProb + i)         = nan,
			*(result->tcpAbruptChange + i) = nan,
			*(result->tcpCI + i)           = nan,
			*(result->tcpCI + prior->trendMaxKnotNum + i) = nan;
	}
	 
}
 
void MV_allocate_result_output(
	MV_RESULT_PTR  result, MV_OPTIONS_PTR  opt,
	xMemPointers * _restrict MEM)
{
	memset(result, 0, sizeof(MV_RESULT));
	rI32 N = opt->N;
	I32  q  = opt->q;
	I32  Nq = q *opt->N;

	MEM_alloc(result->sN,       F32, (*MEM), 1,   0);
	MEM_alloc(result->tN,       F32, (*MEM), 1,   0);
	MEM_alloc(result->marg_lik, F32, (*MEM), 1,   0);
	MEM_alloc(result->sig2,     F32, (*MEM), q*q, 0);

	MEM_alloc(result->sNProb, I32, (*MEM),  (opt->prior.seasonMaxKnotNum + 1), 64 );
	MEM_alloc(result->tNProb, I32, (*MEM),  (opt->prior.trendMaxKnotNum  + 1), 64 );

	MEM_alloc(result->sProb,  I32, (*MEM),  N, 64 );
	MEM_alloc(result->tProb,  I32, (*MEM),  N, 64 );

	MEM_alloc(result->s,   F32, (*MEM), Nq, 64);
	MEM_alloc(result->t,   F32, (*MEM), Nq, 64);
	MEM_alloc(result->b,   F32, (*MEM), Nq, 64);
	MEM_alloc(result->sSD, F32, (*MEM), Nq, 64);
	MEM_alloc(result->tSD, F32, (*MEM), Nq, 64);
	MEM_alloc(result->bSD, F32, (*MEM), Nq, 64);
 

	if (opt->extra.computeCredible) {
		MEM_alloc(result->sCI, F32, (*MEM), 2 * Nq, 64);
		MEM_alloc(result->tCI, F32, (*MEM), 2 * Nq, 64);
		MEM_alloc(result->bCI, F32, (*MEM), 2 * Nq, 64);
	}
	if (opt->extra.computeSlopeSign) 	{
		MEM_alloc(result->bsign, I32, (*MEM), Nq, 64);
	} 
	if (opt->extra.computeSlopeSign) {
		MEM_alloc(result->horder, U32, (*MEM), opt->N, 64);
	}
	if (opt->extra.computeTrendOrder){
		MEM_alloc(result->torder, U32, (*MEM), opt->N, 64);
	} 
	
	if (opt->extra.computeChangepoints)
		MEM_alloc(result->scp,             F32, (*MEM), opt->prior.seasonMaxKnotNum, 64),
		MEM_alloc(result->scpProb,         F32, (*MEM), opt->prior.seasonMaxKnotNum, 64),
		MEM_alloc(result->scpAbruptChange, F32, (*MEM), opt->prior.seasonMaxKnotNum, 64),
		MEM_alloc(result->scpCI          , F32, (*MEM), opt->prior.seasonMaxKnotNum*2, 64),
 
		MEM_alloc(result->tcp,             F32, (*MEM), opt->prior.trendMaxKnotNum, 64),
		MEM_alloc(result->tcpProb,         F32, (*MEM), opt->prior.trendMaxKnotNum, 64),
		MEM_alloc(result->tcpAbruptChange, F32, (*MEM), opt->prior.trendMaxKnotNum, 64),
		MEM_alloc(result->tcpCI          , F32, (*MEM), opt->prior.trendMaxKnotNum*2, 64);	 
}

static I32 GetRawTimeDimension(MV_IO_PTR io){
		
	if ( io->ndim  == 2)           // Not a 3D stack input
		return io->dims[0];			

	if (io->ndim == 3L)    // A 3D stack input
		return io->dims[io->meta.whichDimIsTime - 1];
	
}
I32 MV_GetTotalNumberOfPixels(MV_OPTIONS_PTR opt){
	return opt->io->numOfPixels;
}
static void GetOffsetStride(MV_IO_PTR io, I64 idx, I64 *pStride, I64 *pOffset, I64 *pN)
{
	I64 N = GetRawTimeDimension(io);
	
	I64 stride, offset;	
	if (io->ndim ==  2)       // Not a 3D stack input
		stride = 1L,
		offset = (idx - 1)*N;	
	else if (io->ndim == 3L)  // A 3D stack input
	{
			int  ROW, COL;
			switch (io->meta.whichDimIsTime) {
			case 1:
				ROW    = io->dims[1];
				COL    = io->dims[2];
				stride = 1L;
				offset = (idx - 1)*N;
				break;
			case 2: {
				ROW = io->dims[0];
				COL = io->dims[2];
				int r, c;
				c = (idx-1) / ROW;
				r = idx - c*ROW;
				c = c + 1;
				stride = ROW;
				offset = (c - 1)*(N*ROW) + r - 1;
				break;
			}
			case 3: {
				ROW = io->dims[0],
				COL = io->dims[1];
				int r, c;
				c = (idx - 1) / ROW;
				r = idx - c*ROW;
				c = c + 1;
				stride = ROW*COL;
				offset = (c - 1)*ROW + r - 1;
				break;
			}

		} // switch ( io->timeDimensionIndex)	 

	} // (io->ndim == 3L)
	
	*pStride = stride;
	*pOffset = offset;
	*pN = N;
}
static void fetch_next_timeSeries_MEM(MV_YINFO_PTR  yInfo, int idx, F32PTR GlobalMEMBuf,  MV_OPTIONS_PTR  opt)
{   
	// Fetech the next time series, find missing rows, and normalize it.
	// idx: pixelIdx 
	I64  stride, offset, N;
	GetOffsetStride(opt->io, idx, &stride, &offset, &N);
	/****************************************/
	/* Now offest and stride are appropr-   */
	/*iately set up. Start to extract the ts*/
	/****************************************/
	I32      q = opt->q;
	F32PTR   Y = yInfo->Y;
	for (I32 i = 0; i < q; i++) {
		CopyStrideMEMToF32Arr(Y+i*N/*dst*/, opt->io->pdata[i] /*src*/, N, stride, offset,opt->io->dtype);
		f32_set_nan_by_value(Y+i*N, N, opt->io->meta.missingValue);
	}

	//Normalize Y with NaN ommitted and then
	//Pre-comoute Y'*Y: YtY_plus_Q using gemm
	yInfo->nMissing = f32_normalize_multicols_zeroout_nans(Y, yInfo->rowsMissing, N, N, q, yInfo->mean, yInfo->sd);
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, N, 1.0, Y, N, Y, N, 0.f, yInfo->YtY_plus_Q, q);
	yInfo->n = N - yInfo->nMissing;


}

 
static void fetch_next_timeSeries_MEM_irregular(MV_YINFO_PTR  yInfo, int idx, F32PTR GlobalMEMBuf, MV_OPTIONS_PTR  opt)
{   
		// Fetech the next time series, find missing rows, and normalize it.
	// idx: pixelIdx
	MV_IO_PTR io = opt->io;
	I64  stride, offset, Nold;
	GetOffsetStride(io, idx, &stride, &offset, &Nold);	
	/****************************************/
	/* Now offest and stride are appropr-   */
	/*iately set up. Start to extract the ts*/
	/****************************************/
	I32      q    = opt->q;
	I32      Nnew = opt->N;
	F32PTR   Y    = yInfo->Y;
	for (rI32 i = 0; i < q; i++) {
		CopyStrideMEMToF32Arr(   GlobalMEMBuf  /*dst*/,	opt->io->pdata[i]/*src*/, Nold, stride, offset, opt->io->dtype);
		f32_set_nan_by_value(GlobalMEMBuf, Nold, opt->io->meta.missingValue);
		
		tsAggegrationPerform(Y+Nnew*i, Nnew, GlobalMEMBuf, Nold, opt->io->T.numPtsPerInterval,
						     opt->io->T.sortedTimeIdx + opt->io->T.startIdxOfFirsInterval);
	}

	//Normalize Y with NaN ommitted and then
	//Pre-comoute Y'*Y: YtY_plus_Q using gemm
	yInfo->nMissing = f32_normalize_multicols_zeroout_nans(Y, yInfo->rowsMissing, Nnew, Nnew, q, yInfo->mean, yInfo->sd);
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, q, q, Nnew, 1.0, Y, Nnew, Y, Nnew, 0.f, yInfo->YtY_plus_Q, q);
	yInfo->n = Nnew - yInfo->nMissing;
}

void MV_fetch_next_timeSeries(
	MV_YINFO_PTR  yInfo, int pixelIndex,F32PTR GlobalMEMBuf, MV_OPTIONS_PTR  opt) {

		if (opt->io->isRegularOrdered)
			fetch_next_timeSeries_MEM(yInfo, pixelIndex, GlobalMEMBuf, opt);			
		else
			fetch_next_timeSeries_MEM_irregular(yInfo, pixelIndex, GlobalMEMBuf, opt);
}


#if !defined(NA_REAL)
	#define NA_REAL nan
#endif
static void GetOutputOffsetStride(MV_IO_PTR io, I64 idx,  I64 N, I64 *pStride, I64 *pOffset )
{
	I64 stride, offset;	
	if (io->ndim ==  2)       // Not a 3D stack input
	    stride = 1L, offset = (idx - 1)*N;	
	else if (io->ndim == 3L)  // A 3D stack input
	{
			int  ROW, COL;
			switch (io->meta.whichDimIsTime) {
			    case 1:
					ROW = io->dims[1];
					COL = io->dims[2];
					break;
				case 2:
					ROW = io->dims[0];
					COL = io->dims[2];
					break;
				case 3:
					ROW = io->dims[0],
					COL = io->dims[1];
			}

			switch (io->out.whichDimIsTime) {
			case 1:
				stride = 1L;
				offset = (idx - 1)*N;
				break;
			case 2: {			
				int c = (idx-1) / ROW;
				int r = idx - c*ROW;
				c = c + 1;
				stride = ROW;
				offset = (c - 1)*(N*ROW) + r - 1;
				break;
			}
			case 3: {			
				int c = (idx - 1) / ROW;
				int r = idx - c*ROW;
				c      = c + 1;
				stride = ROW*COL;
				offset = (c - 1)*ROW + r - 1;
				break;
			}

		} // switch ( io->timeDimensionIndex)	 

	} // (io->ndim == 3L)
	
	*pStride = stride;
	*pOffset = offset;	
}
void  MV_WriteOutput(MV_OPTIONS_PTR opt, MV_RESULT_PTR result, I64 pixelIndex) 
{
	MV_IO_PTR       io   = opt->io;
	MV_RESULT_PTR   mat  = io->out.result;

	const I64 N = opt->N;
	const I64 q = opt->q;
	const I64 seasonMaxKnotNum = opt->prior.seasonMaxKnotNum;
	const I64 trendMaxKnotNum  = opt->prior.trendMaxKnotNum;

	DATA_TYPE datType = io->out.dtype;

	I64 len, offset, stride;
	len = 1, offset = pixelIndex - 1, 	stride = 0;
	WriteF32ArraryToStrideMEM(result->sN,       mat->sN,      len,  stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->tN,       mat->tN,      len,  stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->marg_lik, mat->marg_lik,len, stride, offset, datType);

	len = q*q,  offset = pixelIndex - 1, stride = 1;
	WriteF32ArraryToStrideMEM(result->sig2,     mat->sig2,    len,  stride, offset, datType);

	len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
	WriteF32ArraryToStrideMEM(result->sProb, mat->sProb, len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->tProb, mat->tProb, len, stride, offset, datType);

	len = (seasonMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
	WriteF32ArraryToStrideMEM(result->sNProb, mat->sNProb,  len, stride, offset, datType);

	len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
	WriteF32ArraryToStrideMEM(result->tNProb, mat->tNProb, len, stride, offset, datType);

	len = N*q;                       GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
	WriteF32ArraryToStrideMEM(result->s,   mat->s  , len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->sSD, mat->sSD, len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->t,   mat->t  , len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->tSD, mat->tSD, len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->b,   mat->b  , len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->bSD, mat->bSD, len, stride, offset, datType);
	
	if (opt->extra.computeCredible) {
		len = N * q*2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->sCI, mat->sCI, len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->tCI, mat->tCI, len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->bCI, mat->bCI, len, stride, offset, datType);
	}

	if (opt->extra.computeSlopeSign)	  {
		len = N*q, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->bsign, mat->bsign,len, stride, offset, datType);
	}
		
	if (opt->extra.computeHarmonicOrder)	  {
		len    = N,	GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->horder, mat->horder, len, stride, offset, datType);
	}
		

	if (opt->extra.computeTrendOrder)	  {
		len = N, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->torder, mat->torder, len, stride, offset, datType);
	}

 
	if (opt->extra.computeChangepoints)	  {

		len = seasonMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->scp,             mat->scp,            len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->scpProb,         mat->scpProb,        len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->scpAbruptChange, mat->scpAbruptChange,len, stride, offset, datType);

		len = seasonMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->scpCI, mat->scpCI, len, stride, offset, datType);
		
		//--------------Trend----------------------------

		len = trendMaxKnotNum;      GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->tcp,             mat->tcp,             len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->tcpProb,         mat->tcpProb,         len, stride, offset, datType);
		WriteF32ArraryToStrideMEM(result->tcpAbruptChange, mat->tcpAbruptChange, len, stride, offset, datType);

		len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->tcpCI, mat->tcpCI, len, stride, offset, datType);
	} // sProb, tProb, sY, tY, sCI, tCI, tB, tBCI
 
}

void MV_print_options(MV_OPTIONS_PTR  opt)
{
	if (opt->extra.printOptions== 0) return;

	r_printf("......OPTIONS used in the MCMC inference ......\n");

 	MV_METADATA_PTR meta = &(opt->io->meta);
#if M_INTERFACE==1
	char filler = '.';
	r_printf("   metadata = [ ]\n");
#elif R_INTERFACE==1
	char filler = '$';
	r_printf("   metadata = list()\n");
#elif P_INTERFACE==1
	char filler = '.';
	r_printf("   metadata = list()\n");
#endif
	r_printf("   metadata%cmissingValue = %f\n", filler,  meta->missingValue);
	r_printf("   metadata%cwhichDimIsTime = %d\n", filler, meta->whichDimIsTime);
	r_printf("   metadata%cstartTime = %f\n", filler,      meta->startTime);
	r_printf("   metadata%cdeltaTime = %f\n", filler,      meta->deltaTime);
	r_printf("   metadata%cperiod = %d \n",   filler,  (int)meta->period);
	r_printf("   metadata%cmaxMissingRate = %f\n", filler, meta->maxMissingRate);
	r_printf("......End of displaying MetaData ......\n\n");

	MV_PRIOR_PTR prior = &(opt->prior);
#if M_INTERFACE==1
	r_printf("   prior = [ ]\n");
#elif R_INTERFACE==1
	r_printf("   prior = list()\n");
#endif	
	r_printf("   prior%cminSeasonOrder = %d\n", filler, prior->minSeasonOrder);
	r_printf("   prior%cwhichDimIsTime = %d\n", filler, prior->maxSeasonOrder);
	r_printf("   prior%cminTrendOrder = %d\n", filler, prior->minTrendOrder);
	r_printf("   prior%cmaxTrendOrder = %d\n", filler, prior->maxTrendOrder);
	r_printf("   prior%ctrendMinSepDist = %d\n", filler, prior->trendMinSepDist);
	r_printf("   prior%cseasonMinSepDist = %d\n", filler, prior->seasonMinSepDist);
	r_printf("   prior%ctrendMaxKnotNum = %d\n", filler, prior->trendMaxKnotNum);
	r_printf("   prior%cseasonMaxKnotNum = %d\n", filler, prior->seasonMaxKnotNum);
	r_printf("......End of displaying pripr ......\n\n");

 
	MV_MCMC_PTR mcmc = &(opt->mcmc);
#if M_INTERFACE==1
	r_printf("   mcmc = [ ]\n");
#elif R_INTERFACE==1
	r_printf("   mcmc = list()\n");
#endif	
	r_printf("   mcmc%cseed = %d\n", filler, mcmc->seed);
	r_printf("   mcmc%cmaxMoveStepSize = %d\n", filler, mcmc->maxMoveStepSize);
	r_printf("   mcmc%csamples = %d\n", filler, mcmc->samples);
	r_printf("   mcmc%cthinningFactor = %d\n", filler, mcmc->thinningFactor);
	r_printf("   mcmc%cburnin = %d\n",      filler, mcmc->burnin);
	r_printf("   mcmc%cchainNumber = %d\n", filler, mcmc->chainNumber);
	r_printf("   mcmc%ctrendResamplingOrderProb = %f\n", filler,  mcmc->trendResamplingOrderProb);
	r_printf("   mcmc%cseasonResamplingOrderProb = %f\n", filler, mcmc->seasonResamplingOrderProb);
	r_printf("......End of displaying mcmc ......\n\n");

 
	MV_EXTRA_PTR extra = &(opt->extra);
#if M_INTERFACE==1
	r_printf("   extra = [ ]\n");
#elif R_INTERFACE==1
	r_printf("   extra = list()\n");
#endif	
	r_printf("   extra%ccomputeCredible = %d\n", filler, extra->computeCredible);
	r_printf("   extra%cfastCIComputation = %d\n", filler, extra->fastCIComputation);
	r_printf("   extra%ccomputeChangepoints = %d\n", filler, extra->computeChangepoints);
	r_printf("   extra%ccomputeSlopeSign = %d\n", filler, extra->computeSlopeSign);
	r_printf("   extra%ccomputeHarmonicOrder = %d\n", filler, extra->computeHarmonicOrder);
	r_printf("   extra%ccomputeTrendOrder = %d\n", filler, extra->computeTrendOrder);
	r_printf("   extra%cconsoleWidth = %d\n", filler, extra->consoleWidth);
	r_printf("   extra%cprintProgressBar = %d\n", filler, extra->printProgressBar);
	r_printf("   extra%cprintOptions = %d\n", filler, extra->printOptions);
	r_printf("   extra%cnumCPUCoresToUse = %d\n", filler, extra->numCPUCoresToUse);
	r_printf("......End of displaying extra ......\n\n");

	//r_printf("   opt%coutputToDisk = %d\n", filler, opt->outputToDisk);
	//r_printf("   opt%clengthPerTimeSeries_infile = %d\n", filler, opt->N); 

}


ENABLE_MANY_WARNINGS