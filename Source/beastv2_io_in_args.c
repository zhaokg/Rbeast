#include "abc_000_warning.h"

#include "abc_001_config.h"

#include <math.h>
#include <string.h>
#include <time.h>

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_ide_util.h"  //printf
#include "abc_common.h"    //strcicmp
#include "abc_ts_func.h"
#include "abc_date.h"
#include "beastv2_func.h"    
#include "beastv2_io.h"

#include <stdio.h>	               //fprintf fopen FILE #include<stdio.h>  // Need _GNU_SOURCE for manylinux; otherwise report /usr/include/stdio.h:316:6: error: unknown type name '_IO_cookie_io_functions_t'

 
#define CondErrMsgRet0(cond, ...)   if(cond) { r_error(__VA_ARGS__); return 0;}
#define CondErrActionRet0(cond, Action, ...)   if(cond) { (Action) ;r_error(__VA_ARGS__); return 0;}
#define ifelse(cond, a, b)  ((cond)?(a):(b))
static int  GetArg_0th_Data(VOIDPTR prhs[], int nrhs, BEAST2_IO_PTR _OUT_ io) {

	// Check the number of inputs:  the first arg is the algorithm name.
	CondErrMsgRet0( nrhs < 2L, "ERROR: At least one input argument is needed!\n" );

	VOIDPTR DATA  = prhs[1];
	int     numel = GetNumberOfElements(DATA);
	// Check if the input 'data' is numeric and lengthy enough 	
	if (IsEmpty(DATA) || (/*beg*/ !(IsNumeric(DATA) && numel > 2) && !(IsStruct(DATA) && numel >= 1) && !IsCell(DATA) /*end*/) ) {
		/* IsCell  true Only for Matlab*/
		// Yobj has multivariate time series if it is a struct variable
		r_error("ERROR: The input data should be numeric and must be long enough.\n");
		return 0;
	}
 
	/*****************************************************************************************/
	/* If pY is not a string but a matrix or vector  find the dimesions, and configure io-*/
	/*****************************************************************************************/
	I32       q = 0;
	VOID_PTR  Y = NULL;
	if ((IsStruct(DATA) && numel >= 1) || IsCell(DATA)) {
		// For MRBEAST only: Yobj is a struct variable with multivariate time series						
		q = numel;
		// Mem to be dellocated in DeallocatTimeSeriesIO
		io->pdata = malloc(sizeof(VOID_PTR) * q);
		io->dtype = malloc(sizeof(DATA_TYPE) * q);
		//TODO: For MRBEAST,the dimensions of all elments are assumed to be the same
		// A sanity check needs to be done to ensure that is the case.
		for (int i = 0; i < q; i++) {
			Y = GetFieldByIdx(DATA, i);
			io->pdata[i] = GetData(Y);
			io->dtype[i] = GetDataType(Y);
		}		
	} else {
		//DATA is a vector of numeric type, and not a struct variable with multivariate time series			
		q = 1;
		//Mem to be dellocated in DeallocatTimeSeriesIO
		io->pdata   = malloc(sizeof(VOID_PTR) * q);
		io->dtype    = malloc(sizeof(DATA_TYPE) * q);
		io->pdata[0] = GetData(DATA);
		io->dtype[0] = GetDataType(DATA);
		Y = DATA;
	}
	for (int i = 0; i < q; i++) {
		CondErrMsgRet0( io->dtype[i] == DATA_UNKNOWN, "ERROR: The input data has an uknown numeric type!\n");		
	}
	io->q = q;

	/*************************************************/
	/* Get the dimension of the input              */
	/*************************************************/
	io->timedim = UnknownStatus;   // A filler to indicate hte status is undetermined
	I32 ndims = GetNumOfDim(Y);	
	if (ndims == 0 ||	 // the input is a vector: ndims == 0  is possible only for R vectors (not  Matlab)
		ndims == 1)      // the input is a vector: ndims == 1  is possible only for Python vectors (not R or matlab)
	{
		I32 N = GetNumberOfElements(Y);
		io->ndim    = 1L;
		io->dims[0] = N;
		io->dims[1] = 1L;
		io->dims[2] = 1L; 
		io->timedim = 1L; 		
	}
	else if (ndims == 2) {
		// Matlab: a vector or matrix;
		// R:      a matrix or a vector with a dim attribute.
		int N = GetDim1(Y);
		int M = GetDim2(Y);

		if (min(N, M) == 1L) {
			//PY is a vector
			N = max(N, M),
			io->ndim   = 1L;
			io->dims[0] = N;
			io->dims[1] = 1L;
			io->dims[2] = 1L;
			io->timedim = 1L;
		} else {
			//PY is a matrix
			io->ndim = 2L;
			io->dims[0] = N;
			io->dims[1] = M;
			io->dims[2] = 1L; // for 2d matrix input, the last dimenion is always set to 1L 
		}//PY is a matrix
	}
	else if (ndims == 3) {
		// If the input is a 3D array
		io->ndim = 3L;
		GetDimensions(Y, io->dims, 3L);	
	}
	else {
		r_printf("ERROR: The maximum dimension allowed is 3 when data is a 3D stack of images over time, but the input has a dimension of %d.\n", ndims);
		return 0;
	}

	return 1;
}

static float ParsePeriod(BEAST2_IO_PTR _OUT_ io);
static int   Parse_whichDimIsTime(BEAST2_IO_PTR _OUT_ io, int Nrawtime, int userWhichDim);

static int  GetArg_1st_MetaData(VOIDPTR prhs[], int nrhs, BEAST2_IO_PTR _OUT_ io) {
	

	BEAST2_METADATA_PTR meta = &io->meta;
	meta->nrhs = nrhs;

	/**********************************/
	// Get the type of metadata
	/**********************************/
	int METADATA_NONE          = _False_;
	int METADATA_NumericScalar = _False_;
	int METADATA_NumericVector = _False_;
	int METADATA_CharVector    = _False_;
	int METADATA_CharScaler    = _False_;
	int METADATA_Struct        = _False_;
	int METADATA_OTHER         = _False_;

	VOIDPTR pmeta              = (nrhs < 3) ? NULL : prhs[2L];
	VOIDPTR TIMEobj            = NULL;
	int     userWhichDimIsTime = UnknownStatus;
	if ( pmeta==NULL || IsEmpty(pmeta) ) {
		METADATA_NONE       = _True_;
	} else if ( IsNumeric(pmeta) ) {
		int  numel = GetNumberOfElements(pmeta);
		if      (numel == 0) METADATA_NONE          = _True_;
		else if (numel == 1) METADATA_NumericScalar = _True_;
		else if (numel > 1 ) METADATA_NumericVector = _True_, TIMEobj = pmeta;
		else 			     METADATA_OTHER         = _True_;
	} else if ( IsChar(pmeta) ) {
		int numel = GetNumberOfElements(pmeta);
		if      (numel == 0) METADATA_NONE          = _True_;
		else if (numel == 1) METADATA_CharScaler    = _True_;
		else if (numel > 1 ) METADATA_CharVector    = _True_, TIMEobj = pmeta;
		else 			     METADATA_OTHER         = _True_;
	} else if (IsStruct(pmeta) ) {
		METADATA_Struct      = _True_;
		TIMEobj              =  GetField123Check(pmeta, "time", 2) ; // if time is NULL or empty []	VOIDPTR  tmp;
		VOIDPTR  tmp;
		userWhichDimIsTime  = (tmp = GetField123Check(pmeta, "whichDimIsTime", 2)) ? GetScalar(tmp) : UnknownStatus;
	} else {
		METADATA_OTHER = _True_;
	}
	CondErrMsgRet0(METADATA_OTHER == _True_, "ERROR: The 'metadata' parameter given is of unsupported type.\n");

	//*****************************************************************
	// Determnie if it has a seaosn component or not: Moved here
	// becuase period is needed in TSaggragatePrepare()
	//*****************************************************************	
	meta->hasSeasonCmpnt = UnknownStatus;
	I08 ISDATE   = UnknownStatus;
	F32 START    = getNaN();
	F32 DT       = getNaN();
	F32 PERIOD   = getNaN();
	//////////////////////////////////////////////////////////////////////////////////
	if (METADATA_NONE || METADATA_NumericScalar || METADATA_CharVector || METADATA_CharScaler ||  METADATA_NumericVector) {
		//meta->isRegularOrdered  = 1; // We asumme the data is regular/orderred if no meta is supplied
	 
		START  = getNaN();
		DT     = getNaN();
		PERIOD = getNaN();
		if (METADATA_NumericScalar) {
			PERIOD = GetScalar(pmeta);
			if (!IsNaN(PERIOD) ) {
				if (PERIOD <= 0.) {
					q_warning("WARNING: A negative or zero value of period (%g) means no periodic/seasonal component in the input time series!\n", meta->period);
				} else{ //meta->period > 0;
					CondErrMsgRet0(!_IsAlmostInteger(PERIOD),
						"ERROR: When metadata is supplied as a single number %g, it must be an integer to specify the period of the regular time seires!\n", meta->period);
					PERIOD = round(PERIOD);
				}
			}
			meta->hasSeasonCmpnt = PERIOD <= 0.0? _False_ :_True_; // Period may be still a NAN  if meta is a NAN
		}	else if (METADATA_CharScaler) {
		// there is only one possibiliyt: Period nyst be the string 'none"
			char period[10+1];
			GetCharArray(pmeta, period, 10L);
			CondErrMsgRet0( strcicmp(period, "none") != 0, "ERROR: When metadata is supplied as a string to speciify period, it can only be 'none'!\n");
 			meta->hasSeasonCmpnt = _False_;
			PERIOD               = 0.0;
		}	else {
			// has a seasonal component with period still being NAN and to be determined
			meta->hasSeasonCmpnt = _True_;
		}
		meta->seasonForm        = ifelse( meta->hasSeasonCmpnt, 'S','N');

		meta->hasOutlierCmpnt   = _False_;		
		meta->detrend           = _False_;
		meta->deseasonalize     = _False_;
		meta->missingValue      = getNaN();// FLOAT_MAX;
		meta->maxMissingRate    = 0.75f;
	}
	else if (METADATA_Struct) {
	//*****************************************************************
	// Determnie if it has a seaosn component or not when  METADATA_Struct =1
	//*****************************************************************	
		VOIDPTR  tmp;

		//  Deermine the startTime property
		tmp = GetField123Check(pmeta, "startTime", 2); // IsEmpty is checked inside GetField123Check
		TimeScalarInfo timevalue;
		Parse_SingelDateObject(tmp, &timevalue);
		START = timevalue.value;
		if (timevalue.unit == 'Y')  ISDATE = _True_;
		CondErrMsgRet0(timevalue.unit == 'B', "ERROR: cannot interpret the metadata$startTime input,'\n");

		//  Deermine the deltaTime property		
		tmp = GetField123Check(pmeta, "deltaTime", 3);
		Parse_TimeIntervalObject(tmp, &timevalue);
		DT  = timevalue.fyear;
		if (timevalue.unit != 'U' && timevalue.unit != 'B') ISDATE = _True_;
		CondErrMsgRet0(timevalue.unit == 'B', "ERROR: cannot interpret the metadata$deltaTime input,'\n");
		CondErrMsgRet0(DT             <= 0,   "ERROR: metadata$deltaTime must be a positive time interval!\n");

		//  Deermine the period property
		tmp = GetField123Check(pmeta, "period", 2);
		Parse_TimeIntervalObject(tmp, &timevalue);
		PERIOD  = timevalue.fyear;
		if (timevalue.unit != 'U' && timevalue.unit != 'B') ISDATE = _True_;
		CondErrMsgRet0(timevalue.unit == 'B', "ERROR: cannot interpret the metadata$period input,'\n");

		if (PERIOD<=0) {
			meta->hasSeasonCmpnt = _False_;
			q_warning("WARNING: A negative or zero value of period (%g) indicates no periodic/seasonal component in the input time series!\n", PERIOD);
		}
		if (PERIOD > 0 && PERIOD /DT >= 2.0) {
			meta->hasSeasonCmpnt = _True_;			
		}
		if (PERIOD > 0 && PERIOD /DT < 2.0) {
			meta->hasSeasonCmpnt = _False_;
			q_warning("WARNING: The value metadata$period/metadata$deltTime=%g/%g=%g is unreasonalbe for a siginal with seasonal componet, so no seasonal component is "
				"assumed and the trend-only model is fitted. Perios is also set to 0.0\n", PERIOD, DT, PERIOD/DT);
			PERIOD = 0;
		}

		/////////////////////////////////////////////////////////// 
		//   Get the season string to determine the number and types of season components: harnomic or dummy
		////////////////////////////////////////////////////////// 
		tmp = ( tmp=GetField123Check(pmeta, "season", 2) , tmp && !IsChar(tmp) ? NULL : tmp);
		char season[20 + 1];
		GetCharArray(tmp, season, 20); // Will check if tmp is NULL inside the function
		ToUpper(season);
		char a = season[0], b = season[1];
	
	   if (meta->hasSeasonCmpnt == _False_ && a != 'N' && a != '\0' ) {
		   // period has been set and used to determine hasSeasonCMpnt or not
			   q_warning("WARNING: A confilict found between metadata$season=%s and period=%g. '%s' suggests a time series with periodic variations but "
				         "period=%g or period='none' suggests a trend-only time series without any periodic variations. The season parameter is "
				         "ignored and no seasonal component is assumed for the input.\n", season, PERIOD, season, PERIOD);
		} else if  (tmp  && IsChar(tmp) ) {			
			// tmp may be an empty string "" (a=0).
			int  hasSeasonCmpt_bySeasonStr = _True_;			
			if    ((a=='N' && b=='O') || a == '\0')		hasSeasonCmpt_bySeasonStr = _False_;								//none
			else if (a=='H' && b=='A') 					hasSeasonCmpt_bySeasonStr = _True_, meta->seasonForm = 'S';    	//harmonic
			else if (a=='D' && b=='U')					hasSeasonCmpt_bySeasonStr = _True_, meta->seasonForm = 'D';		//dummy
			else if (a=='S' && b=='V')					hasSeasonCmpt_bySeasonStr = _True_, meta->seasonForm = 'V';		//svd			
			else {
				hasSeasonCmpt_bySeasonStr = _True_;
				meta->seasonForm          = 'S';  //the default form is harmonic
				q_warning("WARNING: metadata$season='%s' has an unrecongizable string. The default season='harmonic' is used instead.\n", season);
			}

			if (meta->hasSeasonCmpnt == _True_ && hasSeasonCmpt_bySeasonStr == _False_) {
				hasSeasonCmpt_bySeasonStr = _True_;
				meta->seasonForm          = 'S';
				q_warning("WARNING: A confilict found between metadata$season='none' and period=%g. season='none' suggests a time series with "
					       "no periodic variations but period=%g suggests otherwise. The season='none' parameter is ignored and "
					      "the data is assumed to have a seasonal component.\n", PERIOD, PERIOD);
			}
			meta->hasSeasonCmpnt = hasSeasonCmpt_bySeasonStr;
		} 	else	{
			if (meta->hasSeasonCmpnt == UnknownStatus || meta->hasSeasonCmpnt == _True_) {
				meta->hasSeasonCmpnt = _True_;
				meta->seasonForm     = 'S';
				q_warning("WARNING: metadata$season is either missing or not given as a valid specifier string (e.g., none, harmonic, or dummy). A default season='harmonic' is assumed.\n");
			}
		}	

		if (meta->seasonForm == 'V') {
			meta->svdTerms_Object   = GetFieldCheck(pmeta, "svdTerms");
			meta->svdYseason_Object = GetFieldCheck(pmeta, "svdYseason");
		}

		meta->hasOutlierCmpnt = (tmp = GetField123Check(pmeta, "hasOutlierCmpnt", 2)) ? GetScalar(tmp) : _False_;
		meta->detrend         = (tmp = GetField123Check(pmeta, "detrend"      ,3)) ?  GetScalar(tmp)   : _False_;
		meta->deseasonalize   = (tmp = GetField123Check(pmeta, "deseasonalize",3)) ?  GetScalar(tmp)   : _False_;
		meta->missingValue	  = (tmp = GetField123Check(pmeta, "missingValue",0))   ? GetScalar(tmp)   : getNaN();// FLOAT_MAX;
		meta->maxMissingRate  = (tmp = GetField123Check(pmeta, "maxMissingRate",0)) ? GetScalar(tmp)   : 0.75;
	} //else if (nrhs >= 3)
	else {
		CondErrMsgRet0(  _True_, "ERROR: The 'metadata' parameter given is of unsupported type.\n");
	}

	////////////////////////////////////////////////////////////////////////////////
	// Sanity check
	CondErrMsgRet0(meta->hasSeasonCmpnt == UnknownStatus, "ERROR: Cnnot determine whether the input time series has a seasonal/periodic componnet or net.\n");
	if (meta->hasSeasonCmpnt == 0) PERIOD = 0;        //just double check to make sure it is the case that period=0 for trend-only data
	if (DT  < 0.)                  DT     = getNaN(); //dT should never < 0 at this point, but just double check.
	////////////////////////////////////////////////////////////////////////////////


	TimeVecInfo tvec = { 0 };
	TimeVec_init(&tvec);

	/**********************************/
	// If there is a time object suppiked
	/**********************************/
     if (TIMEobj) { 
		//  Allocated mem that needs to be freeed explicilty for tvec.fyear, but noo mem allocated if TIMEObj
		//  is NULL, which means that the ts is regular/ordered, as dertermined by start and dt only
		int Nrawtime = TimeVec_from_TimeObject(TIMEobj, &tvec);        // isDate may be updated inside based on the TimeObj
	    if(tvec.isDate == 1) ISDATE =_True_;                            // this is the LAST chance to decide whether the time is date
	   
		// Nawtime must be larger than 1L.
		CondErrMsgRet0( Nrawtime<=1, "ERROR: Unable to read and intepret 'time' or 'metadata$time'!\n");

		// io->timedim and meta->whichDim are filled insside
		// Nrawtime>1: we know the time series length
		userWhichDimIsTime=Parse_whichDimIsTime(io, Nrawtime, userWhichDimIsTime);  
		CondErrMsgRet0(userWhichDimIsTime <= 0, "ERROR: Unable to dtermine which dimo of the input data refers to the time'!\n");
	
	} else {
	 /**********************************/
    // TIMEobj=NULLL If there is no time object suppiled; it must be a regular ts
	// Assume it is a regular time series
	/**********************************/
		int Nrawtime = 0;

		// io->timedim and meta->whichDim are filled insside
		// Nrawtime>0: we know the time series length
		userWhichDimIsTime = Parse_whichDimIsTime(io, Nrawtime, userWhichDimIsTime);
		CondErrMsgRet0(userWhichDimIsTime <= 0, "ERROR: Unable to dtermine which dimo of the input data refers to the time'!\n");

		int isDefaultDeltaTime = 0;
		int isDefaultStartTime = 0;
		if (IsNaN(START)) {
			START              = 1.f;
			isDefaultStartTime = 1L;
		}
		if (IsNaN(DT)) {
			DT   = 1.f;
			isDefaultDeltaTime = 1L;
		}

		char *warningMsg= "WARNING: If the input data is regular and ordered in time, the times of individual datapoints are determined fully by 'metadata$startTime' and 'metadata$deltaTime'.";
		if (isDefaultDeltaTime && isDefaultStartTime) {
			q_warning("%s But startTime and deltaTime are missing and a default value 1 is used for both!\n", warningMsg);
			if (ISDATE == _True_) { // as determined from the period field
				q_warning("WARNING! The default sartTime=1 and deltaTime=1  are used, but the 'period' field specifies that the time is date\n");
			}
		} else if (isDefaultStartTime) {
			q_warning("%s But startTime is missing and a default value 1 is used!\n", warningMsg);
			if (ISDATE == _True_) { // as determined from the period field
				q_warning("WARNING! The default startTime=1 is used, but the 'period' field specifies that the time is date\n");
			}
		} else if (isDefaultDeltaTime) {
			q_warning("%s But deltaTime is missing and a default value 1 is used!\n", warningMsg);
			if (ISDATE == _True_) { // as determined from the period field
				q_warning("WARNING! The default deltaTime=1 is used, but the 'period' field specifies that the time is date\n");
			}
		}

		int N = io->dims[userWhichDimIsTime - 1];   // filled in Parse_whichDimIsTime;

		TimeVec_from_StartDeltaTime(&tvec, START, DT, N, ISDATE);         // isDate is not updated inise		
		PERIOD          = tvec.isDateNum == 1 ? PERIOD * 365 : PERIOD; // time unit may be changed to datenum inside the function above
		START           = tvec.data_start;   // time unit might have been changed
		DT              = tvec.data_dt;      // time unit might have been changed
	}
	// No change to update IsData any longer, so if it is still undetermined, then, it is not date
	ISDATE = ISDATE == UnknownStatus ? _False_ : ISDATE; 

	// Sorted_time_indices are allocated here taht need to be de-allocated explicilty
	TimeVec_SortCheckRegularOrder(&tvec); // Update data_start and data_dT

	tvec.isDate       = ISDATE;
	tvec.data_period  = PERIOD;
	tvec.out.start    = START; // may be still NA in the case of TimeObj != NULL
	tvec.out.dT       = DT;  // may be still NA in the case of TimeObj != NULL

	// datt_period may be changed inside when f32time is converted to days
   // out.start may be adjustd slightly to better match the real start of the time series
	TimeVec_UpdateUserStartDt(&tvec); // PERIOD may be needed in this function to get a better estiamte of dT
	                                  
	if (TIMEobj && IsNaN(DT)) {
	// When TIMEObj is not NULL and dT is not available from the input
	//FOr regular inpus: q_warning("WARNING: The input time series is regular and metadata$deltaTime/deltat is missing, the regular time interval %g is used. If not making sense, "
	//      "please specify metadata$deltaTime/deltat explicitly. \n", io->T.dT);
		F32 dT_new = tvec.out.dT;
		if (io->T.isRegular == 0) {
			if (ISDATE && io->T.out.asDailyTS==1) {
				q_warning("WARNING: The input time series is irregular (or may span across leap years) and BEAST needs to aggregate/resample it into regular data "
					"at a user-specified interval 'metadata$deltaTime/deltat' (for faster computation). But deltaTime is missing, a best guess of it %g year = %g months = %g days is used. "
					"If not making sense, please specify metadata$deltaTime/deltat explicitly. \n", dT_new, dT_new * 12, dT_new * 365.0);
			} else {
				q_warning("WARNING: The input data is irregular and for faster computaiton, BEAST needs to aggregate/resample it into regular data "
					"at a user-specified interval 'metadata$deltaTime/deltat' But deltaTime is missing, a best guess of it %g is used. "
					"If not making sense, please specify metadata$deltaTime/deltat explicitly. \n", dT_new);
			}
		}
	}
	PERIOD = tvec.data_period;
	START  = tvec.out.start;
	DT     = tvec.out.dT;;

	io->N = tsAggegrationPrepare(&tvec);

	/////////////////////////////////////////////////////////
	// PERIOD is still uknown. Guess it via auto-correlation
	// Up to this point, if period=NAN, it will be further determined via auto-correlaton, It must be an intger
	/////////////////////////////////////////////////////////	
	meta->IsPeriodEstimated = 0;
	if ( meta->hasSeasonCmpnt && IsNaN(PERIOD)) {		
		io->T = tvec; // TODO: must update io->T because it is needed inside ParsePeriod
		F32  estPeriod = ParsePeriod(io);
		char* season;
		char* msg;
		season =io->meta.seasonForm == 'S' ?      "harmonic" :
			    io->meta.seasonForm == 'V' ?     "svd":
			                                     "dummy";
		msg   = "suggests that the time series has a periodic/seasonal component. \"metadata$period\" is needed but missing.";
		if (estPeriod > 0) {
			q_warning( "WARNING: metadata$season='%s' %s A BEST GUESS of numbers of datapoints per period is %d, giving period = num_sample_per_period * deltaTime "
				"= %d*%g = %g. Please make sure this estimate makes sense; otherwise, the BEAST decomposition result will be incorrect.\n",
				season,  msg,  (int)estPeriod, (int)estPeriod, DT, DT * estPeriod);
		}	else {
			r_error("ERROR: metadata$season='%s' %s BEAST tried to estimate it via an auotcorrelation method but failed to get a reliable estimate. "
				"Please specify the period value EXPLICILTY. Or if your input has no periodic/seasonal component at all, "
				" set metadata$season='none' or period=0, which will fit a trend-only model.\n", season, msg);
			return 0;
		}
		meta->IsPeriodEstimated = 1;
		PERIOD = estPeriod * DT; // Convert period to the unit of deltaTime
	}

 
	// Another check on PERIOD: make sure there are at least 2 point per period
	if (meta->hasSeasonCmpnt == 1) {	
		F32 freq = PERIOD / DT;
		CondErrMsgRet0(freq <= 0, "ERROR: BEAST can't handle a time series with a periodic componnet (\"metadata$season='%s'\") but with metadata$period=%g or period='none' specified. If you mean to "
				    "handle trend-only time series, use the trend-only BEAST version by specifying metadata$season='none'.", meta->seasonForm == 'S' ? "harmonic" : "dummy" ,PERIOD);
		CondErrMsgRet0(freq <2, "ERROR: BEAST can't handle a time series with a periodic componnet (\"metadata$season='%s'\") but with metadata$period=%g and metadata$deltaTime=%g specified, "
			                    "which gives only 'period/deltatime=%g' data points per period, If you mean to "
							    "handle trend-only time series, use the trend-only BEAST version by specifying metadata$season='none'.", meta->seasonForm == 'S' ? "harmonic" : "dummy", PERIOD, DT, PERIOD/DT);
		CondErrMsgRet0(meta->seasonForm == 'D' && !IsNaN(freq) && !_IsAlmostInteger(freq), "ERROR: For a dummy seasonal component (\"metadata$season='dummy'\"), metadata$period=%g must be a multiple of metadata$deltaTime by an INTEGER number. "
			"Your period/deltaTime ratio is %g.", freq* DT, freq);
		CondErrMsgRet0(meta->seasonForm == 'V' && !IsNaN(freq) && !_IsAlmostInteger(freq), "ERROR: For a SVD seasonal component (\"metadata$season='dummy'\"), metadata$period=%g must be a multiple of metadata$deltaTime by an INTEGER number. "
			"Your period/deltaTime ratio is %g.", freq* DT, freq);	 
	}


	meta->startTime   = START;  //shound't be NAN
	meta->deltaTime   = DT;     //shound't be NAN
	meta->period      = PERIOD/DT;

	TimeVec_kill_fyearArray(&tvec); // Deallocate f64time only; the other allocate mem is still needed
	io->T = tvec;

 
	return 1;
}

static int Parse_whichDimIsTime(BEAST2_IO_PTR _OUT_ io, int Nrawtime, int userWhichDim) {
	   // Nrawtime = 0 when TIMEObject is NULL (only start and dt are provided)
	   // Nrawtime >1 when TIMEObjc is not NULL (metadata$time is provided)

		int whichDim_final = userWhichDim;

		if ( io->ndim == 1 && userWhichDim != UnknownStatus && userWhichDim != 1) {
			q_warning("WARNING: metadata$whichDimIsTime = %d is ignored because 'whichDimIsTime' is used only for 2D matrix or 3D array inputs but your input is a 1D vector.\n", userWhichDim);
            // in this case, io->timedim has been assigned 1 in Get_DATA
		}
    

		if (Nrawtime > 0) {
			// TImeobjec is not NULL: We know the the length of the time dim

			int matcheNumDims = (Nrawtime == io->dims[0]) + (Nrawtime == io->dims[1]) + (Nrawtime == io->dims[2]);
			if (matcheNumDims == 0) {
				r_error("ERROR: The input data must have the same length as the time in metadata.\n");
				return -1;
			}
			else if (matcheNumDims == 1) {
				int  timeDimMatched;
				if (Nrawtime == io->dims[0])  timeDimMatched = 1;
				if (Nrawtime == io->dims[1])  timeDimMatched = 2;
				if (Nrawtime == io->dims[2])  timeDimMatched = 3;
				if (userWhichDim != UnknownStatus && userWhichDim != timeDimMatched) {
					q_warning("WARNING: the specified metadata$whichDimIsTime=%d is ignored; 'whichDimIsTime=%d' is instead used based on the match between the input data and time.\n", userWhichDim, timeDimMatched);
				}
				whichDim_final = timeDimMatched;
			}
			else { // matcheNumDims =2 or 3
				if (userWhichDim == UnknownStatus || (io->ndim == 2 && userWhichDim != 1 && userWhichDim != 2)) {
					r_error("ERROR: For a 2D matrix input of size [%d x %d] (i.e., multiple time series), metadata$whichDimIsTime must be given "
						"to tell which dim of the matrix  refers to time. It must take a value out of 1 or 2 only.\n", io->dims[0], io->dims[1]);
					return 0;
				}
				if (userWhichDim == UnknownStatus || (io->ndim == 3 && userWhichDim != 1 && userWhichDim != 2 && userWhichDim != 3)) {
					r_error("ERROR: For a 3D array input of size [%d x %d x %d] (i.e., stacked time series images), metadata$whichDimIsTime must be given "
						"to tell which dim of the 3D array  refers to time. It must take a value out of 1, 2 or 3 only.\n", io->dims[0], io->dims[1], io->dims[2]);
					return 0;
				}

				if (userWhichDim>3 || userWhichDim <1) {
					r_error("ERROR: the input (whichDimIsTime=%d) muust be an integer of 1, 2, or 3.\n", userWhichDim+1);
					return 0;
				} else {
					if (io->dims[userWhichDim - 1] != Nrawtime) {
						r_error("ERROR: The length of the time dimension of the input (whichDimIsTime=%d) doesn't match the length of time/metadata$time (i.e., %d!=%d).\n", userWhichDim, io->dims[userWhichDim], Nrawtime);
						return 0;
					}
				}
		
				whichDim_final = userWhichDim;
			}
		}	
		else {
		// Nraw=0: rregular inputs wihtout the time object

			if (io->timedim == UnknownStatus) {
				// io->timedim has been determined for 1D signals; this branch is only for 2D or 3D data 
				if (userWhichDim == UnknownStatus || (io->ndim == 2 && userWhichDim != 1 && userWhichDim != 2)) {
					r_error("ERROR: For a 2D matrix input of size [%d x %d] (e.g., multiple time series), metadata$whichDimIsTime must be given "
						"to tell which matrix dim refers to time. It must take a value out of 1 or 2 only.\n", io->dims[0], io->dims[1]);
					return 0;
				}
				if (userWhichDim == UnknownStatus || (io->ndim == 3 && userWhichDim != 1 && userWhichDim != 2 && userWhichDim != 3)) {
					r_error("ERROR: For a 3D array input of size [%d x %d x %d] (i.e., stacked time series images), metadata$whichDimIsTime must be given "
						"to tell which aray dim refers to time. It must take a value out of 1, 2 or 3 only.\n", io->dims[0], io->dims[1], io->dims[2]);
					return 0;
				}
				whichDim_final = userWhichDim;
			}	else {
				// it can be only 1
				whichDim_final = io->timedim;
			}
		
		}

		io->timedim = io->meta.whichDimIsTime = whichDim_final;

		if (io->meta.whichDimIsTime == 1) { io->rowdim = 2, io->coldim = 3, io->timedim = 1; }
		if (io->meta.whichDimIsTime == 2) { io->rowdim = 1, io->coldim = 3, io->timedim = 2; }
		if (io->meta.whichDimIsTime == 3) { io->rowdim = 1, io->coldim = 2, io->timedim = 3; }

		io->imgdims[0]  = io->dims[io->rowdim - 1];
		io->imgdims[1]  = io->dims[io->coldim - 1];
		io->numOfPixels = (I64)io->dims[0] * io->dims[1] * io->dims[2] / io->dims[io->timedim - 1L];

		//int N = io->dims[io->timedim - 1L];
		return whichDim_final;		
}

 static float ParsePeriod(BEAST2_IO_PTR _OUT_ io ) {
		
		/*************************************************/
		// There is a seasonal componet, and we need to check and get period	 
		/*************************************************/
		//if ( io->meta.hasSeasonCmpnt==0 ||    // Period is not needed for trend-only data
		//	 io->meta.period        > 0    )  // or if there is a season compnt and the period is already a valid value			 

		/*************************************************/
		//  if period<=0 or IsNan(period)
		//  Determine the period perameter via auto-correlation
		/*************************************************/
		BEAST2_YINFO Yinfo;
		F32PTR       MEMBUF;
		int    N    = io->N;		
		int    q    = io->q;
		int    Nraw = io->dims[io->timedim - 1];
		F32PTR tmp  = malloc(sizeof(F32)*(N*q +N+q +q +q*q +Nraw )); //alocate mem for yInfo and MEMBUF

		Yinfo.Y                 = tmp;
		Yinfo.rowsMissing       = tmp+N*q;
		Yinfo.mean              = tmp + N * q+N;
		Yinfo.sd                = tmp + N * q +N+q;
		Yinfo.YtY_plus_alpha2Q  = tmp + N * q +N+q+q;
		MEMBUF                  = tmp + N * q +N+q+q+q*q; //needed only for irregular time series

		F32 nan                 = getNaN();
		F32 period              = nan;
		I32 goodPixelVisited    = 0;
		I32 MaxNumPixelVisisted = 200;
		// THe pixel index is 1-based and not zero-based
		for (int i = 1; i <= io->numOfPixels; ++i) {

			BEAST2_fetch_timeSeries(&Yinfo, i, MEMBUF, io);
			// yInfo has been fillted above and now compaute mean, std, and YtY						
			// Normalize Y with NaN ommitted and then pre-comoute Y'*Y: YtY_plus_Q using gemm
			Yinfo.nMissing = f32_normalize_multicols_zeroout_nans(Yinfo.Y, Yinfo.rowsMissing, N, N, q, Yinfo.mean, Yinfo.sd);

			U08 skipCurrentPixel = Yinfo.nMissing > (N * io->meta.maxMissingRate);
			if (skipCurrentPixel) {
				continue;
			}		
			
			for (int j = 0; j < q; ++j) {

				F32PTR y = Yinfo.Y + j * N;
				for (int k = 0; k < Yinfo.nMissing; ++k) {
					y[Yinfo.rowsMissing[k]] = nan;
				}

				F32 jthPeriod = DeterminePeriod(y, N);  // return -1 if failing to estimate the period
				if (j == 0) { 
					//if it is the 1st out of the q time series	
					period = jthPeriod;
					if (period < 0)	break;					
				} else {
					if (jthPeriod != period ) {
						//if a later time series dones't give the same period as the first one
						period = nan;
						break;
					}
				}
			}

			if (period > 0   )                           	break;
			if (++goodPixelVisited > MaxNumPixelVisisted) 	break;
		}
		free(tmp);
 
	return period;
}

static int __GetPrecPriorType( VOID_PTR S ) {

	VOID_PTR  tmp = GetField123Check(S, "precPriorType",5);
	
	if (tmp == NULL)
		return UniformPrec;

	// if tmp is not NULL
	if (IsNumeric(tmp)) {
		I32 value = GetScalar(tmp);
		if (value == 0) return ConstPrec;
		if (value == 1) return UniformPrec;
		if (value == 2) return ComponentWise;
		if (value == 3) return OrderWise;

		q_warning("WARNING: The arg prior$precPriorType=(%d) is not a valid value; the default prior$precPriorType='%s' is assumed instread!", value, "uniform");
		return UniformPrec;
	} 
	else if (IsChar(tmp)) {
		char str[60+1];
		GetCharArray(tmp,str, 60);
		ToUpper(str);
		char a = str[0];
		char c = str[2];
		if (a =='U')				return UniformPrec;
		if (a == 'O')				return OrderWise;
		if (a == 'C' && c == 'M')	return ComponentWise;
		if (a == 'C' && c == 'N')	return ConstPrec;	

		q_warning("WARNING: The arg prior$precPriorType=(%s) is not recongizable; the default prior$precPriorType='%s' is assumed instread!", str, "uniform");
		return UniformPrec;
	}
	
	q_warning("WARNING: The arg prior$precPriorType has an supported format or value; the default prior$precPriorType='%s' is assumed instread!",  "uniform");
	return UniformPrec;
}

static int  GetArg_2nd_Prior__(VOIDPTR prhs[], int nrhs, BEAST2_PRIOR_PTR prior, BEAST2_IO_PTR io)
{	 
	// Before running this fuction,  period and N must have been determined

	struct PRIOR_MISSING {	
		U08   seasonMinOrder, seasonMaxOrder, trendMinOrder, trendMaxOrder;
		U08   trendMinSepDist, seasonMinSepDist;
		U08   trendMinKnotNum, seasonMinKnotNum;
		U08   trendMaxKnotNum, seasonMaxKnotNum,  outlierMaxKnotNum;

		U08   trendLeftMargin, trendRightMargin;
		U08   seasonLeftMargin, seasonRightMargin;

		U08   seasonBasisFuncType, trendBasisFuncType,  outlierBasisFuncType;
		U08   modelPriorType;
		U08   precPriorType;

		U08   K_MAX;
		U08   sigFactor;
		U08   outlierSigFactor;
		U08   sig2;
		U08   precValue;
		U08   alpha1, alpha2, delta1, delta2;

 
	} m = {0,};

	#define o  (*prior)

	

	if (nrhs < 4) 	
		memset(&m, 1L, sizeof(struct PRIOR_MISSING));

	if (nrhs >= 4) {		 
		VOIDPTR S = prhs[3L];
		if (!IsStruct(S)) {
			q_warning("WARNING: The arg 'prior' is ignored because it is not a List/Struct variable.");
			memset(&m, 1L, sizeof(struct PRIOR_MISSING));
		}
		else {
			VOIDPTR tmp;
			if (io->meta.hasSeasonCmpnt) {
				o.seasonMinOrder = (tmp = GetField123Check(S, "seasonMinOrder", 10)) ? GetScalar(tmp) : (m.seasonMinOrder = 1);
				o.seasonMaxOrder = (tmp = GetField123Check(S, "seasonMaxOrder", 10)) ? GetScalar(tmp) : (m.seasonMaxOrder = 1);
				o.seasonMinSepDist = (tmp = GetField123Check(S, "seasonMinSepDist", 10)) ? GetScalar(tmp) : (m.seasonMinSepDist = 1);
				o.seasonMinKnotNum = (tmp = GetField123Check(S, "seasonMinKnotNum", 10)) ? GetScalar(tmp) : (m.seasonMinKnotNum = 1);
				o.seasonMaxKnotNum = (tmp = GetField123Check(S, "seasonMaxKnotNum", 10)) ? GetScalar(tmp) : (m.seasonMaxKnotNum = 1);

				o.seasonLeftMargin  = (tmp = GetField123Check(S, "seasonLeftMargin", 10)) ? GetScalar(tmp) : (m.seasonLeftMargin = 1);
				o.seasonRightMargin = (tmp = GetField123Check(S, "seasonRightMargin", 10)) ? GetScalar(tmp) : (m.seasonRightMargin = 1);
			}

			o.trendMinOrder     = (tmp = GetField123Check(S, "trendMinOrder",  10)) ?   GetScalar(tmp) : (m.trendMinOrder = 1);
			o.trendMaxOrder     = (tmp = GetField123Check(S, "trendMaxOrder",  10)) ?   GetScalar(tmp) : (m.trendMaxOrder = 1);
			o.trendMinSepDist   = (tmp = GetField123Check(S, "trendMinSepDist", 10)) ?  GetScalar(tmp) : (m.trendMinSepDist = 1);
			o.trendMinKnotNum   = (tmp = GetField123Check(S, "trendMinKnotNum", 10)) ?  GetScalar(tmp) : (m.trendMinKnotNum = 1);			
			o.trendMaxKnotNum   = (tmp = GetField123Check(S, "trendMaxKnotNum", 10)) ?  GetScalar(tmp) : (m.trendMaxKnotNum = 1);
			
			o.trendLeftMargin  = (tmp = GetField123Check(S, "trendLeftMargin", 10)) ? GetScalar(tmp) : (m.trendLeftMargin = 1);
			o.trendRightMargin = (tmp = GetField123Check(S, "trendRightMargin", 10)) ? GetScalar(tmp) : (m.trendRightMargin = 1);

			if (io->meta.hasOutlierCmpnt) {
				o.outlierMaxKnotNum = (tmp = GetField123Check(S, "outlierMaxKnotNum", 10)) ? GetScalar(tmp) : (m.outlierMaxKnotNum = 1);
				o.outlierSigFactor = (tmp = GetFieldCheck(S, "outlierSigFactor")) ? GetScalar(tmp) : (m.outlierSigFactor = 1);
			}

			o.sigFactor         = (tmp = GetFieldCheck(S,  "sigFactor")) ?			GetScalar(tmp) : (m.sigFactor        = 1);

			o.sig2              = (tmp = GetField123Check(S, "sig2",2)) ?			GetScalar(tmp) : (m.sig2 = 1);
			o.precValue		    = (tmp = GetField123Check(S, "precValue",5)) ?		GetScalar(tmp) : (m.precValue = 1);
			o.alpha1			= (tmp = GetField123Check(S, "alpha1",0)) ?			GetScalar(tmp) : (m.alpha1 = 1);
			o.alpha2			= (tmp = GetField123Check(S, "alpha2",0)) ?			GetScalar(tmp) : (m.alpha2 = 1);
			o.delta1			= (tmp = GetField123Check(S, "delta1",0)) ?			GetScalar(tmp) : (m.delta1 = 1);
			o.delta2			= (tmp = GetField123Check(S, "delta2",0)) ?			GetScalar(tmp) : (m.delta2 = 1);

			o.K_MAX = (tmp = GetField123Check(S, "K_MAX", 1)) ? GetScalar(tmp) : (m.K_MAX = 1);

			if (io->meta.hasSeasonCmpnt)  o.seasonBasisFuncType	  = (tmp = GetField123Check(S, "seasonBasisFuncType",10)) ?  GetScalar(tmp) : (m.seasonBasisFuncType = 1);
			if (1L)                       o.trendBasisFuncType	  = (tmp = GetField123Check(S, "trendBasisFuncType", 10)) ?   GetScalar(tmp) : (m.trendBasisFuncType = 1);
			if (io->meta.hasOutlierCmpnt) o.outlierBasisFuncType  = (tmp = GetField123Check(S, "outlierBasisFuncType", 10)) ? GetScalar(tmp) : (m.outlierBasisFuncType = 1);
			o.modelPriorType		= (tmp = GetField123Check(S, "modelPriorType",  10)) ?		GetScalar(tmp) : (m.modelPriorType = 1);
			
			//o.precPriorType		    = (tmp = GetFieldCheck(S, "precPriorType")) ?	    GetScalar(tmp) : (m.precPriorType = 1);
			o.precPriorType = __GetPrecPriorType(S);

		
		}

	} // if (nrhs >= 4)

	o.numBasis       = 1L + io->meta.hasSeasonCmpnt + io->meta.hasOutlierCmpnt;
	I32  basisIdx    = 0;	
	if (io->meta.hasSeasonCmpnt) {
		I08      seasonFrom = io->meta.seasonForm;
		if      (seasonFrom == 'S')	o.basisType[basisIdx++] = SEASONID;
		else if (seasonFrom == 'D') o.basisType[basisIdx++] = DUMMYID;
		else if (seasonFrom == 'V') o.basisType[basisIdx++] = SVDID;
		else {
			r_error("ERROR: the season character is unrecognized. Valid values are 'none', 'harmonic', 'dummy', and 'svd'. \n");
			return 0;
		}
	}
	o.basisType[basisIdx++] = TRENDID;
	if (io->meta.hasOutlierCmpnt) o.basisType[basisIdx++] = OUTLIERID;
		
	/**************************/
	//Put here because the following setup need the values of period and N (time-series length).
	//These two values have been determined in ParseInputData. 
	I32 period = io->meta.period;
	I32 N      = io->N;
	/**************************/

	if (io->meta.hasSeasonCmpnt) {
		if (io->meta.seasonForm == 'S') {
			if (m.seasonMinOrder)      o.seasonMinOrder = 1L;				   o.seasonMinOrder = min(o.seasonMinOrder, period / 2 - 1);    o.seasonMinOrder = max(o.seasonMinOrder, 1L);
			if (m.seasonMaxOrder)      o.seasonMaxOrder = (period / 2 - 1);    o.seasonMaxOrder = min(o.seasonMaxOrder, (period / 2 - 1));  o.seasonMaxOrder = max(o.seasonMaxOrder, o.seasonMinOrder);
		}	else if (io->meta.seasonForm == 'V') {
			if (m.seasonMinOrder)      o.seasonMinOrder = 1L;				   o.seasonMinOrder = min(o.seasonMinOrder, period - 1);   o.seasonMinOrder = max(o.seasonMinOrder, 1L);
			if (m.seasonMaxOrder)      o.seasonMaxOrder = (period / 2 - 1);    o.seasonMaxOrder = min(o.seasonMaxOrder, period );  o.seasonMaxOrder = max(o.seasonMaxOrder, o.seasonMinOrder);
		}
		
		if (m.seasonMinSepDist || o.seasonMinSepDist <= 0)   o.seasonMinSepDist  = period / 2;         
		o.seasonMinSepDist = max(o.seasonMinSepDist, o.seasonMaxOrder);		 
		o.seasonMinSepDist = min(o.seasonMinSepDist, N / 2 - 1       ); // TODO:N/2-1 can be negtative, and then forced into a lager postive U16 integer

		if (m.seasonLeftMargin  || o.seasonLeftMargin  < 0)   o.seasonLeftMargin  = o.seasonMinSepDist;          
		if (m.seasonRightMargin || o.seasonRightMargin < 0)   o.seasonRightMargin = o.seasonMinSepDist;

		I32 Nleft = N - (1 + o.seasonLeftMargin + o.seasonRightMargin);
		if (Nleft < 1) {
			
			if (N - (1 + o.seasonMinSepDist * 2) < 1) {
				o.seasonLeftMargin  = 0;
				o.seasonRightMargin = 0;
			} else {
				o.seasonLeftMargin  = o.seasonMinSepDist;
				o.seasonRightMargin = o.seasonMinSepDist;
			}
			q_warning("WARNING: prior$seasonLeftMargin and prior$seasonRightMargin are too large, and no remaining data points are "
				      "available as potential changepoints! Their vaules are changed to %d and %d.", o.seasonLeftMargin, o.seasonRightMargin);

			Nleft = N - (1 + o.seasonLeftMargin + o.seasonRightMargin);
		}

		I32 MaxChangePointPossible = ceil((Nleft + 0.0) / (o.seasonMinSepDist + 1.0));

		if (m.seasonMinKnotNum   )   o.seasonMinKnotNum = 0;                  
		if (m.seasonMaxKnotNum==0)   o.seasonMinKnotNum =  min(o.seasonMaxKnotNum, o.seasonMinKnotNum);
		o.seasonMinKnotNum = max(o.seasonMinKnotNum, 0);
		o.seasonMinKnotNum = min(o.seasonMinKnotNum, MaxChangePointPossible);

		if (m.seasonMaxKnotNum)     o.seasonMaxKnotNum  = min(MaxChangePointPossible, 5);
		o.seasonMaxKnotNum = min(o.seasonMaxKnotNum, MaxChangePointPossible);
		o.seasonMaxKnotNum = max(o.seasonMaxKnotNum, o.seasonMinKnotNum);
	}

	{ // the trend component
		if (m.trendMinOrder || o.trendMinOrder < 0)      o.trendMinOrder    = 0L;
		if (m.trendMaxOrder)                             o.trendMaxOrder    = 1L;				                    
		o.trendMaxOrder	  = max(o.trendMaxOrder, o.trendMinOrder);	

		if (m.trendMinSepDist || o.trendMinSepDist <= 0) o.trendMinSepDist = io->meta.hasSeasonCmpnt? period / 2: 3 ;
		o.trendMinSepDist = max(o.trendMinSepDist, o.trendMaxOrder + 1);
		o.trendMinSepDist = min(o.trendMinSepDist, N / 2 - 1);

		if (m.trendLeftMargin  || o.trendLeftMargin < 0)     o.trendLeftMargin = o.trendMinSepDist;
		if (m.trendRightMargin || o.trendRightMargin < 0)    o.trendRightMargin = o.trendMinSepDist;

		I32 Nleft = N - (1 + o.trendLeftMargin + o.trendRightMargin);
		if (Nleft < 1) {

			if (N - (1 + o.trendMinSepDist * 2) < 1) {
				o.trendLeftMargin = 0;
				o.trendRightMargin = 0;
			}
			else {
				o.trendLeftMargin = o.trendMinSepDist;
				o.trendRightMargin = o.trendMinSepDist;
			}
			q_warning("WARNING: prior$trendLeftMargin and prior$trendRightMargin are too large, and no remaining data points are "
				"available as potential changepoints! Their vaules are changed to %d and %d.", o.trendRightMargin, o.trendRightMargin);

			Nleft = N - (1 + o.trendLeftMargin + o.trendRightMargin);
		}

		I32 MaxChangePointPossible = ceil((Nleft + 0.0) / (o.trendMinSepDist + 1.0));


		if (m.trendMinKnotNum)       o.trendMinKnotNum = 0;	   
		o.trendMinKnotNum = max(min(o.trendMaxKnotNum, o.trendMinKnotNum), 0);

		if (m.trendMinKnotNum)        o.trendMinKnotNum = 0;
		if (m.trendMaxKnotNum == 0)   o.trendMinKnotNum = min(o.trendMinKnotNum, o.trendMaxKnotNum);
		o.trendMinKnotNum = max(o.trendMinKnotNum, 0);
		o.trendMinKnotNum = min(o.trendMinKnotNum, MaxChangePointPossible);

		if (m.trendMaxKnotNum)     o.trendMaxKnotNum = min(MaxChangePointPossible, 10);
		o.trendMaxKnotNum = min(o.trendMaxKnotNum, MaxChangePointPossible);
		o.trendMaxKnotNum = max(o.trendMaxKnotNum, o.trendMinKnotNum);	 
	}
	
	if (m.outlierMaxKnotNum) o.outlierMaxKnotNum = o.trendMaxKnotNum;    
	o.outlierMaxKnotNum = max(o.outlierMaxKnotNum, 1L); // at least has one; otherwise, the program crahses if hasOUtliercomponet=1
	
	if (m.K_MAX )            o.K_MAX            = 0;    // exact value of K_max to be determined later if Kmax = 0              
	if (m.sigFactor)         o.sigFactor        = 1.8;            o.sigFactor        = max(o.sigFactor,        1.02);
	if (m.outlierSigFactor)  o.outlierSigFactor = 2.5;            o.outlierSigFactor = max(o.outlierSigFactor, 1.5);

	if (m.sig2 )             o.sig2      = 0.2f;				  o.sig2             = max(o.sig2,      0.01);
	if (m.precValue)         o.precValue = 1.5f;				  o.precValue        = max(o.precValue, 0.01);
	if (m.alpha1)		     o.alpha1	 = 0.00000001f;
	if (m.alpha2)		     o.alpha2	 = 0.00000001f;
	if (m.delta1)		     o.delta1	 = 0.00000001f;
	if (m.delta2)		     o.delta2	 = 0.00000001f;

	if (m.precPriorType)			o.precPriorType = UniformPrec;


	if (m.seasonBasisFuncType) {
		if      (o.precPriorType == UniformPrec)		o.seasonBasisFuncType = 0;
		else if (o.precPriorType == ConstPrec)          o.seasonBasisFuncType = 0;
		else if (o.precPriorType == ComponentWise)      o.seasonBasisFuncType = 1;
		else if (o.precPriorType == OrderWise)          o.seasonBasisFuncType = 1;		
	}
	if (m.trendBasisFuncType) {
		if      (o.precPriorType == UniformPrec)		o.trendBasisFuncType = 0;
		else if (o.precPriorType == ConstPrec)          o.trendBasisFuncType = 0;
		else if (o.precPriorType == ComponentWise)      o.seasonBasisFuncType = 1;
		else if (o.precPriorType == OrderWise)          o.trendBasisFuncType = 1;
	}	 
	if (m.outlierBasisFuncType) {
		if      (o.precPriorType == UniformPrec)		o.outlierBasisFuncType = 0;
		else if (o.precPriorType == ConstPrec)          o.outlierBasisFuncType = 0;
		else if (o.precPriorType == ComponentWise)      o.outlierBasisFuncType = 1;
		else if (o.precPriorType == OrderWise)          o.outlierBasisFuncType = 1;
	}	 
	if (m.modelPriorType)			o.modelPriorType        = 1L;

	return 1;

#undef o
}

static int  GetArg_3rd_MCMC___(VOIDPTR prhs[], int nrhs, BEAST2_MCMC_PTR mcmc,  BEAST2_OPTIONS_PTR opt)
{
	// Before running this function, trendMinSepDist must have been set.

	#define o (*mcmc)
	struct MCMC_MISSING {
		U08   seed;	                  // Unsigned long long seed;
		U08   credIntervalAlphaLevel;
		U08   trendResamplingOrderProb;
		U08   seasonResamplingOrderProb;
		U08   ridgeFactor;
		U08   burnin, samples, chainNumber;
		U08   maxMoveStepSize;
		U08   thinningFactor;
	} m = {0,};


	if (nrhs < 5) 
		memset(&m, 1L, sizeof(struct MCMC_MISSING));

	if (nrhs >= 5) {

		VOIDPTR S = prhs[4L];
		if (!IsStruct(S)) {
			q_warning("WARNING: The arg 'mcmc' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct MCMC_MISSING));
		} else {
			VOIDPTR tmp;
			o.maxMoveStepSize = (tmp = GetField123Check(S, "maxMoveStepSize",2))? GetScalar(tmp) : (m.maxMoveStepSize = 1);
			o.samples         = (tmp = GetField123Check(S, "samples", 3)) ?        GetScalar(tmp) : (m.samples = 1);
			o.thinningFactor  = (tmp = GetField123Check(S, "thinningFactor", 2)) ? GetScalar(tmp) : (m.thinningFactor = 1);
			o.burnin          = (tmp = GetField123Check(S, "burnin", 2)) ?         GetScalar(tmp) : (m.burnin = 1);
			o.chainNumber     = (tmp = GetField123Check(S, "chainNumber", 2)) ?    GetScalar(tmp) : (m.chainNumber = 1);
			o.seed			  = (tmp = GetField123Check(S, "seed", 4)) ?			GetScalar(tmp) : (m.seed = 1);
			o.ridgeFactor	  = (tmp = GetField123Check(S, "ridgeFactor", 2)) ?	GetScalar(tmp) : (m.ridgeFactor = 1);

			o.trendResamplingOrderProb  = (tmp = GetField123Check(S, "trendResamplingOrderProb",  5)) ? GetScalar(tmp) : (m.trendResamplingOrderProb = 1);
			o.seasonResamplingOrderProb = (tmp = GetField123Check(S, "seasonResamplingOrderProb", 5)) ? GetScalar(tmp) : (m.seasonResamplingOrderProb = 1);
			o.credIntervalAlphaLevel    = (tmp = GetField123Check(S, "credIntervalAlphaLevel",    2)) ? GetScalar(tmp) : (m.credIntervalAlphaLevel = 1);
		}

	} // if (nrhs >= 5)
 
	//r_printf("move :%d  %d %d\n", o.maxMoveStepSize, opt->io.meta.hasSeasonCmpnt , opt->prior.trendMinSepDist  );
	if (m.maxMoveStepSize || o.maxMoveStepSize==0) o.maxMoveStepSize = opt->io.meta.hasSeasonCmpnt? opt->io.meta.period: (opt->prior.trendMinSepDist + 1);
	if (m.samples         || o.samples==0)		   o.samples		 = 3000;        o.samples		  = max(o.samples, 800);
	if (m.thinningFactor  || o.thinningFactor==0)  o.thinningFactor  = 1L;			o.thinningFactor = max(o.thinningFactor, 1L);
	if (m.burnin          || o.burnin==0)          o.burnin		     = 150L;		o.burnin		  = max(o.burnin, 150L);
	if (m.chainNumber || o.chainNumber==0)         o.chainNumber	 = 3;			o.chainNumber	  = max(o.chainNumber, 1L);
	if (m.seed)            o.seed			 = 0L;
	if (m.credIntervalAlphaLevel)      o.credIntervalAlphaLevel	 = .95;
	if (m.ridgeFactor)     o.ridgeFactor	 = 0.0001f;

	if (m.trendResamplingOrderProb)  o.trendResamplingOrderProb  = .1f;
	if (m.seasonResamplingOrderProb) o.seasonResamplingOrderProb = .17f;
 
	return 1;

#undef o
}

static int  GetArg_4th_EXTRA__(VOIDPTR prhs[], int nrhs, BEAST2_EXTRA_PTR extra, I32 whichDimIsTime, I32 ndims)
{
	// Before running this function, meta.whichDimIsTime must be first obtained

	#define o (*extra)
	struct OUTFLAGS_MISSING {
		I08   numThreadsPerCPU;
		I08   numParThreads;
		I08   numCPUCoresToUse;		

		I08   consoleWidth;
		I08   whichOutputDimIsTime;
		I08   removeSingletonDims;
		I08   dumpInputData;

		I08   ncpStatMethod;
		I08  smoothCpOccPrCurve;
		I08  useMeanOrRndBeta;
		I08  computeCredible;
		I08  fastCIComputation;
		I08  computeSeasonOrder;
		I08  computeTrendOrder;

		I08  computeSeasonChngpt;
		I08  computeTrendChngpt;
		I08  computeOutlierChngpt;

		I08  computeSeasonAmp;
		I08  computeTrendSlope;

		I08 tallyPosNegSeasonJump;
		I08 tallyPosNegTrendJump;
		I08 tallyIncDecTrendJump;
		I08 tallyPosNegOutliers;

		I08  printOptions;
		I08  printProgressBar;
	} m = {0,};


	if (nrhs < 6) 
		memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));	

	if (nrhs >= 6) {
		VOIDPTR S = prhs[5L];
		if (!IsStruct(S)) {
			q_warning("WARNING: The arg 'extra' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct OUTFLAGS_MISSING));
		}
		else {
			VOIDPTR tmp;
			o.whichOutputDimIsTime = (tmp = GetField123Check(S, "whichOutputDimIsTime",2)) ?	GetScalar(tmp) : (m.whichOutputDimIsTime = 1);
			o.removeSingletonDims  = (tmp = GetField123Check(S, "removeSingletonDims", 8)) ? GetScalar(tmp) : (m.removeSingletonDims = 1);			
			o.numThreadsPerCPU     = (tmp = GetField123Check(S, "numThreadsPerCPU", 4)) ? GetScalar(tmp) : (m.numThreadsPerCPU = 1);
			o.numParThreads        = (tmp = GetField123Check(S, "numParThreads", 4)) ?			GetScalar(tmp) : (m.numParThreads = 1);
			o.numCPUCoresToUse     = (tmp = GetField123Check(S, "numCPUCoresToUse", 4)) ?		GetScalar(tmp) : (m.numCPUCoresToUse = 1);
			o.consoleWidth         = (tmp = GetField123Check(S, "consoleWidth",2)) ?			GetScalar(tmp) : (m.consoleWidth = 1);
			o.dumpInputData        = (tmp = GetField123Check(S, "dumpInputData",2)) ? GetScalar(tmp) : (m.dumpInputData = 1);
			o.smoothCpOccPrCurve   = (tmp = GetField123Check(S, "smoothCpOccPrCurve",2)) ? GetScalar(tmp) : (m.smoothCpOccPrCurve = 1);
			#define _1(x)       o.x = (tmp=GetFieldCheck(S,#x))? GetScalar(tmp): (m.x=1)
			#define _2(x,y)     _1(x);_1(y)
			#define _3(x,y,z)   _1(x);_2(y,z)
			#define _4(x,y,z,w) _2(x,y);_2(z,w)
			
			_2(printProgressBar, printOptions);
			_2(computeCredible,  fastCIComputation);

			_2(computeSeasonOrder,   computeTrendOrder);
			_3(computeSeasonChngpt,  computeTrendChngpt, computeOutlierChngpt);
			_2(computeSeasonAmp,     computeTrendSlope);
			_4(tallyPosNegSeasonJump, tallyPosNegTrendJump, tallyIncDecTrendJump, tallyPosNegOutliers);
			_1(useMeanOrRndBeta);
				 
		} // if (!IsStruct(S)) : S is a struct
	} // if (nrhs >= 5)

	if (m.whichOutputDimIsTime)		o.whichOutputDimIsTime = whichDimIsTime;
	if (m.removeSingletonDims)		o.removeSingletonDims  = 1;

	// the output dim shouldn'tY be larger than the input dims
	if (o.whichOutputDimIsTime > ndims)		o.whichOutputDimIsTime = 1L; 
	 
	if (m.smoothCpOccPrCurve)    o.smoothCpOccPrCurve = 0;
	if (m.dumpInputData)		 o.dumpInputData       = 0;
	if (m.numThreadsPerCPU)      o.numThreadsPerCPU    = 2;
	if (m.numParThreads)         o.numParThreads		= 0;
	if (m.numCPUCoresToUse)      o.numCPUCoresToUse	= 0;	
	if (m.consoleWidth||o.consoleWidth<=0)  o.consoleWidth= GetConsoleWidth(); 	o.consoleWidth = max(o.consoleWidth, 40);
	if (m.printProgressBar)      o.printProgressBar	= 1;
	if (m.printOptions)          o.printOptions		= 1;

	if (m.computeCredible)       o.computeCredible		= 0L;
	if (m.fastCIComputation)     o.fastCIComputation	= 1L;

	if (m.computeSeasonOrder)    o.computeSeasonOrder	= 0L;
	if (m.computeTrendOrder)     o.computeTrendOrder	= 0L;

	if (m.computeSeasonChngpt)   o.computeSeasonChngpt	= 1L;
	if (m.computeTrendChngpt)    o.computeTrendChngpt	= 1L;
	if (m.computeOutlierChngpt)  o.computeOutlierChngpt = 1L;

	if (m.computeSeasonAmp)       o.computeSeasonAmp	= 0L;
	if (m.computeTrendSlope)      o.computeTrendSlope	= 0L;

	if (m.tallyPosNegSeasonJump)     o.tallyPosNegSeasonJump	= 0L;
	if (m.tallyPosNegTrendJump)      o.tallyPosNegTrendJump		= 0L;
	if (m.tallyIncDecTrendJump)      o.tallyIncDecTrendJump		= 0L;
	if (m.tallyPosNegOutliers)       o.tallyPosNegOutliers		= 0L;

	if (o.tallyPosNegSeasonJump) o.computeSeasonChngpt = 1,	o.computeSeasonAmp	= 1;
	if (o.tallyPosNegTrendJump)  o.computeTrendChngpt	= 1,	o.computeTrendSlope = 1;
	if (o.tallyIncDecTrendJump)  o.computeTrendChngpt	= 1,	o.computeTrendSlope = 1;
	if (o.tallyPosNegOutliers)   o.computeOutlierChngpt = 1;

	if (m.useMeanOrRndBeta)      o.useMeanOrRndBeta = 0;


	
	return 1;
#undef o
}

I32 PostCheckArgs(A(OPTIONS_PTR) opt) {

	I08 hasSeasonCmpnt   = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasHarmonicCmpnt = opt->prior.basisType[0] == SEASONID ;
	I08 hasDummyCmpnt    = opt->prior.basisType[0] == DUMMYID;
	I08 hasSVDCmpnt      = opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt  = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt    = 1;
	I08 hasAlways        = 1;

	// Period must be an integer when Dummy is used
	if (hasDummyCmpnt) opt->io.meta.period = ceil(opt->io.meta.period);

	// If the max and  orders are equal, do not re-sample the order
	if (opt->prior.trendMaxOrder  == opt->prior.trendMinOrder)	opt->mcmc.trendResamplingOrderProb = 0;
	if (opt->prior.seasonMaxOrder == opt->prior.seasonMinOrder)	opt->mcmc.seasonResamplingOrderProb = 0;
	if (hasDummyCmpnt)											opt->mcmc.seasonResamplingOrderProb = 0;


	if (opt->io.meta.deseasonalize && !hasSeasonCmpnt) {
		opt->io.meta.deseasonalize = 0;		
	}
	// the data contains only 1 period, and no enough data to fit the seasonal compnt
	if(opt->io.meta.period+3 /* 1:const, 1:slope, 1:extra degree of freedom*/ > opt->io.N) {
		opt->io.meta.deseasonalize = 0;
	}

	if (!hasSeasonCmpnt) {
		opt->extra.computeSeasonOrder    = 0;
		opt->extra.computeSeasonChngpt   = 0;
		opt->extra.computeSeasonAmp      = 0;
		opt->extra.tallyPosNegSeasonJump = 0;
	}
	if (hasDummyCmpnt) {
		opt->extra.computeSeasonOrder	    = 0;	
		opt->extra.computeSeasonAmp         = 0; //TODO: remove this restriction
		opt->prior.seasonMinOrder           = 0;
		opt->prior.seasonMaxOrder           = 0;		
	}
	if (hasSVDCmpnt) { 
		opt->extra.computeSeasonAmp = 0; //TODO: remove this restriction 
	}
	if (!hasOutlierCmpnt) {
		opt->extra.computeOutlierChngpt = 0;
		opt->extra.tallyPosNegOutliers  = 0;
	}

	BEAST2_PRIOR_PTR PRIOR = &opt->prior;
	#define prior (*PRIOR)
	I08  isTrendCmpntFixed  = prior.trendMaxOrder  == prior.trendMinOrder  && prior.trendMaxKnotNum  == 0 && prior.trendMinKnotNum == 0;
	I08  isSeasonCmpntFixed = prior.seasonMaxOrder == prior.seasonMinOrder && prior.seasonMaxKnotNum == 0 && prior.seasonMinKnotNum == 0;
	if (hasDummyCmpnt) isSeasonCmpntFixed = prior.seasonMaxKnotNum == 0 && prior.seasonMinKnotNum == 0;
	#undef prior

	//if numbasis == 2 and hasOutlier == 1, then the first cmpt must not be fixed
	if (opt->prior.numBasis == 2 ) {
		if (hasTrendCmpnt  && hasOutlierCmpnt && isTrendCmpntFixed) {
			q_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0' will"
				     " fix the trend to a global curve.\n");					 
		}
		if (hasHarmonicCmpnt && hasOutlierCmpnt && isSeasonCmpntFixed) {
			q_warning("WARNING: The options 'seasonMaxOrder==seasonMaxOrder && seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
				    " fix the season component to a global curve.\n");		 
		}
		if (hasDummyCmpnt && hasOutlierCmpnt && isSeasonCmpntFixed) {
			// for the dummy cmpnt, maxOrder and minOrder are both set to zeros.
			q_warning("WARNING: The options 'seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
				    " fix the dummy season component to a global curve.\n");			
		}
	}
	if (hasTrendCmpnt && hasDummyCmpnt & isTrendCmpntFixed && isSeasonCmpntFixed) {		
		q_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0 && "
			    " seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
			    " fix the model structures of the trend and dummy season components.\n");
 
	}
	if (hasTrendCmpnt && hasHarmonicCmpnt & isTrendCmpntFixed && isSeasonCmpntFixed) {
		q_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0 && "
			    "seasonMaxOrder==seasonMaxOrder && seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
			    " fix the model structures of the trend and harmonic season components.\n"); 
	}

	I32 KMAX=0;
	for (I32 i = 0; i < PRIOR->numBasis; i++) {		
		I08 type = PRIOR->basisType[i];
		if (type == SEASONID)			KMAX += (PRIOR->seasonMaxOrder*2) * (PRIOR->seasonMaxKnotNum + 1);
		if (type == DUMMYID)			KMAX += (opt->io.meta.period)     * (PRIOR->seasonMaxKnotNum + 1);
		if (type == SVDID)			    KMAX += (opt->io.meta.period)     * (PRIOR->seasonMaxKnotNum + 1);
		if (type == TRENDID)			KMAX += (PRIOR->trendMaxOrder +1) * (PRIOR->trendMaxKnotNum + 1);
		if (type == OUTLIERID)			KMAX += PRIOR->outlierMaxKnotNum;
	}
	if (PRIOR->K_MAX <= 0) {  // Not set yet, so takt the max number of terms as dertermined by the knots
		PRIOR->K_MAX =  KMAX;
	} else {
		PRIOR->K_MAX = min(PRIOR->K_MAX, KMAX);
	}
	
	if (opt->io.N < 5000)
		PRIOR->K_MAX = min(PRIOR->K_MAX, opt->io.N);
	else if (opt->io.N < 15000) {
		int KMAX = min(5000, opt->io.N / 2);
		PRIOR->K_MAX = min(PRIOR->K_MAX, KMAX);
	}
	else  {
		int KMAX = 7500;
		PRIOR->K_MAX = min(PRIOR->K_MAX, KMAX);
	}
	
 

	// The initial model is generated in "basis_genrandombasis". We cacluat the number of terms for it
	// so KMAX must be larger than it
	I32 K_INITIAL_MODEL = 0;
	for (I32 i = 0; i < PRIOR->numBasis; i++) {
		I08 type = PRIOR->basisType[i];
		if (type == SEASONID) {
			I32 order    = ceil((PRIOR->seasonMaxOrder   + PRIOR->seasonMinOrder) / 2.0);
			I32 nKnot    = ceil((PRIOR->seasonMinKnotNum + PRIOR->seasonMaxKnotNum) / 2.0);   //nKnot is zero if maxKnotNum=0
			K_INITIAL_MODEL += (order * 2) * (nKnot + 1);
		}
		if (type == DUMMYID) {
			I32 order   =  ceil(opt->io.meta.period);
			I32 nKnot   =  ceil((PRIOR->seasonMinKnotNum + PRIOR->seasonMaxKnotNum) / 2.0);   //nKnot is zero if maxKnotNum=0
			K_INITIAL_MODEL += (order ) * (nKnot + 1);
		}
		if (type == SVDID) {
			I32 order   = ceil(opt->io.meta.period);
			I32 nKnot   = ceil((PRIOR->seasonMinKnotNum + PRIOR->seasonMaxKnotNum) / 2.0);   //nKnot is zero if maxKnotNum=0
			K_INITIAL_MODEL += (order) * (nKnot + 1);
		}
		if (type == TRENDID)		{
			I32 order    = ceil((PRIOR->trendMaxOrder   + PRIOR->trendMinOrder) / 2.0);
			I32 nKnot    = ceil((PRIOR->trendMinKnotNum + PRIOR->trendMaxKnotNum) / 2.0);   //nKnot is zero if maxKnotNum=0
			K_INITIAL_MODEL += (order +1L) * (nKnot + 1);
		}
		if (type == OUTLIERID) {
			K_INITIAL_MODEL += 0;
		}
	}
	PRIOR->K_MAX = max(PRIOR->K_MAX, K_INITIAL_MODEL);


	opt->io.out.dtype          = IsRinterface() ? DATA_DOUBLE : DATA_FLOAT;
	opt->io.out.whichDimIsTime = opt->extra.whichOutputDimIsTime;

    // no need to use componetwise if there is only one componet
	if (opt->prior.precPriorType == ComponentWise && opt->prior.numBasis == 1) {
		opt->prior.precPriorType = UniformPrec;
		//q_warning("WARNING: prior$precPriorType is changed from 'componentwise' to 'uniform' because the model specified only has a trend component.\n");
	}
	return 1;
}

int BEAST2_GetArgs(VOIDPTR prhs[], int nrhs, A(OPTIONS_PTR) opt) {

  

	int  failed = !GetArg_0th_Data(prhs, nrhs, &opt->io)                 ||
		          !GetArg_1st_MetaData(prhs, nrhs, &opt->io)		      || 				  
			      !GetArg_2nd_Prior__(prhs, nrhs, &opt->prior, &opt->io)  ||
			      !GetArg_3rd_MCMC___(prhs, nrhs, &opt->mcmc,  opt)       ||
			      !GetArg_4th_EXTRA__(prhs, nrhs, &opt->extra, opt->io.meta.whichDimIsTime,opt->io.ndim) ;
	int success = !failed;	
	if (success) 	success=PostCheckArgs(opt); 	
	if (success) 	BEAST2_print_options(opt);	

	return success;
}


void BEAST2_DeallocateTimeSeriesIO(BEAST2_IO_PTR  o) {
	// Free mems allocated in the "tsAggegrationPrepare" function
	if (o->T.out.numPtsPerInterval != NULL) {
		free(o->T.out.numPtsPerInterval);
		o->T.out.numPtsPerInterval = NULL;
	}
	TimeVec_kill(&o->T);

	if (o->pdata != NULL) {
		free(o->pdata);
		o->pdata = NULL;
	}
	if (o->dtype != NULL) {
		free(o->dtype);
		o->dtype = NULL;
	}
	if (o->out.result != NULL) {
		free(o->out.result);
		o->out.result = NULL;
	}
 
}

#include "abc_000_warning.h"