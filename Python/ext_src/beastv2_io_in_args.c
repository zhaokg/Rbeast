#include "abc_000_warning.h"

#include "abc_001_config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_ide_util.h"  //printf
#include "abc_common.h"  //strcicmp
#include "abc_ts_func.h"
#include "abc_date.h"
#include "beastv2_func.h"    
#include "beastv2_io.h"

#define IsAlmostInteger(x)  ( fabs(x-round(x)) <1e-3 )
static int  GetArg_1st_MetaData(VOIDPTR prhs[], int nrhs, BEAST2_METADATA_PTR meta)
{
	
	// Check the number of inputs:  the first arg is the algorithm name.
	if (nrhs < 2L) {	
		r_error("ERROR: At least one input argument is needed!\n");
		return 0;
	}

	meta->rawInput     = prhs[1];	
	
 

	I08 metaProcessed = 0;
	if (nrhs < 3) {
		meta->isMetaStruct		= 0;
		meta->isRegularOrdered  = 1; //We asumme the data is regular if no meta is supplied
		meta->hasSeasonCmpnt    = 1;
		meta->hasOutlierCmpnt   = 0;
		meta->seasonForm        ='S';
		meta->detrend           = 0;
		meta->deseasonalize     = 0;
		meta->period			= getNaN();

		meta->startTime = 1.f;
		meta->deltaTime = 1.f;	
		meta->missingValue = getNaN();// FLOAT_MAX;
		meta->maxMissingRate = 0.75;
		meta->whichDimIsTime = -1; //Not used
		meta->rawTimeVec      = NULL;

		metaProcessed = 1;
	}

	if (nrhs >= 3 && IsNumeric(prhs[2L])) {

		VOIDPTR pmeta = prhs[2L];

		meta->isMetaStruct = 0;
		meta->isRegularOrdered = 1L; //We asumme the data is regular if meta is a scalar numeric value.
		meta->hasSeasonCmpnt   = 1L;
		meta->hasOutlierCmpnt  = 0L;
		meta->seasonForm        ='S';
		meta->detrend           = 0;
		meta->deseasonalize     = 0;

		meta->period    = GetScalar(pmeta);
		meta->startTime = 1.f;
		meta->deltaTime = 1.f;
		meta->missingValue = getNaN();// FLOAT_MAX;
		meta->maxMissingRate = 0.75;
		meta->whichDimIsTime = -1; //Not used
		meta->rawTimeVec     = NULL;

		metaProcessed = 1;
	}

	if (nrhs >= 3 && IsStruct(prhs[2L])  ) {

		meta->isMetaStruct = 1;
	
		VOIDPTR pmeta = prhs[2L];
		VOIDPTR tmp;

		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->isRegularOrdered = (tmp = GetField123Check(pmeta, "isRegularOrdered",2)) ? GetScalar(tmp) : 1L;
		if (tmp == NULL)
			r_warning("WARNING: metadata$isRegularOrdered is not specified. A default 'metadata$isRegularOrdered=TRUE' is assumed.\n");

		//////////////////////////////////////////////////////////////////////////////////////////////
		// Get the season string to determine the number and types of season components: harnomic or dummy
		
		
		if ((tmp = GetField123Check(pmeta, "season",2)) != NULL && IsChar(tmp)) {			
			char season[20 + 1];
			GetCharArray(tmp, season, 20);
			ToUpper(season);

			char a = season[0], b = season[1];
			if      (a=='N' && b=='O') 		meta->hasSeasonCmpnt = 0;							    //none
			else if (a=='H' && b=='A') 		meta->hasSeasonCmpnt = 1, meta->seasonForm = 'S'; 	//harmonic
			else if (a=='D' && b=='U')		meta->hasSeasonCmpnt = 1, meta->seasonForm = 'D';			//dummy
			else if (a=='S' && b=='V')		meta->hasSeasonCmpnt = 1, meta->seasonForm = 'V';			//svd			
			else {
				meta->hasSeasonCmpnt = 1;
				meta->seasonForm     = 'S';  //the default form is harmonic
				r_warning("WARNING: metadata$season=%s is specified but has an unrecongizable string. A default metadata$season='harmonic' is used instead.\n", season);
			}

			if (meta->seasonForm == 'V') meta->svdTerms = GetFieldCheck(pmeta, "svdTerms");

		} else  {
			meta->hasSeasonCmpnt = 1;
			meta->seasonForm     = 'S';
			r_warning("WARNING: metadata$season is either missing or not given as a valid specifier string (e.g., none, harmonic, or dummy). A default "
				      "metadata$season='harmonic' is assumed.\n");
		}	
		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->hasOutlierCmpnt = (tmp = GetField123Check(pmeta, "hasOutlierCmpnt",2)) ? GetScalar(tmp) : 0L;

		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->detrend         = (tmp = GetField123Check(pmeta, "detrend"      ,3)) ?  GetScalar(tmp) : 0L;
		meta->deseasonalize   = (tmp = GetField123Check(pmeta, "deseasonalize",3)) ?  GetScalar(tmp) : 0L;
 
		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->period = (tmp = GetField123Check(pmeta, "period",2)) ? GetScalar(tmp) : getNaN();
		if (!meta->hasSeasonCmpnt && tmp != NULL)
			r_warning("WARNING: For metadata$season='none' (i.e., no periodic component in the time series), metadata$period is ignored!\n");

		//////////////////////////////////////////////////////////////////////////////////////////////
		I08 isDefaultStartTime = 0;
		tmp = GetField123Check(pmeta, "startTime",2);
		if     (tmp && IsClass(tmp,"Date")) 		{
			// this branch is true only for R and Matlab's IsClass always return FALSE
			// in R, tmp is an integer array
			int   days      = GetScalar(tmp);
			meta->startTime = fractional_civil_from_days(days);
		}
		else if (tmp && IsNumeric(tmp)) {
			I32 n = GetNumberOfElements(tmp);
			if (n==1)
				meta->startTime = GetScalar(tmp);
			else if (n == 2) {
				F32 Y = GetNumericElement(tmp, 0);
				F32 M = GetNumericElement(tmp, 1);
				meta->startTime = YMDtoF32time(Y, M, 15);
				r_warning("WARNING: Your metadata$startTime (%f, %f) is interpreted as (Year: %d, Month: %d, Day: 15) and converted to "
					       "a numeric value of %f. If not making sense, make sure to supply a correct start time: "
					       "startTime can take a numeric scalar, a vector of two values (year, month), or a vector of three values (year, month,day).",
					       Y,M, (int)Y, (int)M, meta->startTime);
			}
			else if (n == 3) {
				F32 Y = GetNumericElement(tmp, 0);
				F32 M = GetNumericElement(tmp, 1);
				F32 D = GetNumericElement(tmp, 2);
				meta->startTime = (F32) YMDtoF32time((int)Y, (int)M, (int)D);
				r_warning("WARNING: Your metadata$startTime (%f, %f, %f) is interpreted as (Year: %d, Month: %d, Day: %d) and converted to "
					"a numeric value of %f. If not making sense, make sure to supply a correct start time: "
					"startTime can take a numeric scalar, a vector of two values (year, month), or a vector of three values (year, month,day).",
					Y, M,D, (int)Y, (int)M, (int)D, meta->startTime);
			}
			else { 
				r_error("ERROR: Your metadata$startTime is a vector of more than three elements. A valid 'startTime' should be either a numeric scalar,"
					     "a vector of two values (year, month), or a vector of three values (year, month,day)."	);
				return 0;
			}
			
		} else {
			meta->startTime = 1L;
			isDefaultStartTime = 1L;
		}

		////////////////////////////P//////////////////////////////////////////////////////////////////
		I08 isDefaultDeltaTime = 0;
		meta->deltaTime = (tmp = GetField123Check(pmeta, "deltaTime",3)) ? GetScalar(tmp) : getNaN();			
		if (!meta->isRegularOrdered && IsNaN(meta->deltaTime) ) {
			r_error("ERROR: when metadata$isRegualrOrdered=FALSE, the input data is considered irregular/unordered in time AND "
			 	   "metadata$deltaTime must be specified for BEAST to pre-process and aggregate the irregular input to a regular time series spaced at deltaTime.\n");
			return 0;
		} 
		if (meta->deltaTime <= 0  ) {
			r_error("ERROR: metadata$deltaTime must be a positive time interval!\n");
			return 0;
		}
		if (meta->isRegularOrdered && IsNaN(meta->deltaTime)) {
			meta->deltaTime    = 1.f;
			isDefaultDeltaTime = 1;		
		}
		
		if (meta->isRegularOrdered ) {
			if (isDefaultDeltaTime && isDefaultStartTime) {
				r_warning("WARNING: when metadata$isRegualrOrdered=TRUE, the input data is assumed to be regular and ordered in time AND the times of "
					"individual datapoints are determined fully by 'metadata$startTime' and 'metadata$deltaTime'. But metadata$startTime and deltaTime "
					"are missing and a default value 1 is used for both!\n");
			}
			else if (isDefaultStartTime) {
				r_warning("WARNING: when metadata$isRegualrOrdered=TRUE, the input data is assumed to be regular and ordered in time AND the times of "
					"individual datapoints are determined fully by 'metadata$startTime' and 'metadata$deltaTime'. But startTime is missing "
					"and a default value 1 is used!\n");
			}
			else if (isDefaultDeltaTime) {
				r_warning("WARNING: when metadata$isRegualrOrdered=TRUE, the input data is assumed to be regular and ordered in time AND the times of "
					"individual datapoints are determined fully by 'metadata$startTime' and 'metadata$deltaTime'. But deltaTime is missing and a default"
					" value 1 is used!\n");	
			}			
		}
		

 

		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->rawTimeVec = (tmp = GetField123Check(pmeta, "time",2)) ? tmp : NULL;
		if (meta->isRegularOrdered && tmp)
			r_warning("WARNING: metadata$time is specified but ignored. For regular ordered time series, 'metadata$startTime=%f' and 'metadata$deltaTime=%f' "
					  "alone are enough to determine the times for all data points. \n", meta->startTime, meta->deltaTime);
		if (!meta->isRegularOrdered && tmp == NULL) {
			r_error("ERROR: metadata$isRegualrOrdered=0 (false) indicates that the input data is irregular or unordered in time; 'metadata$time' must be "
				     "also supplied to specify the times at which individual data points are collected.\n");
			return 0;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////
		meta->missingValue			= (tmp = GetField123Check(pmeta, "missingValue",2))   ? GetScalar(tmp) : getNaN();// FLOAT_MAX;
		meta->maxMissingRate 		= (tmp = GetField123Check(pmeta, "maxMissingRate",2)) ? GetScalar(tmp) : 0.75;
		meta->whichDimIsTime		= (tmp = GetField123Check(pmeta, "whichDimIsTime",2)) ? GetScalar(tmp) : -1;

		metaProcessed = 1;

	} //else if (nrhs >= 3)

	if (metaProcessed == 0 ) {
		r_error("ERROR: The 'metadata' parameter given is of unsupported type.\n");
		return 0;
	}

	/*************************************************/
	// Convert period to the unit of deltaTime
	/*************************************************/
	if (meta->hasSeasonCmpnt) {	
		meta->period = meta->period / meta->deltaTime;
		F32 period   = meta->period;
		if (period == 0) {
			r_error("ERROR: Your input parameters are \"metadata$period=0\" and \"metadata$season='%s'\". BEAST can't handle a time series with no perodicity. "
				"If you mean to handle a trend-only time series, please use the trend-only version of the BEAST algorithm by specifying metadata$season='none'."
				, meta->seasonForm == 'S' ? "harmonic" : "dummy" );
			return 0;
		}

		if (meta->seasonForm == 'D' && !IsNaN(period) && !IsAlmostInteger(period)  ) {
			r_error("ERROR: Your input parameters are \"metadata$period=%f\" and \"metadata$season='%s'\". For a dummy seasonal component, "
				"period must be a multiple of deltaTime by an INTEGER number. Your period/deltaTime is %f.",	
				 period*meta->deltaTime, "dummy", period );
			return 0;
		}
	
	}

	// Up to this point, if period=NAN, it will be further determined in ParseInputData



	meta->nrhs = nrhs;
	return 1;
}

static INLINE I32 __GetRawTimeDimension(BEAST2_IO_PTR io) {
	 return io->dims[io->meta.whichDimIsTime - 1];
}

static int  ___PraseMetaData_RegularTS_Dim123(A(METADATA_PTR) meta, BEAST2_IO_PTR io) {
	/*  If the input is regular time series */ 
	if (meta->nrhs < 3 && io->ndim >=2 ) {
		// No metadata input
		r_error("ERROR: For a 2D matrix or 3D stack input, 'metadata$whichDimIsTime' must be given to specify which dim of the input refers to time.\n");
		return 0;
	} 
	if (meta->isMetaStruct==0 && io->ndim >= 2) {
		r_error("ERROR: For a 2D matrix or 3D stack input, the 'metadata' argument must be a LIST (for R) or STRUCT (for Matlab) variable "
			    "(e.g., metedata$period' and 'metadata$whichDimIsTime) providing additional info on the 2D or 3D input.\n");
		return 0;
	}
	io->N = __GetRawTimeDimension(io);
	return 1L;
}



int __ReadRawTime(F32PTR time, VOID_PTR timeField, const I32 Nraw) {		
	// time should have been pre-allocated, with a length of Nraw

	// TimeField is not a struct variable
	if ( IsStruct(timeField)==0 ) {	

		if (IsNumeric(timeField)==0) {
			r_error("ERROR: metadata$time is not numeric. If times are strings, please use metadata$time$dateStr and metadata$time$strFmt to specify data observation times.\n");
			return 0;
		} else {
			I32  Ntime = GetNumberOfElements(timeField);
			if ( Ntime != Nraw ) {
				r_error("ERROR: 'metadata$time' (%d) must be of the same length as the time dim of the input data (%d).\n", Ntime, Nraw);
				return 0;			
			}
					
			if ( IsClass(timeField, "Date")) {
				// this branch is true only for R and Matlab's IsClass always return FALSE
				// in R, tmp is an integer array
				F64PTR days = GetData(timeField);
				for (int i = 0; i < Ntime; ++i) {
					time[i] = fractional_civil_from_days((int)days[i]);
				}	
				return 1L;
			}
			else {
				if (CopyNumericArrToF32Arr(time, timeField, Ntime))
					return 1L;
				else {
					r_error("ERROR: metadata$time has an unsupported data format or type.\n");
					return 0;
				}			
			}


				
		}
		
	}

   // TimeField is a struct variable
	VOIDPTR yr  = GetField123Check(timeField, "year",1);
	VOIDPTR mn  = GetField123Check(timeField, "month",1);
	VOIDPTR day = GetField123Check(timeField, "day",3);
	VOIDPTR doy = GetField123Check(timeField, "doy",3);

	int isTimeProcessed = 0;
	if (!isTimeProcessed && yr && mn && IsNumeric(yr) && IsNumeric(mn)  && GetNumberOfElements(yr)==Nraw && GetNumberOfElements(mn) == Nraw) 	{
	
		F32PTR yr32  = time;
		F32PTR mn32  = malloc(sizeof(F32) *2*Nraw); // one for mn32 and another for day32
		F32PTR day32 = mn32+Nraw;

		if (!CopyNumericArrToF32Arr(yr32, yr, Nraw)) {
			r_error("ERROR: metadata$time$year has an unsupported data format or type.\n");
			free(mn32);
			return 0;
		}
		if (!CopyNumericArrToF32Arr(mn32, mn, Nraw)) {
			r_error("ERROR: metadata$time$month has an unsupported data format or type.\n");
			free(mn32);
			return 0;
		}
 
		if (day && IsNumeric(day) && GetNumberOfElements(day) == Nraw) {
		// time$day is present
			if (!CopyNumericArrToF32Arr(day32,day, Nraw)) {
				r_error("ERROR: metadata$time$day has an unsupported data format or type.\n");
				free(mn32);
				return 0;
			}
			
			for (int i = 0; i < Nraw; ++i) {
				yr32[i] = YMDtoF32time(yr32[i], mn32[i], day32[i]);
				if (yr32[i] < -1e9) {
					r_error("ERROR: The (%d)-ith date (metadata$time$year=%d,metadata$time$month=%d, and metadata$time$day=%d) is not valid.\n",i+1, (int) yr32[i], (int)mn32[i], (int)day32[i] );
					free(mn32);
					return 0;
				}
			}
			
		} else {
			// time$day is not present
			r_warning("WARNING: metadata$time$day is not specified, so only time$year and time$month are used!\n");

			for (int i = 0; i < Nraw; ++i) {
				yr32[i] = YMDtoF32time(yr32[i], mn32[i], 15);
				if (yr32[i] < -1e9) {
					r_error("ERROR: The (%d)-ith date (metadata$time$year=%d,and metadata$time$month=%d) is not valid.\n", i+1, (int)yr32[i], (int)mn32[i]);
					free(mn32);
					return 0;
				}
			}			
		}

		isTimeProcessed = 1;
		free(mn32);
		return 1L;

	}  

	if (!isTimeProcessed && yr && doy && IsNumeric(yr) && IsNumeric(doy) && GetNumberOfElements(yr) == Nraw && GetNumberOfElements(doy) == Nraw)
	{
		F32PTR yr32  = time;
		F32PTR doy32 = malloc(sizeof(F32) * Nraw); 

		if (!CopyNumericArrToF32Arr(yr32, yr, Nraw)) {
			r_error("ERROR: metadata$time$year has an unsupported data format or type.\n");
			free(doy32);
			return 0L;
		}
		if (!CopyNumericArrToF32Arr(doy32, doy, Nraw)) {
			r_error("ERROR: metadata$time$doy has an unsupported data format or type.\n");
			free(doy32);
			return 0L;
		}
 
		for (int i = 0; i < Nraw; ++i) {
			yr32[i] = YDOYtoF32time(yr32[i], doy32[i] );
			if (yr32[i] < -1e9) {
				r_error("ERROR: The (%d)-ith date (metadata$time$year=%d,and metadata$time$doy=%d) is not valid.\n", i + 1, (int)yr32[i], (int)doy32[i] );
				free(doy32);
				return 0;
			}
		}

		isTimeProcessed = 1;
		free(doy32);
		return 1L;
	}


	VOIDPTR datestr = GetField123Check(timeField, "dateStr",3);
	VOIDPTR strfmt  = GetField123Check(timeField, "strFmt",3); 
	if (!isTimeProcessed && datestr && strfmt && IsChar(datestr) && IsChar(strfmt) && GetNumberOfElements(datestr)==Nraw ) 	{		

		char STRFmt[255 + 1];
		GetCharArray(strfmt, STRFmt, 255);
		
		DateFmtPattern1 fmt1;	
		if (GetStrPattern_fmt1(STRFmt, &fmt1)) {
			
			for (int i = 0; i < Nraw; ++i) {	
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					r_error("ERROR: Unable to read the %d-ith date string from metadata$time$dateStr.\n", i + 1);
					return 0L;
				}
				time[i] = Str2F32time_fmt1(TMP, &fmt1);
				if (time[i] < -1e9) {
					r_error("ERROR: The %d-th string ($metadata$time$dateStr=\"%s\") is invalid, incompatiable with the specified "
							" metadata$time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					return 0L;
				}
			}	
		}
	 
		DateFmtPattern2 fmt2;	
		if (GetStrPattern_fmt2(STRFmt, &fmt2)) {
			
			for (int i = 0; i < Nraw; ++i) {	
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					r_error("ERROR: Unable to read the %d-ith date string from metadata$time$dateStr.\n", i + 1);
					return 0L;
				}
				time[i] = Str2F32time_fmt2(TMP, &fmt2);
				if (time[i] < -1e9) {
					r_error("ERROR: The %d-th string (metadata$time$dateStr=\"%s\") is invalid, incompatiable with the specified "
							" metadata$time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					return 0L;
				}
			}	
		}


		DateFmtPattern3 fmt3;	
		if (GetStrPattern_fmt3(STRFmt, &fmt3)) {
			
			for (int i = 0; i < Nraw; ++i) {	
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					r_error("ERROR: Unable to read the %d-ith date string from metadata$time$dateStr.\n", i + 1);
					return 0L;
				}
				time[i] = Str2F32time_fmt3(TMP, &fmt3);
				if (time[i] < -1e9) {
					r_error("ERROR: The %d-th string ($metadata$time$dateStr=\"%s\") is invalid, incompatiable with the specified "
							" metadata$time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					return 0L;
				}
			}	
		}
		return 1L;

	}
	 

	return 0L;
}
static int  ___PraseMetaData_IrregularTS_Dim123(BEAST2_METADATA_PTR meta, BEAST2_IO_PTR io) {
	// If the input time series data are irregular
	if (meta->nrhs < 3 || meta->isMetaStruct==0) {
		r_printf("ERROR: For irregular time series (i.e., metadata$isRegularOdered=FALSE), the argument 'metadata' must be supplied and it should be a LIST (for R) or struct (for Matlab) variable!");
		return 0;
	}
	if (meta->rawTimeVec == NULL) {
		r_printf("ERROR: For irregular time series (i.e., metadata$isRegularOdered=FALSE), metadata$time must be given to "
			     "specify the times at which data points are collected.");
		return 0;
	}

	/************************************************************/
	// REad the raw time and gather the info needed for aggragration
	/************************************************************/
	I32    Nraw = __GetRawTimeDimension(io);
	F32PTR time = malloc(sizeof(F32) * Nraw);	
	// Nraw must is the same as the length of rawTimeVec, to be checked in __ReadRawTime		
	if (__ReadRawTime(time, meta->rawTimeVec, Nraw)) {
		F32  dT = meta->deltaTime;
		io->N   = tsAggegrationPrepare(
					time, Nraw, dT, &io->T.sortedTimeIdx, &io->T.numPtsPerInterval,
					&io->T.startIdxOfFirsInterval, &meta->startTime
				); // meta->stattTime may be adjustd slightly to better match the real start of the time series
		free(time);
	} else {
		free(time);
		r_error("ERROR: Unable to read and intepret 'metadata$time'!\n");
		return 0;
	}

	// Npad:N_extended; the mutiples of 8 closest to N, defined for 32-byte alginment	
	//opt->Npad16 = (N+15)/16 *16;   //(I32)ceil((F32)N / 16.f) * 16;
	return 1L;
}
static int ParseInputData( BEAST2_IO_PTR _OUT_ io){

		// Check if the input 'data' is numeric and lengthy enough
		VOIDPTR DATA  =	 io->meta.rawInput;
		int     numel =  GetNumberOfElements(DATA);
		if (!(IsDouble(DATA) && numel > 2) && 	!(IsSingle(DATA) && numel > 2) &&
			!(IsInt32(DATA)  && numel > 2) &&	!(IsInt16(DATA) && numel > 2) &&
			!(IsInt64(DATA) && numel > 2)  && !(IsStruct(DATA) && numel >= 1 ) &&
			!(IsCell(DATA) ) /*Only true for Matlab*/      // Y=rawInput is a struct variable with multivariate time series
		   ) {   
			r_error("ERROR: The input data should be numeric and must be long enough.\n");
			return 0;
		}

		/*************************************************/
		/* If pY is not a string but a matrix or vector  */
		/* find the dimesions, and configure io-*/
		/*************************************************/
		I32       q;
		VOID_PTR  Y;
		if ( (IsStruct(DATA) && numel >= 1) || IsCell(DATA)) {
			// For MRBEAST only
			// rawInput is a struct variable with multivariate time series						
			q = numel;

			//Mem to be dellocated in DeallocatTimeSeriesIO
			io->pdata = malloc(sizeof(VOID_PTR) * q);
			io->dtype = malloc(sizeof(DATA_TYPE) * q);
			//TODO: For MRBEAST, here the data type an dimensions of all elments are assumed to be the same
			// A sanity check needs to be done to ensure that is the case.
			for (I32 i = 0; i < q; i++) {
				 Y           = GetFieldByIdx(DATA, i);
				io->pdata[i] = GetData(Y);
				io->dtype[i] = GetDataType(Y);
			}		
			//Aded for MRBEAST
			
		} else {
		    //rawIinput is a vector of numeric type, and not a struct variable with multivariate time series			
			q = 1;
			//Mem to be dellocated in DeallocatTimeSeriesIO
			io->pdata    = malloc(sizeof(VOID_PTR) * q);
			io->dtype    = malloc(sizeof(DATA_TYPE) * q);
			io->pdata[0] = GetData(DATA);
			io->dtype[0] = GetDataType(DATA);
			Y = DATA;
		}				
		for (I32 i = 0; i < q; i++) {
			if (io->dtype[i] == DATA_UNKNOWN) {
				r_error("ERROR: The input data has an uknown numeric type!\n");
				return 0;
			}
		}
		io->q = q;
		/*************************************************/
		/* Get the dimension of the input              */		
		/*************************************************/
		I32 ndims = GetNumOfDim(Y);
		if (ndims == 0) {
			// the input is a vector: this branch is possible only for R vectors (not  Matlab)
			I32 N = GetNumberOfElements(Y);

			io->numOfPixels = 1L;
			io->ndim		= 1L;
			io->dims[0]		= N;
			io->dims[1]		= 1L;
			
			if (io->meta.whichDimIsTime != -1) 
				r_warning("WARNING: metadata$whichDimIsTime=%d is ignored because 'whichDimIsTime' is used only if the input "
					      "is a 2D matrix or 3D array but your input is a 1D vector.\n", io->meta.whichDimIsTime);

			//TODO: always set it to 1 for 1D and 2D inputs
			// If hte input is a vector, the time dim is always 1
			io->meta.whichDimIsTime = 1L;
		} 
		else if (ndims == 1) {
			// the input is a vector: this branch is possible only for Python vectors (not R or matlab)
			I32 N = GetNumberOfElements(Y);

			io->numOfPixels = 1L;
			io->ndim     = 1L;
			io->dims[0]  = N;
			io->dims[1]  = 1L;

			if (io->meta.whichDimIsTime != -1)
				r_warning("WARNING: metadata$whichDimIsTime=%d is ignored because 'whichDimIsTime' is used only if the input "
					"is a 2D matrix or 3D array but your input is a 1D vector.\n", io->meta.whichDimIsTime);

			//TODO: always set it to 1 for 1D and 2D inputs
			// If hte input is a vector, the time dim is always 1
			io->meta.whichDimIsTime = 1L;
		}
		//ndims is impossible to be 1L
		else if (ndims == 2) { 
		// Matlab: a vector or matrix;
		// R:      a matrix or a vector with a dim attribute.

			int N = GetDim1(Y);
			int M = GetDim2(Y);

			if (min(N, M) == 1L)  //PY is a vector
			{
				N = max(N, M),
				io->numOfPixels = 1,
				io->ndim      = 1L,
				io->dims[0]   = N,
				io->dims[1]   = 1L,
				io->meta.whichDimIsTime = 1L;

				// If hte input is a vector, the time dim is always 1
				if (io->meta.whichDimIsTime != -1)
					r_warning("WARNING: metadata$whichDimIsTime=%d is ignored because 'whichDimIsTime' is used only if the input "
						"is a 2D matrix or 3D array but your input is a 1D vector.\n", io->meta.whichDimIsTime);
			}
			else                 //PY is a matrix
			{
				io->ndim	= 2L,
				io->dims[0] = N,
				io->dims[1] = M;


				I32 whichDimIsTime = io->meta.whichDimIsTime;
				if (whichDimIsTime == -1 || (whichDimIsTime != 1 && whichDimIsTime != 2 )) {
					r_error("ERROR: For a 2D matrix input of size [%d x %d] (i.e., multiple time series), metadata$whichDimIsTime must be given "
						    "to tell which dim of the matrix  refers to time. It must take a value out of 1 or 2 only.\n", N, M);
					return 0;
				}

				io->numOfPixels = (I64)io->dims[0] * io->dims[1]  / io->dims[io->meta.whichDimIsTime - 1L];
			}//PY is a matrix

			
		}
		else if (ndims == 3) {
			// If the input is a 3D array
			io->ndim = 3L;
			GetDimensions(Y, io->dims, 3L);

			I32 whichDimIsTime = io->meta.whichDimIsTime;
			if (whichDimIsTime == -1 ||	(whichDimIsTime!=1 && whichDimIsTime!=2 && whichDimIsTime!=3 )  ) {
				r_error("ERROR: For a 3D array input of size [%d x %d x %d] (i.e., stacked time series images), metadata$whichDimIsTime must be given "
					     "to tell which dim of the 3D array  refers to time. It must take a value out of 1, 2 or 3 only.\n", 
					     io->dims[0], io->dims[1], io->dims[2]);
				return 0;	}

			io->numOfPixels = (I64)io->dims[0] * io->dims[1] * io->dims[2] / io->dims[io->meta.whichDimIsTime - 1L];
		}
		else {
			r_printf("ERROR: The maximum dimension allowed is 3 when the data is a 3D stack of images over time,"
				    " but the input data has a dimension of %d .\n", ndims);
			return 0;
		}

	

		/*************************************************/
		// Check and get the values for period, startTime, deltaTime
		// or period_in_time
		/*************************************************/
		int res = io->meta.isRegularOrdered? ___PraseMetaData_RegularTS_Dim123(  &io->meta, io):
											 ___PraseMetaData_IrregularTS_Dim123(&io->meta, io); //set io->N
		if (res == 0)
			return 0;

		// Period is not needed for trend-only data
		// or if there is a season compnt and the period is already a valid value
		if ( !io->meta.hasSeasonCmpnt   || io->meta.period > 0)
			return 1;


		// There is a seasonal componet, and we need to check and get period	 
		/*************************************************/
		//  if period<=0 or IsNan(period)
		//  Determine the period perameter via auto-correlation
		/*************************************************/
		BEAST2_YINFO yInfo;
		F32PTR       MEMBUF;
		I32   N    = io->N;
		I32   Nraw = __GetRawTimeDimension(io);
		F32PTR tmp = malloc(sizeof(F32)*(N*q +N+q +q +q*q +Nraw )); //alocate mem for yInfo and MEMBUF

		yInfo.Y               = tmp;
		yInfo.rowsMissing     = tmp+N*q;
		yInfo.mean            = tmp + N * q+N;
		yInfo.sd                = tmp + N * q +N+q;
		yInfo.YtY_plus_alpha2Q  = tmp + N * q +N+q+q;
		MEMBUF                  = tmp + N * q +N+q+q+q*q; //needed only for irregular time series

		F32 period              = -1;
		F32 nan                 = getNaN();
		I32 goodPixelVisited    = 0;
		I32 MaxNumPixelVisisted = 200;

		// THe pixel index is 1-based and not zero-based
		for (int i = 1; i <= io->numOfPixels; ++i) {

			BEAST2_fetch_next_timeSeries(&yInfo, i, MEMBUF, io);
			// yInfo has been fillted above and now compaute mean, std, and YtY			
			int    N = io->N;
			int    q = io->q;
			//Normalize Y with NaN ommitted and then pre-comoute Y'*Y: YtY_plus_Q using gemm
			yInfo.nMissing = f32_normalize_multicols_zeroout_nans(yInfo.Y, yInfo.rowsMissing, N, N, q, yInfo.mean, yInfo.sd);
			U08 skipCurrentPixel = yInfo.nMissing > (N * io->meta.maxMissingRate) ? 1 : 0;
			if (skipCurrentPixel) {
				continue;
			}
		
			
			for (int j = 0; j < q; ++j) {
				F32PTR y = yInfo.Y + j * N;
				for (int k = 0; k < yInfo.nMissing; ++k) {y[yInfo.rowsMissing[k]] = nan;}

				F32 curPeriod = DeterminePeriod(y, N);  // return -1 if failing to estimate the period
				if (j == 0){ //if it is the 1st out of the q time series				
					period = curPeriod;
					if (period < 0) {
						break;
					}
				} else {
					if (period != curPeriod) {
						//if a later time series dones't give the same period as the first one
						period = -1;
						break;
					}
				}
			}

			if (period > 0) {
				break;
			}

			if (++goodPixelVisited > MaxNumPixelVisisted)
				break;
		}
		free(tmp);
		

		if (period > 0) {
			r_warning("WARNING: When metadata$season='%s' (i.e., the time series has a periodic component), \"metadata$period\" "
				"MUST be known in advance and specified by the user but it is missing. A BEST GUESS of it is %f (period=freq*deltaTime=%d*%f) and "
				"will be used. Please make sure this estimate makes sense; otherwise, the BEAST decomposition result will be incorrect.\n",
				io->meta.seasonForm == 'S' ? "harmonic" : "dummy", period* io->meta.deltaTime, (int)period, io->meta.deltaTime);
		} else {
			r_error("ERROR: When metadata$season='%s' (i.e., a periodic component present in the input time series), the \"metadata$period\" parameter "
				"must be known in advance and specified by the user but it is missing. BEAST tried to estimate it via an auotcorrelation method but failed to "
				"get a reliable estimate. Please specify the value for metadata$period EXPLICILTY. Or if your input has no periodic component at all, "
				" set  metadata$season='none', which will fit a trend-only model.\n", io->meta.seasonForm == 'S' ? "harmonic" : "dummy");
			return 0;
		}

		io->meta.period = period;
		return 1;

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

		r_warning("WARNING: The arg prior$precPriorType=(%d) is not a valid value; the default prior$precPriorType='%s' is assumed instread!", value, "uniform");
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

		r_warning("WARNING: The arg prior$precPriorType=(%s) is not recongizable; the default prior$precPriorType='%s' is assumed instread!", str, "uniform");
		return UniformPrec;
	}
	
	r_warning("WARNING: The arg prior$precPriorType has an supported format or value; the default prior$precPriorType='%s' is assumed instread!",  "uniform");
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
			r_warning("WARNING: The arg 'prior' is ignored because it is not a List/Struct variable.");
			memset(&m, 1L, sizeof(struct PRIOR_MISSING));
		}
		else {
			VOIDPTR tmp;
			o.seasonMinOrder    = (tmp = GetField123Check(S, "seasonMinOrder",10)) ?    GetScalar(tmp) : (m.seasonMinOrder = 1);
			o.seasonMaxOrder    = (tmp = GetField123Check(S, "seasonMaxOrder", 10)) ?   GetScalar(tmp) : (m.seasonMaxOrder = 1);
			o.trendMinOrder     = (tmp = GetField123Check(S, "trendMinOrder", 10)) ?    GetScalar(tmp) : (m.trendMinOrder = 1);
			o.trendMaxOrder     = (tmp = GetField123Check(S, "trendMaxOrder", 10)) ?    GetScalar(tmp) : (m.trendMaxOrder = 1);
			o.trendMinSepDist   = (tmp = GetField123Check(S, "trendMinSepDist", 10)) ?  GetScalar(tmp) : (m.trendMinSepDist = 1);
			o.seasonMinSepDist  = (tmp = GetField123Check(S, "seasonMinSepDist", 10)) ? GetScalar(tmp) : (m.seasonMinSepDist = 1);
			o.trendMinKnotNum   = (tmp = GetField123Check(S, "trendMinKnotNum", 10)) ?  GetScalar(tmp) : (m.trendMinKnotNum = 1);
			o.seasonMinKnotNum  = (tmp = GetField123Check(S, "seasonMinKnotNum", 10)) ? GetScalar(tmp) : (m.seasonMinKnotNum = 1);
			o.trendMaxKnotNum   = (tmp = GetField123Check(S, "trendMaxKnotNum", 10)) ?  GetScalar(tmp) : (m.trendMaxKnotNum = 1);
			o.seasonMaxKnotNum  = (tmp = GetField123Check(S, "seasonMaxKnotNum", 10)) ? GetScalar(tmp) : (m.seasonMaxKnotNum = 1);
			o.outlierMaxKnotNum = (tmp = GetField123Check(S, "outlierMaxKnotNum", 10))? GetScalar(tmp) : (m.outlierMaxKnotNum = 1);
			o.K_MAX             = (tmp = GetField123Check(S, "K_MAX",1)) ?			GetScalar(tmp) : (m.K_MAX = 1);

			o.sigFactor         = (tmp = GetFieldCheck(S,  "sigFactor")) ?			GetScalar(tmp) : (m.sigFactor        = 1);
			o.outlierSigFactor  = (tmp = GetFieldCheck(S,  "outlierSigFactor")) ?	GetScalar(tmp) : (m.outlierSigFactor = 1);

			o.sig2              = (tmp = GetField123Check(S, "sig2",2)) ?				GetScalar(tmp) : (m.sig2 = 1);
			o.precValue		    = (tmp = GetField123Check(S, "precValue",5)) ?		GetScalar(tmp) : (m.precValue = 1);
			o.alpha1			= (tmp = GetField123Check(S, "alpha1",0)) ?			GetScalar(tmp) : (m.alpha1 = 1);
			o.alpha2			= (tmp = GetField123Check(S, "alpha2",0)) ?			GetScalar(tmp) : (m.alpha2 = 1);
			o.delta1			= (tmp = GetField123Check(S, "delta1",0)) ?			GetScalar(tmp) : (m.delta1 = 1);
			o.delta2			= (tmp = GetField123Check(S, "delta2",0)) ?			GetScalar(tmp) : (m.delta2 = 1);

		
			o.seasonBasisFuncType	= (tmp = GetField123Check(S, "seasonBasisFuncType",10)) ?  GetScalar(tmp) : (m.seasonBasisFuncType = 1);
			o.trendBasisFuncType	= (tmp = GetField123Check(S, "trendBasisFuncType", 10)) ?   GetScalar(tmp) : (m.trendBasisFuncType = 1);
			o.outlierBasisFuncType	= (tmp = GetField123Check(S, "outlierBasisFuncType", 10)) ? GetScalar(tmp) : (m.outlierBasisFuncType = 1);
			o.modelPriorType		= (tmp = GetField123Check(S, "modelPriorType",  10)) ?		GetScalar(tmp) : (m.modelPriorType = 1);
			
			//o.precPriorType		    = (tmp = GetFieldCheck(S, "precPriorType")) ?	    GetScalar(tmp) : (m.precPriorType = 1);
			o.precPriorType = __GetPrecPriorType(S);

		
		}

	} // if (nrhs >= 4)

	o.numBasis       = 1L + io->meta.hasSeasonCmpnt + io->meta.hasOutlierCmpnt;
	I32  basisIdx    = 0;	
	if (io->meta.hasSeasonCmpnt) {
		I08  seasonFrom = io->meta.seasonForm;
		if      (seasonFrom == 'S')	o.basisType[basisIdx++] = SEASONID;
		else if (seasonFrom == 'D') o.basisType[basisIdx++] = DUMMYID;
		else if (seasonFrom == 'V') o.basisType[basisIdx++] = SVDID;
		else {
			r_error("ERROR: the seasonform character is unrecognized. It must be one of 'S', 'D', or 'V'. \n");
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

	if (m.seasonMinOrder)    o.seasonMinOrder   = 1L;				   o.seasonMinOrder   = min(o.seasonMinOrder, period / 2 - 1);    o.seasonMinOrder = max(o.seasonMinOrder, 1L);
	if (m.seasonMaxOrder)    o.seasonMaxOrder   = (period/2 - 1);      o.seasonMaxOrder   = min(o.seasonMaxOrder, (period / 2 - 1));  o.seasonMaxOrder = max(o.seasonMaxOrder, o.seasonMinOrder);
	if (m.trendMinOrder)     o.trendMinOrder    = 0L;				   o.trendMinOrder	  = max(o.trendMinOrder, 0L);
	if (m.trendMaxOrder)     o.trendMaxOrder    = 1L;				   o.trendMaxOrder	  = max(o.trendMaxOrder, o.trendMinOrder);
	if (m.seasonMinSepDist || o.seasonMinSepDist == 0)				   o.seasonMinSepDist = period / 2;       o.seasonMinSepDist = max(o.seasonMinSepDist, o.seasonMaxOrder);		 o.seasonMinSepDist = min(o.seasonMinSepDist, N / 2 - 1); // TODO:N/2-1 can be negtative, and then forced into a lager postive U16 integer
	
	if (m.trendMinSepDist || o.trendMinSepDist == 0)   o.trendMinSepDist = io->meta.hasSeasonCmpnt? period / 2: 3 ;
	o.trendMinSepDist = max(o.trendMinSepDist, (o.trendMaxOrder+1));
	o.trendMinSepDist = min(o.trendMinSepDist, N / 2 - 1);
	
	
	if (m.seasonMinKnotNum)  o.seasonMinKnotNum = 0;                   o.seasonMinKnotNum = max(min(o.seasonMaxKnotNum, o.seasonMinKnotNum), 0);
	if (m.trendMinKnotNum)   o.trendMinKnotNum  = 0;	               o.trendMinKnotNum  = max( min(o.trendMaxKnotNum,  o.trendMinKnotNum),0);
	
	if (m.seasonMaxKnotNum)  o.seasonMaxKnotNum  = min(floor(N / (o.seasonMinSepDist + 1) - 1.f), 5);  o.seasonMaxKnotNum = min(o.seasonMaxKnotNum, floor(N / (o.seasonMinSepDist + 1) - 1.f));
	if (m.trendMaxKnotNum)   o.trendMaxKnotNum   = min( floor(N/(o.trendMinSepDist  + 1) - 1.f),  10); o.trendMaxKnotNum  = min(o.trendMaxKnotNum,  floor(N / (o.trendMinSepDist + 1) - 1.f));	
	if (m.outlierMaxKnotNum) o.outlierMaxKnotNum = o.trendMaxKnotNum;                                  o.outlierMaxKnotNum = max(o.outlierMaxKnotNum, 1L); // at least has one; otherwise, the program crahses if hasOUtliercomponet=1
	
	if (m.K_MAX )            o.K_MAX            = 500;                  
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
			r_warning("WARNING: The arg 'mcmc' is ignored because it is not a LIST variable.");
			memset(&m, 1L, sizeof(struct MCMC_MISSING));
		} else {
			VOIDPTR tmp;
			o.maxMoveStepSize = (tmp = GetField123Check(S, "maxMoveStepSize",2))? GetScalar(tmp) : (m.maxMoveStepSize = 1);
			o.samples         = (tmp = GetField123Check(S, "samples", 2)) ?        GetScalar(tmp) : (m.samples = 1);
			o.thinningFactor  = (tmp = GetField123Check(S, "thinningFactor", 2)) ? GetScalar(tmp) : (m.thinningFactor = 1);
			o.burnin          = (tmp = GetField123Check(S, "burnin", 2)) ?         GetScalar(tmp) : (m.burnin = 1);
			o.chainNumber     = (tmp = GetField123Check(S, "chainNumber", 2)) ?    GetScalar(tmp) : (m.chainNumber = 1);
			o.seed			  = (tmp = GetField123Check(S, "seed", 2)) ?			GetScalar(tmp) : (m.seed = 1);
			o.ridgeFactor	  = (tmp = GetField123Check(S, "ridgeFactor", 2)) ?	GetScalar(tmp) : (m.ridgeFactor = 1);

			o.trendResamplingOrderProb  = (tmp = GetField123Check(S, "trendResamplingOrderProb",  2)) ? GetScalar(tmp) : (m.trendResamplingOrderProb = 1);
			o.seasonResamplingOrderProb = (tmp = GetField123Check(S, "seasonResamplingOrderProb", 2)) ? GetScalar(tmp) : (m.seasonResamplingOrderProb = 1);
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
			r_warning("WARNING: The arg 'extra' is ignored because it is not a LIST variable.");
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
			#define _4(x,y,z,w) _2(y,z);_2(y,z)
			#define _5(x,y,z,w) _2(y,z);_3(y,z)

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
		opt->mcmc.seasonResamplingOrderProb = 0;
		opt->prior.seasonMinOrder           = 0;
		opt->prior.seasonMaxOrder           = 0;		
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
			r_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0' will"
				     " fix the trend to a global curve.\n");					 
		}
		if (hasHarmonicCmpnt && hasOutlierCmpnt && isSeasonCmpntFixed) {
			r_warning("WARNING: The options 'seasonMaxOrder==seasonMaxOrder && seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
				    " fix the season component to a global curve.\n");		 
		}
		if (hasDummyCmpnt && hasOutlierCmpnt && isSeasonCmpntFixed) {
			// for the dummy cmpnt, maxOrder and minOrder are both set to zeros.
			r_warning("WARNING: The options 'seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
				    " fix the dummy season component to a global curve.\n");			
		}
	}
	if (hasTrendCmpnt && hasDummyCmpnt & isTrendCmpntFixed && isSeasonCmpntFixed) {		
		r_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0 && "
			    " seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
			    " fix the model structures of the trend and dummy season components.\n");
 
	}
	if (hasTrendCmpnt && hasHarmonicCmpnt & isTrendCmpntFixed && isSeasonCmpntFixed) {
		r_warning("WARNING: The options 'trendMaxOrder==trendMaxOrder && trendMaxKnotNum==0 && trendMinKnotNum==0 && "
			    "seasonMaxOrder==seasonMaxOrder && seasonMaxKnotNum==0 && seasonMinKnotNum==0' will"
			    " fix the model structures of the trend and harmonic season components.\n"); 
	}

	I32 KMAX=0;
	for (I32 i = 0; i < PRIOR->numBasis; i++) {		
		I08 type = PRIOR->basisType[i];
		if (type == SEASONID)			KMAX += (PRIOR->seasonMaxOrder*2) * (PRIOR->seasonMaxKnotNum + 1);
		if (type == DUMMYID)			KMAX += (opt->io.meta.period)    * (PRIOR->seasonMaxKnotNum + 1);
		if (type == SVDID)			    KMAX += (opt->io.meta.period) * (PRIOR->seasonMaxKnotNum + 1);
		if (type == TRENDID)			KMAX += (PRIOR->trendMaxOrder +1) * (PRIOR->trendMaxKnotNum + 1);
		if (type == OUTLIERID)			KMAX += PRIOR->outlierMaxKnotNum;
	}
	PRIOR->K_MAX = min(PRIOR->K_MAX, KMAX);

	
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
		r_warning("WARNING: prior$precPriorType is changed from 'componentwise' to 'uniform' because the model specified only has one component.\n");
	}
	return 1;
}

int BEAST2_GetArgs(VOIDPTR prhs[], int nrhs, A(OPTIONS_PTR) opt)
{
	int  failed = !GetArg_1st_MetaData(prhs, nrhs, &opt->io.meta)		 || 
				  !ParseInputData(&opt->io)								 || 
			      !GetArg_2nd_Prior__(prhs, nrhs, &opt->prior, &opt->io)   ||
			      !GetArg_3rd_MCMC___(prhs, nrhs, &opt->mcmc,  opt)     ||
			      !GetArg_4th_EXTRA__(prhs, nrhs, &opt->extra, opt->io.meta.whichDimIsTime,opt->io.ndim) ;
	int success = !failed;	
	if (success) 	success=PostCheckArgs(opt); 	
	if (success) 	BEAST2_print_options(opt);	
	return success;
}


void BEAST2_DeallocateTimeSeriesIO(BEAST2_IO_PTR  o) {
	// Free mems allocated in the "tsAggegrationPrepare" function
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