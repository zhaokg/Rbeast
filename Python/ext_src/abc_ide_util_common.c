#include <math.h>
#include <string.h>
#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"
#include "abc_vec.h"
#include "abc_sort.h"
 

// C++ solution to define a keyword as a function name
// #pragma warning(suppress: 4483)
// extern void __identifier("?ioFlush@@YA_NXZ")(void);

I08 IDE_USER_INTERRUPT;

void* CreateNumVector(DATA_TYPE dtype, int length, VOIDPTR* data_ptr) {	
	return CreateNumVar( dtype, &length, 1, data_ptr);
}

void* CreateNumMatrix(DATA_TYPE dtype, int Nrow, int Ncol, VOIDPTR* data_ptr) {
	int dims[] = {Nrow, Ncol};
	return CreateNumVar(dtype, dims, 2, data_ptr);
}

void* CreateF32NumVector(int length, VOIDPTR* data_ptr) {
	return CreateNumVar(DATA_FLOAT, &length, 1, data_ptr);
}

void* CreateF32NumMatrix(int Nrow, int Ncol, VOIDPTR* data_ptr) {
	int dims[] = { Nrow, Ncol };
	return CreateNumVar(DATA_FLOAT, dims, 2, data_ptr);
}

void* CreateF64NumVector(int length, VOIDPTR* data_ptr) {
	return CreateNumVar(DATA_DOUBLE, &length, 1, data_ptr);
}

void* CreateF64NumMatrix(int Nrow, int Ncol, VOIDPTR* data_ptr) {
	int dims[] = { Nrow, Ncol };
	return CreateNumVar(DATA_DOUBLE, dims, 2, data_ptr);
}

void* CreateI32NumVector(int length, VOIDPTR* data_ptr) {
	return CreateNumVar(DATA_INT32, &length, 1, data_ptr);
}

void* CreateI32NumMatrix(int Nrow, int Ncol, VOIDPTR* data_ptr) {
	int dims[] = { Nrow, Ncol };
	return CreateNumVar(DATA_INT32, dims, 2, data_ptr);
}

void printProgress1(F32 pct, I32 width, char * buf, I32 firstTimeRun) {
	
	//https:// stackoverflow.com/questions/2685435/cooler-ascii-spinners

 	static char spinnerChar[] =  "|/-\\";
	static I32  cnt = 1;
	cnt++;
	cnt = (cnt == 4) ? 0 : cnt;
	width = max(width, 35);  // cause errors when width is 0 or negative
	memset(buf, '*', width); // space = 20

	I32  len = 0;
	buf[len++] = spinnerChar[cnt];

	char prefix[] = "Progress:";
	I32 strLen    = sizeof(prefix)-1L; //the last byte is a Zero
	memcpy(buf+len, prefix, strLen);
	len += strLen;

	snprintf(buf + len, 15L, "%5.1f%% done", pct * 100);
	len += 5+1+5;
	buf[len++] = '[';

	I32 finishedLen = round((width - len - 1)*pct);
	memset(buf + len, '=', finishedLen);
	len += finishedLen;
	buf[len++] = '>';
	//memset(buf + len, ' ', width-len);
	buf[width - 1] = ']';
	buf[width] = 0;

#if R_INTERFACE==1
	r_printf("\r%s", buf);
	//R doesnto allow io operationrs from external libraries
	//fflush(stdout);
#elif P_INTERFACE==1
	r_printf("\r%s", buf);
#elif M_INTERFACE==1
	if (firstTimeRun == 1) 	{
		r_printf("\r\n");
		r_printf("%s", buf);
	} else {
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf);
	}

	StdouFlush();
	//mexEvalString("drawnow");
	//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
	//mexEvalString("pause(0.0000001);drawnow");
#endif	
}

void printProgress2(F32 pct, F64 time, I32 width, char * buf, I32 firstTimeRun)
{//https:// stackoverflow.com/questions/2685435/cooler-ascii-spinners
	
	static char spinnerChar[] = "|/-\\";
	static int  count         = 1;
	count = (++count) == 4 ? 0 : count;
	width = max(width, 40); //cause errors when width is 0 or negative
	memset(buf, '*', width); // space = 20

	I32  len = 0;
	buf[len++] = (pct<1.0) ? spinnerChar[count]: ' ';

	snprintf(buf + len, 10L, "%5.1f%%", pct * 100);
	len += 5 + 1;

	char prefix[] = "done";
	I32 strLen = sizeof(prefix)-1L; // the last byte is a Zero
	memcpy(buf + len, prefix, strLen);
	len += strLen;

	F64 SecsPerDay = 3600 * 24;
	I32 days = time / SecsPerDay;
	F64 tmp = time - days *SecsPerDay;
	I32 hrs = tmp / 3600;
	tmp = (tmp - hrs * 3600);
	I32 mins = tmp / 60;
	tmp = tmp - mins * 60;
	I32 secs = tmp;
	days = days >= 99 ? 99 : days;

	//https://stackoverflow.com/questions/71013186/warning-name0-directive-output-may-be-truncated-writing-10-bytes-into-a-regi
	if (time > SecsPerDay)
		snprintf(buf + len, 100L, "<Remaining%02dday%02dhrs%02dmin>", days, hrs, mins);
	else
		snprintf(buf + len, 100L,"<Remaining%02dhrs%02dmin%02dsec>", hrs, mins, secs);
	len += 26;

	buf[len++] = '[';

	I32 finishedLen = round((width - len - 1)*pct);
	memset(buf + len, '=', finishedLen);
	len += finishedLen;
	buf[len++] = '>';
	//memset(buf + len, ' ', width-len);
	buf[width - 1] = ']';
	buf[width] = 0;

#if R_INTERFACE==1
	r_printf("\r%s", buf);
	//R doesnto allow io operationrs from external libraries
	//fflush(stdout);
#elif P_INTERFACE==1
	r_printf("\r%s", buf);
#elif M_INTERFACE==1
	if (firstTimeRun == 1)	{
		r_printf("\r\n");
		r_printf("%s", buf);
	} else {
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf); 
	}

	//mexEvalString("drawnow");
	//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
	StdouFlush();
#endif	
}

void RemoveField(FIELD_ITEM *fieldList, int nfields, char * fieldName)
{
	for (I64 i = 0; i < nfields && fieldList[i].name[0]!=0; i++) {
		if (strcmp(fieldList[i].name, fieldName) == 0)	{
			// The two strings are equal
			if (fieldList[i].ptr) {
				fieldList[i].ptr[0] = NULL; // the storage pointer is set to NULL
			}
			fieldList[i].ptr = NULL;        // the pointer to the storage pointer is set to NULL
			break;
		}
	}
}


void RemoveSingltonDims(FIELD_ITEM* flist, I32 nlist) {
	// 1 x1    -> 1
	// 2x1x1x2 -> 2x2
	// 1x10    -> 10

	// 10 : in R, this will be created as a vector, and in Matlab 
	// be created as a column vector.

	 #define  MAX_NUM_DIM   10 
	for (int i = 0; i < nlist && flist[i].name[0]!=0; i++) {

		if ( flist[i].ndim == 1) continue;		

		int goodN   = 0;
		int goodDims[MAX_NUM_DIM];
		// FInd all the non-one dims
		for (int j = 0; j < flist[i].ndim; j++) {
			if (flist[i].dims[j] != 1L) {
				goodDims[goodN++] = flist[i].dims[j];
			}
		}

		// Re-assign the non-one dims
		if (goodN == 0) {
			// it is 1x 1 x1 ..,
			flist[i].ndim    = 1;
			flist[i].dims[0] = 1;
		} else {
			flist[i].ndim = goodN;
			for (int j = 0; j < goodN; j++) {
				flist[i].dims[j] = goodDims[j];
			}
		}
	}

}

int CopyNumericObjToF32Arr(F32PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if (IsSingle(infield))     		memcpy(outmem, data, sizeof(F32) * N);
	else if (IsDouble(infield))		for (I32 i = 0; i < N; i++) outmem[i] = (F32)*((double*)data + i);
	else if (IsInt32(infield))		for (I32 i = 0; i < N; i++) outmem[i] = (F32) *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = (F32) *((I64*)data + i);
	else if (IsChar(infield))		return 0;
	else {	
		return 0;
	}
	return 1L;
}

int CopyNumericObjToF64Arr(F64PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if (IsDouble(infield))     		memcpy(outmem, data, sizeof(F64) * N);
	else if (IsSingle(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((F32*)data + i);
	else if (IsInt32(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((I64*)data + i);
	else if (IsChar(infield))		return 0;
	else {	
		return 0;
	}
	return 1L;
}
int CopyNumericObjToI32Arr(I32PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if      (IsInt32(infield))    	memcpy(outmem, data, sizeof(I32) * N);
	else if (IsDouble(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((double*)data + i);
	else if (IsSingle(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((I64*)data + i);
	else if (IsChar(infield))	  return 0;
	else {
		return 0;
	}
	return 1L;
}

int HaveEqualDimesions(const void* p1, const void* p2) {
	int dim1 = GetNumOfDim(p1);
	int dim2 = GetNumOfDim(p2);
	if (dim1 != dim2) return 0;

	I32 dims1[5], dims2[5];
	GetDimensions(p1, dims1, dim1);
	GetDimensions(p2, dims2, dim2);

	I32 equal = 1;
	for (int i = 0; i < dim1; ++i) {
		equal = equal & (dims1[i] == dims2[i]);
	}
	return equal;
}

int GetDataType(VOID_PTR Y) {

	if      (IsInt32(Y)) 						return DATA_INT32;
	else if (IsInt16(Y)) 						return DATA_INT16;
	else if (IsInt64(Y)) 						return DATA_INT64;
	else if (IsDouble(Y))  /* isReal(pY)*/		return DATA_DOUBLE;
	else if (IsSingle(Y))  /* isReal(pY)*/ 		return DATA_FLOAT;
	else                                        return DATA_UNKNOWN;
}

#if R_INTERFACE ==1 || M_INTERFACE ==1
F64 GetNumericElement(const void* Y, I32 idx) {

	if (!IsNumeric(Y)) {
		return getNaN();
	}

	I32 n = GetNumberOfElements(Y);
	if (n == 1) {
		// Y is a scalar
		if (idx==0)		return GetScalar(Y);
		else       		return getNaN();
	} else {
		// Y is a vector of more than one elements
		if (idx < n) {
			VOID_PTR y = GetData(Y);
			if (IsInt32(Y)) 						    return *((I32PTR)y + idx);
			else if (IsInt16(Y)) 						return *((I16PTR)y + idx);
			else if (IsDouble(Y))  /* isReal(pY)*/		return *((F64PTR)y + idx);
			else if (IsSingle(Y))  /* isReal(pY)*/ 		return *((F32PTR)y + idx);
			else                                        return getNaN();
		}	else {
			// idx is out of bounds
			return getNaN();
		}
	}

}

void* CvtToPyArray_NewRef(VOIDPTR Y) { return Y; }
#endif

void* GetField123Check(const void* structVar, char* fname, int nPartial) {
	// Check if the retured field is an empety or R_NilValue
	if (structVar == NULL) return NULL;

	VOID_PTR p = GetField123(structVar, fname, nPartial);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}

}

void* GetFieldCheck(const void* structVar, char* fname) {
	if (structVar == NULL) return NULL;
	// Check if the retured field is an empety or R_NilValue
	VOID_PTR p = GetField(structVar, fname);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}
}

int  GetNumElemTimeObject(  VOID_PTR timeObj ) {		
	// time should have been pre-allocated, with a length of Nraw
	if (timeObj == NULL) return -1L;
 
	/////////////////////////////////////////////////////
	// timeObj is not a struct variable
	/////////////////////////////////////////////////////	
	if (IsNumeric(timeObj) || IsChar(timeObj)) {
		//R's date object is a REALSXP.
		return GetNumberOfElements(timeObj);
	} // 	if (IsNumeric(timeObj)==0) 

	// if it is not numeric or Struct.
	if ( IsStruct(timeObj) == 0 ) {
		return -1;
	} //if ( IsStruct(timeObj) == 0 )

	  /////////////////////////////////////////////////////
	// timeObj  is a struct variable
	/////////////////////////////////////////////////////	
	VOIDPTR yr  = GetField123Check(timeObj, "year",1);
	VOIDPTR mn  = GetField123Check(timeObj, "month",1);
	VOIDPTR day = GetField123Check(timeObj, "day",3);
	VOIDPTR doy = GetField123Check(timeObj, "doy",3);

	int isTimeProcessed = 0;
	if (!isTimeProcessed && yr && mn && IsNumeric(yr) && IsNumeric(mn) ) {
	
		int Nyr = GetNumberOfElements(yr);
		int Nmn = GetNumberOfElements(mn);
		int Ntime = Nyr;
		if (Nyr != Nmn) {		 
			return -1L;
		}	
   
		if (day && IsNumeric(day) && GetNumberOfElements(day) == Ntime) {	} 
		else {
			return -1L;
		}

		isTimeProcessed = 1; 
		return Ntime;
	}  

	if (!isTimeProcessed && yr && doy && IsNumeric(yr) && IsNumeric(doy) ) 	{
 
		int Nyr  = GetNumberOfElements(yr);
		int Ndoy = GetNumberOfElements(doy);
		int Ntime = Nyr;
		if (Nyr != Ndoy) {		 
			return -1;
		}	 
		isTimeProcessed = 1;	 
		return Ntime;
	}

	VOIDPTR datestr = GetField123Check(timeObj, "dateStr",3); 
	if (!isTimeProcessed && datestr && IsChar(datestr) ) 	{
		return GetNumberOfElements(datestr); 
	}
	return -1;
}

#define isEqualOneDay(dT) ( fabs(dT-1./365.)<1e-3 || fabs(dT - 1./366.) < 1e-3) 

 F64 Parse_TimeIntervalObject(VOIDPTR * obj, TimeScalarInfo* tint) {

	 int isProcessed = 0;
	
	if ( obj == NULL || IsEmpty(obj) ) {
		tint->fyear=tint->value = getNaN();
		tint->unit = 'U';
		return 	tint->value;
	}
	else if (obj && IsNumeric(obj) && GetNumberOfElements(obj) == 1) {
		tint->value = GetScalar(obj);
		tint->fyear = tint->value;
		tint->unit = 'U';
		return 	tint->value;	
	}
	else if (obj && IsChar(obj)) {
		char s[40 + 1];
		GetCharArray(obj, s, 40);

		if (strcicmp(s, "none")==0) {
			// This is for metadata$period only
			tint->fyear=tint->value = 0;
			tint->unit  = 'U';
			return tint->value;
		}

		char unit;
		F64 fyear=extract_timeinterval_from_str(s, &tint->value, &unit);
		tint->unit  = unit;
		tint->fyear = fyear;

		if (IsNaN(fyear)) {
			r_error("ERROR: Unable to intepret '%s' as a time interval for metadata$deltaTime or metadata$period !\n", s);
			tint->unit = 'B';
		}
		return tint->value;
	}
	else {

		tint->fyear=tint->value = getNaN();
		tint->unit = 'B';

		r_error("ERROR: metadata$deltaTime pr metadata$period is of unrecognizable type!\n");
		return 	tint->value;
	}

	
}

 F64  Parse_DateOriginObject(VOIDPTR * obj) {

	 // Convert the origin into Julidan Day Number

	int Yorg=-9999, Morg  , Dorg  ;
	F64 origin_datenum;

 
	#if R_INTERFACE ==1
		Yorg = 1970, Morg = 1; Dorg = 1;
	#elif P_INTERFACE==1
		Yorg = 1, Morg = 1; Dorg = 0;
	#elif M_INTERFACE==1
		Yorg = 0, Morg = 1; Dorg = 0;
	#endif

	 int isProcessed = 0;
	
	if ( obj == NULL || IsEmpty(obj) ) {
		// Use the default origin based on the current Platform
		if (Yorg != -9999) { // Unkonw platforms/interfaces
			isProcessed = 1;
		}		
	}
	else if (obj && IsNumeric(obj) && GetNumberOfElements(obj) == 3) {
		Yorg = GetNumericElement(obj, 0);
		Morg = GetNumericElement(obj, 1);
		Dorg = GetNumericElement(obj, 2);

		if (Morg < 0 || Morg>12 || Dorg < 0 || Dorg>31) {
			isProcessed = 0;
		}	else {
			isProcessed = 1;
		}		
	}
	else if (obj && IsChar(obj)) {
		char originstr[40 + 1];
		GetCharArray(obj, originstr, 40);
			
		if (strcicmp(originstr, "matlab") == 0) {
			Yorg = 0, Morg = 1; Dorg = 0;
			isProcessed = 1;
		}
		else if (strcicmp(originstr, "R") == 0 ||
			strcicmp(originstr, "unix") == 0 ||
			strcicmp(originstr, "posix") == 0 ||
			strcicmp(originstr, "pandas") == 0) {
			Yorg = 1970, Morg = 1; Dorg = 1;
			isProcessed = 1;
		}
		else if (strcicmp(originstr, "python") == 0) {
			Yorg = 1, Morg = 1; Dorg = 0;
			isProcessed = 1;
		}
		else {
			int startidx[1] = { 0 };
			F64 ftime;
			int Nout =FracYear_from_Strings(&ftime, originstr, startidx, 1);
			if (Nout>0) {
				FracYear_to_YMD(ftime, &Yorg, &Morg, &Dorg);
				isProcessed = 1;
 
			} 
		}
	}

	if (isProcessed) {
		return  JulianDayNum_from_civil_ag1(Yorg, Morg, Dorg);
	}

	return getNaN(); 
}

 F64  Parse_SingelDateObject(VOIDPTR* obj, TimeScalarInfo* tint) {

	 int isProcessed = 0;

	 if (obj == NULL || IsEmpty(obj)) {
		 tint->value = getNaN();
		 tint->fyear = tint->value;
		 tint->unit = 'U';
		 return 	tint->value;
	 }
	 else if (obj && IsClass(obj, "Date"))  {
		// This branch is true only for R.  Matlab's IsClass always return FALSE
		// in R, tmp is an integer array

		 F64 value = GetScalar(obj);
		 tint->fyear =tint->value= FracYear_from_DateNum(value+R_DATE_ORIGIN);		 
		 tint->unit  = 'Y';

		 return tint->value;
	 } 
	 else if (obj && IsNumeric(obj)  ) {
		 char* msg = "\n"
			 "(1) a numeric scalar (e.g., 2002.33)\n"
			 "(2) a vector of two values (Year, Month) (e.g., c(2002,4) in R, [2002,4] in Matlab/Python) \n"
			 "(3) a vector of three values (Year, Month, Day) (e.g., c(2002,4,15) in R,  [2002,4,15] in Matlab/Python)\n"
			 "(4) a date string (e.g., \"2002-4-15\", \"2002/04/15\", or \"2002/4\")\n"
			 "(5) a datenum as a list/struct variable (e.g., list(datenum=12523, origin='R') or struct('datenum', 731321, 'origin','matlab')\n";
		 I32   n     = GetNumberOfElements(obj);
		 F32   start = getNaN();
		 if (n == 1) {
			 start = GetScalar(obj);
			 tint->unit = 'U';
		 }		
		 else if (n == 2) {
			 F32 Y = GetNumericElement(obj, 0), M = GetNumericElement(obj, 1);
			 start = Y + M / 12 - 1 / 24.0;
			 tint->unit = 'Y';
			 q_warning("INFO: metadata$startTime=[%g,%g] is interpreted as %04d/%02d/15 (Year/Month/Day) and converted to a decimal fractional year "
				       "of %g. If not making sense, supply a correct start time using one of the following for startTime: %s\n", Y, M, (int)Y, (int)M, start, msg);
		 }
		 else if (n == 3) {
			 F32 Y = GetNumericElement(obj, 0), M = GetNumericElement(obj, 1), D = GetNumericElement(obj, 2);
			 start = (F32)FracYear_from_YMD((int)Y, (int)M, (int)D);
			 tint->unit = 'Y';
			 q_warning("INFO: Your metadata$startTime=[%g,%g,%g] is interpreted as %04d/%02d/%02d (Year/Month/Day) and converted to a decimal year "
				       "of %g. If not making sense, supply a correct start time using one of the following for startTime:  %s\n", Y, M, D, (int)Y, (int)M, (int)D, start, msg);
		 }
		 else {
			 F32 Y = GetNumericElement(obj, 0), M = GetNumericElement(obj, 1), D = GetNumericElement(obj, 2);
			 q_warning("ERROR: Your metadata$startTime =[%g, %g, %g, ...] has more than three elements and can't be intepreted as a valid date o time. "
				       " A default value will be used. %s\n",    Y, M, D,	 msg);
			 tint->unit = 'B';
		 }
		 tint->fyear = tint->value = start;
 		 return tint->value;
	 }
	 else  {
		 TimeVecInfo tv = {0,};		 
		 tv.isDate  = UnknownStatus;

		 int N=TimeVec_from_TimeObject(obj, &tv);
		 
		 if (N>0) {
			 tint->fyear = tint->value = tv.f64time[0];
			 tint->unit  = tv.isDate == 1? 'Y': 'U';			 			 
		 }	 else {
			 tint->fyear = tint->value = getNaN();
			 tint->unit = 'B';		  //Bad format
		 }
 
		 TimeVec_kill(&tv);
		 return tint->value;
	 }

	 tint->fyear = tint->value = getNaN();	 
	 tint->unit  = 'U';
	 return 	tint->value;
 }

 static  void* quick_timevec_allocator_f64time(TimeVecInfo* tv, int Nnew) {
	 if (Nnew <= tv->Ncapacity_fyear) {
		 tv->N    = 0;
		 tv->Nbad = 0;
	 }  else {
		 if (tv->f64time) free(tv->f64time);
		 tv->f64time = malloc(sizeof(F64) * Nnew);
		 tv->Ncapacity_fyear = Nnew;
		 tv->N      = 0;
		 tv->Nbad   = 0;		 
	 }
	 return tv->f64time;
 }

 static  void* quick_timevec_allocator_sortidx(TimeVecInfo* tv, int Nnew) {
	 if (Nnew <= tv->Ncapacity_sortidx) {
		 return tv->sorted_time_indices;
	 }
	 else {
		 if (tv->sorted_time_indices) free(tv->sorted_time_indices);
		 tv->sorted_time_indices = malloc(sizeof(F32) * Nnew);
		 tv->Ncapacity_sortidx = Nnew;
		 //tv->N = 0;      DO not touch N and Nbad
		 //tv->Nbad = 0;
		 return tv->sorted_time_indices;
	 }
 }

 void TimeVec_init(TimeVecInfo* tv) {
	 
	 tv->isDaily = UnknownStatus;
	 tv->isDate = UnknownStatus;	 
	 tv->isMonthly = UnknownStatus;
	 tv->isRegular = UnknownStatus;
	 tv->IsOrdered  = UnknownStatus;
	 tv->isStartDeltaOnly = UnknownStatus;
	 tv->isConvertedFromStartDeltaOnly = 0;// Need to be expliclty updated to 1
	 tv->isDateNum                     = 0; // Need to be expliclty updated to 1


	 tv->data_dt     = getNaN();
	 tv->data_start  = getNaN();
	 tv->data_period = getNaN();
	 tv->N    = 0;
	 tv->Nbad = 0;
		 

 }
  void TimeVec_kill(TimeVecInfo* tv) {
	 if (tv->f64time)               free(tv->f64time);
	 if (tv->sorted_time_indices) free(tv->sorted_time_indices);
	 tv->f64time = NULL;
	 tv->sorted_time_indices = NULL;
	 tv->Ncapacity_fyear           = 0;
	 tv->Ncapacity_sortidx         = 0;
 }

  void TimeVec_kill_fyearArray(TimeVecInfo* tv) {
	  if (tv->f64time)               free(tv->f64time);
 	   tv->f64time = NULL;
	  tv->Ncapacity_fyear = 0;
  }


int TimeVec_from_TimeObject(VOID_PTR timeObj, TimeVecInfo * tv ) {
	 
	//  List of UpdateD Fields: fYear (new mem allocated), N, and isDate (only if the time is date. Keep it untouched if uknown).

	tv->isStartDeltaOnly = 0;
	tv->N                = 0;
	tv->Nbad             = 0;

	// timeObj is NULL, then, assume the time series is REGULAR. But do nothing here
	if (timeObj == NULL || IsEmpty(timeObj)) {
		return tv->N;
	}
	int isTimeProcessed = 0;
 	/////////////////////////////////////////////////////
	// timeObj is NOT NULL and it is not  a struct variable
	///////////////////////////////////////////////////////	
	VOIDPTR datestr = NULL; // possibly assigned when timeObj itselt is a char array.
	if (  IsStruct(timeObj) == 0   ) {

		int Ntime = GetNumberOfElements(timeObj);
		if (Ntime < 0) {
			return tv->N;
		}
		
		if ( IsNumeric(timeObj)  ) {

			F64PTR f64time = quick_timevec_allocator_f64time(tv, Ntime);	 
			tv->N          = Ntime;

			if (IsClass(timeObj, "Date")) {
				// This branch is true only for R and Matlab's IsClass always return FALSE
				// in R, tmpObJ is an REAL array
				int status = CopyNumericObjToF64Arr(f64time, timeObj, Ntime);
				F64PTR days = f64time;
				for (int i = 0; i < Ntime; ++i) {
					f64time[i] = FracYear_from_DateNum(f64time[i] + R_DATE_ORIGIN);
				}
				tv->isDate = 1;
			} else if (IsClass(timeObj, "POSIXct") || IsClass(timeObj, "POSIXt")) {
				// This branch is true only for R and Matlab's IsClass always return FALSE
				// in R, tmpObJ is an REAL array
				int status = CopyNumericObjToF64Arr(f64time, timeObj, Ntime);
				F64PTR secs = f64time;
				F64    secsInOneDay = 3600 * 24;
				for (int i = 0; i < Ntime; ++i) {
					f64time[i] = FracYear_from_DateNum(secs[i] / secsInOneDay + R_DATE_ORIGIN);
				}
				tv->isDate = 1;
			}

			else {
				int status = CopyNumericObjToF64Arr(f64time, timeObj, Ntime);
				if ( status == 0 ) {										
					tv->N      = 0;					
					q_warning("ERROR: time has an unsupported data format or type.\n");					 
				}
				// DO not touch tv->Date
			}

			return tv->N;
		}
		
		if (IsChar(timeObj)) {
			// It is a string and direclty jump to the smart section 
			datestr = timeObj;  
			goto __ENTRY_DATESTR_SMART_LOC;
		}

		// Otherwise, if it is not numeric or char. There must be an error	
		tv->N     =  0;
		q_warning("ERROR: time is not numeric or strings. If times are strings, use time$dateStr and time$strFmt to specify data observation times.\n");
		return tv->N;
		
	} //if ( IsStruct(timeObj) == 0 )


	/////////////////////////////////////////////////////
	// timeObj  is a struct variable
	/////////////////////////////////////////////////////	

	VOIDPTR yr  = GetField123Check(timeObj, "year",1);
	VOIDPTR mn  = GetField123Check(timeObj, "month",1);
	VOIDPTR day = GetField123Check(timeObj, "day",3);
	VOIDPTR doy = GetField123Check(timeObj, "doy",3);



	/////////////////////////////////////////////////////
	// From time$year, time$month, and time$day
	/////////////////////////////////////////////////////	 
	if (!isTimeProcessed && yr && mn && IsNumeric(yr) && IsNumeric(mn) ) {
	
		int Nyr = GetNumberOfElements(yr);
		int Nmn = GetNumberOfElements(mn);
		if (Nyr != Nmn) {
			q_warning("WARNING: time$year and time$month should have the same length.\n"); 
			goto __ENTRY_YEAR_DOY_LOC;
		}	
	
		int    Ntime = Nyr;
		F64PTR yr64  = quick_timevec_allocator_f64time(tv, Ntime + (Ntime+1)/2); // Ntime for ur32 and N/2 for mn32
		F32PTR mn32  = (F32PTR)(yr64 + Ntime);         
		F32PTR day32 = quick_timevec_allocator_sortidx(tv, Ntime);

		if (!CopyNumericObjToF64Arr(yr64, yr, Ntime)) {
			q_warning("WARNING: time$year has an unsupported data format or type.\n");
			goto __ENTRY_YEAR_DOY_LOC;
		}
		if (!CopyNumericObjToF32Arr(mn32, mn, Ntime)) {
			q_warning("WARNING: time$month has an unsupported data format or type.\n");
			goto __ENTRY_YEAR_DOY_LOC;
		}
 
		if (day && IsNumeric(day) && GetNumberOfElements(day) == Ntime) {
		/* time$day is present */

			if (!CopyNumericObjToF32Arr(day32,day, Ntime)) {
				q_warning("WARNING: time$day has an unsupported data format or type.\n");
				goto __ENTRY_YEAR_DOY_LOC;
			}
			
			for (int i = 0; i < Ntime; ++i) {
				yr64[i] = FracYear_from_YMD( round(yr64[i]), round(mn32[i]), round(day32[i]) );
				if (yr64[i] < -1e9) {
					q_warning("WARNING: The (%d)-ith date (time$year=%d,time$month=%d, and time$day=%d) is not valid.\n",i+1, (int) yr64[i], (int)mn32[i], (int)day32[i] );
					goto __ENTRY_YEAR_DOY_LOC;
				}
			}		
 
			
		} else {
		/*  time$day is not present */

			q_warning("WARNING: time$day is not specified, so only time$year and time$month are used!\n");
			for (int i = 0; i < Ntime; ++i) {
				yr64[i] = yr64[i] + mn32[i]/12.0 -1./24.0;
				if (yr64[i] < -1e9) {
					q_warning("WARNING: The (%d)-ith date (metadata$time$year=%d,and metadata$time$month=%d) is not valid.\n", i+1, (int)yr64[i], (int)mn32[i]);
					goto __ENTRY_YEAR_DOY_LOC;
				}
			}	

			tv->isMonthly = 1;
		}

		isTimeProcessed = 1;
		tv->isDate = 1;
		tv->N      = Ntime;

		return tv->N;
	}  


	__ENTRY_YEAR_DOY_LOC:
	/////////////////////////////////////////////////////
	// From time$year, time$doy, 
	/////////////////////////////////////////////////////	
	if ( !isTimeProcessed && yr && doy && IsNumeric(yr) && IsNumeric(doy) ) {
 
		int Nyr  = GetNumberOfElements(yr);
		int Ndoy = GetNumberOfElements(doy);	
		if (Nyr != Ndoy) {
			q_warning("WARNING: time$year and time$doy should have the same length.\n");
			goto __ENTRY_DATESTR_LOC;
		}
		int  Ntime = Nyr;

		F64PTR  yr64  = quick_timevec_allocator_f64time(tv, Ntime);
		F32PTR  doy32 = quick_timevec_allocator_sortidx(tv, Ntime);

		if (!CopyNumericObjToF64Arr(yr64, yr, Ntime)) {
			q_warning("WARNING: time$year has an unsupported data format or type.\n");
			goto __ENTRY_DATESTR_LOC;
		}
		if (!CopyNumericObjToF32Arr(doy32, doy, Ntime)) {
			q_warning("WARNING: time$doy has an unsupported data format or type.\n");
			goto __ENTRY_DATESTR_LOC;
		}
 
		for (int i = 0; i < Ntime; ++i) {
			yr64[i] = FracYear_from_intYDOY( round(yr64[i]), round( doy32[i] ) );
			if (yr64[i] < -1e9) {
				q_warning("WARNING: The (%d)-ith date ( time$year=%d,and  time$doy=%d) is not valid.\n", i + 1, (int)yr64[i], (int)doy32[i] );
				goto __ENTRY_DATESTR_LOC;
			}
		}

		isTimeProcessed = 1;
 		tv->isDate = 1;
		tv->N      = Ntime;
		return  tv->N;
	}


	//F32PTR   datestr; // datestr has already been defined above
	VOIDPTR  strfmt;

__ENTRY_DATESTR_LOC:

	/////////////////////////////////////////////////////
	// From datestr, and strfmt
	/////////////////////////////////////////////////////	

	 datestr = GetField123Check(timeObj, "dateStr", 3);
	 strfmt  = GetField123Check(timeObj, "strFmt", 3);

	if (!isTimeProcessed && datestr && strfmt && IsChar(datestr) && IsChar(strfmt)  ) 	{		

		int    Ntime   = GetNumberOfElements(datestr);
		F64PTR f64time = quick_timevec_allocator_f64time(tv, Ntime);
  
		char STRFmt[400 + 1];
		GetCharArray(strfmt, STRFmt, 400);
		
		DateFmtPattern1 fmt1;	
		if (GetStrPattern_fmt1(STRFmt, &fmt1)) {			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[400 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 400)) {			
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				f64time[i] = Str2F32time_fmt1(TMP, &fmt1);
				if (f64time[i] < -1e9) {
					q_warning("WARNING:  The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified "
							" time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}
	 
		DateFmtPattern2 fmt2;	
		if (GetStrPattern_fmt2(STRFmt, &fmt2)) {
			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				f64time[i] = Str2F32time_fmt2(TMP, &fmt2);
				if (f64time[i] < -1e9) {
					q_warning("WARNING: The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}


		DateFmtPattern3 fmt3;	
		if (GetStrPattern_fmt3(STRFmt, &fmt3)) {
			
			for (int i = 0; i < Ntime; ++i) {
				char TMP[255 + 1];
				if (!GetCharVecElem(datestr, i, TMP, 255)) {
					q_warning("WARNING: Unable to read the %d-ith date string from time$dateStr.\n", i + 1);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
				f64time[i] = Str2F32time_fmt3(TMP, &fmt3);
				if (f64time[i] < -1e9) {
					q_warning("WARNING: The %d-th string (time$dateStr=\"%s\") is invalid, incompatiable with the specified time$strFmat=\"%s\".\n", i + 1, TMP, STRFmt);
					goto __ENTRY_DATESTR_SMART_LOC;
				}
			}	
		}

		isTimeProcessed = 1;
 		tv->isDate = 1;
		tv->N      = Ntime;
		return  tv->N;

	}
	 

__ENTRY_DATESTR_SMART_LOC:
	/////////////////////////////////////////////////////
	// From datestr alone and strfmat is missing
	/////////////////////////////////////////////////////	
	if (!isTimeProcessed && datestr &&  IsChar(datestr)  ) {
			DynMemBuf        strarr  = { .max_len = 100 };
			DynAlignedBuf    start   = { .align = 4,.elem_size = 4,.p.raw = NULL, };
			DynAlignedBuf    nchars  = { .align = 4,.elem_size = 4,.p.raw = NULL, };
			CharObj2CharArr(datestr, &strarr, &start, &nchars);

			int    Ntime   = start.cur_len;
			F64PTR f64time = quick_timevec_allocator_f64time(tv, Ntime);
			int    Nout    = FracYear_from_Strings(f64time, strarr.raw, start.p.i32, Ntime);
	
			dynbuf_kill(&strarr);
			adynbuf_kill(&start);
			adynbuf_kill(&nchars);

			if (Nout>0) {
				tv->N        = Nout ;
				tv->isDate   = 1;
			}	else {
				tv->N = 0;
				goto __ENTRY_DATENUM_ORIGIN_LOC;
				// Do not touch isDate
				// COntintue to be processesed to see if there are other possibilites
			}			

			return tv->N;
		 
	}

	/////////////////////////////////////////////////////
	// From datenum and origin
	/////////////////////////////////////////////////////	
	VOIDPTR datenum;
	VOIDPTR origin;

__ENTRY_DATENUM_ORIGIN_LOC:
	datenum = IsStruct(timeObj) ? GetField123Check(timeObj, "datenum", 1) : NULL;
	origin  = IsStruct(timeObj) ? GetField123Check(timeObj, "origin", 1) :NULL;
	if (!isTimeProcessed && datenum && IsNumeric(datenum)  ) {

		int     Ntime   = GetNumberOfElements(datenum);
		F64PTR  f64time = quick_timevec_allocator_f64time(tv,Ntime);

		CopyNumericObjToF64Arr(f64time, datenum, Ntime);

		F64 originJDN = Parse_DateOriginObject(origin);

		if (IsNaN(originJDN)) {

		}	else {
			for (int i = 0; i < Ntime; i++) {	
				f64time[i] = FracYear_from_DateNum(f64time[i]+ originJDN);
			}
			tv->N     =  Ntime;
			tv->isDate = 1;
			return tv->N;
		}	 //	if (IsNaN(originJDN)) 	
	}

	tv->N     = 0;
	return tv->N;
}

void TimeVec_from_StartDeltaTime(TimeVecInfo* tv, F32 start, F32 dt, int N, int isDate) {
	/*
	if (IsNaN(tv->user_dt) && IsNaN(tv->user_start) && (IsNaN(tv->user_period) || tv->user_dt <= 0) && tv->isDate != UnknownStatus) {
		//this cominbation is impossible and there must be something wrong! isDate shouldn't be 1 if none of the variables were checked before
		assert(0);
	}
	*/

	tv->isStartDeltaOnly = 1;

	tv->N          = N;  // Must update N later
	tv->Nbad       = 0;
	tv->isRegular = 1;
	tv->IsOrdered = 1;
	tv->data_dt    = dt; // can be NaNs
	tv->data_start = start;  // can be NaNs

	if (isDate == 1 && isEqualOneDay(dt)) {
	// If this is the case, we will explictly converted it into an allocated array of time
    // as fractioonal years
		tv->isDaily =        1;
		tv->isRegular       = UnknownStatus;
		tv->isStartDeltaOnly             = 0;
		tv->isConvertedFromStartDeltaOnly = 1;

		tv->data_dt =  1. / 365;  // coearce it to 1/365

		F64PTR f64time          = quick_timevec_allocator_f64time(tv, N);  // tv.N is set to zero inside
		tv->sorted_time_indices = quick_timevec_allocator_sortidx(tv, N);
		double  startDateNum    = FracYear_to_DateNum(start);	
		for (int i = 0; i < N; i++) {
			f64time[i] = FracYear_from_DateNum(startDateNum + i);
			tv->sorted_time_indices[i] = i;
		}
		tv->N          =  N;
		tv->isDateNum  =  0;
		tv->data_start = f64time[0];
   }
	return; //No mem is allocated
}

static int IsRegularTS(F64PTR t, int N, F32* dt) {

	/************************************************/
	// Now oldTime is the sorted times, CHeck if the ts is regular
	/************************************************/
	int sRegular = 1;
	F64PTR  SortedTimes  = t;
	F32     dT_estimate = SortedTimes[1] - SortedTimes[0];
	F32     dT_mean = (SortedTimes[N - 1] - SortedTimes[0]) / (N - 1);
	F32     epilson = dT_mean * 1e-3;
	for (int i = 2; i < N; i++) {
		F32 dT = SortedTimes[i] - SortedTimes[i - 1];
		if (fabsf(dT - dT_estimate) > epilson) {
			sRegular = 0;
			break;
		}
	}

	*dt = sRegular ? dT_estimate : dT_mean;
	return sRegular;

	/*
	// Give it one more try if it is date
	if (t.isRegular == 0 && IsDate && t.isDateNum == 0) {

		t.isRegular = 1;
		// CHeck if it is daily data
		for (int i = 1; i < Nold; i++) {
			F32 dt = SortedTimes[i] - SortedTimes[i - 1];
			if (fabsf(dt - 1. / 365) > 1e-3 || fabsf(dt - 1. / 366) > 1e-3) {
				t.isRegular = 0;
				break;
			}
		}
	}

	// Give it one more try if it is date
	if (t.isRegular == 0 && IsDate && t.isDateNum == 0) {
		t.isRegular = 1;
		// CHeck if it is montly data
		for (int i = 1; i < Nold; i++) {
			F32 dt = SortedTimes[i] - SortedTimes[i - 1];
			if (fabsf(dt - 30 / 365) > 3. / 365 || fabsf(dt - 30 / 366) > 3. / 366) {
				t.isRegular = 0;
				break;
			}
		}
	}
	*/
}
 void TimeVec_SortCheckRegularOrder(TimeVecInfo* tv) {

#define t  (*tv)
	 if ( tv->N <= 0) {
		 return;
	 }
	  /**********************************************************************************************/
	  // If t.isStartDeltaOnly = 1 Reular and ordered time series
	  /**********************************************************************************************/
	  if ( tv->isStartDeltaOnly==1  ) { 
		  if (IsNaN(tv->data_dt) || IsNaN(tv->data_start))  {
			  r_printf("ERROR: startTime and deltatTime are NANs.when sorting the regular time series\n");
			  return;
		  }
		  tv->IsOrdered = 1;
		  tv->isRegular = 1;
		  return;
	  }
	  
	  /**********************************************************************************************/
	  // Thanks to Gary Iacobucci to noice this bug: time may also have NaNs, so as a remdy, we fill the 
	  // the NaNs with INF so that in the sorted time, they  will be in the last. We also adjust Nold to 
	  // a new value
	  /**********************************************************************************************/
	  F32 InfValue = (F32)(1.e36) * (F32)(1.e36);

	  F64PTR  oldTime    = t.f64time;
	  int     Nold       = t.N;	 
	  int     Nbadvalues = 0;
	  for (int i = 0; i < Nold; i++) {
		  Nbadvalues = Nbadvalues + IsNaN(oldTime[i]);
		  oldTime[i] = IsNaN(oldTime[i]) ? InfValue : oldTime[i];
	  }
	  t.Nbad = Nbadvalues;

	  /************************************************/
	  t.sorted_time_indices = quick_timevec_allocator_sortidx(tv, Nold);
	  i32_seq(t.sorted_time_indices, 0, 1, Nold); // fill with 0 to N-1.
	  f64a_introSort_index(oldTime, 0, Nold - 1, t.sorted_time_indices);
	  Nold = Nold - Nbadvalues;
	  // Now the sorted values should have all the Infs at the CODE_EOF, if any	
	  /************************************************/

	/************************************************/
	// Determine if hte ts is ordered
	/************************************************/
	 t.IsOrdered = 1L;
	if (tv->isConvertedFromStartDeltaOnly == 0) {		  
		  for (int i = 0; i < Nold; i++) {
			  if (t.sorted_time_indices[i] != i) { t.IsOrdered = 0;	break; }
		  }
	 }
	  

	 /************************************************/
	 // Determine if the ts is regular
	 /************************************************/
	 F32 dT;
	 t.isRegular = IsRegularTS(tv->f64time, Nold, &dT);
	 t.data_start = tv->f64time[0];
	 t.data_dt    = dT;

#undef t

}

 void TimeVec_UpdateUserStartDt(TimeVecInfo* tv) {

	  tv->out.asDailyTS = 0;

	  TimeVecInfo  t;
	  memcpy(&t, tv, sizeof(TimeVecInfo));

	  F32 dT     = t.out.dT;
	  F32 start  = t.out.start;
	  I32 IsDate = t.isDate;
	  F32 period = t.data_period;
	  int Ngood = t.N - t.Nbad;

	  // Sanity check: IsData should have been decieded at this point
	  if (IsDate !=0 && IsDate !=1) {
		  assert(0);  r_printf("ERROR: there must be something wrong with regards to isDate=%d.\n", IsDate);
		  return;
	  }
	  if (IsNaN(t.data_dt)|| IsNaN(t.data_start)) {
		  assert(0); 
		  r_printf("ERROR: there must be something wrong with regards to deltaTime and start.\n");
		  return;
	  }

 
	  /************************************************/
      // If dT is missing (i.e, dT==NAN), determne a best guess of it
	  /************************************************/	
	  if (IsNaN(dT)) {

		  dT = t.data_dt;  // Assign a default value

		  if (t.isRegular == 0) {		 
			  F32 freq = 1. / t.data_dt;
			  if (t.isDate == 1 && t.isDateNum==0 && dT > 0.5 / 366.) {
				  // if it is NOT a sub-daily time series
				  //  'float' which may cause truncation of value [-Wabsolute-value]				
				  F32 dt1 = fabs(freq - 365) / 365, dt2 = fabs(freq - 24) / 24.;
				  F32 dt3 = fabs(freq - 12) / 12.,  dt4 = fabs(freq - 365. / 7) / 52.; //52
				  F32 mindt = min(min(dt1, dt2), min(dt3, dt4));
				  if (dt1 == mindt) dT = 1. / 365;
				  if (dt2 == mindt) dT = 1. / 24;
				  if (dt3 == mindt) dT = 1. / 12;
				  if (dt4 == mindt) dT = 1. / 52;
			  }	  else if (period > 0) {
				  freq   = ceil(period /t.data_dt);
				  dT     = period / freq;
			  }
		  }	 

		  t.out.dT = dT ;
	  }


	  /************************************************/
	  // If start is missing (i.e, dT==NAN), determne a best guess of it
	  /************************************************/

	  F32 Tstart, Tend;	  
	  if ( t.isStartDeltaOnly == 1) {
		  Tstart = t.data_start;
		  Tend = t.data_start + (Ngood - 1) * t.data_dt;
	  }	  else {	 
		  Tstart = t.f64time[0];
		  Tend   = t.f64time[Ngood - 1];
	  }

	  // First of all, if start goes beyond the time series span, pull it back
	  if ( !IsNaN(start) && start >= Tend) {
		  t.out.start=start = t.data_start;
	  }

	  if (IsNaN(start)) {
		  if (t.isStartDeltaOnly == 1) {
			  // Regualr data determined from start and dt only
			  start = t.data_start;
		  }   else {
			  // data that has a time array fyear
			  I32 i0 = round(Tstart / dT), i1 = round(Tend / dT);
			  //F32 i0 = floor(T0/period);
			  //T0     = i0*period + floor( (T0-i0)/dT )*dT;		
			  I32 Nnew          = ((i1 - i0) + 1);
			  F32 start_estiamte = i0 * dT;   // the midpoint of the first internval
			  if ((t.isRegular == 1 && fabs(dT - t.data_dt) < dT * 1e-3) ) {
				  // ifit is regular, and the data dT is the same as the requestt dt
				  start = t.data_start;
			  }	  else {
				  start = start_estiamte;
			  }
		  }
		  t.out.start = start;
	  }
	  // At this point, start cannot be a NAN
	 
	  // Further ajust it if hte start is too far from the start
	   if ( fabs(start - t.data_start) >dT * 1e-2 && start < Tend && fabs(start - Tstart) > dT * 200) {
			  F32 T0 = start;
			  int i0 = round(T0 / dT);
			  int i1 = round(Tend / dT);
			  int Nnew = ((i1 - i0) + 1);
			  F32 start_estiamte = i0 * dT;
			  start = start_estiamte;
			  t.out.start = start;
	   }	   
	   /************************************************/
	   // Now sart and dt shoudl not be NANs
	   /************************************************/

	  /************************************************/
	  // if it is a PERIODIC daily time series with the period less than 300 days:
	  // we will coearce the time into days. No need to handle the leap years
	  /************************************************/

	  // Convert Fractional Yeat to JDN if it is trend-only or the priod is less than 300
		   
	  if ( IsDate == 1 && t.isDateNum==0 && isEqualOneDay(dT) &&
		     (   (period > 0 && period*365 < 300) || // a periodic ts with a less-than-one year period
			      period <= 0                        // a daily trend-only time series
		     )
		 ) 
	  {		   
		  F32  potentialPeriod = period* 365;
		  // force it to an integer if close enough
		  if (fabs(potentialPeriod - round(potentialPeriod)) < 1e-3) {
			  potentialPeriod = round(potentialPeriod);
		  }
 
		  t.out.asDailyTS   = 1;
		  t.out.start   = FracYear_to_DateNum(t.out.start);;
		  t.out.dT      = 1; //1 day		    
		  t.data_dt     = t.data_dt * 365;  // date's dt may be not 1 day
		  t.data_start  = FracYear_to_DateNum(t.data_start);
		  t.data_period = potentialPeriod;
		  t.isDateNum   = 1;  // TTODO: this is the only place wher isDateNum is set to 1
 
		  /***************************************/
		  /* Convert individual times to DateNum and Re-check if it is a regular ts */
		  /***************************************/
		  if (t.isStartDeltaOnly) {
			  // at this point, date_dt is impossible to be one-day. If yes, it has been converted a time array in TImeVec_from_StartDt
			  //We neeed to expliictly convert to a time array
			  t.f64time             = quick_timevec_allocator_f64time(&t, Ngood);
			  t.sorted_time_indices = quick_timevec_allocator_sortidx(&t, Ngood);
			  for (int i = 0; i < Ngood; i++) {
				  t.f64time[i] = t.data_start + t.data_dt * i;
				  t.sorted_time_indices[i] = i;
			  }
			  t.isStartDeltaOnly = 0;
			  t.isConvertedFromStartDeltaOnly = 1;

			  t.isRegular = 1;
		  }	  else {
			  for (int i = 0; i < Ngood; i++) {
				  t.f64time[i] = FracYear_to_DateNum( t.f64time[i] );
			  }
			  F32 dt;
			  t.isRegular = IsRegularTS(t.f64time, Ngood, &dt);
			  t.data_dt   = dt;
		  }

	
	  }

	  memcpy(tv, &t, sizeof(TimeVecInfo));
	  return;
 
}

 

void CharObj2CharArr(VOID_PTR o, DynMemBufPtr str , DynAlignedBufPtr charstart, DynAlignedBufPtr nchars){
	int n = GetNumberOfElements(o);
	dynbuf_init(str, n * 200);
	adynbuf_init(charstart, n);
	adynbuf_init(nchars,    n);
	for (int i = 0; i < n; i++) {
		dynbuf_requestmore(str, 200);
		int  nlen = GetCharVecElem(o, i, str->raw+str->cur_len, 200);
		charstart->p.i32[i] = str->cur_len;
		nchars   ->p.i32[i] = nlen;
		charstart->cur_len++;
		nchars   ->cur_len++;
		str->cur_len += (nlen+1); // 1L is the extra NULL.
	}	 
}

void obj_to_str(VOID_PTR o, DynMemBufPtr s, int leftMargin) {
	int nfld      = GetNumberOfFields(o);

	int maxKeyLen = 0;
	for (int i = 0; i < nfld; ++i) {
		char key[50];
		GetFieldNameByIdx(o, i, key, 50);
		int  keyLen = (int)strlen(key);
		if (keyLen > maxKeyLen) {
			maxKeyLen = keyLen;
		}
	} 

	for (int i = 0; i < nfld; ++i) {
		char key[50];
		GetFieldNameByIdx(o, i, key, 50);
		int  keyLen = (int)strlen(key);
		char tmp[200];
		snprintf(tmp, 199, "%*s%-*.*s : ", leftMargin, "", maxKeyLen, maxKeyLen, key);
		dynbuf_insert_bytes(s, tmp, (int)strlen(tmp)+1);

		VOIDPTR e  = GetFieldByIdx(o, i);
		if (e == NULL) {
			snprintf(tmp, 199, "[]\n");
			dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
		}
		else if (IsNumeric(e)) {
			int size  = GetNumberOfElements(e);
			int dtype = GetDataType(e);
			if (size == 1) {
				F64 value = GetScalar(e);				
				tmp[0] = 0;
				if (dtype == DATA_INT16 || dtype == DATA_INT32 || dtype == DATA_INT64) {
					snprintf(tmp,199, "%d\n",(int) value);
				} else if (dtype == DATA_FLOAT || dtype == DATA_DOUBLE) {
					snprintf(tmp,199, "%g\n", value);
				}	else {
					snprintf(tmp,199, "%g\n", value);
			
				}
				dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
				continue;
			}

			char* etype;
			if      (dtype == DATA_INT16)  etype = "int16";
			else if (dtype == DATA_INT32)  etype = "int32";
			else if (dtype == DATA_INT64)  etype = "int64";
			else if (dtype == DATA_FLOAT)  etype = "float32";
			else if (dtype == DATA_DOUBLE)  etype = "float64";
			else etype = "others";

			int ndims = GetNumOfDim(e);
			int dims[10];			
			GetDimensions(e, dims, ndims);

			snprintf(tmp, 199, "[%d", dims[0]);
			dynbuf_insert_bytes(s, tmp, strlen(tmp) + 1);
			for (int j = 1; j < ndims; j++) {
				snprintf(tmp, 199, "x%d", dims[j]);
				dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
			}
			snprintf(tmp, 199," %s] \n", etype);
			dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
		}
		else if (IsChar(e)) {
			char tmpstr[30];
			GetCharArray(e, tmpstr, 30);
			snprintf(tmp, 199,"'%s'\n", tmpstr);
			dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
		}
		else if (IsStruct(e)) {
			int len = GetNumberOfFields(e);
			snprintf(tmp, 199, "[ 1 object with %d fields] \n", len);
			dynbuf_insert_bytes(s, tmp, (int)strlen(tmp) + 1);
			obj_to_str(e, s, leftMargin + maxKeyLen + 3);
		}
	}
}

int IDEPrintObject(VOID_PTR o) {
	if (!IsStruct(o)) {
		r_printf("Not an object, structure, or list.\n");
	}

	int nfld = GetNumberOfFields(o);
	r_printf("Object of %d field(s): \n\n", nfld);

	DynMemBuf s = { 0, };
	dynbuf_init(&s, 1000);
	int leftMargin = 1;
	obj_to_str(o, &s, leftMargin);
	s.raw[s.cur_len] = 0;
	r_printf("%s", s.raw);
	dynbuf_kill(&s);
	return 0;
}

#include "abc_000_warning.h"