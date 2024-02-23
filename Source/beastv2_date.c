#include "abc_000_warning.h"

#include "abc_001_config.h"

#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>	               //fprintf fopen FILE #include<stdio.h>  // Need _GNU_SOURCE for manylinux; otherwise report /usr/include/stdio.h:316:6: error: unknown type name '_IO_cookie_io_functions_t'

#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_ide_util.h"  //printf
#include "abc_common.h"    //strcicmp
#include "abc_ts_func.h"
#include "abc_date.h"
#include "beastv2_func.h"    
#include "beastv2_io.h"


void* to_fyear(void * TIMEobj) {

	TimeVecInfo tvec = { 0 };
	TimeVec_init(&tvec);


	//  Allocated mem that needs to be freeed explicilty for tvec.fyear, but noo mem allocated if TIMEObj
	//  is NULL, which means that the ts is regular/ordered, as dertermined by start and dt only
	int Nrawtime = TimeVec_from_TimeObject(TIMEobj, &tvec);        // isDate may be updated inside based on the TimeObj

	int   nptr = 0;
	void* ans = NULL;
	if (Nrawtime > 0 ){
		double* data = NULL;
		ans = CreateNumVector(DATA_DOUBLE, Nrawtime, &data);
		PROTECT(ans);
		nptr++;
		memcpy(data, tvec.f64time, sizeof(double) * Nrawtime);
	}

	TimeVec_kill(&tvec); // Deallocate f64time only; the other allocate mem is still needed
	UNPROTECT(nptr);
	return ans;
}

void* from_fyear(void* FYEARobj) {

	if (!IsSingle(FYEARobj) && !IsDouble(FYEARobj)) {
		return NULL;
	}
   
	int N = GetNumberOfElements(FYEARobj);
	int* year, * mon, * day, * hr, * min, *wkday;
	double* sec;
	FIELD_ITEM  fieldList[] = {
	{"year",      DATA_INT32,		1,   {N,},      &year},
	{"mon",      DATA_INT32,		1,   {N,},      &mon}, //Needed to be changed to reflect the output dim
	{"day",       DATA_INT32,		1,   {N,},      &day},
	{"hr",       DATA_INT32,		1,   {N,},      &hr},
    {"min",       DATA_INT32,		1,   {N,},      &min},
	{"sec",       DATA_DOUBLE,		1,   {N,},      &sec},
	{"weekday",       DATA_INT32,		1,   {N,},    &wkday},
	{"",                                                },
	};
 
	int   nptr = 0;
	void * ans = PROTECT(CreateStructVar(fieldList, 99999));
	nptr++;

	YmdHms date;

	if (IsSingle(FYEARobj)){
		float* data = GetData(FYEARobj);
		for (int i = 0; i < N; i++) {
			FracYear_to_YMDHMS(data[i], &date);
			year[i] = date.y;
			day[i] = date.d;
			mon[i] = date.m;
			hr[i] = date.hr;
			min[i] = date.min;
			sec[i] = date.sec;
			wkday[i] = WeekDay(year[i],  mon[i], day[i]);
		}
		
	}
	else if (IsDouble(FYEARobj)) {
		double* data = GetData(FYEARobj);
		for (int i = 0; i < N; i++) {
			FracYear_to_YMDHMS(data[i], &date);
			year[i] = date.y;
			day[i] = date.d;
			mon[i] = date.m;
			hr[i] = date.hr;
			min[i] = date.min;
			sec[i] = date.sec;
			wkday[i] = WeekDay(year[i], mon[i], day[i]);
		}

	}


	UNPROTECT(nptr);
	return ans;
}

#include "abc_000_warning.h"