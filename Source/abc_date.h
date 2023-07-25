#pragma once 

#include<stdio.h>
#include<inttypes.h>


 
int     isValidDate(int year, int mon, int day);
//int64_t countLeapYears(int64_t year);
//int64_t datenum(int year, int mon, int day);

typedef struct  __Date__ {
	int   y;
	short m, d, hr, min;
	float sec;
} DateType;

typedef struct {
	int yearIdx, monIdx, dayIdx;
} DateFmtPattern1;

typedef struct {
	int yearIdx, doyIdx;
} DateFmtPattern2;

typedef struct {
	char order[3];
	char sep1[20];
	char sep2[20];
} DateFmtPattern3;

int     GetStrPattern_fmt1(char* strings, DateFmtPattern1* pattern);
double  Str2F32time_fmt1(char* strings, DateFmtPattern1* pattern);
int     GetStrPattern_fmt2(char* strings, DateFmtPattern2* pattern);
double  Str2F32time_fmt2(char* strings, DateFmtPattern2* pattern);
int     GetStrPattern_fmt3(char* strings, DateFmtPattern3* pattern);
double  Str2F32time_fmt3(char* strings, DateFmtPattern3* pattern);

double extract_timeinterval_from_str(char* s, float* value, char* unit);

int  JulianDayNum_from_civil_ag1(int y, int m, int d);
int JulianDayNum_to_Civil_ag1(int JDN, int* yr, int* mn, int* day);

double FracYear_from_YDOY(int year, int doy);
double FracYear_from_YMD(int year, int mon, int day);
int    FracYear_to_YDOY(double fyear, int* yr);
void   FracYear_to_YMD(double fyear, int* yr, int* mon, int* day);
 
double FracYear_to_DateNum(double fyear);
double FracYear_from_DateNum(double datenum);

void   FracYear_to_DateType(double fyear, DateType* date);
double FracYear_from_DateType(DateType* date);
 

int    FracYear_from_Strings(double * out, char* s, int* strstart, int n);
  

// Matlab datenum: the number of days since midnight on Jan 0st, 0 AD[not Jan 1 but Jan 0].
// 
//R's Origin:     2440588 : 1970, 1,1
//Mat's origin:   1721059 :  0,1,0
//Python's origin: 1721425:  1,1,0
#define NULL_DATE_ORIGIN  0
#define R_DATE_ORIGIN 2440588
#define MATLAB_DATE_ORIGIN 1721059
#define PYTHON_DATE_ORIGIN 1721425