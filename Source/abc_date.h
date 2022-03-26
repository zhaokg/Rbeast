#include<stdio.h>
#include<inttypes.h>


//int     isLeapYear(int year);
//int     date2doy(int year, int mon, int day);
//int     numdays(int year);
int     isValidDate(int year, int mon, int day);
float   fractional_civil_from_days(int days);
//int64_t countLeapYears(int64_t year);
//int64_t datenum(int year, int mon, int day);

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

int    GetStrPattern_fmt1(char* strings, DateFmtPattern1* pattern);
float  Str2F32time_fmt1(char* strings, DateFmtPattern1* pattern);
int    GetStrPattern_fmt2(char* strings, DateFmtPattern2* pattern);
float  Str2F32time_fmt2(char* strings, DateFmtPattern2* pattern);
int    GetStrPattern_fmt3(char* strings, DateFmtPattern3* pattern);
float  Str2F32time_fmt3(char* strings, DateFmtPattern3* pattern);


float YDOYtoF32time(int year, int doy);
float YMDtoF32time(int year, int mon, int day);
