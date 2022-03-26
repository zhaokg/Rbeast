#include<stdio.h>
#include<inttypes.h>
// include stdio.h library
#include <stdlib.h> // atof atoi
#include <string.h> //strchr strrchr strstr str

#include "abc_000_warning.h"

#include "abc_date.h"
#include "abc_common.h"   //ToUpper
#include "abc_ide_util.h" //ToUpper

 

//https://overiq.com/c-examples/c-program-to-calculate-the-day-of-year-from-the-date/
static int IsLeapYear(int year) { return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);}
static int GetNumDays(int year) { return IsLeapYear(year) ? 366 : 365; }
//https://stackoverflow.com/questions/19377396/c-get-day-of-year-from-date
static const int DAYS[2][13] = {
	{ 0, 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
	{ 0, 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 }
};

int Date2Doy(int year, int mon, int day) { 	return DAYS[IsLeapYear(year)][mon] + day; }
int Doy2Date(int doy,  int y, int* d, int* m){

		static int month[13] = { 0, 31, 28, 31, 30, 31, 30,	31, 31, 30, 31, 30, 31 };

		if (IsLeapYear(y))	month[2] = 29;

		int i;
		for (i = 1; i <= 12; i++) 		{
			if (doy <= month[i])
				break;
			doy = doy - month[i];
		}

		*d = doy;
		*m = i;

		return 0;
}


int IsValidDate(int year, int mon, int day )
{
	int is_valid = 1, is_leap = 0;

	if (year >= 1800 && year <= 9999)	{

		//  check whether year is a leap year
		is_leap = (year%4 == 0 && year%100 != 0) || (year % 400 == 0);
		
		// check whether mon is between 1 and 12
		if (mon >= 1 && mon <= 12) 		{
			// check for days in feb
			if (mon == 2) 	{
				if (is_leap && day == 29){
					is_valid = 1;
				} 
				else if (day > 28)
				{
					is_valid = 0;
				}
			}
			else if (mon == 4 || mon == 6 || mon == 9 || mon == 11)
			{	// check for days in April, June, September and November
				if (day > 30) {
					is_valid = 0;
				}
			}
			else if (day > 31)
			{  // check for days in rest of the months:Jan, Mar, May, July, Aug, Oct, Dec
				is_valid = 0;
			}
		}
		else
		{
			is_valid = 0;
		}

	}
	else
	{
		is_valid = 0;
	}

	return is_valid;

}

//codereview.stackexchange.com/questions/152017/simple-days-between-dates-calculator
int64_t CountLeapYears(int64_t year)
{
	// 438 - 17 + 4
	static int64_t fakeLeaps = (1753 / 4) - (1753 / 100) + (1753 / 400);
	// Leaps before 1753
	// We start at 0 in 1753
	int64_t leaps     = year / 4;
	int64_t badLeaps  = year / 100;
	int64_t extraLeaps = year / 400;

	return leaps - badLeaps + extraLeaps - fakeLeaps;
}

float YDOYtoF32time(int year, int doy) {
	return (float)year + ((float)doy - 0.5) / GetNumDays(year);
}

float YMDtoF32time(int year, int mon, int day) { 
	return YDOYtoF32time(year, Date2Doy(year, mon, day)); 
}


int64_t datenum(int year, int mon, int day) {
	int64_t numYears   = year - 1753;
	int64_t leapYears  = CountLeapYears(year);
	int64_t daysInYear = DAYS[IsLeapYear(year)][mon];
	return  numYears * 365LL+ leapYears + daysInYear + (day - 1);
}
//https://overiq.com/c-examples/c-program-to-calculate-the-day-of-year-from-the-date/

//https://stackoverflow.com/questions/14218894/number-of-days-between-two-dates-c
//https://stackoverflow.com/questions/33712685/c-days-between-given-date-and-today%C2%B4s-date
// Civil calendar is the Gregorian calendar: the one we normally use
// Returns number of days since civil 1970-01-01.  Negative values indicate
//    days prior to 1970-01-01.
// Preconditions:  y-m-d represents a date in the civil (Gregorian) calendar
//                 m is in [1, 12]
//                 d is in [1, last_day_of_month(y, m)]
//                 y is "approximately" in
//                   [numeric_limits<Int>::min()/366, numeric_limits<Int>::max()/366]
//                 Exact range of validity is:
//                 [civil_from_days(numeric_limits<Int>::min()),
//                  civil_from_days(numeric_limits<Int>::max()-719468)]
int days_from_civil(int y, unsigned m, unsigned d) {	
	y -= m <= 2;
	const int      era = (y >= 0 ? y : y - 399) / 400;
	const unsigned yoe = (y - era * 400);      // [0, 399]
	const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;  // [0, 365]
	const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;         // [0, 146096]
	return era * 146097 +doe - 719468;
}
//http://howardhinnant.github.io/date_algorithms.html
int civil_from_days(int days, int * yr, int*mn, int* day) 
{ 
	days += 719468;
	const int era = (days >= 0 ? days : days - 146096) / 146097;
	const unsigned doe = (days - era * 146097);          // [0, 146096]
	const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
	const int      y  = (yoe) + era * 400;
	const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                // [0, 365]
	const unsigned mp  = (5 * doy + 2) / 153;                                   // [0, 11]
	const unsigned d  = doy - (153 * mp + 2) / 5 + 1;                             // [1, 31]
	const unsigned m  = mp < 10 ? mp + 3 : mp - 9;                            // [1, 12]
	*yr = y + (m <= 2);
	*mn = m;
	*day = d;
 

	return 0;
}

// //In R, dates are represented as the number of days since 1970-01-01
float fractional_civil_from_days(int days)
{
	days += 719468;
	const int era = (days >= 0 ? days : days - 146096) / 146097;
	const unsigned doe = (days - era * 146097);          // [0, 146096]
	const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;  // [0, 399]
	const int      y = (yoe)+era * 400;
	const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);                // [0, 365]
	const unsigned mp = (5 * doy + 2) / 153;                                   // [0, 11]
	const unsigned d = doy - (153 * mp + 2) / 5 + 1;                             // [1, 31]
	const unsigned m = mp < 10 ? mp + 3 : mp - 9;                            // [1, 12]
	int yr = y + (m <= 2);
	int mn = m;
	 int day = d;
	 
	 return YMDtoF32time(yr, mn, day);

}

void date_jump(int y, int m, int d, int jumpDays, int* y1, int* m1, int* d1) {
	int days=days_from_civil(y, m, d);
	civil_from_days(days + jumpDays, y1, m1, d1);
}



static int __FindPatternStart( char *str, char * token) {

	char * pchar = strstr(str, token);
	if (pchar)	{
		return (int)(pchar - str);
	}
	return -10000;
}

int    GetStrPattern_fmt1(char* fmtstr, DateFmtPattern1* pattern) {
	ToUpper(fmtstr);

	int yearIdx = __FindPatternStart(fmtstr, "YYYY");
	if (yearIdx <0) return 0;
	int monIdx  = __FindPatternStart(fmtstr, "MM");
	if (monIdx < 0) return 0;
	int dayIdx  = __FindPatternStart(fmtstr, "DD");
	if (dayIdx < 0) return 0;

	pattern->yearIdx = yearIdx;
	pattern->monIdx  = monIdx;
	pattern->dayIdx  = dayIdx;
	return 1;
}

float  Str2F32time_fmt1(char* datestr, DateFmtPattern1* pattern) {

	char s[5];
	memcpy(s, datestr + pattern->yearIdx, 4);	s[4] = 0;	int year = atoi(s);
	if (year == 0) { return -1e10; }

	memcpy(s, datestr + pattern->monIdx, 2); 	s[2] = 0;	int mon = atoi(s);
	if (mon < 1 || mon > 12) { return -1e10; }

	memcpy(s, datestr + pattern->dayIdx, 2);  	s[2] = 0;	int day = atoi(s);
	if (day < 1 || day >31) { return -1e10; }

	return YMDtoF32time(year, mon, day);
}


int    GetStrPattern_fmt2(char* fmtstr, DateFmtPattern2* pattern) {
	ToUpper(fmtstr);
	int yearIdx = __FindPatternStart(fmtstr, "YYYY");
	if (yearIdx <0) return 0;
	int doyIdx  = __FindPatternStart(fmtstr, "DOY");
	if (doyIdx < 0) return 0;
	
	pattern->yearIdx = yearIdx;
	pattern->doyIdx  = doyIdx;
	return 1;
}

float  Str2F32time_fmt2(char* datestr, DateFmtPattern2* pattern) {

	char s[5];
	memcpy(s, datestr + pattern->yearIdx, 4);	s[4] = 0;	int year = atoi(s);
	if (year == 0) { return -1e10; }

	memcpy(s, datestr + pattern->doyIdx, 3); 	s[3] = 0;	int doy = atoi(s);	
	if (doy < 0 || doy > 366) { return -1e10; }	 

	return YDOYtoF32time(year, doy);
}

 
static char* _FindCharOccurrence(char* s, char c, int* numTimes) {
	
	*numTimes = 0;	 	
	char *pLast =NULL;
	while (( s=strchr(s, c)) != NULL ) {
		pLast = s++;
		++*numTimes;				
	}	
	return pLast;
	
}

static void insertionSort(void * arr[], char *index, int n)
{
	int i,  j;	
	for (i = 1; i < n; i++)
	{
		void* key = arr[i];
		char  idx = index[i];
		j = i - 1;
		/* Move elements of arr[0..i-1], that are
		greater than key, to one position ahead
		of their current position */
		while (j >= 0 && arr[j] > key)	{
			arr[j + 1]   = arr[j];
			index[j + 1] = index[j];
			j = j - 1;
		}
		arr[j + 1]   = key;
		index[j + 1] = idx;

		//index[j + 1] = index[j];
	}
}

int    GetStrPattern_fmt3(char* fmtstr, DateFmtPattern3* pattern) {
    
	ToUpper(fmtstr);
	int nTimes;
	char* yearPt = _FindCharOccurrence(fmtstr, 'Y', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	char* monPt = _FindCharOccurrence(fmtstr, 'M', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	char* dayPt = _FindCharOccurrence(fmtstr, 'D', &nTimes);
	if (nTimes == 0 || nTimes > 1) return 0;

	pattern->order[0] = 'Y';
	pattern->order[1] = 'M';
	pattern->order[2] = 'D';

	char* pts[] = { yearPt, monPt, dayPt };
	insertionSort(pts, pattern->order, 3);

	int64_t len;
	len =  (pts[1] - 1) - (pts[0] + 1) + 1;
	if (len <= 0) return 0;
	memcpy(pattern->sep1, pts[0] + 1, len); 	pattern->sep1[len] = 0;


	len =  (pts[2] - 1) - (pts[1] + 1) + 1 ;
	if (len <= 0) return 0;
	memcpy(pattern->sep2, pts[1] + 1, len);  	pattern->sep2[len] = 0;
	 
	return 1;
}


float  Str2F32time_fmt3(char* datestr, DateFmtPattern3* pattern) {

	int   N = (int) strlen(datestr);
	char  old;

	char* p0 = datestr;
	char *p1 = strstr(p0, pattern->sep1);
	if (p1 == NULL) return -1e10;
	old =p1[0]; p1[0] = 0;
	int n1=atoi(p0);
	p1[0] = old;

	p0 = p1 + strlen(pattern->sep1);
	p1 = strstr(p0, pattern->sep2);
	if (p1 == NULL) return -1e10;
	old = p1[0]; p1[0] = 0;
	int n2 = atoi(p0);
	p1[0] = old;

	p0 = p1 + strlen(pattern->sep2);
	if (p0 >= datestr+N) return -1e10;
	int n3 = atoi(p0);

	char* p = pattern->order;
	int year = p[0] == 'Y' ? n1 : (p[1] == 'Y' ? n2 : n3);
	int mon  = p[0] == 'M' ? n1 : (p[1] == 'M' ? n2 : n3);
	int day  = p[0] == 'D' ? n1 : (p[1] == 'D' ? n2 : n3);
	 
	return YMDtoF32time(year, mon, day);
}

#include "abc_000_warning.h"