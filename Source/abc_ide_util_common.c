#include <math.h>
#include <string.h>
#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"

 
// https:// stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
extern void matlab_IOflush(void);


// C++ solution to define a keyword as a function name
// #pragma warning(suppress: 4483)
// extern void __identifier("?ioFlush@@YA_NXZ")(void);

I08 IDE_USER_INTERRUPT;

void printProgress(F32 pct, I32 width, char * buf, I32 firstTimeRun)
{//https:// stackoverflow.com/questions/2685435/cooler-ascii-spinners

	
 	static char spinnerChar[] =  "|/-\\";
	static I32  cnt = 1;
	cnt++;
	cnt = cnt == 4 ? 0 : cnt;
	width = max(width, 35); //cause errors when width is 0 or negative
	memset(buf, '*', width); // space = 20

	I32  len = 0;
	buf[len++] = spinnerChar[cnt];

	char prefix[] = "Progress:";
	I32 strLen    = sizeof(prefix)-1L; //the last byte is a Zero
	memcpy(buf+len, prefix, strLen);
	len += strLen;

	sprintf(buf + len, "%5.1f%% done", pct * 100);
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
	Rprintf("\r%s", buf);

#elif P_INTERFACE==1
	r_printf("\r%s", buf);
	//R doesnto allow io operationrs from external libraries
	//fflush(stdout);
#elif M_INTERFACE==1
	if (firstTimeRun == 1)
	{
		r_printf("\r\n");
		r_printf("%s", buf);		
		matlab_IOflush();
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		//mexEvalString("pause(0.0000001);drawnow");
	}
	else
	{
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf);
		matlab_IOflush();
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		//mexEvalString("pause(0.0000001);drawnow");

	}


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

	sprintf(buf + len, "%5.1f%%", pct * 100);
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

	if (time > SecsPerDay)
		sprintf(buf + len, "<Remaining%02dday%02dhrs%02dmin>", days, hrs, mins);
	else
		sprintf(buf + len, "<Remaining%02dhrs%02dmin%02dsec>", hrs, mins, secs);
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

	if (firstTimeRun == 1)
	{
		r_printf("\r\n");
		r_printf("%s", buf);
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		matlab_IOflush();
	}
	else {
		char * back = buf + width + 5;
		memset(back, '\b', width + 2);
		back[width + 2] = 0;

		r_printf(back);
		r_printf("%s\r\n", buf);
		//mexEvalString("drawnow");
		//mexCallMATLAB(0, NULL, 0, NULL, "drawnow");
		matlab_IOflush();
	}
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

int CopyNumericArrToF32Arr(F32PTR outmem, VOID_PTR infield, int N) {

	VOID_PTR data = GetData(infield);

	if (IsSingle(infield))     		memcpy(outmem, data, sizeof(F32) * N);
	else if (IsDouble(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((double*)data + i);
	else if (IsInt32(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((int*)data + i);
	else if (IsInt64(infield))		for (I32 i = 0; i < N; i++) outmem[i] = *((I64*)data + i);
	else if (IsChar(infield))		return 0;
	else {	
		return 0;
	}
	return 1L;
}

int CopyNumericArrToI32Arr(I32PTR outmem, VOID_PTR infield, int N) {

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
	
	VOID_PTR p = GetField123(structVar, fname, nPartial);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}

}

void* GetFieldCheck(const void* structVar, char* fname) {
	// Check if the retured field is an empety or R_NilValue
	VOID_PTR p = GetField(structVar, fname);
	if (p == NULL || IsEmpty(p))
		return NULL;
	else {
		return p;
	}
}

#include "abc_000_warning.h"