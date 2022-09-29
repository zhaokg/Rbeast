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
//#pragma warning(suppress: 4483)
//extern void __identifier("?ioFlush@@YA_NXZ")(void);

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
	else if (IsChar(infield))		 return 0;
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
#if R_INTERFACE==1

#include <string.h>
#include <stdio.h>

//char t[] = "\033[0;35m";
//fflush(stdout);
//R_FlushConsole(); //c function: fflush(stdout)--flush the line buffer to see immediate outputs
//Rf_GetOption

SEXP getListElement(SEXP list, const char *str)
{
	SEXP elmt  = NULL; // R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
		elmt = VECTOR_ELT(list, i);
		break;
	}
	return elmt;
}
static int GetFieldIndex(SEXP list, const char* str)
{
	SEXP elmt  = NULL; // R_NilValue;
	int  index = -1L;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt  = VECTOR_ELT(list, i);
			index = i;
			break;
		}
	return index;
}
SEXP getListElement_CaseIn(SEXP list, const char *str)
{
	SEXP elmt = NULL; // R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
	if (strcicmp(CHAR(STRING_ELT(names, i)), str) == 0) {
		elmt = VECTOR_ELT(list, i);
		break;
	}
	return elmt;
}
 

VOID_PTR GetFieldByIdx(VOID_PTR strucVar, I32 ind) { return VECTOR_ELT(strucVar, ind); }
I32      GetConsoleWidth()  { return (I32) GetOptionWidth(); }

I32 GetCharArray(void *ptr, char * dst, int n) {

	if (TYPEOF((SEXP)ptr)!= STRSXP) return 0;

	char *tmp = CHAR(STRING_ELT(ptr, 0));	
	strncpy(dst, tmp, n);	
	dst[n] = 0;
	return (I32) strlen(dst);
	
}
I32 GetCharVecElem(void* ptr, int idx, char* dst, int n) {
	if (TYPEOF((SEXP)ptr) != STRSXP) return 0;

	char* tmp = CHAR(STRING_ELT((SEXP)ptr, idx));
	strncpy(dst, tmp, n);
	dst[n] = 0;
	return (I32)strlen(dst);
}

I32  GetNumberOfFields(const void* structVar) { return Rf_length (structVar); }

void * GetField(const void * structVar, char *fname) {

	if (structVar==NULL || structVar == R_NilValue )
		return NULL;

	void * elem = (void*)getListElement(structVar, fname);
	if (elem == NULL) {
		elem= (void*)getListElement_CaseIn(structVar, fname);
	}
	return elem;
}

void* GetField123(const void* structVar, char* fname, int nPartial) {
	if (!structVar) return NULL;

	void* elem = (void*)getListElement_CaseIn(structVar, fname);
	if (elem == NULL) {	 
		SEXP names = getAttrib(structVar, R_NamesSymbol);
		for (int i = 0; i < length(structVar); i++)
			if (strcicmp_nfirst(CHAR(STRING_ELT(names, i)), fname,nPartial) == 0) {
				elem = VECTOR_ELT(structVar, i);
				break;
			} 
	}
	return elem;
}
 

#define IS_SCALAR(x, type) (TYPEOF((SEXP)x) == (type) && XLENGTH((SEXP)x) == 1)
F64    GetScalar(const void * ptr) { 
	if (TYPEOF((SEXP)ptr) == INTSXP)
		return (F64)asInteger((SEXP)ptr);
	else if (TYPEOF((SEXP)ptr) == REALSXP)
		return (F64)asReal((SEXP)ptr);
	else if (TYPEOF((SEXP)ptr) == LGLSXP)
		return asLogical((SEXP)ptr);

	return getNaN();
}
void * GetData(const void * ptr) { 
	if (TYPEOF((SEXP)ptr) == INTSXP)
		return INTEGER((SEXP)ptr);
	else if (TYPEOF((SEXP)ptr) == REALSXP)
		return REAL((SEXP)ptr);
	else if (TYPEOF((SEXP)ptr) == LGLSXP)
		return LOGICAL((SEXP)ptr);

	return NULL;
}
int   GetDim1(const void * ptr) {
	SEXP dims = PROTECT( getAttrib((SEXP)ptr, R_DimSymbol) );
	int dim1 = INTEGER(dims)[0];	
	UNPROTECT(1);
	return dim1;
}
int   GetDim2(const void * ptr) {
	SEXP dims = PROTECT(getAttrib((SEXP)ptr, R_DimSymbol));
	int  dim2 = INTEGER(dims)[1];
	UNPROTECT(1);
	return dim2;
}
int   GetNumOfDim(const void * ptr) {
	SEXP dims  = PROTECT(getAttrib(ptr, R_DimSymbol));
	int  ndims = Rf_length(dims);
	UNPROTECT(1);
	return ndims;
}
void GetDimensions(const void * ptr, int dims[], int ndims) {
	int  N    = min(ndims, GetNumOfDim(ptr));
	SEXP DIMS = PROTECT(getAttrib((SEXP)ptr, R_DimSymbol));
	for (int i = 0; i < N; i++)
	{
		dims[i] = INTEGER(DIMS)[i];
	}
	UNPROTECT(1);
}

int GetNumberOfElements(const void * ptr) {	return Rf_length((SEXP)ptr);}

void * GetCellElement(const void * ptr, I32 idx) {		return NULL; }
 
// if (!(TYPEOF(Y) == INTSXP  && XLENGTH(Y) > 2) 
//     &&  !(TYPEOF(Y) == REALSXP && XLENGTH(Y) > 2) 
//     &&  !(TYPEOF(Y) == INTSXP  && isMatrix(Y))     
//     &&  !(TYPEOF(Y) == REALSXP && isMatrix(Y)) 
//     && !(TYPEOF(Y) == STRSXP  && XLENGTH(Y) == 1) 
//     &&  !(TYPEOF(Y) == REALSXP  && isArray(Y)) 		)	

int IsClass(void* ptr, char* class) {
	SEXP klass;	
 
	if (OBJECT(ptr)) {
		klass = getAttrib(ptr, R_ClassSymbol);
		for (int i = 0; i < length(klass); i++) {
			if (strcmp(CHAR(STRING_ELT(klass, i)), class) == 0)
				return 1;
		}		
	}
	return 0;
}
int IsCell(void* ptr)    { return 0L; }
int IsChar(void* ptr)    { return TYPEOF((SEXP)ptr) == STRSXP; }
int IsEmpty(void* ptr)   { return ptr == R_NilValue || GetNumberOfElements(ptr)==0; }
int IsStruct(void* ptr)  { return isNewList((SEXP)ptr) || ptr==R_NilValue;       }
int IsNumeric(void* ptr) { return isNumeric((SEXP)ptr); }

int IsDouble(void* ptr)  { return TYPEOF((SEXP)ptr) == REALSXP; }
int IsSingle(void* ptr)  { return 0; }
int IsInt32(void* ptr)   { return TYPEOF((SEXP)ptr) == INTSXP;  }
int IsInt16(void* ptr) { return 0; }
int IsInt64(void* ptr)   { return 0; }
//http://adv-r.had.co.nz/C-interface.html
int IsLogical(void* ptr) { return TYPEOF((SEXP)ptr) == LGLSXP; };

void *CreateNumVar(DATA_TYPE dtype, int *dims, int ndims, VOIDPTR * data_ptr) {
	         
	        int rtype;
			if (dtype == DATA_INT32) 
				rtype = INTSXP; 
			else if (dtype == DATA_DOUBLE) 
				rtype = REALSXP;			
			else {
				assert(0);
			}


			SEXP     tmpSEXP=NULL;
	 
			SEXP dims4d;
			switch (ndims) {
			case 1:
				//Actually even no need to protect it bcz after alloCVector, no new allocation is called
				PROTECT( tmpSEXP = allocVector(rtype, dims[0]) );
				UNPROTECT(1);
				break;
			case 2:
				PROTECT(tmpSEXP = allocMatrix(rtype,  dims[0], dims[1]) );				
				UNPROTECT(1); 
				break;
			case 3:
				PROTECT(tmpSEXP = alloc3DArray(rtype, dims[0], dims[1], dims[2]));				
				UNPROTECT(1);
				break;
			case 4:
				
				PROTECT(dims4d = allocVector(INTSXP, 4));
				INTEGER(dims4d)[0] = dims[0];
				INTEGER(dims4d)[1] = dims[1];
				INTEGER(dims4d)[2] = dims[2];
				INTEGER(dims4d)[3] = dims[3];
				PROTECT( tmpSEXP = allocArray(rtype, dims4d) );				
				UNPROTECT(2);
				break;
			}

			if (data_ptr != NULL && tmpSEXP !=NULL) {

				if      (rtype == INTSXP) 	data_ptr[0] = INTEGER(tmpSEXP);		 
				else if (rtype == REALSXP) 	data_ptr[0] = REAL(tmpSEXP);
				else {
					assert(0);
					data_ptr[0] = NULL;
				} 
			}

			return tmpSEXP;
}

void *CreateStructVar(FIELD_ITEM *fieldList, int nfields) { 	

	// Find the sentinal element to update nfields
	int nfields_new = 0;
	for (int i = 0; i < nfields; ++i) {
		nfields_new++;
		if (fieldList[i].name[0] == 0) 	break;
	}
	nfields = nfields_new;

// www.mathworks.com/help/matlab/apiref/mxsetfieldbynumber.html	
	
///https://github.com/kalibera/rchk/issues/21
//Function CreateStructVar
//	[UP] protect stack is too deep, unprotecting all variables, results will be incomplete
// [UP] unsupported form of unprotect, unprotecting all variables, results will be incomplete / home / docker / R - svn / packages / build / RWzZRvbc / Rbeast / src / abc_ide_util.c:446

//https://developer.r-project.org/Blog/public/2019/04/18/common-protect-errors/

//https://stackoverflow.com/questions/59050384/protecting-elements-in-vecsxp
//To do the analysis, you need to know every point in your code that does allocations, because an 
//allocation might trigger garbage collection, and that will release any unprotected object.

//..Protecting an R object automatically protects all the R objects pointed to in the
//corresponding SEXPREC, for example all elements of a protected list are automatically
//protected.



	//https://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
	SEXP LIST;
	SEXP NAMES;

	//PROTECT(pY = coerceVector(pY, REALSXP));  nprt++;
	//Writing R extensons: Protection is not needed for objects which R already knows are in use.In particular, this applies to function arguments

	int  nprt = 0L;

	PROTECT(LIST  = allocVector(VECSXP, nfields));   ++nprt;
	PROTECT(NAMES = allocVector(STRSXP, nfields));   ++nprt;

	for (I64 i = 0; i < nfields; i++)
		SET_STRING_ELT(NAMES, i, mkChar(fieldList[i].name));
	
	SEXP tmpSEXP;
	for (I32 i = 0; i < nfields; i++) 	{	 

		if (fieldList[i].ptr == NULL) {
			// no mem will be alloacted if ptr==NULL
			SET_VECTOR_ELT(LIST, i, R_NilValue);
			continue;
		} 
		
		if (fieldList[i].type == DATA_STRUCT) {			
			// ptr should be an ptr to an existing SEXP object
			tmpSEXP = fieldList[i].ptr;
			SET_VECTOR_ELT(LIST, i, tmpSEXP);
		}else  	{		
			// new mem will be allocated and its pointer is saved to *ptr
			tmpSEXP = PROTECT( CreateNumVar(fieldList[i].type, fieldList[i].dims, fieldList[i].ndim, fieldList[i].ptr ));
			SET_VECTOR_ELT(LIST, i, tmpSEXP);
			UNPROTECT(1);		
		}  
		
		 
	}//for (rI64 i = 0; i < nfields; i++) 

	setAttrib(LIST, R_NamesSymbol, NAMES);

	/*
	PROTECT(tmpSEXP = mkString("beast"));  ++nprt;
	setAttrib(LIST, R_ClassSymbol, tmpSEXP);

	PROTECT(tmpSEXP = mkString("beastST")); ++nprt;
	setAttrib(LIST, install("algorithm"), tmpSEXP);
	*/
	UNPROTECT(nprt);
	return (void*)LIST;
}

void ReplaceStructField(VOIDPTR s, char* fname, VOIDPTR newvalue){
	int index= GetFieldIndex(s,fname);
	if (index <0) {
		return;
	}
	SEXP oldvalue = VECTOR_ELT(s, index);

	SET_VECTOR_ELT(s, index, newvalue);
}

void  DestoryStructVar(VOID_PTR strutVar) {

	// Do nothing because in R, varialbes are protected/unprocted. 
	//So, to delete an object, just unprotect it so it can be gc-ed.

}

//https://stackoverflow.com/questions/19688380/inconsistency-when-getting-class-attribute-of-different-r-objects-using-c-functi
void AddStringAttribute(VOID_PTR listVar, const char * field, const char *value) {

	SEXP  tmp=PROTECT(mkString(value)); 
	if (strcmp(field, "class") == 0) {
		setAttrib(listVar, R_ClassSymbol, tmp);
	}else{
		setAttrib(listVar, install(field), tmp);
	}
	UNPROTECT(1);
}
void AddIntegerAttribute(VOID_PTR listVar, const char* field, I32 value) {

	SEXP  tmp = PROTECT(ScalarInteger(value));
	setAttrib(listVar, install(field), tmp);	
	UNPROTECT(1);
}
void RemoveAttribute(VOID_PTR listVar, const char* field) {

	SEXP  sym = PROTECT( install(field)  );

	SEXP  attr = getAttrib(listVar, sym);
	if (!Rf_isNull(attr)) 
		setAttrib(listVar, sym, R_NilValue);

	UNPROTECT(1);
}


//https://stackoverflow.com/questions/40563522/r-how-to-write-interruptible-c-function-and-recover-partial-results
//https://github.com/tidyverse/dplyr/issues/123
static void __chkIntFn(void *dummy) {	R_CheckUserInterrupt();}
// this will call the above in a top-level context so it won't longjmp-out of your context 
I32  CheckInterrupt()        {	return (R_ToplevelExec(__chkIntFn, NULL) == FALSE);}
void ConsumeInterruptSignal() { return ; }
 
#if defined(WIN64_OS)  && 0
// R doesn't allow calling no-API entry points: get_R_HOME
	#define WIN32_LEAN_AND_MEAN
	#include "windows.h"
	#include "Rembedded.h" // for get_R_HOME only

	static void * GetReadConsole(void) 	{
		int(*R_ReadConsole)(const char* prompt, unsigned char* buf, int len, int addtohistory);

		R_ReadConsole = NULL;

		char  dll[1000];
		char *RHOME = get_R_HOME();		
		strncpy(dll, RHOME, 1000);
		strcat(dll, "/bin/x64/R.dll");
	
		HINSTANCE hinstLib = LoadLibrary(dll);
		if (hinstLib != NULL) {			
			R_ReadConsole = GetProcAddress(hinstLib, "R_ReadConsole");
			FreeLibrary(hinstLib);
			//r_printf("%#010x\n", R_ReadConsole);
		}

		/*		
		char s[10];
		while (stricmp(s, "exit")!=0) {	R_ReadConsole("", s, 1, 0);	s[4] = 0;	}		
		*/
		return R_ReadConsole;
	}


	 void* GetR_interrupts_pending(void) {
		char* RHOME = get_R_HOME();
		char  dll[1000];
		strncpy(dll, RHOME, 1000);
		strcat(dll, "/bin/x64/R.dll");

		HINSTANCE hinstLib;
		hinstLib = LoadLibrary(dll);
		int(*R_ReadConsole)(const char* prompt, unsigned char* buf, int len, int addtohistory);
		void * R_interrupts_pending = GetProcAddress(hinstLib, "R_interrupts_pending");
		FreeLibrary(hinstLib);

		return R_interrupts_pending;
	}


	bool GetUserInput(void * nullPtr)
	{
		int (*R_ReadConsole)(const char *prompt, unsigned char *buf, int len, int addtohistory);
		R_ReadConsole = GetReadConsole();
		char str[100];
		str[0]=0;

		while (stricmp(str, "exit") != 0) {
			R_ReadConsole("", str, 10, 0);
			str[4] = 0;
		}
		r_printf("Program is forcely terminated ...");
		IDE_USER_INTERRUPT = 1L;
	}
#endif
 
#elif M_INTERFACE==1


//https://stackoverflow.com/questions/25998442/mex-options-when-compiling-with-visual-studio
//https://undocumentedmatlab.com/articles/mex-ctrl-c-interrupt
// The following 2 function are exported from libut.lib
bool utIsInterruptPending();
void utSetInterruptPending(bool);

I32 CheckInterrupt() { return utIsInterruptPending(); }
void ConsumeInterruptSignal() { utSetInterruptPending(false); }


I32  GetNumberOfFields(const void* structVar) { return mxGetNumberOfFields(structVar); }
//https://www.mathworks.com/help/matlab/apiref/mxgetfieldbynumber.html
VOID_PTR GetFieldByIdx(VOID_PTR ptr, I32 ind) {

	if (IsCell(ptr)) {
		return  mxGetCell(ptr, ind);
	}
	if (IsStruct(ptr)) {
		mxGetFieldByNumber(ptr, 0, ind);
	}
	return NULL;
}

 

I32 GetCharArray(void* ptr, char* dst,   int n) {
	//Dsr needs to be at least of length (n+1);

	int len = 0;
	if (mxIsChar(ptr)){		
		// This func allocates new mem for tmp, which needs to be explicitlly dellocated later
		//http://www.ece.northwestern.edu/local-apps/matlabhelp/techdoc/apiref/mxarraytostring.html#976726
		char* tmp = mxArrayToString(ptr);
		strncpy(dst, tmp, n); dst[n] = 0;	len  = strlen(dst);
		r_free(tmp);
		return len;
	}
	else if (mxIsClass(ptr, "string")) {
		//www.mathworks.com/matlabcentral/answers/330929-how-to-access-matlab-string-data-in-mex-c
		mxArray *tmpMxChar = mxGetProperty(ptr, 0, "data");
		if (tmpMxChar != NULL) {
			char* tmp = mxArrayToString(tmpMxChar);
			strncpy(dst, tmp, n); dst[n] = 0;	len  = strlen(dst);
			r_free(tmp);
			return len;
		}
		else {
			//Matlab's String class is encapsulated,
		//use Matlab call to convert it to char array
			mxArray* string_class[1] = {ptr}, *char_array[1];			
			mexCallMATLAB(1, char_array, 1, string_class, "char");
			char* tmp  = mxArrayToString(char_array[0]);
			mxDestroyArray(char_array[0]);

			strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
			r_free(tmp);
			return len;			
		}
	}
	else {
		return 0L;
	}
	
}
I32 GetCharVecElem(void* ptr, int idx, char* dst, int n) {
	//Dsr needs to be at least of length (n+1);

	int len = 0;

	if (IsCell(ptr) ) {
		
		void* elem = mxGetCell(ptr, idx);

		if (mxIsChar(elem)) {

			// This func allocates new mem for tmp, which needs to be explicitlly dellocated later
			// www.ece.northwestern.edu/local-apps/matlabhelp/techdoc/apiref/mxarraytostring.html#976726
			char* tmp = mxArrayToString(elem);
			strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
			r_free(tmp);
			return len;

		}
		else if (mxIsClass(elem, "string")) {
			//www.mathworks.com/matlabcentral/answers/330929-how-to-access-matlab-string-data-in-mex-c
			mxArray* tmpMxChar = mxGetProperty(elem, 0, "data");
			if (tmpMxChar != NULL) {
				char* tmp = mxArrayToString(tmpMxChar);
				strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
				r_free(tmp);
				return len;
			}
			else {
				//Matlab's String class is encapsulated,
			//use Matlab call to convert it to char array
				mxArray* string_class[1] = { elem }, * char_array[1];
				mexCallMATLAB(1L, char_array, 1L, string_class, "char");
				char* tmp = mxArrayToString(char_array[0]);
				mxDestroyArray(char_array[0]);

				strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
				r_free(tmp);
				return len;
			}
		}
		else {
			return 0L;
		}

	}
	else if (mxIsClass(ptr, "string") && !IsCell(ptr)) {

		//["a","b","c"...]

		void* elem = (void**)mxGetData(ptr) + idx;

		//www.mathworks.com/matlabcentral/answers/330929-how-to-access-matlab-string-data-in-mex-c
		mxArray* tmpMxChar = mxGetProperty(elem, 0, "data");
		if (tmpMxChar != NULL) {
			char* tmp = mxArrayToString(tmpMxChar);
			strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
			r_free(tmp);
			return len;
		}
		else {
			//Matlab's String class is encapsulated,
		//use Matlab call to convert it to char array
			mxArray* string_class[1] = { elem }, * char_array[1];
			mexCallMATLAB(1L, char_array, 1L, string_class, "char");
			char* tmp = mxArrayToString(char_array[0]);
			mxDestroyArray(char_array[0]);

			strncpy(dst, tmp, n); dst[n] = 0;	len = strlen(dst);
			r_free(tmp);
			return len;
		}
	}
	else {
		return 0L;
	}

}
void * GetField(const void * structVar, char *fname) {	
	
	if (structVar==NULL || mxIsEmpty(structVar)) {
		return NULL;
	}

	VOIDPTR ptr = (VOIDPTR) mxGetField(structVar, 0, fname);
	if (ptr != NULL) {
		return ptr;
	}
	
	I32 numFlds = mxGetNumberOfFields(structVar);
	for (I32 idx = 0; idx < numFlds; idx++) {
		char* tmpName = mxGetFieldNameByNumber(structVar, idx);
		if (strcicmp(fname, tmpName) == 0) {
			return mxGetFieldByNumber(structVar, 0, idx);
		}
	}

	return NULL;	
}
void* GetField123(const void* structVar, char* fname, int nPartial) {
	
	if (structVar == NULL || mxIsEmpty(structVar)) {
		return NULL;
	}

	VOIDPTR ptr = (VOIDPTR)mxGetField(structVar, 0, fname);
	if (ptr != NULL) {
		return ptr;
	}

	I32 numFlds = mxGetNumberOfFields(structVar);
	for (I32 idx = 0; idx < numFlds; idx++) {
		char* tmpName = mxGetFieldNameByNumber(structVar, idx);
		if (strcicmp_nfirst(fname, tmpName, nPartial) == 0) {
			return mxGetFieldByNumber(structVar, 0, idx);
		}
	}

	return NULL;
}
F64    GetScalar(const void * ptr) { 	return mxGetScalar((mxArray *)ptr); }
void * GetData(const void * ptr) {   return mxGetData((mxArray *)ptr); }
int    GetDim1(const void * ptr) { return mxGetM((mxArray *)ptr); }
int    GetDim2(const void * ptr) { return mxGetN((mxArray *)ptr); }
int    GetNumOfDim(const void * ptr) { return mxGetNumberOfDimensions((mxArray *)ptr); }
void   GetDimensions(const void * ptr, int dims[], int ndims) {
	int N = min(ndims, GetNumOfDim(ptr) );
	const mwSize *dimArr = mxGetDimensions((mxArray *)ptr);
	for (int i = 0; i < N; i++)
	{
		dims[i] = dimArr[i];
	}	
}
int    GetNumberOfElements(const void * ptr)
{ // Do not work for a vector of opaque objects for some reason 
	return mxGetNumberOfElements(ptr);
}
 
int IsChar(void* ptr)    { 
	if (IsCell(ptr)) {
		// if ptr is a cell, iterate through all elements and return 1L if all the elems are char.
		int n = GetNumberOfElements(ptr);
		for (int i = 0; i < n; ++i) {
			void *tmp = mxGetCell(ptr, i);
			if ( !mxIsChar(tmp) && !mxIsClass(tmp, "string")) return 0; 
		}
		return 1L;		
	}
	else {
		return mxIsChar(ptr) || mxIsClass(ptr, "string");
	}
	
}
int IsClass(void* ptr, char* class) { return 0; }
int IsEmpty(void* ptr) { return mxIsEmpty(ptr); }
int IsStruct(void* ptr)  { return mxIsStruct(ptr) ||  mxIsEmpty(ptr); } // We simply treat an empty mxArray as a struct, though matlab treats it as numeric
int IsCell(void* ptr)    { return mxIsCell(ptr); }
int IsNumeric(void* ptr) { return mxIsNumeric(ptr) && !mxIsEmpty(ptr); } // Matlab treates an empty as Numeric
int IsDouble(void* ptr)  { return mxIsDouble(ptr); }
int IsSingle(void* ptr)  { return mxIsSingle(ptr); }
int IsInt32(void* ptr)   { return mxIsInt32(ptr); }
int IsInt16(void* ptr) { return mxIsInt16(ptr); }
int IsInt64(void* ptr) { return mxIsInt64(ptr); }
int IsLogical(void* ptr) { return mxIsLogical(ptr); }

void * CreateStructVar(FIELD_ITEM *fieldList, int nfields)
{ // www.mathworks.com/help/matlab/apiref/mxsetfieldbynumber.html	
	
		// Find the sentinal element to update nfields
	int nfields_new = 0;
	for (int i = 0; i < nfields; ++i) {
		nfields_new++;
		if (fieldList[i].name[0] == 0) 	break;
	}
	nfields = nfields_new;


	mxArray * _restrict out;
	{
		char * fldNames[100];		
		for (int i = 0; i < nfields; i++) { fldNames[i] = fieldList[i].name; }					
		mwSize dims_2d[2] = { 1,1 };
		out = mxCreateStructArray(2, dims_2d, nfields, fldNames);
	}

	for (int i = 0; i < nfields; i++){	 

		if (fieldList[i].ptr == NULL) continue;

		if (fieldList[i].type == DATA_STRUCT){
			mxSetField(out, 0L, fieldList[i].name, (mxArray *) fieldList[i].ptr);
			continue;
		}

		// Not a struct variable to be created
		mxClassID fieldDataType;
		switch (fieldList[i].type)	{
			case DATA_FLOAT: 	fieldDataType = mxSINGLE_CLASS;		break;
			case DATA_DOUBLE:   fieldDataType = mxDOUBLE_CLASS;		break;
			case DATA_INT32:    fieldDataType = mxINT32_CLASS;		break;
			default:			fieldDataType = mxSINGLE_CLASS;			
		}

		//mxCreateDoubleScalar
		//mxCreateDoubleMatrix
		//mxCreateNumericArray
		//mxCreateNumericArray
		
		mxArray * _restrict mxptr=NULL;
		if (fieldList[i].ndim == 1) { 
			// Make it a column vector if it is a 1D. In Matlab, there is no concpet of dimensonless vector
			// But in R, there is.
			mxptr = mxCreateNumericMatrix(fieldList[i].dims[0], 1L, fieldDataType, mxREAL);
		}else if ( fieldList[i].ndim == 2) {
			mxptr = mxCreateNumericMatrix(fieldList[i].dims[0], fieldList[i].dims[1], fieldDataType, mxREAL);
		} else if (fieldList[i].ndim >= 3)	{
			mwSize DIMS[4] = {fieldList[i].dims[0], fieldList[i].dims[1], fieldList[i].dims[2], fieldList[i].dims[3] };
			mxptr = mxCreateNumericArray(fieldList[i].ndim, DIMS, fieldDataType, mxREAL);
		}
		mxSetField(out, 0L, fieldList[i].name, mxptr);
		*(fieldList[i].ptr) = mxGetData(mxptr);
	}

	return (void*)out;
}

//https://www.mathworks.com/help/matlab/apiref/mxdestroyarray.html
void  DestoryStructVar(VOID_PTR strutVar) {
		mxDestroyArray(strutVar);
}

void AddStringAttribute(VOID_PTR listVar, const char* field, const char* value) {

	mxArray* tmp = mxCreateString(value);

	mxAddField(listVar, field);	
	mxSetField(listVar, 0, field, tmp);

}
void AddIntegerAttribute(VOID_PTR listVar, const char* field, I32 value) {

	mxArray* tmp = mxCreateDoubleScalar(value);

	mxAddField(listVar, field);
	mxSetField(listVar, 0, field, tmp);

}
void RemoveAttribute(VOID_PTR listVar, const char* field) {
	 
}

I32 GetConsoleWidth()
{
	mxArray *pOut[1], *pIn[2];
	pIn[0] = mxCreateDoubleScalar(0);
	pIn[1] = mxCreateString("CommandWindowSize");
	mexCallMATLAB(1, pOut, 2, pIn, "get");
	F64 *ptr = mxGetData(pOut[0]);
	I32 screenWidth = ptr[0];
	mxDestroyArray(pIn[0]);
	mxDestroyArray(pIn[1]);
	mxDestroyArray(pOut[0]);
	return screenWidth;
}
#endif

int  HaveEqualDimesions(const void* p1, const void* p2) {
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
	else if (IsDouble(Y))  /* isReal(pY)*/		return DATA_DOUBLE;
	else if (IsSingle(Y))  /* isReal(pY)*/ 		return DATA_FLOAT;
	else                                        return DATA_UNKNOWN;
}
F64  GetNumericElement(const void* Y, I32 idx) {

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
			if (IsInt32(Y)) 						    return *((I32PTR)y+idx);
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