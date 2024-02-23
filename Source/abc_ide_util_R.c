#include <string.h>
#include "assert.h"

#include "abc_000_warning.h"
#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"

#include<stdio.h>  // Need _GNU_SOURCE for manylinux; otherwise report /usr/include/stdio.h:316:6: error: unknown type name '_IO_cookie_io_functions_t'

//char t[] = "\033[0;35m";
// fflush(stdout);   //: c function: fflush(stdout)--flush the line buffer to see immediate outputs
// Rf_GetOption


#if R_INTERFACE == 1


int  JDN_to_DateNum(int jdn) {
	return jdn - 2440588;
}



void StdouFlush(void) {
// https://stackoverflow.com/questions/13327305/making-rcout-output-appear-on-the-r-console-immediately
//	R_FlushConsole(): a R functin to flush the print
	R_FlushConsole();
}

SEXP getListElement(SEXP list, const char* str) {

	SEXP elmt = NULL; // R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(names); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}


SEXP getListElement_CaseIn(SEXP list, const char* str)
{
	SEXP elmt = NULL; // R_NilValue;
	SEXP names = getAttrib(list, R_NamesSymbol); 
	for (int i = 0; i < length(names); i++)
		if (strcicmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}
	return elmt;
}

I32      GetConsoleWidth() {
	return (I32)GetOptionWidth();
}

// if (!(TYPEOF(Y) == INTSXP  && XLENGTH(Y) > 2) 
//     &&  !(TYPEOF(Y) == REALSXP && XLENGTH(Y) > 2) 
//     &&  !(TYPEOF(Y) == INTSXP  && isMatrix(Y))     
//     &&  !(TYPEOF(Y) == REALSXP && isMatrix(Y)) 
//     && !(TYPEOF(Y) == STRSXP  && XLENGTH(Y) == 1) 
//     &&  !(TYPEOF(Y) == REALSXP  && isArray(Y)) 		)	

int IsClass(void* ptr, char* class) {

	if (OBJECT(ptr)) {
		SEXP klass = getAttrib(ptr, R_ClassSymbol);
		for (int i = 0; i < length(klass); i++) {
			if (strcmp(CHAR(STRING_ELT(klass, i)), class) == 0) {
				return 1;
			}
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
int IsInt16(void* ptr)   { return 0; }
int IsInt64(void* ptr)   { return 0; }
//http://adv-r.had.co.nz/C-interface.html
int IsLogical(void* ptr) { return TYPEOF((SEXP)ptr) == LGLSXP; }

VOID_PTR GetFieldByIdx(VOID_PTR strucVar, I32 ind0) { 
	return VECTOR_ELT(strucVar, ind0); 
}

void GetFieldNameByIdx(VOID_PTR strucVar, I32 ind0, char* str, int buflen) {
 
	SEXP names = getAttrib(strucVar, R_NamesSymbol);
	int  n     = length(names);
	if (ind0 < n) {
		const char* name = CHAR(STRING_ELT(names, ind0));
		strncpy(str, name, buflen);
		str[buflen - 1] = 0;
	}	else {
		str[0] = 0;
	}
 	
}

I32 GetCharArray(void *ptr, char * dst, int n) {
	dst[0] = 0;
	if (ptr==NULL || TYPEOF((SEXP)ptr) != STRSXP) {
		return 0;
	}

	const char* tmpstr= CHAR(STRING_ELT(ptr, 0));
	strncpy(dst, tmpstr, n);
	dst[n] = 0;
	return (I32) strlen(dst);
	
}

I32 GetCharVecElem(void* ptr, int idx, char* dst, int n) {

	if (TYPEOF((SEXP)ptr) != STRSXP) {
		return 0;
	}

	const char* tmpstr = CHAR(STRING_ELT((SEXP)ptr, idx));
	strncpy(dst, tmpstr, n);
	dst[n] = 0;
	return (I32)strlen(dst);
}

I32  GetNumberOfFields(const void* structVar) { 
	return Rf_length (structVar); 
}

void * GetField(const void * structVar, char *fname) {

	if (structVar == NULL || structVar == R_NilValue) {
		return NULL;
	}

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
		SEXP names = getAttrib(structVar, R_NamesSymbol); // names may be NULL
		for (int i = 0; i < length(names); i++)
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
void * SetDimensions(const void* ptr, int dims[], int ndims) {
	
	if (!ptr) return NULL;

	SEXP  tmp = PROTECT(allocVector(INTSXP, ndims));	
	for (int i = 0; i < ndims; i++) 	{
		INTEGER(tmp)[i] = dims[i];
	}
	setAttrib(ptr, R_DimSymbol, tmp);
	UNPROTECT(1);

	return ptr;
}


int GetNumberOfElements(const void * ptr) {	
	return Rf_length((SEXP)ptr);
}


void *CreateNumVar(DATA_TYPE dtype, int *dims, int ndims, VOIDPTR * data_ptr) {
	         
	        int rtype;
			if (dtype == DATA_INT32) 
				rtype = INTSXP; 
			else if (dtype == DATA_DOUBLE) 
				rtype = REALSXP;			
			else {
				assert(0);
			}


			SEXP   tmpSEXP=NULL;
	 
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
	int nfields_actual = 0;
	for (int i = 0; i < nfields; ++i) {
		nfields_actual++;
		if (fieldList[i].name[0] == 0) {
			nfields_actual--;
			break;
		}
	}
	nfields = nfields_actual;
 
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
		} else  	{		
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


static int __GetFieldIndex(SEXP list, const char* str) {
	int  index = -1L;
	SEXP names = getAttrib(list, R_NamesSymbol);
	for (int i = 0; i < length(list); i++)
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			SEXP elmt = VECTOR_ELT(list, i);
			index = i;
			break;
		}
	return index;
}


void ReplaceStructField(VOIDPTR s, char* fname, VOIDPTR newvalue){

	int index= __GetFieldIndex(s,fname);
	if (index <0) {
		return;
	}

	
	// no need to porect old value; we want to let it get destoried
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
static void __chkIntFn(void *dummy) {R_CheckUserInterrupt();}
// this will call the above in a top-level context so it won't longjmp-out of your context 
I32  CheckInterrupt()         {	return (R_ToplevelExec(__chkIntFn, NULL) == FALSE);}
void ConsumeInterruptSignal() { return ; }
 
#if defined(OS_WIN64)  && 0
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
		while (strcicmp(s, "exit")!=0) {	R_ReadConsole("", s, 1, 0);	s[4] = 0;	}		
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


	 Bool GetUserInput(void * nullPtr)
	{
		int (*R_ReadConsole)(const char *prompt, unsigned char *buf, int len, int addtohistory);
		R_ReadConsole = GetReadConsole();
		char str[100];
		str[0]=0;

		while (strcicmp(str, "exit") != 0) {
			R_ReadConsole("", str, 10, 0);
			str[4] = 0;
		}
		r_printf("Program is forcely terminated ...");
		IDE_USER_INTERRUPT = 1L;
	}
#endif

#else

const static char achar = 'a';
#endif

#include "abc_000_warning.h"