#include <math.h>
#include <string.h>
#include "assert.h"
#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_common.h"
#include "abc_date.h"
#include "inttypes.h"  // Possibly needed for StdoutFlush()
// #include "mex.h"

////https://stackoverflow.com/questions/10529500/what-does-this-mean-int-a
// in C++, (int &) is a type punning that can force converstion of a variable to int

#if M_INTERFACE==1  


int  JDN_to_DateNum(int jdn) {
	return jdn - 1721059;
}


void StdouFlush(void) {

// https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running/26271557
// isoFlush is  an undocumented C++ function that resides in libmwservices.dll

// https://blog.csdn.net/saddlesad/article/details/119322795
// https://github.com/gchatelet/gcc_cpp_mangling_documentation
// Name mangling starts with __Z in MacOS and _Z in Linux

#if defined(OS_WIN32) || defined(OS_WIN64)
 #if defined(COMPILER_MSVC)
    // https://stackoverflow.com/questions/53381461/does-visual-c-provide-a-language-construct-with-the-same-functionality-as-a
    #pragma comment(linker, "/alternatename:ioFlush=?ioFlush@@YA_NXZ")
 #else
    extern Bool ioFlush(void)  asm("?ioFlush@@YA_NXZ");
 #endif
#elif defined(OS_MAC)
   extern Bool ioFlush(void)  asm("__Z7ioFlushv");
#elif defined(OS_LINUX)   
   extern Bool ioFlush(void)  asm("_Z7ioFlushv");
#else
   // Do nothing
#endif

//  https://stackoverflow.com/questions/35837694/how-to-manually-mangle-names-in-visual-c
//  https://gcc.gnu.org/onlinedocs/gcc/Asm-Labels.html
//  ON Windows, for MSVC C++ only, use the __identififer macro to give an alias
//  Using gcc, use int a asm("xxx")

 // ioFlush is not available in Octave
 #ifndef O_INTERFACE
	ioFlush();
 #endif

}



 





I32 GetConsoleWidth() {
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

int IsClass(void* ptr, char* class) { return 0; }

int IsCell(void* ptr) { return mxIsCell(ptr); }
int IsChar(void* ptr)    { 

	if (ptr == NULL) return 0;

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
int IsEmpty(void* ptr) { return mxIsEmpty(ptr); }
int IsStruct(void* ptr)  { return mxIsStruct(ptr) ||  mxIsEmpty(ptr); } // We simply treat an empty mxArray as a struct, though matlab treats it as numeric
int IsNumeric(void* ptr) { return mxIsNumeric(ptr) && !mxIsEmpty(ptr); } // Matlab treates an empty as Numeric
int IsDouble(void* ptr)  { return mxIsDouble(ptr); }
int IsSingle(void* ptr)  { return mxIsSingle(ptr); }
int IsInt32(void* ptr)   { return mxIsInt32(ptr); }
int IsInt16(void* ptr) { return mxIsInt16(ptr); }
int IsInt64(void* ptr) { return mxIsInt64(ptr); }
//int IsInt08(void* ptr) { return mxIsInt8(ptr); }
int IsLogical(void* ptr) { return mxIsLogical(ptr); }

VOID_PTR GetFieldByIdx(VOID_PTR ptr, I32 ind) {

	//https://www.mathworks.com/help/matlab/apiref/mxgetfieldbynumber.html
	if (IsCell(ptr)) {
		return  mxGetCell(ptr, ind);
	}
	if (IsStruct(ptr)) {
		return mxGetFieldByNumber(ptr, 0, ind);
	}
	return NULL;
}

void GetFieldNameByIdx(VOID_PTR strucVar, I32 ind0, char *str, int buflen) {
	I32 numFlds   = mxGetNumberOfFields(strucVar); 
	if (ind0 < numFlds) {
		char* name = mxGetFieldNameByNumber(strucVar, ind0);
		strncpy(str, name, buflen);
		str[buflen - 1] = 0;
	}	else {
		str[0] = 0;
	}
	
}

I32 GetCharArray(void* ptr, char* dst,   int n) {
	dst[0] = 0;
	//Dsr needs to be at least of length (n+1);
	if (ptr == NULL) {
		return 0;
	}
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
		dst[0] = 0;
		return 0L;
	}
	
}

I32 GetCharVecElem(void* ptr, int idx, char* dst, int n) {
	//Dsr needs to be at least of length (n+1);

	int len = 0;

	if (IsCell(ptr) ) {
		int  nelem = mxGetNumberOfElements(ptr);
		if (idx >= nelem) {
			dst[0] = 0;
			return 0;
		}

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
			dst[0] = 0;
			return 0L;
		}

	}
	else if (mxIsChar(ptr)) {
		int            ndims = mxGetNumberOfDimensions(ptr);
		const mwSize*  dims  = mxGetDimensions(ptr);
		if (idx >= dims[0] || ndims > 2) {
			dst[0] = 0;
			return 0;
		}
		else {
			int    cols = dims[1];
			int    rows = dims[0];
			//MATLAB uses 16-bit unsigned integer character encoding for Unicode® characters.
			mwSize   newdims[2] = {1,cols };
			mxArray* newrow     = mxCreateCharArray(2, newdims);
			
			mxChar* newpos = mxGetChars(newrow);
			mxChar* oldpos = mxGetChars(ptr) + idx;
			for (int i = 0; i < cols; i++) {
				newpos[i] = *oldpos;
				oldpos   += rows;
			}			

			char* tmp = mxArrayToString(newrow);
			strncpy(dst, tmp, n); dst[n] = 0;
			len = strlen(dst);
			r_free(tmp);

			mxDestroyArray(newrow);

			return len;
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
		dst[0] = 0;
		return 0L;
	}

}

I32  GetNumberOfFields(const void* structVar) {
	return mxGetNumberOfFields(structVar); 
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

F64    GetScalar(const void * ptr) { return mxGetScalar((mxArray *)ptr); }
void * GetData(const void * ptr) {   return mxGetData((mxArray *)ptr); }
int    GetDim1(const void * ptr) {   return mxGetM((mxArray *)ptr); }
int    GetDim2(const void * ptr) {   return mxGetN((mxArray *)ptr); }
int    GetNumOfDim(const void * ptr) { return mxGetNumberOfDimensions((mxArray *)ptr); }
void   GetDimensions(const void * ptr, int dims[], int ndims) {
	int   N            = min(ndims, GetNumOfDim(ptr) );
	const mwSize *mxdims = mxGetDimensions((mxArray *)ptr);
	for (int i = 0; i < N; i++) 	{
		dims[i] = mxdims[i];
	}	
}

void  * SetDimensions(const void* ptr, int dims[], int ndims) {
	if (!ptr) return NULL;
	mwSize mwdims[20];
	for (int i = 0; i < ndims; i++) mwdims[i] = dims[i];
	mwSize ndimension = ndims;
	int res=mxSetDimensions(ptr, mwdims, ndimension);
	mwSize* p = mxGetDimensions(ptr);
	return ptr;
}

int    GetNumberOfElements(const void * ptr)
{ // Do not work for a vector of opaque objects for some reason 
	if (mxIsChar(ptr)) {
		// Suppose it is a 2d mxCharArray
		const mwSize* dims = mxGetDimensions((mxArray*)ptr);
		return dims[0]; 		
	} else {
		return mxGetNumberOfElements(ptr);
	}	
}

void* CreateNumVar(DATA_TYPE dtype, int* dims, int ndims, VOIDPTR* data_ptr) {

		// Not a struct variable to be created
		mxClassID fieldDataType;
		switch (dtype)	{
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
		if (ndims == 1) {
			// Make it a column vector if it is a 1D. In Matlab, there is no concpet of dimensonless vector
			// But in R, there is.
			mxptr = mxCreateNumericMatrix(dims[0], 1L,      fieldDataType, mxREAL);
		}else if (ndims == 2) {
			mxptr = mxCreateNumericMatrix(dims[0], dims[1], fieldDataType, mxREAL);
		} else if (ndims >= 3)	{
			mwSize DIMS[4] = {dims[0], dims[1],dims[2],dims[3] };
			mxptr = mxCreateNumericArray(ndims, DIMS, fieldDataType, mxREAL);
		}

		if (data_ptr  && mxptr ) {
			*data_ptr = mxGetData(mxptr);
		}

		return mxptr;
}
void * CreateStructVar(FIELD_ITEM *fieldList, int nfields)
{ // www.mathworks.com/help/matlab/apiref/mxsetfieldbynumber.html	
	
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

	mxArray * _restrict out;
	{
		char * fldNames[100];		
		for (int i = 0; i < nfields; i++) { fldNames[i] = fieldList[i].name; }					
		mwSize dims_2d[2] = { 1,1 };
		out = mxCreateStructArray(2, dims_2d, nfields, fldNames);
	}

	for (int i = 0; i < nfields; i++) {

		if (fieldList[i].ptr == NULL) {
			continue;
		}

		mxArray* _restrict newFldPtr;
		if (fieldList[i].type == DATA_STRUCT){			
			newFldPtr = (mxArray*) fieldList[i].ptr;
		} else {
			// the data base is stored in fiedlist.ptr
			newFldPtr = CreateNumVar(fieldList[i].type,  fieldList[i].dims, fieldList[i].ndim, fieldList[i].ptr);
		} 

		mxSetField(out, 0L, fieldList[i].name, newFldPtr);		
	}

	return (void*)out;
}

void ReplaceStructField(VOIDPTR s, char* fname, VOIDPTR newvalue){

	mxArray *fld= mxGetField(s, 0, fname);
	if (fld != NULL) {
		mxDestroyArray(fld);
	}
	mxSetField(s, 0, fname, newvalue);
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


//https://stackoverflow.com/questions/25998442/mex-options-when-compiling-with-visual-studio
//https://undocumentedmatlab.com/articles/mex-ctrl-c-interrupt
// The following 2 function are exported from libut.lib
Bool utIsInterruptPending();
void utSetInterruptPending(Bool);

#ifndef O_INTERFACE
I32  CheckInterrupt()         { return utIsInterruptPending(); }
void ConsumeInterruptSignal() { utSetInterruptPending(_False_); }
#else
I32  CheckInterrupt() { return 0; }
void ConsumeInterruptSignal() { ; }
#endif

#else
const static char achar = 'a';
#endif



#include "abc_000_warning.h"