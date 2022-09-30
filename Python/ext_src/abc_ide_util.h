#pragma once
#include "abc_001_config.h"
#include "abc_datatype.h"

extern I08 IDE_USER_INTERRUPT;

typedef enum   IO_TYPE { MEM_IO,	DISK_IO } IO_TYPE;
typedef struct FIELD_ITEM {
	char      name[64 - 1];
	DATA_TYPE type;
	int       ndim;
	int       dims[5];
	void **   ptr;
	int       extra; //added for extension to multivariate ts for mrbeast
} FIELD_ITEM;

VOID_PTR GetFieldByIdx(VOID_PTR strucVar, I32 ind);
void * CreateNumVar(DATA_TYPE dtype, int* dims, int ndims, VOIDPTR* data_ptr);
void   ReplaceStructField(VOIDPTR s, char* fname, VOIDPTR newvalue);
void * CreateStructVar(FIELD_ITEM *fieldList, int nfields);
void   DestoryStructVar(VOID_PTR strutVar);
void   RemoveField(FIELD_ITEM *fieldList, int nfields, char * fieldName);

void AddStringAttribute(VOID_PTR listVar, const char* field, const char* value);
void AddIntegerAttribute(VOID_PTR listVar, const char* field, I32 value);
void RemoveAttribute(VOID_PTR listVar, const char* field);

extern  I32   GetConsoleWidth();
extern  void  printProgress(F32 pct, I32 width, char * buf, I32 firstTimeRun);
extern  void  printProgress2(F32 pct, F64 time, I32 width, char * buf, I32 firstTimeRun);

I32 GetCharArray(void *ptr, char * dst, int n);
I32 GetCharVecElem(void* ptr, int idx, char* dst, int n);

void *GetField123(const void* structVar, char* fname, int nPartial);
void *GetField(const void * structVar, char *fname);

void* GetField123Check(const void* structVar, char* fname, int nPartial);
void* GetFieldCheck(const void* structVar, char* fname);

F64   GetScalar(const void * ptr);
F64   GetNumericElement(const void* Y, I32 idx);
void * GetData(const void * ptr);
int  GetDataType(VOID_PTR Y);
int  GetDim1(const void * ptr);
int  GetDim2(const void * ptr);
int  GetNumOfDim(const void * ptr);
void GetDimensions(const void * ptr, int dims[], int ndims);
int  GetNumberOfElements(const void* ptr);
I32  GetNumberOfFields(const void* structVar);
int IsCell(void* ptr);
int IsChar(void* ptr);
int IsEmpty(void* ptr);
int IsStruct(void* ptr);
int IsNumeric(void* ptr);
int IsDouble(void* ptr);
int IsSingle(void* ptr);
int IsInt32(void* ptr);
int IsInt16(void* ptr);
int IsInt64(void* ptr);
int IsLogical(void* ptr);
int HaveEqualDimesions(const void* p1, const void* p2);
int CopyNumericArrToF32Arr(F32PTR outmem, VOID_PTR infield, int N);
int CopyNumericArrToI32Arr(I32PTR outmem, VOID_PTR infield, int N);

extern I32  CheckInterrupt();
extern void ConsumeInterruptSignal();

static INLINE int IsRinterface() { return R_INTERFACE;}
static INLINE int IsMinterface() { return M_INTERFACE; }
static INLINE int IsPinterface() { return P_INTERFACE; }

/**********************************************/
// Convert the input data into a contiguours Fortran aarray.
// It will genrate a new referrence that must be explicilty
// decref-ed when program exits.

// For R and Matlab, nothing will be done
void* CvtToPyArray_NewRef(VOIDPTR Y);
/**********************************************/

#if M_INTERFACE==1
	#define  r_printf(...)  mexPrintf(__VA_ARGS__)
	//#define  r_error(x)    mexErrMsgTxt(x)  //mexErrMsgTxt dones' take a va_args input.
	#define  r_error(...)   mexPrintf(__VA_ARGS__)
    #define  r_warning(...) mexPrintf(__VA_ARGS__)
	#define  r_malloc(x)    mxMalloc(x) 
	#define  r_free(x)      mxFree(x)
    #define  IDE_NULL       mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxREAL)
#elif R_INTERFACE==1
	//stackoverflow.com/questions/37206118/va-args-not-swallowing-comma-when-zero-args-under-c99
	//#define  mexPrintf(output, ...) Rprintf(output, ##__VA_ARGS__) //"##" is used to swallow the preceding comma if it is empty!
	#define  r_printf(...)   Rprintf(__VA_ARGS__)
	#define  r_error(...)    error(__VA_ARGS__)
	//#define  r_warning(...)  warning(__VA_ARGS__)
    #define  r_warning(...)  Rf_warning(__VA_ARGS__)
	#define  r_malloc(x)     Calloc(x, char)  //from header file R_ext\RS.h
	#define  r_free(x)       Free(x) 
    #define   IDE_NULL          R_NilValue
#elif P_INTERFACE==1
	//stackoverflow.com/questions/37206118/va-args-not-swallowing-comma-when-zero-args-under-c99
	//#define  mexPrintf(output, ...) Rprintf(output, ##__VA_ARGS__) //"##" is used to swallow the preceding comma if it is empty!
	#define  r_printf(...)   printf(__VA_ARGS__)
    #define  r_printf(...)   PySys_WriteStdout(__VA_ARGS__)
	#define  r_error(...)    printf(__VA_ARGS__)
	//#define  r_warning(...)  warning(__VA_ARGS__)
    #define  r_warning(...)  printf(__VA_ARGS__)
	#define  r_malloc(x)     PyMem_RawMalloc(x, char)  //from header file R_ext\RS.h
	#define  r_free(x)       PyMem_RawFree(x) 
    #define   IDE_NULL       Py_None
#endif

#if R_INTERFACE==1

extern  SEXP	getListElement(SEXP list, const char *str);
extern  SEXP    getListElement_CaseIn(SEXP list, const char *str);


// https: //github.com/ryanoasis/public-bash-scripts/blob/master/unix-color-codes.sh
// https: //stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux/5947779

/*

Color_Off = "\[\033[0m\]"       # Text Reset

# Regular Colors
Black = "\[\033[0;30m\]"        # Black
Red = "\[\033[0;31m\]"          # Red
Green = "\[\033[0;32m\]"        # Green
Yellow = "\[\033[0;33m\]"       # Yellow
Blue = "\[\033[0;34m\]"         # Blue
Purple = "\[\033[0;35m\]"       # Purple
Cyan = "\[\033[0;36m\]"         # Cyan
White = "\[\033[0;37m\]"        # White

# Bold
BBlack = "\[\033[1;30m\]"       # Black
BRed = "\[\033[1;31m\]"         # Red
BGreen = "\[\033[1;32m\]"       # Green
BYellow = "\[\033[1;33m\]"      # Yellow
BBlue = "\[\033[1;34m\]"        # Blue
BPurple = "\[\033[1;35m\]"      # Purple
BCyan = "\[\033[1;36m\]"        # Cyan
BWhite = "\[\033[1;37m\]"       # White

# Underline
UBlack = "\[\033[4;30m\]"       # Black
URed = "\[\033[4;31m\]"         # Red
UGreen = "\[\033[4;32m\]"       # Green
UYellow = "\[\033[4;33m\]"      # Yellow
UBlue = "\[\033[4;34m\]"        # Blue
UPurple = "\[\033[4;35m\]"      # Purple
UCyan = "\[\033[4;36m\]"        # Cyan
UWhite = "\[\033[4;37m\]"       # White

# Background
On_Black = "\[\033[40m\]"       # Black
On_Red = "\[\033[41m\]"         # Red
On_Green = "\[\033[42m\]"       # Green
On_Yellow = "\[\033[43m\]"      # Yellow
On_Blue = "\[\033[44m\]"        # Blue
On_Purple = "\[\033[45m\]"      # Purple
On_Cyan = "\[\033[46m\]"        # Cyan
On_White = "\[\033[47m\]"       # White

# High Intensty
IBlack = "\[\033[0;90m\]"       # Black
IRed = "\[\033[0;91m\]"         # Red
IGreen = "\[\033[0;92m\]"       # Green
IYellow = "\[\033[0;93m\]"      # Yellow
IBlue = "\[\033[0;94m\]"        # Blue
IPurple = "\[\033[0;95m\]"      # Purple
ICyan = "\[\033[0;96m\]"        # Cyan
IWhite = "\[\033[0;97m\]"       # White

# Bold High Intensty
BIBlack = "\[\033[1;90m\]"      # Black
BIRed = "\[\033[1;91m\]"        # Red
BIGreen = "\[\033[1;92m\]"      # Green
BIYellow = "\[\033[1;93m\]"     # Yellow
BIBlue = "\[\033[1;94m\]"       # Blue
BIPurple = "\[\033[1;95m\]"     # Purple
BICyan = "\[\033[1;96m\]"       # Cyan
BIWhite = "\[\033[1;97m\]"      # White

# High Intensty backgrounds
On_IBlack = "\[\033[0;100m\]"   # Black
On_IRed = "\[\033[0;101m\]"     # Red
On_IGreen = "\[\033[0;102m\]"   # Green
On_IYellow = "\[\033[0;103m\]"  # Yellow
On_IBlue = "\[\033[0;104m\]"    # Blue
On_IPurple = "\[\033[10;95m\]"  # Purple
On_ICyan = "\[\033[0;106m\]"    # Cyan
On_IWhite = "\[\033[0;107m\]"   # White

BLACK = "\[\e[00;30m\]"
DARY_GRAY = "\[\e[01;30m\]"
RED = "\[\e[00;31m\]"
BRIGHT_RED = "\[\e[01;31m\]"
GREEN = "\[\e[00;32m\]"
BRIGHT_GREEN = "\[\e[01;32m\]"
BROWN = "\[\e[00;33m\]"
YELLOW = "\[\e[01;33m\]"
BLUE = "\[\e[00;34m\]"
BRIGHT_BLUE = "\[\e[01;34m\]"
PURPLE = "\[\e[00;35m\]"
LIGHT_PURPLE = "\[\e[01;35m\]"
CYAN = "\[\e[00;36m\]"
BRIGHT_CYAN = "\[\e[01;36m\]"
LIGHT_GRAY = "\[\e[00;37m\]"
WHITE = "\[\e[01;37m\]"
ENDCOLOR = "\e[m"
*/

#elif M_INTERFACE==1
	#define PROTECT(XXXX)   XXXX
	#define UNPROTECT(XXXX) XXXX
#elif P_INTERFACE == 1 
	#define PROTECT(XXXX)   XXXX
	#define UNPROTECT(XXXX) XXXX
	extern PyTypeObject BarObject_Type;  //tentative definition
	extern PyObject    *currentModule;
	extern PyObject* classOutout; 
	extern PyObject* setClassObjects(PyObject* self, PyObject* args);
#endif