#include "abc_000_macro.h"
#include "abc_001_config.h"
#include "abc_mem.h"
#include "abc_datatype.h"
#include <inttypes.h>
#include <windows.h>
#include <synchapi.h>


#define __in 
#define _Inout_
#define _Out_
#include <synchapi.h> //SleepConditionVariableCS

//extern WINBASEAPI BOOL WINAPI SleepConditionVariableCS(    _Inout_ PCONDITION_VARIABLE ConditionVariable, _Inout_ PCRITICAL_SECTION CriticalSection, _In_ DWORD dwMilliseconds);
//extern WINBASEAPI VOID WINAPI WakeConditionVariable(       _Inout_ PCONDITION_VARIABLE ConditionVariable);
//extern WINBASEAPI VOID WINAPI InitializeConditionVariable( _Out_ PCONDITION_VARIABLE ConditionVariable);



typedef struct style__
{
	F32 hMargin;
	F32 rEdit;	
	F32 hgtButtonBar, vButtonRatio, wButton,sep;
	F32 vScrollRatio;	
	F32 rFig[5]; //the first two for beastTrend
	int  xChar, xCapChar, yChar;
	int labelGap;
	//////////////////////
	int   numRows;
	F32 fractionLabel;
	F32 fractionEdit;
	int widthDialg;
} Style;

typedef struct log_
{
	char *str;
	int  LEN;
	int  len;
} Logger;

typedef struct lpParam
{
#if R_INTERFACE==1
	SEXP pY;
	SEXP pOpt;
#elif M_INTERFACE==1
	mxArray ** plhs, ** prhs;
	int nlhs, nrhs;
#elif P_INTERFACE==1

#endif
	MemPointers *GLOBAL_MEM;
	VOID_PTR     GLOBAL_OPTIONS;
	F32PTR      *GLOBAL_Y;
	VOID_PTR     GLOBAL_RESULT;
} LParam;

typedef struct GlobalStruct
{
	CRITICAL_SECTION   cs;
	CONDITION_VARIABLE cv;
	HWND   hwnd;
	HANDLE threadHandle;
	int    timerInterval;
	int    sleepInterval;
	enum __status { RUN, PAUSE, DONE } status;
	U08  quit;
	int * plotData[5][5];
	F32 yMin, yMax;
	int N,    sample;

	//----------------------
	F32PTR  y;
	F32PTR  t, tCI, curt;
	I32PTR  T;  //must be consisetnw with TKNOT_PTR
	I32PTR  tProb;
	I16     tKnotNum;
	//------------------
	F32PTR  s, sCI, curs;
	I32PTR  S;  //must be consisetnw with TKNOT_PTR
	I32PTR  sProb;
	I16     sKnotNum;
	F32     yMinS, yMaxS;
	F32     yMinT, yMaxT;
	U32     nMissing;
	U32PTR  rowsMissing;
	//-------------
	enum modified { NotAssigned, NoUpdate, NeedUpdate } optStatus;
	//---
	int w[5], h[5], x0[5], y0[5];
	int ite;

	//-------------------
	int curChainNumber;
	F32 sN, tN;
} GlobalStruct;

extern GlobalStruct gData;
extern Style        style;
extern LParam       threadParam;
extern Logger       logger;

extern void LoggerInsert(char *newMsg);
//---------------------------------------------
extern void ResizeControls(int x, int y, int N);
extern void CreateGDIObject(HWND hwnd, int N);
extern void DeleteGDIObject(int N);

extern void EnableButtons(HWND * hButton, HWND *hDiag);
extern void ResizeDialogControls(HWND hwnd);


extern BYTE ANDmaskIcon[];
extern BYTE XORmaskIcon[];

extern HICON   hIcon;
extern HWND    hButton[4], hScroll, hEdit, hStatic[2], hDiag;
extern HBITMAP hBitmap[5], hBufferBitmap[5]; //only the first two used by beastTrend
extern HDC     memDC[5], bufferDC[5];        //only the first two used by beastTrend  
extern HBRUSH  blackBrush, grayBrush;
extern HPEN	   greenPen, bluePen, yellowPen, redPen;

#define OEM_FIXED_FONT      10
#define ANSI_FIXED_FONT     11
#define ANSI_VAR_FONT       12
#define SYSTEM_FONT         13
#define DEVICE_DEFAULT_FONT 14
#define SYSTEM_FIXED_FONT   16
#define DEFAULT_GUI_FONT    17
//hfDefault = GetStockObject(DEFAULT_GUI_FONT);

static const int  ID_TIMER = 1;
#define IDC_MAIN_EDIT	101
