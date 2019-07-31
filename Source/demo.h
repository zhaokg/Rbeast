#include "abc_000_macro.h"
#include "abc_001_config.h"
#include <windows.h>
#include <synchapi.h>
#include <inttypes.h>
#define __in 
#define _Inout_
#define _Out_
#include <synchapi.h> 
extern WINBASEAPI BOOL WINAPI SleepConditionVariableCS(_Inout_ PCONDITION_VARIABLE ConditionVariable,_Inout_ PCRITICAL_SECTION CriticalSection,_In_ DWORD dwMilliseconds);
extern WINBASEAPI VOID WINAPI WakeConditionVariable(_Inout_ PCONDITION_VARIABLE ConditionVariable);
extern WINBASEAPI VOID WINAPI InitializeConditionVariable(_Out_ PCONDITION_VARIABLE ConditionVariable);
#include "beast_common.h"
#include "abc_datatype.h"
typedef struct style__
{
	float hMargin;
	float rEdit;	
	float hgtButtonBar,vButtonRatio,wButton,sep;
	float vScrollRatio;	
	float rFig[5]; 
	int  xChar,xCapChar,yChar;
	int labelGap;
	int   numRows;
	float fractionLabel;
	float fractionEdit;
	int  widthDialg;
}Style;
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
#else
	mxArray ** plhs,** prhs;
	int nlhs,nrhs;
#endif
	MemPointers *GLOBAL_MEM;
	Options     *GLOBAL_OPTIONS;
	F32PTR      *GLOBAL_Y;
	RESULT      *GLOBAL_RESULT;
} LParam;
typedef struct GlobalStruct
{
	CRITICAL_SECTION cs;
	CONDITION_VARIABLE cv;
	HWND   hwnd;
	HANDLE threadHandle;
	int    timerInterval;
	int    sleepInterval;
	enum stat { RUN,PAUSE,DONE } status;
	byte  quit;
	int * plotData[5][5];
	float yMin,yMax;	
	int N,sample;
	float *y;
	float *t,*tCI,*ct;
	uint16_t * T;
	int32_t * tProb;
	int16_t tKnotNum;
	float *s,*sCI,*curs;
	uint16_t * S;
	int32_t * sProb;
	int16_t sKnotNum;
	float yMinS,yMaxS;
	float yMinT,yMaxT;
	uint32_t nMissing;
	U32PTR  rowsMissing;
	Options opt;
	enum modified {NotAssigned,NoUpdate,NeedUpdate } optStatus;
	int w[5],h[5],x0[5],y0[5];
	int ite;
	int curChainNumber;
	float sN,tN;
} GlobalStruct;
extern Style style;
extern GlobalStruct gData;
extern LParam threadParam;
extern Logger logger;
extern void LoggerInsert(char *newMsg);
extern DWORD WINAPI beastST_demo(__in LPVOID lpParameter);
extern DWORD WINAPI beastTrend_demo(__in LPVOID lpParameter);
extern void InitGlobalData();
extern void AllocatePlotData();
extern void ResizeControls(int x,int y,int N);
extern void CreateGDIObject(HWND hwnd,int N);
extern void DeleteGDIObject(int N);
extern void GeneratePlotData();
extern void DrawPlots(HDC hdc);
extern void ResizeDialogControls(HWND hwnd);
extern void EnableButtons(HWND *hButton,HWND *hDiag);
extern DWORD WINAPI beastST_demo(__in LPVOID lpParameter);
extern void InitGlobalData_ST();
extern void AllocatePlotData_ST();
extern void GeneratePlotData_ST();
extern void DrawPlots_ST(HDC hdc);
extern void ResizeDialogControls_ST(HWND hwnd);
extern void EnableButtons_ST(HWND *hButton,HWND *hDiag);
extern BYTE ANDmaskIcon[];
extern BYTE XORmaskIcon[];
extern HICON   hIcon;
extern HWND    hButton[4],hScroll,hEdit,hStatic[2],hDiag;
extern HBITMAP hBitmap[5],hBufferBitmap[5]; 
extern HDC     memDC[5],bufferDC[5];        
extern HBRUSH  blackBrush,grayBrush;
extern HPEN	   greenPen,bluePen,yellowPen,redPen;
#define OEM_FIXED_FONT      10
#define ANSI_FIXED_FONT     11
#define ANSI_VAR_FONT       12
#define SYSTEM_FONT         13
#define DEVICE_DEFAULT_FONT 14
#define SYSTEM_FIXED_FONT   16
#define DEFAULT_GUI_FONT    17
static const int  ID_TIMER=1;
#define IDC_MAIN_EDIT	101
