#include "abc_000_macro.h"
#include "abc_000_warning.h"
 

#if defined(OS_WIN64) 

#include <stdio.h>   // printf
#include <io.h>      // _open_osfhandle, _dup2
#include <math.h>
#include "abc_win32_demo.h"
#include "beastv2_header.h"
#include "globalvars.h"
#include "abc_ide_util.h"
#include "beastv2_io.h"

/*
DWORD WINAPI MY_THREAD(__in LPVOID lpParameter)
{
char s[100];//wsprintf(s, "Thread inside %d \n", GetCurrentThreadId());
for (int i = 0; i < 50; i++)
{
Sleep(gData.timerInterVal);
EnterCriticalSection(&gData.cs);
for (int i = 0; i < (gData.N * 2); i += 2)
{
gData.data[i] = (i / 2.) / 100.*style.wfig1;
gData.data[i + 1] = (int)((F32)rand() / RAND_MAX * style.hfig1);

}
while (gData.status == PAUSE & gData.quit == 0)
{
SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);
}
if (gData.quit)
{
LeaveCriticalSection(&gData.cs);
return;
}
LeaveCriticalSection(&gData.cs);
}
PostMessage(gData.hwnd, WM_USER + 1, 0, 0);
return 0;
}
*/
static LRESULT CALLBACK DialogProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);
static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);
static HICON CreateMyIcon(void) {
	// Yang icon AND bitmask 
	static BYTE ANDmaskIcon[] = { 0xFF, 0xFF, 0xFF, 0xFF,   // line 1 
			0xFF, 0xFF, 0xC3, 0xFF,   // line 2 
			0xFF, 0xFF, 0x00, 0xFF,   // line 3 
			0xFF, 0xFE, 0x00, 0x7F,   // line 4 

			0xFF, 0xFC, 0x00, 0x1F,   // line 5 
			0xFF, 0xF8, 0x00, 0x0F,   // line 6 
			0xFF, 0xF8, 0x00, 0x0F,   // line 7 
			0xFF, 0xF0, 0x00, 0x07,   // line 8 

			0xFF, 0xF0, 0x00, 0x03,   // line 9 
			0xFF, 0xE0, 0x00, 0x03,   // line 10 
			0xFF, 0xE0, 0x00, 0x01,   // line 11 
			0xFF, 0xE0, 0x00, 0x01,   // line 12 

			0xFF, 0xF0, 0x00, 0x01,   // line 13 
			0xFF, 0xF0, 0x00, 0x00,   // line 14 
			0xFF, 0xF8, 0x00, 0x00,   // line 15 
			0xFF, 0xFC, 0x00, 0x00,   // line 16 

			0xFF, 0xFF, 0x00, 0x00,   // line 17 
			0xFF, 0xFF, 0x80, 0x00,   // line 18 
			0xFF, 0xFF, 0xE0, 0x00,   // line 19 
			0xFF, 0xFF, 0xE0, 0x01,   // line 20 

			0xFF, 0xFF, 0xF0, 0x01,   // line 21 
			0xFF, 0xFF, 0xF0, 0x01,   // line 22 
			0xFF, 0xFF, 0xF0, 0x03,   // line 23 
			0xFF, 0xFF, 0xE0, 0x03,   // line 24 

			0xFF, 0xFF, 0xE0, 0x07,   // line 25 
			0xFF, 0xFF, 0xC0, 0x0F,   // line 26 
			0xFF, 0xFF, 0xC0, 0x0F,   // line 27 
			0xFF, 0xFF, 0x80, 0x1F,   // line 28 

			0xFF, 0xFF, 0x00, 0x7F,   // line 29 
			0xFF, 0xFC, 0x00, 0xFF,   // line 30 
			0xFF, 0xF8, 0x03, 0xFF,   // line 31 
			0xFF, 0xFC, 0x3F, 0xFF };  // line 32 
	// Yang icon XOR bitmask 
	static BYTE XORmaskIcon[] = { 0x00, 0x00, 0x00, 0x00,   // line 1 
				0x00, 0x00, 0x00, 0x00,   // line 2 
				0x00, 0x00, 0x00, 0x00,   // line 3 
				0x00, 0x00, 0x00, 0x00,   // line 4 

				0x00, 0x00, 0x00, 0x00,   // line 5 
				0x00, 0x00, 0x00, 0x00,   // line 6 
				0x00, 0x00, 0x00, 0x00,   // line 7 
				0x00, 0x00, 0x38, 0x00,   // line 8 

				0x00, 0x00, 0x7C, 0x00,   // line 9 
				0x00, 0x00, 0x7C, 0x00,   // line 10 
				0x00, 0x00, 0x7C, 0x00,   // line 11 
				0x00, 0x00, 0x38, 0x00,   // line 12 

				0x00, 0x00, 0x00, 0x00,   // line 13 
				0x00, 0x00, 0x00, 0x00,   // line 14 
				0x00, 0x00, 0x00, 0x00,   // line 15 
				0x00, 0x00, 0x00, 0x00,   // line 16 

				0x00, 0x00, 0x00, 0x00,   // line 17 
				0x00, 0x00, 0x00, 0x00,   // line 18 
				0x00, 0x00, 0x00, 0x00,   // line 19 
				0x00, 0x00, 0x00, 0x00,   // line 20 

				0x00, 0x00, 0x00, 0x00,   // line 21 
				0x00, 0x00, 0x00, 0x00,   // line 22 
				0x00, 0x00, 0x00, 0x00,   // line 23 
				0x00, 0x00, 0x00, 0x00,   // line 24 

				0x00, 0x00, 0x00, 0x00,   // line 25 
				0x00, 0x00, 0x00, 0x00,   // line 26 
				0x00, 0x00, 0x00, 0x00,   // line 27 
				0x00, 0x00, 0x00, 0x00,   // line 28 

				0x00, 0x00, 0x00, 0x00,   // line 29 
				0x00, 0x00, 0x00, 0x00,   // line 30 
				0x00, 0x00, 0x00, 0x00,   // line 31 
				0x00, 0x00, 0x00, 0x00 };  // line 32 
	//char szFileName[101];
	//GetModuleFileName(hInstance, szFileName, 100); 
	HINSTANCE   hInstance = GetModuleHandle(NULL);
	HICON hIcon = CreateIcon(
		hInstance,    // application instance  
		32,              // icon width 
		32,              // icon height 
		1,               // number of XOR planes 
		1,               // number of bits per pixel 
		ANDmaskIcon,     // AND bitmask  
		XORmaskIcon);    // XOR bitmask 
	
	return hIcon;
}
static void  ClassRegisterWnd(char * wndClassName, HICON hIcon)
{
	HINSTANCE   hInstance = GetModuleHandle(NULL);

	WNDCLASSEX	wc;
	wc.cbSize	= sizeof(WNDCLASSEX);
	wc.style	= CS_DBLCLKS | CS_HREDRAW | CS_DROPSHADOW;
	//https:// stackoverflow.com/questions/2285860/change-win32-window-style 
	//SetWindowLongPtr SetWindowLong(hWnd, GWL_STYLE, newStyle); ShowWindow(hWnd, SW_SHOW);
	wc.lpfnWndProc = WndProc;
	wc.cbClsExtra	= 0;
	wc.cbWndExtra	= 0;
	wc.hInstance	= hInstance;
	wc.hIcon			= hIcon; //LoadIcon(NULL, IDI_HAND);
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground	= (HBRUSH)(COLOR_WINDOW);
	wc.lpszMenuName		= NULL;
	wc.lpszClassName	= (LPCSTR)wndClassName;
	wc.hIconSm			= hIcon;// LoadIcon(NULL, IDI_APPLICATION);

	if (!RegisterClassEx(&wc)) {
		MessageBox(NULL, "Window Registration Failed!", "Error!", MB_ICONEXCLAMATION | MB_OK);
		return;
	}
}
static void  ClassRegisterDialog(char* wndClassName)
{
 
	WNDCLASSEX wc = { 0 };
	wc.cbSize	 = sizeof(WNDCLASSEX);
	wc.style	 = CS_DBLCLKS | CS_HREDRAW | CS_DROPSHADOW;
	//https:// stackoverflow.com/questions/2285860/change-win32-window-style 
	//SetWindowLongPtr SetWindowLong(hWnd, GWL_STYLE, newStyle); ShowWindow(hWnd, SW_SHOW);
	wc.lpfnWndProc = DialogProc;
	wc.cbClsExtra  = 0;
	wc.cbWndExtra  = 0;
	wc.hInstance  = GetModuleHandle(NULL);
	wc.hIcon		= NULL;// hIcon; //LoadIcon(NULL, IDI_HAND);
	wc.hCursor		= NULL;//LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)(COLOR_WINDOW);
	wc.lpszMenuName  = NULL;
	wc.lpszClassName = (LPCSTR) wndClassName;
	wc.hIconSm = NULL;//hIcon;// LoadIcon(NULL, IDI_APPLICATION);
	if (!RegisterClassEx(&wc))
	{
		MessageBox(NULL, "DialogBox Class Registration Failed!", "Error!", MB_ICONEXCLAMATION | MB_OK);
	}
}

 
/***********************************/ 
// Defined in beastv2_gui_plot.c
/***********************************/
extern void BEAST2_InitGlobalData(void);
extern void BEAST2_AllocatePlotData(void);
extern void BEAST2_GeneratePlotData(void);
extern void BEAST2_DrawPlots(HDC hdc);


DWORD WINAPI beast_thread(__in LPVOID dummy_GLOBAL_OPTIONS) {
	#if R_INTERFACE==1
	//mexPrintf is not thread-safe.
	BEAST2_print_options(GLOBAL_OPTIONS);
	#endif

	VOIDPTR ANS=NULL;
	PROTECT(ANS=BEAST2_Output_AllocMEM(GLOBAL_OPTIONS));
	beast2_main_corev4_gui();
	UNPROTECT(1);
	DestoryStructVar(ANS);

	return 0;
}
void  DllExport BEAST2_WinMain(VOID_PTR  option)
{
	logger.LEN = 2000;
	logger.len = 0;
	logger.str = malloc(logger.LEN + 1);
	

	BEAST2_InitGlobalData();
	{
		DWORD  thread_id;		
		gData.threadHandle = CreateThread(0, 0, beast_thread, option, 0, &thread_id);
		//gData.threadHandle = (HANDLE)_beginthreadex(0, 0, &mythread, 0, 0, 0);
	}
	//First start the thread so gData.N is set in beast2_main_core4. If gData.N is not set in time,
	// CreateWindowEx will conditionally wait for it
	
	HICON       hIcon				  = CreateMyIcon();	
	ClassRegisterWnd( "MyWndClassBeast2", hIcon);
	ClassRegisterDialog("MyDiaglogClass3434");

	HINSTANCE   hInstance = GetModuleHandle(NULL);
	HWND        hwnd      = CreateWindowEx(
								WS_EX_CLIENTEDGE  /*WS_EX_TOOLWINDOW|WS_EX_STATICEDGE*/,
								"MyWndClassBeast2",
								"Bayesian time series decomposition and changepoint detection",
								(WS_OVERLAPPEDWINDOW & ~WS_MAXIMIZEBOX),  // | WS_THICKFRAME
								CW_USEDEFAULT, CW_USEDEFAULT, 520, 740,
								NULL, NULL, hInstance, NULL
							);
	if (!hwnd) {
		MessageBox(NULL, "Window Creation Failed!", "Error!", MB_ICONEXCLAMATION | MB_OK);
		return;
	}
	else {
		// The system printf is thread-safe, but when you compile a mexFunction, but mex.h defines printf as 
		// a macro which is just mexPrintf.  So you get mexPrintf instead, which then dies if used inside the
		//  parallel region.
		// https://www.mathworks.com/matlabcentral/answers/513664-error-message-management-cpp-712-anonymous-namespace-when-compiling-mex-with-openmp
		
		// r_printf("\nWindows created succesfully!\n");
	}
	ShowWindow(hwnd, SW_SHOWDEFAULT);
	UpdateWindow(hwnd);
	SetWindowPos(hwnd, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
	//------------------------------------------------
	
	EnterCriticalSection(&gData.cs);
	gData.hwnd = hwnd;
	LeaveCriticalSection(&gData.cs);
	WakeConditionVariable(&gData.cv); //the beast2_coreve4 threat is waiting on gData.hWnd to be set.

	//WNDPROC fWndProc = (WNDPROC)GetWindowLong(hwnd, GWL_WNDPROC);
	//char str[200];
	//wsprintf(str, "0x%0x  %x", fWndProc, WndProc);
	//---------------------------------------------
	//buf = (char*)GlobalAlloc(GPTR, len + 1);
	//GetWindowTextLength()
	//GetDlgItemText(hwnd, IDC_TEXT, buf, len + 1);
	//GlobalFree((HANDLE)buf);
	//  int index = SendDlgItemMessage(hwnd, IDC_LIST, LB_ADDSTRING, 0, (LPARAM)"Hi there!");
	//  SendDlgItemMessage(hwnd, IDC_LIST, LB_SETITEMDATA, (WPARAM)index, (LPARAM)nTimes);
	//Edits SetDlgItemText(hwnd, IDC_TEXT, "This is a string");

	// HWND    hList  = GetDlgItem(hwnd, IDC_LIST);
	// int     count  = SendMessage(hList, LB_GETSELCOUNT, 0, 0);
	// int     data   = SendMessage(hList, LB_GETITEMDATA, (WPARAM)index, 0);
	// HBITMAP hData  = (HBITMAP)SendMessage(hList, LB_GETITEMDATA, (WPARAM)index, 0);
		

	/*
	{
		SetupConsole();
		HWND htmp = GetConsoleHwnd();
		RECT rc;
		GetWindowRect(hwnd, &rc);
		MoveWindow(htmp, rc.left, rc.bottom, rc.right - rc.left, 500, TRUE);
		//printf("%d, %d\n", GetWindowLongPtr(htmp, GWL_STYLE)& WS_CHILD, hwnd);
		//SetParent(htmp, hwnd);
	}*/

	MSG msg;
	while (GetMessage(&msg, NULL, 0, 0) > 0) {
		if (!IsDialogMessage(hDiag, &msg)) 		{
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}

	WaitForSingleObject(gData.threadHandle, INFINITE);
	CloseHandle(gData.threadHandle);
	gData.threadHandle = NULL;

	DeleteCriticalSection(&gData.cs);
	DestroyIcon(hIcon);
	UnregisterClass("MyWndClassBeast2",  hInstance);
	UnregisterClass("MyDiaglogClass3434", hInstance);

	free(logger.str);
	return ;
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
	switch (msg) {
	case WM_CREATE:
	{
		{
			EnterCriticalSection(&gData.cs);
			while (gData.N <= 0)
				SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);
			LeaveCriticalSection(&gData.cs);

			EnterCriticalSection(&gData.cs);
				BEAST2_AllocatePlotData();	 
			LeaveCriticalSection(&gData.cs);
			WakeConditionVariable(&gData.cv);

		}
		
		if (!SetTimer(hwnd, ID_TIMER, gData.timerInterval, NULL))  
			MessageBox(hwnd, "Could not SetTimer()!", "Error", MB_OK | MB_ICONEXCLAMATION);

		{
			hEdit = CreateWindowEx(WS_EX_CLIENTEDGE, "EDIT", "",
						WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | ES_MULTILINE | ES_AUTOVSCROLL | ES_AUTOHSCROLL,
						0, 0, 100, 100, hwnd, (HMENU)IDC_MAIN_EDIT, GetModuleHandle(NULL), NULL);
			if (hEdit == NULL) 
				MessageBox(hwnd, "Could not create edit box.\0\0\0", "Error\0\0\0", MB_OK | MB_ICONERROR);

			HFONT hfDefault = GetStockObject(SYSTEM_FONT);
			SendMessage(  hEdit, WM_SETFONT, (WPARAM) hfDefault, MAKELPARAM(FALSE, 0));
			SetWindowText(hEdit, "");
			SendDlgItemMessage(hwnd, IDC_MAIN_EDIT, EM_SETREADONLY, (WPARAM)TRUE, 0);
		}
		HFONT hfontDefault = GetStockObject(SYSTEM_FONT);
		{		
			char* caption[] = { "Run\0\0\0", "Pause\0\0\0", "Setting\0\0\0", "Exit\0\0\0" };			
			for (int i = 0; i < 4; i++)
			{
				hButton[i] = CreateWindowEx(WS_EX_CLIENTEDGE, "BUTTON", caption[i], WS_CHILD | WS_VISIBLE /* | BS_PUSHBUTTON   | BS_PUSHLIKE | BS_OWNERDRAW*/,
									        10 + 100 * i, 20, 30 * 9 / 4, 30,   hwnd, (HMENU)i, GetModuleHandle(NULL), NULL);
				if (hButton[i] == NULL)
					MessageBox(hwnd, "Could not create edit box.\0\0\0", "Error", MB_OK | MB_ICONERROR);
				SendMessage(hButton[i], WM_SETFONT, (WPARAM)hfontDefault, MAKELPARAM(FALSE, 0));				
			}
			SendMessage(hButton[0], BM_SETCHECK, (WPARAM)BST_UNCHECKED, MAKELPARAM(FALSE, 0));
			EnableButtons(hButton, NULL);
		}
		{
			SCROLLINFO si;
			hScroll = CreateWindowEx(
				0,                      // no extended styles 
				"SCROLLBAR",           // scroll bar control class 
				(PTSTR)NULL,           // no window text 
				WS_CHILD | WS_VISIBLE   // window styles  
				| SBS_HORZ,         // horizontal scroll bar style 
				400,              // horizontal position 
				50, // vertical position 
				100,             // width of the scroll bar 
				10,               // height of the scroll bar
				hwnd,             // handle to main window 
				(HMENU)10,           // no menu 
				GetModuleHandle(NULL),                // instance owning this window 
				(PVOID)NULL            // pointer not needed 
			);
			si.cbSize = sizeof(si);
			si.nMin = 0;
			si.nMax = 999;
			si.nPos = 500;
			si.nPage = 50;
			si.fMask = SIF_RANGE | SIF_POS | SIF_PAGE;
			SetScrollInfo(hScroll, SB_CTL, &si, TRUE);
			//950=nMax-nPage+1;
			F32 x = (max(950 - si.nPos, 0));
			x = x * x * x * x;
			gData.sleepInterval = 1e-9 * x;
		}
		{
			hStatic[0] = CreateWindowEx(0, "STATIC", "SLOW", WS_CHILD | WS_VISIBLE | SS_CENTER,
											400, 100, 100, 100, hwnd, (HMENU)20, GetModuleHandle(NULL), NULL);
			SendMessage(hStatic[0], WM_SETFONT, (WPARAM)hfontDefault, MAKELPARAM(FALSE, 0));
			SetWindowText(hStatic[0], "Slow");
			//SendDlgItemMessage(hwnd, IDC_MAIN_EDIT, EM_SETREADONLY, (WPARAM)TRUE, 0);

			hStatic[1] = CreateWindowEx(0, "STATIC", "FAST", WS_CHILD | WS_VISIBLE | WS_TABSTOP,
										500, 100, 100, 100, hwnd, (HMENU)21, GetModuleHandle(NULL), NULL);
			SendMessage(hStatic[1], WM_SETFONT, (WPARAM)hfontDefault, MAKELPARAM(FALSE, 0));
			SetWindowText(hStatic[1], "Fast");
		}
		{
			TEXTMETRIC tm;
			HDC hdc;

			hdc = GetDC(hStatic[0]);
			GetTextMetrics(hdc, &tm);
			style.xChar = tm.tmAveCharWidth;
			style.xCapChar = (tm.tmPitchAndFamily & 1 ? 3 : 2) * style.xChar / 2;
			style.yChar = tm.tmHeight + tm.tmExternalLeading;
			ReleaseDC(hStatic[0], hdc);

		}
		{
			EnterCriticalSection(&gData.cs);
				CreateGDIObject(hwnd, 5);
			LeaveCriticalSection(&gData.cs);
		}
		{
		
			hDiag = CreateWindowEx(
				WS_EX_CLIENTEDGE,
				"MyDiaglogClass3434",
				"Settings",
				WS_OVERLAPPEDWINDOW | (WS_CAPTION & ~WS_VISIBLE),
				CW_USEDEFAULT, CW_USEDEFAULT, 320, 240,
				hwnd, (HMENU)NULL, GetModuleHandle(NULL), NULL);
			if (hDiag == 0)
				MessageBox(NULL, "Window Registration Failed!", "Error!", MB_ICONEXCLAMATION | MB_OK);
			ShowWindow(hDiag, SW_HIDE);
			//UpdateWindow(hDiag);
		}



	}
	break;
	case WM_HSCROLL:
	{

					   HWND hScroll;//= (HWND)lParam;
					   WORD action; //= LOWORD(wParam);
					   hScroll = (HWND)lParam;
					   action = LOWORD(wParam);

					   int pos;
					   pos = GetScrollPos(hScroll, SB_CTL);
					   //printf("%d\n", pos);
					   if (action == SB_THUMBPOSITION || action == SB_THUMBTRACK) {
						   pos = HIWORD(wParam);
					   }
					   else if (action == SB_LINERIGHT) {
						   pos = min(pos + 3, 999);
					   }
					   else if (action == SB_LINELEFT) {
						   pos = max(pos - 3, 0);
					   }
					   else if (action == SB_PAGERIGHT) {
						   pos = max(pos + 30, 0);
					   }
					   else if (action == SB_PAGELEFT) {
						   pos = max(pos - 30, 0);
					   }

					   char str[30];
					   SetScrollPos(hScroll, SB_CTL, pos, TRUE);
					   F32 x = (max(950 - pos, 0));
					   x = x*x*x*x;
					   gData.sleepInterval = 1e-9 * x;
					   wsprintf(str, "%d %d", pos, gData.sleepInterval);
					   //SetWindowText(hEdit, str);

	}
		break;
	case WM_COMMAND:
	{
					   switch LOWORD(wParam)
					   {
					   case 0://RUN
						   EnterCriticalSection(&gData.cs);

						   if (gData.status == DONE)
						   {
							   WaitForSingleObject(gData.threadHandle, INFINITE);
							   CloseHandle(gData.threadHandle);
							   DWORD mythreadid;
							   gData.threadHandle = CreateThread(0, 0, beast_thread, NULL, 0, &mythreadid);


							   gData.status = RUN;
							   EnableButtons(hButton, &hDiag);


							   CreateGDIObject(hwnd, 5);
							   InvalidateRect(hwnd, NULL, TRUE);

							   while (gData.y == NULL)
								   SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);

							   {
								   F32PTR  y = gData.y;
								   F32 W, H;
								   W = gData.w[0];
								   H = gData.h[0];
								   F32 a1 = H + H / (gData.yMax - gData.yMin) *gData.yMin;
								   F32 b1 = H / (gData.yMax - gData.yMin);
								   RECT  rect;
								   SetRect(&rect, 0, 0, W, H);
								   FillRect(bufferDC[0], &rect, blackBrush);
								   SelectObject(bufferDC[0], GetStockObject(NULL_BRUSH));
								   int RAIDUS = 4;
								   int curRowIdx = 0;
								   for (int i = 0; i < gData.N; i++)
								   {
									   if ((i + 1) == gData.rowsMissing[curRowIdx] && (curRowIdx + 1) <= gData.nMissing)
									   {
										   curRowIdx++;
										   continue;

									   }

									   F32 xx = ((F32)(i + 1) / gData.N) *W;
									   F32 yy = a1 - y[i] * b1;
									   if (yy == yy)
										   Ellipse(bufferDC[0], xx - RAIDUS, yy - RAIDUS, xx + RAIDUS, yy + RAIDUS);
								   }
							   }

							   LeaveCriticalSection(&gData.cs);

							   /*
							   EnterCriticalSection(&gData.cs);
								   HDC hdc = GetDC(hwnd); 
								   BEAST_GeneratePlotData();
								   BEAST_DrawPlots(hdc);
								   ReleaseDC(hwnd, hdc);
							   LeaveCriticalSection(&gData.cs);
							   */

							   //-------------------------------------------------------
							   {
								   RECT rc;
								   GetWindowRect(hwnd, &rc);
								   MoveWindow(hDiag, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);

								   //Enable resizing the window
								   SetWindowLong(hwnd, GWL_STYLE, GetWindowLong(hwnd, GWL_STYLE) | (WS_SIZEBOX));
							   }

							   char str[100];
							   wsprintf(str, "Re-running....\r\n");
							   LoggerInsert(str);
						   }
						   else if (gData.status == PAUSE)
						   {
							   gData.status = RUN;
							   EnableButtons(hButton, &hDiag);
							   LeaveCriticalSection(&gData.cs);
							   WakeConditionVariable(&gData.cv);

							   char str[100];
							   wsprintf(str, "Resumed....\r\n");
							   LoggerInsert(str);
						   }
						   SetTimer(hwnd, ID_TIMER, gData.timerInterval, NULL);

						   break;
					   case 1://PAUSE
						   KillTimer(hwnd, ID_TIMER);
						   EnterCriticalSection(&gData.cs);
						   gData.status = PAUSE;
						   EnableButtons(hButton, &hDiag);
						   LeaveCriticalSection(&gData.cs);
						   //printf("%d", wParam);
						   char str[100];
						   wsprintf(str, "PAUSED....\r\n");
						   LoggerInsert(str);
						   break;
					   case 2://SETTING
						   if (IsWindowVisible(hDiag))
							   ShowWindow(hDiag, SW_HIDE);
						   else
							   ShowWindow(hDiag, SW_SHOWDEFAULT);

						   UpdateWindow(hDiag);
						   break;
					   case 3://QUIT
						   EnterCriticalSection(&gData.cs);
						   gData.quit = 1;
						   LeaveCriticalSection(&gData.cs);
						   if (gData.status == PAUSE)
							   WakeConditionVariable(&gData.cv);
						   DestroyWindow(hwnd);
						   break;
					   case IDC_MAIN_EDIT:
						   return DefWindowProc(hwnd, msg, wParam, lParam);
						   break;
					   default:
						   return DefWindowProc(hwnd, msg, wParam, lParam);
					   }


	}
		break;
	case WM_GETMINMAXINFO:
	{//http: //www.owlnext.net/old/ttt15.html
		// set the MINMAXINFO structure pointer 
		PMINMAXINFO lpMinMaxInfo = (MINMAXINFO*)lParam;
		lpMinMaxInfo->ptMinTrackSize.x = 600;
		lpMinMaxInfo->ptMinTrackSize.y = 500;
		lpMinMaxInfo->ptMaxTrackSize.x = 5000;
		lpMinMaxInfo->ptMaxTrackSize.y = 5000;
	}
		break;
	case WM_MOVE:
	{
		RECT rc;
		GetWindowRect(hwnd, &rc);
		MoveWindow(hDiag, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);
	}
		break;
	case WM_USER + 1:
		KillTimer(hwnd, ID_TIMER);
		//Disable resizing the window
		SetWindowLong(hwnd, GWL_STYLE, GetWindowLong(hwnd, GWL_STYLE) & ~(WS_SIZEBOX));

		EnterCriticalSection(&gData.cs);
			gData.status = DONE;
			gData.ite = 0;
			EnableButtons(hButton, &hDiag);
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);

		{
			char str[100];
			wsprintf(str, "Finished....\r\n");
			LoggerInsert(str);
		}


		break;
	case WM_USER + 2:

		EnterCriticalSection(&gData.cs);
		{
			char str[100];
			wsprintf(str, "Chain #%d is finished\r\n", gData.curChainNumber);
			LoggerInsert(str);
			snprintf(str, 99, "Mean number of scp is %8.2f\r\n", gData.sN);
			LoggerInsert(str);
			snprintf(str,99, "Mean number of tcp is %8.2f\r\n", gData.tN);
			LoggerInsert(str);
		}
		LeaveCriticalSection(&gData.cs);

		break;
	case WM_SIZE:
	{
					//HWND hEdit = GetDlgItem(hwnd, IDC_MAIN_EDIT);					
					EnterCriticalSection(&gData.cs);

					if (gData.status == DONE)
					{
						LeaveCriticalSection(&gData.cs);
						break;
					}

					CreateGDIObject(hwnd,5);
					InvalidateRect(hwnd, NULL, TRUE);

					while (gData.y == NULL)
						SleepConditionVariableCS(&gData.cv, &gData.cs, INFINITE);

					{
						F32PTR  y = gData.y;
						F32 W, H;
						W = gData.w[0];
						H = gData.h[0];
						F32 a1 = H + H / (gData.yMax - gData.yMin) *gData.yMin;
						F32 b1 = H / (gData.yMax - gData.yMin);
						RECT  rect;
						SetRect(&rect, 0, 0, W, H);
						FillRect(bufferDC[0], &rect, blackBrush);
						SelectObject(bufferDC[0], GetStockObject(NULL_BRUSH));
						int RAIDUS = 4;
						int curRowIdx = 0;
						for (int i = 0; i < gData.N; i++)
						{					
							if ((i + 1) == gData.rowsMissing[curRowIdx] && (curRowIdx+1)<=gData.nMissing)
							{
								curRowIdx++;
								continue;

							}

							F32 xx = ((F32)(i + 1) / gData.N) *W;
							F32 yy = a1 - y[i] * b1;
							if (yy == yy)
								Ellipse(bufferDC[0], xx - RAIDUS, yy - RAIDUS, xx + RAIDUS, yy + RAIDUS);
						}
					}
					LeaveCriticalSection(&gData.cs);

					/*
					{	HDC hdc = GetDC(hwnd);
					EnterCriticalSection(&gData.cs);
					BEAST_GeneratePlotData();
					BEAST_DrawPlots(hdc);
					LeaveCriticalSection(&gData.cs);
					ReleaseDC(hwnd, hdc);
					}
					*/


					//SetWindowPos(hEdit, NULL, 0, rcClient.bottom - height, rcClient.right, height, SWP_NOZORDER);
					//InvalidateRect(hwnd, &rcClient, FALSE);
					//ValidateRect(hwnd, &rcClient);
					{
						RECT rc;
						GetWindowRect(hwnd, &rc);
						MoveWindow(hDiag, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);

					}

	}
		break;
	case WM_SHOWWINDOW:
		if (wParam)
		{
			//SendMessage(hDiag, WM_USER + 1, (WPARAM)NULL, (LPARAM)NULL); //MAKELPARAM(low, high)

		}
		break;
	case WM_CLOSE:
		EnterCriticalSection(&gData.cs);
		gData.quit = 1;
		LeaveCriticalSection(&gData.cs);
		if (gData.status == PAUSE)
			WakeConditionVariable(&gData.cv);
		DestroyWindow(hwnd);
		break;
	case WM_CTLCOLORSTATIC:
	{//https: //stackoverflow.com/questions/2262628/changing-background-of-text-in-edit-control
							  if ((HWND)lParam == hEdit)
							  {
								  HDC hdc = (HDC)wParam;
								  //SetBkMode((HDC)wParam, TRANSPARENT);
								  SetBkColor(hdc, RGB(100, 0, 0));
								  SetTextColor(hdc, RGB(0, 255, 0));
								  return (LRESULT)GetStockObject(GRAY_BRUSH);
							  }
							   else if ((HWND)lParam == hStatic[0] || (HWND)lParam == hStatic[1]) 		  {
								  HDC hdc = (HDC)wParam;
								  SetBkMode((HDC)wParam, TRANSPARENT);
								  SetBkColor(hdc, RGB(100, 0, 0));
								  SetTextColor(hdc, RGB(255, 255, 0));
								  return GetStockObject(GRAY_BRUSH);
							  }

	}
		break;
	case WM_PAINT:
	{

					 RECT rcClient;
					 PAINTSTRUCT ps;
					 HDC hdc;

					 hdc = BeginPaint(hwnd, &ps);
					 //GetClientRect(hwnd, &rcClient);
					 EnterCriticalSection(&gData.cs);
					 if (gData.status == PAUSE || gData.status == DONE) 	 BEAST2_DrawPlots(hdc);
					 LeaveCriticalSection(&gData.cs);
					 EndPaint(hwnd, &ps);
					 /*
					 HBRUSH hBnew, hBold;
					 hdc = GetDC(hwnd);
					 //hdc = BeginPaint(hwnd, &ps);
					 char str[100];

					 SetWindowText(GetDlgItem(hwnd, IDC_MAIN_EDIT), str);

					 TextOut(hdc, 400, 50, str, lstrlen(str));
					 GetClientRect(hwnd, &rcClient);
					 GetUpdateRect(hwnd, &rcClient, FALSE);
					 //Rectangle(hdc, rcClient.left, rcClient.top, rcClient.right, rcClient.bottom);

					 TextOut(hdc, 200, 150, str, lstrlen(str));
					 //DrawIcon(hdc, 400, 20, hIcon);
					 //FrameRect(hdc, &rcClient, hBnew);
					 ValidateRect(hwnd, &rcClient);
					 SelectObject(hdc, hBold);
					 ReleaseDC(hwnd, hdc);
					 //EndPaint(hwnd, &ps);
					 DeleteObject(hBnew);

					 //GetWindowRect(GetDlgItem(hwnd, IDC_MAIN_EDIT), &rcClient);
					 InvalidateRect(GetDlgItem(hwnd, IDC_MAIN_EDIT), (RECT *)NULL,FALSE);
					 InvalidateRect(hButton[0], (RECT *)NULL, FALSE);
					 // GetUpdateRect(GetDlgItem(hwnd, IDC_MAIN_BUTTON), &rcClient, FALSE);
					 //InvalidateRect(GetDlgItem(hwnd, IDC_MAIN_BUTTON), (RECT *)NULL, FALSE);
					 //return DefWindowProc(hwnd, msg, wParam, lParam);
					 */

	}
		break;
	case WM_TIMER:
	{
					 /*
					 RECT rcClient;
					 HDC hdc = GetDC(hwnd);
					 GetClientRect(hwnd, &rcClient);
					 InvalidateRect(hwnd, &rcClient, FALSE);
					 //TextOut(hdc, (F32)rand()/ RAND_MAX *rcClient.right, (F32)rand() / RAND_MAX *rcClient.bottom, "12345test!!", 10);
					 ReleaseDC(hwnd, hdc);
					 */
					 HDC hdc = GetDC(hwnd);
					 EnterCriticalSection(&gData.cs);					 
					 BEAST2_DrawPlots(hdc);
					 LeaveCriticalSection(&gData.cs);
					 ReleaseDC(hwnd, hdc);

					 /*
					 {
						 hdc = GetWindowDC(hwnd);
						 char str[100];
						 wsprintf(str, "%d", gData.ite);
						 TextOut(hdc, 100, 5, str, strlen(str));
					     ReleaseDC(hwnd, hdc);
					 }
					 */
					 
					 // GetClientRect(hwnd, &rcClient);
					 //InvalidateRect(hwnd, &rcClient, FALSE);		 
	}
		break;
	case WM_DRAWITEM:
	{
						//获得绘制结构体，包含绘制的按钮DC和当前按钮状态等  
						LPDRAWITEMSTRUCT pdis = (LPDRAWITEMSTRUCT)lParam;

						if (pdis->CtlType == ODT_BUTTON)//只绘制button类型  
						{
							HDC hdc = pdis->hDC;
							SaveDC(hdc);//保存DC,绘制完必须恢复默认  

							TextOut(hdc, 0, 0, "789", 3);
							//绘制获取焦点时状态  
							if (pdis->itemState & ODS_FOCUS)
							{
								TextOut(hdc, 0, 0, "123", 3);

							}

							//绘制下压状态  
							if (pdis->itemState & ODS_SELECTED)
							{
								TextOut(hdc, 0, 0, "456", 3);

								//SelectObject(hMemDc, pdis->CtlID == IDOK ? hBitmapOK_D : hBitmapCANCEL_D);

							}

							RestoreDC(hdc, -1);
						}
	}
		break;
	case WM_DESTROY:
		KillTimer(hwnd, ID_TIMER);
		DeleteGDIObject(5);

		WaitForSingleObject(gData.threadHandle, INFINITE);
		CloseHandle(gData.threadHandle);
		gData.threadHandle = NULL;

		for (int i = 0; i < sizeof(gData.plotData) / sizeof(gData.plotData[0][0]); i++)
		{
			if (*((int **)gData.plotData + i) != NULL) free(*((int **)gData.plotData + i));
		}
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hwnd, msg, wParam, lParam);
	}
	return 0;
}
static LRESULT CALLBACK DialogProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	static byte  isDrawing = 0;
	static char* names[] = { "seasonMinKnotNum", "seasonMaxKnotNum", "seasonMinOrder",  "seasonMaxOrder",  "seasonMinSepDist",
							  "trendMinKnotNum", "trendMaxKnotNum", "trendMinOrder",  "trendMaxOrder",  "trendMinSepDist",
							   "burnin", "samples", "chainNumber", "thinningFactor","maxMoveStepSize" };
		char   str[100];	
	BEAST2_OPTIONS_PTR opt = GLOBAL_OPTIONS;

	switch (msg) {
	case WM_CREATE:
	{
		int values[] = {
						opt->prior.seasonMinKnotNum,
						opt->prior.seasonMaxKnotNum,
						opt->prior.seasonMinOrder,
						opt->prior.seasonMaxOrder,
						opt->prior.seasonMinSepDist,

						opt->prior.trendMinKnotNum,
						opt->prior.trendMaxKnotNum,
						opt->prior.trendMinOrder,
						opt->prior.trendMaxOrder,
						opt->prior.trendMinSepDist,

						opt->mcmc.burnin,
						opt->mcmc.samples,
						opt->mcmc.chainNumber,
						opt->mcmc.thinningFactor,
						opt->mcmc.maxMoveStepSize,
						
		};


		HWND  htmp;
		HFONT hfDefault = GetStockObject(SYSTEM_FONT);
		
		//HWND  hFrame = CreateWindowEx(0, "STATIC", "XXX", WS_CHILD | WS_VISIBLE | SS_LEFT| SS_BLACKFRAME,	10, 10, 200, 900, hwnd, (HMENU)400, GetModuleHandle(NULL), NULL);
		//SendDlgItemMessage(hwnd, 400, WM_SETTEXT, (WPARAM)NULL, (LPARAM)"Prior");
		for (int i = 0; i < 15; i++)
		{	  //WS_EX_DLGMODALFRAME WS_EX_CLIENTEDGE

			htmp = CreateWindowEx(0, "STATIC", names[i], WS_CHILD | WS_VISIBLE | SS_LEFT,
				10, 10, 40, 40, hwnd, (HMENU)i, GetModuleHandle(NULL), NULL);
			SendMessage(htmp, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM(FALSE, 0));
			SendDlgItemMessage(hwnd, i, WM_SETTEXT, (WPARAM)NULL, (LPARAM)names[i]);

			htmp = CreateWindowEx(0, "EDIT", "434", WS_CHILD | WS_VISIBLE | ES_LEFT | WS_BORDER | ES_NUMBER,
				60, 10, 50, 50, hwnd, (HMENU)(100 + i), GetModuleHandle(NULL), NULL);
			SendMessage(htmp, WM_SETFONT, (WPARAM)hfDefault, MAKELPARAM(FALSE, 0));
			SetDlgItemInt(hwnd, 100 + i, values[i], 1);

		}

		htmp = CreateWindowEx(WS_EX_CLIENTEDGE, "BUTTON", "Apply",
			WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | BS_PUSHLIKE  /*| BS_OWNERDRAW*/,
			10, 70, 30 * 9 / 4, 30, hwnd, (HMENU)200, GetModuleHandle(NULL), NULL);

		htmp = CreateWindowEx(WS_EX_CLIENTEDGE, "BUTTON", "Close",
			WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON | BS_PUSHLIKE  /*| BS_OWNERDRAW*/,
			10, 70, 30 * 9 / 4, 30, hwnd, (HMENU)201, GetModuleHandle(NULL), NULL);

		style.numRows = 15 + 1;

		/*
		U08  minSeasonOrder, maxSeasonOrder, minTrendOrder, maxTrendOrder;
		U16 trendMinSepDist, seasonMinSepDist;
		U16 trendMaxKnotNum, seasonMaxKnotNum;
		U16 maxMoveStepSize;
		F32    trendResamplingOrderProb;
		F32    seasonResamplingOrderProb;
		U32 burnin, samples, chainNumber;
		*/
		wsprintf(str, "Diag Created..%d \r\n", GetTickCount());
		LoggerInsert(str);

	}

		break;
	case WM_USER + 1:
		isDrawing = 1;
		ShowWindow(hwnd, SW_SHOWDEFAULT);
		UpdateWindow(hwnd);
		break;
	case WM_MOVE:
	{
		if (isDrawing)
		{
			RECT rc;
			GetWindowRect(gData.hwnd, &rc);
			MoveWindow(hwnd, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);
			BEAST2_ResizeDialogControls(hwnd);
			wsprintf(str, "WM_MOV_created.Inside.%d \r\n", GetTickCount());
			LoggerInsert(str);
		}

		break;
	}		
	case WM_SIZE:
	{
		if (isDrawing)
		{
			RECT rc;
			GetWindowRect(gData.hwnd, &rc);
			MoveWindow(hwnd, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);

			BEAST2_ResizeDialogControls(hwnd);
			wsprintf(str, "WM_SIZE_created.INsideside.%d \r\n", GetTickCount());
			LoggerInsert(str);
		}

		break;

	}
	
		/*
		case WM_GETMINMAXINFO:
		{//http: //www.owlnext.net/old/ttt15.html
		// set the MINMAXINFO structure pointer
		PMINMAXINFO lpMinMaxInfo = (MINMAXINFO *)lParam;
		lpMinMaxInfo->ptMinTrackSize.x = 600;
		lpMinMaxInfo->ptMinTrackSize.y = 500;
		lpMinMaxInfo->ptMaxTrackSize.x = 5000;
		lpMinMaxInfo->ptMaxTrackSize.y = 5000;
		}
		break;
		*/
	case WM_SHOWWINDOW:
	{

		if (wParam)
		{//the diaglog shows up
			RECT rc;
			GetWindowRect(gData.hwnd, &rc);
			MoveWindow(hwnd, rc.right + 20, rc.top, style.widthDialg, rc.bottom - rc.top, TRUE);

			EnableButtons(NULL, &hDiag);
			BEAST2_ResizeDialogControls(hwnd);
			int values[] = {
					opt->prior.seasonMinKnotNum,
					opt->prior.seasonMaxKnotNum,
					opt->prior.seasonMinOrder,
					opt->prior.seasonMaxOrder,
					opt->prior.seasonMinSepDist,

					opt->prior.trendMinKnotNum,
					opt->prior.trendMaxKnotNum,
					opt->prior.trendMinOrder,
					opt->prior.trendMaxOrder,
					opt->prior.trendMinSepDist,

						opt->mcmc.burnin,
						opt->mcmc.samples,
						opt->mcmc.chainNumber,
						opt->mcmc.thinningFactor,
						opt->mcmc.maxMoveStepSize,
			};

			for (int i = 0; i < 8; i++)
			{	  //WS_EX_DLGMODALFRAME WS_EX_CLIENTEDGE
				SetDlgItemInt(hwnd, 100 + i, values[i], 1);
			}

			BEAST2_ResizeDialogControls(hwnd);

			//SetDlgItemText(hwnd, 2, itoa(opt->burnin, buf, 10));

		}
		else
		{//the diaglog hides
			/*
			// GetDlgItemText(hwnd, 2, buf, 19);
			//opt->burnin = atoi(buf);
			byte translated = 0;
			int i = 100;
			opt->trendMaxKnotNum = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->minTrendOrder = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->maxTrendOrder = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->trendMinSepDist = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->maxMoveStepSize = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->burnin = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->samples = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->chainNumber = GetDlgItemInt(hwnd, i, &translated, 0); i++;

			//----------------------------
			*/
			gData.optStatus = NeedUpdate;
		}



	}
		break;

	case WM_COMMAND:
		switch (LOWORD(wParam))
		{
		case 200://Apply 
		{
			BOOL translated = 0;
			int i = 100;
			opt->prior.seasonMinKnotNum = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.seasonMaxKnotNum = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.seasonMinOrder   = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.seasonMaxOrder   = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.seasonMinSepDist = GetDlgItemInt(hwnd, i, &translated, 0); i++;

			opt->prior.trendMinKnotNum = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.trendMaxKnotNum = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.trendMinOrder   = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.trendMaxOrder   = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->prior.trendMinSepDist = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			
			opt->mcmc.burnin = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->mcmc.samples = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->mcmc.chainNumber = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->mcmc.thinningFactor = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			opt->mcmc.maxMoveStepSize = GetDlgItemInt(hwnd, i, &translated, 0); i++;
			//----------------------------

			gData.optStatus = NeedUpdate;


		}
			break;
		case 201://Close
			ShowWindow(hwnd, SW_HIDE);
			break;
		default:
			//return (DefWindowProc(hwnd, msg, wParam, lParam));
			break;
		}

		break;

	case WM_CLOSE:
		ShowWindow(hwnd, SW_HIDE);
		//DestroyWindow(hwnd);
		break;
	default:
		if (IsWindowUnicode(hwnd))
			return (DefWindowProc(hwnd, msg, wParam, lParam)); //DefWindowProcW
		else
			return (DefWindowProc(hwnd, msg, wParam, lParam)); //DefWindowProcA		

	}
	return 0;

}
#else
static char fileID UNUSED_DECORATOR = 'c';
#endif


 
#include "abc_000_warning.h"