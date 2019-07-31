#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
#if defined(WIN64_OS) 
#include <stdio.h>   
#include <io.h>      
#include <math.h>
#include "demo.h"
static void             DialogClassRegister();
static LRESULT CALLBACK DialogProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam);
static LRESULT CALLBACK WndProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	switch (msg) {
	case WM_CREATE:
	{
		{
			EnterCriticalSection(&gData.cs);
			while (gData.N <=0)
				SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
			LeaveCriticalSection(&gData.cs);
			EnterCriticalSection(&gData.cs);
			AllocatePlotData_ST();
			LeaveCriticalSection(&gData.cs);
			WakeConditionVariable(&gData.cv);
		}
					  HFONT hfDefault;
					  hfDefault=GetStockObject(SYSTEM_FONT);
					  if (!SetTimer(hwnd,ID_TIMER,gData.timerInterval,NULL))   MessageBox(hwnd,"Could not SetTimer()!","Error",MB_OK|MB_ICONEXCLAMATION);
					  {
						  hEdit=CreateWindowEx(WS_EX_CLIENTEDGE,"EDIT","",
							  WS_CHILD|WS_VISIBLE|WS_VSCROLL|WS_HSCROLL|ES_MULTILINE|ES_AUTOVSCROLL|ES_AUTOHSCROLL,
							  0,0,100,100,hwnd,(HMENU)IDC_MAIN_EDIT,GetModuleHandle(NULL),NULL);
						  if (hEdit==NULL)   MessageBox(hwnd,"Could not create edit box.\0\0\0","Error\0\0\0",MB_OK|MB_ICONERROR);
						  SendMessage(hEdit,WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
						  SetWindowText(hEdit,"");
						  SendDlgItemMessage(hwnd,IDC_MAIN_EDIT,EM_SETREADONLY,(WPARAM)TRUE,0);
					  }
					  {
						  char *caption[]={ "Run\0\0\0","Pause\0\0\0","Setting\0\0\0","Exit\0\0\0" };
						  for (int i=0; i < 4; i++)
						  {
							  hButton[i]=CreateWindowEx(WS_EX_CLIENTEDGE,"BUTTON",caption[i],
								  WS_CHILD|WS_VISIBLE|BS_PUSHBUTTON|BS_PUSHLIKE,
								  10+100 * i,20,30 * 9/4,30,hwnd,(HMENU) i,GetModuleHandle(NULL),NULL);
							  if (hButton[i]==NULL)  MessageBox(hwnd,"Could not create edit box.\0\0\0","Error",MB_OK|MB_ICONERROR);
							  SendMessage(hButton[i],WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
							  SendMessage(hButton[0],BM_SETCHECK,(WPARAM)BST_UNCHECKED,MAKELPARAM(FALSE,0));
						  }
						  EnableButtons(hButton,NULL);
					  }
					  {
						  SCROLLINFO si;
						  hScroll=CreateWindowEx(
							  0,                      
							  "SCROLLBAR",           
							  (PTSTR)NULL,           
							  WS_CHILD|WS_VISIBLE   
							  |SBS_HORZ,         
							  400,              
							  50, 
							  100,             
							  10,               
							  hwnd,             
							  (HMENU)10,           
							  GetModuleHandle(NULL),                
							  (PVOID)NULL            
							  );
						  si.cbSize=sizeof(si);
						  si.nMin=0;
						  si.nMax=999;
						  si.nPos=500;
						  si.nPage=50;
						  si.fMask=SIF_RANGE|SIF_POS|SIF_PAGE;
						  SetScrollInfo(hScroll,SB_CTL,&si,TRUE);
						  float x=(max(950 - si.nPos,0));
						  x=x*x*x*x;
						  gData.sleepInterval=1e-9 * x;
					  }
					  {
						  hStatic[0]=CreateWindowEx(0,"STATIC","SLOW",WS_CHILD|WS_VISIBLE|SS_CENTER,
							  400,100,100,100,hwnd,(HMENU)20,GetModuleHandle(NULL),NULL);
						  SendMessage(hStatic[0],WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
						  SetWindowText(hStatic[0],"Slow");
						  hStatic[1]=CreateWindowEx(0,"STATIC","FAST",WS_CHILD|WS_VISIBLE|WS_TABSTOP,
							  500,100,100,100,hwnd,(HMENU)21,GetModuleHandle(NULL),NULL);
						  SendMessage(hStatic[1],WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
						  SetWindowText(hStatic[1],"Fast");
					  }
					  {
						  TEXTMETRIC tm;
						  HDC hdc;
						  hdc=GetDC(hStatic[0]);
						  GetTextMetrics(hdc,&tm);
						  style.xChar=tm.tmAveCharWidth;
						  style.xCapChar=(tm.tmPitchAndFamily & 1 ? 3 : 2) * style.xChar/2;
						  style.yChar=tm.tmHeight+tm.tmExternalLeading;
						  ReleaseDC(hStatic[0],hdc);
					  }
					  {
						  EnterCriticalSection(&gData.cs);
						  CreateGDIObject(hwnd,5);
						  LeaveCriticalSection(&gData.cs);
					  }
					  {
						  DialogClassRegister();
						  hDiag=CreateWindowEx(
							  WS_EX_CLIENTEDGE,
							  "myDiaglogClass3434",
							  "Settings",
							  WS_OVERLAPPEDWINDOW|(WS_CAPTION & ~WS_VISIBLE),
							  CW_USEDEFAULT,CW_USEDEFAULT,320,240,
							  hwnd,(HMENU)NULL,GetModuleHandle(NULL),NULL);
						  if (hDiag==0)
							  MessageBox(NULL,"Window Registration Failed!","Error!",MB_ICONEXCLAMATION|MB_OK);
						  ShowWindow(hDiag,SW_HIDE);
					  }
	}
		break;
	case WM_HSCROLL:
	{
					   HWND hScroll;
					   WORD action; 
					   hScroll=(HWND)lParam;
					   action=LOWORD(wParam);
					   int pos;
					   pos=GetScrollPos(hScroll,SB_CTL);
					   if (action==SB_THUMBPOSITION||action==SB_THUMBTRACK) {
						   pos=HIWORD(wParam);
					   }
					   else if (action==SB_LINERIGHT) {
						   pos=min(pos+3,999);
					   }
					   else if (action==SB_LINELEFT) {
						   pos=max(pos - 3,0);
					   }
					   else if (action==SB_PAGERIGHT) {
						   pos=max(pos+30,0);
					   }
					   else if (action==SB_PAGELEFT) {
						   pos=max(pos - 30,0);
					   }
					   char str[30];
					   SetScrollPos(hScroll,SB_CTL,pos,TRUE);
					   float x=(max(950 - pos,0));
					   x=x*x*x*x;
					   gData.sleepInterval=1e-9 * x;
					   wsprintf(str,"%d%d",pos,gData.sleepInterval);
	}
		break;
	case WM_COMMAND:
	{
					   switch LOWORD(wParam)
					   {
					   case 0:
						   EnterCriticalSection(&gData.cs);
						   if (gData.status==DONE)
						   {
							   WaitForSingleObject(gData.threadHandle,INFINITE);
							   CloseHandle(gData.threadHandle);
							   DWORD mythreadid;
							   gData.threadHandle=CreateThread(0,0,beastST_demo,&threadParam,0,&mythreadid);
							   gData.status=RUN;
							   EnableButtons(hButton,&hDiag);
							   CreateGDIObject(hwnd,5);
							   InvalidateRect(hwnd,NULL,TRUE);
							   while (gData.y==NULL)
								   SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
							   {
								   float * y=gData.y;
								   float W,H;
								   W=gData.w[0];
								   H=gData.h[0];
								   float a1=H+H/(gData.yMax - gData.yMin) *gData.yMin;
								   float b1=H/(gData.yMax - gData.yMin);
								   RECT  rect;
								   SetRect(&rect,0,0,W,H);
								   FillRect(bufferDC[0],&rect,blackBrush);
								   SelectObject(bufferDC[0],GetStockObject(NULL_BRUSH));
								   int RAIDUS=4;
								   int curRowIdx=0;
								   for (int i=0; i < gData.N; i++)
								   {
									   if ((i+1)==gData.rowsMissing[curRowIdx] && (curRowIdx+1) <=gData.nMissing)
									   {
										   curRowIdx++;
										   continue;
									   }
									   float xx=((float)(i+1)/gData.N) *W;
									   float yy=a1 - y[i] * b1;
									   if (yy==yy)
										   Ellipse(bufferDC[0],xx - RAIDUS,yy - RAIDUS,xx+RAIDUS,yy+RAIDUS);
								   }
							   }
							   LeaveCriticalSection(&gData.cs);
							   {
								   RECT rc;
								   GetWindowRect(hwnd,&rc);
								   MoveWindow(hDiag,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
								   SetWindowLong(hwnd,GWL_STYLE,GetWindowLong(hwnd,GWL_STYLE)|(WS_SIZEBOX));
							   }
							   char str[100];
							   wsprintf(str,"Re-running....\r\n");
							   LoggerInsert(str);
						   }
						   else if (gData.status==PAUSE)
						   {
							   gData.status=RUN;
							   EnableButtons(hButton,&hDiag);
							   LeaveCriticalSection(&gData.cs);
							   WakeConditionVariable(&gData.cv);
							   char str[100];
							   wsprintf(str,"Resumed....\r\n");
							   LoggerInsert(str);
						   }
						   SetTimer(hwnd,ID_TIMER,gData.timerInterval,NULL);
						   break;
					   case 1:
						   KillTimer(hwnd,ID_TIMER);
						   EnterCriticalSection(&gData.cs);
						   gData.status=PAUSE;
						   EnableButtons(hButton,&hDiag);
						   LeaveCriticalSection(&gData.cs);
						   char str[100];
						   wsprintf(str,"PAUSED....\r\n");
						   LoggerInsert(str);
						   break;
					   case 2:
						   if (IsWindowVisible(hDiag))
							   ShowWindow(hDiag,SW_HIDE);
						   else
							   ShowWindow(hDiag,SW_SHOWDEFAULT);
						   UpdateWindow(hDiag);
						   break;
					   case 3:
						   EnterCriticalSection(&gData.cs);
						   gData.quit=1;
						   LeaveCriticalSection(&gData.cs);
						   if (gData.status==PAUSE)
							   WakeConditionVariable(&gData.cv);
						   DestroyWindow(hwnd);
						   break;
					   case IDC_MAIN_EDIT:
						   return DefWindowProc(hwnd,msg,wParam,lParam);
						   break;
					   default:
						   return DefWindowProc(hwnd,msg,wParam,lParam);
					   }
	}
		break;
	case WM_GETMINMAXINFO:
	{
							 PMINMAXINFO lpMinMaxInfo=(MINMAXINFO *)lParam;
							 lpMinMaxInfo->ptMinTrackSize.x=600;
							 lpMinMaxInfo->ptMinTrackSize.y=500;
							 lpMinMaxInfo->ptMaxTrackSize.x=5000;
							 lpMinMaxInfo->ptMaxTrackSize.y=5000;
	}
		break;
	case WM_MOVE:
	{
					RECT rc;
					GetWindowRect(hwnd,&rc);
					MoveWindow(hDiag,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
	}
		break;
	case WM_USER+1:
		KillTimer(hwnd,ID_TIMER);
		SetWindowLong(hwnd,GWL_STYLE,GetWindowLong(hwnd,GWL_STYLE) & ~(WS_SIZEBOX));
		EnterCriticalSection(&gData.cs);
			gData.status=DONE;
			gData.ite=0;
			EnableButtons(hButton,&hDiag);
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);
		{
			char str[100];
			wsprintf(str,"Finished....\r\n");
			LoggerInsert(str);
		}
		break;
	case WM_USER+2:
		EnterCriticalSection(&gData.cs);
		{
			char str[100];
			wsprintf(str,"Chain #%d is finished\r\n",gData.curChainNumber);
			LoggerInsert(str);
			sprintf(str,"Mean number of scp is%8.2f\r\n",gData.sN);
			LoggerInsert(str);
			sprintf(str,"Mean number of tcp is%8.2f\r\n",gData.tN);
			LoggerInsert(str);
		}
		LeaveCriticalSection(&gData.cs);
		break;
	case WM_SIZE:
	{
					EnterCriticalSection(&gData.cs);
					if (gData.status==DONE)
					{
						LeaveCriticalSection(&gData.cs);
						break;
					}
					CreateGDIObject(hwnd,5);
					InvalidateRect(hwnd,NULL,TRUE);
					while (gData.y==NULL)
						SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
					{
						float * y=gData.y;
						float W,H;
						W=gData.w[0];
						H=gData.h[0];
						float a1=H+H/(gData.yMax - gData.yMin) *gData.yMin;
						float b1=H/(gData.yMax - gData.yMin);
						RECT  rect;
						SetRect(&rect,0,0,W,H);
						FillRect(bufferDC[0],&rect,blackBrush);
						SelectObject(bufferDC[0],GetStockObject(NULL_BRUSH));
						int RAIDUS=4;
						int curRowIdx=0;
						for (int i=0; i < gData.N; i++)
						{					
							if ((i+1)==gData.rowsMissing[curRowIdx] && (curRowIdx+1)<=gData.nMissing)
							{
								curRowIdx++;
								continue;
							}
							float xx=((float)(i+1)/gData.N) *W;
							float yy=a1 - y[i] * b1;
							if (yy==yy)
								Ellipse(bufferDC[0],xx - RAIDUS,yy - RAIDUS,xx+RAIDUS,yy+RAIDUS);
						}
					}
					LeaveCriticalSection(&gData.cs);
					{
						RECT rc;
						GetWindowRect(hwnd,&rc);
						MoveWindow(hDiag,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
					}
	}
		break;
	case WM_SHOWWINDOW:
		if (wParam)
		{
			SendMessage(hDiag,WM_USER+1,(WPARAM)NULL,(LPARAM)NULL); 
		}
		break;
	case WM_CLOSE:
		EnterCriticalSection(&gData.cs);
		gData.quit=1;
		LeaveCriticalSection(&gData.cs);
		if (gData.status==PAUSE)
			WakeConditionVariable(&gData.cv);
		DestroyWindow(hwnd);
		break;
	case WM_CTLCOLORSTATIC:
	{
							  if ((HWND)lParam==hEdit)
							  {
								  HDC hdc=(HDC)wParam;
								  SetBkColor(hdc,RGB(100,0,0));
								  SetTextColor(hdc,RGB(0,255,0));
								  return (LRESULT)GetStockObject(GRAY_BRUSH);
							  }
							  else if ((HWND)lParam==hStatic[0]||(HWND)lParam==hStatic[1])
							  {
								  HDC hdc=(HDC)wParam;
								  SetBkMode((HDC)wParam,TRANSPARENT);
								  SetBkColor(hdc,RGB(100,0,0));
								  SetTextColor(hdc,RGB(255,255,0));
							  }
	}
		break;
	case WM_PAINT:
	{
					 RECT rcClient;
					 PAINTSTRUCT ps;
					 HDC hdc;
					 hdc=BeginPaint(hwnd,&ps);
					 EnterCriticalSection(&gData.cs);
					 if (gData.status==PAUSE||gData.status==DONE) 	 DrawPlots_ST(hdc);
					 LeaveCriticalSection(&gData.cs);
					 EndPaint(hwnd,&ps);
	}
		break;
	case WM_TIMER:
	{
					 HDC hdc=GetDC(hwnd);
					 DrawPlots_ST(hdc);
					 ReleaseDC(hwnd,hdc);
	}
		break;
	case WM_DRAWITEM:
	{
						LPDRAWITEMSTRUCT pdis=(LPDRAWITEMSTRUCT)lParam;
						if (pdis->CtlType==ODT_BUTTON)
						{
							HDC hdc=pdis->hDC;
							SaveDC(hdc);
							TextOut(hdc,0,0,"789",3);
							if (pdis->itemState & ODS_FOCUS)
							{
								TextOut(hdc,0,0,"123",3);
							}
							if (pdis->itemState & ODS_SELECTED)
							{
								TextOut(hdc,0,0,"456",3);
							}
							RestoreDC(hdc,-1);
						}
	}
		break;
	case WM_DESTROY:
		KillTimer(hwnd,ID_TIMER);
		DeleteGDIObject(5);
		WaitForSingleObject(gData.threadHandle,INFINITE);
		CloseHandle(gData.threadHandle);
		gData.threadHandle=NULL;
		for (int i=0; i < sizeof(gData.plotData)/sizeof(gData.plotData[0][0]); i++)
		{
			if (*((int **)gData.plotData+i) !=NULL) free(*((int **)gData.plotData+i));
		}
		UnregisterClass("myDiaglogClass3434",GetModuleHandle(NULL));
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hwnd,msg,wParam,lParam);
	}
	return 0;
}
void DllExport WinMainDemoST(Options * option,RESULT * result)
{
	HINSTANCE hInstance=GetModuleHandle(NULL);
	WNDCLASSEX wc;
	HWND hwnd;
	MSG  Msg;
	threadParam.GLOBAL_OPTIONS=option;
	threadParam.GLOBAL_RESULT=result;
	hIcon=0;
	InitGlobalData_ST();
	{
		DWORD  mythreadid;
		gData.threadHandle=CreateThread(0,0,beastST_demo,&threadParam,0,&mythreadid);
	}
	logger.LEN=2000;
	logger.str=malloc(logger.LEN+1);
	logger.len=0;
	hIcon=CreateIcon(hInstance,    
		32,              
		32,              
		1,               
		1,               
		ANDmaskIcon,     
		XORmaskIcon);    
	const char g_szClassName[]="myWindowClass";
	wc.cbSize=sizeof(WNDCLASSEX);
	wc.style=CS_DBLCLKS|CS_HREDRAW|CS_DROPSHADOW;
	wc.lpfnWndProc=WndProc;
	wc.cbClsExtra=0;
	wc.cbWndExtra=0;
	wc.hInstance=hInstance;
	wc.hIcon=hIcon; 
	wc.hCursor=LoadCursor(NULL,IDC_ARROW);
	wc.hbrBackground=(HBRUSH)(COLOR_WINDOW);
	wc.lpszMenuName=NULL;
	wc.lpszClassName=(LPCSTR)g_szClassName;
	wc.hIconSm=hIcon;
	if (!RegisterClassEx(&wc))
	{
		MessageBox(NULL,"Window Registration Failed!","Error!",MB_ICONEXCLAMATION|MB_OK);
		return;
	}
	hwnd=CreateWindowEx(
		WS_EX_CLIENTEDGE,
		g_szClassName,
		"Bayesian nonlinear trend analysis",
		(WS_OVERLAPPEDWINDOW & ~WS_MAXIMIZEBOX)|WS_THICKFRAME,
		CW_USEDEFAULT,CW_USEDEFAULT,520,740,
		NULL,NULL,hInstance,NULL);
	gData.hwnd=hwnd;
	if (!hwnd)	{
		MessageBox(NULL,"Window Creation Failed!","Error!",MB_ICONEXCLAMATION|MB_OK);
		return;
	}
	ShowWindow(hwnd,SW_SHOWDEFAULT);
	UpdateWindow(hwnd);
	SetWindowPos(hwnd,HWND_TOPMOST,0,0,0,0,SWP_NOMOVE|SWP_NOSIZE);
	while (GetMessage(&Msg,NULL,0,0) > 0)
	{
		if (!IsDialogMessage(hDiag,&Msg))
		{
			TranslateMessage(&Msg);
			DispatchMessage(&Msg);
		}
	}
	WaitForSingleObject(gData.threadHandle,INFINITE);
	CloseHandle(gData.threadHandle);
	DeleteCriticalSection(&gData.cs);
	DestroyIcon(hIcon);
	UnregisterClass(g_szClassName,hInstance);
	free(logger.str);
	return ;
}
static void DialogClassRegister()
{
	const char g_szClassName[]="myDiaglogClass3434";
	WNDCLASSEX wc={ 0 };
	wc.cbSize=sizeof(WNDCLASSEX);
	wc.style=CS_DBLCLKS|CS_HREDRAW|CS_DROPSHADOW;
	wc.lpfnWndProc=DialogProc;
	wc.cbClsExtra=0;
	wc.cbWndExtra=0;
	wc.hInstance=GetModuleHandle(NULL);
	wc.hIcon=NULL;
	wc.hCursor=NULL;
	wc.hbrBackground=(HBRUSH)(COLOR_WINDOW);
	wc.lpszMenuName=NULL;
	wc.lpszClassName=(LPCSTR)g_szClassName;
	wc.hIconSm=NULL;
	if (!RegisterClassEx(&wc))
	{
		MessageBox(NULL,"DialogBox Class Registration Failed!","Error!",MB_ICONEXCLAMATION|MB_OK);
	}
}
static LRESULT CALLBACK DialogProc(HWND hwnd,UINT msg,WPARAM wParam,LPARAM lParam)
{
	char str[100];
	static byte isDrawing=0;
	switch (msg) {
	case WM_CREATE:
	{
					  HFONT hfDefault=GetStockObject(SYSTEM_FONT);
					  HWND  htmp;
					  char *names[]={ "maxKnotNum_Trend","minTrendOrder","maxTrendOrder",
						  "minSepDist_Trend","maxMoveStepSize","burnin","samples","chainNumber" };
					  int values[]={ gData.opt.maxKnotNum_Trend,
						  gData.opt.minTrendOrder,
						  gData.opt.maxTrendOrder,
						  gData.opt.minSepDist_Trend,
						  gData.opt.maxMoveStepSize,
						  gData.opt.burnin,
						  gData.opt.samples,
						  gData.opt.chainNumber,
					  };
					  for (int i=0; i < 8; i++)
					  {	  
						  htmp=CreateWindowEx(0,"STATIC",names[i],WS_CHILD|WS_VISIBLE|SS_LEFT,
							  10,10,40,40,hwnd,(HMENU)i,GetModuleHandle(NULL),NULL);
						  SendMessage(htmp,WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
						  SendDlgItemMessage(hwnd,i,WM_SETTEXT,(WPARAM)NULL,(LPARAM)names[i]);
						  htmp=CreateWindowEx(0,"EDIT","434",WS_CHILD|WS_VISIBLE|ES_LEFT|WS_BORDER|ES_NUMBER,
							  60,10,50,50,hwnd,(HMENU)(100+i),GetModuleHandle(NULL),NULL);
						  SendMessage(htmp,WM_SETFONT,(WPARAM)hfDefault,MAKELPARAM(FALSE,0));
						  SetDlgItemInt(hwnd,100+i,values[i],1);
					  }
					  htmp=CreateWindowEx(WS_EX_CLIENTEDGE,"BUTTON","Apply",
						  WS_CHILD|WS_VISIBLE|BS_PUSHBUTTON|BS_PUSHLIKE,
						  10,70,30 * 9/4,30,hwnd,(HMENU)200,GetModuleHandle(NULL),NULL);
					  htmp=CreateWindowEx(WS_EX_CLIENTEDGE,"BUTTON","Close",
						  WS_CHILD|WS_VISIBLE|BS_PUSHBUTTON|BS_PUSHLIKE,
						  10,70,30 * 9/4,30,hwnd,(HMENU)201,GetModuleHandle(NULL),NULL);
					  style.numRows=8+1;
					  wsprintf(str,"Diag Created..%d \r\n",GetTickCount());
					  LoggerInsert(str);
	}
		break;
	case WM_USER+1:
		isDrawing=1;
		ShowWindow(hwnd,SW_SHOWDEFAULT);
		UpdateWindow(hwnd);
		break;
	case WM_MOVE:
	{
					if (isDrawing)
					{
						RECT rc;
						GetWindowRect(gData.hwnd,&rc);
						MoveWindow(hwnd,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
						ResizeDialogControls(hwnd);
						wsprintf(str,"WM_MOV_created.Inside.%d \r\n",GetTickCount());
						LoggerInsert(str);
					}
	}
		break;
	case WM_SIZE:
	{
					if (isDrawing)
					{
						RECT rc;
						GetWindowRect(gData.hwnd,&rc);
						MoveWindow(hwnd,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
						ResizeDialogControls(hwnd);
						wsprintf(str,"WM_SIZE_created.INsideside.%d \r\n",GetTickCount());
						LoggerInsert(str);
					}
	}
		break;
	case WM_SHOWWINDOW:
	{
						  if (wParam)
						  {
							  RECT rc;
							  GetWindowRect(gData.hwnd,&rc);
							  MoveWindow(hwnd,rc.right+20,rc.top,style.widthDialg,rc.bottom - rc.top,TRUE);
							  EnableButtons(NULL,&hDiag);
							  ResizeDialogControls(hwnd);
							  int values[8]={ gData.opt.maxKnotNum_Trend,
								  gData.opt.minTrendOrder,
								  gData.opt.maxTrendOrder,
								  gData.opt.minSepDist_Trend,
								  gData.opt.maxMoveStepSize,
								  gData.opt.burnin,
								  gData.opt.samples,
								  gData.opt.chainNumber,
							  };
							  for (int i=0; i < 8; i++)
							  {	  
								  SetDlgItemInt(hwnd,100+i,values[i],1);
							  }
							  ResizeDialogControls(hwnd);
						  }
						  else
						  {
							  gData.optStatus=NeedUpdate;
						  }
	}
		break;
	case WM_COMMAND:
		switch (LOWORD(wParam))
		{
		case 200:
		{
					 BOOL translated=0;
					 int i=100;
					 gData.opt.maxKnotNum_Trend=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.minTrendOrder=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.maxTrendOrder=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.minSepDist_Trend=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.maxMoveStepSize=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.burnin=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.samples=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.opt.chainNumber=GetDlgItemInt(hwnd,i,&translated,0); i++;
					 gData.optStatus=NeedUpdate;
		}
			break;
		case 201:
			ShowWindow(hwnd,SW_HIDE);
			break;
		default:
			break;
		}
		break;
	case WM_CLOSE:
		ShowWindow(hwnd,SW_HIDE);
		break;
	default:
		if (IsWindowUnicode(hwnd))
			return (DefWindowProc(hwnd,msg,wParam,lParam)); 
		else
			return (DefWindowProc(hwnd,msg,wParam,lParam)); 
	}
	return 0;
}
#else
static char fileID UNUSED_DECORATOR='c';
#endif
ENABLE_MANY_WARNINGS
