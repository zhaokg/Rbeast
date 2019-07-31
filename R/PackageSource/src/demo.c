#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
#if defined(WIN64_OS) 
#include "demo.h"
GlobalStruct gData;
LParam threadParam;
Logger logger;
Style style;
HICON   hIcon;
HWND    hButton[4],hScroll,hEdit,hStatic[2],hDiag;
HBITMAP hBitmap[5]={ 0,0 },hBufferBitmap[5]={ 0,0 };
HDC     memDC[5]={ 0,0 },bufferDC[5]={ 0,0 };
HBRUSH  blackBrush=0,grayBrush=0;
HPEN	greenPen=0,bluePen=0,yellowPen=0,redPen=0;
void InitGlobalData()
{
	 style=(Style) {
		.hMargin=3,
		.rEdit=0.25,
		.hgtButtonBar=40,
		.wButton=72,
		.sep=10,
		.vButtonRatio=0.7,
		.vScrollRatio=0.5,
		.labelGap=8,
		.rFig={ 0.6,0.4 },
		.fractionLabel=0.5,
		.fractionEdit=0.5,
		.widthDialg=270,
	};
	 memset(hBitmap,0,sizeof(HBITMAP)* 2);
	 memset(hBufferBitmap,0,sizeof(HBITMAP)* 2);
	 memset(memDC,0,sizeof(HDC)* 2);
	 memset(bufferDC,0,sizeof(HDC)* 2);
	 blackBrush=0;
	 greenPen=0,redPen=0;
	 memset(&gData,0,sizeof(gData));
	 InitializeCriticalSection(&gData.cs);
	 InitializeConditionVariable(&gData.cv);
	 gData.timerInterval=25;
	 gData.sleepInterval=50;
	 gData.status=RUN;
	 gData.optStatus=NotAssigned;
}
void AllocatePlotData()
{
	gData.plotData[0][0]=malloc(sizeof(int)*gData.N * 2);
	gData.plotData[0][1]=malloc(sizeof(int)*gData.N * 2);
	gData.plotData[0][2]=malloc(sizeof(int)*gData.N);
	gData.plotData[0][3]=malloc(sizeof(int)*gData.N);
	gData.plotData[0][4]=malloc(sizeof(int)*gData.N * 4);
	gData.plotData[1][0]=malloc(sizeof(int)*gData.N * 2);
}
void EnableButtons(HWND *hButton,HWND *hDiag)
{
	byte isDone=(gData.status==DONE);
	if (hButton !=NULL)
	{
		EnableWindow(hButton[0],gData.status !=RUN||gData.status==DONE);
		EnableWindow(hButton[1],(gData.status==RUN) && (gData.status !=DONE));
		EnableWindow(hButton[2],isDone);
		EnableWindow(hButton[3],TRUE);
		if ((gData.status==RUN)||(gData.status==PAUSE))
			SetWindowText(hButton[0],"Resume");
		if (gData.status==DONE)
			SetWindowText(hButton[0],"Re-run");
	}
	if (hDiag !=NULL)
	{
		EnableWindow(hDiag[0],isDone);
		for (int i=0; i < 8; i++)
		{
			HWND hwnd=GetDlgItem(hDiag[0],i);
			EnableWindow(hwnd,isDone);
			hwnd=GetDlgItem(hDiag[0],i+100);
			EnableWindow(hwnd,isDone);
		}
		EnableWindow(GetDlgItem(hDiag[0],200),isDone);
		EnableWindow(GetDlgItem(hDiag[0],201),isDone);
	}
}
void ResizeControls(int x,int y,int N)
{
	float sep=style.sep;
	float mgr=style.hMargin;
	float width,height,x0,y0,xExt,yExt;
	width=style.wButton;	
	yExt=style.hgtButtonBar;
	height=yExt*style.vButtonRatio;
	y0=(yExt -height )/2;	
	x0=mgr;
	for (int i=0; i < 4; i++)	{		
		SetWindowPos(hButton[i],NULL,x0,y0,width,height,SWP_NOZORDER);
		x0=x0+(width+sep);
	}
	x0=x0 - sep;
	xExt=x - x0 -mgr;
	int   len;
	len=3 * style.xChar+style.xCapChar;
	width=xExt - 2 * (len+2 * style.labelGap);
	height=yExt*style.vScrollRatio;
	y0=(yExt - height)/2;
	x0=x0+(len+2 * style.labelGap);		
	SetWindowPos(hScroll,NULL,x0,y0,width,height,SWP_NOZORDER);
	x0=mgr+4 * style.wButton+3 * sep+style.labelGap;
	height=style.yChar;
	y0=(yExt - height)/2;
	width=len;
	SetWindowPos(hStatic[0],NULL,x0,y0,width,height,SWP_NOZORDER);
	x0=x - mgr - (len+style.labelGap);
    SetWindowPos(hStatic[1],NULL,x0,y0,width,height,SWP_NOZORDER);
	x0=mgr;
	y0=y*(1 - style.rEdit);
	height=y* style.rEdit;
	width=x - 2 * mgr;
	SetWindowPos(hEdit,NULL,x0,y0,width,height,SWP_NOZORDER);
	int vGap=10;
	width=x - 2 * mgr;
	yExt=y - y*style.rEdit - style.hgtButtonBar - (N+1)*vGap;
	y0=style.hgtButtonBar+vGap;	
	for (int i=0; i < N; i++)
	{
		gData.w[i]=width;
		gData.h[i]=yExt * style.rFig[i];
		gData.x0[i]=mgr;
		gData.y0[i]=y0;
		y0=y0+gData.h[i]+vGap;
	} 
}
void CreateGDIObject(HWND hwnd,int N)
{
	HDC  hdc;
	RECT rect;
	GetClientRect(hwnd,&rect);
	ResizeControls(rect.right-rect.left,rect.bottom - rect.top,N);
	if (blackBrush==0)
	{
		blackBrush=CreateSolidBrush(RGB(30,33,40));
		grayBrush=CreateSolidBrush(RGB(125,125,125));
		greenPen=CreatePen(PS_SOLID,1,RGB(0,205,1));
		bluePen=CreatePen(PS_SOLID,1,RGB(0,5,201));
		redPen=CreatePen(PS_SOLID,1,RGB(181,45,69));
		yellowPen=CreatePen(PS_SOLID,2,RGB(255,255,0));
	}
	hdc=GetDC(hwnd);
	for (int i=0; i < N; i++)
	{
		if (hBitmap[i] !=0)
			DeleteDC(memDC[i]),DeleteObject(hBitmap[i]);
		hBitmap[i]=CreateCompatibleBitmap(hdc,gData.w[i],gData.h[i]);
		memDC[i]=CreateCompatibleDC(hdc);
		SelectObject(memDC[i],hBitmap[i]);
		SelectObject(memDC[i],greenPen);
		SelectObject(memDC[i],blackBrush);
		SetBkMode(memDC[i],TRANSPARENT);
		SetTextColor(memDC[i],RGB(155,155,155));
		if (i==0)
		{
			if (hBufferBitmap[i] !=0)
			{
				DeleteDC(bufferDC[i]);
				DeleteObject(hBufferBitmap[i]);
			}
			hBufferBitmap[i]=CreateCompatibleBitmap(hdc,gData.w[i],gData.h[i]);
			bufferDC[i]=CreateCompatibleDC(hdc);
			SelectObject(bufferDC[i],hBufferBitmap[i]);
			SelectObject(bufferDC[i],redPen);
			SelectObject(bufferDC[i],blackBrush);
		}
	}	
	ReleaseDC(hwnd,hdc);
}
void DeleteGDIObject( int N)
{
	for (int i=0; i < N; i++)
	{
		DeleteDC(memDC[i]);
		DeleteObject(hBitmap[i]);
		if (i==0)
		{
			DeleteDC(bufferDC[i]);
			DeleteObject(hBufferBitmap[i]);
		}
	}
	DeleteObject(blackBrush);
	DeleteObject(grayBrush);
	DeleteObject(greenPen);
	DeleteObject(bluePen);
	DeleteObject(redPen);
	DeleteObject(yellowPen);
}
void GeneratePlotData(  )
{
	float W=gData.w[0];
	float H=gData.h[0];
	int N=gData.N;
	int sample=gData.sample;
	float Ymax=gData.yMax;
	float Ymin=gData.yMin;
	float a=H+H/(Ymax - Ymin) *Ymin;
	float b1=H/(Ymax - Ymin)/sample;
	float b2=H/(Ymax - Ymin);
	float c=W/N;
	float tmpX=0;
	int * _restrict   data1,*_restrict    data2;
	F32PTR Y1,Y2;
	data1=gData.plotData[0][0];
	data2=gData.plotData[0][1];
	Y1=gData.t;
	Y2=gData.ct;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX; 
		data1[2 * i+1]=a - b1*Y1[i];
		data2[2 * i]=tmpX; 
		data2[2 * i+1]=a - b2*Y2[i];
		tmpX+=c;
	}
	data1=gData.plotData[0][2];
	data2=gData.plotData[0][3];
		for (int i=0; i < gData.tKnotNum; i++)
	{
		data1[4 * i]=c * gData.T[i];
		data1[4 * i+1]=0;
		data1[4 * i+2]=c * gData.T[i];
		data1[4 * i+3]=H; 
		data2[i]=2;
	}
	data1=gData.plotData[0][4]; 
	Y2=gData.tCI;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2++);
		tmpX+=c;
	}
	tmpX=c*(N - 1);
	Y2=gData.tCI+2*N-1;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2--);
		tmpX -=c;
	}
	Ymax=0.9;
	Ymin=-0.02;
	W=gData.w[1];
	H=gData.h[1];
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;	
	c=W/N;
	data1=gData.plotData[1][0];	
	int32_t * _restrict Yi=gData.tProb;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*Yi[i];		
		tmpX+=c;
	}
}
void DrawPlots(HDC hdc)
{
	RECT rect;
	int W=gData.w[0];
	int H=gData.h[0];
	SetRect(&rect,0,0,W,H);
	char str[80];
	BitBlt(memDC[0],0,0,W,H,bufferDC[0],0,0,SRCCOPY);
	int mode=SetROP2(memDC[0],R2_MERGEPEN); 
	SelectObject(memDC[0],grayBrush);
	Polygon(memDC[0],(POINT *)gData.plotData[0][4],gData.N * 2);
	SetROP2(memDC[0],mode);
	Polyline(memDC[0],(POINT *)gData.plotData[0][0],gData.N); 
	SelectObject(memDC[0],yellowPen);
	Polyline(memDC[0],(POINT *)gData.plotData[0][1],gData.N);
	PolyPolyline(memDC[0],(POINT *)gData.plotData[0][2],(DWORD*)gData.plotData[0][3],gData.tKnotNum);
	SelectObject(memDC[0],greenPen);
	wsprintf(str,"iteration:%d data+fitted",gData.ite);
	TextOut(memDC[0],0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[0],gData.y0[0],W,H,memDC[0],0,0,SRCCOPY);
	W=gData.w[1];
	H=gData.h[1];
	SetRect(&rect,0,0,W,H);
	FillRect(memDC[1],&rect,blackBrush);
	SelectObject(memDC[1],yellowPen);
	Polyline(memDC[1],(POINT *)gData.plotData[1][0],gData.N); 
	wsprintf(str,"iteration:%d Probability",gData.ite);
	TextOut(memDC[1],0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[1],gData.y0[1],W,H,memDC[1],0,0,SRCCOPY);
}
void ResizeDialogControls(HWND hwnd)
{
	RECT rc;
	GetClientRect(hwnd,&rc);
	int x=rc.right - rc.left;
	int y=rc.bottom - rc.top;
	int rowHeight=y/style.numRows;
	int lastRowHeight=y - rowHeight*(style.numRows - 1);
	float width,height,x0,y0,xExt,yExt;
	float labelWidth=style.xChar*(lstrlen("maxKnotNum_Trend")+2);
	for (int i=0; i < 8; i++)	{
		HWND hc;
		hc=GetDlgItem(hwnd,i);
		x0=5;
		y0=i*rowHeight+rowHeight*(1 - style.fractionLabel)/2;
		width=labelWidth;
		height=rowHeight*style.fractionLabel;
		SetWindowPos(hc,NULL,x0,y0,width,height,SWP_NOZORDER);
		hc=GetDlgItem(hwnd,i+100);
		x0=labelWidth+30;
		y0=i*rowHeight+rowHeight*(1 - style.fractionEdit)/2;
		width=x-x0-5;
		height=rowHeight*style.fractionEdit;
		SetWindowPos(hc,NULL,x0,y0,width,height,SWP_NOZORDER);
	}
	{
		HWND hc;
		hc=GetDlgItem(hwnd,200);
		float r=0.25;
		width=r*x;
		height=lastRowHeight*style.fractionEdit;
		y0=y-lastRowHeight+lastRowHeight*(1 - style.fractionEdit)/2;
		x0=x*0.25/2;
		SetWindowPos(hc,NULL,x0,y0,width,height,SWP_NOZORDER);
		hc=GetDlgItem(hwnd,201);				
		x0=x*(0.5+0.25/2);
		SetWindowPos(hc,NULL,x0,y0,width,height,SWP_NOZORDER);
	}
}
void LoggerInsert(char *newMsg)
{
	int lenMsg=lstrlen(newMsg);
	if ((lenMsg+logger.len) >=logger.LEN)
	{
		char *tmp;
		int   newBytes;
		newBytes=max(lenMsg,1000);
		tmp=malloc(logger.LEN+newBytes+1);
		memcpy(tmp,logger.str,logger.len+1);
		free(logger.str);
		logger.str=tmp;
		logger.LEN=logger.LEN+newBytes;
	}
	memcpy(logger.str+logger.len,newMsg,lenMsg);
	logger.len+=lenMsg;
	*(logger.str+(logger.len+1) - 1)=0;
	SetWindowText(hEdit,logger.str);
	SendMessage(hEdit,EM_LINESCROLL,0,100);
}
HWND GetConsoleHwnd(void)
{ 
#define MY_BUFSIZE 1024 
	HWND hwndFound;         
	char pszNewWindowTitle[MY_BUFSIZE]; 
	char pszOldWindowTitle[MY_BUFSIZE]; 
	GetConsoleTitle(pszOldWindowTitle,MY_BUFSIZE);
	wsprintf(pszNewWindowTitle,"%d/%d",GetTickCount(),GetCurrentProcessId());
	SetConsoleTitle(pszNewWindowTitle);
	Sleep(20);
	hwndFound=FindWindow(NULL,pszNewWindowTitle);
	SetConsoleTitle(pszOldWindowTitle);
	return(hwndFound);
}
#if defined(MSVC_COMPILER)
void SetupConsole()
{
	BOOL bCreated=AllocConsole();
	HANDLE hConOut;
	if (!bCreated)
		return; 
	hConOut=GetStdHandle(STD_OUTPUT_HANDLE);
	int fd=_open_osfhandle((intptr_t)hConOut,0);
	_dup2(fd,1);
	HWND hwnd=GetConsoleHwnd();
	MoveWindow(hwnd,500,500,800,600,TRUE);
}
#endif
BYTE ANDmaskIcon[]={ 0xFF,0xFF,0xFF,0xFF,   
0xFF,0xFF,0xC3,0xFF,   
0xFF,0xFF,0x00,0xFF,   
0xFF,0xFE,0x00,0x7F,   
0xFF,0xFC,0x00,0x1F,   
0xFF,0xF8,0x00,0x0F,   
0xFF,0xF8,0x00,0x0F,   
0xFF,0xF0,0x00,0x07,   
0xFF,0xF0,0x00,0x03,   
0xFF,0xE0,0x00,0x03,   
0xFF,0xE0,0x00,0x01,   
0xFF,0xE0,0x00,0x01,   
0xFF,0xF0,0x00,0x01,   
0xFF,0xF0,0x00,0x00,   
0xFF,0xF8,0x00,0x00,   
0xFF,0xFC,0x00,0x00,   
0xFF,0xFF,0x00,0x00,   
0xFF,0xFF,0x80,0x00,   
0xFF,0xFF,0xE0,0x00,   
0xFF,0xFF,0xE0,0x01,   
0xFF,0xFF,0xF0,0x01,   
0xFF,0xFF,0xF0,0x01,   
0xFF,0xFF,0xF0,0x03,   
0xFF,0xFF,0xE0,0x03,   
0xFF,0xFF,0xE0,0x07,   
0xFF,0xFF,0xC0,0x0F,   
0xFF,0xFF,0xC0,0x0F,   
0xFF,0xFF,0x80,0x1F,   
0xFF,0xFF,0x00,0x7F,   
0xFF,0xFC,0x00,0xFF,   
0xFF,0xF8,0x03,0xFF,   
0xFF,0xFC,0x3F,0xFF };  
BYTE XORmaskIcon[]={ 0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x38,0x00,   
0x00,0x00,0x7C,0x00,   
0x00,0x00,0x7C,0x00,   
0x00,0x00,0x7C,0x00,   
0x00,0x00,0x38,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00,   
0x00,0x00,0x00,0x00 };  
void InitGlobalData_ST()
{
	style=(Style) {
		.hMargin=3,
			.rEdit=0.1,
			.hgtButtonBar=30,
			.wButton=72,
			.sep=10,
			.vButtonRatio=0.7,
			.vScrollRatio=0.5,
			.labelGap=8,
			.rFig={ 0.3,0.7/2. * 2/3.,0.7/2./3,0.7/2. * 2./3,0.7/2./3 },
			.fractionLabel=0.5,
			.fractionEdit=0.5,
			.widthDialg=270,
	};
	memset(hBitmap,0,sizeof(HBITMAP)* 5);
	memset(hBufferBitmap,0,sizeof(HBITMAP)* 5);
	memset(memDC,0,sizeof(HDC)* 5);
	memset(bufferDC,0,sizeof(HDC)* 5);
	blackBrush=0;
	greenPen=0,redPen=0;
	memset(&gData,0,sizeof(gData));
	InitializeCriticalSection(&gData.cs);
	InitializeConditionVariable(&gData.cv);
	gData.timerInterval=25;
	gData.sleepInterval=50;
	gData.status=RUN;
	gData.optStatus=NotAssigned;
}
void AllocatePlotData_ST()
{
	gData.plotData[0][0]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[0][1]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[1][0]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[1][1]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[1][2]=malloc(sizeof(int)*gData.N);     
	gData.plotData[1][3]=malloc(sizeof(int)*gData.N);     
	gData.plotData[1][4]=malloc(sizeof(int)*gData.N * 4); 
	gData.plotData[2][0]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[3][0]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[3][1]=malloc(sizeof(int)*gData.N * 2); 
	gData.plotData[3][2]=malloc(sizeof(int)*gData.N);     
	gData.plotData[3][3]=malloc(sizeof(int)*gData.N);     
	gData.plotData[3][4]=malloc(sizeof(int)*gData.N * 4); 
	gData.plotData[4][0]=malloc(sizeof(int)*gData.N * 2); 
}
void GeneratePlotData_ST()
{
	float W;  
	float H;  
	int N=gData.N;
	int sample=gData.sample;
	float Ymax  ;
	float Ymin  ;
	float a;
	float b1;
	float b2;
	float c;
	float tmpX;
	int * _restrict   data1,*_restrict    data2;
	F32PTR Y1,Y2;
	W=gData.w[0];
	H=gData.h[0];
	Ymax=gData.yMax;
	Ymin=gData.yMin;
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;
	b2=H/(Ymax - Ymin);
	c=W/N;
	tmpX=0;
	data1=gData.plotData[0][0];
	data2=gData.plotData[0][1];
	Y1=gData.s;
	Y2=gData.curs;
	F32PTR Y3,Y4;
	Y3=gData.t;
	Y4=gData.ct;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*(Y1[i]+Y3[i]);
		data2[2 * i]=tmpX;
		data2[2 * i+1]=a - b2*(Y2[i]+Y4[i]);
		tmpX+=c;
	}
	W=gData.w[1];
	H=gData.h[1];
	Ymax=gData.yMaxS;
	Ymin=gData.yMinS;
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;
	b2=H/(Ymax - Ymin);
	c=W/N;
	tmpX=0;
	data1=gData.plotData[1][0];
	data2=gData.plotData[1][1];
	Y1=gData.s;
	Y2=gData.curs;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*Y1[i];
		data2[2 * i]=tmpX;
		data2[2 * i+1]=a - b2*Y2[i];
		tmpX+=c;
	}
	data1=gData.plotData[1][2];
	data2=gData.plotData[1][3];
	for (int i=0; i < gData.sKnotNum; i++)
	{
		data1[4 * i]=c * gData.S[i];
		data1[4 * i+1]=0;
		data1[4 * i+2]=c * gData.S[i];
		data1[4 * i+3]=H;
		data2[i]=2;
	}
	data1=gData.plotData[1][4];
	Y2=gData.sCI;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2++);
		tmpX+=c;
	}
	tmpX=c*(N - 1);
	Y2=gData.sCI+2 * N - 1;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2--);
		tmpX -=c;
	}
	W=gData.w[2];
	H=gData.h[2];
	Ymax=0.9;
	Ymin=-0.02;
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;
	c=W/N;
	data1=gData.plotData[2][0];
	int32_t * _restrict Yi=gData.sProb;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*Yi[i];
		tmpX+=c;
	}
	W=gData.w[3];
	H=gData.h[3];
	Ymax=gData.yMaxT;
	Ymin=gData.yMinT;
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;
	b2=H/(Ymax - Ymin);
	c=W/N;
	tmpX=0;
	data1=gData.plotData[3][0];
	data2=gData.plotData[3][1];
	Y1=gData.t;
	Y2=gData.ct;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*Y1[i];
		data2[2 * i]=tmpX;
		data2[2 * i+1]=a - b2*Y2[i];
		tmpX+=c;
	}
	data1=gData.plotData[3][2];
	data2=gData.plotData[3][3];
	for (int i=0; i < gData.tKnotNum; i++)
	{
		data1[4 * i]=c * gData.T[i];
		data1[4 * i+1]=0;
		data1[4 * i+2]=c * gData.T[i];
		data1[4 * i+3]=H;
		data2[i]=2;
	}
	data1=gData.plotData[3][4];
	Y2=gData.tCI;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2++);
		tmpX+=c;
	}
	tmpX=c*(N - 1);
	Y2=gData.tCI+2 * N - 1;
	for (int i=0; i < N; i++)
	{
		*data1++=tmpX;
		*data1++=a - b2*(*Y2--);
		tmpX -=c;
	}
	W=gData.w[4];
	H=gData.h[4];
	Ymax=0.9;
	Ymin=-0.02;
	a=H+H/(Ymax - Ymin) *Ymin;
	b1=H/(Ymax - Ymin)/sample;
	c=W/N;
	data1=gData.plotData[4][0];
	Yi=gData.tProb;
	tmpX=0;
	for (int i=0; i < N; i++)
	{
		data1[2 * i]=tmpX;
		data1[2 * i+1]=a - b1*Yi[i];
		tmpX+=c;
	}
}
void DrawPlots_ST(HDC hdc)
{
	RECT rect;
	int W; 
	int H; 
	int mode;
	HDC hdcMEM;
	char str[80];
	W=gData.w[1-1];
	H=gData.h[1-1];
	SetRect(&rect,0,0,W,H);
	hdcMEM=memDC[0];
	BitBlt(hdcMEM,0,0,W,H,bufferDC[0],0,0,SRCCOPY);
	Polyline(hdcMEM,(POINT *)gData.plotData[0][0],gData.N); 
	SelectObject(hdcMEM,yellowPen);
	Polyline(hdcMEM,(POINT *)gData.plotData[0][1],gData.N); 
	SelectObject(hdcMEM,greenPen);
	wsprintf(str,"iteration:%d Data+Fitted",gData.ite);
	TextOut(hdcMEM,0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[0],gData.y0[0],W,H,hdcMEM,0,0,SRCCOPY);
	W=gData.w[2 - 1];
	H=gData.h[2 - 1];
	hdcMEM=memDC[1];
	SetRect(&rect,0,0,W,H);
	FillRect(hdcMEM,&rect,blackBrush);
	mode=SetROP2(hdcMEM,R2_MERGEPEN); 
	SelectObject(hdcMEM,grayBrush);
	Polygon(hdcMEM,(POINT *)gData.plotData[1][4],gData.N * 2);
	SetROP2(hdcMEM,mode);
	Polyline(hdcMEM,(POINT *)gData.plotData[1][0],gData.N); 
	SelectObject(hdcMEM,yellowPen);
	Polyline(hdcMEM,(POINT *)gData.plotData[1][1],gData.N);
	PolyPolyline(hdcMEM,(POINT *)gData.plotData[1][2],(DWORD*)gData.plotData[1][3],gData.sKnotNum);
	SelectObject(hdcMEM,greenPen);
	wsprintf(str,"iteration:%d Season",gData.ite);
	TextOut(hdcMEM,0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[1],gData.y0[1],W,H,hdcMEM,0,0,SRCCOPY);
	W=gData.w[2];
	H=gData.h[2];
	hdcMEM=memDC[2];
	SetRect(&rect,0,0,W,H);
	FillRect(hdcMEM,&rect,blackBrush);
	SelectObject(hdcMEM,yellowPen);
	Polyline(hdcMEM,(POINT *)gData.plotData[2][0],gData.N); 
	wsprintf(str,"iteration:%d Probability of scp",gData.ite);
	TextOut(hdcMEM,0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[2],gData.y0[2],W,H,hdcMEM,0,0,SRCCOPY);
	W=gData.w[4 - 1];
	H=gData.h[4 - 1];
	hdcMEM=memDC[3];
	SetRect(&rect,0,0,W,H);
	FillRect(hdcMEM,&rect,blackBrush);
	mode=SetROP2(hdcMEM,R2_MERGEPEN); 
	SelectObject(hdcMEM,grayBrush);
	Polygon(hdcMEM,(POINT *)gData.plotData[3][4],gData.N * 2);
	SetROP2(hdcMEM,mode);
	Polyline(hdcMEM,(POINT *)gData.plotData[3][0],gData.N); 
	SelectObject(hdcMEM,yellowPen);
	Polyline(hdcMEM,(POINT *)gData.plotData[3][1],gData.N);
	PolyPolyline(hdcMEM,(POINT *)gData.plotData[3][2],(DWORD*)gData.plotData[3][3],gData.tKnotNum);
	SelectObject(hdcMEM,greenPen);
	wsprintf(str,"iteration:%d Trend",gData.ite);
	TextOut(hdcMEM,0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[3],gData.y0[3],W,H,hdcMEM,0,0,SRCCOPY);
	W=gData.w[4];
	H=gData.h[4];
	hdcMEM=memDC[4];
	SetRect(&rect,0,0,W,H);
	FillRect(hdcMEM,&rect,blackBrush);
	SelectObject(hdcMEM,yellowPen);
	Polyline(hdcMEM,(POINT *)gData.plotData[4][0],gData.N); 
	wsprintf(str,"iteration:%d Probability of tcp",gData.ite);
	TextOut(hdcMEM,0,0,str,lstrlen(str));
	BitBlt(hdc,gData.x0[4],gData.y0[4],W,H,hdcMEM,0,0,SRCCOPY);
}
#else
static char fileID UNUSED_DECORATOR='c';
#endif
ENABLE_MANY_WARNINGS
