#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if defined(OS_WIN64) 

#include "abc_win32_demo.h"
void BEAST2_InitGlobalData(void) {
	style = (Style) {
			.hMargin		= 3.f,
			.rEdit			= 0.1f,
			.hgtButtonBar	= 30.f,
			.wButton        = 72.f,
			.sep            = 10.f,
			.vButtonRatio   = 0.7f,
			.vScrollRatio   = 0.5f,
			.labelGap       = 8.f,
			.rFig			= { 0.3, 0.7 / 2. * 2 / 3., 0.7 / 2. / 3, 0.7 / 2. * 2. / 3, 0.7 / 2. / 3 },
			.fractionLabel	= 0.6f,
			.fractionEdit	= 0.5f,
			.widthDialg		= 270.f,
	};
	
	memset(hBitmap,       0, sizeof(HBITMAP)* 5);
	memset(hBufferBitmap, 0, sizeof(HBITMAP)* 5);
	memset(memDC,		  0, sizeof(HDC)* 5);
	memset(bufferDC,     0, sizeof(HDC)* 5);

	blackBrush	= 0;
	greenPen	= 0;
	redPen		= 0;
	//--------------------------------
	memset(&gData, 0, sizeof(gData));
	InitializeCriticalSection(&gData.cs);
	InitializeConditionVariable(&gData.cv);

	gData.timerInterval = 25;
	gData.sleepInterval = 50;
	gData.status	= RUN;
	gData.optStatus = NotAssigned;
}
void BEAST2_AllocatePlotData(void)
{
	gData.plotData[0][0] = malloc(sizeof(int)*gData.N * 2); //Yfit_mean
	gData.plotData[0][1] = malloc(sizeof(int)*gData.N * 2); //Yfit_cur
	

	gData.plotData[1][0] = malloc(sizeof(int)*gData.N * 2); //S_mean
	gData.plotData[1][1] = malloc(sizeof(int)*gData.N * 2); //S_cur
	gData.plotData[1][2] = malloc(sizeof(int)*gData.N);     //S_bp
	gData.plotData[1][3] = malloc(sizeof(int)*gData.N);     //S_bp
	gData.plotData[1][4] = malloc(sizeof(int)*gData.N * 4); //S_CI

	gData.plotData[2][0] = malloc(sizeof(int)*gData.N * 2); //S_Porb
 

	gData.plotData[3][0] = malloc(sizeof(int)*gData.N * 2); //T_mean
	gData.plotData[3][1] = malloc(sizeof(int)*gData.N * 2); //T_cur
	gData.plotData[3][2] = malloc(sizeof(int)*gData.N);     //T_bp
	gData.plotData[3][3] = malloc(sizeof(int)*gData.N);     //T_bp
	gData.plotData[3][4] = malloc(sizeof(int)*gData.N * 4); //T_CI

	gData.plotData[4][0] = malloc(sizeof(int)*gData.N * 2); //T_Pprob
 
}
void BEAST2_GeneratePlotData(void)
{
	int N      = gData.N;
	int sample = gData.sample;


	F32 W, H, Ymax, Ymin, dX;
	F32 x;	

	F32 a;// = H + H / (Ymax - Ymin) *Ymin;
	F32 b;// = H / (Ymax - Ymin) / sample;  	//F32 b2;//= H / (Ymax - Ymin);
	I32PTR data;
	I32PTR ptsPerPoly;

	//gData.plotData[0][0] = malloc(sizeof(int)*gData.N * 2); //Yfit_mean
	//gData.plotData[0][1] = malloc(sizeof(int)*gData.N * 2); //Yfit_cur
	W = (F32) gData.w[0];
	H = (F32) gData.h[0];
	Ymax = gData.yMax;
	Ymin = gData.yMin;
	dX = W / N;
	//Figure 1
	if (gData.S != NULL) {

		// The average fitted curve
		data = gData.plotData[0][0];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * (gData.s[i] + gData.t[i]);

		// the current fitted curve
		data = gData.plotData[0][1];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin); 	x = 0;
		for (int i = 0; i < N; i++)	data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * (gData.curs[i] + gData.curt[i]);
	}
	else { // Trend-only 
		// The average fitted curve
		data = gData.plotData[0][0];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * (gData.t[i]);

		// the current fitted curve
		data = gData.plotData[0][1];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin); 	x = 0;
		for (int i = 0; i < N; i++)	data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * (gData.curt[i] );


	}

	  
	//gData.plotData[1][0] = malloc(sizeof(int)*gData.N * 2); //S_mean
	//gData.plotData[1][1] = malloc(sizeof(int)*gData.N * 2); //S_cur
	//gData.plotData[1][2] = malloc(sizeof(int)*gData.N);     //S_bp
	//gData.plotData[1][3] = malloc(sizeof(int)*gData.N);     //S_bp
	//gData.plotData[1][4] = malloc(sizeof(int)*gData.N * 4); //S_CI
	 //Figure 2
	if (gData.S != NULL) {
		W = gData.w[1];
		H = gData.h[1];
		Ymax = gData.yMaxS;
		Ymin = gData.yMinS;
		dX = W / N;

		// the averaged season curve
		data = gData.plotData[1][0];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.s[i];

		// the current season curve
		data = gData.plotData[1][1];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin); 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.curs[i];

		// Season breakpoints

		ptsPerPoly = gData.plotData[1][3];
		data = gData.plotData[1][2];
		for (int i = 0; i < gData.sKnotNum; i++) {
			data[4 * i] = dX * gData.S[i];
			data[4 * i + 1] = 0;
			data[4 * i + 2] = dX * gData.S[i];
			data[4 * i + 3] = H;
			ptsPerPoly[i] = 2;
		}

		// the season CI envelop
		data = gData.plotData[1][4];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin); 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.sCI[i];

		F32PTR Y = gData.sCI + 2 * N - 1;
		data = data + 2 * N;
		x = dX * (N - 1);
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x -= dX, data[2 * i + 1] = a - b * Y[-i];


		//---------Figure 3------------

		W = gData.w[2];
		H = gData.h[2];
		Ymax = 0.9;
		Ymin = -0.02;
		dX = W / N;


		// the seasnn prob curve
		data = gData.plotData[2][0];
		a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
		for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.sProb[i];
	}



	//--------------Figure 4

	//gData.plotData[1][0] = malloc(sizeof(int)*gData.N * 2); //S_mean
	//gData.plotData[1][1] = malloc(sizeof(int)*gData.N * 2); //S_cur
	//gData.plotData[1][2] = malloc(sizeof(int)*gData.N);     //S_bp
	//gData.plotData[1][3] = malloc(sizeof(int)*gData.N);     //S_bp
	//gData.plotData[1][4] = malloc(sizeof(int)*gData.N * 4); //S_CI

	W = gData.w[3];
	H = gData.h[3];
	Ymax = gData.yMaxT;
	Ymin = gData.yMinT;
	dX   = W / N;

	// the aveagre tredn
	data = gData.plotData[3][0];
	a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
	for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.t[i];

	// the current trend
	data = gData.plotData[3][1];
	a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) ; 	x = 0;
	for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.curt[i];

	
	// Trend breakpoints		
	data       = gData.plotData[3][2];
	ptsPerPoly = gData.plotData[3][3];
	for (int i = 0; i < gData.tKnotNum; i++) {
		data[4 * i]     = dX * gData.T[i];
		data[4 * i + 1] = 0;
		data[4 * i + 2] = dX * gData.T[i];
		data[4 * i + 3] = H;
		ptsPerPoly[i] = 2;
	}


	// the trend CI envelop
	data = gData.plotData[3][4];
	a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin); 	x = 0;
	for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.tCI[i];

	F32PTR Y = gData.tCI + 2 * N - 1;
	data = data + 2 * N;
	x    = dX * (N - 1);
	for (I32 i = 0; i < N; i++)  data[2 * i] = x, x -= dX, data[2 * i + 1] = a - b * Y[-i];

 	//---------Figure 5-----------------

	W = gData.w[4];
	H = gData.h[4];
	Ymax = 0.9;
	Ymin = -0.02;
	dX = W / N;
 
	// the seasnn prob curve
	data = gData.plotData[4][0];
	a = H + H / (Ymax - Ymin) * Ymin; 	b = H / (Ymax - Ymin) / sample; 	x = 0;
	for (I32 i = 0; i < N; i++)  data[2 * i] = x, x += dX, data[2 * i + 1] = a - b * gData.tProb[i];

	
}
void BEAST2_DrawPlots(HDC hdc)
{
	
	int  W, H;
	int  mode;
	HDC  hdcMEM;
	RECT rect;
	char str[80];
	

	//FillRect(memDC, &rcClient, (HBRUSH)COLOR_WINDOW + 2);
	//ExtTextOut(memDC, 0, 0, ETO_OPAQUE, &rcClient, "", 0, 0);
	//FillRect(memDC, &rcClient, blackBrush);

	//Figure 1
	W = gData.w[1-1];
	H = gData.h[1-1];
	SetRect(&rect, 0, 0, W, H);
	hdcMEM = memDC[0];

	BitBlt(hdcMEM, 0, 0, W, H, bufferDC[0], 0, 0, SRCCOPY);

	Polyline(hdcMEM, (POINT *)gData.plotData[0][0], gData.N); // Yfit_mean

	SelectObject(hdcMEM, yellowPen);
	Polyline(hdcMEM, (POINT *)gData.plotData[0][1], gData.N); //current Trend

	SelectObject(hdcMEM, greenPen);
	wsprintf(str, "iteration:%d Data+Fitted", gData.ite);
	TextOut(hdcMEM, 0, 0, str, lstrlen(str));

	BitBlt(hdc, gData.x0[0], gData.y0[0], W, H, hdcMEM, 0, 0, SRCCOPY);

	//Figure 2
	W = gData.w[2 - 1];
	H = gData.h[2 - 1];
	hdcMEM = memDC[1];

	if (gData.S != NULL) {
		SetRect(&rect, 0, 0, W, H);
		FillRect(hdcMEM, &rect, blackBrush);

		mode = SetROP2(hdcMEM, R2_MERGEPEN); //R2_MERGENOTPEN
		SelectObject(hdcMEM, grayBrush);
		Polygon(hdcMEM, (POINT*)gData.plotData[1][4], gData.N * 2);
		SetROP2(hdcMEM, mode);

		Polyline(hdcMEM, (POINT*)gData.plotData[1][0], gData.N); // avaraged trend

		SelectObject(hdcMEM, yellowPen);
		Polyline(hdcMEM, (POINT*)gData.plotData[1][1], gData.N);//current Trend
		PolyPolyline(hdcMEM, (POINT*)gData.plotData[1][2], (DWORD*)gData.plotData[1][3], gData.sKnotNum);//curent breakpoints

		SelectObject(hdcMEM, greenPen);
		wsprintf(str, "Season: chain#%d|iteration#%d|#scp%d", gData.curChainNumber, gData.ite, gData.sKnotNum);
		TextOut(hdcMEM, 0, 0, str, lstrlen(str));
		BitBlt(hdc, gData.x0[1], gData.y0[1], W, H, hdcMEM, 0, 0, SRCCOPY);
	}
	else {
		SetRect(&rect, 0, 0, W, H);
		FillRect(hdcMEM, &rect, blackBrush);

		SelectObject(hdcMEM, greenPen);
		wsprintf(str, "%s","No season component for the trend-only version");
		TextOut(hdcMEM, W/2, H/2, str, lstrlen(str));
		BitBlt(hdc, gData.x0[1], gData.y0[1], W, H, hdcMEM, 0, 0, SRCCOPY);
	}

	//-FIGURE3=======
	W = gData.w[2];
	H = gData.h[2];
	hdcMEM = memDC[2];
	if (gData.S != NULL) {
		SetRect(&rect, 0, 0, W, H);
		FillRect(hdcMEM, &rect, blackBrush);

		SelectObject(hdcMEM, yellowPen);
		Polyline(hdcMEM, (POINT*)gData.plotData[2][0], gData.N); // tProb						 
		//SelectObject(memDC[0], greenPen);

		wsprintf(str, "iteration:%d Probability of scp", gData.ite);
		TextOut(hdcMEM, 0, 0, str, lstrlen(str));

		BitBlt(hdc, gData.x0[2], gData.y0[2], W, H, hdcMEM, 0, 0, SRCCOPY);
	}
	else {	
		SetRect(&rect, 0, 0, W, H);
		FillRect(hdcMEM, &rect, blackBrush);

		SelectObject(hdcMEM, greenPen);
		wsprintf(str, "%s", "No season component for the trend-only version");
		TextOut(hdcMEM, W / 2, H / 2, str, lstrlen(str));
		BitBlt(hdc, gData.x0[1], gData.y0[1], W, H, hdcMEM, 0, 0, SRCCOPY);
	
	}

	//Figure 4
	W = gData.w[4 - 1];
	H = gData.h[4 - 1];
	hdcMEM = memDC[3];

	SetRect(&rect, 0, 0, W, H);
	FillRect(hdcMEM, &rect, blackBrush);

	mode = SetROP2(hdcMEM, R2_MERGEPEN); //R2_MERGENOTPEN
	SelectObject(hdcMEM, grayBrush);
	Polygon(hdcMEM, (POINT *)gData.plotData[3][4], gData.N * 2);
	SetROP2(hdcMEM, mode);

	Polyline(hdcMEM, (POINT *)gData.plotData[3][0], gData.N); // avaraged trend

	SelectObject(hdcMEM, yellowPen);
	Polyline(hdcMEM, (POINT *)gData.plotData[3][1], gData.N);//current Trend

	PolyPolyline(hdcMEM, (POINT *)gData.plotData[3][2], (DWORD*)gData.plotData[3][3], gData.tKnotNum);//curent breakpoints

	SelectObject(hdcMEM, greenPen);
	wsprintf(str, "Trend: chain#%d|iteration#%d|#tcp%d", gData.curChainNumber,gData.ite,gData.tKnotNum);
	TextOut(hdcMEM, 0, 0, str, lstrlen(str));
	BitBlt(hdc, gData.x0[3], gData.y0[3], W, H, hdcMEM, 0, 0, SRCCOPY);

	//-FIGURE5=======
	W = gData.w[4];
	H = gData.h[4];
	hdcMEM = memDC[4];
	SetRect(&rect, 0, 0, W, H);
	FillRect(hdcMEM, &rect, blackBrush);

	SelectObject(hdcMEM, yellowPen);
	Polyline(hdcMEM, (POINT *)gData.plotData[4][0], gData.N); // tProb						 
	//SelectObject(memDC[0], greenPen);
	wsprintf(str, "iteration:%d Probability of tcp", gData.ite);
	TextOut(hdcMEM, 0, 0, str, lstrlen(str));
	BitBlt(hdc, gData.x0[4], gData.y0[4], W, H, hdcMEM, 0, 0, SRCCOPY);
	 
}

#else
static char fileID UNUSED_DECORATOR = 'c';
#endif
 
#include "abc_000_warning.h"