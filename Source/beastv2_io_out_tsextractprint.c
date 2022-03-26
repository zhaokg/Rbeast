#include "abc_000_warning.h"

#include "abc_001_config.h"
#include <stdio.h>
#include <math.h>
//#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_common.h"  //WriteF32ArrayToStrideMEM
#include "abc_vec.h"  //WriteF32ArrayToStrideMEM
#include "beastv2_io.h"
 

static void GetOutputOffsetStride_V2(int ROW, int COL, int whichDimIsTime, I64 idx, I64 Nvec, I64* pStride, I64* pOffset)
{
	//index should be 1-based.
	I64 stride, offset;
	if (ROW*COL==1)  
	// A 1d vcetor: io->ndim == 1
		stride = 1L, offset = (idx - 1) * Nvec;
	else if (ROW==1 || COL==1) {
	// A 2 D mat input:io->ndim == 2L
	//TODO: buggy: a edge case may be  a 3D input of 1x1 in dimension

		I32 nPixel = ROW * COL;
		if(whichDimIsTime ==1)
			stride = 1L, offset = (idx - 1) * Nvec;
		else {
			stride = nPixel;
			offset = (idx - 1);
		}
	} else  { 
	// A 3D stack input: io->ndim == 3L

		switch (whichDimIsTime) {
		case 1:
			stride = 1L;
			offset = (idx - 1) * Nvec;
			break;
		case 2: {
			int c = (idx - 1) / ROW;
			int r = idx - c * ROW;
			c = c + 1;
			stride = ROW;
			offset = (c - 1) * (Nvec * ROW) + r - 1;
			break;	}
		case 3: {
			int c = (idx - 1) / ROW;
			int r = idx - c * ROW;
			c = c + 1;
			stride = ROW * COL;
			offset = (c - 1) * ROW + r - 1;
			break;	}
		} // switch ( io->whichDimIsTime)	 

	} // (io->ndim == 3L)

	*pStride = stride;
	*pOffset = offset;
}

void* BEAST2_TsExtract(void *o, void * pindex ) {
 

	VOID_PTR    tmp;	
	
	// Get basic dimension parameters from 0
	DATA_TYPE   dtype;
	int         N;
	tmp   = GetField(o, "time");
	dtype = GetDataType(tmp);
	N     = GetNumberOfElements(tmp);
	
	int        whichDimTim;
	tmp         = GetField(o, "whichOutDimIsTime");
	whichDimTim = (tmp==NULL)?1: GetNumericElement(tmp, 0);
	
	int        ROW, COL;
	tmp = GetField(o, "marg_lik");
	if (GetNumberOfElements(tmp) == 1) {
		ROW = COL = 1;
		r_printf("BEAST2_TsExtract: this branch should never run!\n");
		return NULL;
	} else {		
		// if there are more than one TS, ncrows and ncols must be present.
        ROW = GetNumericElement(GetField(o, "nrows"), 0); 
		COL = GetNumericElement(GetField(o, "ncols"), 0);	
	}

	int index = 1L;
	if (IsNumeric(pindex)) {
		int numel = GetNumberOfElements(pindex);
		if       (numel == 1) {

			index = GetScalar(pindex);
			if (ROW >1 || COL> 1) {
				//r_printf("The input x contains a 3D array of time series and the index should "
				//	      "be a vector of two integers, specifying the row and column of the "
				//	      "selected pixel.\n");
			}

		} else if (numel >= 2) {

			int row = GetNumericElement(pindex, 0);
			int col = GetNumericElement(pindex, 1);
			index = (col - 1) * ROW + row;

			if (ROW == 1 || COL == 1) {
				//r_printf("The input x contains a 2D array of time series and the index should "
				//	"be a single integer, specifying the subscript of the selected time series\n.");
			}

		}	
		
		index = min(index, ROW * COL);
		index = max(index, 1L);

	}	else {
		r_printf("The index should be one or two integers (for 3D array of time series).\n");
		index = 1;
	}


	int maxKnotNumT   =0, maxKnotNumS   = 0,  maxKnotNumO     = 0;
	I08 hasSeasonCmpnt=0, hasTrendCmpnt = 0,  hasOutlierCmpnt = 0;
	tmp = GetField(o, "trend");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumT = size / ROW / COL -1L;
		hasTrendCmpnt = 1;
	}
	tmp = GetField(o, "season");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumS = size / ROW / COL - 1L;
		hasSeasonCmpnt = 1;
	}
	tmp = GetField(o, "outlier");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumO     = size / ROW / COL - 1L;
		hasOutlierCmpnt = 1;
	}

	int         mxKnotNum;
	F32PTR      TmpPointer;
	//https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                 (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                   (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&TmpPointer }
	#define _1(name, ...)                _(name, __VA_ARGS__)  
	#define _2(name1,name2, ...)         _(name1, __VA_ARGS__),  _(name2, __VA_ARGS__)   
	#define _3(n1,n2,n3, ...)            _2(n1, n2, __VA_ARGS__),  _(n3, __VA_ARGS__)  
	#define _4(n1,n2,n3,n4, ...)         _3(n1, n2,n3, __VA_ARGS__),  _(n4, __VA_ARGS__)  
	#define _5(n1,n2,n3,n4,n5, ...)      _4(n1, n2,n3,n4, __VA_ARGS__),  _(n5, __VA_ARGS__)  
	#define _6(n1,n2,n3,n4,n5,n6, ...)   _5(n1, n2,n3,n4,n5, __VA_ARGS__),  _(n6, __VA_ARGS__)  
	//{ "time",    dtype, 2, { N,               1, 0L, 0L }, &mat->time },

	/*********************************************************/
	// For univariate time series, _q() is the same as _()
	/*********************************************************/
	int isMultiVariate = 0L;
	#define _q(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&TmpPointer, isMultiVariate }
	#define _q1(name, ...)                _q(name, __VA_ARGS__)  
	#define _q2(name1,name2, ...)         _q(name1, __VA_ARGS__),  _q(name2, __VA_ARGS__)   
	#define _q3(n1,n2,n3, ...)            _q2(n1, n2, __VA_ARGS__),  _q(n3, __VA_ARGS__)  
	#define _q4(n1,n2,n3,n4, ...)         _q3(n1, n2,n3, __VA_ARGS__),  _q(n4, __VA_ARGS__)  
	#define _q5(n1,n2,n3,n4,n5, ...)      _q4(n1, n2,n3,n4, __VA_ARGS__),  _q(n5, __VA_ARGS__)  
	#define _q6(n1,n2,n3,n4,n5,n6, ...)   _q5(n1, n2,n3,n4,n5, __VA_ARGS__),  _q(n6, __VA_ARGS__) 

 
	mxKnotNum = maxKnotNumT;
	FIELD_ITEM  fieldListT[ ]= {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10,  1),
			_(ncpPr,                     mxKnotNum + 1),
			_2(cpOccPr, order,   N),			
			_2(cp,  cpPr,        mxKnotNum),
            _q(cpAbruptChange,   mxKnotNum),
			_(cpCI,				 mxKnotNum, 2),

			_q2(Y, SD,   N),
			_q(CI,       N,2),

		     _q4(slp, slpSD, slpSgnPosPr,  slpSgnZeroPr, N),

			// tpos_ncp, tneg_ncp, tinc_ncp, tdec_ncp,
			_q2(pos_ncp,    neg_ncp,      1),
			_q2(pos_ncpPr,  neg_ncpPr,    mxKnotNum + 1),
			_q2(pos_cpOccPr,neg_cpOccPr,  N),
			_q6(pos_cp,     neg_cp,   pos_cpPr, neg_cpPr, pos_cpAbruptChange, neg_cpAbruptChange, mxKnotNum),
			_q2(pos_cpCI,   neg_cpCI,     mxKnotNum,2),

			_q2(inc_ncp,      dec_ncp,      1),
			_q2(inc_ncpPr,    dec_ncpPr,    mxKnotNum + 1),
			_q2(inc_cpOccPr,  dec_cpOccPr,  N),
			_q6(inc_cp,    dec_cp, inc_cpPr, dec_cpPr, inc_cpAbruptChange, dec_cpAbruptChange, mxKnotNum),
			_q2(inc_cpCI,  dec_cpCI,        mxKnotNum,2),
	};

	mxKnotNum = maxKnotNumS;
	FIELD_ITEM  fieldListS[] = {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10,  1),
			_(ncpPr,    mxKnotNum + 1),				
			_2(cpOccPr, order, N),													
			_2(cp,  cpPr,  mxKnotNum),
           	_q(cpAbruptChange, mxKnotNum),
			_(cpCI,			   mxKnotNum, 2 ),

			_q2(Y ,SD,   N   ),
			_q(CI,       N,2),
			_q2(amp,  ampSD, N),

			// tpos_ncp, tneg_ncp, tinc_ncp, tdec_ncp,
			_q2(pos_ncp,    neg_ncp,      1 ),
			_q2(pos_ncpPr,  neg_ncpPr,    mxKnotNum+ 1),
			_q2(pos_cpOccPr,neg_cpOccPr,  N),
			_q6(pos_cp,     neg_cp,   pos_cpPr, neg_cpPr, pos_cpAbruptChange, neg_cpAbruptChange, mxKnotNum),
			_q2(pos_cpCI,   neg_cpCI, mxKnotNum,2)				
	};

	mxKnotNum = maxKnotNumO;
	FIELD_ITEM  fieldListO[ ]= {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10,  1),
			_(ncpPr,   mxKnotNum + 1),
			_(cpOccPr, N),

			_2(cp, cpPr,  mxKnotNum),
			_(cpCI,       mxKnotNum, 2),

			_q2(Y, SD, N),
			_q(CI,     N,2),

			_q2(pos_ncp,     neg_ncp,      1),
			_q2(pos_ncpPr,   neg_ncpPr,    mxKnotNum + 1),
			_q2(pos_cpOccPr, neg_cpOccPr,  N),
			_q4(pos_cp,      neg_cp, pos_cpPr, neg_cpPr, mxKnotNum),
			_q2(pos_cpCI,    neg_cpCI,     mxKnotNum,2),
	};
	int nfieldsT = sizeof(fieldListT) / sizeof(FIELD_ITEM);
	int nfieldsS = sizeof(fieldListS) / sizeof(FIELD_ITEM);
	int nfieldsO = sizeof(fieldListO) / sizeof(FIELD_ITEM);

	I32       nprt  = 0;
	VOID_PTR  trend = NULL, season=NULL, outlier=NULL;
	FIELD_ITEM* fieldList;
	int         nfields;

	/***************************************/
	// Create and allocate mem for the three components
	/***************************************/
	fieldList = fieldListT;
	nfields   = nfieldsT;
	if (hasTrendCmpnt)   { 
		VOID_PTR cmpnt = GetField(o, "trend");
		for (int i = 0; i < nfields; i++) {			
			tmp = GetField(cmpnt, fieldList[i].name);
			if (tmp == NULL || IsEmpty(tmp)) { fieldList[i].ptr = NULL; }	// Remove the field if it doesn't exist in o		
		}
		trend   = PROTECT(CreateStructVar(fieldList, nfields));;    nprt++;
	}
	fieldList = fieldListS;
	nfields   = nfieldsS;
	if (hasSeasonCmpnt)  { 
		VOID_PTR cmpnt = GetField(o, "season");
		for (int i = 0; i < nfields; i++) {			
			tmp = GetField(cmpnt, fieldList[i].name);
			if (tmp == NULL || IsEmpty(tmp)) { fieldList[i].ptr = NULL; }	// Remove the field if it doesn't exist in o		
		}
		season = PROTECT(CreateStructVar(fieldList, nfields));;    nprt++;
	}
	fieldList = fieldListO;
	nfields   = nfieldsO;
	if (hasOutlierCmpnt) {
		VOID_PTR cmpnt = GetField(o, "outlier");
		for (int i = 0; i < nfields; i++) {
			tmp = GetField(cmpnt, fieldList[i].name);
			if (tmp == NULL || IsEmpty(tmp)) { fieldList[i].ptr = NULL; }	// Remove the field if it doesn't exist in o		
		}
		outlier = PROTECT(CreateStructVar(fieldList, nfields));;    nprt++;
	}

   // TmpPointer is just a placeholder; the filled value is not actually used 
   // THe data address will be explicilty obtained by qurying with the name
	FIELD_ITEM  fieldListBEAST[ ] ={
		{"time",      dtype,		1,   {N,},        &TmpPointer},		
		{"data",      dtype,		1,   {N, },       &TmpPointer, .extra= isMultiVariate}, //Needed to be changed to reflect the output dim
		{"marg_lik",  dtype,		1,   {1,},        &TmpPointer},
		{"R2",        dtype,		1,   {1,},        &TmpPointer,  .extra = isMultiVariate},
		{"RMSE",      dtype,		1,   {1,},        &TmpPointer,.extra = isMultiVariate},		
		{"sig2",      dtype,		1,   {1,},        &TmpPointer},
		{"trend",     DATA_STRUCT,	0,	 {0,},       (void**)trend},
		{"season",    DATA_STRUCT,	0,	 {0,},       (void**)season},
		{"outlier",   DATA_STRUCT,	0,	 {0,},       (void**)outlier}, 
	};
	I32    nfieldsBEAST = sizeof(fieldListBEAST) / sizeof(FIELD_ITEM);

	// Check if data exists
	tmp = GetField(o, "data");
	if (tmp == NULL || IsEmpty(tmp)) {
		fieldListBEAST[1].ptr = NULL;
	}

	VOID_PTR  out;
	out = PROTECT(CreateStructVar(fieldListBEAST, nfieldsBEAST));          nprt++;
	AddStringAttribute(out,  "class",          "beast");	
	

	/************************************************/
	// Copy elements to the newly created array
	/************************************************/
	I64 stride, offset;

	if (hasTrendCmpnt) {
		fieldList = fieldListT;
		nfields   = nfieldsT;
		VOID_PTR cmpnt = GetField(o, "trend");
		for (int i = 0; i < nfields; i++) {

			if (fieldList[i].ptr == NULL) {
				continue;
			}			
			int Nvec = fieldList[i].ndim == 1 ? fieldList[i].dims[0] : fieldList[i].dims[0] * fieldList[i].dims[1];			
			GetOutputOffsetStride_V2(ROW, COL, whichDimTim,  index, Nvec, &stride, &offset);

			VOID_PTR newData = GetData(  GetField(trend, fieldList[i].name) );
			VOID_PTR oldData = GetData(  GetField(cmpnt, fieldList[i].name));
			WriteStrideMEMToArrMEM(newData, oldData, Nvec, stride, offset, dtype);	
		}
	
	}
	if (hasSeasonCmpnt) {
		fieldList = fieldListS;
		nfields   = nfieldsS;
		VOID_PTR cmpnt = GetField(o, "season");
		for (int i = 0; i < nfields; i++) {

			if (fieldList[i].ptr == NULL) {
				continue;
			}			
			int Nvec = fieldList[i].ndim == 1 ? fieldList[i].dims[0] : fieldList[i].dims[0] * fieldList[i].dims[1];			
			GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);

			VOID_PTR newData = GetData(  GetField(season, fieldList[i].name) );
			VOID_PTR oldData = GetData(  GetField(cmpnt, fieldList[i].name));
			WriteStrideMEMToArrMEM(newData, oldData, Nvec, stride, offset, dtype);
		}
	}
	if (hasOutlierCmpnt) {
		fieldList = fieldListO;
		nfields   = nfieldsO;
		VOID_PTR cmpnt = GetField(o, "outlier");
		for (int i = 0; i < nfields; i++) {

			if (fieldList[i].ptr == NULL) {
				continue;
			}			
			int Nvec = fieldList[i].ndim == 1 ? fieldList[i].dims[0] : fieldList[i].dims[0] * fieldList[i].dims[1];			
			GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);

			VOID_PTR newData = GetData(  GetField(outlier, fieldList[i].name) );
			VOID_PTR oldData = GetData(  GetField(cmpnt, fieldList[i].name));
			WriteStrideMEMToArrMEM(newData, oldData, Nvec, stride, offset, dtype);
		}
	}

	// Ti,e
	{
		VOID_PTR newData = GetData(GetField(out, fieldListBEAST[0].name));
		VOID_PTR oldData = GetData(GetField(o, fieldListBEAST[0].name));
		WriteStrideMEMToArrMEM(newData, oldData, N, 1, 0, dtype);

		for (int i = 1; i < 6; i++) {

			if (fieldListBEAST[i].ptr == NULL) {
				continue;
			}

			int Nvec = fieldListBEAST[i].ndim == 1 ? fieldListBEAST[i].dims[0] : fieldListBEAST[i].dims[0] * fieldListBEAST[i].dims[1];
			GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);

			VOID_PTR newData = GetData(GetField(out, fieldListBEAST[i].name));
			VOID_PTR oldData = GetData(GetField(o, fieldListBEAST[i].name));
			WriteStrideMEMToArrMEM(newData, oldData, Nvec, stride, offset, dtype);
		}

	}



	UNPROTECT(nprt);
	return out;

	#undef NUMARGS
    #undef NARGS
    #undef _
	#undef _2
	#undef _3
	#undef _4
    #undef _5
	#undef _6
	#undef _7
	#undef _q
	#undef _q2
	#undef _q3
	#undef _q4
    #undef _q5
	#undef _q6
	#undef _q7
}
 
 
 void* BEAST2_PrintResult(void *o, void * pindex ) {

	VOID_PTR    tmp;	
	
	// Get basic dimension parameters from 0
	DATA_TYPE   dtype;
	int         N;
	tmp   = GetField(o, "time");
	dtype = GetDataType(tmp);
	N     = GetNumberOfElements(tmp);
	
	int        whichDimTim;
	tmp         = GetField(o, "whichOutDimIsTime");
	whichDimTim = (tmp==NULL)?1: GetNumericElement(tmp, 0);
	
	int        ROW, COL;
	tmp = GetField(o, "marg_lik");
	if (GetNumberOfElements(tmp) == 1) {
		ROW = COL = 1;
	} else {		
		// if there are more than one TS, ncrows and ncols must be present.
        ROW = GetNumericElement(GetField(o, "nrows"), 0); 
		COL = GetNumericElement(GetField(o, "ncols"), 0);	
	}

	int maxKnotNumT   =0, maxKnotNumS   = 0,  maxKnotNumO     = 0;
	I08 hasSeasonCmpnt=0, hasTrendCmpnt = 0,  hasOutlierCmpnt = 0;
	tmp = GetField(o, "trend");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumT = size / ROW / COL -1L;
		hasTrendCmpnt = 1;
	}
	tmp = GetField(o, "season");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumS = size / ROW / COL - 1L;
		hasSeasonCmpnt = 1;
	}
	tmp = GetField(o, "outlier");
	if (tmp != NULL && !IsEmpty(tmp)) {
		int size = GetNumberOfElements(GetField(tmp, "ncpPr"));
		maxKnotNumO     = size / ROW / COL - 1L;
		hasOutlierCmpnt = 1;
	}


	int index = 1L;
	if (IsNumeric(pindex)) {
		int numel = GetNumberOfElements(pindex);
		if       (numel == 1) {

			index = GetScalar(pindex);
			if (ROW >1 || COL> 1) {
				r_printf("The input x contains a 3D array of time series and the index should "
					      "be a vector of two integers, specifying the row and column of the "
					      "selected pixel.\n");
			}

		} else if (numel >= 2) {

			int row = GetNumericElement(pindex, 0);
			int col = GetNumericElement(pindex, 1);
			index = (col - 1) * ROW + row;

			if (ROW == 1 || COL == 1) {
				r_printf("The input x contains a 2D array of time series and the index should "
					"be a single integer, specifying the subscript of the selected time series\n.");
			}

		}	
		
		index = min(index, ROW * COL);
		index = max(index, 1L);

	}	
	else {
		r_printf("The index should be one or two integers (for 3D array of time series).\n");
		index = 1;
	}




	F32    tmpbuf[100];
	int    isallocated = 0;
	F32PTR newData = tmpbuf;
	int    maxLen = max(maxKnotNumT, max(maxKnotNumS, maxKnotNumO)) + 1;
	if (maxLen > 100) {
		newData     = malloc(sizeof(F32) * maxLen);
		isallocated = 1;
	}

	I64  stride, offset;
	char s1[] = "                                                ";
	int  nChar = strlen(s1);

	#define cat r_printf

#if R_INTERFACE==1
	#define boldStart  "\033[1;31m" 
    #define boldEnd    "\033[0m" 
#else
	#define boldStart  "<strong>" 
    #define boldEnd   "</strong>" 
#endif


	cat(boldStart);
	cat("#####################################################################\n");
	cat("#                      Seasonal  Changepoints                       #\n");
	cat("#####################################################################\n");
	cat(boldEnd);

	if (!hasSeasonCmpnt) {
		cat(" No seasonal/periodic component present (i.e., season='none')\n");
	}
		
	if (hasSeasonCmpnt) {

		int      maxKnotNum = maxKnotNumS;
		VOID_PTR cmpnt      = GetField(o, "season");

	
		int      Nvec       = maxKnotNum +1L;
		GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		VOID_PTR oldData  = GetData(GetField(cmpnt, "ncpPr"));
		CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 maxPr;
		I32 maxIx = f32_maxidx(newData, Nvec, &maxPr);
		Nvec = min(Nvec, 99);
		cat(".-------------------------------------------------------------------.\n");
		cat("| Ascii plot of probability distribution for number of chgpts (ncp) |\n");
		cat(".-------------------------------------------------------------------.\n");
		for ( int i=0; i<Nvec; i++) {
			int slen;
			slen = ceil(newData[i] / maxPr * (nChar - 1));
			slen = max(1, slen);
			memset(s1, ' ', nChar);
			memset(s1, '*', slen);
		    cat("|Pr(ncp = %-2d)=%.3f|%s|\n", i, newData[i], s1);
		}

		Nvec    = 1;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "ncp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_median")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_median = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct90")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct90 = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct10")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct10 = newData[0];

		cat(".-------------------------------------------------------------------.\n"); 
		cat("|    Summary for number of Seasonal ChangePoints (scp)              |\n");
		cat(".-------------------------------------------------------------------.\n");
		cat("|ncp_max    = %-4d | MaxSeasonKnotNum: A parameter you set          |\n",     maxKnotNum);
		cat("|ncp_mode   = %-4d | Pr(ncp=%2d)=%3.2f: There is a %3.1f%% probability  |\n", maxIx, min(maxIx, 99), maxPr, maxPr * 100);
		cat("|                  | that the seasonal component has %2d chgnpt(s).  |\n",    min(maxIx, 99));
		cat("|ncp_mean   = %-4.2f | Sum{ncp*Pr(ncp)} for ncp = 0,...,%-4d          |\n",   ncp, maxKnotNum);
		cat("|ncp_pct10  = %-4.2f | 10%% percentile for number of changepoints      |\n", ncp_pct10);
		cat("|ncp_median = %-4.2f | 50%% percentile: Median number of changepoints  |\n", ncp_median);
		cat("|ncp_pct90  = %-4.2f | 90%% percentile for number of changepoints      |\n", ncp_pct90);
		cat(".-------------------------------------------------------------------.\n");
		cat("| List of probable seasonal changepoints ranked by probability of   |\n");
		cat("| occurrence: Please combine the ncp reported above to determine    |\n");
		cat("| which changepoints below are  practically meaningful              |\n");
		cat("'-------------------------------------------------------------------'\n");
		cat("|scp#              |time (cp)                  |prob(cpPr)          |\n");
		cat("|------------------|---------------------------|--------------------|\n");

		F32 cp[200];
		F32 cpPr[200];

		Nvec    = maxKnotNum;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "cp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);		
		memcpy(cp, newData, sizeof(F32)* min(200, Nvec));

		oldData = GetData(GetField(cmpnt, "cpPr"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		memcpy(cpPr, newData, sizeof(F32) * min(200, Nvec));

		Nvec = min(200L, Nvec);

		int ncp_all=0;
		for (int i = 0; i < Nvec; i++) {
			ncp_all = ncp_all + (cp[i] == cp[i]);
		}
		for (int i = 0; i < ncp_all; i++) {
			cat("|%-18d|%-27.6f|%-20.5f|\n", i + 1, cp[i], cpPr[i]);
		}	 
		cat(".-------------------------------------------------------------------.\n\n");
			
	}
	cat("\n\n");
 
 
	cat(boldStart);
	cat("#####################################################################\n");
	cat("#                      Trend  Changepoints                          #\n");
	cat("#####################################################################\n");
	cat(boldEnd);

	if (hasTrendCmpnt) {
		int      maxKnotNum = maxKnotNumT;
		VOID_PTR cmpnt      = GetField(o, "trend");

		int      Nvec       = maxKnotNum +1L;
		GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		VOID_PTR oldData  = GetData(GetField(cmpnt, "ncpPr"));
		CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 maxPr;
		I32 maxIx = f32_maxidx(newData, Nvec, &maxPr);
		Nvec = min(Nvec, 99);
		cat(".-------------------------------------------------------------------.\n");
		cat("| Ascii plot of probability distribution for number of chgpts (ncp) |\n");
		cat(".-------------------------------------------------------------------.\n");
		for ( int i=0; i<Nvec; i++) {
			int slen;
			slen = ceil(newData[i] / maxPr * (F32)(nChar - 1));
			slen = max(1, slen);
			memset(s1, ' ', nChar);
			memset(s1, '*', slen);
		    cat("|Pr(ncp = %-2d)=%.3f|%s|\n", i, newData[i], s1);
		}

		Nvec    = 1;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "ncp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_median")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_median = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct90")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct90 = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct10")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct10 = newData[0];

		cat(".-------------------------------------------------------------------.\n");
		cat("|    Summary for number of Trend ChangePoints (tcp)                 |\n");
		cat(".-------------------------------------------------------------------.\n");
		cat("|ncp_max    = %-4d | MaxTrendKnotNum: A parameter you set           |\n", maxKnotNum);
		cat("|ncp_mode   = %-4d | Pr(ncp=%2d)=%3.2f: There is a %3.1f%% probability  |\n", maxIx, min(maxIx, 99), maxPr, maxPr * 100);
		cat("|                  | that the trend component has %2d changepoint(s).|\n",  min(maxIx, 99));
		cat("|ncp_mean   = %-4.2f | Sum{ncp*Pr(ncp)} for ncp = 0,...,%-4d          |\n", ncp, maxKnotNum);
		cat("|ncp_pct10  = %-4.2f | 10%% percentile for number of changepoints      |\n", ncp_pct10);
		cat("|ncp_median = %-4.2f | 50%% percentile: Median number of changepoints  |\n", ncp_median);
		cat("|ncp_pct90  = %-4.2f | 90%% percentile for number of changepoints      |\n", ncp_pct90);
		cat(".-------------------------------------------------------------------.\n");
		cat("| List of probable trend changepoints ranked by probability of      |\n");
		cat("| occurrence: Please combine the ncp reported above to determine    |\n");
		cat("| which changepoints below are  practically meaningful              |\n");
		cat("'-------------------------------------------------------------------'\n");
		cat("|tcp#              |time (cp)                  |prob(cpPr)          |\n");
		cat("|------------------|---------------------------|--------------------|\n");

		F32 cp[200];
		F32 cpPr[200];

		Nvec    = maxKnotNum;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "cp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);		
		memcpy(cp, newData, sizeof(F32)* min(200, Nvec));

		oldData = GetData(GetField(cmpnt, "cpPr"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		memcpy(cpPr, newData, sizeof(F32) * min(200, Nvec));

		Nvec = min(200L, Nvec);

		int ncp_all=0;
		for (int i = 0; i < Nvec; i++) {
			ncp_all = ncp_all + (cp[i] == cp[i]);
		}
		for (int i = 0; i < ncp_all; i++) {
			cat("|%-18d|%-27.6f|%-20.5f|\n", i + 1, cp[i], cpPr[i]);
		}	 
		cat(".-------------------------------------------------------------------.\n\n");
			
	}
	cat("\n\n");
 
	if (hasOutlierCmpnt) {
		cat(boldStart);
		cat("#####################################################################\n");
		cat("#                      Outlier  Changepoints                        #\n");
		cat("#####################################################################\n");
		cat(boldEnd);
		int      maxKnotNum = maxKnotNumO;
		VOID_PTR cmpnt      = GetField(o, "outlier");


		int      Nvec       = maxKnotNum +1L;
		GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		VOID_PTR oldData = GetData(GetField(cmpnt, "ncpPr"));
		CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 maxPr;
		I32 maxIx = f32_maxidx(newData, Nvec, &maxPr);
		Nvec = min(Nvec, 99);
		cat(".-------------------------------------------------------------------.\n");
		cat("| Ascii plot of probability distribution for number of chgpts (ncp) |\n");
		cat(".-------------------------------------------------------------------.\n");
		for ( int i=0; i<Nvec; i++) {
			int slen;
			slen = ceil(newData[i] / maxPr * (nChar - 1));
			slen = max(1, slen);
			memset(s1, ' ', nChar);
			memset(s1, '*', slen);
		    cat("|Pr(ncp = %-2d)=%.3f|%s|\n", i, newData[i], s1);
		}

		Nvec    = 1;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "ncp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_median")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_median = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct90")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct90 = newData[0];

		oldData = GetData(GetField(cmpnt, "ncp_pct10")); 	CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		F32 ncp_pct10 = newData[0];

		cat(".-------------------------------------------------------------------.\n");
		cat("|    Summary for number of Outlier ChangePoints (ocp)               |\n");
		cat(".-------------------------------------------------------------------.\n");
		cat("|ncp_max    = %-4d | MaxOutlierKnotNum: A parameter you set         |\n", maxKnotNum);
		cat("|ncp_mode   = %-4d | Pr(ncp=%2d)=%3.2f: There is a %3.1f%% probability  |\n", maxIx, min(maxIx, 99), maxPr, maxPr * 100);
		cat("|                  | that the outlier component has %2d chngpnt(s).  |\n", min(maxIx, 99));
		cat("|ncp_mean   = %-4.2f | Sum{ncp*Pr(ncp)} for ncp = 0,...,%-4d          |\n", ncp, maxKnotNum);
		cat("|ncp_pct10  = %-4.2f | 10%% percentile for number of changepoints      |\n", ncp_pct10);
		cat("|ncp_median = %-4.2f | 50%% percentile: Median number of changepoints  |\n", ncp_median);
		cat("|ncp_pct90  = %-4.2f | 90%% percentile for number of changepoints      |\n", ncp_pct90);
		cat(".-------------------------------------------------------------------.\n");
		cat("| List of probable outlier changepoints ranked by probability of    |\n");
		cat("| occurrence: Please combine the ncp reported above to determine    |\n");
		cat("| which changepoints below are  practically meaningful              |\n");
		cat("'-------------------------------------------------------------------'\n");
		cat("|ocp#              |time (cp)                  |prob(cpPr)          |\n");
		cat("|------------------|---------------------------|--------------------|\n");


		F32 cp[200];
		F32 cpPr[200];

		Nvec    = maxKnotNum;	GetOutputOffsetStride_V2(ROW, COL, whichDimTim, index, Nvec, &stride, &offset);
		oldData = GetData(GetField(cmpnt, "cp"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);		
		memcpy(cp, newData, sizeof(F32)* min(200, Nvec));

		oldData = GetData(GetField(cmpnt, "cpPr"));		  CopyStrideMEMToF32Arr(newData, oldData, Nvec, stride, offset, dtype);
		memcpy(cpPr, newData, sizeof(F32) * min(200, Nvec));

		Nvec = min(200L, Nvec);

		int ncp_all=0;
		for (int i = 0; i < Nvec; i++) {
			ncp_all = ncp_all + (cp[i] == cp[i]);
		}
		for (int i = 0; i < ncp_all; i++) {
			cat("|%-18d|%-27.6f|%-20.5f|\n", i + 1, cp[i], cpPr[i]);
		}	 
		cat(".-------------------------------------------------------------------.\n\n");

	}

	#if R_INTERFACE==1
	r_printf("NOTE: the beast output object 'o' is a LIST. Type 'str(o)' to see all \n"
		     "the elements in it. Or use 'plot(o)' or 'plot(o,interactive=TRUE)' to \n"
		     "plot the model output.\n");
	#endif
 

	if (isallocated) {
		free(newData);
	}

	return NULL;

	 
}
#include "abc_000_warning.h"