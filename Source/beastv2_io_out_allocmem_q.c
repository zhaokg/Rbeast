#include "abc_000_warning.h"
#include "abc_001_config.h"

#include <string.h>   // memset and memcpy

#include "abc_ide_util.h"
#include "abc_vec.h"
#include "abc_date.h"
#include "abc_mem.h"
#include "beastv2_io.h"


static void __RemoveFieldsGivenFlags_Trend(A(OPTIONS_PTR)  opt, FIELD_ITEM * fieldList, int nfields) {

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	const A(IO_PTR)      io   = &opt->io;
	const A(RESULT_PTR)  mat  = io->out.result;
	const A(EXTRA_PTR)   flag = &(opt->extra);

	if (!flag->computeCredible)		RemoveField(fieldList, nfields, "CI"),	  mat->tCI = NULL;	
	if (!flag->computeTrendOrder) 	RemoveField(fieldList, nfields, "order"), mat->torder = NULL;	
	if (!flag->computeTrendSlope)
		RemoveField(fieldList, nfields, "slp"),				mat->tslp		= NULL,
		RemoveField(fieldList, nfields, "slpSD"),			mat->tslpSD		= NULL,
		RemoveField(fieldList, nfields, "slpSgnPosPr"),		mat->tslpSgnPosPr	= NULL,
	    RemoveField(fieldList, nfields, "slpSgnZeroPr"),	mat->tslpSgnPosPr	= NULL;	
	if (!flag->computeTrendChngpt)
		RemoveField(fieldList, nfields, "cp"),				mat->tcp = NULL,
		RemoveField(fieldList, nfields, "cpPr"),			mat->tcpPr = NULL,
		RemoveField(fieldList, nfields, "cpCI"),			mat->tcpCI = NULL,
		RemoveField(fieldList, nfields, "cpAbruptChange"),  mat->tcpAbruptChange = NULL;
		

	#define  _(x)						RemoveField(fieldList, nfields, #x), mat->t##x = NULL
	#define _2(x,y)						_(x),    _(y)
	#define _3(x,y,z)					_2(x,y), _(z)
	#define _4(x,y,z,v)					_3(x,y,z), _(v)
	#define _6(x,y,z,v1,v2,v3)			_3(x,y,z),_3(v1,v2,v3)
	#define _7(x,y,z,v1,v2,v3,v4)		_6(x,y,z,v1,v2,v3),_(v4)

	if (!flag->tallyPosNegTrendJump) {
		_4(pos_ncp, neg_ncp, pos_ncpPr, neg_ncpPr),
		_2(pos_cpOccPr, neg_cpOccPr),
		_4(pos_cp, neg_cp, pos_cpPr, neg_cpPr),
		_2(pos_cpAbruptChange, neg_cpAbruptChange),
		_2(pos_cpCI, neg_cpCI);
	}
	if (!flag->tallyIncDecTrendJump) {
		_4(inc_ncp, dec_ncp, inc_ncpPr, dec_ncpPr),
		_2(inc_cpOccPr, dec_cpOccPr),
		_4(inc_cp,  dec_cp, inc_cpPr, dec_cpPr),
		_2(inc_cpAbruptChange, dec_cpAbruptChange),
		_2(inc_cpCI, dec_cpCI);
	}
 
	#undef _
	#undef _2
	#undef _3
	#undef _4
    #undef _5
	#undef _6
	#undef _7


}

static void __RemoveFieldsGivenFlags_Season(A(OPTIONS_PTR)  opt, FIELD_ITEM * fieldList, int nfields) {

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	const A(IO_PTR)      io  = &opt->io;
	const A(RESULT_PTR)  mat = io->out.result;
	const A(EXTRA_PTR)   flag = &(opt->extra);

	#define  _(x)						RemoveField(fieldList, nfields, #x), mat->s##x = NULL
	#define _2(x,y)						_(x),    _(y)
	#define _3(x,y,z)					_2(x,y), _(z)
	#define _4(x,y,z,v)					_3(x,y,z), _(v)
	#define _6(x,y,z,v1,v2,v3)			_3(x,y,z),_3(v1,v2,v3)
	#define _5(x,y,z,v1,v2)				_4(x,y,z,v1), _(v2)
	#define _7(x,y,z,v1,v2,v3,v4)		_6(x,y,z,v1,v2,v3),_(v4)

	if (!hasSeasonCmpnt) {
		_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10);
		_5(  ncpPr,   cpOccPr, Y,  SD, CI);
		_7(order,  amp, ampSD,   cp, cpPr, cpCI,  cpAbruptChange);
		_6(pos_ncp, neg_ncp, pos_ncpPr, neg_ncpPr, pos_cpOccPr, neg_cpOccPr);
		_4(pos_cp,  neg_cp,  pos_cpPr,  neg_cpPr  );
		_4(pos_cpAbruptChange, neg_cpAbruptChange, pos_cpCI, neg_cpCI);
		return;
	}

	if (!flag->computeCredible)		RemoveField(fieldList, nfields, "CI"),		mat->sCI = NULL;		
	if (!flag->computeSeasonOrder) 	RemoveField(fieldList, nfields, "order"),	mat->sorder = NULL;

	if (!flag->computeSeasonAmp)
		RemoveField(fieldList, nfields, "amp"),   mat->samp = NULL,
		RemoveField(fieldList, nfields, "ampSD"), mat->sampSD = NULL;

	if (!flag->computeSeasonChngpt)
		RemoveField(fieldList, nfields, "scp"),   mat->scp = NULL,
		RemoveField(fieldList, nfields, "scpPr"), mat->scpPr = NULL,
		RemoveField(fieldList, nfields, "scpCI"), mat->scpCI = NULL,
		RemoveField(fieldList, nfields, "scpAbruptChange"), mat->scpAbruptChange = NULL;
		

	if (!flag->tallyPosNegSeasonJump) {
		_4(pos_ncp, neg_ncp, pos_ncpPr, neg_ncpPr),
		_2(pos_cpOccPr, neg_cpOccPr),
		_4(pos_cp, neg_cp, pos_cpPr, neg_cpPr),
		_2(pos_cpAbruptChange, neg_cpAbruptChange),
		_2(pos_cpCI, neg_cpCI);
	}

    #undef _
	#undef _2
	#undef _3
	#undef _4
    #undef _5
	#undef _6
	#undef _7
}

static void __RemoveFieldsGivenFlags_Outlier(A(OPTIONS_PTR)  opt, FIELD_ITEM * fieldList, int nfields) {

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	const A(IO_PTR)      io  = &opt->io;
	const A(RESULT_PTR)  mat = io->out.result;
	const A(EXTRA_PTR)   flag = &(opt->extra);

	#define  _(x)						RemoveField(fieldList, nfields, #x), mat->o##x = NULL
	#define _2(x,y)						_(x),    _(y)
	#define _3(x,y,z)					_2(x,y), _(z)
	#define _4(x,y,z,v)					_3(x,y,z), _(v)
	#define _6(x,y,z,v1,v2,v3)			_3(x,y,z),_3(v1,v2,v3)
	#define _5(x,y,z,v1,v2)				_4(x,y,z,v1), _(v2)
	#define _7(x,y,z,v1,v2,v3,v4)		_6(x,y,z,v1,v2,v3),_(v4)

	if (!hasOutlierCmpnt) {				
		_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10);
		_5(ncpPr, cpOccPr, Y, SD, CI);
		_3(cp, cpPr, cpCI);
		_6(pos_ncp, neg_ncp, pos_ncpPr, neg_ncpPr, pos_cpOccPr, neg_cpOccPr);
		_4(pos_cp, neg_cp, pos_cpPr, neg_cpPr);
		_2( pos_cpCI, neg_cpCI);		
		return;
	}
 
	if (!flag->computeCredible)		RemoveField(fieldList, nfields, "CI"), mat->oCI = NULL;

	if (!flag->computeOutlierChngpt)
		RemoveField(fieldList, nfields, "cp"),		mat->ocp = NULL,
		RemoveField(fieldList, nfields, "cpPr"),	mat->ocpPr = NULL,
		RemoveField(fieldList, nfields, "cpCI"),	mat->ocpCI = NULL;	

	if (!flag->tallyPosNegOutliers) {
			_4(pos_ncp, neg_ncp, pos_ncpPr, neg_ncpPr),
			_2(pos_cpOccPr, neg_cpOccPr),
			_4(pos_cp, neg_cp, pos_cpPr, neg_cpPr),
			_2(pos_cpCI, neg_cpCI);
	}
	
    #undef _
	#undef _2
	#undef _3
	#undef _4
    #undef _5
	#undef _6
	#undef _7
}


static void   __AddSpatialDimension(int ROW, int COL, int whichOutDimIsTime, FIELD_ITEM * fieldList, int nfields) {
	  
  	if (whichOutDimIsTime == 1) {
		// Apend the nuMTS at the end
		for (int i = 0; i < nfields; i++) {
			FIELD_ITEM * fld  = fieldList + i;
			int          ndim = fld->ndim;
			fld->dims[ndim]     = ROW;
			fld->dims[ndim + 1] = COL;
			fld->ndim += 2;
		}
		return;
	}

	if (whichOutDimIsTime == 2) {
		// Insert ROW at the start and the COL at the end
		for (int i = 0; i < nfields; i++) {
			FIELD_ITEM * fld  = fieldList + i;
			int          ndim = fld->ndim;
			for (int i = ndim - 1; i >= 0; i--) { 
				fld->dims[i + 1] = fld->dims[i]; 
			}
			fld->dims[0]         = ROW;
			fld->dims[ndim + 1]  = COL;
			fld->ndim           += 2;
		}
		return;
	}

	if (whichOutDimIsTime == 3) {
		// Insert the numTS at the start
		for (int i = 0; i < nfields; i++) {
			FIELD_ITEM * fld  = fieldList + i;
			int          ndim = fld->ndim;
			for (int i = ndim - 1; i >= 0; i--) {
				fld->dims[i + 2] = fld->dims[i]; 
			}
			fld->dims[0] = ROW;
			fld->dims[1] = COL;
			fld->ndim    += 2;
		}
		return;
	}
	 
}


static  I32  __MR_ExtendFieldsToMultiVaraiteTS(FIELD_ITEM *flist, I32 N, I32 q) {

	if (q == 1) 
		return 0;

	//char* nms[] = { "Y1","Y2" ,"Y3", };

	I32 nptr       = 0;
	I32 nptr_dummy = 0;
	for (int i = 0; i < N; i++) {

		if (flist[i].extra == 0 || flist[i].ptr == NULL)
			continue;

		nptr_dummy++;
		FIELD_ITEM qList[100] = { { /*fielditem*/{/*char*/0,},}, };
		for (int j = 0; j < q; j++) {
			snprintf(qList[j].name, 63, "Y%d", j+1); 	//strcpy(qList[j].name, nms[j]);
			qList[j].extra = 0;
			qList[j].type  = flist[i].type;
			qList[j].ndim  = flist[i].ndim;
			memcpy(&qList[j].dims, &flist[i].dims, sizeof(I32) * 5);
	
			if (flist[i].ptr != NULL)
				qList[j].ptr = (char*)(flist[i].ptr) + sizeof(BEAST2_RESULT) * j;
			else 
				// flist[i].ptr==NULL has been skipped in the above, so the else branch will never run
				qList[j].ptr = NULL;
		}
		nptr_dummy--;
		VOID_PTR  out=PROTECT(CreateStructVar(qList, q)); nptr++;
		
		//TODO: buggy: out is not procted after this point, and it may be GCed.

		flist[i].extra = 0;
		flist[i].ndim  = 0;
		flist[i].type  = DATA_STRUCT;
		flist[i].ptr   = out;
	}
 
 
	 UNPROTECT(nptr_dummy); //nptr_dummy=0: added here to trick R's ptr checker
	 return nptr;
}

 
static void* __BEAST2_Output_AllocMEM_Trend(A(OPTIONS_PTR)  opt) {

	const BEAST2_IO_PTR      io = &opt->io;
	const BEAST2_RESULT_PTR  mat = io->out.result;
	DATA_TYPE   dtype     = io->out.dtype; // DATA_FLOAT or DATA_DOUBLE                            
	const int   N         = io->N;
	const int   mxKnotNum = opt->prior.trendMaxKnotNum;
	//https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                 (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                   (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->t##name }
	#define _1(name, ...)                _(name, __VA_ARGS__)  
	#define _2(name1,name2, ...)         _(name1, __VA_ARGS__),  _(name2, __VA_ARGS__)   
	#define _3(n1,n2,n3, ...)            _2(n1, n2, __VA_ARGS__),  _(n3, __VA_ARGS__)  
	#define _4(n1,n2,n3,n4, ...)         _3(n1, n2,n3, __VA_ARGS__),  _(n4, __VA_ARGS__)  
	#define _5(n1,n2,n3,n4,n5, ...)      _4(n1, n2,n3,n4, __VA_ARGS__),  _(n5, __VA_ARGS__)  
	#define _6(n1,n2,n3,n4,n5,n6, ...)   _5(n1, n2,n3,n4,n5, __VA_ARGS__),  _(n6, __VA_ARGS__)  
	//{ "time",    dtype, 2, { N,               1, 0L, 0L }, &mat->time },
	//{ "scp",	   dtype, 2, { mxKnotNumSeason, M, },        & mat->scp },

	int isMultiVariate = (io->q == 1) ? 0L : 1L;
	#define _q(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->t##name, isMultiVariate }
	#define _q1(name, ...)                _q(name, __VA_ARGS__)  
	#define _q2(name1,name2, ...)         _q(name1, __VA_ARGS__),  _q(name2, __VA_ARGS__)   
	#define _q3(n1,n2,n3, ...)            _q2(n1, n2, __VA_ARGS__),  _q(n3, __VA_ARGS__)  
	#define _q4(n1,n2,n3,n4, ...)         _q3(n1, n2,n3, __VA_ARGS__),  _q(n4, __VA_ARGS__)  
	#define _q5(n1,n2,n3,n4,n5, ...)      _q4(n1, n2,n3,n4, __VA_ARGS__),  _q(n5, __VA_ARGS__)  
	#define _q6(n1,n2,n3,n4,n5,n6, ...)   _q5(n1, n2,n3,n4,n5, __VA_ARGS__),  _q(n6, __VA_ARGS__) 

	// For univariate time series, _q() is the same as _()
	FIELD_ITEM  fieldList[] = {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10,    1),
			_1(ncpPr,            mxKnotNum + 1),
			_2(cpOccPr, order,   N),
			_2(cp,  cpPr,        mxKnotNum),
			_q(cpAbruptChange,   mxKnotNum),
			_1(cpCI,			 mxKnotNum, 2),

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

	I32 nfields = sizeof(fieldList) / sizeof(FIELD_ITEM);
	__AddSpatialDimension(io->dims[io->rowdim-1], io->dims[io->coldim - 1], io->out.whichDimIsTime, fieldList, nfields);
	__RemoveFieldsGivenFlags_Trend(opt, fieldList, nfields);
	if (opt->extra.removeSingletonDims) {
		RemoveSingltonDims(fieldList, nfields);
	}

	// Do nothing if io->q is 1L;  within the function, newly created lists are protected
	I32       nptr = __MR_ExtendFieldsToMultiVaraiteTS(fieldList, nfields, io->q); 
	VOID_PTR  out  = PROTECT(CreateStructVar(fieldList, nfields));
					 UNPROTECT(1L);
	UNPROTECT(nptr); 
	return out;

	#undef NUMARGS
    #undef NARGS
    #undef _
	#undef _1
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

static void* __BEAST2_Output_AllocMEM_Season(A(OPTIONS_PTR)  opt)
{
	const BEAST2_IO_PTR      io     = &opt->io;
	const BEAST2_RESULT_PTR  mat    = io->out.result;
	DATA_TYPE   dtype        = io->out.dtype; // DATA_FLOAT or DATA_DOUBLE                            
	const int   N	         = io->N;
	const int   mxKnotNum	 = opt->prior.seasonMaxKnotNum;	 
	//https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                  (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->s##name }
	#define _1(name, ...)               _(name, __VA_ARGS__)  
	#define _2(name1,name2, ...)         _(name1, __VA_ARGS__),  _(name2, __VA_ARGS__)   
	#define _3(n1,n2,n3, ...)            _2(n1, n2, __VA_ARGS__),  _(n3, __VA_ARGS__)  
	#define _4(n1,n2,n3,n4, ...)         _3(n1, n2,n3, __VA_ARGS__),  _(n4, __VA_ARGS__)  
	#define _5(n1,n2,n3,n4,n5, ...)      _4(n1, n2,n3,n4, __VA_ARGS__),  _(n5, __VA_ARGS__)  
	#define _6(n1,n2,n3,n4,n5,n6, ...)   _5(n1, n2,n3,n4,n5, __VA_ARGS__),  _(n6, __VA_ARGS__)  
	//{ "time",    dtype, 2, { N,               1, 0L, 0L }, &mat->time },
	//{ "scp",	   dtype, 2, { mxKnotNumSeason, M, },        & mat->scp },

	int isMultiVariate = (io->q == 1) ? 0L : 1L;
	#define _q(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->s##name, isMultiVariate }
	#define _q1(name, ...)                _q(name, __VA_ARGS__)  
	#define _q2(name1,name2, ...)         _q(name1, __VA_ARGS__),  _q(name2, __VA_ARGS__)   
	#define _q3(n1,n2,n3, ...)            _q2(n1, n2, __VA_ARGS__),  _q(n3, __VA_ARGS__)  
	#define _q4(n1,n2,n3,n4, ...)         _q3(n1, n2,n3, __VA_ARGS__),  _q(n4, __VA_ARGS__)  
	#define _q5(n1,n2,n3,n4,n5, ...)      _q4(n1, n2,n3,n4, __VA_ARGS__),  _q(n5, __VA_ARGS__)  
	#define _q6(n1,n2,n3,n4,n5,n6, ...)   _q5(n1, n2,n3,n4,n5, __VA_ARGS__),  _q(n6, __VA_ARGS__) 

	FIELD_ITEM  fieldList[] = {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90,  ncp_pct10,   1),
			_1(ncpPr,    mxKnotNum + 1),				
			_2(cpOccPr, order, N),													
			_2(cp,  cpPr,  mxKnotNum),
           	_q(cpAbruptChange, mxKnotNum),
			_1(cpCI,			   mxKnotNum, 2 ),

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

	I32 nfields = sizeof(fieldList) / sizeof(FIELD_ITEM);
	__AddSpatialDimension(io->dims[io->rowdim - 1], io->dims[io->coldim - 1], io->out.whichDimIsTime, fieldList, nfields);
	__RemoveFieldsGivenFlags_Season(opt, fieldList, nfields);
	if (opt->extra.removeSingletonDims) {
		RemoveSingltonDims(fieldList, nfields);
	}

	// Do nothing if io->q is 1L;  within the function, newly created lists are protected
	I32       nptr = __MR_ExtendFieldsToMultiVaraiteTS(fieldList, nfields, io->q); 
	VOID_PTR  out  = PROTECT(CreateStructVar(fieldList, nfields));
					 UNPROTECT(1L);
	UNPROTECT(nptr); 
	return out;


	#undef NUMARGS
    #undef NARGS
    #undef _
	#undef _1
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

static void* __BEAST2_Output_AllocMEM_Outlier(A(OPTIONS_PTR)  opt)
{
	const BEAST2_IO_PTR      io     = &opt->io;
	const BEAST2_RESULT_PTR  mat    = io->out.result;
	DATA_TYPE   dtype        = io->out.dtype; // DATA_FLOAT or DATA_DOUBLE                            
	const int   N	         = io->N;
	const int   mxKnotNum	 = opt->prior.outlierMaxKnotNum;	
	//https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                  (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->o##name }
	#define _1(name, ...)               _(name, __VA_ARGS__)  
	#define _2(name1,name2, ...)         _(name1, __VA_ARGS__),  _(name2, __VA_ARGS__)   
	#define _3(n1,n2,n3, ...)            _2(n1, n2, __VA_ARGS__),  _(n3, __VA_ARGS__)  
	#define _4(n1,n2,n3,n4, ...)         _3(n1, n2,n3, __VA_ARGS__),  _(n4, __VA_ARGS__)  
	#define _5(n1,n2,n3,n4,n5, ...)      _4(n1, n2,n3,n4, __VA_ARGS__),  _(n5, __VA_ARGS__)  
	#define _6(n1,n2,n3,n4,n5,n6, ...)   _5(n1, n2,n3,n4,n5, __VA_ARGS__),  _(n6, __VA_ARGS__)  
	//{ "time",    dtype, 2, { N,               1, 0L, 0L }, &mat->time },
	//{ "scp",	   dtype, 2, { mxKnotNumSeason, M, },        & mat->scp },


	int isMultiVariate = (io->q == 1) ? 0L : 1L;
	#define _q(name, ...)                 {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->o##name, isMultiVariate }
	#define _q1(name, ...)                _q(name, __VA_ARGS__)  
	#define _q2(name1,name2, ...)         _q(name1, __VA_ARGS__),  _q(name2, __VA_ARGS__)   
	#define _q3(n1,n2,n3, ...)            _q2(n1, n2, __VA_ARGS__),  _q(n3, __VA_ARGS__)  
	#define _q4(n1,n2,n3,n4, ...)         _q3(n1, n2,n3, __VA_ARGS__),  _q(n4, __VA_ARGS__)  
	#define _q5(n1,n2,n3,n4,n5, ...)      _q4(n1, n2,n3,n4, __VA_ARGS__),  _q(n5, __VA_ARGS__)  
	#define _q6(n1,n2,n3,n4,n5,n6, ...)   _q5(n1, n2,n3,n4,n5, __VA_ARGS__),  _q(n6, __VA_ARGS__) 


	FIELD_ITEM  fieldList[ ]= {
			_5(ncp, ncp_median, ncp_mode, ncp_pct90, ncp_pct10,  1),
			_1(ncpPr,   mxKnotNum + 1),
			_1(cpOccPr, N),

			_2(cp, cpPr,  mxKnotNum),
			_1(cpCI,       mxKnotNum, 2),

			_q2(Y, SD, N),
			_q(CI,     N,2),

			_q2(pos_ncp,     neg_ncp,      1),
			_q2(pos_ncpPr,   neg_ncpPr,    mxKnotNum + 1),
			_q2(pos_cpOccPr, neg_cpOccPr,  N),
			_q4(pos_cp,      neg_cp, pos_cpPr, neg_cpPr, mxKnotNum),
			_q2(pos_cpCI,    neg_cpCI,     mxKnotNum,2),
	};

	I32 nfields = sizeof(fieldList) / sizeof(FIELD_ITEM);
	__AddSpatialDimension(io->dims[io->rowdim - 1], io->dims[io->coldim - 1], io->out.whichDimIsTime, fieldList, nfields);
	__RemoveFieldsGivenFlags_Outlier(opt, fieldList, nfields);
	if (opt->extra.removeSingletonDims) {
		RemoveSingltonDims(fieldList, nfields);
	}

	// Do nothing if io->q is 1L;  within the function, newly created lists are protected
	I32       nptr = __MR_ExtendFieldsToMultiVaraiteTS(fieldList, nfields, io->q); 
	VOID_PTR  out  = PROTECT(CreateStructVar(fieldList, nfields));
					 UNPROTECT(1L);
	UNPROTECT(nptr); 
	return out;


	#undef NUMARGS
    #undef NARGS
    #undef _
    #undef _1
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



void* BEAST2_Output_AllocMEM(A(OPTIONS_PTR)  opt)  {	

	//Moved from glue_code to accomodate multivriate time series
	if (opt->io.out.result) {
		// Added to accomodate the beast_thread() for Win GUI where Output_AllocMem is called each time a thread is created.
		free(opt->io.out.result);
	}
	opt->io.out.result = malloc(sizeof(BEAST2_RESULT) * opt->io.q);
	memset(opt->io.out.result, 0, sizeof(BEAST2_RESULT) * opt->io.q);

	const BEAST2_IO_PTR      io     = &opt->io;
	const BEAST2_RESULT_PTR  mat    = io->out.result;
	DATA_TYPE                dtype  = io->out.dtype; // DATA_FLOAT or DATA_DOUBLE                            
 
	//https://stackoverflow.com/questions/2124339/c-preprocessor-va-args-number-of-arguments
	#define NUMARGS(...)                (sizeof((int[]){__VA_ARGS__})/sizeof(int))
	#define NARGS(...)                  (sizeof((int[]){0, ##__VA_ARGS__})/sizeof(int)-1)
	#define _(name, ...)                {#name, dtype, NUMARGS(__VA_ARGS__),{__VA_ARGS__},(void ** )&mat->tY##name }
	#define _1(name, ...)               _(name, __VA_ARGS__)  

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasDummyCmpnt   = opt->prior.basisType[0] == DUMMYID;
	I08 hasTrendCmpnt   = 1;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;

	I32       nprt  = 0;
	VOID_PTR  trend = NULL, season=NULL, outlier=NULL;
	if (hasTrendCmpnt)   { trend   = PROTECT(__BEAST2_Output_AllocMEM_Trend(opt));    nprt++; }
	if (hasSeasonCmpnt)  { season  = PROTECT(__BEAST2_Output_AllocMEM_Season(opt));   nprt++; }
	if (hasOutlierCmpnt) { outlier = PROTECT(__BEAST2_Output_AllocMEM_Outlier(opt));  nprt++;}
	
	int   ROW, COL;
	if (io->ndim == 1 || io->ndim==2) {
		if (io->ndim==1 || (io->ndim == 2 && io->out.whichDimIsTime==1)) 			
			ROW = 1,               COL = io->numOfPixels;
		else
			ROW = io->numOfPixels, COL = 1;
	} else {
		switch (io->meta.whichDimIsTime) {
		case 1:		ROW = io->dims[1], COL = io->dims[2]; break;
		case 2:		ROW = io->dims[0], COL = io->dims[2]; break;
		case 3:		ROW = io->dims[0], COL = io->dims[1]; break;
		}	
	}

	//#define __DEBUG__  
	I32PTR whichDimIsTime;
	I32PTR nrows, ncols;
	const  int         N  = io->N;
	const  int         q  = io->q;
	int    isMultiVariate = (q == 1) ? 0L : 1L;
	FIELD_ITEM  fieldList[ ] ={
		{"time",      dtype,		2,   {N,1,},      &mat->time},				
		{"data",      dtype,		1,   {N,  },      &mat->data, .extra= isMultiVariate}, //Needed to be changed to reflect the output dim
		{"marg_lik",  dtype,		2,   {ROW,COL,},  &mat->marg_lik},
		#ifdef __DEBUG__
        {"R2",      dtype,		2,   {(N + 7) / 8 * 8,300,},  &mat->R2},
		{"RMSE",    dtype,		2,   {(N + 7) / 8 * 8,300,},  &mat->RMSE},
	    #else
		{"R2",        dtype,		2,   {ROW,COL,},  &mat->R2,  .extra = isMultiVariate},
		{"RMSE",      dtype,		2,   {ROW,COL,},  &mat->RMSE,.extra = isMultiVariate},
		#endif
	    {"sig2",      dtype,		2,   {q,q,},       &mat->sig2}, //Needed to be changed to reflect the output dim
		//  {"sig2",      dtype,		4,   {ROW,COL,q,q},  &mat->sig2},
		{"trend",     DATA_STRUCT,	0,	 {0,},       (void**)trend},
		{"season",    DATA_STRUCT,	0,	 {0,},       (void**)season},
		{"outlier",   DATA_STRUCT,	0,	 {0,},       (void**)outlier},
		{"whichOutDimIsTime", DATA_INT32,  2,  {1,1,},   &whichDimIsTime },
		{"nrows",             DATA_INT32,  2,  {1,1,},   &nrows },
		{"ncols",             DATA_INT32,  2,  {1,1,},   &ncols },
	};
	I32    nfields = sizeof(fieldList) / sizeof(FIELD_ITEM);

	int sig2_index = 6 - 1;
	__AddSpatialDimension(io->dims[io->rowdim - 1], io->dims[io->coldim - 1], io->out.whichDimIsTime, &fieldList[sig2_index], 1L); //1L means only one field to be adjusted
	// Adjust the dimension of out$data based on the input dimesion and the outputtimedim
	// to make it compatible with whichOutputDimIsTime.
	// Change the data field by adding the extra dims and swithcing for the outputTimeDim
	__AddSpatialDimension(io->dims[io->rowdim - 1], io->dims[io->coldim - 1], io->out.whichDimIsTime, &fieldList[1], 1L);
	if (!opt->extra.dumpInputData) {  		 
		RemoveField(fieldList, nfields, "data"); mat->data = NULL;
	}
	if (opt->extra.removeSingletonDims) {
		RemoveSingltonDims(fieldList, nfields);
	}

	VOID_PTR  out;
	// Extend it if and only if q >1
	I32       nptr1 = __MR_ExtendFieldsToMultiVaraiteTS(fieldList, nfields, io->q);  nprt += nptr1;
	if (ROW * COL == 1) {
		// No need to inclulde whichTimeDim and ncols and nrows if there is only one ts
		out = PROTECT(CreateStructVar(fieldList, nfields-3));                    nprt++;
	} else {
		out = PROTECT(CreateStructVar(fieldList, nfields));                      nprt++;		
		whichDimIsTime[0] = io->out.whichDimIsTime;
		nrows[0] = ROW;
		ncols[0] = COL;
	}
	 

	AddStringAttribute(out,  "class",          "beast");	
	if (hasSeasonCmpnt && !hasDummyCmpnt) AddStringAttribute(out,  "season_type",  "harmonic");
	if (hasSeasonCmpnt &&  hasDummyCmpnt) AddStringAttribute(out,  "season_type",  "dummy");
	if (!hasSeasonCmpnt )                 AddStringAttribute(out,  "season_type",  "none");
	//AddIntegerAttribute(out, "hasOutlier", hasOutlierCmpnt);	
 
	if (io->T.out.asDailyTS) {
		for (int i = 0; i < N; i++) {
			mat->time[i]=FracYear_from_DateNum(io->meta.startTime+ io->meta.deltaTime*i);
		} 
	}	else {
		f32_seq(mat->time, io->meta.startTime, io->meta.deltaTime, N);
	}    

	if (dtype == DATA_DOUBLE)  f32_to_f64_inplace(mat->time, N);
	//r_ippsMulC_32f_I(io->meta.deltaTime, mat->time, N);
	//r_ippsSubC_32f_I(-(io->meta.startTime - io->meta.deltaTime), mat->time, N);
 

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

void  BEAST2_Result_AllocMEM(A(RESULT_PTR)  result, A(OPTIONS_PTR)  opt, MemPointers* _restrict MEM)
{	
	const I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	const I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	const I08 hasTrendCmpnt   = 1;
	const I08 hasAlways       = 1;

	const I32 N  = opt->io.N;
	const I32 q  = opt->io.q;
	const I32 Nq = N * q;

	const I32 seasonMaxKnotNum  = opt->prior.seasonMaxKnotNum;
	const I32 trendMaxKnotNum   = opt->prior.trendMaxKnotNum;
	const I32 outlierMaxKnotNum = opt->prior.outlierMaxKnotNum;

	*result = (BEAST2_RESULT){0, }; //memset result to zero

	MemNode nodes[250];
	int     nid    = 0;

	result->time   =  NULL; // No need to allocate mem for time bcz Pits values were directly set in the mat->time

	nodes[nid++] = (MemNode){ &result->marg_lik,	 sizeof(F32) * 1,      .align = 64 };
	nodes[nid++] = (MemNode){ &result->sig2,	     sizeof(F32) * q * q,  .align = 4 };
	nodes[nid++] = (MemNode){ &result->R2,	         sizeof(F32) * q  ,    .align = 4 };
	nodes[nid++] = (MemNode){ &result->RMSE,	     sizeof(F32) * q  ,    .align = 4 };
 
	if (hasSeasonCmpnt) {
		nodes[nid++] = (MemNode){ &result->sncp,	      sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sncp_median,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sncp_mode,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sncp_pct90,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sncp_pct10,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sncpPr,	      sizeof(I32) * (seasonMaxKnotNum + 1), .align = 4 };
		nodes[nid++] = (MemNode){ &result->scpOccPr,	  sizeof(I32) * N,       .align = 4 };
		nodes[nid++] = (MemNode){ &result->sY,	          sizeof(I32) * Nq,      .align = 64 };
		nodes[nid++] = (MemNode){ &result->sSD,	          sizeof(I32) * Nq,      .align = 64 };
 
		if (opt->extra.computeSeasonOrder)
			nodes[nid++] = (MemNode){ &result->sorder,	          sizeof(U32) * N,      .align = 64 };			
		if (opt->extra.computeSeasonAmp) {
			nodes[nid++] = (MemNode){ &result->samp,	          sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
			nodes[nid++] = (MemNode){ &result->sampSD,	          sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
		}
	}
	/////////////////////////////////////////////////////////////

	if (hasTrendCmpnt) {
		nodes[nid++] = (MemNode){ &result->tncp,	      sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tncp_median,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tncp_mode,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tncp_pct90,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tncp_pct10,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tncpPr,	      sizeof(I32) * (trendMaxKnotNum + 1), .align = 4 };
		nodes[nid++] = (MemNode){ &result->tcpOccPr,	  sizeof(I32) * N,       .align = 4 };
		nodes[nid++] = (MemNode){ &result->tY,	          sizeof(I32) * Nq,      .align = 64 };
		nodes[nid++] = (MemNode){ &result->tSD,	          sizeof(I32) * Nq,      .align = 64 };

		if (opt->extra.computeTrendOrder)
			nodes[nid++] = (MemNode){ &result->torder,	          sizeof(U32) * N,      .align = 64 };
		if (opt->extra.computeTrendSlope) {
			nodes[nid++] = (MemNode){ &result->tslp,	          sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
			nodes[nid++] = (MemNode){ &result->tslpSD,	          sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
			nodes[nid++] = (MemNode){ &result->tslpSgnPosPr,	  sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
			nodes[nid++] = (MemNode){ &result->tslpSgnZeroPr,	  sizeof(F32) * Nq,      .align = 64 };  //NEWLY ADDED
		} 
	}
	/////////////////////////////////////////////////////////////

	if (hasOutlierCmpnt) {
		nodes[nid++] = (MemNode){ &result->oncp,	      sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oncp_median,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oncp_mode,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oncp_pct90,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oncp_pct10,	  sizeof(F32) * 1,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oncpPr,	      sizeof(I32) * (outlierMaxKnotNum + 1), .align = 4 };
		nodes[nid++] = (MemNode){ &result->ocpOccPr,	  sizeof(I32) * N,       .align = 4 };
		nodes[nid++] = (MemNode){ &result->oY,	          sizeof(I32) * Nq,      .align = 64 };
		nodes[nid++] = (MemNode){ &result->oSD,	          sizeof(I32) * Nq,      .align = 64 };
	}

	/////////////////////////////////////////////////////////////
	if (opt->extra.computeCredible){
		if (hasSeasonCmpnt)  nodes[nid++] = (MemNode){ &result->sCI,    sizeof(F32) * Nq * 2,      .align = 4 };
		if (hasTrendCmpnt)   nodes[nid++] = (MemNode){ &result->tCI,    sizeof(F32) * Nq * 2,      .align = 4 }; 
		if (hasOutlierCmpnt) nodes[nid++] = (MemNode){ &result->oCI,    sizeof(F32) * Nq * 2,      .align = 4 }; 
	}


	if (opt->extra.computeSeasonChngpt && hasSeasonCmpnt) {
		nodes[nid++] = (MemNode){ &result->scp,				sizeof(U32) * seasonMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->scpPr,			sizeof(U32) * seasonMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->scpAbruptChange,	sizeof(U32) * seasonMaxKnotNum * q,   .align = 4 };
		nodes[nid++] = (MemNode){ &result->scpCI,			sizeof(U32) * seasonMaxKnotNum * 2,   .align = 4 };
  	}	
	if (opt->extra.computeTrendChngpt && hasTrendCmpnt) {
		nodes[nid++] = (MemNode){ &result->tcp,				sizeof(U32) * trendMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tcpPr,			sizeof(U32) * trendMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tcpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,   .align = 4 };
		nodes[nid++] = (MemNode){ &result->tcpCI,			sizeof(U32) * trendMaxKnotNum * 2,   .align = 4 };
	}
	if (opt->extra.computeOutlierChngpt && hasOutlierCmpnt) {
		nodes[nid++] = (MemNode){ &result->ocp,				sizeof(U32) * outlierMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->ocpPr,			sizeof(U32) * outlierMaxKnotNum,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->ocpCI,			sizeof(U32) * outlierMaxKnotNum * 2,   .align = 4 };
	}
		
	/////////////////////////////////////////////////////////////	

	if (opt->extra.tallyPosNegSeasonJump && hasSeasonCmpnt) {
		nodes[nid++] = (MemNode){ &result->spos_ncp,				sizeof(F32) * 1*q,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_ncp,				sizeof(F32) * 1 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->spos_ncpPr,				sizeof(I32) * (seasonMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_ncpPr,				sizeof(I32) * (seasonMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->spos_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };

		nodes[nid++] = (MemNode){ &result->spos_cp,				sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_cp,				sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->spos_cpPr,			sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_cpPr,			sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->spos_cpAbruptChange,	sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_cpAbruptChange,	sizeof(U32) * seasonMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->spos_cpCI,			sizeof(U32) * seasonMaxKnotNum * 2*q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->sneg_cpCI,			sizeof(U32) * seasonMaxKnotNum * 2*q,     .align = 4 };

	}
 

	if (opt->extra.tallyPosNegTrendJump && hasTrendCmpnt) {
		nodes[nid++] = (MemNode){ &result->tpos_ncp,				sizeof(F32) * 1*q,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_ncp,				sizeof(F32) * 1 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tpos_ncpPr,				sizeof(I32) * (trendMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_ncpPr,				sizeof(I32) * (trendMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tpos_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };

		nodes[nid++] = (MemNode){ &result->tpos_cp,				sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_cp,				sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tpos_cpPr,			sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_cpPr,			sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tpos_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tpos_cpCI,			sizeof(U32) * trendMaxKnotNum * 2*q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tneg_cpCI,			sizeof(U32) * trendMaxKnotNum * 2*q,     .align = 4 };
	}
		 

	if (opt->extra.tallyIncDecTrendJump && hasTrendCmpnt) {
		nodes[nid++] = (MemNode){ &result->tinc_ncp,				sizeof(F32) * 1 * q,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_ncp,				sizeof(F32) * 1 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tinc_ncpPr,				sizeof(I32) * (trendMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_ncpPr,				sizeof(I32) * (trendMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tinc_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };

		nodes[nid++] = (MemNode){ &result->tinc_cp,				sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_cp,				sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tinc_cpPr,			sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_cpPr,			sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tinc_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tinc_cpCI,			sizeof(U32) * trendMaxKnotNum * 2 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->tdec_cpCI,			sizeof(U32) * trendMaxKnotNum * 2 * q,     .align = 4 };
	}
  

	if (opt->extra.tallyPosNegOutliers && hasOutlierCmpnt) {
		nodes[nid++] = (MemNode){ &result->opos_ncp,				sizeof(F32) * 1 * q,      .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_ncp,				sizeof(F32) * 1 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->opos_ncpPr,				sizeof(I32) * (outlierMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_ncpPr,				sizeof(I32) * (outlierMaxKnotNum + 1) * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->opos_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_cpOccPr,			sizeof(I32) * Nq,     .align = 4 };

		nodes[nid++] = (MemNode){ &result->opos_cp,				sizeof(U32) * outlierMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_cp,				sizeof(U32) * outlierMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->opos_cpPr,			sizeof(U32) * outlierMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_cpPr,			sizeof(U32) * outlierMaxKnotNum * q,     .align = 4 };
		//nodes[nid++] = (MemNode){ &result->tpos_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		//nodes[nid++] = (MemNode){ &result->tneg_cpAbruptChange,	sizeof(U32) * trendMaxKnotNum * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->opos_cpCI,			sizeof(U32) * outlierMaxKnotNum * 2 * q,     .align = 4 };
		nodes[nid++] = (MemNode){ &result->oneg_cpCI,			sizeof(U32) * outlierMaxKnotNum * 2 * q,     .align = 4 };
	}
		
	if (opt->extra.dumpInputData) {
		nodes[nid++] = (MemNode){ &result->data,				sizeof(F32) * Nq,      .align = 4 };		
	}

	nodes[nid++] = (MemNode){ NULL,		};
 
	MEM->alloclist(MEM, nodes, AggregatedMemAlloc, NULL);
}

void  BEAST2_Result_FillMEM(A(RESULT_PTR)  result, A(OPTIONS_PTR)  opt, const F32 nan)
{
	const I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	const I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	const I08 hasTrendCmpnt   = 1;
	const I08 hasAlways       = 1;

	const I32 N = opt->io.N;
	const I32 q = opt->io.q;
	const I32 Nq = N * q;
	const I32 seasonMaxKnotNum  = opt->prior.seasonMaxKnotNum;
	const I32 trendMaxKnotNum   = opt->prior.trendMaxKnotNum;
	const I32 outlierMaxKnotNum = opt->prior.outlierMaxKnotNum;

	A(EXTRA_PTR) flag = &(opt->extra);

	//No need to fill "data" because its values will be copied directly from yInfo.Y
	// Its value won't be set to nans even if beast fails.
	//f32_fill_val(nan, result->data, N);

	*result->marg_lik = nan;
	f32_fill_val(nan, result->sig2, q*q);
	f32_fill_val(nan, result->R2,   q);
	f32_fill_val(nan, result->RMSE, q);
	if (hasSeasonCmpnt) {
			*result->sncp         = nan;
			*result->sncp_median  = nan;
			*result->sncp_mode    = nan;
			*result->sncp_pct90   = nan;
			*result->sncp_pct10   = nan;
			f32_fill_val(nan, result->sncpPr,  seasonMaxKnotNum + 1);
			f32_fill_val(nan, result->scpOccPr, N);
			f32_fill_val(nan, result->sY,		Nq);
			f32_fill_val(nan, result->sSD,		Nq);
			if (flag->computeSeasonOrder)  //result->sorder is an int32ptr
				f32_fill_val(nan, result->sorder, N); //https://en.wikipedia.org/wiki/Signed_zero#:~:text=In%20IEEE%20754%20binary%20floating,sign%20bit%20set%20to%20one.
		
			if (flag->computeSeasonAmp) {
				f32_fill_val(nan, result->samp,   Nq);
				f32_fill_val(nan, result->sampSD, Nq);
			}
	}

	////////////////////////////////////////////
	if (hasTrendCmpnt) {
			*result->tncp         = nan;
			*result->tncp_median  = nan;
			*result->tncp_mode    = nan;
			*result->tncp_pct90   = nan;
			*result->tncp_pct10   = nan;
			f32_fill_val(nan, result->tncpPr, trendMaxKnotNum + 1);
			f32_fill_val(nan, result->tcpOccPr, N);
			f32_fill_val(nan, result->tY,     Nq);
			f32_fill_val(nan, result->tSD,    Nq);
			if (flag->computeTrendOrder)
				f32_fill_val(nan, result->torder, N);
			if (flag->computeTrendSlope) {
				f32_fill_val(nan, result->tslp,       Nq);
				f32_fill_val(nan, result->tslpSD,     Nq);
				f32_fill_val(nan, result->tslpSgnPosPr, Nq);
				f32_fill_val(nan, result->tslpSgnZeroPr, Nq);
			}
	}

	////////////////////////////////////////////
	if (hasOutlierCmpnt) {
		*result->oncp = nan;
		*result->oncp_median = nan;
		*result->oncp_mode = nan;
		*result->oncp_pct90 = nan;
		*result->oncp_pct10 = nan;
		f32_fill_val(nan, result->oncpPr, outlierMaxKnotNum + 1);
		f32_fill_val(nan, result->ocpOccPr, N);
		f32_fill_val(nan, result->oY,       Nq);
		f32_fill_val(nan, result->oSD,      Nq);
	}

	////////////////////////////////////////////
	if (flag->computeCredible) 	{
		if (hasSeasonCmpnt)  f32_fill_val(nan, result->sCI, 2*Nq);
		if (hasTrendCmpnt)   f32_fill_val(nan, result->tCI, 2*Nq);
		if (hasOutlierCmpnt) f32_fill_val(nan, result->oCI, 2*Nq);
	}

	if (flag->computeSeasonChngpt && hasSeasonCmpnt) 	{
		f32_fill_val(nan, result->scp,		seasonMaxKnotNum);
		f32_fill_val(nan, result->scpPr,	seasonMaxKnotNum);
		f32_fill_val(nan, result->scpAbruptChange, seasonMaxKnotNum*q);
		f32_fill_val(nan, result->scpCI, 2*seasonMaxKnotNum); 
	}

	if (flag->computeTrendChngpt && hasTrendCmpnt) {
		f32_fill_val(nan, result->tcp,		trendMaxKnotNum );
		f32_fill_val(nan, result->tcpPr,	trendMaxKnotNum );
		f32_fill_val(nan, result->tcpAbruptChange, trendMaxKnotNum * q);
		f32_fill_val(nan, result->tcpCI,	2* trendMaxKnotNum );
 
	}

	if (flag->computeOutlierChngpt&& hasOutlierCmpnt) {
		f32_fill_val(nan, result->ocp,		outlierMaxKnotNum );
		f32_fill_val(nan, result->ocpPr,	outlierMaxKnotNum );
		f32_fill_val(nan, result->ocpCI, 2* outlierMaxKnotNum );
	}

	if (flag->tallyPosNegSeasonJump && hasSeasonCmpnt)  {

		f32_fill_val(nan, result->spos_ncp, 1*q);
		f32_fill_val(nan, result->sneg_ncp, 1 * q);

		f32_fill_val(nan, result->spos_ncpPr, (seasonMaxKnotNum + 1) * q);
		f32_fill_val(nan, result->sneg_ncpPr, (seasonMaxKnotNum + 1) * q);

		f32_fill_val(nan, result->spos_cpOccPr, Nq);
		f32_fill_val(nan, result->sneg_cpOccPr, Nq);

		f32_fill_val(nan, result->spos_cp, seasonMaxKnotNum * q);
		f32_fill_val(nan, result->sneg_cp, seasonMaxKnotNum * q);
		f32_fill_val(nan, result->spos_cpPr, seasonMaxKnotNum * q);
		f32_fill_val(nan, result->sneg_cpPr, seasonMaxKnotNum * q);
		f32_fill_val(nan, result->spos_cpAbruptChange, seasonMaxKnotNum * q);
		f32_fill_val(nan, result->sneg_cpAbruptChange, seasonMaxKnotNum * q);

		f32_fill_val(nan, result->spos_cpCI, 2*seasonMaxKnotNum * q);
		f32_fill_val(nan, result->sneg_cpCI, 2 * seasonMaxKnotNum * q);
	}

	if (flag->tallyPosNegTrendJump && hasTrendCmpnt) {
		f32_fill_val(nan, result->tpos_ncp, 1 * q);
		f32_fill_val(nan, result->tneg_ncp, 1 * q);

		f32_fill_val(nan, result->tpos_ncpPr, (trendMaxKnotNum + 1) * q);
		f32_fill_val(nan, result->tneg_ncpPr, (trendMaxKnotNum + 1) * q);

		f32_fill_val(nan, result->tpos_cpOccPr, Nq);
		f32_fill_val(nan, result->tneg_cpOccPr, Nq);

		f32_fill_val(nan, result->tpos_cp, trendMaxKnotNum * q);
		f32_fill_val(nan, result->tneg_cp, trendMaxKnotNum* q);

		f32_fill_val(nan, result->tpos_cpPr, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tneg_cpPr, trendMaxKnotNum* q);

		f32_fill_val(nan, result->tpos_cpAbruptChange, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tneg_cpAbruptChange, trendMaxKnotNum* q);

		f32_fill_val(nan, result->tpos_cpCI, 2* trendMaxKnotNum * q);
		f32_fill_val(nan, result->tneg_cpCI, 2 * trendMaxKnotNum * q);
	}


	if (flag->tallyIncDecTrendJump && hasTrendCmpnt) {

		f32_fill_val(nan, result->tinc_ncp, 1 * q);
		f32_fill_val(nan, result->tdec_ncp, 1 * q);

		f32_fill_val(nan, result->tinc_ncpPr, (trendMaxKnotNum + 1) * q);
		f32_fill_val(nan, result->tdec_ncpPr, (trendMaxKnotNum + 1)* q);

		f32_fill_val(nan, result->tinc_cpOccPr, Nq);
		f32_fill_val(nan, result->tdec_cpOccPr, Nq);

		f32_fill_val(nan, result->tinc_cp, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tdec_cp, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tinc_cpPr, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tdec_cpPr, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tinc_cpAbruptChange, trendMaxKnotNum* q);
		f32_fill_val(nan, result->tdec_cpAbruptChange, trendMaxKnotNum* q);

		f32_fill_val(nan, result->tinc_cpCI, 2* trendMaxKnotNum);
		f32_fill_val(nan, result->tdec_cpCI, 2 * trendMaxKnotNum); 
 
	}

	if (flag->tallyPosNegOutliers && hasOutlierCmpnt) {
		f32_fill_val(nan, result->opos_ncp, 1 * q);
		f32_fill_val(nan, result->oneg_ncp, 1 * q);

		f32_fill_val(nan, result->opos_ncpPr, (outlierMaxKnotNum + 1)* q);
		f32_fill_val(nan, result->opos_ncpPr, (outlierMaxKnotNum + 1) * q);

		f32_fill_val(nan, result->opos_cpOccPr, Nq);
		f32_fill_val(nan, result->oneg_cpOccPr, Nq);

		f32_fill_val(nan, result->opos_cp, outlierMaxKnotNum* q);
		f32_fill_val(nan, result->oneg_cp, outlierMaxKnotNum* q);
		f32_fill_val(nan, result->opos_cpPr, outlierMaxKnotNum* q);
		f32_fill_val(nan, result->oneg_cpPr, outlierMaxKnotNum * q);

		f32_fill_val(nan, result->opos_cpCI, 2 * outlierMaxKnotNum * q);
		f32_fill_val(nan, result->oneg_cpCI, 2 * outlierMaxKnotNum * q);

	}
}


 

#include "abc_000_warning.h"