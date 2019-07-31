#include "abc_001_config.h"
DISABLE_MANY_WARNINGS
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "beast_common.h"
#include "beast_interface.h"
#if !R_RELEASE
#include "trend_interface.h"
#endif
#if defined(WIN64_OS) 
extern void DllExport WinMainDemoST(Options * option,RESULT * result);
#endif
#if R_INTERFACE==1
#include "abc_R_util.h"
SEXP DllExport beast2(SEXP Y,SEXP opt)
{	
	char	missing[31];
	Options beastOption;
	if (!R_check_input(Y,opt)||
		!R_read_input(&beastOption,Y,opt,missing)||
		!check_options(&beastOption,missing))
		return R_NilValue;
	print_options(&beastOption);
	RESULT	result;
	SEXP	ANS;
	PROTECT(ANS=R_allocate_output(&result,&beastOption));
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	beastST2();
	UNPROTECT(1);
	return ANS;
}
#if !R_RELEASE 
SEXP DllExport trend2(SEXP Y,SEXP opt)
{
	char	missing[31];
	Options beastOption;
	if (!R_trend_check_input(Y,opt) &&
		!R_trend_read_input(&beastOption,Y,opt,missing) &&
		!trend_check_options(&beastOption,missing))
		return R_NilValue;
	print_options_trend(&beastOption);
	RESULT	result;
	SEXP	ANS;
	PROTECT(ANS=R_trend_allocate_output(&result,&beastOption));
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	beastTrend2();
	UNPROTECT(1);
	return ANS;
}
#endif
#if defined(WIN64_OS) 
SEXP DllExport GUI_beast(SEXP Y,SEXP opt)
{
	char	missing[31];
	Options beastOption;
	if (!R_check_input(Y,opt)||
		!R_read_input(&beastOption,Y,opt,missing)||
		!check_options(&beastOption,missing))
		return R_NilValue;
	print_options(&beastOption);
	RESULT	result;
	SEXP	ANS;
	PROTECT(ANS=R_allocate_output(&result,&beastOption));
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	WinMainDemoST(GLOBAL_OPTIONS,GLOBAL_RESULT);
	UNPROTECT(1);
	return ANS;
}
#endif
#if !R_RELEASE && defined(WIN64_OS) 
SEXP DllExport GUI_trend(SEXP Y,SEXP opt)
{
	if (!R_trend_check_input(Y,opt))				return R_NilValue;
	char	missing[31];
	Options beastOption;
	R_trend_read_input(&beastOption,Y,opt,missing);
	if (!trend_check_options(&beastOption,missing)) return R_NilValue;
	print_options_trend(&beastOption);
	RESULT	result;
	SEXP	ANS;
	PROTECT(ANS=R_trend_allocate_output(&result,&beastOption));
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	WinMainDemoTrend(GLOBAL_OPTIONS,GLOBAL_RESULT);
	UNPROTECT(1);
	return ANS;
}
SEXP DllExport sbm2(SEXP Y,SEXP opt)
{
	if (!R_sbm_check_input(Y,opt))		 			return R_NilValue;
	char	missing[31];
	Options beastOption;
	R_sbm_read_input(&beastOption,Y,opt,missing);
	if (!sbm_check_options(&beastOption,missing))	return R_NilValue;
	print_options(&beastOption);
	RESULT	result;
	SEXP	ANS;
	PROTECT(ANS=R_allocate_output(&result,&beastOption));
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	sbm();
	UNPROTECT(1);
	return ANS;
}
#endif
static char fileID UNUSED_DECORATOR='R';
#ifdef __GNU__
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xresource.h>
#include <X11/Xlocale.h>
#include <X11/Xatom.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#define CALLDEF(name,n) {#name,(DL_FUNC) &name,n}
#if defined(WIN64_OS) 
#endif
static const R_CallMethodDef CallEntries[]={
			CALLDEF(beast2,2),
#if !R_RELEASE 
			CALLDEF(trend2,2),
#endif
#if defined(WIN64_OS) 
			CALLDEF(GUI_beast,2),
	#if !R_RELEASE 
			CALLDEF(GUI_trend,2),
	#endif
#endif
			{ NULL,NULL,0 } 
};
void  R_init_Rbeast(DllInfo *dll)
{
	R_registerRoutines(dll,NULL,CallEntries,NULL,NULL);
	R_useDynamicSymbols(dll,FALSE);
}
#elif M_INTERFACE==1
void DllExport beastST_multipleChain_fast(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport beastTrend_multipleChain_fast(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport mvST_multipleChain(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport SBM_ST(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport SBM_ST_BIC(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
void DllExport WinMainDemoTrend(Options * option,RESULT * result);
void DllExport WinMainDemoST(Options * option,RESULT * result);
void DllExport mexFunction(int nlhs,mxArray * _restrict plhs[],int nrhs,const mxArray * _restrict prhs[])
{
	const mxArray * S;
	S=prhs[1];
	if (!mxIsStruct(S))
	{
		r_error("Error (mexFunction):The input parameter OPTION should be a struct variable!\n");
		return;
	}
	mxArray *tmp;
	tmp=mxGetField(S,0,"algorithm");
	if (tmp==NULL)
	{
		r_error("Error (mexFunction): In the option,the algorithm parameter must be specified!\n");
		return;
	}
	Options beastOption;
	char    missing[31];
	RESULT  result;
	GLOBAL_OPTIONS=&beastOption;
	GLOBAL_RESULT=&result;
	char *algName=mxIsChar(tmp) ? mxArrayToString(tmp) : NULL;
	if (algName !=NULL)
	{
		if (strcicmp(algName,"beastST_old")==0)
		{
		}
		else if (strcicmp(algName,"beastST_multipleChain_fast")==0)
		{
		}
		else if (strcicmp(algName,"beastTrend_multipleChain_fast")==0)
		{
		}
		else if (strcicmp(algName,"mvST_multipleChain")==0)
		{
			mvST_multipleChain(nlhs,plhs,nrhs,prhs);
		}
		else if (strcicmp(algName,"SBM_ST")==0)
		{
		}
		else if (strcicmp(algName,"SBM_ST_BIC")==0)
		{
		}
		else if (strcicmp(algName,"WinMainDemoTrend")==0)
		{
			if (!M_trend_check_input(nrhs,prhs)) {r_free(algName); return;}
			M_trend_read_input(&beastOption,nrhs,prhs,missing);
			if (!trend_check_options(&beastOption,missing)) {r_free(algName); return;}
			print_options_trend(&beastOption);
			plhs[0]=M_trend_allocate_output(&result,&beastOption);
			WinMainDemoTrend(GLOBAL_OPTIONS,GLOBAL_RESULT);
		}
		else if (strcicmp(algName,"WinMainDemoST")==0)
		{
			if (!M_check_input(nrhs,prhs)) {r_free(algName); return;}
			M_read_input(&beastOption,nrhs,prhs,missing);
			if (!check_options(&beastOption,missing)) {r_free(algName); return;}
			print_options(&beastOption);			
			plhs[0]=M_allocate_output(&result,&beastOption);
			WinMainDemoST(GLOBAL_OPTIONS,GLOBAL_RESULT);
		}
		else if (strcicmp(algName,"beast2")==0)
		{
			if (!M_check_input(nrhs,prhs)||
				!M_read_input(&beastOption,nrhs,prhs,missing)||
				!check_options(&beastOption,missing)
				)
			{r_free(algName); return;}
			print_options(&beastOption);			
			plhs[0]=M_allocate_output(&result,&beastOption);
			beastST2();
		}
		else if (strcicmp(algName,"beastTrend2")==0)
		{
			if (!M_trend_check_input(nrhs,prhs)) {r_free(algName); return;}
			M_trend_read_input(&beastOption,nrhs,prhs,missing);
			if (!trend_check_options(&beastOption,missing)) {r_free(algName); return;}
			print_options_trend(&beastOption);
			plhs[0]=M_trend_allocate_output(&result,&beastOption);
			beastTrend2();
		}
		r_free(algName);
	}
	return;
}
#endif
ENABLE_MANY_WARNINGS
