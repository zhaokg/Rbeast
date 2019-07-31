#pragma once
#include <stdio.h>
#include "abc_001_config.h"
#include "abc_blas_lapack_lib.h"
#include "abc_mem.h"
#include "abc_datatype.h" 
#if M_INTERFACE==1
	#define  r_printf(...) mexPrintf(__VA_ARGS__)
	#define  r_error(x)    mexErrMsgTxt(x)
	#define  r_malloc(x)   mxMalloc(x) 
	#define  r_free(x)     mxFree(x)
#elif R_INTERFACE==1
	#define  r_printf(...)  Rprintf(__VA_ARGS__)
	#define  r_error(x)     error(x)
	#define  r_malloc(x)    Calloc(x,char)  
	#define  r_free(x)      Free(x) 
#endif
extern void    fill_1_to_N(rF32PTR p,int N);
extern void    quickSortD(F32PTR arr,int32_t * _restrict index,int32_t low,int32_t high);
extern void    quickSortA(F32PTR arr,int32_t * _restrict index,int32_t low,int32_t high);
extern void    transpose_inplace(F32PTR m,int w,int h);
extern int32_t strcicmp(char const * _restrict a,char const * _restrict b);
extern float   determine_period(F32PTR data,int32_t N,float  omissionValue);
extern int32_t find_changepoint(F32PTR prob,F32PTR mem,float threshold,I32PTR cpt,F32PTR cptCI,int32_t N,int32_t minSepDist,int32_t maxCptNumber);
extern float fastlog(float x);
extern float sum_log_diag(rF32PTR p,rI32 K);
extern float fastexp(float x);
extern float fast_sqrt(float x);
static INLINE void normalize(rF32PTR ptr,int N)
{
	float mean,std;
	r_ippsMeanStdDev_32f(ptr,N,&mean,&std,ippAlgHintAccurate);
	r_ippsSubC_32f_I(mean,ptr,N);
	r_cblas_sscal(N,1.f/std,ptr,1L);
}
static INLINE void normalize_x_factor(F32PTR ptr,int N,float factor)
{
	float mean,std;
	r_ippsMeanStdDev_32f(ptr,N,&mean,&std,ippAlgHintAccurate);
	r_ippsSubC_32f_I(mean,ptr,N);
	r_cblas_sscal(N,1/std*factor,ptr,1);
}
