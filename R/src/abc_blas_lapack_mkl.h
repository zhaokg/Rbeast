#pragma once
#include "abc_001_config.h"

#define r_mkl_simatcopy(C, T, N, numCISample, alpha, Matrix, N1, numCISample1)  \
		  mkl_simatcopy(C, T, N, numCISample, alpha, Matrix, N1, numCISample1) 

/*******************************************************/
#define r_ippsMaxIndx_32f				ippsMaxIndx_32f 
#define r_ippsMinIndx_32f				ippsMinIndx_32f 
#define r_ippsSumLn_32f					ippsSumLn_32f
#define r_ippsMeanStdDev_32f			ippsMeanStdDev_32f
#define r_vsPowx(n, a,b, r)				vsPowx(n, a, b,r)
#define r_vmsLn(n, a, r, mode)			vmsLn(n, a, r, mode)
#define r_vmsCos(n, a, r, mode)			vmsCos(n, a, r, mode)
#define r_vmsSin(n, a, r, mode)			vmsSin(n, a, r, mode)
#define r_ippsSet_8u					ippsSet_8u  
#define r_ippsLn_32f_I					ippsLn_32f_I
#define r_ippsSqrt_32f_I				ippsSqrt_32f_I
#define r_ippsMul_32f_I					ippsMul_32f_I
#define r_ippsSubC_32f_I				ippsSubC_32f_I
#define r_ippsSubCRev_32f_I				ippsSubCRev_32f_I
#define r_ippsSub_32f_I					ippsSub_32f_I
#define r_ippsSub_32f                   ippsSub_32f
#define ONE_STEP_DIFF					ippsSub_32f
#define r_ippsAdd_32f_I					ippsAdd_32f_I
#define r_ippsMul_32f					ippsMul_32f
#define r_ippsSum_32f					ippsSum_32f
#define r_ippsSum_32s_Sfs				ippsSum_32s_Sfs
#define r_ippsMulC_32f_I				ippsMulC_32f_I
#define r_LAPACKE_strtrs				LAPACKE_strtrs
#define r_LAPACKE_spotrs				LAPACKE_spotrs
#define r_LAPACKE_spotrf				LAPACKE_spotrf
#define r_cblas_strmv 					cblas_strmv 
#define r_cblas_sgemv					cblas_sgemv
#define r_cblas_sgemm					cblas_sgemm
#define r_cblas_ssymv					cblas_ssymv
#define r_ippsSet_32f(value,dst,N)	    ippsSet_32f(value,dst,N)
#define r_ippsSet_32s(value,dst,N)	    ippsSet_32s(value,dst,N)
#define DOT(N,X,Y)						cblas_sdot(N, X,1,Y,1)
#define r_cblas_sscal(N, alpha, X, incX)				cblas_sscal(N, alpha, X, incX);
#define r_cblas_scopy(N, src, incX, dst, incY)			cblas_scopy(N, src, incX, dst, incY)
#define r_ippsAddC_32s_ISfs(val, X, N, scaleFactor)		ippsAddC_32s_ISfs(val, X, N, scaleFactor)  

 
 