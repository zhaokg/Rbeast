#pragma once
#include "abc_000_macro.h"
#if (MYMAT_LIBRARY==1) 
	#include "abc_blas_lapack_myl.h"
#elif (MKL_LIBRARY==1)
	#include "abc_blas_lapack_mkl.h"
#endif
#if (MYMAT_LIBRARY==1)||(MATLAB_LIBRARY==1)
#define r_mkl_simatcopy(C,T,N,numCISample,alpha,Matrix,N1,numCISample1) transpose_inplace(Matrix,N,numCISample )
#endif
