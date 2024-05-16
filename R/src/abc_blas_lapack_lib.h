#pragma once


#include "abc_000_macro.h"
#if (MYMAT_LIBRARY == 1) 
	#include "abc_blas_lapack_myl.h"
#elif (MKL_LIBRARY ==1)
	#include "abc_blas_lapack_mkl.h"
#endif

#if (MYMAT_LIBRARY == 1) || (MATLAB_LIBRARY==1)

#include "abc_vec.h" //tranpose_inplace
#include "abc_tranpose.h" //tranpose_inplace
#define r_mkl_simatcopy(C, T, N, numCISample, alpha, Matrix, N1, numCISample1) \
	    i32_transpose_inplace_prev_two_ends(Matrix, (U64) N,(U64) numCISample )

#endif
