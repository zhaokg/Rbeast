#pragma once

#include "abc_000_macro.h"
#include "abc_vec.h"
#include "abc_mat.h"

#if defined(COMPILER_MSVC)
	#define F77__CALL(x)  x
#elif defined(COMPILER_CLANG)|| defined(COMPILER_GCC) ||defined(COMPILER_SOLARIS)  
	#define PRIMITIVE_CAT(a, ...) a##__VA_ARGS__
	#define F77__CALL(x)		PRIMITIVE_CAT(x, _)

	#define RSGEMV 			rsgemv 
	#define RSGEMM 			rsgemm 
	#define RSTRMV			rstrmv
	#define RSSYMV          rssymv

	#define RSPOTRS			rspotrs
	#define RSTRTRS 		rstrtrs 
	#define RSPOTRF 		rspotrf 
#endif

#define CblasUpper    121
#define CblasLower    122
#define CblasNonUnit  131
#define CblasUnit     132
#define CblasNoTrans  111
#define CblasTrans    112


// T Declarations of FOTRAN BLAS/LAPACK functions 
void  F77__CALL(RSGEMV)(char *trans, int *m, int *n, F32PTR alpha, F32PTR a, const int *lda, F32PTR x, int *incx, F32PTR beta, F32PTR y, int *incy, size_t xxx_Char_Len3);
void  F77__CALL(RSGEMM)(char *, char*, int*, int*, int*, F32PTR , F32PTR , int*, F32PTR , int*, F32PTR , F32PTR , int*, size_t xxx_Char_Len1, size_t xxx_Char_Len2);
void  F77__CALL(RSTRMV)(char *uplo, char *trans, char *diag, int *n, F32PTR a, int *lda, F32PTR x, int *incx, size_t xxx_Char_Len1, size_t xxx_Char_Len2, size_t xxx_Char_Len3);
void  F77__CALL(RSSYMV)(char *UPLO, int* N, F32PTR  ALPHA, F32PTR  A, int * lda, F32PTR X, int * INCX, F32PTR  beta, F32PTR Y, int * INCY, size_t xxx_Char_Len1);
void  F77__CALL(RSTRTRS)(char* uplo, char* trans, char* diag, int* n, int* nrhs, F32PTR  a, int* lda, F32PTR  b, int* ldb, int* info, size_t xxx_Char_Len1, size_t xxx_Char_Len2, size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRF)(char* uplo, int* n, F32PTR  a, int* lda, int* info, size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRS)(char* uplo, int* n, int* nrhs, F32PTR  a, int* lda, F32PTR  b, int* ldb, int* info, size_t xxx_Char_Len1);
/* Declarations of Vector-based Fortran math funcitions Vector-based Fortran math funcitions */


/****************************************************************************************/
//  Declarations of C functions 
/****************************************************************************************/

//https: //en.wikipedia.org/wiki/Inline_function
/*
//warning: 'F77__ssymv' declared 'static' but never defined [-Wunused-function]
static INLINE void  F77__ssymv(int uplo, int n, F32 alpha, F32PTR A, int lda, F32PTR x, int incx, F32 beta, F32PTR y, int incy);
static INLINE void  F77__strmv(int uplo, int trans, int diag, int n, F32PTR a, int lda, F32PTR x, int incx);
//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);

static INLINE void  F77__sgemv(int trans, int m, int n, F32 alpha, F32PTR a, int lda,F32PTR x, int incx, F32  beta, F32PTR y, int incy);
static INLINE void  F77__sgemm(int transa, int transb, int m, int n, int k, F32 alpha, F32PTR a, int lda,	F32PTR b, int ldb, F32 beta, F32PTR c, int ldc);
static INLINE int   F77__strtrs(char uplo, char trans, char diag, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb);
//LAPACKE_strtrs( int matrix_layout, char uplo, char trans, char diag,lapack_int n, lapack_int nrhs, const F32PTR  a,lapack_int lda, F32PTR  b, lapack_int ldb );
static INLINE int   F77__spotrs(char uplo, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb);
//lapack_int LAPACKE_spotrs(int matrix_layout, char uplo, lapack_int n,lapack_int nrhs, const F32PTR  a, lapack_int lda,F32PTR  b, lapack_int ldb);
static INLINE int   F77__spotrf(char uplo, int n, F32PTR  a, int lda);
//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, F32PTR  a, lapack_int lda);
*/

/****************************************************************************************/
// C interfaces to the above Fortran Math functions 
/****************************************************************************************/
 
/*

#define r_cblas_ssymv(layout, uplo, n, alpha, A, lda, x, incx, beta, y, incy) \
		F77__ssymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)
//void cblas_ssymv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const F32 alpha, const F32PTR a, const MKL_INT lda, const F32PTR x, const MKL_INT incx, const F32 beta, F32PTR y, const MKL_INT incy);
static INLINE void  F77__ssymv(int uplo, int n, F32 alpha, F32PTR A, int lda, F32PTR x, int incx, F32 beta, F32PTR y, int incy)
{
	char UP = (uplo == CblasUpper) ? 'U' : 'L';
	//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
	F77__CALL(RSSYMV)(&UP, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}


//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
#define r_cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX)  \
	        F77__strmv(Uplo, TransA, Diag, N, A, lda, X, incX)
static INLINE void  F77__strmv(int uplo, int trans, int diag, int n, F32PTR a, int lda, F32PTR x, int incx)
{
	char UP = (uplo == CblasUpper) ? 'U' : 'L';
	char TRANSPOSE = (trans == CblasTrans) ? 'T' : 'N';
	char DIAGONAL = (diag == CblasUnit) ? 'U' : 'N';
	//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
	F77__CALL(RSTRMV)(&UP, &TRANSPOSE, &DIAGONAL, &n, a, &lda, x, &incx, 1, 1, 1);
}


//LAPACKE_strtrs( int matrix_layout, char uplo, char trans, char diag,lapack_int n, lapack_int nrhs, const F32PTR  a,lapack_int lda, F32PTR  b, lapack_int ldb );
#define r_LAPACKE_strtrs(matrix_layout, uplo, trans, diag,n, nrhs,  a, lda,  b, ldb )  \
      	F77__strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb)
static INLINE int  F77__strtrs(char uplo, char trans, char diag, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb)
{
	int info;
	F77__CALL(RSTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
	return info;
}


//lapack_int LAPACKE_spotrs(int matrix_layout, char uplo, lapack_int n,lapack_int nrhs, const F32PTR  a, lapack_int lda,F32PTR  b, lapack_int ldb);
#define r_LAPACKE_spotrs(matrix_layout, uplo, n,nrhs, a, lda,b, ldb)   \
        F77__spotrs(uplo, n, nrhs, a, lda, b, ldb)
static INLINE int  F77__spotrs(char uplo, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb)
{
	int info;
	F77__CALL(RSPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
	return info;
}


//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, F32PTR  a, lapack_int lda);
#define  r_LAPACKE_spotrf(matrix_layout, uplo, n,  a, lda)            F77__spotrf(uplo, n, a, lda)
static INLINE int  F77__spotrf(char uplo, int n, F32PTR  a, int lda) {
	int info;
	F77__CALL(RSPOTRF)(&uplo, &n, a, &lda, &info, 1);
	return info;
}



#define r_cblas_sgemv(CBLAS_LAYOUT, CBLAS_TRANSPOSE,  M, N, alpha, A, lda, X, incX, beta, Y, incY)  \
	                 F77__sgemv(CBLAS_TRANSPOSE, M, N, alpha, A, lda, X, incX, beta, Y, incY)
static INLINE void  F77__sgemv(int trans, int m, int n, F32 alpha, F32PTR a, int lda,F32PTR x, int incx, F32  beta, F32PTR y, int incy)
{
	char TRANSPOSE = (trans == CblasTrans) ? 'T' : 'N';
	F77__CALL(RSGEMV)(&TRANSPOSE, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}


#define r_cblas_sgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)  \
	             	F77__sgemm(CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
static  INLINE void  F77__sgemm(int transa, int transb, int m, int n, int k, F32 alpha, F32PTR a, int lda,F32PTR b, int ldb, F32 beta, F32PTR c, int ldc)
{
	char TRANSA = (transa == CblasTrans) ? 'T' : 'N';
	char TRANSB = (transb == CblasTrans) ? 'T' : 'N';
	F77__CALL(RSGEMM)(&TRANSA, &TRANSB, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}
*/
 


/***************************************************************************************/
/*                           C interfaces to Fortran vector math fucntions             */
/***************************************************************************************/

 
#define r_cblas_sgemv(CBLAS_LAYOUT, CBLAS_TRANSPOSE,  M, N, alpha, A, lda, X, incX, beta, Y, incY)  \
	       f32__sgemv(CBLAS_TRANSPOSE, M, N, alpha, A, lda, X, incX, beta, Y, incY)
static INLINE void  f32__sgemv(int trans, int m, int n, F32 alpha, F32PTR a, int lda, F32PTR x, int incx, F32  beta, F32PTR y, int incy)
{
	if      (trans == CblasTrans)		f32_gemm_XtY2x2(n, 1, m, a, lda, x, m, y, n);
	else if (trans == CblasNoTrans)		f32_gemv_Xb(m, n, a, lda, x, y);
}

#define r_cblas_sgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)  \
		   f32__sgemm(CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
static  INLINE void  f32__sgemm(int transa, int transb, int m, int n, int k, F32 alpha, F32PTR a, int lda, F32PTR b, int ldb, F32 beta, F32PTR c, int ldc)
{
	if (transa == CblasTrans) {
		if      (transb == CblasNoTrans)	f32_gemm_XtY2x2(m, n, k, a, lda, b, ldb, c, ldc);
		else if (transb == CblasTrans)		f32_gemm_XtYt2x2(m, n, k, a, lda, b, ldb, c, ldc);
	}
	else if (transa == CblasNoTrans) {
		if      (transb == CblasNoTrans)	f32_gemm_XY2x2(m, n, k, a, lda, b, ldb, c, ldc);
		else if (transb == CblasTrans)		f32_gemm_XYt2x1(m, n, k, a, lda, b, ldb, c, ldc);
	}

}


//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, F32PTR  a, lapack_int lda);
#define  r_LAPACKE_spotrf(matrix_layout, uplo, n,  a, lda)        f32_chol(a, lda, n);
static INLINE void  f32_chol(F32PTR  a, int lda, int n )        {
	inplace_chol(a, lda, n);
}


//lapack_int LAPACKE_spotrs(int matrix_layout, char uplo, lapack_int n,lapack_int nrhs, const F32PTR  a, lapack_int lda,F32PTR  b, lapack_int ldb);
#define r_LAPACKE_spotrs(matrix_layout, uplo, n,nrhs, a, lda,b, ldb)   \
             f32__spotrs(uplo, n, nrhs, a, lda, b, ldb)
static INLINE void  f32__spotrs(char uplo, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb) {
	solve_U_as_LU_rectmat_multicols(a, b, b, lda, n, nrhs); 
}


 
//LAPACKE_strtrs( int matrix_layout, char uplo, char trans, char diag,lapack_int n, lapack_int nrhs, const F32PTR  a,lapack_int lda, F32PTR  b, lapack_int ldb );
#define r_LAPACKE_strtrs(matrix_layout, uplo, trans, diag,n, nrhs,  a, lda,  b, ldb )  \
      	f32__strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb)
static INLINE void  f32__strtrs(char uplo, char trans, char diag, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb)
{
 
	if     (uplo == 'U' && trans == 'N') {
		for (int i = 0; i < nrhs; ++i) {
			solve_U_as_U(a, b, lda, n);
			b = b + ldb;
		}	
	}
	else if (uplo == 'U' && trans == 'T') {
		for (int i = 0; i < nrhs; ++i) {
			solve_U_as_L(a, b, lda, n);
			b = b + ldb;
		}
	}
	else if (uplo == 'L' && trans == 'N') {
		for (int i = 0; i < nrhs; ++i) {
			solve_L_as_L(a, b, lda, n);
			b = b + ldb;
		}
	}
	else if (uplo == 'L' && trans == 'T') {
		for (int i = 0; i < nrhs; ++i) {
			solve_L_as_U(a, b, lda, n);
			b = b + ldb;
		}
	}
}

 
 

/////////////////////////////////////////////////////////////////////////////////////////////////////
#define r_ippsMaxIndx_32f(X, N, val, idx)				*(I32PTR)(idx)=f32_maxidx( X,  N, val)
#define r_ippsMinIndx_32f(X, N, val, idx)				*(I32PTR)(idx)=f32_minidx( X,  N, val)
#define r_ippsSumLn_32f(pSrc, len, pSum)				*(F32PTR)(pSum)=f32_sumlog(pSrc, len) //#define r_ippsSumLn_32f(pSrc, len, pSum) ippsSumLnR_32f(pSrc, len, pSum) 

#define r_ippsMeanStdDev_32f(pSrc, N, pMean, pSd, hint) f32_avgstd( pSrc,  N, pMean, pSd) ////ippsMeanStdDev_32f(pSrc, N, pMean, pStdDev, hint)
#define r_vsPowx(n, src,power, res)						f32_pow_vec_inplace(src,power, n) //vmsPowx(const MKL_INT n, const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
#define r_vmsLn(n, src, res, mode)						f32_log_vec_inplace(src, n) ////ippsLn_32f_I(pSrcDst, N)
#define r_ippsLn_32f_I(pSrcDst, N)						f32_log_vec_inplace(pSrcDst, N)
#define r_vmsCos(n, a, r, mode)							f32_cos_vec_inplace(a,n)
#define r_vmsSin(n, a, r, mode)							f32_sin_vec_inplace(a, n)
#define r_ippsSqrt_32f_I( pSrc, len)					f32_sqrt_vec_inplace(pSrc, len)

#define r_ippsSub_32f_I(pSrc, pDst,len)					f32_sub_vec_inplace( pSrc,  pDst, len)// F77__ssub_i(pSrc,pDst, len)
#define r_ippsSub_32f(pSrc1, pSrc2, pDst,len)			f32_sub_vec(pSrc1, pSrc2, pDst,len)
#define r_ippsAdd_32f_I(pSrc, pDst,len)					f32_add_vec_inplace(pSrc,pDst, len)
#define r_ippsSubC_32f_I(val, pSrcDst, N)				f32_add_val_inplace(-(val),pSrcDst,N)//F77__ssubc(val, pSrcDst, N)
#define r_ippsSubCRev_32f_I(val, pSrcDst, N)			f32_subrev_val_inplace(val, pSrcDst, N) 
#define r_ippsMul_32f( pSrc1, pSrc2, pDst, N)			f32_mul_vec(pSrc1,pSrc2,pDst, N)
#define r_ippsMul_32f_I(pSrc, pSrcDst, N)				f32_mul_vec_inplace(pSrc, pSrcDst, N)

//#define ONE_STEP_DIFF(pSrcX,pSrcY,pDst,N)				F77__sdiff(pSrcX,N)
#define ONE_STEP_DIFF(pSrcX,pSrcY,pDst,N)				f32_diff_back(pSrcX, pDst,N)

#define r_ippsSum_32f(Src, N, Sum, PrecisionLevel)		*(F32PTR)(Sum)=f32_sum(Src, N)
#define r_ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor)	*(I32PTR)(Sum)=i32_sum(Src, N)

#define r_cblas_scopy(N, src, incX, dst, incY)			memcpy(dst, src, (N)*sizeof(F32))

//https: //stackoverflow.com/questions/32020854/conditional-inside-c-macro
//https: //github.com/pfultz2/Cloak/wiki/C-Preprocessor-tricks,-tips,-and-idioms
/*
#define r_ippsSet_32f(value, dst, N)					r_ippsSet_32f##value(dst, N)
#define r_ippsSet_32f0(dst, N)							memset(dst, 0, sizeof(F32)*(N))
#define r_ippsSet_32f1(dst, N)							fill_float32(dst, N)
*/
#define r_ippsSet_32f(value, dst, N)					f32_fill_val(value, dst, N)

//Zero-fill an int32 array
#define r_ippsSet_32s(value, dst, N)					r_ippsSet_32s##value(dst,N)
#define r_ippsSet_32s0(dst, N)							memset(dst, 0, sizeof(I32)*(N))

//ippsSet_8u(val, PDst, N)
#define r_ippsSet_8u(val, pDst, N)						memset(pDst, val, N) 

#define r_cblas_sscal(N, alpha, X, incX)				f32_mul_val_inplace( alpha, X,   N)
#define r_ippsMulC_32f_I(alpha, SrcDst, N)				f32_mul_val_inplace( alpha, SrcDst, N) ////ippsMulC_32f_I(Ipp32f  val, Ipp32f*  pSrcDst, int len))
#define r_ippsAddC_32s_ISfs(val, X, N, scaleFactor)     i32_add_val_inplace(val,X,N)
#define DOT(N,X,Y)										f32_dot(X,Y,N)

//This function should be included only once;  otherwise, the redefination error will ocucr. Here we
// use the header guard to ensure this header is included only once in each translation unit.
//static void fill_float32(rF32 value, rF32PTR dst, int N) { for (int i = N; i > 0; i--) *dst++ = value; }

