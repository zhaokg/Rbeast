//#ifndef  BEAST_BLAS_LAPACK_HEADER2
//#define  BEAST_BLAS_LAPACK_HEADER2
#pragma once

#include "abc_000_macro.h"
//#include <inttypes.h> //#include <stdint.h>
#include "abc_datatype.h"

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

/****************************************************************************************/
#define r_cblas_ssymv(layout, uplo, n, alpha, A, lda, x,incx,beta,y, incy) \
		F77__ssymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

#define r_cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX)  \
		F77__strmv(Uplo, TransA, Diag, N, A, lda, X, incX)

#define r_cblas_sgemv(CBLAS_LAYOUT, CBLAS_TRANSPOSE,  M, N, alpha, A, lda, X, incX, beta, Y, incY)  \
		F77__sgemv(CBLAS_TRANSPOSE, M, N, alpha, A, lda, X, incX, beta, Y, incY)

#define r_cblas_sgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B,M, N,K, alpha, A,\
		lda, B, ldb, beta, C, ldc)  \
		F77__sgemm('T', 'N', M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)


#define r_LAPACKE_strtrs(matrix_layout, uplo, trans, diag,n, nrhs,  a, lda,  b, ldb ) \
	F77__strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb)

#define r_LAPACKE_spotrs(matrix_layout, uplo, n,nrhs, a, lda,b, ldb)  \
	F77__spotrs(uplo, n, nrhs, a, lda, b, ldb)

#define  r_LAPACKE_spotrf(matrix_layout, uplo, n,  a, lda)  \
	F77__spotrf(uplo, n, a, lda)

//{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
//https: //en.wikipedia.org/wiki/Inline_function
static INLINE void  F77__ssymv(int uplo, int n, F32 alpha, F32PTR A, int lda, F32PTR x, int incx, F32 beta, F32PTR y, int incy);

static INLINE void  F77__strmv(int uplo, int trans, int diag, int n, F32PTR a, int lda, F32PTR x, int incx);
//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);

static INLINE void  F77__sgemv(int trans, int m, int n, F32 alpha, F32PTR a, int lda,
	F32PTR x, int incx, F32  beta, F32PTR y, int incy);

static INLINE void   F77__sgemm(char transa, char transb, int m, int n, int k, F32 alpha, F32PTR a, int lda,
	F32PTR b, int ldb, F32 beta, F32PTR c, int ldc);

static INLINE int  F77__strtrs(char uplo, char trans, char diag, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb);
//LAPACKE_strtrs( int matrix_layout, char uplo, char trans, char diag,lapack_int n, lapack_int nrhs, const F32PTR  a,lapack_int lda, F32PTR  b, lapack_int ldb );

static INLINE int  F77__spotrs(char uplo, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb);
//lapack_int LAPACKE_spotrs(int matrix_layout, char uplo, lapack_int n,lapack_int nrhs, const F32PTR  a, lapack_int lda,F32PTR  b, lapack_int ldb);

static INLINE int  F77__spotrf(char uplo, int n, F32PTR  a, int lda);
//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, F32PTR  a, lapack_int lda);

//}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
 



/*
 * START Declarations of Fortran Blas or lapack functions
 */
void  F77__CALL(RSGEMV)(char *trans, int *m, int *n, F32PTR alpha, F32PTR a, const int *lda, F32PTR x, int *incx, F32PTR beta, F32PTR y, int *incy, size_t xxx_Char_Len3);
void  F77__CALL(RSGEMM)(char *, char*, int*, int*, int*, F32PTR , F32PTR , int*, F32PTR , int*, F32PTR , F32PTR , int*, size_t xxx_Char_Len1, size_t xxx_Char_Len2);
void  F77__CALL(RSTRMV)(char *uplo, char *trans, char *diag, int *n, F32PTR a, int *lda, F32PTR x, int *incx, size_t xxx_Char_Len1, size_t xxx_Char_Len2, size_t xxx_Char_Len3);
void  F77__CALL(RSSYMV)(char *UPLO, int* N, F32PTR  ALPHA, F32PTR  A, int * lda, F32PTR X, int * INCX, F32PTR  beta, F32PTR Y, int * INCY, size_t xxx_Char_Len1);

void  F77__CALL(RSTRTRS)(char* uplo, char* trans, char* diag, int* n, int* nrhs, F32PTR  a, int* lda, F32PTR  b, int* ldb, int* info, size_t xxx_Char_Len1, size_t xxx_Char_Len2, size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRF)(char* uplo, int* n, F32PTR  a, int* lda, int* info, size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRS)(char* uplo, int* n, int* nrhs, F32PTR  a, int* lda, F32PTR  b, int* ldb, int* info, size_t xxx_Char_Len1);
/* Declarations of Vector-based Fortran math funcitions Vector-based Fortran math funcitions */




/*<<<<<START:  C interfaces to the above Fortran Math functions */
//void cblas_ssymv(const CBLAS_LAYOUT Layout, const CBLAS_UPLO uplo, const MKL_INT n, const F32 alpha, const F32PTR a, const MKL_INT lda, const F32PTR x, const MKL_INT incx, const F32 beta, F32PTR y, const MKL_INT incy);
static INLINE void  F77__ssymv(int uplo, int n, F32 alpha, F32PTR A, int lda, F32PTR x, int incx, F32 beta, F32PTR y, int incy)
{
	char UP = (uplo == CblasUpper) ? 'U' : 'L';
	//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
	F77__CALL(RSSYMV)(&UP, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
#define r_cblas_ssymv(layout, uplo, n, alpha, A, lda, x, incx, beta, y, incy) \
	F77__ssymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy)

static INLINE void  F77__strmv(int uplo, int trans, int diag, int n, F32PTR a, int lda, F32PTR x, int incx)
{
	char UP = (uplo == CblasUpper) ? 'U' : 'L';
	char TRANSPOSE = (trans == CblasTrans) ? 'T' : 'N';
	char DIAGONAL = (diag == CblasUnit) ? 'U' : 'N';
	//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
	F77__CALL(RSTRMV)(&UP, &TRANSPOSE, &DIAGONAL, &n, a, &lda, x, &incx, 1, 1, 1);
}
//cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX);
#define r_cblas_strmv(Layout, Uplo,TransA, Diag, N, A,  lda, X,incX)  \
	F77__strmv(Uplo, TransA, Diag, N, A, lda, X, incX)

static INLINE int  F77__strtrs(char uplo, char trans, char diag, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb)
{
	int info;
	F77__CALL(RSTRTRS)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
	return info;
}
//LAPACKE_strtrs( int matrix_layout, char uplo, char trans, char diag,lapack_int n, lapack_int nrhs, const F32PTR  a,lapack_int lda, F32PTR  b, lapack_int ldb );
#define r_LAPACKE_strtrs(matrix_layout, uplo, trans, diag,n, nrhs,  a, lda,  b, ldb ) \
	F77__strtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb)

static INLINE int  F77__spotrs(char uplo, int n, int nrhs, F32PTR  a, int lda, F32PTR  b, int ldb)
{
	int info;
	F77__CALL(RSPOTRS)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
	return info;
}
//lapack_int LAPACKE_spotrs(int matrix_layout, char uplo, lapack_int n,lapack_int nrhs, const F32PTR  a, lapack_int lda,F32PTR  b, lapack_int ldb);
#define r_LAPACKE_spotrs(matrix_layout, uplo, n,nrhs, a, lda,b, ldb)  \
	F77__spotrs(uplo, n, nrhs, a, lda, b, ldb)


static INLINE int  F77__spotrf(char uplo, int n, F32PTR  a, int lda) {
	int info;
	F77__CALL(RSPOTRF)(&uplo, &n, a, &lda, &info, 1);
	return info;
}
//lapack_int LAPACKE_spotrf(int matrix_layout, char uplo, lapack_int n, F32PTR  a, lapack_int lda);
//#define  r_LAPACKE_spotrf(matrix_layout, uplo, n,  a, lda)  F77__spotrf(uplo,n,a,lda)

static INLINE void  F77__sgemv(int trans, int m, int n, F32 alpha, F32PTR a, int lda,
	F32PTR x, int incx, F32  beta, F32PTR y, int incy)
{
	char TRANSPOSE = (trans == CblasTrans) ? 'T' : 'N';
	F77__CALL(RSGEMV)(&TRANSPOSE, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}
#define r_cblas_sgemv(CBLAS_LAYOUT, CBLAS_TRANSPOSE,  M, N, alpha, A, lda, X, incX, beta, Y, incY)  \
	F77__sgemv(CBLAS_TRANSPOSE, M, N, alpha, A, lda, X, incX, beta, Y, incY)

static  INLINE void  F77__sgemm(int transa, int transb, int m, int n, int k, F32 alpha, F32PTR a, int lda,
	F32PTR b, int ldb, F32 beta, F32PTR c, int ldc)
{
	char TRANSA = (transa == CblasTrans) ? 'T' : 'N';
	char TRANSB = (transb == CblasTrans) ? 'T' : 'N';
	F77__CALL(RSGEMM)(&TRANSA, &TRANSB, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}
#define r_cblas_sgemm(CBLAS_LAYOUT, CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A,\
	                  lda, B, ldb, beta, C, ldc)  \
	F77__sgemm(CBLAS_TRANSPOSE_A, CBLAS_TRANSPOSE_B, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

/* CODE_EOF>>>>>>>:  C interfaces to the Fortran Math functions            */


/***************************************************************************************/
/*                           C interfaces to Fortran vector math fucntions             */
/***************************************************************************************/

#if defined(COMPILER_MSVC)
	//#define F77__CALL(x)  x
#elif defined(COMPILER_CLANG)|| defined(COMPILER_GCC) ||defined(COMPILER_SOLARIS) 

	//#define PRIMITIVE_CAT(a, ...) a##__VA_ARGS__
	//#define F77__CALL(x)		PRIMITIVE_CAT(x, _)

	#define RSSUMLOG 		rssumlog
	#define RSPOWVEC 		rspowvec
	#define RSLOGVEC 		rslogvec
	#define RSCOSVEC 		rscosvec 
	#define RSSINVEC 		rssinvec 
	#define RSSQRTVEC 		rssqrtvec 
	#define RSSUB_I 		rssub_i
	#define RSSUB 			rssub
	#define RSDIFF 			rsdiff 
	#define RSADD 			rsadd         
	#define RSSUBC 			rssubc 
	#define RSCSUB 			rscsub 
	#define RSSQRVEC 		rssqrvec 
	#define RSSUM 			rssum 
	#define RISUM			risum
	#define RSMINIDX 		rsminidx 
	#define RSMAXIDX 		rsmaxidx 
	#define RSMNSD 			rsmnsd 
	#define RSDOT 			rsdot 
	#define RIADDC          riaddc
	#define RSSCAL          rsscal

#endif

/****************************************************************************************/
#define r_ippsMaxIndx_32f(X, N, val, idx)				F77__smaxidx(X, N, val, idx)
#define r_ippsMinIndx_32f(X, N, val, idx)				F77__sminidx(X, N, val, idx)
#define r_ippsSumLn_32f(pSrc, len, pSum)				F77__ssumlog(pSrc, len, pSum) //#define r_ippsSumLn_32f(pSrc, len, pSum) ippsSumLnR_32f(pSrc, len, pSum) 
#define r_ippsMeanStdDev_32f(pSrc, N, pMean, pStdDev, hint) F77__smnsd(pSrc, N, pMean, pStdDev) ////ippsMeanStdDev_32f(pSrc, N, pMean, pStdDev, hint)
#define r_vsPowx(n, a,b, r)								F77__spowvec(a,b, n) //vmsPowx(const MKL_INT n, const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
#define r_vmsLn(n, a, r, mode)							F77__slogvec(a, n) ////ippsLn_32f_I(pSrcDst, N)
#define r_ippsLn_32f_I(pSrcDst, N)						F77__slogvec(pSrcDst, N)
#define r_vmsCos(n, a, r, mode)							F77__scosvec(a, n)
#define r_vmsSin(n, a, r, mode)							F77__ssinvec(a, n)
#define r_ippsSqrt_32f_I( pSrc, len)					F77__ssqrtvec(pSrc ,len)
#define r_ippsSub_32f_I(pSrc, pDst,len)					F77__ssub_i(pSrc,pDst, len)
#define r_ippsSub_32f(pSrc1, pSrc2, pDst,len)			F77__ssub(pSrc1, pSrc2, pDst,len)
#define ONE_STEP_DIFF(pSrcX,pSrcY,pDst,N)				F77__sdiff(pSrcX,N)
#define r_ippsAdd_32f_I(pSrc, pDst,len)					F77__sadd(pSrc,pDst, len)
#define r_ippsSubC_32f_I(val, pSrcDst, N)				F77__ssubc(val, pSrcDst, N)
#define r_ippsSubCRev_32f_I(val, pSrcDst, N)			F77__scsub(val, pSrcDst, N) 
#define r_ippsMul_32f( pSrc1, pSrc2, pDst, N)			F77__ssqrvec(pSrc1,pDst, N)
#define r_ippsMul_32f_I(pSrc, pSrcDst, N)				F77__ssqrvec(pSrc, pSrcDst, N)

#define r_ippsSum_32f(Src, N, Sum, PrecisionLevel)		F77__ssum(Src, N, Sum)
#define r_ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor)	F77__isum(Src, N, Sum)

#define r_cblas_scopy(N, src, incX, dst, incY)			memcpy(dst, src, (N)*sizeof(F32))

//https: //stackoverflow.com/questions/32020854/conditional-inside-c-macro
//https: //github.com/pfultz2/Cloak/wiki/C-Preprocessor-tricks,-tips,-and-idioms
/*
#define r_ippsSet_32f(value, dst, N)					r_ippsSet_32f##value(dst, N)
#define r_ippsSet_32f0(dst, N)							memset(dst, 0, sizeof(F32)*(N))
#define r_ippsSet_32f1(dst, N)							fill_float32(dst, N)
*/
#define r_ippsSet_32f(value, dst, N)					fill_float32(value, dst, N)

//Zero-fill an int32 array
#define r_ippsSet_32s(value, dst, N)					r_ippsSet_32s##value(dst,N)
#define r_ippsSet_32s0(dst, N)							memset(dst, 0, sizeof(I32)*(N))

//ippsSet_8u(val, PDst, N)
#define r_ippsSet_8u(val, pDst, N)						memset(pDst, val, N) 
#define DOT(N,X,Y)										F77__sdot(N,X,1,Y,1)
#define r_cblas_sscal(N, alpha, X, incX)				F77__sscal(N, alpha, X, incX);
#define r_ippsMulC_32f_I(val, SrcDst, N)				F77__sscal(N, val, SrcDst, 1) ////ippsMulC_32f_I(Ipp32f  val, Ipp32f*  pSrcDst, int len))
#define r_ippsAddC_32s_ISfs(val, X, N, scaleFactor)     F77__iaddc(val,X,N)

/****************************************************************************************/

static INLINE void  F77__smaxidx(F32 *X, int N, F32PTR val, int * idx);
static INLINE void  F77__sminidx(F32 *X, int N, F32PTR val, int * idx);
static INLINE void  F77__ssumlog(F32 *X, int N, F32PTR pSum);	//ippsSumLn_32f(const Ipp32f* pSrc, int len, Ipp32f* pSum))
static INLINE void  F77__smnsd(F32 *X, int N, F32PTR MEAN, F32PTR  STD);
static INLINE void  F77__spowvec(F32 *X, F32 pow, int N);  //vmsPowx(const MKL_INT n,  const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
static INLINE void F77__slogvec(F32 *X, int N);	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
static INLINE void F77__scosvec(F32 *X, int N);	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
static INLINE void F77__ssinvec(F32 *X, int N);   //vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
static INLINE void F77__ssqrtvec(F32 *X, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
static INLINE void F77__ssub_i(F32 *X, F32PTR Y, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
static INLINE void F77__ssub(F32 *X, F32PTR Y, F32PTR Z, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
static INLINE void F77__sdiff(F32 *X, int N);	//Z=Y-Z;	
static INLINE void F77__sadd(F32 *X, F32PTR Y, int N);//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
static INLINE void F77__ssubc(F32 C, F32PTR Y, int N); //ippsSubC_32f_I(val, pSrcDst, N)
static INLINE void F77__scsub(F32 C, F32PTR Y, int N);//r_ippsSubCRev_32f_I(val, pSrcDst, N)
static INLINE void F77__ssqrvec(F32 *X, F32PTR Y, int N);//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
static INLINE void F77__ssum(F32 *X, int N, F32PTR ANS);//ippsSum_32f(Src, N, Sum, PrecisionLevel);
static INLINE void F77__isum(int *X, int N, int *ANS);//ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor);
static INLINE F32  F77__sdot(const int N, const F32PTR x, const int incx, const F32PTR y, const int incy);
static INLINE void  F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX);
static INLINE void  F77__iaddc(const int C, int *X, const int N);

/**
**********************************************************
**/
static INLINE void fill_float32(rF32 value, rF32PTR dst, int N) { for (rI32 i = N; i > 0; i--) *dst++ = value; }
//This function should be included only once;  otherwise, the redefination error will ocucr. Here
//we use the header guard to ensure this header is included only once in each translation unit.


//#endif


/*<<<<<START Declarations of Vector-based Fortran math funcitions--*/
void F77__CALL(RSSUMLOG)(F32 *x, int *N, F32PTR ans);
void F77__CALL(RSPOWVEC)(F32 *x, F32PTR  pow, int *N);
void F77__CALL(RSLOGVEC)(F32 *x, int *N);
void F77__CALL(RSCOSVEC)(F32 *x, int *N);
void F77__CALL(RSSINVEC)(F32 *x, int *N);
void F77__CALL(RSSQRTVEC)(F32 *x, int *N);
void F77__CALL(RSSUB_I)(F32 *x, F32PTR y, int *N);
void F77__CALL(RSSUB)(F32 *x, F32PTR y, F32PTR z, int *N);
void F77__CALL(RSDIFF)(F32 *x, int *N);
void F77__CALL(RSADD)(F32 *x, F32PTR y, int *N);
void F77__CALL(RSSUBC)(F32 *C, F32PTR y, int *N);
void F77__CALL(RSCSUB)(F32 *C, F32PTR y, int *N);
void F77__CALL(RSSQRVEC)(F32 *x, F32PTR y, int *N);
void F77__CALL(RSSUM)(F32 *x, int* N, F32PTR  ANS);
void F77__CALL(RISUM)(int *x, int* N, int* ANS);
void F77__CALL(RSMINIDX)(F32 *, int *N, F32PTR val, int * idx);
void F77__CALL(RSMAXIDX)(F32 *, int *N, F32PTR val, int * idx);
void F77__CALL(RSMNSD)(F32 * SX, int * N, F32PTR MEAN, F32PTR  STD);
F32 F77__CALL(RSDOT)(const int *n, const F32 * x, const int *incx, const F32PTR y, const int *incy);
void F77__CALL(RSSCAL)(const int *n, const F32 * alpha, F32PTR x, const int *incx);
void F77__CALL(RIADDC)(const int *C, int *x, const int *N);
/*Declarations of Vector-based Fortran math funcitions Vector-based Fortran math funcitions CODE_EOF>>>>>>
*/

/************************************************************************/
/****** START:  C interfaces to the above Fortran Math functions *******/
/************************************************************************/
static INLINE void   F77__smaxidx(F32 *X, int N, F32PTR val, int * idx) {
	F77__CALL(RSMAXIDX)(X, &N, val, idx);
}

static INLINE void F77__sminidx(F32 *X, int N, F32PTR val, int * idx) {
	F77__CALL(RSMINIDX)(X, &N, val, idx);
}

static INLINE void  F77__ssumlog(F32 *X, int N, F32PTR pSum){	//ippsSumLn_32f(const Ipp32f* pSrc, int len, Ipp32f* pSum))
	F77__CALL(RSSUMLOG)(X, &N, pSum);
}

static INLINE void   F77__smnsd(F32 *X, int N, F32PTR MEAN, F32PTR  STD){
	F77__CALL(RSMNSD)(X, &N, MEAN, STD);
}

static INLINE void  F77__spowvec(F32 *X, F32 pow, int N){//vmsPowx(const MKL_INT n,  const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
	F77__CALL(RSPOWVEC)(X, &pow, &N);

}

static INLINE void  F77__slogvec(F32 *X, int N){	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__CALL(RSLOGVEC)(X, &N);
}//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))

static INLINE void  F77__scosvec(F32 *X, int N){	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__CALL(RSCOSVEC)(X, &N);
}//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))

static INLINE void  F77__ssinvec(F32 *X, int N){//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__CALL(RSSINVEC)(X, &N);
}//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))


static INLINE void F77__ssqrtvec(F32 *X, int N){	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__CALL(RSSQRTVEC)(X, &N);
} //ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))


static INLINE void  F77__ssub_i(F32 *X, F32PTR Y, int N){	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__CALL(RSSUB_I)(X, Y, &N);
}//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))

static INLINE void  F77__ssub(F32 *X, F32PTR Y, F32PTR Z, int N){	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__CALL(RSSUB)(X, Y, Z, &N);
}//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))


static INLINE void  F77__sdiff(F32 *X, int N){	//Z=Y-Z;	
	F77__CALL(RSDIFF)(X, &N);
}

static INLINE void  F77__sadd(F32 *X, F32PTR Y, int N){
	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__CALL(RSADD)(X, Y, &N);
}//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))


static INLINE  void  F77__ssubc(F32 C, F32PTR Y, int N){
	//ippsSubC_32f_I(val, pSrcDst, N)
	F77__CALL(RSSUBC)(&C, Y, &N);
} //ippsSubC_32f_I(val, pSrcDst, N)


static INLINE void  F77__scsub(F32 C, F32PTR Y, int N){
	//ippsSubC_32f_I(val, pSrcDst, N)
	F77__CALL(RSCSUB)(&C, Y, &N);
}//r_ippsSubCRev_32f_I(val, pSrcDst, N)


static INLINE void  F77__ssqrvec(F32 *X, F32PTR Y, int N) {
	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__CALL(RSSQRVEC)(X, Y, &N);
} //ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))

static INLINE void  F77__ssum(F32 *X, int N, F32PTR ANS){
	//ippsSum_32f(Src, N, Sum, PrecisionLevel);
	F77__CALL(RSSUM)(X, &N, ANS);
}//ippsSum_32f(Src, N, Sum, PrecisionLevel);


static INLINE void  F77__isum(int *X, int N, int *ANS){
	//ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor);
	F77__CALL(RISUM)(X, &N, ANS);
}//ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor);


static INLINE F32  F77__sdot(const int N, const F32PTR x, const int incx, const F32PTR y, const int incy)
{
	return F77__CALL(RSDOT)(&N, x, &incx, y, &incy);
}


static INLINE void  F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX){
	F77__CALL(RSSCAL)(&N, &alpha, X, &incX);
}

static INLINE void  F77__iaddc(const int C, int *X, const int N){
	F77__CALL(RIADDC)(&C, X, &N);
}

/************************************************************************/
/****** CODE_EOF : C interfaces to the Fortran Math functions          *******/
/************************************************************************/




