
#pragma once
#include "abc_000_macro.h"
#include "abc_datatype.h"
#if defined(MSVC_COMPILER)
	#define F77__CALL(x)  x
#elif defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER)  
	#define PRIMITIVE_CAT(a,...) a##__VA_ARGS__
	#define F77__CALL(x)		PRIMITIVE_CAT(x,_)
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
#define r_cblas_ssymv(layout,uplo,n,alpha,A,lda,x,incx,beta,y,incy) \
		F77__ssymv(uplo,n,alpha,A,lda,x,incx,beta,y,incy)
#define r_cblas_strmv(Layout,Uplo,TransA,Diag,N,A,lda,X,incX)  \
		F77__strmv(Uplo,TransA,Diag,N,A,lda,X,incX)
#define r_cblas_sgemv(CBLAS_LAYOUT,CBLAS_TRANSPOSE,M,N,alpha,A,lda,X,incX,beta,Y,incY)  \
		F77__sgemv(CBLAS_TRANSPOSE,M,N,alpha,A,lda,X,incX,beta,Y,incY)
#define r_cblas_sgemm(CBLAS_LAYOUT,CBLAS_TRANSPOSE_A,CBLAS_TRANSPOSE_B,M,N,K,alpha,A,\
		lda,B,ldb,beta,C,ldc)  \
		F77__sgemm('T','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc)
#define r_LAPACKE_strtrs(matrix_layout,uplo,trans,diag,n,nrhs,a,lda,b,ldb ) \
	F77__strtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb)
#define r_LAPACKE_spotrs(matrix_layout,uplo,n,nrhs,a,lda,b,ldb)  \
	F77__spotrs(uplo,n,nrhs,a,lda,b,ldb)
#define  r_LAPACKE_spotrf(matrix_layout,uplo,n,a,lda)  \
	F77__spotrf(uplo,n,a,lda)
static INLINE void  F77__ssymv(int uplo,int n,float alpha,float *A,int lda,float *x,int incx,float beta,float *y,int incy);
static INLINE void  F77__strmv(int uplo,int trans,int diag,int n,float *a,int lda,float *x,int incx);
static INLINE void  F77__sgemv(int trans,int m,int n,float alpha,float *a,int lda,
	float *x,int incx,float  beta,float *y,int incy);
static INLINE void   F77__sgemm(char transa,char transb,int m,int n,int k,float alpha,float *a,int lda,
	float *b,int ldb,float beta,float *c,int ldc);
static INLINE int  F77__strtrs(char uplo,char trans,char diag,int n,int nrhs,float * a,int lda,float* b,int ldb);
static INLINE int  F77__spotrs(char uplo,int n,int nrhs,float* a,int lda,float* b,int ldb);
static INLINE int  F77__spotrf(char uplo,int n,float* a,int lda);
void  F77__CALL(RSGEMV)(char *trans,int *m,int *n,float *alpha,float *a,const int *lda,float *x,int *incx,float *beta,float *y,int *incy,size_t xxx_Char_Len3);
void  F77__CALL(RSGEMM)(char *,char*,int*,int*,int*,float*,float*,int*,float*,int*,float*,float*,int*,size_t xxx_Char_Len1,size_t xxx_Char_Len2);
void  F77__CALL(RSTRMV)(char *uplo,char *trans,char *diag,int *n,float *a,int *lda,float *x,int *incx,size_t xxx_Char_Len1,size_t xxx_Char_Len2,size_t xxx_Char_Len3);
void  F77__CALL(RSSYMV)(char *UPLO,int* N,float * ALPHA,float * A,int * lda,float *X,int * INCX,float * beta,float *Y,int * INCY,size_t xxx_Char_Len1);
void  F77__CALL(RSTRTRS)(char* uplo,char* trans,char* diag,int* n,int* nrhs,float * a,int* lda,float* b,int* ldb,int* info,size_t xxx_Char_Len1,size_t xxx_Char_Len2,size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRF)(char* uplo,int* n,float* a,int* lda,int* info,size_t xxx_Char_Len3);
void  F77__CALL(RSPOTRS)(char* uplo,int* n,int* nrhs,float* a,int* lda,float* b,int* ldb,int* info,size_t xxx_Char_Len1);
static INLINE void  F77__ssymv(int uplo,int n,float alpha,float *A,int lda,float *x,int incx,float beta,float *y,int incy)
{
	char UP=(uplo==CblasUpper) ? 'U' : 'L';
	F77__CALL(RSSYMV)(&UP,&n,&alpha,A,&lda,x,&incx,&beta,y,&incy,1);
}
#define r_cblas_ssymv(layout,uplo,n,alpha,A,lda,x,incx,beta,y,incy) \
	F77__ssymv(uplo,n,alpha,A,lda,x,incx,beta,y,incy)
static INLINE void  F77__strmv(int uplo,int trans,int diag,int n,float *a,int lda,float *x,int incx)
{
	char UP=(uplo==CblasUpper) ? 'U' : 'L';
	char TRANSPOSE=(trans==CblasTrans) ? 'T' : 'N';
	char DIAGONAL=(diag==CblasUnit) ? 'U' : 'N';
	F77__CALL(RSTRMV)(&UP,&TRANSPOSE,&DIAGONAL,&n,a,&lda,x,&incx,1,1,1);
}
#define r_cblas_strmv(Layout,Uplo,TransA,Diag,N,A,lda,X,incX)  \
	F77__strmv(Uplo,TransA,Diag,N,A,lda,X,incX)
static INLINE int  F77__strtrs(char uplo,char trans,char diag,int n,int nrhs,float * a,int lda,float* b,int ldb)
{
	int info;
	F77__CALL(RSTRTRS)(&uplo,&trans,&diag,&n,&nrhs,a,&lda,b,&ldb,&info,1,1,1);
	return info;
}
#define r_LAPACKE_strtrs(matrix_layout,uplo,trans,diag,n,nrhs,a,lda,b,ldb ) \
	F77__strtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb)
static INLINE int  F77__spotrs(char uplo,int n,int nrhs,float* a,int lda,float* b,int ldb)
{
	int info;
	F77__CALL(RSPOTRS)(&uplo,&n,&nrhs,a,&lda,b,&ldb,&info,1);
	return info;
}
#define r_LAPACKE_spotrs(matrix_layout,uplo,n,nrhs,a,lda,b,ldb)  \
	F77__spotrs(uplo,n,nrhs,a,lda,b,ldb)
static INLINE int  F77__spotrf(char uplo,int n,float* a,int lda) {
	int info;
	F77__CALL(RSPOTRF)(&uplo,&n,a,&lda,&info,1);
	return info;
}
static INLINE void  F77__sgemv(int trans,int m,int n,float alpha,float *a,int lda,
	float *x,int incx,float  beta,float *y,int incy)
{
	char TRANSPOSE=(trans==CblasTrans) ? 'T' : 'N';
	F77__CALL(RSGEMV)(&TRANSPOSE,&m,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy,1);
}
#define r_cblas_sgemv(CBLAS_LAYOUT,CBLAS_TRANSPOSE,M,N,alpha,A,lda,X,incX,beta,Y,incY)  \
	F77__sgemv(CBLAS_TRANSPOSE,M,N,alpha,A,lda,X,incX,beta,Y,incY)
static  INLINE void  F77__sgemm(char transa,char transb,int m,int n,int k,float alpha,float *a,int lda,
	float *b,int ldb,float beta,float *c,int ldc)
{
	F77__CALL(RSGEMM)(&transa,&transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc,1,1);
}
#define r_cblas_sgemm(CBLAS_LAYOUT,CBLAS_TRANSPOSE_A,CBLAS_TRANSPOSE_B,M,N,K,alpha,A,\
	lda,B,ldb,beta,C,ldc)  \
	F77__sgemm('T','N',M,N,K,alpha,A,lda,B,ldb,beta,C,ldc)
#if defined(MSVC_COMPILER)
#elif defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER) 
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
#define r_ippsMaxIndx_32f(X,N,val,idx)				F77__smaxidx(X,N,val,idx)
#define r_ippsMinIndx_32f(X,N,val,idx)				F77__sminidx(X,N,val,idx)
#define r_ippsSumLn_32f(pSrc,len,pSum)				F77__ssumlog(pSrc,len,pSum) 
#define r_ippsMeanStdDev_32f(pSrc,N,pMean,pStdDev,hint) F77__smnsd(pSrc,N,pMean,pStdDev) 
#define r_vsPowx(n,a,b,r)								F77__spowvec(a,b,n) 
#define r_vmsLn(n,a,r,mode)							F77__slogvec(a,n) 
#define r_ippsLn_32f_I(pSrcDst,N)						F77__slogvec(pSrcDst,N)
#define r_vmsCos(n,a,r,mode)							F77__scosvec(a,n)
#define r_vmsSin(n,a,r,mode)							F77__ssinvec(a,n)
#define r_ippsSqrt_32f_I( pSrc,len)					F77__ssqrtvec(pSrc,len)
#define r_ippsSub_32f_I(pSrc,pDst,len)					F77__ssub_i(pSrc,pDst,len)
#define r_ippsSub_32f(pSrc1,pSrc2,pDst,len)			F77__ssub(pSrc1,pSrc2,pDst,len)
#define ONE_STEP_DIFF(pSrcX,pSrcY,pDst,N)				F77__sdiff(pSrcX,N)
#define r_ippsAdd_32f_I(pSrc,pDst,len)					F77__sadd(pSrc,pDst,len)
#define r_ippsSubC_32f_I(val,pSrcDst,N)				F77__ssubc(val,pSrcDst,N)
#define r_ippsSubCRev_32f_I(val,pSrcDst,N)			F77__scsub(val,pSrcDst,N) 
#define r_ippsMul_32f( pSrc1,pSrc2,pDst,N)			F77__ssqrvec(pSrc1,pDst,N)
#define r_ippsMul_32f_I(pSrc,pSrcDst,N)				F77__ssqrvec(pSrc,pSrcDst,N)
#define r_ippsSum_32f(Src,N,Sum,PrecisionLevel)		F77__ssum(Src,N,Sum)
#define r_ippsSum_32s_Sfs(Src,N,Sum,ScaleFfactor)	F77__isum(Src,N,Sum)
#define r_cblas_scopy(N,src,incX,dst,incY)			memcpy(dst,src,(N)*sizeof(float))
#define r_ippsSet_32f(value,dst,N)					fill_float32(value,dst,N)
#define r_ippsSet_32s(value,dst,N)					r_ippsSet_32s##value(dst,N)
#define r_ippsSet_32s0(dst,N)							memset(dst,0,sizeof(int32_t)*(N))
#define r_ippsSet_8u(val,pDst,N)						memset(pDst,val,N) 
#define DOT(N,X,Y)										F77__sdot(N,X,1,Y,1)
#define r_cblas_sscal(N,alpha,X,incX)				F77__sscal(N,alpha,X,incX);
#define r_ippsMulC_32f_I(val,SrcDst,N)				F77__sscal(N,val,SrcDst,1) 
#define r_ippsAddC_32s_ISfs(val,X,N,scaleFactor)     F77__iaddc(val,X,N)
static INLINE void  F77__smaxidx(float *X,int N,float *val,int * idx);
static INLINE void  F77__sminidx(float *X,int N,float *val,int * idx);
static INLINE void  F77__ssumlog(float *X,int N,float *pSum);	
static INLINE void  F77__smnsd(float *X,int N,float *MEAN,float * STD);
static INLINE void  F77__spowvec(float *X,float pow,int N);  
static INLINE void F77__slogvec(float *X,int N);	
static INLINE void F77__scosvec(float *X,int N);	
static INLINE void F77__ssinvec(float *X,int N);   
static INLINE void F77__ssqrtvec(float *X,int N);	
static INLINE void F77__ssub_i(float *X,float *Y,int N);	
static INLINE void F77__ssub(float *X,float *Y,float *Z,int N);	
static INLINE void F77__sdiff(float *X,int N);	
static INLINE void F77__sadd(float *X,float *Y,int N);
static INLINE void F77__ssubc(float C,float *Y,int N); 
static INLINE void F77__scsub(float C,float *Y,int N);
static INLINE void F77__ssqrvec(float *X,float *Y,int N);
static INLINE void F77__ssum(float *X,int N,float *ANS);
static INLINE void F77__isum(int *X,int N,int *ANS);
static INLINE float  F77__sdot(const int N,const float *x,const int incx,const float *y,const int incy);
static INLINE void  F77__sscal(const int N,const float alpha,float *X,const int incX);
static INLINE void  F77__iaddc(const int C,int *X,const int N);
static INLINE void fill_float32(rF32 value,rF32PTR dst,int N) { for (rI32 i=N; i > 0; i--) *dst++=value; }
void F77__CALL(RSSUMLOG)(float *x,int *N,float *ans);
void F77__CALL(RSPOWVEC)(float *x,float * pow,int *N);
void F77__CALL(RSLOGVEC)(float *x,int *N);
void F77__CALL(RSCOSVEC)(float *x,int *N);
void F77__CALL(RSSINVEC)(float *x,int *N);
void F77__CALL(RSSQRTVEC)(float *x,int *N);
void F77__CALL(RSSUB_I)(float *x,float *y,int *N);
void F77__CALL(RSSUB)(float *x,float *y,float *z,int *N);
void F77__CALL(RSDIFF)(float *x,int *N);
void F77__CALL(RSADD)(float *x,float *y,int *N);
void F77__CALL(RSSUBC)(float *C,float *y,int *N);
void F77__CALL(RSCSUB)(float *C,float *y,int *N);
void F77__CALL(RSSQRVEC)(float *x,float *y,int *N);
void F77__CALL(RSSUM)(float *x,int* N,float* ANS);
void F77__CALL(RISUM)(int *x,int* N,int* ANS);
void F77__CALL(RSMINIDX)(float *,int *N,float *val,int * idx);
void F77__CALL(RSMAXIDX)(float *,int *N,float *val,int * idx);
void F77__CALL(RSMNSD)(float * SX,int * N,float *MEAN,float* STD);
float F77__CALL(RSDOT)(const int *n,const float *x,const int *incx,const float *y,const int *incy);
void F77__CALL(RSSCAL)(const int *n,const float *alpha,float *x,const int *incx);
void F77__CALL(RIADDC)(const int *C,int *x,const int *N);
static INLINE void   F77__smaxidx(float *X,int N,float *val,int * idx) {
	F77__CALL(RSMAXIDX)(X,&N,val,idx);
}
static INLINE void F77__sminidx(float *X,int N,float *val,int * idx) {
	F77__CALL(RSMINIDX)(X,&N,val,idx);
}
static INLINE void  F77__ssumlog(float *X,int N,float *pSum){	
	F77__CALL(RSSUMLOG)(X,&N,pSum);
}
static INLINE void   F77__smnsd(float *X,int N,float *MEAN,float * STD){
	F77__CALL(RSMNSD)(X,&N,MEAN,STD);
}
static INLINE void  F77__spowvec(float *X,float pow,int N){
	F77__CALL(RSPOWVEC)(X,&pow,&N);
}
static INLINE void  F77__slogvec(float *X,int N){	
	F77__CALL(RSLOGVEC)(X,&N);
}
static INLINE void  F77__scosvec(float *X,int N){	
	F77__CALL(RSCOSVEC)(X,&N);
}
static INLINE void  F77__ssinvec(float *X,int N){
	F77__CALL(RSSINVEC)(X,&N);
}
static INLINE void F77__ssqrtvec(float *X,int N){	
	F77__CALL(RSSQRTVEC)(X,&N);
} 
static INLINE void  F77__ssub_i(float *X,float *Y,int N){	
	F77__CALL(RSSUB_I)(X,Y,&N);
}
static INLINE void  F77__ssub(float *X,float *Y,float *Z,int N){	
	F77__CALL(RSSUB)(X,Y,Z,&N);
}
static INLINE void  F77__sdiff(float *X,int N){	
	F77__CALL(RSDIFF)(X,&N);
}
static INLINE void  F77__sadd(float *X,float *Y,int N){
	F77__CALL(RSADD)(X,Y,&N);
}
static INLINE  void  F77__ssubc(float C,float *Y,int N){
	F77__CALL(RSSUBC)(&C,Y,&N);
} 
static INLINE void  F77__scsub(float C,float *Y,int N){
	F77__CALL(RSCSUB)(&C,Y,&N);
}
static INLINE void  F77__ssqrvec(float *X,float *Y,int N) {
	F77__CALL(RSSQRVEC)(X,Y,&N);
} 
static INLINE void  F77__ssum(float *X,int N,float *ANS){
	F77__CALL(RSSUM)(X,&N,ANS);
}
static INLINE void  F77__isum(int *X,int N,int *ANS){
	F77__CALL(RISUM)(X,&N,ANS);
}
static INLINE float  F77__sdot(const int N,const float *x,const int incx,const float *y,const int incy)
{
	return F77__CALL(RSDOT)(&N,x,&incx,y,&incy);
}
static INLINE void  F77__sscal(const int N,const float alpha,float *X,const int incX){
	F77__CALL(RSSCAL)(&N,&alpha,X,&incX);
}
static INLINE void  F77__iaddc(const int C,int *X,const int N){
	F77__CALL(RIADDC)(&C,X,&N);
}
