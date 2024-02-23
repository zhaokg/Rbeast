#pragma once
#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif

	typedef enum { CMP_LT, CMP_LE, CMP_GT, CMP_GE, CMP_EQ }CmpFlag;

	//https://stackoverflow.com/questions/57518041/are-extern-extern-c-and-extern-c-extern-allowed
	//https://stackoverflow.com/questions/61467251/extern-used-twice-in-c

	extern void (*i32_add_val_inplace)(const int C, const I32PTR X, const int N);
	extern I32(*i32_sum)(const I32PTR X, const int N);
	extern void (*f32_fill_val)(const F32 C, F32PTR X, int N);
	extern F32(*f32_sum)(const F32PTR X, int N);
	extern void (*f32_add_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
	extern void (*f32_sub_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
	extern void (*f32_add_vec_inplace)(const F32PTR SRC, const F32PTR DST, const int N);
	extern void (*f32_sub_vec_inplace)(const F32PTR SRC, F32PTR DST, int N);
	extern void (*f32_subrev_val_inplace)(const F32 C, F32PTR X, int N);
	extern void (*f32_add_val_inplace)(const F32 C, F32PTR X, int N);
	extern void (*f32_mul_val_inplace)(const F32 C, F32PTR X, const int N);
	extern void (*f32_mul_vec_inplace)(const F32PTR SRC, F32PTR DST, int N);
	extern void (*f32_mul_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
	extern F32(*f32_dot)(const F32PTR x, const F32PTR y, const int N);

	extern F32(*f32_dot2x1)(const F32PTR x, const F32PTR y, const F32PTR v, const int N, F32PTR res);
	extern void (*f32_dot2x2)(const F32PTR x1, const F32PTR x2, const F32PTR y1, const F32PTR y2, const int N, F32PTR res1, F32PTR res2);
	extern void (*f32_add_v_v2_vec_inplace)(const F32PTR SRC, const F32PTR x, F32PTR x2, int N);
	extern void (*f32_cos_vec_inplace)(const F32PTR X, const int N);
	extern void (*f32_sin_vec_inplace)(const F32PTR X, const int N);
	extern void (*f32_sincos_vec_inplace)(const F32PTR in_outsin, F32PTR outcos, const int N);
	extern void (*f32_pow_vec_inplace)(F32PTR X, F32 pow, int N);
	extern void (*f32_log_vec_inplace)(const F32PTR X, const int N);
	extern void (*f32_exp_vec_inplace)(const F32PTR X, const int N);
	extern void (*f32_sqrt_vec_inplace)(const F32PTR X, const int N);
	extern void (*f32_sqrt_vec)(const F32PTR X, F32PTR Y, int N);

	extern void(*f32_avgstd)(const F32PTR X, int N, F32PTR AVG, F32PTR STD);
	extern void (*f32_sx_sxx_to_avgstd_inplace)(F32PTR SX, F32PTR SXX, I32 Nsample, F32 scale, F32 offset, int N);
	// the index is zero-based
	extern I32(*f32_maxidx_slow)(const F32PTR  X, const int N, F32PTR val);
	// the index is zero-based
	extern I32(*f32_maxidx)(const F32PTR  X, const  int N, F32PTR val);
	extern I32(*f32_minidx)(const F32PTR  X, const int  N, F32PTR val);
	extern void (*f32_diff_back)(const F32PTR  X, F32PTR result, const int N);
	extern void (*f32_seq)(F32PTR p, F32 x0, F32 dX, int N);
	extern void (*i32_seq)(I32PTR p, I32 x0, I32 dX, int N);
	extern void (*f32_to_f64_inplace)(F32PTR data32, int N);
	extern void (*f64_to_f32_inplace)(F64PTR data64, int N);
	extern void (*i32_to_f32_scaleby_inplace)(I32PTR X, int N, F32 scale);
	extern void (*i32_increment_bycond_inplace)(I32PTR x, F32PTR cond, int N);
	extern void (*i32_increment_vec2_bycond_inplace)(I32PTR x, I32PTR y, F32PTR cond, int N);
	extern I32(*i08_sum_binvec)(U08PTR binvec, I32 N);

	extern void f32_cumsum_inplace(const F32PTR X, int N);
	extern void f32_cumsumsqr_inplace(const F32PTR X, int N);
	extern F32 f32_sumlog_slow(const F32PTR  X, const int N);
	extern F32 f32_sumlog(const F32PTR  X, const int N);
	extern I32 i08_find_nth_onebyte_binvec(U08PTR binvec, I32 N, I32 nth);
	extern I32 i08_find_nth_onebyte_binvec_v2(U08PTR binvec, I32 N, I32 numOneBytes, U32 rnd);
	extern I64 i08_sum(I08PTR x, int N);


	extern int  i32_insert_noduplicate(I32PTR x, I32 N, I32PTR Xnew, I32 Nnew);
	extern int  i32_unique_inplace(I32PTR x, int N);
	extern int  i32_exclude_inplace(I32PTR x, int N, I32PTR excludeList, I32 Nexclude);
	extern void f32_sumfilter(const F32PTR X, F32PTR Y, int N, int winSize);
	extern F32  f32_corr_rmse_nan(const F32PTR X, const F32PTR Y, int N, F32PTR rmse);
	extern void f32_truncate_inplace(const F32PTR X, F32 value, int N);


	I32 i32_maxidx(const I32PTR  X, const  int N, I32PTR val);
	I32 i32_minidx(const I32PTR  X, const  int N, I32PTR val);

	void f32_to_strided_f64(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset);
	void f32_to_strided_i64(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset);
	void f32_to_strided_f32(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset);
	void f32_to_strided_i32(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset);
	void f32_to_strided_i16(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset);

	void f32_from_strided_f64(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset);
	void f32_from_strided_i64(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset);
	void f32_from_strided_f32(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset);
	void f32_from_strided_i32(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset);
	void f32_from_strided_i16(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset);

	void f32_to_strided_mem(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset, DATA_TYPE dtype);
	void f32_from_strided_mem(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset, DATA_TYPE srcDataType);
	void arr_from_strided_mem(VOID_PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset, DATA_TYPE srcDstDataType);

	void f32_interp1dvec_cycled_inplace(F32PTR Y, int P, I32PTR goodIndices, int Pgood);
	void f32_rep_vec1d_upto_inplace(F32PTR Y, int P, int N);
	void f32_compute_seasonal_avg(F32PTR y, int N, int P, F32PTR mean, I32PTR NumGoodPtsPerTime);
	void f32_deseasonalize_inplace(F32PTR y, int N, int P,  F32PTR mean_tmp, I32PTR NumGoodPtsPerTime_tmp);

	F32 f32_nansum(F32PTR x, int N) ;
	F32 f32_nanmean(F32PTR x, int N, int* Ngood);
	F32 f32_absmax(F32PTR x, int N);

	I64  sub2ind(int* dims, int ndim, int* subs);
	void ind2sub(int* dims, int ndim, I64 ind, int* subs);
	int  ndarray_get1d_stride_offset(int* dims, int ndim, int* subs, int whichdim, I64* stride, I64* offset);
	void f32_get1d_from_ndarray(F32PTR dst, VOID_PTR src, int* dims, int ndim, int* subs, int whichdim, DATA_TYPE srcDtype);
	void f32_set1d_to_ndarray(F32PTR src, VOID_PTR dst, int* dims, int ndim, int* subs, int whichdim, DATA_TYPE dstDtype);
	void f32_get2d_from_ndarray(F32PTR dst, VOID_PTR src, int* dims, int ndim, int* subs, int d1, int d2, DATA_TYPE srcDtype);
	void f32_set2d_from_ndarray(F32PTR src, VOID_PTR dst, int* dims, int ndim, int* subs, int d1, int d2, DATA_TYPE dstDtype);



	void f32_set_value_to_nan(F32PTR a, I32 N, F32 missingValue);
	int f32_normalize_multicols_zeroout_nans(F32PTR Y, I32PTR BadRowIndices, I32 ldy, I32 N, I32 q, F32PTR mean, F32PTR sd);
	extern void f32_transpose_inplace(F32PTR Mat, I32 ROW, I32 COL);
	extern void i32_transpose_inplace(I32PTR Mat, I32 NROW, I32 NCOL);
	void i32_transpose_inplace_prev(I32PTR Mat, I32 NROW, I32 NCOL);
	void i32_transpose_inplace_prev_two_ends(I32PTR Mat, U64 NROW, U64 NCOL);
	extern void f32_fill_val_matrixdiag(F32PTR mat, const F32 value, I32 N);
	extern void f32_add_val_matrixdiag(F32PTR mat, const F32 value, I32 N);
	extern F32 f32_sum_matrixdiag(F32PTR mat, I32 N);
	extern F32 f32_abs_sum(F32PTR X, I32 N);

	extern void f32_mat_multirows_extract_set_by_scalar(F32PTR X, I32 ROW, I32 COL, F32PTR Xcopy, I32PTR RowIndices, I32 nRows, F32 newValue);
	extern void f32_mat_multirows_set_by_submat(F32PTR X, I32 ROW, I32 COL, F32PTR Xcopy, I32PTR RowIndices, I32 nRows);
	extern void f32_normalize_std_avg_inplace(F32PTR X, I32 N, F32PTR avg, F32PTR std);
	extern void f32_normalize_inplace(F32PTR X, I32 N);
	extern void f32_normalize_x_factor_inplace(F32PTR X, I32 N, F32 factor);
	extern I32  f32_find_nans(const F32PTR X, int N, I32PTR index);

	extern void  (*f32_hinge_neg)(const F32PTR X, const F32PTR Y, const F32 knot, const int N);
	extern void  (*f32_hinge_pos)(const F32PTR X, const F32PTR Y, const F32 knot, const int N);
	extern void  (*f32_step_neg)(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N);
	extern void  (*f32_step_pos)(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N);
	extern void  (*f32_axpy_inplace)(const F32 a, const F32PTR x, F32PTR y, const int N);


   //https://stackoverflow.com/questions/47450718/gcc7-2-argument-range-exceeds-maximum-object-size-9-7-werror-alloc-size-larg
   #define f32_copy(src, dst,N)  memcpy(dst,src, (U32)sizeof(F32)*(U32)(N))


	extern void (*f32_gemm_XtY2x1)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
	extern void (*f32_gemm_XtY2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);

	extern void (*f32_gemm_XY1x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
	extern void (*f32_gemm_XY2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
	extern void (*f32_gemm_XtYt2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
	extern void (*f32_gemm_XYt2x1)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);

	extern void(*f32_gemv_Xb)(int N, int K, F32PTR X, int lda, F32PTR b, F32PTR C);

	extern I32(*f32_findindex)(F32PTR  x, I32PTR indices, F32 value, int N, CmpFlag flag);
	extern void  (*f32_scatter_vec_byindex)(F32PTR  x, I32PTR indices, F32PTR values, int N);
	extern void (*f32_gatherVec_scatterVal_byindex)(F32PTR  x, I32PTR indices, F32PTR values, F32 newValue, int N);
	extern void (*f32_gather2Vec_scatterVal_byindex)(F32PTR  x, F32PTR  y, I32PTR indices, F32PTR values, F32 newValue, int N);
	extern void (*f32_scale_inplace)(const F32 gain, const F32 offset, const F32PTR x, const int N);
	extern void SetupVectorFunction_AVX2(void);
	extern void SetupVectorFunction_AVX512(void);
	extern void SetupVectorFunction_Generic(void);

	void print_funcs(void);

	/////////////////////////////////////////
	// THis is the older header of abc_generic.h, kept here to 
	// give a quick referene to the conversions between the C and fortran function interfaces
	////////////////////////////////////////////////

	/*

	 F77__smaxidx(F32* X, int N, F32PTR val, int* idx);
	 F77__sminidx(F32* X, int N, F32PTR val, int* idx);
	F77__ssumlog(F32* X, int N, F32PTR pSum);	//ippsSumLn_32f(const Ipp32f* pSrc, int len, Ipp32f* pSum))
	 F77__smnsd(F32* X, int N, F32PTR MEAN, F32PTR  STD);
	F77__spowvec(F32* X, F32 pow, int N);  //vmsPowx(const MKL_INT n,  const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
	F77__slogvec(F32* X, int N);	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__scosvec(F32* X, int N);	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__ssinvec(F32* X, int N);   //vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	F77__ssqrtvec(F32* X, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__ssub_i(F32* X, F32PTR Y, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__ssub(F32* X, F32PTR Y, F32PTR Z, int N);	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	xxxxxx void F77__sdiff(F32* X, int N);	//Z=Y-Z;
	F77__sadd(F32* X, F32PTR Y, int N);//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__ssubc(F32 C, F32PTR Y, int N); //ippsSubC_32f_I(val, pSrcDst, N)
	F77__scsub(F32 C, F32PTR Y, int N);//r_ippsSubCRev_32f_I(val, pSrcDst, N)
	F77__ssqrvec(F32* X, F32PTR Y, int N);//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	F77__ssum(F32* X, int N, F32PTR ANS);//ippsSum_32f(Src, N, Sum, PrecisionLevel);
	F77__isum(int* X, int N, int* ANS);//ippsSum_32s_Sfs(Src, N, Sum, ScaleFfactor);
	F77__sdot(const int N, const F32PTR x, const int incx, const F32PTR y, const int incy);
	F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX);
	F77__iaddc(const int C, int* X, const int N);


	void  i32_add_val_inplace(const int C, const I32PTR X, const int N);
	I32   i32_sum(const I32PTR X, const int N);

	void f32_fill_val(const F32 C, F32PTR X, int N);

	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	void f32_add_vec_inplace(const F32PTR SRC, const F32PTR DST, const int N);

	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	void f32_sub_vec_inplace(const F32PTR SRC, F32PTR DST, int N);

	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	void f32_sub_vec(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);

	//X+C
	void f32_add_val_inplace(const F32 C, F32PTR X, int N);

	//F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX);
	void f32_mul_val_inplace(const F32 C, F32PTR X, const int N);

	//F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX);
	void f32_mul_vec_inplace(const F32PTR SRC, F32PTR DST, int N);

	//F77__sscal(const int N, const F32 alpha, F32PTR X, const int incX);
	void f32_mul_vec(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
	//F77__sdot(const int N, const F32PTR x, const int incx, const F32PTR y, const int incy);
	F32 f32_dot(const F32PTR x, const F32PTR y, const int N);

	//ippsSum_32f(Src, N, Sum, PrecisionLevel);
	F32 f32_sum(const F32PTR X, int N);

	void f32_add_v_v2_vec_inplace(const F32PTR SRC, const F32PTR x, F32PTR x2, int N);

	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	void f32_cos_vec_inplace(const F32PTR X, const int N);
	void f32_sin_vec_inplace(const F32PTR X, const int N);
	//vmsPowx(const MKL_INT n,  const F32  a[], const F32   b, F32  r[], MKL_INT64 mode))
	void f32_pow_vec_inplace(F32* X, F32 pow, int N);
	//vmsSin(const MKL_INT n, const F32  a[], F32  r[], MKL_INT64 mode))
	void f32_log_vec_inplace(const F32PTR X, const int N);
	void f32_sqrt_vec_inplace(const F32PTR X, const int N);

	//ippsMul_32f, (const Ipp32f*  pSrc1, const Ipp32f*  pSrc2, Ipp32f*  pDst, int len))
	void f32_sqrt_vec(const F32PTR X, F32PTR Y, int N);
	//ippsSumLn_32f(const Ipp32f* pSrc, int len, Ipp32f* pSum))
	F32 f32_sumlog_slow(const F32PTR  X, const int N);

	F32 f32_sumlog(const F32PTR  X, const int N);

	//r_ippsSubCRev_32f_I(val, pSrcDst, N):  C-X
	void f32_subrev_val_inplace(const F32 C, F32PTR X, int N);

	//void  F77__smnsd(F32* X, int N, F32PTR MEAN, F32PTR  STD);
	F32 f32_avgstd(const F32PTR X, int N, F32PTR STD);
	void f32_sx_sxx_to_avgstd_inplace(F32PTR SX, F32PTR SXX, I32 Nsample, F32 scale, F32 offset, int N);
	// the index is zero-based
	I32 f32_maxidx_slow(const F32PTR  X, const int  N, F32PTR val);

	// the index is zero-based
	I32 f32_maxidx(const F32PTR  X, const int N, F32PTR val);

	// the index is zero-based
	//F77__sminidx(F32* X, int N, F32PTR val, int* idx);
	I32 f32_minidx(const F32PTR  X, const int N, F32PTR val);



	// Compute the 1-step backward difference:
	// First, get the backward difference: diff=[d1=NA,d2=x2-x1,...,d_N=x_n-x_(n-1)]
	// The first elemn of result (d1=y1-y0 where y0 is non-existent) is borrowed from its
	// immediate right neighor: d1 =d2=x2-x1, so that diff=(d1=y2-y1, d2=y2-y1, d3=d3-d2...)
	// diff[0] = diff[1]; //fill the leftmost

	void f32_diff_back(const F32PTR  X, F32PTR result, const  int N);

	void f32_seq(F32PTR p, F32 x0, F32 dX, int N);
	void f32_to_f64_inplace(F32PTR data32, int N);
	void i32_increment_bycond_inplace(I32PTR x, F32PTR cond, int N);
	void i32_to_f32_scaleby_inplace(I32PTR X, int N, F32 scale);
	I32 i08_sum_binvec(U08PTR binvec, I32 N);
	I32 i08_find_nth_onebyte_binvec(U08PTR binvec, I32 N, I32 nth);
	I32  i08_find_nth_onebyte_binvec_v2(U08PTR binvec, I32 N, I32 numOneBytes, U32 rnd);

	*/

#ifdef __cplusplus
}
#endif