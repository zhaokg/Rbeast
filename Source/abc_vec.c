#include <math.h>
#include <string.h>
#include "abc_000_warning.h"

#include "abc_datatype.h"
#include "abc_vec.h"
#include "abc_ide_util.h"

void (*i32_add_val_inplace)(const int C, const I32PTR X, const int N)   = NULL;
I32  (*i32_sum)(const I32PTR X, const int N)                            = NULL;
void (*f32_fill_val)(const F32 C, F32PTR X, int N);
F32  (*f32_sum)(const F32PTR X, int N)=NULL;
void (*f32_add_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
void (*f32_sub_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
void (*f32_add_vec_inplace)(const F32PTR SRC, const F32PTR DST, const int N);
void (*f32_sub_vec_inplace)(const F32PTR SRC, F32PTR DST, int N);
void (*f32_subrev_val_inplace)(const F32 C, F32PTR X, int N);
void (*f32_add_val_inplace)(const F32 C, F32PTR X, int N);
void (*f32_mul_val_inplace)(const F32 C, F32PTR X, const int N);
void (*f32_mul_vec_inplace)(const F32PTR SRC, F32PTR DST, int N);
void (*f32_mul_vec)(const F32PTR SRC1, const F32PTR SRC2, F32PTR DST, int N);
F32 (*f32_dot)(const F32PTR x, const F32PTR y, const int N);

F32 (*f32_dot2x1)(const F32PTR x, const F32PTR y, const F32PTR v, const int N, F32PTR res);

void (*f32_dot2x2)(const F32PTR x1, const F32PTR x2, const F32PTR y1, const F32PTR y2, const int N, F32PTR res1, F32PTR res2);
void (*f32_add_v_v2_vec_inplace)(const F32PTR SRC, const F32PTR x, F32PTR x2, int N);
void (*f32_cos_vec_inplace)(const F32PTR X, const int N);
void (*f32_sin_vec_inplace)(const F32PTR X, const int N);
void (*f32_sincos_vec_inplace)(const F32PTR in_outsin, F32PTR outcos, const int N);
void (*f32_pow_vec_inplace)(F32PTR X, F32 pow, int N);
void (*f32_log_vec_inplace)(const F32PTR X, const int N);
void (*f32_exp_vec_inplace)(const F32PTR X, const int N);
void (*f32_sqrt_vec_inplace)(const F32PTR X, const int N);
void (*f32_sqrt_vec)(const F32PTR X, F32PTR Y, int N);
//F32 (*f32_sumlog_slow)(const F32PTR  X, const int N);
//F32 (*f32_sumlog)(const F32PTR  X, const int N);
void  (*f32_avgstd)(const F32PTR X, int N, F32PTR AVG, F32PTR STD);
void (*f32_sx_sxx_to_avgstd_inplace)(F32PTR SX, F32PTR SXX, I32 Nsample, F32 scale, F32 offset, int N);
// the index is zero-based
I32 (*f32_maxidx_slow)(const F32PTR  X, const int N, F32PTR val);
// the index is zero-based
I32 (*f32_maxidx)(const F32PTR  X, const  int N, F32PTR val);
I32 (*f32_minidx)(const F32PTR  X, const int  N, F32PTR val);
void (*f32_diff_back)(const F32PTR  X, F32PTR result, const int N);
void (*f32_seq)(F32PTR p, F32 x0, F32 dX, int N);
void (*i32_seq)(I32PTR p, I32 x0, I32 dX, int N);
void (*f32_to_f64_inplace)(F32PTR data32, int N);
void (*f64_to_f32_inplace)(F64PTR data64, int N);
void (*i32_to_f32_scaleby_inplace)(I32PTR X, int N, F32 scale);
void (*i32_increment_bycond_inplace)(I32PTR x, F32PTR cond, int N);
void (*i32_increment_vec2_bycond_inplace)(I32PTR x, I32PTR y, F32PTR cond, int N);
//I32  (*i08_find_nth_onebyte_binvec)(U08PTR binvec, I32 N, I32 nth);
//I32  (*i08_find_nth_onebyte_binvec_v2)(U08PTR binvec, I32 N, I32 numOneBytes, U32 rnd);
I32  (*i08_sum_binvec)(U08PTR binvec, I32 N);

void (*f32_gemm_XtY2x1)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
void (*f32_gemm_XtY2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
void (*f32_gemm_XY1x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
void (*f32_gemm_XY2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
void (*f32_gemm_XtYt2x2)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);
void (*f32_gemm_XYt2x1)(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc);

void(*f32_gemv_Xb)(int N, int K, F32PTR X, int lda, F32PTR b, F32PTR C);

I32   (*f32_findindex)(F32PTR  x, I32PTR indices, F32 value, int N, CmpFlag flag);
void  (*f32_scatter_vec_byindex)(F32PTR  x, I32PTR indices, F32PTR values, int N);
void  (*f32_gatherVec_scatterVal_byindex)(F32PTR  x, I32PTR indices, F32PTR values, F32 newValue, int N);
void  (*f32_gather2Vec_scatterVal_byindex)(F32PTR  x, F32PTR  y, I32PTR indices, F32PTR values, F32 newValue, int N);
void  (*f32_scale_inplace)(const F32 gain, const F32 offset, const F32PTR x, const int N);

void  (*f32_hinge_neg)(const F32PTR X, const F32PTR Y, const F32 knot, const int N);
void  (*f32_hinge_pos)(const F32PTR X, const F32PTR Y, const F32 knot, const int N);

void  (*f32_step_neg)(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N);
void  (*f32_step_pos)(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N);
void  (*f32_axpy_inplace)(const F32 a, const F32PTR x, F32PTR y, const int N);

void print_funcs(void) {
	
	// stackoverflow.com/questions/34705680/warning-format-x-expects-argument-of-type-unsigned-int
	// stackoverflow.com/questions/14733761/printf-formatting-for-hexadecimal
	// warning: format ‘%x’ expects argument of type ‘unsigned int’ (-Wformat=)
	// "%#x": # will print a 0x, but this won't work if the value is 0 (e.g., print("%#x", 0))
	// Note that if i = 0;, the versions using %# will not include the 0x prefix.

	r_printf("\n\n"  );
	r_printf("%s:%8p\n", "i32_add_val_inplace", (void *) i32_add_val_inplace);
	r_printf("%s:%8p\n", "i32_sum", (void*)i32_sum);
	r_printf("%s:%8p\n", "f32_fill_val", (void*)f32_fill_val);
	r_printf("%s:%8p\n", "f32_sum", (void*)f32_sum);
	r_printf("%s:%8p\n", "f32_add_vec", (void*)f32_add_vec);
	r_printf("%s:%8p\n", "f32_sub_vec", (void*)f32_sub_vec);
	r_printf("%s:%8p\n", "f32_add_vec_inplace", (void*)f32_add_vec_inplace);
	r_printf("%s:%8p\n", "f32_sub_vec_inplace", (void*)f32_sub_vec_inplace);
	r_printf("%s:%8p\n", "f32_subrev_val_inplace", (void*)f32_subrev_val_inplace);
	r_printf("%s:%8p\n", "f32_add_val_inplace", (void*)f32_add_val_inplace);
	r_printf("%s:%8p\n", "f32_mul_val_inplace", (void*)f32_mul_val_inplace);
	r_printf("%s:%8p\n", "f32_mul_vec_inplace", (void*)f32_mul_vec_inplace);
	r_printf("%s:%8p\n", "f32_mul_vec", (void*)f32_mul_vec);
	r_printf("%s:%8p\n", "f32_dot", (void*)f32_dot);

	r_printf("%s:%8p\n", "*f32_dot2x1", (void*)f32_dot2x1);

	r_printf("%s:%8p\n", "f32_dot2x2", (void*)f32_dot2x2);
	r_printf("%s:%8p\n", "f32_add_v_v2_vec_inplace", (void*)f32_add_v_v2_vec_inplace);
	r_printf("%s:%8p\n", "f32_cos_vec_inplace", (void*)f32_cos_vec_inplace);
	r_printf("%s:%8p\n", "f32_sin_vec_inplace", (void*)f32_sin_vec_inplace);
	r_printf("%s:%8p\n", "f32_sincos_vec_inplace", (void*)f32_sincos_vec_inplace);
	r_printf("%s:%8p\n", "f32_pow_vec_inplace", (void*)f32_pow_vec_inplace);
	r_printf("%s:%8p\n", "f32_log_vec_inplace", (void*)f32_log_vec_inplace);
	r_printf("%s:%8p\n", "f32_exp_vec_inplace", (void*)f32_exp_vec_inplace);
	r_printf("%s:%8p\n", "f32_sqrt_vec_inplace", (void*)f32_sqrt_vec_inplace);
	r_printf("%s:%8p\n", "f32_sqrt_vec", (void*)f32_sqrt_vec);
	//F32 (*f32_sumlog_slow)(const F32PTR  X, const int N);
	//F32 (*f32_sumlog)(const F32PTR  X, const int N);
	r_printf("%s:%8p\n", "f32_avgstd", (void*)f32_avgstd);
	r_printf("%s:%8p\n", "f32_sx_sxx_to_avgstd_inplace", (void*)f32_sx_sxx_to_avgstd_inplace);
	// the index is zero-based
	r_printf("%s:%8p\n", "f32_maxidx_slow", (void*)f32_maxidx_slow);
	// the index is zero-based
	r_printf("%s:%8p\n", "f32_maxidx", (void*)f32_maxidx);
	r_printf("%s:%8p\n", "f32_minidx", (void*)f32_minidx);
	r_printf("%s:%8p\n", "f32_diff_back", (void*)f32_diff_back);
	r_printf("%s:%8p\n", "f32_seq", (void*)f32_seq);
	r_printf("%s:%8p\n", "f32_to_f64_inplace", (void*)f32_to_f64_inplace);
	r_printf("%s:%8p\n", "f64_to_f32_inplace", (void*)f64_to_f32_inplace);
	r_printf("%s:%8p\n", "i32_to_f32_scaleby_inplace", (void*)i32_to_f32_scaleby_inplace);
	r_printf("%s:%8p\n", "i32_increment_bycond_inplace", (void*)i32_increment_bycond_inplace);
	r_printf("%s:%8p\n", "i32_increment_vec2_bycond_inplace", (void*)i32_increment_vec2_bycond_inplace);
	//I32  (*i08_find_nth_onebyte_binvec)(U08PTR binvec, I32 N, I32 nth);
	//I32  (*i08_find_nth_onebyte_binvec_v2)(U08PTR binvec, I32 N, I32 numOneBytes, U32 rnd);
	r_printf("%s:%8p\n", "i08_sum_binvec", (void*)i08_sum_binvec);

	r_printf("%s:%8p\n", "f32_gemm_XtY2x1", (void*)f32_gemm_XtY2x1);
	r_printf("%s:%8p\n", "f32_gemm_XtY2x2", (void*)f32_gemm_XtY2x2);
	r_printf("%s:%8p\n", "f32_gemm_XY1x2", (void*)f32_gemm_XY1x2);
	r_printf("%s:%8p\n", "f32_gemm_XY2x2", (void*)f32_gemm_XY2x2);
	r_printf("%s:%8p\n", "f32_gemm_XtYt2x2", (void*)f32_gemm_XtYt2x2);
	r_printf("%s:%8p\n", "f32_gemm_XYt2x1", (void*)f32_gemm_XYt2x1);

	r_printf("%s:%8p\n", "f32_gemv_Xb", (void*)f32_gemv_Xb);

	r_printf("%s:%8p\n", "f32_findindex", (void*)f32_findindex);
	r_printf("%s:%8p\n", "f32_scatter_vec_byindex", (void*)f32_scatter_vec_byindex);
	r_printf("%s:%8p\n", "f32_gatherVec_scatterVal_byindex", (void*)f32_gatherVec_scatterVal_byindex);
	r_printf("%s:%8p\n", "f32_gather2Vec_scatterVal_byindex", (void*)f32_gather2Vec_scatterVal_byindex);
 

}
void f32_cumsum_inplace(const F32PTR X, int N) {
	#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	F32 csum = 0;
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		csum += X[i];     X[i]     = csum;
		csum += X[i + 1]; X[i + 1] = csum;
		csum += X[i + 2]; X[i + 2] = csum;
		csum += X[i + 3]; X[i + 3] = csum;
	}
	for (; i < N; ++i) {
		csum += X[i];     X[i] = csum;
	}

	#undef UNROLL_NUMBER
}
void f32_cumsumsqr_inplace(const F32PTR X, int N) {

	#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	F32 csum = 0;
	I32 i    = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		csum += X[i]* X[i];     X[i]     = csum;
		csum += X[i + 1]* X[i+1]; X[i + 1] = csum;
		csum += X[i + 2] * X[i+2]; X[i + 2] = csum;
		csum += X[i + 3] * X[i+3]; X[i + 3] = csum;
	}
	for (; i < N; ++i) {
		csum += X[i]*X[i];     X[i] = csum;
	}
}
void f32_sumfilter(const F32PTR X,  F32PTR Y, int N, int winSize) {

	if (winSize <= 1) {
		memcpy(Y, X, sizeof(F32) * N);
		return;
	}

	// Half the inter-changepoint distacne 	
	I32 wLeft  = winSize / 2;          // leftside half
	I32 wRight = (winSize - wLeft)-1;  // rightside half
	//there is an extra midpoint between wLeft and wRight

	F32 csumRIghtSide = 0;
	I32 Nadj          = min(wRight, N-1);
	for (int i = 1; i <= Nadj; ++i) {
		csumRIghtSide += X[i];
	}


	F32 csumLeftSide_Mid = 0;
	//0,1,.. wLeft-1, (Wleft=midpoint), 
	for (int i=0; i<= wLeft;++i){
		csumLeftSide_Mid  += X[i];
		Y[i]               = csumLeftSide_Mid + csumRIghtSide ;

		F32  newRightPoint = (i + wRight + 1) >= N ? 0 : X[i + wRight + 1];
		F32  oldRightPoint = (i +  1)         >= N ? 0 : X[i + 1];
		csumRIghtSide      +=  newRightPoint - oldRightPoint;
	}

	 F32 csumAll = Y[wLeft]; 
	for (int i = wLeft+1; i < N-wRight; ++i) {
		csumAll += X[i + wRight] - X[i-1 -  wLeft]; 
		Y[i]    = csumAll;		
	}

	csumAll = csumAll;
	for (int i = N - wRight; i < N; ++i) {
		csumAll -= X[i-wLeft-1];
		Y[i]     = csumAll;		
	}
 
	#undef UNROLL_NUMBER
}
F32  f32_corr_rmse_nan(const F32PTR X,  const F32PTR Y, int N, F32PTR rmse) {

	#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	F32 sumX  = 0;
	F32 sumXX = 0;
	F32 sumY  = 0;
	F32 sumYY = 0;
	F32 sumXY = 0;
	F32 DXY2   = 0;
	I32 n     = 0;	
	for (int i = 0; i < N; ++i) {
		I32 isGood = (X[i] == X[i]) & (Y[i] == Y[i]);
		n += isGood;
		F32 x = isGood ? X[i]:0;
		F32 y = isGood ? Y[i]:0;
		sumX  += x;   sumY+= y;
		sumXX += x*x; sumYY+= y*y;
		sumXY += x * y; 
		DXY2  += (x-y)*(x-y);
	}
	//https://www.statisticshowto.com/probability-and-statistics/correlation-coefficient-formula/

	F32 r=(n*sumXY-sumX*sumY)/sqrtf((n*sumXX-sumX* sumX)*(n * sumYY - sumY * sumY));
	
	*rmse = sqrtf(DXY2 / n);
 
	return r;

		#undef UNROLL_NUMBER
}
void f32_truncate_inplace(const F32PTR X,  F32 value, int N) {

	#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
 
	I32 i    = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		X[i]   = X[i] > value ? value : X[i];
		X[i+1] = X[i+1] > value ? value : X[i+1];
		X[i+2] = X[i+1] > value ? value : X[i+2];
		X[i+3] = X[i+1] > value ? value : X[i+3];
 
	}
	for (; i < N; ++i) {
		X[i] = X[i] > value ? value : X[i];
	}
 
		#undef UNROLL_NUMBER
}

I32  f32_find_nans(const F32PTR X, int N, I32PTR index ) {

#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)

	I32 nMissing = 0;
	I32 i        = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
				
		index[nMissing] = i;
		nMissing += X[i] != X[i];		
		index[nMissing] = i+1;
		nMissing += X[i + 1] != X[i + 1];		
		index[nMissing] = i+2;
		nMissing += X[i + 2] != X[i + 2];		
		index[nMissing] = i+3;	
		nMissing += X[i + 3] != X[i + 3];

	}
	for (; i < N; ++i) {
		index[nMissing] = i;
		nMissing += X[i] != X[i];
	}
	return nMissing;
   #undef UNROLL_NUMBER
}
//https://stackoverflow.com/questions/20065406/whats-the-largest-denormalized-and-normalized-number64bit-iee-754-1985
// Denormalized floating nubmer are 10x slower to process
//largest Normal: (1 + 1 / 2 + 1 / 2 ^ 2 + ... + 1 / 2 ^ 52) * 2 ^ 1023 = (2 - 2 ^ -52) * 2 ^ 1023 = 1.797...e + 308
//largst Denormal : (0 + 1 / 2 + 1 / 2 ^ 2 + ... + 1 / 2 ^ 52) * 2 ^ -1022 = (1 - 2 ^ -52) * 2 ^ -1022 = 2.225...e - 308

//https://stackoverflow.com/questions/32193791/single-precision-floating-point-format-range
// single precison
//The maximum positive number has the maximum mantissa 1.111111111111111111111112 and maximum non - infinite exponent 127. Thus 1.111111111111111111111112 x 2127 = (2 − 2−23) x 2127 ≈ 3.402 × 10^38 ≈ 2128.
//The minimum positive number has the non - zero mantissa 0.000000000000000000000012 and the minimum exponent −126 for subnormal / denormalized numbers.Thus 0.000000000000000000000012 × 2−126 = 2−23 × 2-126 = 2−149 ≈ 1.401x10−45.

//https://www.agner.org/optimize/blog/read.php?i=142
//Floating point underflow and denormal numbers
//Denormal numbers are floating point numbers that are coded in a non - normal way which is used when the value is close to underflow, according to the official IEEE 754 standard.Most processors are unable to handle floating point underflow, denormal numbers, and other special cases in the general floating point execution units.These special cases are typically handled by microcode exceptions at the cost of 150 - 200 clocks per instruction.The Sandy Bridge can handle many of these special cases in hardware without any penalty.In my tests, the cases of underflowand denormal numbers were handled just as fast as normal floating point numbers for addition, but not for multiplication.

//https://stackoverflow.com/questions/9314534/why-does-changing-0-1f-to-0-slow-down-performance-by-10x
//To demonstrate that this has everything to do with denormalized numbers, if we flush denormals to zero by adding this to the start of the code :
//_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
//https://randomascii.wordpress.com/2012/05/20/thats-not-normalthe-performance-of-odd-floats/

/*
Performance implications on SSE
Intel handles NaNsand infinities much better on their SSE FPUs than on their x87 FPUs.NaNsand infinities have long been handled at full speed on this floating - point unit.However denormals are still a problem.

On Core 2 processors the worst - case I have measured is a 175 times slowdown, on SSE additionand multiplication.

On SandyBridge Intel has fixed this for addition – I was unable to produce any slowdown on ‘addps’ instructions.However SSE multiplication(‘mulps’) on Sandybridge has about a 140 cycle penalty if one of the inputs or results is a denormal.
*/

#include "math.h"

// the index is zero-based
I32 i32_maxidx(const I32PTR  X, const  int N, I32PTR val) {

	#define UNROLL_NUMBER  2 // it would be even slower if changed to 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;

	I32 maxVal = X[0];
	I32 maxIdx = 0;
	I32 i = 0;
	for (; i < regularPart; i+= UNROLL_NUMBER) {
		I32 val;
		I32 idx = X[i + 1] > X[i] ? (val=X[i+1],i + 1) : (val = X[i], i);
		if (maxVal < val) {
			maxVal = val;
			maxIdx = idx;
		}		
	}

	for (; i < N; i++) {
		if (maxVal < X[i]) {
			maxVal = X[i];
			maxIdx = i;
		}
	}

	*val = maxVal;
	return maxIdx;

	#undef UNROLL_NUMBER
}

// the index is zero-based
//F77__sminidx(F32* X, int N, F32PTR val, int* idx);
I32 i32_minidx(const I32PTR  X, const int  N, I32PTR val) {

	#define UNROLL_NUMBER  2 // it would be even slower if changed to 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;

	I32 minVal = X[0];
	I32 minIdx = 0;
	I32 i = 0;
	for (; i < regularPart; i+= UNROLL_NUMBER) {
		I32 val;
		I32 idx = X[i + 1] < X[i] ? (val=X[i+1],i + 1) : (val = X[i], i);
		if (minVal > val) {
			minVal = val;
			minIdx = idx;
		}		
	}

	for (; i < N; i++) {
		if (minVal > X[i]) {
			minVal = X[i];
			minIdx = i;
		}
	}

	*val = minVal;
	return minIdx;

	#undef UNROLL_NUMBER
}


F32 f32_sumlog(const F32PTR  X, const int N) {
// first take the products of multiple elements to reduce the number of calling logf
	F64 sumlog	= 0;
	F64 cumprod = 1.;
	for (I32 i = 0; i < N; i++) {		
		F64 cumprod_new = cumprod * X[i];
		if (cumprod_new > 1e-307 && cumprod_new < 1e308) { /*largest denormal: 2.225...e-308*/
			cumprod = cumprod_new;
		} else {
			sumlog += log(cumprod);
			cumprod = X[i];
		}
	}
	sumlog += logf(cumprod);
	return sumlog;
}

I32  i08_find_nth_onebyte_binvec(U08PTR binvec, I32 N, I32 nth)
{
	//binvec: 0,0,1,1,0,1,0,0,0: an array with elements valued at either 0 or 1
	
	//N must be a multiplifer of 16.
	I32 nthPos;

	//Cumulatively add up to j: IDX_tmp until reaching GoodPos_tmp :i
	I32 count   = 0;
	I32 failed  = 1;
	{
		I32 i=0, N16 = N / 16;		
		I32 deltaSum=0;
		for (i = 0; i < N16; i++) {
			I64 sum = *((I64PTR)binvec) + *((I64PTR)binvec + 1);
			*((I32PTR)&sum) += *((I32PTR)&sum + 1);
			*((I16PTR)&sum) += *((I16PTR)&sum + 1);
			*((I08PTR)&sum) += *((I08PTR)&sum + 1);
			deltaSum = *((I08PTR)&sum);
			count	+= deltaSum;
			if (count >= nth) { 
				failed = 0;
				break; 
			}

			binvec += 16;
		}//newKnot is the newly chosen breakpoint	

		count -= deltaSum;
		nthPos = i * 16;
	}

	if (failed) {
		return 0xffffffff;
	}

	{
		I32 j;
		for (j = 0; j < 16; j++) {
			count += binvec[j];
			if (count == nth) break;
		}
		nthPos += (j + 1);
	}

	//nthPos is a 1-based index
	return nthPos;
}

#include "stdio.h"

I32  i08_find_nth_onebyte_binvec_v2(U08PTR binvec, I32 N, I32 numOneBytes, U32 rnd)
{
	static int I1 = 0, I2 = 0;
	/***************************************/
	// First, cnnvert rnd to an random index within 0 to (N-1)
	// check if the assocaite byte is 1 or not. If yes, return this index
	
	I32 nth = rnd % N;	
	if (binvec[nth]) {
		I1++;
		return nth + 1;
	}

	// If not, go to the more intesnsve part
	nth = (rnd % numOneBytes) + 1;
	I2++;

	if (I2 % 100 == 0) {
		r_printf("%d %d\n",I1,I2);
	}
	/***************************************/


	//binvec: 0,0,1,1,0,1,0,0,0: an array with elements valued at either 0 or 1

	//N must be a multiplifer of 16.
	I32 nthPos;

	//Cumulatively add up to j: IDX_tmp until reaching GoodPos_tmp :i
	I32 count = 0;
	{
		I32 i = 0, N16 = N / 16;
		I32 deltaSum = 0;
		for (i = 0; i < N16; i++) {
			I64 sum = *((I64PTR)binvec) + *((I64PTR)binvec + 1);
			*((I32PTR)&sum) += *((I32PTR)&sum + 1);
			*((I16PTR)&sum) += *((I16PTR)&sum + 1);
			*((I08PTR)&sum) += *((I08PTR)&sum + 1);
			deltaSum = *((I08PTR)&sum);
			count += deltaSum;
			if (count >= nth) break;
			binvec += 16;
		}//newKnot is the newly chosen breakpoint	

		count -= deltaSum;
		nthPos = i * 16;
	}
	{
		I32 j;
		for (j = 0; j < 16; j++) {
			count += binvec[j];
			if (count == nth) break;
		}
		nthPos += (j + 1);
	}

	return nthPos;
}

I64 i08_sum(I08PTR x, int N) {

#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i   = 0;
	I64 sum = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) sum += (I64) x[i] + (I64)x[i+1] + (I64)x[i+2] + (I64)x[i+3];
	for (; i < N; ++i)  sum += (I64) x[i];
	return sum;
#undef UNROLL_NUMBER
}

int i32_insert_noduplicate(I32PTR x, I32 N, I32PTR Xnew, I32 Nnew) {

	for (int i = 0; i < Nnew; i++) {
		int newvalue = Xnew[i];
		int matched = _False_;
		for (int j = 0; j < N; j++) {
			if (x[j] == newvalue) {
				matched = _True_;
				break;
			}
		}
		if (!matched) {
			x[N++] = newvalue;
		}
	}

	// the number of elements of the new vector
	return N;
}

int i32_unique_inplace(I32PTR x, int N) {
	int Nuniq = 0;
	int next  = 0;
	while (next < N) {
		int xcur = x[next];
		// first, move Inext foward to check if consectuve elements are the same
		while (next < (N - 1) && x[next + 1] == xcur) {
			next++;
		}

		int matched = 0;
		for (int i = 0; i < Nuniq; ++i) {
			if (x[i] == xcur) {
				matched = 1;
				break;
			}
		}
		if (!matched) {
			x[Nuniq++] = xcur;
		}
		next++;
	}

	return Nuniq;
}

int i32_exclude_inplace(I32PTR x, int N, I32PTR excludeList, I32 Nexclude) {

	if (x == NULL || excludeList == NULL) return N;

	for (int i = 0; i < Nexclude && N>0; i++) {
		int xcur = excludeList[i];
		for (int j = 0; j < N; j++) {
			if (x[j] == xcur) {
				// find a match; move the last element to fill the pos
				x[j] = x[N - 1];
				N--;
				break;
			}
		}
	}
	return N;
}


void f32_fill_val_matrixdiag(F32PTR mat, const F32 value, I32 N) {
	for (int i = 0; i < N; ++i) {
		*mat = value;
		mat += (N + 1);
	}
}

void f32_add_val_matrixdiag(F32PTR mat, const F32 value, I32 N) {
	for (int i = 0; i < N; ++i) {
		(*mat) += value;
		mat    += (N + 1);
	}
}

F32 f32_sum_matrixdiag(F32PTR mat,  I32 N) {
	F64 sum = 0;
	for (int i = 0; i < N; ++i) {
		sum += *mat;
		mat += (N + 1);		
	}
	return sum;
}

F32 f32_abs_sum(F32PTR X, I32 N) {
	F64 sum = 0;
	for (int i = 0; i < N; ++i) {
		sum += fabs(X[i]);		 
	}
	return sum;
}

void f32_mat_multirows_extract_set_by_scalar(F32PTR X, I32 ROW, I32 COL, F32PTR Xcopy, I32PTR RowIndices, I32 nRows, F32 newValue) {
	// equivalent to Matlab: Xcopy=X(idx,:);X(idx,:)=val;
	int i = 0;
	for (; i < COL - 1; i += 2) f32_gather2Vec_scatterVal_byindex(X + i * ROW, X + i * ROW + ROW, (I32PTR) RowIndices,        Xcopy + i * nRows, newValue, nRows);
	if (i < COL)          		f32_gatherVec_scatterVal_byindex( X + i * ROW, (I32PTR)RowIndices,        Xcopy + i * nRows, newValue, nRows);
}

void f32_mat_multirows_set_by_submat(F32PTR X, I32 ROW, I32 COL, F32PTR Xcopy, I32PTR RowIndices, I32 nRows) {
	// equivalent to Matlab: X(idx,:)=Xcopy	
	for (int i = 0; i < COL; ++i) f32_scatter_vec_byindex(X + i * ROW, (I32PTR)RowIndices, Xcopy + nRows * i, nRows);
}

void f32_to_strided_f64(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset) {
	dst = (F64PTR)dst + dstOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		*(F64PTR)dst = src[i];
		*((F64PTR)dst + stride) = src[i + 1];
		*((F64PTR)dst + 2 * stride) = src[i + 2];
		*((F64PTR)dst + 3 * stride) = src[i + 3];
		dst = (F64PTR)dst + 4 * stride;
	}
	for (; i < N; ++i) {
		*(F64PTR)dst = src[i];
		dst = (F64PTR)dst + stride;
	}
#undef UNROLL_NUMBER
}

void f32_to_strided_i64(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset) {
	dst = (I64PTR)dst + dstOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		*(I64PTR)dst                = src[i];
		*((I64PTR)dst + stride)     = src[i + 1];
		*((I64PTR)dst + 2 * stride) = src[i + 2];
		*((I64PTR)dst + 3 * stride) = src[i + 3];
		dst = (I64PTR)dst + 4 * stride;
	}
	for (; i < N; ++i) {
		*(I64PTR)dst = src[i];
		dst = (I64PTR)dst + stride;
	}
}

void f32_to_strided_f32(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset) {
	dst = (F32PTR)dst + dstOffset;
	if (stride == 1)
		memcpy(dst, src, sizeof(F32) * N);
	else {
		#define UNROLL_NUMBER 4
		const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
		// (AND-ing with -vectorsize will round down to nearest lower multiple of 
		// vectorsize. This works only if vectorsize is a power of 2)
		I32 i = 0;
		for (; i < regularPart; i += UNROLL_NUMBER) {
			*(F32PTR)dst              = src[i];			*( (F32PTR)dst+stride)    = src[i+1];
			*( (F32PTR)dst+2*stride)  = src[i+2];		*( (F32PTR)dst+3*stride)  = src[i+3];
			dst = (F32PTR)dst + 4*stride;
		}
		for (; i < N; ++i) {
			*(F32PTR)dst = src[i];			
			dst = (F32PTR)dst + stride;
		}	 	 
	} //stide ==1 

#undef UNROLL_NUMBER
}

void f32_to_strided_i32(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset) {
	dst = (I32PTR)dst + dstOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		*(I32PTR)dst = src[i];			*((I32PTR)dst + stride) = src[i + 1];
		*((I32PTR)dst + 2 * stride) = src[i + 2];		*((I32PTR)dst + 3 * stride) = src[i + 3];
		dst = (I32PTR)dst + 4 * stride;
	}
	for (; i < N; ++i) {
		*(I32PTR)dst = src[i];
		dst = (I32PTR)dst + stride;
	} 
#undef UNROLL_NUMBER
}
void f32_to_strided_i16(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset) {
	dst = (I16PTR)dst + dstOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		*(I16PTR)dst                = src[i];		
		*((I16PTR)dst + stride) = src[i + 1];
		*((I16PTR)dst + 2 * stride) = src[i + 2];	
		*((I16PTR)dst + 3 * stride) = src[i + 3];
		dst = (I16PTR)dst + 4 * stride;
	}
	for (; i < N; ++i) {
		*(I16PTR)dst = src[i];
		dst = (I16PTR)dst + stride;
	} 
#undef UNROLL_NUMBER
}


 /*
	// If omissonVaue is not an NaN
	rF32 badValue = missingValue;
	rF32 nan = getNaN();
	for (rI32 i = N; i > 0; i--) {
		rF32 y = *data;
		data += srcStride;
		*dst++ = fabsf(y - badValue) > 1e-10 ? y : nan;
	}
 */
 

void f32_from_strided_f64(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset) { 

	src = (F64PTR)src + srcOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		dst[i]     = (F32) *(F64PTR)src;
		dst[i + 1] = (F32) *((F64PTR)src + srcStride);
		dst[i + 2] = (F32) *((F64PTR)src + 2 * srcStride);
		dst[i + 3] = (F32) *((F64PTR)src + 3 * srcStride);

		src = (F64PTR)src + 4 * srcStride;
	}
	for (; i < N; ++i) {
		dst[i] = (F32) *(F64PTR)src;
		src = (F64PTR)src + srcStride;
	}
#undef UNROLL_NUMBER		 
}

void f32_from_strided_i64(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset) { 

	src = (I64PTR)src + srcOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		dst[i]     = *(I64PTR)src;
		dst[i + 1] = *((I64PTR)src + srcStride);
		dst[i + 2] = *((I64PTR)src + 2 * srcStride);
		dst[i + 3] = *((I64PTR)src + 3 * srcStride);

		src = (I64PTR)src + 4 * srcStride;
	}
	for (; i < N; ++i) {
		dst[i] = *(I64PTR)src;
		src = (I64PTR)src + srcStride;
	}
#undef UNROLL_NUMBER	 
}

void f32_from_strided_f32(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset) {

	src = (F32PTR)src + srcOffset;
	if (srcStride == 1) {
		memcpy(dst, src, sizeof(F32) * N);
	}
	else {
#define UNROLL_NUMBER 4
		const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
		// (AND-ing with -vectorsize will round down to nearest lower multiple of 
		// vectorsize. This works only if vectorsize is a power of 2)
		I32 i = 0;
		for (; i < regularPart; i += UNROLL_NUMBER) {
			dst[i] = *(F32PTR)src;
			dst[i + 1] = *((F32PTR)src + srcStride);
			dst[i + 2] = *((F32PTR)src + 2 * srcStride);
			dst[i + 3] = *((F32PTR)src + 3 * srcStride);

			src = (F32PTR)src + 4 * srcStride;
		}
		for (; i < N; ++i) {
			dst[i] = *(F32PTR)src;
			src = (F32PTR)src + srcStride;
		}
	}
#undef UNROLL_NUMBER
}
void f32_from_strided_i32(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset)
{ 
	src = (I32PTR)src + srcOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		dst[i]     = *(I32PTR)src;
		dst[i + 1] = *((I32PTR)src + srcStride);
		dst[i + 2] = *((I32PTR)src + 2 * srcStride);
		dst[i + 3] = *((I32PTR)src + 3 * srcStride);

		src = (I32PTR)src + 4 * srcStride;
	}
	for (; i < N; ++i) {
		dst[i] = *(I32PTR)src;
		src = (I32PTR)src + srcStride;
	}
#undef UNROLL_NUMBER	 
}

void f32_from_strided_i16(F32PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset)
{ 
	src = (I16PTR)src + srcOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		dst[i]     = *(I16PTR)src;
		dst[i + 1] = *((I16PTR)src + srcStride);
		dst[i + 2] = *((I16PTR)src + 2 * srcStride);
		dst[i + 3] = *((I16PTR)src + 3 * srcStride);

		src = (I16PTR)src + 4 * srcStride;
	}
	for (; i < N; ++i) {
		dst[i] = *(I16PTR)src;
		src = (I16PTR)src + srcStride;
	}
#undef UNROLL_NUMBER	 
}	
 

void f64_from_strided_f64(F64PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset) { 

	src = (F64PTR)src + srcOffset;
	#define UNROLL_NUMBER 4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;
	// (AND-ing with -vectorsize will round down to nearest lower multiple of 
	// vectorsize. This works only if vectorsize is a power of 2)
	I32 i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		dst[i]     = *(F64PTR)src;
		dst[i + 1] = *((F64PTR)src + srcStride);
		dst[i + 2] = *((F64PTR)src + 2 * srcStride);
		dst[i + 3] = *((F64PTR)src + 3 * srcStride);
		src = (F64PTR)src + 4 * srcStride;
	}
	for (; i < N; ++i) {
		dst[i] = *(F64PTR)src;
		src = (F64PTR)src + srcStride;
	}
#undef UNROLL_NUMBER	 
}

void f32_to_strided_mem(F32PTR src, VOID_PTR dst, I64 N, I64 stride, I64 dstOffset, DATA_TYPE dstDtype) {
	if (dstDtype == DATA_FLOAT) 		    f32_to_strided_f32(src, dst, N, stride, dstOffset);
	else if (dstDtype == DATA_DOUBLE)    	f32_to_strided_f64(src, dst, N, stride, dstOffset);
	else if (dstDtype == DATA_INT64)    	f32_to_strided_i64(src, dst, N, stride, dstOffset);
	else if (dstDtype == DATA_INT32)    	f32_to_strided_i32(src, dst, N, stride, dstOffset);
	else if (dstDtype == DATA_INT16)    	f32_to_strided_i16(src, dst, N, stride, dstOffset);
}

void f32_from_strided_mem(F32PTR dst,  VOID_PTR src, int N, int srcStride, int srcOffset, DATA_TYPE srcDataType) {

	if      (srcDataType== DATA_FLOAT) 	   f32_from_strided_f32(dst, src, N, srcStride, srcOffset);
	else if (srcDataType == DATA_DOUBLE)   f32_from_strided_f64(dst, src, N, srcStride, srcOffset);
	else if (srcDataType == DATA_INT64)    f32_from_strided_i64(dst, src, N, srcStride, srcOffset);
	else if (srcDataType == DATA_INT32) {
		f32_from_strided_i32(dst, src, N, srcStride, srcOffset);
		#if R_INTERFACE==1
		//https://www.markvanderloo.eu/yaRb/2012/07/08/representation-of-numerical-nas-in-r-and-the-1954-enigma/
		// In R, the smallest negative interger is treated as NA
		I32 INTMIN    = 0x80000001;
		F32 NAinteger = (F32)(INTMIN);  //the most negative integer
		f32_set_value_to_nan(dst, N, NAinteger);
		#endif
	}
	else if (srcDataType == DATA_INT16)	   f32_from_strided_i16(dst, src, N, srcStride, srcOffset);

}

void arr_from_strided_mem(VOID_PTR dst, VOID_PTR src, int N, int srcStride, int srcOffset, DATA_TYPE srcDstDataType) {

	if (srcDstDataType == DATA_FLOAT) {
		f32_from_strided_f32(dst, src, N, srcStride, srcOffset);
	} else if (srcDstDataType == DATA_DOUBLE) {
		f64_from_strided_f64(dst, src, N, srcStride, srcOffset);
	}

}

void f32_set_value_to_nan(F32PTR a, I32 N, F32 missingValue) {

	if (missingValue != missingValue)		
		return;
		
	/* TODO: THis is terribly wrong bcz the first element a[0] is never visited and a[n] is accidently visted
	  for (I32 i = N; i > 0; i--) {	a[i] = fabsf(a[i]- missingValue) < 1e-10 || IsInf(a[i]) ?  nan: a[i] ;} */

	// If omissonVaue is not an NaN
	F32 nan = getNaN();
	for (I32 i = N-1; i >= 0; i--) {
		a[i] = fabsf(a[i] - missingValue) < 1e-10 || IsInf(a[i]) ? nan : a[i];
	}

}
 
int f32_normalize_multicols_zeroout_nans(F32PTR Y, I32PTR BadRowIndices, I32 ldy, I32 N, I32 q, F32PTR mean, F32PTR sd) {
	// BadRowIndices must be at least of length N.

    /***********************************************/
	// If only 1 col needs to be normalized
	/***********************************************/
	if (q == 1) {		
		I32     nMissing = 0;
		F64     sumY     = 0,  sumYY = 0;
		for (I32 i = 0; i < N; i++) {
			BadRowIndices[nMissing] = i;
			nMissing                = nMissing + (Y[i] != Y[i] || IsInf(Y[i])  );

			F32  y = (Y[i] != Y[i]  || IsInf(Y[i]) ) ? 0.f : Y[i]; // Is NAN or Inf
 			sumY  += y;
			sumYY += y*y;		 
			
		}
		I32 n      = N - nMissing;
		F32 MEAN32 = sumY / n;
		F64 SD32   = (sumYY / n - MEAN32 * MEAN32); //https://www.mun.ca/biology/scarr/Simplified_calculation_of_variance.html
		SD32       = SD32 > 0 ? sqrt(SD32) : 1.f;

		I32   jOmit = 0;
		for (I32 i = 0; i < N; ++i) {
			if (jOmit < nMissing && i == BadRowIndices[jOmit]) {
				Y[i] = 0.f; // fill bad points with zeros
				jOmit++;
			}  else {
				Y[i] = (Y[i] - MEAN32) * (1. / SD32);
			}
		}
		mean[0] = MEAN32;
		sd[0]    = SD32;
		return nMissing;
	}


	/***********************************************/
	// If there are multiple colums to be normalized (q>1)
	/***********************************************/

	I32PTR  isNANValue = (I32PTR)BadRowIndices;
	memset(isNANValue,0,sizeof(I32)* N);
	for (I32 i = 0; i < q; i++) {
		for (I32 j = 0; j < N; j++) { 	 
			// Is NAN or Inf
			if ( Y[j]!= Y[j] || IsInf(Y[i]) ) 	isNANValue[j] = 1;
		}
		Y += ldy;
	}
	Y -= ldy * q;

    // Find the indices from the binary vector isNNAValue
	// isNANValue and BadRowindices point to the same buf
	I32 nMissing = 0; 	
	for (I32 j = 0; j < N; j++) {
		I32 keep = BadRowIndices[j];
		BadRowIndices[nMissing] = j;
		nMissing               += keep;
	}


	I32 n = N - nMissing;
	for (I32 col = 0; col < q; col++) {
		F64 sumY  = 0;
		F64 sumYY = 0;
		I32 jOmit = 0;
		for (I32 i = 0; i < N; ++i) {
			if (jOmit < nMissing && i == BadRowIndices[jOmit]) {
				jOmit++;
			} else {
				sumY  += Y[i];
				sumYY += Y[i] * Y[i];
			}
		}

		F32 MEAN32 = sumY / n;
		F64 SD32 = (sumYY / n - MEAN32 * MEAN32); //https://www.mun.ca/biology/scarr/Simplified_calculation_of_variance.html
		SD32 = SD32 > 0 ? sqrt(SD32) : 1.f;

		jOmit = 0;
		for (I32 i = 0; i < N; ++i) {
			if (jOmit < nMissing && i == BadRowIndices[jOmit]) {
				Y[i] = 0.f; // fill bad points with zeros
				jOmit++;
			}	else {
				Y[i] = (Y[i] - MEAN32) * (1. / SD32);
			}
		}
		mean[col] = MEAN32;
		sd[col]   = SD32;


		Y = Y + ldy;

	}

	return nMissing;
}


void f32_normalize_std_avg_inplace(F32PTR X, I32 N, F32PTR avg, F32PTR std) {
	f32_avgstd(X, N, avg, std);
	F32 gain    = 1 / std[0];
	F32 offset = -avg[0] / std[0];  //(x-m)/s= 1/s *m - m/s
	f32_scale_inplace(gain, offset, X, N);
}
void f32_normalize_inplace(F32PTR X, I32 N) {
	F32 avg, std;
	f32_avgstd(X, N, &avg, &std);
	F32 gain   = 1 / std;
	F32 offset = -avg/std;  //(x-m)/s= 1/s *m - m/s
	f32_scale_inplace(gain, offset, X, N);
}
void f32_normalize_x_factor_inplace(F32PTR X, I32 N, F32 factor) {
	F32 avg, std;
	f32_avgstd(X, N, &avg, &std);
	F32 gain   = factor/std;
	F32 offset = -factor*avg/std;  //(x-m)*f/s= f/s *x - mf*/s
	f32_scale_inplace(gain, offset, X, N);
}

void f32_interp1dvec_cycled_inplace(F32PTR Y, int P, I32PTR goodIndices, int Pgood) {
	// Y has all the values filled at the goodIndices. The function predicts the missing 
	// values at the bad points, assuming a periodic boundary condition

	int lastIdx = goodIndices[Pgood - 1];
	int lastGoodIdx = lastIdx - P;
	for (int j = 0; j < Pgood; j++) {
		int curGoodIdx = goodIndices[j];
		if ((lastGoodIdx + 1) < curGoodIdx) {
			// Bad vaules located
			F32 curv = Y[curGoodIdx];
			F32 lastv = Y[lastGoodIdx >= 0 ? lastGoodIdx : lastGoodIdx + P];
			F32 full = curGoodIdx - lastGoodIdx;
			for (int k = lastGoodIdx + 1; k < curGoodIdx; k++) {
				Y[k >= 0 ? k : k + P] = lastv * (curGoodIdx - k) / full + curv * (k - lastGoodIdx) / full;
			}
		}
		lastGoodIdx = curGoodIdx;
	}
}


void f32_rep_vec1d_upto_inplace(F32PTR Y, int P, int N) {
	// Y has P values initially, and repeat it up to a length of N
	// N may be not an integer multiplicity of P
    
	int Ncylce_int = N / P;
	for (int i = 1; i < Ncylce_int; ++i) {
		memcpy(Y + i * P, Y, P * sizeof(F32));
	}

	int j = 0;
	for (int i = Ncylce_int * P; i < N; i++) {
		Y[i] = Y[j++];
	}
}

void f32_compute_seasonal_avg(F32PTR y, int N, int P, F32PTR mean, I32PTR NumGoodPtsPerTime) {
	// mean and NUmgoodPtsPer time must at least have a length of P
	// NumGoodPtsPerTime cannot be NULL
	memset(NumGoodPtsPerTime, 0, P * sizeof(I32));
		
	if ( mean == NULL) {
		int p = 0;
		for (int i = 0; i < N; i++) {
			NumGoodPtsPerTime[p] += (y[i] == y[i]);
			p++;
			p = (p == P) ? 0 : p;
		}
	} else {
		memset(mean, 0, P * sizeof(F32));;
		int p = 0;
		for (int i = 0; i < N; i++) {
			NumGoodPtsPerTime[p] += (y[i] == y[i]);
			mean[p]              += (y[i] == y[i])? y[i]:0.f;
			p++;
			p = (p == P) ? 0 : p;
		}

		for (int i = 0; i < P; i++) {
			mean[i] = (NumGoodPtsPerTime[i] > 0) ? mean[i] / NumGoodPtsPerTime[i] : getNaN();
		}	
	}

}

F32 f32_nansum(F32PTR x, int N) {

	#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;    

	F32  sum = 0;
	int  i   = 0;	
	for (; i < regularPart; i += UNROLL_NUMBER) {
		sum += x[i] == x[i] ? x[i] : 0;
		sum += x[i + 1] == x[i + 1] ? x[i + 1] : 0;
		sum += x[i + 2] == x[i + 2] ? x[i + 2] : 0;
		sum += x[i + 3] == x[i + 3] ? x[i + 3] : 0;
	}
		
	for (; i < N; ++i) {
		sum += x[i] == x[i] ? x[i] : 0;
	}		 
	return sum;
   #undef UNROLL_NUMBER 
}

F32 f32_nanmean(F32PTR x, int N, int * Ngood) {

#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;    

	F32  sum = 0;
	int  NgoodNum = 0;
	int  i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		sum       += x[i] == x[i] ? x[i] : 0;     
		NgoodNum += x[i] == x[i];

		sum += x[i + 1] == x[i + 1] ? x[i + 1] : 0;
		NgoodNum += x[i+1] == x[i+1];

		sum += x[i + 2] == x[i + 2] ? x[i + 2] : 0;
		NgoodNum += x[i + 2] == x[i + 2];

		sum += x[i + 3] == x[i + 3] ? x[i + 3] : 0;
		NgoodNum += x[i + 3] == x[i + 3];
	}

	for (; i < N; ++i) {
		sum       += x[i] == x[i] ? x[i] : 0;
		NgoodNum += x[i] == x[i];
	}
	if (Ngood) Ngood[0] = NgoodNum;
	return sum / NgoodNum;;
#undef UNROLL_NUMBER 
}

F32 f32_absmax(F32PTR x, int N) {

#define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;    

	F32  maxVal  = 0;
	
	int  i = 0;
	for (; i < regularPart; i += UNROLL_NUMBER) {
		F32 a = max(fabsf(x[i]), fabsf(x[i+1]));
		F32 b = max(fabsf(x[i+2]), fabsf(x[i + 3]));
		maxVal = max(maxVal, a);
		maxVal = max(maxVal, b);
	}

	for (; i < N; ++i) {
		maxVal = max(maxVal,fabsf(x[i]));
	}
	
	return maxVal;;
#undef UNROLL_NUMBER 

}

void f32_deseasonalize_inplace(F32PTR y, int N, int P,  F32PTR mean_tmp, I32PTR NumGoodPtsPerTime_tmp) {

	// mean_buf and NUmgoodPtsPertime_Buf have at least a lenght of P
	f32_compute_seasonal_avg(y, N, P, mean_tmp,  NumGoodPtsPerTime_tmp);
	int p = 0;
	for (int i = 0; i < N; i++) {
		y[i] -= mean_tmp[p];
		p++;
		p = (p == P) ? 0 : p;
	}

}

I64  sub2ind(int *dims, int ndim, int * subs) {
	
	if (ndim == 1) {
		return subs[0];
	} else if (ndim == 2) {
		I32 rows = dims[0];		
		return subs[0] + (subs[1] - 1) * rows;
	} else if (ndim == 3) {
		I64 rows = dims[0];
		I64 cols = dims[1];
		return subs[0] + (subs[1] - 1) * rows + (subs[2] - 1) * rows*cols;
	}	else {
		I64 idx    = subs[0];
		I64 stride = 1;
		for (int i = 1; i < ndim; i++) {
			stride = stride *dims[i - 1];
			idx = idx + subs[i] * stride;		
		}		
		return  idx;
	}	
}
void ind2sub(int *dims, int ndim, I64 ind, int * subs) {
	
	if (ndim == 1) {
		subs[0] = ind;
	} else if (ndim == 2) {
		I32 ROW = dims[0];				
		int c   = (ind - 1) / ROW;
		int r   = ind - c * ROW;
		c = c + 1;
		subs[0] = r;
		subs[1] = c;
	} else if (ndim == 3) {
		I64 ROW = dims[0];
		I64 COL = dims[0];
		I64 ind0 = ind - 1; // zero-based index
		int d3 = ind0 / (ROW * COL);
		ind0 = ind0 - d3 * (ROW * COL);
		int d2 = ind0 / ROW;
		int d1 = ind0 - d2 * ROW;
		subs[0] = d1+1;
		subs[1] = d2+1;
		subs[1] = d2+1;
		
	}	else {	 
		I64 ind0 = ind - 1; // zero-based index		

		subs[0] = 1; // subs used here as a temp buf
		for (int i = 1; i < ndim; i++) {
			// cacluate the stride of each dim
			subs[i] = dims[i-1] * subs[i - 1];
		}
		for (int i = ndim-1; i >0; i--) {
			// cacluate the stride of each dim
			I32 dimStride = subs[i];
			subs[i]       = ind0 / dimStride;          // now subs[i] is the index for the ith dim
			ind0          = ind0 - subs[i] * dimStride;
			subs[i]++;
		}
		subs[0] = ind0 + 1;
	}	

}

 
int ndarray_get1d_stride_offset(int *dims, int ndim, int *subs, int whichdim, I64 * stride,  I64 *offset) {
	// the dimensio requested is indicated by 0 or -1 in subs

	int whichdim0 = whichdim - 1;

	I64 STRIDE;
	I64 OFFSET         = 0;
	I64 CulDimStride   = 1;
	for (int i = 0; i < ndim; i++) {
		if (i == whichdim0) {STRIDE = CulDimStride;	} 
		OFFSET       += CulDimStride * (subs[i] - 1);
		CulDimStride *= dims[i];		
	}		
	OFFSET -= STRIDE * (subs[whichdim0] - 1);

	*stride = STRIDE;
	*offset = OFFSET;
	return dims[whichdim0];
}

void f32_get1d_from_ndarray(F32PTR dst, VOID_PTR src, int* dims, int ndim, int* subs,int whichdim, DATA_TYPE srcDtype) {

	I64 stride, offset;
	int N= ndarray_get1d_stride_offset(dims, ndim, subs, whichdim, &stride, &offset);
	 
	f32_from_strided_mem(dst, src,N, stride, offset, srcDtype);

}

void f32_set1d_to_ndarray(F32PTR src, VOID_PTR dst, int* dims, int ndim, int* subs, int whichdim, DATA_TYPE dstDtype) {

	I64 stride, offset;
	int N = ndarray_get1d_stride_offset(dims, ndim, subs, whichdim, &stride, &offset);
	f32_to_strided_mem(src, dst, N, stride, offset, dstDtype);
}

void f32_get2d_from_ndarray(F32PTR dst, VOID_PTR src, int* dims, int ndim, int* subs, int d1, int d2, DATA_TYPE srcDtype) {

	d1--; 
	d2--;
	if (d1 > d2) {
		int tmp = d1;
		d1 = d2;
		d2 = tmp;
	}

	I64 stride1, stride2;
	I64 cumstride = 1;
	I64 offset    = 0;
	for (int i = 0; i < ndim; i++) {
		if (i == d1) { stride1 = cumstride; };
		if (i == d2) { stride2 = cumstride; };
		offset    += cumstride * (subs[i]-1);
		cumstride *= dims[i];
	}
	offset -= stride1 * (subs[d1]-1);
	offset -= stride2 * (subs[d2]-1);

	int n1 = dims[d1];
	int n2 = dims[d2];
	if (d2 - d1 == 1) {		
		f32_from_strided_mem(dst, src, n1 * n2, stride1, offset, srcDtype);
	} else {
		for (int i = 0; i < n2; i++) {
			f32_from_strided_mem(dst, src, n1, stride1, offset, srcDtype);
			dst    += n1;
			offset += stride2;
		}
	}
}

void f32_set2d_from_ndarray(F32PTR src, VOID_PTR dst, int* dims, int ndim, int* subs, int d1, int d2, DATA_TYPE dstDtype) {

	d1--; 
	d2--;
	if (d1 > d2) {
		int tmp = d1;
		d1 = d2;
		d2 = tmp;
	}

	I64 stride1, stride2;
	I64 cumstride = 1;
	I64 offset    = 0;
	for (int i = 0; i < ndim; i++) {
		if (i == d1) { stride1 = cumstride; };
		if (i == d2) { stride2 = cumstride; };
		offset    += cumstride * (subs[i]-1);
		cumstride *= dims[i];
	}
	offset -= stride1 * (subs[d1]-1);
	offset -= stride2 * (subs[d2]-1);

	int n1 = dims[d1];
	int n2 = dims[d2];
	if (d2 - d1 == 1) {		
		f32_to_strided_mem(src,dst, n1 * n2, stride1, offset, dstDtype);
	} else {
		for (int i = 0; i < n2; i++) {
			f32_to_strided_mem(src, dst, n1, stride1, offset, dstDtype);
			src    += n1;
			offset += stride2;
		}
	}
}

#include "abc_000_warning.h"