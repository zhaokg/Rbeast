#pragma once
#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif
		 
	extern void f32_transpose_inplace(F32PTR Mat, I32 ROW, I32 COL);
	extern void i32_transpose_inplace(I32PTR Mat, I32 NROW, I32 NCOL);
	void i32_transpose_inplace_prev(I32PTR Mat, I32 NROW, I32 NCOL);          // this is a beter versio than i32_tranpose_inplace_next
	void i32_transpose_inplace_prev_two_ends(I32PTR Mat, U64 NROW, U64 NCOL); // this is the fatested version
	void i32_permuate_nd(I32PTR mat, int* dims, int* order, int ndim);

#ifdef __cplusplus
}
#endif