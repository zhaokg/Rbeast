#pragma once 

#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void  f32_QuickSortA(F32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern void  f32_QuickSortD(F32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern  void f64_QuickSortD(F64PTR  arr, I32PTR INDEX, I32 low, I32 high);
extern void  f64_QuickSortA(F64PTR  arr, I32PTR INDEX, I32 low, I32 high);
extern void  i32_QuickSortA(I32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern void  i32_QuickSortD(I32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern  void i32_InsertionSort(I32PTR arr, I32PTR index, int n);
extern  void VOIDPTR_InsertionSort(void* arr[], char* index, int n);


extern void i32_sort_d_iterative(I32PTR  arr, int* idx, int* stack, int l, int h);
extern int  i32_find_majority_fast(I32PTR arr, int n, int* status);
extern int i32_find_unique_occurrance_inplace(I32PTR arr, int n, I32PTR counts);
#ifdef __cplusplus
}
#endif