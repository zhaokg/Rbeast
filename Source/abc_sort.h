#pragma once 

#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif

extern  void f64a_introSort(F64PTR arr,  I32 low, I32 high);
extern  void f64d_introSort(F64PTR arr,  I32 low, I32 high);
extern  void f32a_introSort(F32PTR arr,  I32 low, I32 high);
extern  void f32d_introSort(F32PTR arr,  I32 low, I32 high);
extern  void i16a_introSort(I16PTR arr, I32 low, I32 high);
extern  void i16d_introSort(I16PTR arr, I32 low, I32 high);
extern  void i32a_introSort(I32PTR arr,  I32 low, I32 high);
extern  void i32d_introSort(I32PTR arr,  I32 low, I32 high);
extern  void i64a_introSort(I64PTR arr, I32 low, I32 high);
extern  void i64d_introSort(I64PTR arr, I32 low, I32 high);

extern  void f64a_introSort_index(F64PTR arr, I32 low, I32 high, I32PTR index);
extern  void f64d_introSort_index(F64PTR arr, I32 low, I32 high, I32PTR index);
extern  void f32a_introSort_index(F32PTR arr, I32 low, I32 high, I32PTR index);
extern  void f32d_introSort_index(F32PTR arr, I32 low, I32 high, I32PTR index);
extern  void i16a_introSort_index(I16PTR arr, I32 low, I32 high, I32PTR index);
extern  void i16d_introSort_index(I16PTR arr, I32 low, I32 high, I32PTR index);
extern  void i32a_introSort_index(I32PTR arr, I32 low, I32 high, I32PTR index);
extern  void i32d_introSort_index(I32PTR arr, I32 low, I32 high, I32PTR index);
extern  void i64a_introSort_index(I64PTR arr, I32 low, I32 high, I32PTR index);
extern  void i64d_introSort_index(I64PTR arr, I32 low, I32 high, I32PTR index);

extern int  i32_find_majority_fast(I32PTR arr, int n, int* status);
extern int  i32_find_unique_occurrance_inplace(I32PTR arr, int n, I32PTR counts);

#ifdef __cplusplus
}
#endif