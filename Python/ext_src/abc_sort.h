#include "abc_datatype.h"

#ifdef __cplusplus
extern "C" {
#endif


extern void QuickSortD(F32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern void QuickSortA(F32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern void f64_QuickSortA(F64PTR  arr, I32PTR INDEX, I32 low, I32 high);
extern void  i32_QuickSortA(I32PTR arr, I32PTR INDEX, I32 low, I32 high);
extern  void i32_InsertionSort(I32PTR arr, I32PTR index, int n);
extern  void VOIDPTR_InsertionSort(void* arr[], char* index, int n);


#ifdef __cplusplus
}
#endif