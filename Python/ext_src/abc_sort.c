#include "abc_000_warning.h"

#include "abc_sort.h"


/*
* Introsort is used here. Alternative methods like timsort and powersort are also fast, but they 
* are based on mergesort and not in-place.
*/

//https://stackoverflow.com/questions/24991208/expand-a-macro-in-a-macro

#define JOIN1(X,Y)              X ## _ ## Y
#define FUNC_JOIN(fun, name)    JOIN1(fun,name)
#define FUNC(name)              FUNC_JOIN(funprefix,name) 

// FUNC(123)
//  -> FUNC_JOIN(funprefix, 123)
//  -> JOIN1(f32a, 123)
//  ->  f32a_123 

#define JOIN2(X, Y)             X ## _ ## Y ## _index
#define FINDEX_JOIN(fun, name)  JOIN2(fun, name) 
#define FINDEX(name)            FINDEX_JOIN(funprefix, name)

#define SwapIndex( index, i, j)             { int   tmp =index[i];  index[i]=index[j]; index[j]=tmp;}
#define SwapElements( arr, i, j)            { DTYPE tmp =arr[i];    arr[i]  =arr[j];   arr[j]=tmp;}
#define SwapElemIndex( arr, index, i, j)    { SwapElements(arr,i,j); SwapIndex(index,i,j); }

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         float
#define funprefix     f32a
#define LESS(a, b)    (a < b)
#define LESSEQ(a, b)  (a <= b)
#include "abc_sort_template.h"

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         float
#define funprefix     f32d
#define LESS(a, b)    (a > b)
#define LESSEQ(a, b)  (a >= b)
#include "abc_sort_template.h"


//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         double
#define funprefix     f64a
#define LESS(a, b)    (a < b)
#define LESSEQ(a, b)  (a <= b)
#include "abc_sort_template.h"

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         double
#define funprefix     f64d
#define LESS(a, b)    (a > b)
#define LESSEQ(a, b)  (a >= b)
#include "abc_sort_template.h"


//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int16_t
#define funprefix     i16a
#define LESS(a, b)    (a < b)
#define LESSEQ(a, b)  (a <= b)
#include "abc_sort_template.h"

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int16_t
#define funprefix     i16d
#define LESS(a, b)    (a > b)
#define LESSEQ(a, b)  (a >= b)
#include "abc_sort_template.h"


//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int32_t
#define funprefix     i32a
#define LESS(a, b)    (a < b)
#define LESSEQ(a, b)  (a <= b)
#include "abc_sort_template.h"

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int32_t
#define funprefix     i32d
#define LESS(a, b)    (a > b)
#define LESSEQ(a, b)  (a >= b)
#include "abc_sort_template.h"


//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int64_t
#define funprefix     i64a
#define LESS(a, b)    (a < b)
#define LESSEQ(a, b)  (a <= b)
#include "abc_sort_template.h"

//////////////////////////////////////////////////////////////////////////////////////////////
#define DTYPE         int64_t
#define funprefix     i64d
#define LESS(a, b)    (a > b)
#define LESSEQ(a, b)  (a >= b)
#include "abc_sort_template.h"

int i32_find_unique_occurrance_inplace(I32PTR arr, int n, I32PTR counts) {

	I32PTR stack = counts; // Use counts a temp buf for stack in the sorting
		
	i32a_introSort(arr, 0, n-1);

	int Nunique = 0;
	for (int i = 0; i < n; ) {
		int xcur   = arr[i];
		int ncount = 0;
		while (i < n && arr[i] == xcur) {
			ncount++;
			i++;
		}
		arr[Nunique] = xcur;
		counts[Nunique] = ncount;
		Nunique++;
	}

 	i32a_introSort_index(counts, 0, Nunique - 1, arr);
 
	return Nunique;

}


int i32_find_majority_fast( I32PTR arr, int n, int * status) {
// www.geeksforgeeks.org/boyer-moore-majority-voting-algorithm
// This is the Boyer-Moore Majority Voting Algorithm, guarantted to find the majority element with occurrences > N/2 times
// But might fail if the occurrence is not. So use status to indicate failure or success
		
	int  candidate;
	int  votes = 0;
	// Finding majority candidate
	for ( int i = 0; i < n; i++) {
		if (votes == 0) {
			candidate = arr[i];
			votes     = 1;
		} else {
			votes += (arr[i] == candidate) ? 1 : -1;
		}
	}

	// Checking if  themajority candidate occurs more than n/2 times
	int count = 0;	
	for (int i = 0; i < n; i++) {
		count += (arr[i] == candidate);			
	}

	status[0] = (count > n / 2) ? 1L : 0L;		 
	return candidate;
}





/*********************************************************
*          The quick sort algorithm
**********************************************************/

/* This function takes last element as pivot, places the pivot element at its correct
position in sorted array, and places all smaller (smaller than pivot) to left of pivot 
and all greater elements to right of pivot 

 THis default quicksort has serious problemes: The complexity for orderred arrays are n^2, taking too long and alos
 making the swap for all the elments
*/


static INLINE void SwapValueI32(I32PTR a, I32PTR b) { I32 t = *a; *a = *b;	*b = t; }
static INLINE void SwapIntIndex(I32PTR a, I32PTR b) { I32 t = *a;	*a = *b;	*b = t; }


/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low  --> Starting INDEX,
high  --> Ending INDEX */

static I32 i32_PartitionA(I32PTR arr,   I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] <= pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueI32(&arr[i], &arr[j]);
			//SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueI32(&arr[i + 1], &arr[high]);
	//SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void i32_QuickSortA(I32PTR arr, I32 low, I32 high) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now
		at right place */
		I32 pi = i32_PartitionA(arr, low, high);

		// Separately sort elements before
		// partition and after partition
		i32_QuickSortA(arr, low, pi - 1);
		i32_QuickSortA(arr, pi + 1, high);
	}
}
 
  
  
// www.geeksforgeeks.org/iterative-quick-sort/
 
/* A[]    --> Array to be sorted, l --> Starting index,		h --> Ending index */
static void i32_quicksortA_iterative(I32PTR  arr, int *stack, int l, int h) {

 	// Stack should be at least of length [h-l+1]

	// Create an auxiliary stack
    // int stack[h - l + 1];
     
    int top      = -1;     // initialize top of stack    
    stack[++top] = l;      // push initial values of l and h to stack
    stack[++top] = h;
 
    // Keep popping from stack while is not empty
    while (top >= 0) {
        // Pop h and l
        h = stack[top--];
        l = stack[top--];
 
        // Set pivot element at its correct position  in sorted array
        int p = i32_PartitionA(arr, l, h);
 
        // If there are elements on left side of pivot,
        // then push left side to stack
        if (p - 1 > l) {
            stack[++top] = l;
            stack[++top] = p - 1;
        }
 
        // If there are elements on right side of pivot,
        // then push right side to stack
        if (p + 1 < h) {
            stack[++top] = p + 1;
            stack[++top] = h;
        }
    }
	return;
}



//////////////////////////////////////////////////////////////////////////////////////////////

static INLINE void SwapValueF64(F64PTR a, F64PTR b) { F64 t = *a; *a = *b;	*b = t; }



/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low  --> Starting INDEX,
high  --> Ending INDEX */

static I32 f64_PartitionD(F64PTR arr,   I32 low, I32 high, I32PTR INDEX) {
	F64 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] > pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueF64(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF64(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void f64_QuickSortD(F64PTR arr, I32 low, I32 high, I32PTR INDEX) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now
		at right place */
		I32 pi = f64_PartitionD(arr, low, high, INDEX);

		// Separately sort elements before
		// partition and after partition
		f64_QuickSortD(arr, low, pi - 1, INDEX);
		f64_QuickSortD(arr, pi + 1, high, INDEX);
	}
}
 
  

static INLINE void SwapValueI64(I64PTR a, I64PTR b) { 	I64 t = *a; 	*a = *b;	*b = t; }

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low  --> Starting INDEX,
high  --> Ending INDEX */

static I32 i64_PartitionD(I64PTR arr,   I32 low, I32 high, I32PTR INDEX) {
	I64 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] > pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueI64(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueI64(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void i64_QuickSortD(I64PTR arr, I32 low, I32 high, I32PTR INDEX) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now
		at right place */
		I32 pi = i64_PartitionD(arr, low, high, INDEX);

		// Separately sort elements before
		// partition and after partition
		i64_QuickSortD(arr, low, pi - 1, INDEX);
		i64_QuickSortD(arr, pi + 1, high, INDEX);
	}
}

/*
#include "time.h"
#include "stdlib.h"  //srand


void test_sort_alg( void ) {

	double a[256], b[256];
	int    indices[256];

	srand(time(NULL));

	for (int k = 0; k < 91256; k++) {

		for (int i = 0; i < 256; i++) {
			indices[i] = i;
			double tmp = (double)rand() / RAND_MAX;
			double tmp2 = (double)rand() / RAND_MAX;
			if (tmp2 > 0.2) {
				tmp = 0;
			}
			a[i] = b[i] = tmp;
		}

		f64d_introSort_index(a, 0, 255, indices);

		for (int j = 0; j < 256 - 1; j++) {
			if (a[j] < a[j + 1]) {
				printf("Error\n");
			}

			if (a[j] != b[indices[j]]) {
				printf("Error b\n");
			}
		}
	}
}

*/

 
#include "abc_000_warning.h"

