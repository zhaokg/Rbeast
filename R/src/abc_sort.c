#include "abc_000_warning.h"

#include "abc_sort.h"


/*********************************************************
*          The quick sort algorithm
**********************************************************/

// A utility function to swap two elements
static INLINE void SwapValueI32(I32PTR a, I32PTR b) {   I32 t = *a; *a = *b;	*b = t; }
static INLINE void SwapValueF32(F32PTR a, F32PTR b) {	F32 t = *a; *a = *b;	*b = t;}
static INLINE void SwapValueF64(F64PTR a, F64PTR b) {   F64 t = *a; *a = *b;	*b = t; }
static INLINE void SwapIntIndex(I32PTR a, I32PTR b) {	I32 t = *a;	*a = *b;	*b = t;}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */

static I32 f32_PartitionD(F32PTR arr, I32PTR INDEX, I32 low, I32 high) {

	F32 pivot = arr[high];   // pivot
	I32 i     = (low - 1);   // Index of smaller element

	for (rI32 j = low; j <= high - 1; j++) 	{
		// If current element is smaller than or  equal to pivot		
		if (arr[j] > pivot)  	{ //if (arr[j] <= pivot)
			i++;                   // increment INDEX of smaller element
			SwapValueF32(&arr[i],   &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF32(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low  --> Starting INDEX,
high  --> Ending INDEX */
void f32_QuickSortD(F32PTR arr, I32PTR INDEX, I32 low, I32 high)
{
	if (low < high)  {
		/* pi is partitioning INDEX, arr[p] is now at right place */
		I32 pi = f32_PartitionD(arr, INDEX, low, high);

		// Separately sort elements before  partition and after partition
		f32_QuickSortD(arr, INDEX, low, pi - 1);
		f32_QuickSortD(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////

static I32 f32_PartitionA(F32PTR arr, I32PTR INDEX, I32 low, I32 high)
{
	F32 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] <= pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueF32(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF32(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void f32_QuickSortA(F32PTR arr, I32PTR INDEX, I32 low, I32 high)
{
	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now 	at right place */
		I32 pi = f32_PartitionA(arr, INDEX, low, high);

		// Separately sort elements before partition and after partition
		f32_QuickSortA(arr, INDEX, low, pi - 1);
		f32_QuickSortA(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////

static I32 i32_PartitionA(I32PTR arr, I32PTR INDEX, I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] <= pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueF32(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF32(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void i32_QuickSortA(I32PTR arr, I32PTR INDEX, I32 low, I32 high) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now
		at right place */
		I32 pi = i32_PartitionA(arr, INDEX, low, high);

		// Separately sort elements before
		// partition and after partition
		i32_QuickSortA(arr, INDEX, low, pi - 1);
		i32_QuickSortA(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////
static I32 i32_PartitionD_v1(I32PTR arr, I32PTR INDEX, I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	I32 i     = (low - 1);  // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {	
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] > pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueF32(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF32(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

static I32 i32_PartitionD (I32PTR arr, I32PTR INDEX, I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	I32 i     = low;  // Index of the first element immedidalely after smaller element

	//for (I32 j = low; j <= high - 1; j++) {
	for (I32 j = low; j < high; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] > pivot) {
			SwapValueF32(&arr[i],   &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
			i++;    // increment INDEX of smaller element
		}
	}
	SwapValueF32(&arr[i], &arr[high]);
	SwapIntIndex(&INDEX[i], &INDEX[high]);
	return (i);
}

static I32 i32_PartitionD_nullidx (I32PTR arr, I32PTR INDEX, I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	I32 i     = low;  // Index of the first element immedidalely after smaller element

	//for (I32 j = low; j <= high - 1; j++) {
	for (I32 j = low; j < high; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)
		if (arr[j] > pivot) {
			SwapValueI32(&arr[i],   &arr[j]);
			//SwapIntIndex(&INDEX[i], &INDEX[j]);
			i++;    // increment INDEX of smaller element
		}
	}
	SwapValueI32(&arr[i], &arr[high]);
	//SwapIntIndex(&INDEX[i], &INDEX[high]);
	return (i);
}
static I32 i32_PartitionD_v3    (I32PTR arr, I32PTR INDEX, I32 low, I32 high) {
	I32 pivot = arr[high];    // pivot
	//I32 i = (low - 1);  // Index of smaller element
	I32 i     = low;  // Index of the first element immedidalely after smaller element


	//for (I32 j = low; j <= high - 1; j++) {
	for (I32 j = low; j < high; j++) {
		// If current element is smaller than or equal to pivot
		//if (arr[j] <= pivot)

		const int shouldMove = arr[j] > pivot;		
		int newValueI = shouldMove ? arr[j] : arr[i];
		int newValueJ = shouldMove ? arr[i] : arr[j];
		int newIndexI = shouldMove ? INDEX[j] : INDEX[i];
		int newIndexJ = shouldMove ? INDEX[i] : INDEX[j];
		arr[i]   = newValueI;
		arr[j]   = newValueJ;
		INDEX[i] = newIndexI;
		INDEX[j] = newIndexJ; 
		i = i+shouldMove;    // increment INDEX of smaller element
	}


	SwapValueF32(&arr[i], &arr[high]);
	SwapIntIndex(&INDEX[i], &INDEX[high]);
	return (i);
}

 
void i32_QuickSortD(I32PTR arr, I32PTR INDEX, I32 low, I32 high) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now
		at right place */
		I32 pi = i32_PartitionD(arr, INDEX, low, high);

		// Separately sort elements before
		// partition and after partition
		i32_QuickSortD(arr, INDEX, low, pi - 1);
		i32_QuickSortD(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////

static I32 f64_PartitionA(F64PTR arr, I32PTR INDEX, I32 low, I32 high) {
	F64 pivot = arr[high];    // pivot
	I32 i     = (low - 1);    // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot if (arr[j] <= pivot)
		if (arr[j] <= pivot) {
			i++;    // increment INDEX of smaller element
			SwapValueF64(&arr[i], &arr[j]);
			SwapIntIndex(&INDEX[i], &INDEX[j]);
		}
	}
	SwapValueF64(&arr[i + 1], &arr[high]);
	SwapIntIndex(&INDEX[i + 1], &INDEX[high]);
	return (i + 1);
}

void f64_QuickSortA( F64PTR  arr, I32PTR INDEX, I32 low, I32 high) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now	at right place */
		I32 pi = f64_PartitionA(arr, INDEX, low, high);

		// Separately sort elements before partition and after partition
		f64_QuickSortA(arr, INDEX, low,    pi - 1);
		f64_QuickSortA(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////

static I32 f64_PartitionD(F64PTR arr, I32PTR INDEX, I32 low, I32 high) {
	F64 pivot = arr[high];    // pivot
	I32 i     = (low - 1);    // Index of smaller element

	for (I32 j = low; j <= high - 1; j++) {
		// If current element is smaller than or equal to pivot if (arr[j] <= pivot)
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

void f64_QuickSortD( F64PTR  arr, I32PTR INDEX, I32 low, I32 high) {

	if (low < high)	{
		/* pi is partitioning INDEX, arr[p] is now	at right place */
		I32 pi = f64_PartitionD(arr, INDEX, low, high);

		// Separately sort elements before partition and after partition
		f64_QuickSortA(arr, INDEX, low,    pi - 1);
		f64_QuickSortA(arr, INDEX, pi + 1, high);
	}
}

/////////////////////////////////////////////////////////////////////////
 void VOIDPTR_InsertionSort(void* arr[], char* index, int n) {
	int i, j;
	for (i = 1; i < n; i++) {
		void* key = arr[i];
		char  idx = index[i];
		j = i - 1;
		/* Move elements of arr[0..i-1], that are
		greater than key, to one position ahead
		of their current position */
		while (j >= 0 && arr[j] > key) {
			arr[j + 1]   = arr[j];
			index[j + 1] = index[j];
			j = j - 1;
		}
		arr[j + 1]   = key;
		index[j + 1] = idx;
		//index[j + 1] = index[j];
	}
}

void i32_InsertionSort(I32PTR arr, I32PTR index, int n) {
	int i, j;
	for (i = 1; i < n; i++) {
		I32    key = arr[i];
		I32    idx = index[i];
		j = i - 1;
		/* Move elements of arr[0..i-1], that are
		greater than key, to one position ahead
		of their current position */
		while (j >= 0 && arr[j] > key) {
			arr[j + 1]   = arr[j];
			index[j + 1] = index[j];
			j = j - 1;
		}
		arr[j + 1]   = key;
		index[j + 1] = idx;
		//index[j + 1] = index[j];
	}
}


// www.geeksforgeeks.org/iterative-quick-sort/
 
/* A[]    --> Array to be sorted, l --> Starting index,		h --> Ending index */
void i32_sort_d_iterative(I32PTR  arr, int* idx, int *stack, int l, int h) {


	I32  (*i32_PartitionD_fptr)(I32PTR arr, I32PTR INDEX, I32 low, I32 high) = i32_PartitionD;
	if (!idx) {
		i32_PartitionD_fptr = i32_PartitionD_nullidx;
	}

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
        int p = i32_PartitionD_fptr(arr, idx, l, h);
 
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

int i32_find_unique_occurrance_inplace(I32PTR arr, int n, I32PTR counts) {

	I32PTR stack = counts; // Use counts a temp buf for stack in the sorting
	I32PTR idx   = NULL;
	i32_sort_d_iterative(arr, idx, stack, 0, n-1);

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

	// CHeck if there is engough space left in arr as a temp stack
	if (Nunique < n / 2) {
		I32PTR stack = arr + Nunique;
		i32_sort_d_iterative(counts, arr, stack, 0, Nunique - 1);
	} else {
		i32_QuickSortD(counts, arr, 0, Nunique - 1);
	}

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


#include "abc_000_warning.h"

