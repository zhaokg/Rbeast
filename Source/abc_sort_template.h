
static INLINE void  FUNC(sort2)(DTYPE * arr, int i, int j) {
	if (LESS(arr[j], arr[i])) { SwapElements(arr, i, j); }
}

static INLINE void  FINDEX(sort2)(DTYPE* arr, int* index, int a, int b) {
	if (LESS(arr[b], arr[a])) { SwapElemIndex(arr, index,a, b); }
}

static INLINE void  FUNC(sort3)(DTYPE * arr, int a, int  b, int c) {
		
	FUNC(sort2)(arr, a, b);
	if ( LESS(arr[c], arr[b]) ) {
		SwapElements(arr, b, c);
		FUNC(sort2)(arr, a, b);
	}

}

static INLINE void  FINDEX(sort3)(DTYPE* arr, int* index, int a, int  b, int c ) {

	FINDEX(sort2)(arr, index, a, b);
	if (LESS(arr[c], arr[b])) {
		SwapElemIndex(arr, index, b, c);
		FINDEX(sort2)(arr, index, a, b);
	}

}

static int FUNC(binary_search)(DTYPE * arr, DTYPE key, int low, int high) {
	 
	while (high > low) {
		int mid = (low + high) / 2;
		if ( LESS(key,  arr[mid]) )   high = mid - 1;
		else       		              low  = mid + 1;
	}

	if (high < low) {// this happens only if item < arr[low]
		return low;
	} else {
		return LESS(key, arr[low]) ? low : low + 1;
	}
}


// www.geeksforgeeks.org/binary-insertion-sort/

static void FUNC(insertionSort_standard)(DTYPE* arr, int low, int high) {

	for (int i = low + 1; i <= high; ++i) {
		DTYPE   key = arr[i];
		int     cursor = i - 1;
		/* Move elements of arr[0..i-1], that are greater than key,
		 to one position ahead of their current position */
		while (cursor >= low && LESS(key, arr[cursor])) {
			arr[cursor + 1] = arr[cursor];
			--cursor;
		}
		arr[cursor + 1] = key;
		//index[j + 1] = index[j];
	}
}

static  void FUNC(insertionSort_binary_byOne)(DTYPE* arr, int L, int H) {

	for (int i = L + 1; i <= H; i++) {
		DTYPE key = arr[i];
		int   newloc = FUNC(binary_search)(arr, key, L, i - 1);
		for (int j = i - 1; j >= newloc; j--) {
			arr[j + 1] = arr[j];
		}
		arr[newloc] = key;
	}

}

static  void FUNC(insertionSort_binary_byTwo)(DTYPE* arr, int low, int high) {
	int N = high - low + 1;

	if (N == 1) return;

	int skip = 1;
	if (N % 2 == 0) {
		FUNC(sort2)(arr, low, low + 1);
		skip = 2;
	}

	for (int i = low + skip; i <= (high - 1); i += 2) {

		FUNC(sort2)(arr, i, i + 1);

		DTYPE key1 = arr[i];
		DTYPE key2 = arr[i + 1];

		int newloc1 = FUNC(binary_search)(arr, key1, low, i - 1);

		if (newloc1 == i) {
			continue;
		}

		int newloc2 = FUNC(binary_search)(arr, key2, newloc1, i - 1);

		for (int j = i - 1; j >= newloc2; j--) { arr[j + 2] = arr[j]; }
		arr[newloc2 + 1] = key2;

		for (int j = newloc2 - 1; j >= newloc1; j--) { arr[j + 1] = arr[j]; }
		arr[newloc1] = key1;
	}

}


static  void FUNC(insertionSort_binary_byThree)(DTYPE* arr, int low, int high) {

	int N = high - low + 1;

	if (N == 1) return;

	FUNC(sort2)(arr, low, low + 1);

	int skip = N % 3;
    if (skip == 0){
		if (LESS(arr[low+2], arr[low+1])) {
			SwapElements(arr, low+1, low+2);
			FUNC(sort2)(arr, low, low+1);
		}
		skip = 3;
	}

	for (int i = low + skip; i <= (high - 2); i += 3) {

		FUNC(sort3)(arr, i, i + 1, i+2);

		DTYPE key1 = arr[i];
		DTYPE key2 = arr[i + 1];
		DTYPE key3 = arr[i + 2];

		int newloc1 = FUNC(binary_search)(arr, key1, low, i - 1);

		if (newloc1 == i) {
			continue;
		}

		int newloc2 = FUNC(binary_search)(arr, key2, newloc1, i - 1);
		int newloc3;
		if (newloc2 == i) {
			newloc3 = i;
		}  else {
			newloc3 = FUNC(binary_search)(arr, key3, newloc2, i - 1);
		}

		for (int j = i - 1; j >= newloc3; j--) { arr[j + 3] = arr[j]; }
		arr[newloc3 + 2] = key3;

		for (int j = newloc3 - 1; j >= newloc2; j--) { arr[j + 2] = arr[j]; }
		arr[newloc2 + 1] = key2;

		for (int j = newloc2 - 1; j >= newloc1; j--) { arr[j + 1] = arr[j]; }
		arr[newloc1] = key1;

	}

}

static  void FINDEX(insertionSort_binary_byThree)(DTYPE* arr, int low, int high, int *index) {

	int N = high - low + 1;

	if (N == 1) return;

	FINDEX(sort2)(arr, index,low, low + 1);

	int skip = N % 3;
	if (skip == 0) {
		if (LESS(arr[low + 2], arr[low + 1])) {
			SwapElemIndex(arr, index, low + 1, low + 2);
			FINDEX(sort2)(arr, index,  low,     low + 1);
		}
		skip = 3;
	}

	for (int i = low + skip; i <= (high - 2); i += 3) {

		FINDEX(sort3)(arr, index, i, i + 1, i + 2);

		DTYPE key1 = arr[i];
		DTYPE key2 = arr[i + 1];
		DTYPE key3 = arr[i + 2];
	
		int newloc1 = FUNC(binary_search)(arr, key1, low, i - 1);

		if (newloc1 == i) {
			continue;
		}

		int newloc2 = FUNC(binary_search)(arr, key2, newloc1, i - 1);
		int newloc3;
		if (newloc2 == i) {
			newloc3 = i;
		} else {
			newloc3 = FUNC(binary_search)(arr, key3, newloc2, i - 1);
		}

		int index1 = index[i];
		int index2 = index[i+1];
		int index3 = index[i+2];

		for (int j = i - 1; j >= newloc3; j--) { arr[j + 3] = arr[j]; index[j + 3] = index[j];		}
		arr[  newloc3 + 2] = key3;
		index[newloc3 + 2] = index3;

		for (int j = newloc3 - 1; j >= newloc2; j--) { arr[j + 2] = arr[j]; index[j + 2] = index[j];	}
		arr[newloc2 + 1]   = key2;
		index[newloc2 + 1] = index2;

		for (int j = newloc2 - 1; j >= newloc1; j--) { arr[j + 1] = arr[j];  index[j +1] = index[j];	}
		arr[newloc1]       = key1;
		index[newloc1] = index1;

	}

}

// attern-defeating Quicksort: https://arxiv.org/pdf/2106.05123.pdf
static int FUNC(partition_right)(DTYPE* arr, int low, int high, int * no_swap) {

	int  N = high - low + 1;

	// if Ttere are only two elments
	if (N == 2) {
		FUNC(sort2)(arr, low, high);
		return low;
	}

	// Find the mdedian of the three as the pviot
	int mid   = (low + high) / 2;
	int first = low;
	int last;
	if (N > 512) {
		FUNC(sort3)(arr, low,    mid-1,   high-2);
		FUNC(sort3)(arr, low+1,  mid,     high-1);
		FUNC(sort3)(arr, low+2,  mid+1,   high);
		FUNC(sort3)(arr, mid-1,   mid,  mid+1);
		SwapElements(arr, low, mid);
		last = high + 1;   // the last elment may be smaller than the pivot
	} else {		
		FUNC(sort3)(arr, mid, low, high);
		last = high;   // the last elment is known to be larger than the pivot
	}
	

	DTYPE pivot  = arr[low]; // the leftmost element is the privot

	while ( LESS( arr[++first], pivot) );
	if (first == low +1 && !LESS(arr[mid], pivot) && !LESS(arr[low+2], pivot) ) {
		while (last > first && !LESS(arr[--last], pivot) );
	}  else {
		while (!LESS(arr[--last], pivot));
	}

	no_swap[0] = first >= last;

	while (first < last) {
		SwapElements(arr, first, last);
		while ( LESS(arr[++first], pivot) );
		while ( !LESS(arr[--last], pivot) );
	}

	int pivot_finalidx = first - 1;
	SwapElements(arr, low, pivot_finalidx);

	return pivot_finalidx;
}

static int FINDEX(partition_right)(DTYPE* arr, int* index, int low, int high,  int* no_swap) {

	int  N = high - low + 1;

	// if Ttere are only two elments
	if (N == 2) {
		FINDEX(sort2)(arr, index, low, high);
		return low;
	}

	// Find the mdedian of the three as the pviot
	int mid   = (low + high) / 2;
	int first = low;
	int last;
	if (N > 512) {
		FINDEX(sort3)(arr, index, low,    mid-1,   high-2);
		FINDEX(sort3)(arr, index, low+1,  mid,     high-1);
		FINDEX(sort3)(arr, index, low+2,  mid+1,   high);
		FINDEX(sort3)(arr, index, mid-1,   mid,    mid+1);
		SwapElemIndex(arr, index, low, mid);
		last = high + 1;   // the last elment may be smaller than the pivot
	} else {		
		FINDEX(sort3)(arr, index, mid, low, high);
		last = high;   // the last elment is known to be larger than the pivot
	}
	

	DTYPE pivot  = arr[low]; // the leftmost element is the privot

	while ( LESS( arr[++first], pivot) );
	if (first == low +1 && !LESS(arr[mid], pivot) && !LESS(arr[low+2], pivot) ) {
		while (last > first && !LESS(arr[--last], pivot) );
	}  else {
		while (!LESS(arr[--last], pivot));
	}

	no_swap[0] = first >= last;

	while (first < last) {
		SwapElemIndex(arr, index, first, last);
		while ( LESS(arr[++first], pivot) );
		while ( !LESS(arr[--last], pivot) );
	}

	int pivot_finalidx = first - 1;
	SwapElemIndex(arr, index, low, pivot_finalidx);

	return pivot_finalidx;
}

static void FUNC(quickSort)(DTYPE arr[], int low, int high) {

	// when low is less than high
	if (low < high)	{
		// pivotidx is the partition return index of pivot
		int no_swap;
		int pivotidx = FUNC(partition_right)(arr, low, high, &no_swap);

		//Recursion Call:smaller element than pivot goes left and higher element goes right
		FUNC(quickSort)(arr, low,          pivotidx - 1);
		FUNC(quickSort)(arr, pivotidx + 1, high);
	}
}

static void FINDEX(quickSort)(DTYPE arr[], int low, int high, int *index) {

	// when low is less than high
	if (low < high) {
		// pivotidx is the partition return index of pivot
		int no_swap;
		int pivotidx = FINDEX(partition_right)(arr, index, low, high, &no_swap);

		//Recursion Call:smaller element than pivot goes left and higher element goes right
		FINDEX(quickSort)(arr, low, pivotidx - 1 , index);
		FINDEX(quickSort)(arr, pivotidx + 1, high, index);
	}
}

/////////////////////////////////////////////////////////////////////////////////////
// www.geeksforgeeks.org/heap-sort/ 
// To heapify a subtree rooted with node i which is an index in arr[]. n is size of heap

static  void FUNC(heapify)(DTYPE arr[], int N, int i) {

	// Find largest among root, left child and right child
	// Initialize largest as root
	int largest = i;
	int left    = 2 * i + 1;
	int right   = 2 * i + 2;

	// If left child is larger than root
	if (left < N && LESS( arr[largest], arr[left]) )   { largest = left;	}

	// If right child is larger than largest so far
	if (right < N && LESS(arr[largest], arr[right]))   {largest = right;	}
	
	// Swap and continue heapifying::  if root is not largest;
	if (largest != i) {		
		SwapElements(arr, i, largest);
		// Recursively heapify the affected: sub-tree
		FUNC(heapify)(arr, N, largest);
	}
}

// Main function to do heap sort
static void FUNC(heapSort)(DTYPE arr[], int low, int high) {

	DTYPE* arrstart  = arr + low;
	int    N         = high - low + 1;

	// Build max heap
	for (int i = N / 2 - 1; i >= 0; i--) {
		FUNC(heapify)(arrstart, N, i);
	}
	
	// Heap sort
	for (int i = N - 1; i >= 0; i--) {
		SwapElements(arrstart, 0, i);
		// Heapify root element to get highest element atm root again
		FUNC(heapify)(arrstart, i, 0);
	}
}

static void FINDEX(heapify)(DTYPE arr[], int N, int i, int* index) {

	// Find largest among root, left child and right child
	// Initialize largest as root
	int largest = i;
	int left = 2 * i + 1;
	int right = 2 * i + 2;

	// If left child is larger than root
	if (left < N && LESS(arr[largest], arr[left])) { largest = left; }

	// If right child is larger than largest so far
	if (right < N && LESS(arr[largest], arr[right])) { largest = right; }

	// Swap and continue heapifying::  if root is not largest;
	if (largest != i) {
		SwapElemIndex(arr, index, i, largest);
		// Recursively heapify the affected: sub-tree
		FINDEX(heapify)(arr, N, largest, index);
	}
}

static void FINDEX(heapSort)(DTYPE arr[], int low, int high, int *index) {

	DTYPE* arrstart   = arr   + low;
	int  * indexstart = index + low;
	int  N            = high - low + 1;

	// Build max heap
	for (int i = N / 2 - 1; i >= 0; i--) {
		FINDEX(heapify)(arrstart, N, i, indexstart);
	}
	
	// Heap sort
	for (int i = N - 1; i >= 0; i--) {
		SwapElemIndex(arrstart, indexstart, 0, i);
		// Heapify root element to get highest element atm root again
		FINDEX(heapify)(arrstart, i, 0, indexstart);
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

static int FUNC(isorderedarray)(DTYPE arr[], int  low, int  high) {

	int ordered = 1;
	for (int i = low + 1; i <= high; i++) {
		if ( LESS(arr[i],  arr[i - 1]) )  {
			ordered = 0;
			break;
		}
	}

	return ordered; 
}

// A Utility function to perform intro sort
static void FUNC(introsortUtil)(DTYPE arr[], int  low, int  high, int depthLimit) {

	// Count the number of elements
	int size = high - low +1;

	if (size <=1) {
		return;
	}

	// If partition size is low then do insertion sort
	if (size <= 32) {
		FUNC(insertionSort_binary_byThree)(arr, low, high);
		return;
	}

	// If the depth is zero use heapsort
	if ( depthLimit == 0) {
		FUNC(heapSort)(arr, low, high);
		return;
	}

	// Else use a median-of-three concept to
	int no_swap;
	int pivotidx = FUNC(partition_right)(arr, low, high, &no_swap);

	int leftOrdered = 0, rightOrdered = 0;
	if (no_swap) {
		leftOrdered  = FUNC(isorderedarray)(arr, low, pivotidx - 1);
		rightOrdered = FUNC(isorderedarray)(arr, pivotidx + 1, high);
	}


	// Check for a highly unbalanced partition.
	int  l_size = pivotidx - low;
	int  r_size = high     - (pivotidx + 1);
	int  highly_unbalanced = l_size < size / 8 || r_size < size / 8;


	//Recursion Call:smaller element than pivot goes left and higher element goes right
	if (!leftOrdered)   FUNC(introsortUtil)(arr, low,           pivotidx - 1 , depthLimit- highly_unbalanced);
	if (!rightOrdered) 	FUNC(introsortUtil)(arr, pivotidx + 1, high,           depthLimit - highly_unbalanced);
	return;
}



static void FINDEX(introsortUtil)(DTYPE arr[], int  low, int  high, int * index, int depthLimit) {

	// Count the number of elements
	int size = high - low +1;

	if (size <=1) {
		return;
	}

	// If partition size is low then do insertion sort
	if (size <= 32 ){
		FINDEX(insertionSort_binary_byThree)(arr, low, high,index);
		return;
	}

	//ERROR: Fixed error with the following function
	//FINDEX(heapSort)(arr, low, high, index); 

	// If the depth is zero use heapsort
	if ( depthLimit == 0) {
		FINDEX(heapSort)(arr, low, high,index);
		return;
	}

	// Else use a median-of-three concept to
	int no_swap;
	int pivotidx = FINDEX(partition_right)(arr, index, low, high,  &no_swap);

	int leftOrdered = 0, rightOrdered = 0;
	if (no_swap) {
		leftOrdered  = FUNC(isorderedarray)(arr, low, pivotidx - 1);
		rightOrdered = FUNC(isorderedarray)(arr, pivotidx + 1, high);
	}


	// Check for a highly unbalanced partition.
	int  l_size = pivotidx - low;
	int  r_size = high     - (pivotidx + 1);
	int  highly_unbalanced = l_size < size / 8 || r_size < size / 8;


	//Recursion Call:smaller element than pivot goes left and higher element goes right
	if (!leftOrdered)   FINDEX(introsortUtil)(arr, low,           pivotidx - 1 , index, depthLimit- highly_unbalanced);
	if (!rightOrdered) 	FINDEX(introsortUtil)(arr, pivotidx + 1, high,           index, depthLimit - highly_unbalanced);
	return;
}


/* Implementation of introsort*/
void FUNC(introSort)(DTYPE arr[], int low, int high) {

	int N          = high - low + 1;
	int depthLimit = 0;
	while (N >>= 1) ++depthLimit;  // a sloppy way to compute log2(N)

	FUNC(introsortUtil)(arr, low, high, depthLimit);  // Perform a recursive Introsort
	return;
}

/* Implementation of introsort*/
void FINDEX(introSort)(DTYPE arr[], int low, int high, int *index) {

	int N          = high - low + 1;
	int depthLimit = 0;
	while (N >>= 1) ++depthLimit;  // a sloppy way to compute log2(N)

	FINDEX(introsortUtil)(arr, low, high,index, depthLimit);  // Perform a recursive Introsort
	return;
}




/*****************************************************************************/
// Implementation of introsort based on a stack without recursion
// This code is not used and kept here as a d
/*****************************************************************************/

static void FUNC(introSort_iterative)(DTYPE arr[], int low, int high) {

	int N = high - low + 1;
	int depthLimit = 0;
	while (N >>= 1) ++depthLimit;  // a sloppy way to compute log2(N)

	struct { int left, right, depthlimit; } stack[32];

	int top=-1;
	top++;
	stack[top].left       = low;
	stack[top].right      = high;
	stack[top].depthlimit = depthLimit;

	while (top >= 0) {
       	
		int left  = stack[top].left;
		int right = stack[top].right;
		int depth = stack[top].depthlimit;
		top--;

		// Count the number of elements
		int size = right - left + 1;

		if (size <= 1) {
			continue;
		}

		// If partition size is low then do insertion sort
		if (size <= 32) {
			FUNC(insertionSort_binary_byTwo)(arr, left, right);
			continue;
		}

		// If the depth is zero use heapsort
		if (depth == 0) {
			FUNC(heapSort)(arr, left, right);
			continue;
		}

		// Else use a median-of-three concept to
		int no_swap;
		int pivotidx = FUNC(partition_right)(arr, left, right, &no_swap);

		int leftOrdered = 0, rightOrdered = 0;
		if (no_swap) {
			leftOrdered  = FUNC(isorderedarray)(arr, left, pivotidx - 1);
			rightOrdered = FUNC(isorderedarray)(arr, pivotidx + 1, right);
		}

		// Check for a highly unbalanced partition.
		int  l_size = pivotidx - left;
		int  r_size = right - (pivotidx + 1);
		int  highly_unbalanced = l_size < size / 8 || r_size < size / 8;
	
		//Recursion Call:smaller element than pivot goes left and higher element goes right
		if (!leftOrdered) {
			top++;
			stack[top].left  = left;
			stack[top].right = pivotidx - 1;			
			stack[top].depthlimit  =depth - highly_unbalanced;	 
		} 

		if (!rightOrdered) {
			top++;
			stack[top].left       = pivotidx + 1;
			stack[top].right      = right;	
			stack[top].depthlimit = depth - highly_unbalanced;
		}	
	}

	 
	return;
}

#undef DTYPE
#undef funprefix
#undef LESS
#undef LESSEQ
