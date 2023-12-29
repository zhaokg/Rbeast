#include "abc_datatype.h"
#include "string.h"
#include "stdlib.h"
#include "abc_sort.h"

void f32_transpose_inplace(F32PTR Mat, I32 ROWs, I32 COLs) { 
  //w: number of rows, h: number of columns
  //https: //rosettacode.org/wiki/Matrix_transposition#C
	
	I32 totalElement = ROWs * COLs;
	for (I32 start = 0; start < totalElement; start++) {
		I32 next = start;
		I32 i    = 0;
		do {
			++i;
			next = (next % COLs) * ROWs + next / COLs;
		} while (next > start);
		if (next < start || i == 1) continue; // s is a fixed point or s has been previously rorated

		F32 tmp = Mat[next = start];
		do {
			i = (next % COLs) * ROWs + next / COLs;
			Mat[next] = (i == start) ? tmp : Mat[i];
			next = i;
		} while (next > start);
	}
}


// ALGORITHM 513 : Analysis of In - Situ Transposition I - FI - I ;ESKO G.CATE and DAVID W.TWIGG
static int greatest_common_divsor(int A, int B) {	
	int remainder;
	do {
		remainder = A % B;
		A = B;
		B = remainder;
	}  while (remainder != 0);
	return  A;
}

#define NEXT(s)  ((s*NCOL) % K)
#define PREV(s)  ((s*NROW) % K)
#define Nwork    1000

void i32_transpose_inplace_next(I32PTR Mat, I32 NROW, I32 NCOL) {

	if (NROW == 1 || NCOL == 1) { return; }

	I32 NfixedPoints      = greatest_common_divsor(NROW - 1, NCOL - 1) + 1L;
	I32 Nprocessed        = NfixedPoints;
	I08 WORKED[Nwork + 1] = { 0, };
	I32 K                 = NROW * NCOL - 1;
	for (I32 start = 1; start <= K - 1; start++) {
		if (start <= Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                 break;

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			int next         = NEXT(start);			
			if (next == start) continue; // No need to increment Nprocessed

			int  tmpCurrValue = Mat[start];
			while (next != start) {
				int tmpNextValue = Mat[next];
				Mat[next]    = tmpCurrValue;				
				tmpCurrValue = tmpNextValue;
				WORKED[(next <= Nwork)*next] = 1; // avoid the if branch by using the multiplication trick
				Nprocessed++;
				next         = NEXT(next);
			}
			Nprocessed++;      // next must be equal to start
			Mat[next] = tmpCurrValue;		
		}	else {
			// Do not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork

			I32 prev = NEXT(start);
			if (prev <= start) {
				continue;  // it is not the smallest index and must have been processed bofre
				// or if it ia fixed point
			}

			I32 next = start; // start is not a fixed point because it has already been checked above
			do {
				next = NEXT(next);
			} while (next > start); 
			if (next < start) {
				continue;  // it is not the smallest index and must have been processed bofre
			}
			
			     next         = NEXT(start);
			int  tmpCurrValue = Mat[start];
			while (next != start) {
				int tmpNextValue = Mat[next];
				Mat[next] = tmpCurrValue;
				tmpCurrValue = tmpNextValue;	
				Nprocessed++;
				next = NEXT(next);
			}
			Nprocessed++;            // next must be equal to start
			Mat[next] = tmpCurrValue;

			
		}
		

	}
}

void i32_transpose_inplace_prev(I32PTR Mat, I32 NROW, I32 NCOL) {
	//https: //rosettacode.org/wiki/Matrix_transposition#C
	if (NROW == 1 || NCOL == 1) { 	return; }

	I32 NfixedPoints    = greatest_common_divsor(NROW - 1, NCOL - 1) + 1L;
	I32 Nprocessed      = NfixedPoints;
	I08 WORKED[Nwork+1] = { 0, };
	I32 K               = NROW * NCOL - 1;
	for (I32 start = 1; start <= K-1; start++) {
		if (start<=Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                 break;

		int  curr = start;
		int  prev = PREV(start);
		if (prev <= start) continue; // it s a fixed point or has already been processed

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			int  tmpStartValue = Mat[start];			
			while (prev != start) {				
				WORKED[(prev <= Nwork)*prev] = 1;
				Mat[curr] = Mat[prev];				
				curr      = prev;
				prev      = PREV(curr);								
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;
		}	else {
			// Do know not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork
 
			I32 next = start; // start is not a fixed point because it has already been checked above
			do {
				next = NEXT(next);
			} while (next > start); 
			if (next < start) {
				continue;  // it is not the smallest index and must have been processed bofre
			}

			int  tmpStartValue = Mat[start];			
			while (prev != start) {				
				Mat[curr] = Mat[prev];				
				curr  = prev;
				prev  = PREV(curr);								
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;			
		}
		

	}
}

void i32_transpose_inplace_prev_two_ends(I32PTR Mat, U64 NROW, U64 NCOL) {
	// rosettacode.org/wiki/Matrix_transposition#C
	if (NROW == 1 || NCOL == 1) { return; }

	I32 NfixedPoints = greatest_common_divsor(NROW - 1, NCOL - 1) + 1L;
	I32 Nprocessed   = NfixedPoints;
	I08 WORKED[Nwork + 1] = { 0, };
	U64 K            = NROW * NCOL - 1;        // use int64 t- avoid ovwerflow in PREV()
	for (U64 start = 1; start <= K - 1; start++) {
		if (start <= Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                     break;

		U64  curr = start;
		U64  prev = PREV(start);
		if (prev <= start) continue; // it s a fixed point or has already been processed

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			U64  startCompanion = K-start;
			U64  prevCompanion  = K - prev;
			U64  currCompanion  = startCompanion;

			WORKED[(currCompanion <= Nwork) * currCompanion] = 1;

			int  tmpStartValue          = Mat[start];
			int  tmpStartValueCompanion = Mat[startCompanion];
			while (prev != start) {
				if (prev == startCompanion) {
					break;
				}

				WORKED[(prev <= Nwork) * prev] = 1;
				Mat[curr] = Mat[prev];
				Nprocessed++;

				WORKED[(prevCompanion <= Nwork) * prevCompanion] = 1;
				Mat[currCompanion] = Mat[prevCompanion];
				Nprocessed++;

				curr = prev;
				prev = PREV(curr);

				currCompanion = prevCompanion;
				prevCompanion = PREV(currCompanion);
			
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;

			Nprocessed++;            // next must be equal to start
			Mat[currCompanion] = tmpStartValueCompanion;
		}	else {
			// Do know not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork

			U64 next = start; // start is not a fixed point because it has already been checked above
			do {
				next = NEXT(next);
			} while (next > start);
			if (next < start) {
				continue;  // it is not the smallest index and must have been processed bofre
			}

			int  tmpStartValue = Mat[start];
			while (prev != start) {
				Mat[curr] = Mat[prev];
				curr = prev;
				prev = PREV(curr);
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;
		}


	}
}
 
void i32_permute3d_inplace_abc123_acb132(I32PTR Mat, U64  NROW, U64 NCOL, int NZ) {
   // THis is essentilly the 2D tranpose, except each element is the Nrow-elem vector.
   // Must allocate a temp buf to hold two such elmeents: one for the cycle from one end,
  // another for the cycle from the other end

	if (NROW == 1) {
		// this essentially is a 2D matrix of dimenion [NCOL, NZ]
		i32_transpose_inplace_prev_two_ends(Mat, NCOL, NZ);
		return;
	}

	if (NROW == 2) {
		// this essentially is a 2D matrix of dimenion [NCOL, NZ], each elment being 64-bit
		//i64_transpose_inplace_prev_two_ends(Mat, NCOL, NZ);
		return;
	}

	NROW = NCOL;
	NCOL = NZ;
	if (NROW == 1 || NCOL == 1) { return; }

	int    ElemSize               = NROW;
	I32PTR tempbuf                = malloc(ElemSize * sizeof(I32) * 2);
	I32PTR tmpStartValue          = tempbuf;
	I32PTR tmpStartValueCompanion = tempbuf + ElemSize;

	I32 NfixedPoints = greatest_common_divsor(NROW - 1, NCOL - 1) + 1L;
	I32 Nprocessed   = NfixedPoints;
	I08 WORKED[Nwork + 1] = { 0, };
	U64 K             = NROW * NCOL - 1;

	for (U64 start = 1; start <= K - 1; start++) {

		if (start <= Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                     break;

		U64  curr = start;
		U64  prev = PREV(start);
		if (prev <= start) continue; // it s a fixed point or has already been processed

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			U64  startCompanion = K - start;
			U64  prevCompanion  = K - prev;
			U64  currCompanion  = startCompanion;

			WORKED[(currCompanion <= Nwork) * currCompanion] = 1;

			memcpy(tmpStartValue,          Mat +start* ElemSize,ElemSize * sizeof(I32));
			memcpy(tmpStartValueCompanion, Mat + startCompanion * ElemSize, ElemSize * sizeof(I32));

			while (prev != start) {
				if (prev == startCompanion) {
					break;
				}

				WORKED[(prev <= Nwork) * prev] = 1;
				memcpy(Mat + curr * ElemSize, Mat + prev * ElemSize, ElemSize * sizeof(I32));
				Nprocessed++;

				WORKED[(prevCompanion <= Nwork) * prevCompanion] = 1;
				memcpy(Mat + currCompanion + curr * ElemSize, Mat + prevCompanion * ElemSize, ElemSize * sizeof(I32));
				Nprocessed++;

				curr = prev;
				prev = PREV(curr);

				currCompanion = prevCompanion;
				prevCompanion = PREV(currCompanion);

			}
			Nprocessed++;            // next must be equal to start
			memcpy(Mat + curr * ElemSize, tmpStartValue, ElemSize * sizeof(I32));

			Nprocessed++;            // next must be equal to start
			memcpy(Mat + currCompanion * ElemSize, tmpStartValueCompanion, ElemSize * sizeof(I32));
		}
		else {
			// Do know not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork

			U64 next = start; // start is not a fixed point because it has already been checked above
			do {
				next = NEXT(next);
			} while (next > start);
			if (next < start) {
				continue;  // it is not the smallest index and must have been processed bofre
			}

		 
			memcpy(tmpStartValue, Mat + start * ElemSize, ElemSize * sizeof(I32));
			while (prev != start) {
				//Mat[curr] = Mat[prev];
				memcpy(Mat + curr * ElemSize, Mat + prev * ElemSize, ElemSize * sizeof(I32));
				curr = prev;
				prev = PREV(curr);
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			memcpy(Mat + curr * ElemSize, tmpStartValue, ElemSize * sizeof(I32));
		}


	}

	free(tempbuf);
}


void i32_permute3d_inplace_abc123_bac213(I32PTR Mat, I32 NROW, I32 NCOL, int NZ) {
   // THis is essentilly the 2D tranpose, except each element is the Nrow-elem vector.
   // Must allocate a temp buf to hold two such elmeents: one for the cycle from one end,
  // another for the cycle from the other end
 
	for (int i = 0; i < NZ; i++) {
		i32_transpose_inplace_prev_two_ends(Mat + NROW * NCOL * i, NROW, NCOL);
	}
}

void i32_permute3d_inplace_abc123_cab312(I32PTR Mat, I32 NROW, I32 NCOL, int NZ) {
	// THis is essentilly the 2D tranpose, except each element is the Nrow-elem vector.
	// Must allocate a temp buf to hold two such elmeents: one for the cycle from one end,
   // another for the cycle from the other end
	i32_transpose_inplace_prev_two_ends(Mat, NROW*NCOL, NZ);	 
}

void i32_permute3d_inplace_abc123_bca231(I32PTR Mat, I32 NROW, I32 NCOL, int NZ) {
	// THis is essentilly the 2D tranpose, except each element is the Nrow-elem vector.
	// Must allocate a temp buf to hold two such elmeents: one for the cycle from one end,
   // another for the cycle from the other end
	i32_transpose_inplace_prev_two_ends(Mat, NROW , NCOL*NZ);
}

void i32_permute3d_inplace_abc123_cba321_prev(I32PTR Mat, I32 NROW, I32 NCOL, int NZ) {

	//rosettacode.org/wiki/Matrix_transposition#C
	if ( NCOL == 1) { 
		i32_transpose_inplace_prev_two_ends(Mat, NROW, NZ);
		return;
	}
	 
	I32 Nprocessed      = 2 ;// the first and last elements ar fixed points
	I08 WORKED[Nwork+1] = { 0, };
	I32 K               = NROW * NCOL*NZ - 1;

	I32 AB = NROW * NCOL;
	I32 CB = NZ * NCOL;

	for (I32 start = 1; start <= K-1; start++) {

		int x, y, z;  // Temp varialbes used to calculate indices
#ifdef PREV
	#undef PREV
#endif
#ifdef NEXT
	#undef NEXT
#endif
#define PREV(J)  (x=J/CB, y=J/NZ  -NCOL*x,  z=J % NZ,   x+NROW*y+AB*z)
#define NEXT(I)  (z=I/AB, y=I/NROW-NCOL*z, x=I % NROW, z+NZ*y  +CB*x)

		if (start<=Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                 break;

		int  curr = start;
		int  prev = PREV(start);
		if (prev <= start) {
			Nprocessed += (prev == start); // it is a fixed point and we need to increase the counter
			// but we won't increase the counter f prev < start (i.e., the elme has already been processed
			continue; 
		}

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			int  tmpStartValue = Mat[start];			
			while (prev != start) {				
				WORKED[(prev <= Nwork)*prev] = 1;
				Mat[curr] = Mat[prev];				
				curr      = prev;
				prev      = PREV(curr);								
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;
		}	else {
			// Do know not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork

			// This is not a FIXED point because it has been checked above using "prev == start"
 
			I32 next = start; // start is not a fixed point because it has already been checked above
			do {
				next = NEXT(next);
			} while (next > start); 
			if (next < start) {
				// but we won't increaset the counter
				continue;  // it is not the smallest index and must have been processed bofre
			}

			int  tmpStartValue = Mat[start];			
			while (prev != start) {				
				Mat[curr] = Mat[prev];				
				curr  = prev;
				prev  = PREV(curr);								
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;			
		} //	if (start <= Nwork)
		

	}
 
}

void i32_permuate3d(I32PTR mat, int N1, int N2, int N3, int* new_order) {
	int order = new_order[0] * 100 + new_order[1] * 10 + new_order[2];
	if (order == 132 ) {
		i32_permute3d_inplace_abc123_acb132(mat, N1, N2, N3);
	}
	else if (order == 321) {
		  i32_permute3d_inplace_abc123_cba321_prev(mat, N1, N2, N3);
	}
	else if (order == 213) {
		i32_permute3d_inplace_abc123_bac213(mat, N1, N2, N3);
	}
	else if (order == 312) {
		i32_permute3d_inplace_abc123_cab312(mat, N1, N2, N3);
	}
	else if (order == 231) {
		i32_permute3d_inplace_abc123_bca231(mat, N1, N2, N3);
	}
	return;
#undef PREV
#undef NEXT
}

void i32_permute_any_dim(I32PTR Mat, I32 *N, I32 * order, int ndims) {

	I32 Nprocessed        = 2;// the first and last elements ar fixed points
	I08 WORKED[Nwork + 1] = { 0, };
	I32 K  =1;	
	I32 AB[20];
	I32 CB[20];
	I32 prod1=1, prod2 = 1;
	for (int i = 0; i < ndims; i++) {
		K = K * N[i];
		AB[i] = prod1;
		CB[i] = prod2;
		prod1 *=N[i];
		prod2 *=N[order[i]-1];

		order[i]--;  // change it to zero-based indices
	}
	K = K - 1;

	int orderreverse[20];
	for (int i = 0; i < ndims; i++) {
		orderreverse[order[i]] = i;
	}

	for (I32 start = 1; start <= K - 1; start++) {

 
#define NEXT(next, I)  do {  \
			int INDICES[20];  \
			int offset = 0;     \
			for (int i = ndims - 1; i > 0; i--) {  \
				INDICES[i] = I / AB[i] - offset;    \
                offset     =  N[i-1]*(INDICES[i]+offset); \
			}		                                \
		    INDICES[0]= I%N[0];   \
			next = INDICES[order[0] ]; \
			for (int i = 1; i < ndims; i++) { \
				next += INDICES[order[i]] * CB[i]; \
			} \
		} while(0)

#define PREV(prev, J)  do {  \
			int INDICES[20];  \
			int offset = 0;     \
			for (int i = ndims - 1; i > 0; i--) {  \
				INDICES[i] = J / CB[i] - offset;    \
                offset     =  N[order[i-1]]*(INDICES[i]+offset); \
		     }		                      		\
		    INDICES[0]= J%N[order[0]];   \
			prev      = INDICES[orderreverse[0] ]; \
			for (int i = 1; i < ndims; i++) { \
				prev += INDICES[orderreverse[i]] * AB[i]; \
			} \
		} while(0)

		if (start <= Nwork && WORKED[start]) 	continue; 		// has been already rotated		
		if (Nprocessed > K)                 break;

		int  curr = start;
		int  prev; 	PREV(prev,start);
		if (prev <= start) {
			Nprocessed += (prev == start); // it is a fixed point and we need to increase the counter
			// but we won't increase the counter f prev < start (i.e., the elme has already been processed
			continue;
		}

		if (start <= Nwork) {
			// has never been rotrated and start should aslo be the smallest index because all the previous s have been procesed
			// WORKED[start] = 1; // Not really needed bcz this element will never be visited 
			int  tmpStartValue = Mat[start];
			while (prev != start) {
				WORKED[(prev <= Nwork) * prev] = 1;
				Mat[curr] = Mat[prev];
				curr = prev;
				PREV(prev,curr);
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;
		}
		else {
			// Do know not whether or not s has been rotated
			// No need to mark WORKED because it has already pass the tracked region of size Nwork

			// This is not a FIXED point because it has been checked above using "prev == start"

			I32 next = start; // start is not a fixed point because it has already been checked above
			do {
				NEXT(next,next);
			} while (next > start);
			if (next < start) {
				// but we won't increaset the counter
				continue;  // it is not the smallest index and must have been processed bofre
			}

			int  tmpStartValue = Mat[start];
			while (prev != start) {
				Mat[curr] = Mat[prev];
				curr = prev;
				PREV(prev,curr);
				Nprocessed++;
			}
			Nprocessed++;            // next must be equal to start
			Mat[curr] = tmpStartValue;
		} //	if (start <= Nwork)


	}

	for (int i = 0; i < ndims; i++) {
		order[i]++;  // change it back to the one-based indices
	}

}

void i32_permuate_nd(I32PTR mat, int *dims, int* order, int ndim) {

	if (ndim > 20 || ndim < 2) {
		return;
	}

	if (ndim == 2) {
		i32_transpose_inplace_prev_two_ends(mat, dims[0], dims[1]);
		return;
	}
	
	if(ndim == 3) {
		i32_permuate3d(mat, dims[0], dims[1], dims[2], order);
		return;
	}

	int  new_dims[20];
	int  new_groupstart[20];
	int  new_ndim = 0;
	// Re-group the order into consective groups
	int i = 0;
	int culdims_group = 1;
	while (i < ndim-1) {		
		new_groupstart[new_ndim] = order[i];				
		culdims_group = 1;
		while ((order[i + 1] - order[i]) == 1) {
			culdims_group *= order[i];		
			i++;
		}
		culdims_group         *= dims[i];
		new_dims[new_ndim]    = culdims_group;

		new_ndim++;
		i++;
	}
	if (order[ndim - 1] == order[ndim - 2]) {
		// the last dim belongs to the prior grou[
		new_dims[new_ndim - 1] = culdims_group * dims[ndim - 1];
	}	else {
		// the last dim does not belong to the prior grou[
		new_dims[new_ndim]       = dims[ndim - 1];
		new_groupstart[new_ndim] = order[ndim - 1];
		new_ndim++;
	}

	int new_index[20];
	for (int i = 0; i < new_ndim; i++) {
		new_index[i] = i;
	}
		
	i32_QuickSortA(new_groupstart, new_index, 0, new_ndim - 1);
	int new_order[20];
	int new_orgdims[20];
	for (int i = 0; i < new_ndim; i++) {
		new_order[new_index[i]]   = i + 1;
		new_orgdims[i]            = new_dims [new_index[i]];
	}

	// TODO:remove singletons
	for (int i = 0; i < new_ndim; i++) {
 
	}

	if (new_ndim == 2) {
		i32_transpose_inplace_prev_two_ends(mat, new_orgdims[0], new_orgdims[1]);
		return;
	}

	if (new_ndim == 3) {
		i32_permuate3d(mat, new_orgdims[0], new_orgdims[1], new_orgdims[2], new_order);
		return;
	}
	i32_permute_any_dim(mat, new_orgdims, new_order, new_ndim);
	return;
}