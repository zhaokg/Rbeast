#include "abc_000_warning.h"
#include "abc_mat.h"


static void chol(F32PTR XtX, F32PTR U, I32 K, I32 k)
{
	if (k == 0) 	{
		U[0] = sqrtf(XtX[0]);
		U[K] = XtX[K] / U[0];
		U[K + 1] = sqrtf(XtX[K + 1] - U[K] * U[K]);
		return;
	}


	 for (rI32 m = 0; m <= 1; m++)
	{
		I32	   k_K = k*K;
		F32PTR    oldColPtr = U;
		F32PTR    newColPtr;

		memcpy(U + k_K, XtX + k_K, sizeof(F32)*(k + 1));

		F32 sum = 0.f;
		for (I32 i = 0; i < k; i++)
		{
			newColPtr = U + k_K;
			F32   tmp = 0.f;
			for (I32 j = 0; j < i; j++)
			{
				tmp += (*newColPtr++)* (*oldColPtr++);
			}
			*newColPtr = tmp = ((*newColPtr) - tmp) / (*oldColPtr);
			sum = sum + tmp*tmp;
			oldColPtr = oldColPtr + (K - i);
		}
		newColPtr++;
		*newColPtr = sqrtf(*newColPtr - sum);
		k++;
	}



}

void chol_update_U(F32PTR U, F32PTR x, I32 ldu, I32 n) {
 // ldu: the leading deimsion of U
 // n:   the lenght of x, and the number of columns to be procesed
	F32PTR Ubase = U;
	for (I32 row = 1; row <= n; row++) 	{
		U        = Ubase + (row -1)* ldu +(row -1);
		F32 Ukk  = (*U);		
		F32 r    = sqrtf(Ukk*Ukk + (*x)*(*x));
		F32 c    = r / Ukk, cinv= Ukk / r;
		F32 s    = (*x) / Ukk;
		(*U) = r;
	
		// only c, cinv, and s needed in the loop
		for (I32 col = row +1; col <= n; col++)		{
				U = U + ldu;
				x = x +  1;
				*U = (*U + s* (*x))*cinv;
				*x = c*(*x) - s*(*U);		
		}
		x = x - (n-1) + row; //first go to the start of x and then jump to the (k+1)-th element
	}
	// out of the loop, x points to one element after the end of the original x vector
}
void chol_dwdate_U(F32PTR U, F32PTR x, I32 ldu, I32 n) {
 // ldu: the leading deimsion of U
 // n:   the lenght of x, and the number of columns to be procesed
	F32PTR Ubase = U;
	for (I32 row = 1; row <= n; row++) 	{
		U        = Ubase + (row -1)* ldu +(row -1);
		F32 Ukk  = (*U);		
		F32 r    = sqrtf(Ukk*Ukk - (*x)*(*x));
		F32 c    = r / Ukk, cinv= Ukk / r;
		F32 s    = (*x) / Ukk;
		(*U) = r;
	
		// only c, cinv, and s needed in the loop
		for (I32 col = row +1; col <= n; col++)		{
				U = U + ldu;
				x = x +  1;
				*U = (*U - s* (*x))*cinv;
				*x = c*(*x) - s*(*U);		
		}
		x = x - (n-1) + row; //first go to the start of x and then jump to the (k+1)-th element
	}
	// out of the loop, x points to one element after the end of the original x vector
}
void chol_update_L(F32PTR L, F32PTR x, I32 ldu, I32 n) {
 // ldu: the leading deimsion of U
 // n:   the lenght of x, and the number of columns to be procesed
	F32PTR Lbase = L;
	for (I32 col = 1; col <= n;col++) 	{
		L        = Lbase + (col -1)* ldu +(col -1);
		F32 Lkk  = (*L);
		F32 r    = sqrtf(Lkk* Lkk + (*x)*(*x));
		F32 c    = r / Lkk, cinv= Lkk / r;
		F32 s    = (*x) / Lkk;
		(*L) = r;
	
		// only c, cinv, and s needed in the loop
		for (I32 row = col +1; row<= n; row++)		{
				L = L +  1;
				x = x +  1;
				*L = (*L + s* (*x))*cinv;
				*x = c*(*x) - s*(*L);
		}
		x = x - (n-1) + col; //first go to the start of x and then jump to the (k+1)-th element
	}
	// out of the loop, x points to one element after the end of the original x vector
}
void chol_dwdate_L(F32PTR L, F32PTR x, I32 ldu, I32 n) {
  // ldu: the leading deimsion of U
 // n:   the lenght of x, and the number of columns to be procesed
	F32PTR Lbase = L;
	for (I32 col = 1; col <= n;col++) 	{
		L        = Lbase + (col -1)* ldu +(col -1);
		F32 Lkk  = (*L);
		F32 r    = sqrtf(Lkk* Lkk - (*x)*(*x));
		F32 c    = r / Lkk, cinv= Lkk / r;
		F32 s    = (*x) / Lkk;
		(*L) = r;
	
		// only c, cinv, and s needed in the loop
		for (I32 row = col +1; row<= n; row++)		{
				L = L +  1;
				x = x +  1;
				*L = (*L - s* (*x))*cinv;
				*x = c*(*x) - s*(*L);
		}
		x = x - (n-1) + col; //first go to the start of x and then jump to the (k+1)-th element
	}
	// out of the loop, x points to one element after the end of the original x vector
}

void chol_columwise(F32PTR A, F32PTR U, I64  N, I64 K)
{
	// A is am upper triagnel matrix
	// A =U'*U:  A of size NxN, decompose from Col 1 up to Col K
	// the chol decomposition under the guise of a backward solver

	F32PTR A_base = A;
	F32PTR U_base = U;

	for (I64 COL = 1; COL <= K; COL++) 	{
		A  = A_base + (COL - 1)*N;
		U  = U_base;
		F32PTR Ucol    = U_base + (COL - 1) * N;
		F64    SUM_UxU = 0.f;
		for (I64 col = 1; col< COL; col++) 	{
			F64 sum = 0.f;
			for (I32 row = 1; row < col; row++) {
				sum += (*U++)*(*Ucol++);
			}
			*Ucol   = (A[col - 1] - sum) / (*U);
			SUM_UxU += (*Ucol)*(*Ucol);

			Ucol = Ucol - (col - 1);  //go to the start of the current col
			U    = U - (col - 1) + N; //jump to the next col
		}

		Ucol[COL - 1] = sqrt( A[COL-1] - SUM_UxU);
	}

}
void chol_columwise_v2( F32PTR A, F32PTR U, I64  N, I64 K )
{   // A is an upper triangle matrix
	// A =U'*U:  A of size NxN

	F32PTR A_base = A;
	F32PTR U_base = U;
	for (I32 COL = 1; COL <= K; ++COL) {			       // (i,j) for row and col indices, respectively				
		U  = U_base ;								       // U  points to the start of the mat
		A  = A_base + (COL - 1) * N ;				       // Au points to the start of the col-th column		
		F32PTR  Ucol             = U_base + (COL - 1) * N; // Ucol points to the start of the col-th column
		F64		SUM_Ucol_x_Ucol  = 0;				
		for (I32 col = 1; col < COL; ++col) {       // col is a col index
			F64 sum = 0.f;
			for (I32 row = 1; row < col; ++row)		//U points to the start of the j-th col
				sum += U[row - 1] * Ucol[row - 1];

			F64 res		    = (A[col -1]-sum) / U[col -1];
			Ucol[col - 1]	= res;
			SUM_Ucol_x_Ucol += res * res;	

			U += N;								 // Jump to the next col (i.e., the new j-th col)
		}
		Ucol[COL -1] = sqrt(A[COL -1]-SUM_Ucol_x_Ucol);
	}
	    

}
void chol_rowwise( F32PTR A, F32PTR U, I64  N, I64 K ) { 

	// A is am upper triagnel matrix
	// A =U'*U:  A of size NxN

	F32PTR A_base = A;
	F32PTR U_base = U;

	for (I32 ROW = 1; ROW <= K; ++ROW) 	{
		U  = U_base + (ROW-1)*N;         //go to the start of the ROW-th column
		A  = A_base + (ROW-1)*N ;        //go to the start  of the ROW-th column

		F64 sum = 0.0; 	for (I32 row = 1; row < ROW; ++row) {sum += U[row-1]* U[row-1]; }
		F32 Ukk     = sqrt(A[ROW-1]-sum);				
		F32 Ukk_inv = 1.f / Ukk;
		U[ROW - 1] = Ukk;
	
		F32PTR Ucurcol = U;			
		for (I32 col = ROW + 1; col <= K; ++col) {
			A = A + N;
			U = U + N;
			F64 sum = 0.0;	for (I32 row = 1; row < ROW; ++row)	{sum += U[row-1] * Ucurcol[row-1];}
			U[ROW-1] = ( A[ROW-1]-sum)*Ukk_inv;			
		}

	}

}
void chol_addCol(F32PTR A, F32PTR U, I64 N, I64 K0, I64 K1)
{
	// A is the K0-K1-th colums of an upper triangle matrix: A is not  a full matrix
	// A =U'*U:  A of size NxN
	F32PTR Abase = A;
	F32PTR Ubase = U;

	for (I32 COL = K0; COL <= K1; COL++) 	{		
		A			= Abase + (COL-K0)*N;		
		U			= Ubase;
		F32PTR Ucol = Ubase + (COL - 1) * N;
		F64    SUM  = 0.f;
		for (I32 col = 1; col< COL; col++)	{			
			F64 sum = 0.f;
			for (I32 row = 1; row < col; row++)	{sum += (*U++)* (*Ucol++);}			
			F64 Uk = (A[col - 1] - sum) / (*U);
			*Ucol  = Uk;
			SUM   +=Uk * Uk;

			Ucol = Ucol - (col - 1);
			U    = U - (col - 1) + N;
		}

		Ucol[COL - 1] = sqrt(A[COL - 1] - SUM);
	}

}

void inplace_chol(F32PTR A, I64  N, I64 K)
{// A =U'*U:  A of size NxN

	F32PTR Abase = A;	
	for (I64 COL = 1; COL <= K; COL++) 	{
		A = Abase + (COL-1)*N;		
		F32 Ukk_inv;
		{	F64 sum = 0.f;
			for (I64 row = 1; row < COL; row++) { sum += A[row - 1] * A[row - 1]; };			
			F64  Ukk = sqrt(A[COL-1] - sum);			
			A[COL - 1] = Ukk;
			Ukk_inv    = 1.f/Ukk;
		}		

		F32PTR U_curCol_base = A;
		A   =   A + N;
		for (I64 col2 = COL + 1; col2 <= K; col2++) {
			F32PTR U_curCol = U_curCol_base;
			F32    sum  = 0.f;
			for (rI64 row = 1; row < COL; row++) {	sum += (*A++)*(*U_curCol++);	}
			*A = (*A - sum)* Ukk_inv;
			A = A - (COL - 1) + N;
		}

	}

}
void inplace_chol_addCol(F32PTR A, I64 N, I64 K0, I64 K1)
{	// A =U'*U:  A of size NxN
	F32PTR Abase = A;
	for (I64 COL=K0; COL <= K1; COL++)
	{
		F32PTR Ucol = Abase + (COL - 1)*N;

		A = Abase;
		F32  SUM = 0.f;
		for (I64 col = 1; col< COL; col++) 	{
			F64 sum = 0.f;
			for (I64 row = 1; row < col; row++){	sum += (*A++)*(*Ucol++);}
			F64  Uk= (*Ucol - sum) / (*A);
			*Ucol = Uk;
			SUM   += Uk * Uk;

			Ucol = Ucol - (col - 1);
			A    = A - (col - 1) + N;
		}

		Ucol[COL - 1] = sqrt(Ucol[COL-1] - SUM);
	}

}
void solve_U_as_L(F32PTR A, F32PTR x, I64 lda, I64 K) {

	// A =U'*U:  A of size NxN
	for (I64 col = 1; col<= K; col++) 	{
		F64 sum = 0.f;
		for (I64 row = 1; row < col; row++)	{sum += A[row-1]*x[row-1];}
		x[col-1] =  (x[col-1] - sum) / A[col-1];			
		A = A + lda;
	}

}
void solve_U_as_U(F32PTR U, F32PTR x, I64 lda, I64 K)
{	// A =U'*U:  A of size NxN

	x = x + (K-1);

	F32PTR UlastCol_end = U + (K-1)* lda + (K-1);
	for (I64 col = 1; col <= K; col++) {
		U = UlastCol_end - (col - 1);
		F64 sum = 0.f;
		for (I64 row = 1; row < col; row++) {
			sum += (*U)*(*x--);
			U   -= lda;
		}
		*x = (*x - sum) / (*U);
		 x = x + (col - 1);
	}
	//At this point, U points back to the original start
}
void solve_L_as_L(F32PTR A, F32PTR x, I64 lda, I64 K) {

	// A =U'*U:  A of size NxN
	F32PTR Abase = A;
	for (I64 row= 1; row<= K; row++) 	{
		A = Abase+row-1;
		F64 sum = 0.f;
		for (I64 col = 1; col < row; ++col)	{
			sum += (*A)*x[col-1];
			A   +=lda;
		}
		x[row-1] =  (x[row-1] - sum) / (*A);	
	}

}
void solve_L_as_U(F32PTR A, F32PTR x, I64 lda, I64 K) {

	// A =U'*U:  A of size NxN
	x = x + (K-1);
	A = A + (K-1)* lda + (K-1); //got to the end of the last colum
	for (I64 col = K; col >= 1; col--) {		
		F64 sum = 0.f;
		for (I64 row = K; row > col; row--) {
			sum += (*A--)*(*x--);		 
		}
		*x = (*x - sum) / (*A);
		 x = x + (K-col);
		 A = A + (K - col) - lda;
	}

}
void pack_chol_update(F32PTR x, F32PTR  U, I64 K)
{
	for (I64 col = 1; col <= K; col++)
	{// Upon entry in each loop, U points to the diag element of the current col

		F32 c, cinv, s;
		{
			rF32 Ukk = *U;
			rF32 r;
			*U = r = sqrtf(Ukk*Ukk + (*x)*(*x));
			s	= (*x)/ Ukk; 	c	= r   / Ukk;   cinv = Ukk / r;
		}
		
		U = U + col;
		F32PTR U_nextCol_diag = U + 1;
		for (rI64 i = col+1; i <=K; i++)
		{			
			x++;
			*U = (*U + s*(*x))*cinv;
			*x = c*(*x) - s*(*U);

			U += i;
		}
		x = x - (K - 1) + col;


		U = U_nextCol_diag;
	}
}
void pack_chol_dwdate(F32PTR x, F32PTR U, I64 K)
{
	for (I64 col = 1; col <= K; col++)
	{// Upon entry in each loop, U points to the diag element of the current col

		F32 c, cinv, s;
		{
			rF32 Ukk = *U;
			rF32 r;
			*U = r = sqrtf(Ukk*Ukk - (*x)*(*x));
			s = (*x) / Ukk; 	c = r / Ukk;   cinv = Ukk / r;
		}

		U = U + col;
		F32PTR U_nextCol_diag = U + 1;
		for (rI64 i = col + 1; i <= K; i++)
		{
			x++;
			*U = (*U - s*(*x))*cinv;
			*x = c*(*x) - s*(*U);

			U += i;
		}
		x = x - (K - 1) + col;


		U = U_nextCol_diag;
	}
}
void pack_chol(F32PTR Au, F32PTR U, I64  N)
{
	// A =U'*U:  A of size NxN
	for (I64 COL = 1; COL <= N ; COL++)
	{
		// Upon entry in each loop:
		// U  points to the col base
		// Au points to the diag element of the col

		rF32 sum = 0.f;
		for (rI64 row = 1; row < COL; row++)
		{
			sum += U[row-1]*U[row-1];
		}

		rF32 Ukk = sqrt( *Au - sum);
		U[COL-1] = Ukk;
		Ukk = 1.f / Ukk;


		rF32PTR U_curCol_base = U;
		F32PTR Au_curCol_diagElem = Au;

		U  = U_curCol_base      + COL; //u: next col	
		Au = Au_curCol_diagElem + COL ;
		for (rI64 col2 = COL + 1; col2 <= N; col2++)
		{
			rF32PTR U_curCol = U_curCol_base;

			sum = 0.f;
			for (rI64 row = 1; row < COL; row++)
			{
				sum += (*U++)*(*U_curCol++);
			}
			*U = ((*Au) - sum)*Ukk;

			U  = U - (COL-1) + col2;
			Au = Au -(COL-1) + col2 + (COL-1);
		}

		U  = U_curCol_base + COL; //u: next col
		Au = Au_curCol_diagElem + COL + 1;
	}

}
void pack_chol_addCol(F32PTR Au, F32PTR U, I64 K0, I64 K1)
{	// A =U'*U:  A of size NxN

	F32PTR  U_base   = U;
	rF32PTR U_curCol = U + (1+(K0-1))*(K0-1)/2 ;

	for (; K0 <= K1; K0++)
	{
		// Upon entry in each loop:
		// U  points to the col base
		// Au points to the col base of Au
 
		U = U_base;

		rF32  SUM = 0.f;
		for (rI64 col = 1; col< K0; col++)
		{
			rF32 sum = 0.f;
			for (rI64 row = 1; row < col; row++)
			{
				sum += (*U++)*(*U_curCol++);
			}

			sum = ( (*Au++) - sum) / (*U++);
			*U_curCol = sum;
			SUM = SUM + sum*sum;

			U_curCol = U_curCol - (col - 1);
			//U = U - (col - 1) + N;
		}
		//U_curCol[COL - 1] = sqrt( Au[COL - 1] - SUM);
		U_curCol[K0 - 1] = sqrt((*Au++) - SUM);

		U_curCol = U_curCol + K0; //move to the next COL

	}

}

void pack_inplace_chol(F32PTR A, I64  N)
{// A =U'*U:  A of size NxN

	for (I64 COL = 1; COL <= N; COL++)
	{
		// Upon entry in each loop:
		// A  points to the col base

		rF32 sum = 0.f;
		for (rI64 row = 1; row < COL; row++)
		{
			sum += A[row - 1] * A[row - 1];
		}

		rF32 Ukk = sqrt(A[COL-1] - sum);
		A[COL-1] = Ukk;
		Ukk      = 1.f / Ukk;


		rF32PTR U_curCol_base     = A;
		
		A  = U_curCol_base      + COL; //A: next col	
		
		for (rI64 col2 = COL + 1; col2 <= N; col2++)
		{
			rF32PTR U_curCol = U_curCol_base;

			sum = 0.f;
			for (rI64 row = 1; row < COL; row++)
			{
				sum += (*A++)*(*U_curCol++);
			}
			*A = (*A - sum)*Ukk;

			A  = A - (COL - 1) + col2;			
		}

		A = U_curCol_base + COL; //u: next col
		
	}

}
void pack_inplace_chol_addCol(F32PTR A, I64 K0, I64 K1)
{	// A =U'*U:  A of size NxN

	F32PTR  U_base   = A;
	rF32PTR U_curCol = A + (1 + (K0 - 1))*(K0 - 1) / 2;

	for (; K0 <= K1; K0++)
	{
		// Upon entry in each loop:
		// A  points to the matrix base
		// U_curCol points to the col base  

		A = U_base;

		rF32  SUM = 0.f;
		for (rI64 col = 1; col< K0; col++)
		{
			rF32 sum = 0.f;
			for (rI64 row = 1; row < col; row++)
			{
				sum += (*A++) * (*U_curCol++);
			}

			sum = (*U_curCol - sum) / (*A++);
			*U_curCol = sum;
			SUM = SUM + sum*sum;

			U_curCol = U_curCol - (col - 1);
			//U = U - (col - 1) + N;
		}
		//U_curCol[COL - 1] = sqrt( Au[COL - 1] - SUM);
		U_curCol[K0 - 1] = sqrt(U_curCol[K0 - 1] - SUM);

		U_curCol = U_curCol + K0; //move to the next COL
	}

}
void pack_solve_L(F32PTR A, F32PTR x, I64 K)
{	// A =U'*U:  A of size NxN
	 
		for (rI64 col = 1; col<=K; col++)
		{
			rF32 sum = 0.f;
			for (rI64 row = 1; row < col; row++)
			{sum += (*A++)*(*x++);}
			*x = (*x - sum) / (*A++);	 

			x = x - (col - 1);		
		}	
}
void pack_solve_U(F32PTR A, F32PTR x, I64 K)
{	// A =U'*U:  A of size NxN

	rF32PTR A_lastCol_end = A + (K + 1)*K /2L - 1L;
	x = x + (K - 1);

	for (rI64 col = 1; col<= K; col++)
	{
		A = A_lastCol_end - (col - 1);

		rF32 sum = 0.f;
		for (rI64 row = 1; row < col; row++)
		{
			sum += (*A)*(*x--);
			A = A - (K - row);
		}
		*x = (*x - sum) / (*A);

		x = x + (col - 1);
	}
}

///////////////////////////
void chol_addCol_skipleadingzeros(F32PTR Au, F32PTR U, I64 N, I64 K0, I64 K1) {
	// A =U'*U:  A of size NxN
	// F32PTR A_BASE= Au;
	// Au starts at the K0-th column of the square input matrix

	F32PTR  Ubase = U;
	F32PTR  Ucol  = Ubase + (K0-1) * N;
	for (I64 COL = K0; COL <= K1; COL++){
		//Ucol = U_BASE + (COL-1) * N;
		//Au   = A_BASE + (COL-K0)*N;		
		//U    = Ubase;

		// Skip the leading zeros--if any--for the new row
		// TODO: Buggy. Some new rows may be almost zeros, so the new Cols of XtX are all zero
		// the loop below may go far over the legal boundary to overwrite other data; therefore
		// causuing crashes. Adding "&& rIdxFirstNonZero<COL" as a safeguard
		I64 rIdxFirstNonZero = 1;
		for (; Au[rIdxFirstNonZero-1] == 0 && rIdxFirstNonZero<COL;  Ucol[rIdxFirstNonZero-1]=0, rIdxFirstNonZero++ );
		U = Ubase + (rIdxFirstNonZero - 1) * N;

		F64  SUM = 0.f;
		for (I64 col = rIdxFirstNonZero; col< COL; col++) 	{			
			F64 sum = 0.f;
			for (I64 row = rIdxFirstNonZero; row < col; row++) {
				sum += U[row-1]*Ucol[row-1];
			}
			F64  Ucol_curElem =  (Au[col-1] - sum) / U[col - 1];
			Ucol[col-1]       =  Ucol_curElem;
			SUM              += Ucol_curElem * Ucol_curElem;
			U  += N;
		}

		Ucol[COL - 1] = sqrt( Au[COL-1] - SUM);

		Ucol += N;
		Au   += N;
	}

}
void chol_addCol_skipleadingzeros_prec(F32PTR Au, F32PTR U, F32 precPrior, I64 N, I64 K0, I64 K1)
{
	// A =U'*U:  A of size NxN
	// F32PTR A_BASE= Au;
	// Au starts at the K0-th column of the square input matrix

	F32PTR Ubase = U;
	F32PTR Ucol  = Ubase + (K0-1) * N;

	for (I64 COL = K0; COL <= K1; COL++){
		//Ucol = U_BASE + (COL-1) * N;
		//Au   = A_BASE + (COL-K0)*N;		
		//U    = U_BASE;

		// Skip the leading zeros--if any--for the new row
		// TODO: Buggy. Some new rows may be almost zeros, so the new Cols of XtX are all zero
		// the loop below may go far over the legal boundary to overwrite other data; therefore
		// causuing crashes. Adding "&& rIdxFirstNonZero<COL" as a safeguard
		I64 rIdxFirstNonZero = 1;
		for (; Au[rIdxFirstNonZero - 1] == 0 && rIdxFirstNonZero < COL; Ucol[rIdxFirstNonZero - 1] = 0, rIdxFirstNonZero++);
		U = Ubase + (rIdxFirstNonZero - 1) * N;

		F64  SUM = 0.f;
		for (I64 col = rIdxFirstNonZero; col< COL; col++) 	{			

			F64 sum = 0.f;for (I64 row = rIdxFirstNonZero; row < col; row++) {sum += U[row-1]*Ucol[row-1];}

			F64  Ucol_curElem = (Au[col - 1] - sum) / U[col - 1]; //invert diagoal element
			Ucol[col-1]    =  Ucol_curElem;
			SUM           += Ucol_curElem * Ucol_curElem;
			U  += N;
		}
		Ucol[COL - 1] = sqrt( (Au[COL-1]+ precPrior)- SUM); //invert diagoal element

		Ucol += N;
		Au   += N;
	}

}

void chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(F32PTR Au, F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1)
{
	// A =U'*U:  A of size NxN
	// F32PTR A_BASE= Au;
	// Au starts at the K0-th column of the square input matrix


	F32PTR  Ubase = U;
	F32PTR  Ucol = Ubase + (K0 - 1) * N;

	for (I64 COL = K0; COL <= K1; COL++) {
		//Ucol = U_BASE + (COL-1) * N;
		//Au   = A_BASE + (COL-K0)*N;		
		//U      = U_BASE;

		// Skip the leading zeros--if any--for the new row
		// TODO: Buggy. Some new rows may be almost zeros, so the new Cols of XtX are all zero
		// the loop below may go far over the legal boundary to overwrite other data; therefore
		// causuing crashes. Adding "&& rIdxFirstNonZero<COL" as a safeguard
		I64 rIdxFirstNonZero = 1;
		for (; Au[rIdxFirstNonZero - 1] == 0 && rIdxFirstNonZero < COL; Ucol[rIdxFirstNonZero - 1] = 0, rIdxFirstNonZero++);
		U = Ubase + (rIdxFirstNonZero - 1) * N;

		F64  SUM = 0.f;
		for (I64 col = rIdxFirstNonZero; col < COL; col++) {
			F64 sum = 0.f;
			for (I64 row = rIdxFirstNonZero; row < col; row++) { sum += U[row - 1] * Ucol[row - 1]; }

			F32  Ukk_invert = U[col - 1];
			F64  Ucol_curElem = (Au[col - 1] - sum) * Ukk_invert; //invert diagoal element
			Ucol[col - 1] = Ucol_curElem;
			SUM += Ucol_curElem * Ucol_curElem;
			U += N;
		}

		F32 prec       = (COL == 1) ? 0.f : *precPrior;
		Ucol[COL - 1]  = 1.f / sqrt((Au[COL - 1] + prec) - SUM); //invert diagoal element

		Ucol += N;
		Au   += N;
	}

}

void chol_addCol_skipleadingzeros_prec_invdiag(F32PTR Au, F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1)
{
	// A =U'*U:  A of size NxN
	// F32PTR A_BASE= Au;
	// Au starts at the K0-th column of the square input matrix


	F32PTR  Ubase = U;
	F32PTR  Ucol  = Ubase + (K0-1) * N;

	for (I64 COL = K0; COL <= K1; COL++){
		//Ucol = U_BASE + (COL-1) * N;
		//Au   = A_BASE + (COL-K0)*N;		
		//U      = U_BASE;

		// Skip the leading zeros--if any--for the new row
		// TODO: Buggy. Some new rows may be almost zeros, so the new Cols of XtX are all zero
		// the loop below may go far over the legal boundary to overwrite other data; therefore
		// causuing crashes. Adding "&& rIdxFirstNonZero<COL" as a safeguard
		I64 rIdxFirstNonZero = 1;
		for (; Au[rIdxFirstNonZero - 1] == 0 && rIdxFirstNonZero < COL; Ucol[rIdxFirstNonZero - 1] = 0, rIdxFirstNonZero++);
		U = Ubase + (rIdxFirstNonZero - 1) * N;

		F64  SUM = 0.f;
		for (I64 col = rIdxFirstNonZero; col< COL; col++) 	{			
			F64 sum = 0.f;
			for (I64 row = rIdxFirstNonZero; row < col; row++) {sum += U[row-1]*Ucol[row-1];}
			
			F32  Ukk_invert   = U[col - 1];
			F64  Ucol_curElem = (Au[col-1] - sum) * Ukk_invert; //invert diagoal element
			Ucol[col-1]    = Ucol_curElem;
			SUM           += Ucol_curElem * Ucol_curElem;
			U  += N;
		}		
		Ucol[COL - 1] = 1.f/sqrt( (Au[COL-1]+ *precPrior)- SUM); //invert diagoal element

		Ucol += N;
		Au   += N;
	}

}
void chol_addCol_skipleadingzeros_precVec_invdiag(   F32PTR Au, F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1)
{
	// A =U'*U:  A of size NxN
	// F32PTR A_BASE= Au;
	// Au starts at the K0-th column of the square input matrix
	// precPrior starts from the 1-str col and has a length of N or at least up to K1

	F32PTR  Ubase = U;
	F32PTR  Ucol   = Ubase + (K0-1) * N;

	for (I64 COL = K0; COL <= K1; COL++){
		//Ucol = U_BASE + (COL-1) * N;
		//Au   = A_BASE + (COL-K0)*N;		
		// U      = Ubase;

		// Skip the leading zeros--if any--for the new row
		// TODO: Buggy. Some new rows may be almost zeros, so the new Cols of XtX are all zero
		// the loop below may go far over the legal boundary to overwrite other data; therefore
		// causuing crashes. Adding "&& rIdxFirstNonZero<COL" as a safeguard
		I64 rIdxFirstNonZero = 1;
		for (; Au[rIdxFirstNonZero - 1] == 0 && rIdxFirstNonZero < COL; Ucol[rIdxFirstNonZero - 1] = 0, rIdxFirstNonZero++);
		U = Ubase + (rIdxFirstNonZero - 1) * N;

		F64  SUM = 0.f;
		for (rI64 col = rIdxFirstNonZero; col< COL; col++) 
		{			
			F64 sum = 0.f;	for (rI64 row = rIdxFirstNonZero; row < col; row++) {sum += U[row-1]*Ucol[row-1];}
						
			F64  Ukk_invert   = U[col - 1];
			F64  Ucol_curElem = (Au[col-1] - sum) * Ukk_invert; //invert diagoal element
			Ucol[col-1]      =  Ucol_curElem;
			SUM             += Ucol_curElem * Ucol_curElem;
			U  += N;
		}

		Ucol[COL - 1] = 1.f/sqrt( (Au[COL-1]+ precPrior[COL-1])- SUM); //invert diagoal element

		Ucol += N;
		Au   += N;
	}

}


void solve_U_as_LU(F32PTR U, F32PTR y, F32PTR x, I64 N, I64 K) {	

	// A =U'*U:  A of size NxN

	for (I64 col = 1; col<= K; col++) {
		F64 sum = 0.f;	
		for (I64 row = 1; row < col; row++)	{sum += U[row-1]*x[row-1];	}
		x[col-1] = (y[col-1] - sum) / U[col-1];			
		U += N;
	}


 	// A =U'*U:  A of size NxN	
	F32PTR U_lastCol_end = U-N + (K-1);
	for (I64 nCol = 1; nCol <= K; nCol++) 	{
		U = U_lastCol_end - (nCol - 1);

		F64 sum = 0.f;
		for (I64 col = K; col> (K- nCol)+1; col--) {
			sum += (*U)*x[col-1];
			U   -= N;
		}
		x[K- nCol] = (x[K- nCol]-sum) / (*U);
	}
	//At this point, U points back to the original start

}

void solve_U_as_LU_rectmat_multicols(F32PTR U, F32PTR y, F32PTR x, I64 ldu, I64 K, I64 nCols) {	
	// A =U'*U:  A of size NxN

	for (I32 I = 1; I <= nCols; ++I) {

		for (I64 col = 1; col<= K; ++col) {
			F32 sum = 0.f;
			for (I64 row = 1; row < col; ++row)	{sum += U[row-1]*x[row-1];}
			x[col-1] = ( y[col-1]-sum ) / U[col-1];	
			U += ldu;
		}

 		// A = U'*U:  A of size NxN	
		F32PTR U_LastCol_End = U- ldu + (K-1);
		for (I64 nCol = 1; nCol <= K; ++nCol)	{
			U = U_LastCol_End - (nCol-1);
			F32 sum = 0.f;
			for (I64 col = K; col > (K- nCol)+1; --col) {
				sum += (*U)*x[col-1];
				U   -= ldu;
			}
			x[K- nCol] = (x[K- nCol]-sum) /(*U);
		}	
		x += K;
		y += K;
		//At this point, U points back to the original start
	}

}
void solve_U_as_LU_invdiag_rectmat(F32PTR U, F32PTR y, F32PTR x, I64 ldu, I64 K) {	
	// A =U'*U:  A of size NxN
	for (I64 col = 1; col<= K; ++col) {
		F32 sum = 0.f;
		for (I64 row = 1; row < col; ++row)	{
			sum += U[row-1]*x[row-1];
		}
		x[col-1] = ( y[col-1]-sum ) * U[col-1];		 //invert diagoal element	
		U += ldu;
	}

 	// A = U'*U:  A of size NxN	
	F32PTR U_LastCol_End = U- ldu + (K-1);
	for (I64 nCol = 1; nCol <= K; ++nCol)	{
		U = U_LastCol_End - (nCol-1);
		F32 sum = 0.f;
		for (I64 col = K; col > (K- nCol)+1; --col) {
			sum += (*U)*x[col-1];
			U   -= ldu;
		}
		x[K- nCol] = (x[K- nCol]-sum) * (*U);  	 //invert diagoal element
	}
	//At this point, U points back to the original start
}
void solve_U_as_LU_invdiag_sqrmat(F32PTR U, F32PTR y, F32PTR x, I64 K) {	
	// Slove A*x=y
	// A =U'*U:  A of size NxN
	for (I64 col = 1; col<= K; ++col) {
		F32 sum = 0.f;
		for (I64 row = 1; row < col; ++row)	{
			sum += U[row-1]*x[row-1];
		}
		x[col-1] = (y[col-1] - sum) * U[col-1];		 //invert diagoal element	
		U += K;
	}

 	// A =U'*U:  A of size NxN	
	F32PTR U_LastCol_End = U-K + (K-1);
	for (I64 nCol = 1; nCol <= K; ++nCol) {
		U = U_LastCol_End - (nCol-1);
		F32 sum = 0.f;
		for (I64 col = K; col> (K- nCol)+1; --col) {
			sum += (*U)*x[col-1];
			U   -= K;
		}
		x[K- nCol] = (x[K- nCol]-sum) * (*U);  	 //invert diagoal element
	}
	//At this point, U points back to the original start
}
void solve_U_as_LU_invdiag_sqrmat_multicols(F32PTR U, F32PTR y, F32PTR x, I64 K, I64 nColY) {	
	// Slove A*x=y
	// A =U'*U:  A of size NxN

	for (int I = 0; I < nColY; ++I) {

		for (I64 col = 1; col <= K; ++col) {
			F64 sum = 0.f;
			for (I64 row = 1; row < col; ++row) {sum += U[row - 1] * x[row - 1];}
			x[col - 1] = (y[col - 1] - sum) * U[col - 1];		 //invert diagoal element	
			U += K;
		}

	
		F32PTR U_LastCol_End = U - K + (K - 1);
		for (I64 nCol = 1; nCol <= K; ++nCol) {
			U = U_LastCol_End - (nCol - 1);
			F64 sum = 0.f;
			for (I64 col = K; col > (K - nCol) + 1; --col) {
				sum += (*U) * x[col - 1];
				U   -= K;
			}
			x[K - nCol] = (x[K - nCol] - sum) * (*U);  	 //invert diagoal element
		}
		// at this point, U should point to the start of the Ubase
		x += K;
		y += K;
	}


}
void solve_U_as_U_invdiag(F32PTR U, F32PTR x, I64 ldu, I64 K) 
{	
 	// A =U'*U:  A of size KxX
	
	F32PTR U_lastCol_end = U + (K-1)* ldu + (K-1);
	for (I64 nCol = 1; nCol <= K; nCol++) 	{
		U = U_lastCol_end - (nCol - 1);
		F32 sum = 0.f;
		for (I64 col = K; col> (K- nCol)+1; col--) {
			sum += (*U)*x[col-1];
			U   -= ldu;
		}
		x[K- nCol] = (x[K- nCol]-sum) * (*U);  //invert diagoal element
	}
	//At this point, U points back to the original start
 
}
void solve_U_as_U_invdiag_multicols(F32PTR U, F32PTR x, I64 ldu, I64 K, I32 nColx)
{	
 	// A = U'*U:  A of size KxK
	
	for (I32 I = 0; I < nColx; ++I) {

		F32PTR U_lastCol_end = U + (K - 1) * ldu + (K - 1);
		for (I64 nCol = 1; nCol <= K; ++nCol) {
			U = U_lastCol_end - (nCol - 1);
			F64 sum = 0.f;
			for (I64 col = K; col > (K - nCol) + 1; --col) {
				sum += (*U) * x[col - 1];
				U -= ldu;
			}
			x[K - nCol] = (x[K - nCol] - sum) * (*U);  //invert diagoal element
		}
		//At this point, U points back to the original start

		x += K;
	}

}


#include "abc_blas_lapack_lib.h"


void linear_regression(F32PTR Y, F32PTR X, int ldx, int N, int K, F32PTR B,F32PTR Yfit, F32PTR Yerror, F32PTR TMP) {
	
	// Get XtY
	r_cblas_sgemv(CblasColMajor, CblasTrans, N, K, 1.f, X, ldx, Y, 1L, 0.f, B, 1L);

	F32PTR XtX = TMP;
	r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, 1.0f, X, ldx, X, ldx, 0., XtX, K);
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR, 'U',  K, XtX, K);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR, 'U',  K, 1L,  XtX, K, B, K);

	r_cblas_sgemv(CblasColMajor, CblasNoTrans, N, K, 1.f, X, ldx, B, 1L, 0.f, Yfit, 1L);

	if (Yerror)
		r_ippsSub_32f(Yfit, Y, Yerror,N); //Yerror=Y-Yfit 

}

void simple_linear_regression_nan(F32PTR Y, F32PTR X, int N, F32PTR Yfit, F32PTR Yerror) {
	F64 Xsum, Ysum, XY, XX;
	F64 b0, b1;
	int Ngood;

	Xsum = Ysum = XX = XY = Ngood = 0;
	
	if (X) {
		for (int i = 0; i < N; i++) {
			int isGoodPt = (X[i] == X[i]) && (Y[i] == Y[i]);
			Xsum += isGoodPt ? X[i] : 0;
			Ysum += isGoodPt ? Y[i] : 0;
			XY += isGoodPt ? X[i] * Y[i] : 0;
			XX += isGoodPt ? X[i] * X[i] : 0;
			Ngood += isGoodPt;
		}

	}
	else {
		//Assume X is 0, 1/N, 2/N, (N-1)/N
		for (int i = 0; i < N; i++) {
			int isGoodPt =  (Y[i] == Y[i]);
			F64 x = (double) i / (double)N;
			Xsum += isGoodPt ? x : 0;
			Ysum += isGoodPt ? Y[i] : 0;
			XY += isGoodPt ? x * Y[i] : 0;
			XX += isGoodPt ? x *x : 0;
			Ngood += isGoodPt;
		}
	}
	F64 SXY = XY - Xsum * Ysum / Ngood;
	F64 SXx = XX - Xsum * Xsum / Ngood;

	b1 = SXY / SXx;
	b0 = Ysum / Ngood - b1 * Xsum / Ngood;

	if (Yfit) {
		if (X) {
			for (int i = 0; i < N; i++) {	 
				Yfit[i] = X[i] * b1 + b0;
			}
		}else {
			for (int i = 0; i < N; i++) {
				Yfit[i] =  b1*i/N + b0;
			}
		}
	}

	if (Yerror) {
		if (X) {
			for (int i = 0; i < N; i++) {
				Yerror[i] = Y[i]-(X[i] * b1 + b0);
			}
		}
		else {
			for (int i = 0; i < N; i++) {
				Yerror[i] = Y[i] - (b1 * i / N + b0);
			}
		}
	}
}



/*
void f32_gemm_XtY1(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc) 
{
 	for (int col = 0; col < N; ++col) {
		int dir = 1;
		int row = 0;		
		for (; row < M-(2-1); row +=2) {
			//f32_dot3x1(A + ROW * lda, A + (ROW + 1) * lda, A + (ROW + 2) * lda, B, K, C + ROW);
			C[row +1] = f32_dot2x1(A + row *lda, A + (row + 1)*lda, B, K, C+ row,dir);
			//dir = !dir;
		}
		for (;  row <M; ++row) {
			C[row]=f32_dot(A + row * lda, B, K);
		}
		B += ldb;
		C += ldc;
	}
}
*/

void update_XtX_from_Xnewterm(F32PTR X, F32PTR Xnewterm, F32PTR XtX, F32PTR XtXnew, NEWCOLINFO * new ) {

	I32 k1       = new->k1;
	I32 k2_old   = new->k2_old;
	I32 k2_new   = new->k2_new;
	I32 Knewterm = new->Knewterm; // k2_new - k1 + 1L;
	I32 KOLD     = new->KOLD;
	I32 KNEW     = new->KNEW;

	I32 N    = new->N;
	I32 Nlda = new->Nlda;
	/*************************************************************************/
	//               The FIRST component:	
	/*************************************************************************/
	// There'sY no first component if k1_old/k1_new=1 for SEASON	
	for (I32 i = 1; i < k1; i++) SCPY(i, XtX + (i - 1L) * KOLD, XtXnew + (i - 1L) * KNEW);

	/*************************************************************************/
	//              The SECOND component
	/*************************************************************************/
	// No new cols/terms if flag=ChORDER && isInsert=0:the resampled basis has a higher order than the old basis
	if (Knewterm != 0) {

		FILL0(XtXnew + (k1 - 1) * KNEW, (KNEW - k1 + 1) * KNEW); // zero out the cols from k1-th to the end
		if (k1 > 1) {
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1 - 1, Knewterm, N, 1.0f,
				X, Nlda,
				Xnewterm, Nlda, 0.f,
				XtXnew + (k1 - 1L) * KNEW, KNEW);
		}


		// Three alternative ways to compute Xnewterm'*XnewTerm. Note that the resulting matrix is symmetric

		// THE FIRST WAY:
		r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
			Knewterm, Knewterm, N, 1.0,
			Xnewterm, Nlda,
			Xnewterm, Nlda, 0.f,
			XtXnew + (k1 - 1) * KNEW + k1 - 1, KNEW);


		//THE SECOND WAY: 
		//sgemmt only updates the upper triangular part of the resulting matrix, which is supposed to be faster than sgemm, but it is NOT
		//cblas_sgemmt(CblasColMajor, CblasUpper, CblasTrans, CblasNoTrans, K_newTerm, Npad, 1.0f, Xnewterm, Npad, Xnewterm, Npad, 0.f, GlobalMEMBuf_2nd, K_newTerm);

		//THE THIRD WAY: 
		//This is the fastest way when using Intel'sY MKL
		/*
		{  for (int i = 1; i <= K_newTerm; i++)
		   for (int j = 1; j <= i; j++)
		   GlobalMEMBuf_2nd[K_newTerm*(i - 1) + j - 1] = DOT(N, Xnewterm + (j - 1)*Npad, Xnewterm + (i - 1)*Npad);
		} */

		/*
		//After obtaining Xnewterm'*Xnewterm, insert it into XtX_prop at appropriate locations
		for (rI32 i = k1_new, j = 1; i <= k2_new; i++, j++) {
			if (k1_new != 1) r_cblas_scopy(k1_new - 1, MEMBUF1 + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
			r_cblas_scopy(j, MEMBUF2 + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);
		}*/
	}
	/*{
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, k1_new - 1, K_newTerm, N, 1, X_mars, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, k1_new - 1);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
	r_cblas_scopy(k1_new - 1, GlobalMEMBuf_1st + (j - 1)*(k1_new - 1), 1, XtX_prop + (i - 1)*KNEW, 1);
	cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, K_newTerm, K_newTerm, N, 1, X_mars_prop + (k1_new - 1)*N, N, X_mars_prop + (k1_new - 1)*N, N, 0, GlobalMEMBuf_1st, K_newTerm);
	for (int i = k1_new, j = 1; i <= k2_new; i++, j++)
	r_cblas_scopy(j, GlobalMEMBuf_1st + (j - 1)*K_newTerm, 1, XtX_prop + (i - 1)*KNEW + k1_new - 1, 1);
	}*/

	/*************************************************************************/
	//                  The THRID component: 
	/*************************************************************************/
	//There is no third componet if k2_old=KOLD 
	if (k2_old != KOLD) {
		/*for (rI32  j = 1; i <= KOLD; i++, j++) {r_cblas_scopy(K_newTerm,  MEMBUF1 + (j - 1)*K_newTerm, 1, XtX_prop + (k - 1)*KNEW + k1_new - 1, 1),					*/
		for (I32 kold = k2_old + 1, knew = k2_new + 1; kold <= KOLD; kold++, knew++) {
			F32PTR ColStart_old = XtX + (kold - 1) * KOLD;
			F32PTR ColStart_new = XtXnew + (knew - 1) * KNEW;
			SCPY(k1 - 1,       ColStart_old, ColStart_new); //the upper part of the third componet
			SCPY(kold - k2_old, ColStart_old + (k2_old + 1) - 1, ColStart_new + (k2_new + 1) - 1); // the bottom part of the 3rd cmpnt
		}

		// If there is a MIDDLE part of the componet (i.e, Knewterm>0); this part 
		// will be missing if flag is resmaplingOder and isInsert = 0.
		if (Knewterm != 0) {
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
				Knewterm, (KOLD - k2_old), N, 1.0,
				Xnewterm, Nlda,
				X + (k2_old + 1 - 1) * Nlda, Nlda, 0.0,
				XtXnew + (k2_new + 1 - 1) * KNEW + k1 - 1, KNEW);
		}

	}

			 
}

void update_XtY_from_Xnewterm(F32PTR Y, F32PTR Xnewterm, F32PTR XtY, F32PTR XtYnew, NEWCOLINFO* new, I32 q) {

	// X and Xnewterm has a leading dimesnion of new.Nlada
	// Y has a leading dimension of new.N

	I32 k1       = new->k1;
	I32 k2_old   = new->k2_old;
	I32 k2_new   = new->k2_new;
	I32 Knewterm = new->Knewterm;
	I32 N        = new->N;
	I32 Nlda     = new->Nlda;
	I32 KOLD     = new->KOLD;
	I32 KNEW     = new->KNEW;
/*********************************************************************************/
//                 Compute XtY_prop from XtY
/********************************************************************************/
	if (q == 1) {
		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1)       SCPY(k1 - 1, XtY,XtYnew);
		// New components : XnewTemrm*Y
		if (Knewterm > 0) { 
				r_cblas_sgemv(CblasColMajor, CblasTrans, N, Knewterm, 1.f,
						Xnewterm,   Nlda,
						Y,         1L, 0.f,
					    XtYnew + k1 - 1, 1L);
		}
		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD) SCPY(KNEW - k2_new, XtY + (k2_old + 1L) - 1L, XtYnew + (k2_new + 1) - 1);

	}
	else {
		// FOR MrBEAST

		// Skipped if k1_old=1 when dealing with SEASON
		if (k1 > 1) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(k1 - 1, XtY + KOLD * c, XtYnew + KNEW * c);
			}
		}
		// New components : XnewTemrm*Y
		if (Knewterm > 0) {
			//MatxVec(NEW.SEG, NEW.numSeg, Xnewterm, yInfo.Y, MODEL.prop.XtY + NEW.k1 - 1, N);
			r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
				Knewterm, q, N, 1.f, 
				Xnewterm, Nlda,
				Y,        N, 0.f,
				XtYnew + k1 -1, KNEW);
		}

		//this part will be skipped if k2_old=KOLD when dealing with TREND(Istrend==1)
		if (k2_old != KOLD) {
			for (I32 c = 0; c < q; ++c) {
				SCPY(KNEW - k2_new, XtY + (k2_old + 1L) - 1L + KOLD * c, XtYnew + (k2_new + 1) - 1 + KNEW * c);
			}
		}


	}
}



void get_parts_for_newinfo(NEWCOLINFOv2* new) {
	// x1, xnew1, x2, xnew2, x3, xnew3, x4

	int Knewterm    = 0;
	int Kbase_dst   = 1;

	int jParts      = 0;
	for (int i = 0; i <  new->nbands; ++i) {

		new->parts[jParts].X      = new->X;
		new->parts[jParts].ks_dst = Kbase_dst;
		
		if (i == 0) {
			new->parts[jParts].ks_src = 1;
		} else {
			new->parts[jParts].ks_src = new->ks_x[i-1] + new->kterms_x[i-1];
		}
		new->parts[jParts].kterms = new->ks_x[i]- new->parts[jParts].ks_src;

		Kbase_dst += new->parts[jParts].kterms;
		jParts++;
		///////////////////////////////////////////////////////	


		//////////////////////////////////////////////////////
		new->parts[jParts].X      = new->Xnewterm;
		new->parts[jParts].ks_dst = Kbase_dst;
		
		new->parts[jParts].ks_src = new->ks_xnewterm[i];
		new->parts[jParts].kterms = new->kterms_xnewterm[i];
		Kbase_dst    += new->parts[jParts].kterms;
		Knewterm     += new->kterms_xnewterm[i];

		jParts++;
		
	}
	// The last part is the X's remaining part
	new->parts[jParts].X       = new->X;
	new->parts[jParts].ks_dst = Kbase_dst;	
	new->parts[jParts].ks_src = new->ks_x[new->nbands - 1] + new->kterms_x[new->nbands - 1];
	new->parts[jParts].kterms = (new->K - new->parts[jParts].ks_src)+1L;
    
	new->Knew     = Kbase_dst + (new->parts[jParts].kterms - 1L);
	new->Knewterm = Knewterm;
	new->Kchol    = new->ks_x[0];

	new->isEqualSwap = 0;
	if (new->K == new->Knew) {
		new->isEqualSwap = 1;
		for (int i = 0; i < new->nbands; i++) {
			if (new->kterms_x[i] != new->kterms_xnewterm[i]) {
				new->isEqualSwap = 0;
				break;
			}
		}
	}	
}


void update_XtX_from_Xnewterm_v2(F32PTR XtX, F32PTR XtXnew, NEWCOLINFOv2* new ) {

 
	I32 N    = new->N;
	I32 Nlda = new->Nlda;

	I32 KOLD = new->K;
	I32 KNEW = new->Knew;

	if (new->isEqualSwap) {
		SCPY(KOLD * KOLD, XtX, XtXnew);
	}

	I32 colbase = 1;	
	for (int i = 0; i < (new->nbands * 2 + 1); i++) {		
		I32 Kcol = new->parts[i].kterms;
		if (Kcol == 0) {
			continue;
		}

		I32     rowbase = 1;
		F32PTR  Xcol0 = new->parts[i].X;
		F32PTR  Xcol  = Xcol0 + (new->parts[i].ks_src - 1) * Nlda;		
		for (int j = 0; j <= i; j++) {
			I32  Krow = new->parts[j].kterms;
			if (Krow == 0) {
				continue;
			}

			F32PTR  Xrow0  = new->parts[j].X;
			F32PTR  Xrow   = Xrow0 + (new->parts[j].ks_src - 1) * Nlda;
			if (Xcol0 != new->X || Xrow0 != new->X) {
			 // New XtX values needed to be computed
				r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans, 
					Krow, Kcol, N, 1.0f,
					Xrow, Nlda,	Xcol, Nlda, 0.f,
					XtXnew + (colbase - 1L) * KNEW + rowbase - 1, KNEW);
			}
			else if (new->isEqualSwap==0) {
			// Need to copy old values from XtX to XtXnew
			
				I32 row_old = new->parts[j].ks_src;
				I32 col_old = new->parts[i].ks_src;

				F32PTR ColStart_old = XtX    + (col_old - 1) * KOLD + row_old -1L;
				F32PTR ColStart_new = XtXnew + (colbase - 1) * KNEW + rowbase- 1L;
				if (i == j) {
					// Copy the existing values of the upper triangle from XtX to XtXnew
					for (int k = 1; k <= Kcol; ++k) {
						SCPY(k, ColStart_old + (k - 1) * KOLD, ColStart_new + (k - 1) * KNEW);
					}
				} else {
					// Copy the existing values of the full rectangle from XtX to XtXnew
					for (int k = 1; k <= Kcol; ++k) {						
						SCPY(Krow, ColStart_old + (k - 1) * KOLD, ColStart_new + (k - 1) * KNEW);
					}
				}
				

	 	   } // if (Xcol0 != new->X || Xrow0 != new->X) 	

		   rowbase += Krow;

		} //for (int j = 0; j <= i; j++)

		colbase += Kcol;
	}
	

			 
}

void update_XtY_from_Xnewterm_v2(F32PTR XtY, F32PTR XtYnew, F32PTR Y, NEWCOLINFOv2* new, I32 q) {

	// X and Xnewterm has a leading dimesnion of new.Nlada
	// Y has a leading dimension of new.N

	I32 N    = new->N;
	I32 Nlda = new->Nlda;

	I32 KOLD = new->K;
	I32 KNEW = new->Knew;

	if (new->isEqualSwap) {
		SCPY(KOLD*q, XtY, XtYnew);
	}

	if (q == 1) {
		I32     rowbase = 1;
		for (int i = 0; i < (new->nbands * 2 + 1); i++) {
			I32 Krow = new->parts[i].kterms;
			I32 row_old = new->parts[i].ks_src;
			if (Krow == 0) {
				continue;
			}

			F32PTR  Xrow0 = new->parts[i].X;
			F32PTR  Xrow = Xrow0 + (row_old - 1) * Nlda;
			if (Xrow0 != new->X) {
				// New XtX values needed to be computed
				r_cblas_sgemv(CblasColMajor, CblasTrans,
					N, Krow, 1.f,
					Xrow, Nlda,
					Y, 1L, 0.f,
					XtYnew + rowbase - 1, 1L);
			}
			else if (new->isEqualSwap == 0) {
				// Need to copy old values from XtX to XtXnew			
				SCPY(Krow, XtY + row_old - 1L, XtYnew + rowbase - 1);

			} // if (Xcol0 != new->X || Xrow0 != new->X) 	

			rowbase += Krow;
		}
	}
	else {	// q > 1

		I32     rowbase = 1;
		for (int i = 0; i < (new->nbands * 2 + 1); i++) {
			I32 Krow = new->parts[i].kterms;
			I32 row_old = new->parts[i].ks_src;
			if (Krow == 0) {
				continue;
			}

			F32PTR  Xrow0 = new->parts[i].X;
			F32PTR  Xrow = Xrow0 + (row_old - 1) * Nlda;
			if (Xrow0 != new->X) {
				// New XtX values needed to be computed
				r_cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
					Krow, q, N, 1.f,
					Xrow, Nlda,
					Y, N, 0.f,
					XtYnew + rowbase - 1, KNEW);
			}
			else if (new->isEqualSwap == 0) {
				// Need to copy old values from XtX to XtXnew							
				for (I32 c = 0; c < q; ++c) {
					SCPY(Krow, XtY +c*KOLD+ row_old - 1L, XtYnew + c * KNEW + rowbase - 1);
				}
			} // if (Xcol0 != new->X || Xrow0 != new->X) 	

			rowbase += Krow;
		}
	
	} //else {	// q > 1

 
}

void shift_last_elems(void * X, I32 Kstart, I32 Kend, I32 Knewstart, I32 elemSize) {
	// K are all one-based indices
	
	I08PTR x = X;
	if (Knewstart == Kstart) { 	// a. No need to move at all
		return;
	}	else if (Knewstart < Kstart || Knewstart > Kend) {	// b.Move forwward or no overlapping
		// dst(k2_new):-----123455----
		// src(k2_old):----------123455----

		// dst(k2_new):----------------123455----
		// src(k2_old):-----123455----		
		memcpy(x + (Knewstart - 1) * elemSize, x+(Kstart - 1) *elemSize,  (Kend - Kstart + 1) * elemSize);
		//memcpy(dst, src, size)
	}	else {	// c. Movee backward with some overlapping
		// dst(k2_new):-------------123456----
		// src(k2_old):---------123456----
		I32 segStartIdx = Kend + 1;	
		int nElemsCopy  = Knewstart-Kstart;  //j>0
		while (_True_) {
			segStartIdx = segStartIdx - nElemsCopy;
			if (segStartIdx > Kstart) {
				//memcpy(dst, src, size)
				memcpy(x + ((segStartIdx + nElemsCopy) - 1) * elemSize, x + (segStartIdx - 1) * elemSize, nElemsCopy * elemSize) ;
			} else {
				nElemsCopy = (segStartIdx + nElemsCopy) - Kstart;
				memcpy(x + (Knewstart - 1) * elemSize,  x+ (Kstart - 1) * elemSize, nElemsCopy * elemSize);
				break;
			}
		}//while (_True_)
	} //if (j < 0 || k2_new + 1 > Kold)

}

void swap_elem_bands(NEWCOLINFOv2* new, void *x, void *xnew, I32 elemSize) {

	U08PTR X    = x;
	U08PTR Xnew = xnew;

	if (new->isEqualSwap == 0) {
		// First, shift X's col bands
		// X1, Xnewterm1, X2, Xnewterm2, X3, Xnewtem4, Xn4
		// Need to shift X2, X3, and X4
		int    Kend    = new->K;
		int    Kadjust = 0;
		for (int i = 3; i <= (new->nbands * 2 + 1); i += 2) {
			if (new->parts[i - 1].kterms == 0) {
				continue;
			}
			int Kstart    = new->parts[i - 1].ks_src + Kadjust;
			int Knewstart = new->parts[i - 1].ks_dst;
			shift_last_elems(X,Kstart, Kend, Knewstart, elemSize);

			// then the original cols have been shifted
			Kadjust += (Knewstart - Kstart);
			Kend    += Kadjust;
		}

	}

	// Insert the Xnewterm nto the shiftted X 
	for (int i = 2; i <= (new->nbands * 2 + 1); i += 2) {
		if (new->parts[i - 1].kterms == 0) {
			continue;
		}
		int Knewterm = new->parts[i - 1].kterms;
		// memcpy(dst, src, size)
		memcpy(X + (new->parts[i - 1].ks_dst - 1) * elemSize, Xnew + (new->parts[i - 1].ks_src - 1) * elemSize, Knewterm * elemSize);		
	}

}

void shift_lastcols_within_matrix(F32PTR X, I32 N, I32 Kstart, I32 Kend, I32 Knewstart) {

	int offset = Knewstart - Kstart;
	if (offset == 0) {
		return;
	}

	int Kseg = Kend - Kstart + 1;
	if ( offset>=Kseg || offset <= -Kseg) { // Not overlapping	
		// dst(k2_new):--------------------{Knewstart]123455----  [CASE 1]
		// src(k2_old):-----123455[Kend]----

		// dst(k2_new):-----123455----                            [CASE 2] 
		// src(k2_old):---------------[Ksart]123455----
		r_cblas_scopy(Kseg * N, X + (Kstart - 1) * N, 1, X + (Knewstart - 1) * N, 1);
	}	
	else if (offset < 0) { // Overlaping [ -Kseg <offset <0]
		// Use of memcpy in this case is not safe, having undefined behavior accroding to the C standard
		// dst(k2_new):-----123455----
		// src(k2_old):----------123455----		
		memmove( X+(Knewstart - 1) * N, X + (Kstart - 1) * N, Kseg * N * sizeof(F32));
	}
    else {  // Overlaping [ 0 <offset < Kseg]
		// dst(k2_new):-------------123456----
		// src(k2_old):---------123456----
		I32 segStartIdx = Kend + 1;
		while (_True_) {
			segStartIdx = segStartIdx -offset;
			if (segStartIdx > Kstart) {
				SCPY(offset * N, X + (segStartIdx - 1) * N, X + ((segStartIdx + offset) - 1) * N);
			} else {
				int Kremaining = (segStartIdx + offset) - Kstart;
				SCPY(Kremaining * N, X + (Kstart - 1) * N, X + (Knewstart - 1) * N);
				break;
			}
		}//while (_True_)
	}//if (j < 0 || k2_new + 1 > Kold)

}

void swap_cols_bands_within_matrx(NEWCOLINFOv2* new) {
	
	F32PTR X = new->X;

	if (new->isEqualSwap == 0) {
		// First, shift X's col bands
		// X1, Xnewterm1, X2, Xnewterm2, X3, Xnewtem4, Xn4
		// Need to shift X2, X3, and X4
		int    Kend    = new->K;
		int    Kadjust = 0;		
		for (int i = 3; i <= (new->nbands * 2 + 1); i += 2) {
			if (new->parts[i - 1].kterms == 0) {
				continue;
			}
			int Kstart     = new->parts[i - 1].ks_src + Kadjust;
			int Knewstart  = new->parts[i - 1].ks_dst;
			shift_lastcols_within_matrix(X, new->Nlda, Kstart, Kend, Knewstart);
						
			// then the original cols have been shifted
			Kadjust += (Knewstart - Kstart);
			Kend    += Kadjust;		
		}

	}

	// Insert the Xnewterm nto the shiftted X
	int N = new->Nlda;
	for (int i = 2; i <= (new->nbands * 2 + 1); i += 2) {
		if (new->parts[i - 1].kterms == 0) {
			continue;
		}
		int Knewterm = new->parts[i - 1].kterms;
		SCPY(Knewterm*N, new->Xnewterm + (new->parts[i-1].ks_src- 1) * N, X + (new->parts[i - 1].ks_dst - 1) * N);
	}
		
}
#include "abc_000_warning.h"
 

