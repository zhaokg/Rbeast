#pragma once
#include <math.h>  
#include <string.h>  

#include "abc_datatype.h" 
extern void chol_update_U(F32PTR U, F32PTR x, I32 ldu, I32 n);
extern void chol_dwdate_U(F32PTR U, F32PTR x, I32 ldu, I32 n);
extern void chol_columwise(F32PTR Au, F32PTR U, I64  N, I64 K);
extern void chol_columwise_v2(F32PTR Au, F32PTR U, I64  N, I64 K);
extern void chol_rowwise(F32PTR Au, F32PTR U, I64  N, I64 K);
extern void chol_addCol(F32PTR Au, F32PTR U, I64  N, I64 K0, I64 K1);
extern void inplace_chol(F32PTR A, I64  N, I64 K);
extern void inplace_chol_addCol(F32PTR A, I64 N, I64 K0, I64 K1);
extern void chol_full_update(F32PTR x, F32PTR  U, I64 N, I64 K);
extern void chol_full_downdate(F32PTR x, F32PTR  U, I64 N, I64 K);
extern void solve_U_as_U(F32PTR A, F32PTR x, I64 N, I64 K);
extern void solve_U_as_L(F32PTR A, F32PTR x, I64 N, I64 K);
extern void solve_L_as_L(F32PTR A, F32PTR x, I64 lda, I64 K);
extern void solve_L_as_U(F32PTR A, F32PTR x, I64 lda, I64 K);

extern void pack_chol(F32PTR Au, F32PTR U, I64  N);
extern void pack_chol_addCol(F32PTR Au, F32PTR U, I64 K0, I64 K1);
extern void pack_chol_update(F32PTR x, F32PTR  U, I64 K);
extern void pack_chol_dwdate(F32PTR x, F32PTR  U, I64 K);

extern void pack_inplace_chol(F32PTR A, I64  N);
extern void pack_inplace_chol_addCol(F32PTR A, I64 K0, I64 K1);

extern void pack_solve_U(F32PTR A, F32PTR x, I64 K);
extern void pack_solve_L(F32PTR A, F32PTR x, I64 K);
///////////////////////////
void chol_addCol_skipleadingzeros(F32PTR Au, F32PTR U, I64 N, I64 K0, I64 K1);
void chol_addCol_skipleadingzeros_prec(F32PTR Au, F32PTR U, F32 precPrior, I64 N, I64 K0, I64 K1);
void solve_U_as_LU(F32PTR U, F32PTR y, F32PTR x, I64 N, I64 K);
void solve_U_as_LU_rectmat_multicols(F32PTR U, F32PTR y, F32PTR x, I64 ldu, I64 K, I64 nCols);
void solve_U_as_LU_invdiag_rectmat(F32PTR U, F32PTR y, F32PTR x, I64 ldu, I64 K);
void solve_U_as_LU_invdiag_sqrmat(F32PTR U, F32PTR y, F32PTR x,  I64 K);
void solve_U_as_LU_invdiag_sqrmat_multicols(F32PTR U, F32PTR y, F32PTR x, I64 K, I64 nColY);
void solve_U_as_U_invdiag(F32PTR U, F32PTR x, I64 N, I64 K);
void solve_U_as_U_invdiag_multicols(F32PTR U, F32PTR x, I64 ldu, I64 K, I32 nColx);

void chol_addCol_skipleadingzeros_prec_nostartprec_invdiag(F32PTR Au, F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1);
void chol_addCol_skipleadingzeros_prec_invdiag(F32PTR Au, F32PTR U,  F32PTR precPrior, I64 N, I64 K0, I64 K1);
void chol_addCol_skipleadingzeros_precVec_invdiag(   F32PTR Au,  F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1);



void linear_regression(F32PTR Y, F32PTR X, int ldx, int N, int K, F32PTR B, F32PTR Yfit, F32PTR Yerror, F32PTR TMP);
void simple_linear_regression_nan(F32PTR Y, F32PTR X, int N, F32PTR Yfit, F32PTR Yerror);

typedef struct {
   I32 N;
   I32 Nlda; // the leading dimensio of the X

   I16 k1;   // k1_old=k1_new;
   I16 k2_old, k2_new;
   I16 Knewterm;
   I16 KOLD,  KNEW;
   //I15 KNew =KOLD+K2_new+k2_old;
} NEWCOLINFO, * _restrict NEWCOLINFO_PTR;


extern void update_XtX_from_Xnewterm(F32PTR X, F32PTR Xnewterm, F32PTR XtX, F32PTR XtXnew, NEWCOLINFO* new);

extern void update_XtY_from_Xnewterm(F32PTR Y, F32PTR Xnewterm, F32PTR XtY, F32PTR XtYnew, NEWCOLINFO* new, I32 q);



typedef struct {
	I32 N;
	I32 Nlda; // the leading dimensio of the X

	F32PTR X;
	F32PTR Xnewterm;

	I16 nbands;
	I16 ks_x[5];
	I16 kterms_x[5];
	I16 ks_xnewterm[5];
	I16 kterms_xnewterm[5];

	struct {
		F32PTR X;
		I16    ks_src;
		I16    kterms;
		I16    ks_dst;
	} parts[11];

	I16 K;       // the col dim of X
	I16 Knewterm;// the col dim of Xnewterm
	I16 Knew;    // the col dim of Xnew (after inserting Xnewterm into X)
	I16 Kchol;   // the start col index for the chol decomposition

	I16 isEqualSwap;

} NEWCOLINFOv2, * _restrict NEWCOLINFOvs_PTR;
extern void get_parts_for_newinfo(NEWCOLINFOv2* new);
extern void update_XtX_from_Xnewterm_v2(F32PTR XtX, F32PTR XtXnew, NEWCOLINFOv2* new); 
extern void update_XtY_from_Xnewterm_v2(F32PTR XtY, F32PTR XtYnew, F32PTR Y, NEWCOLINFOv2* new, I32 q);


extern void shift_last_elems(void* X, I32 Kstart, I32 Kend, I32 Knewstart, I32 elemSize);
extern void swap_elem_bands(NEWCOLINFOv2* new, void* x, void* xnew, I32 elemSize);

extern void shift_lastcols_within_matrix(F32PTR X, I32 N, I32 Kstart, I32 Kend, I32 Knewstart);
extern void swap_cols_bands_within_matrx(NEWCOLINFOv2* new);