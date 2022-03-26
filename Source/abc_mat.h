#pragma once
#include <math.h>  
#include <string.h>  

#include "abc_datatype.h" 
void chol_update_U(F32PTR U, F32PTR x, I32 ldu, I32 n);
void chol_dwdate_U(F32PTR U, F32PTR x, I32 ldu, I32 n);
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

void chol_addCol_skipleadingzeros_prec_invdiag(F32PTR Au, F32PTR U,  F32PTR precPrior, I64 N, I64 K0, I64 K1);
void chol_addCol_skipleadingzeros_precVec_invdiag(   F32PTR Au,  F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1);



void linear_regression(F32PTR Y, F32PTR X, int ldx, int N, int K, F32PTR B, F32PTR Yfit, F32PTR Yerror, F32PTR TMP);
