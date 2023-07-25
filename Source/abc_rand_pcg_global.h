#pragma once
#include "abc_datatype.h"

 


//LCG:
// X     =  X    * multiplier + shift(plus/state.inc) mod 2^64
// State =  State* multiplier + shift(plus/INC)       mod 2^64
extern void pcg_get_lcg_multiplier_shift_multistep(I32 delta, U64 cur_mult, U64 cur_plus, U64* acc_mult, U64* acc_shift);

void pcg_set_seed(U64 initstate,     U64 initseq);
void pcg_random(U32 * _restrict rnd, I32 N);
void pcg_gauss(F32PTR RND,			I32 N);
void pcg_gamma(F32PTR rnd,			F32 a, I32 N);
extern void pcg_wishart_unit_lowtriangle_zeroout(F32PTR rnd, F32PTR tmp, I32 m, F32 df);
extern void pcg_wishart_unit_lowtriangle_zeroout_notmp(F32PTR wishrnd, I32 m, F32 df);
extern void pcg_invwishart_upper(F32PTR iwrnd_upper, F32PTR iwrnd_upper_inv, F32PTR tmp, I32 m, F32PTR Qup, F32 df);



/*****************************************************************************************/
/*                                                                                        */
/*****************************************************************************************/
typedef struct GAUSS_CONSTANT {
	F32    x[64];
	F32    yRatio[63];	
	I16    indices[125];
	F32    amax;
	F32    exp_lamda;
	I32    inflectionId;
} GAUSS_CONSTANT;
extern GAUSS_CONSTANT GAUSS;