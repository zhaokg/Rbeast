#pragma once
#include "abc_datatype.h"

 
// LCG:
// X     =  X    * multiplier + shift(plus/state.inc) mod 2^64
// State =  State* multiplier + shift(plus/INC)       mod 2^64

typedef union local_pcg32_random_struct {
	//https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.h
	/**************/
	//for pcg_avx
	/**************/
	struct {
	U64 ALIGN32_BEG
		state[4]
		ALIGN32_END;  // RNG state. All values are possible.
	U64	increment;	  // Controls which RNG sequence (stream) is
					  // selected. Must *always* be odd.
	U64 MULTIPLIER_4steps;
	U64 INCREMENT_4steps;
	
	U32 INTERNAL_RNDBUF[4];  // Used in pcg_random_with_internalbuf to save unused random numbers
	int BUF_PTR;
	};

	/**************/
	//for pcg_avx512
	/**************/
	struct {
		U64 ALIGN32_BEG
			state512[8]
			ALIGN32_END;    // RNG state. All values are possible.
		U64	increment512;	// Controls which RNG sequence (stream) is
						    // selected. Must *always* be odd.
		U64 MULTIPLIER_8steps;
		U64 INCREMENT_8steps;
	};

	// forl pcg_local_generic
	struct {
		U64 STATE;
		U64	INCREMENT;	// Controls which RNG sequence (stream) is
	};
} local_pcg32_random_t;

extern void (*local_pcg_set_seed)(local_pcg32_random_t* rng, U64 initstate, U64 initseq);
extern void (*local_pcg_random)(local_pcg32_random_t* rng, U32PTR rnd, I32 N);
extern void (*local_pcg_print_state)(local_pcg32_random_t* rng);

extern F32  local_pcg_trgauss_oneside_scalar(local_pcg32_random_t* rng, F32 a, int whichside);
extern  void local_pcg_randint(local_pcg32_random_t* rng, int bnd, I32PTR RND, int N);
extern void local_pcg_gauss(local_pcg32_random_t* rng, F32PTR RND, int N);
extern void local_pcg_gamma(local_pcg32_random_t* rng, F32PTR rnd, F32 a, int N);
extern void local_pcg_wishart_unit_lowtriangle_zeroout(local_pcg32_random_t* rng, F32PTR rnd, F32PTR tmp, I32 m, F32 df);
extern void local_pcg_wishart_unit_lowtriangle_zeroout_notmp(local_pcg32_random_t* rng, F32PTR wishrnd, I32 m, F32 df);
extern void local_pcg_invwishart_upper(local_pcg32_random_t* rng, F32PTR iwrnd_upper, F32PTR iwrnd_upper_inv, F32PTR tmp, I32 m, F32PTR Qup, F32 df);
extern void SetupPCG_GENERIC(void);
extern void SetupPCG_AVX2(void);
extern void SetupPCG_AVX512(void);


// r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, k, beta, 0, 1);
#define VSLStreamStatePtr													  local_pcg32_random_t														 
#define r_vslNewStream(stream, METHOD_dummy, seed)							  local_pcg_set_seed(stream, 0x853c49e6748fea9bULL, seed)  //	pcg_set_seed(seed, (size_t) &globalStream)
#define r_viRngUniformBits32(  METHOD_dummy, stream, N, rnd)				  local_pcg_random(&stream,rnd, N)
#define r_vsRngGamma(		   METHOD_dummy, stream, N, rnd, a, b, beta_dummy) local_pcg_gamma(&stream, rnd, a,N)
#define r_vsRngGaussian(	   MTEHDO_dummy, stream,N, rnd, a, b)			  local_pcg_gauss(&stream, rnd, N)
#define r_vslDeleteStream(stream_dummy)  

/*
#define VSLStreamStatePtr													 local_pcg32_random_t														 
#define r_vslNewStream(stream, METHOD_dummy, seed)							 pcg_set_seed( 0x853c49e6748fea9bULL,seed)  //	pcg_set_seed(seed, (size_t) &globalStream)
#define r_viRngUniformBits32(METHOD_dummy, stream, N, rnd)					 pcg_random(rnd, N)
#define r_vsRngGamma(		 METHOD_dummy, stream, N, rnd, a, b, beta_dummy) pcg_gamma( rnd, a,N)
#define r_vsRngGaussian(	 MTEHDO_dummy,stream,N, rnd, a, b)				 pcg_gauss( rnd, N)
#define r_vslDeleteStream(stream_dummy)  
*/