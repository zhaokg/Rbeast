#include <math.h>  
#include "abc_000_warning.h"

#include "abc_datatype.h" 
#include "abc_ide_util.h"
#include "abc_rand_pcg_global.h" 
#include "abc_rand_pcg_local.h" 


//https://hal.archives-ouvertes.fr/hal-02700791/document: a nice intro to the pcg alg
/*Practical seed-recovery for the PCG Pseudo-Random Numbr Generator*/

/*
Internally, this RNG uses two 64-bit integers for its internal state, consisting of:
the current state — the RNG iterates through all 264 possible internal states.
the RNG-sequence constant — a value that defines which of 263 possible random sequences
the current state is iterating through; it holds the same value over the lifetime of the RNG.
*/

/*
For this generator, there are 2^63 possible sequences of pseudorandom numbers. Each sequence is
entirely distinct and has a period of 2^64. The initseq argument selects which stream you will use. 
The initstate argument specifies where you are in that 264 period.
*/

#define PCG_DEFAULT_MULTIPLIER_64  6364136223846793005ULL //https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.c
#define PCG_DEFAULT_INCREMENT_64   1442695040888963407ULL ////https://github.com/imneme/pcg-c/blob/master/include/pcg_variants.h


//https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.h
#define PCG_DEFAULT_GLOBAL_STATE_64     0x853c49e6748fea9bULL
#define PCG_DEFAULT_GLOBAL_INCREMENT_64 0xda3e39cb94b95bdbULL

// Linear Congruential Generator (LCG):
// X     =  X    * MULTIPLIER + INCREMENT(shift/plus/state.inc) mod 2^64
// State =  State* MULTIPLIER + INCREMENT(shift/plus/INC)       mod 2^64


/*
struct { 
	U64 STATE;     // RNG state.All values are possible.
	U64	INCREMENT;	// Controls which RNG sequence (stream) is 
					// selected. Must *always* be odd.
};
 */

void gen_pcg_print_state(local_pcg32_random_t* rng) {

	r_printf("PCG State: \n");
	r_printf("State: %"  PRIx64 "\n", rng->STATE);	
	r_printf("Increment: %"  PRIx64  "\n", rng->INCREMENT);
 
}


void gen_pcg_random(local_pcg32_random_t* rng, U32PTR rnd, I32 N);
void gen_pcg_set_seed(local_pcg32_random_t* rng, U64 initstate, U64 initseq)
{

	initstate = PCG_DEFAULT_GLOBAL_STATE_64 ^ initseq; //Added bcz only initseq is supplied as a seed. We run initseq to randomize initstate a little bit

	initstate = initstate == 0 ? PCG_DEFAULT_GLOBAL_STATE_64     : initstate;
	initseq   = initseq == 0   ? PCG_DEFAULT_GLOBAL_INCREMENT_64 : initseq;

	rng->STATE		 = 0U;
	rng->INCREMENT   = (initseq << 1u) | 1u;  	//inc must be an odd number

	U32 rnd;
	gen_pcg_random(rng, &rnd, 1);
	rng->STATE += initstate;
	gen_pcg_random(rng, &rnd, 1);

	extern void init_gauss_rnd(void);
	init_gauss_rnd(); //Indepedent of PCG, used to initialize the GAUSS structure
}
 
void gen_pcg_random(local_pcg32_random_t* rng, U32PTR rnd, I32 N) 
{
	U64 oldstate = rng->STATE;
	U64 shift    = rng->INCREMENT;

	for (int i = 0; i < N; i++)
	{
		U32 xorshifted =     ((oldstate >> 18u) ^ oldstate) >> 27u   ;
		U32 rot = oldstate >> 59u;
		//*rnd++ = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
		//https://en.wikipedia.org/wiki/Circular_shift
		*rnd++ = (xorshifted >> rot) | (xorshifted << ((- *((I32*)(&rot))  ) & 31)); 
		oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;
	}
	rng->STATE = oldstate;
}

void SetupPCG_GENERIC() {

	local_pcg_set_seed = gen_pcg_set_seed;
	local_pcg_random   = gen_pcg_random;
	local_pcg_print_state = gen_pcg_print_state;
}
/**********************************************/
#include "abc_000_warning.h"
