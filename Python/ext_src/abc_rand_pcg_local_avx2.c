#include <stdio.h> //fprintf FILE fopen
#include <string.h>
#include <math.h>  
//#include <stdlib.h> //malloc and free

#include "abc_000_warning.h"


#include "abc_datatype.h"
#include "abc_ide_util.h"

#include "abc_vec.h"

#include "abc_rand_pcg_global.h"
#include "abc_rand_pcg_local.h"
 
#include "assert.h"

#ifdef COMPILER_MSVC
#define __attribute__(x)
//https://stackoverflow.com/questions/21116270/gcc-attributes-with-c-methods
// In a function definition, the attribute field should go before the function name
//, but it can be put after the function decleration in the decleration
#endif


///////////////////////////////////////////////////////////////////////////
//stackoverflow.com/questions/2622017/suppressing-deprecated-warnings-in-xcode
/*
#ifdef COMPILER_CLANG
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wdeprecated-declarations"
			//do something////
	#pragma clang diagnostic pop
#endif
#ifdef COMPILER_GCC
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
			//do something////
	#pragma GCC diagnostic pop
#endif
*/
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//https://clickhouse.tech/codebrowser/html_report/ClickHouse/src/Functions/TargetSpecific.h.html
#if  defined(COMPILER_CLANG) && !defined(cpu_ARM64)  

    //https://stackoverflow.com/questions/31373885/how-to-change-optimization-level-of-one-function/49353441
    #pragma clang optimize on
    //https://stackoverflow.com/questions/46165752/does-clang-have-something-like-pragma-gcc-target
    #pragma clang attribute push (__attribute__((target("sse,sse2,sse3,ssse3,sse4,popcnt,avx,fma,avx2"))), apply_to=function)
	//#pragma clang attribute push (__attribute__((target("avx,avx2,avx512f,avx512dq"))), apply_to=function)
    //#pragma clang attribute pop
#endif

#if  defined(COMPILER_GCC) && !defined(cpu_ARM64)  
    //https://www.geeksforgeeks.org/speed-up-naive-algorithms-in-competitive-coding-in-c-cpp/
    //https://codeforces.com/blog/entry/78897
    //https://stackoverflow.com/questions/61759552/why-some-top-level-competitive-programmers-use-pragma    
    //#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
    //#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
    #pragma optimization_level 3
    #pragma GCC optimize("O3,Ofast,inline,omit-frame-pointer,no-asynchronous-unwind-tables")  //Comment optimisations for interactive problems (use endl)
    //#pragma GCC target("avx,avx2,fma")
     #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,avx,avx2,fma,tune=haswell")
	 //#pragma GCC target("avx,avx2,avx512f,avx512dq") //https://en.wikipedia.org/wiki/AVX-512
    // -mavx256-split-unaligned-load is spit out by default, hurting new CPUS
    // an altertive is to set -march=hashwell, but it doesntot work in pragram for all gcc versions
    // #pragma GCC optimization ("unroll-loops")
    //#pragma GCC target "arch=core-avx2,tune=core-avx2"
     
    // https://stackoverflow.com/questions/51003218/gcc-target-for-avx2-disabling-sse-instruction-set
    // #pragma GCC target "arch=haswell" : // A buf for some gcc versions, not working
    
#endif
///////////////////////////////////////////////////////////////////////////


#if !defined(COMPILER_SOLARIS) && defined(TARGET_64) && !defined(cpu_ARM64)

#include "abc_math_avx.h"
// LCG  X  = X * Multiplier + Shift(plus/state.inc) MOD 2^64
#define PCG_DEFAULT_MULTIPLIER_64  6364136223846793005ULL //https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.c
#define PCG_DEFAULT_INCREMENT_64   1442695040888963407ULL ////https://github.com/imneme/pcg-c/blob/master/include/pcg_variants.h


//https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.h
#define PCG_DEFAULT_GLOBAL_STATE_64     0x853c49e6748fea9bULL
#define PCG_DEFAULT_GLOBAL_INCREMENT_64 0xda3e39cb94b95bdbULL

 
void avx_pcg_print_state(local_pcg32_random_t* rng) {

	r_printf("PCG State: \n");
	r_printf("State: %"  PRIx64 " %" PRIx64 " %" PRIx64 " %" PRIx64 "\n", rng->state[0], rng->state[1], rng->state[2], rng->state[3] );
	r_printf("Increment: %"  PRIx64  "\n", rng->increment);
	r_printf("Mutiplier4: %"  PRIx64  "\n", rng->MULTIPLIER_4steps);
	r_printf("INcrementr4: %"  PRIx64  "\n\n", rng->INCREMENT_4steps);
}

 void avx_pcg_set_seed(local_pcg32_random_t *rng, U64 initstate, U64 initseq)
{
	
	 initstate = PCG_DEFAULT_GLOBAL_STATE_64 ^ initseq; //Added bcz only initseq is supplied as a seed. We run initseq to randomize initstate a little bit

	initstate = initstate == 0 ? PCG_DEFAULT_GLOBAL_STATE_64 : initstate;
	initseq   = initseq == 0 ?   PCG_DEFAULT_GLOBAL_INCREMENT_64 : initseq;

	rng->increment = (initseq << 1u) | 1u; //inc must be an odd number
	U64 state = 0U;						
	state   = state * PCG_DEFAULT_MULTIPLIER_64 + rng->increment;	// pcg_random(&rnd, 1);
	state  += initstate;	
	rng->state[0] = state = state * PCG_DEFAULT_MULTIPLIER_64 + rng->increment;// pcg_random(&rnd, 1);
	rng->state[1] = state = state * PCG_DEFAULT_MULTIPLIER_64 + rng->increment;// pcg_random(&rnd, 1);
	rng->state[2] = state = state * PCG_DEFAULT_MULTIPLIER_64 + rng->increment;// pcg_random(&rnd, 1);
	rng->state[3] = state = state * PCG_DEFAULT_MULTIPLIER_64 + rng->increment;// pcg_random(&rnd, 1);
	__m256i state256 = _mm256_set_epi64x(rng->state[3], rng->state[2], rng->state[1], rng->state[0]);;

	//r_printf("state: %30" PRIu64 " inc: %30"  PRIu64 "\n", rng->state[0], rng->increment);

	// Compute the multiper and shift constants for a four-step advance
	pcg_get_lcg_multiplier_shift_multistep(4L, PCG_DEFAULT_MULTIPLIER_64, rng->increment, &rng->MULTIPLIER_4steps, &rng->INCREMENT_4steps);

	rng->BUF_PTR = 4;// used in the pcg_random_internalbuf (4 meaans the buf is empty: 0,1,2,3: there are still data)

	extern void init_gauss_rnd(void);
	init_gauss_rnd(); //Indepedent of PCG, used to initialize the GAUSS structure
}

 static INLINE  __m256i  __attribute__((always_inline)) GetMoveMask(int n)    {
    __m128i maskIdx = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i maskNum = _mm256_cvtepu8_epi32(maskIdx);
    __m256i nvec    = _mm256_set1_epi32(n);
    //__m256i maskmov = _mm256_sub_epi32(maskNum, nvec);
    __m256i maskmov = _mm256_cmpgt_epi32(nvec, maskNum);
    return maskmov;
}

#if defined(COMPILER_MSVC)

 // replace hadd -> shuffle (4 uops) with shift/add/and (3 uops)
// The constant takes 2 insns to generate outside a loop.
 static INLINE __m256i __attribute__((always_inline)) __mul64_haswell(__m256i a, __m256i b)
 {
	 //https://stackoverflow.com/questions/37296289/fastest-way-to-multiply-an-array-of-int64-t

	 /*
	 // instruction does not exist. Split into 32-bit multiplies
	 __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);           // swap H<->L
	 __m256i prodlh = _mm256_mullo_epi32(a, bswap);            // 32 bit L*H products

	 // or use pshufb instead of psrlq to reduce port0 pressure on Haswell
	 __m256i prodlh2 = _mm256_srli_epi64(prodlh, 32);          // 0  , a0Hb0L,          0, a1Hb1L
	 __m256i prodlh3 = _mm256_add_epi32(prodlh2, prodlh);      // xxx, a0Lb0H+a0Hb0L, xxx, a1Lb1H+a1Hb1L
	 __m256i prodlh4 = _mm256_and_si256(prodlh3, _mm256_set1_epi64x(0x00000000FFFFFFFF)); // zero high halves

	 __m256i prodll = _mm256_mul_epu32(a, b);                  // a0Lb0L,a1Lb1L, 64 bit unsigned products
	 __m256i prod = _mm256_add_epi64(prodll, prodlh4);       // a0Lb0L+(a0Lb0H+a0Hb0L)<<32, a1Lb1L+(a1Lb1H+a1Hb1L)<<32
	 return  prod;
	 */
	 __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);           // swap H<->L
	 __m256i prodlh = _mm256_mullo_epi32(a, bswap);            // 32 bit L*H products
	 __m256i zero = _mm256_setzero_si256();                 // 0
	 __m256i prodlh2 = _mm256_hadd_epi32(prodlh, zero);         // a0Lb0H+a0Hb0L,a1Lb1H+a1Hb1L,0,0
	 __m256i prodlh3 = _mm256_shuffle_epi32(prodlh2, 0x73);     // 0, a0Lb0H+a0Hb0L, 0, a1Lb1H+a1Hb1L
	 __m256i prodll = _mm256_mul_epu32(a, b);                  // a0Lb0L,a1Lb1L, 64 bit unsigned products
	 __m256i prod = _mm256_add_epi64(prodll, prodlh3);       // a0Lb0L+(a0Lb0H+a0Hb0L)<<32, a1Lb1L+(a1Lb1H+a1Hb1L)<<32
	 return  prod;
 }
#else
 
 // replace hadd -> shuffle (4 uops) with shift/add/and (3 uops)
// The constant takes 2 insns to generate outside a loop.
 static INLINE __m256i __mul64_haswell_ptr(__m256i * pa, __m256i *pb)
 {
	 __m256i a = _mm256_loadu_si256(pa);
	 __m256i b = _mm256_loadu_si256(pb);
	 //https://stackoverflow.com/questions/37296289/fastest-way-to-multiply-an-array-of-int64-t

	 /*
	 // instruction does not exist. Split into 32-bit multiplies
	 __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);           // swap H<->L
	 __m256i prodlh = _mm256_mullo_epi32(a, bswap);            // 32 bit L*H products

	 // or use pshufb instead of psrlq to reduce port0 pressure on Haswell
	 __m256i prodlh2 = _mm256_srli_epi64(prodlh, 32);          // 0  , a0Hb0L,          0, a1Hb1L
	 __m256i prodlh3 = _mm256_add_epi32(prodlh2, prodlh);      // xxx, a0Lb0H+a0Hb0L, xxx, a1Lb1H+a1Hb1L
	 __m256i prodlh4 = _mm256_and_si256(prodlh3, _mm256_set1_epi64x(0x00000000FFFFFFFF)); // zero high halves

	 __m256i prodll = _mm256_mul_epu32(a, b);                  // a0Lb0L,a1Lb1L, 64 bit unsigned products
	 __m256i prod = _mm256_add_epi64(prodll, prodlh4);       // a0Lb0L+(a0Lb0H+a0Hb0L)<<32, a1Lb1L+(a1Lb1H+a1Hb1L)<<32
	 return  prod;
	 */
	 __m256i bswap = _mm256_shuffle_epi32(b, 0xB1);           // swap H<->L
	 __m256i prodlh = _mm256_mullo_epi32(a, bswap);            // 32 bit L*H products
	 __m256i zero = _mm256_setzero_si256();                 // 0
	 __m256i prodlh2 = _mm256_hadd_epi32(prodlh, zero);         // a0Lb0H+a0Hb0L,a1Lb1H+a1Hb1L,0,0
	 __m256i prodlh3 = _mm256_shuffle_epi32(prodlh2, 0x73);     // 0, a0Lb0H+a0Hb0L, 0, a1Lb1H+a1Hb1L
	 __m256i prodll = _mm256_mul_epu32(a, b);                  // a0Lb0L,a1Lb1L, 64 bit unsigned products
	 __m256i prod = _mm256_add_epi64(prodll, prodlh3);       // a0Lb0L+(a0Lb0H+a0Hb0L)<<32, a1Lb1L+(a1Lb1H+a1Hb1L)<<32
	 return  prod;
 }

#define __mul64_haswell(a,b) __mul64_haswell_ptr(&(a), &(b))
#endif
void avx_pcg_random(local_pcg32_random_t* rng, U32PTR rnd, I32 N) {

	const __m256i	INCREMENT_SHIFT		= _mm256_set1_epi64x(rng->INCREMENT_4steps);
	const __m256i	MULITPLIER	        = _mm256_set1_epi64x(rng->MULTIPLIER_4steps);
 
	#define srl	_mm256_srli_epi64
	#define xor _mm256_xor_si256

	__m256i			oldstate = _mm256_loadu_si256(rng->state);

	I32 N4 = (N + 3) / 4 * 4;		
	for (int i = 0; i < N4; i+=4)	{
		__m256i xorshifted =  srl ( xor(srl(oldstate, 18u), oldstate), 27u) ;
		__m256i rot		   =  srl(oldstate, 59u);	
		oldstate		   =  _mm256_add_epi64(__mul64_haswell(oldstate, MULITPLIER), INCREMENT_SHIFT); //oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;
	 
		__m256i result= _mm256_or_si256(
							_mm256_srlv_epi32(xorshifted, rot),							
							_mm256_sllv_epi32(xorshifted,  _mm256_sub_epi32(_mm256_set1_epi32(32), rot)) //	//https://en.wikipedia.org/wiki/Circular_shift
						);	//*rnd++ = (xorshifted >> rot) | (xorshifted << ((- *((I32*)(&rot))  ) & 31)); //*rnd++ = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));

		__m128i r1 = _mm256_castsi256_si128(result);
		__m128i r2 = _mm256_extracti128_si256(result, 1);		
		__m128  r  = _mm_shuffle_ps(_mm_castsi128_ps(r1), _mm_castsi128_ps(r2), _MM_SHUFFLE(2, 0, 2, 0));

		 
		if (i < N - 3) {
			_mm_storeu_ps(rnd + i, r);
		} else {
			int n = N - i;
			 
			__m256i mask = GetMoveMask(n);
			_mm_maskstore_ps(rnd + i, _mm256_castsi256_si128(mask), r);
			 
		}

	}
	_mm256_storeu_si256(rng->state, oldstate);

	_mm256_zeroupper();
}
 
void avx_pcg_random_with_internalbuf(local_pcg32_random_t* rng, U32PTR rnd, I32 N) {

// Move these two local variable  to the rng structure so they can be initizeled during pcg_set_seed
//	static  __m128 INTERNAL_RNDBUF;
//	static  I32    BUF_PTR = 4;   // 4 means theare are no data: (0,1,2,3: there are data left)

	__m256i oldstate = _mm256_loadu_si256(rng->state);

	// N is changing, indicating the remaining number
	while (N > 0) {
		
		if (rng->BUF_PTR < 4) {
        // there is still data in the buffer; first consume the buffer
			U32PTR rndbuf = (U32PTR) rng->INTERNAL_RNDBUF + rng->BUF_PTR;
			if (N == 1L || rng->BUF_PTR ==3) {
				rnd[0] = rndbuf[0];
				++rng->BUF_PTR;
				++rnd;
				--N;				
			}
			else if (rng->BUF_PTR == 0 && N >= 4) {
				_mm_storeu_ps(rnd, _mm_loadu_ps(rng->INTERNAL_RNDBUF));
				rng->BUF_PTR = 4L;
				rnd     += 4;
				N       -= 4;
				// THen run to the bottom part to replish the buffer
			} else {
			// nCopy must be one of 2 and 3.
				int nAvailable = (4 - rng->BUF_PTR);
				int nCopy      = min(nAvailable, N);
				rnd[0] = rndbuf[0];
				rnd[1] = rndbuf[1];
                if (nCopy == 3) {
					rnd[2] = rndbuf[2];
				}
				rng->BUF_PTR += nCopy;
				rnd          += nCopy;
				N            -= nCopy;
			}
		}
	
		if (rng->BUF_PTR != 4) {
			//assert(N == 0);
			continue;
		}
 
		// Jump here to sample if BUF_PTR == 4

		const __m256i	INCREMENT_SHIFT = _mm256_set1_epi64x(rng->INCREMENT_4steps);
		const __m256i	MULITPLIER     = _mm256_set1_epi64x(rng->MULTIPLIER_4steps);

		#define srl	_mm256_srli_epi64
		#define xor _mm256_xor_si256

		__m256i xorshifted = srl(xor (srl(oldstate, 18u), oldstate), 27u);
		__m256i rot = srl(oldstate, 59u);
		oldstate = _mm256_add_epi64(__mul64_haswell(oldstate, MULITPLIER), INCREMENT_SHIFT); //oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;

		__m256i result = _mm256_or_si256(
			_mm256_srlv_epi32(xorshifted, rot),
			_mm256_sllv_epi32(xorshifted, _mm256_sub_epi32(_mm256_set1_epi32(32), rot)) //	//https://en.wikipedia.org/wiki/Circular_shift
		);	//*rnd++ = (xorshifted >> rot) | (xorshifted << ((- *((I32*)(&rot))  ) & 31)); //*rnd++ = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));

		__m128i r1 = _mm256_castsi256_si128(result);
		__m128i r2 = _mm256_extracti128_si256(result, 1);
		__m128  r = _mm_shuffle_ps(_mm_castsi128_ps(r1), _mm_castsi128_ps(r2), _MM_SHUFFLE(2, 0, 2, 0));

		_mm_storeu_ps(rng->INTERNAL_RNDBUF, r);

		rng->BUF_PTR = 0;
	 

	}
	_mm256_storeu_si256(rng->state, oldstate);
	_mm256_zeroupper();
}

void avx_pcg_random_vec8_slow(local_pcg32_random_t* rng,U32PTR rnd, I32 N) {

	__m256i			oldstate	= _mm256_loadu_si256(rng->state );
	const __m256i	SHIFT		= _mm256_set1_epi64x(rng->INCREMENT_4steps);
	const __m256i	MULITPLIER	= _mm256_set1_epi64x(rng->MULTIPLIER_4steps);
 
	#define srl	_mm256_srli_epi64
	#define xor _mm256_xor_si256

	I32 N8 = (N + 7) / 8 * 8;
	for (int i = 0; i < N8; i+=8)	{
		__m256i xorshifted =  srl ( xor(srl(oldstate, 18u), oldstate), 27u) ;
		__m256i rot		   =  srl(oldstate, 59u);	

		oldstate =  _mm256_add_epi64(__mul64_haswell(oldstate, MULITPLIER), SHIFT); //oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;

		__m256i xorshifted1 = srl(xor (srl(oldstate, 18u), oldstate), 27u);
		__m256i rot1		= srl(oldstate, 59u);

		xorshifted1 = _mm256_slli_epi64(xorshifted1, 32);
		rot1		 = _mm256_slli_epi64(rot1, 32);

		//https://stackoverflow.com/questions/18298585/how-to-disable-warning-binary-constants-are-a-gcc-extension
		//
		xorshifted = _mm256_blend_epi32(xorshifted, xorshifted1, 0xAA); // 0b10101010);
		rot		   = _mm256_blend_epi32(rot,		rot1,        0xAA); // 0b10101010);

		__m256i result= _mm256_or_si256(
							_mm256_srlv_epi32(xorshifted, rot),							
							_mm256_sllv_epi32(xorshifted,  _mm256_sub_epi32(_mm256_set1_epi32(32), rot))
						);	//*rnd++ = (xorshifted >> rot) | (xorshifted << ((- *((I32*)(&rot))  ) & 31)); //*rnd++ = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));

		oldstate = _mm256_add_epi64(__mul64_haswell(oldstate, MULITPLIER), SHIFT); //oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;

		if (i < N - 7) {
			_mm256_storeu_si256(rnd + i, result);
		} else {
			int n = N - i; 
			_mm256_maskstore_epi32((I32PTR)rnd + i, GetMoveMask(n), result);
		}

	}
	_mm256_storeu_si256(rng->state, oldstate);	
}



void SetupPCG_AVX2(void){

	 local_pcg_set_seed = avx_pcg_set_seed;
	 //local_pcg_random = avx_pcg_random;
	 local_pcg_random=avx_pcg_random_with_internalbuf;
	 local_pcg_print_state = avx_pcg_print_state;
 
}

#endif


///////////////////////////////////////////////////////////////////////////
#if defined(COMPILER_CLANG) && !defined(cpu_ARM64)
    //pragma clang attribute push (__attribute__((target("avx,avx2"))), apply_to=function)
    #pragma clang attribute pop
#endif
///////////////////////////////////////////////////////////////////////////



#include "abc_000_warning.h"
