
#include <math.h> //sqrtf
#include <stdio.h>
#include <string.h>

#include "abc_000_macro.h"
#include "abc_000_warning.h"
#include "abc_vec.h"

#define STATIC static

#if defined (OS_WIN64) || defined(OS_WIN32)
    #include <malloc.h> //alloca
#else
    #include <alloca.h> //alloca
#endif



#ifdef COMPILER_MSVC
    #define __attribute__(x)
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
	// #pragma GCC target("avx,avx2,avx512f,avx512dq") //https://en.wikipedia.org/wiki/AVX-512
    // -mavx256-split-unaligned-load is spit out by default, hurting new CPUS
    // an altertive is to set -march=hashwell, but it doesntot work in pragram for all gcc versions
    // #pragma GCC optimization ("unroll-loops")
    //#pragma GCC target "arch=core-avx2,tune=core-avx2"
     
    // https://stackoverflow.com/questions/51003218/gcc-target-for-avx2-disabling-sse-instruction-set
    // #pragma GCC target "arch=haswell" : // A buf for some gcc versions, not working
    
#endif
///////////////////////////////////////////////////////////////////////////


/*
// This is an exmaple to illustarate the special processing of the start and teh end of the atray
void  sadd0(float* x, float s, int N) {
    __m256  S = set1(s);
    int i = 0;
    for (; i < N - (NF - 1); i += NF)       _mm256_store_ps(x, add(*((__m256*) (x + i)), S));
    for (; i < N; i++)                      x[i] += s;
}
void  sadd(float* x, float s, int N) {
    __m256  S = set1(s);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        _mm256_store_ps(x+i,    add(*((__m256*) (x + i)),   S));
        _mm256_store_ps(x+i+8,  add(*((__m256*) (x + i+8)), S));
        _mm256_store_ps(x+i+16, add(*((__m256*) (x+i+16)),  S));
        _mm256_store_ps(x+i+24, add(*((__m256*) (x+i+24)),  S));
    }
    for (; i < N - (NF - 1); i += NF) {
        _mm256_store_ps(x + i, add(*((__m256*) (x + i)), S));
    }
    for (; i < N; i++) {    x[i] *= s;   }
    _mm256_zeroupper();
}
float saddu(float* x, float s, int N) {
    __m256  S = set1(s);
    int i = 0;
    for (; i < N && (((intptr_t)&x[i] & 31ULL) != 0); ++i)        x[i] *= s;
    for (; i < N - (NF - 1); i += NF)                     _mm256_store_ps(x, add(*((__m256*) (x + i)), S));
    for (; i < N; i++)                                   x[i] *= s;
}
*/

#define NF          8
#define NF2         (NF*2)
#define NF3         (NF*3)
#define NF4         (NF*4)
#define load        _mm256_loadu_ps
#define store       _mm256_storeu_ps
#define loadi       _mm256_loadu_si256
#define storei      _mm256_storeu_si256
#define maskload    _mm256_maskload_ps
#define maskstore   _mm256_maskstore_ps
#define maskloadi    _mm256_maskload_epi32
#define maskstorei  _mm256_maskstore_epi32
#define set1        _mm256_set1_ps
#define set1_i32    _mm256_set1_epi32
#define set0        _mm256_setzero_ps
#define set0i       _mm256_setzero_si256
#define mul         _mm256_mul_ps
#define add         _mm256_add_ps
#define sub         _mm256_sub_ps
#define addi32      _mm256_add_epi32

#if !defined(COMPILER_SOLARIS) && defined(TARGET_64) && !defined(cpu_ARM64)
#include <immintrin.h>
#include "abc_math_avx.h"

// GCC will impiclty convert all the x+i read/write using aligned load/store while
// MSVC and ICC use the unligned load/store.  So for GCC, we have to expliclty use 
// the unaligned load/store to avoid seg faults. That is,
// addi32(*((__m256i*) (x + i)), C)--> addi32 (loadu(x+i),C)

#if defined(COMPILER_MSVC)

static INLINE __m256i GetMoveMask(int n) {
    __m128i maskIdx = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i maskNum = _mm256_cvtepu8_epi32(maskIdx);
    __m256i nvec    = set1_i32(n);
    __m256i maskmov = _mm256_cmpgt_epi32(nvec, maskNum);
    //__m256i maskmov = _mm256_sub_epi32(maskNum, nvec);    
    return maskmov;
}

static U64   masktemplate[8];
static void  FillMaskTempalte(void) {
    for (int i=0; i<8;i++)     masktemplate[i] = (1ULL << (i * 8)) - 1;
    //{0x00|FF, 0x|FF|FF,0x|FF|FF|FF....}    
}

static INLINE __m256i GetMoveMask1(int n) {
    __m128i maskIdx = _mm_cvtsi64_si128(masktemplate[n]);
    __m256i mask    = _mm256_cvtepi8_epi32(maskIdx);
    return  mask;
}

#else

static U64   masktemplate[8];
static void  FillMaskTempalte(void) {
    for (int i=0; i<8;i++)     masktemplate[i] = (1ULL << (i * 8)) - 1;
    //{0x00|FF, 0x|FF|FF,0x|FF|FF|FF....}
}

#define GetMoveMask(n) ({__m128i maskIdx = _mm_cvtsi64_si128(masktemplate[n]); __m256i mask    = _mm256_cvtepi8_epi32(maskIdx); mask;})

#endif
//https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction
//https://stackoverflow.com/questions/4120681/how-to-calculate-single-vector-dot-product-using-sse-intrinsic-functions-in-c/4121295#4121295


// _mm512_reduce_add_ps: no a single instrution but an emulated sequence of instruction
//https://stackoverflow.com/questions/26896432/horizontal-add-with-m512-avx512
//https://stackoverflow.com/questions/55296777/summing-8-bit-integers-in-m512i-with-avx-intrinsics
static INLINE F32  __attribute__((always_inline))  f32_hsum_slow(__m256 r) { 
    ///horizontal summing 
    //https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction   
    __m128 vlow  = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128  = _mm_add_ps(vlow, vhigh);     // add the low 128
                                               //float hsum_ps_sse3(__m128 v) {
    __m128 sums = _mm_hadd_ps(v128, v128);
    sums = _mm_hadd_ps(sums, sums);
    F32 sum = _mm_cvtss_f32(sums);    
    return sum;
}
 
static INLINE F32  __attribute__((always_inline)) f32_hsum ( __m256 r) {
    //https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction   
    __m128 vlow  = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128  = _mm_add_ps(vlow, vhigh);     // add the low 128

    //float hsum_ps_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v128, shuf);
    shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums = _mm_add_ss(sums, shuf);
    F32 sum = _mm_cvtss_f32(sums);
    return sum;
}

static INLINE I32  __attribute__((always_inline)) i32_hsum(__m256i r)   {
    //https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction   

    __m128i vlow  = _mm256_castsi256_si128(r);
    __m128i vhigh = _mm256_extracti128_si256(r, 1); // high 128
    __m128i v128  = _mm_add_epi32(vlow, vhigh);     // add the low 128

    __m128i hi64 = _mm_unpackhi_epi64(v128, v128);           // 3-operand non-destructive AVX lets us save a byte without needing a mov
   //__m128i hi64 = _mm_shuffle_epi32(x, _MM_SHUFFLE(1, 0, 3, 2));

    __m128i sum64 = _mm_add_epi32(hi64, v128);
    __m128i hi32 = _mm_shufflelo_epi16(sum64, _MM_SHUFFLE(1, 0, 3, 2));    // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);       // SSE2 movd    
}
 

/*
void avx2_f32_memcpy(const I32PTR dst, const I32PTR src, const int N) {    
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(dst + i,       load(src + i)),
        store(dst + i + NF,  load(src + i + NF)),
        store(dst + i + NF2, load(src + i + NF2)),
        store(dst + i + NF3, load(src + i + NF3));
    for (; i < N - (NF - 1); i += NF)
        store(dst + i, load(src + i));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), load(src + i));
    _mm256_zeroupper();
}
*/
STATIC void avx2_i32_add_val_inplace (const int c, const I32PTR x, const int N) {
    __m256i  C = set1_i32(c);
    int i = 0;
    for (; i < N - (NF4-1); i += NF4) 
        storei(x + i,       addi32(loadi(x + i),       C) ),
        storei(x + i + NF,  addi32(loadi(x + i + NF),  C) ),
        storei(x + i + NF2, addi32(loadi(x + i + NF2), C) ),
        storei(x + i + NF3, addi32(loadi(x + i + NF3), C) );
    for (; i < N - (NF-1); i += NF)  
        storei(x + i, addi32(loadi(x + i), C) );
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstorei(x + i,  GetMoveMask(n),  addi32(loadi(x + i), C)  );
    _mm256_zeroupper();    
}

STATIC I32  avx2_i32_sum(const I32PTR x, const int N) {
    __m256i  S = set0i();
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256i s12   = addi32(loadi(x + i),        loadi(x + i + NF));
        __m256i s34   = addi32(loadi(x + i + NF2),  loadi(x + i + NF3));
         S = addi32(S, addi32(s12, s34));
    }
    for (; i < N - (NF - 1); i += NF)  
        S = addi32(S, loadi(x + i));        
    int n = N - i;
    if (n > 0) 
        S= addi32( S,  maskloadi(x + i, GetMoveMask(n) ) );

    I32 sum = i32_hsum(S);
     _mm256_zeroupper();
     return sum;
}

STATIC F32  avx2_f32_sum(const F32PTR x, const int N) {
    __m256  S = set0();
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 s12      = add(load(x + i),        load(x + i + NF) );
        __m256 s34      = add(load(x + i + NF2),  load(x + i + NF3));        
         S = add(S, add(s12, s34));
    }
    for (; i < N - (NF - 1); i += NF)      S = add(S, load(x + i));        
    int n = N - i;
    if (n > 0)        S = add( S,  maskload(x + i, GetMoveMask(n) ) );
    F32 sum = f32_hsum(S);
     _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
     return sum;
}

STATIC void avx2_f32_fill_val(const F32 c, F32PTR x, int N) {
    __m256  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store((x + i),    C);
        store((x + i+NF),  C);
        store((x + i+NF2), C);
        store((x + i+NF3), C);
    }
    for (; i < N - (NF - 1); i += NF)
        store((x+i), C);
    int n = N - i;
    if (n > 0) 
        maskstore(x + i, GetMoveMask(n), C);    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_add_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       //dst=dst-src
        store(dst + i,       add(load(src2 + i),     load(src1 + i)));
        store(dst + i + NF,   add(load(src2 + i + NF), load(src1 + i  + NF)));
        store(dst + i + NF2,  add(load(src2 + i + NF2), load(src1 + i + NF2)));
        store(dst + i + NF3,  add(load(src2 + i + NF3), load(src1 + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) 
        store(dst + i,     add(load(src2 + i), load(src1 + i))  );    
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), add(load(src2 + i), load(src1 + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_sub_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
    //dst=src2-src1
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       //dst=dst-src
        store(dst + i,       sub(load(src2 + i), load(src1 + i)));
        store(dst + i + NF,  sub(load(src2 + i + NF), load(src1 + i + NF)));
        store(dst + i + NF2, sub(load(src2 + i + NF2), load(src1 + i + NF2)));
        store(dst + i + NF3, sub(load(src2 + i + NF3), load(src1 + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) 
        store(dst + i, sub( load(src2 + i), load(src1 + i))  ); 
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), sub(load(src2 + i), load(src1 + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_add_vec_inplace(const F32PTR src, const F32PTR dst, const int N) {
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store(dst + i,           add(load(dst + i),       load(src + i)));
        store(dst + i + NF,      add(load(dst + i + NF),  load(src + i + NF)));
        store(dst + i + NF2,     add(load(dst + i + NF2), load(src + i + NF2)));
        store(dst + i + NF3,     add(load(dst + i + NF3), load(src + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF)  {
        store(dst + i,      add(load(dst + i), load(src + i))  );
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), add(load(dst + i), load(src + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_sub_vec_inplace(const F32PTR src, F32PTR dst, int N) {
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       //dst=dst-src
        store(dst + i,        sub(load(dst + i),      load(src + i)));
        store(dst + i + NF,   sub(load(dst + i + NF),  load(src + i + NF)));
        store(dst + i + NF2,  sub(load(dst + i + NF2), load(src + i + NF2)));
        store(dst + i + NF3,  sub(load(dst + i + NF3), load(src + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) {
        store(dst + i,      sub(load(dst + i), load(src + i))  );
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), sub(load(dst + i), load(src + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_add_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m256  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) 
        store(x + i,          add(load( x + i),     C)   ),
        store(x + i + NF,     add(load(x + i + NF),  C) ),
        store(x + i + NF2,    add(load(x + i + NF2), C) ),
        store(x + i + NF3,    add(load(x + i + NF3), C) );
    for (; i < N - (NF - 1); i += NF)
        store(x + i, add(load(x + i), C) );
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i,  GetMoveMask(n), add(load(x + i), C)  );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_f32_subrev_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m256  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) 
        store(x + i,      sub(C, load( x + i))   ),
        store(x + i + NF,  sub(C, load(x + i + NF)) ),
        store(x + i + NF2, sub(C, load(x + i + NF2)) ),
        store(x + i + NF3, sub(C, load(x + i + NF3)) );
    for (; i < N - (NF - 1); i += NF)
        store(x + i, sub(C, load(x + i)) );
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i,  GetMoveMask(n), sub(C, load(x + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_f32_mul_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m256  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(x + i,      mul(load(x + i), C)),
        store(x + i + NF, mul(load(x + i + NF), C)),
        store(x + i + NF2, mul(load(x + i + NF2), C)),
        store(x + i + NF3, mul(load(x + i + NF3), C));
    for (; i < N - (NF - 1); i += NF)
        store(x + i, mul(load(x + i), C));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i, GetMoveMask(n), mul(load(x + i), C));
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_f32_mul_vec_inplace(const F32PTR src, F32PTR dst, int N) {
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store(dst + i,       mul(load(dst + i),      load(src + i)));
        store(dst + i + NF,   mul(load(dst + i + NF),  load(src + i + NF)));
        store(dst + i + NF2,  mul(load(dst + i + NF2), load(src + i + NF2)));
        store(dst + i + NF3,  mul(load(dst + i + NF3), load(src + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) {
        store(dst + i, mul(load(dst + i), load(src + i))  );
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), mul(load(dst + i), load(src + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_mul_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       
        store(dst + i,       mul(load(src2 + i),      load(src1 + i)));
        store(dst + i + NF,  mul(load(src2 + i + NF),  load(src1 + i + NF)));
        store(dst + i + NF2,  mul(load(src2 + i + NF2), load(src1 + i + NF2)));
        store(dst + i + NF3,  mul(load(src2 + i + NF3), load(src1 + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) {
        store(dst + i,     mul(load(src2 + i), load(src1 + i))  );
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), mul(load(src2 + i), load(src1 + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}

STATIC void avx2_f32_add_v_v2_vec_inplace(const F32PTR src, const F32PTR  v, const F32PTR  v2, int N) {

    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store(v + i,  add(load(v + i), load(src + i)));
        store(v2 + i, add(load(v2 + i), mul(load(src + i), load(src + i))));

        store(v + i+NF,  add(load(v + i + NF),  load(src + i + NF)));
        store(v2 +i+ NF, add(load(v2 + i + NF), mul(load(src+i + NF), load(src + i + NF))));

        store(v + i + NF2,  add(load(v + i + NF2), load(src + i + NF2)) );
        store(v2 + i + NF2, add(load(v2 + i + NF2), mul(load(src + i + NF2), load(src + i + NF2))));

        store(v +  i+NF3,  add(load(v + i + NF3), load(src + i + NF3)));
        store(v2 + i+ NF3, add(load(v2 + i + NF3), mul(load(src+i+ NF3), load(src + i + NF3) )));
    }
    for (; i < N - (NF - 1); i += NF) {
        store(v  + i, add(load(v + i), load(src + i)));
        store(v2 + i, add(load(v2 + i), mul(load(src + i), load(src + i))));
    }

    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __m256i mask = GetMoveMask(n);
        maskstore(v + i, mask, add(    load(v + i),   load(src + i)  )                           );
        maskstore(v2 + i, mask, add(   load(v2 + i), mul(load(src + i), load(src + i)))     );
    }
    _mm256_zeroupper();

 
}

STATIC F32 avx2_f32_dot(F32PTR  x, F32PTR y, int N) {

    __m256 r = set0();    
    int i = 0;    
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 r0 = mul(load(x + i),     load(y + i)) ;
        __m256 r1 = mul(load(x + i+NF),   load(y + i+ NF));
        __m256 r2 = mul(load(x + i+ NF2),  load(y + i + NF2));
        __m256 r3 = mul(load(x + i + NF3), load(y + i + NF3));
        __m256 s1   = add(r0, r1);
        __m256 s2   = add(r2, r3);
        __m256 s12  = add(s1, s2);
        r=add(r, s12);
    }        
    for (; i < N - (NF - 1); i += NF) 
         r =add(r,  mul(load(x + i), load(y + i))  ); 
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        r = add(r,  mul(   maskload(x + i, mask), maskload(y + i, mask)   )   );
    }    
    F32 sum = f32_hsum(r);
     _mm256_zeroupper();    
    return sum;
}

STATIC F32 avx2_f32_dot2x1(F32PTR  x, F32PTR y, F32PTR v, int N, F32PTR res) {
    __m256 RX = set0();
    __m256 RY = set0();
    int i = 0;    
    for (; i < N - (NF2-1); i += NF2) {
        __m256 x0 = mul(load(x + i),      load(v + i)) ;
        __m256 x1 = mul(load(x + i+NF),   load(v + i+ NF));
        __m256 y0 = mul(load(y + i),      load(v + i)) ;
        __m256 y1 = mul(load(y + i+ NF),  load(v + i+NF));
   
        RX =add(RX, add(x0, x1)); 
        RY =add(RY, add(y0, y1));
    }        
    for (; i < N - (NF - 1); i += NF) {
         RX =add(RX,  mul(load(x + i), load(v + i))  );
         RY = add(RY, mul(load(y + i), load(v + i))  );
    }

    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        __m256  vvec = maskload(v + i, mask); 
        RX = add(RX, mul(maskload(x + i, mask), vvec) );        
        RY = add(RY, mul(maskload(y + i, mask), vvec) );
    }

    float sumX = f32_hsum(RX);
    float sumY = f32_hsum(RY);
   
     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    res[0] = sumX;
    return sumY;
}

STATIC void avx2_f32_dot2x2(F32PTR  x, F32PTR y, F32PTR v, F32PTR w, int N, F32PTR res1, F32PTR res2) {
    __m256 RX1 = set0();
    __m256 RX2 = set0();
    __m256 RY1 = set0();
    __m256 RY2 = set0();
    int i = 0;    
    
    for (; i  < N - (NF2-1); i += NF2) {
        __m256 x0 = mul( load(x+i),      load(v+i)    ) ;
        __m256 x1 = mul( load(x+i+NF),   load(v+i+NF))  ;
        RX1 = add(RX1, add(x0, x1));

        __m256 y0 = mul( load(y+i),      load(v+i)    ) ;
        __m256 y1 = mul( load(y+i+NF),   load(v+i+NF) );
        RY1  = add(RY1, add(y0, y1));

        __m256 xx0 = mul( load(x+i),       load(w+i)     );
        __m256 xx1 = mul( load(x+i+NF),    load(w+i+NF)  );
        RX2 = add(RX2, add(xx0, xx1));

        __m256 yy0 = mul( load(y + i),      load(w + i)     );
        __m256 yy1 = mul( load(y + i + NF), load(w + i + NF));
        RY2 = add(RY2, add(yy0, yy1));
    }      
    for (; i < N - (NF - 1); i += NF) {
         RX1 = add(RX1,  mul(load(x + i), load(v + i))  );
         RY1 = add(RY1, mul(load(y + i), load(v + i))  );

         RX2 = add(RX2, mul(load(x + i), load(w + i)));
         RY2 = add(RY2, mul(load(y + i), load(w + i)));
    }
    
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        __m256  vvec = maskload(v + i, mask);

        __m256  r0   =  mul(maskload(x + i, mask), vvec);       RX1 = add(RX1, r0);
        __m256  r1 = mul(maskload(y + i, mask), vvec);          RY1 = add(RY1, r1);

        __m256  wvec = maskload(w + i, mask);
        __m256  rr0 = mul(maskload(x + i, mask), wvec);        RX2 = add(RX2, rr0);
        __m256  rr1 = mul(maskload(y + i, mask), wvec);        RY2 = add(RY2, rr1);
    }

    float sumX1 = f32_hsum(RX1);
    float sumY1 = f32_hsum(RY1);
    float sumX2 = f32_hsum(RX2);
    float sumY2 = f32_hsum(RY2);


     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    res1[0] = sumX1;
    res1[1] = sumY1;
    res2[0] = sumX2;
    res2[1] = sumY2;
    return;
}

STATIC F32 avx2_fma_f32_dot(F32PTR  x, F32PTR y, int N) {

    __m256 r0 = set0();
    __m256 r1 = set0();
    __m256 r2 = set0();
    __m256 r3 = set0();
    
    int i = 0;    
    for (; i < N - (NF4 - 1); i += NF4) {
        r0 = _mm256_fmadd_ps(load(x + i),   load(y + i),   r0) ;
        r1 = _mm256_fmadd_ps(load(x + i+8), load(y + i+8), r1);
        r2 = _mm256_fmadd_ps(load(x + i+16), load(y + i+16),r2);
        r3 = _mm256_fmadd_ps(load(x + i+24), load(y + i+24), r3);        
    }        
    for (; i < N - (NF - 1); i += NF) 
         r0 = _mm256_fmadd_ps(load(x + i), load(y + i),  r0  );    
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        r1   = _mm256_fmadd_ps( maskload(x + i, mask), maskload(y + i, mask),r1);
    }
    __m256 s1 = add(r0, r1);
    __m256 s2 = add(r2, r3);
    __m256 r  = add(s1, s2);
    float sum = f32_hsum(r);                      
     _mm256_zeroupper();
    return sum;
}

//https://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
//https://stackoverflow.com/questions/54642663/how-sqrt-of-gcc-works-after-compiled-which-method-of-root-is-used-newton-rap
//https://stackoverflow.com/questions/31555260/fast-vectorized-rsqrt-and-reciprocal-with-sse-avx-depending-on-precision
STATIC void avx2_f32_sqrt_vec_inplace(const F32PTR x, const int N)
{
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 r1 = load(x + i);          store(x + i,      _mm256_sqrt_ps(r1)  );
        __m256 r2 = load(x + i+NF);       store(x + i+ NF,  _mm256_sqrt_ps(r2) );
        __m256 r3 = load(x + i + NF2);    store(x + i+ NF2, _mm256_sqrt_ps(r3) );
        __m256 r4 = load(x + i + NF3);    store(x + i+ NF3, _mm256_sqrt_ps(r4) );
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);           store(x + i,       _mm256_sqrt_ps(r)  );
    }
    
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);
        maskstore(x + i, GetMoveMask(n), _mm256_sqrt_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_sqrt_vec(const F32PTR x, const F32PTR y, const int N)
{
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 r1 = load(x + i);          store(y + i,      _mm256_sqrt_ps(r1)  );
        __m256 r2 = load(x + i+NF);       store(y + i+ NF,  _mm256_sqrt_ps(r2) );
        __m256 r3 = load(x + i + NF2);    store(y + i+ NF2, _mm256_sqrt_ps(r3) );
        __m256 r4 = load(x + i + NF3);    store(y + i+ NF3, _mm256_sqrt_ps(r4) );     
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);    store(y + i, _mm256_sqrt_ps(r));
    }
    
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);
        maskstore(y + i, GetMoveMask(n), _mm256_sqrt_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

#ifdef COMPILER_MSVC
STATIC void avx2_f32_sin_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);    store(x + i,   _mm256_sin_ps(r));
    }    
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm256_sin_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_cos_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);    store(x + i,   _mm256_cos_ps(r));
    }    
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm256_cos_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_sincos_vec_inplace_MSVC(const F32PTR in_outsin, F32PTR outcos, const int N)
{
//Compute the sine and cosine of packed single-precision (32-bit) floating-point elements in a expressed in radians,
//store the sine in dst, and store the cosine into memory at mem_addr.
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(in_outsin + i);
        store(in_outsin + i, _mm256_sincos_ps(outcos+i,r));
    }
    int n = N - i;
    if (n > 0) {
        __m256i mask   = GetMoveMask(n);
        __m256  tmpcos;
        __m256  r = load(in_outsin + i);   
        maskstore(in_outsin + i, mask, _mm256_sincos_ps(&tmpcos, r) );
        maskstore(outcos  + i,   mask, tmpcos);
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_pow_vec_inplace_MSVC(F32PTR x, F32 pow, int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    
    __m256 C = set1(pow);

    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);  store(x + i,  _mm256_pow_ps(r,C) );
    }    
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);  maskstore(x + i, GetMoveMask(n), _mm256_pow_ps(r, C));
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_log_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);    store(x + i, _mm256_log_ps(r));
    }
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm256_log_ps(r));
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}

STATIC void avx2_f32_exp_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 r = load(x + i);    store(x + i, _mm256_exp_ps(r));
    }
    int n = N - i;
    if (n > 0) {
        __m256 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm256_exp_ps(r));
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
#endif

STATIC void avx2_f32_sin_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i,  sin256_ps(x + i) );    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), sin256_ps(x + i) );    
    _mm256_zeroupper();    
}

STATIC void avx2_f32_cos_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, cos256_ps(x + i));    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), cos256_ps(x + i) );    
    _mm256_zeroupper();    
}

STATIC void avx2_f32_sincos_vec_inplace(const F32PTR in_outsin, F32PTR outcos, const int N)
{
//Compute the sine and cosine of packed single-precision (32-bit) floating-point elements in a expressed in radians,
//store the sine in dst, and store the cosine into memory at mem_addr.
    int i = 0;
    for (; i < N - (NF - 1); i += NF) 
        sincos256_ps(in_outsin + i, in_outsin + i, outcos + i);    
    int n = N - i;
    if (n > 0) {
        __m256i mask   = GetMoveMask(n);
        __m256  tmpcos, tmpsin;
        sincos256_ps(in_outsin + i, &tmpsin, &tmpcos);
        maskstore(in_outsin + i,   mask,  tmpsin);
        maskstore(outcos    + i,   mask,  tmpcos);
    }
    _mm256_zeroupper();    
}

static INLINE __m256 pow256(F32PTR x, float n) {
    //https://community.intel.com/t5/Intel-ISA-Extensions/AVX512-reciprocal-approximations/td-p/1068416
    __m256 res;
    res = _mm256_mul_ps(log256_ps(x), _mm256_set1_ps(n));
    res = exp256_ps(&res);
    return res;
}
static INLINE __m256 pow256_int(__m256 x, int n)
{
    __m256 res = _mm256_set1_ps(1.0f);
    int npositive = n >= 0 ? n : -n;
    while (npositive) {
        if (npositive & 1L)
            res = _mm256_mul_ps(res, x);
        x = _mm256_mul_ps(x, x);
        npositive = npositive >> 1;
    }
    if (n < 0) {
        res = _mm256_rcp_ps(res);
    }
    return res;
}

STATIC void avx2_f32_pow_vec_inplace(F32PTR x, F32 pow, int N) {
         int   nInteger   = pow;
        float  nRemainder = pow - nInteger;

        if (nRemainder != 0) {
            int i = 0;
            for (; i < N - (NF - 1); i += NF) {
                store(x + i, pow256(x+i, pow ) );
            }            
            int n = N - i;
            if (n > 0) 
                maskstore(x + i, GetMoveMask(n), pow256(x + i, pow) );              
        }
        else {
            int i = 0;
            for (; i < N - (NF - 1); i += NF) {
                store(x + i, pow256_int(load(x + i), nInteger));
            }
            int n = N - i;
            if (n > 0)
                maskstore(x + i, GetMoveMask(n), pow256_int(load(x + i), nInteger));         
        }
        _mm256_zeroupper();
}

STATIC void avx2_f32_log_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, log256_ps(x + i) );    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), log256_ps(x + i) );    
    _mm256_zeroupper();    
}

STATIC void avx2_f32_exp_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, exp256_ps(x + i) );    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), exp256_ps(x + i) );
    _mm256_zeroupper();    
}

STATIC void  avx2_f32_avgstd(const F32PTR x, int N, F32PTR avg, F32PTR std) {
    __m256  S    = set0();
    __m256  SS   = set0();
    int i = 0;
    for (; i < N - 15; i += 16) {
        __m256 s1  =  load(x + i);
        __m256 s2  =  load(x + i+8);

        __m256 s12  =  add(s1, s2);                      
        __m256 ss12 = add(mul(s1, s1), mul(s2, s2));

         S  = add(S,  s12);
         SS = add(SS, ss12);
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 s1 = load(x + i);
        S  = add(S,  s1);
        SS = add(SS, mul(s1, s1));
    }
    
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        __m256  s1  =  maskload(x + i,mask);
        S   = add(S, s1);
        SS  = add(SS, mul(s1, s1));
    }
    
    F32 sumx = f32_hsum(S);
    F32 sumxx = f32_hsum(SS);
     _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
     F64 AVG = sumx / N;
     std[0] = sqrtf((sumxx - AVG * sumx) / (N - 1));// sqrtf((sumxx - AVG * AVG * N) / (N - 1));
     avg[0] =AVG;
}

STATIC void avx2_f32_sx_sxx_to_avgstd_inplace(F32PTR SX, F32PTR SXX, I32 Nsample, F32 scale, F32 offset, int N) {
    /*
                   r_ippsMulC_32f_I(inv_sample* sd, resultChain.tY, N);
                   r_ippsMul_32f(resultChain.tY, resultChain.tY, MEMBUF1, N);
                   r_ippsMulC_32f_I(inv_sample* sd* sd, resultChain.tSD, N);
                   r_ippsSub_32f_I(MEMBUF1, resultChain.tSD, N);
                   r_ippsSqrt_32f_I(resultChain.tSD, N);
    */

    __m256 invsample_scale  = set1(1. / (F64)Nsample * scale);
    __m256 invsample_scale2 = set1(1. / (F64)Nsample * scale*scale);
    __m256 offset_vec       = set1(offset);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 avg   = mul(load(SX + i),    invsample_scale);
        __m256 sxx_n = mul(load(SXX + i),   invsample_scale2);

        __m256 sd2  = sub(sxx_n, mul(avg, avg));
        __m256 sd   = mul(_mm256_rsqrt_ps(sd2), sd2); //sqrt(sd2)

        avg = add(avg, offset_vec);
        store(SXX + i, sd);
        store(SX  + i, avg); 

    }

    int n = N - i;
    if (n >0){
        __m256i mask  = GetMoveMask(n);
        __m256  avg   = mul(maskload(SX + i,mask), invsample_scale);
        __m256  sxx_n = mul(maskload(SXX + i,mask), invsample_scale2);

        __m256 sd2 = sub(sxx_n, mul(avg, avg));
        __m256 sd = mul(_mm256_rsqrt_ps(sd2), sd2); //sqrt(sd2)

        avg = add(avg, offset_vec);
        maskstore(SXX + i, mask,sd);
        maskstore(SX + i, mask, avg);

    }

}


#ifdef OS_WIN64 
    #define WIN32_LEAN_AND_MEAN
    #include "windows.h"

STATIC I32 avx2_f32_maxidx1(const F32PTR x, const int N, F32PTR val) {
//https://stackoverflow.com/questions/23590610/find-index-of-maximum-element-in-x86-simd-vector    
    if (N < 8) {        
        int maxIdx = 0;
        F32 maxVal = x[0];        
        for (int i=1; i < N; i++) {
            if (maxVal < x[i])  maxVal = x[i], maxIdx = i;            
        }
        val[0] = maxVal;
        return maxIdx;
    }


    //__m256 r1 = set0();    
    __m256  maxVec      = load(x);
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i maxIdx      = _mm256_cvtepu8_epi32(_maskIdx);
    __m256i EIGHT       = set1_i32(8);

    __m256i idx         = maxIdx;

    int i = 8;    
    for (; i < N - (NF - 1); i += NF) {        
        idx   = addi32(idx, EIGHT);

        __m256  vec     = load(x + i);
        __m256  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);

         maxVec  = _mm256_blendv_ps(maxVec,  vec, cmpmask); //0: take a, 1: take from vec        
         maxIdx  = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));        
    }        
 
    int n = N - i;
    if (n > 0) {
        __m256i mask    = GetMoveMask(n);
        __m256  vec     = maskload(x + i,mask);
        idx = addi32(idx, EIGHT);
        
        idx = _mm256_blendv_epi8(maxIdx,idx, mask);
        vec = _mm256_blendv_ps(  maxVec,vec, _mm256_castsi256_ps(mask)); //o: take a, 1: take from vec        

        __m256  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);
        maxVec = _mm256_blendv_ps(maxVec, vec, cmpmask); //o: take a, 1: take from vec        
        maxIdx = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));
    }

    __m128 vlow     = _mm256_castps256_ps128(maxVec);
    __m128 vhigh    = _mm256_extractf128_ps(maxVec, 1); // high 128    
    __m128 tmp      = _mm_max_ps(vlow, vhigh);
    tmp = _mm_max_ps(tmp, _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(1, 0, 3, 2) ) );
    tmp = _mm_max_ps(tmp, _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 1))  );    
    F32 max = _mm_cvtss_f32(tmp);
    __m256  max256 = set1(max);

    __m256   vcmp      = _mm256_cmp_ps(maxVec, max256, _CMP_EQ_OQ);
    uint32_t finalmask = _mm256_movemask_epi8(_mm256_castps_si256(vcmp));
 
    DWORD index;
    _BitScanReverse(&index, finalmask);     //tmp = _mm_packs_epi16(_mm_packs_epi32(_mm_cmpeq_epi32(_mm_set1_epi32(max), lo),    _mm_cmpeq_epi32(_mm_set1_epi32(max), hi)),       _mm_setzero_si128());
    union {
        I32    idx[8];
        __m256i dummy;
    }ID;
    _mm256_store_si256(&ID.dummy, maxIdx);
    _mm256_zeroupper();     
    val[0] = max;
    return ID.idx[index>>2] ;
}
#endif

STATIC I32 avx2_f32_maxidx(const F32PTR x, const int N, F32PTR val) {
       
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i idx         = _mm256_cvtepu8_epi32(_maskIdx);
    __m256i EIGHT       = set1_i32(8);

    __m256  maxVec = set1(x[0]);
    __m256i maxIdx = set0i();
    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {        
        __m256  vec     = load(x + i);
        __m256  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);

         maxVec  = _mm256_blendv_ps(maxVec,  vec, cmpmask); //o: take a, 1: take from vec        
         maxIdx  = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));        
         idx = addi32(idx, EIGHT);
    }         
    int n = N - i;
    if (n > 0) {
        __m256i mask    = GetMoveMask(n);
        __m256  vec     = maskload(x + i,mask);
        //idx = addi32(idx, EIGHT);
        
        vec = _mm256_blendv_ps(maxVec, vec, _mm256_castsi256_ps(mask)); //o: take a, 1: take from vec        
        idx = _mm256_blendv_epi8(maxIdx,idx, mask);
        
        __m256  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);
        maxVec = _mm256_blendv_ps(maxVec, vec, cmpmask); //o: take a, 1: take from vec        
        maxIdx = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));
    }
    __m128 vlow     = _mm256_castps256_ps128(maxVec);
    __m128 vhigh    = _mm256_extractf128_ps(maxVec, 1); // high 128
    __m128i idxlow = _mm256_castsi256_si128(maxIdx);
    __m128i idxhigh = _mm256_extracti128_si256(maxIdx, 1); // high 128

    __m128 cmpmask  = _mm_cmp_ps(vlow, vhigh, _CMP_LE_OQ);  
    __m128 max128   = _mm_blendv_ps(vlow, vhigh, cmpmask); //o: take a, 1: take from vec     
    __m128i idx128  = _mm_blendv_epi8(idxlow, idxhigh, _mm_castps_si128(cmpmask)); //o: take a, 1: take from vec     

     F32 VAL[4];
     union {
        I32 idx[4];
        __m128i dummy;
     }ID;
     _mm_storeu_ps(VAL,   max128);
     _mm_store_si128(&ID.dummy, idx128);
     _mm256_zeroupper();
     I32 i1 = VAL[0] > VAL[1] ? 0 : 1;
     I32 i2 = VAL[2] > VAL[3] ? 2 : 3;
     I32 IDX = VAL[i1] > VAL[i2] ? i1 : i2;

     val[0] = VAL[IDX];
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return ID.idx[IDX];
}

STATIC I32 avx2_f32_minidx(const F32PTR x, const int N, F32PTR val) {
    
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i idx         = _mm256_cvtepu8_epi32(_maskIdx);
    __m256i EIGHT       = set1_i32(8);

    __m256   minVec = set1(x[0]);
    __m256i  minIdx = set0i();

    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {        
        __m256  vec     = load(x + i);
        __m256  cmpmask = _mm256_cmp_ps(minVec, vec, _CMP_GE_OQ);

        minVec = _mm256_blendv_ps(minVec,  vec, cmpmask); //o: take a, 1: take from vec        
        minIdx = _mm256_blendv_epi8(minIdx, idx, _mm256_castps_si256(cmpmask));

        idx = addi32(idx, EIGHT);
    }        
 
    int n = N - i;
    if (n > 0) {
        __m256i mask    = GetMoveMask(n);
        __m256  vec     = maskload(x + i,mask);
        //idx = addi32(idx, EIGHT);
        
        vec = _mm256_blendv_ps(minVec, vec, _mm256_castsi256_ps(mask)); //o: take a, 1: take from vec  
        idx = _mm256_blendv_epi8(minIdx,idx, mask);
             
        __m256  cmpmask = _mm256_cmp_ps(minVec, vec, _CMP_GE_OQ);
        minVec = _mm256_blendv_ps(minVec, vec, cmpmask); //o: take a, 1: take from vec        
        minIdx = _mm256_blendv_epi8(minIdx, idx, _mm256_castps_si256(cmpmask));
    }
    __m128 vlow     = _mm256_castps256_ps128(minVec);
    __m128 vhigh    = _mm256_extractf128_ps(minVec, 1); // high 128
    __m128 cmpmask  = _mm_cmp_ps(vlow, vhigh, _CMP_GE_OQ);
  
    __m128 max128 = _mm_blendv_ps(vlow, vhigh, cmpmask); //o: take a, 1: take from vec     

    __m128i idxlow  = _mm256_castsi256_si128(minIdx);
    __m128i idxhigh = _mm256_extracti128_si256(minIdx, 1); // high 128
    __m128i idx128  = _mm_blendv_epi8(idxlow, idxhigh, _mm_castps_si128(cmpmask)); //o: take a, 1: take from vec     

     F32 VAL[4];
     union {
        I32 idx[4];
        __m128i dummy;
     }ID;
     _mm_storeu_ps(VAL,   max128);
     _mm_store_si128(&ID.dummy, idx128);
     _mm256_zeroupper();
     I32 i1 = VAL[0] < VAL[1] ? 0 : 1;
     I32 i2 = VAL[2] < VAL[3] ? 2 : 3;
     I32 IDX = VAL[i1] < VAL[i2] ? i1 : i2;

     val[0] = VAL[IDX];
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return ID.idx[IDX];
}

STATIC void avx2_f32_diff_back(const F32PTR  x, F32PTR result, const int N) {
    //skip the first elements because evaluation of diff0-x[0]-x[-1] needs an nonexistent element at -1
    //res[i]=x[i]-x[i-1]

    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store(result + i,        sub(load(x + i),         load(x + i-1)));
        store(result + i + NF,   sub(load(x + i + NF),    load(x + i + NF-1)));
        store(result + i + NF2,  sub(load(x + i + NF2),    load(x + i + NF2 -1)));
        store(result + i + NF3,  sub(load(x + i + NF3),    load(x + i + NF3 -1)));
    }
    for (; i < N - (NF - 1); i += NF) {
        store(result + i, sub(load(x + i), load(x + i-1))  );
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(result + i, GetMoveMask(n), sub(load(x + i), load(x + i-1))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];

    result[0] = result[1];

}

STATIC void avx2_f32_seq(F32PTR p, F32 x0, F32 dX, int N) {
   
    __m256 X;
    {    
        __m128i tmp     = _mm_cvtsi64_si128(0x0706050403020100);
        __m256i tmp256  = _mm256_cvtepu8_epi32(tmp);
        __m256  v1to7   = _mm256_cvtepi32_ps(tmp256);        
        __m256  DX      = mul(v1to7, set1(dX));
          
        X = add(set1(x0),DX);    
    }
    __m256 DX8 = set1(dX*8);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        store(p + i, X);
        X = add(X, DX8);        
    } 
    int n = N - i;
    if (n > 0) {
        maskstore(p + i, GetMoveMask(n), X);
    }     
    _mm256_zeroupper();    
}

STATIC void avx2_i32_seq(I32PTR p, I32 x0, I32 dX, int N) {

    __m256i X;
    {
        __m128i tmp    = _mm_cvtsi64_si128(0x0706050403020100);
        __m256i v1to7  = _mm256_cvtepu8_epi32(tmp);
        //__m256  v1to7  = _mm256_cvtepi32_ps(tmp256);
        __m256i  DX     = _mm256_mullo_epi32(v1to7, set1_i32(dX)); //https://stackoverflow.com/questions/28479429/multiply-two-vectors-of-32bit-integers-producing-a-vector-of-32bit-result-eleme

        X = addi32(set1_i32(x0), DX);
    }
    __m256i DX8 = set1_i32(dX * 8);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        storei(p + i, X);
        X = addi32(X, DX8);
    }
    int n = N - i;
    if (n > 0) {
        maskstorei(p + i, GetMoveMask(n), X);
    }
    _mm256_zeroupper();
}

STATIC void avx2_f32_to_f64_inplace(F32PTR data32, int N) {

    F64PTR data64 = data32;
    int i = N - 8;
    for (; i >= 0; i -= 8) {
        __m256 r = load(data32 + i);
        __m128 lowps    = _mm256_castps256_ps128(r);
        __m128 highps   = _mm256_extractf128_ps(r, 1);         // high 128
        _mm256_storeu_pd(data64 + i, _mm256_cvtps_pd(lowps));
        _mm256_storeu_pd(data64 + i+4, _mm256_cvtps_pd(highps));
    }
    i = i + 8;
    if (i >= 4) {
        i = i - 4;
        __m128 r = _mm_loadu_ps(data32 + i);        
        _mm256_storeu_pd(data64 + i, _mm256_cvtps_pd(r));
    }
    i = i - 1;
    for (; i >= 0; --i)   data64[i] = data32[i];    
    _mm256_zeroupper();    
}

STATIC void avx2_f64_to_f32_inplace(F64PTR data64, int N) {

    F32PTR data32 = data64;
    int i = 0;
    for (; i < N-7; i += 8) {    
        __m128  R1 = _mm256_cvtpd_ps(_mm256_loadu_pd(data64 + i));
        __m128  R2 = _mm256_cvtpd_ps(_mm256_loadu_pd(data64 + i + 4));
        _mm_storeu_ps(data32 + i,   R1);
        _mm_storeu_ps(data32 + i+4, R2);
    }
    for (; i < N - 3; i +=4) {
        __m128  R1 = _mm256_cvtpd_ps(_mm256_loadu_pd(data64 + i));
        _mm_storeu_ps(data32 + i, R1);
    }    
    _mm256_zeroupper();

    for (; i <N; ++i)      data32[i] = data64[i];         
}

STATIC void avx2_i32_to_f32_scaleby_inplace(I32PTR x, int N, F32 scale) {
    __m256  C = set1(scale);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(x + i,      mul(_mm256_cvtepi32_ps(loadi(x + i)), C)),
        store(x + i + NF, mul(_mm256_cvtepi32_ps(loadi(x + i+NF)), C)),
        store(x + i + NF2, mul(_mm256_cvtepi32_ps(loadi(x + i+ NF2)), C)),
        store(x + i + NF3, mul(_mm256_cvtepi32_ps(loadi(x + i+ NF3)), C));
    for (; i < N - (NF - 1); i += NF)
        store(x + i, mul(_mm256_cvtepi32_ps(loadi(x + i)), C));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i, GetMoveMask(n), mul(_mm256_cvtepi32_ps(loadi(x + i)), C));
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_i32_increment_bycond_inplace(I32PTR x,  F32PTR cond, int N) {
    __m256   C0 = set0();
    __m256i  C1 = set1_i32(1);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {            
        __m256i mask = _mm256_castps_si256(  _mm256_cmp_ps( load(cond+i), C0, _CMP_GT_OQ) );
        __m256i vec = maskloadi(x + i, mask);
        vec = addi32(vec, C1);
        maskstorei(x + i, mask, vec);
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __m256i loadmask = GetMoveMask(n);
        __m256i mask = _mm256_castps_si256(_mm256_cmp_ps(maskload(cond + i,loadmask), C0, _CMP_GT_OQ));
        mask = _mm256_and_si256(mask, loadmask);
        __m256i vec = maskloadi(x + i, mask);
        vec = addi32(vec, C1);
        maskstorei(x + i, mask, vec);

    }
    _mm256_zeroupper();
}

STATIC void avx2_i32_increment_vec2_bycond_inplace(I32PTR x, I32PTR y, F32PTR cond, int N) {
    __m256   C0       = set0();
    __m256   Ceplison1 = set1(1e-10);
    __m256   Ceplison2 = set1(-1e-10);
    __m256i  C1       = set1_i32(1);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {         
        __m256  data  = load(cond + i);
        __m256i mask1 = _mm256_castps_si256(  _mm256_cmp_ps(data, Ceplison1, _CMP_GT_OQ) );
        __m256i vec1   = maskloadi(x + i, mask1);
        vec1 = addi32(vec1, C1);     maskstorei(x + i, mask1, vec1);

        __m256i mask2= _mm256_castps_si256( _mm256_and_ps( _mm256_cmp_ps(Ceplison1, data, _CMP_GT_OQ), _mm256_cmp_ps(data, Ceplison2, _CMP_GT_OQ)) );
        __m256i vec2 = maskloadi(y + i, mask2);
        vec2 = addi32(vec2, C1);     maskstorei(y + i, mask2, vec2);
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __m256i loadmask = GetMoveMask(n);

        __m256  data  = maskload(cond + i, loadmask);
        __m256i mask1 = _mm256_castps_si256(  _mm256_cmp_ps(data, Ceplison1, _CMP_GT_OQ) );
        mask1 = _mm256_and_si256(mask1, loadmask);
        __m256i vec1   = maskloadi(x + i, mask1);
        vec1 = addi32(vec1, C1);     maskstorei(x + i, mask1, vec1);

        __m256i mask2= _mm256_castps_si256( _mm256_and_ps( _mm256_cmp_ps(Ceplison1, data, _CMP_GT_OQ), _mm256_cmp_ps(data, Ceplison2, _CMP_GT_OQ)) );
        mask2 = _mm256_and_si256(mask2, loadmask);
        __m256i vec2 = maskloadi(y + i, mask2);
        vec2 = addi32(vec2, C1);     maskstorei(y + i, mask2, vec2);        

    }
    _mm256_zeroupper();
}

STATIC I32  avx2_i08_sum_binvec(U08PTR binvec, I32 N) {
    //stackoverflow.com/questions/36998538/fastest-way-to-horizontally-sum-sse-unsigned-byte-vector
	I32   SUM = 0;
	I32	  i   = 0;
	//|8bytees|8bytees|8bytees|8bytees|
    __m256i r   = set0i();
    I32  ite_read_for_sad = 0;
	for (; i < N - (32*4-1); i+=32*4) {
        __m256i sum1  = _mm256_adds_epu16(loadi(binvec + i),    loadi(binvec + i + 32));
        __m256i sum2  = _mm256_adds_epu16(loadi(binvec + i+64), loadi(binvec + i + 96));      
        __m256i sum12 = _mm256_adds_epu16(sum1, sum2);
                 r    = _mm256_adds_epu16(r, sum12);                 
                ++ite_read_for_sad;
                if (ite_read_for_sad > 60) {
                    ite_read_for_sad = 0;
                    //stackoverflow.com/questions/36998538/fastest-way-to-horizontally-sum-sse-unsigned-byte-vector
                    __m256i  sum= _mm256_sad_epu8(r, set0i());  

                    __m128i sumLow  = _mm256_castsi256_si128(sum);
                    __m128i sumHigh = _mm256_extracti128_si256(sum, 1); // high 128

                    __m128i sum128 = _mm_adds_epu16(sumLow, sumHigh);
                    I32    localSum = _mm_extract_epi16(sum128, 0) + _mm_extract_epi16(sum128, 4);
                    SUM = SUM + localSum;
                    r = set0i();
                }
	}
    for (; i < N - (32 - 1); i += 32) {
            r = _mm256_adds_epu16(r, loadi(binvec + i));
    }   
    {
        r = _mm256_sad_epu8(r, set0i());
        __m128i sumLow = _mm256_castsi256_si128(r);
        __m128i sumHigh = _mm256_extracti128_si256(r, 1); // high 128
        __m128i sum = _mm_adds_epu16(sumLow, sumHigh);
        I32    localSum = _mm_extract_epi16(sum, 0) + _mm_extract_epi16(sum, 4);
        SUM = SUM + localSum;
    }
    
    I32  SUM_remaining = 0;
    binvec += i;
	for (; i < N - (8 - 1); i += 8) {
        SUM_remaining   += *((I32PTR)binvec) + *((I32PTR)binvec +1);		
        binvec += 8;
    }
	for (; i < N-1; i+=2 )  {
        SUM_remaining   += *((I16PTR)binvec);   
        binvec += 2;
    }
    if (i < N) 
        SUM_remaining += *binvec;    
    if(SUM_remaining){
        *((I16PTR)&SUM_remaining) = *((I16PTR)&SUM_remaining) + *((I16PTR)&SUM_remaining + 1);
        *((I08PTR)&SUM_remaining) = *((I08PTR)&SUM_remaining) + *((I08PTR)&SUM_remaining + 1);
        SUM = SUM + *((I08PTR)&SUM_remaining);
    }

	return SUM;
}



/*
SGEMM  performs one of the matrix - matrix operations
* > C : = alpha * op(A) * op(B) + beta * C,
* > where  op(X) is one of
* > op(X) = X or op(X) = X * *T,
* > alpha and beta are scalars, and A, Band C are matrices, with op(A)
* > an m by k matrix, op(B)  a  k by n matrixand C an m by n matrix.
*/

STATIC void avx2_f32_gemm_XtY2x1(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{ 
	for (int col = 0; col < N; ++col) {
		int row = 0;
		for (; row < M - (2 - 1); row += 2) 
            C[row + 1] = avx2_f32_dot2x1(A + row * lda, A + (row + 1) * lda, B, K, C + row);
			//f32_dot3x1(A + ROW * lda, A + (ROW + 1) * lda, A + (ROW + 2) * lda, B, K, C + ROW);			
			//C[row] = f32_dot(A + row * lda, B, K);		
		if (row < M) 
			C[row] = avx2_f32_dot(A + row * lda, B, K);
		B += ldb;
		C += ldc;
	}
 
}
 
STATIC void avx2_f32_gemm_XtY2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{
	int COL;
 	for (COL = 0; COL < N-(2-1); COL+=2) {
		int ROW;
		for (ROW=0; ROW < M-(2-1); ROW+=2) {			
            avx2_f32_dot2x2(A+ROW*lda, A + (ROW+1)*lda, B, B+ldb, K, C+ROW, C +ldc+ROW);
		}
	    if(ROW<M){
			*(C + ldc + ROW)= avx2_f32_dot2x1(B, B + ldb, A + ROW* lda, K, C+ROW);
		}
		B += ldb+ ldb;
		C += ldc+ldc;
	}

	if (COL < N) {
		int ROW;
		for (ROW = 0; ROW < M - (2 - 1); ROW += 2) 
			*(C+ROW+1)= avx2_f32_dot2x1(A + ROW * lda, A + (ROW + 1) * lda, B, K, C + ROW);
		if (ROW < M) 
			C[ROW] = avx2_f32_dot(A + ROW * lda, B,K);
	}
}
////////////////////////////////////////////////////////////////////////////////////////// 


static INLINE void GetOneRowFromMatrix(F32PTR row, F32PTR mat, I32 lda, I32 K,I32 kRemaining, __m256i offset, __m256i mask) 
{
    int i = 0;
    for (; i < K - 7; i += 8)    store(row + i, _mm256_i32gather_ps(mat +i* lda, offset, 4));
    if (kRemaining) {
        __m256 subrow = _mm256_mask_i32gather_ps(set0(), mat+ i*lda, offset, _mm256_castsi256_ps(mask), 4);
        maskstore(row+ i, mask, subrow);
    }
}

static INLINE void GetTwoRowFromMatrix(F32PTR row1, F32PTR row2, F32PTR mat, I32 lda, I32 K, I32 kRemaining, __m256i offset, __m256i mask)
{
    int i = 0;
    for (; i < K - 7; i += 8) {
        store(row1 + i, _mm256_i32gather_ps(mat   + i * lda, offset, 4));
        store(row2 + i, _mm256_i32gather_ps(mat+1 + i * lda, offset, 4));
    }
    if (kRemaining) {
        __m256 subrow1 = _mm256_mask_i32gather_ps(set0(), mat   + i * lda, offset, _mm256_castsi256_ps(mask), 4);
        __m256 subrow2 = _mm256_mask_i32gather_ps(set0(), mat+1 + i * lda, offset, _mm256_castsi256_ps(mask), 4);
        maskstore(row1 + i, mask, subrow1);
        maskstore(row2 + i, mask, subrow2);
    }
}
static INLINE __m256i GetRowOffset( I32 lda ) { 
    __m256i vecLDA = set1_i32(lda);
    __m256i offset = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(0x0706050403020100));
    __m256i offset1 = _mm256_mullo_epi32(offset, vecLDA);
    return offset1;
}
/*
void f32_gemm_XY2x1_old(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{
    // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
    // MSVC seems not to support it

   
    //https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
   // Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
 
    F32     XROW_FIXEXD[4096];
    F32PTR  Xrow = (K <= 4096) ? XROW_FIXEXD : alloca(K * sizeof(F32));
    I32     nRemainder = K % 8;
    __m256i mask = GetMoveMask(nRemainder);
    __m256i offset;
    {
        __m256i vecLDA = set1_i32(lda);
        offset = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(0x0706050403020100));
        offset = _mm256_mul_epi32(offset, vecLDA);
    }

    for (int row = 0; row < M; row++) {
 
        //  Get one row
 
        int i = 0;
        for (; i < K - 7; i += 8)
            store(Xrow + i, _mm256_i32gather_ps(A + row + i * lda, offset, 4));
        if (nRemainder) {
            __m256 subrow = _mm256_mask_i32gather_ps(set0(), A + row + i * lda, offset, _mm256_castsi256_ps(mask), 4);
            maskstore(Xrow + i, mask, subrow);
        }
   
        // End of getting one row
  

        F32PTR Crow = C + row;
        int col = 0;
        for (; col < N - 1; col += 2) {
            *(Crow + ldc * (col + 1)) = avx2_f32_dot2x1(B + col * ldb, B + col * ldb + ldb, Xrow, K, Crow + ldc * col);
        }
        if (col < N) {
            *(Crow + ldc * col) = avx2_f32_dot(B + col * ldb, Xrow, K);
        }

    }

}
*/
STATIC void avx2_f32_gemm_XY1x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{
    // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
    // MSVC seems not to support it

    /*
    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
    Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
    */
    F32     XROW_FIXEXD[4096];
    F32PTR  Xrow         = (K<=4096)? XROW_FIXEXD: alloca(K*sizeof(F32));
    I32     nRemainder   = K % 8;
    __m256i mask         =  GetMoveMask(nRemainder);
    __m256i offset       =  GetRowOffset(lda);
  

    for (int row = 0; row < M; row++) {
    //  Get one row of A
        GetOneRowFromMatrix(Xrow, A + row, lda, K, nRemainder, offset, mask);        
        F32PTR Crow = C + row;
        int    col  = 0;
        for (; col < N - 1; col += 2) 
            *(Crow + ldc * (col+1))=avx2_f32_dot2x1(B+col*ldb,B+col*ldb+ldb,Xrow,K, Crow+ldc*col);        
        if (col < N) 
            *(Crow+ldc*col)=avx2_f32_dot(B + col * ldb, Xrow, K);
    }  

}
STATIC void avx2_f32_gemm_XY2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{ // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
  // MSVC seems not to support it
  //    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
  //   Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
   
    F32     XROW_FIXEXD[4096];
    F32PTR  Xrow = (2*K <= 4096) ? XROW_FIXEXD : alloca(2*K * sizeof(F32));
    F32PTR  Xrow1 = Xrow, Xrow2 = Xrow1 + K;
    I32     nRemainder = K % 8;
    __m256i mask   =  GetMoveMask(nRemainder);
    __m256i offset = GetRowOffset(lda);

    int row = 0;
    for (; row < M-1; row+=2) {
        //  Get one row of A
        GetTwoRowFromMatrix(Xrow1,Xrow2, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;        
        int    col = 0;
        for (; col < N - 1; col += 2)
            avx2_f32_dot2x2(Xrow1, Xrow2, B + col * ldb, B + col * ldb + ldb, K, Crow + ldc * col, Crow+ldc*(col+1));
        if (col < N)
            *(Crow + ldc * col + 1) = avx2_f32_dot2x1(Xrow1,Xrow2,B + col * ldb,  K, Crow + ldc * col);;
   
    }

	if (row < M) {
        //  Get one row of A
        GetOneRowFromMatrix(Xrow1, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;
        int    col = 0;
        for (; col < N - 1; col += 2)
             *(Crow + ldc * (col+1))=avx2_f32_dot2x1(B+col*ldb,B+(col+1)*ldb, Xrow1,K, Crow+col*ldc);              
        if (col < N)
            *(Crow + ldc * col ) = avx2_f32_dot(Xrow1, B+col*ldb,  K);            
	}
}

STATIC void avx2_f32_gemm_XtYt2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{ // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
  // MSVC seems not to support it
  //    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
  //   Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
   
    F32     YROW_FIXEXD[4096];
    F32PTR  Yrow = (2*K <= 4096) ? YROW_FIXEXD : alloca(2*K * sizeof(F32));
    F32PTR  Yrow1 = Yrow, Yrow2 = Yrow1 + K;
    I32     nRemainder = K % 8;
    __m256i mask   =  GetMoveMask(nRemainder);
    __m256i Boffset = GetRowOffset(ldb);
 
    int col = 0;
    for (; col < N-1; col+=2) {
        //  Get one row of B
        GetTwoRowFromMatrix(Yrow1,Yrow2, B + col, ldb, K, nRemainder, Boffset, mask);        
        int    row = 0;
        for (; row < M - 1; row += 2)
            avx2_f32_dot2x2(A+row*lda, A + (row+1)*lda, Yrow1, Yrow2, K, C+row, C+ldc+row);
        if (row < M)
            *(C + ldc + row) = avx2_f32_dot2x1(Yrow1, Yrow2, A + row * lda,  K, C+row);
   
        C += ldc+ldc;
    }
	if (col < N) {
        //  Get one row of B
        GetOneRowFromMatrix(Yrow1, B + col, ldb, K, nRemainder, Boffset, mask);        
        int    row = 0;
        for (; row < M - 1; row += 2)
             *(C+row + 1)=avx2_f32_dot2x1(A+row*lda,A+(row+1)*lda, Yrow1,K, C+row);              
        if (row < M)
            *(C + row) = avx2_f32_dot(A+row*lda, Yrow1, K);          
	}
}


/////////////////////////////////////////////////////////////////////////////
static INLINE F32 __avx2_f32_dot2x1stride(F32PTR  x, F32PTR y, F32PTR v, int N, __m256i mask, __m256i voffset, I32 stride,  F32PTR res) {
    __m256 RX = set0();
    __m256 RY = set0();
    int i = 0;
    for (; i < N - (NF2 - 1); i += NF2) {
        __m256 v0 = _mm256_i32gather_ps(v + i * stride,     voffset, 4);
        __m256 v1 = _mm256_i32gather_ps(v + (i+8) * stride, voffset, 4);
        __m256 x0 = mul(load(x + i),     v0);
        __m256 x1 = mul(load(x + i + NF), v1);
        __m256 y0 = mul(load(y + i),     v0);
        __m256 y1 = mul(load(y + i + NF), v1);

        RX = add(RX, add(x0, x1));
        RY = add(RY, add(y0, y1));
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 V = _mm256_i32gather_ps(v + i * stride, voffset, 4);
        RX = add(RX, mul(load(x + i), V));
        RY = add(RY, mul(load(y + i), V));
    }

    int n = N - i;
    if (n > 0) {
        __m256  vvec = _mm256_mask_i32gather_ps(set0(), v + i*stride, voffset, _mm256_castsi256_ps(mask), 4);
        RX = add(RX, mul(maskload(x + i, mask), vvec));
        RY = add(RY, mul(maskload(y + i, mask), vvec));
    }

    float sumX = f32_hsum(RX);
    float sumY = f32_hsum(RY);

    _mm256_zeroupper();
    res[0] = sumX;
    return sumY;
}

static INLINE F32 __avx2_f32_dot_stride(F32PTR  x, F32PTR y, int N, __m256i mask, __m256i yoffset, I32 stride) {

    __m256 r = set0();    
    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {
        __m256 Y = _mm256_i32gather_ps(y + i * stride, yoffset, 4);
        r = add(  r, mul( load(x + i), Y ) );
    }    
    int n = N - i;
    if (n > 0) {        
        __m256 Y = _mm256_mask_i32gather_ps(set0(), y + i * stride, yoffset, _mm256_castsi256_ps(mask), 4);
        r = add(r,  mul(   maskload(x + i, mask), Y  )   );
    }    
    F32 sum = f32_hsum(r);
     _mm256_zeroupper();    
    return sum;
}

STATIC void avx2_f32_gemm_XYt2x1(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{
    // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
    // MSVC seems not to support it
    /*
    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
    Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
    */
    F32     XROW_FIXEXD[4096];
    F32PTR  Xrow  = (2 * K <= 4096) ? XROW_FIXEXD : alloca( 2*K * sizeof(F32) );
    F32PTR  Xrow1 = Xrow, Xrow2 = Xrow1 + K;
    I32     nRemainder   = K % 8;
    __m256i mask         = GetMoveMask(nRemainder);
    __m256i offset       = GetRowOffset(lda);
    __m256i BOffset      = GetRowOffset(ldb);
 
    int row = 0;
    for (; row < M - 1; row += 2) {
        //  Get one row of A
        GetTwoRowFromMatrix(Xrow1, Xrow2, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;
        int    col = 0;
        for (; col < N; col++) {
            *(Crow + ldc * col + 1) = __avx2_f32_dot2x1stride(
                Xrow1, Xrow2, B + col, K, mask,BOffset, ldb, Crow + ldc * col);;
        }

    }
    if (row < M) {
        //  Get one row of A
        GetOneRowFromMatrix(Xrow1, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;
        int    col = 0;
        for (; col < N; col++)
            *(Crow + ldc * col) = __avx2_f32_dot_stride(Xrow1, B+col, K,mask,BOffset,ldb);
    }
}


///////////////////////////////////////////////////////////////////////////

STATIC void avx2_f32_gemv_Xy1x1_slow(int N, int K, F32PTR X, int lda, F32PTR y, F32PTR C)
{     
    I32     nRemainder = K % 8;
    __m256i mask       = GetMoveMask(nRemainder);
    __m256i Xoffset    = GetRowOffset(lda);    
    for (int row = 0; row < N ; ++row){
       C[row] =__avx2_f32_dot_stride(y, X + row, K, mask, Xoffset, lda);
    }

}


static INLINE void fma_f32_axpy_inplace( const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m256  C = set1(a);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(y + i,       _mm256_fmadd_ps(load(x + i), C, load(y + i))),
        store(y + i + NF,  _mm256_fmadd_ps(load(x + i + NF),  C, load(y + i + NF))),
        store(y + i + NF2, _mm256_fmadd_ps(load(x + i + NF2), C, load(y + i + NF2))),
        store(y + i + NF3, _mm256_fmadd_ps(load(x + i + NF3), C, load(y + i + NF3)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, _mm256_fmadd_ps(load(x + i), C, load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), _mm256_fmadd_ps(load(x + i), C, load(y + i)) );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

static  void avx2_f32_axpy_inplace(const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m256  C = set1(a);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3)
        store(y + i, add(mul(load(x + i), C), load(y + i))),
        store(y + i + 8, add(mul(load(x + i + 8), C), load(y + i + 8))),
        store(y + i + 16, add(mul(load(x + i + 16), C), load(y + i + 16)));
    //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, add(mul(load(x + i), C), load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), add(mul(load(x + i), C), load(y + i)));
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}


static INLINE void __avx2_f32_axpy_inplace(const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m256  C = set1(a);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3)
        store(y + i, add(mul(load(x + i), C),       load(y + i))),
        store(y + i+8, add(mul(load(x + i+8), C),   load(y + i+8))),
        store(y + i+16, add(mul(load(x + i+16), C), load(y + i+16)));
        //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, add(mul(load(x + i), C), load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), add(mul(load(x + i), C), load(y + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

static INLINE void __avx2_f32_ax1ax2py_inplace(const F32PTR a, const F32PTR x1, const F32PTR x2, F32PTR y, const int N) {
    //y=a*x+y;
    __m256  C1 = set1(a[0]);
    __m256  C2 = set1(a[1]);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3) {
        __m256 t1 = add(mul(load(x1 + i), C1), mul(load(x2 + i), C2)); store(y + i, add(t1, load(y + i)));
        __m256 t2 = add(mul(load(x1 + i + 8), C1), mul(load(x2 + i + 8), C2)); store(y + i + 8, add(t2, load(y + i + 8)));
        __m256 t3 = add(mul(load(x1 + i+16), C1), mul(load(x2 + i+16), C2)); store(y + i+16, add(t3, load(y + i+16)));    
    }
        //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF) {
        __m256 t1 = add(mul(load(x1 + i), C1), mul(load(x2 + i), C2)); store(y + i, add(t1, load(y + i)));
    }    
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __m256 t1 = add(mul(load(x1 + i), C1), mul(load(x2 + i), C2)); 
        maskstore(y + i, GetMoveMask(n), add(t1, load(y + i)));      
    }  
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_f32_gemv_Xb(int N, int K, F32PTR X, int lda, F32PTR b, F32PTR C)
{
    //I32     nRemainder = K % 8;
    //__m256i mask    = GetMoveMask(nRemainder);
    //__m256i Xoffset = GetRowOffset(lda);
    memset(C, 0, sizeof(F32) * N);
    int row = 0;
    for (; row < N - (256 * 2 -1); row += 256*2) {
        int col = 0;
        for (; col < K-1; col+=2) {
            __avx2_f32_ax1ax2py_inplace(b+col, X + col * lda + row, X + (col+1) * lda + row, C + row, 256 * 2);
            //avx2_f32_axpy_inplace_v2(b[col], X + col * lda+row, C+row, 256);
        }
        if(col<K)
            __avx2_f32_axpy_inplace(b[col], X + col * lda + row, C + row, 256 * 2);
    }
    int n = N - row;
    if (n > 0) {
        int col = 0;
        for (; col < K-1; col+=2) {
            __avx2_f32_ax1ax2py_inplace(b+col, X + col * lda + row, X + (col+1) * lda + row, C + row, n);
        }
        if(col<K)
            __avx2_f32_axpy_inplace(b[col], X + col * lda + row, C + row, n);
    }
    

}


////////////////////////////////////////////////////////////////
static INLINE I32 _avx2_f32_findindex_LT(F32PTR  x, I32PTR indices, F32 value, int N) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m256  VALUE = set1(value);

    int cnt = 0;
    int i   = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 cmpmask = _mm256_cmp_ps(load(x + i), VALUE, _CMP_LT_OQ);
        int    mask    = _mm256_movemask_ps(cmpmask);
        int    segIdx  = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed uo: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask >>= 1L;
        }
    }
    int n = N - i;
    if (n > 0) {
        __m256i movemask = GetMoveMask(n);
        __m256  cmpmask = _mm256_cmp_ps(maskload(x + i, movemask), VALUE, _CMP_LT_OQ);
        int     mask     = _mm256_movemask_ps(_mm256_and_ps(cmpmask, _mm256_castsi256_ps(movemask)));
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask = mask >> 1L;
        }
    }

    _mm256_zeroupper();

    return cnt;


    

}
static INLINE I32 _avx2_f32_findindex_LE(F32PTR  x, I32PTR indices, F32 value, int N) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m256  VALUE = set1(value);

    int  cnt = 0;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 cmpmask = _mm256_cmp_ps(load(x + i), VALUE, _CMP_LE_OQ);
        int    mask = _mm256_movemask_ps(cmpmask);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed uo: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask >>= 1L;
        }
    }
    int n = N - i;
    if (n > 0) {
        __m256i movemask = GetMoveMask(n);
        __m256  cmpmask = _mm256_cmp_ps(maskload(x + i, movemask), VALUE, _CMP_LE_OQ);
        int    mask = _mm256_movemask_ps(_mm256_and_ps(cmpmask, _mm256_castsi256_ps(movemask)));
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask = mask >> 1L;
        }
    }

    _mm256_zeroupper();

    return cnt;


    

}
static INLINE I32 _avx2_f32_findindex_GT(F32PTR  x, I32PTR indices, F32 value, int N) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m256  VALUE = set1(value);

    int  cnt = 0;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 cmpmask = _mm256_cmp_ps(load(x + i), VALUE, _CMP_GT_OQ);
        int    mask = _mm256_movemask_ps(cmpmask);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed uo: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask >>= 1L;
        }
    }
    int n = N - i;
    if (n > 0) {
        __m256i movemask = GetMoveMask(n);
        __m256  cmpmask = _mm256_cmp_ps(maskload(x + i, movemask), VALUE, _CMP_GT_OQ);
        int    mask = _mm256_movemask_ps(_mm256_and_ps(cmpmask, _mm256_castsi256_ps(movemask)));
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask = mask >> 1L;
        }
    }

    _mm256_zeroupper();

    return cnt;


    

}
static INLINE I32 _avx2_f32_findindex_GE(F32PTR  x, I32PTR indices, F32 value, int N) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m256  VALUE = set1(value);

    int  cnt = 0;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 cmpmask = _mm256_cmp_ps(load(x + i), VALUE, _CMP_GE_OQ);
        int    mask = _mm256_movemask_ps(cmpmask);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed uo: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask >>= 1L;
        }
    }
    int n = N - i;
    if (n > 0) {
        __m256i movemask = GetMoveMask(n);
        __m256  cmpmask = _mm256_cmp_ps(maskload(x + i, movemask), VALUE, _CMP_GE_OQ);
        int    mask = _mm256_movemask_ps(_mm256_and_ps(cmpmask, _mm256_castsi256_ps(movemask)));
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask = mask >> 1L;
        }
    }

    _mm256_zeroupper();

    return cnt;


    

}
static INLINE I32 _avx2_f32_findindex_EQ(F32PTR  x, I32PTR indices, F32 value, int N) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m256  VALUE = set1(value);
    __m256  EPS   = set1(1e-20);
    __m256  ABSMASK = _mm256_castsi256_ps(set1_i32(0x7FFFFFFF));

    int  cnt = 0;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m256 cmpmask = _mm256_cmp_ps(_mm256_and_ps(sub(load(x + i), VALUE), ABSMASK), EPS, _CMP_LT_OQ);
        int    mask    = _mm256_movemask_ps(cmpmask);
        int    segIdx  = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed uo: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask >>= 1L;
        }
    }
    int n = N - i;
    if (n > 0) {
        __m256i movemask = GetMoveMask(n);
        __m256 cmpmask = _mm256_cmp_ps(_mm256_and_ps(sub(maskload(x + i, movemask), VALUE), ABSMASK), EPS, _CMP_LT_OQ);
        int    mask = _mm256_movemask_ps(_mm256_and_ps(cmpmask, _mm256_castsi256_ps(movemask)));
        int    segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt += keep;
            mask = mask >> 1L;
        }
    }

    _mm256_zeroupper();
    return cnt;


    

}
STATIC I32 avx2_f32_findindex(F32PTR  x, I32PTR indices, F32 value, int N, CmpFlag flag) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */

    switch (flag) {
    case CMP_LT:
        return _avx2_f32_findindex_LT(x, indices, value, N);
    case CMP_LE:
        return _avx2_f32_findindex_LE(x, indices, value, N);
    case CMP_GT:
        return _avx2_f32_findindex_GT(x, indices, value, N);
    case CMP_GE:
        return _avx2_f32_findindex_GE(x, indices, value, N);
    case CMP_EQ:
        return _avx2_f32_findindex_EQ(x, indices, value, N);

    }
    return 0;
}


STATIC void  avx2_f32_scatter_val_byindex(F32PTR  x, I32PTR indices, F32 value, int N) {
// No efficent avx2 version avaialble, so the generic version is used here
   
    #define UNROLL_NUMBER  4
	const int regularPart = N & (-UNROLL_NUMBER); // 4: 100 ; -4: ffff00;    
    int  i   = 0;
    for (; i < regularPart; i += UNROLL_NUMBER) {
        x[indices[i]] = value;
        x[indices[i+1]] = value;
        x[indices[i+2]] = value;
        x[indices[i+3]] = value;
    }   
    for (; i < N; ++i)    x[indices[i]] = value; 
}
static void avx2_f32_scatter_val_byindex_slow_not_used(F32PTR  x, I32PTR indices, F32 value, int N) {
    int i = 0;
    __m128  vec = _mm_set1_ps(value);
    for (; i < N - (NF - 1); i += NF) {
        __m256i index  = loadi(indices+i);
        __m128i indlow = _mm256_castsi256_si128(index);
        int i0 = _mm_extract_epi32(indlow, 0);       _mm_store_ss(x+i0,   vec);
        int i1 = _mm_extract_epi32(indlow, 1);       _mm_store_ss(x + i1, vec);
        int i2 = _mm_extract_epi32(indlow, 2);       _mm_store_ss(x + i2, vec);
        int i3 = _mm_extract_epi32(indlow, 3);       _mm_store_ss(x + i3, vec);
        __m128i indhigh = _mm256_extracti128_si256(index, 1);
        int i4 = _mm_extract_epi32(indhigh, 0);      _mm_store_ss(x + i4, vec);
        int i5 = _mm_extract_epi32(indhigh, 1);      _mm_store_ss(x + i5, vec);
        int i6 = _mm_extract_epi32(indhigh, 2);     _mm_store_ss(x + i6, vec);
        int i7 = _mm_extract_epi32(indhigh, 3);      _mm_store_ss(x + i7, vec);
    }
    int n = N - i;
    if (n >= 4) {
        __m128i indlow = _mm_loadu_si128(indices+i);        
        int i0 = _mm_extract_epi32(indlow, 0);       _mm_store_ss(x + i0, vec);
        int i1 = _mm_extract_epi32(indlow, 1);      _mm_store_ss(x + i1, vec);
        int i2 = _mm_extract_epi32(indlow, 2);       _mm_store_ss(x + i2, vec);
        int i3 = _mm_extract_epi32(indlow, 3);      _mm_store_ss(x + i3, vec);
        i += 4;
    }
    _mm256_zeroupper();
    for (; i < N; ++i)    x[indices[i]] = value;
    
}



STATIC void avx2_f32_scatter_vec_byindex(F32PTR  x, I32PTR indices, F32PTR values, int N) {

     int i = 0;
     #define f2i _mm_castps_si128
     #define i2f _mm_castsi128_ps
     
    for (; i < N - (NF - 1); i += NF) {
        __m256  value    = load(values + i); 
        __m128  valueLow = _mm256_castps256_ps128(value);
         _mm_store_ss(x + indices[i],   valueLow );
         _mm_store_ss(x + indices[i+1], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),4) ));
         _mm_store_ss(x + indices[i+2], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),8) ));
         _mm_store_ss(x + indices[i+3], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),12) ));
         __m128  valueHigh = _mm256_extractf128_ps(value, 1);
         _mm_store_ss(x + indices[i+4], valueHigh);
         _mm_store_ss(x + indices[i+5], i2f(  _mm_alignr_epi8(f2i(valueHigh), f2i(valueHigh),4) ));
         _mm_store_ss(x + indices[i+6], i2f(  _mm_alignr_epi8(f2i(valueHigh), f2i(valueHigh),8) ));
         _mm_store_ss(x + indices[i+7], i2f(  _mm_alignr_epi8(f2i(valueHigh), f2i(valueHigh),12) ));
    }
    int n = N - i;
    if (n >= 4) {
        __m128 valueLow = _mm_loadu_ps(values +i);
         _mm_store_ss(x + indices[i],   valueLow );
         _mm_store_ss(x + indices[i+1], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),4) ));
         _mm_store_ss(x + indices[i+2], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),8) ));
         _mm_store_ss(x + indices[i+3], i2f(  _mm_alignr_epi8(f2i(valueLow), f2i(valueLow),12) ));
    }
    _mm256_zeroupper();
    for (; i < N; ++i)    x[indices[i]] = values[i];
    
} 
STATIC void avx2_f32_gatherVec_scatterVal_byindex(F32PTR  x, I32PTR indices, F32PTR values, F32 newValue, int N) {

     int i = 0;   
    for (; i < N - (NF - 1); i += NF) {
        __m256i index  = loadi(indices + i);   //   __m128i indlow = _mm256_castsi256_si128(index);
        v8sf  oldValues = _mm256_i32gather_ps(x, index, 4);
        store(values + i, oldValues);
        x[indices[i]]   = newValue;
        x[indices[i+1]] = newValue;
        x[indices[i+2]] = newValue;
        x[indices[i+3]] = newValue;
        x[indices[i+4]] = newValue;
        x[indices[i+5]] = newValue;
        x[indices[i+6]] = newValue;
        x[indices[i+7]] = newValue;
    }
    _mm256_zeroupper();
    int n = N - i;
    if (n >= 4) {
        __m128i index     = _mm_loadu_si128(indices + i);   //   __m128i indlow = _mm256_castsi256_si128(index);
        v4sf    oldValues = _mm_i32gather_ps(x, index, 4);
        _mm_storeu_ps(values + i, oldValues);
        x[indices[i]]   = newValue;
        x[indices[i+1]] = newValue;
        x[indices[i+2]] = newValue;
        x[indices[i+3]] = newValue;        
        i = i + 4;
    }    
    for (; i < N; ++i) {
        values[i] = x[indices[i]];
        x[indices[i]] = newValue;
    }
    
} 
STATIC void avx2_f32_gather2Vec_scatterVal_byindex(F32PTR  x, F32PTR  y, I32PTR indices, F32PTR values, F32 newValue, int N) {

     int i = 0;   
    for (; i < N - (NF - 1); i += NF) {
        __m256i index   = loadi(indices + i);   //   __m128i indlow = _mm256_castsi256_si128(index);
        v8sf  oldXValues = _mm256_i32gather_ps(x, index, 4);
        store(values + i, oldXValues);
        v8sf  oldYValues = _mm256_i32gather_ps(y, index, 4);
        store(values +N + i, oldYValues);
        y[indices[i]]=  x[indices[i]]  = newValue;
        y[indices[i+1]] = x[indices[i+1]] = newValue;
        y[indices[i+2]] = x[indices[i+2]] = newValue;
        y[indices[i+3]] = x[indices[i+3]] = newValue;
        y[indices[i+4]] = x[indices[i+4]] = newValue;
        y[indices[i+5]] = x[indices[i+5]] = newValue;
        y[indices[i+6]] = x[indices[i+6]] = newValue;
        y[indices[i+7]] = x[indices[i+7]] = newValue;
    }
    _mm256_zeroupper();
    int n = N - i;
    if (n >= 4) {
        __m128i index     = _mm_loadu_si128(indices + i);   //   __m128i indlow = _mm256_castsi256_si128(index);
        v4sf    oldValues = _mm_i32gather_ps(x, index, 4);
        _mm_storeu_ps(values + i, oldValues);
        v4sf    oldYValues = _mm_i32gather_ps(y, index, 4);
        _mm_storeu_ps(values + N + i, oldYValues);
        y[indices[i]] =x[indices[i]]   = newValue;
        y[indices[i]] =x[indices[i+1]] = newValue;
        y[indices[i]] =x[indices[i+2]] = newValue;
        y[indices[i]] =x[indices[i+3]] = newValue;
        i = i + 4;
    }    
    for (; i < N; ++i) {
        values[i] = x[indices[i]];
        values[N+i] = y[indices[i]];
        x[indices[i]] = newValue;
        y[indices[i]] = newValue;
    }
    
} 
STATIC void avx2_f32_scale_inplace(const F32 gain, const F32 offset,const F32PTR x, const int N) {
    __m256  G = set1(gain);
    __m256  O = set1(offset);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(x + i,       add(mul(load(x + i),       G), O) ),
        store(x + i + NF,  add(mul(load(x + i + NF),  G), O)),
        store(x + i + NF2, add(mul(load(x + i + NF2), G), O)),
        store(x + i + NF3, add(mul(load(x + i + NF3), G), O));
    for (; i < N - (NF - 1); i += NF)
        store(x + i, add(mul(load(x + i), G), O) );
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i, GetMoveMask(n), add(mul(load(x + i), G), O) );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx2_f32_hinge_pos(const F32PTR X, const F32PTR Y, const F32 knot, const int N){
    __m256  O = set0();
    __m256  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 d1 = sub(load(X + i),    C);   store(Y + i,       _mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d2 = sub(load(X + i+NF), C);   store(Y + i + NF,  _mm256_blendv_ps(O, d2, _mm256_cmp_ps(d2, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d3 = sub(load(X + i+NF2), C);  store(Y + i + NF2, _mm256_blendv_ps(O, d3, _mm256_cmp_ps(d3, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d4 = sub(load(X + i+NF3), C);  store(Y + i + NF3, _mm256_blendv_ps(O, d4, _mm256_cmp_ps(d4, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 d1 = sub(load(X + i), C);   store(Y + i, _mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m256 d1 = sub(load(X + i), C);   maskstore(Y + i, GetMoveMask(n),_mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}
STATIC void avx2_f32_hinge_neg(const F32PTR X, const F32PTR Y, const F32 knot, const int N){
    __m256  O = set0();
    __m256  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 d1 = sub(C, load(X + i));      store(Y + i,       _mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d2 = sub(C, load(X + i+NF));   store(Y + i + NF,  _mm256_blendv_ps(O, d2, _mm256_cmp_ps(d2, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d3 = sub(C, load(X + i+NF2));  store(Y + i + NF2, _mm256_blendv_ps(O, d3, _mm256_cmp_ps(d3, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d4 = sub(C, load(X + i+NF3));  store(Y + i + NF3, _mm256_blendv_ps(O, d4, _mm256_cmp_ps(d4, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 d1 = sub(C,load(X + i));   store(Y + i, _mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m256 d1 = sub(C,load(X + i));   maskstore(Y + i, GetMoveMask(n),_mm256_blendv_ps(O, d1, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}

STATIC void avx2_f32_step_pos(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N){
    __m256  O = set0();
    __m256  I = set1(1.0);
    __m256  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 d1 = sub(load(X + i),    C);   store(Z + i,       _mm256_blendv_ps(O, load(Y + i),       _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d2 = sub(load(X + i+NF), C);   store(Z + i + NF,  _mm256_blendv_ps(O, load(Y + i+NF),    _mm256_cmp_ps(d2, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d3 = sub(load(X + i+NF2), C);  store(Z + i + NF2, _mm256_blendv_ps(O, load(Y + i+NF2),   _mm256_cmp_ps(d3, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d4 = sub(load(X + i+NF3), C);  store(Z + i + NF3, _mm256_blendv_ps(O, load(Y + i+NF3),   _mm256_cmp_ps(d4, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 d1 = sub(load(X + i), C);     store(Z + i,       _mm256_blendv_ps(O, load(Y + i),        _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m256 d1 = sub(load(X + i), C);   maskstore(Z + i, GetMoveMask(n),_mm256_blendv_ps(O, load(Y + i), _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}
STATIC  void avx2_f32_step_neg(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N){
    __m256  O = set0();
    __m256  I = set1(1.0);
    __m256  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m256 d1 = sub(load(X + i),    C);   store(Z + i,       _mm256_blendv_ps(load(Y + i ),      O, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d2 = sub(load(X + i+NF), C);   store(Z + i + NF,  _mm256_blendv_ps(load(Y + i + NF),  O, _mm256_cmp_ps(d2, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d3 = sub(load(X + i+NF2), C);  store(Z + i + NF2, _mm256_blendv_ps(load(Y + i + NF2), O, _mm256_cmp_ps(d3, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
        __m256 d4 = sub(load(X + i+NF3), C);  store(Z + i + NF3, _mm256_blendv_ps(load(Y + i + NF3), O, _mm256_cmp_ps(d4, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m256 d1 = sub(load(X + i), C);     store(Z + i,       _mm256_blendv_ps(load(Y + i),        O,  _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m256 d1 = sub(load(X + i), C);   maskstore(Z + i, GetMoveMask(n),_mm256_blendv_ps(load(Y + i), O, _mm256_cmp_ps(d1, O, _CMP_GE_OQ))); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();     
}
 void SetupVectorFunction_AVX2() {

    FillMaskTempalte();

    i32_add_val_inplace = &avx2_i32_add_val_inplace;;
    i32_sum   = &avx2_i32_sum;
    f32_fill_val = &avx2_f32_fill_val;
    f32_sum  = &avx2_f32_sum;
    f32_add_vec   = &avx2_f32_add_vec;
    f32_sub_vec    = &avx2_f32_sub_vec;
    f32_add_vec_inplace    = &avx2_f32_add_vec_inplace;
    f32_sub_vec_inplace    = &avx2_f32_sub_vec_inplace;
    f32_subrev_val_inplace    = &avx2_f32_subrev_val_inplace;
    f32_add_val_inplace    = &avx2_f32_add_val_inplace;
    f32_mul_val_inplace    = &avx2_f32_mul_val_inplace;
    f32_mul_vec_inplace    = &avx2_f32_mul_vec_inplace;
    f32_mul_vec    = &avx2_f32_mul_vec;
    f32_dot    = &avx2_f32_dot;

    
    f32_dot2x1    = &avx2_f32_dot2x1;
    f32_dot2x2    = &avx2_f32_dot2x2;
    f32_add_v_v2_vec_inplace    = &avx2_f32_add_v_v2_vec_inplace;
  
    f32_cos_vec_inplace    = &avx2_f32_cos_vec_inplace;  
    f32_sin_vec_inplace    = &avx2_f32_sin_vec_inplace;   
    f32_sincos_vec_inplace   = & avx2_f32_sincos_vec_inplace;
 
    f32_pow_vec_inplace   = & avx2_f32_pow_vec_inplace;
    f32_log_vec_inplace    = &avx2_f32_log_vec_inplace;
    f32_exp_vec_inplace     = &avx2_f32_exp_vec_inplace;
    f32_sqrt_vec_inplace    = &avx2_f32_sqrt_vec_inplace;
    f32_sqrt_vec    = &avx2_f32_sqrt_vec;
    //f32_sumlog_slow    = &avx2_f32_sumlog_slow;
    //f32_sumlog    = &avx2_f32_sumlog;
  
    f32_avgstd                      = &avx2_f32_avgstd;
    f32_sx_sxx_to_avgstd_inplace   = & avx2_f32_sx_sxx_to_avgstd_inplace;
  
    f32_maxidx    = &avx2_f32_maxidx;
    f32_minidx    = &avx2_f32_minidx;
    
    f32_diff_back   = & avx2_f32_diff_back;
  
    f32_seq   = & avx2_f32_seq;
    i32_seq = &avx2_i32_seq;

    /////////////////2
    f32_to_f64_inplace   = &avx2_f32_to_f64_inplace;
    f64_to_f32_inplace   = &avx2_f64_to_f32_inplace;
    
    i32_to_f32_scaleby_inplace    = &avx2_i32_to_f32_scaleby_inplace;
    i32_increment_bycond_inplace  = &  avx2_i32_increment_bycond_inplace;
    i32_increment_vec2_bycond_inplace = &avx2_i32_increment_vec2_bycond_inplace;
    //i08_find_nth_onebyte_binvec     = & avx2_i08_find_nth_onebyte_binvec;
    //i08_find_nth_onebyte_binvec_v2  = &  avx2_i08_find_nth_onebyte_binvec_v2;
    i08_sum_binvec   = &avx2_i08_sum_binvec;
   
    f32_gemm_XtY2x1  = avx2_f32_gemm_XtY2x1;
    f32_gemm_XtY2x2  = avx2_f32_gemm_XtY2x2;
    f32_gemm_XY1x2   = avx2_f32_gemm_XY1x2;
    f32_gemm_XY2x2   = avx2_f32_gemm_XY2x2;
    f32_gemm_XtYt2x2 = avx2_f32_gemm_XtYt2x2;
    f32_gemm_XYt2x1  = avx2_f32_gemm_XYt2x1;
    
    f32_gemv_Xb = &avx2_f32_gemv_Xb;

    f32_findindex = avx2_f32_findindex;
    f32_scatter_vec_byindex = avx2_f32_scatter_vec_byindex;
    f32_gatherVec_scatterVal_byindex = avx2_f32_gatherVec_scatterVal_byindex;
    f32_gather2Vec_scatterVal_byindex = &avx2_f32_gather2Vec_scatterVal_byindex;

    f32_scale_inplace = avx2_f32_scale_inplace;
    f32_hinge_pos = avx2_f32_hinge_pos;
    f32_hinge_neg = avx2_f32_hinge_neg;
    f32_step_pos = avx2_f32_step_pos;
	f32_step_neg  = avx2_f32_step_neg;
    f32_axpy_inplace = avx2_f32_axpy_inplace;
}

#else
static char a = 'a';
#endif



///////////////////////////////////////////////////////////////////////////
#if defined(COMPILER_CLANG) && !defined(cpu_ARM64)
    //pragma clang attribute push (__attribute__((target("avx,avx2"))), apply_to=function)
    #pragma clang attribute pop
#endif
///////////////////////////////////////////////////////////////////////////


#include "abc_000_warning.h"