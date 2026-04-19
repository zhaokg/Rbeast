#include <stdio.h>
#include <string.h>
#include <math.h>         //sqrtf
#include "abc_000_macro.h"
#include "abc_000_warning.h"

#define STATIC   static

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
    #pragma clang optimize on
    //https://stackoverflow.com/questions/46165752/does-clang-have-something-like-pragma-gcc-target
   // #pragma clang attribute push (__attribute__((target("sse,sse2,sse3,ssse3,sse4,popcnt,avx,avx2"))), apply_to=function)
	#pragma clang attribute push (__attribute__((target("avx,avx2,avx512f,avx512dq,avx512bw,avx512vl"))), apply_to=function)
    //#pragma clang attribute pop
#endif

#if  defined(COMPILER_GCC) && !defined(cpu_ARM64) 
    //https://www.geeksforgeeks.org/speed-up-naive-algorithms-in-competitive-coding-in-c-cpp/
    //https://codeforces.com/blog/entry/78897
    //https://stackoverflow.com/questions/61759552/why-some-top-level-competitive-programmers-use-pragma    
    //#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
    //#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
    #pragma optimization_level 3

    //https://stackoverflow.com/questions/43152633/invalid-register-for-seh-savexmm-in-cygwin
    #pragma GCC optimize("O3,Ofast,inline,omit-frame-pointer,no-asynchronous-unwind-tables")  //Comment optimisations for interactive problems (use endl)
    //#pragma GCC target("avx,avx2,fma")
     //#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,avx,avx2,fma,tune=haswell")
	 #pragma GCC target("avx,avx2,avx512f,avx512dq,avx512bw,avx512vl") //https://en.wikipedia.org/wiki/AVX-512
    // -mavx256-split-unaligned-load is spit out by default, hurting new CPUS
    // an altertive is to set -march=hashwell, but it doesntot work in pragram for all gcc versions
    // #pragma GCC optimization ("unroll-loops")
    //#pragma GCC target "arch=core-avx2,tune=core-avx2"
     
    // https://stackoverflow.com/questions/51003218/gcc-target-for-avx2-disabling-sse-instruction-set
    // #pragma GCC target "arch=haswell" : // A buf for some gcc versions, not working
    
#endif
///////////////////////////////////////////////////////////////////////////

//Move the pragmas before abc_math_avx512 to enable the instruction sets for the REMAINING CODE

#include "abc_vec.h"

#if defined(OS_WIN64) || defined(OS_WIN32)
    #include <malloc.h> //alloca
#else
    #include <alloca.h> //alloca
#endif

#ifdef COMPILER_MSVC
    #define __attribute__(x)
#endif



#define NF        16
#define NF2       (NF*2)
#define NF3       (NF*3)
#define NF4       (NF*4)
#define load        _mm512_loadu_ps
#define store       _mm512_storeu_ps
#define loadi       _mm512_loadu_si512
#define storei      _mm512_storeu_si512
#define maskload(src,mask)    _mm512_mask_loadu_ps(  _mm512_setzero_ps(),mask,src)
#define maskloadi(src,mask)   _mm512_mask_loadu_epi32(_mm512_setzero_epi32(),mask,src)
#define maskstore             _mm512_mask_storeu_ps
#define maskstorei            _mm512_mask_storeu_epi32
#define set1        _mm512_set1_ps
#define set1_i32    _mm512_set1_epi32
#define set0        _mm512_setzero_ps
#define set0i       _mm512_setzero_si512
#define mul         _mm512_mul_ps
#define add         _mm512_add_ps
#define sub         _mm512_sub_ps
#define addi32      _mm512_add_epi32

#if !defined(COMPILER_SOLARIS) && defined(TARGET_64) && !defined(cpu_ARM64)
#include <immintrin.h>
#include "abc_math_avx512.h"

// GCC will impiclty convert all the x+i read/write using aligned load/store while
// MSVC and ICC use the unligned load/store.  So for GCC, we have to expliclty use 
// the unaligned load/store to avoid seg faults. That is,
// addi32(*((__m512i*) (x + i)), C)--> addi32 (loadu(x+i),C)

static __mmask16 masktemplate[16];
static INLINE void      FillMaskTemplate(void) {   for (I32 i = 0; i < 16; i++)    masktemplate[i] = (1UL << i) - 1UL;  }
static INLINE __mmask16  __attribute__((always_inline)) GetMoveMask(int n)  {
    return masktemplate[n];
}


static INLINE F32 __attribute__((always_inline)) f32_hsum_m256( __m256 r)  {
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
static INLINE I32 __attribute__((always_inline)) i32_hsum_m256(__m256i r)  {
    __m128i vlow = _mm256_castsi256_si128(r);
    __m128i vhigh = _mm256_extracti128_si256(r, 1); // high 128
    __m128i v128 = _mm_add_epi32(vlow, vhigh);     // add the low 128

    __m128i hi64 = _mm_unpackhi_epi64(v128, v128);           // 3-operand non-destructive AVX lets us save a byte without needing a mov
   //__m128i hi64 = _mm_shuffle_epi32(x, _MM_SHUFFLE(1, 0, 3, 2));

    __m128i sum64 = _mm_add_epi32(hi64, v128);
    __m128i hi32 = _mm_shufflelo_epi16(sum64, _MM_SHUFFLE(1, 0, 3, 2));    // Swap the low two elements
    __m128i sum32 = _mm_add_epi32(sum64, hi32);
    return _mm_cvtsi128_si32(sum32);       // SSE2 movd    
}

static INLINE F32 __attribute__((always_inline)) f32_hsum( __m512 r)   {

    __m256 low  =  _mm512_castps512_ps256(r);
    __m256 high =  _mm512_extractf32x8_ps(r, 1);
    __m256 R    =  _mm256_add_ps(low, high);
    return f32_hsum_m256(R);
}
static INLINE I32 __attribute__((always_inline)) i32_hsum(__m512i r) {
    __m256i low = _mm512_castsi512_si256(r);
    __m256i high = _mm512_extracti32x8_epi32(r, 1);
    __m256i R = _mm256_add_epi32(low, high);
    return i32_hsum_m256(R);// SSE2 movd    
}


STATIC void avx512_i32_add_val_inplace(const int c, const I32PTR x, const int N) {
    __m512i  C = set1_i32(c);
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
STATIC I32  avx512_i32_sum(const I32PTR x, const int N) {
    __m512i  S = set0i();
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512i s12   = addi32(loadi(x + i),        loadi(x + i + NF));
        __m512i s34   = addi32(loadi(x + i + NF2),  loadi(x + i + NF3));
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
STATIC F32  avx512_f32_sum(const I32PTR x, const int N) {
    __m512  S = set0();
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 s12      = add(load(x + i),        load(x + i + NF) );
        __m512 s34      = add(load(x + i + NF2),  load(x + i + NF3));
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
STATIC void avx512_f32_fill_val(const F32 c, F32PTR x, int N) {
    __m512  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store((x + i),     C);
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
STATIC void avx512_f32_add_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       //dst=dst-src
        store(dst + i,        add(load(src2 + i),       load(src1 + i)));
        store(dst + i + NF,   add(load(src2 + i + NF),  load(src1 + i  + NF)));
        store(dst + i + NF2,  add(load(src2 + i + NF2), load(src1 + i + NF2)));
        store(dst + i + NF3,  add(load(src2 + i + NF3), load(src1 + i + NF3)));
    }
    for (; i < N - (NF - 1); i += NF) 
        store(dst + i,   add(load(src2 + i), load(src1 + i))  );    
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(dst + i, GetMoveMask(n), add(load(src2 + i), load(src1 + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)        dst[i] += src[i];
}
STATIC void avx512_f32_sub_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
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


STATIC void avx512_f32_add_vec_inplace(const F32PTR src, const F32PTR dst, const int N) {
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        store(dst + i,           add(load(dst + i),          load(src + i)));
        store(dst + i + NF,      add(load(dst + i + NF),     load(src + i + NF)));
        store(dst + i + NF2,     add(load(dst + i + NF2),    load(src + i + NF2)));
        store(dst + i + NF3,     add(load(dst + i + NF3),    load(src + i + NF3)));
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
STATIC void avx512_f32_sub_vec_inplace(const F32PTR src, F32PTR dst, int N) {
   int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {       //dst=dst-src
        store(dst + i,       sub(load(dst + i),      load(src + i)));
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
STATIC void avx512_f32_add_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m512  C = set1(c);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) 
        store(x + i,        add(load( x + i),     C)   ),
        store(x + i + NF,    add(load(x + i + NF),  C) ),
        store(x + i + NF2,   add(load(x + i + NF2), C) ),
        store(x + i + NF3,    add(load(x + i + NF3), C) );
    for (; i < N - (NF - 1); i += NF)
        store(x + i, add(load(x + i), C) );
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i,  GetMoveMask(n), add(load(x + i), C)  );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}
STATIC void avx512_f32_subrev_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m512  C = set1(c);
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

STATIC void avx512_f32_mul_val_inplace(const F32 c, const F32PTR x, const int N) {
    __m512  C = set1(c);
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
STATIC void avx512_f32_mul_vec_inplace(const F32PTR src, F32PTR dst, int N) {
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

STATIC void avx512_f32_mul_vec(const F32PTR src1, const F32PTR src2, F32PTR dst, int N) {
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

STATIC void avx512_f32_add_v_v2_vec_inplace(const F32PTR src, const F32PTR  v, const F32PTR  v2, int N) {

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
        __mmask16 mask = GetMoveMask(n);
        maskstore(v + i,  mask, add(   load(v + i),   load(src + i)  )                           );
        maskstore(v2 + i, mask, add(   load(v2 + i),  mul(load(src + i), load(src + i)))     );
    }
    _mm256_zeroupper();

 
}
STATIC F32 avx512_f32_dot(F32PTR  x, F32PTR y, int N) {

    __m512 r = set0();    
    int i = 0;    
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 r0   = mul(load(x + i),      load(y + i )) ;
        __m512 r1   = mul(load(x + i+NF),   load(y + i+NF));
        __m512 r2   = mul(load(x + i+NF2),  load(y + i+NF2));
        __m512 r3   = mul(load(x + i+NF3),  load(y + i+NF3));
        __m512 s1   = add(r0, r1);
        __m512 s2   = add(r2, r3);
        __m512 s12  = add(s1, s2);
        r=add(r, s12);
    }        
    for (; i < N - (NF - 1); i += NF) 
         r =add(r,  mul(load(x + i), load(y + i))  ); 
    int n = N - i;
    if (n > 0) {
        __mmask16  mask = GetMoveMask(n);
        r = add(r,  mul(   maskload(x + i, mask), maskload(y + i, mask)   )   );
    }    
    F32 sum = f32_hsum(r);
     _mm256_zeroupper();    
    return sum;
}
STATIC F32 avx512_f32_dot2x1(F32PTR  x, F32PTR y, F32PTR v, int N, F32PTR res) {
    __m512 RX = set0();
    __m512 RY = set0();
    int i = 0;    
    for (; i < N - (NF2-1); i += NF2) {
        __m512 x0 = mul(load(x + i),    load(v + i)) ;
        __m512 x1 = mul(load(x + i+NF),  load(v + i+NF));
        __m512 y0 = mul(load(y + i),    load(v + i)) ;
        __m512 y1 = mul(load(y + i+NF),  load(v + i+NF));
   
        RX =add(RX, add(x0, x1)); 
        RY =add(RY, add(y0, y1));
    }        
    for (; i < N - (NF - 1); i += NF) {
         RX = add(RX,  mul(load(x + i), load(v + i))  );
         RY = add(RY, mul(load(y + i), load(v + i))  );
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 mask = GetMoveMask(n);
        __m512    vvec = maskload(v + i, mask); 
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
STATIC void avx512_f32_dot2x2(F32PTR  x, F32PTR y, F32PTR v, F32PTR w, int N, F32PTR res1, F32PTR res2) {
    __m512 RX1 = set0();
    __m512 RX2 = set0();
    __m512 RY1 = set0();
    __m512 RY2 = set0();

    int i = 0;    
    for (; i  < N - (NF2-1); i += NF2) {
        __m512 x0 = mul( load(x+i),      load(v+i)    ) ;
        __m512 x1 = mul( load(x+i+NF),   load(v+i+NF))  ;
        RX1 = add(RX1, add(x0, x1));

        __m512 y0 = mul( load(y+i),      load(v+i)    ) ;
        __m512 y1 = mul( load(y+i+NF),   load(v+i+NF) );      
        RY1  = add(RY1, add(y0, y1));

        __m512 xx0 = mul( load(x+i),       load(w+i)     );
        __m512 xx1 = mul( load(x+i+NF),    load(w+i+NF)  );
        RX2 = add(RX2, add(xx0, xx1));

        __m512 yy0 = mul( load(y + i),      load(w + i)     );
        __m512 yy1 = mul( load(y + i + NF), load(w + i + NF));        
        RY2 = add(RY2, add(yy0, yy1));
    }  
    for (; i < N - (NF - 1); i += NF) {
         RX1 = add(RX1, mul(load(x + i), load(v + i))  );
         RY1 = add(RY1, mul(load(y + i), load(v + i))  );

         RX2 = add(RX2, mul(load(x + i), load(w + i)));
         RY2 = add(RY2, mul(load(y + i), load(w + i)));
    }
  
    int n = N - i;
    if (n > 0) {
        __mmask16 mask = GetMoveMask(n);

        __m512  vvec = maskload(v + i, mask);
        __m512  r0  =  mul(maskload(x + i, mask), vvec);       RX1 = add(RX1, r0);
        __m512  r1  =  mul(maskload(y + i, mask), vvec);       RY1 = add(RY1, r1);

        __m512  wvec = maskload(w + i, mask);
        __m512  rr0  = mul(maskload(x + i, mask), wvec);        RX2 = add(RX2, rr0);
        __m512  rr1  = mul(maskload(y + i, mask), wvec);        RY2 = add(RY2, rr1);
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
STATIC F32 avx512_fma_f32_dot(F32PTR  x, F32PTR y, int N) {

    __m512 r0 = set0();
    __m512 r1 = set0();
    __m512 r2 = set0();
    __m512 r3 = set0();
    
    int i = 0;    
    for (; i < N - (NF4 - 1); i += NF4) {
        r0 = _mm512_fmadd_ps(load(x + i),   load(y + i),   r0) ;
        r1 = _mm512_fmadd_ps(load(x + i+NF), load(y + i+ NF), r1);
        r2 = _mm512_fmadd_ps(load(x + i+ NF2), load(y + i+ NF2),r2);
        r3 = _mm512_fmadd_ps(load(x + i+ NF3), load(y + i+ NF3), r3);
    }        
    for (; i < N - (NF - 1); i += NF) 
         r0 = _mm512_fmadd_ps(load(x + i), load(y + i),  r0  );
    int n = N - i;
    if (n > 0) {
        __mmask16 mask = GetMoveMask(n);
        r1   = _mm512_fmadd_ps( maskload(x + i, mask), maskload(y + i, mask),r1);
    }
    __m512 s1 = add(r0, r1);
    __m512 s2 = add(r2, r3);
    __m512 r  = add(s1, s2);
    float sum = f32_hsum(r);                      
     _mm256_zeroupper();
    return sum;
}

//https://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
//https://stackoverflow.com/questions/54642663/how-sqrt-of-gcc-works-after-compiled-which-method-of-root-is-used-newton-rap
//https://stackoverflow.com/questions/31555260/fast-vectorized-rsqrt-and-reciprocal-with-sse-avx-depending-on-precision
STATIC void avx512_f32_sqrt_vec_inplace(const F32PTR x, const int N)
{
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        //__m512 r1 = load(x + i);          store( x + i,  mul(r1, _mm512_sqrt_ps(r1) ) );
        __m512 r1 = load(x + i);          store(x + i,      _mm512_sqrt_ps(r1) );
        __m512 r2 = load(x + i+NF);       store(x + i+NF,   _mm512_sqrt_ps(r2));
        __m512 r3 = load(x + i + NF2);    store(x + i+NF2,  _mm512_sqrt_ps(r3));
        __m512 r4 = load(x + i + NF3);    store(x + i+NF3,  _mm512_sqrt_ps(r4));
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(x + i, _mm512_sqrt_ps(r));
    }
    
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);
        maskstore(x + i, GetMoveMask(n), _mm512_sqrt_ps(r));
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
//https://stackoverflow.com/questions/1528727/why-is-sse-scalar-sqrtx-slower-than-rsqrtx-x
//https://stackoverflow.com/questions/35885170/handling-zeroes-in-mm256-rsqrt-ps
//https://stackoverflow.com/questions/8924729/using-avx-intrinsics-instead-of-sse-does-not-improve-speed-why
//x0 = vrsqrtps(a)
//x1 = 0.5 * x0 * (3 - (a * x0) * x0)
STATIC void avx512_f32_sqrt_vec(const F32PTR x, const F32PTR y, const int N)
{
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 r1 = load(x + i);          store(y + i,      _mm512_sqrt_ps(r1)  );
        __m512 r2 = load(x + i+NF);       store(y + i+ NF,  _mm512_sqrt_ps(r2)  );
        __m512 r3 = load(x + i + NF2);    store(y + i+ NF2, _mm512_sqrt_ps(r3)  );
        __m512 r4 = load(x + i + NF3);    store(y + i+ NF3, _mm512_sqrt_ps(r4)  );
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(y + i, _mm512_sqrt_ps(r) );
    }
    
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);
        maskstore(y + i, GetMoveMask(n), _mm512_sqrt_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}


#ifdef COMPILER_MSVC
STATIC void avx512_f32_sin_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(x + i,   _mm512_sin_ps(r));
    }    
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm512_sin_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
STATIC void avx512_f32_cos_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(x + i, _mm512_cos_ps(r));
    }    
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm512_cos_ps(r) );
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
STATIC void avx512_f32_sincos_vec_inplace_MSVC(const F32PTR in_outsin, F32PTR outcos, const int N)
{
//Compute the sine and cosine of packed single-precision (32-bit) floating-point elements in a expressed in radians,
//store the sine in dst, and store the cosine into memory at mem_addr.
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(in_outsin + i);
        store(in_outsin + i, _mm512_sincos_ps(outcos+i,r));
    }
    int n = N - i;
    if (n > 0) {
        __mmask16 mask   = GetMoveMask(n);
        __m512  tmpcos;
        __m512  r = load(in_outsin + i);   
        maskstore(in_outsin + i, mask, _mm512_sincos_ps(&tmpcos, r) );
        maskstore(outcos  + i,   mask, tmpcos);
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
STATIC void avx512_f32_pow_vec_inplace_MSVC(F32PTR x, F32 pow, int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    
    __m512 C = set1(pow);

    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);  store(x + i,  _mm512_pow_ps(r,C) );
    }    
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);  maskstore(x + i, GetMoveMask(n), _mm512_pow_ps(r, C));
    }    
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
STATIC void avx512_f32_log_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(x + i, _mm512_log_ps(r));
    }
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm512_log_ps(r));
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
STATIC void avx512_f32_exp_vec_inplace_MSVC(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 r = load(x + i);    store(x + i, _mm512_exp_ps(r));
    }
    int n = N - i;
    if (n > 0) {
        __m512 r = load(x + i);    maskstore(x + i, GetMoveMask(n), _mm512_exp_ps(r));
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;
}
#endif


STATIC void avx512_f32_sin_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
    int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i,  sin512_ps(load(x + i)));    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), sin512_ps(load(x + i)) );
    _mm256_zeroupper();    
}
STATIC void avx512_f32_cos_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, cos512_ps(load(x + i)));    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), cos512_ps(load(x + i)) );
    _mm256_zeroupper();    
}
STATIC void avx512_f32_sincos_vec_inplace(const F32PTR in_outsin, F32PTR outcos, const int N)
{
//Compute the sine and cosine of packed single-precision (32-bit) floating-point elements in a expressed in radians,
//store the sine in dst, and store the cosine into memory at mem_addr.
    int i = 0;
    for (; i < N - (NF - 1); i += NF) 
        sincos512_ps(load(in_outsin + i), in_outsin + i, outcos + i);    
    int n = N - i;
    if (n > 0) {
        __mmask16 mask   = GetMoveMask(n);
        __m512  tmpcos;
        __m512  tmpsin;
        sincos512_ps(load(in_outsin + i), &tmpsin, &tmpcos);
        maskstore(in_outsin + i, mask,  tmpsin);
        maskstore(outcos   + i,   mask, tmpcos);
    }
    _mm256_zeroupper();    
}

static INLINE __m512 pow512(__m512 x, float n) {
    //https://community.intel.com/t5/Intel-ISA-Extensions/AVX512-reciprocal-approximations/td-p/1068416
    __m512 res =  _mm512_mul_ps(log512_ps(x), _mm512_set1_ps(n));
    res = exp512_ps(res);
    return res;
}
static INLINE __m512 pow512_int(__m512 x, int n)
{
    __m512 res = _mm512_set1_ps(1.0f);
    int npositive = n >= 0 ? n : -n;
    while (npositive) {
        if (npositive & 1L)
            res = _mm512_mul_ps(res, x);
        x          = _mm512_mul_ps(x, x);
        npositive  = npositive >> 1;
    }
    if (n < 0) {
        res = _mm512_rcp14_ps(res);
    }
    return res;
}
STATIC void avx512_f32_pow_vec_inplace(F32PTR x, F32 pow, int N)
{
         int   nInteger   = pow;
        float  nRemainder = pow - nInteger;

        if (nRemainder != 0) {
            int i = 0;
            for (; i < N - (NF - 1); i += NF) {
                store(x + i, pow512(load(x+i), pow ) );
            }            
            int n = N - i;
            if (n > 0) 
                maskstore(x + i, GetMoveMask(n), pow512(load(x + i), pow) );              
        }
        else {
            int i = 0;
            for (; i < N - (NF - 1); i += NF) {
                store(x + i, pow512_int(load(x + i), nInteger));
            }
            int n = N - i;
            if (n > 0)
                maskstore(x + i, GetMoveMask(n), pow512_int(load(x + i), nInteger));         
        }
        _mm256_zeroupper();
}
STATIC void avx512_f32_log_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, log512_ps(load(x + i)));    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), log512_ps(load(x + i)) );
    _mm256_zeroupper();    
}
STATIC void avx512_f32_exp_vec_inplace(const F32PTR x, const int N)
{
    //https://github.com/reyoung/avx_mathfun/blob/master/avx_mathfun.h
    //http://gruntthepeon.free.fr/ssemath/sse_mathfun.h
    //https://web.archive.org/web/20191019010229/http://software-lisc.fbk.eu/avx_mathfun/avx_mathfun.h
     int i = 0;
    for (; i < N - (NF - 1); i += NF)    
        store(x + i, exp512_ps(load(x + i)));    
    int n = N - i;
    if (n > 0)
        maskstore(x + i, GetMoveMask(n), exp512_ps(load(x + i)) );
    _mm256_zeroupper();    
}



STATIC void avx512_f32_avgstd(const F32PTR x, int N, F32PTR avg, F32PTR std) {
    __m512  S    = set0();
    __m512  SS   = set0();
    int i = 0;
    for (; i < N - (NF3-1); i += NF4) {
        __m512 s1  =  load(x + i);
        __m512 s2  =  load(x + i+NF);
        __m512 s3  =  load(x + i+NF2);
        __m512 s4   = load(x + i+NF3);

        __m512 s12  = add(s1, s2);                      
        __m512 s34  = add(s3, s4);                ;
        __m512 ss12 = add(mul(s1, s1), mul(s2, s2));
        __m512 ss34 = add(mul(s3, s3), mul(s4, s4));               

        __m512 s1234  = add( s12,s34);
        __m512 ss1234 = add(ss12, ss34);
         S  = add(S,  s1234);
         SS = add(SS, ss1234);
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 s1 = load(x + i);
        S  = add(S,  s1);
        SS = add(SS, mul(s1, s1));
    }
    
    int n = N - i;
    if (n > 0) {
        __mmask16 mask = GetMoveMask(n);
        __m512  s1  =  maskload(x + i,mask);
        S   = add(S, s1);
        SS  = add(SS, mul(s1, s1));
    }
    
    F32 sumx  = f32_hsum(S);
    F32 sumxx = f32_hsum(SS);
     _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
     F64 AVG = sumx / N;
     std[0] = sqrtf((sumxx - AVG * sumx) / (N - 1));// sqrtf((sumxx - AVG * AVG * N) / (N - 1));
     avg[0] = AVG;
}
STATIC void avx512_f32_sx_sxx_to_avgstd_inplace(F32PTR SX, F32PTR SXX, I32 Nsample, F32 scale, F32 offset, int N) {
    /*
                   r_ippsMulC_32f_I(inv_sample* sd, resultChain.tY, N);
                   r_ippsMul_32f(resultChain.tY, resultChain.tY, MEMBUF1, N);
                   r_ippsMulC_32f_I(inv_sample* sd* sd, resultChain.tSD, N);
                   r_ippsSub_32f_I(MEMBUF1, resultChain.tSD, N);
                   r_ippsSqrt_32f_I(resultChain.tSD, N);
    */

    __m512 invsample_scale  = set1(1. / (F64) Nsample * scale);
    __m512 invsample_scale2 = set1(1. / (F64) Nsample * scale*scale);
    __m512 offset_vec       = set1(offset);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        __m512 avg   = mul(load(SX + i),    invsample_scale);
        __m512 sxx_n = mul(load(SXX + i),   invsample_scale2);

        __m512 sd2  = sub(sxx_n, mul(avg, avg));
        __m512 sd   = mul(_mm512_rsqrt14_ps(sd2), sd2); //sqrt(sd2)

        avg = add(avg, offset_vec);
        store(SXX + i, sd);
        store(SX  + i, avg); 

    }

    int n = N - i;
    if (n >0){
        __mmask16 mask  = GetMoveMask(n);
        __m512  avg   = mul(maskload(SX + i,mask), invsample_scale);
        __m512  sxx_n = mul(maskload(SXX + i,mask), invsample_scale2);

        __m512 sd2 = sub(sxx_n, mul(avg, avg));
        __m512 sd = mul(_mm512_rsqrt14_ps(sd2), sd2); //sqrt(sd2)

        avg = add(avg, offset_vec);
        maskstore(SXX + i, mask,sd);
        maskstore(SX + i, mask, avg);

    }

}


#ifdef OS_WIN64 
    #define WIN32_LEAN_AND_MEAN
    #include "windows.h"
#endif
/*
I32 avx512_f32_maxidx1(const F32PTR x, const int N, F32PTR val) {
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


    //__m512 r1 = set0();    
    __m512  maxVec      = load(x);
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);
    __m512i maxIdx      = _mm256_cvtepu8_epi32(_maskIdx);
    __m512i EIGHT       = set1_i32(8);

    __m512i idx         = maxIdx;

    int i = 8;    
    for (; i < N - (NF - 1); i += NF) {        
        idx   = addi32(idx, EIGHT);

        __m512  vec     = load(x + i);
        __m512  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);

         maxVec  = _mm256_blendv_ps(maxVec,  vec, cmpmask); //o: take a, 1: take from vec        
         maxIdx  = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));        
    }        
 
    int n = N - i;
    if (n > 0) {
        __mmask16 mask    = GetMoveMask(n);
        __m512  vec     = maskload(x + i,mask);
        idx = addi32(idx, EIGHT);
        
        idx = _mm256_blendv_epi8(maxIdx,idx, mask);
        vec = _mm256_blendv_ps(  maxVec,vec, _mm256_castsi256_ps(mask)); //o: take a, 1: take from vec        

        __m512  cmpmask = _mm256_cmp_ps(maxVec, vec, _CMP_LE_OQ);
        maxVec = _mm256_blendv_ps(maxVec, vec, cmpmask); //o: take a, 1: take from vec        
        maxIdx = _mm256_blendv_epi8(maxIdx, idx, _mm256_castps_si256(cmpmask));
    }

    __m128 vlow     = _mm256_castps256_ps128(maxVec);
    __m128 vhigh    = _mm256_extractf128_ps(maxVec, 1); // high 128    
    __m128 tmp      = _mm_max_ps(vlow, vhigh);
    tmp = _mm_max_ps(tmp, _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(1, 0, 3, 2) ) );
    tmp = _mm_max_ps(tmp, _mm_shuffle_ps(tmp, tmp, _MM_SHUFFLE(0, 0, 0, 1))  );    
    F32 max = _mm_cvtss_f32(tmp);
    __m512  max256 = set1(max);

    __m512   vcmp      = _mm256_cmp_ps(maxVec, max256, _CMP_EQ_OQ);
    uint32_t finalmask = _mm256_movemask_epi8(_mm256_castps_si256(vcmp));
 
    DWORD index;
    _BitScanReverse(&index, finalmask);     //tmp = _mm_packs_epi16(_mm_packs_epi32(_mm_cmpeq_epi32(_mm_set1_epi32(max), lo),    _mm_cmpeq_epi32(_mm_set1_epi32(max), hi)),       _mm_setzero_si128());
    union {
        I32    idx[8];
        __m512i dummy;
    }ID;
    _mm256_store_si256(&ID.dummy, maxIdx);
    _mm256_zeroupper();     
    val[0] = max;
    return ID.idx[index>>2] ;
}
*/

STATIC I32 avx512_f32_maxidx(const F32PTR x, const int N, F32PTR val) {
    //https://stackoverflow.com/questions/23590610/find-index-of-maximum-element-in-x86-simd-vector    
      
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);    
            _maskIdx    = _mm_insert_epi64(_maskIdx, 0x0F0E0D0C0B0A0908, 1);
    __m512i idx         = _mm512_cvtepu8_epi32(_maskIdx);
    __m512i Sixteen     = set1_i32(16);

    __m512   maxVec = set1(x[0]);
    __m512i  maxIdx = set0i();

    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {            
        __m512     vec     = load(x + i);
        __mmask16  cmpmask = _mm512_cmp_ps_mask(maxVec, vec, _CMP_LE_OQ);

         maxVec  = _mm512_mask_blend_ps( cmpmask,   maxVec, vec ); //o: take a, 1: take from vec        
         maxIdx  = _mm512_mask_blend_epi32(cmpmask, maxIdx, idx );

         idx = addi32(idx, Sixteen);
    }         
    int n = N - i;
    if (n > 0) {
        __mmask16 mask    = GetMoveMask(n);
        __m512    vec     = maskload(x + i,mask);        
        
        vec = _mm512_mask_blend_ps(mask, maxVec, vec); //o: take a, 1: take from vec        
        idx = _mm512_mask_blend_epi32(mask, maxIdx,idx);
        
        __mmask16  cmpmask = _mm512_cmp_ps_mask(maxVec, vec, _CMP_LE_OQ);
         maxVec  = _mm512_mask_blend_ps( cmpmask,   maxVec, vec ); //o: take a, 1: take from vec        
         maxIdx  = _mm512_mask_blend_epi32(cmpmask, maxIdx, idx );
    }    
    __m256  low256     = _mm512_castps512_ps256(maxVec);
    __m256  high256    = _mm512_extractf32x8_ps(maxVec, 1); // high 128 
    __m256i idxlow256  = _mm512_castsi512_si256(maxIdx);
    __m256i idxhigh256 = _mm512_extracti32x8_epi32(maxIdx, 1); // high 128

    __m256 cmpmask    = _mm256_cmp_ps(low256, high256, _CMP_LE_OQ);      
    __m256  maxVec256 = _mm256_blendv_ps(low256, high256, cmpmask); //o: take a, 1: take from vec     
    __m256i maxIdx256 = _mm256_blendv_epi8(idxlow256, idxhigh256, _mm256_castps_si256(cmpmask)); //o: take a, 1: take from vec     

    //////////////////////////////////////////////////////////////////

    __m128  vlow     = _mm256_castps256_ps128(maxVec256);
    __m128  vhigh    = _mm256_extractf128_ps(maxVec256, 1); // high 128
    __m128i idxlow   = _mm256_castsi256_si128(maxIdx256);
    __m128i idxhigh  = _mm256_extracti128_si256(maxIdx256, 1); // high 128

    __m128  cmpmask128  = _mm_cmp_ps(vlow, vhigh, _CMP_LE_OQ);  
    __m128  max128      = _mm_blendv_ps(vlow, vhigh, cmpmask128); //o: take a, 1: take from vec     
    __m128i idx128      = _mm_blendv_epi8(idxlow, idxhigh, _mm_castps_si128(cmpmask128)); //o: take a, 1: take from vec     

     F32 VAL[4];
     union {
        I32     idx[4];
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
STATIC I32 avx512_f32_minidx(const F32PTR x, const int N, F32PTR val) {
    
    //https://stackoverflow.com/questions/23590610/find-index-of-maximum-element-in-x86-simd-vector    
      
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);    
            _maskIdx    = _mm_insert_epi64(_maskIdx, 0x0F0E0D0C0B0A0908, 1);
    __m512i idx         = _mm512_cvtepu8_epi32(_maskIdx);
    __m512i Sixteen     = set1_i32(16);

    __m512   minVec = set1(x[0]);
    __m512i  minIdx = set0i();

    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {            
        __m512     vec     = load(x + i);
        __mmask16  cmpmask = _mm512_cmp_ps_mask(minVec, vec, _CMP_LE_OQ);

        minVec = _mm512_mask_blend_ps( cmpmask, minVec, vec ); //o: take a, 1: take from vec        
        minIdx = _mm512_mask_blend_epi32(cmpmask, minIdx, idx );

        idx    = addi32(idx, Sixteen);
    }         
    int n = N - i;
    if (n > 0) {
        __mmask16 mask    = GetMoveMask(n);
        __m512    vec     = maskload(x + i,mask);        
        
        vec = _mm512_mask_blend_ps(  mask,  minVec, vec); //o: take minVec, 1: take from vec        
        idx = _mm512_mask_blend_epi32(mask, minIdx, idx);

        __mmask16  cmpmask = _mm512_cmp_ps_mask(minVec, vec, _CMP_LE_OQ);
         minVec  = _mm512_mask_blend_ps( cmpmask,   minVec, vec ); //o: take a, 1: take from vec        
         minIdx  = _mm512_mask_blend_epi32(cmpmask, minIdx, idx );
    }    
    __m256 low256      = _mm512_castps512_ps256(minVec);
    __m256 high256     = _mm512_extractf32x8_ps(minVec, 1); // high 256
    __m256i idxlow256  = _mm512_castsi512_si256(minIdx);
    __m256i idxhigh256 = _mm512_extracti32x8_epi32(minIdx, 1); // high 128

    __m256 cmpmask    = _mm256_cmp_ps(low256, high256, _CMP_LE_OQ);  
    __m256  minVec256  = _mm256_blendv_ps(low256, high256, cmpmask); //o: take a, 1: take from vec     
    __m256i minIdx256 = _mm256_blendv_epi8(idxlow256, idxhigh256, _mm256_castps_si256(cmpmask)); //o: take a, 1: take from vec     

    //////////////////////////////////////////////////////////////////

    __m128 vlow        = _mm256_castps256_ps128(minVec256);
    __m128 vhigh       = _mm256_extractf128_ps(minVec256, 1); // high 128
    __m128i idxlow     = _mm256_castsi256_si128(minIdx256);
    __m128i idxhigh    = _mm256_extracti128_si256(minIdx256, 1); // high 128

    __m128 cmpmask128  = _mm_cmp_ps(vlow, vhigh, _CMP_LE_OQ);  
    __m128 max128      = _mm_blendv_ps(vlow, vhigh, cmpmask128); //o: take a, 1: take from vec     
    __m128i idx128     = _mm_blendv_epi8(idxlow, idxhigh, _mm_castps_si128(cmpmask128)); //o: take a, 1: take from vec     

     F32 VAL[4];
     union {
        I32     idx[4];
        __m128i dummy;
     }ID;
     _mm_storeu_ps(VAL,          max128);
     _mm_storeu_si128(&ID.dummy, idx128);
     _mm256_zeroupper();
     I32 i1 = VAL[0] < VAL[1] ? 0 : 1;
     I32 i2 = VAL[2] < VAL[3] ? 2 : 3;
     I32 IDX = VAL[i1] < VAL[i2] ? i1 : i2;

     val[0] = VAL[IDX];
    return ID.idx[IDX];
}

STATIC void avx512_f32_diff_back(const F32PTR  x, F32PTR result, const int N) {
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
STATIC void avx512_f32_seq(F32PTR p, F32 x0, F32 dX, int N) {
    
    __m512 X;
    {    
        /*
        __m128i tmp1 = _mm_cvtsi64_si128(0x0706050403020100);
        __m128i tmp2 = _mm_cvtsi64_si128(0x0F0E0D0C0B0A0908);
        __m256i tmp2561 = _mm256_cvtepu8_epi32(tmp1);
        __m256i tmp2562 = _mm256_cvtepu8_epi32(tmp2);
        __m512i v1to15 = _mm512_inserti64x4(_mm512_castsi256_si512(tmp2561), tmp2562, 1);
        __m512  vf1to15 = _mm512_cvtepi32_ps(v1to15);
        __m512  DX = mul(vf1to15, set1(dX));   
        */
        __m128i tmp;
        tmp  = _mm_cvtsi64_si128(0x0706050403020100);
        tmp  = _mm_insert_epi64(tmp, 0x0F0E0D0C0B0A0908, 1);        
        __m512i tmp512 = _mm512_cvtepu8_epi32(tmp);
        __m512  v1to15 = _mm512_cvtepi32_ps(tmp512);
        __m512  DX     = mul(v1to15, set1(dX));
        X = add(set1(x0),DX);        
    }    
    __m512 dX16 = set1(dX*16)  ;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        store(p+i, X);  
        X = add(X, dX16);        
    } 
    int n = N - i;
    if (n > 0) {
        maskstore(p + i, GetMoveMask(n), X);
    }     
    _mm256_zeroupper();    
}
STATIC void avx512_i32_seq(I32PTR p, I32 x0, I32 dX, int N) {
    
    __m512i X;
    {    
        /*
        __m128i tmp1 = _mm_cvtsi64_si128(0x0706050403020100);
        __m128i tmp2 = _mm_cvtsi64_si128(0x0F0E0D0C0B0A0908);
        __m256i tmp2561 = _mm256_cvtepu8_epi32(tmp1);
        __m256i tmp2562 = _mm256_cvtepu8_epi32(tmp2);
        __m512i v1to15 = _mm512_inserti64x4(_mm512_castsi256_si512(tmp2561), tmp2562, 1);
        __m512  vf1to15 = _mm512_cvtepi32_ps(v1to15);
        __m512  DX = mul(vf1to15, set1(dX));   
        */
        __m128i tmp;
        tmp  = _mm_cvtsi64_si128(0x0706050403020100);
        tmp  = _mm_insert_epi64(tmp, 0x0F0E0D0C0B0A0908, 1);        
        __m512i v1to15 = _mm512_cvtepu8_epi32(tmp);
        //__m512  v1to15 = _mm512_cvtepi32_ps(tmp512);
        __m512i  DX     = _mm512_mullo_epi32(v1to15, set1_i32(dX));
        X = addi32(set1_i32(x0),DX);        
    }    
    __m512i dX16 = set1_i32(dX*16)  ;
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {
        storei(p+i, X);  
        X = addi32(X, dX16);        
    } 
    int n = N - i;
    if (n > 0) {
        maskstorei(p + i, GetMoveMask(n), X);
    }     
    _mm256_zeroupper();    
}
STATIC void avx512_f32_to_f64_inplace(F32PTR data32, int N) {

    F64PTR data64 = data32;
    int i = N - 16;
    for (; i >= 0; i -= 16) {
        __m512 r        = load(data32 + i);
        __m256 lowps    = _mm512_castps512_ps256(r);
        __m256 highps   = _mm512_extractf32x8_ps(r, 1);         // high 128
        _mm512_storeu_pd(data64 + i,    _mm512_cvtps_pd(lowps));
        _mm512_storeu_pd(data64 + i+8,  _mm512_cvtps_pd(highps));
    }
    i = i + 16; // i is the remaining number and is also the last index processed above
    if (i >= 8) {
        i = i - 8;
        __m256 r        = _mm256_loadu_ps(data32 + i);
        __m128 lowps    = _mm256_castps256_ps128(r);
        __m128 highps   = _mm256_extractf128_ps(r, 1);         // high 128
        _mm256_storeu_pd(data64 + i,   _mm256_cvtps_pd(lowps));
        _mm256_storeu_pd(data64 + i+4, _mm256_cvtps_pd(highps));
    }
    if (i >= 4) {
        i = i - 4;
        __m128 r = _mm_loadu_ps(data32 + i);        
        _mm256_storeu_pd(data64 + i, _mm256_cvtps_pd(r));
    }
    i = i - 1;
    for (; i >= 0; --i)   data64[i] = data32[i];    
    _mm256_zeroupper();    
}
STATIC void avx512_f64_to_f32_inplace(F64PTR data64, int N) {

    F32PTR data32 = data64;
    int i = 0;
    for (; i < N-15; i += 16) {    
        __m256  R1 = _mm512_cvtpd_ps( _mm512_loadu_pd(data64 + i  )   );
        __m256  R2 = _mm512_cvtpd_ps( _mm512_loadu_pd(data64 + i+8)   );
        _mm256_storeu_ps(data32 + i,   R1);
        _mm256_storeu_ps(data32 + i+8, R2);
    }
    for (; i < N - 7; i += 8) {
        __m128  R1 = _mm256_cvtpd_ps( _mm256_loadu_pd(data64 + i)     );
        __m128  R2 = _mm256_cvtpd_ps( _mm256_loadu_pd(data64 + i + 4) );
        _mm_storeu_ps(data32 + i,       R1);
        _mm_storeu_ps(data32 + i + 4,   R2);
    }
    for (; i < N - 3; i +=4) {
        __m128  R1 = _mm256_cvtpd_ps(_mm256_loadu_pd(data64 + i));
        _mm_storeu_ps(data32 + i, R1);
    }    
    _mm256_zeroupper();
    for (; i <N; ++i)      data32[i] = data64[i];         
}
STATIC void avx512_i32_to_f32_scaleby_inplace(I32PTR x, int N, F32 scale) {
    __m512  C = set1(scale);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(x + i,       mul(_mm512_cvtepi32_ps(loadi(x + i    )), C)),
        store(x + i + NF,  mul(_mm512_cvtepi32_ps(loadi(x + i+NF )), C)),
        store(x + i + NF2, mul(_mm512_cvtepi32_ps(loadi(x + i+NF2)), C)),
        store(x + i + NF3, mul(_mm512_cvtepi32_ps(loadi(x + i+NF3)), C));
    for (; i < N - (NF - 1); i += NF)
        store(x + i, mul(_mm512_cvtepi32_ps(loadi(x + i)), C));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(x + i, GetMoveMask(n), mul(_mm512_cvtepi32_ps(loadi(x + i)), C));
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

STATIC void avx512_i32_increment_bycond_inplace(I32PTR x,  F32PTR cond, int N) {
    __m512   C0 = set0();
    __m512i  C1 = set1_i32(1);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {            
        __mmask16 mask   = _mm512_cmp_ps_mask( load(cond+i), C0, _CMP_GT_OQ) ;
        __m512i   vec    =  maskloadi(x + i, mask);
        vec = addi32(vec, C1);
        maskstorei(x + i, mask, vec);
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __mmask16 loadmask = GetMoveMask(n);
        __mmask16 mask     = _mm512_cmp_ps_mask(maskload(cond + i,loadmask), C0, _CMP_GT_OQ);
         mask = mask & loadmask;
        __m512i vec = maskloadi(x + i, mask);
        vec = addi32(vec, C1);
        maskstorei(x + i, mask, vec);
    }
    _mm256_zeroupper();
}
STATIC void avx512_i32_increment_vec2_bycond_inplace(I32PTR x, I32PTR y, F32PTR cond, int N) {
    __m512   C0        = set0();
    __m512   Ceplison1 = set1(1e-10);
    __m512   Ceplison2 = set1(-1e-10);
    __m512i  C1        = set1_i32(1);
    int i = 0;
    for (; i < N - (NF - 1); i += NF) {   
        __m512    data    = load(cond + i);
        __mmask16 mask1   = _mm512_cmp_ps_mask(data, Ceplison1, _CMP_GT_OQ) ;
        __m512i   vec1    =  maskloadi(x + i, mask1);
        vec1 = addi32(vec1, C1);    maskstorei(x + i, mask1, vec1);

        __mmask16 mask2 = _mm512_cmp_ps_mask(Ceplison1, data, _CMP_GT_OQ) & _mm512_cmp_ps_mask(data, Ceplison2, _CMP_GT_OQ);
        __m512i   vec2  = maskloadi(y + i, mask2);
        vec2 = addi32(vec2, C1);    maskstorei(y + i, mask2, vec2);
    }
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
    {
        __mmask16 loadmask = GetMoveMask(n);
        __m512    data     = maskload(cond + i, loadmask);
        __mmask16 mask1     = _mm512_cmp_ps_mask(data, Ceplison1, _CMP_GT_OQ);
         mask1 = mask1 & loadmask;
        __m512i vec1 = maskloadi(x + i, mask1);
        vec1 = addi32(vec1, C1);   maskstorei(x + i, mask1, vec1);

        __mmask16 mask2 = _mm512_cmp_ps_mask(Ceplison1, data, _CMP_GT_OQ) & _mm512_cmp_ps_mask(data, Ceplison2, _CMP_GT_OQ);
         mask2 = mask2 & loadmask;
        __m512i   vec2  = maskloadi(y + i, mask2);
        vec2 = addi32(vec2, C1);    maskstorei(y + i, mask2, vec2);

    }
    _mm256_zeroupper();
}
STATIC I32  avx512_i08_sum_binvec(U08PTR binvec, I32 N) {
    //stackoverflow.com/questions/36998538/fastest-way-to-horizontally-sum-sse-unsigned-byte-vector
	I32   SUM = 0;
	I32	  i   = 0;
	//|8bytees|8bytees|8bytees|8bytees|
    __m512i r   = set0i();
    I32  ite_read_for_sad = 0;
	for (; i < N - (64*4-1); i+=64*4) {
        __m512i sum1  = _mm512_adds_epu16(loadi(binvec + i),    loadi(binvec  + i + 64));
        __m512i sum2  = _mm512_adds_epu16(loadi(binvec + i+128), loadi(binvec + i + 192));
        __m512i sum12 = _mm512_adds_epu16(sum1, sum2);
                 r    = _mm512_adds_epu16(r, sum12);
                ++ite_read_for_sad;
                if (ite_read_for_sad > 25) { /*25*4=100  100*2=200*/
                    ite_read_for_sad = 0;
                    
                    __m256i low     = _mm512_castsi512_si256(r);
                    __m256i high    = _mm512_extracti32x8_epi32(r, 1);

                    //https://stackoverflow.com/questions/14962349/adding-two-signed-or-unsigned-integers
                    //https://stackoverflow.com/questions/25609091/what-happens-when-i-mix-signed-and-unsigned-types
                    __m256i R       = _mm256_add_epi32(low, high);

                    //stackoverflow.com/questions/36998538/fastest-way-to-horizontally-sum-sse-unsigned-byte-vector
                    __m256i  sum   = _mm256_sad_epu8(R, _mm256_setzero_si256());  

                    __m128i sumLow  = _mm256_castsi256_si128(sum);
                    __m128i sumHigh = _mm256_extracti128_si256(sum, 1); // high 128

                    __m128i sum128  = _mm_adds_epu16(sumLow, sumHigh);
                    I32    localSum = _mm_extract_epi16(sum128, 0) + _mm_extract_epi16(sum128, 4);
                    SUM = SUM + localSum;
                    r = set0i();
                }
	}
    for (; i < N - (64 - 1); i += 64) {
        r = _mm512_adds_epu16(r, loadi(binvec + i));
    }
    for (; i < N - (32 - 1); i += 32) {
            __m256i tmp256   = _mm256_loadu_si256(binvec + i);
            __m512i tmp512   = set0i();
            tmp512 = _mm512_inserti64x4(tmp512, tmp256, 0);
            r = _mm512_adds_epu16(r, tmp512);
    }   
    {
        __m256i low = _mm512_castsi512_si256(r);
        __m256i high = _mm512_extracti32x8_epi32(r, 1);
        __m256i R = _mm256_add_epi32(low, high);

        R = _mm256_sad_epu8(R, _mm256_setzero_si256());
        
        __m128i sumLow  = _mm256_castsi256_si128(R);
        __m128i sumHigh = _mm256_extracti128_si256(R, 1); // high 128
        __m128i sum     = _mm_adds_epu16(sumLow, sumHigh);
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

STATIC void avx512_f32_gemm_XtY2x1(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{ 
	for (int col = 0; col < N; ++col) {
		int row = 0;
		for (; row < M - (2 - 1); row += 2) 
            C[row + 1] = avx512_f32_dot2x1(A + row * lda, A + (row + 1) * lda, B, K, C + row);
			//f32_dot3x1(A + ROW * lda, A + (ROW + 1) * lda, A + (ROW + 2) * lda, B, K, C + ROW);			
			//C[row] = f32_dot(A + row * lda, B, K);		
		if (row < M) 
			C[row] = avx512_f32_dot(A + row * lda, B, K);
		B += ldb;
		C += ldc;
	}
 
}
  
STATIC void avx512_f32_gemm_XtY2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{
	int COL;
 	for (COL = 0; COL < N-(2-1); COL+=2) {
		int ROW;
		for (ROW=0; ROW < M-(2-1); ROW+=2) {			
            avx512_f32_dot2x2(A+ROW*lda, A + (ROW+1)*lda, B, B+ldb, K, C+ROW, C +ldc+ROW);
		}
	    if(ROW<M){
			*(C + ldc + ROW)= avx512_f32_dot2x1(B, B + ldb, A + ROW* lda, K, C+ROW);
		}
		B += ldb+ldb;
		C += ldc+ldc;
	}

	if (COL < N) {
		int ROW;
		for (ROW = 0; ROW < M - (2-1); ROW += 2) 
			*(C+ROW+1)= avx512_f32_dot2x1(A + ROW * lda, A + (ROW + 1) * lda, B, K, C + ROW);
		if (ROW < M) 
			C[ROW]   = avx512_f32_dot(A + ROW * lda, B,K);
	}
}
////////////////////////////////////////////////////////////////////////////////////////// 


////////////////////////////////////////////////////////////////////////////////////////// 
static INLINE void GetOneRowFromMatrix(F32PTR row, F32PTR mat, I32 lda, I32 K,I32 kRemainder, __m512i offset, __mmask16 maskRemainder) 
{
    int i = 0;
    for (; i < K-(NF-1); i += NF)  
        store(row + i, _mm512_i32gather_ps( offset, mat + i * lda, 4));
    if (kRemainder) {
        __m512 subrow = _mm512_mask_i32gather_ps(set0(), maskRemainder, offset, mat + i * lda, 4);
        maskstore(row+i, maskRemainder, subrow);
    }
}

static INLINE void GetTwoRowFromMatrix(F32PTR row1, F32PTR row2, F32PTR mat, I32 lda, I32 K, I32 kRemainder, __m512i offset, __mmask16 maskRemainder)
{
    int i = 0;
    for (; i < K - (NF - 1); i += NF) {
        store(row1 + i, _mm512_i32gather_ps(offset,mat   + i * lda,  4));
        store(row2 + i, _mm512_i32gather_ps(offset,mat+1 + i * lda,  4));
    }
    if (kRemainder) {
        __m512 subrow1 = _mm512_mask_i32gather_ps(set0(), maskRemainder, offset, mat   + i * lda,  4);
        __m512 subrow2 = _mm512_mask_i32gather_ps(set0(), maskRemainder, offset, mat+1 + i * lda,  4);
        maskstore(row1 + i, maskRemainder, subrow1);
        maskstore(row2 + i, maskRemainder, subrow2);
    }
}
static INLINE __m512i GetRowOffset( I32 lda ) {     
    __m128i tmp = _mm_cvtsi64_si128(0x0706050403020100);
            tmp = _mm_insert_epi64(tmp, 0x0F0E0D0C0B0A0908, 1);
    __m512i vec0to15  = _mm512_cvtepu8_epi32(tmp);
    __m512i vecLDA    = set1_i32(lda);
    __m512i offset    = _mm512_mullo_epi32(vec0to15, vecLDA);
    return offset;
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
    __mmask16 mask = GetMoveMask(nRemainder);
    __m512i offset;
    {
        __m512i vecLDA = set1_i32(lda);
        offset = _mm256_cvtepu8_epi32(_mm_cvtsi64_si128(0x0706050403020100));
        offset = _mm256_mul_epi32(offset, vecLDA);
    }

    for (int row = 0; row < M; row++) {
 
        //  Get one row
 
        int i = 0;
        for (; i < K - 7; i += 8)
            store(Xrow + i, _mm256_i32gather_ps(A + row + i * lda, offset, 4));
        if (nRemainder) {
            __m512 subrow = _mm256_mask_i32gather_ps(set0(), A + row + i * lda, offset, _mm256_castsi256_ps(mask), 4);
            maskstore(Xrow + i, mask, subrow);
        }
   
        // End of getting one row
  

        F32PTR Crow = C + row;
        int col = 0;
        for (; col < N - 1; col += 2) {
            *(Crow + ldc * (col + 1)) = avx512_f32_dot2x1(B + col * ldb, B + col * ldb + ldb, Xrow, K, Crow + ldc * col);
        }
        if (col < N) {
            *(Crow + ldc * col) = avx512_f32_dot(B + col * ldb, Xrow, K);
        }

    }

}
*/
STATIC void avx512_f32_gemm_XY1x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{
    // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
    // MSVC seems not to support it

    /*
    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
    Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
    */
    F32       XROW_FIXEXD[4096];
    F32PTR    Xrow         = (K<=4096)? XROW_FIXEXD: alloca(K*sizeof(F32));
    I32       nRemainder = K % 16L;
    __mmask16 mask       =  GetMoveMask(nRemainder);
    __m512i   offset      =  GetRowOffset(lda);
  

    for (int row = 0; row < M; row++) {
    //  Get one row of A
        GetOneRowFromMatrix(Xrow, A + row, lda, K, nRemainder, offset, mask);        
        F32PTR Crow = C + row;
        int    col  = 0;
        for (; col < N - 1; col += 2) 
            *(Crow + ldc * (col+1))=avx512_f32_dot2x1(B+col*ldb,B+col*ldb+ldb,Xrow,K, Crow+ldc*col);        
        if (col < N) 
            *(Crow+ldc*col)=avx512_f32_dot(B + col * ldb, Xrow, K);
    }  

}
STATIC void avx512_f32_gemm_XY2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{ // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
  // MSVC seems not to support it
  //    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
  //   Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
   
    F32     XROW_FIXEXD[4096*2];
    F32PTR  Xrow1 = (2*K <= 2*4096) ? XROW_FIXEXD : alloca(2*K*sizeof(F32));
    F32PTR  Xrow2 = Xrow1 + K;
    I32     nRemainder = K % 16L;
    __mmask16 mask   =  GetMoveMask(nRemainder);
    __m512i   offset = GetRowOffset(lda);

    int row = 0;
    for (; row < M-1; row+=2) {
        //  Get one row of A
        GetTwoRowFromMatrix(Xrow1,Xrow2, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;        
        int    col = 0;
        for (; col < N - 1; col += 2)
            avx512_f32_dot2x2(Xrow1, Xrow2, B + col * ldb, B + col * ldb + ldb, K, Crow + ldc * col, Crow+ldc*(col+1));
        if (col < N)
            *(Crow + ldc * col + 1) = avx512_f32_dot2x1(Xrow1,Xrow2,B + col * ldb,  K, Crow + ldc * col);;
   
    }

	if (row < M) {
        //  Get one row of A
        GetOneRowFromMatrix(Xrow1, A + row, lda, K, nRemainder, offset, mask);
        F32PTR Crow = C + row;
        int    col = 0;
        for (; col < N - 1; col += 2)
             *(Crow + ldc * (col+1))=avx512_f32_dot2x1(B+col*ldb,B+(col+1)*ldb, Xrow1,K, Crow+col*ldc);              
        if (col < N)
            *(Crow + ldc * col ) = avx512_f32_dot(Xrow1, B+col*ldb,  K);            
	}
}

STATIC void avx512_f32_gemm_XtYt2x2(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb,  F32PTR C, int ldc)
{ // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
  // MSVC seems not to support it
  //    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
  //   Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
   
    F32     YROW_FIXEXD[4096*2];
    F32PTR  Yrow1 = (K <= 4096) ? YROW_FIXEXD : alloca(2*K * sizeof(F32));
    F32PTR  Yrow2 = Yrow1 + K;
    I32     nRemainder = K % 16L;
    __mmask16 mask     =  GetMoveMask(nRemainder);
    __m512i   Boffset  = GetRowOffset(ldb);
 
    int col = 0;
    for (; col < N-1; col+=2) {
        //  Get one row of B
        GetTwoRowFromMatrix(Yrow1,Yrow2, B + col, ldb, K, nRemainder, Boffset, mask);        
        int    row = 0;
        for (; row < M - 1; row += 2)
            avx512_f32_dot2x2(A+row*lda, A + (row+1)*lda, Yrow1, Yrow2, K, C+row, C+ldc+row);
        if (row < M)
            *(C + ldc + row) = avx512_f32_dot2x1(Yrow1, Yrow2, A + row * lda,  K, C+row);
   
        C += ldc+ldc;
    }
	if (col < N) {
        //  Get one row of B
        GetOneRowFromMatrix(Yrow1, B + col, ldb, K, nRemainder, Boffset, mask);        
        int    row = 0;
        for (; row < M - 1; row += 2)
             *(C+row + 1)=avx512_f32_dot2x1(A+row*lda,A+(row+1)*lda, Yrow1,K, C+row);              
        if (row < M)
            *(C + row) = avx512_f32_dot(A+row*lda, Yrow1, K);          
	}
}


/////////////////////////////////////////////////////////////////////////////
static INLINE F32 __avx512_f32_dot2x1stride(F32PTR  x, F32PTR y, F32PTR v, int N, __mmask16 mask, __m512i voffset, I32 vstride,  F32PTR res) {
    __m512 RX = set0();
    __m512 RY = set0();
    int i = 0;
    for (; i < N - (NF2 - 1); i += NF2) {
        __m512 v0 = _mm512_i32gather_ps(voffset, v +     i  * vstride,    4);
        __m512 v1 = _mm512_i32gather_ps(voffset, v + (i+NF) * vstride,    4);
        __m512 x0 = mul(load(x + i),     v0);
        __m512 x1 = mul(load(x + i + NF), v1);
        __m512 y0 = mul(load(y + i),     v0);
        __m512 y1 = mul(load(y + i + NF), v1);

        RX = add(RX, add(x0, x1));
        RY = add(RY, add(y0, y1));
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 V = _mm512_i32gather_ps(voffset, v + i * vstride, 4);
        RX = add(RX, mul(load(x + i), V));
        RY = add(RY, mul(load(y + i), V));
    }

    int n = N - i;
    if (n > 0) {
        __m512  vec = _mm512_mask_i32gather_ps(set0(), mask, voffset, v + i*vstride, 4);
        RX = add(RX, mul(maskload(x + i, mask), vec));
        RY = add(RY, mul(maskload(y + i, mask), vec));
    }

    float sumX = f32_hsum(RX);
    float sumY = f32_hsum(RY);

    _mm256_zeroupper();
    res[0] = sumX;
    return sumY;
}
static INLINE F32 __avx512_f32_dot_stride(F32PTR  x, F32PTR y, int N, __mmask16 mask, __m512i yoffset, I32 ystride) {

    __m512 r = set0();    
    int i = 0;    
    for (; i < N - (NF - 1); i += NF) {
        __m512 Y = _mm512_i32gather_ps(yoffset, y + i * ystride, 4);
        r = add(  r, mul( load(x + i), Y ) );
    }    
    int n = N - i;
    if (n > 0) {        
        __m512 Y = _mm512_mask_i32gather_ps(set0(), mask, yoffset, y + i * ystride, 4);
        r = add(r,  mul(   maskload(x + i, mask), Y  )   );
    }    
    F32 sum = f32_hsum(r);
     _mm256_zeroupper();    
    return sum;
}
STATIC void avx512_f32_gemm_XYt2x1(int M, int N, int K, F32PTR A, int lda, F32PTR B, int ldb, F32PTR C, int ldc)
{
    // F32 Xrow[K] : variable-length array defined in C99 but optional in C11
    // MSVC seems not to support it
    /*
    https://stackoverflow.com/questions/14666665/trying-to-understand-gcc-option-fomit-frame-pointer
    Typically the compiler can keep track of stack depth on its own and does not need a frame pointer. The exception is if the function uses alloca which moves the stack pointer by a variable amount. Frame pointer omission does make debugging significantly harder. Local variables are harder to locate and stack traces are much harder to reconstruct without a frame pointer to help out. Also, accessing parameters can get more expensive since they are far away from the top of the stack and may require more expensive addressing mode
    */
    F32     XROW_FIXEXD[4096*2];
    F32PTR  Xrow1  = (K <= 4096) ? XROW_FIXEXD : alloca( 2*K * sizeof(F32) );
    F32PTR  Xrow2  = Xrow1 + K;
    I32     nRemainder   = K % 16L;
    __mmask16 mask       = GetMoveMask(nRemainder);
    __m512i   Aoffset    = GetRowOffset(lda);
    __m512i   BOffset    = GetRowOffset(ldb);
 
    int row = 0;
    for (; row < M - 1; row += 2) {
        //  Get one row of A
        GetTwoRowFromMatrix(Xrow1, Xrow2, A + row, lda, K, nRemainder, Aoffset, mask);
        F32PTR Crow = C + row;
        int    col = 0;
        for (; col < N; col++) {
            Crow[ldc*col+1] = __avx512_f32_dot2x1stride(Xrow1, Xrow2, B+col, K, mask,BOffset, ldb, Crow + ldc*col);;
        }

    }
    if (row < M) {
        //  Get one row of A
        GetOneRowFromMatrix(Xrow1, A + row, lda, K, nRemainder, Aoffset, mask);
        F32PTR Crow = C + row;
        int    col  = 0;
        for (; col < N; col++)
            Crow[ldc * col] = __avx512_f32_dot_stride(Xrow1, B+col, K,mask,BOffset,ldb);
    }
}


///////////////////////////////////////////////////////////////////////////

STATIC void avx512_f32_gemv_Xy1x1_slow(int N, int K, F32PTR X, int lda, F32PTR y, F32PTR C)
{     
    I32       nRemainder = K % NF;
    __mmask16 mask       = GetMoveMask(nRemainder);
    __m512i   Xoffset    = GetRowOffset(lda);    
    for (int row = 0; row < N ; ++row){
            C[row] =__avx512_f32_dot_stride(y, X+row, K, mask, Xoffset, lda);
    }
}


static INLINE void fma_f32_axpy_inplace( const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m512  C = set1(a);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4)
        store(y + i,       _mm512_fmadd_ps(load(x + i),       C, load(y + i))),
        store(y + i + NF,  _mm512_fmadd_ps(load(x + i + NF),  C, load(y + i + NF))),
        store(y + i + NF2, _mm512_fmadd_ps(load(x + i + NF2), C, load(y + i + NF2))),
        store(y + i + NF3, _mm512_fmadd_ps(load(x + i + NF3), C, load(y + i + NF3)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, _mm512_fmadd_ps(load(x + i), C, load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), _mm512_fmadd_ps(load(x + i), C, load(y + i)) );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

static  void avx512_f32_axpy_inplace(const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m512  C = set1(a);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3)
        store(y + i, add(mul(load(x + i), C), load(y + i))),
        store(y + i + NF, add(mul(load(x + i + NF), C), load(y + i + NF))),
        store(y + i + NF2, add(mul(load(x + i + NF2), C), load(y + i + NF2)));
    //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, add(mul(load(x + i), C), load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), add(mul(load(x + i), C), load(y + i)));
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

static INLINE void __avx512_f32_axpy_inplace(const F32 a, const F32PTR x, F32PTR y, const int N) {
    //y=a*x+y;
    __m512  C = set1(a);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3)
        store(y + i,      add(mul(load(x + i), C),       load(y + i))),
        store(y + i+NF,   add(mul(load(x + i+ NF), C),   load(y + i+ NF))),
        store(y + i+NF2,  add(mul(load(x + i+ NF2), C),  load(y + i+ NF2)));
        //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF)
        store(y + i, add(mul(load(x + i), C), load(y + i)));
    int n = N - i;
    if (n > 0) //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd
        maskstore(y + i, GetMoveMask(n), add(mul(load(x + i), C), load(y + i))  );
    _mm256_zeroupper();
    //for (; i < N; i++)         x[i] += c;
}

#define fmadd _mm512_fmadd_ps  
/*
TODO: for uknown reasons, fmadd in the following function does not give the correct results. It could be
a gcc-specific compiler bug.
*/

static INLINE void __fma_f32_ax1bx2py_inplace(const F32PTR a, const F32PTR x1, const F32PTR x2, F32PTR y, const int N) {
    //y=a1*x1+a2x2+y;
    __m512  A = set1(a[0]);
    __m512  B = set1(a[1]);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3) {
        __m512 Y1 = load(y + i);      Y1 = fmadd(load(x1 + i),     A, Y1); Y1 = fmadd(load(x2 + i),     B, Y1);   store(y + i,     Y1);
        __m512 Y2 = load(y + i+NF);   Y2 = fmadd(load(x1 + i+NF ), A, Y2); Y2 = fmadd(load(x2 + i+NF),  B, Y2);   store(y + i+NF,  Y2);
        __m512 Y3 = load(y + i+NF2);  Y3 = fmadd(load(x1 + i+NF2), A, Y3); Y3 = fmadd(load(x2 + i+NF2), B, Y3);   store(y + i+NF2, Y3);
    }
        //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF) {
        __m512 Y1 = load(y + i);      Y1 = fmadd(load(x1 + i), A, Y1); Y1 = fmadd(load(x2 + i), B, Y1);   store(y + i, Y1);
    }    
    int n = N - i;
    if (n > 0) { //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd 
        __m512 Y1 = load(y + i);      Y1 = fmadd(load(x1 + i), A, Y1); Y1 = fmadd(load(x2 + i), B, Y1);   store(y + i, Y1);
        maskstore(y+i, GetMoveMask(n), Y1);      
    }  
    _mm256_zeroupper(); 
}
static INLINE void __avx512_f32_ax1bx2py_inplace(const F32PTR a, const F32PTR x1, const F32PTR x2, F32PTR y, const int N) {
    //y=a1*x1+a2x2+y;
    __m512  C1 = set1(a[0]);
    __m512  C2 = set1(a[1]);
    int i = 0;
    for (; i < N - (NF3 - 1); i += NF3) {
        __m512 t1 = add(mul(load(x1 + i), C1),      mul(load(x2 + i), C2));          store(y + i,     add(t1, load(y + i)));
        __m512 t2 = add(mul(load(x1 + i +NF), C1),  mul(load(x2 + i + NF), C2));     store(y + i +NF, add(t2, load(y + i + NF)));
        __m512 t3 = add(mul(load(x1 + i+ NF2), C1),   mul(load(x2 + i+ NF2), C2));   store(y + i+ NF2,  add(t3, load(y + i+ NF2)));
    }
        //store(y + i+24, add(mul(load(x + i+24), C), load(y + i+24)));
    for (; i < N - (NF - 1); i += NF) {
        __m512 t1 = add(mul(load(x1 + i), C1), mul(load(x2 + i), C2)); store(y + i, add(t1, load(y + i)));
    }    
    int n = N - i;
    if (n > 0) { //TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd 
        __m512 t1 = add(mul(load(x1 + i), C1), mul(load(x2 + i), C2)); 
        maskstore(y + i, GetMoveMask(n), add(t1, load(y + i)));      
    }  
    _mm256_zeroupper(); 
}
STATIC void avx512_f32_gemv_Xb(int N, int K, F32PTR X, int lda, F32PTR b, F32PTR C)
{
    //I32     nRemainder = K % 8;
    //__mmask16 mask    = GetMoveMask(nRemainder);
    //__m512i Xoffset = GetRowOffset(lda);
    memset(C, 0, sizeof(F32) * N);
    int row = 0;
    for (; row < N - (256 * 2 -1); row += 256*2) {
        int col = 0;
        for (; col < K-1; col+=2) {
            __avx512_f32_ax1bx2py_inplace(b+col, X + col * lda + row, X + (col+1) * lda + row, C + row, 256 * 2);
            //fma_f32_ax1bx2py_inplace(b + col, X + col * lda + row, X + (col + 1) * lda + row, C + row, 256 * 2);
            //avx512_f32_axpy_inplace_v2(b[col], X + col * lda+row, C+row, 256);
        }
        if(col<K)
            //fma_f32_axpy_inplace(b[col], X + col * lda + row, C + row, 256 * 2);
            __avx512_f32_axpy_inplace(b[col], X + col * lda + row, C + row, 256 * 2);
    }
    int n = N - row;
    if (n > 0) {
        int col = 0;
        for (; col < K-1; col+=2) {
            __avx512_f32_ax1bx2py_inplace(b+col, X + col * lda + row, X + (col+1) * lda + row, C + row, n);
        }
        if(col<K)
            __avx512_f32_axpy_inplace(b[col], X + col * lda + row, C + row, n);
    }
    

}
 


////////////////////////////////////////////////////////////////
static INLINE I32 _avx512_f32_findindex_cmp_lt(F32PTR  x, I32PTR indices, F32 value, int N ) {
    // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m512  Values = set1(value);
    int     cnt    = 0;

    int     i      = 0;
    for (; i < N - (NF - 1); i += NF) {
        __mmask16 mask = _mm512_cmp_ps_mask( load(x + i), Values, _CMP_LT_OQ);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed up: if (keep)            
            int keep     =  mask & 0x1;
            indices[cnt] =  segIdx++;     //if keep=false, the assigned sgeIdx value is invaliad and will be upldated later (cnt is unchanged)        
            cnt          = cnt  + keep;
            mask         = mask >> 1L;
        }
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 movemask = GetMoveMask(n);
        __mmask16  mask    = _mm512_mask_cmp_ps_mask(movemask, load(x + i), Values, _CMP_LT_OQ);
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep     = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt  = cnt  + keep;
            mask = mask >> 1L;
        }
    }
    _mm256_zeroupper();

    return cnt;
}
static INLINE I32 _avx512_f32_findindex_cmp_le(F32PTR  x, I32PTR indices, F32 value, int N ) {
    // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m512  Values = set1(value);
    int     cnt    = 0;

    int     i      = 0;
    for (; i < N - (NF - 1); i += NF) {
        __mmask16 mask = _mm512_cmp_ps_mask( load(x + i), Values, _CMP_LE_OQ);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed up: if (keep)            
            int keep     =  mask & 0x1;
            indices[cnt] =  segIdx++;     //if keep=false, the assigned sgeIdx value is invaliad and will be upldated later (cnt is unchanged)        
            cnt          = cnt  + keep;
            mask         = mask >> 1L;
        }
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 movemask = GetMoveMask(n);
        __mmask16  mask    = _mm512_mask_cmp_ps_mask(movemask, load(x + i), Values, _CMP_LE_OQ);
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep     = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt  = cnt  + keep;
            mask = mask >> 1L;
        }
    }
    _mm256_zeroupper();

    return cnt;
}
static INLINE I32 _avx512_f32_findindex_cmp_gt(F32PTR  x, I32PTR indices, F32 value, int N ) {
    // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m512  Values = set1(value);
    int     cnt    = 0;

    int     i      = 0;
    for (; i < N - (NF - 1); i += NF) {
        __mmask16 mask = _mm512_cmp_ps_mask( load(x + i), Values, _CMP_GT_OQ);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed up: if (keep)            
            int keep     =  mask & 0x1;
            indices[cnt] =  segIdx++;     //if keep=false, the assigned sgeIdx value is invaliad and will be upldated later (cnt is unchanged)        
            cnt          = cnt  + keep;
            mask         = mask >> 1L;
        }
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 movemask = GetMoveMask(n);
        __mmask16  mask    = _mm512_mask_cmp_ps_mask(movemask, load(x + i), Values, _CMP_GT_OQ);
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep     = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt  = cnt  + keep;
            mask = mask >> 1L;
        }
    }
    _mm256_zeroupper();

    return cnt;
}
static INLINE I32 _avx512_f32_findindex_cmp_ge(F32PTR  x, I32PTR indices, F32 value, int N ) {
    // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m512  Values = set1(value);
    int     cnt    = 0;

    int     i      = 0;
    for (; i < N - (NF - 1); i += NF) {
        __mmask16 mask = _mm512_cmp_ps_mask( load(x + i), Values, _CMP_GE_OQ);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed up: if (keep)            
            int keep     =  mask & 0x1;
            indices[cnt] =  segIdx++;     //if keep=false, the assigned sgeIdx value is invaliad and will be upldated later (cnt is unchanged)        
            cnt          = cnt  + keep;
            mask         = mask >> 1L;
        }
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 movemask = GetMoveMask(n);
        __mmask16  mask    = _mm512_mask_cmp_ps_mask(movemask, load(x + i), Values, _CMP_GE_OQ);
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep     = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt  = cnt  + keep;
            mask = mask >> 1L;
        }
    }
    _mm256_zeroupper();

    return cnt;
}
static INLINE I32 _avx512_f32_findindex_cmp_eq(F32PTR  x, I32PTR indices, F32 value, int N ) {
    // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    __m512  Values = set1(value);
    int     cnt    = 0;

    int     i      = 0;
    for (; i < N - (NF - 1); i += NF) {
        __mmask16 mask = _mm512_cmp_ps_mask( load(x + i), Values, _CMP_EQ_OQ);
        int  segIdx = i;
        while (mask) {
            //https://stackoverflow.com/questions/20927710/quickly-count-number-of-zero-valued-bytes-in-an-array/20933337#20933337
            //Remove the branches to speed up: if (keep)            
            int keep     =  mask & 0x1;
            indices[cnt] =  segIdx++;     //if keep=false, the assigned sgeIdx value is invaliad and will be upldated later (cnt is unchanged)        
            cnt          = cnt  + keep;
            mask         = mask >> 1L;
        }
    }

    int n = N - i;
    if (n > 0) {
        __mmask16 movemask = GetMoveMask(n);
        __mmask16  mask    = _mm512_mask_cmp_ps_mask(movemask, load(x + i), Values, _CMP_EQ_OQ);
        int  segIdx = i;
        while (mask) {
            //avoid branches: if (keep)
            int keep     = mask & 0x1L;
            indices[cnt] = segIdx++;
            cnt  = cnt  + keep;
            mask = mask >> 1L;
        }
    }
    _mm256_zeroupper();

    return cnt;
}
STATIC I32 avx512_f32_findindex(F32PTR  x, I32PTR indices, F32 value, int N, CmpFlag flag) {
    /*
    https://stackoverflow.com/questions/18971401/sparse-array-compression-using-simd-avx2
    https://stackoverflow.com/questions/36932240/avx2-what-is-the-most-efficient-way-to-pack-left-based-on-a-mask
    AVX512 has the compress insturctions that allow copying scattered elements into a contiguous memory accroding to a mask
    https://stackoverflow.com/questions/25074197/compact-avx2-register-so-selected-integers-are-contiguous-according-to-mask
    */
    I32 cnt;
    int cmpflag;
    //   // cmpflag: _CMP_LT_OQ _CMP_LE_OQ _CMP_GT_OQ _CMP_GE_OQ _CMP_LT_OQ
    switch (flag) {
    case CMP_LT:
        cnt=_avx512_f32_findindex_cmp_lt(x, indices, value, N);  break;
    case CMP_LE:
        cnt = _avx512_f32_findindex_cmp_le(x, indices, value, N);  break; //cmpflag = _CMP_LE_OQ; break;
    case CMP_GT:
        cnt = _avx512_f32_findindex_cmp_gt(x, indices, value, N);  break; //cmpflag = _CMP_GT_OQ; break;
    case CMP_GE:
        cnt = _avx512_f32_findindex_cmp_ge(x, indices, value, N);  break; //cmpflag = _CMP_GE_OQ; break;
    case CMP_EQ:
        cnt = _avx512_f32_findindex_cmp_eq(x, indices, value, N);  break; // cmpflag = _CMP_EQ_OQ; break;
    }
  
    return cnt;
}


STATIC void avx512_f32_scatter_val_byindex(F32PTR  x, I32PTR index, F32 value, int N) {
    
    __m512  vec = set1(value);
    int     i   = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512i indices1  = loadi(index +i);           _mm512_i32scatter_ps(x, indices1, vec , 4);        
        __m512i indices2  = loadi(index + i+NF);       _mm512_i32scatter_ps(x, indices2, vec, 4);
        __m512i indices3  = loadi(index + i+NF2);      _mm512_i32scatter_ps(x, indices3, vec, 4);
        __m512i indices4  = loadi(index + i + NF3);     _mm512_i32scatter_ps(x, indices4, vec, 4);
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512i indices1 = loadi(index + i);           _mm512_i32scatter_ps(x, indices1, vec, 4);        
    }
    int n = N - i;
    if (n >0) {
        __mmask16 loadmask = GetMoveMask(n);
        __m512i   indices1 = maskloadi(index+i, loadmask);
        _mm512_mask_i32scatter_ps(x, loadmask, indices1, vec, 4);  
    }
    _mm256_zeroupper();
  
    
}

STATIC void avx512_f32_scatter_vec_byindex(F32PTR  x, I32PTR index, F32PTR value, int N) {

    int     i   = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512i indices1 = loadi(index +i);          __m512 vec1  = load(value + i);         _mm512_i32scatter_ps(x, indices1, vec1, 4);        
        __m512i indices2 = loadi(index + i+NF);      __m512 vec2 = load(value + i + NF);     _mm512_i32scatter_ps(x, indices2, vec2, 4);
        __m512i indices3 = loadi(index + i+NF2);     __m512 vec3 = load(value + i + NF2);    _mm512_i32scatter_ps(x, indices3, vec3, 4);
        __m512i indices4 = loadi(index + i+NF3);     __m512 vec4 = load(value + i + NF3);    _mm512_i32scatter_ps(x, indices4, vec4, 4);
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512i indices1 = loadi(index + i);   __m512 vec1 = load(value + i);     _mm512_i32scatter_ps(x, indices1, vec1, 4);
    }
    int n = N - i;
    if (n >0) {
        __mmask16 loadmask = GetMoveMask(n);
        __m512i   indices1 = maskloadi(index+i, loadmask);
        __m512    vec1     = maskload(value+i, loadmask);
        _mm512_mask_i32scatter_ps(x, loadmask, indices1, vec1, 4);
  
    }
    _mm256_zeroupper();
  
    
} 

#define vf5 __m512
#define vi5 __m512i

STATIC void avx512_f32_gatherVec_scatterVal_byindex(F32PTR  x, I32PTR index, F32PTR value, F32 newValue, int N) {

    __m512  vec = set1(newValue);

    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        vi5  indices1  = loadi(index + i);       vf5  oldValue1 = _mm512_i32gather_ps(indices1, x, 4);
        store(value + i, oldValue1);             _mm512_i32scatter_ps(x, indices1, vec, 4);

        vi5  indices2 = loadi(index + i+NF);     vf5  oldValue2 = _mm512_i32gather_ps(indices2, x, 4);
        store(value +i+NF, oldValue2);           _mm512_i32scatter_ps(x, indices2, vec, 4);

        vi5  indices3 = loadi(index + i+NF2);     vf5  oldValue3 = _mm512_i32gather_ps(indices3, x, 4);
        store(value + i + NF2, oldValue3);        _mm512_i32scatter_ps(x, indices3, vec, 4);

        vi5  indices4 = loadi(index + i + NF3);   vf5  oldValue4 = _mm512_i32gather_ps(indices4, x, 4);
        store(value + i + NF3, oldValue4);        _mm512_i32scatter_ps(x, indices4, vec, 4);
    }    
    for (; i < N - (NF - 1); i += NF) {
        vi5  indices1 = loadi(index + i);       vf5  oldValue1 = _mm512_i32gather_ps(indices1, x, 4);
        store(value + i, oldValue1);            _mm512_i32scatter_ps(x, indices1, vec, 4); 
    }

    int n = N - i;
    if (n >0) {
        __mmask16 loadmask = GetMoveMask(n);
        vi5  indices1   = maskloadi(index + i,loadmask);   
        vf5  oldValue1 = _mm512_mask_i32gather_ps(set0(),loadmask, indices1, x, 4);
        maskstore(value + i, loadmask,oldValue1);            
        _mm512_mask_i32scatter_ps(x, loadmask, indices1, vec, 4);
    }
    _mm256_zeroupper();

} 
STATIC void avx512_f32_gather2Vec_scatterVal_byindex(F32PTR  x, F32PTR  y, I32PTR index, F32PTR value, F32 newValue, int N) {

     __m512  vec = set1(newValue);

    int i = 0;
    for (; i < N - (NF2 - 1); i += NF2) {
        vi5  indices1  = loadi(index + i);       vf5  oldX1 = _mm512_i32gather_ps(indices1, x, 4);   vf5  oldY1 = _mm512_i32gather_ps(indices1, y, 4);
        store(value +   i, oldX1);               _mm512_i32scatter_ps(x, indices1, vec, 4);
        store(value +N+ i, oldY1);               _mm512_i32scatter_ps(y, indices1, vec, 4);

        vi5  indices2  = loadi(index + i+NF);    vf5  oldX2 = _mm512_i32gather_ps(indices2, x, 4);   vf5  oldY2 = _mm512_i32gather_ps(indices2, y, 4);
        store(value +   i + NF, oldX1);          _mm512_i32scatter_ps(x, indices2, vec, 4);
        store(value +N+ i + NF, oldY1);          _mm512_i32scatter_ps(y, indices2, vec, 4);

    }    
    for (; i < N - (NF - 1); i += NF) {
        vi5  indices1  = loadi(index + i);       vf5  oldX1 = _mm512_i32gather_ps(indices1, x, 4);   vf5  oldY1 = _mm512_i32gather_ps(indices1, y, 4);
        store(value +   i, oldX1);               _mm512_i32scatter_ps(x, indices1, vec, 4);
        store(value +N+ i, oldY1);               _mm512_i32scatter_ps(y, indices1, vec, 4);
    }

    int n = N - i;
    if (n >0) {
        __mmask16 loadmask = GetMoveMask(n);
        vi5  indices1   = maskloadi(index + i,loadmask);   
        vf5  oldX1 = _mm512_mask_i32gather_ps(set0(),loadmask, indices1, x, 4);
        vf5  oldY1 = _mm512_mask_i32gather_ps(set0(), loadmask, indices1, x, 4);
        maskstore(value    +i, loadmask, oldX1);            
        maskstore(value +N +i, loadmask, oldY1);
        _mm512_mask_i32scatter_ps(x, loadmask, indices1, vec, 4);
        _mm512_mask_i32scatter_ps(y, loadmask, indices1, vec, 4);
    }
    _mm256_zeroupper();

} 
STATIC void avx512_f32_scale_inplace(const F32 gain, const F32 offset,const F32PTR x, const int N) {
    __m512  G = set1(gain);
    __m512  O = set1(offset);
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

STATIC void avx512_f32_hinge_pos(const F32PTR X, const F32PTR Y, const F32 knot, const int N){
    __m512  O = set0();
    __m512  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 d1 = sub(load(X + i),    C);   store(Y + i,       _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1)); //0: take a, 1: take from vec        
        __m512 d2 = sub(load(X + i+NF), C);   store(Y + i + NF,  _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d2, O, _CMP_GE_OQ), O, d2)); //0: take a, 1: take from vec        
        __m512 d3 = sub(load(X + i+NF2), C);  store(Y + i + NF2, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d3, O, _CMP_GE_OQ), O, d3)); //0: take a, 1: take from vec        
        __m512 d4 = sub(load(X + i+NF3), C);  store(Y + i + NF3, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d4, O, _CMP_GE_OQ), O, d4)); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 d1 = sub(load(X + i), C);   store(Y + i, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1)); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m512 d1 = sub(load(X + i), C);  
        maskstore(Y + i, GetMoveMask(n), _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1) ); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}
STATIC void avx512_f32_hinge_neg(const F32PTR X, const F32PTR Y, const F32 knot, const int N){
    __m512  O = set0();
    __m512  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 d1 = sub(C, load(X + i)    );   store(Y + i,       _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1)); //0: take a, 1: take from vec        
        __m512 d2 = sub(C, load(X + i+NF) );   store(Y + i + NF,  _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d2, O, _CMP_GE_OQ), O, d2)); //0: take a, 1: take from vec        
        __m512 d3 = sub(C, load(X + i+NF2));   store(Y + i + NF2, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d3, O, _CMP_GE_OQ), O, d3)); //0: take a, 1: take from vec        
        __m512 d4 = sub(C, load(X + i+NF3));   store(Y + i + NF3, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d4, O, _CMP_GE_OQ), O, d4)); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 d1 = sub(C, load(X + i) );      store(Y + i,       _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1)); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m512 d1 = sub(C, load(X + i));  
        maskstore(Y + i, GetMoveMask(n), _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, d1) ); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}

STATIC void avx512_f32_step_pos(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N){
    __m512  O = set0();
    __m512  I = set1(1.0);
    __m512  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 d1 = sub(load(X + i),    C);   store(Z + i,       _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, load(Y + i ))); //0: take a, 1: take from vec        
        __m512 d2 = sub(load(X + i+NF), C);   store(Z + i + NF,  _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d2, O, _CMP_GE_OQ), O, load(Y + i + NF))); //0: take a, 1: take from vec        
        __m512 d3 = sub(load(X + i+NF2), C);  store(Z + i + NF2, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d3, O, _CMP_GE_OQ), O, load(Y + i + NF2))); //0: take a, 1: take from vec        
        __m512 d4 = sub(load(X + i+NF3), C);  store(Z + i + NF3, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d4, O, _CMP_GE_OQ), O, load(Y + i + NF3))); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 d1 = sub(load(X + i), C);      store(Z + i,      _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ),  O, load(Y+ i))); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m512 d1 = sub(load(X + i), C);  
        maskstore(Z + i, GetMoveMask(n), _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), O, load(Y + i  )) ); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}
STATIC void avx512_f32_step_neg(const F32PTR X, const F32PTR Y, const F32PTR Z, const F32 knot, const int N){
    __m512  O = set0();
    __m512  I = set1(1.0);
    __m512  C = set1(knot);
    int i = 0;
    for (; i < N - (NF4 - 1); i += NF4) {
        __m512 d1 = sub(load(X + i),    C);   store(Z + i,       _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), load(Y + i), O)); //0: take a, 1: take from vec        
        __m512 d2 = sub(load(X + i+NF), C);   store(Z + i + NF,  _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d2, O, _CMP_GE_OQ), load(Y + i+NF), O)); //0: take a, 1: take from vec        
        __m512 d3 = sub(load(X + i+NF2), C);  store(Z + i + NF2, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d3, O, _CMP_GE_OQ), load(Y + i+NF2), O)); //0: take a, 1: take from vec        
        __m512 d4 = sub(load(X + i+NF3), C);  store(Z + i + NF3, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d4, O, _CMP_GE_OQ), load(Y + i+NF3), O)); //0: take a, 1: take from vec        
    }
    for (; i < N - (NF - 1); i += NF) {
        __m512 d1 = sub(load(X + i), C);   store(Z + i, _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), load(Y + i), O)); //0: take a, 1: take from vec        
    }
    
    int n = N - i;
    if (n > 0) {//TODO:BUGGY, loadi should also been a masked load. If not, it may fail during the page bnd        
        __m512 d1 = sub(load(X + i), C);  
        maskstore(Z + i, GetMoveMask(n), _mm512_mask_blend_ps(_mm512_cmp_ps_mask(d1, O, _CMP_GE_OQ), load(Y + i), O) ); //0: take a, 1: take from vec        
    }
    _mm256_zeroupper();    
}

void SetupVectorFunction_AVX512() {

    FillMaskTemplate();
    i32_add_val_inplace = &avx512_i32_add_val_inplace;;
  
    i32_sum   = &avx512_i32_sum;                  
    f32_fill_val = &avx512_f32_fill_val;
    f32_sum  = &avx512_f32_sum;
    f32_add_vec   = &avx512_f32_add_vec;
    f32_sub_vec    = &avx512_f32_sub_vec;  
    f32_add_vec_inplace    = &avx512_f32_add_vec_inplace;
    f32_sub_vec_inplace    = &avx512_f32_sub_vec_inplace;
    f32_subrev_val_inplace    = &avx512_f32_subrev_val_inplace;

 
    f32_add_val_inplace    = &avx512_f32_add_val_inplace;
    f32_mul_val_inplace    = &avx512_f32_mul_val_inplace;
    f32_mul_vec_inplace    = &avx512_f32_mul_vec_inplace;
    f32_mul_vec            = &avx512_f32_mul_vec;
    f32_dot                = &avx512_f32_dot;

  
    f32_dot2x1    = &avx512_f32_dot2x1;
    f32_dot2x2    = &avx512_f32_dot2x2;
    f32_add_v_v2_vec_inplace    = &avx512_f32_add_v_v2_vec_inplace;
    f32_cos_vec_inplace    = &avx512_f32_cos_vec_inplace;
    f32_sin_vec_inplace    = &avx512_f32_sin_vec_inplace;
    f32_sincos_vec_inplace   = & avx512_f32_sincos_vec_inplace;
    f32_pow_vec_inplace   = & avx512_f32_pow_vec_inplace;
    f32_log_vec_inplace    = &avx512_f32_log_vec_inplace;
    f32_exp_vec_inplace     = &avx512_f32_exp_vec_inplace;
    f32_sqrt_vec_inplace    = &avx512_f32_sqrt_vec_inplace;
    f32_sqrt_vec    = &avx512_f32_sqrt_vec;
 
    //f32_sumlog_slow    = &avx512_f32_sumlog_slow;
    //f32_sumlog    = &avx512_f32_sumlog;
    f32_avgstd                      = &avx512_f32_avgstd;
    f32_sx_sxx_to_avgstd_inplace   = & avx512_f32_sx_sxx_to_avgstd_inplace;
    
    //f32_maxidx_slow     = &avx512_f32_maxidx_slow;    
    
 
    f32_maxidx     = &avx512_f32_maxidx;
    f32_minidx     = &avx512_f32_minidx;    
    f32_diff_back  = & avx512_f32_diff_back;
  
    f32_seq   = & avx512_f32_seq;
    i32_seq    =&avx512_i32_seq;

    f32_to_f64_inplace   = &avx512_f32_to_f64_inplace;
    f64_to_f32_inplace   = &avx512_f64_to_f32_inplace;
    
    i32_to_f32_scaleby_inplace        = &avx512_i32_to_f32_scaleby_inplace;
    i32_increment_bycond_inplace      = &avx512_i32_increment_bycond_inplace; 
    i32_increment_vec2_bycond_inplace = &avx512_i32_increment_vec2_bycond_inplace;

    //i08_find_nth_onebyte_binvec     = & avx512_i08_find_nth_onebyte_binvec;
    //i08_find_nth_onebyte_binvec_v2  = &  avx512_i08_find_nth_onebyte_binvec_v2;
    i08_sum_binvec   = &avx512_i08_sum_binvec;

    f32_gemm_XtY2x1 = avx512_f32_gemm_XtY2x1;
    f32_gemm_XtY2x2 = avx512_f32_gemm_XtY2x2;
  
    f32_gemm_XY1x2  = avx512_f32_gemm_XY1x2;
    f32_gemm_XY2x2   = avx512_f32_gemm_XY2x2;
    f32_gemm_XtYt2x2 = avx512_f32_gemm_XtYt2x2;
    f32_gemm_XYt2x1 = avx512_f32_gemm_XYt2x1;
    
    f32_gemv_Xb = &avx512_f32_gemv_Xb;
 
  
    f32_findindex = avx512_f32_findindex;
    f32_scatter_vec_byindex = avx512_f32_scatter_vec_byindex;
    f32_gatherVec_scatterVal_byindex = avx512_f32_gatherVec_scatterVal_byindex;
    f32_gather2Vec_scatterVal_byindex = &avx512_f32_gather2Vec_scatterVal_byindex;

    f32_scale_inplace = avx512_f32_scale_inplace;
    f32_hinge_pos = avx512_f32_hinge_pos;
    f32_hinge_neg = avx512_f32_hinge_neg;
    f32_step_pos = avx512_f32_step_pos;
    f32_step_neg = avx512_f32_step_neg;
    f32_axpy_inplace = avx512_f32_axpy_inplace;
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