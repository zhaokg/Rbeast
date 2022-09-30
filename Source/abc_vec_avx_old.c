#include "abc_datatype.h"
#include <stdio.h>
#include <immintrin.h>
#include <string.h>

#ifdef TARGET_64
//https://stackoverflow.com/questions/2622017/suppressing-deprecated-warnings-in-xcode
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
 //do something////
#pragma clang diagnostic pop

////////////////////////////////////////
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//do something////
#pragma GCC diagnostic pop
////////////////////////////////////////

#if defined (GCC_COMPILER)
    ////https://www.geeksforgeeks.org/speed-up-naive-algorithms-in-competitive-coding-in-c-cpp/
    ////https://codeforces.com/blog/entry/78897
    //#pragma GCC optimize("Ofast")
    //#pragma GCC target("avx,avx2")
    //https://stackoverflow.com/questions/61759552/why-some-top-level-competitive-programmers-use-pragma
    #pragma optimization_level 3
    //#pragma GCC optimize("Ofast,no-stack-protector,unroll-loops,fast-math,O3")
    //#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
    #pragma GCC optimize("Ofast")  //Comment optimisations for interactive problems (use endl)
    #pragma GCC target("avx,avx2,fma")
   // #pragma GCC optimization ("unroll-loops")
    //#pragma GCC target "arch=core-avx2,tune=core-avx2"
#endif

static INLINE __m256i GetMoveMask(int n) {
    __m128i maskIdx = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i maskNum = _mm256_cvtepu8_epi32(maskIdx);
    __m256i nvec    = _mm256_set1_epi32(n);
    __m256i maskmov = _mm256_sub_epi32(maskNum, nvec);
    return maskmov;
}

void avxu_mulC_inplace(F32PTR x, F32 c, int N) {

    __m256  C  = _mm256_set1_ps(c);    

    int i = 0;
    for (; i < N && (((intptr_t)&x[i] & 31ULL) != 0); ++i)
        x[i] *= c;

    for (; i < N - 31; i += 32) {
        _mm256_store_ps(x+i,    _mm256_mul_ps(*((__m256*) (x+ i)),    C));
        _mm256_store_ps(x+i+8,  _mm256_mul_ps(*((__m256*) (x+ i+8)),  C));
        _mm256_store_ps(x+i+16, _mm256_mul_ps(*((__m256*) (x+ i+16)), C));
        _mm256_store_ps(x+i+24, _mm256_mul_ps(*((__m256*) (x+ i+24)), C));
    }

    for (; i<N-7; i+=8)     
        _mm256_store_ps(x+i, _mm256_mul_ps( *((__m256*) (x + i)), C)  );    

    _mm256_zeroupper();

    for (; i < N; i++)  x[i] *= c;
}
void avxa_mulC_inplace(float* x, float c, int N) {

     __m256  C  = _mm256_set1_ps(c);    

    int i = 0;
    for (; i < N - 31; i += 32) {
        _mm256_store_ps(x+i,    _mm256_mul_ps(*((__m256*) (x+ i)), C));
        _mm256_store_ps(x+i+8,  _mm256_mul_ps(*((__m256*) (x+ i+8)), C));
        _mm256_store_ps(x+i+16, _mm256_mul_ps(*((__m256*) (x+ i+16)), C));
        _mm256_store_ps(x+i+24, _mm256_mul_ps(*((__m256*) (x+ i+24)), C));
    }

    for (; i<N-7; i+=8)     
        _mm256_store_ps(x+i, _mm256_mul_ps( *((__m256*) (x + i)), C)  );    


    _mm256_zeroupper();

    for (; i < N; i++)           x[i] *= c;
    
}

void avxu_addCi32_inplace(int* x, int c, int N) {

    __m256i  C  = _mm256_set1_epi32(c);    

    int i = 0;
    for (; i < N && (((intptr_t)&x[i] & 31ULL) != 0); ++i)
        x[i] += c;

    for (; i < N - 31; i += 32) {
        _mm256_store_si256(x+i,    _mm256_add_epi32(*((__m256i*) (x+ i)), C));
        _mm256_store_si256(x+i+8,  _mm256_add_epi32(*((__m256i*) (x+ i+8)), C));
        _mm256_store_si256(x+i+16, _mm256_add_epi32(*((__m256i*) (x+ i+16)), C));
        _mm256_store_si256(x+i+24, _mm256_add_epi32(*((__m256i*) (x+ i+24)), C));
    }

    for (; i<N-7; i+=8)     
        _mm256_store_si256(x+i, _mm256_add_epi32( *((__m256i*) (x + i)), C)  );


    _mm256_zeroupper();

    for (; i < N; i++)   
        x[i] += c;
    
}
void avxa_addCi32_inplace(int* x, int c, int N) {

    __m256i  C  = _mm256_set1_epi32(c);    

    int i = 0;

    for (; i < N - 31; i += 32) {
        _mm256_store_si256(x+i,    _mm256_add_epi32(*((__m256i*) (x+ i)), C));
        _mm256_store_si256(x+i+8,  _mm256_add_epi32(*((__m256i*) (x+ i+8)), C));
        _mm256_store_si256(x+i+16, _mm256_add_epi32(*((__m256i*) (x+ i+16)), C));
        _mm256_store_si256(x+i+24, _mm256_add_epi32(*((__m256i*) (x+ i+24)), C));
    }

    for (; i<N-7; i+=8)     
        _mm256_store_si256(x+i, _mm256_add_epi32( *((__m256i*) (x + i)), C)  );


    _mm256_zeroupper();

    for (; i < N; i++)   
        x[i] += c;
}

void  avxu_sfill( float s, float* x, int N) {
    
    __m256  S = _mm256_set1_ps(s);

    int i = 0;

    for (; i < N && (((intptr_t)&x[i] & 31ULL) != 0); ++i)
        x[i] = s;

    for (; i < N - 31; i += 32) {
        _mm256_store_ps((x + i),    S);
        _mm256_store_ps((x + i+8),  S);
        _mm256_store_ps((x + i+16), S);
        _mm256_store_ps((x + i+24), S);
    }

    for (; i < N - 7; i += 8)
        _mm256_store_ps((x+i), S);

    int n = N - i;
    if (n > 0) {
        _mm256_maskstore_ps(x + i, GetMoveMask(n), S);
    }    
    _mm256_zeroupper();
    //for (; i < N; i++)   x[i] = s;

}

void  avxa_sfill(float s, float* x, int N) {
    
 
    __m256  S = _mm256_set1_ps(s);

    int i = 0;

    for (; i < N - 31; i += 32) {
        _mm256_store_ps((x + i),    S);
        _mm256_store_ps((x + i+8),  S);
        _mm256_store_ps((x + i+16), S);
        _mm256_store_ps((x + i+24), S);
    }

    for (; i < N - 7; i += 8)
        _mm256_store_ps((x+i), S);

   

    int n = N - i;
    if (n > 0) {
        _mm256_maskstore_ps(x + i, GetMoveMask(n), S);
    }
    _mm256_zeroupper();
    // for (; i < N; i++)       x[i] = s;

}

/*
void  sadd0(float* x, float s, int N) {

    __m256  S = _mm256_set1_ps(s);

    int i = 0;
    for (; i < N - 7; i += 8)
        _mm256_store_ps(x, _mm256_add_ps(*((__m256*) (x + i)), S));

    for (; i < N; i++)
        x[i] += s;

}
void  sadd(float* x, float s, int N) {
    __m256  S = _mm256_set1_ps(s);

    int i = 0;
    for (; i < N - 31; i += 32) {
        _mm256_store_ps(x+i,    _mm256_add_ps(*((__m256*) (x + i)),   S));
        _mm256_store_ps(x+i+8,  _mm256_add_ps(*((__m256*) (x + i+8)), S));
        _mm256_store_ps(x+i+16, _mm256_add_ps(*((__m256*) (x+i+16)),  S));
        _mm256_store_ps(x+i+24, _mm256_add_ps(*((__m256*) (x+i+24)),  S));
    }

    for (; i < N - 7; i += 8) {
        _mm256_store_ps(x + i, _mm256_add_ps(*((__m256*) (x + i)), S));        
    }
 
    for (; i < N; i++) {
        x[i] *= s;
    }
    

    _mm256_zeroupper();
}
float saddu(float* x, float s, int N) {

    __m256  S = _mm256_set1_ps(s);

    int i = 0;
    for (; i < N && (((intptr_t)&x[i] & 31ULL) != 0); ++i)
        x[i] *= s;

    for (; i < N - 7; i += 8)
        _mm256_store_ps(x, _mm256_add_ps(*((__m256*) (x + i)), S));

    for (; i < N; i++)
        x[i] *= s;

}
*/

void  avxa_addv_inplace(float* __restrict src, float * __restrict dst, int N) {

    int i = 0;
    for (; i < N - 31; i += 32) {
        _mm256_store_ps(dst +i,    _mm256_add_ps(   *(__m256*) (dst+i),    *(__m256*) (src+i)     )   );
        _mm256_store_ps(dst +i+8,  _mm256_add_ps(   *(__m256*) (dst+i+8),  *(__m256*) (src+i+ 8)  )   );
        _mm256_store_ps(dst +i+16, _mm256_add_ps(   *(__m256*) (dst+i+16), *(__m256*) (src+i+ 16) )   );
        _mm256_store_ps(dst +i+24, _mm256_add_ps(   *(__m256*) (dst+i+24), *(__m256*) (src+i+ 24) )   );
    }

    for (; i < N - 7; i += 8) {
        _mm256_store_ps(dst +i,    _mm256_add_ps(*(__m256*) (dst+i),    *(__m256*) (src+i)     )   );
    }
 
    _mm256_zeroupper();

    for (; i < N; i++) {
        dst[i] += src[i];
    }
      
}

void  avxu_addv_inplace(float* __restrict src, float* __restrict dst, int N) {

    int i = 0;
    for (; i < N && (((intptr_t)&dst[i] & 31ULL) != 0); ++i)
        dst[i] += src[i];

    for (; i < N - 31; i += 32) {
        _mm256_store_ps(dst + i, _mm256_add_ps(*(__m256*) (dst + i), *(__m256*) (src + i)));
        _mm256_store_ps(dst + i + 8, _mm256_add_ps(*(__m256*) (dst + i + 8), *(__m256*) (src + i + 8)));
        _mm256_store_ps(dst + i + 16, _mm256_add_ps(*(__m256*) (dst + i + 16), *(__m256*) (src + i + 16)));
        _mm256_store_ps(dst + i + 24, _mm256_add_ps(*(__m256*) (dst + i + 24), *(__m256*) (src + i + 24)));
    }

    for (; i < N - 7; i += 8) {
        _mm256_store_ps(dst + i, _mm256_add_ps(*(__m256*) (dst + i), *(__m256*) (src + i)));
    }

    _mm256_zeroupper();

    for (; i < N; i++) {
        dst[i] += src[i];
    }

}




//https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction

F32 avxa_dot(F32PTR  x, F32PTR y, int N) {

    __m256 r = _mm256_setzero_ps();
    
    int i = 0;    
    for (; i < N - 31; i += 32) {
        __m256 r0 = _mm256_mul_ps(*(__m256*)(x + i),        *(__m256*)(y + i)) ;
        __m256 r1 = _mm256_mul_ps(*(__m256*)(x + i+8),      *(__m256*)(y + i+8));
        __m256 r2 = _mm256_mul_ps(*(__m256*)(x + i+16),     *(__m256*)(y + i +16));
        __m256 r3 = _mm256_mul_ps(*(__m256*)(x + i +24),    *(__m256*)(y + i + 24));

        __m256 s1   = _mm256_add_ps(r0, r1);
        __m256 s2   = _mm256_add_ps(r2, r3);
        __m256 s12  = _mm256_add_ps(s1, s2);
        r=_mm256_add_ps(r, s12);
    }        
    for (; i < N - 7; i += 8) {
         r =_mm256_add_ps(r,  _mm256_mul_ps(*(__m256*)(x + i), *(__m256*)(y + i))  );        
    }

    float sum;
    ///horizontal summing
    __m128 vlow  = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128  = _mm_add_ps(vlow, vhigh);     // add the low 128

    //float hsum_ps_sse3(__m128 v) {
     __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
     __m128 sums = _mm_add_ps(v128, shuf);
     shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
     sums = _mm_add_ss(sums, shuf);
     sum  = _mm_cvtss_f32(sums);
    //}
                      
     _mm256_zeroupper();
    for (; i < N; i++)     sum+=x[i] *y[i];
    return sum;
}
F32 avxa_dot1(F32PTR  x, F32PTR y, int N) {

    __m256 r = _mm256_setzero_ps();
    
    int i = 0;    
    for (; i < N - 31; i += 32) {
        __m256 r0 = _mm256_mul_ps(*(__m256*)(x + i),        *(__m256*)(y + i)) ;
        __m256 r1 = _mm256_mul_ps(*(__m256*)(x + i+8),      *(__m256*)(y + i+8));
        __m256 r2 = _mm256_mul_ps(*(__m256*)(x + i+16),     *(__m256*)(y + i+16));
        __m256 r3 = _mm256_mul_ps(*(__m256*)(x + i+24),    *(__m256*)(y + i+24));

        __m256 s1   = _mm256_add_ps(r0, r1);
        __m256 s2   = _mm256_add_ps(r2, r3);
        __m256 s12  = _mm256_add_ps(s1, s2);
        r=_mm256_add_ps(r, s12);
    }        
    for (; i < N - 7; i += 8) {
         r =_mm256_add_ps(r,  _mm256_mul_ps(*(__m256*)(x + i), *(__m256*)(y + i))  );        
    }

    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        __m256  r0   =  _mm256_mul_ps( _mm256_maskload_ps(x + i, mask), _mm256_maskload_ps(y + i, mask));
        r = _mm256_add_ps(r, r0);
    }

    float sum;
    ///horizontal summing
    __m128 vlow  = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128  = _mm_add_ps(vlow, vhigh);     // add the low 128

    //float hsum_ps_sse3(__m128 v) {
     __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
     __m128 sums = _mm_add_ps(v128, shuf);
     shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
     sums = _mm_add_ss(sums, shuf);
     sum  = _mm_cvtss_f32(sums);
    //}
                      
     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return sum;
}
F32 avxa_dot2x1(F32PTR  x, F32PTR y, F32PTR v, int N, F32PTR res) {

    __m256 RX = _mm256_setzero_ps();
    __m256 RY = _mm256_setzero_ps();

    int i = 0;    
    for (; i < N - 15; i += 16) {
        __m256 r0 = _mm256_mul_ps(*(__m256*)(x + i),        *(__m256*)(v + i)) ;
        __m256 r1 = _mm256_mul_ps(*(__m256*)(x + i+8),      *(__m256*)(v + i+8));
        
        __m256 y0 = _mm256_mul_ps(*(__m256*)(y + i),        *(__m256*)(v + i)) ;
        __m256 y1 = _mm256_mul_ps(*(__m256*)(y + i+8),      *(__m256*)(v + i+8));
        
        __m256 s1   = _mm256_add_ps(r0, r1);                
        RX =_mm256_add_ps(RX, s1);

        __m256 ss1 = _mm256_add_ps(y0, y1);        
         RY   = _mm256_add_ps(RY, ss1);
    }        
    for (; i < N - 7; i += 8) {
         RX =_mm256_add_ps(RX,  _mm256_mul_ps(*(__m256*)(x + i), *(__m256*)(v + i))  );        
         RY = _mm256_add_ps(RY, _mm256_mul_ps(*(__m256*)(y + i), *(__m256*)(v + i)));
    }

    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        __m256  vvec = _mm256_maskload_ps(v + i, mask);

        __m256  r0   =  _mm256_mul_ps( _mm256_maskload_ps(x + i, mask), vvec);
        RX = _mm256_add_ps(RX, r0);

        __m256  r1 = _mm256_mul_ps(_mm256_maskload_ps(y + i, mask), vvec);
        RY = _mm256_add_ps(RY, r1);
    }

    float sumX;
    {
        ///horizontal summing
        __m128 vlow = _mm256_castps256_ps128(RX);
        __m128 vhigh = _mm256_extractf128_ps(RX, 1); // high 128
        __m128 v128 = _mm_add_ps(vlow, vhigh);     // add the low 128

        //float hsum_ps_sse3(__m128 v) {
        __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
        __m128 sums = _mm_add_ps(v128, shuf);
        shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
        sums = _mm_add_ss(sums, shuf);
        sumX = _mm_cvtss_f32(sums);
        //}
    }      
    res[0] = sumX;

    float sumY; 
    {
        ///horizontal summing
        __m128 vlow = _mm256_castps256_ps128(RY);
        __m128 vhigh = _mm256_extractf128_ps(RY, 1); // high 128
        __m128 v128 = _mm_add_ps(vlow, vhigh);     // add the low 128

        //float hsum_ps_sse3(__m128 v) {
        __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
        __m128 sums = _mm_add_ps(v128, shuf);
        shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
        sums = _mm_add_ss(sums, shuf);
        sumY = _mm_cvtss_f32(sums);
        //}
    }


     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return sumY;
}
F32 avxa_dot2(F32PTR  x, F32PTR y, int N) {

    __m256 r0 = _mm256_setzero_ps();
    __m256 r1 = _mm256_setzero_ps();
    __m256 r2 = _mm256_setzero_ps();
    __m256 r3 = _mm256_setzero_ps();
    
    int i = 0;    
    for (; i < N - 31; i += 32) {
        r0 = _mm256_fmadd_ps(*(__m256*)(x + i),        *(__m256*)(y + i),   r0) ;
        r1 = _mm256_fmadd_ps(*(__m256*)(x + i+8),      *(__m256*)(y + i+8), r1);
        r2 = _mm256_fmadd_ps(*(__m256*)(x + i+16),     *(__m256*)(y + i+16),r2);
        r3 = _mm256_fmadd_ps(*(__m256*)(x + i+24),    *(__m256*)(y + i+24), r3);
        
    }        
    for (; i < N - 7; i += 8) {
         r0 = _mm256_fmadd_ps( *(__m256*)(x + i), *(__m256*)(y + i),  r0  );
    }
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
        r1   = _mm256_fmadd_ps( _mm256_maskload_ps(x + i, mask), _mm256_maskload_ps(y + i, mask),r1);
    }

    __m256 s1 = _mm256_add_ps(r0, r1);
    __m256 s2 = _mm256_add_ps(r2, r3);
    __m256 r = _mm256_add_ps(s1, s2);

    float sum;
    ///horizontal summing
    __m128 vlow  = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128  = _mm_add_ps(vlow, vhigh);     // add the low 128

    //float hsum_ps_sse3(__m128 v) {
     __m128 shuf = _mm_movehdup_ps(v128);        // broadcast elements 3,1 to 2,0
     __m128 sums = _mm_add_ps(v128, shuf);
     shuf = _mm_movehl_ps(shuf, sums); // high half -> low half
     sums = _mm_add_ss(sums, shuf);
     sum  = _mm_cvtss_f32(sums);
    //}
                      
     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return sum;
}
F32 avxa_max(const F32PTR x, const int N, F32PTR val) {

    
    __m256 r1 = _mm256_setzero_ps();
    
    __m256  maxVec      = _mm256_loadu_ps(x);
    __m128i _maskIdx    = _mm_cvtsi64_si128(0x0706050403020100);
    __m256i maxIdx      = _mm256_cvtepu8_epi32(_maskIdx);
    __m256i eight       = _mm256_set1_epi32(8);

    int i = 0;    
    for (; i < N - 7; i += 8) {

        __m256 cmpmask = _mm256_cmp_ps(maxVec, *(__m256*)(x + i), _CMP_LE_OQ);
         maxVec  = _mm256_blendv_ps(maxVec, *(__m256*)(x + i), cmpmask);

        __m256i newIdx = _mm256_add_epi32(maxIdx, eight);
        maxIdx  = _mm256_blendv_epi8( maxIdx,  newIdx, _mm256_castps_si256(cmpmask));        
    }        
    for (; i < N - 7; i += 8) {
        // r0 = _mm256_fmadd_ps( *(__m256*)(x + i), *(__m256*)(y + i),  r0  );
    }
    int n = N - i;
    if (n > 0) {
        __m256i mask = GetMoveMask(n);
       // r1   = _mm256_fmadd_ps( _mm256_maskload_ps(x + i, mask), _mm256_maskload_ps(y + i, mask),r1);
    }

   
                      
     _mm256_zeroupper();
    //for (; i < N; i++)     sum+=x[i] *y[i];
    return 1;
}
F32 sdot1(F32PTR x, F32PTR y, int N) {


    float  sum;

    __m256 r = _mm256_setzero_ps();

    int i = 0;
    for (; i < N - 7; i += 8)
        r = _mm256_add_ps(r, _mm256_mul_ps(*(__m256*)(x + i), *(__m256*)(y + i)));

    ///horizontal summing
    __m128 vlow = _mm256_castps256_ps128(r);
    __m128 vhigh = _mm256_extractf128_ps(r, 1); // high 128
    __m128 v128 = _mm_add_ps(vlow, vhigh);     // add the low 128


    //float hsum_ps_sse3(__m128 v) {
    __m128 sums = _mm_hadd_ps(v128, v128);
    sums = _mm_hadd_ps(sums, sums);
    sum = _mm_cvtss_f32(sums);
    //}


    for (; i < N; i++)
        sum += x[i] * y[i];
    return sum;
}
 
#else
static char a = 'a';
#endif