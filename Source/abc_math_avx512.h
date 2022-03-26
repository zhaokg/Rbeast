#pragma once
#include <immintrin.h>
typedef __m512  v16sf; // vector of 8 float (avx)
typedef __m512i v16si; // vector of 8 int   (avx)
typedef __m256  v8sf; // vector of 8 float (avx)
typedef __m256i v8si; // vector of 8 int   (avx)
typedef __m128i v4si; // vector of 8 int   (avx)

extern v16sf log512_ps(v16sf x);
extern v16sf exp512_ps(v16sf x);
extern v16sf sin512_ps(v16sf x);
extern v16sf cos512_ps(v16sf x);
extern void  sincos512_ps(v16sf x, v16sf* s, v16sf* c);
extern v16sf pow512_ps(v16sf x, float n);

 