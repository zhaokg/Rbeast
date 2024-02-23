#pragma once
#include <immintrin.h>
#include "abc_000_macro.h"
#include "abc_datatype.h"

typedef __m256  v8sf; // vector of 8 float (avx)
typedef __m256i v8si; // vector of 8 int   (avx)
typedef __m128i v4si; // vector of 8 int   (avx) 
typedef __m128  v4sf; // vector of 8 int   (avx) 

#if defined(COMPILER_MSVC)
extern v8sf log256_ps(F32PTR x);
extern v8sf exp256_ps(F32PTR x);
extern v8sf sin256_ps(F32PTR x);
extern v8sf cos256_ps(F32PTR x);
extern void sincos256_ps(F32PTR x, F32PTR s, F32PTR c);
extern v8sf pow256_ps(F32PTR x, float n);
#else
extern void log256_ps_ptr(F32PTR x, F32PTR out);
extern void exp256_ps_ptr(F32PTR x, F32PTR out);
extern void sin256_ps_ptr(F32PTR x, F32PTR out);
extern void cos256_ps_ptr(F32PTR x, F32PTR out);
extern void sincos256_ps_ptr(F32PTR x, F32PTR s, F32PTR c);
extern void pow256_ps_ptr(F32PTR x, float n, F32PTR out);


//https://stackoverflow.com/questions/6440021/compiler-support-of-gnu-statement-expression
#define log256_ps( x)                ({ float y[8]; log256_ps_ptr(x, y); _mm256_loadu_ps(y);})
#define exp256_ps( x)                ({ float y[8]; exp256_ps_ptr(x, y); _mm256_loadu_ps(y);})
#define sin256_ps(x)				 ({ float y[8]; sin256_ps_ptr(x, y); _mm256_loadu_ps(y);})
#define cos256_ps( x)				 ({ float y[8]; cos256_ps_ptr(x, y); _mm256_loadu_ps(y);})
#define sincos256_ps( x,  s,  c)	 ({ float y[8]; sincos256_ps_ptr(x,s,c); _mm256_loadu_ps(y);})
#define pow256_ps( x,  n)			 ({ float y[8]; pow256_ps_ptr(x,n, y); _mm256_loadu_ps(y);})
#endif

