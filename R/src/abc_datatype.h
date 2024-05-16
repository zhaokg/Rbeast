#pragma once

//#include <stdint.h>
#include <inttypes.h> /* contains stdint.h*/
#include <float.h>
#include "abc_000_macro.h"

//#define FLOAT_TYPE 8  //FLOAT is double
#define FLOAT_TYPE 4    //FLOAT is F32

#ifdef FLOAT
#undef FLOAT
#endif

#if FLOAT_TYPE==4
	//#define FLOAT F32
	typedef float FLOAT;
	#define FLOAT_MAX  FLT_MAX
	#define FLOAT_MIN (-FLT_MAX)
	#define FLOAT_EPSILON FLT_EPSILON
#elif FLOAT_TYPE==8
	//#define FLOAT double
	typedef double FLOAT;
	#define FLOAT_MAX  DBL_MAX
	#define FLOAT_MIN (-DBL_MAX)
	#define FLOAT_EPSILON DBL_EPSILON
#endif


#if defined(IsNaN)
	#undef IsNaN	
#endif
#define IsNaN(x)  (  (x) != (x))



#if defined(IsInf)
    #undef IsInf	
#endif
#define IsInf(f)  (  ((float)f) >FLOAT_MAX || ((float)f)< FLOAT_MIN)

/*------------------------------------------------------------*/

#ifdef __cplusplus
	extern "C" {
#endif 

	typedef float     F32;
	typedef double    F64;
	typedef int64_t   I64;
	typedef uint64_t  U64;
	typedef int32_t   I32;
	typedef uint32_t  U32;
	typedef int16_t   I16;
	typedef uint16_t  U16;
	typedef int8_t    I08;
	typedef uint8_t   U08;

	#define rFLOAT		register FLOAT
	#define rF32		register float
	#define rF64		register double
	#define rI64		register int64_t
	#define rU64		register uint64_t
	#define rI32		register int32_t
	#define rU32		register uint32_t
	#define rI16		register int16_t
	#define rU16		register uint16_t
	#define rI08		register int8_t
	#define rU08		register uint8_t

	typedef FLOAT* _restrict FLOATPTR;
	typedef float* _restrict F32PTR;
	typedef double* _restrict F64PTR;
	typedef int64_t* _restrict I64PTR;
	typedef uint64_t* _restrict U64PTR;
	typedef int32_t* _restrict I32PTR;
	typedef uint32_t* _restrict U32PTR;
	typedef int16_t* _restrict I16PTR;
	typedef uint16_t* _restrict U16PTR;
	typedef int8_t* _restrict I08PTR;
	typedef uint8_t* _restrict U08PTR;

	typedef  I08PTR I08PTR08; //8-byte alignment
	typedef  U08PTR U08PTR08; //8-byte alignment

	typedef void* _restrict  VOID_PTR;
	typedef void* _restrict  VOIDPTR;

	#define rFLOATPTR	register  FLOATPTR
	#define rF32PTR		register  F32PTR	
	#define rF64PTR		register  F64PTR
	#define rI64PTR		register  I64PTR	
	#define rU64PTR		register  U64PTR
	#define rI32PTR		register  I32PTR	
	#define rU32PTR		register  U32PTR
	#define rI16PTR		register  I16PTR	
	#define rU16PTR		register  U16PTR
	#define rI08PTR		register  I08PTR	
	#define rU08PTR		register  U08PTR
	#define rVOIDPTR	register  VOIDPTR

	typedef enum   DATA_TYPE {
		DATA_FLOAT = 0,
		DATA_DOUBLE,
		DATA_INT32,
		DATA_INT16,
		DATA_INT64,
		DATA_CHAR,
		DATA_STRUCT,
		DATA_UNKNOWN
	} DATA_TYPE;



// In C++, true and false are keywords and you can't have variables with that name. 
// C doesn't have Bool but C99 has _Bool.
// In C23,, Bool, true, and false become keywords:: https://groups.google.com/g/comp.lang.c/c/uTlZY_qKQmE

enum { _False_ = 0, _True_ = 1 };


//static INLINE F32 getNaN(void) { F32 nan = 1e38f;	nan = nan * nan * 0.f;	return nan;}
#if defined (COMPILER_MSVC)
	static INLINE float  getNaN(void) { return (float)1e38f * (float)1e38f * (float)0.f; }
#else
	static INLINE float  getNaN(void) { return (float)(0.0 / 0.0); }
#endif


#include <string.h>   // memset and memcpy

#define  free0(p)     if(p){ free(p);p=NULL;}
#define  malloc0(n)   calloc(1L, n)



#ifdef __cplusplus
}
#endif

#define max(a,b)			(((a) > (b)) ? (a) : (b))
#define min(a,b)			(((a) < (b)) ? (a) : (b))
// The outermost parentheses are DEFINITELY needed. If not, subtile errors
// can creep. Here is a failing example from beastv2_func.c
// min(Kbase + b->ke[j], Klastcol) - (Kbase + b->ks[j]) + 1

#define max3(a,b,c)         max( max(a,b), c)
#define max4(a,b,c,d)       max( max3(a,b,c),d)
#define _IsAlmostInteger(x)  ( fabs(x-round(x)) <1e-3 )

#define RoundTo64(N)       ((N+63)/64*64)
#define RoundTo8(N)        ((N+7)/8*8)
 


#define PostiveMod(i, n)  (  i%(n) + (i<0)*(n)  )
 
// A/B is a fix division: the rseult goes toward zero
// #define FLOORdiv(a,b) (int)floor((float)(a)/(b))
#define FLOORdiv(y,N) ((y >= 0 ? y : y - (N-1)) /N) 
#define  UnknownStatus 99
	
	