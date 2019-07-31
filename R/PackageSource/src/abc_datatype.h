#pragma once
#include "abc_000_macro.h"
#include "float.h"
#define FLOAT_TYPE 4   
#ifdef FLOAT
#undef FLOAT
#endif
#if FLOAT_TYPE==4
typedef float FLOAT;
#define FLOAT_MAX  FLT_MAX
#define FLOAT_MIN (-FLT_MAX)
#define FLOAT_EPSILON FLT_EPSILON
#elif FLOAT_TYPE==8
typedef double FLOAT;
#define FLOAT_MAX  DBL_MAX
#define FLOAT_MIN (-DBL_MAX)
#define FLOAT_EPSILON DBL_EPSILON
#endif
#if !defined(PI)
	#define PI   3.141592653589793
#endif
#if !defined(ISNAN)
#define ISNAN(f) (f) !=(f)
#endif
#define ISINF(f) (f) >FLOAT_MAX||(f)< FLOAT_MIN
#include <inttypes.h>
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
typedef FLOAT     * _restrict FLOATPTR;
typedef float   *   _restrict F32PTR;
typedef double   *  _restrict F64PTR;
typedef int64_t *   _restrict I64PTR;
typedef uint64_t*   _restrict U64PTR;
typedef int32_t *   _restrict I32PTR;
typedef uint32_t*   _restrict U32PTR;
typedef int16_t  *  _restrict I16PTR;
typedef uint16_t  * _restrict U16PTR;
typedef int8_t  *   _restrict I08PTR;
typedef uint8_t *   _restrict U08PTR;
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
