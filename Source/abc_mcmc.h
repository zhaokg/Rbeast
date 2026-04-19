#pragma once

#include "abc_datatype.h"
#include "abc_mem.h"

//clang complains "warning: redefinition of typedef 'MemPointers' is a C11 feature [-Wtypedef-redefinition]"
//typedef struct MemPointers MemPointers;

typedef struct _CRED_INTERVAL {
	F32PTR CI95;
	F32PTR minValPerStrip;
	I32PTR minIdxPerStrip;
	I32PTR whichStripHasMin;

	F32PTR CI05;                // = r_malloc(sizeof(F32)*N*numCISample);
	F32PTR maxValPerStrip;      //= r_malloc(sizeof(F32)*N * numStrips);
	I32PTR maxIdxPerStrip;      // = r_malloc(sizeof(int)*N * numStrips)	
	I32PTR whichStripHasMax;    //= r_malloc(sizeof(int)*N);

	F32PTR result;     //size(result)=2*N;
	F32PTR newDataRow; //size(newDataRow)=N;
	U32    N;
	U32    samplesInserted; // Must be reset to zero at the start of cacluation
} CI_RESULT;

typedef struct _CI {
	//U32       subsampleFraction_x_INT32MAX;// = 4294967295UL;
	U16       subsampleFraction_x_INT16MAX;// = 65535UL;
	U32       nSamples;
	U32       nStrips;

	I32PTR    SamplesPerStrip;
	I32PTR    OffsetsPerStrip;
} CI_PARAM;

 


#define _inout_ 
void ConstructCIStruct(
	F32 alpahLevel, I32 MCMC_SAMPLES, I32 N, I32 numCIVars, MemPointers* MEM,
	U08PTR _inout_  fastCIComputation, CI_PARAM* _out_ ciInfo, CI_RESULT* _out_ CI);

extern void InsertNewRowToUpdateCI(CI_PARAM* _restrict info, CI_RESULT* _restrict ci);
#define  RANDINTEGER(x0,x1,SEED)             ((x0)+ (SEED) % ( (x1)- (x0)+1))
#define  RANDI08(a,b)						 RANDINTEGER(a,b,*rnd08++)
#define  RANDI16(a,b)					  	 RANDINTEGER(a,b,*rnd16++)
#define  RANDI32(a,b)						 RANDINTEGER(a,b,*rnd32++)
#define  RANDINT(x0,x1,SEED)                ( x0==x1 ? x1 : x0 + (I32) ( (U32)(SEED) % (U32) ((x1)-(x0)+1) ) )

#define  RANDI_Skip_MIDNUM(rndInt, L, M, U)  out =( out=RANDI32(L,U-1),       out<(M)? out:(out+1) )
#define  RANDINT_SKIPONE(out,L,M,U,SEED)     out =( out=RANDINT(L,U-1, SEED), out<(M)? out:(out+1) )
 