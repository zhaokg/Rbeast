#pragma once
#include "abc_001_config.h"
#include <stdio.h>
#include <inttypes.h>
#include "abc_datatype.h"
#define MAX_K			300
#define MAX_SEG			30
#define MAX_RAND_NUM	5000
 typedef struct YINFO {
	float   yMean;
	float   yStd;
	F32PTR   Y;
	float    YtY;
	uint32_t n;
	uint32_t nMissing;
	U32PTR  rowsMissing;
} YINFO;
struct ModelPar {
	float alpha_1,alpha_2,del_1,del_2,sig2,prec[3];
};
typedef struct RESULT {
	F32PTR  time;
	F32PTR  sN,tN;
	I32PTR sNProb,tNProb;
	I32PTR sProb,tProb;
	F32PTR s,t,b;
	F32PTR sCI,tCI,bCI;
	F32PTR  sSD,tSD,bSD;
	F32PTR  marg_lik;	
	F32PTR  sig2;
	I32PTR bsign;
	I32PTR  horder;
	I32PTR torder;
	F32PTR  tcp,scp;
	F32PTR  tcpCI,scpCI;
} RESULT;
struct BASIS {
	float marg_lik;
	float alpha_star;
	F32PTR XtX,XtY,post_P_U;
	F32PTR beta_mean,beta;
	U16PTR S;
	U08PTR sOrder;
	I16PTR sks,ske;
	U16PTR T;
	U08PTR tOrder;
	I16PTR tks,tke;
	U08PTR termType;
	int16_t  sKnotNum;
	int16_t  tKnotNum;
	int16_t  k,k_SN,k_const;
} ;
typedef struct BASIS BASIS;
typedef struct Options {
	void * input,* output;
	char  inputType,outputType;
	int   M,L;
	int8_t   timeDimensionIndex;
	int8_t   isInput3DStack;
	char  *inputFile;
	char  *outputFolder;
	float *yInputData;
	uint64_t seed; 	
	bool   isSingleyInput;
	char   separator;
	uint32_t totalPixelNum;
	float    startTime,timeInterval;
	float    period;
	float    ridgeFactor;
	float    omissionValue;
	int32_t  N,Npad,Npad16;
	uint8_t  minSeasonOrder,maxSeasonOrder,minTrendOrder,maxTrendOrder;
	uint16_t minSepDist_Trend,minSepDist_Season;
	uint16_t maxKnotNum_Trend,maxKnotNum_Season;
	uint16_t maxMoveStepSize;	
	float    resamplingTrendOrderProb;
	float    resamplingSeasonOrderProb;
	uint32_t burnin,samples,chainNumber;
	uint16_t thinningFactor;
	float    alphaLevel;
	bool     computeSlopeSign;
	bool     computeHarmonicOrder;
	bool     computeTrendOrder;
	bool     computeChangepoints;
	bool	 inputFromDisk,outputToDisk;
	bool	 computeCredible;
	bool	 fastCIComputation;
	bool	 printToScreen;
	int16_t	 printCharLen;
} Options;
typedef struct FILE_LIST
{
	FILE *sN,*tN,*sNProb,*tNProb,*sProb,*tProb,*s,*sCI,*sSD;
	FILE *t,*tCI,*tSD,*b,*bCI,*bSD;
	FILE * bsign;
	FILE *infp;
} FILE_LIST;
enum  MOVETYPE { birth,death,merge,move,chorder};
static INLINE void zeroOut_Xmars_zero(F32PTR Xt_mars,F32PTR Xt_zeroBackup,
	U32PTR rowsMissing,uint32_t nMissing,int32_t N,int32_t Npad,int k)
{
	register  ptrdiff_t tmpidx;
	for (rI32 j=k; j>0; j--)
	{
		for (rI32 i=nMissing; i>0; i--)		{
			tmpidx=(*rowsMissing++);
			tmpidx--;
			*Xt_zeroBackup++=Xt_mars[tmpidx];
			Xt_mars[tmpidx]=0.f;
		}
		rowsMissing=rowsMissing - nMissing;
		Xt_mars=Xt_mars+Npad;
	}
}
static INLINE void zeroOut_Xmars_fill(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,
	uint32_t nMissing,int32_t N,int32_t Npad,int k)
{
	register size_t i,j,tmpidx;
	for (j=1; j <=k; j++)
	{
		for (i=1; i <=nMissing; i++)
		{
			tmpidx=rowsMissing[i - 1] - 1;
			Xt_mars[tmpidx]=*Xt_zeroBackup++;
		}
		Xt_mars+=Npad;
	}
}
extern void zeroOut_Xmars_zero0(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,int32_t nMissing,int32_t N,int32_t Npad,int32_t k);
extern void zeroOut_Xmars_zero1(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,int32_t nMissing,int32_t N,int32_t Npad,int32_t k);
extern void zeroOut_Xmars_zero2(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,int32_t nMissing,int32_t N,int32_t Npad,int32_t k);
extern void zeroOut_Xmars_zero3(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,int32_t nMissing,int32_t N,int32_t Npad,int32_t k);
extern int32_t find_index_by_csum(rU08PTR good,rI32 N,rI32 randIdx);
extern int32_t int8_arr_sum(rU08PTR good,rI32 N);
#include "abc_mem_ext.h"
extern  void print_error(int code,MemPointers *MEM);
extern Options     *GLOBAL_OPTIONS;
extern RESULT      *GLOBAL_RESULT;
#if M_INTERFACE==1
#define NULL_RET 
#elif R_INTERFACE==1
#define NULL_RET R_NilValue
#endif
