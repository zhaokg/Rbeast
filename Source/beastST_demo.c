#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
#if defined(WIN64_OS) 
#if defined(MSVC_COMPILER)
#include <windows.h>               
#include "intrin.h"                
#endif
#include <stdio.h>	               
#include <string.h>	               
#include <time.h>
#include <math.h>
#include "abc_001_config.h"
#include "abc_mem.h"              
#include "abc_common.h"           
#include "abc_blas_lapack_lib.h"  
#include "beast_lib.h"  
#include "beast_common.h"
#if MYRAND_LIBRARY==1
#include "abc_rand_pcg.h"
#elif MKLRAND_LIBRARY==1
#include "abc_rand_mkl.h"
VSLStreamStatePtr stream;  
#endif
#ifdef __MACH__
#include <mach/mach_time.h>
#endif
static F32PTR GlobalMEMBuf_1st=NULL;
static F32PTR GlobalMEMBuf_2nd=NULL;
#if PTHREAD_INOUT==1
struct THREADPAR_WRITE {
	int		curIdx;
	struct RESULT output;
	Options *opt;
};
static pthread_mutex_t MUTEX_WRITE;
static pthread_cond_t  CONDITION_WRITE;
static bool			   DATA_AVAILABLE_WRITE=false;
static void *WRITE(void* arg)
{
	struct THREADPAR_WRITE *par=(struct THREADPAR_WRITE *) arg;
	Options *opt=par->opt;
	int M=opt->M;
	if (M==0)	 return NULL;
	r_printf("Entering the WRITE thread...... \n");
	char *outPath=opt->outputFolder;
	char fn[230];
	strcpy(fn,outPath); 	strcat(fn,"/sN");		FILE *f_sN=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/tN");  	FILE *f_tN=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/sNProb");  FILE *f_sNProb=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/tNProb");  FILE *f_tNProb=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/sProb");  	FILE *f_sProb=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/tProb");  	FILE *f_tProb=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/s");  		FILE *f_s=fopen(fn,"wb+");
	FILE *f_sCI; if (opt->computeCredible)  strcpy(fn,outPath),strcat(fn,"/sCI"),f_sCI=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/sSD");  	FILE *f_sSD=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/t");  		FILE *f_t=fopen(fn,"wb+");
	FILE *f_tCI; if (opt->computeCredible)  strcpy(fn,outPath),strcat(fn,"/tCI"),f_tCI=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/tSD");  	FILE *f_tSD=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/b");  		FILE *f_b=fopen(fn,"wb+");
	FILE *f_bCI; 	if (opt->computeCredible)  strcpy(fn,outPath),strcat(fn,"/bCI"),f_bCI=fopen(fn,"wb+");
	strcpy(fn,outPath); 	strcat(fn,"/bSD");  	FILE *f_bSD=fopen(fn,"wb+");
	FILE *f_bsign; 	if (opt->computeSlopeSign) strcpy(fn,outPath),strcat(fn,"/bsign"),f_bsign=fopen(fn,"wb+");
	while (true)
	{
		pthread_mutex_lock(&MUTEX_WRITE);
		if (!DATA_AVAILABLE_WRITE)
		{
			while (!DATA_AVAILABLE_WRITE)
				pthread_cond_wait(&CONDITION_WRITE,&MUTEX_WRITE);
		}
		else if (DATA_AVAILABLE_WRITE)
		{
			int N=opt->N;
			r_printf("Starting to write output ..%d\n",par->curIdx);
			fwrite(par->output.sN,sizeof(float),1,f_sN);
			fwrite(par->output.tN,sizeof(float),1,f_tN);
			fwrite(par->output.sNProb,sizeof(int32_t),opt->maxKnotNum_Season+1,f_sNProb);
			fwrite(par->output.tNProb,sizeof(int32_t),opt->maxKnotNum_Trend+1,f_tNProb);
			fwrite(par->output.sProb,sizeof(int32_t),N,f_sProb);
			fwrite(par->output.tProb,sizeof(int32_t),N,f_tProb);
			fwrite(par->output.s,sizeof(float),N,f_s);
			if (opt->computeCredible)  fwrite(par->output.sCI,sizeof(float),N * 2,f_sCI);
			fwrite(par->output.sSD,sizeof(float),N,f_sSD);
			fwrite(par->output.t,sizeof(float),N,f_t);
			if (opt->computeCredible)  fwrite(par->output.tCI,sizeof(float),N * 2,f_tCI);
			fwrite(par->output.tSD,sizeof(float),N,f_tSD);
			fwrite(par->output.b,sizeof(float),N,f_b);
			if (opt->computeCredible)  fwrite(par->output.bCI,sizeof(float),N * 2,f_bCI);
			fwrite(par->output.bSD,sizeof(float),N,f_bSD);
			if (opt->computeSlopeSign)  fwrite(par->output.bsign,sizeof(float),N,f_bsign);
			r_printf("Finished writing output ..%d\n",par->curIdx);
			if (par->curIdx==M)
			{
				DATA_AVAILABLE_WRITE=false;
				pthread_mutex_unlock(&MUTEX_WRITE);
				pthread_cond_signal(&CONDITION_WRITE);
				break;
			}
			DATA_AVAILABLE_WRITE=false;
			pthread_cond_signal(&CONDITION_WRITE);
		}
		pthread_mutex_unlock(&MUTEX_WRITE);
	}
	fclose(f_sN);
	fclose(f_tN);
	fclose(f_sNProb);
	fclose(f_tNProb);
	fclose(f_sProb);
	fclose(f_tProb);
	fclose(f_s);
	fclose(f_t);
	fclose(f_b);
	if (opt->computeCredible) fclose(f_sCI);
	if (opt->computeCredible) fclose(f_tCI);
	if (opt->computeCredible) fclose(f_bCI);
	if (opt->computeSlopeSign) fclose(f_bsign);
	fclose(f_sSD);
	fclose(f_tSD);
	fclose(f_bSD);
	r_printf("Quitting the WRITE thread...... \n");
	return NULL;
}
struct THREADPAR_READ {
	int			curIdx;
	Options		*opt;
	YINFO input;
};
static pthread_mutex_t MUTEX_READ;
static pthread_cond_t  CONDITION_READ;
static bool  DATA_AVAILABLE_READ=false;
static void *READ(void* arg)
{
	struct THREADPAR_READ *par=(struct THREADPAR_READ *) arg;
	Options *opt=par->opt;
	int		M=opt->M;
	int		N=opt->N;
	float  omissionValue=opt->omissionValue;
	if (M==0)
	{
		return NULL;
	}
	r_printf("Entering the READ thread...... \n");
	FILE * infp=fopen(opt->inputFile,"rb");
	if (infp==NULL)  return NULL;
	float *tmp_MEMBuf=r_malloc(sizeof(float)*N);
	while (true)
	{
		pthread_mutex_lock(&MUTEX_READ);
		if (DATA_AVAILABLE_READ)
		{
			while (DATA_AVAILABLE_READ)
				pthread_cond_wait(&CONDITION_READ,&MUTEX_READ);
		}
		else if (!DATA_AVAILABLE_READ)
		{
			if (par->curIdx > M)
			{
				DATA_AVAILABLE_READ=false;
				pthread_mutex_unlock(&MUTEX_READ);
				pthread_cond_signal(&CONDITION_READ);
				break;
			}
			r_printf("Starting to read a time series from the file ...%d \n",par->curIdx);
			fseek(infp,(par->curIdx - 1)*N*sizeof(float),SEEK_SET);
			fread(par->input.Y,sizeof(float),N,infp);
			float *yInput=par->input.Y;
			float *buf_start=tmp_MEMBuf; 
			int nMissing=0;
			for (int32_t i=1; i <=N; i++)
			{
				if (*yInput !=*yInput||fabs(*yInput - omissionValue) < 1e-5)
					par->input.rowsMissing[nMissing++]=i;
				else
					*buf_start++=*yInput;
				yInput++;
			}
			par->input.nMissing=nMissing;
			par->input.n=N - nMissing;
			r_ippsMeanStdDev_32f(tmp_MEMBuf,par->input.n,&(par->input.yMean),&(par->input.yStd),ippAlgHintAccurate);
			normalize(tmp_MEMBuf,par->input.n);
			yInput=par->input.Y;
			buf_start=tmp_MEMBuf;
			for (int32_t i=1; i <=N; i++)
			{
				if (*yInput !=*yInput||fabs(*yInput - omissionValue) < 1e-5)
					*yInput++=0;
				else
					*yInput++=*buf_start++;
			}
			if (par->curIdx==M)
			{
				DATA_AVAILABLE_READ=1;
				pthread_mutex_unlock(&MUTEX_READ);
				pthread_cond_signal(&CONDITION_READ);
				break;
			}
			DATA_AVAILABLE_READ=true;
			pthread_cond_signal(&CONDITION_READ);
		}
		pthread_mutex_unlock(&MUTEX_READ);
	}
	r_free(tmp_MEMBuf);
	fclose(infp);
	r_printf("Quiting the READ thread...... \n");
	return NULL;
}
#endif
static uint64_t elapsedTime;
#if defined(MSVC_COMPILER)
static LARGE_INTEGER t1,t2;
static LARGE_INTEGER Frequency;
#elif ( defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER) ) && !(defined(__APPLE__)||defined(__MACH__))
static struct timespec t1,t2;
#elif defined(__MACH__)
static uint64_t t1,t2;
#endif
#include "demo.h" 
#undef NULL_RET
#define NULL_RET 0
extern Options     *GLOBAL_OPTIONS;
extern RESULT      *GLOBAL_RESULT;
DWORD WINAPI beastST_demo(__in LPVOID lpParameter)
{
	GLOBAL_OPTIONS=((LParam*)lpParameter)->GLOBAL_OPTIONS;
	GLOBAL_RESULT=((LParam*)lpParameter)->GLOBAL_RESULT;
#if MKL_LIBRARY==1
#endif
#if defined(MSVC_COMPILER)
#else	
#endif
	float nan;
	nan=1.00000e300;
	nan=nan*nan*0.f;
	MemPointers MEM=(MemPointers) { .init=mem_init };
	MEM.init(&MEM);
	Options	opt;
	memcpy(&opt,GLOBAL_OPTIONS,sizeof(Options));
	if (gData.N <=0)
	{
		EnterCriticalSection(&gData.cs);
		gData.N=opt.N;
		memcpy(&gData.opt,&opt,sizeof(opt));
		gData.optStatus=NoUpdate;
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);
	}
	if (gData.optStatus==NeedUpdate)
	{
		EnterCriticalSection(&gData.cs);
		float * tmp=opt.yInputData;
		memcpy(&opt,&gData.opt,sizeof(opt));
		opt.yInputData=tmp;
		LeaveCriticalSection(&gData.cs);
	}
	if (gData.plotData[0][0]==NULL)
	{
		EnterCriticalSection(&gData.cs);
		while (gData.plotData[0][0]==NULL)
			SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
		LeaveCriticalSection(&gData.cs);
	}
	uint32_t  numCISample=0;
	uint32_t  numStrips=0;
	I32PTR    rowsPerStrip=NULL;
	I32PTR    startIdxOfStrip=NULL;
	uint32_t  subsampleFraction_x_INT32MAX=4294967295;
	if (opt.computeCredible)
	{
		rF32 alpha=(1 - opt.alphaLevel)/2.f;
		if (opt.fastCIComputation)
		{
			numCISample=100;
			float subsampleFraction=(float)numCISample/alpha/opt.samples;
			if (subsampleFraction >=0.99f)
			{
				opt.fastCIComputation=false;
			}
			else
			{
				subsampleFraction_x_INT32MAX=(uint32_t)((double)subsampleFraction * (double)4294967295LL);
			}
		}
		if (!opt.fastCIComputation)
		{
			numCISample=(uint32_t)((float)opt.samples*alpha);
		}
		uint32_t stripWidth=(uint32_t)ceil(fast_sqrt((float)numCISample));
		numStrips=(uint32_t)floor(numCISample/stripWidth);
		rowsPerStrip=MyALLOC(MEM,numStrips,int,0);
		startIdxOfStrip=MyALLOC(MEM,numStrips,int,0);
		uint32_t stripTmpNum=0;
		for (rU32 i=1; i <=numStrips; i++)
		{
			rowsPerStrip[i - 1]=(i !=numStrips) ? stripWidth : numCISample - stripTmpNum;
			startIdxOfStrip[i - 1]=stripTmpNum;
			stripTmpNum+=stripWidth;
		}
	}
	F32PTR credT_upper=NULL;
	I32PTR min_credT_upper_IdxInStrip=NULL;
	F32PTR min_credT_upperStrip=NULL;
	I32PTR whichMin_credT_upper=NULL;
	F32PTR credT_lower=NULL;
	I32PTR max_credT_lower_IdxInStrip=NULL;
	F32PTR max_credT_lowerStrip=NULL;
	I32PTR whichMax_credT_lower=NULL;
	F32PTR credS_upper=NULL;
	I32PTR min_credS_upper_IdxInStrip=NULL;
	F32PTR min_credS_upperStrip=NULL;
	I32PTR whichMin_credS_upper=NULL;
	F32PTR credS_lower=NULL;
	I32PTR max_credS_lower_IdxInStrip=NULL;
	F32PTR max_credS_lowerStrip=NULL;
	I32PTR whichMax_credS_lower=NULL;
	F32PTR credB_upper=NULL;
	I32PTR min_credB_upper_IdxInStrip=NULL;
	F32PTR min_credB_upperStrip=NULL;
	I32PTR whichMin_credB_upper=NULL;
	F32PTR credB_lower=NULL; 
	I32PTR max_credB_lower_IdxInStrip=NULL;
	F32PTR max_credB_lowerStrip=NULL;
	I32PTR whichMax_credB_lower=NULL;
	if (opt.computeCredible)
	{
		credT_upper=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		min_credT_upper_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N*numStrips,0);
		min_credT_upperStrip=MEM.alloc(&MEM,sizeof(float)*opt.N*numStrips,0);
		whichMin_credT_upper=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
		credT_lower=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		max_credT_lower_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N * numStrips,0);
		max_credT_lowerStrip=MEM.alloc(&MEM,sizeof(float)*opt.N * numStrips,0);
		whichMax_credT_lower=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
		credS_upper=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		min_credS_upper_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N*numStrips,0);
		min_credS_upperStrip=MEM.alloc(&MEM,sizeof(float)*opt.N*numStrips,0);
		whichMin_credS_upper=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
		credS_lower=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		max_credS_lower_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N * numStrips,0);
		max_credS_lowerStrip=MEM.alloc(&MEM,sizeof(float)*opt.N * numStrips,0);
		whichMax_credS_lower=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
		credB_upper=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		min_credB_upper_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N*numStrips,0);
		min_credB_upperStrip=MEM.alloc(&MEM,sizeof(float)*opt.N*numStrips,0);
		whichMin_credB_upper=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
		credB_lower=MEM.alloc(&MEM,sizeof(float)*opt.N*numCISample,0);
		max_credB_lower_IdxInStrip=MEM.alloc(&MEM,sizeof(int)*opt.N * numStrips,0);
		max_credB_lowerStrip=MEM.alloc(&MEM,sizeof(float)*opt.N * numStrips,0);
		whichMax_credB_lower=MEM.alloc(&MEM,sizeof(int)*opt.N,0);
	}
	BASIS   BASISTEMP1,BASISTEMP2;
	BASIS * _restrict basis=&BASISTEMP1;
	BASIS * _restrict basis_prop=&BASISTEMP2;
	{
		BASISTEMP1.XtX=MyALLOC(MEM,MAX_K*MAX_K,float,64);				
		BASISTEMP2.XtX=MyALLOC(MEM,MAX_K*MAX_K,float,64);
		BASISTEMP1.XtY=MyALLOC(MEM,MAX_K,float,64);
		BASISTEMP2.XtY=MyALLOC(MEM,MAX_K,float,64);
		BASISTEMP1.post_P_U=MyALLOC(MEM,MAX_K*MAX_K,float,64);
		BASISTEMP2.post_P_U=MyALLOC(MEM,MAX_K*MAX_K,float,64);
		BASISTEMP1.beta_mean=MyALLOC(MEM,MAX_K,float,64);
		BASISTEMP2.beta_mean=MyALLOC(MEM,MAX_K,float,64);
		BASISTEMP1.beta=MyALLOC(MEM,MAX_K,float,64);
		BASISTEMP2.beta=MyALLOC(MEM,MAX_K,float,64);
		if ((opt.maxKnotNum_Season+2L) * 2 * sizeof(uint16_t) <=64)
		{
			BASISTEMP1.S=MyALLOC(MEM,64,uint16_t,64);			
			*BASISTEMP1.S++=1L;
			BASISTEMP2.S=(uint16_t *)BASISTEMP1.S+(opt.maxKnotNum_Season+1);
			*BASISTEMP2.S++=1L;
		}
		else
		{
			BASISTEMP1.S=MyALLOC(MEM,opt.maxKnotNum_Season+2,uint16_t,64);
			*BASISTEMP1.S++=1L;
			BASISTEMP2.S=MyALLOC(MEM,opt.maxKnotNum_Season+2,uint16_t,64); 
			*BASISTEMP2.S++=1L;
		}
		if ((opt.maxKnotNum_Season+1L) * 2 * sizeof(char) <=64)
		{
			BASISTEMP1.sOrder=MyALLOC(MEM,64,uint8_t,64);
			BASISTEMP2.sOrder=(uint8_t *)BASISTEMP1.sOrder+(opt.maxKnotNum_Season+1);
		}
		else
		{
			BASISTEMP1.sOrder=MyALLOC(MEM,opt.maxKnotNum_Season+1,char,64); 
			BASISTEMP2.sOrder=MyALLOC(MEM,opt.maxKnotNum_Season+1,char,64); 
		}
		BASISTEMP1.sks=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP2.sks=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP1.ske=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP2.ske=MyALLOC(MEM,MAX_K,int16_t,64);
		if ((opt.maxKnotNum_Trend+2L) * 2 * sizeof(uint16_t) <=64)
		{
			BASISTEMP1.T=MyALLOC(MEM,64,uint16_t,64);
			*BASISTEMP1.T++=1L;
			BASISTEMP2.T=(uint16_t *)BASISTEMP1.T+opt.maxKnotNum_Trend+1;
			*BASISTEMP2.T++=1L;
		}
		else
		{
			BASISTEMP1.T=MyALLOC(MEM,opt.maxKnotNum_Trend+2,uint16_t,64);  
			*BASISTEMP1.T++=1L;
			BASISTEMP2.T=MyALLOC(MEM,opt.maxKnotNum_Trend+2,uint16_t,64);
			*BASISTEMP2.T++=1L;
		}
		if ((opt.maxKnotNum_Trend+1) * 2 * sizeof(char) <=64)
		{
			BASISTEMP1.tOrder=MyALLOC(MEM,64,uint8_t,64);
			BASISTEMP2.tOrder=(uint8_t *)BASISTEMP1.tOrder+(opt.maxKnotNum_Trend+1);
		}
		else
		{
			BASISTEMP1.tOrder=MyALLOC(MEM,opt.maxKnotNum_Trend+1,uint8_t,64);
			BASISTEMP2.tOrder=MyALLOC(MEM,opt.maxKnotNum_Trend+1,uint8_t,64);
		}
		BASISTEMP1.tks=MyALLOC(MEM,MAX_K,int16_t,64); 
		BASISTEMP2.tks=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP1.tke=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP2.tke=MyALLOC(MEM,MAX_K,int16_t,64);
		BASISTEMP1.termType=MyALLOC(MEM,MAX_K,uint8_t,64);
		BASISTEMP2.termType=MyALLOC(MEM,MAX_K,uint8_t,64);
	}
	F32PTR Xt_mars=MyALLOC(MEM,opt.Npad*MAX_K,float,64); 
	F32PTR Xnewterm=MyALLOC(MEM,opt.Npad*MAX_K,float,64); 
	F32PTR Xt_zeroBackup=MyALLOC(MEM,opt.Npad*MAX_K,float,64);
	U08PTR goodT=NULL,goodS=NULL;
	uint32_t goodT_num,goodS_num;
	goodT=MyALLOC(MEM,opt.Npad16,uint8_t,64);  
	goodS=MyALLOC(MEM,opt.Npad16,uint8_t,64);
	memset(goodT+opt.N,0L,opt.Npad16 - opt.N);
	memset(goodS+opt.N,0L,opt.Npad16 - opt.N);
	F32PTR SEASON_TERMS=NULL;
	F32PTR TREND_TERMS=NULL;
	F32PTR INV_SQR=NULL;
	F32PTR COEFF_A=NULL;
	F32PTR COEFF_B=NULL;
	F32PTR SEASON_SQR_CSUM=NULL;
	SEASON_TERMS=MyALLOC(MEM,opt.N*2* opt.maxSeasonOrder,float,64);  
	TREND_TERMS=MyALLOC(MEM,opt.N*(opt.maxTrendOrder+1L),float,64);
	INV_SQR=MyALLOC(MEM,opt.N,float,0);   
	COEFF_A=MyALLOC(MEM,opt.N,float,0);
	COEFF_B=MyALLOC(MEM,opt.N,float,0);
	SEASON_SQR_CSUM=MyALLOC(MEM,(opt.N+1L) * 2 * opt.maxSeasonOrder,float,64);  
	preCompute_Xmars_terms_fast( 
		SEASON_TERMS,SEASON_SQR_CSUM,TREND_TERMS, 
		INV_SQR,COEFF_A,COEFF_B,
		opt.N,opt.period,opt.maxSeasonOrder,opt.maxTrendOrder);
	RESULT resultChain={ NULL,};
	RESULT result={ NULL,};
	RESULT matOutput=(RESULT) { NULL,};
	allocate_single_output(&resultChain,&opt,&MEM);
	allocate_single_output(&result,&opt,&MEM);
	memcpy(&matOutput,GLOBAL_RESULT,sizeof(RESULT));
#if PTHREAD_INOUT==1
	pthread_t		THREADID_WRITE;
	pthread_attr_t	attr;
	struct THREADPAR_WRITE threadParWrite;
	if (opt.outputType=='F')
	{
		threadParWrite.opt=&opt;
		threadParWrite.output.sN=MEM.alloc(&MEM,sizeof(float)* 1,0);
		threadParWrite.output.tN=MEM.alloc(&MEM,sizeof(float)* 1,0);
		threadParWrite.output.sProb=MEM.alloc(&MEM,sizeof(int32_t)*opt.N,0);
		threadParWrite.output.tProb=MEM.alloc(&MEM,sizeof(int32_t)*opt.N,0);
		threadParWrite.output.sNProb=MEM.alloc(&MEM,sizeof(int32_t)*(opt.maxKnotNum_Season+1),0);
		threadParWrite.output.tNProb=MEM.alloc(&MEM,sizeof(int32_t)*(opt.maxKnotNum_Trend+1),0);
		threadParWrite.output.s=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		if (opt.computeCredible) { threadParWrite.output.sCI=MEM.alloc(&MEM,sizeof(float)*opt.N * 2,0); }
		threadParWrite.output.sSD=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		threadParWrite.output.t=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		if (opt.computeCredible) { threadParWrite.output.tCI=MEM.alloc(&MEM,sizeof(float)*opt.N * 2,0); }
		threadParWrite.output.tSD=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		threadParWrite.output.b=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		if (opt.computeCredible) { threadParWrite.output.bCI=MEM.alloc(&MEM,sizeof(float)*opt.N * 2,0); }
		threadParWrite.output.bSD=MEM.alloc(&MEM,sizeof(float)*opt.N,0);
		if (opt.computeSlopeSign) { threadParWrite.output.bsign=MEM.alloc(&MEM,sizeof(float)*opt.N,0); }
		pthread_mutex_init(&MUTEX_WRITE,NULL);
		pthread_cond_init(&CONDITION_WRITE,NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
		DATA_AVAILABLE_WRITE=false;
		pthread_create(&THREADID_WRITE,&attr,WRITE,(void *)&threadParWrite);
		pthread_attr_destroy(&attr);
	}
#else
	struct FILE_LIST file;
	memset(&file,0,sizeof(FILE_LIST));
	if (opt.outputType=='F')
	{
		char fn[230];
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/sN"); 	    file.sN=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/tN");  	file.tN=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/sNProb");  file.sNProb=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/tNProb");  file.tNProb=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/sProb");  	file.sProb=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/tProb");  	file.tProb=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/s");  	    file.s=fopen(fn,"wb+");
		if (opt.computeCredible)  		strcpy(fn,opt.outputFolder),strcat(fn,"/sCI"),file.sCI=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/sSD");  	file.sSD=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/t");  	    file.t=fopen(fn,"wb+");
		if (opt.computeCredible) 		strcpy(fn,opt.outputFolder),strcat(fn,"/tCI"),file.tCI=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/tSD");  	file.tSD=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/b");  	    file.b=fopen(fn,"wb+");
		if (opt.computeCredible)	    strcpy(fn,opt.outputFolder),strcat(fn,"/bCI"),file.bCI=fopen(fn,"wb+");
		strcpy(fn,opt.outputFolder); 	strcat(fn,"/bSD");  	file.bSD=fopen(fn,"wb+");
		if (opt.computeSlopeSign)		strcpy(fn,opt.outputFolder),strcat(fn,"/bsign"),file.bsign=fopen(fn,"wb+");
	}
#endif
	F32PTR scaleFactorSeason;
	F32PTR scaleFactorTrend;
	scaleFactorSeason=MyALLOC(MEM,opt.maxKnotNum_Season+1,float,0);
	scaleFactorTrend=MyALLOC(MEM,opt.maxKnotNum_Trend+1,float,0); 
	GlobalMEMBuf_1st=Xnewterm; 
	GlobalMEMBuf_2nd=Xnewterm+(opt.maxKnotNum_Season+1)+(opt.maxKnotNum_Trend+1);
	preCompute_scale_factor(scaleFactorSeason,scaleFactorTrend,opt.N,opt.maxKnotNum_Season,\
		opt.maxKnotNum_Trend,opt.minSepDist_Season,opt.minSepDist_Trend,GlobalMEMBuf_1st,GlobalMEMBuf_2nd);
	YINFO	yInfo;
	yInfo.Y=MyALLOC(MEM,opt.Npad,float,64); 
	yInfo.rowsMissing=MyALLOC(MEM,opt.N,uint32_t,64);
#if PTHREAD_INOUT==1
	pthread_t				THREADID_READ;
	struct THREADPAR_READ	threadParRead;
	if (opt.inputFromDisk)
	{
		threadParRead.opt=&opt;
		threadParRead.input.Y=MEM.alloc(&MEM,sizeof(float)*opt.N,0);         
		threadParRead.input.rowsMissing=MEM.alloc(&MEM,sizeof(int)*opt.N,0);; 
		threadParRead.curIdx=1;
		pthread_mutex_init(&MUTEX_READ,NULL);
		pthread_cond_init(&CONDITION_READ,NULL);
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
		DATA_AVAILABLE_READ=0;
		pthread_create(&THREADID_READ,&attr,READ,(void *)&threadParRead);
		pthread_attr_destroy(&attr);
	}
#else
	if (opt.inputType=='F')
	{
		file.infp=fopen(opt.inputFile,"rb");		
	}
#endif
	struct {
		uint8_t birth;
		uint8_t death;
		uint8_t merge;
		uint8_t move;
	} T_propProb,S_propProb;
	if (opt.maxTrendOrder==opt.minTrendOrder)
		opt.resamplingTrendOrderProb=0;
	opt.resamplingTrendOrderProb=1 - opt.resamplingTrendOrderProb;
	T_propProb.birth=(uint8_t)(opt.resamplingTrendOrderProb *  0.33 * 127);
	T_propProb.death=(uint8_t)(opt.resamplingTrendOrderProb * (0.33+0.165) * 127);
	T_propProb.merge=(uint8_t)(opt.resamplingTrendOrderProb * (0.33+0.33) * 127);
	T_propProb.move=(uint8_t)(opt.resamplingTrendOrderProb * 1.0 * 127);
	if (opt.maxSeasonOrder==opt.minSeasonOrder)
		opt.resamplingSeasonOrderProb=0;
	opt.resamplingSeasonOrderProb=1 - opt.resamplingSeasonOrderProb;
	S_propProb.birth=(uint8_t)(opt.resamplingSeasonOrderProb *  0.33 * 127);
	S_propProb.death=(uint8_t)(opt.resamplingSeasonOrderProb * (0.33+0.165) * 127);
	S_propProb.merge=(uint8_t)(opt.resamplingSeasonOrderProb * (0.33+0.33) * 127);
	S_propProb.move=(uint8_t)(opt.resamplingSeasonOrderProb * 1.0 * 127);
#if defined(MSVC_COMPILER)
	QueryPerformanceFrequency(&Frequency); 
#endif
#if	defined(MSVC_COMPILER)
	if (opt.seed==0) 	opt.seed=GetTickCount64();
#elif (defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER)) && ! (defined(__APPLE__)||defined(__MACH__))
	struct timespec tmpTimer;
	clock_gettime(CLOCK_REALTIME,&tmpTimer);
	if (opt.seed==0) 	opt.seed=tmpTimer.tv_sec * 1000000000LL+tmpTimer.tv_nsec;
#elif defined(__MACH__) 
	if (opt.seed==0) opt.seed=mach_absolute_time();
#endif
	r_vslNewStream(&stream,VSL_BRNG_MT19937,opt.seed);
	for (uint32_t pixelIndex=1; pixelIndex <=opt.M; pixelIndex++)
	{
		if (pixelIndex%200==0)
			r_printf("Processing the%d -th time series out of%d \n ..",pixelIndex,opt.M);
		elapsedTime=0;
#if defined(MSVC_COMPILER)
		QueryPerformanceCounter(&t1);
#elif (defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER))  && ! (defined(__APPLE__)||defined(__MACH__))
		clock_gettime(CLOCK_REALTIME,&t1);
#endif
		int32_t N=opt.N;
		int32_t Npad=opt.Npad;
		if (opt.inputType=='F')
		{
#if PTHREAD_INOUT==1
			pthread_mutex_lock(&MUTEX_READ);
			if (DATA_AVAILABLE_READ)
			{
				memcpy(yInfo.rowsMissing,threadParRead.input.rowsMissing,sizeof(int)* threadParRead.input.nMissing);
				r_cblas_scopy(N,threadParRead.input.Y,1,yInfo.Y,1);
				yInfo.n=threadParRead.input.n;
				yInfo.nMissing=threadParRead.input.nMissing;
				yInfo.yMean=threadParRead.input.yMean;
				yInfo.yStd=threadParRead.input.yStd;
				DATA_AVAILABLE_READ=false;
				threadParRead.curIdx=pixelIndex+1;
				pthread_mutex_unlock(&MUTEX_READ);
				pthread_cond_signal(&CONDITION_READ);
			}
			else if (!DATA_AVAILABLE_READ)
			{
				while (!DATA_AVAILABLE_READ)
					pthread_cond_wait(&CONDITION_READ,&MUTEX_READ);
				memcpy(yInfo.rowsMissing,threadParRead.input.rowsMissing,sizeof(int)* threadParRead.input.nMissing);
				r_cblas_scopy(N,threadParRead.input.Y,1,yInfo.Y,1);
				yInfo.n=threadParRead.input.n;
				yInfo.nMissing=threadParRead.input.nMissing;
				yInfo.yMean=threadParRead.input.yMean;
				yInfo.yStd=threadParRead.input.yStd;
				DATA_AVAILABLE_READ=0;
				threadParRead.curIdx=pixelIndex+1;
				pthread_mutex_unlock(&MUTEX_READ);
				pthread_cond_signal(&CONDITION_READ);
			}
#else
			GlobalMEMBuf_1st=Xnewterm;
			r_printf("Starting to read a time series from the file ...%d \n",pixelIndex);
			fseek(file.infp,(pixelIndex - 1)*opt.N*sizeof(float),SEEK_SET);
			rI32 readBytes=(int32_t) fread(GlobalMEMBuf_1st,sizeof(float),opt.N,file.infp);
			opt.yInputData=GlobalMEMBuf_1st;
			fetch_next_time_series(&yInfo,opt.yInputData,1,GlobalMEMBuf_1st+N,opt.isSingleyInput,N,opt.omissionValue);
#endif
		}
		else 
		{
			GlobalMEMBuf_1st=Xnewterm;
			fetch_next_time_series1(&yInfo,opt.input,pixelIndex,GlobalMEMBuf_1st,opt.inputType,N,opt.omissionValue);
		}
		void(*p_zeroOut_Xmars_zero)(F32PTR,F32PTR,U32PTR,int32_t,int32_t,int32_t,int32_t);
		{
			void(*funArray[4])(F32PTR,F32PTR,U32PTR,int32_t,int32_t,int32_t,int32_t);
			funArray[0]=&zeroOut_Xmars_zero0;
			funArray[1]=&zeroOut_Xmars_zero1;
			funArray[2]=&zeroOut_Xmars_zero2;
			funArray[3]=&zeroOut_Xmars_zero3;
			p_zeroOut_Xmars_zero=funArray[yInfo.nMissing%4];
		}
		EnterCriticalSection(&gData.cs);
		{
			int idx;
			r_ippsMaxIndx_32f(yInfo.Y,N,&gData.yMax,&idx);
			r_ippsMinIndx_32f(yInfo.Y,N,&gData.yMin,&idx);
			gData.yMin=gData.yMin - (gData.yMax - gData.yMin)/10;
			gData.yMax=gData.yMax+(gData.yMax - gData.yMin)/10;
			gData.t=resultChain.t;
			gData.s=resultChain.s;
			gData.y=yInfo.Y;
			gData.rowsMissing=yInfo.rowsMissing;
			gData.nMissing=yInfo.nMissing;
		}
		LeaveCriticalSection(&gData.cs);
		WakeConditionVariable(&gData.cv);
		struct ModelPar modelPar;
		modelPar.alpha_1=0.00000001f;
		modelPar.alpha_2=0.00000001f;
		modelPar.del_1=0.00000001f;
		modelPar.del_2=0.00000001f;
		modelPar.sig2=1.f;
		modelPar.prec[0]=modelPar.prec[1]=modelPar.prec[2]=10.f;
		int32_t accS[5]={ 0,0,0,0,0 };
		int32_t accT[5]={ 0,0,0,0,0 };
		uint32_t RND[MAX_RAND_NUM];
		float	 RNDGAMMA[MAX_RAND_NUM];
		uint8_t  RNDCHAR[MAX_RAND_NUM * 4];
		U32PTR rnd=RND;
		U32PTR RND_END=RND+MAX_RAND_NUM - 50;
		F32PTR rndgamma=RNDGAMMA;
		F32PTR RNDGAMMA_END=RNDGAMMA+MAX_RAND_NUM - 50;
		U08PTR rndchar=RNDCHAR;
		U08PTR RNDCHAR_END=RNDCHAR+MAX_RAND_NUM * 4 - 10;
		r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD,stream,MAX_RAND_NUM,(uint32_t *)RND);
		r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,stream,MAX_RAND_NUM,rndgamma,(modelPar.del_1+yInfo.n *0.5f),0,1);
		r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD,stream,MAX_RAND_NUM,(uint32_t *)RNDCHAR);
		zero_out_result_output(&opt,&result);
		uint32_t chainNumber=0;
		while (chainNumber< opt.chainNumber)
		{
			LARGE_INTEGER tStart,tEnd,tFrequency;
			QueryPerformanceCounter(&tStart);
			QueryPerformanceFrequency(&tFrequency);
			zero_out_result_output(&opt,&resultChain);
			basis->sKnotNum=0; 
			basis->sOrder[0]=(unsigned char)opt.minSeasonOrder+(unsigned char)(*rnd++)%(opt.maxSeasonOrder - opt.minSeasonOrder+1);
			basis->S[basis->sKnotNum]=N+1;
			basis->tKnotNum=0; 
			basis->tOrder[0]=(unsigned char)opt.minTrendOrder+(unsigned char)(*rnd++)%(opt.maxTrendOrder - opt.minTrendOrder+1);
			basis->T[basis->tKnotNum]=N+1;
			{
				convert_basis_both(basis);
				GlobalMEMBuf_1st=Xnewterm;
				evaluate_basis_both_fast(basis,opt.N,Xt_mars,Xt_zeroBackup,&yInfo,SEASON_TERMS,TREND_TERMS,&modelPar,GlobalMEMBuf_1st);
				findGoodKnotPositionFromBasis(basis,goodS,goodT,opt.N,opt.minSepDist_Season,opt.minSepDist_Trend);
				goodT_num=goodS_num=0;
				for (rI32 i=N; i >0; i--)
					goodS_num+=(*goodS++),goodT_num+=(*goodT++);
				goodS=goodS - N;
				goodT=goodT - N;
			}
			if (opt.printToScreen)
			{
				register char * _restrict string;
				string=(char*)Xnewterm;
				rI32 i=opt.printCharLen+9; 
				r_ippsSet_8u('#',string,i); 
				string[i+1 - 1]=0;
				r_printf("%s\n",string);
				r_ippsSet_8u('-',string,i); 
				char  *tmpCharEnd=string+i+1 - 1;
				sprintf(tmpCharEnd,\
					"TIMESEIRES #%04d being processed: "
					"A total of%05d rows omitted.",\
					pixelIndex,yInfo.nMissing);
				rI32 j=(int)strlen(tmpCharEnd);
				memcpy(string+(int)(i - j)/2 - 1,tmpCharEnd,j);
				string[0]='#';
				string[i - 1]='#';
				string[i+1 - 1]=0;
				r_printf("%s\n",string);
				r_ippsSet_8u(35,string,i);
				string[i+1 - 1]=0;
				r_printf("%s\n\n",string);
			}
			uint16_t MAX_K_STOP_ADDING_NEWseasonTERMS=MAX_K - (2 * opt.maxSeasonOrder);
			uint16_t MAX_K_STOP_ADDING_NEWtrendTERMS=MAX_K - (opt.maxTrendOrder+1);
			uint32_t ite=0;
			uint32_t sample=0;
			uint32_t subSampleIndex=0;
			F32PTR XtX_prop,XtY_prop,post_P_U_prop,beta_prop,beta_mean_prop;
			F32PTR XtX,XtY,post_P_U,beta,beta_mean;
			int32_t K,K_SN;
			int16_t tKnotNum,sKnotNum;
			{
				XtX_prop=basis_prop->XtX;
				XtY_prop=basis_prop->XtY;
				beta_mean_prop=basis_prop->beta_mean;
				beta_prop=basis_prop->beta;
				post_P_U_prop=basis_prop->post_P_U;
				XtX=basis->XtX;
				XtY=basis->XtY;
				beta_mean=basis->beta_mean;
				beta=basis->beta;
				post_P_U=basis->post_P_U;
				K=basis->k;
				K_SN=basis->k_SN;
				tKnotNum=basis->tKnotNum;
				sKnotNum=basis->sKnotNum;
			}
			while (sample < opt.samples)
			{
				ite++;
				if (rnd >=RND_END){
					rnd=RND;
					r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD,stream,MAX_RAND_NUM,(uint32_t *)rnd);
				}
				if (rndchar >=RNDCHAR_END)	{
					rndchar=RNDCHAR;
					r_viRngUniformBits32(VSL_RNG_METHOD_UNIFORMBITS32_STD,stream,MAX_RAND_NUM,(uint32_t *)rndchar);
				}
				if (rndgamma >=RNDGAMMA_END){
					rndgamma=RNDGAMMA;
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,stream,MAX_RAND_NUM,rndgamma,(modelPar.del_1+yInfo.n *0.5f),0.f,1.f);
				}
				int16_t k;
				int16_t newIdx,endIdx;
				int16_t k1_old,k2_old;
				int16_t k1_new,k2_new;
				int16_t K_newTerm;
				int32_t newTerm_startidx,newTerm_endidx;
				char     oldOrder,newOrder;
				uint16_t newKnot,oldKnot;
				uint8_t unifRnd;
				bool IsTrend;
				enum MOVETYPE flag;
				unifRnd=*rndchar++;
				IsTrend=unifRnd  >=128;
				unifRnd=unifRnd & 0x7F; 
				if (IsTrend)
				{
					if (unifRnd <=T_propProb.birth) 
					{
						flag=birth;
						if (tKnotNum >=opt.maxKnotNum_Trend||goodT_num==0)
							flag=move;
						if (K > MAX_K_STOP_ADDING_NEWtrendTERMS)
							flag=(tKnotNum==0) ? birth : move;
					}
					else if (unifRnd<=T_propProb.death) 
					{
						flag=tKnotNum==0 ? birth : death;
					}
					else if (unifRnd<=T_propProb.merge) 
					{
						if (tKnotNum >=2)
							flag=merge;
						else
							flag=tKnotNum==0 ? birth : death;
					}
					else if (unifRnd<=T_propProb.move)
					{
						flag=(tKnotNum==0) ? birth : move;
					}
					else
					{
						flag=chorder;
					} 
					basis_prop->sKnotNum=sKnotNum;
					memcpy(basis_prop->S,basis->S,sizeof(uint16_t)*(sKnotNum+1));
					memcpy(basis_prop->sOrder,basis->sOrder,sizeof(char)    *(sKnotNum+1));
					rU16PTR knotListNew=basis_prop->T;
					rU08PTR orderListNew=basis_prop->tOrder;
					rU16PTR knotListOld=basis->T;
					rU08PTR orderListOld=basis->tOrder;
					switch (flag)
					{
					case birth:
					{
								  rI32 randIdx;
								  randIdx=(*rnd++)%goodT_num+1;
								  newKnot=find_index_by_csum(goodT,opt.Npad16,randIdx);
								  rU16PTR knotList=basis->T;
								  for (newIdx=1; newKnot > *knotList++; newIdx++);
								  basis_prop->tKnotNum=tKnotNum+1;
								  memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
								  *(knotListNew+newIdx - 1)=newKnot;
								  memcpy(knotListNew+(newIdx+1) - 1,knotListOld+newIdx - 1,sizeof(uint16_t)*(tKnotNum - (newIdx - 1)));
								  memcpy(orderListNew,orderListOld,sizeof(char)*newIdx);
								  memcpy(orderListNew+(newIdx+1) - 1,orderListOld+newIdx - 1,sizeof(char)*(tKnotNum+1 - newIdx+1));
								  k1_old=basis->tks[(newIdx)-1];
								  k2_old=basis->tke[(newIdx)-1];
								  endIdx=newIdx+1;
								  break;
					}
					case death:
					{
								  newIdx=(*rnd++)%tKnotNum+1;
								  basis_prop->tKnotNum=tKnotNum - 1;
								  memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
								  memcpy(knotListNew+newIdx - 1,knotListOld+(newIdx+1) - 1,sizeof(uint16_t)*(tKnotNum - newIdx));
								  memcpy(orderListNew,orderListOld,sizeof(char)*(newIdx - 1));
								  memcpy(orderListNew+newIdx - 1,orderListOld+(newIdx+1) - 1,sizeof(char)*(tKnotNum+1 - (newIdx+1)+1));
								  k1_old=basis->tks[(newIdx)-1];
								  k2_old=basis->tke[(newIdx+1) - 1];
								  endIdx=newIdx;
								  break;
					}
					case merge:
					{
								  newIdx=(*rnd++)%(tKnotNum - 1)+1; 
								  rI16  r1=knotListOld[(newIdx)-1];
								  rI16  r2=knotListOld[(newIdx+1) - 1];
								  rI16 count=(r2 - r1)+1L - 2L;
								  if (count==0L)
								  { 
									  newKnot=*(*(int32_t **)&rnd)++> 0 ? r1 : r2;
								  }
								  else
								  {  
									  newKnot=r1+((*rnd++)%count+1);
								  }
								  basis_prop->tKnotNum=tKnotNum - 1;
								  memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
								  knotListNew[newIdx - 1]=newKnot;
								  memcpy(knotListNew+(newIdx+1) - 1,knotListOld+(newIdx+2) - 1,sizeof(uint16_t)*(tKnotNum - (newIdx+1)));
								  memcpy(orderListNew,orderListOld,sizeof(char)* newIdx);
								  memcpy(orderListNew+(newIdx+1) - 1,orderListOld+(newIdx+2) - 1,sizeof(char)*(tKnotNum+1 - (newIdx+2)+1));
								  k1_old=basis->tks[(newIdx)-1];
								  k2_old=basis->tke[(newIdx+2) - 1];
								  endIdx=newIdx+1;
								  break;
					}
					case move: 
					{
								   newIdx=(*rnd++)%tKnotNum+1; 
								   oldKnot=knotListOld[newIdx - 1];
								   rI16 r1=knotListOld[(newIdx - 1) - 1];
								   rI16 r2=knotListOld[(newIdx+1) - 1];
								   r1=max(r1+opt.minSepDist_Trend+1,oldKnot - opt.maxMoveStepSize);
								   r2=min(r2 - opt.minSepDist_Trend - 1,oldKnot+opt.maxMoveStepSize);
								   if (r2==r1)  {
									   newKnot=oldKnot;
								   }
								   else if (r2 > r1)  {
									   rI32 idx=(*rnd++)%((r2 - r1+1) - 1)+1;
									   newKnot=(r1 - 1)+idx;
									   newKnot=newKnot < oldKnot ? newKnot : (newKnot+1L);
								   }
								   else  {
									   MEM.free_all(&MEM);
									   print_error(3,&MEM);
									   return 0;
								   }
								   basis_prop->tKnotNum=tKnotNum;
								   memcpy(knotListNew,knotListOld,sizeof(uint16_t)*tKnotNum);
								   knotListNew[newIdx - 1]=newKnot;
								   memcpy(orderListNew,orderListOld,sizeof(char)*(tKnotNum+1));
								   k1_old=basis->tks[(newIdx)-1];
								   k2_old=basis->tke[(newIdx+1) - 1];
								   endIdx=newIdx+1;
								   break;
					}
					case chorder:
					{
									newIdx=(*rnd++)%(tKnotNum+1)+1;
									oldOrder=basis->tOrder[newIdx - 1];
									if (oldOrder==opt.minTrendOrder)
										newOrder=oldOrder+1;
									else if (oldOrder==opt.maxTrendOrder)
										newOrder=oldOrder - 1;
									else
										newOrder=*rndchar++>=128 ? oldOrder - 1 : oldOrder+1;
									basis_prop->tKnotNum=tKnotNum;
									memcpy(knotListNew,knotListOld,sizeof(uint16_t)*tKnotNum);
									memcpy(orderListNew,orderListOld,sizeof(char)*(tKnotNum+1));
									orderListNew[newIdx - 1]=newOrder;
									endIdx=newIdx;
									break;
					}
					}
					knotListNew[(basis_prop->tKnotNum+1L)-1L]=N+1L; 
					convert_basis_both(basis_prop);  
					newTerm_startidx=knotListNew[(newIdx - 1) - 1];
					newTerm_endidx=knotListNew[(endIdx)-1] - 1;
					if (flag !=chorder)
					{
						k1_new=k1_old;
						k2_new=k1_new+(basis_prop->tke[endIdx - 1] - basis_prop->tks[(newIdx)-1]);
						K_newTerm=k2_new - k1_new+1;
						rI32 numOfSeg=basis_prop->tKnotNum+1;
						rF32PTR   Xt=Xnewterm;
						r_ippsSet_32f(0,Xt,K_newTerm*Npad);
						k=k1_new;
						for (rI32 i=newIdx; i <=endIdx; i++)
						{
							rI16 r1=knotListNew[(i-1) -1] ;
							rI16 r2=knotListNew[ (i)  -1] -1;
							rI32 segLength=r2 - r1+1;
							rI32   ORDER=orderListNew[i - 1]+1;
							F32PTR tmpTERM=TREND_TERMS+(1L - 1L)*N+r1 - 1;
							rF32   scale=INV_SQR[(segLength)-1];
							for (rI32 j=1; j <=ORDER; j++)
							{
#if BASIS_METHODS==1
								r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
								if (numOfSeg !=1||j !=1) normalize(tmpXt,N);
#elif BASIS_METHODS==2
								if (j==1) 
								{							 
									r_ippsSet_32f(scale,Xt+(r1)-1,segLength);			 
								} 
								else if (j==2)
								{
									rF32PTR tmpFlt=Xt+r1 - 1;
									rF32         b=COEFF_B[segLength - 1];
									rF32         a=COEFF_A[segLength - 1];
									for (rI32 n=0; n < segLength; n++) {
										*tmpFlt++=a;
										a+=b;
									}
								} 
								else 
								{
									r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
									normalize_x_factor(Xt+r1 - 1,segLength,scale);						  
								}
#elif BASIS_METHODS==3
								r_cblas_scopy(segLength,tmpTERM,1L,Xt+r1 - 1,1L);
								if (j !=1)
								{
									float tmpSum;
									r_ippsSum_32f(Xt+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
									tmpSum=tmpSum/segLength;
									r_ippsSubC_32f_I(tmpSum,Xt+r1 - 1,segLength);
								}
								else
								{
								}
#endif
								k++;
								Xt+=Npad;
								tmpTERM+=N;
							} 
						} 
#if MY_DEBUG==1
						k--;
						if (k !=k2_new)
						{
							r_error("The two K's differ1; there must be something wrong!");
						}
#endif
					}
					else 
					{
						if (newOrder > oldOrder)
						{
							k2_old=basis->tke[newIdx - 1]; 
							k1_old=k2_old+1; 
							k1_new=k1_old;
							k2_new=k1_new;
							K_newTerm=k2_new - k1_new+1;
							rI32 numOfSeg=basis_prop->tKnotNum+1;
							rF32PTR Xt=Xnewterm;
							r_ippsSet_32f(0,Xt,K_newTerm*Npad);
							k=k1_new;
							for (rI32 i=newIdx; i <=newIdx; i++)
							{
								rI16 r1=knotListNew[(i - 1) - 1];
								rI16 r2=knotListNew[(i)-1] - 1;
								rI32   segLength=r2 - r1+1;
								rF32   scale=INV_SQR[(segLength)-1];
								rI32   ORDER=orderListNew[i - 1]+1;
								F32PTR tmpTERM=TREND_TERMS+(ORDER - 1)*N+r1 - 1;
								for (rI32 j=ORDER; j <=ORDER; j++)
								{
#if BASIS_METHODS==1
									if (j==1)
									{
										tmpTERM=tmpTERM+N;
										continue;
									}
									r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
									if (numOfSeg !=1||j !=1) normalize(Xt_mars_prop+(k - 1)*N,N);
#elif BASIS_METHODS==2
									if (j==1)
									{
										tmpTERM=tmpTERM+N;
										print_error(4,&MEM);
										MEM.free_all(&MEM);
										return 0;
										continue;
									}
									else if (j==2)
									{
										rF32PTR tmpFlt=Xt+r1 - 1;
										rF32         b=COEFF_B[segLength - 1];
										rF32         a=COEFF_A[segLength - 1];
										for (rI32 n=0; n < segLength; n++) {
											*tmpFlt++=a;
											a+=b;
										}
									}
									else {
										r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
										normalize_x_factor(Xt+r1 - 1,segLength,scale);
									}
#elif BASIS_METHODS==3
									if (j==1)
									{
										tmpTERM=tmpTERM+N;
										MEM.free_all(&MEM);
										r_error("j should be not at least 2. There must be something wrong! \n");
										continue;
									}
									r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
									if (j !=1)
									{
										float tmpSum;
										r_ippsSum_32f(Xt+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
										tmpSum=tmpSum/segLength;
										r_ippsSubC_32f_I(tmpSum,Xt+r1 - 1,segLength);
									}
									else
									{
										MEM.free_all(&MEM);
										r_error("j should be at least 2. There must be something wrong! \n");
									}
#endif
									k++;
									Xt=Xt+Npad;
								}
							} 
							k--;
							if (k !=k2_new)
							{
								print_error(5,&MEM);
								MEM.free_all(&MEM);
								return 0;
							}
						}
						else 
						{
							newTerm_endidx=newTerm_startidx=-999;
							k2_old=k1_old=basis->tke[newIdx - 1];
							k1_new=k1_old;
							k2_new=k1_old - 1;
							K_newTerm=k2_new - k1_new+1;
						}
					}
				}
				if (!IsTrend)
				{
					if (unifRnd <=S_propProb.birth) {
						flag=birth;
						if (sKnotNum==opt.maxKnotNum_Season||goodS_num==0)
							flag=move;
						if (K >MAX_K_STOP_ADDING_NEWseasonTERMS) {
							flag=(sKnotNum==0) ? birth : move;
						}
					}
					else if (unifRnd <=S_propProb.death) {
						flag=sKnotNum==0 ? birth : death;
					}
					else if (unifRnd <=S_propProb.merge){
						if (sKnotNum >=2) 					
							flag=merge;
						else
							flag=sKnotNum==0 ? birth : death;
					}
					else if (unifRnd <=S_propProb.move) { 
						flag=sKnotNum==0 ? birth : move;
					}
					else { 
						flag=chorder;
					}
					basis_prop->tKnotNum=tKnotNum;
					memcpy(basis_prop->T,basis->T,sizeof(uint16_t)*(tKnotNum+1));
					memcpy(basis_prop->tOrder,basis->tOrder,sizeof(char)  * (tKnotNum+1));
					rU16PTR knotListNew=basis_prop->S;
					rU08PTR  orderListNew=basis_prop->sOrder;
					rU16PTR knotListOld=basis->S;
					rU08PTR  orderListOld=basis->sOrder;
					switch (flag)
					{
					case birth: {
									{																				
										rI32 randIdx=(*rnd++)%goodS_num+1;
										newKnot=find_index_by_csum(goodS,opt.Npad16,randIdx);
									}	
									{
										rU16PTR knotList=basis->S;
										for (newIdx=1; newKnot > *knotList++; newIdx++);
					                }
									basis_prop->sKnotNum=sKnotNum+1;
									memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
									knotListNew[newIdx - 1]=newKnot;
									memcpy(knotListNew+(newIdx+1) - 1,knotListOld+newIdx - 1,sizeof(uint16_t)*(sKnotNum - (newIdx - 1)));
									memcpy(orderListNew,orderListOld,sizeof(char)*newIdx);
									memcpy(orderListNew+(newIdx+1) - 1,orderListOld+newIdx - 1,sizeof(char)*(sKnotNum+1 - newIdx+1));
									k1_old=basis->sks[(newIdx)-1];
									k2_old=basis->ske[(newIdx)-1];
									endIdx=newIdx+1;
									break;
					}
					case  death:
					{
								   newIdx=(*rnd++)%sKnotNum+1;
								   basis_prop->sKnotNum=sKnotNum - 1;
								   memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
								   memcpy(knotListNew+newIdx - 1,knotListOld+(newIdx+1) - 1,sizeof(uint16_t)*(sKnotNum - newIdx));
								   memcpy(orderListNew,orderListOld,sizeof(char)*(newIdx - 1));
								   memcpy(orderListNew+newIdx - 1,orderListOld+(newIdx+1) - 1,sizeof(char)*(sKnotNum+1 - (newIdx+1)+1));
								   k1_old=basis->sks[(newIdx)-1];
								   k2_old=basis->ske[(newIdx+1) - 1];
								   endIdx=newIdx;
								   break;
					}
					case  merge:
					{
								   newIdx=(*rnd++)%(sKnotNum - 1)+1; 
								   rI16  r1=knotListOld[(newIdx)-1];
								   rI16  r2=knotListOld[(newIdx+1) - 1];
								   rI16 count=(r2 - r1)+1L - 2L;
								   if (count==0L)
								   { 
									   newKnot=*(*(int32_t **)&rnd)++> 0 ? r1 : r2;
								   }
								   else
								   {  
									   newKnot=r1+((*rnd++)%count+1);
								   }
								   basis_prop->sKnotNum=sKnotNum - 1;
								   memcpy(knotListNew,knotListOld,sizeof(uint16_t)*(newIdx - 1));
								   knotListNew[newIdx - 1]=newKnot;
								   memcpy(knotListNew+(newIdx+1) - 1,knotListOld+(newIdx+2) - 1,sizeof(uint16_t)*(sKnotNum - (newIdx+1)));
								   memcpy(orderListNew,orderListOld,sizeof(char)* newIdx);
								   memcpy(orderListNew+(newIdx+1) - 1,orderListOld+(newIdx+2) - 1,sizeof(char)*(sKnotNum+1 - (newIdx+2)+1));
								   k1_old=basis->sks[(newIdx)-1];
								   k2_old=basis->ske[(newIdx+2) - 1];
								   endIdx=newIdx+1;
								   break;
					}
					case move:
					{
								 newIdx=(*rnd++)%sKnotNum+1;
								 oldKnot=knotListOld[newIdx - 1];
								 rI16  r1=knotListOld[(newIdx - 1) - 1];
								 rI16  r2=knotListOld[(newIdx+1) - 1];
								 r1=max(r1+opt.minSepDist_Season+1,oldKnot - opt.maxMoveStepSize);
								 r2=min(r2 - opt.minSepDist_Season - 1,oldKnot+opt.maxMoveStepSize);
								 if (r2==r1) {
									 newKnot=oldKnot;
								 }
								 else if (r2 > r1)
								 {
									 rI32 idx=(*rnd++)%((r2 - r1+1) - 1)+1;
									 newKnot=(r1 - 1)+idx;
									 newKnot=newKnot < oldKnot ? newKnot : (newKnot+1L);
								 }
								 else
								 {
									 MEM.free_all(&MEM); 
									 r_error("r2<r1(season move): there must be something wrong!\n");
									 return 0;
								 }
								 basis_prop->sKnotNum=sKnotNum;
								 memcpy(knotListNew,knotListOld,sizeof(uint16_t)*sKnotNum);
								 knotListNew[newIdx - 1]=newKnot;
								 memcpy(orderListNew,orderListOld,sizeof(char)*(sKnotNum+1));
								 k1_old=basis->sks[(newIdx)-1];
								 k2_old=basis->ske[(newIdx+1) - 1];
								 endIdx=newIdx+1;
								 break;
					}
					case chorder:
					{
									newIdx=(*rnd++)%(sKnotNum+1)+1;
									oldOrder=basis->sOrder[newIdx - 1];
									if (oldOrder==opt.minSeasonOrder)
										newOrder=oldOrder+1;
									else if (oldOrder==opt.maxSeasonOrder)
										newOrder=oldOrder - 1;
									else
										newOrder=*rndchar++>=128 ? oldOrder - 1 : oldOrder+1;
									basis_prop->sKnotNum=sKnotNum;
									memcpy(knotListNew,knotListOld,sizeof(uint16_t)*sKnotNum);
									memcpy(orderListNew,orderListOld,sizeof(char)*(sKnotNum+1));
									orderListNew[newIdx - 1]=newOrder;
									endIdx=newIdx;
									break;
					}
					}
					knotListNew[basis_prop->sKnotNum]=(N+1);
					convert_basis_both(basis_prop);
					newTerm_startidx=knotListNew[(newIdx - 1) - 1];
					newTerm_endidx=knotListNew[(endIdx)-1] - 1;
					if (flag !=chorder)
					{
						k1_new=k1_old;
						k2_new=k1_new+(basis_prop->ske[endIdx - 1] - basis_prop->sks[(newIdx)-1]);
						K_newTerm=k2_new - k1_new+1;
						k=k1_new;
						rF32PTR Xt=Xnewterm;
						r_ippsSet_32f(0,Xt,K_newTerm*Npad);
						for (rI32 i=newIdx; i <=endIdx; i++)
						{
							rI16 r1=knotListNew[(i - 1) - 1] ;
							rI16 r2=knotListNew[(i)-1]-1;
							rI32 segLength=r2 - r1+1;
							rF32 scale=INV_SQR[(segLength)-1];
							rI32   ORDER=2 * orderListNew[i - 1];
							rF32PTR tmpTERM=SEASON_TERMS+(1 - 1)*N+r1 - 1;
							rF32PTR season_csum=SEASON_SQR_CSUM+1L;
							for (rI32 j=1; j <=ORDER; j++)
							{
#if BASIS_METHODS==1
								r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
								normalize(tmpXt,N);
#elif BASIS_METHODS==2
								r_cblas_scopy(segLength,tmpTERM,1L,Xt+r1 - 1,1L);
								rF32 scalingFactor=sqrtf( N/( season_csum[r2-1] - season_csum[(r1-1)-1]) );
								r_cblas_sscal(segLength,scalingFactor,Xt+r1 - 1,1L);
#elif BASIS_METHODS==3
								r_cblas_scopy(segLength,tmpTERM,1,Xt+r1 - 1,1);
								float tmpSum;
								r_ippsSum_32f(Xt+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
								tmpSum=tmpSum/segLength;
								r_ippsSubC_32f_I(tmpSum,Xt+r1 - 1,segLength);
#endif
								k++;
								tmpTERM+=N;
								Xt+=Npad;
								season_csum+=(N+1);
							}
						} 
#if MY_DEBUG==1
						k--;
						if (k !=k2_new)
						{
							r_error("The two K's differ5; there must be something wrong!");
						}
#endif
					} 
					else 
					{
						if (newOrder > oldOrder)
						{
							k2_old=basis->ske[newIdx - 1];
							k1_old=k2_old+1;
							k1_new=k1_old;
							k2_new=k1_new+1;
							K_newTerm=k2_new - k1_new+1;
							k=k1_new;
							rF32PTR Xt=Xnewterm;
							r_ippsSet_32f(0,Xt,K_newTerm*Npad);
							for (rI32 i=newIdx; i <=newIdx; i++)
							{
								rI16 r1=knotListNew[(i - 1) - 1];
								rI16 r2=knotListNew[(i)-1] - 1;
								rI32 segLength=r2 - r1+1;
								rF32 scale=INV_SQR[(segLength)-1];
								rI32 ORDER=orderListNew[i - 1] * 2;
								for (rI32 j=1; j <=ORDER; j++)
								{
#if BASIS_METHODS==1
									if (j <=(ORDER - 2)) continue;
									r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1;,1,Xt+r1 - 1,1);
									normalize(tmpXt,N);
#elif BASIS_METHODS==2
									if (j <=(ORDER - 2)) continue;
									r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1L,Xt+r1 - 1,1L);
									rF32PTR season_csum=SEASON_SQR_CSUM+(j - 1)*(N+1);
									season_csum++;
									rF32 scalingFactor=season_csum[r2 - 1] - season_csum[(r1 - 1) - 1];
									scalingFactor=sqrtf(N/scalingFactor);
									r_cblas_sscal(segLength,scalingFactor,Xt+r1 - 1,1L);
#elif BASIS_METHODS==3
									if (j <=(ORDER - 2)) continue;
									r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt+r1 - 1,1);
									float tmpSum;
									r_ippsSum_32f(Xt+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
									tmpSum=tmpSum/segLength;
									r_ippsSubC_32f_I(tmpSum,Xt+r1 - 1,segLength);
#endif
									k++;
									Xt=Xt+Npad;
								}
							} 
							k--;
							if (k !=k2_new)
							{
								MEM.free_all(&MEM);
								r_error("The two K's differ8; there must be something wrong!");
								return 0;
							}
						}
						else
						{
							newTerm_endidx=newTerm_startidx=-999;
							k2_old=basis->ske[newIdx - 1];
							k1_old=k2_old - 1;
							k1_new=k1_old;
							k2_new=k1_old - 1;
							K_newTerm=k2_new - k1_new+1;
						}
					}
				}
				if (yInfo.nMissing > 0 && K_newTerm !=0) {
					p_zeroOut_Xmars_zero(Xnewterm,Xt_zeroBackup,yInfo.rowsMissing,yInfo.nMissing,N,Npad,K_newTerm);
				}
				int Kold=K;             
				int Knew=basis_prop->k; 
				GlobalMEMBuf_1st=basis_prop->post_P_U;
				GlobalMEMBuf_2nd=Xnewterm+K_newTerm*Npad;
				{
					for (rI32 i=1; i <=(k1_new - 1); i++)
						r_cblas_scopy(i,XtX+(i - 1)*Kold,1L,XtX_prop+(i - 1)*Knew,1L);
				}
				if (K_newTerm !=0)
				{
					rI32 Nseg=newTerm_endidx - newTerm_startidx+1;
					if (k1_new !=1)
					{
						r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k1_new - 1,K_newTerm,Nseg,1.0f,Xt_mars+(newTerm_startidx)-1,Npad, 
							Xnewterm+(newTerm_startidx)-1,Npad,0.f,GlobalMEMBuf_1st,k1_new - 1);
					}
					{
						r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,K_newTerm,K_newTerm,Nseg,1.0,Xnewterm+newTerm_startidx - 1, 
							Npad,Xnewterm+newTerm_startidx - 1,Npad,0.f,GlobalMEMBuf_2nd,K_newTerm);
					}
					for (rI32 i=k1_new,j=1; i <=k2_new; i++,j++)
					{
						if (k1_new !=1)
							r_cblas_scopy(k1_new - 1,GlobalMEMBuf_1st+(j - 1)*(k1_new - 1),1,XtX_prop+(i - 1)*Knew,1);
						r_cblas_scopy(j,GlobalMEMBuf_2nd+(j - 1)*K_newTerm,1,XtX_prop+(i - 1)*Knew+k1_new - 1,1);
					}
				}
				if (k2_old !=Kold)
				{
					rI32 Nseg=newTerm_endidx - newTerm_startidx+1;
					if (K_newTerm !=0)
					{
						r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,K_newTerm,(Kold - k2_old),Nseg,1,Xnewterm+newTerm_startidx-1,
							Npad,Xt_mars+(k2_old+1 - 1)*Npad+newTerm_startidx - 1,Npad,0,GlobalMEMBuf_1st,K_newTerm);
						k=k2_new+1;
						for (rI32 i=k2_old+1,j=1; i <=Kold; i++,j++)
							r_cblas_scopy(k1_old - 1,XtX+(i - 1)*Kold,1,XtX_prop+(k - 1)*Knew,1),
							r_cblas_scopy(K_newTerm,GlobalMEMBuf_1st+(j - 1)*K_newTerm,1,XtX_prop+(k - 1)*Knew+k1_new - 1,1),
							r_cblas_scopy(i - k2_old,XtX+(i - 1)*Kold+(k2_old+1) - 1,1,XtX_prop+(k - 1)*Knew+(k2_new+1) - 1,1),
							k++;
					}
					else
					{
						k=k2_new+1;
						for (rI32 i=k2_old+1,j=1; i <=Kold; i++,j++)
							r_cblas_scopy(k1_old - 1,XtX+(i - 1)*Kold,1,XtX_prop+(k - 1)*Knew,1),
							r_cblas_scopy(i - k2_old,XtX+(i - 1)*Kold+(k2_old+1) - 1,1,XtX_prop+(k - 1)*Knew+(k2_new+1) - 1,1),
							k++;
					}
				}
				if (k1_old !=1)
				{
					r_cblas_scopy(k1_old - 1,XtY,1,XtY_prop,1);
				}
				if (K_newTerm !=0)
				{
					rI32 Nseg=newTerm_endidx - newTerm_startidx+1;
					r_cblas_sgemv(CblasColMajor,CblasTrans,Nseg,K_newTerm,1.f,Xnewterm+newTerm_startidx-1,
						Npad,yInfo.Y+newTerm_startidx-1,1L,0.f,GlobalMEMBuf_1st,1L);
					r_cblas_scopy(K_newTerm,GlobalMEMBuf_1st,1,XtY_prop+k1_new - 1,1);
				}
				if (k2_old !=Kold)
				{
					r_cblas_scopy(Knew - k2_new,XtY+(k2_old+1) - 1,1L,XtY_prop+(k2_new+1) - 1,1L);
				}
				{
					cp(Knew*Knew,XtX_prop,post_P_U_prop);
					rU08PTR termType=basis_prop->termType;
					for (rI32 i=1,j=0; i <=Knew; i++)
					{
						post_P_U_prop[j+(i)-1]+=modelPar.prec[*termType++];
						j+=Knew;
					}
				}
				r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',Knew,post_P_U_prop,Knew); 
				r_cblas_scopy(Knew,XtY_prop,1,beta_mean_prop,1);
				r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',Knew,1,post_P_U_prop,Knew,beta_mean_prop,Knew);
				basis_prop->alpha_star=yInfo.YtY - DOT(Knew,XtY_prop,beta_mean_prop);
				{
					rF32  half_log_det_prior;
					rF32  half_log_det_post;
					half_log_det_post=-sum_log_diag(post_P_U_prop,Knew);
					half_log_det_prior=-0.5f*(basis_prop->k_SN * fastlog(modelPar.prec[0])+basis_prop->k_const*fastlog(modelPar.prec[1])+(Knew - basis_prop->k_SN - basis_prop->k_const)*fastlog(modelPar.prec[2]));
					basis_prop->marg_lik=half_log_det_post - half_log_det_prior - (yInfo.n * 0.5f+modelPar.alpha_2) *fastlog(modelPar.alpha_1+basis_prop->alpha_star  * 0.5f);
				}
				float scale=1;
				{					
					if (flag==birth)
						scale=IsTrend ? scaleFactorTrend[basis->tKnotNum] : scaleFactorSeason[basis->sKnotNum];
					else if (flag==death)
						scale=IsTrend ? 1.0f/scaleFactorTrend[basis->tKnotNum - 1] : 1.0f/scaleFactorSeason[basis->sKnotNum - 1];
				}
				float factor;
				int delta_k1=IsTrend ? basis->tKnotNum : basis->sKnotNum;
				delta_k1++;
				int delta_k2=IsTrend ? basis_prop->tKnotNum : basis_prop->sKnotNum;
				delta_k2++;
				if (delta_k1==delta_k2 && basis->k==basis_prop->k)
					factor=0;
				else
				{
					int k=min(delta_k1,delta_k2);
					int j=IsTrend ? max(basis->k - basis->k_SN,basis_prop->k - basis_prop->k_SN) : max(basis->k_SN,basis_prop->k_SN)/2;
					int j0=IsTrend ? min(basis->k - basis->k_SN,basis_prop->k - basis_prop->k_SN) : min(basis->k_SN,basis_prop->k_SN)/2;
					factor=1;
					for (rI32 i=1; i <=k - 1; i++)
						factor=factor*(j - 1 - i+1)/(j0 - 1 - i+1);
					factor=factor * (j - 1 - (k - 1))/k;
					factor=fastlog(factor);
					factor=delta_k1 < delta_k2 ? -factor : factor;
				}
				float delta_lik=basis_prop->marg_lik - basis->marg_lik+factor;
				if (delta_lik > 0||*rnd++<fastexp(delta_lik) * scale * 4.294967296000000e+09)
				{
					if (IsTrend)
						++(accT[flag]) ;
					else
						++(accS[flag]);
					if (yInfo.nMissing >0 && K_newTerm !=0)
					{
						zeroOut_Xmars_fill(Xnewterm,Xt_zeroBackup,yInfo.rowsMissing,yInfo.nMissing,N,Npad,K_newTerm);
					}
					if (Kold !=k2_old  && k2_new !=k2_old)
					{
						rI32 j=k2_new - k2_old;
						if (j < 0||(k2_new+1) > Kold)
							r_cblas_scopy((Kold - k2_old)*Npad,Xt_mars+(k2_old+1 - 1)*Npad,1,Xt_mars+(k2_new+1 - 1)*Npad,1);
						else
						{
							rI32 i=Kold+1;
							while (true)
							{
								i=i - j;
								if (i > (k2_old+1))
								{
									r_cblas_scopy(j*Npad,Xt_mars+(i - 1)*Npad,1,Xt_mars+((i+j) - 1)*Npad,1);
								}
								else
								{
									j=(i+j) - (k2_old+1);
									r_cblas_scopy(j*Npad,Xt_mars+((k2_old+1) - 1)*Npad,1,Xt_mars+((k2_new+1) - 1)*Npad,1);
									break;
								}
							}
						}
					}
					if (K_newTerm !=0)
						r_cblas_scopy(K_newTerm*Npad,Xnewterm,1,Xt_mars+(k1_new - 1)*Npad,1);
					if (IsTrend)
					{
						if (flag==birth)
						{
							r_ippsSet_8u(0,goodT+(newKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
						}
						else if (flag==death)
						{
							oldKnot=basis->T[newIdx - 1];
							rI16 r1=basis->T[(newIdx - 1) - 1] ;
							rI16 r2=basis->T[(newIdx+1) - 1] -1;
							r_ippsSet_8u(1,goodT+(oldKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r1)-1,opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r2-opt.minSepDist_Trend+1) - 1,opt.minSepDist_Trend);
						}
						else if (flag==move)
						{
							rI16 r1=basis->T[(newIdx - 1) - 1] ;
							rI16 r2=basis->T[(newIdx+1) - 1]-1;
							r_ippsSet_8u(1,goodT+(oldKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(newKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r1)-1,opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r2 - opt.minSepDist_Trend+1) - 1,opt.minSepDist_Trend );
						}
						else if (flag==merge)
						{
							rI16 r1=basis->T[(newIdx - 1) - 1] ;
							rI16 r2=basis->T[(newIdx+2) - 1]-1;
							oldKnot=basis->T[newIdx - 1];
							r_ippsSet_8u(1,goodT+(oldKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							oldKnot=basis->T[(newIdx+1) - 1];
							r_ippsSet_8u(1,goodT+(oldKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(newKnot - opt.minSepDist_Trend) - 1,2 * opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r1)-1,opt.minSepDist_Trend+1);
							r_ippsSet_8u(0,goodT+(r2 - opt.minSepDist_Trend+1) - 1,opt.minSepDist_Trend );
						}
						if (flag !=chorder)
						{
							goodT_num=int8_arr_sum(goodT,opt.Npad16);
						}
					}
					else 
					{
						if (flag==birth)
						{
							r_ippsSet_8u(0,goodS+(newKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
						}
						else if (flag==death)
						{
							oldKnot=basis->S[newIdx - 1];
							rI16 r1=basis->S[(newIdx - 1) - 1] ;
							rI16 r2=basis->S[(newIdx+1) - 1] -1;
							r_ippsSet_8u(1,goodS+(oldKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r1)-1,opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r2 - opt.minSepDist_Season+1) - 1,opt.minSepDist_Season );
						}
						else if (flag==move)
						{
							rI16 r1=basis->S[(newIdx - 1) - 1] ;
							rI16 r2=basis->S[(newIdx+1) - 1]-1;
							r_ippsSet_8u(1,goodS+(oldKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(newKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r1)-1,opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r2 - opt.minSepDist_Season+1) - 1,opt.minSepDist_Season );
						}
						else if (flag==merge)
						{
							rI16 r1=basis->S[(newIdx - 1) - 1] ;
							rI16 r2=basis->S[(newIdx+2) - 1] -1;
							oldKnot=basis->S[newIdx - 1];
							r_ippsSet_8u(1,goodS+(oldKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							oldKnot=basis->S[(newIdx+1) - 1];
							r_ippsSet_8u(1,goodS+(oldKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(newKnot - opt.minSepDist_Season) - 1,2 * opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r1)-1,opt.minSepDist_Season+1);
							r_ippsSet_8u(0,goodS+(r2 - opt.minSepDist_Season+1) - 1,opt.minSepDist_Season );
						}
						if (flag !=chorder)
						{
							goodS_num=int8_arr_sum(goodS,opt.Npad16);
						}
					}
					{
						register struct BASIS * tmp=basis;
						basis=basis_prop;
						basis_prop=tmp;
						XtX_prop=basis_prop->XtX;
						XtY_prop=basis_prop->XtY;
						beta_mean_prop=basis_prop->beta_mean;
						beta_prop=basis_prop->beta;
						post_P_U_prop=basis_prop->post_P_U;
						XtX=basis->XtX;
						XtY=basis->XtY;
						beta_mean=basis->beta_mean;
						beta=basis->beta;
						post_P_U=basis->post_P_U;
						K=basis->k;
						K_SN=basis->k_SN;
						tKnotNum=basis->tKnotNum;
						sKnotNum=basis->sKnotNum;
					}
				}
				if (ite < opt.burnin) continue;
				if (ite%20==0||(ite%opt.thinningFactor==0 && ite>opt.burnin))
				{
					modelPar.sig2=(*rndgamma++)*1.f/(modelPar.alpha_1+basis->alpha_star *0.5f);
					modelPar.sig2=1.f/modelPar.sig2;
					r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,K,beta,0,1);
					r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',K,1,post_P_U,K,beta,K);
					r_ippsMulC_32f_I(fast_sqrt(modelPar.sig2),beta,K);
					r_ippsAdd_32f_I(beta_mean,beta,K);
				}
				if (ite%20==0)
				{
					GlobalMEMBuf_1st=Xnewterm;
					rF32 sumq=DOT(K,beta,beta);
					r_vsRngGamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE,stream,1,modelPar.prec,(modelPar.del_1+K *0.5f),0,1.f);
					modelPar.prec[2]=modelPar.prec[1]=modelPar.prec[0]=(*modelPar.prec)/(modelPar.del_2+0.5f*sumq/modelPar.sig2);
					r_cblas_scopy(K*K,XtX,1L,post_P_U,1L);
					{
						rU08PTR termType=basis->termType;
						for (rI32 i=1,j=0; i <=K; i++)
						{
							post_P_U[j+(i)-1]+=modelPar.prec[*termType++];
							j+=K;
						}
					}
					r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',K,post_P_U,K); 
					r_cblas_scopy(K,XtY,1,beta_mean,1);
					r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',K,1,post_P_U,K,beta_mean,K);
					basis->alpha_star=yInfo.YtY - DOT(K,beta_mean,XtY);
					rF32 half_log_det_post ;
					half_log_det_post=-sum_log_diag(post_P_U,K);
					rF32 half_log_det_prior=-.5f*(K_SN * fastlog(modelPar.prec[0])+basis->k_const*fastlog(modelPar.prec[1])+(K - K_SN - basis->k_const)*fastlog(modelPar.prec[2]));
					basis->marg_lik=half_log_det_post - half_log_det_prior - (yInfo.n *0.5f+modelPar.alpha_2) *fastlog(modelPar.alpha_1+basis->alpha_star * 0.5f);
				}
				if (ite <=opt.burnin||ite%opt.thinningFactor !=0) continue;
				sample++;
				if (opt.printToScreen && sample%5000==0)
				{
					register char * _restrict string;
					string=(char*)Xnewterm;
					r_printf("[Ite:%06d---Sample:%06d----MargLik:%7.3E]\n",ite,sample,basis->marg_lik);
					r_ippsSet_8u(32,string,opt.printCharLen+1);
					rI32 numOfSeg=sKnotNum+1;
					rI16  r2=-9999;
					rI16  r1;
					for (rI32 i=1; i <=numOfSeg; i++)
					{
						rI32 r1=(int)((basis->S[(i - 1) - 1]+1.f)/(N+0.f) *opt.printCharLen);
						r1=max(r1,r2+1);
						r2=(int)((basis->S[i - 1])/(N+0.f) *opt.printCharLen);
						r1=max(r1,1);
						r2=max(r2,1);
						r_ippsSet_8u(48+basis->sOrder[i - 1],string+r1 - 1,max(r2 - r1+1,0));
						if (i !=numOfSeg) 	string[r2 - 1]='*';
					}
					string[opt.printCharLen+1 - 1]=0;
					r_printf("S[#bp%02d]:%s \n",sKnotNum,string);
					r_ippsSet_8u(32,string,opt.printCharLen+1);
					numOfSeg=tKnotNum+1;
					r2=-9999;
					for (rI32 i=1; i <=numOfSeg; i++)
					{
						r1=(i==1) ? 1 : (int)((basis->T[(i - 1) - 1]+1.f)/(N+0.f) *opt.printCharLen);
						r1=max(r1,r2+1);
						r2=(i==numOfSeg) ? opt.printCharLen : (int)((basis->T[i - 1])/(N+0.f) *opt.printCharLen);
						r1=max(r1,1); 
						r2=max(r2,1);
						r_ippsSet_8u(48+basis->tOrder[i - 1],string+r1 - 1,max(r2 - r1+1,0));
						if (i !=numOfSeg) 	string[r2 - 1]='*';
					}
					string[opt.printCharLen+1 - 1]=0;
					r_printf("T[#bp%02d]:%s \n\n",tKnotNum,string);
				}
				*resultChain.sig2=*resultChain.sig2+modelPar.sig2;
				GlobalMEMBuf_1st=Xnewterm;
				GlobalMEMBuf_2nd=Xnewterm+Npad;
				{
					resultChain.sNProb[sKnotNum]+=1;
					resultChain.tNProb[tKnotNum]+=1;
					*resultChain.marg_lik+=basis->marg_lik;
					rU16PTR knotList=basis->T;
					for (rI32 i=1; i <=tKnotNum; i++)
						resultChain.tProb[(*knotList++) - 1]+=1;
					knotList=basis->S;
					for (rI32 i=1; i <=sKnotNum; i++)
						resultChain.sProb[(*knotList++) - 1]+=1;
					if (opt.computeHarmonicOrder)
					{
						uint8_t * _restrict sOrderList=basis->sOrder;
						knotList=basis->S;
						for (rI32 i=0; i <=sKnotNum; i++)
						{
							rI16 r1=knotList[i - 1] ;
							rI16 r2=knotList[i]-1;
							r_ippsAddC_32s_ISfs(*sOrderList++,resultChain.horder+r1 - 1,r2 - r1+1,0);
						}
					}
					if (opt.computeTrendOrder)
					{
						uint8_t * _restrict  tOrderList=basis->tOrder;
						knotList=basis->T;
						for (rI32 i=0; i <=tKnotNum; i++)
						{
							rI16 r1=knotList[i - 1];
							rI16 r2=knotList[i] -1L;
							r_ippsAddC_32s_ISfs(*tOrderList++,resultChain.torder+r1 - 1,r2 - r1+1,0);
						}
					}
					r_cblas_sgemv(CblasColMajor,CblasNoTrans,Npad,K_SN,1,Xt_mars,Npad,beta_mean,1L,0.f,GlobalMEMBuf_1st,1L);
					r_ippsAdd_32f_I(GlobalMEMBuf_1st,resultChain.s,N); 
					r_ippsMul_32f_I(GlobalMEMBuf_1st,GlobalMEMBuf_1st,N);
					r_ippsAdd_32f_I(GlobalMEMBuf_1st,resultChain.sSD,N); 
					r_cblas_sgemv(CblasColMajor,CblasNoTrans,Npad,(K - K_SN),1,Xt_mars+(K_SN+1 - 1)*Npad,Npad,beta_mean+(K_SN+1) - 1,1L,0.f,GlobalMEMBuf_1st,1L);
					r_cblas_scopy(N,GlobalMEMBuf_1st,1,GlobalMEMBuf_2nd,1);
					r_ippsAdd_32f_I(GlobalMEMBuf_1st,resultChain.t,N);
					r_ippsMul_32f_I(GlobalMEMBuf_1st,GlobalMEMBuf_1st,N);
					r_ippsAdd_32f_I(GlobalMEMBuf_1st,resultChain.tSD,N); 
					ONE_STEP_DIFF(GlobalMEMBuf_2nd,GlobalMEMBuf_2nd+1,GlobalMEMBuf_2nd,N - 1); 
					*GlobalMEMBuf_1st=*GlobalMEMBuf_2nd; 
					r_cblas_scopy(N - 1,GlobalMEMBuf_2nd,1,GlobalMEMBuf_1st+1,1);
					*resultChain.b+=*GlobalMEMBuf_2nd;
					r_ippsAdd_32f_I(GlobalMEMBuf_2nd,resultChain.b+1,N - 1);
					r_ippsMul_32f_I(GlobalMEMBuf_2nd,GlobalMEMBuf_2nd,N - 1);
					*resultChain.bSD+=*GlobalMEMBuf_2nd;
					r_ippsAdd_32f_I(GlobalMEMBuf_2nd,resultChain.bSD+1,N - 1); 
					if (opt.computeSlopeSign)
					{
						for (int i=0; i<N; i++)
						{
							if (*GlobalMEMBuf_1st++>0)++*resultChain.bsign;
							resultChain.bsign++;
						}
						GlobalMEMBuf_1st  -=N;
						resultChain.bsign -=N;
					}
					if (!opt.computeCredible)
						continue;
					if (opt.fastCIComputation)
					{
						if (*rnd++< subsampleFraction_x_INT32MAX)
							++subSampleIndex;
						else
							continue;
					}
					else
					{
						subSampleIndex=sample;
					}
					F32PTR tmpCredUpper;
					int		* _restrict tmpMinIdx;
					F32PTR tmpMinAll;
					F32PTR tmpMinStrip;
					int		* _restrict whichStripMin;
					F32PTR tmpCredLower;
					int		* _restrict tmpMaxIdx;
					F32PTR tmpMaxAll;
					F32PTR tmpMaxStrip;
					int		* _restrict whichStripMax;
					r_cblas_sgemv(CblasColMajor,CblasNoTrans,Npad,K_SN,1,Xt_mars,Npad,beta,1,0.f,GlobalMEMBuf_1st+Npad,1);
					r_cblas_sgemv(CblasColMajor,CblasNoTrans,Npad,(K - K_SN),1,Xt_mars+(K_SN+1 - 1)*Npad,Npad,beta+K_SN+1 - 1,1,0.f,GlobalMEMBuf_1st+Npad+Npad,1);
					if (subSampleIndex <=numCISample)
					{
						int64_t offset=(subSampleIndex - 1)*N;
						r_cblas_scopy(N,GlobalMEMBuf_1st,1,credB_upper+offset,1);
						r_cblas_scopy(N,GlobalMEMBuf_1st+Npad,1,credS_upper+offset,1);
						r_cblas_scopy(N,GlobalMEMBuf_1st+Npad+Npad,1,credT_upper+offset,1);
						if (subSampleIndex==numCISample)
						{
							r_mkl_simatcopy('C','T',N,numCISample,1,credT_upper,N,numCISample);
							r_cblas_scopy(N*numCISample,credT_upper,1,credT_lower,1);
							tmpCredUpper=credT_upper;
							tmpMinIdx=min_credT_upper_IdxInStrip;
							tmpMinAll=resultChain.tCI;
							tmpMinStrip=min_credT_upperStrip;
							whichStripMin=whichMin_credT_upper;
							tmpCredLower=credT_lower;
							tmpMaxIdx=max_credT_lower_IdxInStrip;
							tmpMaxAll=resultChain.tCI+N;
							tmpMaxStrip=max_credT_lowerStrip;
							whichStripMax=whichMax_credT_lower;
							for (rI32 i=1; i <=N; i++)
							{ 
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMinIndx_32f(tmpCredUpper,rowsPerStrip[j - 1],tmpMinStrip,tmpMinIdx);
									tmpMinStrip++;
									tmpMinIdx++;
									tmpCredUpper+=rowsPerStrip[j - 1];
								}
								r_ippsMinIndx_32f(tmpMinStrip - numStrips,numStrips,tmpMinAll,whichStripMin);
								tmpMinAll++;
								whichStripMin++;
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMaxIndx_32f(tmpCredLower,rowsPerStrip[j - 1],tmpMaxStrip,tmpMaxIdx);
									tmpMaxStrip++;
									tmpMaxIdx++;
									tmpCredLower+=rowsPerStrip[j - 1];
								}
								r_ippsMaxIndx_32f(tmpMaxStrip - numStrips,numStrips,tmpMaxAll,whichStripMax);
								tmpMaxAll++;
								whichStripMax++;
							}
							r_mkl_simatcopy('C','T',N,numCISample,1,credS_upper,N,numCISample);
							r_cblas_scopy(N*numCISample,credS_upper,1,credS_lower,1);
							tmpCredUpper=credS_upper;
							tmpMinIdx=min_credS_upper_IdxInStrip;
							tmpMinAll=resultChain.sCI;
							tmpMinStrip=min_credS_upperStrip;
							whichStripMin=whichMin_credS_upper;
							tmpCredLower=credS_lower;
							tmpMaxIdx=max_credS_lower_IdxInStrip;
							tmpMaxAll=resultChain.sCI+N;
							tmpMaxStrip=max_credS_lowerStrip;
							whichStripMax=whichMax_credS_lower;
							for (rI32 i=1; i <=N; i++)
							{ 
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMinIndx_32f(tmpCredUpper,rowsPerStrip[j - 1],tmpMinStrip,tmpMinIdx);
									tmpMinStrip++;
									tmpMinIdx++;
									tmpCredUpper+=rowsPerStrip[j - 1];
								}
								r_ippsMinIndx_32f(tmpMinStrip - numStrips,numStrips,tmpMinAll,whichStripMin);
								tmpMinAll++;
								whichStripMin++;
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMaxIndx_32f(tmpCredLower,rowsPerStrip[j - 1],tmpMaxStrip,tmpMaxIdx);
									tmpMaxStrip++;
									tmpMaxIdx++;
									tmpCredLower+=rowsPerStrip[j - 1];
								}
								r_ippsMaxIndx_32f(tmpMaxStrip - numStrips,numStrips,tmpMaxAll,whichStripMax);
								tmpMaxAll++;
								whichStripMax++;
							}
							r_mkl_simatcopy('C','T',N,numCISample,1,credB_upper,N,numCISample);
							r_cblas_scopy(N*numCISample,credB_upper,1,credB_lower,1);
							tmpCredUpper=credB_upper;
							tmpMinIdx=min_credB_upper_IdxInStrip;
							tmpMinAll=resultChain.bCI;
							tmpMinStrip=min_credB_upperStrip;
							whichStripMin=whichMin_credB_upper;
							tmpCredLower=credB_lower;
							tmpMaxIdx=max_credB_lower_IdxInStrip;
							tmpMaxAll=resultChain.bCI+N;
							tmpMaxStrip=max_credB_lowerStrip;
							whichStripMax=whichMax_credB_lower;
							for (rI32 i=1; i <=N; i++)
							{ 
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMinIndx_32f(tmpCredUpper,rowsPerStrip[j - 1],tmpMinStrip,tmpMinIdx);
									tmpMinStrip++;
									tmpMinIdx++;
									tmpCredUpper+=rowsPerStrip[j - 1];
								}
								r_ippsMinIndx_32f(tmpMinStrip - numStrips,numStrips,tmpMinAll,whichStripMin);
								tmpMinAll++;
								whichStripMin++;
								for (rU32 j=1; j <=numStrips; j++)
								{
									r_ippsMaxIndx_32f(tmpCredLower,rowsPerStrip[j - 1],tmpMaxStrip,tmpMaxIdx);
									tmpMaxStrip++;
									tmpMaxIdx++;
									tmpCredLower+=rowsPerStrip[j - 1];
								}
								r_ippsMaxIndx_32f(tmpMaxStrip - numStrips,numStrips,tmpMaxAll,whichStripMax);
								tmpMaxAll++;
								whichStripMax++;
							}
						}
					}
					else
					{
						tmpCredUpper=credB_upper;
						tmpMinIdx=min_credB_upper_IdxInStrip;
						tmpMinAll=resultChain.bCI;
						tmpMinStrip=min_credB_upperStrip;
						whichStripMin=whichMin_credB_upper;
						tmpCredLower=credB_lower;
						tmpMaxIdx=max_credB_lower_IdxInStrip;
						tmpMaxAll=resultChain.bCI+N;
						tmpMaxStrip=max_credB_lowerStrip;
						whichStripMax=whichMax_credB_lower;
						float * tmpY=GlobalMEMBuf_1st;
						int64_t which;
						float * tmp;
						for (rI32 i=1; i <=N; i++)
						{
							rF32 CurrentYValue=*tmpY++;
							if (CurrentYValue > *tmpMinAll)
							{
								which=(int64_t)*whichStripMin;
								tmp=tmpCredUpper+startIdxOfStrip[which];
								*(tmp+tmpMinIdx[which])=CurrentYValue;
								r_ippsMinIndx_32f(tmp,rowsPerStrip[which],tmpMinStrip+which,tmpMinIdx+which);
								r_ippsMinIndx_32f(tmpMinStrip,numStrips,tmpMinAll,whichStripMin);
							}
							tmpCredUpper+=numCISample;
							tmpMinStrip+=numStrips;
							tmpMinIdx+=numStrips;
							tmpMinAll++;
							whichStripMin++;
							if (CurrentYValue < *tmpMaxAll)
							{
								which=(int64_t)*whichStripMax;
								tmp=tmpCredLower+startIdxOfStrip[which];
								*(tmp+tmpMaxIdx[which])=CurrentYValue;
								r_ippsMaxIndx_32f(tmp,rowsPerStrip[which],tmpMaxStrip+which,tmpMaxIdx+which);
								r_ippsMaxIndx_32f(tmpMaxStrip,numStrips,tmpMaxAll,whichStripMax);
							}
							tmpCredLower+=numCISample;
							tmpMaxStrip+=numStrips;
							tmpMaxIdx+=numStrips;
							tmpMaxAll++;
							whichStripMax++;
						}
						tmpCredUpper=credS_upper;
						tmpMinIdx=min_credS_upper_IdxInStrip;
						tmpMinAll=resultChain.sCI;
						tmpMinStrip=min_credS_upperStrip;
						whichStripMin=whichMin_credS_upper;
						tmpCredLower=credS_lower;
						tmpMaxIdx=max_credS_lower_IdxInStrip;
						tmpMaxAll=resultChain.sCI+N;
						tmpMaxStrip=max_credS_lowerStrip;
						whichStripMax=whichMax_credS_lower;
						tmpY=GlobalMEMBuf_1st+Npad;
						for (rI32 i=1; i <=N; i++)
						{
							rF32 CurrentYValue=*tmpY++;
							if (CurrentYValue > *tmpMinAll)
							{
								which=(int64_t)*whichStripMin;
								tmp=tmpCredUpper+startIdxOfStrip[which];
								*(tmp+tmpMinIdx[which])=CurrentYValue;
								r_ippsMinIndx_32f(tmp,rowsPerStrip[which],tmpMinStrip+which,tmpMinIdx+which);
								r_ippsMinIndx_32f(tmpMinStrip,numStrips,tmpMinAll,whichStripMin);
							}
							tmpCredUpper+=numCISample;
							tmpMinStrip+=numStrips;
							tmpMinIdx+=numStrips;
							tmpMinAll++;
							whichStripMin++;
							if (CurrentYValue < *tmpMaxAll)
							{
								which=(int64_t)*whichStripMax;
								tmp=tmpCredLower+startIdxOfStrip[which];
								*(tmp+tmpMaxIdx[which])=CurrentYValue;
								r_ippsMaxIndx_32f(tmp,rowsPerStrip[which],tmpMaxStrip+which,tmpMaxIdx+which);
								r_ippsMaxIndx_32f(tmpMaxStrip,numStrips,tmpMaxAll,whichStripMax);
							}
							tmpCredLower+=numCISample;
							tmpMaxStrip+=numStrips;
							tmpMaxIdx+=numStrips;
							tmpMaxAll++;
							whichStripMax++;
						}
						tmpCredUpper=credT_upper;
						tmpMinIdx=min_credT_upper_IdxInStrip;
						tmpMinAll=resultChain.tCI;
						tmpMinStrip=min_credT_upperStrip;
						whichStripMin=whichMin_credT_upper;
						tmpCredLower=credT_lower;
						tmpMaxIdx=max_credT_lower_IdxInStrip;
						tmpMaxAll=resultChain.tCI+N;
						tmpMaxStrip=max_credT_lowerStrip;
						whichStripMax=whichMax_credT_lower;
						tmpY=GlobalMEMBuf_1st+Npad+Npad;
						for (rI32 i=1; i <=N; i++)
						{
							rF32 CurrentYValue=*tmpY++;
							if (CurrentYValue > *tmpMinAll)
							{
								which=(int64_t)*whichStripMin;
								tmp=tmpCredUpper+startIdxOfStrip[which];
								*(tmp+tmpMinIdx[which])=CurrentYValue;
								r_ippsMinIndx_32f(tmp,rowsPerStrip[which],tmpMinStrip+which,tmpMinIdx+which);
								r_ippsMinIndx_32f(tmpMinStrip,numStrips,tmpMinAll,whichStripMin);
							}
							tmpCredUpper+=numCISample;
							tmpMinStrip+=numStrips;
							tmpMinIdx+=numStrips;
							tmpMinAll++;
							whichStripMin++;
							if (CurrentYValue < *tmpMaxAll)
							{
								which=(int64_t)*whichStripMax;
								tmp=tmpCredLower+startIdxOfStrip[which];
								*(tmp+tmpMaxIdx[which])=CurrentYValue;
								r_ippsMaxIndx_32f(tmp,rowsPerStrip[which],tmpMaxStrip+which,tmpMaxIdx+which);
								r_ippsMaxIndx_32f(tmpMaxStrip,numStrips,tmpMaxAll,whichStripMax);
							}
							tmpCredLower+=numCISample;
							tmpMaxStrip+=numStrips;
							tmpMaxIdx+=numStrips;
							tmpMaxAll++;
							whichStripMax++;
						}
					}
				}
				EnterCriticalSection(&gData.cs);
				gData.ite=ite;
				gData.sample=sample;
				gData.tKnotNum=basis->tKnotNum;
				gData.T=basis->T;
				gData.ct=GlobalMEMBuf_1st+Npad+Npad;
				gData.tCI=resultChain.tCI;
				gData.tProb=resultChain.tProb;
				gData.sKnotNum=basis->sKnotNum;
				gData.S=basis->S;
				gData.curs=GlobalMEMBuf_1st+Npad;
				gData.sCI=resultChain.sCI;
				gData.sProb=resultChain.sProb;
				if (gData.yMaxT==gData.yMinT)
				{
					int idx;
					r_ippsMaxIndx_32f(gData.ct,N,&gData.yMaxT,&idx);
					r_ippsMinIndx_32f(gData.ct,N,&gData.yMinT,&idx);
					gData.yMinT=gData.yMinT - (gData.yMaxT - gData.yMinT)/10;
					gData.yMaxT=gData.yMaxT+(gData.yMaxT - gData.yMinT)/10;
					r_ippsMaxIndx_32f(gData.curs,N,&gData.yMaxS,&idx);
					r_ippsMinIndx_32f(gData.curs,N,&gData.yMinS,&idx);
					gData.yMinS=gData.yMinS - (gData.yMaxS - gData.yMinS)/10;
					gData.yMaxS=gData.yMaxS+(gData.yMaxS - gData.yMinS)/10;
				}
				if (gData.sleepInterval > 5)
				{
					Sleep(gData.sleepInterval);
					GeneratePlotData_ST();
				}
				else
				{
					QueryPerformanceCounter(&tEnd);
					if ((tEnd.QuadPart - tStart.QuadPart) * 1000/tFrequency.QuadPart > gData.timerInterval)
					{
						GeneratePlotData_ST();
						tStart=tEnd;
					}
				}
				while (gData.status==PAUSE && gData.quit==0)
					SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
				if (gData.quit)
				{
					LeaveCriticalSection(&gData.cs);
					MEM.free_all(&MEM);
					r_vslDeleteStream(&stream);
					return NULL_RET;
				}
				LeaveCriticalSection(&gData.cs);
			}
			r_printf("Model precison parameter is%8.4f; sig2 is%8.4f .\n",modelPar.prec[0],modelPar.sig2);
			GlobalMEMBuf_1st=Xnewterm;
			if (1L||opt.outputType=='F')
			{
				rF32 inv_sample=1.f/sample;
				int sum;
				r_ippsSum_32s_Sfs(resultChain.sProb,opt.N,&sum,0);
				*resultChain.sN=inv_sample * sum;
				r_ippsSum_32s_Sfs(resultChain.tProb,opt.N,&sum,0);
				*resultChain.tN=inv_sample *sum;
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd,resultChain.s,opt.N);
				r_ippsMul_32f(resultChain.s,resultChain.s,GlobalMEMBuf_1st,opt.N);
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd * yInfo.yStd,resultChain.sSD,opt.N);
				r_ippsSub_32f_I(GlobalMEMBuf_1st,resultChain.sSD,opt.N);
				r_ippsSqrt_32f_I(resultChain.sSD,opt.N);
				if (opt.computeCredible)
					r_ippsMulC_32f_I(yInfo.yStd,resultChain.sCI,opt.N+opt.N);
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd,resultChain.t,opt.N);
				r_ippsMul_32f(resultChain.t,resultChain.t,GlobalMEMBuf_1st,opt.N);
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd * yInfo.yStd,resultChain.tSD,opt.N);
				r_ippsSub_32f_I(GlobalMEMBuf_1st,resultChain.tSD,opt.N);
				r_ippsSqrt_32f_I(resultChain.tSD,opt.N);
				r_ippsSubC_32f_I(-yInfo.yMean,resultChain.t,opt.N);
				if (opt.computeCredible)
					r_ippsMulC_32f_I(yInfo.yStd,resultChain.tCI,opt.N+opt.N),
					r_ippsSubC_32f_I(-yInfo.yMean,resultChain.tCI,opt.N+opt.N); 
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd/opt.timeInterval,resultChain.b,opt.N);
				r_ippsMul_32f(resultChain.b,resultChain.b,GlobalMEMBuf_1st,opt.N);
				r_ippsMulC_32f_I(inv_sample  * yInfo.yStd * yInfo.yStd/(opt.timeInterval* opt.timeInterval),resultChain.bSD,opt.N);
				r_ippsSub_32f_I(GlobalMEMBuf_1st,resultChain.bSD,opt.N);
				r_ippsSqrt_32f_I(resultChain.bSD,opt.N);
				if (opt.computeCredible)
					r_ippsMulC_32f_I(yInfo.yStd/opt.timeInterval,resultChain.bCI,opt.N+opt.N);
				*resultChain.marg_lik=*resultChain.marg_lik *inv_sample ;
				*resultChain.sig2=*resultChain.sig2 *inv_sample *yInfo.yStd * yInfo.yStd;
				for (int i=0; i < opt.N; i++)
				{
					*((float *)resultChain.sProb)=inv_sample *(*resultChain.sProb);
					*((float *)resultChain.tProb)=inv_sample *(*resultChain.tProb);
					resultChain.sProb++;  resultChain.tProb++;
				}
				resultChain.sProb=resultChain.sProb - opt.N;
				resultChain.tProb=resultChain.tProb - opt.N;
				if (opt.computeSlopeSign)
				{
					for (int i=0; i < opt.N; i++)
					{
						*((float *)resultChain.bsign)=inv_sample *(*resultChain.bsign);
						resultChain.bsign++;
					}
					resultChain.bsign=resultChain.bsign - opt.N;
				}
				for (int i=0; i < (opt.maxKnotNum_Season+1); i++)
				{
					*((float *)resultChain.sNProb)=inv_sample *(*resultChain.sNProb);
					resultChain.sNProb++;
				}
				resultChain.sNProb=resultChain.sNProb - (opt.maxKnotNum_Season+1);
				for (int i=0; i < (opt.maxKnotNum_Trend+1); i++)
				{
					*((float *)resultChain.tNProb)=inv_sample *(*resultChain.tNProb);
					resultChain.tNProb++;
				}
				resultChain.tNProb=resultChain.tNProb - (opt.maxKnotNum_Trend+1);
				if (opt.computeHarmonicOrder)
				{
					for (int i=0; i < opt.N; i++)
					{
						*((float *)resultChain.horder)=inv_sample *(*resultChain.horder);
						resultChain.horder++;
					}
					resultChain.horder=resultChain.horder - opt.N;
				}
				if (opt.computeTrendOrder)
				{
					for (int i=0; i < opt.N; i++)
					{
						*((float *)resultChain.torder)=inv_sample *(*resultChain.torder);
						resultChain.torder++;
					}
					resultChain.torder=resultChain.torder - opt.N;
				}
			}
			{
				*result.sN+=*resultChain.sN;
				*result.tN+=*resultChain.tN;
				r_ippsAdd_32f_I((float *)resultChain.sProb,(float *)result.sProb,N);
				r_ippsAdd_32f_I((float *)resultChain.tProb,(float *)result.tProb,N);
				r_ippsAdd_32f_I(resultChain.s,result.s,N);
				r_ippsAdd_32f_I(resultChain.t,result.t,N);
				r_ippsAdd_32f_I(resultChain.b,result.b,N);
				if (opt.computeCredible)
					r_ippsAdd_32f_I(resultChain.sCI,result.sCI,N+N),
					r_ippsAdd_32f_I(resultChain.tCI,result.tCI,N+N),
					r_ippsAdd_32f_I(resultChain.bCI,result.bCI,N+N);
				r_ippsAdd_32f_I(resultChain.sSD,result.sSD,N);
				r_ippsAdd_32f_I(resultChain.tSD,result.tSD,N);
				r_ippsAdd_32f_I(resultChain.bSD,result.bSD,N);
				r_ippsAdd_32f_I((float *)resultChain.sNProb,(float *)result.sNProb,opt.maxKnotNum_Season+1);
				r_ippsAdd_32f_I((float *)resultChain.tNProb,(float *)result.tNProb,opt.maxKnotNum_Trend+1);
				*result.marg_lik+=*resultChain.marg_lik;
				*result.sig2+=*resultChain.sig2;
				if (opt.computeSlopeSign)
					r_ippsAdd_32f_I((float *)resultChain.bsign,(float *)result.bsign,N);
				if (opt.computeHarmonicOrder)
					r_ippsAdd_32f_I((float *)resultChain.horder,(float *)result.horder,N);
				if (opt.computeTrendOrder)
					r_ippsAdd_32f_I((float *)resultChain.torder,(float *)result.torder,N);
			}
			chainNumber++;
			EnterCriticalSection(&gData.cs);
			gData.curChainNumber=chainNumber;
			gData.sN=*resultChain.sN;
			gData.tN=*resultChain.tN;
			PostMessage(gData.hwnd,WM_USER+2,0,0);
			LeaveCriticalSection(&gData.cs);
		}
		if (opt.chainNumber >=2)
		{
			float invChainNumber=1.f/(float)opt.chainNumber;
			*result.sN=*result.sN*invChainNumber;
			*result.tN=*result.tN*invChainNumber;
			r_ippsMulC_32f_I(invChainNumber,(float *)result.sProb,N);
			r_ippsMulC_32f_I(invChainNumber,(float *)result.tProb,N);
			r_ippsMulC_32f_I(invChainNumber,result.s,N);
			r_ippsMulC_32f_I(invChainNumber,result.t,N);
			r_ippsMulC_32f_I(invChainNumber,result.b,N);
			if (opt.computeCredible)
				r_ippsMulC_32f_I(invChainNumber,result.sCI,N+N),
				r_ippsMulC_32f_I(invChainNumber,result.tCI,N+N),
				r_ippsMulC_32f_I(invChainNumber,result.bCI,N+N);
			r_ippsMulC_32f_I(invChainNumber,result.sSD,N);
			r_ippsMulC_32f_I(invChainNumber,result.tSD,N);
			r_ippsMulC_32f_I(invChainNumber,result.bSD,N);
			r_ippsMulC_32f_I(invChainNumber,(float *)result.sNProb,opt.maxKnotNum_Season+1);
			r_ippsMulC_32f_I(invChainNumber,(float *)result.tNProb,opt.maxKnotNum_Trend+1);
			*result.marg_lik=*result.marg_lik*invChainNumber;
			*result.sig2=*result.sig2*invChainNumber;
			if (opt.computeSlopeSign)
				r_ippsMulC_32f_I(invChainNumber,(float *)result.bsign,N);
			if (opt.computeHarmonicOrder)
				r_ippsMulC_32f_I(invChainNumber,(float *)result.horder,N);
			if (opt.computeTrendOrder)
				r_ippsMulC_32f_I(invChainNumber,(float *)result.torder,N);
		}
#if defined(MSVC_COMPILER)
		QueryPerformanceCounter(&t2);
		elapsedTime=elapsedTime+(long long)((double)(t2.QuadPart - t1.QuadPart)* 1000.0/(double)Frequency.QuadPart);
#elif (defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER))  && ! (defined(__APPLE__)||defined(__MACH__))
		clock_gettime(CLOCK_REALTIME,&t2);
		elapsedTime=elapsedTime+(t2.tv_sec - t1.tv_sec) * 1000000000LL+(t2.tv_nsec - t1.tv_nsec);
#endif
		r_printf("Time spent is%d \n",elapsedTime);
		if (opt.outputType=='s')
		{
			*(matOutput.sN+(pixelIndex - 1))=*result.sN;
			memcpy(matOutput.sProb+(pixelIndex - 1)*opt.N,result.sProb,sizeof(int32_t)*opt.N);
			*(matOutput.tN+(pixelIndex - 1))=*result.tN;
			memcpy(matOutput.tProb+(pixelIndex - 1)*opt.N,result.tProb,sizeof(int32_t)*opt.N);
			memcpy(matOutput.sNProb+(pixelIndex - 1)*(opt.maxKnotNum_Season+1),result.sNProb,sizeof(int32_t)* (opt.maxKnotNum_Season+1));
			memcpy(matOutput.tNProb+(pixelIndex - 1)*(opt.maxKnotNum_Trend+1),result.tNProb,sizeof(int32_t)* (opt.maxKnotNum_Trend+1));
			matOutput.marg_lik[pixelIndex - 1]=*result.marg_lik;
			memcpy(matOutput.s+(pixelIndex - 1)*opt.N,result.s,sizeof(float)*opt.N);
			if (opt.computeCredible)	memcpy(matOutput.sCI+(pixelIndex - 1)*(opt.N+opt.N),result.sCI,sizeof(float)*(opt.N+opt.N));
			memcpy(matOutput.sSD+(pixelIndex - 1)*opt.N,result.sSD,sizeof(float)*opt.N);
			memcpy(matOutput.t+(pixelIndex - 1)*opt.N,result.t,sizeof(float)*opt.N);
			if (opt.computeCredible) memcpy(matOutput.tCI+(pixelIndex - 1)*(opt.N+opt.N),result.tCI,sizeof(float)*(opt.N+opt.N));
			memcpy(matOutput.tSD+(pixelIndex - 1)*opt.N,result.tSD,sizeof(float)*opt.N);
			memcpy(matOutput.b+(pixelIndex - 1)*opt.N,result.b,sizeof(float)*opt.N);
			if (opt.computeCredible) memcpy(matOutput.bCI+(pixelIndex - 1)*(opt.N+opt.N),result.bCI,sizeof(float)*(opt.N+opt.N));
			memcpy(matOutput.bSD+(pixelIndex - 1)*opt.N,result.bSD,sizeof(float)*opt.N);
			matOutput.sig2[pixelIndex - 1]=*result.sig2;
			if (opt.computeSlopeSign)	   memcpy(matOutput.bsign+(pixelIndex - 1)*opt.N,result.bsign,sizeof(float)*opt.N);
			if (opt.computeHarmonicOrder)  memcpy(matOutput.horder+(pixelIndex - 1)*opt.N,result.horder,sizeof(float)*opt.N);
			if (opt.computeTrendOrder)     memcpy(matOutput.torder+(pixelIndex - 1)*opt.N,result.torder,sizeof(float)*opt.N);
			if (opt.computeChangepoints)
			{
				F32PTR	mem=Xnewterm; 
				float		        threshold=0.001f;
				int32_t * _restrict cptList=(int32_t * ) mem+5 * N;
				F32PTR  cptCIList=mem+6 * N;
				int32_t             cptNumber=(int32_t)round((double)*result.sN);
				int32_t             trueCptNumber;
				trueCptNumber=find_changepoint((float  *)result.sProb,mem,threshold,cptList,cptCIList,opt.N,opt.minSepDist_Season,cptNumber);
				float T0=(float)opt.startTime;
				float dT=(float)opt.timeInterval;
				int32_t baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Season;
				for (int i=0; i <trueCptNumber; i++)
				{
					*(matOutput.scp+baseIdx+i)=(float)(*(cptList+i))* dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Season; i++)
				{
					*(matOutput.scp+baseIdx+i)=nan;
				}
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Season * 2;
				for (int i=0; i <trueCptNumber; i++)
				{
					*(matOutput.scpCI+baseIdx+i)=(float)(*(cptCIList+i))* dT+T0;
					*(matOutput.scpCI+baseIdx+opt.maxKnotNum_Season+i)=(float)(*(cptCIList+trueCptNumber+i))* dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Season; i++)
				{
					*(matOutput.scpCI+baseIdx+i)=nan;
					*(matOutput.scpCI+baseIdx+opt.maxKnotNum_Season+i)=nan;
				}
				cptNumber=(int32_t)round((double)*result.tN);
				trueCptNumber=find_changepoint((float  *)result.tProb,mem,threshold,cptList,cptCIList,opt.N,opt.minSepDist_Trend,cptNumber);
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Trend;
				for (int i=0; i <trueCptNumber; i++)
				{
					*(matOutput.tcp+baseIdx+i)=(float)(*(cptList+i)) * dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Trend; i++)
				{
					*(matOutput.tcp+baseIdx+i)=nan;
				}
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Trend * 2;
				for (int i=0; i <trueCptNumber; i++)
				{
					*(matOutput.tcpCI+baseIdx+i)=(float)(*(cptCIList+i))* dT+T0;
					*(matOutput.tcpCI+baseIdx+opt.maxKnotNum_Trend+i)=(float)(*(cptCIList+trueCptNumber+i))* dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Trend; i++)
				{
					*(matOutput.tcpCI+baseIdx+i)=nan;
					*(matOutput.tcpCI+baseIdx+opt.maxKnotNum_Trend+i)=nan;
				}
			}
		}
		else if (opt.outputType=='d')
		{
			double T0=(double)opt.startTime;
			double dT=(double)opt.timeInterval;
			for (rI32 i=0; i < N; i++)
			{ 
				*((double *)matOutput.time+i)=(double)i*dT+T0;
			}
			int64_t baseIdx=pixelIndex - 1;
			*((double*)matOutput.sN+baseIdx)=(double)*result.sN;
			*((double*)matOutput.tN+baseIdx)=(double)*result.tN;
			*((double*)matOutput.marg_lik+baseIdx)=(double)*result.marg_lik;
			*((double*)matOutput.sig2+baseIdx)=(double)*result.sig2;
			baseIdx=(pixelIndex - 1)*(opt.maxKnotNum_Season+1);
			for (int i=0; i < (opt.maxKnotNum_Season+1); i++)
			{ 
				*((double *)matOutput.sNProb+baseIdx+i)=(double)(*(float *)(result.sNProb+i));
			}
			baseIdx=(pixelIndex - 1)*(opt.maxKnotNum_Trend+1);
			for (int i=0; i < (opt.maxKnotNum_Trend+1); i++)
			{ 
				*((double *)matOutput.tNProb+baseIdx+i)=(double)(*(float *)(result.tNProb+i));
			}
			baseIdx=(pixelIndex - 1)*N;
			for (int i=0; i < N; i++)
			{
				*((double *)matOutput.sProb+baseIdx+i)=(double)(*(float *)(result.sProb+i));;
				*((double *)matOutput.tProb+baseIdx+i)=(double)(*(float *)(result.tProb+i));;
				*((double *)matOutput.s+baseIdx+i)=(double)*(result.s+i);
				*((double *)matOutput.sSD+baseIdx+i)=(double)*(result.sSD+i);
				*((double *)matOutput.t+baseIdx+i)=(double)*(result.t+i);
				*((double *)matOutput.tSD+baseIdx+i)=(double)*(result.tSD+i);
				*((double *)matOutput.b+baseIdx+i)=(double)*(result.b+i);
				*((double *)matOutput.bSD+baseIdx+i)=(double)*(result.bSD+i);
			}
			if (opt.computeCredible)
			{
				int N2=N+N;
				baseIdx=(pixelIndex - 1)*N2;
				for (int i=0; i < N2; i++)
				{
					*((double *)matOutput.sCI+baseIdx+i)=(double)*(result.sCI+i);
					*((double *)matOutput.tCI+baseIdx+i)=(double)*(result.tCI+i);
					*((double *)matOutput.bCI+baseIdx+i)=(double)*(result.bCI+i);
				}
			}
			if (opt.computeSlopeSign)
			{
				baseIdx=(pixelIndex - 1)*N;
				for (int i=0; i < N; i++)
				{
					*((double *)matOutput.bsign+baseIdx+i)=(double)* ((float *)(result.bsign+i));
				}
			}
			if (opt.computeHarmonicOrder)
			{
				baseIdx=(pixelIndex - 1)*N;
				for (int i=0; i < N; i++)
				{
					*((double *)matOutput.horder+baseIdx+i)=(double)* ((float *)(result.horder+i));
				}
			}
			if (opt.computeTrendOrder)
			{
				baseIdx=(pixelIndex - 1)*N;
				for (int i=0; i < N; i++)
				{
					*((double *)matOutput.torder+baseIdx+i)=(double)* ((float *)(result.torder+i));
				}
			}
			if (opt.computeChangepoints)
			{
				F32PTR mem=Xnewterm; 
				float threshold=0.001f;
				int32_t * _restrict cptList=(int32_t *)mem+5 * N;
				F32PTR cptCIList=mem+6 * N;
				int32_t   cptNumber=(int32_t)round((double)*result.sN);
				int32_t   trueCptNumber;
				trueCptNumber=find_changepoint((float  *)result.sProb,mem,threshold,cptList,cptCIList,opt.N,opt.minSepDist_Season,cptNumber);
				double T0=(double)opt.startTime;
				double dT=(double)opt.timeInterval;
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Season;
				for (int i=0; i <trueCptNumber; i++)
				{
					*((double *)matOutput.scp+baseIdx+i)=(double)(*(int32_t *)(cptList+i)) *dT+T0;
				}
#if !defined(NA_REAL)
#define NA_REAL nan
#endif
				for (int i=trueCptNumber; i < opt.maxKnotNum_Season; i++)
				{
					*((double *)matOutput.scp+baseIdx+i)=NA_REAL;
				}
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Season * 2;
				for (int i=0; i <trueCptNumber; i++)
				{
					*((double *)matOutput.scpCI+baseIdx+i)=(double)(*(cptCIList+i)) *dT+T0;
					*((double *)matOutput.scpCI+baseIdx+opt.maxKnotNum_Season+i)=(double)(*(cptCIList+trueCptNumber+i)) *dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Season; i++)
				{
					*((double *)matOutput.scpCI+baseIdx+i)=NA_REAL;
					*((double *)matOutput.scpCI+baseIdx+opt.maxKnotNum_Season+i)=NA_REAL;
				}
				cptNumber=(int32_t)round((double)*result.tN);
				trueCptNumber=find_changepoint((float  *)result.tProb,mem,threshold,cptList,cptCIList,opt.N,opt.minSepDist_Trend,cptNumber);
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Trend;
				for (int i=0; i <trueCptNumber; i++)
				{
					*((double *)matOutput.tcp+baseIdx+i)=(double)(*(int32_t *)(cptList+i)) *dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Trend; i++)
				{
					*((double *)matOutput.tcp+baseIdx+i)=NA_REAL;
				}
				baseIdx=(pixelIndex - 1)*opt.maxKnotNum_Trend * 2;
				for (int i=0; i <trueCptNumber; i++)
				{
					*((double *)matOutput.tcpCI+baseIdx+i)=(double)(*(cptCIList+i)) *dT+T0;
					*((double *)matOutput.tcpCI+baseIdx+opt.maxKnotNum_Trend+i)=(double)(*(cptCIList+trueCptNumber+i)) *dT+T0;
				}
				for (int i=trueCptNumber; i < opt.maxKnotNum_Trend; i++)
				{
					*((double *)matOutput.tcpCI+baseIdx+i)=NA_REAL;
					*((double *)matOutput.tcpCI+baseIdx+opt.maxKnotNum_Trend+i)=NA_REAL;
				}
			}
		}
		r_printf("SEED is%d\n",opt.seed);
		r_printf("TREND is%d%d%d%d%d\n",accT[0],accT[1],accT[2],accT[3],accT[4]);
		r_printf("SEASO is%d%d%d%d%d\n",accS[0],accS[1],accS[2],accS[3],accS[4]);
		if (opt.outputType !='F') continue;
#if PTHREAD_INOUT==1
		pthread_mutex_lock(&MUTEX_WRITE);
		if (!DATA_AVAILABLE_WRITE)
		{
			threadParWrite.curIdx=pixelIndex;
			{
				*threadParWrite.output.sN=*result.sN;
				memcpy(threadParWrite.output.sProb,result.sProb,sizeof(int32_t)*N);
				*threadParWrite.output.tN=*result.tN;
				memcpy(threadParWrite.output.tProb,result.tProb,sizeof(int32_t)*N);
				memcpy(threadParWrite.output.sNProb,result.sNProb,sizeof(int32_t)* (opt.maxKnotNum_Season+1));
				memcpy(threadParWrite.output.tNProb,result.tNProb,sizeof(int32_t)* (opt.maxKnotNum_Trend+1));
				r_cblas_scopy(N,result.s,1,threadParWrite.output.s,1);
				if (opt.computeCredible)
					r_cblas_scopy(N+N,result.sCI,1,threadParWrite.output.sCI,1);
				r_cblas_scopy(N,result.sSD,1,threadParWrite.output.sSD,1);
				r_cblas_scopy(N,result.t,1,threadParWrite.output.t,1);
				if (opt.computeCredible)
					r_cblas_scopy(N+N,result.tCI,1,threadParWrite.output.tCI,1);
				r_cblas_scopy(N,result.tSD,1,threadParWrite.output.tSD,1);
				r_cblas_scopy(N,result.b,1,threadParWrite.output.b,1);
				if (opt.computeCredible)
					r_cblas_scopy(N+N,result.bCI,1,threadParWrite.output.bCI,1);
				if (opt.computeSlopeSign)
					r_cblas_scopy(N,result.bsign,1,threadParWrite.output.bsign,1);
				r_cblas_scopy(N,result.bSD,1,threadParWrite.output.bSD,1);
			}
			DATA_AVAILABLE_WRITE=1;
			pthread_mutex_unlock(&MUTEX_WRITE);
			pthread_cond_signal(&CONDITION_WRITE);
		}
		else if (DATA_AVAILABLE_WRITE)
		{
			while (DATA_AVAILABLE_WRITE)
				pthread_cond_wait(&CONDITION_WRITE,&MUTEX_WRITE);
			threadParWrite.curIdx=pixelIndex;
			*threadParWrite.output.sN=*result.sN;
			memcpy(threadParWrite.output.sProb,result.sProb,sizeof(int32_t)*N);
			*threadParWrite.output.tN=*result.tN;
			memcpy(threadParWrite.output.tProb,result.tProb,sizeof(int32_t)*N);
			memcpy(threadParWrite.output.sNProb,result.sNProb,sizeof(int32_t)* (opt.maxKnotNum_Season+1));
			memcpy(threadParWrite.output.tNProb,result.tNProb,sizeof(int32_t)* (opt.maxKnotNum_Trend+1));
			r_cblas_scopy(N,result.s,1,threadParWrite.output.s,1);
			if (opt.computeCredible)
				r_cblas_scopy(N+N,result.sCI,1,threadParWrite.output.sCI,1);
			r_cblas_scopy(N,result.sSD,1,threadParWrite.output.sSD,1);
			r_cblas_scopy(N,result.t,1,threadParWrite.output.t,1);
			if (opt.computeCredible)
				r_cblas_scopy(N+N,result.tCI,1,threadParWrite.output.tCI,1);
			r_cblas_scopy(N,result.tSD,1,threadParWrite.output.tSD,1);
			r_cblas_scopy(N,result.b,1,threadParWrite.output.b,1);
			if (opt.computeCredible)
				r_cblas_scopy(N+N,result.bCI,1,threadParWrite.output.bCI,1);
			r_cblas_scopy(N,result.bSD,1,threadParWrite.output.bSD,1);
			if (opt.computeSlopeSign)
				r_cblas_scopy(N,result.bsign,1,threadParWrite.output.bsign,1);
			DATA_AVAILABLE_WRITE=1;
			pthread_mutex_unlock(&MUTEX_WRITE);
			pthread_cond_signal(&CONDITION_WRITE);
		}
#else
		r_printf("Writing Output..%d\n",pixelIndex);
		fwrite(result.sN,sizeof(float),1,file.sN);
		fwrite(result.tN,sizeof(float),1,file.tN);
		fwrite(result.sNProb,sizeof(int32_t),opt.maxKnotNum_Season+1,file.sNProb);
		fwrite(result.tNProb,sizeof(int32_t),opt.maxKnotNum_Trend+1,file.tNProb);
		fwrite(result.sProb,sizeof(int32_t),opt.N,file.sProb);
		fwrite(result.tProb,sizeof(int32_t),opt.N,file.tProb);
		fwrite(result.s,sizeof(float),opt.N,file.s);
		if (opt.computeCredible)  fwrite(result.sCI,sizeof(float),opt.N * 2,file.sCI);
		fwrite(result.sSD,sizeof(float),opt.N,file.sSD);
		fwrite(result.t,sizeof(float),opt.N,file.t);
		if (opt.computeCredible)  fwrite(result.tCI,sizeof(float),opt.N * 2,file.tCI);
		fwrite(result.tSD,sizeof(float),opt.N,file.tSD);
		fwrite(result.b,sizeof(float),opt.N,file.b);
		if (opt.computeCredible)  fwrite(result.bCI,sizeof(float),opt.N * 2,file.bCI);
		fwrite(result.bSD,sizeof(float),opt.N,file.bSD);
		if (opt.computeSlopeSign)  fwrite(result.bsign,sizeof(float),opt.N,file.bsign);
#endif
	}	
	r_vslDeleteStream(&stream);
#if PTHREAD_INOUT==1
	if (opt.outputType=='F')
	{
		pthread_join(THREADID_WRITE,NULL);
		pthread_mutex_destroy(&MUTEX_WRITE);
		pthread_cond_destroy(&CONDITION_WRITE);
	}
	if (opt.inputType=='F')
	{
		pthread_join(THREADID_READ,NULL);
		pthread_mutex_destroy(&MUTEX_READ);
		pthread_cond_destroy(&CONDITION_READ);
	}
#else
	if (opt.outputType=='F')
	{
		fclose(file.sN);		fclose(file.tN);		fclose(file.sNProb);		fclose(file.tNProb);		fclose(file.sProb);
		fclose(file.tProb);		fclose(file.s);		fclose(file.t);		fclose(file.b);		if (opt.computeCredible)  fclose(file.sCI);
		if (opt.computeCredible)  
			fclose(file.tCI); 
		if (opt.computeCredible)  
			fclose(file.bCI);
		if (opt.computeSlopeSign)
			fclose(file.bsign);
		fclose(file.sSD);
		fclose(file.tSD); 	fclose(file.bSD);
	}
	if (opt.inputType=='F')
		fclose(file.infp);
#endif
	MEM.free_all(&MEM);
	PostMessage(gData.hwnd,WM_USER+1,0,0);
	EnterCriticalSection(&gData.cs);
	gData.t=NULL;
	gData.y=NULL;
	gData.s=NULL;
	gData.rowsMissing=NULL;
	while (gData.status !=DONE)
		SleepConditionVariableCS(&gData.cv,&gData.cs,INFINITE);
	LeaveCriticalSection(&gData.cs);
	return 1;
} 
#else
static char fileID  UNUSED_DECORATOR='c';
#endif
ENABLE_MANY_WARNINGS
