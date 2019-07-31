#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "abc_001_config.h"
#include "beast_lib.h"
#if  MYRAND_LIBRARY==1
	#include "abc_rand_pcg.h"
#elif MKLRAND_LIBRARY==1
	#include "abc_rand_mkl.h"
	extern VSLStreamStatePtr stream;
#endif
void preCompute_Xmars_terms(F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,int N,float PERIOD,int maxSeasonOrder,int maxTrendOrder)
{
	register F32PTR ptr;
	if (SEASON_TERMS !=NULL)
	{
		ptr=SEASON_TERMS;
		float freq_factor=2.0f *  3.141592653589793f/PERIOD;
		for (int i=1; i <=maxSeasonOrder; i++)
		{
			fill_1_to_N(ptr,N);
			r_cblas_sscal(N,freq_factor*i,ptr,1);
			r_cblas_scopy(N,ptr,1,ptr+N,1);
			r_vmsSin(N,ptr,ptr,VML_HA); 
			normalize(ptr,N);
			r_vmsCos(N,ptr+N,ptr+N,VML_HA);
			normalize(ptr+N,N);
			ptr=ptr+N+N;
		}
	}
	if (TREND_TERMS !=NULL)
	{
		ptr=TREND_TERMS;
		for (int i=1; i <=(maxTrendOrder+1); i++)
		{
			if (i==1)
			{
				r_ippsSet_32f(1,ptr,N);
			}
			else
			{
				fill_1_to_N(ptr,N);
				r_vsPowx(N,ptr,(float)(i - 1),ptr);
				normalize(ptr,N);
			}
			ptr+=N;
		}
	}
}
void preCompute_Xmars_terms_fast(F32PTR SEASON_TERMS,F32PTR SEASON_SQR_CSUM,F32PTR TREND_TERMS,F32PTR INV_SQR,F32PTR COEFF_A,F32PTR COEFF_B,
	int N,float PERIOD,int maxSeasonOrder,int maxTrendOrder)
{
	rF32PTR ptr;
	rF32PTR ptr1;
	rF32PTR ptr2;
	if (SEASON_TERMS !=NULL)
	{
		rF32 freq_factor=2.0f *  3.141592653589793f/PERIOD;
		ptr=SEASON_TERMS;
		ptr1=SEASON_SQR_CSUM;
		ptr2=SEASON_SQR_CSUM+(N+1);
		for (rI32 order=1; order <=maxSeasonOrder; order++)
		{
			fill_1_to_N(ptr,N);
			r_cblas_sscal(N,freq_factor*order,ptr,1);
			r_cblas_scopy(N,ptr,1,ptr+N,1);
			r_vmsSin(N,ptr,ptr,VML_HA);
			rF32 dotVal=DOT(N,ptr,ptr);
			r_ippsMulC_32f_I(1/sqrtf(dotVal/N),ptr,N);
			r_vmsCos(N,ptr+N,ptr+N,VML_HA);
			dotVal=DOT(N,ptr+N,ptr+N);
			r_ippsMulC_32f_I(1/sqrtf(dotVal/N),ptr+N,N);
			*ptr1++=0.f;
			*ptr2++=0.f;
			rF32 csum1=0.f;
			rF32 csum2=0.f;
			rF32PTR ptr_N=ptr+N;
			for (rI32 i=0; i < N; i++)
			{
				rF32 x=*ptr++;
				csum1+=x*x;
				*ptr1++=csum1;
				rF32 y=*ptr_N++;
				csum2+=y*y;
				*ptr2++=csum2;
			}
			ptr=ptr_N;
			ptr1=ptr2;
			ptr2=ptr1+(N+1);
		}
	}
	if (TREND_TERMS !=NULL)
	{
		ptr=TREND_TERMS;
		for (rI32 i=1; i <=(maxTrendOrder+1); i++)
		{
			if (i==1)
				r_ippsSet_32f(1,ptr,N);
			else			
				fill_1_to_N(ptr,N),
				r_vsPowx(N,ptr,(float)(i - 1),ptr),
				normalize(ptr,N);
			ptr+=N;
		}
	}
	if (INV_SQR !=NULL)
	{
		rF32 sqrt_N=fast_sqrt((float) N);
		for (rI32 i=0; i < N; i++)
			INV_SQR[i]=sqrt_N/fast_sqrt((float)(i+1));
	}
	if (COEFF_A !=NULL &&  COEFF_B !=NULL)
	{
		COEFF_B[1 - 1]=0;
		COEFF_A[1 - 1]=fast_sqrt(N);
		for (rI32 n=2; n <=N; n++)
		{
			rF32 sum=(1L+n)/2.f;
			rF32 b=1.f/(  (n+1L)*(2L*n+1)/6.f - sum*sum);
			b=fast_sqrt(b*N/n);
			COEFF_B[n - 1]=b;
			COEFF_A[n - 1]=-b*sum;
		}
	}
}
void preCompute_scale_factor(F32PTR scaleFactorSeason,F32PTR scaleFactorTrend,int32_t N,int maxKnotNum_Season,
	int maxKnotNum_Trend,int minSepDist_Season,int minSepDist_Trend,F32PTR mem1,F32PTR mem2)
{
	float N_tmp,tmp1,tmp2;
	if (scaleFactorSeason !=NULL)
	{
		for (int k=0; k <=maxKnotNum_Season; k++)
		{
			N_tmp=N - (k+1)*(minSepDist_Season - 1) - 1.f;
			if (k==0)
			{
				*mem1=1.0f;
				tmp1=logf(1.0f);
			}
			else
			{
				fill_1_to_N(mem1,k);
				r_ippsSubC_32f_I(1.f,mem1,k);
				r_ippsSubCRev_32f_I(N_tmp,mem1,k);
				r_ippsLn_32f_I(mem1,k);
				r_ippsSum_32f(mem1,k,&tmp1,ippAlgHintAccurate);
			}
			N_tmp=N - (k+2)*(minSepDist_Season - 1) - 1.f;
			fill_1_to_N(mem2,k+1);
			r_ippsSubC_32f_I(1.f,mem2,k+1);
			r_ippsSubCRev_32f_I(N_tmp,mem2,k+1);
			r_ippsLn_32f_I(mem2,k+1);
			r_ippsSum_32f(mem2,k+1,&tmp2,ippAlgHintAccurate);
			scaleFactorSeason[k]=(N - (k+2)*minSepDist_Season+1) *expf(tmp1 - tmp2);
		}
	}
	if (scaleFactorTrend !=NULL)
	{
		for (int k=0; k <=maxKnotNum_Trend; k++)
		{
			N_tmp=N - (k+1)*(minSepDist_Trend - 1) - 1.f;
			if (k==0)
			{
				*mem1=1.0f;
				tmp1=logf(1.0f);
			}
			else
			{
				fill_1_to_N(mem1,k);
				r_ippsSubC_32f_I(1.f,mem1,k);
				r_ippsSubCRev_32f_I(N_tmp,mem1,k);
				r_ippsLn_32f_I(mem1,k);
				r_ippsSum_32f(mem1,k,&tmp1,ippAlgHintAccurate);
			}
			N_tmp=N - (k+2)*(minSepDist_Trend - 1) - 1.f;
			fill_1_to_N(mem2,k+1);
			r_ippsSubC_32f_I(1.f,mem2,k+1);
			r_ippsSubCRev_32f_I(N_tmp,mem2,k+1);
			r_ippsLn_32f_I(mem2,k+1);
			r_ippsSum_32f(mem2,k+1,&tmp2,ippAlgHintAccurate);
			scaleFactorTrend[k]=(N - (k+2)*minSepDist_Trend+1) *expf(tmp1 - tmp2);
		}
	}
}
void fetch_next_time_series(YINFO * _restrict yInfo,void * _restrict yInputData,int idx,F32PTR GlobalMEMBuf_1st,bool isSingleyInput,int32_t N,float omissionValue)
{
	uint32_t Npad=(uint32_t)ceil((float)N/8.0f) * 8; 
	rU32PTR rowsMissing=yInfo->rowsMissing;
	rF32PTR buf=GlobalMEMBuf_1st;
	int     nMissing=0;
	if (isSingleyInput)
	{
		F32PTR yInput=(F32PTR) yInputData+(idx - 1)*N;
		for (rI32 i=1; i <=N; i++)
		{
			rF32 yValue=*yInput++;
			if (fabs(yValue - omissionValue) < 1e-10||yValue !=yValue)
				rowsMissing[nMissing++]=i;
			else
				*buf++=yValue;
		}
	}
	else
	{
		double *yInput=(double *)yInputData+(idx - 1)*N;
		for (rI32 i=1; i <=N; i++)
		{
			rF32 yValue=(float) *yInput++;
			if (fabsf(yValue - omissionValue) < 1e-5f||yValue !=yValue)
				rowsMissing[nMissing++]=i;
			else
				*buf++=yValue;
		}
	}
	yInfo->nMissing=nMissing;
	rI32 n=yInfo->n=N - nMissing;
	r_ippsMeanStdDev_32f(GlobalMEMBuf_1st,n,&yInfo->yMean,&yInfo->yStd,ippAlgHintFast);
	normalize(GlobalMEMBuf_1st,n);
	buf=GlobalMEMBuf_1st;
	{
		rF32PTR  Y=yInfo->Y;
		rI32     jOmit=0;
		for (int i=1; i <=N; i++)
		{
			if (jOmit < nMissing && i==rowsMissing[jOmit])
				Y[i - 1]=0,jOmit++;			
			else
				Y[i - 1]=*buf++;			
		}		
		r_ippsSet_32f(0,Y+N,Npad - N);
		yInfo->YtY=DOT(N,Y,Y);
	}
}
void fetch_next_time_series1(YINFO * _restrict yInfo,void * _restrict yInputData,int idx,F32PTR GlobalMEMBuf_1st,char inputType,int32_t N,float omissionValue)
{
	uint32_t Npad=(uint32_t)ceil((float)N/8.0f) * 8;
	rU32PTR rowsMissing=yInfo->rowsMissing;
	rF32PTR buf=GlobalMEMBuf_1st;
	int     nMissing=0;
	if (inputType=='s')
	{
		F32PTR yInput=(F32PTR)yInputData+(idx - 1)*N;
		for (rI32 i=1; i <=N; i++)
		{
			rF32 yValue=*yInput++;
			if (fabs(yValue - omissionValue) < 1e-10f||yValue !=yValue)
				rowsMissing[nMissing++]=i;
			else
				*buf++=yValue;
		}
	}
	else if (inputType=='d')
	{
		rF64PTR yInput=(double *)yInputData+(idx - 1)*N;
		for (rI32 i=1; i <=N; i++)
		{
			rF32 yValue=(float)*yInput++;
			if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
				rowsMissing[nMissing++]=i;
			else
				*buf++=yValue;
		}
	}
	else if (inputType=='4')
	{
		rI32PTR yInput=(int32_t *) yInputData+(idx - 1)*N;
		for (rI32 i=1; i <=N; i++)
		{
			rF32 yValue=(float)*yInput++;
			if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
				rowsMissing[nMissing++]=i;
			else
				*buf++=yValue;
		}
	}
	yInfo->nMissing=nMissing;
	rI32 n=yInfo->n=N - nMissing;
	r_ippsMeanStdDev_32f(GlobalMEMBuf_1st,n,&yInfo->yMean,&yInfo->yStd,ippAlgHintFast);
	normalize(GlobalMEMBuf_1st,n);
	buf=GlobalMEMBuf_1st;
	{
		rF32PTR  Y=yInfo->Y;
		rI32     jOmit=0;
		for (int i=1; i <=N; i++)
		{
			if (jOmit < nMissing && i==rowsMissing[jOmit])
				Y[i - 1]=0,jOmit++;
			else
				Y[i - 1]=*buf++;
		}
		r_ippsSet_32f(0,Y+N,Npad - N);
		yInfo->YtY=DOT(N,Y,Y);
	}
}
void fetch_next_time_series2(YINFO * _restrict yInfo,void * _restrict yInputData,int idx,F32PTR GlobalMEMBuf_1st,char inputType,int32_t N,float omissionValue,Options * _restrict opt)
{
	uint32_t Npad=(uint32_t)ceil((float)N/8.0f) * 8;
	rU32PTR rowsMissing=yInfo->rowsMissing;
	rF32PTR buf=GlobalMEMBuf_1st;
	int     nMissing=0;
	if (opt->isInput3DStack !=1L)
	{
		if (inputType=='s')
		{
			F32PTR yInput=(F32PTR)yInputData+(idx - 1)*N;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=*yInput++;
				if (fabs(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		else if (inputType=='d')
		{
			rF64PTR yInput=(double *)yInputData+(idx - 1)*N;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=(float)*yInput++;
				if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		else if (inputType=='4')
		{
			rI32PTR yInput=(int32_t *)yInputData+(idx - 1)*N;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=(float)*yInput++;
				if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		yInfo->nMissing=nMissing;
		rI32 n=yInfo->n=N - nMissing;
		r_ippsMeanStdDev_32f(GlobalMEMBuf_1st,n,&yInfo->yMean,&yInfo->yStd,ippAlgHintFast);
		normalize(GlobalMEMBuf_1st,n);
		buf=GlobalMEMBuf_1st;
		{
			rF32PTR  Y=yInfo->Y;
			rI32     jOmit=0;
			for (int i=1; i <=N; i++)
			{
				if (jOmit < nMissing && i==rowsMissing[jOmit])
					Y[i - 1]=0,jOmit++;
				else
					Y[i - 1]=*buf++;
			}
			r_ippsSet_32f(0,Y+N,Npad - N);
			yInfo->YtY=DOT(N,Y,Y);
		}
	}
	if (opt->isInput3DStack==1L)
	{
		int M=opt->M;
		int L=opt->L;
		int     stride;
		int64_t start;
		if (opt->timeDimensionIndex==1)
		{
			stride=1L;
			start=(idx - 1)*N;
		}
		else if (opt->timeDimensionIndex==2)
		{
			stride=M;
			int m,l;
			l=(idx - 1)/M;
			m=idx - l*M;
			l=l+1;
			start=(l - 1)*(N*M)+m - 1;
		}
		else if (opt->timeDimensionIndex==3)
		{
			stride=M*L;
			int m,l;
			l=(idx - 1)/M;
			m=idx - l*M;
			l=l+1;
			start=(l-1)*M+m-1;
		}
		if (inputType=='s')
		{
			F32PTR yInput=(F32PTR) yInputData+start;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=*yInput;
				yInput+=stride;
				if (fabs(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		else if (inputType=='d')
		{
			rF64PTR yInput=(double *)yInputData+start;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=(float)*yInput;
				yInput+=stride;
				if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		else if (inputType=='4')
		{
			rI32PTR yInput=(int32_t *)yInputData+start;
			for (rI32 i=1; i <=N; i++)
			{
				rF32 yValue=(float)*yInput;
				yInput+=stride;
				if (fabsf(yValue - omissionValue) < 1e-10f||yValue !=yValue)
					rowsMissing[nMissing++]=i;
				else
					*buf++=yValue;
			}
		}
		yInfo->nMissing=nMissing;
		rI32 n=yInfo->n=N - nMissing;
		r_ippsMeanStdDev_32f(GlobalMEMBuf_1st,n,&yInfo->yMean,&yInfo->yStd,ippAlgHintFast);
		normalize(GlobalMEMBuf_1st,n);
		buf=GlobalMEMBuf_1st;
		{
			rF32PTR  Y=yInfo->Y;
			rI32     jOmit=0;
			for (int i=1; i <=N; i++)
			{
				if (jOmit < nMissing && i==rowsMissing[jOmit])
					Y[i - 1]=0,jOmit++;
				else
					Y[i - 1]=*buf++;
			}
			r_ippsSet_32f(0,Y+N,Npad - N);
			yInfo->YtY=DOT(N,Y,Y);
		}
	}
}
void convert_basis_both(struct BASIS * _restrict basis)
{
	int numOfSeg;
	int k;
	int order;
	int i,j;
	register int16_t *  _restrict KS,*_restrict KE;  
	register uint8_t * _restrict vecTermType;
	vecTermType=basis->termType;
	numOfSeg=basis->sKnotNum+1; 
	KS=basis->sks;
	KE=basis->ske;
	k=1; 
	uint8_t * _restrict vecOrder;
	vecOrder=basis->sOrder;
	for (i=1; i <=numOfSeg; i++)
	{
		*KS++=k;
		order=(int)vecOrder[i - 1];
		for (j=1; j <=order; j++)
		{
			*vecTermType++=0; 
			k++;
			*vecTermType++=0; 
			k++;
		}
		*KE++=k - 1;
	}
	basis->k_SN=k - 1;
	numOfSeg=basis->tKnotNum+1;
	basis->k_const=numOfSeg;
	KS=basis->tks;
	KE=basis->tke;
	vecOrder=basis->tOrder;
	for (i=1; i <=numOfSeg; i++)
	{
		*KS++=k;
		order=vecOrder[i - 1];
		{
			*vecTermType++=1;
			k++;
		}
		for (j=1; j <=order; j++) 
		{
			*vecTermType++=2; 
			k++;
		}
		*KE++=k - 1;
	}
	basis->k=k - 1;
}
void evaluate_basis_both(struct BASIS * _restrict  basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st)
{
	int r1,r2;
	int32_t Npad=(int32_t)ceil((float)N/8.0f) * 8;
	r_ippsSet_32f(0,Xt_mars,Npad*basis->k);
	int numOfSeg=basis->sKnotNum+1;
	uint16_t *bkPointsList=basis->S;
	int k=1;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int segLength=r2 - r1+1;
		int TWO_sORDERS=basis->sOrder[(i)-1] * 2;
		for (int j=1; j <=TWO_sORDERS; j++) 
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1;,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			float tmpSum;
			r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
			tmpSum=tmpSum/segLength;
			r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			k++;
#endif
		}
	}
	numOfSeg=basis->tKnotNum+1;
	bkPointsList=basis->T;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int TORDER=basis->tOrder[i - 1]+1;
		int segLength=r2 - r1+1;
		for (int j=1; j <=TORDER; j++)
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			if (numOfSeg !=1||j !=1) normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
				r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			}
			else
			{
				if (numOfSeg !=1) normalize(Xt_mars+(k - 1)*Npad,N);
			}
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				float tmpSum;
				r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
				tmpSum=tmpSum/segLength;
				r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			else
			{
			}
			k++;
#endif
		}
	}
	k--; 
	if (k !=basis->k)
			r_printf("Evaluate_basis_both:The two k's don't match; something wrong!%d%d    \n ",k,basis->k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_zero(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	float *XtX=basis->XtX;
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k,k,N,1,Xt_mars,Npad,Xt_mars,Npad,0,XtX,k);
	float *XtY=basis->XtY;
	r_cblas_sgemv(CblasColMajor,CblasTrans,Npad,k,1,Xt_mars,Npad,pyInfo->Y,1,0,XtY,1);
	float *post_P_U=basis->post_P_U;
	r_cblas_scopy(k*k,XtX,1,post_P_U,1);
	uint8_t *TERMTYPE=basis->termType;
	for (int i=1,ind=0; i <=k; i++)
	{
		post_P_U[ind+(i)-1]+=pmodelPar->prec[*TERMTYPE++];
		ind+=k;
	}
	float *beta_mean=basis->beta_mean;
	float *beta=basis->beta;
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',k,post_P_U,k); 
	r_cblas_scopy(k,XtY,1,beta_mean,1);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',k,1,post_P_U,k,beta_mean,k);
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,k,beta,0,1);
	r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',k,1,post_P_U,k,beta,k);
	r_ippsMulC_32f_I(sqrtf(pmodelPar->sig2),beta,k);
	r_ippsAdd_32f_I(beta_mean,beta,k);
	r_cblas_scopy(k,beta_mean,1,GlobalMEMBuf_1st,1);
	r_cblas_strmv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,k,post_P_U,k,GlobalMEMBuf_1st,1);
	basis->alpha_star=pyInfo->YtY - DOT(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st);
	float half_log_det_post;
	float half_log_det_prior;
	for (int i=1,ind=0; i <=k; i++)
	{
		GlobalMEMBuf_1st[(i)-1]=post_P_U[ind+(i)-1];
		ind+=k;
	}
	r_vmsLn(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st,VML_HA);
	r_ippsSum_32f(GlobalMEMBuf_1st,k,&half_log_det_post,ippAlgHintAccurate);
	half_log_det_post=k* logf(1.0f) - half_log_det_post;
	half_log_det_prior=-.5f*(basis->k_SN * logf(pmodelPar->prec[0])+basis->k_const*logf(pmodelPar->prec[1])+(k - basis->k_SN - basis->k_const)*logf(pmodelPar->prec[2]));
	basis->marg_lik=half_log_det_post - half_log_det_prior - (pyInfo->n *0.5f+pmodelPar->alpha_2) *logf(pmodelPar->alpha_1+basis->alpha_star * 0.5f);
	zeroOut_Xmars_fill(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	return;
}
void evaluate_basis_both_fast(struct BASIS * _restrict  basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st)
{
	int r1,r2;
	int32_t Npad=(int32_t)ceil((float)N/8.0f) * 8;
	r_ippsSet_32f(0,Xt_mars,Npad*basis->k);
	int numOfSeg=basis->sKnotNum+1;
	uint16_t *bkPointsList=basis->S;
	int k=1;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1] ;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1]-1;
		int segLength=r2 - r1+1;
		int TWO_sORDERS=basis->sOrder[(i)-1] * 2;
		for (int j=1; j <=TWO_sORDERS; j++) 
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1;,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			r_ippsMulC_32f_I(sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			float tmpSum;
			r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
			tmpSum=tmpSum/segLength;
			r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			k++;
#endif
		}
	}
	numOfSeg=basis->tKnotNum+1;
	bkPointsList=basis->T;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ?        1 : bkPointsList[(i - 1) - 1];
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1] -1;
		int TORDER=basis->tOrder[i - 1]+1;
		int segLength=r2 - r1+1;
		for (int j=1; j <=TORDER; j++)
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			if (numOfSeg !=1||j !=1) normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			if (j !=1)
			{
				r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
				normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
				r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			}
			else
			{
				r_ippsSet_32f(sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				float tmpSum;
				r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
				tmpSum=tmpSum/segLength;
				r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			else
			{
			}
			k++;
#endif
		}
	}
	k--; 
	if (k !=basis->k)
		r_printf("Evaluate_basis_both:The two k's don't match; something wrong!%d%d    \n ",k,basis->k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_zero(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	float *XtX=basis->XtX;
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k,k,N,1,Xt_mars,Npad,Xt_mars,Npad,0,XtX,k);
	float *XtY=basis->XtY;
	r_cblas_sgemv(CblasColMajor,CblasTrans,Npad,k,1,Xt_mars,Npad,pyInfo->Y,1,0,XtY,1);
	float *post_P_U=basis->post_P_U;
	r_cblas_scopy(k*k,XtX,1,post_P_U,1);
	uint8_t *TERMTYPE=basis->termType;
	for (int i=1,ind=0; i <=k; i++)
	{
		post_P_U[ind+(i)-1]+=pmodelPar->prec[*TERMTYPE++];
		ind+=k;
	}
	float *beta_mean=basis->beta_mean;
	float *beta=basis->beta;
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',k,post_P_U,k); 
	r_cblas_scopy(k,XtY,1,beta_mean,1);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',k,1,post_P_U,k,beta_mean,k);
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,k,beta,0,1);
	r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',k,1,post_P_U,k,beta,k);
	r_ippsMulC_32f_I(sqrtf(pmodelPar->sig2),beta,k);
	r_ippsAdd_32f_I(beta_mean,beta,k);
	r_cblas_scopy(k,beta_mean,1,GlobalMEMBuf_1st,1);
	r_cblas_strmv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,k,post_P_U,k,GlobalMEMBuf_1st,1);
	basis->alpha_star=pyInfo->YtY - DOT(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st);
	float half_log_det_post;
	float half_log_det_prior;
	for (int i=1,ind=0; i <=k; i++)
	{
		GlobalMEMBuf_1st[(i)-1]=post_P_U[ind+(i)-1];
		ind+=k;
	}
	r_vmsLn(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st,VML_HA);
	r_ippsSum_32f(GlobalMEMBuf_1st,k,&half_log_det_post,ippAlgHintAccurate);
	half_log_det_post=k* logf(1.0f) - half_log_det_post;
	half_log_det_prior=-.5f*(basis->k_SN * logf(pmodelPar->prec[0])+basis->k_const*logf(pmodelPar->prec[1])+(k - basis->k_SN - basis->k_const)*logf(pmodelPar->prec[2]));
	basis->marg_lik=half_log_det_post - half_log_det_prior - (pyInfo->n *0.5f+pmodelPar->alpha_2) *logf(pmodelPar->alpha_1+basis->alpha_star * 0.5f);
	zeroOut_Xmars_fill(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	return;
}
void findGoodKnotPositionFromBasis(struct BASIS * _restrict basis,uint8_t * _restrict goodS,uint8_t *_restrict goodT,int32_t N,uint16_t minSepDist_Season,uint16_t minSepDist_Trend)
{
	int posNum,n;
	if (goodS !=NULL)
	{
		r_ippsSet_8u(1,goodS,N);
		posNum=2 * minSepDist_Season+1;
		n=basis->sKnotNum;
		for (int i=1; i <=n; i++)
		{
			uint8_t *tmp=goodS+(basis->S[i - 1] - minSepDist_Season) - 1;
			for (int j=1; j <=posNum; j++)
				*tmp++=0;
		}
		r_ippsSet_8u(0,goodS,(minSepDist_Season+1) );
		r_ippsSet_8u(0,goodS+(N-minSepDist_Season+1) - 1,minSepDist_Season );
	}
	if (goodT !=NULL)
	{
		r_ippsSet_8u(1,goodT,N);
		posNum=2 * minSepDist_Trend+1;
		n=basis->tKnotNum;
		for (int i=1; i <=n; i++)
		{
			uint8_t *tmp=goodT+(basis->T[i - 1] - minSepDist_Trend) - 1;
			for (int j=1; j <=posNum; j++)
				*tmp++=0;
		}
		r_ippsSet_8u(0,goodT,minSepDist_Trend+1 );
		r_ippsSet_8u(0,goodT+(N -minSepDist_Trend+1) - 1,minSepDist_Trend);
	}
}
void allocate_single_output(RESULT * _restrict result,Options * _restrict opt,MemPointers * _restrict MEM)
{
	result->sN=MEM->alloc(MEM,sizeof(float)* 1,0);
	result->tN=MEM->alloc(MEM,sizeof(float)* 1,0);
	result->sNProb=MEM->alloc(MEM,sizeof(int32_t)*(opt->maxKnotNum_Season+1),64);
	result->tNProb=MEM->alloc(MEM,sizeof(int32_t)*(opt->maxKnotNum_Trend+1),64);
	result->sProb=MEM->alloc(MEM,sizeof(int32_t)*opt->N,64);
	result->tProb=MEM->alloc(MEM,sizeof(int32_t)*opt->N,64);
	result->s=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->t=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->b=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->sSD=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->tSD=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->bSD=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->sCI=result->tCI=result->bCI=NULL;
	if (opt->computeCredible)
	{
		result->sCI=MEM->alloc(MEM,sizeof(float)*opt->N * 2,64);
		result->tCI=MEM->alloc(MEM,sizeof(float)*opt->N * 2,64);
		result->bCI=MEM->alloc(MEM,sizeof(float)*opt->N * 2,64);
	}
	result->marg_lik=MEM->alloc(MEM,sizeof(float)* 1,0);
	if (opt->computeSlopeSign)
	{
		result->bsign=MEM->alloc(MEM,sizeof(int32_t)*opt->N,64);
	}
	result->sig2=MEM->alloc(MEM,sizeof(float)* 1,0);
	if (opt->computeHarmonicOrder)
	{
		result->horder=MEM->alloc(MEM,sizeof(uint32_t)*opt->N,64);
	}
	if (opt->computeTrendOrder)
	{
		result->torder=MEM->alloc(MEM,sizeof(uint32_t)*opt->N,64);
	}
}
void print_options(Options * _restrict opt)
{
	char filler=opt->separator;
	r_printf("......Options used in the simulation ......\n");
	r_printf("   opt%cperiod=%d \n",filler,(int) opt->period);
	r_printf("   opt%cstartTime=%f\n",filler,opt->startTime);
	r_printf("   opt%ctimeInterval=%f\n",filler,opt->timeInterval);
	r_printf("   opt%cminSeasonOrder=%d\n",filler,opt->minSeasonOrder);
	r_printf("   opt%cmaxSeasonOrder=%d\n",filler,opt->maxSeasonOrder);
	r_printf("   opt%cminTrendOrder=%d\n",filler,opt->minTrendOrder);
	r_printf("   opt%cmaxTrendOrder=%d\n",filler,opt->maxTrendOrder);
	r_printf("   opt%cminSepDist_Trend=%d\n",filler,opt->minSepDist_Trend);
	r_printf("   opt%cminSepDist_Season=%d\n",filler,opt->minSepDist_Season);
	r_printf("   opt%cmaxKnotNum_Trend=%d\n",filler,opt->maxKnotNum_Trend);
	r_printf("   opt%cmaxKnotNum_Season=%d\n",filler,opt->maxKnotNum_Season);
	r_printf("   opt%cmaxMoveStepSize=%d\n",filler,opt->maxMoveStepSize);
	r_printf("   opt%csamples=%d\n",filler,opt->samples);
	r_printf("   opt%cthinningFactor=%d\n",filler,opt->thinningFactor);
	r_printf("   opt%cburnin=%d\n",filler,opt->burnin);
	r_printf("   opt%cchainNumber=%d\n",filler,opt->chainNumber);
	r_printf("   opt%cresamplingTrendOrderProb=%f\n",filler,opt->resamplingTrendOrderProb);
	r_printf("   opt%cresamplingSeasonOrderProb=%f\n",filler,opt->resamplingSeasonOrderProb);
	r_printf("   opt%comissionValue=%f\n",filler,opt->omissionValue);
	r_printf("   opt%cseed=%d\n",filler,opt->seed);
	r_printf("   opt%coutputToDisk=%d\n",filler,opt->outputToDisk);
	if (opt->outputToDisk)
		r_printf("   opt%coutputFolder=%s\n",filler,opt->outputFolder);
	else
		r_printf("   opt%coutputFolder=%s\n",filler,"NOT USED");
	r_printf("   opt%clengthPerTimeSeries_infile=%d\n",filler,opt->N);
	r_printf("   opt%cprintToScreen=%d\n",filler,opt->printToScreen);
	r_printf("   opt%cprintCharLen=%d\n",filler,opt->printCharLen);
	r_printf("   opt%ccomputeCredible=%d\n",filler,opt->computeCredible);
	r_printf("   opt%cfastCIComputation=%d\n",filler,opt->fastCIComputation);
	r_printf("   opt%ccomputeChangepoints=%d\n",filler,opt->computeChangepoints);
	r_printf("   opt%ccomputeSlopeSign=%d\n",filler,opt->computeSlopeSign);
	r_printf("   opt%ccomputeHarmonicOrder=%d\n",filler,opt->computeHarmonicOrder);
	r_printf("   opt%ccomputeTrendOrder=%d\n",filler,opt->computeTrendOrder);
	r_printf("......End of displaying Options ......\n\n");
}
void evaluate_basis_both_BIC(struct BASIS * _restrict basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st)
{
	int r1,r2;
	int32_t Npad=(int32_t)ceil((float)N/8.0f) * 8;
	r_ippsSet_32f(0,Xt_mars,Npad*basis->k);
	int numOfSeg=basis->sKnotNum+1;
	uint16_t *bkPointsList=basis->S;
	int k=1;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int segLength=r2 - r1+1;
		int TWO_sORDERS=basis->sOrder[(i)-1] * 2;
		for (int j=1; j <=TWO_sORDERS; j++) 
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1;,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,SEASON_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			float tmpSum;
			r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
			tmpSum=tmpSum/segLength;
			r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			k++;
#endif
		}
	}
	numOfSeg=basis->tKnotNum+1;
	bkPointsList=basis->T;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int TORDER=basis->tOrder[i - 1]+1;
		int segLength=r2 - r1+1;
		for (int j=1; j <=TORDER; j++)
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			if (numOfSeg !=1||j !=1) normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
				r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			}
			else
			{
				if (numOfSeg !=1) normalize(Xt_mars+(k - 1)*Npad,N);
			}
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				float tmpSum;
				r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
				tmpSum=tmpSum/segLength;
				r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			else
			{
			}
			k++;
#endif
		}
	}
	k--; 
	if (k !=basis->k)
		r_printf("Evaluate_basis_both:The two k's don't match; something wrong!%d%d    \n ",k,basis->k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_zero(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	float *XtX=basis->XtX;
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k,k,N,1,Xt_mars,Npad,Xt_mars,Npad,0,XtX,k);
	float *XtY=basis->XtY;
	r_cblas_sgemv(CblasColMajor,CblasTrans,Npad,k,1,Xt_mars,Npad,pyInfo->Y,1,0,XtY,1);
	float *post_P_U=basis->post_P_U;
	r_cblas_scopy(k*k,XtX,1,post_P_U,1);
	float *beta_mean=basis->beta_mean;
	float *beta=basis->beta;
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',k,post_P_U,k); 
	r_cblas_scopy(k,XtY,1,beta_mean,1);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',k,1,post_P_U,k,beta_mean,k);
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,k,beta,0,1);
	r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',k,1,post_P_U,k,beta,k);
	r_ippsMulC_32f_I(sqrtf(pmodelPar->sig2),beta,k);
	r_ippsAdd_32f_I(beta_mean,beta,k);
	r_cblas_ssymv(CblasColMajor,CblasUpper,k,1.0,XtX,k,beta_mean,1,0,GlobalMEMBuf_1st,1);
	float RSS=(pyInfo->YtY - 2 * DOT(k,beta_mean,XtY)+DOT(k,beta_mean,GlobalMEMBuf_1st));
	basis->marg_lik=pyInfo->n *logf(RSS )+(k+1)*logf((float)pyInfo->n);
	basis->marg_lik=basis->marg_lik *(-0.5f);
	basis->alpha_star=RSS/(pyInfo->n - k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_fill(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	return;
}
void zero_out_result_output(Options * _restrict opt,RESULT *_restrict result)
{
	int32_t  N=opt->N;
    r_ippsSet_32s(0,(I32PTR) result->sN,1);
	r_ippsSet_32s(0,(I32PTR) result->tN,1);
	r_ippsSet_32s(0,result->sProb,N);
	r_ippsSet_32s(0,result->tProb,N);
	r_ippsSet_32f(0,result->s,N);
	r_ippsSet_32f(0,result->t,N);
	r_ippsSet_32f(0,result->b,N);
	if (opt->computeCredible)
		r_ippsSet_32f(0,result->sCI,N+N),
		r_ippsSet_32f(0,result->tCI,N+N),
		r_ippsSet_32f(0,result->bCI,N+N);
	r_ippsSet_32f(0,result->sSD,N);
	r_ippsSet_32f(0,result->tSD,N);
	r_ippsSet_32f(0,result->bSD,N);
	r_ippsSet_32s(0,result->sNProb,opt->maxKnotNum_Season+1);
	r_ippsSet_32s(0,result->tNProb,opt->maxKnotNum_Trend+1);
	*result->marg_lik=0.f;
	*result->sig2=0.f;
	if (opt->computeSlopeSign)
		r_ippsSet_32s(0,result->bsign,N);
	if (opt->computeHarmonicOrder)
		r_ippsSet_32s(0,result->horder,N);
	if (opt->computeTrendOrder)
		r_ippsSet_32s(0,result->torder,N);	 
}
void nan_fill_result_output(Options * _restrict opt,RESULT *_restrict result)
{
	int32_t  N=opt->N;
	rF32    nan=(1e300*1e300)*0.f;
	*result->sN=nan;
	*result->tN=nan;
	*result->marg_lik=nan;
	*result->sig2=nan;
	for (rI32 i=0; i < N; i++)
	{
		*(result->sProb+i)=nan;
		*(result->tProb+i)=nan;
		*(result->s+i)=nan;
		*(result->t+i)=nan;
		*(result->b+i)=nan;
		*(result->sSD+i)=nan;
		*(result->tSD+i)=nan;
		*(result->bSD+i)=nan;
	}
	if (opt->computeCredible)
	{
		for (rI32 i=0; i < N; i++)
		{
			*(result->sCI+i)=nan;
			*(result->tCI+i)=nan;
			*(result->bCI+i)=nan;
			*(result->sCI+N+i)=nan;
			*(result->tCI+N+i)=nan;
			*(result->bCI+N+i)=nan;
		}
	}
	for (rI32 i=0; i < (opt->maxKnotNum_Season+1); i++)
		*(result->sNProb+i)=nan;
	for (rI32 i=0; i < (opt->maxKnotNum_Trend+1); i++)
		*(result->tNProb+i)=nan;
	if (opt->computeSlopeSign)
	{
		for (rI32 i=0; i < N; i++)
			*(result->bsign+i)=nan;	
	}
	if (opt->computeHarmonicOrder)
	{
		for (rI32 i=0; i < N; i++)
			*(result->horder+i)=nan;
	}
	if (opt->computeTrendOrder)
	{
		for (rI32 i=0; i < N; i++)
			*(result->torder+i)=nan;
	}
}
void convert_basis_both_trend(struct BASIS * _restrict basis)
{
	basis->k_SN=0;
	int numOfSeg;
	int k;
	int order;
	int i,j;
	register int16_t *  _restrict KS,*_restrict KE;  
	register uint8_t * _restrict vecTermType;
	vecTermType=basis->termType;
 	uint8_t * _restrict vecOrder;
	numOfSeg=basis->tKnotNum+1; 
	basis->k_const=numOfSeg;
	KS=basis->tks;
	KE=basis->tke;
	k=1; 
	vecOrder=basis->tOrder;
	for (i=1; i <=numOfSeg; i++)
	{
		*KS++=k;
		order=vecOrder[i - 1];
		{
			*vecTermType++=1;
			k++;
		}
		for (j=1; j <=order; j++) 
		{
			*vecTermType++=2; 
			k++;
		}
		*KE++=k - 1;
	}
	basis->k=k - 1;
}
void evaluate_basis_both_trend(struct BASIS * _restrict  basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st)
{
	int r1,r2;
	int32_t Npad=(int32_t)ceil((float)N/8.0f) * 8;
	r_ippsSet_32f(0,Xt_mars,Npad*basis->k);
	int		 numOfSeg=basis->tKnotNum+1;
	uint16_t *bkPointsList=basis->T;
	int		 k=1;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int TORDER=basis->tOrder[i - 1]+1;
		int segLength=r2 - r1+1;
		for (int j=1; j <=TORDER; j++)
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			if (numOfSeg !=1||j !=1) normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
				r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			}
			else
			{
				if (numOfSeg !=1) normalize(Xt_mars+(k - 1)*Npad,N);
			}
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				float tmpSum;
				r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
				tmpSum=tmpSum/segLength;
				r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			else
			{
			}
			k++;
#endif
		}
	}
	k--; 
	if (k !=basis->k)
		r_printf("Evaluate_basis_both:The two k's don't match; something wrong!%d%d    \n ",k,basis->k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_zero(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	float *XtX=basis->XtX;
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k,k,N,1,Xt_mars,Npad,Xt_mars,Npad,0,XtX,k);
	float *XtY=basis->XtY;
	r_cblas_sgemv(CblasColMajor,CblasTrans,Npad,k,1,Xt_mars,Npad,pyInfo->Y,1,0,XtY,1);
	float *post_P_U=basis->post_P_U;
	r_cblas_scopy(k*k,XtX,1,post_P_U,1);
	uint8_t *TERMTYPE=basis->termType;
	for (int i=1,ind=0; i <=k; i++)
	{
		post_P_U[ind+(i)-1]+=pmodelPar->prec[*TERMTYPE++];
		ind+=k;
	}
	float *beta_mean=basis->beta_mean;
	float *beta=basis->beta;
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',k,post_P_U,k); 
	r_cblas_scopy(k,XtY,1,beta_mean,1);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',k,1,post_P_U,k,beta_mean,k);
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,k,beta,0,1);
	r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',k,1,post_P_U,k,beta,k);
	r_ippsMulC_32f_I(sqrtf(pmodelPar->sig2),beta,k);
	r_ippsAdd_32f_I(beta_mean,beta,k);
	r_cblas_scopy(k,beta_mean,1,GlobalMEMBuf_1st,1);
	r_cblas_strmv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,k,post_P_U,k,GlobalMEMBuf_1st,1);
	basis->alpha_star=pyInfo->YtY - DOT(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st);
	float half_log_det_post;
	float half_log_det_prior;
	for (int i=1,ind=0; i <=k; i++)
	{
		GlobalMEMBuf_1st[(i)-1]=post_P_U[ind+(i)-1];
		ind+=k;
	}
	r_vmsLn(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st,VML_HA);
	r_ippsSum_32f(GlobalMEMBuf_1st,k,&half_log_det_post,ippAlgHintAccurate);
	half_log_det_post=k* logf(1.0f) - half_log_det_post;
	half_log_det_prior=-.5f*(basis->k_SN * logf(pmodelPar->prec[0])+basis->k_const*logf(pmodelPar->prec[1])+(k - basis->k_SN - basis->k_const)*logf(pmodelPar->prec[2]));
	basis->marg_lik=half_log_det_post - half_log_det_prior - (pyInfo->n *0.5f+pmodelPar->alpha_2) *logf(pmodelPar->alpha_1+basis->alpha_star * 0.5f);
	zeroOut_Xmars_fill(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	return;
}
void evaluate_basis_both_trend_fast(struct BASIS * _restrict  basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st)
{
	int r1,r2;
	int32_t Npad=(int32_t)ceil((float)N/8.0f) * 8;
	r_ippsSet_32f(0,Xt_mars,Npad*basis->k);
	int       numOfSeg=basis->tKnotNum+1;
	uint16_t *bkPointsList=basis->T;
	int       k=1;
	for (int i=1; i <=numOfSeg; i++)
	{
		r1=(i==1) ? 1 : bkPointsList[(i - 1) - 1]+1;
		r2=(i==numOfSeg) ? N : bkPointsList[(i)-1];
		int TORDER=basis->tOrder[i - 1]+1;
		int segLength=r2 - r1+1;
		for (int j=1; j <=TORDER; j++)
		{
#if BASIS_METHODS==1
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*N+r1 - 1,1);
			if (numOfSeg !=1||j !=1) normalize(Xt_mars+(k - 1)*N,N);
			k++;
#elif BASIS_METHODS==2
			if (j !=1)
			{
				r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
				normalize(Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
				r_cblas_sscal(segLength,sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,1);
			}
			else
			{
				r_ippsSet_32f(sqrtf(N/(segLength+0.f)),Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			k++;
#elif BASIS_METHODS==3
			r_cblas_scopy(segLength,TREND_TERMS+(j - 1)*N+r1 - 1,1,Xt_mars+(k - 1)*Npad+r1 - 1,1);
			if (j !=1)
			{
				float tmpSum;
				r_ippsSum_32f(Xt_mars+(k - 1)*Npad+r1 - 1,segLength,&tmpSum,ippAlgHintAccurate);
				tmpSum=tmpSum/segLength;
				r_ippsSubC_32f_I(tmpSum,Xt_mars+(k - 1)*Npad+r1 - 1,segLength);
			}
			else
			{
			}
			k++;
#endif
		}
	}
	k--; 
	if (k !=basis->k)
		r_printf("Evaluate_basis_both:The two k's don't match; something wrong!%d%d    \n ",k,basis->k);
	if (pyInfo->nMissing>0)
		zeroOut_Xmars_zero(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	float *XtX=basis->XtX;
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,k,k,N,1,Xt_mars,Npad,Xt_mars,Npad,0,XtX,k);
	float *XtY=basis->XtY;
	r_cblas_sgemv(CblasColMajor,CblasTrans,Npad,k,1,Xt_mars,Npad,pyInfo->Y,1,0,XtY,1);
	float *post_P_U=basis->post_P_U;
	r_cblas_scopy(k*k,XtX,1,post_P_U,1);
	uint8_t *TERMTYPE=basis->termType;
	for (int i=1,ind=0; i <=k; i++)
	{
		post_P_U[ind+(i)-1]+=pmodelPar->prec[*TERMTYPE++];
		ind+=k;
	}
	float *beta_mean=basis->beta_mean;
	float *beta=basis->beta;
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',k,post_P_U,k); 
	r_cblas_scopy(k,XtY,1,beta_mean,1);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',k,1,post_P_U,k,beta_mean,k);
	r_vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,k,beta,0,1);
	r_LAPACKE_strtrs(LAPACK_COL_MAJOR,'U','N','N',k,1,post_P_U,k,beta,k);
	r_ippsMulC_32f_I(sqrtf(pmodelPar->sig2),beta,k);
	r_ippsAdd_32f_I(beta_mean,beta,k);
	r_cblas_scopy(k,beta_mean,1,GlobalMEMBuf_1st,1);
	r_cblas_strmv(CblasColMajor,CblasUpper,CblasNoTrans,CblasNonUnit,k,post_P_U,k,GlobalMEMBuf_1st,1);
	basis->alpha_star=pyInfo->YtY - DOT(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st);
	float half_log_det_post;
	float half_log_det_prior;
	for (int i=1,ind=0; i <=k; i++)
	{
		GlobalMEMBuf_1st[(i)-1]=post_P_U[ind+(i)-1];
		ind+=k;
	}
	r_vmsLn(k,GlobalMEMBuf_1st,GlobalMEMBuf_1st,VML_HA);
	r_ippsSum_32f(GlobalMEMBuf_1st,k,&half_log_det_post,ippAlgHintAccurate);
	half_log_det_post=k* logf(1.0f) - half_log_det_post;
	half_log_det_prior=-.5f*(basis->k_SN * logf(pmodelPar->prec[0])+basis->k_const*logf(pmodelPar->prec[1])+(k - basis->k_SN - basis->k_const)*logf(pmodelPar->prec[2]));
	basis->marg_lik=half_log_det_post - half_log_det_prior - (pyInfo->n *0.5f+pmodelPar->alpha_2) *logf(pmodelPar->alpha_1+basis->alpha_star * 0.5f);
	zeroOut_Xmars_fill(Xt_mars,Xt_zeroBackup,pyInfo->rowsMissing,pyInfo->nMissing,N,Npad,k);
	return;
}
void print_options_trend(Options * _restrict opt )
{
	char filler=opt->separator;
	r_printf("......Options used in the simulation ......n");
	r_printf("   opt%cstartTime=%f\n",filler,opt->startTime);
	r_printf("   opt%ctimeInterval=%f\n",filler,opt->timeInterval);
	r_printf("   opt%cminTrendOrder=%d\n",filler,opt->minTrendOrder);
	r_printf("   opt%cmaxTrendOrder=%d\n",filler,opt->maxTrendOrder);
	r_printf("   opt%cminSepDist_Trend=%d\n",filler,opt->minSepDist_Trend);
	r_printf("   opt%cmaxKnotNum_Trend=%d\n",filler,opt->maxKnotNum_Trend);
	r_printf("   opt%cmaxMoveStepSize=%d\n",filler,opt->maxMoveStepSize);
	r_printf("   opt%csamples=%d\n",filler,opt->samples);
	r_printf("   opt%cthinningFactor=%d\n",filler,opt->thinningFactor);
	r_printf("   opt%cburnin=%d\n",filler,opt->burnin);
	r_printf("   opt%cchainNumber=%d\n",filler,opt->chainNumber);
	r_printf("   opt%cresamplingTrendOrderProb=%f\n",filler,opt->resamplingTrendOrderProb);
	r_printf("   opt%comissionValue=%f\n",filler,opt->omissionValue);
	r_printf("   opt%cseed=%d\n",filler,opt->seed);
	r_printf("   opt%coutputToDisk=%d\n",filler,opt->outputToDisk);
	if (opt->outputToDisk)
		r_printf("   opt%coutputFolder=%s\n",filler,opt->outputFolder);
	else
		r_printf("   opt%coutputFolder=%s\n",filler,"NOT USED");
	r_printf("   opt%clengthPerTimeSeries_infile=%d\n",filler,opt->N);
	r_printf("   opt%cprintToScreen=%d\n",filler,opt->printToScreen);
	r_printf("   opt%cprintCharLen=%d\n",filler,opt->printCharLen);
	r_printf("   opt%ccomputeChangepoints=%d\n",filler,opt->computeChangepoints);
	r_printf("   opt%ccomputeCredible=%d\n",filler,opt->computeCredible);
	r_printf("   opt%cfastCIComputation=%d\n",filler,opt->fastCIComputation);
	r_printf("   opt%ccomputeSlopeSign=%d\n",filler,opt->computeSlopeSign);
	r_printf("   opt%ccomputeTrendOrder=%d\n",filler,opt->computeTrendOrder);
	r_printf("......End of displaying Options ......\n\n");
}
void allocate_single_output_trend(RESULT * _restrict result,Options * _restrict opt,MemPointers * _restrict MEM)
{
	result->tN=MEM->alloc(MEM,sizeof(float)* 1,0);	
	result->tNProb=MEM->alloc(MEM,sizeof(int32_t)*(opt->maxKnotNum_Trend+1),64);	
	result->tProb=MEM->alloc(MEM,sizeof(int32_t)*opt->N,64);
	result->t=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->b=MEM->alloc(MEM,sizeof(float)*opt->N,64);	
	result->tSD=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->bSD=MEM->alloc(MEM,sizeof(float)*opt->N,64);
	result->tCI=result->bCI=(F32PTR) (result->bsign=result->torder=NULL);
	if (opt->computeCredible)
	{
		result->tCI=MEM->alloc(MEM,sizeof(float)*opt->N * 2,64);
		result->bCI=MEM->alloc(MEM,sizeof(float)*opt->N * 2,64);
	}
	if (opt->computeSlopeSign)
	{
		result->bsign=MEM->alloc(MEM,sizeof(int32_t)*opt->N,64);
	}
	if (opt->computeTrendOrder)
	{
		result->torder=MEM->alloc(MEM,sizeof(uint32_t)*opt->N,64);
	}
	result->marg_lik=MEM->alloc(MEM,sizeof(float)* 1,0);
	result->sig2=MEM->alloc(MEM,sizeof(float)* 1,0);
}
void zero_out_result_output_trend(Options * _restrict opt,RESULT *_restrict result)
{
	int32_t  N=opt->N;	
	r_ippsSet_32s(0,(I32PTR) result->tN,1);
	r_ippsSet_32s(0,(I32PTR) result->tProb,N);
	r_ippsSet_32f(0,result->t,N);
	r_ippsSet_32f(0,result->b,N);
	if (opt->computeCredible)		
		r_ippsSet_32f(0,result->tCI,N+N),
		r_ippsSet_32f(0,result->bCI,N+N);	
	r_ippsSet_32f(0,result->tSD,N);
	r_ippsSet_32f(0,result->bSD,N);
	r_ippsSet_32s(0,result->tNProb,opt->maxKnotNum_Trend+1);
	*result->marg_lik=0.f;
	*result->sig2=0.f;
	if (opt->computeSlopeSign)
		r_ippsSet_32s(0,result->bsign,N);
	if (opt->computeTrendOrder)
		r_ippsSet_32s(0,result->torder,N);
}
ENABLE_MANY_WARNINGS
