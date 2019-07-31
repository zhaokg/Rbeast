#include "abc_common.h"
#include "math.h"
#include "string.h"
void transpose_inplace(rF32PTR m,int32_t w,int32_t h ) 
{
	int32_t start,next,i;
	float tmp;
	int32_t totalElement=w * h - 1;
	for (start=0; start <=totalElement; start++) {
		next=start;
		i=0;
		do {
			i++;
			next=(next%h) * w+next/h;
		} while (next > start);
		if (next < start||i==1) continue;
		tmp=m[next=start];
		do {
			i=(next%h) * w+next/h;
			m[next]=(i==start) ? tmp : m[i];
			next=i;
		} while (next > start);
	}
}
void  fill_1_to_N(rF32PTR p,int N)
{
	for (rI32 i=1; i <=N; i++) *p++=(float)i;
	p -=N;
}
int32_t strcicmp(char const * _restrict a,char const * _restrict b)
{
	for (;; a++,b++) {
		rI32 d=((*a)|(uint8_t)32) - ((*b)|(uint8_t)32);
		if (d !=0||!*a)
			return d;
	}
}
float determine_period(rF32PTR data,int32_t N,float  omissionValue)
{
	rF32PTR xData=(float *)r_malloc(sizeof(float)*N * 4);
	float delta=1.f/N;
	float curValue=0.f;
	for (register int32_t i=0; i < N; i++)
	{
		float tmp;
		*(xData+i)=1.0;
		*(xData+N+i)=curValue;
		tmp=curValue;
		tmp=tmp*curValue;
		*(xData+N+N+i)=tmp;
		tmp=tmp*curValue;
		*(xData+N+N+N+i)=tmp;
		curValue=curValue+delta;
	}
	uint8_t * _restrict isNA=(uint8_t * ) r_malloc(sizeof(char)*N);
	memset(isNA,0,sizeof(char)*N);
	rF32PTR y=data;
	for (register int32_t i=0; i < N; i++,y++)
	{
		if (*y !=*y||fabs(*y - omissionValue) < 1e-5)
		{
			isNA[i]=1;
			*(xData+i)=0.f;
			*(xData+N+i)=0.f;
			*(xData+N+N+i)=0.f;
			*(xData+N+N+N+i)=0.f;
			*y=0.f;
		}
	}
	rF32PTR XtX=(float *)r_malloc(sizeof(float)* 4 * 4);
	rF32PTR XtY=(float *)r_malloc(sizeof(float)* 4);
	r_cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,4,4,N,1.0f,xData,N,xData,N,0,XtX,4);
	r_LAPACKE_spotrf(LAPACK_COL_MAJOR,'U',4,XtX,4);
	r_cblas_sgemv(CblasColMajor,CblasTrans,N,4,1.f,xData,N,data,1L,0.f,XtY,1L);
	r_LAPACKE_spotrs(LAPACK_COL_MAJOR,'U',4,1,XtX,4,XtY,4);
	rF32PTR yFit=(float *)r_malloc(sizeof(float)* N);
	r_cblas_sgemv(CblasColMajor,CblasNoTrans,N,4,1.f,xData,N,XtY,1L,0.f,yFit,1L);
	r_ippsSub_32f_I(yFit,data,N);
	r_free(yFit);
	r_free(XtY);
	r_free(XtX);
	r_free(xData);
	int32_t M=(int)ceil(N/2);
	rF32PTR ans=(float *)r_malloc(sizeof(float)*M);
	for (register int32_t i=1; i <=M; i++)
	{
		int32_t len=N - i;
		int32_t start=i+1;
		float XY=0,XX=0,YY=0;
		float MX=0,MY=0;
		int32_t NUM=0;
		for (int32_t j=1; j <=len; j++)
		{
			int32_t Ix=j - 1;
			int32_t Iy=start+(j - 1) - 1;
			if ((isNA[Ix]+isNA[Iy])==0)
			{
				NUM++;
				float x=data[Ix];
				float y=data[Iy];
				MX=MX+x;				MY=MY+y;
				XY=XY+x*y; 				XX=XX+x*x;
				YY=YY+y*y;
			}
		}
		MX=MX/(float)NUM;
		MY=MY/(float)NUM;
		ans[i - 1]=(XY/NUM - MX*MY)/sqrtf((XX/N - MX*MX)*(YY/N - MY*MY));
	}
	memset(isNA,0,M);
	I32PTR index=(int32_t *)r_malloc(sizeof(int32_t)*M);
	size_t  totalLMnum=0;
	for (register int32_t i=2; i <=(M - 1); i++)
	{
		if (ans[(i)-1] > ans[(i - 1) - 1] && ans[(i)-1] > ans[(i+1) - 1])
		{
			isNA[i - 1]=1;
			index[totalLMnum++]=i;
		}
	}
	int32_t period=0;
	if (totalLMnum !=0)
	{
		for (int32_t curIdx=1; curIdx <=max(1,(int)floorf(totalLMnum/3.f)); curIdx++)
		{
			period=index[curIdx - 1];
			int32_t  goodTimes=0;
			int32_t  numOfPeriod=(int)floorf((float)(M - 1)/period);
			for (int32_t i=1; i <=numOfPeriod; i++)
			{
				if (isNA[period*i - 1]==1||isNA[period*i+1 - 1]==1||isNA[period*i - 1 - 1]==1)
				{
					goodTimes++;
				}
			}
			if (goodTimes >=min(3,numOfPeriod))
			{
				break;
			}
			period=0;
		}
	}
	r_free(index);
	r_free(ans);
	r_free(isNA);
	return (float)period;
}
static float confidenceInterval(rF32PTR half,int32_t n,char leftOrRight);
int32_t find_changepoint1(rF32PTR prob,rF32PTR mem,float threshold,I32PTR cpt,rF32PTR cptCI,int32_t N,int32_t minSepDist,int32_t maxCptNumber)
{
	if (maxCptNumber==0)	{return maxCptNumber;}
	int32_t w=(int32_t) round((minSepDist - 1)/2);
	int32_t w2=w * 2+1;
	r_ippsSet_32f(0,mem,N);			
	I32PTR cpfromSum_Pos=(int32_t*) mem+N;
	F32PTR cpfromSum_Val=(float *)mem+N * 2;
	I32PTR cpfromProb_Pos=(int32_t *) mem+N * 3;
	F32PTR cpfromProb_Val=(float *)mem+N * 4;
	for (int32_t i=-w; i <=w; i++)
	{
		int32_t len=i > 0 ? i : -i;
		int32_t startIdx_mem=i <=0 ? 0 : i;
		int32_t startIdx_prob=i <=0 ? -i : 0;
		r_ippsAdd_32f_I(prob+startIdx_prob,mem+startIdx_mem,N - len);
	}
	int32_t  UPPERIDX=N - w;
	int32_t  cptNumber=0;
	for (int32_t i=w; i < UPPERIDX; i++)
	{
		if (mem[i] < threshold) continue;
		bool isLargeThanNearestNeighor=(mem[i] >=mem[i - 1]) && (mem[i] >=mem[i+1]);
		bool isLargeThanNearestTwoNeighors=(mem[i] * 4.0) > (mem[i+1]+mem[i+2]+mem[i - 1]+mem[i - 2]);
		if (!(isLargeThanNearestNeighor && isLargeThanNearestTwoNeighors)) continue;
		int32_t		upperIdx_1=i+w;
		int32_t		maxIdx=-999;
		float		maxVal=-999;
		for (int32_t j=i - w; j <=upperIdx_1; j++)
		{
			if ((prob[j] > prob[j - 1] && prob[j] >=prob[j+1])||(prob[j] >=prob[j - 1] && prob[j] > prob[j+1]))
			{
				if (prob[j] > maxVal) 	maxIdx=j,maxVal=prob[j];
			}			
		}		
		if ( maxVal < 0.f	)	continue;
		int32_t diff_btw_twoNeighbors=maxIdx-cpfromProb_Pos[cptNumber - 1]; 
		if ((cptNumber==0)||diff_btw_twoNeighbors >=w2||diff_btw_twoNeighbors <=-w2)
		{
			cpfromSum_Pos[cptNumber]=i;
			cpfromSum_Val[cptNumber]=mem[i];
			cpfromProb_Pos[cptNumber]=maxIdx;
			cpfromProb_Val[cptNumber]=maxVal;
			cptNumber++;
			continue;
		}
		else
		{
			if (maxVal >=cpfromProb_Val[cptNumber - 1])
			{
				cpfromSum_Pos[cptNumber - 1]=i;
				cpfromSum_Val[cptNumber - 1]=mem[i];
				cpfromProb_Pos[cptNumber - 1]=maxIdx;
				cpfromProb_Val[cptNumber - 1]=maxVal;
				continue;
			}
		}		
	}
	if (cptNumber==0) { return cptNumber; }
	quickSortD(cpfromProb_Val,cpfromProb_Pos,0,cptNumber - 1);
	cptNumber=min(cptNumber,maxCptNumber);
	r_cblas_scopy(cptNumber,(float *)cpfromProb_Pos,1,(float *) cpt,1);
	I32PTR index=(int32_t *) mem;
	float *cpt_float=mem+N;
	for (int32_t i=0; i < cptNumber; i++)
	{
		*index++=i;
		*cpt_float++=(float)cpt[i];
	}
	index=index - cptNumber;
	cpt_float=cpt_float - cptNumber;	
	quickSortA(cpt_float,index,0,cptNumber - 1);
	for (int32_t i=0; i < cptNumber; i++)
	{
		cptCI[i]=-9999.f;
		cptCI[cptNumber+i]=-9999.f;
	}
	float delta;
	delta=confidenceInterval(prob,((int32_t) cpt_float[0]-0+1),'L');
	cptCI[0]=delta;
	delta=confidenceInterval(prob+(int32_t)cpt_float[cptNumber - 1],(N - (int32_t)cpt_float[cptNumber - 1]+1),'R');
	cptCI[cptNumber+cptNumber - 1]=delta;
	if (cptNumber==1) {
		cptCI[0]=cpt_float[0] - cptCI[0];
		cptCI[1]=cpt_float[0]+cptCI[1];
		return cptNumber; 
	}
	for (int32_t i=0; i < (cptNumber-1); i++)
	{ 
		float del1,del2,del;
		del1=cptCI[cptNumber+i] > 0 ? cptCI[cptNumber+i] : cptCI[i];
		del2=cptCI[i+1] > 0 ? cptCI[i+1] : ((cptCI[cptNumber+i+1] > 0) ? cptCI[cptNumber+i+1] : -9999.f);
    	del=cpt_float[i+1] - cpt_float[i];
		if (del2 <=0)
		{
				del1=del1 * 2.f;
				del=(del1 > del) ? del/2 : del1;
		}else
		{
			del=del*del1/(del1+del2);
		}
		delta=confidenceInterval(prob+(int32_t)cpt_float[i],(int32_t) ceil(del),'R');
		cptCI[cptNumber+i]=delta;
		del=cpt_float[i+1] - cpt_float[i];
		if (del2 <=0)
		{
			delta=del - delta * 2;
			del=delta <=0 ? del/2 : delta;
		}
		else
		{
			del2=del2 * 2.f;
			del=(del2 >=del) ? del/2 : del2;
		}
		int32_t len=(int32_t)ceil(del);
		delta=confidenceInterval(prob+(int32_t)cpt_float[i+1]-(len-1),len,'L');
		cptCI[i+1]=delta;
	}
	float *temp=mem+2 * N;
	r_cblas_scopy(2*cptNumber,cptCI,1,temp,1);
	for (int32_t i=0; i < cptNumber; i++)
	{
		int32_t idx  ;
		idx=index[i];
		cptCI[idx]=cpt_float[i] - temp[i];
		cptCI[cptNumber+idx]=cpt_float[i]+temp[cptNumber+i];
	}
	return cptNumber;
}
static float confidenceInterval( rF32PTR half,int32_t n,char leftOrRight)
{
	float delta;
	float sum,inv_sum,cumSum;
	r_ippsSum_32f(      half,n,&sum,ippAlgHintAccurate);
	inv_sum=1.f/sum;
	int32_t j;
	cumSum=0;	
	if (leftOrRight=='R')
	{
		for (j=0; j < n; j++)
		{
			cumSum=cumSum+half[j];
			if (cumSum *inv_sum >=0.95) break;
		}
		float J=j+1.f;
		delta=J - (cumSum - 0.95f*sum)/half[j];	
	}
	else
	{
		for (j=n-1; j >=0; j--)
		{
			cumSum=cumSum+half[j];
			if (cumSum*inv_sum  >=0.95) break;
		}
		float J=(float)(n-j);
		delta=J - (cumSum - 0.95f*sum)/half[j];
	}
	return delta;
}
int32_t find_changepoint(rF32PTR prob,rF32PTR mem,float threshold,I32PTR cpt,rF32PTR cptCI,int32_t N,int32_t minSepDist,int32_t maxCptNumber)
{
	if (maxCptNumber==0)	{ return maxCptNumber; }
	int32_t w0=minSepDist/2;
	int32_t w1=minSepDist - w0;
	r_ippsSet_32f(0,mem,N);
	I32PTR cpfromSum_Pos=(int32_t *) mem+N;
	F32PTR cpfromSum_Val=(float *) mem+N * 2;
	I32PTR cpfromProb_Pos=(int32_t *) mem+N * 3;
	F32PTR cpfromProb_Val=(float *)mem+N * 4;
	for (int32_t i=-w1; i <=w0; i++)
	{
		int32_t len=i > 0 ? i : -i;
		int32_t startIdx_mem=i <=0 ? 0 : i;
		int32_t startIdx_prob=i <=0 ? -i : 0;
		r_ippsAdd_32f_I(prob+startIdx_prob,mem+startIdx_mem,N - len);
	}	
	int32_t  UPPERIDX=N - (minSepDist+1);
	int32_t  cptNumber=0;
	for (int32_t i=(minSepDist+1); i < UPPERIDX; i++)
	{
		if (mem[i] < threshold) continue;
		bool isLargeThanNearestNeighor=(mem[i] >=mem[i - 1]) && (mem[i] >=mem[i+1]);
		bool isLargeThanNearestTwoNeighors=(mem[i] * 4.0) > (mem[i+1]+mem[i+2]+mem[i - 1]+mem[i - 2]);
		if (!(isLargeThanNearestNeighor && isLargeThanNearestTwoNeighors)) continue;
		int32_t		UPPERIDX_1=i+w1;
		int32_t		maxIdx=-999;
		float		maxVal=-999;
		for (int32_t j=i - w0; j <=UPPERIDX_1; j++)
		{
			if ((prob[j] > prob[j - 1] && prob[j] >=prob[j+1])||(prob[j] >=prob[j - 1] && prob[j] > prob[j+1]))
			{
				if (prob[j] > maxVal) 	maxIdx=j,maxVal=prob[j];
			}
		}
		if (maxVal < 0.f)	continue;
		int32_t diff_btw_twoNeighbors=maxIdx - cpfromProb_Pos[cptNumber - 1];
		if ((cptNumber==0)||diff_btw_twoNeighbors > minSepDist||diff_btw_twoNeighbors < -minSepDist)
		{
			cpfromSum_Pos[cptNumber]=i;
			cpfromSum_Val[cptNumber]=mem[i];
			cpfromProb_Pos[cptNumber]=maxIdx;
			cpfromProb_Val[cptNumber]=maxVal;
			cptNumber++;
			continue;
		}
		else
		{
			if (maxVal >=cpfromProb_Val[cptNumber - 1])
			{
				cpfromSum_Pos[cptNumber - 1]=i;
				cpfromSum_Val[cptNumber - 1]=mem[i];
				cpfromProb_Pos[cptNumber - 1]=maxIdx;
				cpfromProb_Val[cptNumber - 1]=maxVal;
				continue;
			}
		}
	}
	if (cptNumber==0) { return cptNumber; }
	quickSortD(cpfromSum_Val,cpfromProb_Pos,0,cptNumber - 1);
	cptNumber=min(cptNumber,maxCptNumber);
	r_cblas_scopy(cptNumber,(float *)cpfromProb_Pos,1,(float *) cpt,1);
	I32PTR index=(int32_t * ) mem;
	F32PTR cpt_float=(float	* ) mem+N;
	for (int32_t i=0; i < cptNumber; i++)
	{
		*index++=i;
		*cpt_float++=(float)cpt[i];
	}
	index=index - cptNumber;
	cpt_float=cpt_float - cptNumber;
	quickSortA(cpt_float,index,0,cptNumber - 1);
	for (int32_t i=0; i < cptNumber; i++)
	{
		cptCI[i]=-9999.f;
		cptCI[cptNumber+i]=-9999.f;
	}
	float	delta;
	F32PTR tmpSeg=(float	* )  mem+3 * N;
	I32PTR nullSeg=(int32_t	*) mem+4 * N;
	for (int32_t i=0; i < cptNumber; i++)
	{
		int32_t startIdx,endIdx;
		endIdx=(int32_t) cpt_float[i];
		startIdx=i==0 ? 0 : (int32_t) cpt_float[i-1];
		startIdx=(startIdx+endIdx)/2;
		int32_t len=endIdx - startIdx+1;
		r_cblas_scopy(len,prob+startIdx,1,tmpSeg,1);
		quickSortA(tmpSeg,nullSeg,0,len - 1);
		delta=confidenceInterval(tmpSeg,len,'L');
		cptCI[i]=delta;
		startIdx=(int32_t) cpt_float[i];
		endIdx=i==(cptNumber - 1) ? (N - 1) : (int32_t) cpt_float[i+1];
		endIdx=(startIdx+endIdx)/2;
	    len=endIdx - startIdx+1;
		r_cblas_scopy(len,prob+startIdx,1,tmpSeg,1);
		quickSortD(tmpSeg,nullSeg,0,len - 1);
		delta=confidenceInterval(tmpSeg,len,'R');		
		cptCI[cptNumber+i]=delta;
	}
	float *temp=mem+2 * N;
	r_cblas_scopy(2 * cptNumber,cptCI,1,temp,1);
	for (int32_t i=0; i < cptNumber; i++)
	{
		int32_t idx;
		idx=index[i];
		cptCI[idx]=cpt_float[i] - temp[i];
		cptCI[cptNumber+idx]=cpt_float[i]+temp[cptNumber+i];
	}
	return cptNumber;
}
static INLINE void swapValue(rF32PTR a,rF32PTR b)
{
	rF32 t=*a;
	*a=*b;
	*b=t;
}
static INLINE void swapIndex(I32PTR a,I32PTR b)
{
	rI32 t=*a;
	*a=*b;
	*b=t;
}
int32_t partitionD(rF32PTR arr,I32PTR index,int32_t low,int32_t high)
{
	rF32 pivot=arr[high];    
	rI32 i=(low - 1);  
	for (rI32 j=low; j <=high - 1; j++)
	{
		if (arr[j] > pivot)
		{
			i++;    
			swapValue(&arr[i],&arr[j]);
			swapIndex(&index[i],&index[j]);
		}
	}
	swapValue(&arr[i+1],&arr[high]);
	swapIndex(&index[i+1],&index[high]);
	return (i+1);
}
void quickSortD(rF32PTR arr,I32PTR index,int32_t low,int32_t high)
{
	if (low < high)
	{
		int32_t pi=partitionD(arr,index,low,high);
		quickSortD(arr,index,low,pi - 1);
		quickSortD(arr,index,pi+1,high);
	}
}
int partitionA(rF32PTR arr,I32PTR index,int32_t low,int32_t high)
{
	float pivot=arr[high];    
	int32_t i=(low - 1);  
	for (int32_t j=low; j <=high - 1; j++)
	{
		if (arr[j] <=pivot)
		{
			i++;    
			swapValue(&arr[i],&arr[j]);
			swapIndex(&index[i],&index[j]);
		}
	}
	swapValue(&arr[i+1],&arr[high]);
	swapIndex(&index[i+1],&index[high]);
	return (i+1);
}
void quickSortA(rF32PTR arr,I32PTR index,int32_t low,int32_t high)
{
	if (low < high)
	{
		int32_t pi=partitionA(arr,index,low,high);
		quickSortA(arr,index,low,pi - 1);
		quickSortA(arr,index,pi+1,high);
	}
}
float fastlog1(float x)
{
	register union { float f; uint32_t i; } vx={ x };
	register union{ uint32_t i; float f; } mx={ (vx.i & 0x007FFFFF)|0x3f000000 };
	vx.f=(float)vx.i* 1.1920928955078125e-7f*0.69314718f;
	vx.f=vx.f - 124.22551499f*0.69314718f
		- 1.498030302f*0.69314718f * mx.f
		- 1.72587999f*0.69314718f/(0.3520887068f+mx.f);
	return vx.f;
}
float fastlog2(float x)
{
	register union { float f; uint32_t i; } vx={ x };
	register union{ uint32_t i; float f; } mx={ (vx.i & 0x007FFFFF)|0x3f000000 };
	vx.f=(float) ( (double)vx.i* (double)(1.1920928955078125e-7f*0.69314718f) );
	vx.f=vx.f - 124.2098722217282f*0.69314718f
		- 1.502704726439149f*0.69314718f * mx.f
		- 1.746553042329640f*0.69314718f/(0.356745518976678f+mx.f);
	return vx.f;
}
float fastlog(float x)
{
	register union { float f; uint32_t i; } vx={ x };
	register union { uint32_t i; float f; } mx={ (vx.i & 0x007FFFFF)|0x3f000000 };
	vx.f=(float)vx.i* (1.1920928955078125e-7f*0.69314718f);
	vx.f=vx.f - 125.5233166734556f*0.69314718f+mx.f*(-0.413356886671142+mx.f*(-0.472721975352920+0.078018528401178*mx.f))*0.69314718f+
		-0.787757784962750f*0.69314718f/(0.1781810261970705f+mx.f);
	return vx.f;
}
float sum_log_diag(rF32PTR p,rI32 K)
{
	rF32 x=0;
	for (rI32 i=0; i < K; i++)
	{
		register union { float f; uint32_t i; } vx={ *p };
		register union{ uint32_t i; float f; } mx={ (vx.i & 0x007FFFFF)|0x3f000000 };
		vx.f=(float)vx.i* (1.1920928955078125e-7f*0.69314718f);
		vx.f=vx.f - 125.5233166734556f*0.69314718f+mx.f*(-0.413356886671142+mx.f*(-0.472721975352920+0.078018528401178*mx.f))*0.69314718f+
			-0.787757784962750f*0.69314718f/(0.1781810261970705f+mx.f);
		x+=vx.f;
		p=p+K+1;
	}
	return x;
}
float fastexp(float x){
	x=1.442695040f*x;
	x=(x < -126) ? -126.0f : x;
	register float z=x - (float)((int)x)+((x < 0) ? 1.0f : 0.0f);
	register  union { uint32_t i; float f; } v;
	v.i=(uint32_t)     (    8388608.f * (x+121.2740575f+27.7280233f/(4.84252568f - z) - 1.49012907f * z)          );
	return v.f;
}
float fast_sqrt(float x)
{
	register union {
		uint32_t i;
		float   f;
	} v;
	v.f=x; 
	v.i -=1 << 23; 
	v.i >>=1; 
	v.i+=1 << 29; 
	v.f=(v.f+x/v.f);
	v.f=(v.f*0.25f+x/v.f);
	return v.f;
}
