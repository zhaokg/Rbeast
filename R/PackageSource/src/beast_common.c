#include "abc_000_macro.h"
DISABLE_MANY_WARNINGS
#include "beast_common.h"
Options     *GLOBAL_OPTIONS;
RESULT      *GLOBAL_RESULT;
void zeroOut_Xmars_zero0(F32PTR Xt_mars,F32PTR Xt_zeroBackup,
	U32PTR rowsMissing,int32_t nMissing,int32_t N,int32_t Npad,int32_t k)
{
	rI32 numPart=nMissing/4;
	for (rI32 j=k; j>0; j--)	{
		rI32 i=0;
		while (i<numPart){
			register  ptrdiff_t tmpidx0,tmpidx1,tmpidx2,tmpidx3;
			tmpidx0=*(rowsMissing);;
			*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
			Xt_mars[tmpidx0 - 1]=0.f;
			tmpidx1=*(rowsMissing+1);
			*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
			Xt_mars[tmpidx1 - 1]=0.f;
			tmpidx2=*(rowsMissing+2);
			*(Xt_zeroBackup+2)=Xt_mars[tmpidx2 - 1];
			Xt_mars[tmpidx2 - 1]=0.f;
			tmpidx3=*(rowsMissing+3);
			*(Xt_zeroBackup+3)=Xt_mars[tmpidx3 - 1];
			Xt_mars[tmpidx3 - 1]=0.f;
			rowsMissing+=4;
			Xt_zeroBackup+=4;
			i++;
		}
		rowsMissing=rowsMissing - nMissing;
		Xt_mars=Xt_mars+Npad;
	}
}
void zeroOut_Xmars_zero1(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,
	int32_t nMissing,int32_t N,int32_t Npad,int32_t k)
{
	rI32 numPart=nMissing/4;
	for (rI32 j=k; j>0; j--)	{
		rI32 i=0;
		register  ptrdiff_t tmpidx0,tmpidx1,tmpidx2,tmpidx3;
		while (i<numPart){
			tmpidx0=*(rowsMissing);;
			*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
			Xt_mars[tmpidx0 - 1]=0.f;
			tmpidx1=*(rowsMissing+1);
			*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
			Xt_mars[tmpidx1 - 1]=0.f;
			tmpidx2=*(rowsMissing+2);
			*(Xt_zeroBackup+2)=Xt_mars[tmpidx2 - 1];
			Xt_mars[tmpidx2 - 1]=0.f;
			tmpidx3=*(rowsMissing+3);
			*(Xt_zeroBackup+3)=Xt_mars[tmpidx3 - 1];
			Xt_mars[tmpidx3 - 1]=0.f;
			rowsMissing+=4;
			Xt_zeroBackup+=4;
			i++;
		}
		tmpidx0=*rowsMissing++;
		*Xt_zeroBackup++=Xt_mars[tmpidx0 - 1];
		Xt_mars[tmpidx0 - 1]=0.f;
		rowsMissing=rowsMissing - nMissing;
		Xt_mars=Xt_mars+Npad;
	}
}
void zeroOut_Xmars_zero2(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,
	int32_t nMissing,int32_t N,int32_t Npad,int32_t k)
{
	rI32 numPart=nMissing/4;
	for (rI32 j=k; j>0; j--)	{
		rI32 i=0;
		register  ptrdiff_t tmpidx0,tmpidx1,tmpidx2,tmpidx3;
		while (i<numPart){
			tmpidx0=*(rowsMissing);;
			*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
			Xt_mars[tmpidx0 - 1]=0.f;
			tmpidx1=*(rowsMissing+1);
			*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
			Xt_mars[tmpidx1 - 1]=0.f;
			tmpidx2=*(rowsMissing+2);
			*(Xt_zeroBackup+2)=Xt_mars[tmpidx2 - 1];
			Xt_mars[tmpidx2 - 1]=0.f;
			tmpidx3=*(rowsMissing+3);
			*(Xt_zeroBackup+3)=Xt_mars[tmpidx3 - 1];
			Xt_mars[tmpidx3 - 1]=0.f;
			rowsMissing+=4;
			Xt_zeroBackup+=4;
			i++;
		}
		tmpidx0=*rowsMissing;
		*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
		Xt_mars[tmpidx0 - 1]=0.f;
		tmpidx1=*(rowsMissing+1);
		*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
		Xt_mars[tmpidx1 - 1]=0.f;
		rowsMissing+=2;
		Xt_zeroBackup+=2;
		rowsMissing=rowsMissing - nMissing;
		Xt_mars=Xt_mars+Npad;
	}
}
void zeroOut_Xmars_zero3(F32PTR Xt_mars,F32PTR Xt_zeroBackup,U32PTR rowsMissing,
	int32_t nMissing,int32_t N,int32_t Npad,int32_t k)
{
	rI32 numPart=nMissing/4;
	for (rI32 j=k; j>0; j--)	{
		rI32 i=0;
		register  ptrdiff_t tmpidx0,tmpidx1,tmpidx2,tmpidx3;
		while (i<numPart){
			tmpidx0=*(rowsMissing);;
			*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
			Xt_mars[tmpidx0 - 1]=0.f;
			tmpidx1=*(rowsMissing+1);
			*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
			Xt_mars[tmpidx1 - 1]=0.f;
			tmpidx2=*(rowsMissing+2);
			*(Xt_zeroBackup+2)=Xt_mars[tmpidx2 - 1];
			Xt_mars[tmpidx2 - 1]=0.f;
			tmpidx3=*(rowsMissing+3);
			*(Xt_zeroBackup+3)=Xt_mars[tmpidx3 - 1];
			Xt_mars[tmpidx3 - 1]=0.f;
			rowsMissing+=4;
			Xt_zeroBackup+=4;
			i++;
		}
		tmpidx0=*rowsMissing;
		*Xt_zeroBackup=Xt_mars[tmpidx0 - 1];
		Xt_mars[tmpidx0 - 1]=0.f;
		tmpidx1=*(rowsMissing+1);
		*(Xt_zeroBackup+1)=Xt_mars[tmpidx1 - 1];
		Xt_mars[tmpidx1 - 1]=0.f;
		tmpidx2=*(rowsMissing+2);
		*(Xt_zeroBackup+2)=Xt_mars[tmpidx2 - 1];
		Xt_mars[tmpidx2 - 1]=0.f;
		rowsMissing+=3;
		Xt_zeroBackup+=3;
		rowsMissing=rowsMissing - nMissing;
		Xt_mars=Xt_mars+Npad;
	}
}
int32_t find_index_by_csum(rU08PTR good,rI32 N,rI32 randIdx) 
{
	  rI32 newKnot;
	  rI32 count=0;
	  {
		  rI32 i,delta,N16;
		  N16=N/16;
		  for (i=0; i<N16; i++) {
			  int64_t sum;
			  sum=*((int64_t*)good)+*((int64_t*)good+1);
			  *((int32_t*)&sum)+=*((int32_t*)&sum+1);
			  *((int16_t*)&sum)+=*((int16_t*)&sum+1);
			  *((int8_t*)&sum)+=*((int8_t *)&sum+1);
			  delta=*((int8_t*)&sum);
			  count+=delta;
			  if (count >=randIdx) break;
			  good+=16;
		  }
		  count -=delta;
		  newKnot=i * 16;
	  }
	  {
		  rI32 j;
		  for (j=0; j < 16; j++) {
			  count+=good[j];
			  if (count==randIdx) break;
		  }
		  newKnot+=(j+1);
	  }
	  return newKnot;
  }
 int32_t int8_arr_sum(rU08PTR good,rI32 N)
 {
	  rI32	 good_num=0;
	  rI32      N16=N/16;
	  for (rI32 i=0; i<N16; i++) {
		  int64_t sum;
		  sum=*((int64_t*)good)+*((int64_t*)good+1);
		  *((int32_t*)&sum)=*((int32_t*)&sum)+*((int32_t*)&sum+1);
		  *((int16_t*)&sum)=*((int16_t*)&sum)+*((int16_t*)&sum+1);
		  *((int8_t*)&sum)=*((int8_t*)&sum)+*((int8_t*)&sum+1);
		  good_num=good_num+*((int8_t*)&sum);
		  good+=16;
	  }
	  return good_num;
  }
 void print_error(int code,MemPointers *MEM)
 {
	 if (!MEM)
		 MEM->free_all(MEM);
	 switch (code)
	 {
	 case 1:
		 r_error("Error (beastST):To run,you must save the output to either " 
			     "a lefthand-side variable (e.g.,out=func(...) ) or specify the output path \n"
			     "in the Option.outputFolder parameter. \n");
		 break;
	 case 2:
		 r_error("Can't open the input data file!");
		 break;
	 case 3:
		 r_error("r2<r1:There must be something wrong!\n");
		 break;
	 case 4:
		 r_error("j should be not at least 2. There must be something wrong! \n");
		 break;
	 case 5:
		 r_error("The two K's differ4; there must be something wrong!");
		 break;
	 default:	
		 ;
	 }
 }
 ENABLE_MANY_WARNINGS
