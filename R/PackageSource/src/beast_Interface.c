#include "abc_001_config.h"
DISABLE_MANY_WARNINGS
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include "abc_datatype.h"
#include "abc_blas_lapack_lib.h"
#include "abc_common.h"
#include "beast_common.h"
#include  "beast_lib.h" 
int check_options(Options * _restrict opt,char *missing)
{
	Options o;
	memcpy(&o,opt,sizeof(Options));
	if (o.outputType=='F' && o.output==NULL)
	{
		r_error("To run,you must save the output to either a lefthand-side"
			"variable (e.g.,out=beast(...) ) or specify the output path \n"
			"in the option$outputFolder parameter. \n");
		return 0;
	}
	rF32PTR yData=NULL;
	if (o.inputType=='F')
	{
		if (o.N<0)
		{			
			r_printf("For input from an external file,the time series length must be" 
				" specified by opt$lengthPerTimeSeries_infile");
			return 0;
		}
		FILE * fileID=fopen(o.input,"rb");
		if (fileID==NULL)
		{
			r_error("The specified file cannot be opened!");
			return 0;
		}
		fseek(fileID,0,SEEK_END);
		rI64 fileSizeinByte=ftell((FILE *)fileID);
		if (fileSizeinByte%(4 * o.N) !=0)
		{
			fclose(fileID);
			r_printf("The file size is%d bytes. It must be a multiplier of the time series length.",fileSizeinByte);
			return 0;
		}
		o.M=(uint32_t) (fileSizeinByte/(4 * o.N));
		r_printf("File Size:%d ; Number of time series/pixels:%d \n",fileSizeinByte,o.M);
		fseek(fileID,0,SEEK_SET);
		yData=malloc(sizeof(float)*o.N);
		fread(yData,sizeof(float),o.N,fileID);
		fclose(fileID);
	}
	if (missing[1]) o.omissionValue=-9999;
	if (missing[2]||o.period <=0)
	{ 
		if (yData==NULL) yData=malloc(sizeof(float)*o.N);
		if (o.inputType=='s')
		{
			rF32PTR  fltPtr=(float*)opt->input;
			for (rI32 i=o.N; i>0; i--)
				*yData++=(float)(*fltPtr++);
			yData=yData - o.N;
		}		
		else if (o.inputType=='d')
		{
			rF64PTR  dblPtr=(double*) opt->input;
			for (rI32 i=o.N; i>0; i--)
				*yData++=(float)(*dblPtr++);
			yData=yData - o.N;
		}
		else if (o.inputType=='4')
		{
			rI32PTR  intPtr=(int32_t *)opt->input;
			for (rI32 i=o.N; i>0; i--)
				*yData++=(float)(*intPtr++);
			yData=yData - o.N;
		}
		else if (o.inputType=='F')
		{
		}
		o.period=determine_period(yData,o.N,o.omissionValue);
		r_printf(
			"\nNote: The \"opt$period\" parameter should be known in advance and supplied by the user but it is missing. A best guess of it is%5d AND will be used in the decomposition. "
			"Please make sure this estimate makes sense; otherwise,the BEAST decomposition result will be incorrect.\n\n",(int) o.period );
	}
	else 
	{
		if (o.period<=0||fabsf(o.period - (int32_t)o.period) > 1e-5f)
		{
			if (yData !=NULL) free(yData);
			r_error("The value of period must be an positive integer!");
			return 0;
		}
	}
	if (yData !=NULL) free(yData);
	if (missing[3]) o.startTime=1.f;
	if (missing[4]) o.timeInterval=1.f;
	if (missing[5]) o.minSeasonOrder=1L;				 o.minSeasonOrder=min(o.minSeasonOrder,o.period/2 - 1); o.minSeasonOrder=max(o.minSeasonOrder,1L);
	if (missing[6]) o.maxSeasonOrder=(o.period/2 - 1);  o.maxSeasonOrder=min(o.maxSeasonOrder,(o.period/2 - 1));  o.maxSeasonOrder=max(o.maxSeasonOrder,o.minSeasonOrder);
	if (missing[7]) o.minTrendOrder=0L;				 o.minTrendOrder=max(o.minTrendOrder,0L);
	if (missing[8]) o.maxTrendOrder=1L;				 o.maxTrendOrder=max(o.maxTrendOrder,o.minTrendOrder);
	if (missing[9]) o.minSepDist_Trend=o.period/2;			 o.minSepDist_Trend=max(o.minSepDist_Trend,0L); o.minSepDist_Trend=min(o.minSepDist_Trend,o.N/2-1);
	if (missing[10]) o.minSepDist_Season=o.period/2;			 o.minSepDist_Season=max(o.minSepDist_Season,2L); o.minSepDist_Season=min(o.minSepDist_Season,o.N/2 - 1);
	if (missing[11]) o.maxKnotNum_Trend=floor(o.N/(o.minSepDist_Trend+1) - 1);	o.maxKnotNum_Trend=min(o.maxKnotNum_Trend,floor(o.N/(o.minSepDist_Trend+1) - 1));
	if (missing[12]) o.maxKnotNum_Season=floor(o.N/(o.minSepDist_Season+1) - 1);	o.maxKnotNum_Season=min(o.maxKnotNum_Season,floor(o.N/(o.minSepDist_Season+1) - 1));
	if (missing[13]) o.maxMoveStepSize=(o.minSepDist_Trend/3+1);	
	if (missing[14]) o.samples=3000;             o.samples=max(o.samples,800);
	if (missing[15]) o.thinningFactor=1L;        o.thinningFactor=max(o.thinningFactor,1L);
	if (missing[16]) o.burnin=200L;			   o.burnin=max(o.burnin,200L);
	if (missing[17]) o.chainNumber=3;			   o.chainNumber=max(o.chainNumber,1L);
	if (missing[18]) o.resamplingTrendOrderProb=.1f;		
	if (missing[19]) o.resamplingSeasonOrderProb=.17f;
	if (missing[20]) o.seed=0L;
	if (missing[21]) o.printToScreen=0L;
	if (missing[22]) o.printCharLen=80L;
	if (missing[23]) o.computeCredible=0L;
	if (missing[24])  o.fastCIComputation=1L;
	if (missing[25])  o.alphaLevel=.95;
	if (missing[26])  o.computeSlopeSign=0L;
	if (missing[27])  o.computeHarmonicOrder=0L;
	if (missing[28])  o.computeTrendOrder=0L;
	if (missing[29])  o.computeChangepoints=1L;
	if (missing[30])  o.ridgeFactor=0.0001f;
	memcpy(opt,&o,sizeof(Options));
	return 1;
}
#if R_INTERFACE==1
#include "abc_R_util.h"
#define IS_SCALAR(x,type) (TYPEOF(x)==(type) && XLENGTH(x)==1)
int R_check_input(SEXP Y,SEXP opt)
{
	if (   !(TYPEOF(Y)==INTSXP  && XLENGTH(Y) > 2)  &&
		   !(TYPEOF(Y)==REALSXP && XLENGTH(Y) > 2)  &&
		   !(TYPEOF(Y)==INTSXP  && isMatrix (Y))    && 
		   !(TYPEOF(Y)==REALSXP && isMatrix(Y) )    &&
		   !(TYPEOF(Y)==STRSXP  && XLENGTH(Y)==1) &&
		   !(TYPEOF(Y)==REALSXP  && Rf_isArray(Y) )
		)
	{
		r_error("Error: The input data should be a numeric type and must be long enough.\n");
		return 0;
	}
	if (!isNewList(opt) && !IS_SCALAR(opt,INTSXP) && !IS_SCALAR(opt,REALSXP))
	{
		r_error("Error: The input parameter OPTION should be a List variable or an integer!\n");
		return 0;
	}
	return 1;	
}
int R_read_input(Options * _restrict pOpt,SEXP Y,SEXP opt,char *missing)
{
	Options o;
	SEXP	tmpSEXP=NULL;
	int		nprt=0;
	if (TYPEOF(Y)==STRSXP)
	{	
		o.inputType='F';
		o.input=(void*)CHAR(STRING_ELT(Y,0));
		tmpSEXP=getListElement(opt,"lengthPerTimeSeries_infile");
		o.N=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,INTEGER_VALUE(tmpSEXP)) : -9999;
	}
	else if (isNumeric(Y))
	{ 
		if (TYPEOF(Y)==INTSXP)
		{
			o.inputType='4';
			o.input=(void *)INTEGER(Y);
		}
		else if (TYPEOF(Y)==REALSXP) 
		{
			o.inputType='d';
			o.input=(void *)REAL(Y);
		}
		SEXP dims=PROTECT(getAttrib(Y,R_DimSymbol)); nprt++;
		if (isMatrix(Y) && Rf_length(dims)==2L )
		{
			o.isInput3DStack=0;
			o.N=(int)INTEGER(dims)[0];
			o.M=(int)INTEGER(dims)[1];
		}
		else if (Rf_isArray(Y) && Rf_length(dims)==3L)
		{
			o.isInput3DStack=1;
			if (getListElement(opt,"period")==NULL||
				getListElement(opt,"timeDimensionIndex")==NULL)
			{
				UNPROTECT(nprt);
				nprt=0;
				r_printf("Error: If the input is a 3D stack,the opt.periord and "
					"opt.timeDimensionIndex must be supplied!\n");
				return 0;
			}
			tmpSEXP=getListElement(opt,"timeDimensionIndex");	
			o.timeDimensionIndex=(PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,asInteger(tmpSEXP));
			o.N=INTEGER(dims)[o.timeDimensionIndex - 1];
			if (o.timeDimensionIndex==1)
			{
				o.M=INTEGER(dims)[1];
				o.L=INTEGER(dims)[2];
			}
			else if (o.timeDimensionIndex==2)
			{
				o.M=INTEGER(dims)[0];
				o.L=INTEGER(dims)[2];
			}
			else
			{
				o.M=INTEGER(dims)[0];
				o.L=INTEGER(dims)[1];
			}
		}
		else if (Rf_length(dims)==0L)
		{
			o.isInput3DStack=0;
			o.N=Rf_length(Y);
			o.M=1;
		}
	}
	if (isNewList(opt))
	{
		tmpSEXP=getListElement(opt,"outputToDisk");
		o.outputToDisk=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : 0;
		if (o.outputToDisk)
		{
			o.outputType='F';
			tmpSEXP=getListElement(opt,"outputFolder");
			o.output=(tmpSEXP !=NULL) ? (void *) CHAR(STRING_ELT(tmpSEXP,0)) : NULL;
		}
		else
		{
			o.outputType='d';
			o.output=NULL;
		}
		memset(missing,0L,31);
		tmpSEXP=getListElement(opt,"omissionValue");
		o.omissionValue=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[1]=1);
		tmpSEXP=getListElement(opt,"period");	
		o.period=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[2]=1);
		tmpSEXP=getListElement(opt,"startTime");
		o.startTime=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[3]=1);
		tmpSEXP=getListElement(opt,"timeInterval");
		o.timeInterval=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : -(missing[4]=1);
		tmpSEXP=getListElement(opt,"minSeasonOrder");
		o.minSeasonOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : (missing[5]=1);
		tmpSEXP=getListElement(opt,"maxSeasonOrder");
		o.maxSeasonOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : (missing[6]=1);
		tmpSEXP=getListElement(opt,"minTrendOrder");
		o.minTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : (missing[7]=1);
		tmpSEXP=getListElement(opt,"maxTrendOrder");
		o.maxTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(unsigned char)asInteger(tmpSEXP)) : (missing[8]=1);
		tmpSEXP=getListElement(opt,"minSepDist_Trend");
		o.minSepDist_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[9]=1);
		tmpSEXP=getListElement(opt,"minSepDist_Season");
		o.minSepDist_Season=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[10]=1);
		tmpSEXP=getListElement(opt,"maxKnotNum_Trend");
		o.maxKnotNum_Trend=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[11]=1);
		tmpSEXP=getListElement(opt,"maxKnotNum_Season");
		o.maxKnotNum_Season=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[12]=1);
		tmpSEXP=getListElement(opt,"maxMoveStepSize");
		o.maxMoveStepSize=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[13]=1);
		tmpSEXP=getListElement(opt,"samples");
		o.samples=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : (missing[14]=1);
		tmpSEXP=getListElement(opt,"thinningFactor");
		o.thinningFactor=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint16_t)asInteger(tmpSEXP)) : (missing[15]=1);
		tmpSEXP=getListElement(opt,"burnin");
		o.burnin=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : (missing[16]=1);
		tmpSEXP=getListElement(opt,"chainNumber");
		o.chainNumber=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(mwSize)asInteger(tmpSEXP)) : (missing[17]=1);
		tmpSEXP=getListElement(opt,"resamplingTrendOrderProb");
		o.resamplingTrendOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[18]=1);
		tmpSEXP=getListElement(opt,"resamplingSeasonOrderProb");
		o.resamplingSeasonOrderProb=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[19]=1);
		tmpSEXP=getListElement(opt,"seed");
		o.seed=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(uint64_t)asInteger(tmpSEXP)) : (missing[20]=1);
		tmpSEXP=getListElement(opt,"printToScreen");
		o.printToScreen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[21]=1);
		tmpSEXP=getListElement(opt,"printCharLen");
		o.printCharLen=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(int16_t)asInteger(tmpSEXP)) : (missing[22]=1);
		tmpSEXP=getListElement(opt,"computeCredible");
		o.computeCredible=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[23]=1);
		tmpSEXP=getListElement(opt,"fastCIComputation");
		o.fastCIComputation=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[24]=1);
		tmpSEXP=getListElement(opt,"alphaLevel");
		o.alphaLevel=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[25]=1);
		tmpSEXP=getListElement(opt,"computeSlopeSign");
		o.computeSlopeSign=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[26]=1);
		tmpSEXP=getListElement(opt,"computeHarmonicOrder");
		o.computeHarmonicOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[27]=1);
		tmpSEXP=getListElement(opt,"computeTrendOrder");
		o.computeTrendOrder=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[28]=1);
		tmpSEXP=getListElement(opt,"computeChangepoints");
		o.computeChangepoints=(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,INTSXP)),nprt++,(bool)asInteger(tmpSEXP)) : (missing[29]=1);
		tmpSEXP=getListElement(opt,"ridgeFactor");
		o.ridgeFactor=(float)(tmpSEXP !=NULL) ? (PROTECT(tmpSEXP=coerceVector(tmpSEXP,REALSXP)),nprt++,(float)asReal(tmpSEXP)) : (missing[30]=1);
	}
	else if (isNumeric(opt))
	{
		o.outputToDisk=0;
		o.output=NULL;
		o.outputType='d';
		memset(missing,1L,31);
		missing[2]=0;
		o.period=(PROTECT(tmpSEXP=coerceVector(opt,REALSXP)),nprt++,(float)asReal(tmpSEXP));
	}  
	UNPROTECT(nprt);
	nprt=0;
	o.Npad=(int32_t)ceil((float)o.N/8.0f) * 8;
	o.Npad16=(int32_t)ceil((float)o.N/16.f) * 16;
	o.separator='$';
	memcpy(pOpt,&o,sizeof(Options));
	return 1;
}
SEXP R_allocate_output(RESULT * _restrict matOutput,Options * _restrict opt)
{
	SEXP ANS=NULL;
	SEXP ansListNames=NULL;
	SEXP tmpSEXP;
	int idx;
	int nprt=0L;
	if (opt->isInput3DStack==0)
	{
		int N=opt->N;
		int M=opt->M;
		int32_t elementNum=15+(int32_t)opt->computeCredible * 3+(int32_t)opt->computeSlopeSign  \
			+(int32_t)opt->computeHarmonicOrder+(int32_t)opt->computeTrendOrder * 1+(int32_t)opt->computeChangepoints * 4;
		PROTECT(ANS=allocVector(VECSXP,elementNum));++nprt;
		PROTECT(ansListNames=allocVector(STRSXP,elementNum));++nprt;
		idx=0; PROTECT(tmpSEXP=allocVector(REALSXP,N));++nprt;
		matOutput->time=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("time"));
		idx++; PROTECT(tmpSEXP=allocVector(REALSXP,M));++nprt;
		matOutput->sN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sN"));
		idx++; PROTECT(tmpSEXP=allocVector(REALSXP,M));++nprt;
		matOutput->tN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tN"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,(opt->maxKnotNum_Season+1),M));++nprt;
		matOutput->sNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sNProb"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,(opt->maxKnotNum_Trend+1),M));++nprt;
		matOutput->tNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tNProb"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->sProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sProb"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->tProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tProb"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->s=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("s"));
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,M));++nprt;
			matOutput->sCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sCI"));
		}
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->sSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sSD"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->t=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("t"));
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,M));++nprt;
			matOutput->tCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tCI"));
		}
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->tSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tSD"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->b=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("b"));
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,2,M));++nprt;
			matOutput->bCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bCI"));
		}
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
		matOutput->bSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bSD"));
		idx++; PROTECT(tmpSEXP=allocVector(REALSXP,M));++nprt;
		matOutput->marg_lik=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("marg_lik"));
		idx++; PROTECT(tmpSEXP=allocVector(REALSXP,M));++nprt;
		matOutput->sig2=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sig2"));
		if (opt->computeSlopeSign)
		{
			idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
			matOutput->bsign=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bsign"));
		}
		if (opt->computeHarmonicOrder)
		{
			idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
			matOutput->horder=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("horder"));
		}
		if (opt->computeTrendOrder)
		{
			idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,N,M));++nprt;
			matOutput->torder=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("torder"));
		}
		if (opt->computeChangepoints)
		{
			idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,opt->maxKnotNum_Season,M));++nprt;
			matOutput->scp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scp"));
			idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,opt->maxKnotNum_Trend,M));++nprt;
			matOutput->tcp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcp"));
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Season,2,M));++nprt;
			matOutput->scpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scpCI"));
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Trend,2,M));++nprt;
			matOutput->tcpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcpCI"));
		}
		setAttrib(ANS,R_NamesSymbol,ansListNames);
		PROTECT(tmpSEXP=mkString("beast"));++nprt;
		setAttrib(ANS,R_ClassSymbol,tmpSEXP);
		PROTECT(tmpSEXP=mkString("beastST"));++nprt;
		setAttrib(ANS,install("algorithm"),tmpSEXP);
		double *doublePointer=(double*)matOutput->time;
		for (ptrdiff_t i=(ptrdiff_t)(N - 1); i >=0; i--)
			*(doublePointer+i)=(double)*(matOutput->time+i);
		UNPROTECT(nprt);
		return  ANS;
	}
	else
	{
		int N=opt->N;
		int M=opt->M;
		int L=opt->L;
		SEXP dim4;
		PROTECT(dim4=allocVector(INTSXP,4));++nprt;
		int32_t elementNum=15+(int32_t)opt->computeCredible * 3+(int32_t)opt->computeSlopeSign  \
			+(int32_t)opt->computeHarmonicOrder+(int32_t)opt->computeTrendOrder * 1+(int32_t)opt->computeChangepoints * 4;
		PROTECT(ANS=allocVector(VECSXP,elementNum));++nprt;
		PROTECT(ansListNames=allocVector(STRSXP,elementNum));++nprt;
		idx=0; PROTECT(tmpSEXP=allocVector(REALSXP,N));++nprt;
		matOutput->time=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("time"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,M,L));++nprt;
		matOutput->sN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sN"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,M,L));++nprt;
		matOutput->tN=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tN"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,(opt->maxKnotNum_Season+1),M,L));++nprt;
		matOutput->sNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sNProb"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,(opt->maxKnotNum_Trend+1),M,L));++nprt;
		matOutput->tNProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tNProb"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->sProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sProb"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->tProb=(int32_t *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tProb"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->s=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("s"));
		INTEGER(dim4)[0]=N; INTEGER(dim4)[1]=2; INTEGER(dim4)[2]=M; INTEGER(dim4)[3]=L;
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=Rf_allocArray(REALSXP,dim4));++nprt;
			matOutput->sCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sCI"));
		}
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->sSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sSD"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->t=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("t"));
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=Rf_allocArray(REALSXP,dim4));++nprt;
			matOutput->tCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tCI"));
		}
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->tSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tSD"));
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->b=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("b"));
		if (opt->computeCredible) {
			idx++; PROTECT(tmpSEXP=Rf_allocArray(REALSXP,dim4));++nprt;
			matOutput->bCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bCI"));
		}
		idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
		matOutput->bSD=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bSD"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,M,L));++nprt;
		matOutput->marg_lik=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("marg_lik"));
		idx++; PROTECT(tmpSEXP=allocMatrix(REALSXP,M,L));++nprt;
		matOutput->sig2=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("sig2"));
		if (opt->computeSlopeSign)
		{
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
			matOutput->bsign=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("bsign"));
		}
		if (opt->computeHarmonicOrder)
		{
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
			matOutput->horder=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("horder"));
		}
		if (opt->computeTrendOrder)
		{
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,N,M,L));++nprt;
			matOutput->torder=(I32PTR)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("torder"));
		}
		if (opt->computeChangepoints)
		{
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Season,M,L));++nprt;
			matOutput->scp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scp"));
			idx++; PROTECT(tmpSEXP=alloc3DArray(REALSXP,opt->maxKnotNum_Trend,M,L));++nprt;
			matOutput->tcp=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcp"));
			INTEGER(dim4)[0]=opt->maxKnotNum_Season; INTEGER(dim4)[1]=2; INTEGER(dim4)[2]=M; INTEGER(dim4)[3]=L;
			idx++; PROTECT(tmpSEXP=allocArray(REALSXP,dim4));++nprt;
			matOutput->scpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("scpCI"));
			INTEGER(dim4)[0]=opt->maxKnotNum_Trend; 
			idx++; PROTECT(tmpSEXP=allocArray(REALSXP,dim4));++nprt;
			matOutput->tcpCI=(float *)REAL(tmpSEXP); SET_VECTOR_ELT(ANS,idx,tmpSEXP); SET_STRING_ELT(ansListNames,idx,mkChar("tcpCI"));
		}
		setAttrib(ANS,R_NamesSymbol,ansListNames);
		PROTECT(tmpSEXP=mkString("beast"));++nprt;
		setAttrib(ANS,R_ClassSymbol,tmpSEXP);
		PROTECT(tmpSEXP=mkString("beastST"));++nprt;
		setAttrib(ANS,install("algorithm"),tmpSEXP);
		double *doublePointer=(double*)matOutput->time;
		for (ptrdiff_t i=(ptrdiff_t)(N - 1); i >=0; i--)
			*(doublePointer+i)=(double)*(matOutput->time+i);
		UNPROTECT(nprt);
		return  ANS;
	}
}
#elif M_INTERFACE==1
int  M_check_input(int nrhs,const mxArray * _restrict prhs[])
{
	if (nrhs==0)
	{
		r_error("Error: There should be at least one input argument!.\n");
		return 0;
	}
	const mxArray * S;
	S=prhs[0];
	int numel=mxGetNumberOfElements(S);
	if (!( mxIsDouble(S) && numel > 2) &&
		!( mxIsSingle(S) && numel > 2) &&
		!( mxIsInt32(S) && numel > 2) &&
		!( mxIsChar(S)  )
		)
	{
		r_error("Error: The input data should be a numeric type and must be long enough.\n");
		return 0;
	}
	if (2==nrhs)
	{
		S=prhs[1];
		numel=mxGetNumberOfElements(S);
		if (!mxIsStruct(S) && !(mxIsNumeric(S) && numel==1))
		{
			r_error("Error: The input parameter OPTION should be a struct variable or an integer!\n");
			return 0;
		}
	}
	return 1;
}
int     M_read_input(Options * _restrict pOpt,int nrhs,const mxArray * _restrict prhs[],char *missing)
{
	Options        o;
	mxArray       *tmp;
	const mxArray *S;
	S=prhs[0];
	if (mxIsChar(S))
	{	
		o.inputType='F';
		o.input=S;
		o.N=-9999;
		if (nrhs >=2 && mxIsStruct(prhs[1]))
		{
			tmp=mxGetField(prhs[1],0,"lengthPerTimeSeries_infile");
			o.N=(tmp !=NULL) ? mxGetScalar((mxArray*)tmp) : -9999;
		}
	}
	else if (mxIsNumeric(S))
	{ 
		o.input=mxGetData(S);
		o.inputType='u';
		if (mxIsInt32(S))
		{
			o.inputType='4';		
		}
		else if (mxIsDouble(S))
		{
			o.inputType='d';	
		}
		else if (mxIsSingle(S))
		{
			o.inputType='s';
		}
		if (mxGetNumberOfDimensions(S) !=3)
		{
			o.isInput3DStack=0;
			int N=mxGetM(S);
			int M=mxGetN(S);
			if (min(N,M)>1L)
			{
				o.N=N;
				o.M=M;
			}
			else
			{
				o.N=max(N,M);
				o.M=1L;
			}
		}
		else
		{
			o.isInput3DStack=1;
			if (   mxGetField(prhs[1],0,"period")==NULL||
				   mxGetField(prhs[1],0,"timeDimensionIndex")==NULL )
			{
				r_printf("Error: If the input is a 3D stack,the opt.periord and "
					"opt.timeDimensionIndex must be supplied!\n");
				return 0;
			}
			o.timeDimensionIndex=mxGetScalar(mxGetField(prhs[1],0,"timeDimensionIndex"));
			const mwSize *dims;
			dims=mxGetDimensions(S);
			o.N=dims[o.timeDimensionIndex - 1];
			if (o.timeDimensionIndex==1)
			{
				o.M=dims[1];
				o.L=dims[2];
			}
			else if (o.timeDimensionIndex==2)
			{
				o.M=dims[0];
				o.L=dims[2];
			}
			else
			{
				o.M=dims[0];
				o.L=dims[1];
			}
		}
	}
	S=prhs[1];
	if (mxIsStruct(S))
	{
		o.outputToDisk=((tmp=mxGetField(S,0,"outputToDisk")) !=NULL) ? mxGetScalar((mxArray*)tmp) : 0;
		if (o.outputToDisk)
		{
			o.outputType='F';			
			o.output=((tmp=mxGetField(S,0,"outputFolder")) !=NULL) ? tmp : NULL;
		}
		else
		{
			o.outputType='s';
			o.output=NULL;
		}
		memset(missing,0L,31);
		o.omissionValue=((tmp=mxGetField(S,0,"omissionValue")) !=NULL) ? mxGetScalar(tmp) :  (missing[1]=1);
		o.period=((tmp=mxGetField(S,0,"period")) !=NULL) ? mxGetScalar(tmp) : (missing[2]=1);
		o.startTime=((tmp=mxGetField(S,0,"startTime")) !=NULL) ? mxGetScalar(tmp) : (missing[3]=1);
		o.timeInterval=((tmp=mxGetField(S,0,"timeInterval")) !=NULL) ? mxGetScalar(tmp) : -(missing[4]=1);
		o.minSeasonOrder=((tmp=mxGetField(S,0,"minSeasonOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[5]=1);
		o.maxSeasonOrder=((tmp=mxGetField(S,0,"maxSeasonOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[6]=1);
		o.minTrendOrder=((tmp=mxGetField(S,0,"minTrendOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[7]=1);
		o.maxTrendOrder=((tmp=mxGetField(S,0,"maxTrendOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[8]=1);
		o.minSepDist_Trend=((tmp=mxGetField(S,0,"minSepDist_Trend")) !=NULL) ? mxGetScalar(tmp) : (missing[9]=1);
		o.minSepDist_Season=((tmp=mxGetField(S,0,"minSepDist_Season")) !=NULL) ? mxGetScalar(tmp) : (missing[10]=1);
		o.maxKnotNum_Trend=((tmp=mxGetField(S,0,"maxKnotNum_Trend")) !=NULL) ? mxGetScalar(tmp) : (missing[11]=1);
		o.maxKnotNum_Season=((tmp=mxGetField(S,0,"maxKnotNum_Season")) !=NULL) ? mxGetScalar(tmp) : (missing[12]=1);
		o.maxMoveStepSize=((tmp=mxGetField(S,0,"maxMoveStepSize")) !=NULL) ? mxGetScalar(tmp) : (missing[13]=1);
		o.samples=((tmp=mxGetField(S,0,"samples")) !=NULL) ? mxGetScalar(tmp) : (missing[14]=1);
		o.thinningFactor=((tmp=mxGetField(S,0,"thinningFactor")) !=NULL) ? mxGetScalar(tmp) : (missing[15]=1);
		o.burnin=((tmp=mxGetField(S,0,"burnin")) !=NULL) ? mxGetScalar(tmp) : (missing[16]=1);
		o.chainNumber=((tmp=mxGetField(S,0,"chainNumber")) !=NULL) ? mxGetScalar(tmp) : (missing[17]=1);
		o.resamplingTrendOrderProb=((tmp=mxGetField(S,0,"resamplingTrendOrderProb")) !=NULL) ? mxGetScalar(tmp) : (missing[18]=1);
		o.resamplingSeasonOrderProb=((tmp=mxGetField(S,0,"resamplingSeasonOrderProb")) !=NULL) ? mxGetScalar(tmp) : (missing[19]=1);
		o.seed=((tmp=mxGetField(S,0,"seed")) !=NULL) ? mxGetScalar(tmp) : (missing[20]=1);
		o.printToScreen=((tmp=mxGetField(S,0,"printToScreen")) !=NULL) ? mxGetScalar(tmp) : (missing[21]=1);
		o.printCharLen=((tmp=mxGetField(S,0,"printCharLen")) !=NULL) ? mxGetScalar(tmp) : (missing[22]=1);
		o.computeCredible=((tmp=mxGetField(S,0,"computeCredible")) !=NULL) ? mxGetScalar(tmp) : (missing[23]=1);
		o.fastCIComputation=((tmp=mxGetField(S,0,"fastCIComputation")) !=NULL) ? mxGetScalar(tmp) : (missing[24]=1);
		o.alphaLevel=((tmp=mxGetField(S,0,"alphaLevel")) !=NULL) ? mxGetScalar(tmp) : (missing[25]=1);
		o.computeSlopeSign=((tmp=mxGetField(S,0,"computeSlopeSign")) !=NULL) ? mxGetScalar(tmp) : (missing[26]=1);
		o.computeHarmonicOrder=((tmp=mxGetField(S,0,"computeHarmonicOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[27]=1);
		o.computeTrendOrder=((tmp=mxGetField(S,0,"computeTrendOrder")) !=NULL) ? mxGetScalar(tmp) : (missing[28]=1);
		o.computeChangepoints=((tmp=mxGetField(S,0,"computeChangepoints")) !=NULL) ? mxGetScalar(tmp) : (missing[29]=1);
		o.ridgeFactor=((tmp=mxGetField(S,0,"ridgeFactor")) !=NULL) ? mxGetScalar(tmp) : (missing[30]=1);
	}
	else if (mxIsNumeric(S))
	{
		o.outputToDisk=0;
		o.output=NULL;
		o.outputType='s';
		memset(missing,1L,31);
		missing[2]=0;
		o.period=mxGetScalar(S);
	}
	o.Npad=(int32_t)ceil((float)o.N/8.0f) * 8;
	o.Npad16=(int32_t)ceil((float)o.N/16.f) * 16;
	o.separator='.';
	memcpy(pOpt,&o,sizeof(Options));
	return 1;
}
mxArray *M_allocate_output(RESULT * _restrict matOutput,Options * _restrict opt)
{
	mxArray * _restrict out;
	mxArray * _restrict mxPointer;
	int  fieldNumber=25;
	char *fieldNames[]={
		"time","sN","tN","sNProb","tNProb","sProb","tProb","s","sCI","sSD",\
		"t","tCI","tSD","b","bCI","bSD","marg_lik","sig2","bsign","horder","torder","tcp","scp","tcpCI","scpCI" };
	if (opt->isInput3DStack==0)
	{
		mwSize	 dim2[2]={ 1,1 };
		mwSize   N=opt->N;
		mwSize   M=opt->M;
		mwSize 	 dim3[3]={ N,2,M };
		out=mxCreateStructArray(2,dim2,fieldNumber,fieldNames);
		mxPointer=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"time",mxPointer); 	matOutput->time=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(1,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sN",mxPointer); 	matOutput->sN=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(1,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tN",mxPointer);	matOutput->tN=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Season+1,M,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"sNProb",mxPointer); 	matOutput->sNProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend+1,M,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"tNProb",mxPointer);	matOutput->tNProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sProb",mxPointer);	matOutput->sProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tProb",mxPointer); matOutput->tProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"s",mxPointer); 	matOutput->s=mxGetData(mxPointer);
		if (opt->computeCredible) mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"sCI",mxPointer),matOutput->sCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sSD",mxPointer); 	matOutput->sSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"t",mxPointer); 	matOutput->t=mxGetData(mxPointer);
		if (opt->computeCredible) 	mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"tCI",mxPointer),matOutput->tCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tSD",mxPointer); 	matOutput->tSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"b",mxPointer); 	matOutput->b=mxGetData(mxPointer);
		if (opt->computeCredible)  mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bCI",mxPointer),matOutput->bCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"bSD",mxPointer); 	matOutput->bSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(1,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"marg_lik",mxPointer); matOutput->marg_lik=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(1,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sig2",mxPointer);    matOutput->sig2=mxGetData(mxPointer);
		if (opt->computeSlopeSign) 	   mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bsign",mxPointer),matOutput->bsign=mxGetData(mxPointer);
		if (opt->computeHarmonicOrder) mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"horder",mxPointer),matOutput->horder=mxGetData(mxPointer);
		if (opt->computeTrendOrder) mxPointer=mxCreateNumericMatrix(N,M,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"torder",mxPointer),matOutput->torder=mxGetData(mxPointer);
		if (opt->computeChangepoints)
		{
			mwSize 	 dim1_3[3]={ opt->maxKnotNum_Season,2,M };
			mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Season,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"scp",mxPointer);    matOutput->scp=mxGetData(mxPointer);
			mxPointer=mxCreateNumericArray(3,dim1_3,mxSINGLE_CLASS,mxREAL); 						mxSetField(out,0,"scpCI",mxPointer);  matOutput->scpCI=mxGetData(mxPointer);
			mwSize 	 dim2_3[3]={ opt->maxKnotNum_Trend,2,M };
			mxPointer=mxCreateNumericMatrix(opt->maxKnotNum_Trend,M,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tcp",mxPointer);    matOutput->tcp=mxGetData(mxPointer);
			mxPointer=mxCreateNumericArray(3,dim2_3,mxSINGLE_CLASS,mxREAL); 	 					mxSetField(out,0,"tcpCI",mxPointer);  matOutput->tcpCI=mxGetData(mxPointer);
		}
		fill_1_to_N(matOutput->time,N);
		r_ippsMulC_32f_I(opt->timeInterval,matOutput->time,(const int)N);
		r_ippsSubC_32f_I(-(opt->startTime - opt->timeInterval),matOutput->time,(const int)N);
		return out;
	}
	else
	{
		mwSize	 dim2[2]={ 1,1 };
		mwSize   N=opt->N;
		mwSize   M=opt->M;
		mwSize   L=opt->L;
		mwSize 	 dim3[3]={ N,2,M };
		mwSize 	 dim4[4]={ N,2,M };
		out=mxCreateStructArray(2,dim2,fieldNumber,fieldNames);
		mxPointer=mxCreateNumericMatrix(N,1,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"time",mxPointer); 	matOutput->time=mxGetData(mxPointer);
		dim2[0]=M; dim2[1]=L;
		mxPointer=mxCreateNumericArray(2,dim2,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sN",mxPointer); 	matOutput->sN=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(2,dim2,mxSINGLE_CLASS,mxREAL);		mxSetField(out,0,"tN",mxPointer);	matOutput->tN=mxGetData(mxPointer);
		dim3[0]=opt->maxKnotNum_Season+1; dim3[1]=M; dim3[2]=L;
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL); 	mxSetField(out,0,"sNProb",mxPointer); 	matOutput->sNProb=mxGetData(mxPointer);
		dim3[0]=opt->maxKnotNum_Trend+1; dim3[1]=M; dim3[2]=L;
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);	mxSetField(out,0,"tNProb",mxPointer);	matOutput->tNProb=mxGetData(mxPointer);
		dim3[0]=N; dim3[1]=M; dim3[2]=L;
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);		mxSetField(out,0,"sProb",mxPointer);	matOutput->sProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);		mxSetField(out,0,"tProb",mxPointer); matOutput->tProb=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);		mxSetField(out,0,"s",mxPointer); 	matOutput->s=mxGetData(mxPointer);
		dim4[0]=N; dim4[1]=2; dim4[2]=M; dim4[3]=L;
		if (opt->computeCredible) mxPointer=mxCreateNumericArray(4,dim4,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"sCI",mxPointer),matOutput->sCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);		mxSetField(out,0,"sSD",mxPointer); 	matOutput->sSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"t",mxPointer); 	matOutput->t=mxGetData(mxPointer);
		if (opt->computeCredible) 	mxPointer=mxCreateNumericArray(4,dim4,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"tCI",mxPointer),matOutput->tCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"tSD",mxPointer); 	matOutput->tSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);  		mxSetField(out,0,"b",mxPointer); 	matOutput->b=mxGetData(mxPointer);
		if (opt->computeCredible)  mxPointer=mxCreateNumericArray(4,dim4,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bCI",mxPointer),matOutput->bCI=mxGetData(mxPointer);
		mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"bSD",mxPointer); 	matOutput->bSD=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(M,L,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"marg_lik",mxPointer); matOutput->marg_lik=mxGetData(mxPointer);
		mxPointer=mxCreateNumericMatrix(M,L,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"sig2",mxPointer);    matOutput->sig2=mxGetData(mxPointer);
		if (opt->computeSlopeSign) 	    mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"bsign",mxPointer),matOutput->bsign=mxGetData(mxPointer);
		if (opt->computeHarmonicOrder)  mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"horder",mxPointer),matOutput->horder=mxGetData(mxPointer);
		if (opt->computeTrendOrder)     mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL),mxSetField(out,0,"torder",mxPointer),matOutput->torder=mxGetData(mxPointer);
		if (opt->computeChangepoints)
		{
			dim3[0]=opt->maxKnotNum_Season; dim3[1]=M; dim3[2]=L;
			mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL); 		mxSetField(out,0,"scp",mxPointer);    matOutput->scp=mxGetData(mxPointer);
			dim4[0]=opt->maxKnotNum_Season; dim4[1]=2; dim4[2]=M; dim4[3]=L;
			mxPointer=mxCreateNumericArray(4,dim4,mxSINGLE_CLASS,mxREAL); 						mxSetField(out,0,"scpCI",mxPointer);  matOutput->scpCI=mxGetData(mxPointer);
			dim3[0]=opt->maxKnotNum_Trend; dim3[1]=M; dim3[2]=L;
			mxPointer=mxCreateNumericArray(3,dim3,mxSINGLE_CLASS,mxREAL);  		mxSetField(out,0,"tcp",mxPointer);    matOutput->tcp=mxGetData(mxPointer);
			dim4[0]=opt->maxKnotNum_Trend; dim4[1]=2; dim4[2]=M; dim4[3]=L;
			mxPointer=mxCreateNumericArray(4,dim4,mxSINGLE_CLASS,mxREAL); 			mxSetField(out,0,"tcpCI",mxPointer);  matOutput->tcpCI=mxGetData(mxPointer);
		}
		fill_1_to_N(matOutput->time,N);
		r_ippsMulC_32f_I(opt->timeInterval,matOutput->time,(const int)N);
		r_ippsSubC_32f_I(-(opt->startTime - opt->timeInterval),matOutput->time,(const int)N);
		return out;
	}
}
#endif
ENABLE_MANY_WARNINGS
