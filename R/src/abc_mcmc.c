#include <math.h>
#include "abc_000_warning.h"

#include "abc_mcmc.h"
#include "abc_blas_lapack_lib.h"
#include "abc_mat.h" //f32_transpose_inplace


static void PrepareCiInfo(CI_PARAM *cinfo, U32 nSamples,  MemPointers * MEM)
{
 	    U32 stripWidth;
		U32 nStrips;
 
		stripWidth     = (U32)ceil(sqrt((F64)nSamples));
		nStrips        = (U32)floor( nSamples / stripWidth);

		cinfo->SamplesPerStrip = MyALLOC((*MEM), nStrips, int, 0);
		cinfo->OffsetsPerStrip = MyALLOC((*MEM), nStrips, int, 0);

		U32 stripCumNum = 0;
		for (rU32 i = 1; i <= nStrips; i++)
			cinfo->SamplesPerStrip[i-1] = (i != nStrips) ? stripWidth : (nSamples-stripCumNum),
			cinfo->OffsetsPerStrip[i-1] = stripCumNum,
			stripCumNum                += stripWidth;

		cinfo->nSamples = nSamples;
		cinfo->nStrips  = nStrips;

}
static void AllocateMemForCISorting(CI_PARAM * cinfo, CI_RESULT *CI, int Nrowlength, int NUM_VARS,  MemPointers* MEM)
{
	
	U32  nSamples  = cinfo->nSamples;
	U32  nStrips   = cinfo->nStrips;
		
	for (I32 i = 0; i < NUM_VARS; i++) {
		
		U32  N = Nrowlength;

		CI[i].N = N;

		CI[i].CI95              = MyALLOC((*MEM), N*nSamples, F32, 0);
		CI[i].minIdxPerStrip    = MyALLOC((*MEM), N*nStrips, int, 0);
		CI[i].minValPerStrip    = MyALLOC((*MEM), N*nStrips, F32, 0);
		CI[i].whichStripHasMin  = MyALLOC((*MEM), N, int, 0);

		CI[i].CI05              = MyALLOC((*MEM), N*nSamples, F32, 0);
		CI[i].maxIdxPerStrip    = MyALLOC((*MEM), N*nStrips, int, 0);
		CI[i].maxValPerStrip    = MyALLOC((*MEM), N*nStrips, F32, 0);
		CI[i].whichStripHasMax  = MyALLOC((*MEM), N, int, 0);
	}
}

#define _inout_ 
void ConstructCIStruct(F32 alpahLevel, I32 MCMC_SAMPLES,  I32 N, I32 numCIVars, MemPointers* MEM,
	                  U08PTR _inout_  fastCIComputation,  CI_PARAM * _out_ ciInfo, CI_RESULT * _out_ CI)
{
		//subsampleFraction_x_INT32MAX=4294967295UL;

		F32 alpha   = (1.0 - alpahLevel) / 2.0;
		U32 nSamples;
		if (*fastCIComputation) {
			nSamples = 100;
			F64 totalSamplesFor100nSamples = (nSamples / alpha);
			F64 fractionOfSubSample        = totalSamplesFor100nSamples / (F32)MCMC_SAMPLES;
			if (fractionOfSubSample < 0.99f)
				//ciInfo->subsampleFraction_x_INT32MAX = (U32)(fractionOfSubSample *(F64)4294967295UL);
			    ciInfo->subsampleFraction_x_INT16MAX = (U16)(fractionOfSubSample*(F64)65535UL);
			else
				*fastCIComputation = 0;
		}
		if (!(*fastCIComputation)) {
			nSamples = (F32) MCMC_SAMPLES * alpha;
		}

		// Allocate the memory and calculte stripNumRows and   // stripOffets	   
		PrepareCiInfo(ciInfo, nSamples, MEM);

		// Allocate the memory for saving CI samples
		AllocateMemForCISorting(ciInfo, CI, N, numCIVars, MEM); 

 
}




static void CalcInitialCI(CI_PARAM* _restrict info, CI_RESULT* _restrict ci) {
// When enough samples are inserted, cacluate the inital min and max; this is 
// the first time that CI are computed
	
	I64     nStrips  = info->nStrips;
	I64     nSamples = info->nSamples;

	F32PTR CI95           = ci->CI95;
	F32PTR minValPerStrip = ci->minValPerStrip;
	I32PTR minIdxPerStrip = ci->minIdxPerStrip;


	F32PTR CI05 = ci->CI05;
	F32PTR maxValPerStrip = ci->maxValPerStrip;
	I32PTR maxIdxPerStrip = ci->maxIdxPerStrip;

	I64     N             = ci->N;
	for (I64 i = 0; i < N; i++)
	{ //ippsMaxIndx_32f(const Ipp64f* pSrc, int len, Ipp64f* pMax, int* pIndx);

		for (U32 j = 0; j < nStrips; j++) {
			r_ippsMinIndx_32f(CI95, info->SamplesPerStrip[j], minValPerStrip + j, minIdxPerStrip + j);
			CI95 += info->SamplesPerStrip[j];
		}
		r_ippsMinIndx_32f(minValPerStrip, nStrips, ci->result + i, ci->whichStripHasMin + i);

		minValPerStrip += nStrips;
		minIdxPerStrip += nStrips;

		///////////////////////////////////////////////////////////////////////////////
		for (U32 j = 0; j < nStrips; j++) {
			r_ippsMaxIndx_32f(CI05, info->SamplesPerStrip[j], maxValPerStrip + j, maxIdxPerStrip + j);
			CI05 += info->SamplesPerStrip[j];
		}
		r_ippsMaxIndx_32f(maxValPerStrip, nStrips, ci->result + N + i, ci->whichStripHasMax + i);

		maxValPerStrip += nStrips;
		maxIdxPerStrip += nStrips;
	}
}

 

void InsertNewRowToUpdateCI(CI_PARAM* _restrict info, CI_RESULT* _restrict ci)
{
	I64     N          = ci->N;
	F32PTR  newDataRow = ci->newDataRow;

	if (ci->samplesInserted < info->nSamples) {
		//Insert the inital nSamples rows 
		I64 offset = ci->samplesInserted * N;
		r_cblas_scopy(N, newDataRow, 1, ci->CI95 + offset, 1); 
		ci->samplesInserted++;

		U32 nSamples = info->nSamples;
		if (ci->samplesInserted == nSamples) {
			r_mkl_simatcopy('C', 'T', N, nSamples, 1, ci->CI95, N, nSamples);
			//make a copy of cred_upper into "cred_lower"
			r_cblas_scopy(N * nSamples, ci->CI95, 1L, ci->CI05, 1L);
			CalcInitialCI(info, ci);
		}

		return;
	}

	// No need to increate samplesInserted further because it won't be needed for the update code below
	// ci->samplesInserted++;

	F32PTR CI95           = ci->CI95;
	I32PTR minIdxPerStrip = ci->minIdxPerStrip;
	F32PTR minValPerStrip = ci->minValPerStrip;

	F32PTR CI05           = ci->CI05;
	I32PTR maxIdxPerStrip = ci->maxIdxPerStrip;
	F32PTR maxValPerStrip = ci->maxValPerStrip;

	
	I64     nStrips = info->nStrips;

	for (I64 i = 0; i < N; i++) {

		F32 newY = newDataRow[i];

		/******************************/
		//Upper credible limit	
		/******************************/
		if (newY > ci->result[i]) {

			I64     which;
			F32PTR  data;

			which = ci->whichStripHasMin[i];
			data = CI95 + info->OffsetsPerStrip[which];
			data[minIdxPerStrip[which]] = newY;

			/*Find the min value and index within the chosen strip*/
			rF32 minVal = data[0];
			rI32 minIdx = 0L;
			for (I64 j = 1L; j < (*info).SamplesPerStrip[which]; j++) {
				if (data[j] < minVal) { minVal = data[j]; minIdx = j; }
			}
			minValPerStrip[which] = minVal;
			minIdxPerStrip[which] = minIdx;

			/*Find the max value and index across all the strips*/
			minVal = minValPerStrip[0];
			minIdx = 0;
			for (I64 j = 1L; j < nStrips; j++) {
				if (minValPerStrip[j] < minVal) { minVal = minValPerStrip[j]; minIdx = j; }
			}
			ci->result[i] = minVal;
			ci->whichStripHasMin[i] = minIdx;
		} // if (newY > ci->result[i])

		CI95 += info->nSamples;
		minValPerStrip += info->nStrips;
		minIdxPerStrip += info->nStrips;


		/******************************/
		//Lower credible limit	
		/******************************/
		if (newY < ci->result[N + i]) {

			I64     which;
			F32PTR  data;

			which = (I64)ci->whichStripHasMax[i];
			data = CI05 + info->OffsetsPerStrip[which];
			data[maxIdxPerStrip[which]] = newY;

			/*Find the max value and index within the chosen strip*/
			F32 maxVal = data[0];
			I32 maxIdx = 0L;
			for (I64 j = 1L; j < (*info).SamplesPerStrip[which]; j++) {
				if (data[j] > maxVal) { maxVal = data[j]; maxIdx = j; }
			}
			maxValPerStrip[which] = maxVal;
			maxIdxPerStrip[which] = maxIdx;


			/*Find the max value and index across all the strips*/
			maxVal = maxValPerStrip[0];
			maxIdx = 0;
			for (I64 j = 1L; j < nStrips; j++) {
				if (maxValPerStrip[j] > maxVal) { maxVal = maxValPerStrip[j]; maxIdx = j; }
			}
			ci->result[i + N] = maxVal;
			ci->whichStripHasMax[i] = maxIdx;
		} //if (newY < ci->result[N + i]) 

		CI05           += info->nSamples;
		maxValPerStrip += info->nStrips;
		maxIdxPerStrip += info->nStrips;
	}
}



 
#include "abc_000_warning.h"

