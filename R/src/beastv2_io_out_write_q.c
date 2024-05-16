#include "abc_000_warning.h"

#include "abc_001_config.h"


#include <stdio.h>

//#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_common.h"  //WriteF32ArrayToStrideMEM
#include "abc_vec.h"
#include "beastv2_io.h"

static void __convert_index_to_outdatasubs3(BEAST2_IO* io, int index, int dims[],int subs3[]) {

	int subs2[2];
	ind2sub(io->imgdims, 2L, index, subs2);

	int rows = io->dims[io->rowdim - 1];
	int cols = io->dims[io->coldim - 1];

	if (io->out.whichDimIsTime == 1) {
		subs3[0] = 1;  // can fill any value
		subs3[1] = subs2[0]; 
		subs3[2] = subs2[1]; 
		
		dims[1] = rows;
		dims[2] = cols;
	}
	else if (io->out.whichDimIsTime == 2) {
		subs3[1] = 1;  // can fill any value
		subs3[0] = subs2[0]; 
		subs3[2] = subs2[1]; 

		dims[0] = rows;
		dims[2] = cols;
	}
	else if (io->out.whichDimIsTime == 3) {
		subs3[2] = 1;  // can fill any value
		subs3[0] = subs2[0];
		subs3[1] = subs2[1];

		dims[0] = rows;
		dims[1] = cols;
	}

}

 

void  BEAST2_WriteOutput(A(OPTIONS_PTR) opt, A(RESULT_PTR) result, I64 pixelIndex) {

	BEAST2_IO_PTR       io  =  &opt->io;
	BEAST2_RESULT_PTR   mat = io->out.result;

	DATA_TYPE datType = io->out.dtype;
	const I64 N		  = io->N;
	const I64 q       = io->q;
	const I64 seasonMaxKnotNum  = opt->prior.seasonMaxKnotNum;
	const I64 trendMaxKnotNum	= opt->prior.trendMaxKnotNum;
	const I64 outlierMaxKnotNum = opt->prior.outlierMaxKnotNum;

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	int outTimeDim0 = io->out.whichDimIsTime-1L;
	int subs3[3];
	int dims[3];
	__convert_index_to_outdatasubs3(io, pixelIndex, dims, subs3);
	

	#define  _(x)           f32_to_strided_mem((F32PTR)result->x,  mat->x,  len,  stride, offset, datType)
	#define  _2(x,y)        _(x), _(y)
	#define  _3(x,y,z)      _2(x,y), _(z)
	#define  _4(x,y,z,v)    _3(x,y,z), _(v)
	#define  _5(x,y,z,v,w)  _4(x,y,z,v), _(w)

#define GET_OFFSET_STRIDE(Ntime)  do {                                                       \
			  dims[outTimeDim0] = Ntime;                                                  \
			  len=ndarray_get1d_stride_offset(dims, 3L, subs3, outTimeDim0 + 1L, &stride, &offset); \
			} while(0);

	I64 len, offset, stride;

	len = 1, offset = pixelIndex - 1, stride = 1; //TODO: changed from 0 to 1
	f32_to_strided_mem(result->marg_lik, mat->marg_lik,	len, stride, offset, datType);
 
	GET_OFFSET_STRIDE(q * q);
	f32_to_strided_mem(result->sig2,		mat->sig2,		len, stride, offset, datType);

    len = 1, offset = (pixelIndex - 1) * 1, stride = 1; 
 	for (I32 i=0; i < q; ++i) {		
		f32_to_strided_mem(result->R2   + i,  mat[i].R2,   len, stride, offset, datType);
        f32_to_strided_mem(result->RMSE + i,  mat[i].RMSE, len, stride, offset, datType);		
	}

	if (opt->extra.dumpInputData) {
		GET_OFFSET_STRIDE(N);
        for (I32 i = 0; i < q; ++i) {
			f32_to_strided_mem(result->data+i*N, mat[i].data, len, stride, offset, datType);
		}
	}


	/*****************************************************************/
	//         Write outputs for SEASON
	/*****************************************************************/
	if (hasSeasonCmpnt) {

		len = 1, offset = pixelIndex - 1, stride = 0;
		_5(sncp, sncp_median, sncp_mode,sncp_pct90,sncp_pct10);
				
		GET_OFFSET_STRIDE(seasonMaxKnotNum + 1);
		_(sncpPr);
 
		GET_OFFSET_STRIDE(N); 
		_(scpOccPr);
        for (I32 i = 0; i < q; ++i) {			
			f32_to_strided_mem(result->sY  + N*i,  mat[i].sY,  len, stride, offset, datType);
			f32_to_strided_mem(result->sSD + N*i,  mat[i].sSD, len, stride, offset, datType);
		}
		if (opt->extra.computeSeasonOrder) {
			_(sorder);
		}
		if (opt->extra.computeSeasonAmp)   {
        	for (I32 i = 0; i < q; ++i) {
				f32_to_strided_mem(result->samp + N * i,   mat[i].samp, len, stride, offset, datType);
				f32_to_strided_mem(result->sampSD + N * i, mat[i].sampSD, len, stride, offset, datType);
			}
        }

		if (opt->extra.computeCredible) {
			GET_OFFSET_STRIDE(N*2); 
			for (I32 i = 0; i < q; ++i) {				
				//f32_to_strided_mem(result->sCI + len*i, mat[i].sCI, len, stride, offset, datType);				
				f32_to_strided_mem(result->sCI +      N*i, mat[i].sCI, N, stride, offset, datType);				
				f32_to_strided_mem(result->sCI + N*q+ N*i, mat[i].sCI, N, stride, offset+N*stride, datType);
			}			
		}
		if (opt->extra.computeSeasonChngpt) {
			GET_OFFSET_STRIDE(seasonMaxKnotNum); 
			_2(scp, scpPr);
			for (int i = 0; i < q; ++i) {
				f32_to_strided_mem(result->scpAbruptChange + len * i, mat[i].scpAbruptChange, len, stride, offset, datType);
			}
			GET_OFFSET_STRIDE(seasonMaxKnotNum*2); 
			_(scpCI);
		}

	}

	/*****************************************************************/
	//         Write outputs for TREND
	/*****************************************************************/

		len = 1, offset = pixelIndex - 1, stride = 0;
		_(tncp);
		_(tncp_median);
		_(tncp_mode);
		_(tncp_pct90);
		_(tncp_pct10);

		GET_OFFSET_STRIDE(trendMaxKnotNum+1); 
		_(tncpPr);
				
		GET_OFFSET_STRIDE(N); 
		_(tcpOccPr);
        for (I32 i = 0; i < q; ++i) {			
			f32_to_strided_mem(result->tY  + N * i,  mat[i].tY, len, stride, offset, datType);
			f32_to_strided_mem(result->tSD + N * i, mat[i].tSD, len, stride, offset, datType);
		}	
		if (opt->extra.computeTrendOrder)   _(torder);
		if (opt->extra.computeTrendSlope)  {
			for (I32 i = 0; i < q; ++i) {
				f32_to_strided_mem(result->tslp + N * i,          mat[i].tslp, len, stride, offset, datType);
				f32_to_strided_mem(result->tslpSD + N * i,        mat[i].tslpSD, len, stride, offset, datType);
				f32_to_strided_mem(result->tslpSgnPosPr + N * i,  mat[i].tslpSgnPosPr, len, stride, offset, datType);		
				f32_to_strided_mem(result->tslpSgnZeroPr + N * i, mat[i].tslpSgnZeroPr, len, stride, offset, datType);
			}	
        }

		if (opt->extra.computeCredible) {
			GET_OFFSET_STRIDE(N*2);
			for (I32 i = 0; i < q; ++i) {
				//f32_to_strided_mem(result->tCI + len * i, mat[i].tCI, len, stride, offset, datType);
				f32_to_strided_mem(result->tCI + N * i,         mat[i].tCI, N, stride, offset, datType);
				f32_to_strided_mem(result->tCI + N * q + N * i, mat[i].tCI, N, stride, offset + N * stride, datType);
			}

			if (opt->extra.computeTrendSlope) {
				for (I32 i = 0; i < q; ++i) {
					//f32_to_strided_mem(result->tCI + len * i, mat[i].tCI, len, stride, offset, datType);
					f32_to_strided_mem(result->tslpCI + N * i,         mat[i].tslpCI, N, stride, offset, datType);
					f32_to_strided_mem(result->tslpCI + N * q + N * i, mat[i].tslpCI, N, stride, offset + N * stride, datType);
				}			
			}

		}

		if (opt->extra.computeTrendChngpt) {
			GET_OFFSET_STRIDE(trendMaxKnotNum); 
			_2(tcp, tcpPr);
			for (int i = 0; i < q; ++i) {
				f32_to_strided_mem(result->tcpAbruptChange + len * i, mat[i].tcpAbruptChange, len, stride, offset, datType);
			}
 
			GET_OFFSET_STRIDE(trendMaxKnotNum*2); 
			_(tcpCI);
		}


	//////////////////////////////////////////////////////////////////////////////////////////////
	if (hasOutlierCmpnt) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_(oncp);
		_(oncp_median);
		_(oncp_mode);
		_(oncp_pct90);
		_(oncp_pct10);

		GET_OFFSET_STRIDE(outlierMaxKnotNum + 1); 
		_(oncpPr);

		GET_OFFSET_STRIDE(N); 
		_(ocpOccPr);
		for (I32 i = 0; i < q; ++i) {
			f32_to_strided_mem(result->oY + N * i, mat[i].oY, len, stride, offset, datType);
			f32_to_strided_mem(result->oSD + N * i, mat[i].oSD, len, stride, offset, datType);
		}

		if (opt->extra.computeCredible) {
			GET_OFFSET_STRIDE(N*2);
			for (I32 i = 0; i < q; ++i) {
				//f32_to_strided_mem(result->oCI + len * i, mat[i].oCI, len, stride, offset, datType);
				f32_to_strided_mem(result->oCI + N * i,         mat[i].oCI, N, stride, offset, datType);
				f32_to_strided_mem(result->oCI + N * q + N * i, mat[i].oCI, N, stride, offset + N * stride, datType);
			}
		}
		if (opt->extra.computeOutlierChngpt) {
			GET_OFFSET_STRIDE(outlierMaxKnotNum);  
			_2(ocp, ocpPr);

			GET_OFFSET_STRIDE(outlierMaxKnotNum*2);
			_(ocpCI);
		}
	}

	//////////////////////////////////////////////////////////////
	if (opt->extra.tallyPosNegSeasonJump && hasSeasonCmpnt ) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(spos_ncp, sneg_ncp);
 
		GET_OFFSET_STRIDE(seasonMaxKnotNum+1); 
		_2(spos_ncpPr, sneg_ncpPr);

		GET_OFFSET_STRIDE(N); 
		_2(spos_cpOccPr, sneg_cpOccPr);

		GET_OFFSET_STRIDE(seasonMaxKnotNum); 
		_3(spos_cp, spos_cpPr, spos_cpAbruptChange);
		_3(sneg_cp, sneg_cpPr, sneg_cpAbruptChange);

		GET_OFFSET_STRIDE(seasonMaxKnotNum*2); 
		_2(spos_cpCI, sneg_cpCI);
	}

	if (opt->extra.tallyPosNegTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tpos_ncp, tneg_ncp);

		GET_OFFSET_STRIDE(trendMaxKnotNum+1); 
		_2(tpos_ncpPr, tneg_ncpPr);

		GET_OFFSET_STRIDE(N);
		_2(tpos_cpOccPr, tneg_cpOccPr);

		GET_OFFSET_STRIDE(trendMaxKnotNum); 
		_3(tpos_cp, tpos_cpPr, tpos_cpAbruptChange);
		_3(tneg_cp, tneg_cpPr, tneg_cpAbruptChange);

		GET_OFFSET_STRIDE(trendMaxKnotNum*2);
		_2(tpos_cpCI, tneg_cpCI);
	}

	if (opt->extra.tallyIncDecTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tinc_ncp, tdec_ncp);

		GET_OFFSET_STRIDE(trendMaxKnotNum+1); 
		_2(tinc_ncpPr, tdec_ncpPr);

		GET_OFFSET_STRIDE(N); 
		_2(tinc_cpOccPr, tdec_cpOccPr);

		GET_OFFSET_STRIDE(trendMaxKnotNum); 
		_3(tinc_cp, tinc_cpPr, tinc_cpAbruptChange);
		_3(tdec_cp, tdec_cpPr, tdec_cpAbruptChange);

		GET_OFFSET_STRIDE(trendMaxKnotNum*2);
		_2(tinc_cpCI, tdec_cpCI);

	}

	if (opt->extra.tallyPosNegOutliers && hasOutlierCmpnt) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(opos_ncp, oneg_ncp);

		GET_OFFSET_STRIDE(outlierMaxKnotNum + 1);
		_2(opos_ncpPr, oneg_ncpPr);

		GET_OFFSET_STRIDE(N); 
		_2(opos_cpOccPr, oneg_cpOccPr);

		GET_OFFSET_STRIDE(outlierMaxKnotNum); 
		_2(opos_cp, opos_cpPr);
		_2(oneg_cp, oneg_cpPr);


		GET_OFFSET_STRIDE(outlierMaxKnotNum * 2); 
		_2(opos_cpCI, oneg_cpCI);
	}

}
 
 
#include "abc_000_warning.h"