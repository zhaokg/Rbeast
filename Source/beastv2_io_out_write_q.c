#include "abc_000_warning.h"

#include "abc_001_config.h"


#include <stdio.h>

//#include "abc_common.h"
#include "abc_ide_util.h"
#include "abc_ts_func.h"
#include "abc_common.h"  //WriteF32ArrayToStrideMEM
#include "beastv2_io.h"

static void GetOutputOffsetStride(A(IO_PTR) io, I64 idx, I64 Nvec, I64* pStride, I64* pOffset)
{
	//index should be 1-based.
	I64 stride, offset;
	if (io->ndim == 1)       // A 1d vcetor
		stride = 1L, offset = (idx - 1) * Nvec;
	else if (io->ndim == 2L)  // A 2 D mat input
	{
		I32 nPixel = io->meta.whichDimIsTime==1? 
							io->dims[1]:  /*whichDim=1*/
			                io->dims[0]; /*whichDim=2*/

		if(io->out.whichDimIsTime==1)
			stride = 1L, offset = (idx - 1) * Nvec;
		else {
			stride = nPixel;
			offset = (idx - 1);
		}
	}
	else if (io->ndim == 3L)  // A 3D stack input
	{
		int  ROW, COL;
		switch (io->meta.whichDimIsTime) {
		case 1:
			ROW = io->dims[1];
			COL = io->dims[2];
			break;
		case 2:
			ROW = io->dims[0];
			COL = io->dims[2];
			break;
		case 3:
			ROW = io->dims[0],
			COL = io->dims[1];
		}

		switch (io->out.whichDimIsTime) {
		case 1:
			stride = 1L;
			offset = (idx - 1) * Nvec;
			break;
		case 2: {
			int c = (idx - 1) / ROW;
			int r = idx - c * ROW;
			c = c + 1;
			stride = ROW;
			offset = (c - 1) * (Nvec * ROW) + r - 1;
			break;
		}
		case 3: {
			int c = (idx - 1) / ROW;
			int r = idx - c * ROW;
			c = c + 1;
			stride = ROW * COL;
			offset = (c - 1) * ROW + r - 1;
			break;
		}

		} // switch ( io->whichDimIsTime)	 

	} // (io->ndim == 3L)

	*pStride = stride;
	*pOffset = offset;
}

void  BEAST2_WriteOutput(A(OPTIONS_PTR) opt, A(RESULT_PTR) result, I64 pixelIndex)
{
	A(IO_PTR)       io =  &opt->io;
	A(RESULT_PTR)   mat = io->out.result;

	DATA_TYPE datType = io->out.dtype;
	const I64 N		  = io->N;
	const I64 seasonMaxKnotNum  = opt->prior.seasonMaxKnotNum;
	const I64 trendMaxKnotNum	= opt->prior.trendMaxKnotNum;
	const I64 outlierMaxKnotNum = opt->prior.outlierMaxKnotNum;

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	I64 len, offset, stride;
	len = 1, offset = pixelIndex - 1, stride = 1; //TODO: changed from 0 to 1
	WriteF32ArraryToStrideMEM(result->marg_lik, mat->marg_lik,	len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->sig2,		mat->sig2,		len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->R2,		mat->R2,		len, stride, offset, datType);
	WriteF32ArraryToStrideMEM(result->RMSE,		mat->RMSE,		len, stride, offset, datType);

	if (opt->extra.dumpInputData) {
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		WriteF32ArraryToStrideMEM(result->data, mat->data, len, stride, offset, datType);
	}

	#define  _(x)             WriteF32ArraryToStrideMEM((F32PTR)result->x,  mat->x,  len,  stride, offset, datType)
	#define  _2(x,y)          _(x), _(y)
	#define  _3(x,y,z)        _2(x,y), _(z)
	#define  _4(x,y,z,v)      _3(x,y,z), _(v)
	#define  _5(x,y,z,v,v1)   _4(x,y,z,v), _(v1)

	if (hasSeasonCmpnt) {

		len = 1, offset = pixelIndex - 1, stride = 0;
		_(sncp);
		_(sncp_median);
		_(sncp_mode);
		_(sncp_pct90);
		_(sncp_pct10);
		len = (seasonMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(sncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(scpOccPr, sY, sSD);
		if (opt->extra.computeSeasonOrder) _(sorder);
		if (opt->extra.computeSeasonAmp)   _2(samp, sampSD);
		if (opt->extra.computeCredible) {
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(sCI);
		}
		if (opt->extra.computeSeasonChngpt) {
			len = seasonMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_3(scp, scpPr, scpAbruptChange);
			len = seasonMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(scpCI);
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////

		len = 1, offset = pixelIndex - 1, stride = 0;
		_(tncp);
		_(tncp_median);
		_(tncp_mode);
		_(tncp_pct90);
		_(tncp_pct10);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(tncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(tcpOccPr, tY, tSD);
		if (opt->extra.computeTrendOrder)   _(torder);
		if (opt->extra.computeTrendSlope)   _4(tslp, tslpSD, tslpSgnPosPr, tslpSgnZeroPr);
		if (opt->extra.computeCredible) {
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(tCI);
		}
		if (opt->extra.computeTrendChngpt) {
			len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_3(tcp, tcpPr, tcpAbruptChange);

			len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
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
		len = (outlierMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(oncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(ocpOccPr, oY, oSD);
		if (opt->extra.computeCredible) {
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(oCI);
		}
		if (opt->extra.computeOutlierChngpt) {
			len = outlierMaxKnotNum;    GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_2(ocp, ocpPr);

			len = outlierMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(ocpCI);
		}
	}

	//////////////////////////////////////////////////////////////
	if (opt->extra.tallyPosNegSeasonJump && hasSeasonCmpnt ) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(spos_ncp, sneg_ncp);
		len = (seasonMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_ncpPr, sneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_cpOccPr, sneg_cpOccPr);

		len = seasonMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(spos_cp, spos_cpPr, spos_cpAbruptChange);
		_3(sneg_cp, sneg_cpPr, sneg_cpAbruptChange);
		len = seasonMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_cpCI, sneg_cpCI);
	}

	if (opt->extra.tallyPosNegTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tpos_ncp, tneg_ncp);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_ncpPr, tneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_cpOccPr, tneg_cpOccPr);

		len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(tpos_cp, tpos_cpPr, tpos_cpAbruptChange);
		_3(tneg_cp, tneg_cpPr, tneg_cpAbruptChange);
		len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_cpCI, tneg_cpCI);
	}

	if (opt->extra.tallyIncDecTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tinc_ncp, tdec_ncp);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_ncpPr, tdec_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_cpOccPr, tdec_cpOccPr);

		len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(tinc_cp, tinc_cpPr, tinc_cpAbruptChange);
		_3(tdec_cp, tdec_cpPr, tdec_cpAbruptChange);
		len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_cpCI, tdec_cpCI);

	}

	if (opt->extra.tallyPosNegOutliers && hasOutlierCmpnt) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(opos_ncp, oneg_ncp);
		len = (outlierMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_ncpPr, oneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cpOccPr, oneg_cpOccPr);

		len = outlierMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cp, opos_cpPr);
		_2(oneg_cp, oneg_cpPr);
		len = outlierMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cpCI, oneg_cpCI);
	}

}

void  MR_WriteOutput(A(OPTIONS_PTR) opt, A(RESULT_PTR) result, I64 pixelIndex)
{
	A(IO_PTR)       io =  &opt->io;
	A(RESULT_PTR)   mat = io->out.result;
 
	DATA_TYPE datType = io->out.dtype;
	const I64 N		  = io->N;
	const I64 q        = io->q;
	const I64 seasonMaxKnotNum = opt->prior.seasonMaxKnotNum;
	const I64 trendMaxKnotNum	= opt->prior.trendMaxKnotNum;
	const I64 outlierMaxKnotNum = opt->prior.outlierMaxKnotNum;

	I08 hasSeasonCmpnt  = opt->prior.basisType[0] == SEASONID || opt->prior.basisType[0] == DUMMYID || opt->prior.basisType[0] == SVDID;
	I08 hasOutlierCmpnt = opt->prior.basisType[opt->prior.numBasis - 1] == OUTLIERID;
	I08 hasTrendCmpnt   = 1;
	I08 hasAlways       = 1;

	I64 len, offset, stride;
	len = 1, offset = pixelIndex - 1, stride = 1; 
	WriteF32ArraryToStrideMEM(result->marg_lik, mat->marg_lik,	len, stride, offset, datType);
 
	len = q*q, offset = (pixelIndex - 1)*q*q, stride = 1;
	WriteF32ArraryToStrideMEM(result->sig2,		mat->sig2,		len, stride, offset, datType);

	for (I32 i=0; i < q; ++i) {
		len = 1, offset = (pixelIndex - 1) * q, stride = 1;
		WriteF32ArraryToStrideMEM(result->R2+i, mat[i].R2, len, stride, offset, datType);
	}

	for (I32 i = 0; i < q; ++i) {
		len = 1, offset = (pixelIndex - 1) * q, stride = 1;
		WriteF32ArraryToStrideMEM(result->RMSE + i, mat[i].RMSE, len, stride, offset, datType);
	}
 
	if (opt->extra.dumpInputData) {	
		for (I32 i = 0; i < q; ++i) {
			len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			WriteF32ArraryToStrideMEM(result->data+i*N, mat[i].data, len, stride, offset, datType);
		}
	}

	#define  _(x)          WriteF32ArraryToStrideMEM((F32PTR)result->x,  mat->x,  len,  stride, offset, datType)
	#define  _2(x,y)        _(x), _(y)
	#define  _3(x,y,z)      _2(x,y), _(z)
	#define  _4(x,y,z,v)    _3(x,y,z), _(v)
	#define  _5(x,y,z,v,v1) _4(x,y,z,v), _(v1)
 
	/*****************************************************************/
	//         Write outputs for SEASON
	/*****************************************************************/
	if (hasSeasonCmpnt) {

		len = 1, offset = pixelIndex - 1, stride = 0;
		_(sncp); 
		_(sncp_median);
		_(sncp_mode);
		_(sncp_pct90);
		_(sncp_pct10);
		len = (seasonMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(sncpPr);    //#############
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(scpOccPr);   //#############

		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		for (I32 i = 0; i < q; ++i) {			
			WriteF32ArraryToStrideMEM(result->sY  + N*i,  mat[i].sY,  len, stride, offset, datType);
			WriteF32ArraryToStrideMEM(result->sSD + N*i,  mat[i].sSD, len, stride, offset, datType);
		}
		if (opt->extra.computeSeasonOrder) {
			len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(sorder);
		}
		if (opt->extra.computeSeasonAmp) {
			len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			for (I32 i = 0; i < q; ++i) {
				WriteF32ArraryToStrideMEM(result->samp + N * i, mat[i].samp, len, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->sampSD + N * i, mat[i].sampSD, len, stride, offset, datType);
			}
		}   
		if (opt->extra.computeCredible) {
			// MRBEAST
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			for (I32 i = 0; i < q; ++i) {				
				//WriteF32ArraryToStrideMEM(result->sCI + len*i, mat[i].sCI, len, stride, offset, datType);				
				WriteF32ArraryToStrideMEM(result->sCI +      N*i, mat[i].sCI, N, stride, offset, datType);				
				WriteF32ArraryToStrideMEM(result->sCI + N*q+ N*i, mat[i].sCI, N, stride, offset+N*stride, datType);
			}			
		}
		if (opt->extra.computeSeasonChngpt) {
			len = seasonMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_2(scp, scpPr);
			for (int i = 0; i < q; ++i) {
				WriteF32ArraryToStrideMEM(result->scpAbruptChange + len * i, mat[i].scpAbruptChange, len, stride, offset, datType);
			}
			len = seasonMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(scpCI);
		}

	}

	/*****************************************************************/
	//         Write outputs for TREND
	/*****************************************************************/
	if (hasTrendCmpnt) {

		len = 1, offset = pixelIndex - 1, stride = 0;
		_(tncp);
		_(tncp_median);
		_(tncp_mode);
		_(tncp_pct90);
		_(tncp_pct10);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(tncpPr);         //#############
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(tcpOccPr);       //#############  

		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		for (I32 i = 0; i < q; ++i) {			
			WriteF32ArraryToStrideMEM(result->tY + N * i,  mat[i].tY, len, stride, offset, datType);
			WriteF32ArraryToStrideMEM(result->tSD + N * i, mat[i].tSD, len, stride, offset, datType);
		}		
		if (opt->extra.computeTrendOrder)   _(torder);
		if (opt->extra.computeTrendSlope) {
			len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			for (I32 i = 0; i < q; ++i) {
				WriteF32ArraryToStrideMEM(result->tslp + N * i,       mat[i].tslp, len, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->tslpSD + N * i,     mat[i].tslpSD, len, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->tslpSgnPosPr + N * i, mat[i].tslpSgnPosPr, len, stride, offset, datType);		
				WriteF32ArraryToStrideMEM(result->tslpSgnZeroPr + N * i, mat[i].tslpSgnZeroPr, len, stride, offset, datType);
			}			
		}   		
		if (opt->extra.computeCredible) {
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			for (I32 i = 0; i < q; ++i) {
				//WriteF32ArraryToStrideMEM(result->tCI + len * i, mat[i].tCI, len, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->tCI + N * i,         mat[i].tCI, N, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->tCI + N * q + N * i, mat[i].tCI, N, stride, offset + N * stride, datType);
			}	
		}
		
		if (opt->extra.computeTrendChngpt) {
			len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_2(tcp, tcpPr);
			for (int i = 0; i < q; ++i) {
				WriteF32ArraryToStrideMEM(result->tcpAbruptChange + len * i, mat[i].tcpAbruptChange, len, stride, offset, datType);
			}
			

			len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(tcpCI);
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	if (hasOutlierCmpnt) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_(oncp); 
		_(oncp_median);
		_(oncp_mode); 
		_(oncp_pct90);
		_(oncp_pct10);
		len = (outlierMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(oncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_(ocpOccPr);

		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		for (I32 i = 0; i < q; ++i) {
			WriteF32ArraryToStrideMEM(result->oY + N * i, mat[i].oY, len, stride, offset, datType);
			WriteF32ArraryToStrideMEM(result->oSD + N * i, mat[i].oSD, len, stride, offset, datType);
		}

		if (opt->extra.computeCredible) {
			len = N * 2, GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			for (I32 i = 0; i < q; ++i) {
				//WriteF32ArraryToStrideMEM(result->oCI + len * i, mat[i].oCI, len, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->oCI + N * i,         mat[i].oCI, N, stride, offset, datType);
				WriteF32ArraryToStrideMEM(result->oCI + N * q + N * i, mat[i].oCI, N, stride, offset + N * stride, datType);
			}
		}
		if (opt->extra.computeOutlierChngpt) {
			len = outlierMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_2(ocp, ocpPr);

			len = outlierMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
			_(ocpCI);
		}
	}

	//////////////////////////////////////////////////////////////
	if (opt->extra.tallyPosNegSeasonJump && hasSeasonCmpnt ) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(spos_ncp, sneg_ncp);
		len = (seasonMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_ncpPr, sneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_cpOccPr, sneg_cpOccPr);

		len = seasonMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(spos_cp, spos_cpPr, spos_cpAbruptChange);
		_3(sneg_cp, sneg_cpPr, sneg_cpAbruptChange);
		len = seasonMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(spos_cpCI, sneg_cpCI);
	}

	if (opt->extra.tallyPosNegTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tpos_ncp, tneg_ncp);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_ncpPr, tneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_cpOccPr, tneg_cpOccPr);

		len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(tpos_cp, tpos_cpPr, tpos_cpAbruptChange);
		_3(tneg_cp, tneg_cpPr, tneg_cpAbruptChange);
		len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tpos_cpCI, tneg_cpCI);
	}

	if (opt->extra.tallyIncDecTrendJump) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(tinc_ncp, tdec_ncp);
		len = (trendMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_ncpPr, tdec_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_cpOccPr, tdec_cpOccPr);

		len = trendMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_3(tinc_cp, tinc_cpPr, tinc_cpAbruptChange);
		_3(tdec_cp, tdec_cpPr, tdec_cpAbruptChange);
		len = trendMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(tinc_cpCI, tdec_cpCI);

	}

	if (opt->extra.tallyPosNegOutliers && hasOutlierCmpnt) {
		len = 1, offset = pixelIndex - 1, stride = 0;
		_2(opos_ncp, oneg_ncp);
		len = (outlierMaxKnotNum + 1);  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_ncpPr, oneg_ncpPr);
		len = N;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cpOccPr, oneg_cpOccPr);

		len = outlierMaxKnotNum;  GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cp, opos_cpPr);
		_2(oneg_cp, oneg_cpPr);
		len = outlierMaxKnotNum * 2; GetOutputOffsetStride(io, pixelIndex, len, &stride, &offset);
		_2(opos_cpCI, oneg_cpCI);

	}

}
 
 
#include "abc_000_warning.h"