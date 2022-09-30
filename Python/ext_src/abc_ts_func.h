#pragma once
#include "abc_datatype.h"

extern void preCalc_terms_season(F32PTR SEASON_TERMS, F32PTR SEASON_SQR_CSUM, F32PTR SCALE_FACTOR, int N, F32 PERIOD, int maxSeasonOrder);
extern void preCalc_terms_trend(F32PTR TREND_TERMS, F32PTR INV_SQR, int N, int maxTrendOrder);
extern void preCalc_XmarsTerms_extra(F32PTR COEFF_A, F32PTR COEFF_B, I32 N);
extern void preCalc_XmarsTerms_extra_fmt3(F32PTR COEFF_A, F32PTR COEFF_B, I32 N);
extern void preCalc_scale_factor(F32PTR sclFactor, I32 N, I32 maxKnotNum, I32 minSepDist, F32PTR mem1, F32PTR mem2);
extern void cvt_knotList_to_bincode(U08PTR  good, I32 N, U16 minSepDist, U16PTR knotList, I64 knotNum);

I32 tsAggegrationPrepare(
	F32PTR oldTime, I32 Nold, F32 dT, I32PTR* SortedTimeIdx, I32PTR* NumPtsPerSeg,
	I32* startIdxOfFirsInterval, F32* startTime);

void tsAggegrationPerform(F32PTR RegularTS, I32 Nnew, F32PTR IrregularTS, 
	I32 Nold, I32PTR NumPerSeg, I32PTR SorteTimeIdx);
void tsRemoveNaNs(F32PTR x, int N);