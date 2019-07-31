
#pragma once
#include "abc_common.h"
#include "abc_blas_lapack_lib.h"
#include "abc_mem.h"
#include <inttypes.h>  
#include "beast_common.h"
extern  void preCompute_Xmars_terms(F32PTR SEASON_TERMS,F32PTR TREND_TERMS,int N,float PERIOD,int maxSeasonOrder,int maxTrendOrder);
extern  void preCompute_Xmars_terms_fast(F32PTR SEASON_TERMS,F32PTR SEASON_SQR_CSUM,F32PTR TREND_TERMS,F32PTR INV_SQR,F32PTR COEFF_A,F32PTR COEFF_B,
	int N,float PERIOD,int maxSeasonOrder,int maxTrendOrder);
extern  void preCompute_scale_factor(F32PTR scaleFactorSeason,F32PTR scaleFactorTrend,int32_t N,int maxKnotNum_Season,int maxKnotNum_Trend,int minSepDist_Season,int minSepDist_Trend,F32PTR mem1,F32PTR mem2);
extern  void fetch_next_time_series( YINFO * _restrict yInfo,void * _restrict yInputData,int idx,float *_restrict GlobalMEMBuf_1st,bool isSingleyInput,int32_t N,float omissionValue);
extern  void fetch_next_time_series1(YINFO * _restrict yInfo,void * _restrict yInputData,int idx,float *_restrict GlobalMEMBuf_1st,char inputType,int32_t N,float omissionValue);
extern  void fetch_next_time_series2(YINFO * _restrict yInfo,void * _restrict yInputData,int idx,F32PTR GlobalMEMBuf_1st,char inputType,int32_t N,float omissionValue,Options * _restrict opt);
extern  void convert_basis_both(struct BASIS * _restrict basis);
extern  void evaluate_basis_both(struct BASIS * _restrict basis,int32_t N,F32PTR  Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,F32PTR TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st);
extern  void evaluate_basis_both_fast(struct BASIS * _restrict basis,int32_t N,float *_restrict Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,F32PTR TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st);
extern  void evaluate_basis_both_trend_fast(struct BASIS * _restrict  basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,float* _restrict TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st);
extern  void findGoodKnotPositionFromBasis(struct BASIS * _restrict basis,uint8_t * _restrict goodS,uint8_t * _restrict goodT,int32_t N,uint16_t minSepDist_Season,uint16_t minSepDist_Trend);
extern  void allocate_single_output(RESULT * _restrict result,Options * _restrict opt,MemPointers * _restrict MEM);
extern  void print_options(Options * _restrict opt);
extern  void zero_out_result_output(Options * _restrict,RESULT *_restrict result);
extern  void nan_fill_result_output(Options * _restrict opt,RESULT *_restrict result);
extern void evaluate_basis_both_BIC(struct BASIS * _restrict basis,int32_t N,F32PTR Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,F32PTR TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st);
extern void convert_basis_both_trend(struct BASIS * _restrict basis);
extern void evaluate_basis_both_trend(struct BASIS * _restrict basis,int32_t N,float *_restrict Xt_mars,F32PTR Xt_zeroBackup,YINFO * _restrict pyInfo,F32PTR SEASON_TERMS,F32PTR TREND_TERMS,struct ModelPar * _restrict pmodelPar,F32PTR GlobalMEMBuf_1st);
extern void print_options_trend(Options * _restrict opt);
extern void allocate_single_output_trend(RESULT * _restrict result,Options * _restrict opt,MemPointers * _restrict MEM);
extern void zero_out_result_output_trend(Options * _restrict,RESULT *_restrict result);
