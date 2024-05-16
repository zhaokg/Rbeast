#pragma once
#include "abc_datatype.h"
#include "abc_common.h"
#include "beastv2_header.h"

/****************beastv2_in_args******************/
extern int  BEAST2_GetArgs(VOID_PTR prhs[], int nrhs, A(OPTIONS_PTR) opt);
extern void BEAST2_DeallocateTimeSeriesIO(BEAST2_IO_PTR T);


/****************beastv2_in_readts******************/
//Need to check to make sure Xnewterm/GlobalMEMBuf is big enough to handle ts aggregration. 
//Raw ts may be longer than the aggregated ts but Xnewterm is allocated based on the length 
// of the aggegrate ts			
extern void  BEAST2_fetch_timeSeries(A(YINFO_PTR)  yInfo, int pixelIndex, F32PTR GlobalMEMBuf, A(IO_PTR)  io);
extern I08 BEAST2_preprocess_timeSeries(A(YINFO_PTR)  yInfo, BEAST2_BASIS_PTR basis, F32PTR Xtmp, BEAST2_OPTIONS_PTR opt);

/****************beastv2_out_allocmem******************/
extern void* BEAST2_Output_AllocMEM(BEAST2_OPTIONS_PTR  opt);
extern void  BEAST2_Result_FillMEM(BEAST2_RESULT_PTR  result, BEAST2_OPTIONS_PTR  opt, const F32 nan);
extern void  BEAST2_Result_AllocMEM(BEAST2_RESULT_PTR  result, BEAST2_OPTIONS_PTR  opt, MemPointers* _restrict MEM);

/****************beastv2_out_print******************/
extern void BEAST2_print_options(BEAST2_OPTIONS_PTR  opt);

/****************beastv2_out_write******************/
extern void  BEAST2_WriteOutput(BEAST2_OPTIONS_PTR opt, BEAST2_RESULT_PTR result, I64 pixelIndex);

