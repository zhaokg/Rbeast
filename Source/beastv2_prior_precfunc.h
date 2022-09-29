#include "abc_000_macro.h"

#include <math.h>
#include <string.h>
#include "beastv2_header.h"

#include "abc_rand.h"

#include "globalvars.h"
#include "abc_mcmc.h"
#include "abc_ts_func.h"

typedef struct PREC_FUNC {
	void (*IncreasePrecValues)(        BEAST2_MODEL_PTR model);
	void (*GetNumTermsPerPrecGrp)(     BEAST2_MODEL_PTR model);
	void (*GetXtXPrecDiag)(            BEAST2_MODEL_PTR model);
	void (*UpdateXtXPrec_nTermsPerGrp)(BEAST2_MODEL_PTR model, BEAST2_BASIS_PTR basis, NEWTERM_PTR new, NEWCOLINFO_PTR newcol);
	//F32  (*ComputeLog_det_prior)(BEAST2_MODEL_PTR model, I08PTR nTermsPerPrecGrp, I32 K);
	void (*ComputeMargLik)(            BEAST2_MODELDATA_PTR data, BEAST2_MODEL_PTR model,BEAST2_YINFO_PTR yInfo, BEAST2_HyperPar_PTR hyperPar);
	void (*ResamplePrecValues)(        BEAST2_MODEL_PTR model,    BEAST2_HyperPar* hyperPar,VOID_PTR stream);
	void (*chol_addCol)(               F32PTR Au,                 F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1);
} PREC_FUNCS;
void SetUpPrecFunctions(I08 precPriorType, I32 q, PREC_FUNCS* funcs);
