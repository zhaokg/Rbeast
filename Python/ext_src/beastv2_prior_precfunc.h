#include "abc_000_macro.h"

#include <math.h>
#include <string.h>
#include "beastv2_header.h"

#include "abc_rand.h"

#include "globalvars.h"
#include "abc_mcmc.h"
#include "abc_ts_func.h"

typedef struct PREC_FUNCS {
	void (*SetNtermsPerPrecGrp)(I16PTR nTermsPerGrp, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState);
	void (*SetPrecXtXDiag)(F32PTR precXtXDiag, BEAST2_BASIS_PTR b, int NUMBASIS, PRECSTATE_PTR precState);
	void (*IncreasePrecValues)(        BEAST2_MODEL_PTR model);	
	void (*SetPropPrecXtXDiag_NtermsPerGrp)(BEAST2_MODEL_PTR model, BEAST2_BASIS_PTR basis, NEWTERM_PTR new );
	//F32  (*ComputeLog_det_prior)(BEAST2_MODEL_PTR model, I08PTR nTermsPerPrecGrp, I32 K);
	void (*ComputeMargLik)(            BEAST2_MODELDATA_PTR data, PRECSTATE_PTR precState,BEAST2_YINFO_PTR yInfo, BEAST2_HyperPar_PTR hyperPar);
	void (*ResamplePrecValues)(        BEAST2_MODEL_PTR model,    BEAST2_HyperPar* hyperPar,VOID_PTR stream);
	void (*chol_addCol)(               F32PTR Au,                 F32PTR U, F32PTR precPrior, I64 N, I64 K0, I64 K1);
} PREC_FUNCS;

void SetUpPrecFunctions(I08 precPriorType, I32 q, PREC_FUNCS* funcs);
