#pragma once

#include "abc_common.h"
#include "abc_blas_lapack_lib.h"
#include "abc_ide_util.h"
#include "abc_mem.h"
#include <inttypes.h>  //#include <stdint.h>
#include "beastv2_header.h"


void AllocateXXXMEM(F32PTR* Xt_mars, F32PTR* Xnewterm, F32PTR* Xt_zeroBackup, BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt, MemPointers* MEM);
void AllocateYinfoMEM(BEAST2_YINFO_PTR yInfo, BEAST2_OPTIONS_PTR opt, MemPointers* MEM);
