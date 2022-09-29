#pragma once

#include "beastv2_header.h"


extern void AllocInitModelMEM(BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt, MemPointers* MEM);
extern void  ReInit_PrecValues(BEAST2_MODEL_PTR model, BEAST2_OPTIONS_PTR opt);