#pragma once
#include "abc_000_macro.h"
#if  PCGRAND_LIBRARY == 1
	#include "abc_rand_pcg_global.h"
	#include "abc_rand_pcg_local.h"
#elif MKLRAND_LIBRARY==1
	#include "abc_rand_mkl.h"
#endif