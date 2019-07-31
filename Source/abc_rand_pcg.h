#ifndef  BEAST_RAND_HEADER
#define  BEAST_RAND_HEADER
#include "abc_000_macro.h"
#include "abc_datatype.h"
#define r_vslNewStream(stream_dummy,METHOD_dummy,seed)					pcg_set_seed(seed,56567456ULL) 
#define r_viRngUniformBits32(METHOD_dummy,stream_dummy,N,rnd)			pcg_random(rnd,N)
#define r_vsRngGamma(METHOD_dummy,stream_dummy,N,rnd,a,b,beta_dummy)  pcg_gamma(rnd,a,N)
#define r_vsRngGaussian(MTEHDO_dummy,stream_dummy,N,rnd,a,b)				pcg_gauss(rnd,N)
#define r_vslDeleteStream(stream_dummy)  
void pcg_set_seed(uint64_t initstate,uint64_t initseq);
void pcg_random(uint32_t * _restrict rnd,int32_t N);
void pcg_gauss(F32PTR RND,int N);
void pcg_gamma(F32PTR rnd,float a,int N);
void pcg_wishart_unit(F32PTR rnd,F32PTR tmp,int32_t m,int32_t df,char upperOrLower,char zerofill);
#endif
