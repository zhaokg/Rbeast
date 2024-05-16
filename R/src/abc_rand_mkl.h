#pragma once
#define r_vslNewStream(stream, METHOD, seed)			  vslNewStream(stream, METHOD, seed)		
#define r_viRngUniformBits32(METHOD, stream, N, rnd)	  viRngUniformBits32(METHOD, stream, N, rnd)    
#define r_vsRngGamma(METHOD, stream, N, rnd, a, b, beta)  vsRngGamma(METHOD, stream, N, rnd, a, b, beta)  
#define r_vsRngGaussian(MTEHDO,stream,N, rnd, a, b)       vsRngGaussian(MTEHDO,stream,N, rnd, a, b)    
#define r_vslDeleteStream(stream)						  vslDeleteStream(stream) 
//-----------------------------------------------------------------------------------------