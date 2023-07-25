#include "math.h"  //expf logf sqrtf

#include "abc_000_warning.h"

#include "abc_rand_pcg_global.h"
#include "abc_rand_pcg_local.h"
#include "abc_mat.h"  

 void (*local_pcg_print_state)(local_pcg32_random_t* rng);
 void (*local_pcg_set_seed)(local_pcg32_random_t* rng, U64 initstate, U64 initseq);
 void (*local_pcg_random)(local_pcg32_random_t* rng, U32PTR rnd, I32 N);

#define INV_2p24  ((F64)1./(1LL<<24))
#define INV_2p25  ((F64)1./(1LL<<25))
#define INV_2p26  ((F64)1./(1LL<<26))
#define INV_2p32  ((F64)1./(1LL<<32))
 
 void local_pcg_randint(local_pcg32_random_t* rng, int bnd, I32PTR rndunif, int N) {

	 // Generate a random integer between 0 and bnd-1 : 0, 1,2,3,..., bnd-1L.
	 // or 1, 2, 3, bnd
	 // github.com/imneme/pcg-c-basic/blob/master/pcg_basic.c

	 if (bnd <= 256) {
		 U08 RAND[128];
		 I32 SIZE       = (min(N, 128) +3)/4*4;		 
		 I32 SIZEminus1 = SIZE - 1L;
		 I32 CURSOR     = SIZE;     // force drawing RAND numbers inside the loop

		 U08 threashold = 256L % bnd;
		 for (int i = 1; i <= N; ++i) {
			 while (1) {
				 if (CURSOR > SIZEminus1) {
					 local_pcg_random(rng, RAND, SIZE/4); CURSOR = 0;
				 }
				 U08 u = RAND[CURSOR++];
				 if (u >= threashold) {
					 *rndunif++ = (u % bnd) + 1;
					 break;
				 }
			 }
		 } // Completion of the For loop


	 } 
	 else if ( bnd  <= (1L<<16) ) {
		 U16 RAND[128];
		 I32 SIZE       = ( min(N, 128) +1)/2L*2L;		 
		 I32 SIZEminus1 = SIZE - 1L;
		 I32 CURSOR     = SIZE;     // force drawing RAND numbers inside the loop

		 U08 threashold = (1L<<16) % bnd;
		 for (int i = 1; i <= N; ++i) {
			 while (1) {
				 if (CURSOR > SIZEminus1) {
					 local_pcg_random(rng, RAND, SIZE/2); CURSOR = 0;
				 }
				 U16 u = RAND[CURSOR++];
				 if (u >= threashold) {
					 *rndunif++ = (u % bnd) + 1;
					 break;
				 }
			 }
		 } // Completion of the For loop

	 }
	 else  {
		 U32 RAND[128];
		 I32 SIZE       = min(N, 128) ;
		 I32 SIZEminus1 = SIZE - 1L;
		 I32 CURSOR     = SIZE;     // force drawing RAND numbers inside the loop

		 U08 threashold = (1ULL << 32) % bnd;
		 for (int i = 1; i <= N; ++i) {
			 while (1) {
				 if (CURSOR > SIZEminus1) {
					 local_pcg_random(rng, RAND, SIZE  ); CURSOR = 0;
				 }
				 U32 u = RAND[CURSOR++];
				 if (u >= threashold) {
					 *rndunif++ = (u % bnd) + 1;
					 break;
				 }
			 }
		 } // Completion of the For loop

	 }

 
}
 
 static F32 local_pcg_gauss_binwise(local_pcg32_random_t* rng, int BinIdx) {
	 
	 F32 x;
	 if (BinIdx < 63) {
		 F32 delta  = GAUSS.x[BinIdx + 1] - GAUSS.x[BinIdx];
		 F32 yRatio = GAUSS.yRatio[BinIdx];
		 while (1) {
			 U32 RAND;   //  Need up to 2 random numbers
			 local_pcg_random(rng, &RAND, 1);
			 F32 U0 = RAND * INV_2p32;
			 if (U0 <= yRatio) {
				 x = GAUSS.x[BinIdx] + delta * U0 / yRatio;
				 break;
			 }

			 local_pcg_random(rng, &RAND, 1);
			 F64 U1 = RAND * INV_2p32;
			 int check = U0 <= yRatio + (1.f - yRatio) * U1;
			 x = GAUSS.x[BinIdx + 1] - delta * U1;
			 if (BinIdx  < GAUSS.inflectionId && check)                          { break; }
			 if (log(U0) <= -0.5f * (x * x - GAUSS.x[BinIdx] * GAUSS.x[BinIdx])) { break; }
		 }
		 return x;
	 }  else {
		 // Use the exponential generator for accept-reject sampling
		 while (1) {
			 U32 RAND[2];
			 local_pcg_random(rng, RAND, 2);
			 F32 U0 = RAND[0] * INV_2p32;
			 F32 U1 = RAND[1] * INV_2p32;
			 x = GAUSS.x[63] - log(U1) / GAUSS.exp_lamda;        // x is inf if U1==0
			 //if (expf(-GAUSS.PARAM_R * (x - 0.5f * GAUSS.PARAM_R)) * U2  < expf(-0.5f * x * x))
			 if (log(U0) < -0.5f * ((x - GAUSS.exp_lamda) * (x - GAUSS.exp_lamda)))
				 break;
		 }
		 return x;
	 }
	 

 }

 F32 local_pcg_trgauss_oneside_scalar(local_pcg32_random_t* rng, F32 a, int whichside) {
	 // whichsize=1:  rightside
	 // whichside=-1; leftside
	 a = a * whichside;

	 F32 x;
	 if (a >= GAUSS.amax) {		 
		 F32 lamda = (a + sqrtf(a * a + 4)) / 2; //www.rob-mcculloch.org/2019_cs/webpage/Geweke.pdf
		 while (1) {			 
			 U32 RAND[2];  
			 local_pcg_random(rng, RAND, 2);
			 F32 U0 = RAND[0] * INV_2p32;
			 F32 U1 = RAND[1] * INV_2p32;
			 x = a - log(U1) /lamda;        // x is inf if U1==0
			 //if (expf(-GAUSS.PARAM_R * (x - 0.5f * GAUSS.PARAM_R)) * U2  < expf(-0.5f * x * x))
			 if (log(U0) < -0.5f * ((x - lamda) * (x - lamda)))
				 break;
		 }
		 return x*whichside;
	 }

	 if (a >= 0) {
		 F32 dx     = GAUSS.x[1];
		 int kidx   = (int)(a / dx);
		 int BinIdx = GAUSS.indices[kidx];
		 //-------i-#-i+1--------------------
		 // ----|-i----|-----|-----|
		 // A dx bin may be linked to two BINs. kidx refers to the leftside one
		 // Check if a falls within the left Bin or the right bin
		 int BinStart = a >= GAUSS.x[BinIdx] && a <= GAUSS.x[BinIdx + 1] ? BinIdx : BinIdx + 1;
		 int NumBins  = (63 - BinStart) + 1;

		 U32 threashold = (1ULL << 32) % NumBins;
		 for (;;) {
			 int BinSel;
			 while (1) {
				 U32 RAND;  
				 local_pcg_random(rng, &RAND, 1);
				 if (RAND >= threashold) {
					 BinSel = (RAND % NumBins);
					 break;
				 }
			 }

			 BinSel = BinSel +BinStart;
			 x = local_pcg_gauss_binwise(rng, BinSel);
			 if (x >= a) {
				 break;
			 }
		 }

		 return  x * whichside;
	 
	 }
 
    // a <0

	
	 int BinEnd = 63 + 64;  // the first 64 are the positive bins, followed by the negative bins
	 if (a > -GAUSS.x[63]) {
		 // if not the leftmost bin
		 F32 apos    = -a;
		 F32 dx     = GAUSS.x[1];
		 int kidx   = (int)(apos / dx);
		 int BinIdx = GAUSS.indices[kidx];
		 BinIdx = apos >= GAUSS.x[BinIdx] && apos <= GAUSS.x[BinIdx + 1] ? BinIdx : BinIdx + 1;
		 BinEnd = 63+ 1+BinIdx;
	 }

	  // ----|-i----|-----|--63-|--0|1|2|
	 int NumBins    = BinEnd + 1;
	 U32 threashold = (1ULL << 32) % NumBins;

	 for (;;) {
		 int BinSel;
		 while (1) {
			 U32 RAND;  local_pcg_random(rng, &RAND, 1);
			 if (RAND >= threashold) {
				 BinSel = RAND % NumBins;
				 break;
			 }
		 }

		 int BinSelPos=BinSel;
		 if (BinSel >= 64) {
			 BinSelPos = BinSel - 64;
		 }
		 x = local_pcg_gauss_binwise(rng, BinSelPos);
		 if (BinSel >= 64) {
			 x = -x;
		 }
		 if (x >= a) {
			 break;
		 }
	 }

	 return  x * whichside;

}

 F32 local_pcg_inverse_gauss(local_pcg32_random_t* rng, F32 u, F32 lamda) {
 // en.wikipedia.org/wiki/Inverse_Gaussian_distribution
	 U32 RAND;
	 local_pcg_random(rng, &RAND, 1L);  
	 F64 U1     = (RAND >> 6) * INV_2p26; // divuded by 2^26;
	 I32 BinIdx =  RAND & 0x3f;			  // the lower 6 bits

	 F32 v  = local_pcg_gauss_binwise(rng, BinIdx);
	 F32 y  = v * v;
	 F32 uy = u * y;
	 F32 x  = u + u / (lamda + lamda) * (uy - sqrtf(4 * lamda * uy + uy * uy));

	 if (U1 < u / (u + x)) {
		 return x;
	 } else {
		 return u* u / x;
	 }

  

 }
void local_pcg_gauss(local_pcg32_random_t* rng, F32PTR RND, int N) 
{
	U32 RAND[128];
	I32 SIZEminus2 = min(2 * N, 128) - 2;
	I32 CURSOR     = SIZEminus2 + 1;     // force drawing RAND numbers inside the loop

	for (int i = 1; i <= N; ++i) {

		//  Need at least two rand numbers
		if (CURSOR > SIZEminus2) {
			local_pcg_random(rng, RAND, SIZEminus2 + 2); CURSOR = 0;
		}

		F64 U24_32  = (RAND[CURSOR] >> 7) * INV_2p25; // divuded by 2^25;
		I64 BinIdx  = RAND[CURSOR]  & 0x3f;			   // the lower 6 bits
		I32 sign    = (RAND[CURSOR] & 0x40)  ;  // the highest bit         
		CURSOR++;

		F32 x;
		if (BinIdx < 63) {

			F32 delta  = (GAUSS.x[BinIdx + 1] - GAUSS.x[BinIdx]);	
			F32 yRatio = GAUSS.yRatio[BinIdx];
			while (1) {
			
				if (U24_32 <= yRatio) {
					x = GAUSS.x[BinIdx] + delta * U24_32 / yRatio;
					break;
				}

				F64 U1    = RAND[CURSOR++] * INV_2p32;
				int check = U24_32 <= yRatio + (1.f - yRatio) * U1;	
				x         = GAUSS.x[BinIdx+1] -  delta * U1;
			 
				if (BinIdx < GAUSS.inflectionId && check)      { break; }
				if (BinIdx > GAUSS.inflectionId && check == 0) { goto UpdateSamplerLabel; }

				if ( log(U24_32) <= -0.5f* (x * x- GAUSS.x[BinIdx]* GAUSS.x[BinIdx])  ) {
					break;
				}
				
				// Need at least two rand numbers
				UpdateSamplerLabel:
				if (CURSOR > SIZEminus2) {
					local_pcg_random(rng, RAND, SIZEminus2 + 2); CURSOR = 0;
				}					 
				U24_32 = (F64) RAND[CURSOR++] * INV_2p32; //32-bit random number
				
			}

		}
		else {
			// Use the exponential generator for accept-reject sampling
				
			while (1) {
				F64 U1 =  RAND[CURSOR++] * INV_2p32;                   // divideed by 2^32			
				x = GAUSS.x[63] - log(U1) / GAUSS.exp_lamda;        // x is inf if U1==0
				//if (expf(-GAUSS.PARAM_R * (x - 0.5f * GAUSS.PARAM_R)) * U2  < expf(-0.5f * x * x))
				if ( log(U24_32) < -0.5f * ( (x- GAUSS.exp_lamda)* (x - GAUSS.exp_lamda))  )
					break;

				// Need at least two rand numbers
				if (CURSOR > SIZEminus2) {
					local_pcg_random(rng, RAND, SIZEminus2 + 2); CURSOR = 0;
				}
				U24_32 = RAND[CURSOR++] * INV_2p32; //divideed by 2^32
			}
		}

		*RND++ =  sign ? x:-x;
	} // Completion of the For loop


}
 
void local_pcg_gamma(local_pcg32_random_t* rng, F32PTR rnd, F32 a, int N) {

	/*
	% This function is used to draw gamma random variables.
	% I can't remember who I nicked this off????
	% if you let me know then I'll credit them, please email me.
	%
	% It was addapted from........
	%
	% RANDGAMM Generates N gamma random deviates.
	%     RANDGAMM(A, N) is a random deviate from the standard gamma
	%     distribution with shape parameter A.
	%
	%     B*RANDGAMM(A, N) is a random deviate from the gamma distribution
	%     with shape parameter A and scale parameter B.The distribution
	%     then has mean A*B and variance A*B^2.
	%
	%     See RAND.
	% GKS 31 July 93
	% Algorithm for A >= 1 is Best's rejection algorithm XG
	% Adapted from L.Devroye, "Non-uniform random variate
	% generation", Springer-Verlag, New York, 1986, p. 410.
	% Algorithm for A < 1 is rejection algorithm GS from
	% Ahrens, J.H.and Dieter, U.Computer methods for sampling
	% from gamma, beta, Poisson and binomial distributions.
	% Computing, 12 (1974), 223 - 246.  Adapted from Netlib
	% Fortran routine.

	gamma_store = zeros(N, M);
	*/

	const F32 INV_MAX = 2.328306436538696e-10f;

	U32 RAND[200];
	I32 SIZEminus2 = ((2 * N < 200) ? (2 * N + 2) : 200) - 2;
	I32 CURSOR;
	local_pcg_random(rng, RAND, SIZEminus2 + 2); CURSOR = 0;

	if (a > 1.f) {

		F32 b = a - 1.0f;             //b>0
		F32 c = a+a+a - 0.75f;

		for (int i = 0; i < N; i++)	{
			F32 gam;
			while (1) {

//#if defined(__GNUC__) || defined(__CLANG__) 
//	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				if (CURSOR > SIZEminus2) { //Need at least two rand numbers
					local_pcg_random(rng, RAND, SIZEminus2 + 2);	CURSOR = 0;
				}
			 
				F32  u0 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;
				F32  u1 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;

				F32 w = u0 * (1.f - u0);
				F32 y = sqrtf( c/w ) * (u0 - 0.5f);

				gam = b + y;
				if (gam >= 0)  {
					F32 z = 64.0f * (w*w*w) * (u1 * u1);
					if (z <= (1 - 2 * (y*y) / gam)) break;
					F32 logZ = logf(z);
					if (logZ <= 2 * (b * logf(gam / b) - y)) break;
				} //CODE_EOF of  (gam >= 0)

			}//CODE_EOF WHILE

			*rnd++ = gam;

		}//CODE_EOF FOR
		return;
	}

	if (a == 1.f)	{

		F32 b = a - 1.0f;         // b=0
		F32 c = a + a + a - 0.75f;

		for (int i = 0; i < N; i++) {
			F32 gam;
			while (1){

//#if defined(__GNUC__) || defined(__CLANG__) 
//	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				if (CURSOR > SIZEminus2) { //Need at least two rand numbers
					local_pcg_random(rng, RAND, SIZEminus2 + 2);	CURSOR = 0;
				}
			 
				F32  u0 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;
				F32  u1 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;

				F32 w = u0 * (1.f - u0);
				F32 y = sqrtf(c / w) * (u0 - 0.5f);

				gam = b + y;
				if (gam >= 0)	{
					F32 z = 64.0f * (w*w*w) * (u1 * u1);
					if (z <= (1 - 2 * (y*y) / gam)) break;
					F32 logZ = logf(z);					
   				    if (logZ <= -2 * y) break;
				} //CODE_EOF of  (gam >= 0)

			}//CODE_EOF WHILE

			*rnd++ = gam;

		}//CODE_EOF FOR
		return;
	}

	if (a > 0) // 0<a<1
	{
		F32 b = 1 + .3678794f * a;
		for (int i = 0; i < N; i++)	{

			F32 gam;
			while (1) 	{

//#if defined(__GNUC__) || defined(__CLANG__) 
//	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				if (CURSOR > SIZEminus2) { //Need at least two rand numbers
					local_pcg_random(rng, RAND, SIZEminus2 + 2);	CURSOR = 0;
				}
				 
				F32 u0 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;
				F32 u1 = (F32)(*(U32PTR)&RAND[CURSOR++]) * INV_MAX;
//#if defined(__GNUC__) || defined(__CLANG__) 
//	ENABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				F32 c = b * u0;
				if (c < 1) {
					gam = expf(logf(c) / a);
					if (-logf(1.f - u1) >= gam) break;
				} else	{
					gam = -logf((b - c) / a);
					if (-logf(1.f - u1) >= (1 - a) * logf(gam))
						break;
				}
			}//CODE_EOF WHILE
			*rnd++ = gam;
		}//CODE_EOF FOR
		return;
	}

	if (a < 0.f){
		F32 nan = getNaN();
		for (int i = 0; i < N; i++)	*rnd++ = nan;		
		return;
	}
	if (a == 0.f){
		for (int i = 0; i < N; i++)		*rnd++ = 0.f;		
		return;
	}
}

void local_pcg_wishart_unit_lowtriangle_zeroout(local_pcg32_random_t* rng, F32PTR wishrnd, F32PTR tmp, I32 m, F32 df)
{   
	// This is an implmentation of the Barllet algortihm
	memset(wishrnd, 0, sizeof(F32) * m * m);

	// the number of elemnts below the diagonal: (0+ (m=1))*m/2
	I32 numGaussRnd = (m-1)*m / 2;
	local_pcg_gauss(rng,tmp, numGaussRnd);

	for (int col = 1; col< m; ++ col) {				
		for (int row = col+1; row <= m; row++) 	{
			wishrnd[row-1] = *tmp++;
		}		
		wishrnd += m;
	}
	wishrnd -= (m-1) * m;  //go back to the start of rnd

   //  Fill the diagonl with gamma random variables
	for (int i = 1; i <=m; i++) 	{
		F32 chisqaure;
		local_pcg_gamma(rng, &chisqaure, (F32)(df-i+1)/ 2.f, 1);
		*wishrnd=chisqaure = sqrtf(chisqaure * 2.f);				
		wishrnd += (m + 1);
	}

}

 
#include "abc_vec.h"
void local_pcg_wishart_unit_lowtriangle_zeroout_notmp(local_pcg32_random_t* rng, F32PTR wishrnd, I32 m, F32 df)
{
	// This is an implmentation of the Barllet algortihm
 
	/*
	// the number of elemnts below the diagonal: (0+ (m=1))*m/2	 
	I32     numGaussRnd = (m - 1) * m / 2;
	F32PTR  gauss =  wishrnd + m * m - numGaussRnd;
	local_pcg_gauss(rng, gauss , numGaussRnd);
	
	for (int col = 1; col < m; ++col) {
		for (int row = col + 1; row <= m; row++) {
			wishrnd[row - 1] = *gauss++;
		}
		wishrnd += m;
	}
	wishrnd -= (m - 1) * m;  //go back to the start of rnd
	 */

 
	// the number of elemnts below the diagonal: (0+ (m=1))*m/2	
	I32     numGaussRnd = (m - 1) * m / 2;
	F32PTR  gauss       =  wishrnd + 1;
	local_pcg_gauss(rng, gauss, numGaussRnd);
	gauss = gauss + numGaussRnd-1; //go to the CODE_EOF of the guass

	wishrnd += (m-2) * m; // go to the start of the second row from the last
	for (int col =m-1; col >=2; col--) {
		for (int row = m; row >= col+1; row--) {
			wishrnd[row - 1] = *gauss--;
		}
		wishrnd -= m;
	}
	 wishrnd= wishrnd;  //go back to the start of rnd
	 


	 // Zeroout the upper triangle
	 for (int col = 0; col <m; col++) {
		 for (int row = 0; row<col; row++) {
			 wishrnd[row] = 0;
		 }
		 wishrnd += m;
	 }
	 wishrnd = wishrnd -m*m;  //go back to the start of rnd	

   //  Fill the diagonl with gamma random variables
	for (int i = 1; i <= m; i++) {
		F32 chisqaure;
		local_pcg_gamma(rng, &chisqaure, (F32)(df - i + 1) / 2.f, 1);
		*wishrnd = chisqaure = sqrtf(chisqaure * 2.f);
		wishrnd += (m + 1);
	}

}

void local_pcg_invwishart_upper(local_pcg32_random_t* rng, F32PTR iwrnd_upper, F32PTR iwrnd_upper_inv, F32PTR tmp, I32 m, F32PTR Qup, F32 df)
{   
	//IW(Q,v): Q=Qu'*Qup
	//The IW matrix should be R'*R where R=iwrnd_upper= inv(W)*Qup
	// 
	// irwnd_upper_inv- inv(R)=inv(Qup)*W.
	 
	F32PTR Wlower = tmp;

	// This is an implmentation of the Barllet algortihm
	//local_pcg_wishart_unit_lowtriangle_nozeroout(rng, Wlower, Wlower +m*m, m, df);
	local_pcg_wishart_unit_lowtriangle_zeroout_notmp(rng, Wlower, m, df);
	

	// Solve Wl*iwrnd_upper=Qup
	f32_copy(Qup, iwrnd_upper, m * m);
	for (int i = 0; i < m; ++i) {
		solve_L_as_L(Wlower, iwrnd_upper, m, m);
		iwrnd_upper = iwrnd_upper + m;
	}

	
	//Solve Qup*irwnd_upper_inv=W
	if (iwrnd_upper_inv) {
		// irwnd_upper_inv-inv(R)=inv(Qup)*W.		
		f32_copy(Wlower, iwrnd_upper_inv, m * m);
		for (int i = 0; i < m; ++i) {
			solve_U_as_U(Qup, iwrnd_upper_inv, m, m);
			iwrnd_upper_inv = iwrnd_upper_inv + m;
		}
	}
	
}


#include "abc_000_warning.h"