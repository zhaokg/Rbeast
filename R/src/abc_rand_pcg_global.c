#include <math.h>  
#include <string.h>     //memset
#include "abc_000_warning.h"

#include "abc_datatype.h" 
#include "abc_rand_pcg_global.h" 
#include "abc_ide_util.h"    //r_printf
#include "abc_vec.h"  

#include "assert.h"  
 

//https://hal.archives-ouvertes.fr/hal-02700791/document: a nice intro to the pcg alg
/*Practical seed-recovery for the PCG Pseudo-Random Numbr Generator*/

/*
Internally, this RNG uses two 64-bit integers for its internal state, consisting of:
the current state — the RNG iterates through all 264 possible internal states.
the RNG-sequence constant — a value that defines which of 263 possible random sequences
the current state is iterating through; it holds the same value over the lifetime of the RNG.
*/

/*
For this generator, there are 2^63 possible sequences of pseudorandom numbers. Each sequence is
entirely distinct and has a period of 2^64. The initseq argument selects which stream you will use. 
The initstate argument specifies where you are in that 264 period.
*/

#define PCG_DEFAULT_MULTIPLIER_64  6364136223846793005ULL //https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.c
#define PCG_DEFAULT_INCREMENT_64   1442695040888963407ULL ////https://github.com/imneme/pcg-c/blob/master/include/pcg_variants.h


//https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.h
#define PCG_DEFAULT_GLOBAL_STATE_64     0x853c49e6748fea9bULL
#define PCG_DEFAULT_GLOBAL_INCREMENT_64 0xda3e39cb94b95bdbULL

// Linear Congruential Generator (LCG):
// X     =  X    * MULTIPLIER + INCREMENT(shift/plus/state.inc) mod 2^64
// State =  State* MULTIPLIER + INCREMENT(shift/plus/INC)       mod 2^64


typedef struct pcg32_random_struct { 	//https://github.com/imneme/pcg-c-basic/blob/master/pcg_basic.h
	U64 state; /// RNG state.All values are possible.
	U64 inc;   // Controls which RNG sequence (stream) is
			   // selected. Must *always* be odd.
} pcg32_random_t;

static pcg32_random_t global_state = {
	.state= PCG_DEFAULT_GLOBAL_STATE_64,
	.inc  = PCG_DEFAULT_GLOBAL_INCREMENT_64 };

void pcg_set_seed(U64 initstate, U64 initseq)
{	
	extern void init_gauss_rnd(void);
	init_gauss_rnd();

	if (initstate == 0 ) {
	// Used to just initialzie the Guass structure.
		return;
	}

	initstate = initstate & initseq; //Added bcz only initseq is supplied as a seed. We run initseq to randomize initstate a little bit

	initstate = initstate == 0 ? PCG_DEFAULT_GLOBAL_STATE_64    : initstate;
	initseq   = initseq == 0   ? PCG_DEFAULT_GLOBAL_INCREMENT_64 : initseq;

	global_state.state = 0U;
	global_state.inc   = (initseq << 1u) | 1u;  	//inc must be an odd number

	U32 rnd;
	pcg_random(&rnd, 1);
	global_state.state += initstate;
	pcg_random(&rnd, 1);
}
void pcg_print_state(void) {r_printf("state: %30" PRIu64 " inc: %30"  PRIu64 "\n", global_state.state, global_state.inc);}

// pcg32_random_r(rng)
//   Generate a uniformly distributed 32-bit random number
void pcg_random(U32PTR rnd, I32 N)
{
	U64 oldstate = global_state.state;
	U64 shift    = global_state.inc;

	for (int i = 0; i < N; i++) {
		U32 xorshifted =     ((oldstate >> 18u) ^ oldstate) >> 27u   ;
		U32 rot = oldstate >> 59u;
		//*rnd++ = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
		*rnd++ = (xorshifted >> rot) | (xorshifted << ((- *((I32*)(&rot))  ) & 31));
		oldstate = oldstate * PCG_DEFAULT_MULTIPLIER_64 + shift;
	}
	global_state.state = oldstate;
}

/**********************************************/
// No need to export this func
U64 pcg_advance_lcg_64(U64 state, U64 delta, U64 cur_mult, U64 cur_plus)
{   
	//Advance the state by a total of delta steps
	//https://github.com/imneme/pcg-c/blob/master/src/pcg-advance-64.c
	U64 acc_mult = 1u;
	U64 acc_plus = 0u;
	while (delta > 0) {
		if (delta & 1) {
			acc_mult *= cur_mult;
			acc_plus = acc_plus * cur_mult + cur_plus;
		}
		cur_plus = (cur_mult + 1) * cur_plus;
		cur_mult *= cur_mult;
		delta /= 2;
	}
	return acc_mult * state + acc_plus;
}
void pcg_get_lcg_multiplier_shift_multistep( I32 delta, U64 cur_mult, U64 cur_plus, U64 * mult, U64 *shift) 
{
	// Find the multliplier and shift (plus) for a multi-step advance

	//I32 delta = 4; //Move forward by 4 steps
	//U64 cur_mult = PCG_DEFAULT_MULTIPLIER_64;   //Fiexed to a const
	//U64 cur_plus = global_state.INC;            //Supplied by the user
	
	U64 acc_mult = 1u;
	U64 acc_plus = 0u;
	while (delta > 0) {
		if (delta & 1) {
			acc_mult *= cur_mult;
			acc_plus = acc_plus * cur_mult + cur_plus;
		}
		cur_plus = (cur_mult + 1) * cur_plus;
		cur_mult *= cur_mult;
		delta /= 2;
	}	
	//r_printf("%" PRIu64 " %"  PRIu64 "\n",  acc_mult, acc_plus);
	mult[0]  = acc_mult;
	shift[0] = acc_plus;
	//return acc_mult * state + acc_plus;
}

 

/**************************************************************************/
//  Pre-caclcuated constants for the Gaussain radom generator based on 
//  a hybrid rejection sampling method
/**************************************************************************/
GAUSS_CONSTANT GAUSS
= { 
	
.x = { 0, 0.019584285230127, 0.039176085503098, 0.058782936068943, 0.078412412733112,
0.098072152488661, 0.117769874579095, 0.137513402144336, 0.157310684610171, 0.177169820991740, 0.197099084294312,
0.217106947210130, 0.237202109328788, 0.257393526100938, 0.277690439821577, 0.298102412930487, 0.318639363964375,
0.339311606538817, 0.360129891789569, 0.381105454763557, 0.402250065321725, 0.423576084201200, 0.445096524985516,
0.466825122852590, 0.488776411114670, 0.510965806738247, 0.533409706241281, 0.556125593618691, 0.579132162255556,
0.602449453164424, 0.626099012346421, 0.650104070647995, 0.674489750196082, 0.699283302383220, 0.724514383492366,
0.750215375467941, 0.776421761147928, 0.803172565597918, 0.830510878205399, 0.858484474141832, 0.887146559018876,
0.916556667533113, 0.946781756301046, 0.977897543940542, 1.009990169249582, 1.043158263318454, 1.077515567040280,
1.113194277160929, 1.150349380376008, 1.189164350199337, 1.229858759216589, 1.272698641190536, 1.318010897303537,
1.366203816372098, 1.417797137996268, 1.473467577947101, 1.534120544352547, 1.601008664886076, 1.675939722773444,
1.761670410363067, 1.862731867421652, 1.987427885929896, 2.153874694061456, 2.41755901623650 },

//y[i+1]/y[i]*2^24
.yRatio = { 16773999, 16767562, 16761112, 16754640, 16748136, 16741589, 16734989, 16728325, 16721587, 16714763,
16707840, 16700807, 16693651, 16686358, 16678913, 16671302, 16663507, 16655512, 16647297, 16638843, 16630128, 16621128,
16611818, 16602170, 16592154, 16581736, 16570879,16559544, 16547684, 16535250, 16522186, 16508431, 16493913, 16478554,
16462265, 16444942, 16426470, 16406715, 16385522, 16362712, 16338074,16311363, 16282287, 16250497, 16215576, 16177016,
16134195, 16086343, 16032497, 15971429, 15901559, 15820813, 15726414, 15614561, 15479906,15314686, 15107183, 14838856,
14478538, 13969524, 13196768, 11886086, 9182634 },

//1/y[i]*2^24
/*
.yInverse = { 0.167772160000000e8, 0.167804337107029e8, 0.167900954887396e8, 0.168062273321792e8, 0.168288727728840e8,
0.168580931819422e8, 0.168939682029590e8, 0.169365963190389e8, 0.169860955611753e8,0.170426043678181e8, 0.171062826076636e8,
0.171773127802660e8, 0.172559014119663e8, 0.173422806679528e8, 0.174367102051039e8, 0.175394792947299e8, 0.176509092495599e8,
0.177713561954893e8, 0.179012142359128e8,0.180409190651708e8, 0.181909520980666e8, 0.183518451949447e8, 0.185241860769783e8,
0.187086245447135e8, 0.189058796353774e8, 0.191167478819898e8, 0.193421128712637e8, 0.195829563393308e8, 0.198403710967247e8,
0.201155761397016e8, 0.204099343877364e8, 0.207249735920063e8, 0.210624110937242e8, 0.214241832835665e8, 0.218124808367550e8,
0.222297910899452e8, 0.226789493099905e8, 0.231632011146443e8, 0.236862789891804e8,0.242524967693829e8, 0.248668672301306e8,
0.255352496766927e8, 0.262645369022485e8, 0.270628943829048e8, 0.279400696443366e8, 0.289077971595910e8, 0.299803352204569e8,
0.311751880874552e8, 0.325140929293103e8,0.340243927554388e8, 0.357409846346236e8, 0.377091470037256e8, 0.399887489452612e8,
0.426607037407803e8, 0.458372068395134e8, 0.496786432318855e8, 0.544228824001546e8, 0.604390940468508e8, 0.683340892537681e8,
0.791831158870118e8, 0.950978898378280e8, 1.208991317887296e8, 1.706491807080511 },

 .PARAM_R     = 2.41755901623650,
 .INV_PARAM_R = 0.413640367529367
 */
};

#define INV_2p24  ((F64)1./(1LL<<24))
#define INV_2p25  ((F64)1./(1LL<<25))
#define INV_2p32  ((F64)1./(1LL<<32))


extern void init_gauss_rnd(void);
void init_gauss_rnd(void) {

	static I08 isInitialized = 0;
	if (isInitialized) {
		return;
	}

	for (int i = 0; i < 63; i++) {
		GAUSS.yRatio[i] = exp(-0.5 * (GAUSS.x[i + 1] * GAUSS.x[i + 1] - GAUSS.x[i] * GAUSS.x[i]));
	}

	for (int i = 0; i < 63; i++) {
		if (GAUSS.x[i + 1] >= 1.0) {
			GAUSS.inflectionId = i;
			break;
		}
	}
	
	// this is the optimal lamda maximizing the acceptance ratio.
	F32 a = GAUSS.x[63];
	GAUSS.exp_lamda = (a + sqrt(a * a + 4)) / 2;

	F32PTR x    = GAUSS.x;
	F32    del  = GAUSS.x[1];
	I32    Imax = ceil(GAUSS.x[63] / del); 

	// GUASS.indices saves the left-side bin index
	for (int i = 0; i < Imax; i++) {
		F32 low = i * del;
		F32 up  = (i + 1) * del;
		
		GAUSS.indices[i] = -9999;

		int bingo = 0;
		for (int k = 0; k < 64; k++) {

			if (x[k] >= low && x[k] <= up) {
				if (x[k] == low) {
					GAUSS.indices[i] = k;					
				}	else {
					GAUSS.indices[i] = k - 1;
				}	

				bingo = 1;
				break;		
			}


			if (x[k] <= low && up <= x[k + 1]) {
				GAUSS.indices[i] = k;
				bingo = 1;
				break;
			}

		}
	
		assert(bingo);
	}


	isInitialized = 1;
}

void pcg_gauss(F32PTR RND, int N)
{
	for (int i = 1; i <= N; i++) {

		U32 rnd[2];
		pcg_random(rnd, 2);

		U32  U24_32 = (rnd[0] >> 7) * INV_2p25;
		U32  BinIdx  = rnd[0] & 0x3f;
		I32  sign    = rnd[0] & 0x40;  //U24 & 0x80		

		F32 x; 
		if (BinIdx < 63) 	{
			F32 delta   = GAUSS.x[BinIdx + 1] - GAUSS.x[BinIdx];
			F32 yRatio  = GAUSS.yRatio[BinIdx] ;
			while (1) {
 
				if (U24_32 <= yRatio) {
					x = GAUSS.x[BinIdx] + delta * U24_32 / yRatio;
					break;
				}

				F64 U1    = rnd[1] * INV_2p32;
				int check = U24_32 <= yRatio + (1.f - yRatio) * U1;	
				x         = GAUSS.x[BinIdx+1] -  delta * U1;
			 
				if ((I32) BinIdx < GAUSS.inflectionId && check)      { break; }
				if ((I32)BinIdx > GAUSS.inflectionId && check == 0) { goto UpdateSamplerLabel; }
				if ( log(U24_32) <= -0.5f* (x * x- GAUSS.x[BinIdx]* GAUSS.x[BinIdx])  ) {
					break;
				}				

				// Need at least two rand numbers
				UpdateSamplerLabel:
				pcg_random(rnd, 2);
				U24_32 = (F64)rnd[0] * INV_2p32; //32-bit random number
			}
		} else	{
		while (1) {
				F64 U1 = rnd[1] * INV_2p32;                   // divideed by 2^32			
				x = GAUSS.x[63] - log(U1) / GAUSS.exp_lamda;        // x is inf if U1==0
				//if (expf(-GAUSS.PARAM_R * (x - 0.5f * GAUSS.PARAM_R)) * U2  < expf(-0.5f * x * x))
				if (log(U24_32) < -0.5f * ((x - GAUSS.exp_lamda) * (x - GAUSS.exp_lamda)))
					break;

				// Need at least two rand numbers
				pcg_random(rnd, 2);
				U24_32 = rnd[0] * INV_2p32; //divideed by 2^32
			}
		}

		*RND++ = sign? x:-x;
	}//Completion of the For loop
}
 
void pcg_gamma(F32PTR rnd, F32 a, int N)
{
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
	F32 INV_MAX = 2.328306436538696e-10f;
	if (a >= 1)
	{
		F32 b = a - 1.0f;
		F32 c = a + a + a - 0.75f;

		for (int i = 0; i < N; i++)
		{
			F32 gam;
			while (1)
			{

				
//#if defined(__GNUC__) || defined(__CLANG__) 
//	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				F32 u[2];
				pcg_random((U32PTR)u, 2);
				u[0] = (F32)(*(U32PTR)&u[0]) * INV_MAX;
				u[1] = (F32)(*(U32PTR)&u[1]) * INV_MAX;

				F32 w = u[0] * (1 - u[0]);
				F32 y = sqrtf(c / w) * (u[0] - 0.5f);
				gam = b + y;
				if (gam >= 0)
				{
					F32 z = 64.0f * (w*w*w) * (u[1] * u[1]);
					if (z <= (1 - 2 * (y*y) / gam))
						break;

					F32 logZ = logf(z);
					if (b == 0)					{
						if (logZ <= -2 * y) break;
					}		else					{
						if (logZ <= 2 * (b * logf(gam / b) - y)) break;
					}
				} //CODE_EOF of  (gam >= 0)

			}//CODE_EOF WHILE
			*rnd++ = gam;

		}//CODE_EOF FOR
		return;
	}

	if (a > 0) // 0<a<1
	{
		F32 b = 1 + .3678794f * a;
		for (int i = 0; i < N; i++)
		{
			F32 gam;
			while (1)
			{
//#if defined(__GNUC__) || defined(__CLANG__) 
//	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				F32 u[2];
				pcg_random((U32PTR)u, 2);
				u[0] = (F32)(*(U32PTR)&u[0]) * INV_MAX;
				u[1] = (F32)(*(U32PTR)&u[1]) * INV_MAX;
//#if defined(__GNUC__) || defined(__CLANG__) 
//	ENABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED)
//#endif
				F32 c = b * u[0];
				if (c < 1)
				{
					gam = expf(logf(c) / a);
					if (-logf(1 - u[1]) >= gam) break;
				}
				else
				{
					gam = -logf((b - c) / a);
					if (-logf(1 - u[1]) >= (1 - a) * logf(gam))
						break;
				}
			}//CODE_EOF WHILE
			*rnd++ = gam;
		}//CODE_EOF FOR
		return;
	}


	if (a < 0.f) 	{
		F32 nan = getNaN();
		for (int i = 0; i < N; i++)		*rnd++ = nan;
		return;
	}


	if (a == 0.f)	{
		for (int i = 0; i < N; i++)	*rnd++ = 0.f;
		return;
	}


}


void  pcg_wishart_unit_lowtriangle_nozeroout(F32PTR rnd, F32PTR tmp, I32 m, F32 df )
{
	I32 N = (m - 1)*m / 2;
     pcg_gauss(tmp, N);

	F32PTR tpt;
	tpt = rnd + 1;

	for (int i = 1; i<m  ; i++)
	{			
		I32  movElem = m - i;
		for (int j = 0; j < movElem; j++)
		{
			*tpt++ = *tmp++;
		}
		tpt = tpt + i+1;
	}
	
	tmp = rnd;
	for (int i = 1; i <=m; i++)
	{
		pcg_gamma(tmp, (F32)(df-i+1.)/ 2.f, 1);
		*tmp = sqrtf((*tmp) * 2.f);
		tmp = tmp + (m + 1);
	}

}

void pcg_wishart_unit_lowtriangle_nozeroout_notmp(F32PTR wishrnd, I32 m, F32 df)
{
	// This is an implmentation of the Barllet algortihm
	memset(wishrnd, 0, sizeof(F32) * m * m);

	// the number of elemnts below the diagonal: (0+ (m=1))*m/2
	I32     numGaussRnd = (m - 1) * m / 2;
	F32PTR  gauss       = wishrnd + m * m  - numGaussRnd;
	pcg_gauss(gauss, numGaussRnd);

	for (int col = 1; col < m; ++col) {
		for (int row = col + 1; row <= m; row++) {
			wishrnd[row - 1] = *gauss++;
		}
		wishrnd += m;
	}
	wishrnd -= (m - 1) * m;  //go back to the start of rnd

   //  Fill the diagonl with gamma random variables
	for (int i = 1; i <= m; i++) {
		F32 chisqaure;
		pcg_gamma(&chisqaure, (F32)(df - i + 1) / 2.f, 1);
		*wishrnd = chisqaure = sqrtf(chisqaure * 2.f);
		wishrnd += (m + 1);
	}

}

 

 

void  pcg_wishart_unit_lowtriangle_zeroout_notmp(F32PTR wishrnd, I32 m, F32 df)
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
	pcg_gauss(gauss, numGaussRnd);
	gauss = gauss + numGaussRnd-1; //go to the CODE_EOF of the guass

	wishrnd += (m-2) * m; // go to the start of the second row from the last
	for (int col =m-1; col >=2; col--) {
		for (int row = m; row >= col+1; row--) {
			wishrnd[row - 1] = *gauss--;
		}
		wishrnd -= m;
	}
	 wishrnd = wishrnd;  //go back to the start of rnd
	 


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
		pcg_gamma(&chisqaure, (F32)(df - i + 1) / 2.f, 1);
		*wishrnd = chisqaure = sqrtf(chisqaure * 2.f);
		wishrnd += (m + 1);
	}

}

#include "abc_mat.h"
void  pcg_invwishart_upper(F32PTR iwrnd_upper, F32PTR iwrnd_upper_inv, F32PTR tmp, I32 m, F32PTR Qup, F32 df)
{   
	//IW(Q,v): Q=Qu'*Qup
	//The IW matrix should be R'*R where R=iwrnd_upper= inv(W)*Qup
	// 
	// irwnd_upper_inv- inv(R)=inv(Qup)*W.
	 
	F32PTR Wlower = tmp;

	// This is an implmentation of the Barllet algortihm
	//local_pcg_wishart_unit_lowtriangle_nozeroout(rng, Wlower, Wlower +m*m, m, df);
	pcg_wishart_unit_lowtriangle_zeroout_notmp(Wlower, m, df);
	

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
