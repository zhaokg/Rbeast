#include <string.h> //memset
#include <math.h>   //log
#include "abc_000_warning.h"

#include "beastv2_header.h"
#include "abc_mcmc.h"
#include "abc_math.h"  // for fastlog only
#include "globalvars.h"
 

static F64 K_OUT_OF_N(I32 N, I32 k)     {
	if ( k < 0 || k > N) return 0;
	k = min(k, N-k); 
	F64 y = 1, Ndbl = N, Kdbl = k;
	for (I32 i = 0; i <k; i++)	y *= Ndbl-- / Kdbl--;
	return y;
}

static F64 GetNumModelsGivenK(I32 N, I32 minSep, I32 K) {
	return K_OUT_OF_N(N - minSep * K + K- 1, K - 1);
}

static F64 GetGroupingNumv1(I32 n, I32 ng,  I32 k) {
	// math.stackexchange.com/questions/900828/number-of-groups-containing-at-least-1-and-at-most-k-elements/902458#902458
    // this is the recursion algorithim
	if (k * ng < n || n < ng) 	  return 0;	
	if (n==ng)			          return 1;	
	if (n < k+ng)			      return K_OUT_OF_N(n - 1, ng - 1);	
	F64 r= 0;
	for (I32 i = 0; i <= (n/k); ++i) {	r += K_OUT_OF_N(ng, i) * GetGroupingNumv1(n-i*k, ng-i, k-1);}
	return r;
}

static  F64 GetGroupingNumv2(I32 n, I32 ng, I32 k) {
	// math.stackexchange.com/questions/900828/number-of-groups-containing-at-least-1-and-at-most-k-elements/902458#902458
	// this is the fast algorithm

	if (n == ng*k) return 1;
	if (k == 2   ) return K_OUT_OF_N(ng, n - ng);

	F64 res      = 0, flipsign = -1.0;
	int ilimit   = min( (n - ng) / k, ng); 
	for (int i = 0; i <= ng; i++) {
		flipsign = -flipsign;
		res      += flipsign * K_OUT_OF_N(ng, i) * K_OUT_OF_N(n - i * k - 1, ng - 1);
	}	
	return res < 0.5 ? 0 : res; 
}
  
void PreCaclModelNumber(I32 minOrder, I32 maxOrder, I32 maxNumseg, I32 N, I32 minSep, F64PTR numMat, F64PTR totalNum, F64PTR NmodelsPerK, int priorType) {
 	
	if (priorType == 0 || NmodelsPerK == NULL) {
		return;
	}

	// Compute number of modles given the number of segments
	for (I32 ng = 1; ng <= maxNumseg; ng++) {
		NmodelsPerK[ng - 1] = GetNumModelsGivenK(N, minSep, ng);
		if (NmodelsPerK[ng - 1] > 0) NmodelsPerK[ng - 1] = log(NmodelsPerK[ng - 1]);
	}

	I32 KMAX = maxNumseg * maxOrder;
	// Approximate  number of modles given the number of segments and the number of total terms
	memset(numMat, 0, sizeof(F64) * KMAX * maxNumseg);
	for (int ng = 1; ng <= maxNumseg; ng++) {		
		for (int n = minOrder * ng; n <= maxOrder * ng; ++n) {
			numMat[(ng - 1) * KMAX + n - 1] = K_OUT_OF_N( n-ng* minOrder+ng-1, ng-1);
			if (numMat[(ng - 1) * KMAX + n - 1]  > 0) numMat[(ng - 1) * KMAX + n - 1] = log(numMat[(ng - 1) * KMAX + n - 1]);
		}	
	}
	return;


	if (priorType == 2) {
		// Compute number of modles given the number of segments fir all the possible cobinations of orders
		for (I32 ng = 1; ng <= maxNumseg; ng++) {
			NmodelsPerK[ng - 1] = GetNumModelsGivenK(N, minSep, ng) * pow( maxOrder-minOrder + 1, ng);
		}
		return;
	}
	   
	
	// Compute number of modles given the number of segments and the number of total terms
	memset(numMat, 0, sizeof(F64) * KMAX * maxNumseg);
	for (int ng = 1; ng <= maxNumseg; ng++) {
		NmodelsPerK[ng - 1] = GetNumModelsGivenK(N, minSep, ng);	
		for (int n = minOrder * ng; n <= maxOrder * ng; ++n) {
			numMat[(ng-1) * KMAX + n-1] = GetGroupingNumv2(n - (minOrder-1) * ng, ng, maxOrder-(minOrder-1) );
		}
		
	}

		
	if (priorType == 3) {
		memset(totalNum, 0, sizeof(F64) * KMAX);
		for (I32 ng = 1; ng <= maxNumseg; ++ng) {
			for (I32 n = minOrder * ng; n <= maxOrder * ng; ++n) {
				totalNum[n - 1] += numMat[(ng - 1L) * KMAX + n - 1] * NmodelsPerK[ng-1];
			}
		}
		return;
	}

	if (priorType == 4) {
		for (I32 ng = 1; ng <= maxNumseg; ++ng) {
			for (I32 n = minOrder * ng; n <= maxOrder * ng; ++n) {
				numMat[(ng - 1) * KMAX + n - 1] *= NmodelsPerK[ng - 1];
				//r_printf("%d: %f  %f\n", ng, numMat[(ng - 1) * KMAX + n - 1], Nmodels[ng - 1]);
			}
		}
		return;
	}
	
    	
	//for (I32 ng = 1; ng <= maxNumseg; ++ng) 	    r_printf("%d: %f \n", ng, Nmodels[ng - 1]);
	//for (I32 n = minOrder * 1L; n <= KMAX; ++n)   r_printf("totalNum %d: %f \n", n, totalNum[n - 1]);
	//for (I32 n = minOrder * 1L; n <= KMAX; ++n)   totalNum[n - 1] = -log(totalNum[n - 1]);
}

static F32 ST_ModelPriorFactor0(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {
	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  O1       = basis->prior.minOrder + (basis->type == TRENDID);
 	
	//I32  N        = newcol->N;
	//I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
	//I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	F64 factor0  ;         //  *(F32)(NUMORDER * nSegOld + 1) / (NUMORDER * nSegNew + 1);	                

	factor0 = 0 + basis->prior.modelComplexity * ( (nSegNew - nSegOld)+ (Knew- nSegNew*O1)- (Kold-nSegOld * O1));
		
	return (F32)factor0;

}

static F32 ST_ModelPriorFactor1(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {
		
	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  O1       = basis->prior.minOrder + (basis->type == TRENDID);
 	
	//I32  N        = newcol->N;
	//I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
	I32  KMAX = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	//F64 factor0  = VEC[nSegOld - 1] / VEC[nSegNew - 1] * MAT[(nSegOld - 1) * KMAX + Kold - 1] / MAT[(nSegNew - 1) * KMAX + Knew - 1];
	                //  *(F32)(NUMORDER * nSegOld + 1) / (NUMORDER * nSegNew + 1);	                

	F64 factor0  = VEC[nSegOld - 1] - VEC[nSegNew - 1] + MAT[(nSegOld - 1) * KMAX + Kold - 1] - MAT[(nSegNew - 1) * KMAX + Knew - 1];
	factor0 = factor0 + basis->prior.modelComplexity * ( (nSegNew - nSegOld)+ (Knew- nSegNew*O1)- (Kold-nSegOld * O1));
	
	F32 factor;
	if (new->jumpType == ChORDER)  
		factor=factor0; 
	else if (new->jumpType == BIRTH) {
		F32 propprobRatio =  (F32)basis->goodNum / ( basis->nKnot+1L);
		factor            = factor0 + log(propprobRatio);
	} else {
		F32 propprobRatio = (F32)(basis->goodNum + basis->prior.minSepDist*2)/basis->nKnot;
		factor            = factor0  - log(propprobRatio);
	}
	return (F32) factor ;

	//r_printf("%5.2f  %d %d %d \n", OLD/ MAT[(nSegOld - 1) * KMAX + Kold - 1], NUMORDER, nSegOld, Kold);
}

static F32 ST_ModelPriorFactor2(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {
		
	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  O1       = basis->prior.minOrder + (basis->type == TRENDID);

	//I32  N        = newcol->N;
	//I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
	I32  KMAX = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	//F64 factor0  = VEC[nSegOld - 1] / VEC[nSegNew - 1] * MAT[(nSegOld - 1) * KMAX + Kold - 1] / MAT[(nSegNew - 1) * KMAX + Knew - 1];
	                //  *(F32)(NUMORDER * nSegOld + 1) / (NUMORDER * nSegNew + 1);	               
	F64 factor0  = VEC[nSegOld - 1] - VEC[nSegNew - 1] + MAT[(nSegOld - 1) * KMAX + Kold - 1] - MAT[(nSegNew - 1) * KMAX + Knew - 1];

	factor0 = factor0 + basis->prior.modelComplexity * log((F32)nSegNew / nSegOld *(F32)(Knew - nSegNew * O1+1)  / (Kold - nSegOld * O1 + 1));

	F32 factor;
	if (new->jumpType == ChORDER)  
		factor=factor0; 
	else if (new->jumpType == BIRTH) {
		F32 propprobRatio =  (F32)basis->goodNum / ( basis->nKnot+1L);
		factor            = factor0 + log(propprobRatio);
	} else {
		F32 propprobRatio = (F32)(basis->goodNum + basis->prior.minSepDist*2)/basis->nKnot;
		factor            = factor0  - log(propprobRatio);
	}
	return (F32) factor ;

	//r_printf("%5.2f  %d %d %d \n", OLD/ MAT[(nSegOld - 1) * KMAX + Kold - 1], NUMORDER, nSegOld, Kold);
}

static F32 ST_ModelPriorFactor3(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {

	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	
	//I32  N        = newcol->N;
	//I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorVec;

	F64 factor0 = VEC[Kold - 1] / VEC[Knew - 1]; //* (F32)(NUMORDER * nSegOld + 1) / (NUMORDER * nSegNew + 1);

	factor0  =  log(factor0) + basis->prior.modelComplexity * log((F32)nSegNew / nSegOld);

	F32 factor;
	if (new->jumpType == ChORDER)  
		factor=factor0; 
	else if (new->jumpType == BIRTH) {
		F32 propprobRatio =  (F32)basis->goodNum / ( basis->nKnot+1L);
		factor            = factor0 + log(propprobRatio);
	} else {
		F32 propprobRatio = (F32)(basis->goodNum + basis->prior.minSepDist*2)/basis->nKnot;
		factor            = factor0  - log(propprobRatio);
	}
	return (F32) factor ;

	//r_printf("%5.2f  %d %d %d \n", OLD/ MAT[(nSegOld - 1) * KMAX + Kold - 1], NUMORDER, nSegOld, Kold);
}

static F32 ST_ModelPriorFactor4(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {

	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;

	//I32  N        = newcol->N;	
	I32  KMAX = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));

	F64PTR MAT  = basis->priorMat;
	F64 factor0 = MAT[(nSegOld - 1) * KMAX + Kold - 1] / MAT[(nSegNew - 1) * KMAX + Knew - 1] ; // *(F32)(NUMORDER * nSegOld + 1) / (NUMORDER * nSegNew + 1);

	factor0 = basis->prior.modelComplexity * log(factor0);

	//r_printf("%5.2f  %d %d %d \n", OLD/ MAT[(nSegOld - 1) * KMAX + Kold - 1], NUMORDER, nSegOld, Kold);
	F32 factor;
	if      (new->jumpType == ChORDER)
		factor = factor0;
	else if (new->jumpType == BIRTH) {
		F32 propprobRatio = (F32)basis->goodNum / ( basis->nKnot+1L);
		factor            = factor0 + log(propprobRatio);
		//if (basis->nKnot < 1) r_printf("\nbirth %d %f   %f  %d  %d %f %f\n", basis->nKnot, factor0, factor, Kold, Knew, MAT[(nSegOld - 1) * KMAX + Kold - 1], MAT[(nSegNew - 1) * KMAX + Knew - 1]);
	} else {
		F32 propprobRatio = (F32)(basis->goodNum + basis->prior.minSepDist*2)/basis->nKnot;
		factor            = factor0  - log(propprobRatio);
	}
	return (F32) factor;
}


static F32 ST_ModelPriorFactor3_old(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {   

 // The factor GLOBAL_MODEL_PRIOR_FACTOR is added if assuming the k (number of segments) follows  a geometric distrubtion lamda^k
 // GLOBAL_MODEL_PRIOR_FACTOR * (2 * basis->nKnot * 0 + 1) is used if assuming k follows a lamda^(k^2)
	I32  Kold     = basis ->K                        / (1L+(basis->type == SEASONID));
	I32  Knew     = (newcol->k2_new - newcol->k2_old)/ (1L+(basis->type == SEASONID)) + Kold;
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;

	I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	F32  factor1  = VEC[nSegOld - 1]/ VEC[nSegNew - 1];	
	F64  factor2  = MAT[(nSegOld - 1) * KMAX + Kold - 1] / (MAT[(nSegNew - 1) * KMAX + Knew - 1]) 
		           *  (nSegOld* NUMORDER+1.0)/ (nSegNew * NUMORDER + 1.0);
	 
	F32 factor;
	if (new->jumpType == ChORDER) {
		factor = log(factor2);
		//r_printf("%d %f\n", Knew - Kold, factor),
	} else if (new->jumpType == BIRTH) {
		F32 factor3 = log((F32)basis->goodNum / (basis->nKnot + 1));
		factor1 = log(factor1);
		factor2 = log(factor2);
		//factor  = GLOBAL_MODEL_PRIOR_FACTOR * (2 * basis->nKnot * 0 + 1)+ factor1 + factor2 + factor3;
		// if (basis->nKnot <1) r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);

	} else {
		F32 factor3 = log((F32)(basis->goodNum) / basis->nKnot );  // For death, basis->nKnot >=1
		factor1 = log(factor1);
		factor2 = log(factor2);
		//factor  = -GLOBAL_MODEL_PRIOR_FACTOR * (2 * basis->nKnot * 0 + 1)+ factor1 + factor2 - factor3;
		//r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);
    }

	return (factor);

}
 
static F32 ST_ModelPriorFactor4_old(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {

	I32  N        = newcol->N;
	I32  Kold     = basis ->K / (1 + (basis->type == SEASONID));
	I32  Knew     = Kold + (newcol->k2_new - newcol->k2_old)/ (1 + (basis->type == SEASONID));
	I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
 

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	F32  factor1  = VEC[nSegOld - 1]/ VEC[nSegNew - 1];	
	F64  factor2  = MAT[(nSegOld - 1) * KMAX + Kold - 1] / (MAT[(nSegNew - 1) * KMAX + Knew - 1]) 
		           *  (F32)(nSegOld* NUMORDER+1)/ (nSegNew * NUMORDER + 1);
	 
	F32 factor;
	if (new->jumpType == ChORDER) {
		factor = log(factor2);
		//r_printf("%d %f\n", Knew - Kold, factor),
	} else if (new->jumpType == BIRTH) {
		F32 factor3 = log((F32)basis->goodNum / (basis->nKnot + 1));
		factor1 = log(factor1);
		factor2 = log(factor2);
		factor  = basis->prior.modelComplexity * factor1 + factor2 + factor3;
		//if (basis->nKnot <1) 	r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);

	} else {
		F32 factor3 = log((F32)(basis->goodNum + basis->prior.minSepDist * 2) / basis->nKnot );  // For death, basis->nKnot >=1
		factor1 = log(factor1);
		factor2 = log(factor2);
		factor  = basis->prior.modelComplexity * factor1 + factor2 - factor3;
		//r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);
    }

	return (factor);
}

static F32 ST_ModelPriorFactor5(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {

	I32  N        = newcol->N;
	I32  Kold     = basis ->K / (1 + (basis->type == SEASONID));
	I32  Knew     = Kold + (newcol->k2_new - newcol->k2_old)/ (1 + (basis->type == SEASONID));
	I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
 

	F64PTR MAT = basis->priorMat;
	F64PTR VEC = basis->priorNmodelsPerNseg;

	F32  factor1  = VEC[nSegOld - 1]/ VEC[nSegNew - 1];	
	F64  factor2  = -(nSegNew- nSegOld) *log(NUMORDER+1);
	 
	F32 factor;
	if (new->jumpType == ChORDER) {
		factor = factor2;
		//r_printf("%d %f\n", Knew - Kold, factor),
	} else if (new->jumpType == BIRTH) {
		F32 factor3 = log((F32)basis->goodNum / (basis->nKnot + 1));
		factor1 = log(factor1);
		factor  = basis->prior.modelComplexity * factor1 + factor2 + factor3;
		if (basis->nKnot <1)
		r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);

	} else {
		F32 factor3 = log((F32)(basis->goodNum + basis->prior.minSepDist * 2) / basis->nKnot );  // For death, basis->nKnot >=1
		factor1 = log(factor1);
		factor  = basis->prior.modelComplexity * factor1 + factor2 - factor3;
		//r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);
    }

	return (factor);
}


static F32 ST_ModelPriorFactor6(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {

	I32  N        = newcol->N;
	I32  Kold     = basis ->K / (1 + (basis->type == SEASONID));
	I32  Knew     = Kold + (newcol->k2_new - newcol->k2_old)/ (1 + (basis->type == SEASONID));
	I32  KMAX     = (basis->prior.maxKnotNum + 1L) * (basis->prior.maxOrder + (basis->type == TRENDID));
	I32  nSegOld  = basis->nKnot     + 1L;
	I32  nSegNew  = new  ->nKnot_new + 1L;
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
 
	F64PTR DimVec   = basis->priorVec;
	F64PTR ModelVec = basis->priorNmodelsPerNseg;

	F32  factor1  = DimVec[Kold- 1]/ DimVec[Knew - 1];
	F64  factor2  = ModelVec[nSegOld - 1] / ModelVec[nSegNew - 1];
	
	factor1 = basis->prior.modelComplexity * log(factor1);
	factor2 = log(factor2);

	F32 factor;
	if (new->jumpType == ChORDER) {
		factor = factor1 + factor2;
		//r_printf("%d %f\n", Knew - Kold, factor),
	} else if (new->jumpType == BIRTH) {
		F32 factor3 = log((F32)basis->goodNum / (basis->nKnot + 1));		
		factor  =  factor1 + factor2 + factor3;
		if (basis->nKnot <1)
		r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);

	} else {
		F32 factor3 = log((F32)(basis->goodNum + basis->prior.minSepDist * 2) / basis->nKnot );  // For death, basis->nKnot >=1
		factor  =  factor1 + factor2 - factor3;
		//r_printf("\nbirth %d %f  %f  %f %f\n", basis->nKnot, factor1, factor2, factor3, factor);
    }

	return (factor);
}


void* Get_ModelPrior (I08 id) {
	switch (id) {
	case 0: return ST_ModelPriorFactor0;
	case 1: return ST_ModelPriorFactor1;
	case 2: return ST_ModelPriorFactor1;
	//case 3: return ST_ModelPriorFactor3;
	//case 4: return ST_ModelPriorFactor4;
	//case 5: return ST_ModelPriorFactor5;
	//case 6: return ST_ModelPriorFactor6;
	}
	return NULL;
}

#include "abc_000_warning.h"