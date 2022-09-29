#include <string.h> //memset
#include <math.h>   //log
#include "abc_000_warning.h"

#include "beastv2_header.h"
#include "abc_mcmc.h"
#include "abc_math.h"  // for fastlog only

 

static F64 K_OUT_OF_N(I32 N, I32 k)                     {F64 y = 1;for (I32 i = 0; i <k; i++)	y *= (F64) (N-i) / (F64)(k-i) ;	return y;}
static F64 GetNumModelsGivenK(I32 N, I32 minSep, I32 K) {return K_OUT_OF_N(N - minSep * K - 1, K - 1);}
static F64 GetGroupingNum(I32 n, I32 ng,  I32 k) {
	// math.stackexchange.com/questions/900828/number-of-groups-containing-at-least-1-and-at-most-k-elements/902458#902458

	if (k * ng < n || n < ng) 	   return 0;	
	if (n==ng)					   return 1;	
	if (n < k+ng)				   return K_OUT_OF_N(n - 1, ng - 1);
	
	F64 r= 0;
	for (I32 i = 0; i <= (n/k); ++i) {	r += K_OUT_OF_N(ng, i) * GetGroupingNum(n-i*k, ng-i, k-1);}
	return r;
}

void PreCaclModelNumber(I32 minOrder, I32 maxOrder, I32 maxNumseg, I32 N, I32 minSep,F64PTR TNUM, F64PTR totalNum) {
 
	I32 KMAX = maxNumseg * maxOrder;
	memset(TNUM, 0, KMAX*maxNumseg* sizeof(F64));
	for (I32 ng = 1; ng <= maxNumseg; ng++) {
		for (I32 n = minOrder*ng; n <= maxOrder*ng;  ++n) {
			TNUM[(ng-1)* KMAX +n - 1]=GetGroupingNum(n - (minOrder-1)*ng, ng, maxOrder - (minOrder-1));
		}
	}

	memset(totalNum, 0, sizeof(F64) * KMAX);
	for (I32 ng = 1; ng <= maxNumseg; ++ng) {
		F64 NumModelsPerK = GetNumModelsGivenK(N, minSep, ng);
		for (I32 n = minOrder*ng; n <= maxOrder*ng; ++n) {
			totalNum[n-1] += TNUM[(ng-1L)* KMAX + n - 1] * NumModelsPerK;
		}
	}
	for (I32 n = minOrder*1L; n <= KMAX; ++n) totalNum[n - 1]=-log(totalNum[n - 1]);	
}


static F32 ST_ModelPriorFactor0(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {
	return 0.f; 
}
static F32 ST_ModelPriorFactor1(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new) {
	I32  N    = newcol->N;
	I32  Kold = basis->K;
	I32  Knew = Kold + (newcol->k2_new- newcol->k2_old);
	if (basis->type == SEASONID) { Knew /= 2;	Kold /= 2; }

	I32  Sold = basis->nKnot+1;
	I32  Snew = new->nKnot_new +1;

	I32 O1       = basis->prior.minOrder + (basis->type==TRENDID);
	I32 O2       = basis->prior.maxOrder + (basis->type==TRENDID);
	I32 NUMORDER = basis->prior.maxOrder - basis->prior.minOrder;
 
	F64 OLD = K_OUT_OF_N(Kold - Sold * (O1 - 1) -1, Sold-1);
	F64 NEW = K_OUT_OF_N(Knew - Snew * (O1 - 1) -1, Snew-1);
	F32 factor0 = OLD/NEW  * (F32)(NUMORDER*Sold+1) / (NUMORDER*Snew+1);	

	F32 factor;
	if (new->jumpType == ChORDER)  
		factor=factor0; 
	else if (new->jumpType == BIRTH) {
		factor = (F32)basis->goodNum / (N - basis->nKnot);
		factor = factor0 * factor;
	} else {
		factor = (F32)(basis->goodNum + basis->prior.minSepDist*2)/ (N-new->nKnot_new);
		factor =factor0/factor;
	}
	return (F32) log(factor);
	
}
static F32 ST_ModelPriorFactor2(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol,  NEWTERM_PTR new) {

	I32 N   = newcol->N;

	I32  Kold = basis->K;
	I32  Knew = Kold + (newcol->k2_new-newcol->k2_old);
	if (basis->type == SEASONID) { Knew /= 2;	Kold /= 2; }

	I32  Sold = basis->nKnot+1;
	I32  Snew = new->nKnot_new +1;
 
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder; 
	I32  KMAX     = (basis->prior.maxKnotNum+1L) * (basis->prior.maxOrder+(basis->type ==TRENDID));
 
	F64PTR MAT      = basis->priorMat;
	F32    factor0  =     MAT[(Sold-1) * KMAX + Kold - 1] * (NUMORDER*Sold+1L)
		              / ( MAT[(Snew-1) * KMAX + Knew - 1] * (NUMORDER*Snew+1L));
 

	F32 factor;
	if (new->jumpType == ChORDER)  
		factor=factor0; 
	else if (new->jumpType == BIRTH) {
		factor = (F32)basis->goodNum / (N - basis->nKnot);
		factor = factor0 * factor;
	}
	else {
		factor = (F32)(basis->goodNum + basis->prior.minSepDist*2)/ (N-new->nKnot_new);
		factor =factor0/factor;
	}
	return log(factor);
}


static F32 ST_ModelPriorFactor3(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new)
{   
	    I32  N    = newcol->N ;

		I32  Kold = basis->K;
		I32  Knew = Kold + (newcol->k2_new - newcol->k2_old);
		if (basis->type == SEASONID) { Knew /= 2;	Kold /= 2; }
	    F32 factor0 = basis->priorVec[Knew-1] - basis->priorVec[Kold - 1];

		F32 factor;
		if (new->jumpType == ChORDER)
			factor= factor0;
		else if (new->jumpType == BIRTH) {
			factor = (F32)basis->goodNum / (N - basis->nKnot+1);
			factor = factor0 +log( factor);
		}	else {
			factor = (F32)(basis->goodNum + basis->prior.minSepDist * 2) /(new->nKnot_new+1);
			factor = factor0 -log(factor);
		}
		return (factor);	  
}

static F32 ST_ModelPriorFactor4(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new)
{	
	I32  Kold = basis->K;
	I32  Knew = Kold + (newcol->k2_new-newcol->k2_old);
	if (basis->type == SEASONID) { Knew /= 2;	Kold /= 2; }
	I32  Sold = basis->nKnot+1;
	I32  Snew = new->nKnot_new +1; 
	I32  NUMORDER = basis->prior.maxOrder - basis->prior.minOrder; 
	I32  KMAX     = (basis->prior.maxKnotNum+1L) * (basis->prior.maxOrder+(basis->type ==TRENDID));
	
	F64PTR MAT     = basis->priorMat;
	F32    factor = MAT[(Sold-1)* KMAX + Kold- 1] * (NUMORDER * Sold + 1) / (MAT[(Snew-1) *KMAX + Knew - 1] *(NUMORDER*Snew+1));
	return (F32) logf(factor);	
}

static F32 ST_ModelPriorFactor5(BEAST2_BASIS_PTR basis, NEWCOLINFO_PTR newcol, NEWTERM_PTR new)
{
	F32 factor = 0;
	int delta_k1 = basis->nKnot;
	delta_k1++;
	int delta_k2 = new->nKnot_new;
	delta_k2++;
	I32 Kold = basis->K;
	I32 Knew = basis->K + newcol->k2_new - newcol->k2_old;
	if (delta_k1 == delta_k2 && Kold == Knew)
		factor = 0;
	else {
		int k = min(delta_k1, delta_k2);
		int j;
		int j0;
		j = max(Kold, Knew);
		j0 = min(Kold, Knew);
		if (basis->type == SEASONID) {
			j = j / 2;
			j0 = j0 / 2;
		}

		factor = 1;
		for (rI32 i = 1; i <= k - 1; i++)
			factor = factor * (j - 1 - i + 1) / (j0 - 1 - i + 1);

		factor = factor * (j - 1 - (k - 1)) / k;
		factor = fastlog(factor);
		factor = delta_k1 < delta_k2 ? -factor : factor;
	}
	return factor;
}
/*
static F32 ST_ModelPriorFactor1_old(BEAST2_BASIS_PTR basis, NEWTERM_PTR new, I32 N) {


	F32 factor;
	if (new->jumpType == ChORDER) {
		factor = log( (F32)(new->oldOrder-basis->prior.minOrder + 1)/(new->newOrder-basis->prior.minOrder+1) );
	}
	else {
		I32 order = (new->SEG[0].ORDER2 - new->SEG[0].ORDER1) + 1;
		if (new->jumpType == BIRTH) {

			factor = (F32)basis->goodNum / (N - basis->nKnot) / order;
			factor = log(factor);
		}
		else {
			factor = (F32)(basis->goodNum + basis->prior.minSepDist*2)/ (N - new->nKnot_new) / order;
			factor = -log(factor);
		}

	}

	return factor;
}
static F32 ST_ModelPriorFactor2(BEAST2_BASIS_PTR basis, NEWTERM_PTR new) {
	MOVETYPE flag = new->jumpType;
	F32 factor;
	if (flag == ChORDER || flag == MOVE)
		factor = 0;
	else {
		I32 order = (basis->prior.maxOrder - basis->prior.minOrder + 1);
		factor = order == 1 ? 0 : log(order);
		factor = flag == BIRTH ? -factor : factor;
	}

	return factor;
}


*/

void* Get_ModelPrior (I08 id) {
	switch (id) {
	case 0: return ST_ModelPriorFactor0;
	case 1: return ST_ModelPriorFactor1;
	case 2: return ST_ModelPriorFactor2;
	case 3: return ST_ModelPriorFactor3;
	case 4: return ST_ModelPriorFactor4;
	case 5: return ST_ModelPriorFactor5;
	}
	return NULL;
}

#include "abc_000_warning.h"