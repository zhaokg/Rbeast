#include <string.h>
#include <math.h>
#include "abc_000_warning.h"
#include "abc_mcmc.h"  //RANDINT RANDINT_SKIPINE
#include "beastv2_header.h"

#include "abc_vec.h"   // for i08_sum_binvec only

 

static INLINE void  __CalcAbsDeviation(F32PTR  deviation, F32PTR avgDeviation, PROP_DATA_PTR info, I32 NumBasis) {
		
	F32 invsample = 1.f / info->samples[0];
	I32 N         = info->N;
	I32 q         = info->yInfo->q; //Added for MRBEAST

	F32PTR  Y           = info->yInfo->Y;
	I32PTR  rowsMissing = info->yInfo->rowsMissing;
	I32     nMissing    = info->yInfo->nMissing;		
	F32     nan         = getNaN();

	
	if (NumBasis == 1) {		
		F32PTR  Ypred1 = info->keyresult[0].x;	
		for (int col = 0; col < q; col++) {
				I32     idxMissing = 0;
				F32     sumError = 0;
				for (int i = 0; i < N; i++) {
					if (idxMissing < nMissing && i == rowsMissing[idxMissing]) {  //rowMising is 0-based
						deviation[i] = nan;
						idxMissing++;
					}	else {
						F32 error = fabsf(Y[i] - Ypred1[i] * invsample);
						deviation[i] = error;
						sumError     += error;
					}
				}
				avgDeviation[col] = sumError / info->yInfo->n;
				Y      += N;  //The leading dimenison of Y is N
				Ypred1 += N;
				deviation += N;
			
		}
;
	} 
	else if (NumBasis == 2) {
 
		F32PTR  Ypred1 = info->keyresult[0].x;
		F32PTR  Ypred2 = info->keyresult[1].x;
		for (int col = 0; col < q; col++) {
			I32     idxMissing = 0;
			F32     sumError = 0;
			for (int i = 0; i < N; i++) {
				if (idxMissing < nMissing && (i) == rowsMissing[idxMissing])  //rwoMissing is 0-based
					deviation[i] = nan, idxMissing++;
				else {
					F32 error = fabsf(Y[i] - (Ypred1[i] + Ypred2[i]) * invsample);
					deviation[i] = error;
					sumError += error;
				}
			}
			avgDeviation[col] = sumError / info->yInfo->n;
			Y      += N;  //The leading dimenison of Y is N
			Ypred1 += N;
			Ypred2 += N;
			deviation += N;
		}
		
	} 
	else  {
 
		F32PTR  Ypred1 = info->keyresult[0].x;
		F32PTR  Ypred2 = info->keyresult[1].x;
		F32PTR  Ypred3 = info->keyresult[2].x;
		for (int col = 0; col < q; col++) {
			I32     idxMissing = 0;
			F32     sumError   = 0;
			for (int i = 0; i < N; i++) {
				if (idxMissing<nMissing && (i) == rowsMissing[idxMissing])  ////rowMissing is 0-based
					deviation[i] = nan, 	idxMissing++;
				else {
					F32 error    = fabsf(  Y[i] - (Ypred1[i]+Ypred2[i]+ Ypred3[i]*0) * invsample );
					deviation[i] = error;
					sumError     += error;
				}
			}	
			avgDeviation[col] = sumError / info->yInfo->n;
			Y      += N;  //The leading dimenison of Y is N
			Ypred1 += N;
			Ypred2 += N;
			Ypred3 += N;
			deviation += N;
		}
	}

	// For MBBEAST: compute the relative errors and then sum over the multiple time series
	if (q > 1) {
		deviation = deviation - N * q;
		for (int i = 0; i < q; ++i) {
		// Compute the relative error
			f32_mul_val_inplace(1. / avgDeviation[i], (deviation + i * N), N);
		}
		for (int i = 1; i < q; ++i) {
		// Sum over the multiple time series
			//f32_add_vec_inplace(deviation + i * N, deviation, N);
			F32PTR deviation2 = deviation + i * N;
			for (int j = 0; j < N; ++j) {
				deviation[j] = max(deviation[j], deviation2[j]);
			}
		}

	}
	
}

static INLINE void __CalcExtremKnotPos_ST_BirthOnly(I08PTR extremePosVec, F32PTR deviation, I32 N, F32 threshold) {
	 
	//memset(model->extremePosVec, 0, N); // N+1:Npad16 are set to zero when alllocated
	int i = 0;
	for (; i < N - 3; i += 4) {
		extremePosVec[i]   = deviation[i]   > threshold;
		extremePosVec[i+1] = deviation[i+1] > threshold;
		extremePosVec[i+2] = deviation[i+2] > threshold;
		extremePosVec[i+3] = deviation[i+3] > threshold;
	}	
	for (; i < N ; ++i) 
		extremePosVec[i] = deviation[i] > threshold;

	// No need to compute the number of good extremeKnot Position 
	// because goodVecNum will be recalcuated in the propse function,
	 // model->extremPosNum=i08_sum_binvec(model->extremePosVec, info->Npad16);
}


static INLINE void _CalcDevExtremPos(PROP_DATA_PTR info ) {
	BEAST2_MODEL_PTR model    = info->model;
	I32              NumBasis = model->NUMBASIS;
	//Deviation is needed for sampling Outlier knots but not needed for ST. For ST, only extremPos is needed
	__CalcAbsDeviation( model->deviation, model->avgDeviation,info, NumBasis);
	F32 threshold = info->yInfo->q == 1
		? (model->avgDeviation[0] * info->sigFactor)  // for univariate ts, deviation is the absolute errors
		: ( info->sigFactor);                         // for MRBEAST,      deivation is the summed relative errors
	__CalcExtremKnotPos_ST_BirthOnly(model->extremePosVec, model->deviation, info->N, threshold);
}

static void DSVT_Propose( BEAST2_BASIS_PTR basis, NEWTERM_PTR new, NEWCOLINFO_PTR newcol, PROP_DATA_PTR info)
{	
	I32					goodNum = basis->goodNum;
	I16					nKnot   = basis->nKnot;
	BEAST2_RANDSEEDPTR	PRND	= info->pRND;
	/*********************************************************************/
	//  DETERMINE WHICH PROPOSAL STEPTO MOVE: MAKE A NEW PROPOSAL		
	/**********************************************************************/
	MOVETYPE flag;
	{
		I32  Ktotal                   = info->model->curr.K;
		I32  MINKNOTNUM               = basis->prior.minKnotNum;
		I32  MAXKNOTNUM               = basis->prior.maxKnotNum;
		I32  MAX_K_StopAddingNewTerms = basis->mcmc_Kstopping;
		U08  rnd					  = *(PRND->rnd08)++;

		if (MINKNOTNUM != MAXKNOTNUM) {
				if (rnd < basis->propprob.birth) { 
				// birth:  must be "<" not "<=", so if birht=0, this branch will never chosen
					flag = BIRTH;
					if (nKnot   >= MAXKNOTNUM || goodNum == 0)	flag = MOVE;
					if (Ktotal  >= MAX_K_StopAddingNewTerms  )  flag = (nKnot == 0) ? BIRTH : MOVE;
				}
				else if (rnd < basis->propprob.move)   //move					
					flag = nKnot == 0        ? BIRTH : MOVE;
				else if (rnd < basis->propprob.death)  //death
					flag = nKnot==MINKNOTNUM ? BIRTH : DEATH;
				else if (rnd <= basis->propprob.merge) { //merge  
				// must be "<=" not "<", so if ChOrder=0 (i.e., merge=255), this ChORDER will never chosen					
					if (nKnot == MINKNOTNUM)
						flag = BIRTH;
					else {
						if (nKnot >= 2) flag = MERGE;
						else  			flag = nKnot == 0 ? BIRTH : DEATH;
					}				
				}
				else {                                 //Re-sampling orders					
					if (Ktotal < MAX_K_StopAddingNewTerms)	flag = ChORDER;
					else                                    flag = (nKnot == MINKNOTNUM) ? BIRTH : MOVE;
				}
		}
		else { 	
		// (MINKNOTNUM == MAXKNOTNUM): No birth, merge, and death are allowed at all

			if (MINKNOTNUM == 0) {
			// No changepoint allowed at all, so only ChORDER is possible. 
				if (basis->propprob.merge == 255)   // merge = 255 if maxOrder = minOrder (always the case for DUMMY)
				{
				 // As ensured by PickBasisID, this branch will be visted ONLY if all the non-OUTLIER bases are fixed
			     // to a GLOBAL curve. That is, minKnotN=maxKnotN=0, and minOrder=maxOrder for all non-outlier bases.
			     // For example, if NUMBasis=2 but only one basis is fixed, this branch won't be visited
				// The basis is a global fitting, and so no change is made during the proposal move
				// NoChangeFixGlobal is added to deal with a global fitting with no changepoints at a fixed trend
				// or seasonal order.							 	
					flag = NoChangeFixGlobal; // 				
				}
				else
					// Never visited for DUMMY BASIS
					flag = ChORDER; // minOrder is not equal maxOrder, so we can resample the order				
				
			} else {
			// Only MOVE and CHORDEr allowed for ST, only Move allowed for DUMMY
				if (basis->propprob.merge == 255) 	 // merge = 255 if maxOrder = minOrder (always the case for DUMMY) 			
					flag = MOVE;	 // No ChangeOrder			
				else
					flag = rnd > 128 ? MOVE : ChORDER;
			}
		}

	} // DONE WITH MAKING A NEW PROPOSAL	 

	/*****************************************************************************************/
	// IMPLEMENT THE NEW PROPOSED BASIS: Define local variables here to avoid allocating duplicate stack variables	
	//*****************************************************************************************/
	TKNOT_PTR  knotList  = basis->KNOT;
	TORDER_PTR orderList = basis->ORDER; //DUMMY: Not used for the dummy basis

	
	I32 newIdx, endIdx;
	switch (flag)
	{
	case BIRTH:
	{
		I32    Npad16  = info->Npad16;		
		U64PTR goodVec = basis->goodvec; // Aliased is used here, goodvec must be 8-byte algined
		U64PTR tmpGoodVec;
		I32    tmpGoodNum;
		
		if ( *(PRND->rnd08)++ < 255* PROB_SAMPLE_EXTREME_VECTOR ) {				             

			//  When samples=0, extremPosVec should be all 1s'; it was pre-set when creating the initla model before entering the mcmc loop.			
			I32 samples = info->samples[0];
			//  Do not recompute extremPosVec that often; only do it if the last update is 40 samples apart
			//  Do not compuute extremPos until at least sample is 1: nSample_ExtremVecNeedUpdate must be set to 1L or higher at the start of each chain
			if ( samples >= info->nSample_ExtremVecNeedUpdate) {				
				// Re-calcluate MODEL.extremePosVec, MODEL.deviation, & MODEL.avgDeviation
				_CalcDevExtremPos(info);  
				info->nSample_ExtremVecNeedUpdate = samples+100;
			}
			
			// Aliased is used here, extremPosVec must be 8-byte algined
			U64PTR extremeVec = info->model->extremePosVec;
			tmpGoodVec        = info->mem;   // a temp mem buffer pointed to Xnewterms
			for (I32 i = 0; i < (Npad16/8)-1; i+=2) {
				tmpGoodVec[i]   = extremeVec[i]   & goodVec[i];
				tmpGoodVec[i+1] = extremeVec[i+1] & goodVec[i+1];
			}

			// Number of the good positions is recalcuated here, so no need to cmpute extremPosNum
			// in the __CalcExtremKnotPos_for_ST_BIRTH_ONLY function
			tmpGoodNum = i08_sum_binvec(tmpGoodVec, Npad16);
			if (tmpGoodNum == 0) { tmpGoodVec = goodVec, tmpGoodNum = goodNum; }
			
		}	else {
			tmpGoodVec = basis->goodvec;
			tmpGoodNum = goodNum;
		}

		// Randomly choose a good loc: newKnot is the newly chosen breakpoint	
		I32  randLoc = RANDINT(1, (I32)tmpGoodNum, *(PRND->rnd32)++);
		new->newKnot = i08_find_nth_onebyte_binvec(tmpGoodVec, (I32)Npad16, randLoc);
		//new->newKnot = i08_find_nth_onebyte_binvec_v2(tmpGoodVec, Npad16, tmpGoodNum, *(PRND->rnd32)++ );
		newIdx = 1;	for (TKNOT_PTR tmp=knotList; *tmp++ < new->newKnot; newIdx++);
		//Now, newIdx is the position of the newKnot relative to the old bks.

		//Next, find the index of the newly generated bk-relative to existing bks: newIdx is initiated as the very last bk, but it will be changed in the following loop if not being the last one.	
		/* newIdx = tKnotNum + 1;  rU16PTR knotList = basis->tY.KNOT;  for (rI16 i = 1; i <= tKnotNum; i++)  {  if (newKnot < *knotList++) { newIdx = i; break; } }*/  
		//Now, newIdx is the position of the newKnot relative to the old bks.
		new->numSeg		= 2;
		new->SEG[0].R1	= knotList[(newIdx - 1) - 1];
		new->SEG[0].R2	= new->newKnot - 1;
		new->SEG[1].R1	= new->newKnot;
		new->SEG[1].R2	= knotList[(newIdx)-1] - 1;
		
        // DUMMY: orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[1].ORDER2 = new->SEG[0].ORDER2 = orderList[(newIdx)-1];
	
		endIdx            = newIdx;
		new->newIdx       = newIdx;
		new->nKnot_new    = nKnot + 1;		
		break;
	}
	case DEATH:
	{
		newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);
		new->numSeg     = 1;
		new->SEG[0].R1  = knotList[(newIdx - 1) - 1];
		new->SEG[0].R2  = knotList[(newIdx + 1) - 1] - 1L;
		
        //DUMMY:  orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[0].ORDER2 = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		endIdx            = newIdx + 1L;
		new->newIdx       = newIdx;
		new->nKnot_new = nKnot - 1;
		break;
	}//flag==2/2
	case MERGE:
	{
		newIdx = RANDINT(1, (U16)nKnot - 1, *(PRND->rnd16)++);  //(*rnd32++) % (tKnotNum-1) + 1; // (int)ceilfunc((*rnd32++) *(tKnotNum-1));	
		// count is the number of samples bettween the two MERGEd knots
		I32  r1 = knotList[(newIdx)-1];
		I32  r2 = knotList[(newIdx + 1) - 1];
		I32  count = (r2 - r1) + 1L - 2L;
		if (count == 0L) {
			// count is zero if the two chosen points are next to each other: This possibly happens only if minSeptDist=0
			//https://stackoverflow.com/questions/2733960/pointer-address-type-casting
			new->newKnot = *(*(I08**)&(PRND->rnd08))++ > 0 ? r1 : r2;
		}
		else {
			new->newKnot = RANDINT(r1 + 1, r2 - 1, *(PRND->rnd32)++);  // newKnot = r1 + ((*rnd32++) % count + 1);
		}

		new->numSeg = 2;
		new->SEG[0].R1 = knotList[newIdx - 1L - 1L];
		new->SEG[0].R2 = new->newKnot - 1L;
		new->SEG[1].R1 = new->newKnot;
		new->SEG[1].R2 = knotList[newIdx + 2L - 1L] - 1L;

        //orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[0].ORDER2 = orderList[(newIdx)-1L];
		new->SEG[1].ORDER2 = orderList[newIdx + 2L - 1L];


		endIdx            = newIdx + 2L;
		new->newIdx       = newIdx;
		new->nKnot_new = nKnot - 1;	 		
	 
		break;
	}//flag==merge
	case MOVE: //MOVE
	{
		newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	

		I32 oldKnot = knotList[newIdx - 1];
		I32 r1      = newIdx==1     ?  knotList[INDEX_FakeStart]: knotList[(newIdx - 1) - 1];
		I32 r2      = newIdx==nKnot ?  knotList[INDEX_FakeEnd]  : knotList[(newIdx + 1) - 1];

		I32 minSepDist           = basis->prior.minSepDist;
		I32 MCMC_maxMoveStepSize = basis->mcmc_MoveStep;
		r1 = max(r1 + minSepDist + 1, oldKnot - MCMC_maxMoveStepSize);
		r2 = min(r2 - minSepDist - 1, oldKnot + MCMC_maxMoveStepSize);
		/*??????????????????????????????????*/
		// oldKNot, r1 and r2 should be sigend (though knotList are unsigned); otherwise
		// oldKnot-maxMoveSize could be a negative number and aliased to a super large postiive number
		//, which will cause big r1 or r2. Another possible solution is to use signed intergers for knotList
		/*??????????????????????????????????*/
		if (r2 == r1) {
			new->newKnot = oldKnot;
		} else if (r2 > r1) {
			//rI32 idx = (*rnd32++) % ((r2 - r1 + 1) - 1) + 1;// ceilfunc((*rnd32++) *((r2 - r1 + 1) - 1));
			//newKnot = (r1 - 1) + idx;   newKnot = newKnot < oldKnot ? newKnot : (newKnot + 1L);						   
			RANDINT_SKIPONE(new->newKnot, r1, oldKnot, r2, *(PRND->rnd32)++);
			// r1,..., oldKnot-1,oldKnot(skipped),oldKnot+1,...,r2
		} else {
			r_error("ERROR: r1 < r2; there must be something wrong!\n");
			return ;
		}

		new->numSeg    = 2;
		new->SEG[0].R1 = knotList[newIdx - 1L - 1L];
		new->SEG[0].R2 = new->newKnot - 1L;
		new->SEG[1].R1 = new->newKnot;
		new->SEG[1].R2 = knotList[newIdx + 1L - 1L] - 1L;
		
		
        //orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[0].ORDER2 = orderList[newIdx - 1L];
		new->SEG[1].ORDER2 = orderList[newIdx + 1L - 1L];

		//endIdx = newIdx + 1;
		endIdx         = newIdx + 1L;
		new->newIdx    = newIdx;
		new->nKnot_new = nKnot;

		break;
	}//flag=move MOVE MOVE
	case ChORDER:
	{
	   // this branch should be NEVER visited by DUMMY bases
	   // if (basis->type==DUMMYID) r_error("ChORDER should never been chosen for DUMMYID; there must be something wrong!\n");
	   
		newIdx = RANDINT(1, (U16)nKnot + 1, *(PRND->rnd16)++);  //// (int)ceilfunc((*rnd32++) *(tKnotNum + 1));

		I32 newOrder;
		I32 oldOrder = orderList[newIdx - 1];
		{
			I32 minOrder = basis->prior.minOrder;
			I32 maxOrder = basis->prior.maxOrder;
			if (oldOrder == minOrder)		newOrder = oldOrder + 1;
			else if (oldOrder == maxOrder)	newOrder = oldOrder - 1;
			else           			        newOrder = *(*(I08**)&(PRND->rnd08))++ > 0 ? oldOrder - 1 : oldOrder + 1;
		}

		new->newOrder  = newOrder;
		new->oldOrder  = oldOrder;
		new->SEG[0].R1 = knotList[(newIdx - 1) - 1];
		new->SEG[0].R2 = knotList[(newIdx)-1] - 1L;
		new->SEG[0].ORDER2 = newOrder;//Duplicate there to be used for generating basis function below: to maintan consistent APIs
		//new->orders[0] = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		new->SEG[0].ORDER1 = newOrder; //the starting of the order: not used if decreasing the order

		new->numSeg = newOrder > oldOrder ? 1/*Add new terms*/ : 0/*Rmove old terms*/;
		endIdx = newIdx;// used only if newOrder > oldOrder (add new terms_)	//used for computing newTerm_endIdx

		new->newIdx    = newIdx;
		new->nKnot_new = nKnot;
		break;
	}//flag==4//f
	case NoChangeFixGlobal: {
	   //Added to handle a global model without any changepoint or order-resampling
		new->numSeg  = 0;

		new->newKnot = -9999;
		new->SEG[0].R1 = 0x0fffffff;
		new->SEG[0].R2 = 0x0fffffff;
		new->SEG[1].R1 = 0x0fffffff;
		new->SEG[1].R2 = 0x0fffffff;

		//orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[0].ORDER2 = 0x0fff;
		new->SEG[1].ORDER2 = 0x0fff;

		//endIdx = newIdx + 1;
		//endIdx = newIdx + 1L;
		endIdx = newIdx = 1;
		new->newIdx     = 0x0fff;
		new->nKnot_new  = nKnot;

		break;
	} // flag == 5//f

	}

 


	if (flag != ChORDER) {
		// Even if DUMMY bases can visited this branch, the order values are not used at all
		//orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		
		// basisID is needed to find the starting orders of the newly added segments
		TORDER startOrder = (basis->type == TRENDID) ? 0 : 1;
		new->SEG[1].ORDER1 = new->SEG[0].ORDER1 = startOrder;
	}

	/*******************************************/
	//Get k1_old, k2_old, k1_new and k2_new		
	/*******************************************/
	//I16 k1_old, k2_old;	//The range of terms to be removed from Xt_mars		   			
	//I16 k1_new, k2_new; //The range of terms to be added to Xt_mars_prop			

	I16PTR  KS_old = basis->ks;
	I16PTR  KE_old = basis->ke;

	if (flag != ChORDER) {

		if (flag == NoChangeFixGlobal) {
			//The range of terms to be removed from Xt_mars		   
			// No new term will be added or removed for NoChangeFixGlobal
			// so new->KnewTerm=0;
			newcol->k1     = basis->K+1;
			newcol->k2_old = basis->K;
			//The range of terms to be added to Xt_mars_prop		
			//new->k2_new = new->k1_new + new->Knewterm - 1L;
		}	else {
			//The range of terms to be removed from Xt_mars		   
			newcol->k1     = KS_old[(newIdx)-1];
			newcol->k2_old = KE_old[(endIdx)-1];
			//The range of terms to be added to Xt_mars_prop		
			//new->k2_new = new->k1_new + new->Knewterm - 1L;
		}
		
	}
	else { // if (flag == ChORDER)
	
	   // this branch should be never visited by DUMMY bases
	  // if (basis->type ==DUMMYID) r_error("ChORDER should never been chosen for DUMMYID; there must be something wrong!\n");
		
		if (new->newOrder <= new->oldOrder) { // if (newOrder < oldOrder): remove one or two existing terms
			//NEW.numSeg = 0;	NEW.r[0] = -999;
			newcol->k2_old  = KE_old[newIdx - 1];
			newcol->k1      = basis->type == SEASONID ? /*season*/(newcol->k2_old - 1) :/*trend*/ newcol->k2_old;
			//new->k2_new = new->k1_old - 1;
			//Knewterm = k2_new - k1_new + 1; //Equal to 0, meaning no new terms will be added
		}
		else {// (NEW. newOrder > oldOrder): Add new terms					
			newcol->k2_old = KE_old[newIdx - 1];    //the term immediately preceding the start of the next segment
			newcol->k1     = newcol->k2_old + 1;       //the term immediately following the end of the selected segment

			//new->k2_new = basis->type == SEASONID ? (new->k1_new + 1) : new->k1_new;
			//Knewterm    = k2_new - k1_new + 1; //Knewterm is either 1 or 2, depending on the new term being a trend or season cmpnt
		} // (NEW. newOrder > oldOrder)

	} // if (flag == ChORDER)

	new->jumpType = flag;
   
	//new->SEG[0:1].K is not calculated here and will be obtained by GenTerms()/
}

/*
static void OO_ReAdjustGoodVec(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info) {

	I32 sample = *(info->samples);
	I32 N      = info->N;
	if (sample > 1 && sample % 20 == 0) {

		F32PTR  mem       = info->mem;
		
		F32PTR  Y           = info->yInfo->Y;
		U32PTR  rowsMissing = info->yInfo->rowsMissing;
		I32     nMissing    = info->yInfo->nMissing;
		I32     idxMissing  = 0;
		F32     nan         = getNaN();
		CORESULT* result = info->keyresult;

		F32   invsample = 1. / sample;
		F32   errorSum  = 0;
		
		for (int i = 1; i <= N; i++) {

			if (idxMissing < nMissing && i == rowsMissing[idxMissing]) {
				mem[i - 1] = nan;
				idxMissing++;
			}else	{
				F32 error = fabs(Y[i - 1] - (result[0].x[i - 1] + result[1].x[i - 1] + result[2].x[i - 1]) * invsample);
				mem[i - 1] = error;
				errorSum   += error;
			}
		}


		F32 sig = errorSum / info->yInfo->n;
		sig = 1. * sig;
		memset(basis->goodvec, 0, N);
		for (int i = 1; i <= N; i++) {
			if (mem[i - 1] > sig) {
				basis->goodvec[i - 1] = 1;
			}
		}
	
	}

	BEAST2_MODEL_PTR model = info->model;
	for (int J = 0; J < model->NUMBASIS; J++) {

		TKNOT_PTR KNOT = model->b[J].KNOT;
		for (int i = 0; i < model->b[J].nKnot; i++) {
			I32 idx = KNOT[i] - 1;
			basis->goodvec[idx] = 0;

			basis->goodvec[max(idx-1,0)] = 0;
			basis->goodvec[min(idx+1, N-1)] = 0;
		}
	}
	
	basis->goodNum = i08_sum_binvec(basis->goodvec, info->Npad16);

}
*/


static int __OO_NewKnot_BirthMove_old(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info) {
		
	I32 N                   = info->N;
	I32 Npad16              = info->Npad16;
	BEAST2_MODEL_PTR model  = info->model;
	BEAST2_RANDSEEDPTR PRND = info->pRND;

	// TODO: reconstruct the goodVec array on-the-fly, which means the existing goodVec is dicarded.
	I08PTR goodvec = (I08PTR) basis->goodvec; 
	memset(goodvec, 1, N);

	for (int J = 0; J < model->NUMBASIS; J++) {

		TKNOT_PTR KNOT    = model->b[J].KNOT;
		I32       nKnot   = model->b[J].nKnot;

		if (model->b[J].type == OUTLIERID) {
			// Exclude the current outlier chgt
			for (int i = 0; i < nKnot; ++i) goodvec[KNOT[i] - 1] = 0;		
		} else {
		// TODO: do we really need this extra step?
		// Exclude the current season/trend changepoints and their immediate neighgors
	     	/*
			for (int i = 0; i < nKnot; i++) {
				I32 idx = KNOT[i] - 1;
				goodvec[idx]                 = 0;
				goodvec[max(idx - 1, 0)]     = 0;
				goodvec[min(idx + 1, N - 1)] = 0;
			}
			*/
		// 2022/3/5: The above turns out to be a very bad strategy as the tcps/scps will steal
		 // the chances for a real ocp. Here we go to the contraty and istead include all of thems
	    // in the canidate ocps.
			if ((*(PRND->rnd08)++) > 100) {  
				// Wo do this randomly
				for (I32 i = 0; i < nKnot; i++) {
					I32 idx = KNOT[i] - 1;
					goodvec[idx]                 = 2L;
					goodvec[max(idx - 1, 0)]     = 2L;
					goodvec[min(idx + 1, N - 1)] = 2L;
				}
			}//if ((*(PRND->rnd08)++) > 130)
			
		}		 

	} // for (int J = 0; J < model->NUMBASIS; J++)

	I32PTR IndicesLargeDeviation = info->mem;   // mem points to Xnewterm, used here a tmp buff
	I32    numLargeDev       = 0;
	F32PTR deviation         = model->deviation;
	F32    threshold         = info->yInfo->q==1
							   ? model->avgDeviation[0]* info->outlierSigFactor
							   : info->outlierSigFactor;

	F32 maxValue = 0;
	I32 maxIdx   =-1;	
	
	for (I32 i = 0; i < N; i++) {
		F32 value = deviation[i];
		if ( !goodvec[i] || IsNaN (value) ) {
			continue;
		}

		value = fabsf(value);
		if ( value > maxValue) {
			// if no points have DEV larger than the threshold, then choose the point with the largest DEV
			maxValue = value;
			maxIdx   = i;
		}

		if (value > threshold || goodvec[i]==2) {
			// if large than the threshold or already being a TCP or SCP
			IndicesLargeDeviation[numLargeDev++] = i;
		}
	}

	if (numLargeDev > 1) {	
		if ((*(PRND->rnd08)++) > 50) {
			I32 rndIdx = RANDINT(1, (U16)numLargeDev, *(PRND->rnd16)++);
			maxIdx = IndicesLargeDeviation[rndIdx - 1];
		} else { // a fixed probability to choose the largest deviation
			maxIdx = maxIdx;			
		}
	}

	if (maxIdx < 0) {
		r_printf("__OO_NewKnot_BirthMove: maxIdx=-1, and there must be something wrong!");
	}
	return maxIdx + 1;
}

static int __OO_NewKnot_BirthMove(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info, I32PTR maxIndex) {
		
	I32 N                   = info->N;
	I32 Npad16              = info->Npad16;
	BEAST2_MODEL_PTR model  = info->model;
	BEAST2_RANDSEEDPTR PRND = info->pRND;

	// TODO: reconstruct the goodVec array on-the-fly, which means the existing goodVec is dicarded.
	I08PTR goodvec = (I08PTR) basis->goodvec; 
	memset(goodvec, 1, N);

	for (int J = 0; J < model->NUMBASIS; J++) {

		TKNOT_PTR KNOT    = model->b[J].KNOT;
		I32       nKnot   = model->b[J].nKnot;

		if (model->b[J].type == OUTLIERID) {
			// Exclude the current outlier chgt
			for (int i = 0; i < nKnot; ++i) goodvec[KNOT[i] - 1] = 0;		
		}  
	} // for (int J = 0; J < model->NUMBASIS; J++)

	I32PTR IndicesLargeDeviation = info->mem;   // mem points to Xnewterm, used here a tmp buff
	I32    numLargeDev       = 0;
	F32PTR deviation         = model->deviation;
	F32    threshold         = info->yInfo->q==1
							   ? model->avgDeviation[0]* info->outlierSigFactor
							   : info->outlierSigFactor;

	F32 maxValue = 0;
	I32 maxIdx   =-1;	
	
	for (I32 i = 0; i < N; i++) {
		F32 value = deviation[i];
		if ( !goodvec[i] || IsNaN (value) ) {
			continue;
		}

		value = fabsf(value);
		if ( value > maxValue) {
			// if no points have DEV larger than the threshold, then choose the point with the largest DEV
			maxValue = value;
			maxIdx   = i;
		}

		if (value > threshold) {
			// if large than the threshold or already being a TCP or SCP
			IndicesLargeDeviation[numLargeDev++] = i;
		}
	}

	int newKnot=-1L;
	if (numLargeDev > 1) {	
		I32 rndIdx = RANDINT(1, (U16)numLargeDev, *(PRND->rnd16)++);
		newKnot = IndicesLargeDeviation[rndIdx - 1];
	} else if (numLargeDev == 1) {
		newKnot = IndicesLargeDeviation[0];
	}

	if (maxIdx < 0) {
		r_printf("__OO_NewKnot_BirthMove: maxIdx=-1, and there must be something wrong!");
	}
	*maxIndex = maxIdx+1L;
	return newKnot+1L;

}

static int __OO_NewIdx_MoveDeath(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info) {
		
	I32 N                  = info->N;
	I32 Npad16             = info->Npad16;
	BEAST2_MODEL_PTR model = info->model;
	 
	F32PTR deviation = model->deviation;
	F32    minValue  = 1e34;
	I32    minIdx    = -1;

	I32       nKnot = basis->nKnot;
	TKNOT_PTR KNOT  = basis->KNOT;
	for (int i = 0; i < nKnot; i++) {
	// Choose the ocp with the leasst deviation
		I32 idx   = KNOT[i] - 1;
		F32 value = fabsf(deviation[idx]);
 	 
		if (minValue > value) {
			minValue = value;
			minIdx   = i;
		}
	}
	if (minIdx < 0) {
		r_printf("__OO_NewKnot_BirthMove: maxIdx=-1, and there must be something wrong!");
	}
	return minIdx + 1;
}

static void OO_Propose_01_old(	BEAST2_BASIS_PTR basis, NEWTERM_PTR new, NEWCOLINFO_PTR newcol,PROP_DATA_PTR info)
{	
	I32  goodNum = basis->goodNum; // not used in this function
	I16  nKnot = basis->nKnot;

	BEAST2_RANDSEEDPTR PRND = info->pRND;
	/*********************************************************************/
	//  DETERMINE WHICH PROPOSAL STEPTO MOVE: MAKE A NEW PROPOSAL		
	/**********************************************************************/
	MOVETYPE flag;
	{
		I32  Ktotal                   = info->model->curr.K;
		I32  MAXKNOTNUM               = basis->prior.maxKnotNum;
		I32  MAX_K_StopAddingNewTerms = basis->mcmc_Kstopping;
		U08  rnd                      = *(PRND->rnd08)++;
		if (rnd < basis->propprob.birth) { //birth					
			flag = BIRTH;
			if (nKnot >= MAXKNOTNUM )	           flag = MOVE;			
			if (Ktotal > MAX_K_StopAddingNewTerms) flag = (nKnot == 0) ? BIRTH : MOVE;			
		} else if (rnd < basis->propprob.move)   //move							
			flag = nKnot == 0 ? BIRTH : MOVE;		
		else //if (unifRnd <= basis->propprob.death)  //death
			flag = nKnot == 0 ? BIRTH : DEATH;
 
	}  // DONE WITH MAKING A NEW PROPOSAL	 

	////////////////////////////////////////////////////////////////////
	I32              samples = info->samples[0];	
	if (samples > 0) { // this should be always the case bcz OID is picked up only if samples >0
	 //  Here, we need only the deviation for the outlire proposal but as a side effect,
	 //  we also update the extremKnotPos vector.
		_CalcDevExtremPos(info);
		info->nSample_ExtremVecNeedUpdate = samples + 40L;
	}
	////////////////////////////////////////////////////////////////////


	/*****************************************************************************************/
	// IMPLEMENT THE NEW PROPOSED BASIS
	// Define local variables here to avoid allocating duplicate stack variables	
	//*****************************************************************************************/
	TKNOT_PTR  knotList = basis->KNOT;
	//TORDER_PTR orderList = basis->ORDER;  

	I32 newIdx;      // , endIdx;
 		
	/*****************************************************************************************/
	
	switch (flag)
	{
	case BIRTH:
	{
		
		new->newKnot   = __OO_NewKnot_BirthMove_old(basis, info);
		new->numSeg    = 1;
		new->SEG[0].R1 = new->newKnot;
		new->SEG[0].R2 = new->newKnot;
		new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

		//new->orders2[1] = new->orders2[0] = orderList[(newIdx)-1];
	
		//implicit conversion from 'int' to 'I16' (aka 'short') changes value from -9999999 to 27009 [-Wconstant-conversion]
		//new->newIdx       = -9999999;            // not used at alll
		new->newIdx       = -9999;            // not used at alll
		new->nKnot_new = nKnot + 1;

		break;
	}
	case DEATH:
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);
		newIdx  = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot = knotList[newIdx - 1];
		new->numSeg  = 0;
		
		//new->r1[0] = knotList[(newIdx - 1) - 1];
		//new->r2[0] = knotList[(newIdx + 1) - 1]  -1L;
		//new->orders2[0] = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		//endIdx = newIdx;
		//endIdx = newIdx + 1L;
		new->newIdx       = newIdx;
		new->nKnot_new = nKnot - 1;
		break;
	}//flag==2/2
	case MOVE: //MOVE
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	

		// Randomly choose a good loc: newKnot is the newly chosen breakpoint		
		//rI32 randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++); // will fail if goodNum=0
		//I32  Npad16  = info->Npad16;
		//new->newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, (I32)Npad16, randLoc);

		newIdx         = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot   = __OO_NewKnot_BirthMove_old(basis, info);;
		new->numSeg = 1;
		new->SEG[0].R1         = new->newKnot;
		new->SEG[0].R2         = new->newKnot;
		new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

		//endIdx = newIdx + 1;
		//endIdx = newIdx + 1L;
		new->newIdx = newIdx;
		new->nKnot_new = nKnot;
		break;
	}//flag=move MOVE MOVE

	}

 
	/*******************************************/
	//Get k1_old, k2_old, k1_new and k2_new		
	/*******************************************/
 
	//I16 k1_old, k2_old;	//The range of terms to be removed from Xt_mars		   			
	//I16 k1_new, k2_new; //The range of terms to be added to Xt_mars_prop			

	I16PTR  KS_old = basis->ks;
	I16PTR  KE_old = basis->ke;

	if (flag == BIRTH) {
		//The range of terms to be removed from Xt_mars		
		I32 nKnot = basis->nKnot;
		newcol->k2_old = KE_old[nKnot - 1];
		newcol->k1     = newcol->k2_old + 1;
	}
	else if (flag == DEATH) {
		newcol->k2_old = KE_old[newIdx - 1];
		newcol->k1     = KS_old[newIdx - 1];
	}
	else if (flag == MOVE) {
		newcol->k2_old = KE_old[newIdx - 1];
		newcol->k1     = KS_old[newIdx - 1];
	}

	new->jumpType = flag;
}

static void OO_Propose_01(	BEAST2_BASIS_PTR basis, NEWTERM_PTR new, NEWCOLINFO_PTR newcol, PROP_DATA_PTR info)
{	
	I32  goodNum = basis->goodNum; // not used in this function
	I16  nKnot = basis->nKnot;

	BEAST2_RANDSEEDPTR PRND = info->pRND;
	/*********************************************************************/
	//  DETERMINE WHICH PROPOSAL STEPTO MOVE: MAKE A NEW PROPOSAL		
	/**********************************************************************/
	MOVETYPE flag;
	{
		I32  Ktotal                   = info->model->curr.K;
		I32  MAXKNOTNUM               = basis->prior.maxKnotNum;
		I32  MAX_K_StopAddingNewTerms = basis->mcmc_Kstopping;
		U08  rnd                      = *(PRND->rnd08)++;
		if (rnd < basis->propprob.birth) { //birth					
			flag = BIRTH;
			if (nKnot >= MAXKNOTNUM )	           flag = MOVE;			
			if (Ktotal > MAX_K_StopAddingNewTerms) flag = (nKnot == 0) ? BIRTH : MOVE;			
		} else if (rnd < basis->propprob.move)   //move							
			flag = nKnot == 0 ? BIRTH : MOVE;		
		else //if (unifRnd <= basis->propprob.death)  //death
			flag = nKnot == 0 ? BIRTH : DEATH;
 
	}  // DONE WITH MAKING A NEW PROPOSAL	 

	////////////////////////////////////////////////////////////////////
	I32              samples = info->samples[0];	
	if (samples > 0) { // this should be always the case bcz OID is picked up only if samples >0
	 //  Here, we need only the deviation for the outlire proposal but as a side effect,
	 //  we also update the extremKnotPos vector.
		_CalcDevExtremPos(info);
		info->nSample_ExtremVecNeedUpdate = samples + 40L;
	}
	////////////////////////////////////////////////////////////////////


	/*****************************************************************************************/
	// IMPLEMENT THE NEW PROPOSED BASIS
	// Define local variables here to avoid allocating duplicate stack variables	
	//*****************************************************************************************/
	TKNOT_PTR  knotList = basis->KNOT;
	//TORDER_PTR orderList = basis->ORDER;  

	I32 newIdx;      // , endIdx;
 		
	/*****************************************************************************************/
	I32 maxIdx;
	switch (flag)
	{
	case BIRTH:
	{
		
		new->newKnot   = __OO_NewKnot_BirthMove(basis, info,&maxIdx);
		if (new->newKnot == 0 && nKnot==0) {
			new->newKnot = maxIdx;
		}
		if (new->newKnot > 0) {
			new->numSeg = 1;
			new->SEG[0].R1 = new->newKnot;
			new->SEG[0].R2 = new->newKnot;
			new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

			//new->orders2[1] = new->orders2[0] = orderList[(newIdx)-1];

			//implicit conversion from 'int' to 'I16' (aka 'short') changes value from -9999999 to 27009 [-Wconstant-conversion]
			//new->newIdx       = -9999999;            // not used at alll
			new->newIdx = -9999;            // not used at alll
			new->nKnot_new = nKnot + 1;
		}	else {
			// Death
			flag = DEATH;			 
		}
		
		break;
	}

	case MOVE: //MOVE
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	

		// Randomly choose a good loc: newKnot is the newly chosen breakpoint		
		//rI32 randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++); // will fail if goodNum=0
		//I32  Npad16  = info->Npad16;
		//new->newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, (I32)Npad16, randLoc);

		newIdx         = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot   = __OO_NewKnot_BirthMove(basis, info,&maxIdx);
		if (new->newKnot > 0) {
			new->numSeg = 1;
			new->SEG[0].R1 = new->newKnot;
			new->SEG[0].R2 = new->newKnot;
			new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

			//endIdx = newIdx + 1;
			//endIdx = newIdx + 1L;
			new->newIdx = newIdx;
			new->nKnot_new = nKnot;
			
		}	else {
			flag = DEATH;
		}
		break;
	}//flag=move MOVE MOVE

	}

 	if (flag== DEATH)
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);
		newIdx = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot = knotList[newIdx - 1];
		new->numSeg = 0;

		//new->r1[0] = knotList[(newIdx - 1) - 1];
		//new->r2[0] = knotList[(newIdx + 1) - 1]  -1L;
		//new->orders2[0] = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		//endIdx = newIdx;
		//endIdx = newIdx + 1L;
		new->newIdx = newIdx;
		new->nKnot_new = nKnot - 1;
	 
	}//flag==2/2

	/*******************************************/
	//Get k1_old, k2_old, k1_new and k2_new		
	/*******************************************/
 
	//I16 k1_old, k2_old;	//The range of terms to be removed from Xt_mars		   			
	//I16 k1_new, k2_new; //The range of terms to be added to Xt_mars_prop			

	I16PTR  KS_old = basis->ks;
	I16PTR  KE_old = basis->ke;

	if (flag == BIRTH) {
		//The range of terms to be removed from Xt_mars		
		I32 nKnot = basis->nKnot;
		newcol->k2_old = KE_old[nKnot - 1];
		newcol->k1     = newcol->k2_old + 1;
	}
	else if (flag == DEATH) {
		newcol->k2_old  = KE_old[newIdx - 1];
		newcol->k1      = KS_old[newIdx - 1];
	}
	else if (flag == MOVE) {
		newcol->k2_old = KE_old[newIdx - 1];
		newcol->k1     = KS_old[newIdx - 1];
	}

	new->jumpType = flag;
}

static void OO_Propose_2(	BEAST2_BASIS_PTR basis, NEWTERM_PTR new, NEWCOLINFO_PTR newcol, PROP_DATA_PTR info)
{	
	I32  goodNum = basis->goodNum; // not used in this function
	I16  nKnot = basis->nKnot;

	BEAST2_RANDSEEDPTR PRND = info->pRND;
	/*********************************************************************/
	//  DETERMINE WHICH PROPOSAL STEPTO MOVE: MAKE A NEW PROPOSAL		
	/**********************************************************************/
	MOVETYPE flag;
	{
		I32  Ktotal                   = info->model->curr.K;
		I32  MAXKNOTNUM               = basis->prior.maxKnotNum;
		I32  MAX_K_StopAddingNewTerms = basis->mcmc_Kstopping;
		U08  rnd                      = *(PRND->rnd08)++;
		if (rnd < basis->propprob.birth) { //birth					
			flag = BIRTH;
			if (nKnot >= MAXKNOTNUM )	           flag = MOVE;			
			if (Ktotal > MAX_K_StopAddingNewTerms) flag = (nKnot == 0) ? BIRTH : MOVE;			
		} else if (rnd < basis->propprob.move)   //move							
			flag = nKnot == 0 ? BIRTH : MOVE;		
		else //if (unifRnd <= basis->propprob.death)  //death
			flag = nKnot == 0 ? BIRTH : DEATH;
 
	}  // DONE WITH MAKING A NEW PROPOSAL	 


	////////////////////////////////////////////////////////////////////
	I32              samples = info->samples[0];	
	if (samples > 0) { // this should be always the case bcz OID is picked up only if samples >0
	 //  Here, we need only the deviation for the outlire proposal but as a side effect,
	 //  we also update the extremKnotPos vector.
		_CalcDevExtremPos(info);
		info->nSample_ExtremVecNeedUpdate = samples + 40L;
	}
	////////////////////////////////////////////////////////////////////


	/*****************************************************************************************/
	// IMPLEMENT THE NEW PROPOSED BASIS
	// Define local variables here to avoid allocating duplicate stack variables	
	//*****************************************************************************************/
	TKNOT_PTR  knotList = basis->KNOT;
	//rTORDER_PTR orderList = basis->ORDER;  

	I16 newIdx;      
	switch (flag)
	{
	case BIRTH:
	{
		// Randomly choose a good loc: newKnot is the newly chosen breakpoint		
		//rI32 randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++);
		//I32  Npad16  = info->Npad16;
		//new->newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, (I32)Npad16, randLoc);

		newIdx       = -9999;
		new->newKnot = __OO_NewKnot_BirthMove_old(basis, info);
		
		//newIdx = 1;	for (rTKNOT_PTR tmpList = knotList; *tmpList++ < new->newKnot; newIdx++);
		//Now, newIdx is the position of the newKnot relative to the old bks.
		new->numSeg    = 1;
		new->SEG[0].R1 = 1;       //new->newKnot;
		new->SEG[0].R2 = info->N; // new->newKnot;
		new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

		//new->orders2[1] = new->orders2[0] = orderList[(newIdx)-1];

		new->newIdx       = newIdx;
		new->nKnot_new = nKnot + 1;
		break;
	}
	case DEATH:
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);

		newIdx = __OO_NewIdx_MoveDeath(basis, info);		
		new->newKnot = knotList[newIdx - 1];
		new->numSeg = 0;
		//new->r1[0] = knotList[(newIdx - 1) - 1];
		//new->r2[0] = knotList[(newIdx + 1) - 1]  -1L;
		//new->orders2[0] = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		//endIdx = newIdx;
		//endIdx = newIdx + 1L;
		new->newIdx       =newIdx;
		new->nKnot_new = nKnot - 1;
		break;
	}//flag==2/2
	case MOVE: //MOVE
	{
		//newIdx = RANDINT(1, (U16)nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	


		// Randomly choose a good loc: newKnot is the newly chosen breakpoint		
		//rI32 randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++); // will fail if goodNum=0
		//I32  Npad16  = info->Npad16;
		//new->newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, (I32)Npad16, randLoc);
		newIdx       = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot = __OO_NewKnot_BirthMove_old(basis, info); 

		new->numSeg    = 1;
		new->SEG[0].R1 = 1;       // new->newKnot;
		new->SEG[0].R2 = info->N; // new->newKnot;
		new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

		//endIdx = newIdx + 1;
		//endIdx = newIdx + 1L;
		new->newIdx      =  newIdx;
		new->nKnot_new = nKnot;
		break;
	}//flag=move MOVE MOVE

	}
	/*******************************************/
	//Get k1_old, k2_old, k1_new and k2_new		
	/*******************************************/

	//I16 k1_old, k2_old;	//The range of terms to be removed from Xt_mars		   			
	//I16 k1_new, k2_new; //The range of terms to be added to Xt_mars_prop			

	I16PTR  KS_old = basis->ks;
	I16PTR  KE_old = basis->ke;

	if (flag == BIRTH) {
		//The range of terms to be removed from Xt_mars		
		I32 nKnot = basis->nKnot;
		newcol->k2_old = KE_old[nKnot - 1];
		newcol->k1     = newcol->k2_old + 1;
	}
	else if (flag == DEATH) {
		newcol->k2_old  = KE_old[newIdx - 1];
		newcol->k1      = KS_old[newIdx - 1];
	}
	else if (flag == MOVE) {
		newcol->k2_old = KE_old[newIdx - 1];
		newcol->k1     = KS_old[newIdx - 1];
	}

	new->jumpType = flag;
}

void * Get_Propose(I08 id, BEAST2_OPTIONS_PTR opt) {

	switch (id) {
	case SVDID:
	case DUMMYID:
	case SEASONID: 
	case TRENDID: 
		return DSVT_Propose;
	case OUTLIERID: {
		if      (opt->prior.outlierBasisFuncType==0) 			return OO_Propose_01;
		else if (opt->prior.outlierBasisFuncType == 1)			return OO_Propose_01;
	}
	}
	return NULL;
}

#include "abc_000_warning.h"