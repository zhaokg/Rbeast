#include <string.h>
#include <math.h>
#include "abc_000_warning.h"
#include "abc_mcmc.h"  //RANDINT RANDINT_SKIPINE
#include "abc_ide_util.h"  
#include "beastv2_header.h"

#include "abc_vec.h"   // for i08_sum_binvec only

static  void  __CalcAbsDeviation(F32PTR  deviation, F32PTR avgDeviation, PROP_DATA_PTR info, I32 NumBasis) {
		
	F32     invsample   = 1.f / info->samples[0];
	I32     N           = info->N;
	I32     q           = info->yInfo->q; //Added for MRBEAST

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
				Y         += N;  //The leading dimenison of Y is N
				Ypred1    += N;
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
			f32_mul_val_inplace(1. / avgDeviation[i], (deviation + i * N), N); 	// Compute the relative error
		}
		for (int i = 1; i < q; ++i) {
		    // Sum over the multiple time series
			// f32_add_vec_inplace(deviation + i * N, deviation, N);
			F32PTR relDeviation2 = deviation + i * N;
			for (int j = 0; j < N; ++j) {
				deviation[j] = max(deviation[j], relDeviation2[j]);
			}
		}

	}
	
}

static INLINE void __CalcExtremKnotPos(I08PTR extremePosVec, F32PTR deviation, I32 N, F32 threshold) {
	 
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

	// Note; No need to compute the number of good extremeKnot Position  because 
	// goodVecNum will be recalcuated in the propse function,
	// model->extremPosNum=i08_sum_binvec(model->extremePosVec, info->Npad16);
}
 
static void DSVT_Propose( BEAST2_BASIS_PTR basis, NEWTERM_PTR new,  PROP_DATA_PTR info)
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
		I32    Npad16  = info ->Npad16;		
		U64PTR goodVec = basis->goodvec; // Aliased is used here, goodvec must be 8-byte algined
		U64PTR tmpGoodVec;
		I32    tmpGoodNum;
		
		if ( *(PRND->rnd08)++ < 255 * PROB_SAMPLE_EXTREME_VECTOR ) {				             

			//  When samples=0, extremPosVec should be all 1s'; it was pre-set when creating the initla model before entering the mcmc loop.			
			I32 samples = info->samples[0];

			//  Re-calcluate MODEL.extremePosVec, MODEL.deviation, & MODEL.avgDeviation	
			//  but do not recompute them that often; only if the last update is 40 samples apart	
			//  deviation may be also updated in the outlier proposal function		
			if ( samples >= info->nSample_DeviationNeedUpdate) {
							
				BEAST2_MODEL_PTR model    = info->model;		
				// Deviation is needed for sampling Outlier knots but not needed for ST. For ST, only extremPos is needed
				__CalcAbsDeviation(model->deviation, model->avgDeviation, info, info->numBasisWithoutOutlier);
	
				// Do not compuute extremPos until at least sample is 1: nSample_ExtremVecNeedUpdate must be 
				// set to 1L or higher at the start of each chain
				I32  extraSamples = min( (10 + samples / 8), 200);
				info->nSample_DeviationNeedUpdate = samples+ extraSamples;
				info->shallUpdateExtremVec        =1L;
			}
			
			// Update extremVec using shallUpdateExtremVec==1 because Deviation can be updated here or in the outlier proposal
			if (info->shallUpdateExtremVec) {
				BEAST2_MODEL_PTR model = info->model;
				// For univariate ts, deviation is the absolute errors. For MRBEAST, it ts abd is the summed relative errors
				F32 threshold = (info->yInfo->q == 1) ? (model->avgDeviation[0] * info->sigFactor) : info->sigFactor;
				__CalcExtremKnotPos(model->extremePosVec, model->deviation, info->N, threshold);

				info->shallUpdateExtremVec = 0L;
			}

			
			// Aliased is used here, extremPosVec must be 8-byte algined
			U64PTR extremeVec = info->model->extremePosVec;
			tmpGoodVec        = info->mem;   // a temp mem buffer pointed to Xnewterms
			for (I32 i = 0; i < (Npad16/8)-1; i+=2) {
				tmpGoodVec[i]   = extremeVec[i]   & goodVec[i];
				tmpGoodVec[i+1] = extremeVec[i+1] & goodVec[i+1];
			}

			// Number of the good positions is recalcuated here, so no need to cmpute extremPosNum in  __CalcExtremKnotPos_for_ST_BIRTH_ONLY
			tmpGoodNum = i08_sum_binvec(tmpGoodVec, Npad16);

			if (tmpGoodNum == 0) { 
				tmpGoodVec = goodVec, tmpGoodNum = goodNum; 
			}

			
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
		newIdx = RANDINT(1,  nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);
		new->numSeg     = 1;
		new->SEG[0].R1  = knotList[(newIdx - 1) - 1];
		new->SEG[0].R2  = knotList[(newIdx + 1) - 1] - 1L;
		
        //DUMMY:  orders are not used for DUMMY at all. "if (basis->type==DUMMYID)" is not used here, istead
		// garbage values from orderList are copied for DUMMY basis. THis require allocating mem for
		// basis->ORDER in ST_BASIS_InitializeAllocateMEM, though the mem is ever used
		new->SEG[0].ORDER2 = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		endIdx            = newIdx + 1L;
		new->newIdx       = newIdx;
		new->nKnot_new    = nKnot - 1;
		break;
	}//flag==2/2
	case MERGE:
	{
		newIdx = RANDINT(1, nKnot - 1, *(PRND->rnd16)++);  //(*rnd32++) % (tKnotNum-1) + 1; // (int)ceilfunc((*rnd32++) *(tKnotNum-1));	
		// count is the number of samples bettween the two MERGEd knots
		I32  r1 = knotList[(newIdx)-1];
		I32  r2 = knotList[(newIdx + 1) - 1];
		I32  count = (r2 - r1) + 1L - 2L;
		if (count == 0L) {
			// count is zero if the two chosen points are next to each other: This possibly happens only if minSeptDist=0
			//https://stackoverflow.com/questions/2733960/pointer-address-type-casting
			new->newKnot = *(*(I08**)&(PRND->rnd08))++ > 0 ? r1 : r2;
		}	else {
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
		newIdx      = RANDINT(1, nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	

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
		// Jan 27-2004: Change knotList from u32 (unsigned) to signed (i32) because KNOT[INDEX_FakeStart]
		// may be negative. (this shouldn't be a problem bcz int r1 = (uint) -2 is still -2
		/*??????????????????????????????????*/
		 

		if (r2 == r1) {
			new->newKnot = oldKnot;
		} else if (r2 > r1) {
			//rI32 idx = (*rnd32++) % ((r2 - r1 + 1) - 1) + 1;// ceilfunc((*rnd32++) *((r2 - r1 + 1) - 1));
			//newKnot = (r1 - 1) + idx;   newKnot = newKnot < oldKnot ? newKnot : (newKnot + 1L);						   
			RANDINT_SKIPONE(new->newKnot, r1, oldKnot, r2, *(PRND->rnd32)++);
			// r1,..., oldKnot-1,oldKnot(skipped),oldKnot+1,...,r2
		} else {
			/*
			r_printf("ERROR: r1 < r2; this should never happen and there must be something wrong!\n");
			r_printf("More info for diagnostic purpuse:\n");
			r_printf("r1=%d r2=%d minSepDist=%d, maxMoveStepSize=%d \n", r1, r2, minSepDist, MCMC_maxMoveStepSize);
			r_printf("INDEX_FakeStart=%d  INDEX_FakeEnd=%d newIdx=%d  oldknot=%d \n", INDEX_FakeStart, INDEX_FakeEnd, newIdx, oldKnot);
			r_printf("KnotList (nKnot=%d): \n", nKnot);
			for (int i = INDEX_FakeStart; i <= nKnot; i++) {
				r_printf("  [%d, %d] ", i, knotList[i]);
				if ((i - INDEX_FakeStart + 1) % 10 == 0) {
					r_printf("\n");
				}
			}
			r_printf("\n"); 
			StdouFlush();
			*/
			r_error("ERROR: r1 < r2; this should never happen and there must be something wrong!\n");			
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
	   
		newIdx = RANDINT(1, nKnot + 1, *(PRND->rnd16)++);  //// (int)ceilfunc((*rnd32++) *(tKnotNum + 1));

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

		new->newKnot   = -9999;
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
	// Get k1_old, k2_old, k1_new and k2_new		
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
			new->newcols.k1     = basis->K+1;
			new->newcols.k2_old = basis->K;
			//The range of terms to be added to Xt_mars_prop		
			//new->k2_new = new->k1_new + new->Knewterm - 1L;
		}	else {
			//The range of terms to be removed from Xt_mars		   
			new->newcols.k1     = KS_old[(newIdx)-1];
			new->newcols.k2_old = KE_old[(endIdx)-1];
			//The range of terms to be added to Xt_mars_prop		
			//new->k2_new = new->k1_new + new->Knewterm - 1L;
		}
		
	}
	else { // if (flag == ChORDER)
	
	   // this branch should be never visited by DUMMY bases
	  // if (basis->type ==DUMMYID) r_error("ChORDER should never been chosen for DUMMYID; there must be something wrong!\n");
		
		if (new->newOrder <= new->oldOrder) { // if (newOrder < oldOrder): remove one or two existing terms
			//NEW.numSeg = 0;	NEW.r[0] = -999;
			new->newcols.k2_old  = KE_old[newIdx - 1];
			new->newcols.k1      = basis->type == SEASONID ? /*season*/(new->newcols.k2_old - 1) : /*trend*/ new->newcols.k2_old;
			//new->k2_new = new->k1_old - 1;
			//Knewterm = k2_new - k1_new + 1; //Equal to 0, meaning no new terms will be added
		}
		else {// (NEW. newOrder > oldOrder): Add new terms					
			new->newcols.k2_old = KE_old[newIdx - 1];          //the term immediately preceding the start of the next segment
			new->newcols.k1     = new->newcols.k2_old + 1;     //the term immediately following the end of the selected segment

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

 
static int __OO_NewKnot_BirthMove(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info, I32PTR maxIndex) {
		
	I32 N                   = info->N;
	I32 Npad16              = info->Npad16;
	BEAST2_MODEL_PTR model  = info->model;
	BEAST2_RANDSEEDPTR PRND = info->pRND;

	// TODO: reconstruct the goodVec array on-the-fly, which means the existing goodVec is dicarded.
	I08PTR goodvec = (I08PTR) basis->goodvec; 
	memset(goodvec, 1, N);

 
	// The current basis shoould be the OUTLIER component
	TKNOT_PTR KNOT  = basis->KNOT;
	I32       nKnot = basis->nKnot; 
	// Exclude the current outlier chgt
	for (int i = 0; i < nKnot; ++i) goodvec[KNOT[i] - 1] = 0;
	 

	I32PTR IndicesLargeDeviation = info->mem;   // mem points to Xnewterm, used here a tmp buff
	I32    numLargeDev       = 0;
	F32PTR deviation         = model->deviation;
	F32    threshold         = info->yInfo->q==1
							   ? model->avgDeviation[0]* info->outlierSigFactor
							   : info->outlierSigFactor;

	F32 maxValue = 0;
	I32 maxIdx   =-1;	
	
	for (I32 i = 0; i < N; i++) {
		F32 value = deviation[i];  // value has already been made positive in CalAbsDev
		if ( !goodvec[i] || IsNaN (value) ) {
			continue;
		}
				
		// Find ther maximum absolute deviation 		
		// If no points have DEV larger than the threshold,  then choose the point with the largest DEV	
		if ( value > maxValue) {						
			maxValue = value;
			maxIdx   = i;
		}

		// If large than the threshold or already being a TCP or SCP		
		// This is a branchless assignment
		IndicesLargeDeviation[numLargeDev] = i;
		numLargeDev += (value > threshold);
	}

	int newKnot = -1L;
	if (numLargeDev > 1) {	
		I32 rndIdx = RANDINT(1, (U16)numLargeDev, *(PRND->rnd16)++);
		newKnot = IndicesLargeDeviation[rndIdx - 1];
	} else if (numLargeDev == 1) {
		newKnot = IndicesLargeDeviation[0];
	}

	if (maxIdx < 0) {
		r_printf("ERROR: __OO_NewKnot_BirthMove: maxIdx=-1, and there must be something wrong!");
	}

	*maxIndex = maxIdx+1L;
	return newKnot+1L;

}

static int __OO_NewIdx_MoveDeath(BEAST2_BASIS_PTR basis, PROP_DATA_PTR info) {
		 
	F32PTR deviation = info->model->deviation; // deviation has already been postive 
	F32    minValue  = 1e34;
	I32    minIdx    = -1;

	// Choose the ocp with the leasst deviation 
	I32       nKnot = basis->nKnot;
	TKNOT_PTR KNOT  = basis->KNOT;
	for (int i = 0; i < nKnot; i++) {
		F32 value = deviation[KNOT[i] - 1]; 	 
		if (value < minValue ) {
			minValue = value;
			minIdx   = i;
		}
	}

	if (minIdx < 0) {
		r_printf("__OO_NewIdx_MoveDeath: maxIdx=-1, and there must be something wrong!");
	}

	return minIdx + 1;
}

static void OO_Propose_01(	BEAST2_BASIS_PTR basis, NEWTERM_PTR new,  PROP_DATA_PTR info)
{	
	I32  goodNum = basis->goodNum; // not used in this function
	I16  nKnot   = basis->nKnot;

	BEAST2_RANDSEEDPTR PRND = info->pRND;
	/*********************************************************************/
	//  DETERMINE WHICH PROPOSAL STEPTO MOVE: MAKE A NEW PROPOSAL		
	/**********************************************************************/
	MOVETYPE flag;
	{
		I32  Ktotal                   = info->model->curr.K;
		I32  MINKNOTNUM               = basis->prior.minKnotNum;
		I32  MAXKNOTNUM               = basis->prior.maxKnotNum;
		I32  MAX_K_StopAddingNewTerms = basis->mcmc_Kstopping;

		if (MINKNOTNUM != MAXKNOTNUM) {
			U08  rnd                      = *(PRND->rnd08)++;
			if (rnd < basis->propprob.birth) { //birth					
				flag = BIRTH;
				if (nKnot >= MAXKNOTNUM )	           flag = MOVE;			
				if (Ktotal > MAX_K_StopAddingNewTerms) flag = (nKnot == 0) ? BIRTH : MOVE;			
			} else if (rnd < basis->propprob.move)   //move							
				flag = nKnot == 0 ? BIRTH : MOVE;		
			else //if (unifRnd <= basis->propprob.death)  //death
				flag = nKnot == 0 ? BIRTH : DEATH;
		} else {
			//MINKNOTNUM == MAXKNOTNUM; they must be larger than zero
			flag = MOVE;
		}
 
	}  // DONE WITH MAKING A NEW PROPOSAL	 

	////////////////////////////////////////////////////////////////////
	I32  samples = info->samples[0];	
	if (samples >= info->nSample_DeviationNeedUpdate) { // this should be always the case bcz OID is picked up only if samples >0
	 	BEAST2_MODEL_PTR model = info->model;
		// Deviation is needed for sampling Outlier knots but not needed for ST. For ST, only extremPos is needed
		__CalcAbsDeviation(model->deviation, model->avgDeviation, info, info->numBasisWithoutOutlier);
		I32  extraSamples = min((10 + samples / 8), 200);
		info->nSample_DeviationNeedUpdate = samples + extraSamples;
		info->shallUpdateExtremVec        = 1L;
	}
	////////////////////////////////////////////////////////////////////

	/*****************************************************************************************/
	// IMPLEMENT THE NEW PROPOSED BASIS
	//*****************************************************************************************/
	TKNOT_PTR  knotList = basis->KNOT;
		
	switch (flag)
	{
	case BIRTH:
	{		
		I32 maxIdx;
		new->newKnot   = __OO_NewKnot_BirthMove(basis, info,&maxIdx);

		if (new->newKnot == 0 && nKnot==0) {
			new->newKnot = maxIdx;
		}

		if (new->newKnot > 0) {
			new->numSeg    = 1;
			new->SEG[0].R1 = new->newKnot;
			new->SEG[0].R2 = new->newKnot;
			new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

			//new->orders2[1] = new->orders2[0] = orderList[(newIdx)-1];

			//implicit conversion from 'int' to 'I16' (aka 'short') changes value from -9999999 to 27009 [-Wconstant-conversion]
			//new->newIdx       = -9999999;            // not used at alll
			new->newIdx = -9999;                       // not used at alll
			new->nKnot_new = nKnot + 1;
		}	else {
			flag = DEATH; 			// Death			  
		}
		
		break;
	}

	case MOVE: 
	{
		//newIdx = RANDINT(1, nKnot, *(PRND->rnd16)++);  //(*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) *tKnotNum);	

		// Randomly choose a good loc: newKnot is the newly chosen breakpoint		
		// I32 randLoc = RANDINT(1, (I32)goodNum, *(PRND->rnd32)++); // will fail if goodNum=0
		// new->newKnot = i08_find_nth_onebyte_binvec(basis->goodvec, (I32)Npad16, randLoc);

		I32 maxIdx;
		new->newKnot   = __OO_NewKnot_BirthMove(basis, info,&maxIdx);
		if (new->newKnot > 0) {
			new->numSeg = 1;
			new->SEG[0].R1 = new->newKnot;
			new->SEG[0].R2 = new->newKnot;
			new->SEG[0].outlierKnot = new->newKnot; //knotBirth is only used for the outlier compnent

			I32 newIdx = __OO_NewIdx_MoveDeath(basis, info);
			//endIdx = newIdx + 1;
			//endIdx = newIdx + 1L;
			new->newIdx    = newIdx;
			new->nKnot_new = nKnot;
			
		}	else {
			flag = DEATH;
		}
		break;
	}//flag=move MOVE MOVE

	}

 	if (flag== DEATH)
	{
		//newIdx = RANDINT(1, nKnot, *(PRND->rnd16)++);       // (*rnd32++) % tKnotNum + 1; // (int)ceilfunc((*rnd32++) * tKnotNum);
		I32 newIdx   = __OO_NewIdx_MoveDeath(basis, info);
		new->newKnot = knotList[newIdx - 1];
		new->numSeg  = 0;

		//new->r1[0] = knotList[(newIdx - 1) - 1];
		//new->r2[0] = knotList[(newIdx + 1) - 1]  -1L;
		//new->orders2[0] = orderList[(newIdx + 1L) - 1]; // ------1|-----------2(newIdx:deleted)|----order kept----(3)|-----

		//endIdx = newIdx;
		//endIdx = newIdx + 1L;
		new->newIdx    = newIdx;
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
		new->newcols.k2_old = KE_old[nKnot - 1];
		new->newcols.k1     = new->newcols.k2_old + 1;
	} else if (flag == DEATH) {
		new->newcols.k2_old  = KE_old[new->newIdx - 1];
		new->newcols.k1      = KS_old[new->newIdx - 1];
	} else if (flag == MOVE) {
		new->newcols.k2_old  = KE_old[new->newIdx - 1];
		new->newcols.k1      = KS_old[new->newIdx - 1];
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