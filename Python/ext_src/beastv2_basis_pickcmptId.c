#include "abc_000_warning.h"
#include "beastv2_header.h"
 
//none: fix none
//fix1: fix the 1st cmpnt

// If all the non-outlier bases are fixed, we will already return the 1st basis
// even if both season and trend are present. This won't matter because no terms
// will be changed at all during the proposal MOVE.
static I32 _PickBaisID_1none_2fix1_2fix01(PROP_DATA_PTR PROPINFO) {
	// for 2fix01: both bases are fixed; then the basis will be a NoChangeFixGloBAL
	return 0;
}

static I32 _PickBaisID_______2fix0______(PROP_DATA_PTR PROPINFO) { return 1;}
static I32 _PickBaisID_______2none_3fix2(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM *PRAND = PROPINFO->pRND;
	return (*PRAND->rnd08++ > 128);
}
static I32 _PickBaisID_____________3none(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM*	PRAND	= PROPINFO->pRND;
	U08					unifRnd = *PRAND->rnd08++;
	
	/*if		(unifRnd <  (U08)(255*0.33) )		return 0;   //R1=1 R2=1 2
	else if (unifRnd < (U08)(255*0.66) )		return 1;   //R1=0 R2=1 1
	else										return 2;   //R1=0 R2=0 0 
	*/
	/*****************/
	//Remove the if branches
	/*****************/
	I32 R1 = unifRnd < (U08)(255 * 0.33);
	I32 R2 = unifRnd < (U08)(255 * 0.66);
	return 2L-(R1 + R2);
}
static I32 _PickBaisID_____________3fix1(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM* PRAND  = PROPINFO->pRND;
	U08				  unifRnd = *PRAND->rnd08++;
	/*	if		(unifRnd <  (U08)(255*0.33) )		return 0;   //R1=1 R2=1 2
		else if (unifRnd < (U08)(255*0.66) )		return 1;   //R1=0 R2=1 1
		else										return 2;   //R1=0 R2=0 0	
	*/
	/*****************/
	//Remove the if branches
	/*****************/ 
	return unifRnd < 128? 0:2; //gives only 0 or 2
}


//THe outlier is always the last component
static I32 _PickBaisID_hasOutlier_2fix0(PROP_DATA_PTR PROPINFO) {
	// Do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) {
		//this basis will be a NoChangeFIXGlobal
		return 0; 
	}
	return 1;
}
static I32 _PickBaisID_hasOutlier_2none(PROP_DATA_PTR PROPINFO) {
	// Do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) {return 0;}
	
    BEAST2_RNDSTREAM* PRAND = PROPINFO->pRND;
	return (*PRAND->rnd08++ > 128)?0:1;
}
static I32 _PickBaisID_hasOutlier_3none(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM*	PRAND	= PROPINFO->pRND;
	U08					unifRnd = *PRAND->rnd08++;

	// do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) {	return (unifRnd > 128);	}

	if		(unifRnd <  (U08)(255*0.33) )		return 0;
	else if (unifRnd < (U08)(255*0.66) )		return 1;
	else										return 2; 
}
static I32 _PickBaisID_hasOutlier_3fix0(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM* PRAND   = PROPINFO->pRND;
	U08				  unifRnd = *PRAND->rnd08++;

	// do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) {return 1;}
	return (unifRnd > 128) ? 1 : 2;

}
static I32 _PickBaisID_hasOutlier_3fix1(PROP_DATA_PTR PROPINFO) {
	BEAST2_RNDSTREAM* PRAND   = PROPINFO->pRND;
	U08				  unifRnd = *PRAND->rnd08++;

	// do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) {	return 0;	}
	return (unifRnd > 128) ? 0 : 2;
}
static I32 _PickBaisID_hasOutlier_3fix01(PROP_DATA_PTR PROPINFO) {
	// Do not pick oid if sample == 0; pick oid only if sammple > 0
	if (*PROPINFO->samples == 0) { 
		//this basis will be a NoChangeFIXGlobal
		return 0; 
	}
	return 2;
}
void * Get_PickBasisID(I08 numBasis, I08 hasOutlier, I32PTR isComponentFixed)
{
	if		(numBasis == 1) 
		return _PickBaisID_1none_2fix1_2fix01;  //NoChangeFIXGLOBAL is possible to be generated, if the basis is fixed
	else if (numBasis == 2) {
		if (!hasOutlier) {
			// Four possibiites 
			if (isComponentFixed[0] && isComponentFixed[1]) //both are fixed: NoChangeFIXGLOBAL will be generated
				return _PickBaisID_1none_2fix1_2fix01;			 
			else if (isComponentFixed[0])                   //1st is fixed
				return _PickBaisID_______2fix0______;
			else if (isComponentFixed[1])				    //2nd is fixed
				return _PickBaisID_1none_2fix1_2fix01;
			else                                            // neither is fixed
				return _PickBaisID_______2none_3fix2;
		} else { 
			//hasOulier==1
			if (isComponentFixed[0])
				return _PickBaisID_hasOutlier_2fix0;		//NoChangeFIXGLOBAL will be generated
			else //the second cmpt must be outliers
				return _PickBaisID_hasOutlier_2none;
		}
	}
	else if (numBasis == 3) {
		if (!hasOutlier) {
			// If numbasis==3, the third compnt must be the oultiier
		}
		else { //hasOulier==1
			// Four possibilities 
			if   (isComponentFixed[0] && isComponentFixed[1])  //both bases are FIXED: NoChangeFIXGLOBAL will be generated
				return _PickBaisID_hasOutlier_3fix01;
			else if  (isComponentFixed[0])                     //1st is fixed
				return _PickBaisID_hasOutlier_3fix0;
			else if (isComponentFixed[1])                      //2st is fixed
				return _PickBaisID_hasOutlier_3fix1;
			else                                               //none is fixed
				return _PickBaisID_hasOutlier_3none;
		}

	}
 
	return NULL;
}

#include "abc_000_warning.h"