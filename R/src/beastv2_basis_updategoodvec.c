#include "abc_000_warning.h"  //put before abc_blas_lapack_lib bcz teh header contains static code
#include "abc_mcmc.h"
#include "abc_vec.h"               //for i08_sum_binvec only
 
#include "abc_blas_lapack_lib.h"   //r_ippsSet_8u
#include "beastv2_header.h"

static void DSVT_UpdateGoodVecForNewTerm(BEAST2_BASIS_PTR basis, NEWTERM_PTR new, I32 Npad16)
{
	I32     newKnot = new->newKnot;
	I32     newIdx  = new->newIdx;

	U08PTR    goodVec  = basis->goodvec;
	I32       MINSEP   = basis->prior.minSepDist;
	TKNOT_PTR knotList = basis->KNOT;

	MOVETYPE flag = new->jumpType;
	if (flag == BIRTH) 
		r_ippsSet_8u(0, goodVec + (newKnot - MINSEP) - 1, 2 * MINSEP + 1);	
	else if (flag == DEATH) {//death						
		I32 oldKnot = knotList[newIdx - 1];         //If flag=death, oldKnot was not assgined before until now
		I32 r1      = knotList[(newIdx - 1) - 1];
		I32 r2      = knotList[(newIdx + 1) - 1] - 1;

		r_ippsSet_8u(1, goodVec + (oldKnot - MINSEP) - 1,   2 * MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (r1)-1,                   MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (r2 - MINSEP + 1) - 1,    MINSEP);
		//GOOD.tY.num = countOnes(N, GOOD.tY.binvec);//ippsSum_16s32s_Sfs(GOOD.tY.binvec, N, &GOOD.tY.num, 0);
	}
	else if (flag == MOVE) {//move
		I32 oldKnot = knotList[newIdx - 1];
		I32 r1      = knotList[(newIdx - 1) - 1];
		I32 r2      = knotList[(newIdx + 1) - 1] - 1;

		//If flag=move, oldKnot has been already assgined before.
		r_ippsSet_8u(1, goodVec + (oldKnot - MINSEP) - 1, 2 * MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (newKnot - MINSEP) - 1, 2 * MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (r1)-1,                 MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (r2 - MINSEP + 1) - 1,  MINSEP);
	}
	else if (flag == MERGE) {//merge
		I32 r1 = knotList[(newIdx - 1) - 1];
		I32 r2 = knotList[(newIdx + 2) - 1] - 1;

		I32 oldKnot = knotList[newIdx - 1]; //If flag=merge, oldKnot was not assgined before until now
		r_ippsSet_8u(1, goodVec + (oldKnot - MINSEP) - 1, 2 * MINSEP + 1);

		oldKnot = knotList[(newIdx + 1) - 1];
		r_ippsSet_8u(1, goodVec + (oldKnot - MINSEP) - 1, 2 * MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (newKnot - MINSEP) - 1, 2 * MINSEP + 1);

		r_ippsSet_8u(0, goodVec + (r1)-1, MINSEP + 1);
		r_ippsSet_8u(0, goodVec + (r2 - MINSEP + 1) - 1, MINSEP);
	}

	//Compute the total number of good positions in Trend if flag is not change-order
	if (flag != ChORDER) {
		//ippsSum_16s32s_Sfs(GOOD.tY.binvec, N, &GOOD.tY.num, 0);
		//GOOD.tY.num = countOnes(N, GOOD.tY.binvec);
		basis->goodNum = i08_sum_binvec(goodVec, Npad16);
	}

	/*
	Feb 19: this is a tricky bug that  literally took me 70+ hrs to figure out, thanks to Michael Faran at Tel aviv university
	for the resources to identify the bug.
	Even though mempcy can handle overlapping memory if moving data forward (e.g., arr[i]=arr[i+1]), this may not be honored on 
	some systems (e.g., avx512). So it is better to just use memmove.
	https://stackoverflow.com/questions/25629736/memcpy-of-overlapping-buffers
	https://sourceware.org/bugzilla/show_bug.cgi?id=12518

	A comment from Linus: 
	And now applications will randomly do different things depending on the phase of the moon (technically, depending on which CPU 
	they have and what particular version of memcpy() glibc happens to choose).	
	*/

#define InsertNewElem(dst,n,newIdx, newValue, T)            for(I32 i=(n); i>=(newIdx);i--) *((T*)dst+i)=*((T*)dst+i-1); \
                                                                 *((T*)dst+newIdx-1)=newValue;
#define RepeatElem(dst, n,newIdx, T)                        for(I32 i=(n); i>=(newIdx); i--) *((T*)dst+i)=*((T*)dst+i-1); 
#define ReplaceElem(dst,n,newIdx, newValue, T)              *((T*)dst+newIdx-1)=newValue;
#define DeleteElem(dst,n,newIdx, T)                          memmove((T*)dst+newIdx-1,(T*)dst+newIdx+1-1,sizeof(T)*(n-(newIdx)));
#define MergeTwoElemWithNewValue(dst,n,newIdx,newValue, T)   *((T*)dst+newIdx-1)=newValue;\
			                                                      memmove((T*)dst+(newIdx+1) - 1, (T*)dst + (newIdx +2) - 1, sizeof(T)* (n - (newIdx+2L)+1L));
#define DoNothing(dst,n, T) ;


	// DUMMY: not used for the DUMMY basis at all, but kept here to ensure re-using the code for DST
	// this requires mem allocated for ORDER, though the values are garbage
	TORDER_PTR orderList = basis->ORDER; 	
    I32        nKnot   = basis->nKnot;
	 
	switch (flag)
	{
	case BIRTH: {
		//Insert the newly chosen bk into the trend knot list of basis_prop
		//The last knot (i.e., an extra knot given after the end of the TS) is fixed at (N+1).	
		//That is why the length given herei is numKNot+1 rather than nKnot.
		InsertNewElem(knotList, nKnot + 1L, newIdx, newKnot, TKNOT);
		RepeatElem(orderList, nKnot + 1, newIdx, TORDER);  //DUMMY: orderList not used for DUMMY
		break;
	}
	case DEATH: {
		DeleteElem(knotList, nKnot + 1, newIdx, TKNOT);
		DeleteElem(orderList, nKnot + 1, newIdx, TORDER);   //DUMMY: orderList not used for DUMMY
		break;
	}//flag==2/2
	case MERGE: {
		MergeTwoElemWithNewValue(knotList, nKnot + 1L, newIdx, newKnot, TKNOT);
		DeleteElem(orderList, nKnot + 1, newIdx + 1, TORDER);  //DUMMY: orderList not used for DUMMY
		break;
	}//flag==merge
	case MOVE: //MOVE
	{
		ReplaceElem(knotList, nKnot + 1, newIdx, newKnot, TKNOT)
		DoNothing(orderList, nKnot + 1, TORDER);        //DUMMY: orderList not used for DUMMY
		break;
	}//flag=move MOVE MOVE
	case ChORDER:
	{
		DoNothing(knotList, nKnot + 1, TKNOT);
		ReplaceElem(orderList, nKnot + 1, newIdx, new->newOrder, TORDER);  //DUMMY: orderList not used for DUMMY
		break;
	}//flag==4//f

	}
	basis->nKnot = new->nKnot_new;
}

static void OO_UpdateGoodVecForNewTerm(BEAST2_BASIS_PTR basis, NEWTERM_PTR new, I32 Npad16_not_used)
{
	I32     newKnot = new->newKnot;
	U08PTR  goodVec = basis->goodvec;

	MOVETYPE flag = new->jumpType;
	if (flag == BIRTH) {
		goodVec[newKnot - 1] = 0;
		basis->goodNum--;
	}
	else if (flag == DEATH) {//death	
		goodVec[newKnot - 1] = 1;
		basis->goodNum++;
	}
	else if (flag == MOVE) {//move
		I32 newIdx  = new->newIdx;
		I32 oldKnot = basis->KNOT[newIdx - 1];
		goodVec[oldKnot - 1] = 1;
		goodVec[newKnot - 1] = 0;
	}

	TKNOT_PTR knotList = basis->KNOT;
	I32       nKnot    = basis->nKnot;
	switch (flag)
	{
	case BIRTH: {
		//Insert the newly chosen bk into the trend knot list of basis_prop
		//The last knot (i.e., an extra knot given after the end of the TS) is fixed at (N+1).	
		//That is why the length given herei is numKNot+1 rather than nKnot.
		knotList[nKnot + 1 - 1] = newKnot;
		break;
	}
	case DEATH: {
		I32 newIdx = new->newIdx;
		DeleteElem(knotList, nKnot, newIdx, TKNOT);
		break;
	}//flag==2/2 
	case MOVE: //MOVE
	{
		I32 newIdx = new->newIdx;
		ReplaceElem(knotList, nKnot, newIdx, newKnot, TKNOT)
			break;
	}//flag=move MOVE MOVE

	}//switch (flag)

	basis->nKnot = new->nKnot_new;
}

void* Get_UpdateGoodVec_KnotList(I08 id) {
	if      (id == SEASONID)    return DSVT_UpdateGoodVecForNewTerm;
	else if (id == SVDID)       return DSVT_UpdateGoodVecForNewTerm;
	else if (id == DUMMYID)     return DSVT_UpdateGoodVecForNewTerm;
	else if (id == TRENDID)     return DSVT_UpdateGoodVecForNewTerm;
	else if (id == DUMMYID)     return DSVT_UpdateGoodVecForNewTerm;
	else if (id == OUTLIERID)   return OO_UpdateGoodVecForNewTerm;

	return NULL;
}

#include "abc_000_warning.h"