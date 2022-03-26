#include <string.h>
#include "abc_000_warning.h"
#include "beastv2_header.h"
#include "abc_mem.h"
 
 

static void DSVT_AllocInitBasis(BEAST2_BASIS_PTR basis, I32 N, I32 K_MAX, MemPointers* MEM)
{
	//K_MAX is needed only for allocating termsType

	I32 MAX_KNOTNUM = basis->prior.maxKnotNum;
	I32 MAX_NUM_SEG = (MAX_KNOTNUM + 1L);

    // there is one extra chpt fixed at the first element, and another extra fixed at the very end (N+1)
	// For the dummy basis, ORDER is not needed...		
	I32 nBytes		= sizeof(TKNOT) * (1 + MAX_KNOTNUM + 1)/*knot len*/
					+ sizeof(TORDER) * MAX_NUM_SEG		   /*order len*/;

	basis->KNOT    = (TKNOT_PTR) MyALLOC(*MEM, nBytes, char, 64);
	*basis->KNOT++ = 1L;                              //The first changepont is fixed at the ts start, and move the pointer one step foward, so the fixed br has a index of -1.
	basis->ORDER   = (TORDER_PTR)(basis->KNOT + (MAX_KNOTNUM + 1)); //there is one extra chngpt fixed at N+1, so the skip lenghth here is MAX_KNOTNUM+1
    // basis->ORDER is not used at all for DUMMY basis but mem is still allocated here in order to 
	// re-use the same propose() routinue for DST bases. But the values in basis->ORDER are garbages for DUMMY
	
	I32 nElem =  MAX_NUM_SEG * 2;
	// ks and ke are local to the basis and they must be added with Kbase to get the true col idx
	basis->ks = MyALLOC(*MEM, nElem, I16, 64);
	basis->ke = basis->ks + MAX_NUM_SEG;

	/*************************************************************************************************/
	// goodVec is a N-element binary vector indicating the times/points available as potential knots*/
	// its length is extended to a multiple of 16 bytes for use in i08_find_nth_onebyte_binvec; 
	/*************************************************************************************************/
	I32 Npad16 = 16 * ((N + 15) / 16);
	basis->goodvec = MyALLOC(*MEM, Npad16, U08, 8);
	memset(basis->goodvec + N, 0L, Npad16 - N);    //the extra padded bytes are zeroed out

	//TODO: No need for the dummy basis, except for precType=3 (Orderwise, not imeplemented yet)
	basis->termType = MyALLOC(*MEM, K_MAX, U08, 64);// Added for handling term-specific precisions

	/*??????????????????????????????????????*/
	//termType is not used here and there is a problem with it
	//convert_basis_trend/season will update to get the correct termType for a proposed model
	//but convert_basis_trend/season is called only if the proposed model is accepted....
	//so, if future implemenation will consider different prior precisison parameters for different
	//term types. This must be corrected: convert_basis is needed to be called before doing the
	//the Chol decomposition
	/*??????????????????????????????????????*/
	/**************************************/
	//The same issue is with MODEL.k_const: needed for computing marg_lik but is evaluated
	//only after the marg_Lik is finished and the model is accepted
	//: K*a+(Knew-K)*a is not always equal to (K+x)*a +(Knew-(K+x))*a
	/**************************************/
}
static void OO_AllocInitBasis(BEAST2_BASIS_PTR basis, I32 N, I32 K_MAX, MemPointers* MEM)
{
	//K_MAX is needed only for allocating termsType

	I32 MAX_KNOTNUM = basis->prior.maxKnotNum;
	I32 MAX_NUM_SEG = (MAX_KNOTNUM); // For the outlier part, one knot gives one term so that NUM_SEG=KNOT_NUM, which differs from trend and seasonal parts.
	I32 nBytes;


	//nBytes = sizeof(TKNOT) * (1+ MAX_KNOTNUM + 1)/*knot len*/ 
			//+ sizeof(TORDER) * MAX_NUM_SEG /*order len*/;
	//basis->KNOT  = MyALLOC(*MEM, nBytes, char, 64);
	//*basis->KNOT++ = 1L;                            //The first changepont is fixed at the ts start, and move the pointer one step foward, so the fixed br has a index of -1.
	//basis->ORDER = basis->KNOT + (MAX_KNOTNUM + 1); //there is one extra chngpt fixed at N+1, so the skip lenghth here is MAX_KNOTNUM+1

	nBytes       = sizeof(TKNOT) * MAX_KNOTNUM; //*knot len*/			
	basis->KNOT = (TKNOT_PTR) MyALLOC(*MEM, nBytes, char, 64);  // There is no ORDER needed for the outlier component
	//*basis->KNOT++ = 1L;                            //The first changepont is fixed at the ts start, and move the pointer one step foward, so the fixed br has a index of -1.
	//basis->ORDER = basis->KNOT + (MAX_KNOTNUM + 1); //there is one extra chngpt fixed at N+1, so the skip lenghth here is MAX_KNOTNUM+1


	//The outlier cmpnt is optional, so it can be empty. One extra element (i.e., 1+max_num_seg) is allocated
	// for saving the starts (ks[-1]=1, ke[-1]=0), which is needed to insert new terms into an empty
	// outlier component.
	nBytes         = sizeof(I16) * (1L + MAX_NUM_SEG) * 2;
	basis->ks      = (I16PTR) MyALLOC(*MEM, nBytes, U08, 64);
	*basis->ks++   = 1;
	basis->ke      = basis->ks + MAX_NUM_SEG;
	*basis->ke++   = 0; 
	// the current seg width is 0-1+1=0;

	//basis->termType = NULL; //No term types are needed here bcz all terms are the same
	basis->termType = MyALLOC(*MEM, K_MAX, U08, 64);// Added for handling term-specific precisions
	/*************************************************************************************************/
	// goodVec is a N-element binary vector indicating the times/points available as potential knots*/
	// their length is extended to be a multiple of 16 bytes for use in i08_find_nth_onebyte_binvec; the extra padded
	// bytes are zeroed out
	/*************************************************************************************************/
	I32 Npad16 = 16 * ((N + 15) / 16);
	basis->goodvec = MyALLOC(*MEM, Npad16, U08, 64);
	memset(basis->goodvec + N, 0L, Npad16 - N);
}

void* Get_AllocInitBasis(I08 id) {
	switch (id) {
	case SVDID:
	case DUMMYID:	
	case SEASONID:  
	case TRENDID:  
		return DSVT_AllocInitBasis;
		break;
	case OUTLIERID:
		return OO_AllocInitBasis;
	}

	return NULL;
}

#include "abc_000_warning.h"