#pragma once

#include <stdio.h>
#include "abc_000_macro.h"
#include "abc_001_config.h"
#include "abc_mem.h"
#include "abc_datatype.h"       //#include <inttypes.h>//#include <stdint.h>

char* word_wrap(char* str, int LINE_LENGTH);
char* word_wrap_indented(char* str, int LINE_LENGTH, int nspace);

extern void    ToUpper(char* s);
extern void    QuickSortD(F32PTR arr, I32PTR  index, I32 low, I32 high);
extern void    QuickSortA(F32PTR arr, I32PTR  index, I32 low, I32 high);
extern I32	   strcicmp(char const * _restrict a, char const * _restrict b);
extern I32     strcicmp_nfirst(char const* _restrict a, char const* _restrict b, int nfirst);
extern I32     strmatch(char const* _restrict full, char const* _restrict part);
extern F32     DeterminePeriod(F32PTR Y, I32 N);
extern I32     FindChangepoint(F32PTR prob, F32PTR mem, F32 threshold, I32PTR cpt, F32PTR cptCI, I32 N, I32 minSepDist, I32 maxCptNumber);


extern void SetupRoutines_AutoByCPU(Bool quiet);
extern void SetupRoutines_UserChoice(int avxOption);

//SUM: IppStatus ippsSum_32f(const Ipp64f* pSrc, int len, Ipp64f* pSum); 
//r_ippsSum_32f(ptr, N, &mValue, ippAlgHintAccurate);
//r_ippsMeanStdDev_32f(ptr, N, &mean, &std, ippAlgHintAccurate);
//ippsSubC_32f_I(Ipp64f val, Ipp64f* pSrcDst, int len);
//r_ippsSubC_32f_I(mean, ptr, N);
//IppStatus ippsStdDev_32f(const Ipp64f* pSrc, int len, Ipp64f* pStdDev);
//ippsStdDev_32f(ptr, N, &mValue, ippAlgHintAccurate);
 

#define    CopyInsertNewElement(SRC, NSRC, DST, newPos, newValue, T)                                     \
                         memcpy(DST,                SRC,            (newPos-1)*sizeof(T));        \
	                     *((T*)DST+newPos-1L)  =  newValue;                                             \
                         memcpy((T*)DST+newPos,  (T*)SRC+newPos-1,  (NSRC-(newPos-1))*sizeof(T));   

#define   CopyRepeatChosenElement(SRC, NSRC, DST, newIdx, T)                                                  \
                         memcpy( (T*)DST,               (T*)SRC,           (newIdx)*sizeof(T)          );  \
                         memcpy( (T*)DST+(newIdx+1) -1, (T*)SRC+newIdx-1,  (NSRC-(newIdx-1))*sizeof(T) );

#define   CopyDeleteChosenElement(SRC, NSRC, DST, newIdx,T)                                     \
              memcpy( (T*)DST,           (T*)SRC,                 sizeof(T)*(newIdx -  1      )    );  \
              memcpy( (T*)DST+newIdx-1,  (T*)SRC+(newIdx+1)-1,    sizeof(T)*(NSRC   - (newIdx) )   );  

#define   CopyMergeTwoValuesWithNewValue(SRC, NSRC, DST, newIdx, newKnot, T)                            \
              memcpy( DST,    SRC,    sizeof(T)*(newIdx-1)      );                                      \
              *((T*)DST+newIdx-1) = newKnot;                                                            \
              memcpy( (T*)DST+(newIdx+1)-1,  (T*)SRC+(newIdx+2)-1,    sizeof(T)*(NSRC-(newIdx+1))   );  

#define   CopyReplaceChosenElemWithNewValue( SRC, NSRC, DST, newIdx, newKnot, T )                          \
              memcpy( DST,   SRC,   (NSRC)*sizeof(T) );  \
              *((T*)DST+newIdx-1) = newKnot; 

#define   CopyNoChange( SRC, NSRC, DST, T )       memcpy(DST,SRC, (NSRC)*sizeof(T) ); 
              
