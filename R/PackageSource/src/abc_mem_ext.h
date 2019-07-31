
#pragma once
#include <inttypes.h>
#include "abc_000_macro.h"
#include "abc_common.h"
typedef struct xMemPointers xMemPointers;
struct xMemPointers
{
	void   ** memPointer;
	int8_t  * mem64Aligned;
	char    **memNames;
	int16_t   memNum;
	int16_t   maxNumOfPointers;
	int16_t   printInfo;
	void	(*init)(xMemPointers*self);
	void * (*alloc)(xMemPointers*self,int64_t size,uint8_t alignment,char * name);
	void   (*free_all)(xMemPointers*self);
};
extern void  mem_init_x(xMemPointers* _restrict self);
extern void *mem_alloc_x(xMemPointers * _restrict self,int64_t sizeInByte,uint8_t alignment,char * name);
extern void  mem_free_all_x(xMemPointers * _restrict self);
#define MEM_alloc(name,type,MEM,numElement,align) name=(type *)(MEM).alloc(&(MEM),sizeof(type)*(numElement),align,#name )
