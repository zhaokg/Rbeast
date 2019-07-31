#pragma once
#include <inttypes.h>
#include "abc_000_macro.h"
typedef struct MemPointers MemPointers;
struct MemPointers
{
	void   ** memPointer;
	int8_t  *mem64Aligned;
	int16_t  memNum;
	int16_t  maxNumOfPointers;
	void(*init)(MemPointers * _restrict  self);
	void * (*alloc)(MemPointers * _restrict  self,int64_t size,uint8_t alignment);
	void(*free_all)(MemPointers * _restrict  self);
};
 extern void  mem_init(MemPointers* _restrict self);
#define MyALLOC(MEM,bytes,type,alignment) (type *)(MEM).alloc(&(MEM),sizeof(type)*(bytes),alignment)
