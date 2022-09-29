#pragma once

#include <inttypes.h>
#include "abc_000_macro.h"
#include "abc_datatype.h"

typedef struct {
	void** addr;
	int    size;
	int    align;
} MemNode;

typedef struct MemPointers MemPointers;
struct MemPointers
{
	I64        bytesAllocated;
	VOID_PTR * memPointer;
	I08PTR     memAlignOffset;
	I16        npts;
	I16        nptsMax;
	U64      * memHeaderBackup;
	I08        checkHeader;

// https: //stackoverflow.com/questions/1403890/how-do-you-implement-a-class-in-c
//Object Oriented Programming in ANSI-C : http: //www.planetpdf.com/codecuts/pdfs/ooc.pdf
	void       (*init)(     MemPointers *  self);
	VOID_PTR   (*alloc)(    MemPointers *  self, I64 size, U08 alignment);
	VOID_PTR   (*alloc0)(   MemPointers *  self, I64 size, U08 alignment);
	void       (*alloclist) (MemPointers* self, MemNode* list, int aggregatedAllocation, VOIDPTR* nodesRemove);
	void       (*free_all)( MemPointers *  self);
	I32        (*verify_header)(MemPointers* self);
};
 extern void  mem_init(MemPointers* _restrict self);

#define MyALLOC(MEM,  numElem,type,alignment) (type *)(MEM).alloc(&(MEM), (I64) sizeof(type)* (I64)(numElem),alignment)
#define MyALLOC0(MEM, numElem,type,alignment) (type *)(MEM).alloc0(&(MEM), (I64) sizeof(type)* (I64)(numElem),alignment)
//#endif