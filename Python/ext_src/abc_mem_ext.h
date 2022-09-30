//#ifndef  BEAST_MEM_HEADER
//#define def  BEAST_MEM_HEADER
#pragma once
#include "abc_datatype.h"
//https: //stackoverflow.com/questions/5061392/aligned-memory-management

typedef struct xMemPointers xMemPointers;
struct xMemPointers
{
	VOID_PTR * memPointer;
	I08PTR  mem64Aligned;
	char    **memNames;
	I16   memNum;
	I16   maxNumOfPointers;
	I16   printInfo;


// https: //stackoverflow.com/questions/1403890/how-do-you-implement-a-class-in-c
//Object Oriented Programming in ANSI-C : http: //www.planetpdf.com/codecuts/pdfs/ooc.pdf
	void	(*init)(xMemPointers*self, int printInfo);
	VOID_PTR  (*alloc)(xMemPointers*self, I64 size, U08 alignment, char * name);
	void   (*free_all)(xMemPointers*self);

};
extern void mem_init_x(xMemPointers* _restrict self, int printInfo);
extern VOID_PTR mem_alloc_x(xMemPointers * _restrict self, I64 sizeInByte, U08 alignment, char * name);
extern void  mem_free_all_x(xMemPointers * _restrict self);
//#endif
//#define MEM_alloc(name,type, MEM,numElement,align) name##=( type *)##(MEM).alloc((&##(MEM)), sizeof(type)*(numElement), align,#name )
#define MEM_alloc(name,type, MEM,numElement,align) name=(type *)(MEM).alloc(&(MEM), sizeof(type)*(numElement), align,#name )