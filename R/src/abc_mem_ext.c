#include <stdlib.h>
#include <string.h>

#include "abc_000_warning.h"

#include "abc_mem_ext.h"
#include "abc_001_config.h"
#include "abc_ide_util.h" //r_printf
 

 

//https: //stackoverflow.com/questions/5061392/aligned-memory-management
static VOID_PTR  malloc_64(size_t N)	{
	VOID_PTR  mem = malloc(N + 64);
	VOID_PTR  ptr = (VOID_PTR )(((uintptr_t)mem + 64) & ~(uintptr_t)0x3F);
	*((char *)((char*)ptr - 1)) = (char)((char *)ptr - (char *)mem);
	return ptr;
}

static  void  free_64(VOID_PTR  p)	{
	char * _restrict porig = (char *)p;
	porig = porig - *(porig - 1);
	free(porig);
}

void mem_init_x(xMemPointers* _restrict self, int printInfo)
{
	self->alloc = mem_alloc_x;
	self->free_all = mem_free_all_x,
	self->maxNumOfPointers = 250,
	self->memNames = NULL,
	self->printInfo = printInfo,

	self->memPointer = (VOID_PTR *)malloc(sizeof(VOID_PTR )*	self->maxNumOfPointers);  //mem64Aligned[MAX_NUM_OF_POINTERS];
	self->mem64Aligned = (I08PTR)malloc(sizeof(I08)*self->maxNumOfPointers); //mem64Aligned[MAX_NUM_OF_POINTERS];
	self->memNames = (char **)malloc(sizeof(char *)*self->maxNumOfPointers);
	self->memNum = 0;
}
VOID_PTR  mem_alloc_x(xMemPointers * _restrict self, I64 sizeInByte, U08 alignment, char * name)
{
	VOID_PTR  newPointer;

	if (alignment == 0)
		newPointer = malloc(sizeInByte);
	else
		newPointer = malloc_64(sizeInByte);

	self->memPointer[self->memNum]  = newPointer;
	self->mem64Aligned[self->memNum] = alignment;
	
	//---------------------------
	self->memNames[self->memNum] = malloc(strlen(name)+1);
	strcpy(self->memNames[self->memNum], name);
	if (self->printInfo )
		// abc_mem_ext.c:42:26: warning: format ‘%x’ expects argument of type ‘unsigned int’, but argument 2 has type ‘VOID_PTR’ 
		//r_printf("%#012x: %d bytes of MEM allocated for '%s' \n", newPointer,sizeInByte, self->memNames[self->memNum]);
		//  abc_mem_ext.c:58:26: warning: '0' flag used with '%p' gnu_printf format [-Wformat=]
		r_printf("%12p: %" PRId64 "bytes of MEM allocated for '%s' \n", newPointer, sizeInByte, self->memNames[self->memNum]);
	 
	self->memNum++;
	return newPointer;
}

void mem_free_all_x(xMemPointers * _restrict self)
{
	for (int i = 0; i < self->memNum; i++)
	{
		if (self->mem64Aligned[i] == 0)
			free(self->memPointer[i]);
		else
			free_64(self->memPointer[i]);
				
		if (self->printInfo)
			//r_printf("%p: Memory de-allocated for %s  \n", self->memPointer[i], self->memNames[i]);
			r_printf("%12p: Memory de-allocated for '%s' \n", self->memPointer[i], self->memNames[i]);

		free(self->memNames[i]);
	}
	if (self->memPointer != NULL)
	{
		free(self->memPointer);
		self->memPointer = NULL;
	}
	if (self->mem64Aligned != NULL)
	{
		free(self->mem64Aligned);
		self->mem64Aligned = NULL;
	}
	if (self->memNames != NULL)
	{
		free(self->memNames);
		self->memNames = NULL;
	}
}



#include "abc_000_warning.h"