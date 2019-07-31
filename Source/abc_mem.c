#include "abc_001_config.h"
#include "abc_mem.h"
#include <stdlib.h>
#if M_INTERFACE==1
	#define  r_malloc(x)	mxMalloc(x) 
	#define  r_free(x)		mxFree(x)
#elif R_INTERFACE==1
	#define  r_malloc(x) Calloc(x,char)  
	#define  r_free(x)   Free(x) 
#endif
 static void* malloc_64(size_t N)	{
	void * mem=malloc(N+64);
	void * ptr=(void *)(((uintptr_t)mem+64) & ~(uintptr_t)0x3F);
	*((char *)((char*)ptr - 1))=(char)((char *)ptr - (char *)mem);
	return ptr;
}
 static void  free_64(void * _restrict p)	{
	char * _restrict porig=(char *)p;
	porig=porig - *(porig - 1);
	free(porig);
}
static void * mem_alloc(MemPointers * _restrict self,int64_t sizeInByte,uint8_t alignment)
{
	void * newPointer;
	if (alignment==0)
		newPointer=malloc(sizeInByte);
	else
		newPointer=malloc_64(sizeInByte);
	self->memPointer[self->memNum]=newPointer;
	self->mem64Aligned[self->memNum]=alignment;
	self->memNum++;
	return newPointer;
}
static void mem_free_all(MemPointers * _restrict self)
{
	for (int i=0; i < self->memNum; i++)
	{
		if (self->mem64Aligned[i]==0)
			free(self->memPointer[i]);
		else
			free_64(self->memPointer[i]);
	}
	if (self->memPointer !=NULL)
	{
		free(self->memPointer);
		self->memPointer=NULL;
	}
	if (self->mem64Aligned !=NULL)
	{
		free(self->mem64Aligned);
		self->mem64Aligned=NULL;
	}
}
void mem_init(MemPointers* _restrict self)
{
	*self=(MemPointers) {
			.alloc=mem_alloc,
			.init=mem_init,
			.free_all=mem_free_all,
			.MaxNumOfPointers=250
			};	
	self->memPointer=(void **)malloc(sizeof(void *)* self->MaxNumOfPointers);
	self->mem64Aligned=(int8_t *)malloc(sizeof(int8_t) * self->MaxNumOfPointers);
	self->memNum=0;
}
