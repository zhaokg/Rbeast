#include "abc_mem_ext.h"
#include "abc_001_config.h"
#include <stdlib.h>
#include <string.h>
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
static  void  free_64(void * _restrict p)	{
	char * _restrict porig=(char *)p;
	porig=porig - *(porig - 1);
	free(porig);
}
void mem_init_x(xMemPointers* _restrict self)
{
	self->memPointer=(void **)malloc(sizeof(void *)*	self->maxNumOfPointers);  
	self->mem64Aligned=(int8_t *)malloc(sizeof(int8_t)*self->maxNumOfPointers); 
	self->memNames=(char **)malloc(sizeof(char *)*self->maxNumOfPointers);
	self->memNum=0;
}
void * mem_alloc_x(xMemPointers * _restrict self,int64_t sizeInByte,uint8_t alignment,char * name)
{
	void * newPointer;
	if (alignment==0)
		newPointer=malloc(sizeInByte);
	else
		newPointer=malloc_64(sizeInByte);
	self->memPointer[self->memNum]=newPointer;
	self->mem64Aligned[self->memNum]=alignment;
	self->memNames[self->memNum]=malloc(strlen(name)+1);
	strcpy(self->memNames[self->memNum],name);
	if (self->printInfo !=0)
		r_printf("%d bytes of memory allocated for%s  \n",sizeInByte,self->memNames[self->memNum]);
	self->memNum++;
	return newPointer;
}
void mem_free_all_x(xMemPointers * _restrict self)
{
	for (int i=0; i < self->memNum; i++)
	{
		if (self->mem64Aligned[i]==0)
			free(self->memPointer[i]);
		else
			free_64(self->memPointer[i]);
		if (self->printInfo !=0)
			r_printf("Memory de-allocated for%s  \n",self->memNames[i]);
		free(self->memNames[i]);
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
	if (self->memNames !=NULL)
	{
		free(self->memNames);
		self->memNames=NULL;
	}
}
