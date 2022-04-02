#include <stdlib.h>
#include <string.h>  // for memcpy

#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_001_config.h"
#include "abc_mem.h"

//https://stackoverflow.com/questions/227897/how-to-allocate-aligned-memory-only-using-the-standard-library
//stackoverflow.com/questions/5061392/aligned-memory-management

/************************************************************************/
// malloc_64 and free_63 use the preceding byte to save potential offset.
// We won't use them here bcz of the wate of mem: we have to allocate
// N+64
/************************************************************************/
 static VOID_PTR  malloc_64(size_t N)	{ 
	VOID_PTR  mem = malloc(N + 64);
	VOID_PTR  ptr = (VOID_PTR )(((uintptr_t)mem + 64) & ~(uintptr_t)0x3F);
	*((char *)((char*)ptr - 1)) = (char)((char *)ptr - (char *)mem);
	return ptr;
}
 static void  free_64(VOID_PTR  p)	{
	char * porig = (char*)p - *((char*)p - 1);
	free(porig);
}
 /************************************************************************/
 // malloc_64 and free_63 use the preceding byte to save potential offset.
 // We won't use them here bcz of the wate of mem: we have to allocate
 // N+64
 /************************************************************************/


 static void  ExpandInternelBuf(MemPointers* self ) {
 
	 if (self->npts >= self->nptsMax) {
		 int       oldMax      = self->nptsMax;
		 VOID_PTR* oldPointers = self->memPointer;
		 I08PTR    oldAlign    = self->memAlignOffset;
		 
		 self->nptsMax       = oldMax + 200;
		 self->memPointer    = (VOID_PTR*)malloc(sizeof(VOID_PTR) * self->nptsMax);
		 self->memAlignOffset = (I08PTR)  malloc(sizeof(I08) *      self->nptsMax);

		 if (oldPointers) {
			 memcpy((const void* )self->memPointer,  (const void*)oldPointers, sizeof(VOID_PTR) * oldMax);
			 memcpy((const void*)self->memAlignOffset, (const void*)oldAlign,    sizeof(I08) * oldMax);
			 free(( void*)oldPointers);
			 free(( void*)oldAlign);
		 } 

		 if (self->checkHeader) {	 
			 U64PTR   oldHeader     = self->memHeaderBackup; 
			 self->memHeaderBackup = (U64PTR)malloc(sizeof(64) * self->nptsMax);

			 if (oldHeader) {
				 memcpy((const void*)self->memHeaderBackup, (const void*)oldHeader, sizeof(64) * oldMax);				 
				 free((void*)oldHeader);
			 }
		 }
	 }
 }

static VOID_PTR  MemAlloc(MemPointers * self, I64 N, U08 alignment)
{

	ExpandInternelBuf(self);

	// 0 means arbitray locations, which corresponds to an alignment of 1L
	alignment = alignment == 0 ? 1 : alignment;

	int       isSuccess = 0;
	
	VOID_PTR  ptr  = NULL;
	VOID_PTR  ptrAligned;
	// For small-byte alignent, we try to directy allocate it and try our luck for aligment
	if (alignment <= 8) {
		
		ptr       = malloc(N);
		// 64: (64-1)=0x3F=0b11 1111 ~(64-1)=11111111 00 0000
		ptrAligned = (VOID_PTR)((uintptr_t)ptr & ~(uintptr_t) (alignment-1) );  
		isSuccess = (ptr == ptrAligned);
		self->bytesAllocated += isSuccess?N:0;
	}

	// the allocation above for align<8 failled or if alignment is > 8
	if (!isSuccess) { 
		// Deallocate the preivously allocated ptr
		if (ptr) free(ptr);

		// No need to align N+64 because the worst scenario is that pts is just one byte off from the
		// the alignment boundary
		ptr        = malloc(N  + (alignment-1)  );
		// 64: (64-1)=0x3F=0b11 1111 ~(64-1)=11111111 00 0000
		ptrAligned = (VOID_PTR)( ( (uintptr_t)ptr+alignment-1) &  ~(uintptr_t)(alignment - 1) );		
		self->bytesAllocated += N + (alignment - 1);
	}

	
	self->memPointer[   self->npts]    = ptrAligned;
	self->memAlignOffset[self->npts]   = (uintptr_t)ptrAligned- (uintptr_t)ptr;
	

	if (self->checkHeader) {
		self->memHeaderBackup[self->npts] = *(U64PTR) ((uintptr_t)ptr - 8); 
	}

    self->npts++;
	return ptrAligned;
}

static VOID_PTR  MemAlloc0(MemPointers* _restrict self, I64 sizeInByte, U08 alignment){
	VOID_PTR  ptr = MemAlloc(self, sizeInByte, alignment);
	memset(ptr, 0, sizeInByte);
	return ptr;
}
static void mem_free_all(MemPointers * _restrict self)
{
	for (int i = 0; i < self->npts; i++) 	{
		free(  (char*) self->memPointer[i] - self->memAlignOffset[i]);	
	}

	if (self->memPointer) 	{
		free(self->memPointer);
		self->memPointer = NULL;
	}
	if (self->memAlignOffset) 	{
		free(self->memAlignOffset);
		self->memAlignOffset = NULL;
	}
	if (self->memHeaderBackup) {
		free(self->memHeaderBackup);
		self->memHeaderBackup = NULL;
	}
	self->bytesAllocated = 0;
	self->npts           =0;
	self->nptsMax        = 0;
}
static I32  verify_header(MemPointers* _restrict self) {
	if (!self->checkHeader || self->npts==0) {
		return 0;
	}

	int badHeaderNum = 0;
	for (int i = 0; i < self->npts; ++i) {
		U64 curheader = *(U64PTR)((uintptr_t)self->memPointer[i] - self->memAlignOffset[i] - 8);
		if (curheader != self->memHeaderBackup[i]) {
			badHeaderNum++;
		}
	}
	return badHeaderNum;
}
void mem_init(MemPointers* self) {	 
	*self = (MemPointers) {
			.alloc		        = MemAlloc,
			.alloc0             = MemAlloc0,
			.init		        = mem_init,
			.free_all	        = mem_free_all,
			.nptsMax            =0,
			.npts               =0,
			.bytesAllocated		=0,
			.memAlignOffset     =NULL,
			.memPointer         = NULL,
			.memHeaderBackup    =NULL,
			.checkHeader        =0,
			.verify_header   = verify_header,
			};			
}

#include "abc_000_warning.h"
