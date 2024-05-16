#include <stdlib.h>  // malloc
#include <string.h>  // memcpy
#include <stdarg.h>

#include "abc_000_warning.h"

#include "abc_ide_util.h"
#include "abc_001_config.h"
#include "abc_mem.h"

// Two's complement also gives the correct mask.
// #define AlignedPointer(ptr, Alignment)   (VOID_PTR) (  ( (ptr) + Alignment-1) &  (int64_t) -Alignment  )

// EXAMPLE for alginmengt=64: (64-1)=0x3F=0b11 1111 ~(64-1)=11111111 00 0000
// (VOID_PTR)(((uintptr_t)ptr + buf->align - 1) & ~(uintptr_t)(buf->align - 1));
#define AlignedPointer(ptr, Alignment)  (VOID_PTR)  ( ((uintptr_t) (ptr) + Alignment-1) &  ~(uintptr_t)(Alignment-1) )

#define FreeNullify(ptr) {if (ptr!=NULL){ free(ptr); ptr=NULL;} }
// https:// stackoverflow.com/questions/227897/how-to-allocate-aligned-memory-only-using-the-standard-library
//          stackoverflow.com/questions/5061392/aligned-memory-management

/************************************************************************/
// malloc_64 and free_63 use the preceding byte to save potential offset.
// We won't use them here bcz of the wate of mem: we have to allocate
// N+64
/************************************************************************/
 static VOID_PTR  malloc_64(size_t N)	{ 
	VOID_PTR  mem = malloc(N + 64);
	VOID_PTR  ptr  = (VOID_PTR )(((uintptr_t)mem + 64) & ~(uintptr_t)0x3F);
	*((char *)((char*)ptr - 1)) = (char)((char *)ptr - (char *)mem);
	return ptr;
}

 static void  free_64(VOID_PTR  p)	{
	char * porig = (char*)p - *((char*)p - 1);
	free(porig);
}

 typedef struct {
	 VOIDPTR  ptr;
	 VOIDPTR  alignedptr;
	 int      byte_allocated;
 } MemAlignedPtr;

 static MemAlignedPtr malloc_aligned(size_t N, int alignment) {

	 alignment = max(1, alignment);

	 int       isSuccess      = 0;
	 size_t    bytesAllocated = 0;
	 VOID_PTR  ptr       = NULL;
	 VOID_PTR  ptrAligned;
	 // For small-byte alignent, we try to directy allocate it and try our luck for aligment
	 if (alignment <= 8) {
		 ptr            = malloc(N);
		 ptrAligned     = AlignedPointer(ptr, alignment);
		 isSuccess      = (ptr == ptrAligned);
		 bytesAllocated = isSuccess ? N : 0;
	 }

	 // The allocation above for align<8 failled or if alignment is > 8
	 if (!isSuccess) {
		 // Deallocate the preivously allocated ptr
		 if (ptr) free(ptr);

		 // No need to request for a total of N+64 because the worst scenario is that pts is just 
		 // one byte off from the the alignment boundary
		 ptr = malloc(N + (alignment - 1));
		 // 64: (64-1)=0x3F=0b11 1111 ~(64-1)=11111111 00 0000
		 ptrAligned = AlignedPointer(ptr, alignment);
		 bytesAllocated  = N + (alignment - 1);
	 }

	 MemAlignedPtr result;
	 result.ptr        = ptr;
	 result.alignedptr = ptrAligned;
	 result.byte_allocated = (int) bytesAllocated;
	 return result;
 }


 static void  ExpandInternelBuf(MemPointers* self ) {
 
	 if (self->npts >= self->nptsMax) {
		 int       oldMax      = self->nptsMax;
		 VOID_PTR* oldPointers = self->memPointer;
		 I08PTR    oldAlign    = self->memAlignOffset;
		 
		 self->nptsMax       = oldMax + 200;
		 self->memPointer    = (VOID_PTR*)malloc(sizeof(VOID_PTR) * self->nptsMax);
		 self->memAlignOffset = (I08PTR)  malloc(sizeof(I08) *      self->nptsMax);

		 if (oldPointers) {
			 memcpy((  void* )self->memPointer,  (const void*)oldPointers, sizeof(VOID_PTR) * oldMax);
			 memcpy(( void*)self->memAlignOffset, (const void*)oldAlign,    sizeof(I08) * oldMax);
			 free(( void*)oldPointers);
			 free(( void*)oldAlign);
		 } 

		 if (self->checkHeader) {	 
			 U64PTR   oldHeader     = self->memHeaderBackup; 
			 self->memHeaderBackup = (U64PTR)malloc(sizeof(U64) * self->nptsMax);

			 if (oldHeader) {
				 memcpy((  void*)self->memHeaderBackup, (const void*)oldHeader, sizeof(U64) * oldMax);
				 free((void*)oldHeader);
			 }
		 }
	 }
 }

static VOID_PTR  mem_alloc(MemPointers * self, I64 N, U08 alignment) {

	if (N <= 0) {
		return NULL;
	}

	ExpandInternelBuf(self);

	// 0 means arbitray locations, which corresponds to an alignment of 1L
	alignment = max(1, alignment);

	MemAlignedPtr resultPtr = malloc_aligned(N, alignment);	
	self->memPointer[   self->npts]    = resultPtr.alignedptr;
	self->memAlignOffset[self->npts]   = (int8_t) ((uintptr_t)resultPtr.alignedptr - (uintptr_t)resultPtr.ptr);
	self->bytesAllocated              += resultPtr.byte_allocated;

	if (self->checkHeader) {
		self->memHeaderBackup[self->npts] = *(U64PTR) ((uintptr_t)resultPtr.ptr - 8);
	}

    self->npts++;
	return resultPtr.alignedptr;
}

static VOID_PTR  mem_alloc0(MemPointers* _restrict self, I64 sizeInByte, U08 alignment){
	VOID_PTR  ptr = mem_alloc(self, sizeInByte, alignment);
	if (ptr) memset(ptr, 0, sizeInByte);
	return ptr;
}

static void mem_free_all(MemPointers * _restrict self) {

	for (int i = 0; i < self->npts; i++) 	{
		free( (char*) self->memPointer[i] - self->memAlignOffset[i]);	
	}

	FreeNullify(self->memPointer);
	FreeNullify(self->memAlignOffset);
	FreeNullify(self->memHeaderBackup);
 
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
		if ((curheader & 0xffff0000ffffffff) != (self->memHeaderBackup[i] & 0xffff0000ffffffff)) {
			// 0x117e306c323496b3 vs 	0x117e687b323496b3
			// For some reason, the header of >memHeaderBackup[1] is changed slightly 
			// when new blocks are allocated. So instead check the full 8 bytes and we check
			// only partially according to the stable bytes 
			badHeaderNum++;
		}
	}
	return badHeaderNum;
}

static void  memnodes_nullify_badnodes(MemNode* list, VOIDPTR * nodesRemove) {

	if (nodesRemove == NULL) {
		return;
	}

	int i = 0;
	while (nodesRemove[i]) {	        
		int j = 0;
		while (list[j].addr) {
			if ( list[j].addr == nodesRemove[i] ) {
				list[j].size   = 0;
				break;
			}
			j++;
		}
		i++;
	}
}

 I64  memnodes_calc_offsets(MemNode* list, int * maxAlignment) {

	// Updat the offset_from_origin fields of the list
    int maxAlign    = 1;
	I64 curr_offset = 0; // in the loop below, cur_offset is 0 for the first iteration 	i=0
	int i = 0;
	while (list[i].addr) {
		int curalign               = ( list[i].size == 0 || list[i].align <=0) ? 1 : list[i].align; 
		list[i].offset_from_origin = (list[i].size == 0)   ? curr_offset : (int) (intptr_t) AlignedPointer(curr_offset, curalign);
		curr_offset                = list[i].offset_from_origin + list[i].size;
		i++;

		maxAlign = max(maxAlign, curalign);
	}
	if (maxAlignment) maxAlignment[0] = maxAlign;
	I64 totalSize = curr_offset;

	list[0].offset_from_origin   =   i;  // i is the total number of valid notes
	list[i].size                 = (int) totalSize;
	list[i].align                = maxAlign;
	
	return totalSize;
}

 void   memnodes_assign_from_alignedbase(MemNode* list, VOIDPTR pbase) {

	 // Must call memnode_get_offset first to fill list.offset_from_orign

	 int nNodes   = (int) list[0].offset_from_origin;     // the offset of the first node is used to record the total number of nodes
	 int maxAlign = nNodes>0? list[nNodes].align:1L;      // the sentinal node contain the maxalginment and total size
	 if (nNodes <= 0) {
		 memnodes_calc_offsets(list, NULL); // offset_from_orign has been updated inside
		 nNodes   = (int) list[0].offset_from_origin;
		 maxAlign = (int) list[nNodes].align;
	 }

	 if (AlignedPointer(pbase, maxAlign) != pbase) {
		 r_printf("Error: the input base pointer is not aligned!\n");
		 return;
	 }

	 int j = 0;
	 list[0].offset_from_origin = 0;  // TEmporally change it to zero
	 while (list[j].addr) {
		 list[j].addr[0] = (list[j].size == 0) ? NULL : (char*)pbase + list[j].offset_from_origin;
		 j++;
	 }

	 list[0].offset_from_origin = nNodes; // Re-save the number of nodes
 }

 void   memnodes_assign_from_unalignedbase(MemNode* list, VOIDPTR pbase, int bufsize) {
	 // Must call memnode_get_offset first to fill list.offset_from_orign
	 int nNodes    = (int)  list[0].offset_from_origin;    // the offset of the first node is used to record the total number of nodes
	 int totalsize = (nNodes > 0) ? list[nNodes].size : 0; // the sentinal node contain the maxalginment and total size
	 int maxAlign  = (nNodes > 0) ? list[nNodes].align : 1L; 
	 if (nNodes == 0) {
		 memnodes_calc_offsets(list, NULL); // offset_from_orign has been updated inside
		 nNodes    = (int) list[0].offset_from_origin;
		 totalsize = (int) list[nNodes].size;
		 maxAlign  = (int) list[nNodes].align;
	 }

	 char* palgiend = AlignedPointer(pbase, maxAlign);
	 int offset     = (int) (palgiend-(char*)pbase);
	 
	if (offset + totalsize > bufsize) {
		 r_printf("Error: the buf has no enough space!\n");
		 return;
	 }
 
	 int j = 0;
	 list[0].offset_from_origin = 0;  // TEmporally change it to zero
	 while (list[j].addr) {
		 list[j].addr[0] = (list[j].size == 0) ? NULL : palgiend + list[j].offset_from_origin;
		 j++;
	 }

	 list[0].offset_from_origin = nNodes; // Re-save the number of nodes
 }

 I64    memnodes_find_max_common_size(MemNode* list, ...) {
// https:// www.cprogramming.com/tutorial/c/lesson17.html
	 va_list arguments;
	 va_start(arguments, list); 	 /* Initializing arguments to store all values after num */
	 
	 int nlist = 0;
	 MemNode * listVec[100];	 
	 listVec[nlist++] = list;

	  MemNode* nextList;
	 do {
		 nextList         = va_arg(arguments, MemNode*);
		 listVec[nlist++] = nextList;
	 } while (nextList !=NULL);
	 va_end(arguments);             
	 
	 nlist--;
	 // nlist has to be at least one

	 I64 maxTotalSize = 0;
	 for (int i = 0; i < nlist; i++) {
		 MemNode* LIST= listVec[i];
		 int nNodes = (int) LIST[0].offset_from_origin;
		 if (nNodes <= 0) {
			 memnodes_calc_offsets(LIST, NULL);
			 nNodes = (int) LIST[0].offset_from_origin;
		 }
		 int totalSize = LIST[nNodes].size;
		 int maxAlgin  = LIST[nNodes].align;

		 maxTotalSize = max(maxTotalSize, totalSize + maxAlgin - 1);
	 }

	 return maxTotalSize;
 }

static void  mem_alloc_list(MemPointers* self, MemNode* list, int aggregatedAllocation, VOIDPTR * nodesRemove) {
 
	memnodes_nullify_badnodes(list, nodesRemove);

	/*************************************/
	// AggregatedAllocation == FALSE
	/*************************************/ 
	if (!aggregatedAllocation) {
		int i = 0;
		while (list[i].addr ) {									 
			*list[i].addr = mem_alloc(self, list[i].size, list[i].align);
			i++;
		}    
		return;
	}

	/*************************************/
	// AggregatedAllocation == TRUE
	/*************************************/ 
	int maxAlign;
	I64 totalSize = memnodes_calc_offsets(list, &maxAlign);  // Update the offset_from_orign field

	// Very imortant to align the start of the mem at 64, so the aligments of other blocks are honored
	I08 *paligned  = mem_alloc(self, totalSize, maxAlign);   // REturn NULL if ootalSize <=0
	memnodes_assign_from_alignedbase(list, paligned);
}


void mem_init(MemPointers* self) {	 

	*self = (MemPointers) {
			.alloc		        = mem_alloc,
			.alloc0             = mem_alloc0,
			.alloclist          = mem_alloc_list,
			.init		        = mem_init,
			.free_all	        = mem_free_all,
			.nptsMax            =0,
			.npts               =0,
			.bytesAllocated		=0,
			.memAlignOffset     =NULL,
			.memPointer         = NULL,
			.memHeaderBackup    =NULL,
			.checkHeader        =0,
			.verify_header     = verify_header,
			};			
}


/**************************************/
// Dynanic Buffer
/**************************************/

void dynbuf_init(DynMemBufPtr buf, int init_max_len) {

	buf->cur_len      = 0;

	if ( (size_t) init_max_len > buf->max_len) {		
		if (buf->raw) {
			free(buf->raw);			
			buf->raw = NULL;
		}
		buf->max_len = init_max_len;		
	}

	if (buf->raw == NULL) {		
		buf->raw = malloc(buf->max_len);
		return;
	}
}
	

void dynbuf_kill(DynMemBufPtr buf) {
	if (buf->raw) {
		free(buf->raw);
	}
	memset(buf, 0, sizeof(DynMemBuf));
}

void dynbuf_requestmore(DynMemBufPtr buf, int moreBytes) {

	// Request more bytes relative to the current size (cur_len)

	size_t newLength = (size_t) (buf->cur_len + moreBytes) ;

	if (newLength <= buf->max_len) {
		if (buf->raw == NULL) {
			buf->raw = malloc(buf->max_len);
			buf->cur_len = 0;
		}
		return;
	}
	else {
    // The data bufer need to be expaned
		newLength = (int) max(newLength, buf->max_len + buf->max_len / 2);
		int8_t* newptr;
		if (buf->cur_len == 0) {
			//f there is no exisintg data, so no need to keep a copy
			FreeNullify(buf->raw);
			newptr = malloc(newLength);
		}	else {
			// this exist data to be kept
			newptr = realloc(buf->raw, newLength);
		}
	
		if (newptr) {
			buf->max_len = newLength;
			buf->raw     = newptr;
		} else {
			buf->max_len = 0;
			buf->raw     = NULL;
		}
	}
 
}

void dynbuf_insert_bytes(DynMemBufPtr buf, char * eewData, int nbytes) {	 
	 dynbuf_requestmore(buf, nbytes);                  
	 memcpy(buf->raw + buf->cur_len, eewData, nbytes); // Do not append the NUL at the CODE_EOF
	 buf->cur_len += nbytes;
}

void dynbuf_insert_str(DynMemBufPtr buf, char* newstr) {
	int newstr_len = (int) strlen(newstr)+1;    // 1L for the extra NUL
	dynbuf_requestmore(buf, newstr_len);// 
	memcpy(buf->raw + buf->cur_len, newstr, newstr_len); // Do not append the NUL at the CODE_EOF
	buf->cur_len += newstr_len;
}


void dynbuf_alloc_list(DynMemBufPtr buf, MemNode* list) {

	// ALl the data in buf will be lost
	buf->cur_len = 0;
	buf->max_len = buf->raw == NULL ? 0 : buf->max_len; // Not needed but enforched here to make sure max_len=-0 when raw=NULL

	int maxALign;
	I64 totalSize = memnodes_calc_offsets(list, &maxALign);  // Update the offset_from_orign field

	VOIDPTR paligned      = AlignedPointer(buf->raw, maxALign);
	int     offset        = (int) ((char*)paligned - (char*)buf->raw);
	int     bytesAvailble = (int) (buf->max_len - offset);

	if (bytesAvailble< totalSize) {
		// No enough soace available, need to allocate new MEM
		FreeNullify(buf->raw);
		MemAlignedPtr resultPtr = malloc_aligned(totalSize, maxALign);
		buf->raw     = resultPtr.ptr;
		buf->max_len = resultPtr.byte_allocated;
		paligned     = resultPtr.alignedptr;		
	}  


	buf->cur_len = totalSize + ((char*)paligned - (char*)buf->raw);
	memnodes_assign_from_alignedbase(list, paligned);
}

/**************************************/
// Dynanic Algined Buffer
/*************************************/

 void adynbuf_init(DynAlignedBufPtr buf, int init_max_len) {

	 buf->cur_len = 0;

	 if (buf->elem_size == 0 || buf->align == 0) {
		 r_printf("ERROR: elem_size and algin should not be zeros (in abynbuf_nit).\n");
		 return;
	 }

	 
	 if ((size_t)init_max_len > buf->max_len) {
		 buf->max_len = init_max_len;
		 if (buf->p.raw) {
			 free(buf->p.raw - buf->offset);
			 buf->p.raw = NULL;
		 }	 
	 }
 
 
	 if (buf->max_len>0 && buf->p.raw == NULL) {
		 //int bufsize = buf->max_len * buf->elem_size + (buf->align - 1);
		 // for efficeny, we just request one extra byte to make sure the total size is a nice even number
		 // but in theory, only a number 'align-1' of extra bytes are needed
		 int bufsize = buf->max_len * buf->elem_size + (buf->align);

		 // 64: (64-1)=0x3F=0b11 1111 ~(64-1)=11111111 00 0000
		 char* ptr    = malloc(bufsize);
		 char* palign = AlignedPointer(ptr, buf->align);

		 buf->p.raw  = palign;
		 buf->offset = (int) (palign - ptr);	 
	 }

 }
 void adynbuf_kill(DynAlignedBufPtr buf) {
	 if ((*buf).p.raw) {
		 free(buf->p.raw - buf->offset);
	 }
	 memset(buf, 0, sizeof(DynAlignedBuf));
 }

 void adynbuf_requestmore(DynAlignedBufPtr buf, int moreElements) {

	 size_t newLength = moreElements + buf->cur_len;
	 if (newLength <= buf->max_len) {
		 return;
	 }
	 newLength         = max(newLength, buf->max_len + buf->max_len / 2);

	 int newbufsize    = newLength * buf->elem_size + (buf->align);
	 int8_t* newptr    = realloc(buf->p.raw - buf->offset, newbufsize);
	 int8_t* newpalign = (VOID_PTR)(((uintptr_t)newptr + buf->align - 1) & ~(uintptr_t)(buf->align - 1));	 
	 if (newptr) {
		 		 
		 int  newoffset = (int) (newpalign - newptr);

		 if (newoffset == buf->offset) {			 
			 buf->p.raw     = newpalign;
			 buf->max_len = newLength;
		 }
		 else if (newoffset < buf->offset) {			 
			 // new |  --------------------
			 // old |       ----------------- 
			 memcpy(newpalign,newptr+buf->offset, buf->elem_size * buf->cur_len);
			 buf->p.raw    = newpalign;
			 buf->offset = newoffset;
			 buf->max_len = newLength;
		 }  else if (newoffset > buf->offset) {		 
			 int8_t* newptr1    =  malloc( newbufsize);
			 int8_t* newpalign1 = (VOID_PTR)(((uintptr_t)newptr1 + buf->align - 1) & ~(uintptr_t)(buf->align - 1));
		 
			 memcpy(newpalign1, newptr + buf->offset, buf->elem_size * buf->cur_len);
			 buf->p.raw   = newpalign1;
			 buf->offset  = (int) (newpalign1 - newptr1);
			 buf->max_len = newLength;
			 free(newptr); 
		 }
	 }
	 

 }
 
#include "abc_000_warning.h"
