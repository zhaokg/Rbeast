#pragma once

#include <inttypes.h>
#include "abc_000_macro.h"
#include "abc_datatype.h"

typedef struct {
	void**  addr;
	int     size;
	int     align;
	I64     offset_from_origin;                        
} MemNode;
/*
MemNode * list;
offset_from_origin:  No need to mannualy set and will be autommatically computed by memnodes_calc_offsets
list[0].offset_from_origin is always 0, so it is used to store the total number of nodes 
The last list note is a sentinel one,  so its size and algin fields are used to save the toteal required mem size and the max
alignment required for the whoe list
*/

I64    memnodes_calc_offsets(MemNode* list, int* maxAlignment);
void   memnodes_assign_from_alignedbase(MemNode* list, VOIDPTR prawaligned);  // Assume that the buf has enough space
void   memnodes_assign_from_unalignedbase(MemNode* list, VOIDPTR pbase, int bufsize);
I64    memnodes_find_max_common_size(MemNode* list, ...);

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


 /**************************************/
 // Dynanic Buffer
 /**************************************/

 typedef struct {
	 int8_t* raw;
	 size_t  max_len;
	 size_t  cur_len;
 } DynMemBuf, * _restrict  DynMemBufPtr;

 void dynbuf_init(DynMemBufPtr buf, int init_max_len);
 void dynbuf_kill(DynMemBufPtr buf);
 void dynbuf_requestmore(DynMemBufPtr buf, int moreBytes);
 void dynbuf_insert_bytes(DynMemBufPtr buf, char* newstr, int nbytes);
 void dynbuf_insert_str(DynMemBufPtr buf, char* newstr);

 void dynbuf_alloc_list(DynMemBufPtr buf, MemNode *nodes);

 /**************************************/
 // Dynanic Algined Buffer
 /**************************************/

 typedef struct {
	 union {
		 char   * raw;
		 int8_t * i08;
		 int16_t* i16;
		 int32_t* i32;
		 int64_t* i64;
		 float*   f32;
		 double*  f64;
	 } p;
	 size_t  max_len;
	 size_t  cur_len;
	 int     elem_size;
	 int     align;
	 int     offset;
 } DynAlignedBuf, * _restrict DynAlignedBufPtr;

 void adynbuf_init(DynAlignedBufPtr buf, int init_max_len);
 void adynbuf_kill(DynAlignedBufPtr buf);
 void adynbuf_requestmore(DynAlignedBufPtr buf, int moreBytes);
 void adynbuf_insert_bytes(DynAlignedBufPtr buf, char* newstr);
 void adynbuf_insert_str(DynAlignedBufPtr buf, char* newstr);

