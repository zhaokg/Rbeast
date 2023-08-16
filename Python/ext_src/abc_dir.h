#pragma once

#include "abc_dirent.h"
#include "abc_datatype.h"
typedef struct DIR_FILELIST {	
	I32     num;
	char *  base;
	ptrdiff_t  *  offset;
} FILELIST, *_restrict FILELIST_PTR ;

FILELIST_PTR GetFlist(const char *path, const char * ext);
void FreeFlist(FILELIST_PTR flist);

void PrintFlist(FILELIST_PTR flist);