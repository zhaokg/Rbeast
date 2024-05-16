#include "abc_000_macro.h"
#include "abc_000_warning.h"

#if defined(_MSC_VER)
	#ifndef _CRT_SECURE_NO_DEPRECATE
		#define _CRT_SECURE_NO_DEPRECATE (1)
		#define _CRT_NONSTDC_NO_DEPRECATE
	#endif
	//#pragma warning(disable : 4996)
#endif


#if defined(OS_WIN64) || defined (OS_WIN32)

#include "abc_dir.h"
#include "abc_ide_util.h" //r_printf
// github.com/tronkko/dirent/tree/master/include
// codeforwin.org/2018/03/c-program-to-list-all-files-in-a-directory-recursively.html
void listFiles(const char *path, const char * ext)
{
	DIR *dir = opendir(path);
	// Unable to open directory stream
	if (dir == NULL) return ; 

	struct dirent *dp;
	while ((dp = readdir(dir)) != NULL)
	{
		char  fullpath[1000];
		strcpy(fullpath, path);
		strcat(fullpath, "/");
		strcat(fullpath, dp->d_name);

		// stackoverflow.com/questions/22886290/c-get-all-files-with-certain-extension
		if (dp->d_type == DT_REG /*DT_DIR*/) {
			char * const pext = strrchr(dp->d_name, '.');
			if (pext != NULL && pext != dp->d_name)		{
				if (strcmp(pext, ".tif") == 0){
					r_printf("%s\n", dp->d_name);
				}				
			}
		}

		/*
		struct stat s;
		if (stat(fullpath, &s) == 0)
		{
		if (s.st_mode & S_IFDIR) 		;	//it's a directory;
		else if (s.st_mode & S_IFREG)	;//it's a file
		}*/

	} //while ( (dp = readdir(dir)) != NULL)

	closedir(dir); // Close directory stream
}

FILELIST_PTR GetFlist(const char *path, const char * ext)
{
	DIR *dir = opendir(path);
	// Unable to open directory stream
	if (dir == NULL)return NULL;

	int fNum    = 0;
	int memSize = 0;

	struct dirent *dp;
	while ((dp = readdir(dir)) != NULL)
	{
	  // stackoverflow.com/questions/22886290/c-get-all-files-with-certain-extension
		if (dp->d_type == DT_REG /*DT_DIR*/) {
			char * const pEXT = strrchr(dp->d_name, '.');
			if (pEXT != NULL && pEXT != dp->d_name)
			{
				if (stricmp(pEXT+1, ext) == 0){
					fNum++;
					memSize += (int)strlen(dp->d_name)+ 1L;
				}
			}
		}

		/*
		struct stat s;
		if (stat(fullpath, &s) == 0)
		{
		if (s.st_mode & S_IFDIR) 		;	//it's a directory;
		else if (s.st_mode & S_IFREG)	;//it's a file
		}*/

	} //while ( (dp = readdir(dir)) != NULL)


	/*************************************************/
	FILELIST_PTR flist = malloc(sizeof(FILELIST));
	memset(flist, 0, sizeof(FILELIST));

	flist->num         = fNum;

	flist->base        = malloc(memSize);
	memset(flist->base, 0, memSize);
	flist->offset      = malloc(fNum*sizeof(ptrdiff_t));

	fNum = 0;
	char * ptr = flist->base;	

	rewinddir(dir);
	while ((dp = readdir(dir)) != NULL) {
		// stackoverflow.com/questions/22886290/c-get-all-files-with-certain-extension
		if (dp->d_type == DT_REG /*DT_DIR*/) {
			char * const pEXT = strrchr(dp->d_name, '.');
			if (pEXT != NULL && pEXT != dp->d_name)
			{
				if (stricmp(pEXT+1, ext) == 0){
					int size = (int)strlen(dp->d_name);
					flist->offset[fNum] =  ptr; //-Wint-conversion
					memcpy(ptr, dp->d_name, size);
					ptr[size] = 0;

					fNum++;
					ptr += (size + 1);
				}
			}
		}

		/*
		struct stat s;
		if (stat(fullpath, &s) == 0)
		{
		if (s.st_mode & S_IFDIR) 		;	//it's a directory;
		else if (s.st_mode & S_IFREG)	;//it's a file
		}*/

	} //whil
	/*************************************************/



	closedir(dir); // Close directory stream

	return flist;
}
 
void FreeFlist(FILELIST_PTR flist) {
	if (flist == NULL) return;
	if (flist->base != NULL)
		free(flist->base);
	if (flist->offset != NULL)
		free(flist->offset);
	free(flist);
}

void PrintFlist(FILELIST_PTR flist) {
	if (flist == NULL) return;
	for (int i = 0; i < flist->num; i++){
		r_printf("%s \n", (char*) flist->offset[i]);
	}
}

#else
static char _dummy = 'c';
#endif


#include "abc_000_warning.h"