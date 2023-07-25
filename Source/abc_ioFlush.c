#include "abc_000_macro.h"
#include "abc_000_warning.h"

#include "abc_001_config.h"
#include "abc_datatype.h"

 #if M_INTERFACE ==1 &&  !defined(MSVC_COMPILER)

	#include "inttypes.h"
	#include "mex.h"
    //https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running/26271557
    //isoFlush is  an undocumented C++ function that resides in libmwservices.dll
	#if defined(WIN_OS)
		extern Bool ioFlush(void)  asm("?ioFlush@@YA_NXZ");
	#elif defined(MAC_OS)
      // Linux		
		extern Bool ioFlush(void)  asm("__Z7ioFlushv");
	#else 
		extern Bool ioFlush(void)  asm("_Z7ioFlushv");
	#endif

//https://stackoverflow.com/questions/35837694/how-to-manually-mangle-names-in-visual-c
		//https://gcc.gnu.org/onlinedocs/gcc/Asm-Labels.html
		// ON Windows, for MSVC C++ only, use the __identififer macro to give an alias
		// Using gcc, use  int a asm("xxx")

    void matlab_IOflush(void)	{
		#ifndef O_INTERFACE
			ioFlush();
		#endif
	}

	////https://stackoverflow.com/questions/10529500/what-does-this-mean-int-a
	/*
		in C++, (int &) is a type punning that can force converstion of a variable to int
	*/
	 
#else
static char achar UNUSED_DECORATOR ='c';
#endif

#include "abc_000_warning.h"



