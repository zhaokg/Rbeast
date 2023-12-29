#pragma once

//#define  PLAY_MODE
//#define  M_RELEASE
//#define  P_RELEASE
//#define  R_RELEASE

// https://social.msdn.microsoft.com/Forums/vstudio/en-US/355ed7af-4037-4587-8614-34d51d865f03/missing-prototype-warning?forum=vclanguage
// In MSVC, set the warning level to Level 3 to get warnings of unfctions without prototypes
//ERROR: function returning a value
//ERROR: undefined; assuming extern returning int

//gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 
// file externsion: .mexw64 for Matlab, .pyd for Pyhton

#define R_INTERFACE   0
#define M_INTERFACE   0
#define P_INTERFACE    0
/*------------------------------------------------------------*/
#define MYMAT_LIBRARY   1
#define MKL_LIBRARY     0
#define MATLAB_LIBRARY  0 // For some uknown reason, the Matlab version of Blas and Lapack function doesn't work.

/*------------------------------------------------------------*/
#define PCGRAND_LIBRARY 1
#define MKLRAND_LIBRARY 0
 
//#define  R_RELEASE
//#define  P_RELEASE
//#define  M_RELEASE

#ifdef R_RELEASE
        // For R interface
		#undef   R_INTERFACE
		#undef   M_INTERFACE
		#undef   P_INTERFACE
		#define  R_INTERFACE 1
		#define  M_INTERFACE 0
        #define  P_INTERFACE 0

		#undef   MYMAT_LIBRARY
	    #undef   MKL_LIBRARY

		#define MYMAT_LIBRARY 1
	    #define MKL_LIBRARY   0

	    #define PCGRAND_LIBRARY 1
        #define MKLRAND_LIBRARY 0
#endif
#ifdef P_RELEASE
        // For R interface
		#undef   R_INTERFACE
		#undef   M_INTERFACE
		#undef   P_INTERFACE
		#define  R_INTERFACE 0
		#define  M_INTERFACE 0
        #define  P_INTERFACE 1

		#undef   MYMAT_LIBRARY
	    #undef   MKL_LIBRARY

		#define MYMAT_LIBRARY 1
	    #define MKL_LIBRARY   0

	    #define PCGRAND_LIBRARY 1
        #define MKLRAND_LIBRARY 0
#endif

/*
*  Note for compling mex using max in Matlab:
*  mex can't handle the mixed inputs of C and CPP. Check these
*  https://stackoverflow.com/questions/57579194/matlab-r2016b-mex-fails-to-compile-c-code
*  https://www.mathworks.com/matlabcentral/answers/827145-how-can-i-split-c-c-code-in-multiple-files-and-still-use-mex
*  mex accept a wildcard for source files (*.c) but not a wildcard for object files (*.o)
* 
* mex -v -c CFLAGS='$CFLAGS -DM_RELEASE  -Wall -v -Wl,-v'      *.c     // compiling all the C source
* mex -v -c CXXFLAGS='$CXXFLAGS -DM_RELEASE -Wall -v -Wl,-v'  *.cpp   // compiling all the C source
* mex -v  CFLAGS='$CFLAGS  -Wall -Wl,-v' -lmwservices -lut abc_ioFlush.o *.o -output Rbeast  // for some reason, the output is bad
* 
* Below is my solution:
  system("gcc -c -fPIC  -pthread -DNDEBUG -DM_RELEASE   -DMATLAB_DEFAULT_RELEASE=R2017b  -DMATLAB_MEX_FILE  -I/MATLAB/extern/include   -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  *.c")
  system("g++ -c -fPIC  -pthread -DNDEBUG -DM_RELEASE   -DMATLAB_DEFAULT_RELEASE=R2017b  -I/MATLAB/extern/include   -O2 -Wall    -mfpmath=sse -msse2 -mstackrealign  *.cpp")
  system("gcc -shared -pthread  -L/MATLAB/bin/glnxa64 -lmx -lmex -lmat -lm  -lut -lmwservices *.o -o Rbeast.mexa64")
*/

#ifdef M_RELEASE
        // For Matlab interface
		#undef   R_INTERFACE
		#undef   M_INTERFACE
		#undef   P_INTERFACE
		#define  R_INTERFACE 0
		#define  M_INTERFACE 1
        #define  P_INTERFACE 0

		#undef   MYMAT_LIBRARY
	    #undef   MKL_LIBRARY

		#define MYMAT_LIBRARY 1
	    #define MKL_LIBRARY   0

	    #define PCGRAND_LIBRARY 1
        #define MKLRAND_LIBRARY 0
        
        #ifndef MATLAB_MEX_FILE
			#define MATLAB_MEX_FILE
        #endif
        #define MATLAB_DEFAULT_RELEASE  R2017b
#endif
/*------------------------------------------------------------*/

//http:// nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
//https:// stackoverflow.com/questions/2166483/which-macro-to-wrap-mac-os-x-specific-code-in-c-c
//https:// blog.kowalczyk.info/article/j/guide-to-predefined-macros-in-c-compilers-gcc-clang-msvc-etc..html

#ifdef _MSC_VER
	#define MSVC_COMPILER
//#elif defined(__GNUC__) || defined(__clang__)  || defined(__APPLE__) || defined(__linux__) || defined(__MINGW32__) || defined(__MINGW64__) ||defined(__MACH__)||defined(__SUNPRO_C)||defined(__SUNPRO_CC)
#elif defined(__clang__)
	#define CLANG_COMPILER
#elif (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
	#define GCC_COMPILER
#elif defined(__SUNPRO_C)||defined(__SUNPRO_CC)
	#define SOLARIS_COMPILER
#endif


//nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system
#if defined(_WIN64) || defined (__MINGW64__) && defined(_WIN32) && !defined(__i386)  && !defined(__i686) && !defined(i386) && !defined(__i686)
	#define WIN64_OS
#elif defined(_WIN32) && !defined(_WIN64)
	#define WIN32_OS
#endif

#if defined(__APPLE__ ) && defined (__MACH__) 
	#define MAC_OS
#endif

#if defined(__linux__) 
	//https://stackoverflow.com/questions/142508/how-do-i-check-os-with-a-preprocessor-directive
	#define LINUX_OS
    #ifndef _GNU_SOURCE
		#define _GNU_SOURCE // for including CPU_ZERO in sched.h
    #endif
#endif

#if (defined(unix) || defined(__unix__) || defined(__unix) ) && !defined(__APPLE__)
	#define UNIX_OS
	#ifndef _GNU_SOURCE
		#define _GNU_SOURCE // for including CPU_ZERO in sched.h
	#endif
#endif

#if defined(sun) && defined(__sun) && defined(__SVR4) 
	#define SOLARIS_OS
	
	//https://web.archive.org/web/20191012035921/http://nadeausoftware.com/articles/2012/01/c_c_tip_how_use_compiler_predefined_macros_detect_operating_system#Solaris
#endif




//https://stackoverflow.com/questions/735647/ifdef-for-32-bit-platform
#if _WIN64 || __amd64__ || defined(__LP64__) || (defined(__x86_64__) &&    !defined(__ILP32__) ) || defined(_M_X64) || defined(__ia64) || defined (_M_IA64) || defined(__aarch64__) || defined(__powerpc64__)
	#define TARGET_64
#else
	#define TARGET_32
#endif

#if __GNUC__
	#if __x86_64__ || __ppc64__
		#define TARGET_64
	#else
		#define TARGET_32
	#endif
#endif

#if defined(__aarch64__)
	//https://stackoverflow.com/questions/60588765/how-to-get-cpu-brand-information-in-arm64
    //https://stackoverflow.com/questions/23934862/what-predefined-macro-can-i-use-to-detect-the-target-architecture-in-clang

    #define  ARM64_OS
#endif


//https://stackoverflow.com/questions/152016/detecting-cpu-architecture-compile-time
#if defined(__powerpc) || defined(__powerpc__) || defined(__powerpc64__)  || defined(__POWERPC__) || defined(__ppc__) || defined(__PPC__) || defined(_ARCH_PPC) \
     || defined(__PPC64__) || defined(__ppc64__) || defined(_ARCH_PPC64)
	#define POWERPC_OS
#endif

#if defined(TARGET_32) && defined (MSVC_COMPILER)
    //'strcpy': This function or variable may be unsafe. Consider using strcpy_s instead. To disable deprecation, use _CRT_SECURE_NO_WARNING
	#define _CRT_SECURE_NO_WARNINGS
	#pragma warning (disable: 4703) //potentially uninitialized local pointer variable

#endif


/*------------------------------------------------------------*/
#if MYMAT_LIBRARY==1
		#define PCGRAND_LIBRARY 1
		#define MKLRAND_LIBRARY 0
#endif


/*------------------------------------------------------------*/
#if defined(MSVC_COMPILER)
		#define INLINE    __inline
		#define _restrict __restrict
        #define UNUSED_DECORATOR 
#elif defined(SOLARIS_COMPILER) 
	//#if defined(__SUNPRO_C)||defined(__SUNPRO_CC)
    //https: //docs.oracle.com/cd/E24457_01/html/E21990/gipgw.html
		#define INLINE     inline // __inline__
		#define _restrict _Restrict //available for both -xc99=none and -xc99=all
        #define UNUSED_DECORATOR 
#elif defined(GCC_COMPILER) || defined(CLANG_COMPILER)
		#define INLINE     inline
		#define _restrict __restrict__		
        #define UNUSED_DECORATOR  __attribute__((unused))
#endif

 
// CHeck if it is the MUSL-linux C library
// https://stackoverflow.com/questions/58177815/how-to-actually-detect-musl-libc
#if   defined(LINUX_OS) && ( defined(GCC_COMPILER) || defined(CLANG_COMPILER) )

	#ifdef _GNU_SOURCE 
		#include <features.h>
		#ifndef __USE_GNU
			#define __MUSL__ 
		#endif
	#else
		#define _GNU_SOURCE
		#include <features.h>
		#ifndef __USE_GNU
		#define __MUSL__ 
		#endif
		#undef _GNU_SOURCE /* don't contaminate other includes unnecessarily */
	#endif

#endif




/* yes I know, the top of this file is quite ugly */
#ifdef MSVC_COMPILER
    # define ALIGN32_BEG __declspec(align(32))
    # define ALIGN32_END 
	//https://stackoverflow.com/questions/4750880/can-i-treat-a-specific-warning-as-an-error
	// #pragma warning (error: 4013)
#else
    # define ALIGN32_BEG
    # define ALIGN32_END __attribute__((aligned(32)))
#endif


/*------------------------------------------------------------*/
	#define DIAG_STR(s) #s
	#define DIAG_JOINSTR(x,y) DIAG_STR(x ## y)
	#ifdef MSVC_COMPILER
		#define DIAG_DO_PRAGMA(x) __pragma (#x)
		#define DIAG_PRAGMA(compiler,x) DIAG_DO_PRAGMA(warning(x))
	#else
		#define DIAG_DO_PRAGMA(x)       _Pragma (#x)
		#define DIAG_PRAGMA(compiler,x) DIAG_DO_PRAGMA(compiler diagnostic x)
	#endif

	#if defined(CLANG_COMPILER)
		# define DISABLE_WARNING(gcc_unused,clang_option,msvc_unused) DIAG_PRAGMA(clang,push) DIAG_PRAGMA(clang,ignored DIAG_JOINSTR(-W,clang_option))
		# define ENABLE_WARNING(gcc_unused,clang_option,msvc_unused) DIAG_PRAGMA(clang,pop)
	#elif defined(MSVC_COMPILER)
		# define DISABLE_WARNING(gcc_unused,clang_unused,msvc_errorcode) DIAG_PRAGMA(msvc,push) DIAG_DO_PRAGMA(warning(disable:##msvc_errorcode))
		# define ENABLE_WARNING(gcc_unused,clang_unused,msvc_errorcode) DIAG_PRAGMA(msvc,pop)
	#elif defined(GCC_COMPILER)
		#if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,push) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,pop)
		#else
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_option,msvc_unused) DIAG_PRAGMA(GCC,warning DIAG_JOINSTR(-W,gcc_option))
		#endif
	#endif
	//https:// stackoverflow.com/questions/3378560/how-to-disable-gcc-warnings-for-a-few-lines-of-code
 

//https: //stackoverflow.com/questions/39821164/how-do-you-define-a-multiline-macro-in-c
//https: //www.geeksforgeeks.org/multiline-macros-in-c/
#if defined(GCC_COMPILER)  

//https: //github.com/BioMedIA/MIRTK/issues/230
	#define  DISABLE_MANY_WARNINGS   \
	DISABLE_WARNING(unknown-pragmas, unknown-pragmas, NOT_USED) \
	DISABLE_WARNING(pragmas, pragmas, NOT_USED) \
	DISABLE_WARNING(unused-variable, unused-variable, NOT_USED) \
	DISABLE_WARNING(unused-function, unused-function, NOT_USED) \
	DISABLE_WARNING(pointer-sign, pointer-sign, NOT_USED) \
	DISABLE_WARNING(implicit-function-declaration, implicit-function-declaration, NOT_USED) \
	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED) \
	DISABLE_WARNING(unused-but-set-variable, unused-but-set-variable, NOT_USED) \
	DISABLE_WARNING(maybe-uninitialized, maybe-uninitialized, NOT_USED)\
	DISABLE_WARNING(pointer-to-int-cast, pointer-to-int-cast, NOT_USED)\
	DISABLE_WARNING(misleading-indentation, NOT_USED, NOT_USED)\
	DISABLE_WARNING(discarded-qualifiers, discarded-qualifiers, NOT_USED)\
	DISABLE_WARNING(int-to-pointer-cast, int-to-pointer-cast, NOT_USED)\
	DISABLE_WARNING(unused-result, unused-result, NOT_USED) \
    DISABLE_WARNING(unused-const-variable, unused-const-variable, NOT_USED)\
    DISABLE_WARNING(incompatible-pointer-types-discards-qualifiers, incompatible-pointer-types-discards-qualifiers, NOT_USED)\
	DISABLE_WARNING(incompatible-pointer-types,incompatible-pointer-types, NOT_USED)\
	DISABLE_WARNING(self-assign,self-assign, NOT_USED) \
    DISABLE_WARNING(unused-value,unused-value, NOT_USED) \
    DISABLE_WARNING(int-conversion,int-conversion, NOT_USED) \
    DISABLE_WARNING(restrict,restrict, NOT_USED)\
    DISABLE_WARNING(switch, switch, NOT_USED) \
    DISABLE_WARNING(uninitialized, uninitialized, NOT_USED)\
    DISABLE_WARNING(pedantic, pedantic, NOT_USED) \
    DISABLE_WARNING(div-by-zero,div-by-zero, NOT_USED)

	#define  ENABLE_MANY_WARNINGS   \
    ENABLE_WARNING(div-by-zero,div-by-zero, NOT_USED) \
    ENABLE_WARNING(pedantic, pedantic, NOT_USED)\
    ENABLE_WARNING(uninitialized, uninitialized, NOT_USED)\
    ENABLE_WARNING(switch, switch, NOT_USED)\
    ENABLE_WARNING(restrict, restrict, NOT_USED)\
    ENABLE_WARNING(int-conversion,int-conversion, NOT_USED) \
    ENABLE_WARNING(unused-value,unused-value, NOT_USED) \
	ENABLE_WARNING(self-assign,self-assign, NOT_USED) \
	ENABLE_WARNING(incompatible-pointer-types,incompatible-pointer-types, NOT_USED) \
	ENABLE_WARNING(incompatible-pointer-types-discards-qualifiers, incompatible-pointer-types-discards-qualifiers, NOT_USED)\
    ENABLE_WARNING(unused-const-variable, unused-const-variable, NOT_USED)\
	ENABLE_WARNING(unused-result, unused-result, NOT_USED)\
	ENABLE_WARNING(int-to-pointer-cast, int-to-pointer-cast, NOT_USED)\
	ENABLE_WARNING(discarded-qualifiers, discarded-qualifiers, NOT_USED)\
	ENABLE_WARNING(misleading-indentation, NOT_USED, NOT_USED)\
	ENABLE_WARNING(pointer-to-int-cast, pointer-to-int-cast, NOT_USED)\
	ENABLE_WARNING(maybe-uninitialized, maybe-uninitialized, NOT_USED)\
	ENABLE_WARNING(unused-but-set-variable, unused-but-set-variable, NOT_USED) \
	ENABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED) \
	ENABLE_WARNING(implicit-function-declaration, implicit-function-declaration, NOT_USED) \
	ENABLE_WARNING(pointer-sign, pointer-sign, NOT_USED) \
	ENABLE_WARNING(unused-function, unused-function, NOT_USED) \
	ENABLE_WARNING(unused-variable, unused-variable, NOT_USED)  \
	ENABLE_WARNING(pragmas, pragmas, NOT_USED) \
	ENABLE_WARNING(unknown-pragmas,unknown-pragmas, NOT_USED) 

#elif defined(CLANG_COMPILER)  
//https://stackoverflow.com/questions/14261534/temporarily-overwrite-a-macro-in-c-preprocessor
//https://clang.llvm.org/doxygen/classclang_1_1Preprocessor.html#a04dec9fbfa220dfea433bcbeffa270c3
//https://gcc.gnu.org/onlinedocs/gcc-5.4.0/gcc/Push_002fPop-Macro-Pragmas.html

	#define  DISABLE_MANY_WARNINGS  \
	DISABLE_WARNING(unknown-pragmas, unknown-pragmas, NOT_USED)  \
	DISABLE_WARNING(pragmas, pragmas, NOT_USED) \
	DISABLE_WARNING(unused-variable, unused-variable, NOT_USED)  \
	DISABLE_WARNING(unused-function, unused-function, NOT_USED)  \
	DISABLE_WARNING(pointer-sign, pointer-sign, NOT_USED)  \
	DISABLE_WARNING(implicit-function-declaration, implicit-function-declaration, NOT_USED) \
	DISABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED) \
	DISABLE_WARNING(xxxxxxx, unknown-warning-option, xxxxxxxxx) \
    DISABLE_WARNING(pointer-to-int-cast, pointer-to-int-cast, NOT_USED)\
	DISABLE_WARNING(int-to-pointer-cast, int-to-pointer-cast, NOT_USED) \
	DISABLE_WARNING(unused-result, unused-result, NOT_USED)             \
    DISABLE_WARNING(unused-const-variable, unused-const-variable, NOT_USED) \
	DISABLE_WARNING(incompatible-pointer-types-discards-qualifiers, incompatible-pointer-types-discards-qualifiers, NOT_USED)\
	DISABLE_WARNING(incompatible-pointer-types,incompatible-pointer-types, NOT_USED)\
	DISABLE_WARNING(self-assign,self-assign, NOT_USED)  \
    DISABLE_WARNING(unused-value,unused-value, NOT_USED) \
    DISABLE_WARNING(int-conversion,int-conversion, NOT_USED) \
    DISABLE_WARNING(switch, switch, NOT_USED) \
    DISABLE_WARNING(uninitialized, uninitialized, NOT_USED)\
    DISABLE_WARNING(pedantic, pedantic, NOT_USED) \
    DISABLE_WARNING(typedef-redefinition, typedef-redefinition, NOT_USED) \
    DISABLE_WARNING(div-by-zero,div-by-zero, NOT_USED)

	/*DISABLE_WARNING(restrict, restrict, NOT_USED)\*/
	/*ENABLE_WARNING(restrict, restrict, NOT_USED)\*/

	//https: //github.com/nasa/trick/issues/600
	//Clang no longer supports -Wno-unused-but-set-variable

	//https://clang.llvm.org/docs/DiagnosticsReference.html#wunknown-warning-option
	//Supress unknown warming by disabling "unknown-warning-option"
	//https://stackoverflow.com/questions/41673546/clang-warning-warning-unknown-warning-option-wno-maybe-uninitialized

	//DISABLE_WARNING(unused-but-set-variable, unused-but-set-variable, NOT_USED) 
	//DISABLE_WARNING(maybe-uninitialized, maybe-uninitialized, NOT_USED)

	//https://clang.llvm.org/docs/DiagnosticsReference.html#wpragmas
	//https://clang.llvm.org/docs/DiagnosticsReference.html#wrestrict-expansion
	#define  ENABLE_MANY_WARNINGS  \
    ENABLE_WARNING(div-by-zero,div-by-zero, NOT_USED) \
    ENABLE_WARNING(typedef-redefinition, typedef-redefinition, NOT_USED)\
    ENABLE_WARNING(pedantic, pedantic, NOT_USED)\
    ENABLE_WARNING(uninitialized, uninitialized, NOT_USED)\
    ENABLE_WARNING(switch, switch, NOT_USED)\
    ENABLE_WARNING(int-conversion,int-conversion, NOT_USED) \
    ENABLE_WARNING(unused-value,unused-value, NOT_USED) \
	ENABLE_WARNING(self-assign,self-assign, NOT_USED) \
	ENABLE_WARNING(incompatible-pointer-types,incompatible-pointer-types, NOT_USED) \
	ENABLE_WARNING(incompatible-pointer-types-discards-qualifiers, incompatible-pointer-types-discards-qualifiers, NOT_USED)\
    ENABLE_WARNING(unused-const-variable, unused-const-variable, NOT_USED)\
	ENABLE_WARNING(unused-result, unused-result, NOT_USED)\
	ENABLE_WARNING(int-to-pointer-cast, int-to-pointer-cast, NOT_USED)\
	ENABLE_WARNING(pointer-to-int-cast, pointer-to-int-cast, NOT_USED)\
	ENABLE_WARNING(xxxxxxx, unknown-warning-option, xxxxxxxxx) \
	ENABLE_WARNING(strict-aliasing, strict-aliasing, NOT_USED) \
	ENABLE_WARNING(implicit-function-declaration, implicit-function-declaration, NOT_USED) \
	ENABLE_WARNING(pointer-sign, pointer-sign, NOT_USED) \
	ENABLE_WARNING(unused-function, unused-function, NOT_USED) \
	ENABLE_WARNING(unused-variable, unused-variable, NOT_USED) \
	ENABLE_WARNING(pragmas, pragmas, NOT_USED) \
	ENABLE_WARNING(unknown-pragmas,unknown-pragmas, NOT_USED) 

#else
	#define  DISABLE_MANY_WARNINGS 
	#define  ENABLE_MANY_WARNINGS 	
#endif

#ifdef _WIN32_WINNT
#undef  _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif

/*
//Needed to put after the definiton of WIN_WINNT
In file included from C:/rtools40/mingw32/i686-w64-mingw32/include/crtdefs.h:10,
				 from C:/rtools40/mingw32/i686-w64-mingw32/include/stdint.h:28,
				 from C:/rtools40/mingw32/lib/gcc/i686-w64-mingw32/8.3.0/include/stdint.h:9,
				 from abc_000_macro.h:110,
				 from _beastv2_gui_plot.c:1:
C:/rtools40/mingw32/i686-w64-mingw32/include/_mingw.h:225: note: this is the location of the previous definition
 #define _WIN32_WINNT 0x502

*/
//https://stackoverflow.com/questions/1505582/determining-32-vs-64-bit-in-c
#include <stdint.h>
#if   INTPTR_MAX == INT32_MAX
	#define TARGET_32
#elif INTPTR_MAX == INT64_MAX
	#define TARGET_64
#else
	#error "Environment not 32 or 64-bit."
#endif

//https://stackoverflow.com/questions/45477355/difference-between-pragma-and-pragma-in-c
/*
//https://mikejsavage.co.uk/blog/cpp-tricks-disable-optimisations-macro.html
#  define DISABLE_OPTIMISATIONS() \
        _Pragma( "GCC push_options" ) \
        _Pragma( "GCC optimize (\"O0\")" )
*/

//https://stackoverflow.com/questions/45477355/difference-between-pragma-and-pragma-in-c
//_Pragma ("GCC dependency \"parse.y\"")

#define  CHANGE_TO_AVX_GCC  \
          DIAG_DO_PRAGMA(GCC optimization_level 3) \
          DIAG_DO_PRAGMA(GCC optimize("O3,Ofast,inline,omit-frame-pointer,no-asynchronous-unwind-tables")) \
          DIAG_DO_PRAGMA(GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,avx,avx2,fma,tune=haswell"))  



#define _in_
#define _inout_
#define _out_




#define mv(n, src, dest)	r_cblas_scopy( n,src, 1L, dest, 1L) 
#define cp(n, src, dest)    memcpy(dest, src, sizeof(F32)*(size_t)(n))
#define SCPY(n, src, dest)  memcpy(dest, src, sizeof(F32)*(size_t)(n))
#define FILL0(dest,n)       memset(dest, 0L,  sizeof(F32)*(size_t)(n))