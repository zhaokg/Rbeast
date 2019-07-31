
#pragma once
#ifdef _MSC_VER
	#define MSVC_COMPILER
#elif defined(__clang__)
	#define CLANG_COMPILER
#elif (defined(__GNUC__)||defined(__GNUG__)) && !(defined(__clang__)||defined(__INTEL_COMPILER))
	#define GCC_COMPILER
#elif defined(__SUNPRO_C)||defined(__SUNPRO_CC)
	#define SOLARIS_COMPILER
#endif
#if defined(_WIN64)||defined (__MINGW64__) && defined(_WIN32) && !defined(__i386)  && !defined(__i686) && !defined(i386) && !defined(__i686)
	#define WIN64_OS
#endif
	#define MYMAT_LIBRARY 1
	#define MKL_LIBRARY   0
	#define MATLAB_LIBRARY 0 
	#define R_INTERFACE 0
	#define M_INTERFACE 1
	#define MYRAND_LIBRARY  1
	#define MKLRAND_LIBRARY 0
	#define R_RELEASE   0
	#if R_RELEASE==1
		#undef   R_INTERFACE
		#undef   M_INTERFACE
		#define R_INTERFACE 1
		#define M_INTERFACE 0
		#undef   MYMAT_LIBRARY
	    #undef   MKL_LIBRARY
		#define MYMAT_LIBRARY 1
	    #define MKL_LIBRARY   0
    #endif
	#define BASIS_METHODS 2
	#define PTHREAD_INOUT 0
	#define MY_DEBUG 0
	#if MYMAT_LIBRARY==1
		#define MYRAND_LIBRARY 1
		#define MKLRAND_LIBRARY 0
	#endif
	#if defined(MSVC_COMPILER)
		#define INLINE    __inline
		#define _restrict __restrict
	#elif defined(SOLARIS_COMPILER) 
		#define INLINE     inline 
		#define _restrict _Restrict 
	#elif defined(GCC_COMPILER)||defined(CLANG_COMPILER)
		#define INLINE     inline
		#define _restrict __restrict__		
	#endif
	#define max(a,b)    (((a) > (b)) ? (a) : (b))
	#define min(a,b)    (((a) < (b)) ? (a) : (b))
	#define mv(n,src,dest)	r_cblas_scopy( n,src,1L,dest,1L) 
	#define cp(n,src,dest)    memcpy(dest,src,sizeof(float)*(n))
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
		#if ((__GNUC__ * 100)+__GNUC_MINOR__) >=406
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,push) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,pop)
		#else
			# define DISABLE_WARNING(gcc_option,clang_unused,msvc_unused) DIAG_PRAGMA(GCC,ignored DIAG_JOINSTR(-W,gcc_option))
			# define ENABLE_WARNING(gcc_option,clang_option,msvc_unused) DIAG_PRAGMA(GCC,warning DIAG_JOINSTR(-W,gcc_option))
		#endif
	#endif
#if defined(GCC_COMPILER)  
	#define  DISABLE_MANY_WARNINGS   \
	DISABLE_WARNING(unknown-pragmas,unknown-pragmas,NOT_USED) \
	DISABLE_WARNING(pragmas,pragmas,NOT_USED) \
	DISABLE_WARNING(unused-variable,unused-variable,NOT_USED) \
	DISABLE_WARNING(unused-function,unused-function,NOT_USED) \
	DISABLE_WARNING(pointer-sign,pointer-sign,NOT_USED) \
	DISABLE_WARNING(implicit-function-declaration,implicit-function-declaration,NOT_USED) \
	DISABLE_WARNING(strict-aliasing,strict-aliasing,NOT_USED) \
	DISABLE_WARNING(unused-but-set-variable,unused-but-set-variable,NOT_USED) \
	DISABLE_WARNING(maybe-uninitialized,maybe-uninitialized,NOT_USED)\
	DISABLE_WARNING(pointer-to-int-cast,pointer-to-int-cast,NOT_USED)\
	DISABLE_WARNING(misleading-indentation,NOT_USED,NOT_USED)\
	DISABLE_WARNING(discarded-qualifiers,NOT_USED,NOT_USED)\
	DISABLE_WARNING(int-to-pointer-cast,int-to-pointer-cast,NOT_USED)\
	DISABLE_WARNING(unused-result,unused-result,NOT_USED)
	#define  ENABLE_MANY_WARNINGS   \
	ENABLE_WARNING(unused-result,unused-result,NOT_USED)\
	ENABLE_WARNING(int-to-pointer-cast,int-to-pointer-cast,NOT_USED)\
	ENABLE_WARNING(discarded-qualifiers,NOT_USED,NOT_USED)\
	ENABLE_WARNING(misleading-indentation,NOT_USED,NOT_USED)\
	ENABLE_WARNING(pointer-to-int-cast,pointer-to-int-cast,NOT_USED)\
	ENABLE_WARNING(maybe-uninitialized,maybe-uninitialized,NOT_USED)\
	ENABLE_WARNING(unused-but-set-variable,unused-but-set-variable,NOT_USED) \
	ENABLE_WARNING(strict-aliasing,strict-aliasing,NOT_USED) \
	ENABLE_WARNING(implicit-function-declaration,implicit-function-declaration,NOT_USED) \
	ENABLE_WARNING(pointer-sign,pointer-sign,NOT_USED) \
	ENABLE_WARNING(unused-function,unused-function,NOT_USED) \
	ENABLE_WARNING(unused-variable,unused-variable,NOT_USED)  \
	ENABLE_WARNING(pragmas,pragmas,NOT_USED) \
	ENABLE_WARNING(unknown-pragmas,unknown-pragmas,NOT_USED) 
#elif defined(CLANG_COMPILER)  
	#define  DISABLE_MANY_WARNINGS  \
	DISABLE_WARNING(unknown-pragmas,unknown-pragmas,NOT_USED)  \
	DISABLE_WARNING(pragmas,pragmas,NOT_USED) \
	DISABLE_WARNING(unused-variable,unused-variable,NOT_USED)  \
	DISABLE_WARNING(unused-function,unused-function,NOT_USED)  \
	DISABLE_WARNING(pointer-sign,pointer-sign,NOT_USED)  \
	DISABLE_WARNING(implicit-function-declaration,implicit-function-declaration,NOT_USED) \
	DISABLE_WARNING(strict-aliasing,strict-aliasing,NOT_USED) \
	DISABLE_WARNING(xxxxxxx,unknown-warning-option,xxxxxxxxx) \
    DISABLE_WARNING(pointer-to-int-cast,pointer-to-int-cast,NOT_USED)\
	DISABLE_WARNING(int-to-pointer-cast,int-to-pointer-cast,NOT_USED) \
	DISABLE_WARNING(unused-result,unused-result,NOT_USED)
	#define  ENABLE_MANY_WARNINGS  \
	ENABLE_WARNING(unused-result,unused-result,NOT_USED)\
	ENABLE_WARNING(int-to-pointer-cast,int-to-pointer-cast,NOT_USED)\
	ENABLE_WARNING(pointer-to-int-cast,pointer-to-int-cast,NOT_USED)\
	ENABLE_WARNING(xxxxxxx,unknown-warning-option,xxxxxxxxx) \
	ENABLE_WARNING(strict-aliasing,strict-aliasing,NOT_USED) \
	ENABLE_WARNING(implicit-function-declaration,implicit-function-declaration,NOT_USED) \
	ENABLE_WARNING(pointer-sign,pointer-sign,NOT_USED) \
	ENABLE_WARNING(unused-function,unused-function,NOT_USED) \
	ENABLE_WARNING(unused-variable,unused-variable,NOT_USED) \
	ENABLE_WARNING(pragmas,pragmas,NOT_USED) \
	ENABLE_WARNING(unknown-pragmas,unknown-pragmas,NOT_USED) 
#else
	#define  DISABLE_MANY_WARNINGS 
	#define  ENABLE_MANY_WARNINGS 	
#endif
#if defined(GCC_COMPILER)||defined(CLANG_COMPILER) 
	#define UNUSED_DECORATOR  __attribute__((unused))
#else
	#define UNUSED_DECORATOR 
#endif
