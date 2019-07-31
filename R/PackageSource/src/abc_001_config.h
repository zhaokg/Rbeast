#pragma once
#include "abc_000_macro.h"
#ifdef MSVC_COMPILER
	#if R_INTERFACE==1
		#define DllExport   __declspec( dllexport ) 
	#elif M_INTERFACE==1
		#define DllExport 
		#pragma comment(linker,"/export:mexFunction")
    #endif
	#define QUOTE_IT(x) #x
	#define pthreadPath(__X) QUOTE_IT(C:/USERS/zhaokg/Documents/Visual Studio 2013/Projects/Matlab_Mex_Test/Pthread_IncludeLib/__X )
	#define matlabPath(__X) QUOTE_IT(D:/Program Files/MATLAB/MATLAB Production Server/R2015a/extern/include/__X )
	#define mklPath(__X) QUOTE_IT(C:/Program Files(x86)/IntelSWTools/compilers_and_libraries/windows/mkl/include/__X )
	#define ippPath(__X) QUOTE_IT(C:/Program Files(x86)/IntelSWTools/compilers_and_libraries/windows/ipp/include/__X )
	#define RPath(__X) QUOTE_IT(C:/Program Files/R/R-3.4.2/include/__X )
	#define pthreadLib(__X) QUOTE_IT(C:/USERS/zhaokg/Documents/Visual Studio 2013/Projects/Matlab_Mex_Test/Pthread_IncludeLib/__X )
	#define matlabLib(__X) QUOTE_IT(C:/Program Files/MATLAB/R2013b/extern/lib/win64/microsoft/__X )
	#define tbbLib(__X) QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/tbb/lib/intel64_win/vc_mt/__X )
	#define tbbLibstatic(__X) QUOTE_IT(G:/Intel_Library/tbb-2019_U3/build/vs2013/x64/Release-MT/__X )
	#define mklLib(__X) QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64_win/__X )
	#define ippLib(__X) QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/ipp/lib/intel64_win/__X )
   #define openMPLib(__X) QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64_win/__X )
	#define myMatLib(__X) QUOTE_IT(E:/BEAST_Code/blas/x64/Debug/__X )
	#define RLib(__X) QUOTE_IT(E:/BEAST_Code/beast_st_pthread_rw_Merge_tOrder_std_for_R_ST/R_import_lib/__X )
	#define FortranLib(__X) QUOTE_IT(C:/Program Files (x86)/Intel/Compiler/11.1/038/lib/intel64/__X )
	#if R_INTERFACE==1
		#pragma comment( lib,RLib(R.lib)       )
		#pragma comment( lib,RLib(Rblas.lib)   )
		#pragma comment( lib,RLib(Rlapack.lib) )
	#endif
#elif defined(CLANG_COMPILER)||defined(GCC_COMPILER)||defined(SOLARIS_COMPILER)
	#define  DllExport  
	#define _GNU_SOURCE
	#include <fenv.h>	
#endif
#if PTHREAD_INOUT==1
	#include pthreadPath(pthread.h) 
	#pragma comment(lib,pthreadLib(pthreadVC2.lib) )
#endif
#if MKL_LIBRARY==1
	#include "mkl.h"  
	#include "ipps.h"
	#pragma comment(lib,mklLib(mkl_intel_ilp64.lib) )
	#pragma comment(lib,mklLib(mkl_core.lib))
    #pragma comment(lib,mklLib(mkl_tbb_thread.lib) )
	#pragma comment(lib,tbbLib(tbb.lib))
	#pragma comment(lib,ippLib(ippcoremt.lib))
	#pragma comment(lib,ippLib(ippvmmt.lib))
	#pragma comment(lib,ippLib(ippsmt.lib))
#endif
#if (MYMAT_LIBRARY==1) && defined(MSVC_COMPILER)
	#pragma comment(lib,myMatLib(blas.lib))
#endif
#if  MATLAB_LIBRARY==1
	#include "blas.h" 
	#include "lapack.h"
	#pragma comment(lib,matlabLib(libmwblas.lib) )
	#pragma comment(lib,matlabLib(libmwlapack.lib) )
	#pragma comment(lib,myMatLib(blas.lib))
	#pragma comment(lib,FortranLib(libifcoremdd.lib))
#endif
#if M_INTERFACE==1
	#include  "mex.h"
	#pragma comment(lib,matlabLib(libmx.lib) )
	#pragma comment(lib,matlabLib(libmex.lib) )
	#pragma comment(lib,matlabLib(libmat.lib) )
#endif
#if R_INTERFACE==1
	#include <R.h>
	#include <Rinternals.h>
	#include <Rdefines.h>
	#include "Rmath.h"
	#if MYMAT_LIBRARY==1
	#endif
	#ifndef bool
		#define  bool  uint8_t
	#endif
	enum { false=0,true=1 };
	#define mwSize size_t
	#ifdef CLANG_COMPILER
	   #undef sign
	   #undef warning
	#endif
#endif
