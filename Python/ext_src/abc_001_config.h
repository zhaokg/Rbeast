#pragma once

#include "abc_000_macro.h"

//#include "intrin.h"     // _rdstc
//#include <stdio.h>	  // fprintf fopen FILE
//#include <stdlib.h>     // malloc and free
//include  <string.h>	  // memset memcpy
//#include <io.h>

// Dym Lib file externsion: 
// .dll for R
// .mexw64 for Matlab, 
// .pyd for Pyhton


#define R_INTERFACE    0
#define M_INTERFACE    0
#define P_INTERFACE    0
//#define  PLAY_MODE

//#define  R_RELEASE
//#define  P_RELEASE
//#define  M_RELEASE

#define MYMAT_LIBRARY   1
#define MKL_LIBRARY     0
#define MATLAB_LIBRARY  0 // For some uknown reason, the Matlab version of Blas and Lapack function doesn't work.

#define PCGRAND_LIBRARY 1
#define MKLRAND_LIBRARY 0

/**************************************************************************************************/
//  Define which sofware platform is for: R, Python or Matlab
/**************************************************************************************************/
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

#elif  defined(P_RELEASE)
        // For P interface
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

#elif defined(M_RELEASE)
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
 
#if MYMAT_LIBRARY == 1
		#define PCGRAND_LIBRARY 1
		#define MKLRAND_LIBRARY 0
#endif


#ifdef COMPILER_MSVC
	#if   R_INTERFACE==1
		#define DllExport   __declspec( dllexport ) 
	#elif P_INTERFACE==1
		#define DllExport   __declspec( dllexport ) 
	#elif M_INTERFACE==1
		#define DllExport 
		#pragma comment(linker, "/export:mexFunction")  //(Linker option) export:mexFunction   
    #endif
#elif defined(COMPILER_CLANG)|| defined(COMPILER_GCC) || defined(COMPILER_SOLARIS)   
   #define  DllExport   	   
#endif
  
/**************************************************************************************************/
//  Define the include and lib paths
/**************************************************************************************************/
  
#define QUOTE_IT(x) #x

#ifdef COMPILER_MSVC

	#define INCLUDE_PTHREAD(_X_) QUOTE_IT(C:/USERS/zhaokg/Documents/Visual Studio 2013/Projects/Matlab_Mex_Test/Pthread_IncludeLib/_X_ )
	#define INCLUDE_MATLAB(_X_)  QUOTE_IT(C:/Program Files/MATLAB/R2019a/extern/include/_X_ )
	#define INCLUDE_MKL(_X_)     QUOTE_IT(C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/_X_ )
	#define INCLUDE_IPP(_X_)     QUOTE_IT(C:/Program Files (x86)/Intel/oneAPI/ipp/latest/include/_X_ )
	#define INCLUDE_R(_X_)       QUOTE_IT(C:/Program Files/R/R-4.2.2/include/_X_ )

    #define LIB_FORTRAN(_X_)     QUOTE_IT(C:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/compiler/lib/intel64_win/_X_ )
	#define LIB_PTHREAD(_X_)     QUOTE_IT(C:/USERS/zhaokg/Documents/Visual Studio 2013/Projects/Matlab_Mex_Test/Pthread_IncludeLib/_X_ )
	#define LIB_MATLAB(_X_)      QUOTE_IT(C:/Program Files/MATLAB/R2019a/extern/lib/win64/microsoft/_X_ )
	#define LIB_TBB(_X_)         QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/tbb/lib/intel64_win/vc_mt/_X_ )
	//#define LIB_TBBstatic(_X_) QUOTE_IT(G:/Intel_Library/tbb-2019_U3/build/vs2013/x64/Release-MT/_X_ )
	#define LIB_MKL(_X_)         QUOTE_IT(C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/intel64/_X_ )
	#define LIB_IPP(_X_)         QUOTE_IT(C:/Program Files (x86)/Intel/oneAPI/ipp/latest/lib/intel64/_X_ )
    #define LIB_OpenMP(_X_)      QUOTE_IT(C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64_win/_X_ )
	#define LIB_MyLIB(_X_)       QUOTE_IT(D:/Share/Fortran_blas_lib/_X_ )
    #define LIB_Python(_X_)      QUOTE_IT(C:/Anaconda3/libs/_X_ )
   
	#ifdef TARGET_64
        #define LIB_R(_X_)   QUOTE_IT(C:/Program Files/R/R-4.5.2/implib/_X_ )
	#else
        #define LIB_R(_X_)   QUOTE_IT(C:/Program Files/R/R-4.2.2/implib/i386/_X_ )				   
	#endif
  
#endif

#define PTHREAD_INOUT 0
#if PTHREAD_INOUT == 1
    //#include "pthread.h"
	#include INCLUDE_PTHREAD(pthread.h) 
	#pragma comment(lib , LIB_PTHREAD(pthreadVC2.lib) )
#endif

#if MKL_LIBRARY == 1
	//#include "mkl.h"   //#include "mkl_vsl.h"  
	//#include "ipps.h"

	#include INCLUDE_IPP(ipps.h)
	#include INCLUDE_MKL(mkl.h)
/*
	#pragma comment(lib , LIB_MKL(mkl_intel_ilp64.lib) )
	#pragma comment(lib , LIB_MKL(mkl_core.lib))
    #pragma comment(lib , LIB_MKL(mkl_tbb_thread.lib) )
	#pragma comment(lib , LIB_TBB(tbb.lib))
*/	

/*
    #pragma comment(lib , LIB_TBBstatic(tbb.lib))
    #pragma comment(lib , LIB_TBBstatic(tbbmalloc.lib))
    #pragma comment(lib , LIB_TBBstatic(tbbmalloc_proxy.lib))
	#pragma comment(lib , LIB_MKL(mkl_intel_thread.lib))
	#pragma comment(lib , LIB_OpenMP(libiomp5md.lib))
*/	

    // Non-threaded sequential library 
	#pragma comment(lib , LIB_MKL(mkl_intel_ilp64.lib) )
	#pragma comment(lib , LIB_MKL(mkl_core.lib))
	#pragma comment(lib , LIB_MKL(mkl_sequential.lib) )

	#pragma comment(lib , LIB_IPP(ippcoremt.lib))
	#pragma comment(lib , LIB_IPP(ippvmmt.lib))  // Vector math
	#pragma comment(lib , LIB_IPP(ippsmt.lib))   // Signal processing 

#endif


#if (MYMAT_LIBRARY ==1) && defined(COMPILER_MSVC) && 0
	#pragma comment(lib , LIB_MyLIB(blas_oneAPI.lib))

     /* 
	 Liblifcormdd.lib is an import lib for the Multithread runtime DLL. If used, 
	 libifcoremdd.dll must be also supplised. To use the static version instead,
	 compile blas.lib by setting 
	  (1) "Fortran->Libary->Multithreaded" 
	  (2) "Librian->Additional depencies->libifcoremt.lib" 
	  (3) Librian->Link Libray dependencies (yes). 
	 This way, "libifcoremt.lib" will be inserted into blas.lib. 
	 */
	 
	//#pragma comment(lib , LIB_FORTRAN(libifcoremdd.lib))

    // The libs below are needed for linking the static lib "blas.lib"; they
    // are Intel Fortran's runtime-like librareies. These should be all static lib
    // (e.g., libifcoremt: multithreaded static lib).
    #pragma comment(lib , LIB_FORTRAN(ifconsol.lib))
    #pragma comment(lib , LIB_FORTRAN(libifcoremt.lib))
	#pragma comment(lib , LIB_FORTRAN(libifport.lib))
	#pragma comment(lib , LIB_FORTRAN(libirc.lib))
	#pragma comment(lib , LIB_FORTRAN(libircmt.lib)) //Needed for the oneAPI version of Fortran
	#pragma comment(lib , LIB_FORTRAN(libmmt.lib))	
    #pragma comment(lib , LIB_FORTRAN(svml_dispmt.lib))	
	#pragma comment(lib , LIB_FORTRAN(ifmodintr.lib))	
#endif

#if  MATLAB_LIBRARY == 1
	#include "blas.h" 
	#include "lapack.h"
	#pragma comment(lib , LIB_MATLAB(libmwblas.lib) )
	#pragma comment(lib , LIB_MATLAB(libmwlapack.lib) )
	#pragma comment(lib , LIB_MyLIB(blas.lib))
	#pragma comment(lib , LIB_FORTRAN(libifcoremdd.lib))
#endif


/**************************************************************************************************/
//  the R statistical software
/**************************************************************************************************/

#if R_INTERFACE == 1
    
    #ifdef ERROR 
       #undef ERROR  // ERROR is defined in both windwows.h/wingdi.h and R.h
    #endif

    //https://stackoverflow.com/questions/69672625/how-to-disable-r-no-remap-when-calling-rs-c-rinterals-h 
    //#define R_NO_REMAP  // mapping r api to aliaes (e.g., sign <> Rf_sign)   
	
	#include <R.h>
	#include <Rinternals.h>
	#include <Rdefines.h>
	#include <Rmath.h>
	#include <R_ext/GraphicsEngine.h>
	#include <R_ext/GraphicsDevice.h>	
		//#include <R_ext\Rdynload.h>
	#if MYMAT_LIBRARY ==1
		//#include "R_ext\Lapack.h"
		//#include "R_ext\BLAS.h"
	#endif


	// This must appear after the inclusion of R.h because without the "R_NO_REMAP" macro, sign and warning will be mapped
	#undef sign    // In R 'sign'   is  defined as Rf_sign if R_IO_REMAP"
	#undef warning // In R 'warning' is defined as Rf_warning; warning is also a C++ keyword

	#define mwSize size_t
		
	#ifdef COMPILER_MSVC
        // These R import libs are generated by
        //  >  gendef  Rlapack.dll
        //  >  lib /def:Rlapack.def /machine:x64 /out:Rlapack.lib
		#pragma comment( lib , LIB_R(R.lib)          )
		#pragma comment( lib , LIB_R(Rblas.lib)      )
		#pragma comment( lib , LIB_R(Rlapack.lib)    )
	#endif
	
#endif

/**************************************************************************************************/
//  Matlab
/**************************************************************************************************/
#if M_INTERFACE == 1

	#define    MEX_DOUBLE_HANDLE
	#include  "mex.h"
	
	#ifdef COMPILER_MSVC	
		#pragma comment(lib , LIB_MATLAB(libmx.lib) )
		#pragma comment(lib , LIB_MATLAB(libmex.lib) )
		#pragma comment(lib , LIB_MATLAB(libmat.lib) )

		// See: stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running/26271557
		//  - libmwservices for ioFlush (this is a C++ function with naming managling)
		// See: www.advanpix.com/2016/07/02/devnotes-3-proper-handling-of-ctrl-c-in-mex-module/
		//  - libut.lib     for utIsInterruptPending and utSetInterruptPending(Bool)
		
		#pragma comment(lib,  LIB_MATLAB(libmwservices.lib) )	 
		#pragma comment(lib , LIB_MATLAB(libut.lib) ) 
   #endif
	

/*  
    In Matlab,  mex can't handle the mixed inputs of C and CPP. Check these:
    https://stackoverflow.com/questions/57579194/matlab-r2016b-mex-fails-to-compile-c-code
    https://www.mathworks.com/matlabcentral/answers/827145-how-can-i-split-c-c-code-in-multiple-files-and-still-use-mex
    mex accept a wildcard for source files (*.c) but not a wildcard for object files (*.o)

    mex -v -c CFLAGS='$CFLAGS -DM_RELEASE     -Wall -v -Wl,-v'  *.c     // compiling all the C source
    mex -v -c CXXFLAGS='$CXXFLAGS -DM_RELEASE -Wall -v -Wl,-v'  *.cpp   // compiling all the C source
    mex -v    CFLAGS='$CFLAGS  -Wall -Wl,-v' -lmwservices -lut abc_ioFlush.o *.o -output Rbeast  // for some reason, the output is bad

    Below is my solution:
    > system("gcc -c -fPIC  -pthread -DNDEBUG -DM_RELEASE   -DMATLAB_DEFAULT_RELEASE=R2017b  -DMATLAB_MEX_FILE  -I/MATLAB/extern/include   -O2 -Wall  -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  *.c")
    > system("g++ -c -fPIC  -pthread -DNDEBUG -DM_RELEASE   -DMATLAB_DEFAULT_RELEASE=R2017b  -I/MATLAB/extern/include   -O2 -Wall    -mfpmath=sse -msse2 -mstackrealign  *.cpp")
    > system("gcc -shared -pthread  -L/MATLAB/bin/glnxa64 -lmx -lmex -lmat -lm  -lut -lmwservices *.o -o Rbeast.mexa64")
*/	
	
#endif


/**************************************************************************************************/
//  Python
/**************************************************************************************************/

#if P_INTERFACE == 1

// num[y 1.x:  C:\Anaconda3\pkgs\numpy-base-1.16.5-py37hc3f5095_0\Lib\site-packages\numpy\core\include
// num[y 2.x:  C:\Anaconda3\envs\py311\Lib\site-packages\numpy\core\include

// Python3.8:   C:\Anaconda3\include
// Python3.11:  C:\Anaconda3\envs\py311\include

//Python 3.8 lib:  C:\Anaconda3\pkgs\python-3.8.8-hdbf39b2_5\libs
//Python 3.11 lib:  C:\Anaconda3\envs\py311\libs

	#include "Python.h"                 // inside, it checks/defines "_GNU_SOURCE" 
	#include "structmember.h"

// https:// stuff.mit.edu/afs/sipb/project/python/src/python-numeric-22.0/doc/www.pfdubois.com/numpy/html2/numpy-13.html
// * If not defining it, PyArray_API will have local copies in each file.
// https:// stackoverflow.com/questions/32899621/numpy-capi-error-with-import-array-when-compiling-multiple-modules
// https:// docs.scipy.org/doc/numpy-1.10.1/reference/c-api.array.html#miscellaneous

   #define  NPY_NO_DEPRECATED_API   NPY_1_7_API_VERSION  //suppress the warning "Using deprecated NumPy API, disable it with"

/*
	Suppose I have two files coolmodule.c and coolhelper.c which need to be compiled and linked into a
	single extension module. Suppose coolmodule.c contains the required initcool module initialization 
	function (with the import_array() function called). Then, coolmodule.c would have at the top:

	#define PY_ARRAY_UNIQUE_SYMBOL cool_ARRAY_API
	#include "numpy/arrayobject.h"
	
	On the other hand, coolhelper.c would contain at the top:

	#define NO_IMPORT_ARRAY
	#define PY_ARRAY_UNIQUE_SYMBOL cool_ARRAY_API
	#include "numpy/arrayobject.h"
	
	You can also put the common two last lines into an extension-local header file as long as you make
	sure that NO_IMPORT_ARRAY is #defined before #including that file.
*/

// We will not use the standard method any longer and instead use the explicit loading

	#ifdef  USE_STANDARD_METHOD_IMPORT_NUMPY
		#ifdef  SHOULD_IMPORT_NUMPY
		   #define PY_ARRAY_UNIQUE_SYMBOL NumpyAPIList
		   #include "numpy/arrayobject.h"
		#else
			#define NO_IMPORT_ARRAY
			#define PY_ARRAY_UNIQUE_SYMBOL NumpyAPIList
			#include "numpy/arrayobject.h"
		#endif	   
	#endif


/*
* How to import Nummpy and set up the C API function list
* (1) Numpy 1.x: Call import_arraY()
* (2) Num[y 2.x  Call PyArray_ImportNumPyAPI(), which will indirectly call _import_array()
 PyArray_ImportNumPyAPI()
    -> if (PyArray_API == NULL)
           import_array1(-1)
               -> _import_array()
                    -> import NumPy module
                    -> fetch _ARRAY_API capsule
                    -> set PyArray_API
                    -> check ABI/C-API versions
*/

  // #pragma comment(lib , LIB_Python(python37_msvc.lib) ) 
  // No need to set up the lib explicitly 
  // (1) pynconfig.h has a pragma line to add python.lib:  pragma comment(lib,"python37.lib")   
  // (2) Numpy uses the import_array to load the (np.core._multiarray_umath._ARRAY_API) variable to fill the api arrays

// Improt the pyd extension module
// import sys
// sys.path.append('y:/testold')
// import Rbeast

#endif


/**************************************************************************************************/
//  Define Bool
//  * bool is defined in matlab's tmwtypes.h
//  * bool is a keyword in C24
/**************************************************************************************************/
#if M_INTERFACE == 1
	#ifndef Bool
		 #include "mex.h"
		 #define  Bool  bool          // typedef uint8_t Bool;
    #endif
#elif R_INTERFACE==1 || P_INTERFACE==1
	#ifndef Bool		
		#define  Bool  unsigned char  // typedef uint8_t Bool;
	#endif
#endif

 

 