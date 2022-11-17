## Compilation BEAST from source code (for developers and experts only)

The BEAST source code appears more complicated than necessary, mainly because the same source is used for R, Matlab, and Python interfaces (also for Julia and IDL) as well as for various compilations settings (e.g., compiler variants, alternative library dependencies, cross-platform compatibility, CPU instruction sets, mixed language interfaces, and Win32 API native interfaces). The complication control variables are defined as MACROs in abc_macro.h. Of the soure code files, there are dozens of "abc_xxxx.c" files, which are some auxiliary files; the BEAST algorithm itself is coded in beastv2_COREV4.c; and the R and Matlab interfaces are coded in glue.c and abc_ide_util.c.

We tested our source code under many common compliers (e.g., MSVC, gcc, clang, icc, mingw gcc, and Solaris) and all successfully passed (e.g., see the [Rbeast package status report](https://cran.r-project.org/web/checks/check_results_Rbeast.html)). To compile for R, you need to make sure your machine has a C compiler appropriately set up. For example, see [Package Development Prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) for the tools needed for your operating system. In particular, on Windows platforms, the most convenient option is to go with the Rtools toolkit. To compile for Matlab, the appropriate C/C++ header files (e.g., mex.h) have to be correctly specified. Below are some compilation schemes using the gnu compilers as an example.

To compile from the source, first download all the C/C++ files in the Source folder to your local folder, and go to your local folder and make it as the current working directory.

* To create the Matlab library, run the following steps.
     1.  Compile the C/C++ sources into object files and link the object file as a mex lib
     ```C
     gcc -shared -fPIC -pthread -DM_RELEASE -I/MATLAB/extern/include -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign -L/MATLAB/bin/glnxa64  -lpthread -lmx -lmex -lmat -lm -lut -lmwservices  *.c -o Rbeast.mexa64
     ```
     > `/MATLAB/extern/include` is the Matlab's path for the include header files such as mex.h. Replace it with the correct one for your machine. On Windows, the path is typically `"C:/Program Files/MATLAB/R2019a/extern/include"`.
     > `/MATLAB/bin/glnxa64` is the Matlab's path for the static/import libraries such as libmex.lib and libmat.lib. Replace it with the correct one for your machine. On Windows, the path is typically `C:\Program Files\MATLAB\R2019a\extern\lib\win64\microsoft` for the Visual studio compiler and `C:\Program Files\MATLAB\R2019a\extern\lib\win64\mingw64` for the MinGW gcc compiler. Also, for Windows, the output should be `Rbeast.mexw64`.
     2.  Alternatively, if your Matlab has the mex command correctly set up, the mex library can be compiled from
     ```C
     mex -v CFLAGS='-DM_RELEASE -UUSE_MEX_CMD -fPIC -O2 -Wall -std=gnu99 -march=native' -lmwservices -lut *.c -output Rbeast.mexa64
     ```
     3. Put the resulting Rbeast.mex library together with other m scripts (e.g., beast.m) to call Rbeast via beast or beast123; if needed, Rbeast.mex can be called directly as follows:
     `Rbeast('beastv4',Y,metadata, prior,mcmc, extra)`      
    
* An R dynamic lib (which is the dll/so/dynliab file--part of the R package but not the whole R package itself) probably never needs to be created mannually. But just in case that it is needed, run the following steps.

     1. Compile the C/C++ sources into object files and link them as a shared lib

        ```C
        gcc  -shared  -fPIC -pthread -DR_RELEASE  -I/opt/R-devel/lib/R/include -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  -L/opt/R-devel/lib/R/lib -lpthread -lm -lR *.c -o Rbeast.dll
        ```
        > `/opt/R-devel/lib/R/include` is the R's path for the include header files such as R.h. Replace it with the correct one for your machine. On Windows, the path is typically `C:\Program Files\R\R-4.1.0\include`.
        
        > `/opt/R-devel/lib/R/lib` is the R's path for the static/import libraries such as libR.so. Replace it with the correct one for your machine. On Windows, the path is typically `C:\Program Files\R\R-4.1.0\bin\x64`. On Windows, MinGW compilers should be able to link with the R.dll file directly, but for MSVC, R.dll has to be first exported as an import library to be linked.
        
    2. The R dll library is useful ONLY if you intend to call the dll directly, as shown below.
    ```
      dyn.load("Rbeast.dll");  
      o =.Call('rexFunction',list('beastv4', co2, metadata=metadata,prior,mcmc,extra,1 ),12345 );
      dyn.unload("Rbeast.dll")
    ```
 
