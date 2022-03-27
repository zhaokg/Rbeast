#  BEAST:  A Bayesian Ensemble Algorithm for Change-Point Detection and Time Series Decomposition
###  BEAST is a Bayesian model averaging algorithm to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations, as described in [Zhao et al. (2019)](https://go.osu.edu/beast2019). BEAST is useful for changepoint detection (i.e., breakpoints or structural breaks), nonlinear trend analysis, time series decomposition, and time series segmentation

> **BEAST** was impemented in C/C++. Check the "Source" folder for the source code. R and Matlab interfaces are also provided and can be found under the "R" and "Matlab" folders above. Or follow the instructions below to install and run BEAST in R or Matlab.

# Installation
## R


   1. from **CRAN**: An R package **`Rbeast`** has been deposited at [CRAN](https://CRAN.R-project.org/package=Rbeast). (On CRAN, there is another Bayesian time-series package named "beast", which has nothing to do with the BEAST algorithim. Our package name is `Rbeast`.) Install `Rbeast` in your R console using
      ```R
       install.packages("Rbeast")
      ```

   2. from **GitHub**: The latest versions of package files for **Rbeast** are available here at [GitHub](https://github.com/zhaokg/Rbeast). Alternative ways to install Rbeast are:

      ```R
        # Windows x86 or x64 (install from binary)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Windows/Rbeast_0.9.4.zip" ,repos=NULL)
         
        # Linux/Mac (install from source)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Rbeast_0.9.4.tar.gz", repos = NULL, type="source")
      ```
  

#### Run and test Rbeast

The main functions in Rbeast are `beast(Y, ...)`, `beast.irreg(Y, ...)`, and `beast123(Y, metadata, prior, mcmc, extra)`. The code snippet below  provides a starting point for the basic usage.

  ```R
      library(Rbeast)
      data(Nile)                       #  annual streamflow of the Nile River    
      out = beast(Nile, season='none') #  'none': trend-only data without seasonlaity   
      print(out)                   
      plot(out)
      ?Rbeast          # See more details about the usage of `beast`                 
 ```
---- 
## Matlab

#### Installation and usage (Windows and Linux)

The C code of BEAST has also been compiled into a Matlab mex library (i.e., Rbeast.mexw64 for Windows and Rbeast.mexa64 for Linux) with some wrapper matlab functions (e.g.,beast.m) similar to the R interface, all available at the Rbeast\Matlab folder above. Download the files to your local drive and run BEAST. Alternatively, run the following Matlab code to automatically download the files and install BEAST to your local drive:

  ```Matlab
  % Installation path of your choice; Write permission needed; the var name has to be 'beastPath'
  beastPath = 'C:\rbeast\';                    
  eval( webread('http://go.osu.edu/rbeast', weboptions('cert','')) );
  ```

The Matlab API is similar to those of R. Below is a quick example:
  ```Matlab
   help beast
   help beast123  
   load('Nile.mat')                                   % annual streamflow of the Nile River startin from year 1871
   out = beast(Nile, 'season', 'none','start', 1871)  % trend-only data without seasonality
   printbeast(out)
   plotbeast(out)
  ```
We generated the Matlab mex binary library only for Windows and Linux only. Mex libraries for other OS systems such as Mac can be compiled from the source code files under "\Rbeast\Source". If needed, we are happy to work with you to compile for your specific OS or machines. Additional informaiton on compliation from the C source is also given below.

---- 
## Python

A wrapper in Python is being developed: We wecolme contributions and help from interested developers. If interested, contact Kaiguang Zhao at zhao.1423@osu.edu.

---- 
## Julia

A wrapper in Julia is also being developed: We wecolme contributions and help from interested developers. If interested, contact Kaiguang Zhao at zhao.1423@osu.edu

## Description
Interpretation of time series data is affected by model choices. Different models can give different or even contradicting estimates of patterns, trends, and mechanisms for the same data–a limitation alleviated by the Bayesian estimator of abrupt change,seasonality, and trend (BEAST) of this package. BEAST seeks to improve time series decomposition by forgoing the "single-best-model" concept and embracing all competing models into the inference via a Bayesian model averaging scheme. It is a flexible tool to uncover abrupt changes (i.e., change-points), cyclic variations (e.g., seasonality), and nonlinear trends in time-series observations. BEAST not just tells when changes occur but also quantifies how likely the detected changes are true. It detects not just piecewise linear trends but also arbitrary nonlinear trends. BEAST is applicable to real-valued time series data of all kinds, be it for remote sensing, economics, climate sciences, ecology, and hydrology. Example applications include its use to identify regime shifts in ecological data, map forest disturbance and land degradation from satellite imagery, detect market trends in economic data, pinpoint anomaly and extreme events in climate data, and unravel system dynamics in biological data. Details on BEAST are reported in [Zhao et al. (2019)](https://go.osu.edu/beast2019). The paper is available at https://go.osu.edu/beast2019.

## Reference
>Zhao, K., Wulder, M. A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick B., Zhang, X., & Brown, M. (2019). [Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm.](https://go.osu.edu/beast2019) Remote Sensing of Environment, 232, 111181. (the BEAST paper) 
>
>Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. [Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection](https://www.academia.edu/download/55199778/Hyperspectral-biochemical-Bayesian-chlorophyll-carotenoid-LAI-water-content-foliar-pigment.pdf). Remote Sensing of Environment, 132, pp.102-119. (the mcmc sampler used for BEAST)
>
>Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. [Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine](https://pages.charlotte.edu/gang-chen/wp-content/uploads/sites/184/2021/05/Hu_2021_BEAST-HF-s.pdf). ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261.

---- 
## Additonal Notes

1. **Computation**

As a Bayesian algorithm, BEAST is fast and is possibly among the fastest implementations of Bayesian time-series analysis algorithms of the same nature. For applications dealing with a few to thousands of time series, the computation won't be an practical concern at all. But for remote sensing applications that may easily involve millions or billions of time series, compuation will be a big challenge for Desktop computer users. We suggest first testing BEAST on a single time series or small image chips first to determine whether BEAST is appropriate for your applications and, if yes, estimate how long it may take to process the whole image. In any case, for stacked time-series images, do not use `beast` or `beast.irreg` and use **`beast123`** instead, which can hanlde 3D data cubes and allow parallel computing.  We also welcome consulation with Kaiguang Zhao (zhao.1423@osu.edu) to give specific suggestions if you see some value of BEAST for your applications.

2. **Complifation from source code**

The BEAST source code appears more complicated than neccessary, mainly because the same source is used for both R and Matlab interfaces (also for Python and Julia) as well as for various different compliation settings (e.g., compiler variants, alternative library dependencies, cross-platform compatibility, mixed language interfaces, and Win32 API native interfaces). The compiliation control variables are defined as MARCOs in abc_marco.h. Of the soure code files, there are dozens of "abc_xxxx.c" files, which are some auxliary files; the BEAST aglrotihim itself is coded in beastv2_COREV4.c; and the R and Matlab inferfaces are coded in glue.c and abc_ide_util.c.

We tested our source code under many common compliers (e.g., MSVC, gcc, clang, and Oracle Developer Studio) and all succesfully passed (e.g., see the [Rbeast package status report](https://cran.r-project.org/web/checks/check_results_Rbeast.html)). To complie for R, you need to make sure your machine has a C compiler appropriately set up. For example, see [Package Development Prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) for the tools needed for your operating system. In particular, on Windows platforms, the most convenient option is to go with the Rtools toolkit. To complile for Matlab, the appropriate C/C++ header files (e.g., mex.h) have to be correctly specified. Below are some compliation schemes using the gnu compilers as an example.

* To create the Matlab libray, run the following steps.

     1. Download all the files in the Source folder to your local folder
     2. Go to your local folder and set it as the current dictory. Complie the C/C++ sources into object files

        ```C
        gcc -c -fPIC -pthread -DM_RELEASE -DMATLAB_DEFAULT_RELEASE=R2017b -I/MATLAB/extern/include -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  *.c
        g++ -c -fPIC -pthread -DM_RELEASE -DMATLAB_DEFAULT_RELEASE=R2017b -I/MATLAB/extern/include -O2 -Wall            -mfpmath=sse -msse2 -mstackrealign  *.cpp
        ```
        > "**/MATLAB/extern/include**" is the Matlab's path for the include header files such as mex.h. Replace it with the correct one for your machine. On Windows, the path is typically "**C:/Program Files/MATLAB/R2019a/extern/include**".

    3. Link all the objects into the Matlab mex library
     
        `gcc -shared -L/MATLAB/bin/glnxa64  -lpthread -lmx -lmex -lmat -lm -lut -lmwservices *.o -o Rbeast.mexa64`
        > "**/MATLAB/bin/glnxa64**" is the Matlab's path for the static/import libraries such as libmex.lib and libmat.lib. Replace it with the correct one for your machine. On Windows, the path is typically "**C:\Program Files\MATLAB\R2019a\extern\lib\win64\microsoft**" for the Visual studio compiler and "**C:\Program Files\MATLAB\R2019a\extern\lib\win64\mingw64**" for the MinGW gcc compiler. Also, for Windows, the output should be `Rbeast.mexw64`.
        
    4. Put the Rbeast.mex library together with other m scripts (e.g., beast.m) to call Rbeast via beast or beast123; alternatively Rbeast.mex can be called directly as follows:
       `Rbeast('beastv4',Y,metadata, prior,mcmc, extra)`
       
       
* To create the R dynamic libray (which is part of the R pacakge but not the whole R pacakge itself), run the following steps.

     1. Download all the files in the Source folder to your local folder
     2. Go to your local folder and set it as the current dictory. Complie the C/C++ sources into object files

        ```C
        gcc -c -fPIC -pthread -DR_RELEASE  -I/opt/R-devel/lib/R/include -O2 -Wall -std=gnu99 -mfpmath=sse -msse2 -mstackrealign  *.c
        g++ -c -fPIC -pthread -DR_RELEASE  -I/opt/R-devel/lib/R/include -O2 -Wall            -mfpmath=sse -msse2 -mstackrealign  *.cpp
        ```
        > "**/opt/R-devel/lib/R/include**" is the R's path for the include header files such as R.h. Replace it with the correct one for your machine. On Windows, the path is typically "**C:\Program Files\R\R-4.1.0\include**".

    3. Link all the objects into the R dll library
     
        `gcc -shared  -L/opt/R-devel/lib/R/lib  -lpthread -lm -lR *.o -o Rbeast.dll`
        > "**/opt/R-devel/lib/R/lib**" is the R's path for the static/import libraries such as libR.so. Replace it with the correct one for your machine. On Windows, the path is typically "**C:\Program Files\R\R-4.1.0\bin\x64**". On Windows, MinGW compilers should be able to link with the R.dll file directly, but for MSVC, R.dll has to be first exported as an import library to be linked.
        
    4. The R dll libray is useful ONLY if you intend to call the dll directly, as shown below.
 
      ```   
      dyn.load("Rbeast.dll");  
      o =.Call('rexFunction',list('beastv4', co2, metadata=metadata,prior,mcmc,extra,1 ),12345 );
      dyn.unload("Rbeast.dll")
      ```
 
 

## Reporting Bugs

BEAST is distributed as is and without warranty of suitability for application. The one distribubuted above is a beta version, with potentail room for further improvement. If you encounter flaws with the software (i.e. bugs) please report the issue. Providing a detailed description of the conditions under which the bug occurred will help to identify the bug. *Use the [Issues tracker](https://github.com/zhaokg/Rbeast/issues) on GitHub to report issues with the software and to request feature enchancements. Alternatively, you can directly email its maintainer Dr. Kaiguang Zhao at zhao.1423@osu.edu
