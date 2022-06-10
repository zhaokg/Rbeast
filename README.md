

##  BEAST:  A Bayesian Ensemble Algorithm for Change-Point Detection and Time Series Decomposition

<img align="right"  height="500" src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/beach.png">

####  BEAST (Bayesian Estimator of Abrupt change, Seasonality, and Trend) is a fast, generic Bayesian model averaging algorithm to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations, as described in <ins>[Zhao et al. (2019)](https://go.osu.edu/beast2019)</ins>. BEAST is useful for changepoint detection (i.e., breakpoints or structural breaks), nonlinear trend analysis, time series decomposition, and time series segmentation
> **BEAST** was impemented in C/C++ but accessible from R and Matlab. Check the `Source`, `R`, and `Matlab` folders at [Github](https://github.com/zhaokg/Rbeast) for the C, R, and Matlab code.

**Quick installation**:
   * In Matlab, run **`eval(webread('http://b.link/beast',weboptions('cert','')))`**  
   * In R,     run **`install.packages("Rbeast")`**
   * Or follow the more detailed instructions below to install and run BEAST

## Installation for R

[![](https://www.r-pkg.org/badges/version/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast)
[![](http://cranlogs.r-pkg.org/badges/grand-total/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast)

In CRAN-Task-View: [[Time Series Analysis]](https://cran.r-project.org/web/views/TimeSeries.html#forecasting-and-univariate-modeling) and [[Bayesian inference]](https://cran.r-project.org/web/views/Bayesian.html#time-series-models)


   1. from **CRAN**: An R package **`Rbeast`** has been deposited at [CRAN](https://CRAN.R-project.org/package=Rbeast). (On CRAN, there is another Bayesian time-series package named "beast", which has nothing to do with the BEAST algorithim. Our package is `Rbeast`.) Install it in R using
      ```R
       install.packages("Rbeast")
      ```

   2. from **GitHub**: The latest R versions are also available here at [GitHub](https://github.com/zhaokg/Rbeast) and can be installed using

      ```R
        # Windows x86 or x64 (install from binary)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Windows/Rbeast_0.9.4.zip" ,repos=NULL)         
        # Linux/Mac (install from source)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Rbeast_0.9.4.tar.gz", repos = NULL, type="source")
      ```


   #### Run and test Rbeast in R
   
<img align="left"  height="330" src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/Nile.png">

The main functions in Rbeast are `beast(Y, ...)`, `beast.irreg(Y, ...)`, and `beast123(Y, metadata, prior, mcmc, extra)`. The code snippet below  provides a starting point for the basic usage.  

  ```R
      library(Rbeast)
      data(Nile)                       #  annual streamflow of the Nile River    
      out = beast(Nile, season='none') #  'none': trend-only data without seasonlaity   
      print(out)                   
      plot(out)
      ?Rbeast                          # See more details about the usage of `beast`    
      
      tetris()                         # if you dare to waste a few moments of your life 
      minesweeper()                    # if you dare to waste a few more moments of your life 
  ```


<br/>
<br/>
<br/>

----

## Installation for Matlab  [![View Rbeast on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72515-rbeast)
 

Install the Matlab version of BEAST automatically to a local folder of your choice by running 
  ```Matlab  
  beastPath = 'C:\beast\'                   
  eval( webread('http://b.link/beast') )  
  
  % NOTE:
  % 1. Write permission needed for your chosen path; the variable name must be 'beastPath'
  % 2. If webread has a certificate error, run the following line instead:
  %     eval(  webread( 'http://b.link/beast', weboptions('cert','') )  );
  ```
The above will download all the files in the [Rbeast\Matlab folder at Github](https://github.com/zhaokg/Rbeast) to the chosen folder: if `beastPath` is missing, a default temporary folder (e.g., `C:\Users\$user_name$\AppData\Local\Temp\Rbeast for Windows 10`) will be used. If the automatic script fails, please download the Matlab files from [Github](https://github.com/zhaokg/Rbeast) manually. These files include a Matlab mex library compiled from the C soure code (e.g., `Rbeast.mexw64` for Windows, `Rbeast.mexa64` for Linux, `Rbeast.mexmaci64` for MacOS) and some Matlab wrapper functions (e.g.,`beast.m`, and `beast123.m`) similar to the R interface, as well as some test datasets (e.g., Nile.mat, and co2.mat).

#### Usage 
The Matlab API is similar to those of R. Below is a quick example:
  ```Matlab
   help beast
   help beast123  
   load('Nile.mat')                                   % annual streamflow of the Nile River startin from year 1871
   out = beast(Nile, 'season', 'none','start', 1871)  % trend-only data without seasonality
   printbeast(out)
   plotbeast(out)
  ```
We generated the Matlab mex binary library on our own machines with Win10, Ubuntu 22.04, and macOS High Sierra. If they fail on your machine, the mex library can be compiled  from the C source code files under `Rbeast\Source`. If needed, we are happy to work with you to compile for your specific OS or machines. Additional information on compilations from the C source is also given below.

----
### Python/Julia

Wrappers in Python and Julia are being developed: We welcome contributions and help from interested developers. If interested, contact Kaiguang Zhao at zhao.1423@osu.edu.

## Description
Interpretation of time series data is affected by model choices. Different models can give different or even contradicting estimates of patterns, trends, and mechanisms for the same data–a limitation alleviated by the Bayesian estimator of abrupt change,seasonality, and trend (BEAST) of this package. BEAST seeks to improve time series decomposition by forgoing the "single-best-model" concept and embracing all competing models into the inference via a Bayesian model averaging scheme. It is a flexible tool to uncover abrupt changes (i.e., change-points), cyclic variations (e.g., seasonality), and nonlinear trends in time-series observations. BEAST not just tells when changes occur but also quantifies how likely the detected changes are true. It detects not just piecewise linear trends but also arbitrary nonlinear trends. BEAST is applicable to real-valued time series data of all kinds, be it for remote sensing, finance, public health, economics, climate sciences, ecology, and hydrology. Example applications include its use to identify regime shifts in ecological data, map forest disturbance and land degradation from satellite imagery, detect market trends in economic data, pinpoint anomaly and extreme events in climate data, and unravel system dynamics in biological data. Details on BEAST are reported in [Zhao et al. (2019)](https://go.osu.edu/beast2019). The paper is available at https://go.osu.edu/beast2019.

## Reference
>Zhao, K., Wulder, M. A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick B., Zhang, X., & Brown, M. (2019). [Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm.](https://go.osu.edu/beast2019) Remote Sensing of Environment, 232, 111181. (the BEAST paper) 
>
>Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. [Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection](https://www.academia.edu/download/55199778/Hyperspectral-biochemical-Bayesian-chlorophyll-carotenoid-LAI-water-content-foliar-pigment.pdf). Remote Sensing of Environment, 132, pp.102-119. (the mcmc sampler used for BEAST)
>
>Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. [Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine](https://pages.charlotte.edu/gang-chen/wp-content/uploads/sites/184/2021/05/Hu_2021_BEAST-HF-s.pdf). ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261. (an application paper)

----
## Selected publications using BEAST/Rbeast
| Discipline | Publication Title |
| --- | --- |
|  Population Ecology | *Henderson, P. A. (2021). Southwood's Ecological Methods (5th edition). Oxford University Press., page 475-476*|
| Finance|*Candelaria, Christopher A., Shelby M. McNeill, and Kenneth A. Shores. (2022). What is a School Finance Reform? Uncovering the ubiquity and diversity of school finance reforms using a Bayesian changepoint estimator.(EdWorkingPaper: 22-587). Retrieved from Annenberg Institute at Brown University: https://doi.org/10.26300/4vey-3w10*|
| Social Science|*Linnell, K., Fudolig, M., Schwartz, A., Ricketts, T.H., O'Neil-Dunne, J.P., Dodds, P.S. and Danforth, C.M., 2022. Spatial changes in park visitation at the onset of the pandemic. arXiv preprint arXiv:2205.15937.*|
| Hydrology | *Zohaib, M. and Choi, M., 2020. Satellite-based global-scale irrigation water use and its contemporary trends. Science of The Total Environment, 714, p.136719.* |
| Energy Engineering |*Lindig, S., Theristis, M. and Moser, D., 2022. Best practices for photovoltaic performance loss rate calculations. Progress in Energy, 4(2), p.022003.*|
|Virology|*Shen, L., Sun, M., Song, S., Hu, Q., Wang, N., Ou, G., Guo, Z., Du, J., Shao, Z., Bai, Y. and Liu, K., 2022. The impact of anti‐COVID‐19 nonpharmaceutical interventions on hand, foot, and mouth disease—A spatiotemporal perspective in Xi'an, northwestern China. Journal of medical virology.*|
| Pharmaceutical Sciences|*Patzkowski, M.S., Costantino, R.C., Kane, T.M., Nghiem, V.T., Kroma, R.B. and Highland, K.B., 2022. Military Health System Opioid, Tramadol, and Gabapentinoid Prescription Volumes Before and After a Defense Health Agency Policy Release. Clinical Drug Investigation, pp.1-8.*|
| Geography|*Cai, Y., Liu, S. and Lin, H., 2020. Monitoring the vegetation dynamics in the Dongting Lake Wetland from 2000 to 2019 using the BEAST algorithm based on dense Landsat time series. Applied Sciences, 10(12), p.4209.*|
| Oceanography|*Pitarch, J., Bellacicco, M., Marullo, S. and Van Der Woerd, H.J., 2021. Global maps of Forel–Ule index, hue angle and Secchi disk depth derived from 21 years of monthly ESA Ocean Colour Climate Change Initiative data. Earth System Science Data, 13(2), pp.481-490.*|
|Photovoltaics|*Micheli, L., Theristis, M., Livera, A., Stein, J.S., Georghiou, G.E., Muller, M., Almonacid, F. and Fernández, E.F., 2021. Improved PV soiling extraction through the detection of cleanings and change points. IEEE Journal of Photovoltaics, 11(2), pp.519-526.*|
|Climate Sciences|*White, J.H., Walsh, J.E. and Thoman Jr, R.L., 2021. Using Bayesian statistics to detect trends in Alaskan precipitation. International Journal of Climatology, 41(3), pp.2045-2059.*|
|Field Hydrology|*Merk, M., Goeppert, N. and Goldscheider, N., 2021. Deep desiccation of soils observed by long-term high-resolution measurements on a large inclined lysimeter. Hydrology and Earth System Sciences, 25(6), pp.3519-3538.*|
|Remote Sensing|*Mardian, J., Berg, A. and Daneshfar, B., 2021. Evaluating the temporal accuracy of grassland to cropland change detection using multitemporal image analysis. Remote Sensing of Environment, 255, p.112292.*|
|Forest Ecology|*Moreno-Fernández, D., Viana-Soto, A., Camarero, J.J., Zavala, M.A., Tijerín, J. and García, M., 2021. Using spectral indices as early warning signals of forest dieback: The case of drought-prone Pinus pinaster forests. Science of The Total Environment, 793, p.148578.*|
|Atmospheric Sciences|*Tingwei, C., Tingxuan, H., Bing, M., Fei, G., Yanfang, X., Rongjie, L., Yi, M. and Jie, Z., 2021. Spatiotemporal pattern of aerosol types over the Bohai and Yellow Seas observed by CALIOP. Infrared and Laser Engineering, 50(6), p.20211030.*|
|Terrestrial ecology|*Dashti, H., Pandit, K., Glenn, N.F., Shinneman, D.J., Flerchinger, G.N., Hudak, A.T., de Graaf, M.A., Flores, A., Ustin, S., Ilangakoon, N. and Fellows, A.W., 2021. Performance of the ecosystem demography model (EDv2. 2) in simulating gross primary production capacity and activity in a dryland study area. Agricultural and Forest Meteorology, 297, p.108270.*|
|Environmental Engineering|*Bainbridge, R., Lim, M., Dunning, S., Winter, M.G., Diaz-Moreno, A., Martin, J., Torun, H., Sparkes, B., Khan, M.W. and Jin, N., 2022. Detection and forecasting of shallow landslides: lessons from a natural laboratory. Geomatics, Natural Hazards and Risk, 13(1), pp.686-704.*|
|Hydrology|*Yang, X., Tian, S., You, W. and Jiang, Z., 2021. Reconstruction of continuous GRACE/GRACE-FO terrestrial water storage anomalies based on time series decomposition. Journal of Hydrology, 603, p.127018.*|
|Landscape Ecology|*Adams, B.T., Matthews, S.N., Iverson, L.R., Prasad, A.M., Peters, M.P. and Zhao, K., 2021. Spring phenological variability promoted by topography and vegetation assembly processes in a temperate forest landscape. Agricultural and Forest Meteorology, 308, p.108578.*|

## Additional Notes

1. **Computation**

As a Bayesian algorithm, BEAST is fast and is possibly among the fastest implementations of Bayesian time-series analysis algorithms of the same nature. (But it is slower, compared to nonBayesian methods.) For applications dealing with a few to thousands of time series, the computation won't be an practical concern at all. But for remote sensing applications that may easily involve millions or billions of time series, computation will be a big challenge for Desktop computer users. We suggest first testing BEAST on a single time series or small image chips first to determine whether BEAST is appropriate for your applications and, if yes, estimate how long it may take to process the whole image. In any case, for stacked time-series images, do not use `beast` or `beast.irreg` and use **`beast123`** instead, which can handle 3D data cubes and allow parallel computing.  We also welcome consultation with Kaiguang Zhao (zhao.1423@osu.edu) to give specific suggestions if you see some value of BEAST for your applications.

2. **Compilation from source code**

The BEAST source code appears more complicated than necessary, mainly because the same source is used for both R and Matlab interfaces (also for Python and Julia) as well as for various compilations settings (e.g., compiler variants, alternative library dependencies, cross-platform compatibility, mixed language interfaces, and Win32 API native interfaces). The complication control variables are defined as MACROs in abc_macro.h. Of the soure code files, there are dozens of "abc_xxxx.c" files, which are some auxiliary files; the BEAST algorithm itself is coded in beastv2_COREV4.c; and the R and Matlab interfaces are coded in glue.c and abc_ide_util.c.

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
 

## Reporting Bugs or getting help

BEAST is distributed as is and without warranty of suitability for application. The one distributed above is still a beta version, with potential room for further improvement. If you encounter flaws with the software (i.e. bugs) please report the issue. Providing a detailed description of the conditions under which the bug occurred will help to identify the bug, you can directly email its maintainer Dr. Kaiguang Zhao at zhao.1423@osu.edu. Alternatively, Use the [Issues tracker](https://github.com/zhaokg/Rbeast/issues) on GitHub to report issues with the software and to request feature enhancements. 
