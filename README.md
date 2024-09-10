


##    BEAST:  A Bayesian Ensemble Algorithm for Change-Point Detection and Time Series Decomposition

**BEAST (Bayesian Estimator of Abrupt change, Seasonality, and Trend) is a fast, generic Bayesian model averaging algorithm to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations, as described in <ins>[Zhao et al. (2019)](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view)</ins>. BEAST is useful for *changepoint detection (e.g., breakpoints, structural breaks,joinpoints, regime shifts, or anomalies), trend analysis, time series decomposition (e.g., trend vs seasonality), time series segmentation, and interrupted time series analysis*. See a list of <a href="#publicationid"> selected studies using BEAST </a>.**


**Quick Installation**
> BEAST was impemented in C/C++ but accessible from  R, Python, Matlab, and Octave.  Install it as follows:

* [Python](#python):   **`pip install Rbeast`**   
* [Matlab](#matlab):  **`eval(webread('http://b.link/rbeast',weboptions('cert','')))`**
* [Octave](#octave):  **`eval(webread('http://b.link/rbeast'))`**  
* [R lang](#r):  **`install.packages("Rbeast")`** 


   
**Quick Usage**

> One-liner code for Python, Matlab and R:
```
# Python example
import Rbeast as rb; (Nile, Year)=rb.load_example('nile'); o=rb.beast(Nile,season='none'); rb.plot(o)

# Matlab/Octave example
load('Nile'); o = beast(Nile, 'season','none'); plotbeast(o)

# R example
library(Rbeast); data(Nile); o = beast(Nile); plot(o)
```



## Installation for R <a name=r> </a>

 <!---
 =====================================================================================
 =====================================================================================
 =====================================================================================
 This is our comments:

<table border="0"  style='border:none;'  bordercolor="#ffffff"  > <tr style='border:none;'  >   
   <td valign="center" style='border:none;'   > 
           <img src="https://www.r-pkg.org/badges/version/Rbeast?color=green">  
           <img src="http://cranlogs.r-pkg.org/badges/grand-total/Rbeast?color=green">
   </td>
   <td valign="center" style='border:none;'  > 
   In CRAN-Task-View: 
   <a href="https://cran.r-project.org/web/views/TimeSeries.html#forecasting-and-univariate-modeling">[Time Series Analysis]</a> 
   <a href="https://cran.r-project.org/web/views/Bayesian.html#time-series-models">[Bayesian inference]</a> 
   <a href="https://cran.r-project.org/web/views/Environmetrics.html#environmental-time-series">[Environmetrics]</a> 
    </td>  
  
 </tr></table>
 
  [![](https://www.r-pkg.org/badges/version/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast)  [![](http://cranlogs.r-pkg.org/badges/grand-total/Rbeast?color=green)](https://cran.r-project.org/package=Rbeast) 

 In CRAN-Task-View: [[**Time Series Analysis**]](https://cran.r-project.org/web/views/TimeSeries.html#forecasting-and-univariate-modeling), [[**Bayesian inference**]](https://cran.r-project.org/web/views/Bayesian.html#time-series-models), and [[**Environmetrics**]](https://cran.r-project.org/web/views/Environmetrics.html#environmental-time-series)  

from **GitHub**: The latest R versions are also available here at [GitHub](https://github.com/zhaokg/Rbeast) and can be installed using

      ```R
        # Windows x86 or x64 (install from binary)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Windows/Rbeast_0.9.5.zip" ,repos=NULL)         
        # Linux/Mac (install from source)
        install.packages("https://github.com/zhaokg/Rbeast/raw/master/R/Rbeast_0.9.5.tar.gz", repos = NULL, type="source")
      ```
=====================================================================
=====================================================================================================
=====================================================================================

 -->
 
 
<p  align="left">   
<a href= "https://cran.r-project.org/package=Rbeast">
<img src="https://www.r-pkg.org/badges/version/Rbeast?color=green"          height="20">  
<img src="https://cranlogs.r-pkg.org/badges/grand-total/Rbeast?color=green" height="20">  
<img src="https://img.shields.io/static/v1?style=plastic&logo=r&label=Rbeast%20%20&message= In CRAN-Task-View&color=brightgreen" height="20">
</a>
</p>   
  
> **Rbeast in CRAN-TASK-VIEW**:   <a href="https://cran.r-project.org/web/views/TimeSeries.html#forecasting-and-univariate-modeling">[Time Series Analysis]</a>    <a href="https://cran.r-project.org/web/views/Bayesian.html#time-series-models">[Bayesian inference]</a>     <a href="https://cran.r-project.org/web/views/Environmetrics.html#environmental-time-series">[Environmetrics]</a>  

An R package **`Rbeast`** has been deposited at [CRAN](https://CRAN.R-project.org/package=Rbeast). ( On CRAN, there is another Bayesian time-series package named "beast", which has nothing to do with the BEAST algorithim. Our package is `Rbeast`. Also, `Rbeast` has nothing to do with the famous "Bayesian evolutionary analysis by sampling trees" aglorithm.) Install  `Rbeast` in R using

```R
install.packages("Rbeast")
 ```
 



#### Run and test Rbeast in R
   
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

<table border="0"  style='border:none;'  bordercolor="#ffffff"  width=100%  >
<tr style='border:none;'  >   
   <td valign="center" style='border:none;'  > 
        <img  height="300" align="left" src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/beach.png">
   </td>
   <td valign="center" style='border:none;'  > 
       <img  height="300" align="center"  src="https://github.com/zhaokg/Rbeast/raw/master/R/Images/Nile.png">
   </td>   
 </tr>
 </table>
 
 <!---
 [![](https://github.com/zhaokg/Rbeast/raw/master/R/Images/beach.png)](https://cran.r-project.org/package=Rbeast)
 
 ----
--->

<br/>

## Installation for Matlab  <a name=matlab>    </a>

[![View Rbeast on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/72515-bayesian-changepoint-detection-time-series-decomposition)
 

Install the Matlab version of **BEAST** automatically to a local folder of your choice by running 
  ```Matlab  
  beastPath = 'C:\beast\'                   
  eval( webread('http://b.link/rbeast') )  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%% Note on Automatic Installtion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  % 1. Write permission needed for your chosen path; the variable name must be 'beastPath'.   %
  % 2. If webread has a certificate error, run the following line instead:                    %
      eval(  webread( 'http://b.link/rbeast', weboptions('cert','') )  )                       %
  % 3. If the automatic installation fails, please manually download all the files (see blow) %       
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ```
> The above will download all the files in the [Rbeast\Matlab folder at Github](https://github.com/zhaokg/Rbeast) to the chosen folder: if `beastPath` is missing, a default temporary folder (e.g., `C:\Users\$user_name$\AppData\Local\Temp\Rbeast for Windows 10`) will be used. If the automatic script fails, please download the Matlab files from [Github](https://github.com/zhaokg/Rbeast) manually. These files include a Matlab mex library compiled from the C soure code (e.g., `Rbeast.mexw64` for Windows, `Rbeast.mexa64` for Linux, `Rbeast.mexmaci64` for MacOS) and some Matlab wrapper functions (e.g.,`beast.m`, and `beast123.m`) similar to the R interface, as well as some test datasets (e.g., Nile.mat, and co2.mat). 

> We generated the Matlab mex binary library on our own machines with Win10, Ubuntu 22.04, and macOS High Sierra. If they fail on your machine, the mex library can be compiled  from the C source code files under `Rbeast\Source`. If needed, we are happy to work with you to compile for your specific OS or machines. Additional information on compilations from the C source is also given below.

#### Run and test Rbeast in Matlab
The Matlab API is similar to those of R. Below is a quick example:
  ```Matlab
   help beast
   help beast123  
   load('Nile.mat')                                   % annual streamflow of the Nile River startin from year 1871
   out = beast(Nile, 'season', 'none','start', 1871)  % trend-only data without seasonality
   printbeast(out)
   plotbeast(out)
  ```
## Installation for Octave <a name=octave>    </a>

The same as for Matlab. Now, only Windows platforms are supported. If needed for other platforms (e.g., Octave in Linux and Mac), please contact the author at zhao.1423@osu.edu for support.

----
##  Installation for Python  <a name=python>   </a>


<p  align="left">   
 <a href= "https://pypi.org/project/Rbeast/"> <img src="https://img.shields.io/static/v1?style=plastic&logo=python&label=Python%20%20&message=pip install Rbeast&color=brightgreen" height="20"></a>
</p> 

A package **`Rbeast`** has been deposited at PyPI: https://pypi.org/project/Rbeast/. Run the command below in a console to install:
  
  ```python
    pip install Rbeast
  ```
> Binary wheel files were built on Windows, MacOS, and Linux for Python version 3.7 to 3.11 (either x86_64 or arm64 CPU). If the installation fails, please install from the source to build the package using `pip install Rbeast --no-binary :none:`, which requies a C/C++ compliler (e.g., requiring gcc on Linux or xcode on Mac). If needed, contact Kaiguang Zhao (zhao.1423@osu.edu) to help build the package for your OS platform and Python version.

#### Run and test Rbeast in Python

 
 
`Nile` is annual streamflow of the River Nile, starting from Year 1871. As annual observations, it has no periodic component (i.e., `season='none'`).

  ```python
import Rbeast as rb                                       # Import the Rbeast package as `rb`
nile, year = rb.load_example('nile')                      # a sample time series
o          = rb.beast( nile, start=1871, season='none')
rb.plot(o)
rb.print(o)
o  # see a list of output fields in the output variable o
```

The second example `googletrend` is a monthly time series of the Google Search popularity of the word ***beach*** over the US. This monthly time series is reguarly-spaced (i.e., deltat=`1 month` =`1/12 year`); it has a cyclyic component with a period of 1 year. That is, the number of data points per period is `period` / `deltat` = 1 year / 1 month = 1/(1/12) = 12.

  ```python
beach, year = rb.load_example('googletrend')
o = rb.beast(beach, start= 2004.0, deltat=1/12, period = 1.0)       # the time unit is unknown or arbitrary
o = rb.beast(beach, start= 2004.0, deltat=1/12, period ='1.0 year') # the time unit is fractional year
o = rb.beast(beach, start= 2004.0, deltat='1 month', period =1.0)   # the time unit is fractional year
rb.plot(o)
rb.print(o)
  ```

## Installation for Julia/IDL (yet to come)

Wrappers in Julia and IDL are being developed: We welcome contributions and help from interested developers. If interested, contact Kaiguang Zhao at zhao.1423@osu.edu.

## Description of BEAST
Interpretation of time series data is affected by model choices. Different models can give different or even contradicting estimates of patterns, trends, and mechanisms for the same data–a limitation alleviated by the Bayesian estimator of abrupt change,seasonality, and trend (BEAST) of this package. BEAST seeks to improve time series decomposition by forgoing the "single-best-model" concept and embracing all competing models into the inference via a Bayesian model averaging scheme. It is a flexible tool to uncover abrupt changes (i.e., change-points), cyclic variations (e.g., seasonality), and nonlinear trends in time-series observations. BEAST not just tells when changes occur but also quantifies how likely the detected changes are true. It detects not just piecewise linear trends but also arbitrary nonlinear trends. BEAST is applicable to real-valued time series data of all kinds, be it for remote sensing, finance, public health, economics, climate sciences, ecology, and hydrology. Example applications include its use to identify regime shifts in ecological data, map forest disturbance and land degradation from satellite imagery, detect market trends in economic data, pinpoint anomaly and extreme events in climate data, and unravel system dynamics in biological data. Details on BEAST are reported in [Zhao et al. (2019)](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view). The paper is available at https://go.osu.edu/beast2019.

## Note on computation

As a Bayesian algorithm, BEAST is fast and is possibly among the fastest implementations of Bayesian time-series analysis algorithms of the same nature. (But it is still slower, compared to nonBayesian methods.) For applications dealing with a few to thousands of time series, the computation won't be an practical concern. But for remote sensing/geospatial applications that may easily involve millions or billions of time series, computation will be a big challenge for Desktop computer users. We suggest first testing BEAST on a single time series or small image chips first to determine whether BEAST is appropriate for your applications and, if yes, estimate how long it may take to process the whole image. 

In any case, for those users handling stacked time-series images, do not use `beast` or `beast.irreg`. Use **`beast123`** instead, which can handle 3D data cubes and allow parallel computing.  We also welcome consultation with Kaiguang Zhao (zhao.1423@osu.edu) to give specific suggestions if you see some value of BEAST for your applications.

## Reference
* Zhao, K., Wulder, M. A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick B., Zhang, X., & Brown, M. (2019). [Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm.](https://drive.google.com/file/d/1MFZ0FpK1NwTieVSAf5jicLgl85Lm48uh/view) Remote Sensing of Environment, 232, 111181. (the BEAST paper) 

* Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. [Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection](https://drive.google.com/file/d/1WGOAbH_h5f8ptQvsJsZqwlykO2JevsgT/view?usp=sharing). Remote Sensing of Environment, 132, pp.102-119. (the mcmc sampler used for BEAST)

* Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. [Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine](https://pages.charlotte.edu/gang-chen/wp-content/uploads/sites/184/2021/05/Hu_2021_BEAST-HF-s.pdf). ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261. (an application paper)


## Compilation from C source code (for developers and experts only)

 Though not needed but if preferred, the code can be compiled for your specific machines. Check the [Rbeast\Source](https://github.com/zhaokg/Rbeast/tree/master/Source) folder at GitHub for details.
 

## Reporting Bugs or getting help

BEAST is distributed as is and without warranty of suitability for application. The one distributed above is still a beta version, with potential room for further improvement. If you encounter flaws with the software (i.e. bugs) please report the issue. Providing a detailed description of the conditions under which the bug occurred will help to identify the bug, you can directly email its maintainer Dr. Kaiguang Zhao at zhao.1423@osu.edu. Alternatively, Use the [Issues tracker](https://github.com/zhaokg/Rbeast/issues) on GitHub to report issues with the software and to request feature enhancements. 


## Acknowledgement:

BEAST is developed by Yang Li, Tongxi Hu, Xuesong Zhang, and Kaiguang Zhao. The development of BEAST received supported through Microsoft Azure for Research (CRM0518513) and a USGS 104B grant and a Harmful Algal Bloom Research Initiative grant from the Ohio Department of Higher Education. The contribution of Xuesong Zhang was supported through USDA-ARS.

<h2 id="publicationid" name="publicationid"> Selected publications using BEAST/Rbeast  </h2>   <a name=publications></a>

Despite being published originally for ecological and enviornmental applications, BEAST is developed as a generic tool applicable to time series or time-series-like data arising from all disciplines. BEAST is not a heuristic algorithm but a rigorous statistical model. Below is a list of selected pulications that used BEAST for statistical data analysis.

| Discipline | Publication Title |
| --- | --- |
| Remote Sensing|*Li, J., Li, Z., Wu, H., and You, N., 2022. [Trend, seasonality, and abrupt change detection method for land surface temperature time-series analysis: Evaluation and improvement](https://doi.org/10.1016/j.rse.2022.113222). Remote Sensing of Environment, 10.1016/j.rse.2022.113222*|
|   Population Ecology  | *Henderson, P. A. (2021). [Southwood's Ecological Methods (5th edition)](https://www.google.com/books/edition/Southwood_s_Ecological_Methods/snEhEAAAQBAJ?hl=en&gbpv=1&pg=PA473&printsec=frontcover). Oxford University Press., page 475-476*|
|Climatology|*Webster, M.A., Riihelä, A., Kacimi, S., Ballinger, T.J., Blanchard-Wrigglesworth, E., Parker, C.L. and Boisvert, L., 2024. Summer snow on Arctic sea ice modulated by the Arctic Oscillation. Nature Geoscience, pp.1-8.*|
|Cardiology|*Ozier, D., Rafiq, T., de Souza, R. and Singh, S.M., 2023. [Use of Sacubitril/Valsartan Prior to Primary Prevention Implantable Cardioverter Defibrillator Implantation](https://doi.org/10.1016/j.cjco.2022.10.005). CJC Open.*|
|Political Science|*Reuning, K., Whitesell, A. and Hannah, A.L., 2022. [Facebook algorithm changes may have amplified local republican parties](https://doi.org/10.1177/20531680221103809). Research & Politics, 9(2), p.20531680221103809.*|
|Paleoclimatology|*Anastasia Zhuravleva et al., 2023. Caribbean salinity anomalies contributed to variable North Atlantic circulation and climate during the Common Era.  Science Advances,   DOI:10.1126/sciadv.adg2639*|
|Epidemiology|*Perofsky, A.C., Huddleston, J., Hansen, C., Barnes, J.R., Rowe, T., Xu, X., Kondor, R., Wentworth, D.E., Lewis, N., Whittaker, L. and Ermetal, B., 2023. Antigenic drift and subtype interference shape A (H3N2) epidemic dynamics in the United States. medRxiv, pp.2023-10.*|
|Mechanic Engineering|*Bao, X., Chen, L., Zhong, J., Wu, D. and Zheng, Y., 2024. A self-supervised contrastive change point detection method for industrial time series. Engineering Applications of Artificial Intelligence, 133, p.108217.*|
|Biomedical Engineering|*Saghbiny, E., Da Silva, J., Leblanc, L., Bobbio, C., Morel, G.G., Vialle, R. and Tamadazte, B., 2023, September. Breach detection in spine surgery based on cutting torque with ex-vivo experimental validation. In Conference on New Technologies for Computer and Robot Assisted Surgery.*|
|Food Science|*Zaytsev, V., Tutukina, M.N., Chetyrkina, M.R., Shelyakin, P.V., Ovchinnikov, G., Satybaldina, D., Kondrashov, V.A., Bandurist, M.S., Seilov, S., Gorin, D.A. and Fedorov, F.S., 2024. Monitoring of meat quality and change-point detection by a sensor array and profiling of bacterial communities. Analytica Chimica Acta, p.343022.*|
|Ecology|*Dashti, H., Chen, M., Smith, B., Zhao, K. and Moore, D., 2024. Rethinking Ecosystems Disturbance Recovery: what it was or what it could have been?. Geophysical Research Letters*|
|Atmospheric Science|*Ganji, A., Saeedi, M., Lloyd, M., Xu, J., Weichenthal, S. and Hatzopoulou, M., 2024. Air pollution prediction and backcasting through a combination of mobile monitoring and historical on-road traffic emission inventories. Science of the Total Environment, 915, p.170075.*|
|Coastal Ecology|*Xiang, J., Cui, T., Li, X., Zhang, Q., Mu, B., Liu, R. and Zhao, W., 2023. Evaluating the effectiveness of coastal environmental management policies in China: The case of Bohai Sea. Journal of Environmental Management, 338, p.117812.*|
|Vegetation Science|*Fan, X., Qu, Y., Zhang, J. and Bai, E., 2024. China's vegetation restoration programs accelerated vegetation greening on the Loess Plateau. Agricultural and Forest Meteorology, 350, p.109994.*|
|Geophysics|*Mu, Y., Biggs, T.W. and Jones, C., 2023. Importance in shifting circulation patterns for dry season moisture sources in the Brazilian Amazon. Geophysical Research Letters, 50(9), p.e2023GL103167.*|
|Spatial Ecology|*Laurin, G.V., Cotrina-Sanchez, A., Belelli-Marchesini, L., Tomelleri, E., Battipaglia, G., Cocozza, C., Niccoli, F., Kabala, J.P., Gianelle, D., Vescovo, L. and Da Ros, L., 2024. Comparing ground below-canopy and satellite spectral data for an improved and integrated forest phenology monitoring system. Ecological Indicators, 158, p.111328.*|
|Anthropocene Science|*Thomas, E.R., Vladimirova, D.O., Tetzner, D.R., Emanuelsson, D.B., Humby, J., Turner, S.D., Rose, N.L., Roberts, S.L., Gaca, P. and Cundy, A.B., 2023. The Palmer ice core as a candidate Global boundary Stratotype Section and Point for the Anthropocene series. The Anthropocene Review, p.20530196231155191.*|
|Computer Science|*Xue, X., Mao, H. and Li, Q., 2023. Dag-acfl: Asynchronous clustered federated learning based on dag-dlt. arXiv preprint arXiv:2308.13158.*|
|Spatial Hydrology|*Wang, Yiming, Xuesong Zhang, Kaiguang Zhao, and Debjani Singh. "Streamflow in the United States: Characteristics, trends, regime shifts, and extremes." Scientific Data 11, no. 1 (2024): 788.*|
|Business Science|*Li, Z. and Tian, Y., 2024. Skewed multifractal cross-correlation between price and volume during the COVID-19 pandemic: Evidence from China and European carbon markets. Applied Energy, 371, p.123716.*|
|Ecography|*Smith, M.M. and Pauli, J.N., 2024. Small but connected islands can maintain populations and genetic diversity under climate change. Ecography, p.e07119.*|
|Forestry|*Mulverhill, C., Coops, N.C., White, J.C., Tompalski, P. and Achim, A., 2024. Evaluating the potential for continuous update of enhanced forest inventory attributes using optical satellite data. Forestry: An International Journal of Forest Research, p.cpae029.*|
|Ecosystem Science|*Moreno-Fernández, D., Camarero, J.J., García, M., Lines, E.R., Sánchez-Dávila, J., Tijerín, J., Valeriano, C., Viana-Soto, A., Zavala, M.A. and Ruiz-Benito, P., 2022. The interplay of the tree and stand-level processes mediate drought-induced forest dieback: evidence from complementary remote sensing and tree-ring approaches. Ecosystems, 25(8), pp.1738-1753.*|
|Economics|*Sapkota, B.P., 2024. Analysis of Climate Policy and Monetary Policy Nexus in the Norwegian Context: A DSGE Approach (Master's thesis, Norwegian University of Life Sciences).*|
|Cognitive Science|*Prein, J.C., Maurits, L., Werwach, A., Haun, D.B. and Bohn, M., 2024. Variation in gaze following across the life span: A process‐level perspective. Developmental Science, p.e13546.*|
|Neuroscience|*Aqel, K., Wang, Z., Peng, Y.B. and Maia, P.D., 2024. Reconstructing rodent brain signals during euthanasia with eigensystem realization algorithm (ERA). Scientific Reports, 14(1), p.12261.*|
|Energies|*Livera, A., Tziolis, G., Theristis, M., Stein, J.S. and Georghiou, G.E., 2023. Estimating the performance loss rate of photovoltaic systems using time series change point analysis. Energies, 16(9), p.3724.*|
|Glaciology|*Ramón, C.L., Rueda, F.J., Priet‐Mahéo, M.C. and Andradóttir, H., 2024. The impact of deep glacial water diversions from a hydroelectric reservoir in the thermal dynamics of a sub-arctic lake. Journal of Hydrology, 635, p.131081.*|
|Biomedical Engineering|*Ling, Y.T., 2022. Ultrafast ultrasound imaging for muscle onset mapping during its contraction.*|
|Education Science|*McNeill, S.M., 2022. Three Essays Examining Increases in State Elementary-Secondary Education Spending (Doctoral dissertation).*|
|Cryospheric Sciences|*Wauthy, S., Tison, J.L., Inoue, M., El Amri, S., Sun, S., Fripiat, F., Claeys, P. and Pattyn, F., 2024. Spatial and temporal variability of environmental proxies from the top 120 m of two ice cores in Dronning Maud Land (East Antarctica). Earth System Science Data, 16(1), pp.35-58.*|
|Quaternary Science|*Gibson, D.K., Bird, B.W., Finney, B.P. and Steinman, B.A., 2024. Holocene insolation and sea surface temperature influences on the polar front jet stream and precipitation in the midcontinental United States. Quaternary Science Reviews, 340, p.108865.*|
|GIScience|*Morresi, D., Maeng, H., Marzano, R., Lingua, E., Motta, R. and Garbarino, M., 2024. High-dimensional detection of Landscape Dynamics: a Landsat time series-based algorithm for forest disturbance mapping and beyond. GIScience & Remote Sensing, 61(1), p.2365001.*|
|Geography|*Lyu, R., Pang, J., Zhang, J. and Zhang, J., 2024. The impacts of disturbances on mountain ecosystem services: Insights from BEAST and Bayesian network. Applied Geography, 162, p.103143.*|
|Watershed Hydrology|*Sakizadeh, M., Milewski, A. and Sattari, M.T., 2023. Analysis of Long-Term Trend of Stream Flow and Interaction Effect of Land Use and Land Cover on Water Yield by SWAT Model and Statistical Learning in Part of Urmia Lake Basin, Northwest of Iran. Water, 15(4), p.690.*|
|Limnology|*Wang, M., Chen, X., Cao, L., Kurban, A., Shi, H., Wu, N., Eziz, A., Yuan, X. and De Maeyer, P., 2023. Correlation analysis between the Aral Sea shrinkage and the Amu Darya River. Journal of Arid Land, 15(7), pp.757-778.*|
|Phenology|*Marien, B., Papadimitriou, D., Kotilainen, T., Zuccarini, P., Dox, I., Verlinden, M., Heinecke, T., Marien, J., Willems, P., Decoster, M. and Gasco, A., 2022. Timing leaf senescence: A generalized additive models for location, scale and shape approach. Agricultural and Forest Meteorology, 315, p.108823.*|
|Oceanography|*Huang, W., Zhu, X., Jin, Y. and Shen, X., 2024. Nonstationary modelling of significant wave height using time series decomposition method. Ocean Engineering, 310, p.118731.*|
|Oceanography|*Oehlert, A.M., Hunter, H., Riopelle, C. and Purkis, S.J., 2023. Perturbation to North Atlantic Ocean‐Climate Dynamics Tripled Whitings Mud Production in the Bahamas. Journal of Geophysical Research: Oceans, 128(11), p.e2023JC020021.*|
| Hydraulic Engineering | *Xu, X., Yang, J., Ma, C., Qu, X., Chen, J. and Cheng, L., 2022. Segmented modeling method of dam displacement based on BEAST time series decomposition. Measurement, 202, p.111811.* |
|Social Media|*Barrie, C., Ketchley, N., Siegel, A. and Bagdouri, M., 2023. Measuring Media Freedom.*|
|Political Economy|*Benchimol, J. and Palumbo, L., 2023. Sanctions and Russian Online Prices.*|
|Physiology|*Shakeel, M., Brockmann, A. Temporal effects of sugar intake on fly local search and honey bee dance behaviour. J Comp Physiol A (2023). https://doi.org/10.1007/s00359-023-01670-6*|
|Injuries & Hazards|*Delavary, M., Kalantari, A.H., Mohammadzadeh Moghaddam, A., Fakoor, V., Lavallière, M. and Wilhelm Siebert, F., 2024. Road traffic mortality in Iran: longitudinal trend and seasonal analysis, March 2011-February 2020. International journal of injury control and safety promotion, 31(1), pp.125-137.*|
|Grassland Ecology|*Wang, F., Men, R., Lai, H., Feng, K., Yan, S., Gao, S., Wang, Z., Tian, Q., Guo, W. and Yang, H., 2024. The response of vegetation dynamics to drought and its driving factors identification in Inner Mongolia of China. Ecological Indicators, 164, p.112125.*|
|Soil Sciences|*Gökpinar, T. and Heinze, T., 2024. Algorithm-based segmentation of temperature-depth profiles: Examples from a mine. Environmental Modelling & Software, 180, p.106143.*|
|Regional Studies|*Vo, T.H. and Liou, Y.A., 2024. Four-decade spring droughts in Taiwan. Journal of Hydrology: Regional Studies, 54, p.101849.*|
|Sociopolitical Science|*Savolainen, S., Saarinen, V.P. and Chen, T.H.Y., 2024. Social Media Affordances Sustain Social Movements Facing Repression: Evidence from Climate Activism (No. p4yvk). Center for Open Science.*|
|Civil Engineering|*Langtry, M., Wichitwechkarn, V., Ward, R., Zhuang, C., Kreitmair, M.J., Makasis, N., Conti, Z.X. and Choudhary, R., 2024. Impact of data for forecasting on performance of model predictive control in buildings with smart energy storage. Energy and Buildings, p.114605.*|
|Ichthyology|*Kaeding, L.R., 2023. Climate-change and nonnative-piscivore impacts on a renowned Oncorhynchus metapopulation, requirements for metapopulation recovery, and similarities to anadromous salmonid metapopulations. Aquatic Sciences, 85(4), p.88.*|
|Earth Science|*Yang, Z., Peng, J., Liu, Y., Jiang, S., Cheng, X., Liu, X., Dong, J., Hua, T. and Yu, X., 2024. GloUTCI-M: a global monthly 1 km Universal Thermal Climate Index dataset from 2000 to 2022. Earth System Science Data, 16(5), pp.2407-2424.*|
|Electrical Engineering|*Nguyen, T.K., Ahmad, Z. and Kim, J.M., 2023. Leak State Detection and Size Identification for Fluid Pipelines with a Novel Acoustic Emission Intensity Index and Random Forest. Sensors, 23(22), p.9087.*|
|Remote Sensing|*Mulverhill, C., Coops, N.C. and Achim, A., 2023. Continuous monitoring and sub-annual change detection in high-latitude forests using Harmonized Landsat Sentinel-2 data. ISPRS Journal of Photogrammetry and Remote Sensing, 197, pp.309-319.*|
|Remote Sensing|*Huang, Y., Li, G., Zhao, Y., Yang, J. and Li, Y., 2023. Analysis of the Characteristics and Causes of Land Degradation and Development in Coastal China (1982–2015). Remote Sensing, 15(9), p.2249.*|
|Environmental Science|*Gore, B.T.O., Aman, A., Kouadio, Y. and Duclos, O.M., 2023. Recent Vegetation Cover Dynamics and Climatic Parameters Evolution Study in the Great Green Wall of Senegal. Journal of Environmental Protection, 14(4), pp.254-284.*|
|Marine Science|*Stamieszkin, K., Millette, N.C., Luo, J.Y., Follett, E., Record, N.R. and Johns, D.G., 2024. Large protistan mixotrophs in the North Atlantic Continuous Plankton Recorder time series: associated environmental conditions and trends. Frontiers in Marine Science, 11, p.1320046.*|
|Marine Science|*Haditiar, Y., Ikhwan, M., Mahdi, S., Siregar, A.N., Haridhi, H.A., Setiawan, I., Nanda, M., Prajaputra, V. and Irham, M., 2024. Oceanographic characteristics in the North of Aceh waters. Regional Studies in Marine Science, 71, p.103408.*|
|Physical Chemistry|*Faran, M. and Bisker, G., 2023. Nonequilibrium Self-Assembly Time Forecasting by the Stochastic Landscape Method. The Journal of Physical Chemistry B.*|
|Biogeochemistry|*Dahl, M., Gullström, M., Bernabeu, I., Serrano, O., Leiva‐Dueñas, C., Linderholm, H.W., Asplund, M.E., Björk, M., Ou, T., Svensson, J.R. and Andrén, E., 2024. A 2,000‐year record of eelgrass (Zostera marina L.) colonization shows substantial gains in blue carbon storage and nutrient retention. Global Biogeochemical Cycles, 38(3), p.e2023GB008039.*|
|Data Science|*Livina, V., 2023. Tipping point analysis and its applications.*|
|Data Science|*Roudiere, E., 2024. Data analytics and machine learning for railway track degradation: Using Bothnia Line track measurements for maintenance forecasting.*|
|Soil Science|*Woo, D.K. and Seo, Y., 2022. Effects of elevated temperature and abnormal precipitation on soil carbon and nitrogen dynamics in a Pinus densiflora forest. Frontiers in Forests and Global Change, 5, p.1051210.*|
|Mechanobiology|*Faran, M., Ray, D., Nag, S., Raucci, U., Parrinello, M. and Bisker, G., 2024. A Stochastic Landscape Approach for Protein Folding State Classification. Journal of Chemical Theory and Computation.*|
|Analytical Chemistry|*Simic, M., Neuper, C., Hohenester, U. and Hill, C., 2023. Optofluidic force induction as a process analytical technology. Analytical and Bioanalytical Chemistry, pp.1-11.*|
| Ecosystem Sciences | *Lyu, R., Zhao, W., Pang, J., Tian, X., Zhang, J. and Wang, N., 2022. [Towards a sustainable nature reserve management: Using Bayesian network to quantify the threat of disturbance to ecosystem services](https://doi.org/10.1016/j.ecoser.2022.101483). Ecosystem Services, 58, p.101483.* |
| Environmental Sciences| *Nickerson, S., Chen, G., Fearnside, P., Allan, C.J., Hu, T., de Carvalho, L.M. and Zhao, K., 2022. [Forest loss is significantly higher near clustered small dams than single large dams per megawatt of hydroelectricity installed in the Brazilian Amazon](https://iopscience.iop.org/article/10.1088/1748-9326/ac8236/meta). Environmental Research Letters.*|
|Geology| *Fan, X., Goeppert, N. and Goldscheider, N., 2023. Quantifying the historic and future response of karst spring discharge to climate variability and change at a snow-influenced temperate catchment in central Europe. Hydrogeology Journal, pp.1-17.*|
|Hydro-climatology|*Nair, A., 2022. Trend analysis of hydro-climatological factors using a bayesian ensemble algorithm with reasoning from dynamic and static variables. Atmosphere, 13(12), p.1961.*|
|Wildlife|*Smith, Matthew M., and Jonathan N. Pauli. "Connectivity maintains genetic diversity and population persistence within an archipelagic refugia even under declining lake ice." Mechanisms of species recovery for a forest carnivore in a changing landscape: 173.*|
|  Climate Sciences|*Duke, N.C., Mackenzie, J.R., Canning, A.D., Hutley, L.B., Bourke, A.J., Kovacs, J.M., Cormier, R., Staben, G., Lymburner, L. and Ai, E., 2022. [ENSO-driven extreme oscillations in mean sea level destabilise critical shoreline mangroves—An emerging threat](https://journals.plos.org/climate/article?id=10.1371/journal.pclm.0000037). PLOS Climate, 1(8), p.e000003*|
| Finance|*Candelaria, Christopher A., Shelby M. McNeill, and Kenneth A. Shores. (2022). What is a School Finance Reform? Uncovering the ubiquity and diversity of school finance reforms using a Bayesian changepoint estimator.(EdWorkingPaper: 22-587). Retrieved from Annenberg Institute at Brown University: https://doi.org/10.26300/4vey-3w10*|
| Public health|*Linnell, K., Fudolig, M., Schwartz, A., Ricketts, T.H., O'Neil-Dunne, J.P., Dodds, P.S. and Danforth, C.M., 2022. [Spatial changes in park visitation at the onset of the pandemic](https://journals.plos.org/globalpublichealth/article?id=10.1371/journal.pgph.0000766). arXiv preprint arXiv:2205.15937.*|
|Biometerology|*Li, Y., Liu, Y., Bohrer, G., Cai, Y., Wilson, A., Hu, T., Wang, Z. and Zhao, K., 2022. [Impacts of forest loss on local climate across the conterminous United States: Evidence from satellite time-series observations](https://doi.org/10.1016/j.scitotenv.2021.149651). Science of The Total Environment, 802, p.149651.*|
| Applied Math|*Ferguson, Daniel, and Francois G. Meyer. [Probability density estimation for sets of large graphs with respect to spectral information using stochastic block models](https://arxiv.org/abs/2207.02168). arXiv preprint arXiv:2207.02168 (2022).*|
|Transportation Science| *Delavary, M., Kalantari, A.H., Mohammadzadeh Moghaddam, A., Fakoor, V., Lavalliere, M. and Wilhelm Siebert, F., 2023. Road traffic mortality in Iran: longitudinal trend and seasonal analysis, March 2011-February 2020. International Journal of Injury Control and Safety Promotion, pp.1-12.*|
|Water quality|*He, Ziming, Jiayu Yao, Yancen Lu, and Danlu Guo. "Detecting and explaining long-term changes in river water quality in south eastern Australia." Hydrological Processes: e14741.*|
|Air quality|*Wu, S., Yao, J., Wang, Y. and Zhao, W., 2023. Influencing factors of PM2. 5 concentrations in the typical urban agglomerations in China based on wavelet perspective. Environmental Research, p.116641.*|
| Hydrology | *Zohaib, M. and Choi, M., 2020. [Satellite-based global-scale irrigation water use and its contemporary trends](https://doi.org/10.1016/j.rse.2022.113222). Science of The Total Environment, 714, p.136719.* |
| Energy Engineering |*Lindig, S., Theristis, M. and Moser, D., 2022. Best practices for photovoltaic performance loss rate calculations. Progress in Energy, 4(2), p.022003.*|
|Virology|*Shen, L., Sun, M., Song, S., Hu, Q., Wang, N., Ou, G., Guo, Z., Du, J., Shao, Z., Bai, Y. and Liu, K., 2022. The impact of anti-COVID19 nonpharmaceutical interventions on hand, foot, and mouth disease—A spatiotemporal perspective in Xi'an, northwestern China. Journal of medical virology.*|
| Pharmaceutical Sciences|*Patzkowski, M.S., Costantino, R.C., Kane, T.M., Nghiem, V.T., Kroma, R.B. and Highland, K.B., 2022. Military Health System Opioid, Tramadol, and Gabapentinoid Prescription Volumes Before and After a Defense Health Agency Policy Release. Clinical Drug Investigation, pp.1-8.*|
| Geography|*Cai, Y., Liu, S. and Lin, H., 2020. Monitoring the vegetation dynamics in the Dongting Lake Wetland from 2000 to 2019 using the BEAST algorithm based on dense Landsat time series. Applied Sciences, 10(12), p.4209.*|
| Oceanography|*Pitarch, J., Bellacicco, M., Marullo, S. and Van Der Woerd, H.J., 2021. [Global maps of Forel-Ule index, hue angle and Secchi disk depth derived from 21 years of monthly ESA Ocean Colour Climate Change Initiative data](https://essd.copernicus.org/articles/13/481/2021/). Earth System Science Data, 13(2), pp.481-490.*|
|Photovoltaics|*Micheli, L., Theristis, M., Livera, A., Stein, J.S., Georghiou, G.E., Muller, M., Almonacid, F. and Fernadez, E.F., 2021. Improved PV soiling extraction through the detection of cleanings and change points. IEEE Journal of Photovoltaics, 11(2), pp.519-526.*|
|Climate Sciences|*White, J.H., Walsh, J.E. and Thoman Jr, R.L., 2021. [Using Bayesian statistics to detect trends in Alaskan precipitation](https://doi.org/10.1002/joc.6946). International Journal of Climatology, 41(3), pp.2045-2059.*|
|Glaciology|*Li, G., Lv, M., Quincey, D.J., Taylor, L.S., Li, X., Yan, S., Sun, Y. and Guo, H., 2023. Characterizing the surge behaviour and associated ice-dammed lake evolution of the Kyagar Glacier in the Karakoram. The Cryosphere, 17(7), pp.2891-2907.*|
|Hydrogeology|*Fan, X., Goeppert, N. and Goldscheider, N., 2023. Quantifying the historic and future response of karst spring discharge to climate variability and change at a snow-influenced temperate catchment in central Europe. Hydrogeology Journal, 31(8), pp.2213-2229.*|
|Siviculture|*Kacic, P., Gessner, U., Holzwarth, S., Thonfeld, F. and Kuenzer, C., 2024. Assessing experimental silvicultural treatments enhancing structural complexity in a central European forest–BEAST time‐series analysis based on Sentinel‐1 and Sentinel‐2. Remote Sensing in Ecology and Conservation.*|
|Field Hydrology|*Merk, M., Goeppert, N. and Goldscheider, N., 2021. [Deep desiccation of soils observed by long-term high-resolution measurements on a large inclined lysimeter](https://hess.copernicus.org/articles/25/3519/2021/). Hydrology and Earth System Sciences, 25(6), pp.3519-3538.*|
|Sports Science|*Roccetti, M. and Piconi, F., Minutes played by Under 21 players in Serie A: a descriptive analytics study.*|
|Remote Sensing|*Chen, X., Zhao, X., Zhao, Y., Wang, R., Lu, J., Zhuang, H. and Bai, L., 2023. Interaction of Climate Change and Anthropogenic Activity on the Spatiotemporal Changes of Surface Water Area in Horqin Sandy Land, China. Remote Sensing, 15(7), p.1918.*|
|Remote Sensing|*Many, G., Escoffier, N., Ferrari, M., Jacquet, P., Odermatt, D., Mariethoz, G., Perolo, P. and Perga, M.E., 2022. Long-term spatiotemporal variability of whitings in Lake Geneva from multispectral remote sensing and machine learning. Remote Sensing, 14(23), p.6175.*|
|Remote Sensing|*Ma, Q., Liu, Y., Qiu, T., Huang, T., Deng, T., Hu, Z. and Cui, T., 2022. Satellite-Observed Four-Dimensional Spatiotemporal Characteristics of Maritime Aerosol Types over the Coastal Waters of the Guangdong–Hong Kong–Macao Greater Bay Area and the Northern South China Sea. Remote Sensing, 14(21), p.5464.*|
|Water Science|*Yıldız, M.B., Di Nunno, F., Đurin, B., Pham, Q.B., de Marinis, G. and Granata, F., 2024. A Combined Seasonal Mann–Kendall and Innovative Approach for the Trend Analysis of Streamflow Rate in Two Croatian Rivers. Water, 16(10), p.1422.*|
|Agriculture|*Gao, Y.X., Leng, P., Li, J., Shang, G.F., Zhang, X. and Li, Z.L., 2024. Identification of irrigation events using Bayesian statistics-based change detection and soil moisture measurements. Agricultural Water Management, 302, p.108999.*|
|Agriculture|*Chen, L., He, Z., Gu, X., Xu, M., Pan, S., Tan, H. and Yang, S., 2023. Construction of an Agricultural Drought Monitoring Model for Karst with Coupled Climate and Substratum Factors—A Case Study of Guizhou Province, China. Water, 15(9), p.1795.*|
|Climate Science|*Shazil, M.S., Mahmood, S.A., Ahmad, S., Haseeb, M., Masood, A., Qureshi, J. and Batool, S., 2024. Assessing long-term variability and trends in temperature and precipitation in Gilgit and Hunza river basins. Environmental Earth Sciences, 83(8), p.248.*|
|Forest Ecology|*Moreno-Fernandez, D., Viana-Soto, A., Camarero, J.J., Zavala, M.A., Tijerin, J. and Garcia, M., 2021. Using spectral indices as early warning signals of forest dieback: The case of drought-prone Pinus pinaster forests. Science of The Total Environment, 793, p.148578.*|
|Petroleum  Engineering|*Pan, Y., Bi, R., Yang, S., Lyu, Z. and Ju, X., 2024, February. Application of a Bayesian Ensemble Algorithm for Automated Production Diagnostic of Gas Wells with Plunger-Lift. In International Petroleum Technology Conference (p. D011S029R008). IPTC.*|
|Atmospheric Sciences|*Tingwei, C., Tingxuan, H., Bing, M., Fei, G., Yanfang, X., Rongjie, L., Yi, M. and Jie, Z., 2021. Spatiotemporal pattern of aerosol types over the Bohai and Yellow Seas observed by CALIOP. Infrared and Laser Engineering, 50(6), p.20211030.*|
|Terrestrial ecology|*Dashti, H., Pandit, K., Glenn, N.F., Shinneman, D.J., Flerchinger, G.N., Hudak, A.T., de Graaf, M.A., Flores, A., Ustin, S., Ilangakoon, N. and Fellows, A.W., 2021. [Performance of the ecosystem demography model (EDv2. 2) in simulating gross primary production capacity and activity in a dryland study area](https://doi.org/10.1016/j.agrformet.2020.108270). Agricultural and Forest Meteorology, 297, p.108270.*|
|Statistics|*Storath, M. and Weinmann, A., 2023. Smoothing splines for discontinuous signals. Journal of Computational and Graphical Statistics, (just-accepted), pp.1-26
|Global Ecology|*Seo, Minji, and Hyun-Cheol Kim. "Arctic Greening Trends: Change Points in Satellite-Derived Normalized Difference Vegetation Indexes and Their Correlation with Climate Variables over the Last Two Decades." Remote Sensing 16, no. 7 (2024): 1160.*|
|Sustainability|*Baade, J., Gessner, U., Hahndiek, E., Harmse, C., Hill, S., Hirner, A., Maruping-Mzileni, N., Otte, I., Pathe, C., Renner, P. and Schellenberg, K., 2024. Observational support for regional policy implementation: land surface change under anthropogenic and climate pressure in SALDi study sites. In Sustainability of Southern African Ecosystems under Global Change: Science for Management and Policy Interventions (pp. 845-877). Cham: Springer International Publishing.*|
|Biogeochemistry|*Yang, R., Wu, J., Zhao, J., Guan, Q., Fan, X. and Zhao, L., 2024. Marked interannual variability in the relative dominance of phytoplankton over submerged macrophytes rather than regime shifts in a shallow eutrophic lake: Evidence from long-term observations. Ecological Indicators, 166, p.112301.*|
|Environmental Engineering|*Bainbridge, R., Lim, M., Dunning, S., Winter, M.G., Diaz-Moreno, A., Martin, J., Torun, H., Sparkes, B., Khan, M.W. and Jin, N., 2022. Detection and forecasting of shallow landslides: lessons from a natural laboratory. Geomatics, Natural Hazards and Risk, 13(1), pp.686-704.*|
|Fishery|*Theis, S., Wallace, A.,, Poesch, M.,  Portiss, R. & Ruppert, J.. (2024). Balancing boat-electrofishing sampling effort against costs for nearshore fish communities in the Toronto waterfront, Lake Ontario. Fisheries Management and Ecology. 10.1111/fme.12733.*|
|Water Science|*Sakizadeh, M., Milewski, A. and Sattari, M.T., 2023. Analysis of long-term trend of stream flow and interaction effect of land use and land cover on water yield by SWAT model and statistical learning in part of Urmia Lake Basin, Northwest of Iran. Water, 15(4), p.690.*|
|Water Science|*Wang, F., Men, R., Yan, S., Wang, Z., Lai, H., Feng, K., Gao, S., Li, Y., Guo, W. and Tian, Q., 2024. Identification of the Runoff Evolutions and Driving Forces during the Dry Season in the Xijiang River Basin. Water, 16(16), p.2317.*|
|Water Science|*Sikora, M., Canteri, E., Fernandez-Guerra, A., Oskolkov, N., Agren, R., Hansson, L., Irving-Pease, E.K., Mühlemann, B., Holtsmark Nielsen, S., Scorrano, G. and Allentoft, M.E., 2023. The landscape of ancient human pathogens in Eurasia from the Stone Age to historical times. bioRxiv, pp.2023-10.*|
|Urban Science|*Niazmardi, S., Sadrykia, M. and Rezazadeh, M., 2023. Analysis of spatiotemporal household water consumption patterns and their relationship with meteorological variables. Urban Climate, 52, p.101707.*|
|Hydrology|*Di Nunno, F., de Marinis, G. and Granata, F., 2024. Analysis of SPI index trend variations in the United Kingdom-A cluster-based and bayesian ensemble algorithms approach. Journal of Hydrology: Regional Studies, 52, p.101717.*|
|Hydrology|*Song, X., Chen, H., Chen, T., Qin, Z., Chen, S., Yang, N. and Deng, S., 2024. GRACE-based groundwater drought in the Indochina Peninsula during 1979–2020: Changing properties and possible teleconnection mechanisms. Science of The Total Environment, 908, p.168423.*|
|Hydrology|*Sahu, R.T., Verma, M.K. and Ahmad, I., 2023. Density-based spatial clustering of application with noise approach for regionalisation and its effect on hierarchical clustering. International Journal of Hydrology Science and Technology, 16(3), pp.240-269.*|
|Hydrology|*Yang, X., Tian, S., You, W. and Jiang, Z., 2021. Reconstruction of continuous GRACE/GRACE-FO terrestrial water storage anomalies based on time series decomposition. Journal of Hydrology, 603, p.127018.*|
|Hydrology|*Di Nunno, F. and Granata, F., 2024. Analysis of trends and abrupt changes in groundwater and meteorological droughts in the United Kingdom. Journal of Hydrology, p.131430.*|
|Hydrology|*Kong, L., Ma, L., Li, Y., Abuduwaili, J. and Zhang, J., 2024. Assessing the intensity of the water cycle utilizing a Bayesian estimator algorithm and wavelet coherence analysis in the Issyk-Kul Basin of Central Asia. Journal of Hydrology: Regional Studies, 52, p.101680.*|
|Hydrology|*Ramón, C.L., Rueda, F.J., Priet‐Mahéo, M.C. and Andradóttir, H., 2024. The impact of deep glacial water diversions from a hydroelectric reservoir in the thermal dynamics of a sub-arctic lake. Journal of Hydrology, 635, p.131081.*|
|Landscape Ecology|*Adams, B.T., Matthews, S.N., Iverson, L.R., Prasad, A.M., Peters, M.P. and Zhao, K., 2021. Spring phenological variability promoted by topography and vegetation assembly processes in a temperate forest landscape. Agricultural and Forest Meteorology, 308, p.108578.*|

