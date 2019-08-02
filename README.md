# BEAST:  A Bayesian Ensemble Algorithm for Change-Point Detection and Time Series Decomposition

## Description
Interpretation of time series data is affected by model choices. Different models can give different or even contradicting estimates of patterns, trends, and mechanisms for the same data–a limitation alleviated by the Bayesian estimator of abrupt change,seasonality, and trend (BEAST) of this package. BEAST seeks to improve time series decomposition by forgoing the "single-best-model" concept and embracing all competing models into the inference via a Bayesian model averaging scheme. It is a flexible tool to uncover abrupt changes (i.e., change-points), cyclic variations (e.g., seasonality), and nonlinear trends in time-series observations. BEAST not just tells when changes occur but also quantifies how likely the detected changes are true. It detects not just piecewise linear trends but also arbitrary nonlinear trends. BEAST is applicable to real-valued time series data of all kinds, be it for remote sensing, economics, climate sciences, ecology, and hydrology. Example applications include its use to identify regime shifts in ecological data, map forest disturbance and land degradation from satellite imagery, detect market trends in economic data, pinpoint anomaly and extreme events in climate data, and unravel system dynamics in biological data. Details on BEAST are reported in [Zhao et al. (2019)](http://authors.elsevier.com/c/1ZTiS7qzSnIRT). A 50-days free access to the paper is available at http://authors.elsevier.com/c/1ZTiS7qzSnIRT.


>The core of **BEAST** was impemented through the mixed use of C and Fortran, and the source code is under "Source" folder. To run BEAST, R and Matlab interfaces are provided and can be found under the "R" and "Matlab" folders, respectively. Follow the instructions below to install and test BEAST in R or Matlab.


 ---- 
## R

#### Installation

1. **CRAN**: The R package, named **Rbeast**, has already been deposited at [CRAN](https://CRAN.R-project.org/package=Rbeast). (On CRAN, there is a similar Bayesian time-series package called "beast", which has nothing to do with the BEAST algorithim. Our package name is "Rbeast".) The easist way to intall it is to run below in your R console:

```R
install.packages("Rbeast")
```

2. **GitHub**: Pre-compiled binary package files for **Rbeast** are available at [GitHub](https://github.com/zhaokg/Rbeast). Alternative ways to install Rbeast are:

* Windows x86 or x64

```R
install.packages("https://raw.github.com/zhaokg/Rbeast/master/R/CompiledPackage/Windows/Rbeast_0.2.1.zip" ,repos=NULL)
```

* Mac
```R
install.packages("https://raw.github.com/zhaokg/Rbeast/master/R/CompiledPackage/Mac/Rbeast_0.2.1.tgz" ,repos=NULL)
```

#### Run and test Rbeast

The main function in Rbeast is `beast(data, option=list(),demoGUI=FALSE,...)`. The following code snippt provides a starting point for the basic usage of Rbeast.

```R
library(Rbeast)
data(modis_ohio) # a MODIS NDVI time series for a forest pixel in Ohio
plot(modis_ohio) # Rbeast/BEAST can only process uniformly-spaced time 
                 # series. Irregularly-spaced data have to be aggregated
                 #/resampled to regular time intervals before running beast
out=beast(modis_ohio)
plot(out)
?Rbeast          # See more details about the usage of `beast`                 
```


 ---- 
## Matlab

#### Installation and usage (Windows x64 only)

We generated the Matlab mex library only for the Windows 64 OS:  `beast_default.mexw64` and `beast_mkl.mexw64` under the "Rbeast\Matlab" folder. Mex libraries for other OS systems such as Linux and Mac can be compiled from the source code files under "\Rbeast\Source". If needed, we are happy to work with you to compile for your specific OS or machines.

The two Win64 Matlab mex libraries ( `beast_default.mexw64` and `beast_mkl.mexw64`) are the same BEAST algorithm but linked against different math libs. `beast_default` is generated using the standard Lapack lib (http://www.netlib.org/lapack/) which is the 'Rbeast\Source\sfloatMath.f" Fortran file.  `beast_mkl` is generated using Intel's Math Kernel Library (mkl) (https://software.intel.com/en-us/mkl). On average, `beast_mkl` is slightly faster than `beast_default`.

Download "beast_mkl.mexw64" or "beast_default.mexw64" to your local folder, say, "C:\BEAST\".  Use BEAST as follows:

```Matlab
addpath('C:\BEAST\');
beast_default(YOUR_DATA, YOUR_OPTION_PARAMETER)
```

> Detailed examples on how to use BEAST in Matlab are given in the Matlab script files under "Rbeast\Matlab".


---- 
## Additonal Notes

1. Computation
As a Bayesian algorithm, BEAST is fast and is possibly among the fastest implementations of Bayesian time-series analysis algorithms of the same nature. For applications dealing with a few to thousands of time series, the computation won't be an practical concern at all. But for remote sensing applications that may easily involve millions or billions of time series, compuation will be a big challenge for Desktop computer users. We suggest first testing BEAST on singles or small image chips first to determine whether BEAST is appropriate for your applications and, if yes, estimate how long it may take to process the whole image.

2. Complifation from source code
The BEAST source code is complicated than neccessary, mainly because the same source is used for both R and Matlab interfaces as well as for various different compliation settigns (e.g., compilers and library dependencies). The compiliation control variables are defined as MARCOs in abc_marco.h. 
Note that the avove will install "hemiphoto2LAI" from source. Becaues hemiphoto2LAI was written in the mixed use of C/C++ and Fortran. You need to make sure your machine is able to have a C and a Fotran compiler appropriately set up. For example, see [Package Development Prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) for the tools needed for your operating system. In particular, on Windows platforms, the most convenient option is to go with the Rtools toolkit.
 

## Reporting Bugs

BEAST is distributed as is and without warranty of suitability for application. The one distribubuted above is a beta version, with potentail room for further improvement. If you encounter flaws with the software (i.e. bugs) please report the issue. Providing a detailed description of the conditions under which the bug occurred will help to identify the bug. *Use the [Issues tracker](https://github.com/zhaokg/Rbeast/issues) on GitHub to report issues with the software and to request feature enchancements. Alternatively, you can directly email its maintainer Dr. Kaiguang Zhao at lidar.rs@gmail.com.
