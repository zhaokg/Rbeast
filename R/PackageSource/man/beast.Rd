\name{beast}
\docType{methods}
\alias{beast}
\alias{Beast}
\alias{BEAST}
\alias{rbeast}
\alias{Rbeast}
\alias{RBEAST}
\title{Bayesian estimation of abrupt changepoint, nonlinear trend, and periodicity}
\description{
  Apply a Bayesian model averaging algorithm called "BEAST" to decompose time series data into three contrasting components: Abrupt changes, trends, and cyclic/seasonal variations.
  }
\usage{
    beast(data, option=list(),demoGUI=FALSE,...)
}
\arguments{
  \item{data}{a vector or matrix of input data. Missing values are allowed and can be indicated by \code{NA}, \code{NaN}, or a customized value specified in the 2nd argument \code{option} (e.g., \code{option$omittedValue=-9999} ). If \code{data} is a vector of dimension \code{Nx1} or \code{1xN}, it is treated as a single time series of length \code{N}. If \code{data} is a matrix of size \code{NxM}, it is considered as multiple time series: The row dim "\code{N}" is the time-series length; the col dim "\code{M}" is the number of time series. For earth sciences or remote sensing applications, \code{data} can be a 3D array made of stacked time-series images.}
  \item{option}{(optional). If present, \code{option} can be either an INTEGER specifying the known period of the cyclic/seasonal component or a LIST specifying various paramaters for the BEAST algorithm. The "\code{period}" parameter must be an INTEGER specifying the number of samples/values per cycle (e.g, a monthly-sampled time series with an annual period has a period of 12, that is, \code{option$period=12}). Other possible parameters are demonstrated below in Example 3 of the Examples Secction. If \code{option} is absent, BEAST will use default model parameters; in particular, the \code{period} of the cyclic component of time series will be best guessed via auto-correlation before running BEAST.}
  \item{demoGUI}{a boolean indicator. If set to \code{TRUE}, BEAST will be run in a GUI demostration mode, with a GUI window to show an animation of the MCMC sampling in the model space step by step. Note that "\code{demoGUI=TRUE}" works only for Windows x64 systems not Windows 32 or Linux/Mac systems.}
  \item{\dots}{additional parameters, not used currently but reserved for future extension}
}
\value{
  The output is an object of class "beast". It is a list, consisting of the following components. In the explanations below, the input \code{data} is supposed to be of size \code{NxM}: \code{N} is the time series length and \code{M} is the number of time series (\code{M>=1}): 
  
  \item{time}{a vector of size \code{1xN}: the times at the \code{N} sampled locatons. By default, it is simply set to \code{1:N}}
  \item{sN}{a vector of size \code{1xM}. \code{sN} gives the mean number of seasonal changepoints for each of the M time series. If \code{data} is a single time series (i.e., \code{M=1}), \code{sN} will be a scalar.}
  \item{tN}{a vector of size \code{1xM}. \code{tN} gives the mean number of trend changepoints for each of the M time series. If \code{data} is a single time series, \code{tN}  will be a scalar.}
  \item{sNProb}{a matrix of size (\code{option$maxKnotNum_Season+1)xM}. For the i-th time series (i.e.,\code{1=<i<=M}),\code{sNProb[,i] }  gives a probability distribution of having a certain number of changepoints over the range of [0,option$maxKnotNum_Season]; for example, \code{sNProb[1,i]}  is the probability of having no seasonal changepoint in the ith time series; \code{sNProb[3,i]}  is the probability of having 2 changepoints (i.e, 3-1=2).}
  \item{tNProb}{a matrix of size \code{(option$maxKnotNum_Trend+1)xM}. For the i-th time series (i.e.,\code{1=<i<=M}), tNProb[,i] gives a probability distribution of having a certain number of trend changepoints over the range of [0,option$maxKnotNum_Trend]; for example, \code{tNprob[1,i]} is the probability of having no trend changepoint in the ith time series; \code{tNprob[4,i]} is the probability of having 3 changepoints (i.e, 4-1=3).}
  \item{sProb}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{sProb[,i]} gives a probability distribution of having a seasonal changepoint at a certain time over the time range of [1,N]. Plotting \code{sProb[,i]} will depict a curve of probability-of-being-changepoint over the time for the i-th time series.}
  \item{tProb}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{tProb[,i]} gives a probability distribution of having a trend changepoint at a certain time over the time range of [1,N]. Plotting \code{tProb[,i]} will depict a curve of probability-of-being-trend-changepoint over the time for the i-th time series.}
  \item{s}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{s[,i]} gives the best fitted seasonal component.}
  \item{t}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{t[,i]} gives the best fitted trend component.}
  \item{b}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{b[,i]} gives the estimated slope in the fitted trend component over time.}
  \item{sSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{sSD[,i]} gives the standard deviations of the fitted SEASONAL component.}
  \item{tSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{tSD[,i]} gives the standard deviations of the fitted TREND component.}
  \item{bSD}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{bSD[,i]} gives the standard deviations of estimated SLOPE.}  
  \item{sCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{sCI[,i]} gives the 95\% credible intervals of the fitted SEASONAL component.}  
  \item{tCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{tCI[,i]} gives the 95\% credible intervals of the fitted TREND component.} 
  \item{bCI}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{bCI[,i]}  gives the the 95\% credible intervals of the estimated SLOPE.}
  \item{horder}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{horder[,i]}  gives the estimated harmonic order used to approximate the SEASONAL component over time.}
  \item{torder}{a matrix of size \code{NxM}. For the i-th time series (i.e., \code{1=<i<=M}), \code{torder[,i]}  gives the estimated polynomial order used to approximate the TREND component over time.}
  \item{scp}{a matrix of size \code{(option$maxKnotNum_Season)xM}. For the i-th time series (i.e., \code{1=<i<=M}), scp[,i] gives the most possible locations of changepoints in the  the SEASONAL component.}
  \item{tcp}{a matrix of size \code{(option$maxKnotNum_Trend)xM}. For the i-th time series (i.e., \code{1=<i<=M}), tcp[,i] gives gives the most possible locations of changepoints in the TREND component.}
  

}

\examples{

library(Rbeast)

# A MODIS time series of NDVI for a forest plot in Ohio. It has 23 samples 
# per year (i.e., period=23). Note that the input time series to "beast" must
# be spaced/observed at regular time intervals, with missing values indicated
# by NAs, NaNs, or a custimized value (see Example 3). Iregular-sampled time 
# series data need to be first aggegrated at a regular time interval of your
# choice before running beast; if not, beast may give meaningless results or, 
# even worse, terminate abnormally.
data(modis_ohio)
plot(modis_ohio)

\donttest{
#--------------------------------Example 1--------------------------------#
# No "option" argument supplied below so default parameters are used. The period 
# (i.e., 23) will be estimated via auto-correlation. Letting the program compute the
# period of a cyclic time series is not always reliable, so it is always suggested 
# to directly supply the period as in Example 2 and Example 3.
out=beast(modis_ohio)
plot(out)
#*****************************End of Example 1****************************#


#--------------------------------Example 2--------------------------------#
# "option" is set to 23, specicfying the period of modis_ohio as 23
out=beast(modis_ohio,23)
plot(out)
plot(out$s)       #The same as plot(out$s[,1]): plot the seasonal curve
plot(out$sProb)   #Plot the probability of observing seasonal changepoints
plot(out$t)       #The same as plot(out$t[,1]): plot the trend
plot(out$sProb)   #Plot the probability of observing trend changepoints
}
#*****************************End of Example 2****************************#

#--------------------------------Example 3--------------------------------#
# Specify the option parameters explicilty 
opt=list()              #Create an empty list to append individual model parameters

opt$period=23           #Period of the cyclic/seasonal component of the modis time series
opt$minSeasonOrder=2    #Min harmonic order allowed in fitting season component
opt$maxSeasonOrder=8    #Max harmonic order allowed in fititing season component
opt$minTrendOrder=0     #Min polynomial order allowed to fit trend (0 for constant)
opt$maxTrendOrder=1     #Max polynomial order allowed to fit trend (1 for linear term)
opt$minSepDist_Season=20#Min seperation time btw neighboring season change-pts(mustbe >=0)
opt$minSepDist_Trend=20 #Min seperation time btw neighboring trend change-pts(must be >=0)
opt$maxKnotNum_Season=4 #Max number of season changepoints allowed 
opt$maxKnotNum_Trend=10 #Max number of trend changepoints allowed  
opt$omittedValue=-999   #A customized value to indicate bad/missing values in the time 
                        #series, in additon to those NA or NaN values.

opt$printToScreen=1     #If set to 1, display some progress status while running
opt$printCharLen=150    #The length of chars in each status line when printToScreen=1
						
# The following parameters used to configure the reverisible-jump MCMC (RJMCC) sampler
opt$chainNumber=2       #Number of parallel MCMC chains 
opt$sample=1000         #Number of samples to be collected per chain
opt$thinningFactor=3    #A factor to thin chains (e.g., samples taken every 3 iterations)
opt$burnin=500          #Number of burn-in samples discarded at the start of each chain
opt$maxMoveStepSize=30  #For the move proposal, the max window allowed in jumping from 
                        #the current changepoint
opt$resamplingSeasonOrderProb=0.2 #The probability of selecting a re-sampling proposal 
                                  #(e.g., resample seasonal harmonic order)
opt$resamplingTrendOrderProb=0.2  #The probability of selecting a re-sampling proposal 
                                  #(e.g., resample trend polynomial order)
								   
opt$seed=65654          #A seed for the random generator: If seed=0,random numbers differ
                        #for different BEAST runs. Setting seed to a chosen non-zero integer 
                        #will allow reproducing the same result for different BEAST runs.
opt$computeCredible=0   #If set to 1, compute 95% credible intervals: The results will be 
                        #saved as sCI, tCI, and bCI in the output variable.
opt$fastCIComputation=0 #If set to 1, employ a fast algorithm to compute credible intervals 
opt$computeSlopeSign=1  #If set to 1, compute the probability of having a postive slope in 
                        #the estimated trend. The result will be saved as bsign in the output
                        #variable.  

opt$computeHarmonicOrder=1 #If set to 1, compute the mean harmonic order of the fitted 
                           #seasonal component. The result will be saved as "horder" in 
                           #the output variable.
opt$computeTrendOrder=1    #If set to 1, compute the mean polynomial order of the fitted
                           #trend component. The result will be saved as "torder" in 
                           #the output variable.


#opt$outputToDisk=0        #(if set to 1, results will be written to files in a folder)
#opt$outputFolder ='c:/out'#Specify the output folder when outputToDisk=1     
#opt$lengthPerTimeSeries_infile=300#the time series length if input data come from a binary file

\donttest{
# Use "opt" defined above in the beast function. Note that to run beast(), not all the individual 
# parameters in option need to be explicitly specified. If an parameter is not given in option, its 
# default value will be used.
out=beast(modis_ohio, opt)
plot(out)
}
#*****************************End of Example 3****************************#

\donttest{
#--------------------------------Example 4--------------------------------#
# Run an interactive GUI to visualize how BEAST is samplinig from
# the possible model spaces in terms of the numbers and timings of 
# seasonal and trend changepoints. The GUI inferface allows changing
# the option parameters interactively. This GUI is only available on 
# Win x64 machines, not Mac or Linux.
beast(modis_ohio, 23, demoGUI=TRUE)
#*****************************End of Example 4****************************#

#--------------------------------Example 5--------------------------------#
# 'simdata' is a 300x3 matrix, consisting of three time series
data(simdata)

# Plot individual time series. As a toy example, all the three time series
# are the same.
plot(simdata[,1]) 
plot(simdata[,2]) 

# Below, the option is defined in the command line as a temporary list.
out=beast( simdata, list(period=24, chainNumber=3, sample=1000, burnin=200) )

# "out" contains results for the three time series. Plot the result for the second one
plot(out,2)
#*****************************End of Example 5****************************#
}

}
\references{
Zhao K, Wulder M, Hu, T, Bright R, Wu Q, Qin H, Li Y, Toman E, Mallick B, Zhang X, and Brown M (2019).  Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment (accepted manuscript) 
}
\keyword{misc}
