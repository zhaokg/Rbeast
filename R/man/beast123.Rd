\name{beast123}
\docType{methods}
\alias{beast123}
\alias{BEAST123}
\alias{beast.123}
 
\title{Bayesian time series decomposition for changepoint, trend, and periodicity or seasonality}
\description{
   A Bayesian model averaging algorithm called BEAST to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations. BEAST is useful for changepoint detection (e.g., breakpoints or structural breaks), nonlinear trend analysis, time series decomposition, and time series segmentation.
  }
\usage{
    beast123( Y, 
              metadata=list(),
              prior   =list(),
              mcmc    =list(),
              extra   =list(),
              season  = c('harmonic','svd','dummy','none'),
              method  = c( 'bayes', 'bic', 'aic', 'aicc','hic'),
              ...)
}
\arguments{

  \item{Y}{ 
  	a 1D vector, 2D matrix, or 3D array of numeric data. Missing values are allowed and can be indicated by \code{NA}, \code{NaN}, or a value customized in the 2nd argument \code{metadata} (e.g., \code{metadata$missingValue=-9999}).	 
	\itemize{
    	\item If \code{Y} is a vector of size \code{Nx1} or \code{1xN}, it is treated as a single time series of length \code{N}.
		
		\item If \code{Y} is a 2D matrix or 3D array of dimension \code{N1xN2} or \code{N1xN2xN3} (e.g., stacked images of geospatial data), it includes multiple time series of equal length: Which dimension is time has to be specified in the 2nd argument using \code{metadata$whichDimIsTime}. For example, \code{metadata$whichDimIsTime} = 1 for a 190x35 2D input indicates 35 time series of length 190 each; \code{metadata$whichDimIsTime} = 2 for a 100x200x300 3D input indicates 30000=100*300 time series of length 200 each.
    } 	
    \code{Y} can be either regular (i.e., evenly-spaced in time) or irregular/unordered in time.
  	\itemize{
    	\item If regular, individual times are determined from the time of the 1st data point \code{startTime} and the time span between consecutive points \code{deltaTime}, which are specified in the 2nd arg through \code{metadata$startTime and metadata$deltaTime}; if not given, \code{startTime} and \code{deltaTime} take a default 1.0.
		\item  If irregular or regular but unordered, the times have to be explicitly given through \code{metadata$time}. The BEAST model is currently formulated for regular data only, so internally, the \code{beast123} function will aggregate/re-bin irregular data into regular ones; for the aggregation, the  \code{metadata$deltaTime} parameter should also be also provided to specify the desired bin size or time interval.
	}

    \code{Y} can have a periodic component or have a trend component only. Use the argument \code{season} to specify the cases.
  	\itemize{
   		\item  \code{season='none'}: \code{Y} is treated as trend-only; no periodic components are present in the time series.
		\item  \code{season='harmonic'}: \code{Y} has a periodic/seasonal component. The term 'season' is a misnomer being used here to broad refer to any periodic variations present in \code{Y}. The periodicity is not a statistical parameter estimated by BEAST but a known constant given by the user through \code{metadata$freq}. The periodic component is modeled as a harmonic curve--a combination of sins and cosines.
		\item  \code{season='dummy'}: the same as 'harmonic' except that the periodic/seasonal component is modeled as a non-parametric curve.
		\item  \code{season='svd'}: (experimental feature) the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a linear combination of function bases derived from a Single-value decomposition. The SVD-based basis functions are more parsimonious than the harmonic sin/cos bases in parameterizing the seasonal variations; therefore, more subtle changepoints are likely to be detected.	
		
    }
  }

  \item{metadata}{

	(optional). If present, \code{metadata} may (1) a scalar value to specify the period of the input Y, (2) a vector of numbers, strings, or R Dates to specify the times of Y, or (3) more often, a LIST object specifying various parameters to describe the 1st argument \code{Y}. If missing, default values will be used. But \code{metadata} should be explicitly provided if the input \code{Y} is a 2D matrix or 3D array to avoid mis-interpreting the input Y. \code{metadata} is not part of BEAST's Bayesian formulation but just some additional info to interpret \code{Y}. If \code{metadata} is provided as a LIST, below are possible fields; not all of them are always needed, depending on the types of inputs (e.g., 1D, 2D or 3D; regular or irregular). 
	
         	\itemize{
    		\item \code{metadata$whichDimIsTime}:  integer (<=3). Needed to specify which dimension of \code{Y} is time for a matrix or 3D array input. Ignored if the input \code{Y} is a vector.
    		\item \code{metadata$isRegularOrdered}: logical. Obsolete and no longer used in this version.  Now, \code{metadata$time} is analyzed to determine whether the input is irregular or not; if \code{metadata$time} is missing,  \code{Y} is assumed to be regular.
    		\item \code{metadata$time}: a vector of the same length as \code{Y}'s time dimension to provide the times for datapoints. It can be a vector of numbers, Dates, or date strings; it can also be a list of vectors of year, months, and days. Possible formats include:
      		 
			   \enumerate{
					\item a vector of numerical values [e.g., c(1984.23, 1984.27, 1984.36, ...)]. The unit of the times is irrelevant to BEAST as long as it is consistent with the unit used for specifying \code{startTime}, \code{deltaTime}, and \code{period}.
					
					\item a vector of R Dates [e.g., \code{as.Date}( c("1984-03-27", "1984-04-10", "1984-05-12",... )].
					
					\item a vector of char strings. Examples are:
					   \itemize{
						 \item  c("1984-03-27", "1984-04-10", "1984-05-12")
						 \item  c("1984/03/27", "1984,04,10", "1984 05 12") (i.e., the delimiters differ as long as the YMD order is consistent)
						 \item  c("LT4-1984-03-27",  "LT4-1984-04-10",  "LT4-1984+05,12")
						 \item  c("LT4-1984087ndvi", "LT4-1984101ndvi", "LT4-1984133ndvi")
						 \item  c("1984,,abc 3/ 27", "1984,,ddxfdd 4/ 10"  "ggd1984,, 5/ ttt 12")
					   }
					   BEAST uses several heuristics to automatically parse the date strings without a format specifier but may fail due to ambiguity (e.g., in "LC8-2020-09-20-1984", no way to tell if 2020 or 1984 is the year). To ensure correctness, use a list object as explained below to provide a date format specifier.
					
					\item a list object \code{time=list(datestr=..., strfmat='...')} consisting of a vector of date strings (\code{time$datestr}) and a format specifier (\code{time$strFmt}). The string \code{time$strFmt} specifies how to parse \code{dateStr}. Three formats are currently supported:
					
							\itemize{
								\item (a). All the date strings have a fixed pattern in terms of the relative positions of Year, Month, and Day. For example, to extract 2001/12/02 etc from \code{time$dateStr} = c('P23R34-2001.1202333xd', 'O93X94-2002.1108133fd', 'TP3R34-2009.0122333td') use \code{time$strFmt}=\code{'P23R34-yyyy.mmdd333xd'} where \code{yyyy}, \code{mm}, and \code{dd} are the specifiers and other positions are wildcards and can be filled with any other letters different from \code{yyyy, mm and dd}.
								\item (b). All the date strings have a fixed pattern in terms of the relative positions of year and doy. For example, to extract 2001/045(day of year) from 'P23R342001888045', use strFmt='123123yyyy888doy' where \code{yyyy} and \code{doy} are the specifiers and other positions are wildcards and can be filled with any other letters different from yyyy, and doy. 'doy' must be three digit in length.
								\item (c). All the date strings have a fixed pattern in terms of the separation characters between year, month, and day. For example, to extract 2002/12/02 from '2002,12/02', '  2002 , 12/2', '2002,12  /02 ', use strFmt='Y,M/D' where the whitespaces are ignored. To get 2002/12/02 from '2--12, 2012 ', use strmFmt='D--M,Y'.	
							 }  
							 
				   \item a list object of vectors to specify individual dates of the time series. Use \code{time$year},\code{time$month},and \code{time$day} to give the dates; or alternatively use \code{time$year} and \code{time$doy} where each value of the \code{doy} vector is a number within 1 and 365/366. Each vector must have the same length as the time dimension of \code{Y}.

			} %\enumerate
		  
			
           \item \code{metadata$startTime}: numeric (default to 1.0 if missing). It gives the time of the 1st data point. It can be specified as a scalar (e.g., 2021.23) or a vector of three values in the order of year, month, and day (e.g., \code{metadata$startTime} = \code{c(2021,1,24)}). \code{metadata$startTime} is needed for regular input data but optional for irregular data: If missing, startTime will be computed from \code{metadata$time} for irregular \code{Y}.
           \item \code{metadata$deltaTime}: numeric or string. It specifies the time interval between consecutive data points. It is optional for regular data (default to 1.0 if not supplied), but should be specified for irregular data because deltaTime is needed to aggregate/resample the irregular time series into regular ones.  The unit of \code{deltaTime} needs to be consistent with \code{metadata$time}. If \code{metadata$time} takes a numeric vector, the unit of \code{deltaTime} is arbitrary and irrelevant to BEAST. If \code{time} takes a vector of Dates or date strings, the unit for \code{deltaTime} is assumed to Fractional YEAR. If needed, use a string instead of a number to specify whether the unit of deltaTime is day, month, or year. Examples include '7 days', '7d', '1/2 months', '1mn', '1.0 year', and '1y'.
			   
    	   \item \code{metadata$period}:  numeric or string. Specify the period for the periodic/seasonal component in \code{Y}. Needed only for data with a periodic/cyclic component (i.e., season=\code{'harmonic'} or \code{'dummy'}) and not used for trend-only data (i.e., \code{season='none'}).  The period of the cyclic component should have a unit consisent with the unit of \code{deltaTime}. It holds that \code{period=deltaTime*freq} where \code{freq} is the number of data samples per period. (Note that the \code{freq} argument in earlier versions becomes obsolete and now is replaced by \code{period}.)  \code{period} or the number of data points per period is not a BEAST model parameter and it has to be specified by the user. But if \code{period} is missing, BEAST first attempts to guess its value via auto-correlation before fitting the model. If \code{period} <= 0, no seasonal/cyclic component is assumed (i.e, \code{season='none'}) and the trend-only model is used.  If needed, use a string to specify whether the unit of period is day, month, or year. Examples are '1.0 year', '12 months', '365d', '366 days'.
    	   \item \code{metadata$missingValue}: numeric; a customized value to indicate bad/missing values in the time series, in addition to those NA or NaN values.
    	   \item \code{metadata$maxMissingRate} a fractional number within [0, 1] as the maximum percentage of missing values, above which the time series will be skipped and won't be fitted by BEAST.
           } %\itemize{


   }
  

 \item{prior}{
	(optional). a list object consisting of the hyperprior parameters in the Bayesian formulation of the BEAST model. Because they are part of the model, the fitting result may be sensitive to the choices of these hyperparameters. If \code{prior} is missing, a set of default values will be used and the exact values used will be printed to the console at the start of the BEAST run. Below are possible parameters:
    \itemize{
    	\item \code{prior$seasonMinOrder}: integer (>=1)
    	\item \code{prior$seasonMaxOrder}: integer (>=1); the min and max harmonic orders considered to fit the seasonal component. \code{seasonMinOrder} and \code{seasonMaxOrder} are only used if the time series has a seasonal component (i.e., \code{season='harmonic'}) and ignored for trend-only data or when \code{season='dummy'}. If \code{seasonMinOrder}=\code{seasonMaxOrder}, BEAST assumes a constant harmonic order used and won't infer the posterior probability of harmonic orders.
  	    \item \code{prior$seasonMinKnotNum}: integer (>=0)
    	\item \code{prior$seasonMaxKnotNum}: integer (>=0); the min and max number of seasonal changepoints allowed in segmenting and fitting the seasonal component. \code{seasonMinKnotNum} and \code{seasonMaxKnotNum} are only used if the time series has a seasonal component (i.e., \code{season='harmonic'} or \code{season='dummy'}) and ignored for trend-only data. If \code{seasonMinOrder}=\code{seasonMaxOrder}, BEAST assumes a constant number of changepoints and won't infer the posterior probability of the number of changepoints, but it will still estimate the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur). If \code{seasonMinOrder}=\code{seasonMaxOrder}=0, no changepoints are allowed in the seasonal component; then a global harmonic model is used to fit the seasonal component.
	    \item \code{prior$seasonMinSepDist}: integer (>0). the min separation time between two neighboring season changepoints. That is, when fitting a piecewise harmonic seasonal model, no two changepoints are allowed to occur within a time window of seasonMinSepDist. \code{seasonMinSepDist} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is seasonMinSepDist*metadata$deltaTime.
        \item \code{prior$seasonLeftMargin}:  integer (>=0); the number of leftmost data points excluded for seasonal changepoint detection. That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the starting window/segment of length \code{seasonLeftMargin}. \code{seasonLeftMargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{seasonLeftMargin*deltat}. If missing, \code{seasonLeftMargin} defaults to \code{seasonMinSepDist}.
        \item \code{prior$seasonRightMargin}: integer (>=0); the number of rightmost data points excluded for seasonal changepoint detection. That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the ending window/segment of length \code{seasonRightMargin}. \code{seasonRightMargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{seasonRightMargin*deltat}. If missing, \code{seasonRightMargin} defaults to \code{seasonMinSepDist}.  
  	    \item \code{prior$trendMinOrder}:  integer (>=0)
    	\item \code{prior$trendMaxOrder}:  integer (>=0); the min and max orders of the polynomials considered to fit the trend component. The zero-th order corresponds to a constant term/ a flat line and the 1st order is a line. If \code{trendMinOrder}=\code{trendMaxOrder}, BEAST assumes a constant polynomial order used and won't infer the posterior probability of polynomial orders.
  	    \item \code{prior$trendMinKnotNum}: 
    	\item \code{prior$trendMaxKnotNum}:  integer (>=0); the min and max number of trend changepoints allowed in segmenting and fitting the trend component. If \code{trendMinOrder}=\code{trendMaxOrder}, BEAST assumes a constant number of changepoints in the fitted trend and won't infer the posterior probability of the number of trend changepoints, but it will still estimate the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur). If \code{trendMinOrder}=\code{trendMaxOrder}=0, no changepoints are allowed in the trend component; then a global polynomial model is used to fit the trend.
	    \item \code{prior$trendMinSepDist}:  integer (>0). the min separation time between two neighboring trend changepoints.
		\item \code{prior$trendLeftMargin}:  integer (>=0); the number of leftmost data points excluded for trend changepoint detection. That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the starting window/segment of length \code{trendLeftMargin}. \code{trendLeftMargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{trendLeftMargin*deltat}. If missing, \code{trendLeftMargin} defaults to \code{trendMinSepDist}. 
        \item \code{prior$trendRightMargin}: integer (>=0); the number of rightmost data points excluded for trend changepoint detection. That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the ending window/segment of length \code{trendRightMargin}. \code{trendRightMargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{trendRightMargin*deltat}. If missing, \code{trendRightMargin} defaults to \code{trendMinSepDist}.
		\item \code{prior$precValue}:     numeric (>0); the default value is 10. Useful only if prior$precPriorType='constant'
	    \item \code{prior$precPriorType}: characters. It takes one of 'constant', 'uniform' (the default),  'componentwise', and 'orderwise'. Below are the differences between them.
			\enumerate{
                    \item \code{precPriorType='constant'}:  the precision parameter used to parameterize the model coefficients is fixed to a constant specified by \code{prior$precValue}. In other words, \code{prior$precValue} is a user-defined hyperparameter and the fitting result may be sensitive to the chosen values of \code{prior$precValue}.
                    \item \code{precPriorType='uniform'}:  the precision parameter used to parameterize the model coefficients is a random variable; its initial value is specified by \code{prior$precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result is insensitive to the choice in \code{prior$precValue}.
                    \item \code{precPriorType='componentwise'}: multiple precision parameters are used to parameterize the model coefficients for individual components (e.g., one for season and another for trend); their initial values is specified by \code{prior$precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result is insensitive to the choice in \code{prior$precValue}.
                    \item \code{precPriorType='orderwise'}: multiple precision parameters are used to parameterize the model coefficients not just for individual components but also for individual orders of each component; their initial values is specified by \code{prior$precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result is insensitive to the choice in \code{prior$precValue}.
             } % enumerate

    } %\itemize

   }

  

   \item{mcmc}{

	(optional). a list object consisting of parameters to configure the MCMC inference. These parameter are not part of the Bayesian formulation of the BEAST model but are the settings for the reversible-jump MCMC to generate MCMC chains. Due to the MCMC nature, the longer the simulation chain is, the better the fitting result. Below are possible parameters:

           	\itemize{
    			\item \code{mcmc$seed}: integer (>=0); the seed for the random number generator.  If \code{mcmc$seed=0}, an arbitrary seed will be picked up and the fitting result will var across runs. If fixed to the same on-zero integer, the results can be re-produced for different runs. Note that the results may still vary if run on different computers with the same seed because the random generator library depends on CPU's instruction sets.
    			\item \code{mcmc$samples}: integer (>0); the number of samples collected per MCMC chain.
    			\item \code{mcmc$chainNumber}: integer (>0); the number of parallel MCMC chains.
    			\item \code{mcmc$thinningFactor}: integer (>0); a factor to thin chains (e.g., if thinningFactor=5, samples will be taken every 3 iterations).
    			\item \code{mcmc$burnin}: integer (>0); the number of burn-in samples discarded at the start of each chain.
    			\item \code{mcmc$maxMoveStepSize}: integer (>0). The RJMCMC sampler employs a move proposal when traversing the model space or proposing new positions of changepoints. 'maxMoveStepSize' is used in the move proposal to specify the max window allowed in jumping from the current changepoint.
    			\item \code{mcmc$seasonResamplingOrderProb}: a fractional number less than 1.0; the probability of selecting a re-sampling proposal (e.g., resample seasonal harmonic order).
    			\item \code{mcmc$trendResamplingOrderProb}:  a fractional number less than 1.0; the probability of selecting a re-sampling proposal (e.g., resample trend polynomial order)
  			\item \code{mcmc$credIntervalAlphaLevel}: a fractional number less than 1.0 (default to 0.95); the level of confidence used to compute credible intervals.  			
                 }
  }  

 
   \item{extra}{

	(optional). a list object consisting of flags to control the outputs from the BEAST runs or configure other program setting. Below are possible parameters:

        \itemize{
		    \item \code{extra$quiet}: logical (default to FALSE). If TRUE, no warning messages will be printed out. 
    		\item \code{extra$dumpInputData}: logical (default to FALSE). If TRUE, the input time series will be copied into the output. When the input \code{Y} is irregular (i.e., \code{metadata$isRegularOrdered=FALSE}), the dumped copies will be the aggregated regular time series. 
    		\item \code{extra$whichOutputDimIsTime}: integer (<=3). If the input \code{Y} is a 2D or 3D array (i.e., multiple time series such as stacked images), the \code{whichOutputDimIsTime} specifies which dimension is the time in the output variables. \code{whichOutputDimIsTime} defaults to 3 for 3D inputs and is ignored if the input is a vector (i.e., a single time series).	
			
			\item \code{extra$ncpStatMethod}: character (deprecated). A string to specify which statistic is used to determine the Number of ChangePoint (ncp) when computing the most likely changepoint locations (e.g., out$trend$cp, and out$season$cp). Three values are possible: 'mode', 'mean', and 'median'; the default is 'mode'. Individual models sampled by BEAST has a varying dimension (e.g., number of changepoints or knots). For example, if mcmc$samples=10, the numbers of changepoints for the 10 sampled models are assumed to be c(0, 2, 4, 1, 1, 2, 7, 6, 6, 1). The mean ncp is 3.1 (rounded to 3), the median is 2.5 (2), and the mode is 1. This argument is deprecated; now all the possible changepoints are outputted, together with several versions of ncp, including \code{ncp}, \code{ncp_median}, \code{ncp_mode}, and \code{ncp_pct90}. A similar parameter \code{ncpStat} is added to the \code{\link[Rbeast]{plot.beast}} function to specify which ncp is used when plotting. 
			
    		\item \code{extra$computeCredible}: logical (default to TRUE). Credible intervals will be computed and outputted only if set to TRUE.
    		\item \code{extra$fastCIComputation}:  logical (default to TRUE). If TRUE, a fast method is used to compute credible intervals (CI). Computation of CI is one of the most computational parts and fastCIComputation should be set to TRUE unless more accurate CI estimation is desired.  
    		\item \code{extra$computeSeasonOrder}: logical (default to TRUE). If TRUE, a posterior estimate of the seasonal harmonic order will be outputted; this flag is only valid if the time series has a seasonal component (i.e., season='harmonic' and prior$seasonMinOrder is not equal to prior$seasonMaxOrder).
    		\item \code{extra$computeTrendOrder}: logical (default to TRUE). If TRUE, a posterior estimate of the tend polynomial order will be outputted; this flag is only valid when prior$trendMinOrder is not equal to prior$trendMaxOrder).
  			\item \code{extra$computeTrendOrder}: logical (default to TRUE). If TRUE, a posterior estimate of the tend polynomial order will be outputted; this flag is only valid when prior$trendMinOrder is not equal to prior$trendMaxOrder).
  			\item \code{extra$computeSeasonChngpt}: logical (default to TRUE). If TRUE, compute the most likely times/positions where changepoints occur in the seasonal component. This flag is not valid if there is a seasonal component in the time series  (i.e., season='harmonic' or season='dummy' and prior$seasonMaxKnotNum is non-zero).
  			\item \code{extra$computeTrendChngpt}: logical (default to TRUE). If TRUE, compute the most likely times/positions where changepoints occur in the trend component.
  			\item \code{extra$computeSeasonAmp}: logical (default to FALSE). If TRUE, compute and output the time-varying amplitude of the seasonality.
  			\item \code{extra$computeTrendSlope}: logical (default to FALSE). If TRUE, compute and output the time-varying slope of the estimated trend.
  			\item \code{extra$tallyPosNegSeasonJump}: logical (default to FALSE). If TRUE, compute and differentiate seasonal changepoints in terms of the direction of the jumps in the estimated seasonal signal. Those changepoints with a positive jump will be outputted separately from those with a negative jump. A series of output variables (some for positive-jump changepoints, and others for negative-jump changepoints will be dumped).
  			\item \code{extra$tallyPosNegTrendJump}: logical (default to FALSE). If TRUE, compute and differentiate trend changepoints in terms of the direction of the jumps in the estimated trend. Those changepoints with a positive jump will be outputted separately from those with a negative jump. A series of output variables (some for positive-jump changepoints, and others for negative-jump changepoints will be dumped).
  	                \item \code{extra$tallyIncDecTrendJump}: logical (default to FALSE). If TRUE, compute and differentiate trend changepoints in terms of the direction of the jumps in the estimated slope of the trend signal. Those changepoints with a increase in the slope will be outputted separately from those with a decrease in the slope. A series of output variables (some for increase-jump changepoints, and others for decrease-jump changepoints will be dumped).
  			\item \code{extra$printProgressBar}: logical (default to FALSE). If TRUE, a progress bar will be displayed to show the status of the running. When running on multiple time series (e.g. stacked image time series), the progress bar will also report an estimate of the remaining time for completion.
  			\item \code{extra$consoleWidth}: integer (default to 0); the length of chars in each status line when setting printProgressBar=TRUE. If 0, the current width of the console will be used.
  			\item \code{extra$printOptions}: logical (default to FALSE). If TRUE, the values used in the arguments \code{metadata}, \code{prior}, \code{mcmc}, and \code{extra} will be printed to the console at the start of the run.
     			\item \code{extra$numThreadsPerCPU}: integer (default to 2); the number of threads to be scheduled for each CPU core. 
     			\item \code{extra$numParThreads}: integer (default to 0). When handling many time series, BEAST can use multiple concurrent threads. \code{extra$numParThreads} specifies how many concurrent threads will be used in total. If numParThreads=0, the actual number of threads will be \code{numThreadsPerCPU * cpuCoreNumber}; that is, each CPU core will generate a number 'numThreadsPerCPU' of threads. On Windows 64, ,BEAST is group-aware and will affine or distribute the threads to all the NUMA node. But currently, up to 256 CPU cores are supported.
			}
  }  


  \item{season}{ characters (default to 'harmonic'); specify if \code{y} has a periodic component or not. Three strings are possible.
				\itemize{
					\item  \code{'none'}: \code{y} is trend-only; no periodic components are present in the time series. The args for the seasonal component (i.e.,\code{sorder.minmax}, \code{scp.minmax} and \code{sseg.max}) will be irrelevant and ignored. 
					\item  \code{'harmonic'}: \code{y} has a periodic/seasonal component. The term \code{season} is a misnomer, being used here to broadly refer to any periodic variations present in \code{y}. The periodicity is NOT a model parameter estimated by BEAST but a known constant given by the user through \code{freq}. By default, the periodic component is modeled as a harmonic curve--a combination of sins and cosines.
					\item  \code{'dummy'}: the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a non-parametric curve. The harmonic order arg \code{sorder.minmax} is irrelevant and is ignored.
					\item  \code{'svd'}: (experimental feature) the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a linear combination of function bases derived from a Single-value decomposition. The SVD-based basis functions are more parsimonious than the harmonic sin/cos bases in parameterizing the seasonal variations; therefore, more subtle changepoints are likely to be detected.
				}  
  }				
 
	
  \item{method}{ a string (default to 'bayes'); specify the method for formulating model posterior probability.
         \itemize{
				  \item  \code{'bayes'}: the full Bayesian formulation as described in Zhao et al. (2019). 					
				  \item  \code{'bic'}:  approximation of posterior probability using the Bayesian information criterion (bic).					
				  \item  \code{'aic'}:  approximation of posterior probability using the Akaike information criterion (aic).					
 				  \item  \code{'aicc'}:  approximation of posterior probability using the corrected Akaike information criterion (aicc).					
 				  \item  \code{'hic'}:  approximation of  posterior probability using the Hannan-Quinn information criterion  (hic)					
				}    		   
 }
 				

  \item{\dots}{additional parameters, not used currently but reserved for future extension}

}

\value{
  The output is an object of class "beast". It is a list, consisting of the following variables. Exact sizes of the variables depend on the types of the input \code{Y} as well as the specified output time dimension \code{extra$whichOutputDimIsTime}. In the explanations below, we assume the input \code{Y} is a single time series of length \code{N};  the dimensions for 2D or 2D inputs may be interpreted accordingly:

  \item{time}{a vector of size \code{1xN}: the times at the \code{N} sampled locations. By default, it is simply set to \code{1:N} if the input arguments \code{metadata$deltaTime}, \code{metadata$startTime}, or \code{metadata$time} are missing.}
  \item{data}{a vector, matrix, or 3D array; this is a copy of the input \code{Y} if \code{extra$dumpInputData} = \code{TRUE}. If extra$dumpInputData=FALSE, it is set to NULL. If the original input \code{Y} is irregular, the copy here is the regular version aggregated from the original at the time interval specified by metadata$deltaTime. }
  \item{marg_lik}{numeric; the average of the model marginal likelihood; the larger marg_lik, the better the fitting for a given time series.}
  \item{R2}{numeric; the R-square of the model fitting.}
  \item{RMSE}{numeric; the RMSE of the model fitting.}
  \item{sig2}{numeric; the estimated variance of the model error.}
  \item{trend}{a list object numeric consisting of various outputs related to the estimated trend component:
  	\itemize{
		 \item \code{ncp}: [Number of ChangePoints]. a numeric scalar; the mean number of trend changepoints. Individual models sampled by BEAST has a varying dimension (e.g., number of changepoints or knots), so several alternative statistics (e.g., \code{ncp_mode}, \code{ncp_median}, and \code{ncp_pct90}) are also given to summarize the number of changepoints. For example, if \code{mcmc$samples=10}, the numbers of changepoints for the 10 sampled models are assumed to be c(0, 2, 4, 1, 1, 2, 7, 6, 6, 1). The mean ncp is 3.1 (rounded to 3), the median is 2.5 (2),  the mode is 1, and the 90th percentile (ncp_pct90) is 6.5. 
		 \item \code{ncp_mode}: [Number of ChangePoints]. a numeric scalar; the mode for number of changepoints. See the above for explanations.
		 \item \code{ncp_median}: [Number of ChangePoints]. a numeric scalar; the median for number of changepoints. See the above for explanations.
		 \item \code{ncp_pct90}: [Number of ChangePoints]. a numeric scalar; the 90th percentile for number of changepoints. See the above for explanations.
		 \item \code{ncpPr}: [Probability of the Number of ChangePoints]. A vector of length \code{(prior$trendMaxKnotNum+1)}. It gives a probability distribution of having a certain number of trend changepoints over the range of [0,prior$trendMaxKnotNum]; for example, \code{ncpPr[1]} is the probability of having no trend changepoint; \code{ncpPr[i]} is the probability of having (i-1) changepoints: Note that it is \code{ncpPr[i]} not \code{ncpPr[i-1]} because ncpPr[1] is used for having zero changepoint.
		 \item \code{cpOccPr}: [ChangePoint OCCurence PRobability]. a vector of length N; it gives a probability distribution of having a changepoint in the trend at each point of time. Plotting \code{cpOccPr} will depict a continious curve of probability-of-being-changepoint. Of particular note, in the curve, a higher peak indicates a higher chance of being a changepoint only at that particular SINGLE point in time and does not necessarily mean a higher chance of observing a changepoint AROUND that time. For example, a window of cpOccPr values \code{c(0,0,0.5,0,0)} (i.e., the peak prob is 0.5 and the summed prob is 0.5) is less likely to be a changepoint compared to another window \code{c(0.1,0.2,0.21,0.2,0.1)} (i.e., the peak prob is 0.21 but the summed prob is 0.71).
		 \item \code{order}: a vector of length N; the average polynomial order needed to approximate the fitted trend. As an average over many sampled individual piece-wise polynomial trends, \code{order} is not necessarily an integer.

		 \item \code{cp}: [Changepoints] a vector of length \code{tcp.max=tcp.minmax[2]}; the most possible changepoint locations in the trend component. The locations are obtained by first applying a sum-filtering to the \code{cpOccPr} curve with a filter window size of \code{tseg.min} and then picking up to a total \code{prior$MaxKnotNum}/\code{tcp.max} of the highest peaks in the filtered curve.  NaNs are possible if no enough changepoints are identified. \code{cp} records all the possible changepoints identified and many of them are bound to be false positives. Do not blindly treat all of them as actual changepoints.
		 
		 \item \code{cpPr}: [Changepoints PRobability] a vector of length \code{metadata$trendMaxKnotNum}; the probabilities associated with the changepoints \code{cp}. Filled with NaNs for the remaining elements if \code{ncp<trendMaxKnotNum}.
		 \item \code{cpCI}: [Changepoints Credible Interval] a matrix of dimension \code{metadata$trendMaxKnotNum x 2}; the credibable intervals for the detected changepoints \code{cp}.
		 \item \code{cpAbruptChange}: [Abrupt change at Changepoints]  a vector of length \code{metadata$trendMaxKnotNum}; the jumps in the fitted trend curves at the detected changepoints \code{cp}.
		 \item \code{Y}: a vector of length N; the estimated trend component. It is the Bayesian model averaging of all the individual sampled trend.
		 \item \code{SD}: [Standard Deviation] a vector of length N; the estimated standard deviation of the estimated trend component.
		 \item \code{CI}: [Standard Deviation] a matrix of dimension \code{N x 2}; the estimated credible interval of the estimated trend. One vector of the matrix is for the upper envelope and another for the lower envelope.
		 \item \code{slp}: [Slope] a vector of length N;  the time-varying slope of the fitted trend component .
		 \item \code{slpSD}: [Standar Deviation of Slope] a vector of length N;  the SD of the slope for the trend component.
		 \item \code{slpSgnPosPr}: [PRobability of slope having a positive sign] a vector of length N;  the probability of the slope being positive (i.e., increasing trend) for the trend component. For example, if \code{slpSgnPosPr=0.80} at a given point in time, it means that 80\% of the individual trend models sampled in the MCMC chain has a positive slope at that point.
		 \item \code{slpSgnZeroPr}: [PRobability of slope being zero] a vector of length N;  the probability of the slope being zero (i.e., a flat constant line) for the trend component. For example, if \code{slpSgnZeroPr=0.10} at a given point in time, it means that 10\% of the individual trend models sampled in the MCMC chain has a zero slope at that point.  The probability of slope being negative can be obtained from \code{1}-\code{slpSgnZeroPr}-\code{slpSgnPosPr}.

		 \item \code{pos_ncp}: 
		 \item \code{neg_ncp}:  
		 \item \code{pos_ncpPr}:
		 \item \code{neg_ncpPr}:
		 \item \code{pos_cpOccPr}:
		 \item \code{neg_cpOccPr}:
		 \item \code{pos_cp}:
		 \item \code{neg_cp}:
		 \item \code{pos_cpPr}:
		 \item \code{neg_cpPr}:
		 \item \code{pos_cpAbruptChange}:
		 \item \code{neg_cpAbruptChange}:
		 \item \code{pos_cpCI}:
		 \item \code{neg_cpCI}:  The above variables have the same outputs as those variables without the prefix 'pos' and 'neg', except that we differentiate the changepoints with a POStive jump in the trend from those changepoints with a NEGative jump. For example, \code{pos_ncp} refers to the average number of trend changepoints that jump up (i.e., positively) in the trend.
		 \item \code{inc_ncp}: 
		 \item \code{dec_ncp}:  
		 \item \code{inc_ncpPr}:
		 \item \code{dec_ncpPr}:
		 \item \code{inc_cpOccPr}:
		 \item \code{dec_cpOccPr}:
		 \item \code{inc_cp}:
		 \item \code{dec_cp}:
		 \item \code{inc_cpPr}:
		 \item \code{dec_cpPr}:
		 \item \code{inc_cpAbruptChange}:
		 \item \code{dec_cpAbruptChange}:
		 \item \code{inc_cpCI}:
		 \item \code{dec_cpCI}:  The above variables have the same outputs as those variables without the prefix 'inc' and 'dec', except that we differentiate the changepoints at which the trend slope increases from those changepoints at which the trend slope decreases. For example, if the trend slopes before and after a chngpt is 0.4 and 2.5, then the changepoint is counted toward \code{inc_ncp}.

 
	} %\itemize
		
  }%\item{trend}
 
 

 \item{season}{a list object numeric consisting of various outputs related to the estimated seasonal/periodic component:
  	\itemize{
		 \item \code{ncp}: [Number of ChangePoints]. a numeric scalar; the mean number of seasonal changepoints.
		 \item \code{ncpPr}: [Probability of the Number of ChangePoints]. A vector of length \code{(prior$seasonMaxKnotNum+1)}. It gives a probability distribution of having a certain number of seasonal changepoints over the range of [0,prior$seasonMaxKnotNum]; for example, \code{ncpPr[1]} is the probability of having no seasonal changepoint; \code{ncpPr[i]} is the probability of having (i-1) changepoints: Note that  the index is i rather than (i-1) because ncpPr[1] is used for having zero changepoint.
		 \item \code{cpOccPr}: [ChangePoint OCCurence PRobability]. a vector of length N; it gives a probability distribution of having a changepoint in the seasonal component at each point of time. Plotting \code{cpOccPr} will depict a continious curve of probability-of-being-changepoint over the time. Of particular note, in the curve, a higher value at a peak indicates a higher chance of being a changepoint only at that particular SINGLE point in time, and does not necessarily mean a higher chance of observing a changepoint AROUND that time. For example, a window of cpOccPr values \code{c(0,0,0.5,0,0)} (i.e., the peak prob is 0.5 and the summed prob is 0.5) is less likely to be a changepoint compared to another window values \code{c(0.1,0.2,0.3,0.2,0.1)} (i.e., the peak prob is 0.3 but the summed prob is 0.8).
		 \item \code{order}: a vector of length N; the average harmonic order needed to approximate the seasonal component. As an average over many sampled individual piece-wise harmonic curves, \code{order} is not necessarily an integer.
		 \item \code{cp}: [Changepoints] a vector of length \code{metadata$seasonMaxKnotNum}; the most possible  changepoint locations in the seasonal component. The locations are obtained by first applying a sum-filtering to the \code{cpOccPr} curve with a filter window size of \code{prior$trendMinSeptDist} and then picking up to a total \code{ncp} of the highest peaks in the filtered curve. If \code{ncp<seasonMaxKnotNum}, the remaining of the vector is filled with NaNs.
		 \item \code{cpPr}: [Changepoints PRobability] a vector of length \code{metadata$seasonMaxKnotNum}; the probabilities associated with the changepoints \code{cp}. Filled with NaNs for the remaining elements if \code{ncp<seasonMaxKnotNum}.
		 \item \code{cpCI}: [Changepoints Credible Interval] a matrix of dimension \code{metadata$seasonMaxKnotNum x 2}; the credibable intervals for the detected changepoints \code{cp}.
		 \item \code{cpAbruptChange}: [Abrupt change at Changepoints]  a vector of length \code{metadata$seasonMaxKnotNum}; the jumps in the fitted seasonal curves at the detected changepoints \code{cp}.
		 \item \code{Y}: a vector of length N; the estimated seasonal component. It is the Bayesian model averaging of all the individual sampled seasonal curves.
		 \item \code{SD}: [Standard Deviation] a vector of length N; the estimated standard deviation of the estimated seasonal component.
		 \item \code{CI}: [Standard Deviation] a matrix of dimension \code{N x 2}; the estimated credible interval of the estimated seasonal component. One vector of the matrix is for the upper envelope and another for the lower envelope.
		 \item \code{amp}: [AMPlitude] a vector of length N;  the time-varying amplitude of the estimated seasonality.
		 \item \code{ampSD}: [Standar Deviation of AMPlitude] a vector of length N; , the SD of the amplitude of the seasonality.
 

		 \item \code{pos_ncp}: 
		 \item \code{neg_ncp}:  
		 \item \code{pos_ncpPr}:
		 \item \code{neg_ncpPr}:
		 \item \code{pos_cpOccPr}:
		 \item \code{neg_cpOccPr}:
		 \item \code{pos_cp}:
		 \item \code{neg_cp}:
		 \item \code{pos_cpPr}:
		 \item \code{neg_cpPr}:
		 \item \code{pos_cpAbruptChange}:
		 \item \code{neg_cpAbruptChange}:
		 \item \code{pos_cpCI}:
		 \item \code{neg_cpCI}:  The above variables have the same outputs as those variables without the prefix 'pos' and 'neg', except that we differentiate the changepoints with a POStive jump in the trend from those changepoints with a NEGative jump. For example, \code{pos_ncp} refers to the average number of trend changepoints that jump up (i.e., positively) in the trend.
 

 
	} %\itemize
		
  }%\item{trend}  
 

} %value
 
\seealso{
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{tetris}}, \code{\link[Rbeast]{geeLandsat}} 
}

\examples{

#--------------------------------NOTE-------------------------------------------------#
# beast123() is an all-inclusive function that duplicates the functionalities of beast 
# and beast.irreg. It can handle a single, multiple, or 3D of stacked time series, being
# either regular or irregular. It allows for customization through four LIST arguments:
#   metadata -- additional info about the input Y
#   prior    -- prior parameters for the beast model 
#   mcmc     -- MCMC simulation setting
#   extra    -- misc parameters turning on/off outputs and setting up parallel computations
#
# Despite being essentially the same as beast and beast.irreg, beast123 is provided mainly
# to support concurrent handling of multiple time series (e.g., stacked satellite images)
# via parallel computing: When processing stacked raster layers, DO NOT iterate pixel by pixel
# using beast() or beast.irreg() via an external parallel caller (e.g., doParallel or foreach).
# Instread, please use beast123(), which supports mulithreading internally. 

 
#------------------------------Example 1: one time series with seasonalty-------------#
# Yellowstone is a half-monthly time series of  774 NDVI measurmments at a Yellowstone
# site starting from July 1-15,1981(i.e., start=c(1981,7,7). It has 24 data points per 
# year (freq=24).        

 library(Rbeast)
 data(Yellowstone)
 plot(Yellowstone)

# Below, the four option args are missing, so defalut values will be used, with some 
# warning messages given to altert this. By default, the input Y is assumed to be regular
# with a seasonal component. The default arg values used will be printed out and they can
# serve as a template to customize the parameters.
\donttest{ 
 o = beast123(Yellowstone)
 plot(o)
}

#------------------------------Example 2: a trend-only time series-------------------#
# Nile is an annual river flow time series (i.e., no periodic variation). So, season 
# is set to 'none' to indicate trend-only analysis. Default values are used for other 
# missing options. Unlike the beast() function, beast123 does NOT use the time attributes 
# of a 'ts' object. For example, Nile is treated as a pure data number; its (start=1871,
# end=1970, freq=1) attributes are ignored. The default times 1:length(Nile) are used 
# instead. The true time info need to be specified by the 'metadata' parameter, as shown
# in the next example.
\donttest{ 
 o = beast123(Nile,season='none')
 plot(o)
}

#------------------------------Example 3: call via the full API interface-----------#
# Specify metadata, prior, mcmc, and extra explicitly. Only 'prior' is the true statistical
# model parameters of BEAST; the other three are just options to configure the input/ouput
# or the computation process.

\dontrun{ 

 # metadata is NOT part of BEAST itself, but some extra info to describe the input
 # time series Y. Below, the input Y is the 'Yellowstone' ts. 
 
 metadata                  = list()                                                 
 #metadata$isRegularOrdered = TRUE        # This arg not used any longer in this version
 metadata$whichDimIsTime   = 1            # Which dim of the input refer to time for 
                                          # 2D/3D inputs? Ignored for a single time 
                                          # series input.
 metadata$startTime        = c(1981,7,7)  # Or startTime=1981.5137
                                          #    startTime=as.Date('1981-7-7')   
 metadata$deltaTime        = 1/24         # Half-monthly regular ts: 0.5/12=1/24      
 metadata$period           = 1.0          # The period is 1 year: 
                                          # freq x deltaTime = period
                                          # 24    x  1/24    = 1.0
 metadata$omissionValue    = NaN          # By default, NaNs are ignored
 metadata$maxMissingRateAllowed = 0.7500  # If missingness is higher than .75, the ts
                                          #  is  skipped and not fitted
 metadata$deseasonalize    = FALSE        # Do not remove the global seasonal pattern 
                                          #  before fitting the beast model
 metadata$detrend          = FALSE        # Do not remove the global trend before
                                          # the fitting
 
 # prior is the ONLY true parameters of the beast model,used to specify the priors 
 # in the Bayesian formulation  
 prior = list()
 prior$seasonMinOrder   = 1         #min harmonic order allowed to fit seasonal cmpnt
 prior$seasonMaxOrder   = 5         #max harmonic order allowed to fit seasonal cmpnt
 prior$seasonMinKnotNum = 0         #min number of changepnts in seasonal cmpnt
 prior$seasonMaxKnotNum = 3         #max number of changepnts in seasonal cmpnt
 prior$seasonMinSepDist = 10        #min inter-chngpts separation for seasonal cmpnt
 prior$trendMinOrder	= 0         #min polynomial order allowed to fit trend cmpnt 
 prior$trendMaxOrder	= 1         #max polynomial order allowed to fit trend cmpnt   
 prior$trendMinKnotNum  = 0         #min number of changepnts in trend cmpnt
 prior$trendMaxKnotNum  = 15        #max number of changepnts in trend cmpnt
 prior$trendMinSepDist  = 5         #min inter-chngpts separation for trend cmpnt
 prior$precValue        = 10.0      #Initial value of the precision parameter (no 
                                    # need to change it unless for precPrioType='const')
 prior$precPriorType    = 'uniform' # Possible values: const, uniform, and componentwise

 # mcmc is NOT part of the beast model itself, but some parameters to configure the 
 # MCMC inference. 
 mcmc = list() 
 mcmc$seed                      = 9543434# an arbitray seed for random number generator
 mcmc$samples                   = 3000   # samples collected per chain
 mcmc$thinningFactor            = 3      # take every 3rd sample and discard others
 mcmc$burnin                    = 150    # discard the initial 150 samples per chain
 mcmc$chainNumber               = 3      # number of chains 
 mcmc$maxMoveStepSize           = 4      # max random jump step when proposing new chngpts
 mcmc$trendResamplingOrderProb  = 0.100  # prob of choosing to resample polynomial order
 mcmc$seasonResamplingOrderProb = 0.100  # prob of choosing to resample harmonic order
 mcmc$credIntervalAlphaLevel    = 0.950  # the significance level for credible interval


 # extra is NOT part of the beast model itself, but some parameters to configure the 
 # output and computation process
 extra = list()
 extra$dumpInputData        = FALSE #If true, a copy of input time series is outputted 
 extra$whichOutputDimIsTime = 1     #For 2D or 3D inputs, which dim of the output refers to 
                                    # time? Ignored if the input is a single time series
 extra$computeCredible      = FALSE #If true, compute CI: computing CI is time-intensive. 
 extra$fastCIComputation    = TRUE  #If true, a faster way is used to get CI, but it is 
                                    # still time-intensive. That is why the function beast()
                                    # is slow because it always compute CI.
 extra$computeSeasonOrder   = FALSE #If true, dump the estimated harmonic order over time
 extra$computeTrendOrder    = FALSE #If true, dump the estimated polynomial order over time
 extra$computeSeasonChngpt  = TRUE  #If true, get the most likely locations of s chgnpts
 extra$computeTrendChngpt   = TRUE  #If true, get the most likely locations of t chgnpts
 extra$computeSeasonAmp     = FALSE #If true, get time-varying amplitude of seasonality
 extra$computeTrendSlope    = FALSE #If true, get time-varying slope of trend
 extra$tallyPosNegSeasonJump= FALSE #If true, get those changpts with +/- jumps in season
 extra$tallyPosNegTrendJump = FALSE #If true, get those changpts with +/- jumps in trend 
 extra$tallyIncDecTrendJump = FALSE #If true, get those changpts with increasing/
                                    # decreasing trend slopes 
 extra$printProgressBar     = TRUE
 extra$printOptions         = TRUE
 extra$quiet                = FALSE # print warning messages, if any
 extra$consoleWidth         = 0     # If 0, the console width is from the current console
 extra$numThreadsPerCPU     = 2     # 'numThreadsPerCPU' and 'numParThreads' are used to 
 extra$numParThreads        = 0     # configure multithreading runs; they're used only if
                                    # Y has multiple time series (e.g.,stacked images)

 o = beast123(Yellowstone,metadata,prior,mcmc,extra, season='harmonic')
 plot(o)
}
  


#------------------------------Example 4: Handle irregular time series-----------------#
#  Handle irregular time series: ohio is a data frame of a Landsat NDVI series observed
#  at unevely-spaced times
\donttest{ 
\dontrun{ 

 data(ohio)
 str(ohio)
  
 metadata                  = list()                                                 
 metadata$time             = ohio$time # Must supply individual times for irregular inputs
 metadata$deltaTime        = 1/12      # Must supply the desired time interval for aggregation
 metadata$period           = 1.0
 
 o=beast123(ohio$ndvi, metadata)       # Default values used for those missing parameters

 ##################################################################################### 
 class(ohio$rdate)                      # Another accepted time format for beast123
 
 metadata = list()                                                 
 metadata$deltaTime        = 1/12       # Must supply the desired time interval for aggregation
 metadata$time             = ohio$rdate # Must supply individual times for irregular inputs
 
 o=beast123(ohio$ndvi, metadata)        # Default values used for those missing parameters
 
 ##################################################################################### 
 ohio$Y                              # Another accepted time format for beast123
 ohio$M
 ohio$M
 
 metadata = list()                                                 
 metadata$deltaTime        = 1/12    # Must supply the desired time interval for aggregation
 metadata$time$year        = ohio$Y
 metadata$time$month       = ohio$M
 metadata$time$day         = ohio$D 
 o=beast123(ohio$ndvi, metadata)     # Default values used for those missing parameters
 
 ##################################################################################### 
 ohio$Y                              # Another accepted time format for beast123
 ohio$doy
  
 metadata = list()                                                 
 metadata$deltaTime        = 1/12    # Must supply the desired time interval for aggregation
 metadata$time$year        = ohio$Y
 metadata$time$doy         = ohio$doy
 o=beast123(ohio$ndvi, metadata)     # Default values used for those missing parameters
 
 ##################################################################################### 
 ohio$time                           # Another accepted time format for beast123
  
 metadata = list()                                                 
 metadata$deltaTime        = 1/12      # Must supply the desired time interval for aggregation
 metadata$time             = ohio$time # Fractional year
 
 o=beast123(ohio$ndvi, metadata)       # Default values used for those missing parameters 
 
 ##################################################################################### 
 ohio$datestr1                       # Another accepted time format for beast123
  
 metadata = list()                                                 
 metadata$deltaTime        = 1/12         # Must supply the time interval for aggregation
 metadata$time$datestr     = ohio$datestr1  
 metadata$time$strfmt      = '????yyyy?mm?dd'  
 
 o=beast123(ohio$ndvi, metadata)          # Default values used for those missing parameters
 
 
 ##################################################################################### 
 ohio$datestr2                            # Another accepted time format for beast123  
 metadata = list()                                                 
 metadata$deltaTime        = 1/12         # Must supply a desired time interval for aggregation
 metadata$time$datestr     = ohio$datestr2  
 metadata$time$strfmt      = '????yyyydoy????'  
 
 o=beast123(ohio$ndvi, metadata)          # Default values used for those missing parameters  
  
 ##################################################################################### 
 ohio$datestr3                            # Another accepted time format for beast123  
 metadata = list()                                                 
 metadata$deltaTime        = 1/12         # Must supply the desired time interval for aggregation
 metadata$time$datestr     = ohio$datestr3  
 metadata$time$strfmt      = 'Y,,M/D'  
 
 o=beast123(ohio$ndvi, metadata)          # Default values used for those missing parameters 
}
}


#------------------Example 4: Handle multiple time series (i.e., matrix input)-----------#
# Handle multiple time series: 'simdata' is a 2D matrix of dim 300x3; it consits of 3
# time series of length 300 each. For this toy example, I decide to be lazy and use the same 
# time series for the three columns. 
\dontrun{ 
 data(simdata)                         # dim of simdata: 300 x 3 (time x num_of_time_series)
 dim(simdata)                          # the first dimenion refer to time (i.e, 300)
 
 metadata                  = list() 
 metadata$whichDimIsTime   = 1         # Which dim of the input refer to time for 2D inputs?
                                       # 300 is the ts length, so dim is set to '1' here.
 metadata$period           = 24        # By default, we assume startTime=1 and deltaTime=1 
 
 extra=list()
 extra$whichOutputDimIsTime = 2            # Which dim of the output arrays refers to  time?
 o=beast123(simdata, metadata,extra=extra) # Default values used for those missing parameters
 
 # The lists of arg parameters can also be directly provided inline within the command
 o=beast123( simdata, metadata=list(whichDimIsTime=1,period=24), extra=list(whichOutput=2) )
 
 # The field names of the lists can be shortened as long as no ambiguitity is caused.
 o=beast123( simdata, metadata=list(whichDim=1,per=24), extra=list(whichOut=2) )
 
 #------------------Example 4: Another run by transposing simdata--------------------------#
  
 simdata1=t(simdata)                   # dim of simdata1: 3 x 300 (num of ts  x time )
 
 metadata                  = list() 
 metadata$whichDimIsTime   = 2         # Which dim of the input refer to time for 2D inputs?
                                       # 300 is the ts length, so dim is set to '2' here.
 metadata$period           = 24        # By default, we assume startTime=1 and deltaTime=1 
 o=beast123(simdata1, metadata)        # Default values used for those missing parameters 
 
 o=beast123( simdata1, metadata=list(whichDim=2, per=24) )
 }



#------------------Example 5: Handle stacked time series images (e.g., 3d input)--------#
# Handle 3D stacked images of irregular and unordered time-series: imagestack is a 3D 
# array of size 12x9x1066, each pixel being a time series of length 1066
\dontrun{ 
 data(imagestack)
 dim(imagestack$ndvi)               # Dim: 12 x 9 X 1066 (row x col x time)
 imagestack$datestr                 # A character vector of 1066 date strings
 
 metadata                  = list()        
 metadata$whichDimIsTime   = 3      # Which dim of the input refer to time for 3D inputs?
                                    # 1066 is the ts length, so dim is set to '3' here.
                                    # In this example, this arg is not needed because 
                                    # the time$datestr can also help to match and pick up
                                    # the right time dimesion of imagestack$ndvi.
 metadata$time$datestr     =  imagestack$datestr 
 metadata$time$strfmt      =  'LT05_018032_20080311.yyyy-mm-dd'
 metadata$deltaTime        =  1/12  # Aggregate the irregular ts at a monthly interval:1/12 Yr
 metadata$period           =  1.0   # The period is 1 year: deltaTime*freq=1/12*12=1.0

 extra = list()
 extra$dumpInputData        = TRUE  # Get a copy of aggregated input ts 
 extra$numThreadsPerCPU     = 2     # Each cpu core will be assigned 2 threads
 extra$numParThreads        = 0     # If 0, total_num_threads=numThreadsPerCPU*num_of_cpu_core
                                    # if >0, used to specify the total number of threads
  
 # Default values for missing parameters 
 o=beast123(imagestack$ndvi, metadata=metadata,extra=extra) 
 
 print(o,c(5,3))    # print the result for the pixel at Row 5 and Col 3
 plot(o,c(5,3))     # plot the result for the pixel at Row 5 and Col 3
 image(o$trend$ncp) # number of trend changepoints over space
} 
 
#---------Example 6: Handle stacked GeoTiff image files imported with the raster package------#
# Handle 3D stacked images of irregular time-series : 'ndvi.zip' is  a zip file of 
# 437 NDIV tiff image files, each having a dim of 12 x 9.
# Code availlable at https://github.com/zhaokg/Rbeast/blob/master/R/beast123_raster_example.txt
 
}

\references{
\enumerate{
\item Zhao, K., Wulder, M.A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick, B., Zhang, X. and Brown, M., 2019. Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment, 232, p.111181 (the beast algorithm paper).
\item Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection. Remote Sensing of Environment, 132, pp.102-119 (the Bayesian MCMC scheme used in beast). 
\item Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine. ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261(a beast application paper). 
}
}
\keyword{misc}
