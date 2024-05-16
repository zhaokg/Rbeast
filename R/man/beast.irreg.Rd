\name{beast.irreg}
\docType{methods}
\alias{beast.irreg}
\alias{BEAST.irreg} 

\title{Bayesian time series decomposition for changepoint, trend, and periodicity or seasonality}

\description{
   A Bayesian model averaging algorithm called BEAST to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations. BEAST is useful for changepoint detection (e.g., breakpoints or structural breaks), nonlinear trend analysis, time series decomposition, and time series segmentation.
}

\usage{

  beast.irreg(\
       y, 
       time,        
       deltat         = NULL,
       period         = NULL,  		   
       season         = c("harmonic", "svd", "dummy", "none"),  
       scp.minmax     = c(0,10), sorder.minmax = c(0,5),  	  
       tcp.minmax     = c(0,10), torder.minmax = c(0,1), 	   
       sseg.min       = NULL,    sseg.leftmargin = NULL, sseg.rightmargin = NULL,
       tseg.min       = NULL,    tseg.leftmargin = NULL, tseg.rightmargin = NULL, 
       method         = c( 'bayes', 'bic', 'aic', 'aicc','hic'),
       detrend        = FALSE, 
       deseasonalize  = FALSE,
       mcmc.seed      = 0,      
       mcmc.burnin    = 200, 
       mcmc.chains    = 3,
       mcmc.thin      = 5,
       mcmc.samples   = 8000,
       ci             = FALSE,	
       precValue      = 1.5,
       precPriorType  = c('componentwise','uniform','constant','orderwise'),
       print.options  = TRUE,
       print.progress = TRUE,
       quiet          = FALSE,
       gui            = FALSE,
       ...
      )
}
\note{
x
} 

\arguments{

  \item{y}{ a vector for an irregular or unordered time series. Missing values such as NA and NaN are allowed.
			\itemize{
				\item  If \code{y} is regular and evenly-spaced in time, use the \code{\link[Rbeast]{beast}}function instead.
				\item  If \code{y} is a matrix or 3D array (e.g., stacked images) consisting of multiple regular or irregular time series, use \code{\link[Rbeast]{beast123}} instead. 
            }  
			If \code{y} is a list of multiple time series, the multivariate version of the BEAST algorithm is invoked to decompose the multiple time series and detect common changepoints altogether. This feature is \code{experimental} and under further development. Check \code{\link[Rbeast]{ohio}} for a working example.
  }

  \item{time}{ a vector of the same length as \code{y}'s time dimension to provide the times for datapoints. It can be a vector of numbers, Dates, or date strings; it can also be a list of vectors of year, months, and days. Possible formats include:
      		 
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
  
  }
  
  \item{deltat}{ a number or a string to specify a time interval for aggregating the irregular \code{y} into a regular time series. The BEAST model is currently formulated for regular data only for fast computational, so internally, the \code{beast.irreg} function will aggregate/re-bin irregular data into regular ones. For the aggregation, \code{deltat} is needed to specify the desired bin size or time interval; if missing, a best guess will be used. The unit of \code{deltat} needs to be consistent with \code{time}. If \code{time} takes a numeric vector, the unit of \code{deltat} is arbitrary and irrelevant to beast. If \code{time} takes a vector of Dates or date strings, the unit for \code{deltat} is assumed to Fractional YEAR. If needed, use a string instead of a number to specify whether the unit of deltat is day, month, or year. Examples include '7 days', '7d', '1/2 months', '1mn', '1.0 year', and '1y'.
  
  }  

  \item{period}{ numeric or string. Specify the period for the periodic/seasonal component in \code{y}. Needed only for data with a periodic/cyclic component (i.e., season=\code{'harmonic'} or \code{'dummy'}) and not used for trend-only data (i.e., \code{season='none'}).  The period of the cyclic component should have a unit consisent with the unit of \code{deltat}. It holds that \code{period=deltat*freq} where \code{freq} is the number of data samples per period. (Note that the \code{freq} argument in earlier versions becomes obsolete and now is replaced by \code{period}. \code{freq} is still supported but\code{period} takes precedence if both are provided.)  \code{period} or the number of data points per period is not a BEAST model parameter and it has to be specified by the user. But if \code{period} is missing, BEAST first attempts to guess its value via auto-correlation before fitting the model. If \code{period} <= 0, \code{season='none'} is assumed, and  the trend-only model is fitted without a seasonal/cyclic component. If needed, use a string to specify whether the unit of period is day, month, or year. Examples are '1.0 year', '12 months', '365d', '366 days'.
  }
  
  \item{season}{ characters (default to 'harmonic'); specify if \code{y} has a periodic component or not. Three strings are possible.
				\itemize{
					\item  \code{'none'}: \code{y} is trend-only; no periodic components are present in the time series. The args for the seasonal component (i.e.,sorder.minmax, scp.minmax and sseg.max) will be ignored. 					
					\item  \code{'harmonic'}: \code{y} has a periodic/seasonal component. The term 'season' is a misnomer, being used here to broadly refer to any periodic variations present in \code{y}. The periodicity is NOT a model parameter estimated by beast but a known constant given by the user through \code{freq}. By default, the periodic component is modeled as a harmonic curve--a combination of sins and cosines.
					\item  \code{'dummy'}: the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a non-parametric curve. The arg \code{sorder.minmax} is irrelevant and is ignored.
					\item  \code{'svd'}: (experimental feature) the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a linear combination of function bases derived from a Single-value decomposition. The SVD-based basis functions are more parsimonious than the harmonic sin/cos bases in parameterizing the seasonal variations; therefore, more subtle changepoints are likely to be detected.
				}  
  
  }
  
   \item{scp.minmax}{ a vector of 2 integers (>=0); the min and max number of seasonal changepoints (scp) allowed in segmenting the seasonal component. \code{scp.minmax} is used only if \code{y} has a seasonal component (i.e., \code{season='harmonic' or 'dummy' }) and ignored for trend-only data. If the min and max changepoint numbers are equal, BEAST assumes a constant number of scp and won't infer the posterior probability of the number of changepoints, but it still estimates the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur). If both the min and max numbers are set to 0, no changepoints are allowed; then a global harmonic model is used to fit the seasonal component, but still, the most likely harmonic order will be inferred if sorder.minmax[1] is not equal to sorder.minmax[2].}
	
  \item{sorder.minmax}{ a vector of 2 integers (>=1); the min and max harmonic orders considered to fit the seasonal component. \code{sorder.minmax} is used only used if the time series has a seasonal component (i.e., \code{season='harmonic'}) and ignored for trend-only data or when \code{season='dummy'}. If the min and max orders are equal (\code{sorder.minmax[1]=sorder.minmax[2]}), BEAST assumes a constant harmonic order used and won't infer the posterior probability of harmonic orders.
  }
  

    
 
	   
  %%%%%%%%%%%%%%%%%%%%%	
  
  \item{torder.minmax}{ a vector of 2 integers (>=0); the min and max orders of the polynomials considered to fit the trend component. The 0-th order corresponds to a constant term/a flat line and the 1st order is a line. If \code{torder.minmax[1]=torder.minmax[2]}, BEAST assumes a constant polynomial order used and won't infer the posterior probability of polynomial orders.  }
  
  \item{tcp.minmax}{ a vector of 2 integers; the min and max number of trend changepoints (tcp) allowed in segmenting the trend component. If the min and max changepoint numbers are equal, BEAST assumes a constant number of changepoints and won't infer the posterior probability of the number of changepoints for the trend, but it still estimates the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur in the trend). If both the min and max numbers are set to 0, no changepoints are allowed; then a global polynomial trend is used to fit the trend component, but still, the most likely polynomial order will be inferred if torder.minmax[1] is not equal to torder.minmax[2].}
    
 
  %%%%%%%%%%%%%%%%%%%%%	
  
  \item{sseg.min}{ an integer (>0); the min segment length allowed between two neighboring season changepoints. That is, when fitting a piecewise harmonic seasonal model, two changepoints are not allowed to occur within a time window of length \code{sseg.min}. \code{sseg.min} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{sseg.min*deltat}. \code{sseg.min} defaults to NULL and its value will be given a default value in reference to \code{freq}. }	
  
  \item{sseg.leftmargin}{ an integer (>=0); the number of leftmost data points excluded for seasonal changepoint detection. That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the starting window/segment of length \code{tseg.leftmargin}. \code{sseg.leftmargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{sseg.leftmargin*deltat}. If missing, \code{sseg.leftmargin} defaults to \code{sseg.min}.}
  
  \item{sseg.rightmargin}{ an integer (>=0); the number of rightmost data points excluded for seasonal changepoint detection. That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the ending window/segment of length \code{sseg.rightmargin}. \code{sseg.rightmargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{sseg.rightmargin*deltat}. If missing, \code{sseg.rightmargin} defaults to \code{sseg.min}.}  
    	
  \item{tseg.min}{ an integer (>0); the min segment length allowed between two neighboring trend changepoints. That is, when fitting a piecewise polynomial trend model, two changepoints are not allowed to occur within a time window of length \code{tseg.min}. \code{tseg.min} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{tseg.min*deltat}. \code{tseg.min} defaults to NULL and its value will be given a default value in reference to \code{freq} if the time series has a cyclic component.}
  
  \item{tseg.leftmargin}{ an integer (>=0); the number of leftmost data points excluded for trend changepoint detection. That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the starting window/segment of length \code{tseg.leftmargin}. \code{tseg.leftmargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{tseg.leftmargin*deltat}. If missing, \code{tseg.leftmargin} defaults to \code{tseg.min}.}
  
  \item{tseg.rightmargin}{ an integer (>=0); the number of rightmost data points excluded for trend changepoint detection. That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the ending window/segment of length \code{tseg.rightmargin}. \code{tseg.rightmargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{tseg.rightmargin*deltat}. If missing, \code{tseg.rightmargin} defaults to \code{tseg.min}.}  
  
  %%%%%%%%%%%%%%%%%%%%%	
  
  \item{method}{ a string (default to 'bayes'); specify the method for formulating model posterior probability.
         \itemize{
				  \item  \code{'bayes'}: the full Bayesian formulation as described in Zhao et al. (2019). 					
				  \item  \code{'bic'}:  approximation of posterior probability using the Bayesian information criterion (bic).					
				  \item  \code{'aic'}:  approximation of posterior probability using the Akaike information criterion (aic).					
 				  \item  \code{'aicc'}:  approximation of posterior probability using the corrected Akaike information criterion (aicc).					
 				  \item  \code{'hic'}:  approximation of  posterior probability using the Hannan-Quinn information criterion  (hic)					
				}    		   
    }
 
   
  \item{detrend}{ logical; If \code{TRUE}, a global trend is first fitted and removed from the time series before running BEAST; after BEAST finishes, the global trend is added back to the BEAST result. }
  
  \item{deseasonalize}{ logical; If \code{TRUE}, a global seasonal model is first fitted and removed from the time series before running BEAST; after BEAST finishes, the global seasonal curve is added back to the BEAST result. \code{deseasonalize} is ignored if \code{season='none'} (i.e., trend-only data). }	
  
  \item{mcmc.seed}{ integer (>=0); the seed for the random number generator used for Monte Carlo Markov Chain (mcmc).  If \code{mcmc.seed=0}, an arbitrary seed is picked and the fitting results vary across runs. If fixed to the same non-zero integer, the result can be re-produced for different runs. But the results from the same seed may still vary if run on different computers because the random generator library depends on CPU's instruction sets.}	 
 
  \item{mcmc.chains}{ integer (>0); the number of  MCMC chains.}
  
  \item{mcmc.thin}{ integer (>0); a factor to thin chains (e.g., if thinningFactor=5, samples will be taken every 3 iterations)}
  
  \item{mcmc.burnin}{ integer (>0); the number of burn-in samples discarded at the start of each chain}
  
  \item{mcmc.samples}{ integer (>=0); the number of samples collected per MCMC chain. The total number of iterations is \code{(burnin+samples*thin)*chains}.}		  
  
  \item{ci}{ boolean; If \code{TRUE}, credible intervals (i.e., out$season$CI or out$trend$CI) will be computed for the estimated seasonal and trend components. Computing CI is time-consuming, due to sorting, so set \code{ci} to FALSE if a symmetric credible interval (i.e., out$trend$SD and out$season$SD) suffices.}	

  \item{precValue}{  numeric (>0); the hyperparameter of the precision prior; the default value is 1.5. \code{precValue} is useful only when \code{precPriorType}='constant', as further explained below}
  \item{precPriorType}{ characters. It takes one of 'constant', 'uniform',  'componentwise' (the default), and 'orderwise'. Below are the differences between them.
            \enumerate{
                 \item \code{'constant'}:  the precision parameter used to parameterize the model coefficients is fixed to a constant specified by \code{precValue}. In other words, \code{precValue} is a user-defined hyperparameter and the fitting result may be sensitive to the chosen values of \code{precValue}.
                 \item \code{'uniform'}:  the precision parameter used to parameterize the model coefficients is a random variable; its initial value is specified by \code{precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result will be insensitive to the choice in \code{precValue}.
                 \item \code{'componentwise'}: multiple precision parameters are used to parameterize the model coefficients for individual components (e.g., one for season and another for trend); their initial values is specified by \code{precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result will be insensitive to the choice in \code{precValue}.
                 \item \code{'orderwise'}: multiple precision parameters are used to parameterize the model coefficients not just for individual components but also for individual orders of each component; their initial values is specified by \code{precValue}. In other words, \code{precValue} will be inferred by the MCMC, so the fitting result will be insensitive to the choice in \code{precValue}.
          } %componentwise
  }

  
  \item{print.options}{ boolean. If \code{TRUE},the full list of input parameters to BEAST  will be printed out prior to the MCMC inference; the naming for this list (e.g., metadata, prior, and mcmc) differs slightly from the input to \code{beast}, but there is a one-to-one correspondence (e.g., prior$trendMinSepDist=tseg.min). Internally, beast converts the input parameters to the forms of metadata, prior,and mcmc. Type 'View(beast)' to see the details or check the \code{beast123} function. }	
  
  \item{print.progress}{ boolean;If \code{TRUE}, a progressbar will be displayed. }	
  
  \item{quiet}{ boolean. If \code{TRUE}, warning messages are suppressed and not printed.}	
 	
  
  \item{gui}{ boolean. If \code{TRUE}, BEAST will be run in a GUI demonstration mode, with a GUI window to show an animation of the MCMC sampling in the model space step by step. Note that "\code{gui=TRUE}" works only for Windows x64 systems not Windows 32 or Linux/Mac systems.}	
 	
  \item{\dots}{additional parameters. There are many more settings for the implementation but not made available in the beast() interface; please use the function \code{beast123()} instead }
 
}

\value{
  The output is an object of class "beast". It is a list, consisting of the following variables. In the explanations below, we assume the input \code{y} is a single time series of length \code{N}:

  \item{time}{a vector of size \code{1xN}: the times at the \code{N} sampled locations. By default, it is simply set to \code{1:N} if the input arguments \code{delta}, \code{start}, or \code{time} are missing.}
  \item{data}{a vector, matrix, or 3D array; this is a copy of the input \code{data} if extra$dumpInputData = TRUE. If extra$dumpInputData=FALSE, it is set to NULL. If the original input \code{data} is irregular, the copy here is the regular version aggregated from the original at the time interval specified by metadata$deltaTime. }
  \item{marg_lik}{numeric; the average of the model marginal likelihood; the larger marg_lik, the better the fitting for a given time series.}
  \item{R2}{numeric; the R-square of the model fitting.}
  \item{RMSE}{numeric; the RMSE of the model fitting.}
  \item{sig2}{numeric; the estimated variance of the model error.}
  \item{trend}{a list object consisting of various outputs related to the estimated trend component:
  	\itemize{
		 \item \code{ncp}: [Number of ChangePoints]. a numeric scalar; the mean number of trend changepoints. Individual models sampled by BEAST has a varying dimension (e.g., number of changepoints or knots), so several alternative statistics (e.g., \code{ncp_mode}, \code{ncp_median}, and \code{ncp_pct90}) are also given to summarize the number of changepoints. For example, if \code{mcmc$samples=10}, the numbers of changepoints for the 10 sampled models are assumed to be c(0, 2, 4, 1, 1, 2, 7, 6, 6, 1). The mean ncp is 3.1 (rounded to 3), the median is 2.5 (2),  the mode is 1, and the 90th percentile (ncp_pct90) is 6.5. 
		 \item \code{ncp_mode}: [Number of ChangePoints]. a numeric scalar; the mode for number of changepoints. See the above for explanations.
		 \item \code{ncp_median}: [Number of ChangePoints]. a numeric scalar; the median for number of changepoints. See the above for explanations.
		 \item \code{ncp_pct90}: [Number of ChangePoints]. a numeric scalar; the 90th percentile for number of changepoints. See the above for explanations.
		 \item \code{ncpPr}: [Probability of the Number of ChangePoints]. A vector of length \code{(tcp.minmax[2]+1)=tcp.max+1}. It gives a probability distribution of having a certain number of trend changepoints over the range of [0,tcp.max]; for example, \code{ncpPr[1]} is the probability of having no trend changepoint; \code{ncpPr[i]} is the probability of having (i-1) changepoints: Note that it is \code{ncpPr[i]} not \code{ncpPr[i-1]} because ncpPr[1] is used for having zero changepoint.
		 \item \code{cpOccPr}: [ChangePoint OCCurence PRobability]. a vector of length N; it gives a probability distribution of having a changepoint in the trend at each point of time. Plotting \code{cpOccPr} will depict a continious curve of probability-of-being-changepoint. Of particular note, in the curve, a higher peak indicates a higher chance of being a changepoint only at that particular SINGLE point in time and does not necessarily mean a higher chance of observing a changepoint AROUND that time. For example, a window of cpOccPr values \code{c(0,0,0.5,0,0)} (i.e., the peak prob is 0.5 and the summed prob is 0.5) is less likely to be a changepoint compared to another window \code{c(0.1,0.2,0.21,0.2,0.1)} (i.e., the peak prob is 0.21 but the summed prob is 0.71).
		 \item \code{order}: a vector of length N; the average polynomial order needed to approximate the fitted trend. As an average over many sampled individual piece-wise polynomial trends, \code{order} is not necessarily an integer.
		 \item \code{cp}: [Changepoints] a vector of length \code{tcp.max=tcp.minmax[2]}; the most possible changepoint locations in the trend component. The locations are obtained by first applying a sum-filtering to the \code{cpOccPr} curve with a filter window size of \code{tseg.min} and then picking up to a total \code{prior$MaxKnotNum}/\code{tcp.max} of the highest peaks in the filtered curve.  NaNs are possible if no enough changepoints are identified. \code{cp} records all the possible changepoints identified and many of them are bound to be false positives. Do not blindly treat all of them as actual changepoints.
		 
		 \item \code{cpPr}: [Changepoints PRobability] a vector of length \code{tcp.max=tcp.minmax[2]}; the probabilities associated with the changepoints \code{cp}. Filled with NaNs for the remaining elements if \code{ncp<tcp.max}.
		 \item \code{cpCI}: [Changepoints Credible Interval] a matrix of dimension \code{tcp.max x 2}; the credibable intervals for the detected changepoints \code{cp}.
		 \item \code{cpAbruptChange}: [Abrupt change at Changepoints]  a vector of length \code{tcp.max}; the jumps in the fitted trend curves at the detected changepoints \code{cp}.
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
 

 \item{season}{a list object consisting of various outputs related to the estimated seasonal/periodic component:
  	\itemize{
		 \item \code{ncp}: [Number of ChangePoints]. a numeric scalar; the mean number of seasonal changepoints.
		 \item \code{ncpPr}: [Probability of the Number of ChangePoints]. A vector of length \code{(scp.minmax[2]+1)=scp.max+1}. It gives a probability distribution of having a certain number of seasonal changepoints over the range of [0,scp.max]; for example, \code{ncpPr[1]} is the probability of having no seasonal changepoint; \code{ncpPr[i]} is the probability of having (i-1) changepoints: Note that  the index is i rather than (i-1) because ncpPr[1] is used for having zero changepoint.
		 \item \code{cpOccPr}: [ChangePoint OCCurence PRobability]. a vector of length N; it gives a probability distribution of having a changepoint in the seasonal component at each point of time. Plotting \code{cpOccPr} will depict a continious curve of probability-of-being-changepoint over the time. Of particular note, in the curve, a higher value at a peak indicates a higher chance of being a changepoint only at that particular SINGLE point in time, and does not necessarily mean a higher chance of observing a changepoint AROUND that time. For example, a window of cpOccPr values \code{c(0,0,0.5,0,0)} (i.e., the peak prob is 0.5 and the summed prob is 0.5) is less likely to be a changepoint compared to another window values \code{c(0.1,0.2,0.3,0.2,0.1)} (i.e., the peak prob is 0.3 but the summed prob is 0.8).
		 \item \code{order}: a vector of length N; the average harmonic order needed to approximate the seasonal component. As an average over many sampled individual piece-wise harmonic curves, \code{order} is not necessarily an integer.
		 \item \code{cp}: [Changepoints] a vector of length \code{scp.max=scp.minmax[2]}; the most possible  changepoint locations in the seasonal component. The locations are obtained by first applying a sum-filtering to the \code{cpOccPr} curve with a filter window size of \code{sseg.min} and then picking up to a total \code{ncp} of the highest peaks in the filtered curve. If \code{ncp<scp.max}, the remaining of the vector is filled with NaNs.
		 \item \code{cpPr}: [Changepoints PRobability] a vector of length \code{scp.max}; the probabilities associated with the changepoints \code{cp}. Filled with NaNs for the remaining elements if \code{ncp<scp.max}.
		 \item \code{cpCI}: [Changepoints Credible Interval] a matrix of dimension \code{scp.max x 2}; the credibable intervals for the detected changepoints \code{cp}.
		 \item \code{cpAbruptChange}: [Abrupt change at Changepoints]  a vector of length \code{scp.max}; the jumps in the fitted seasonal curves at the detected changepoints \code{cp}.
		 \item \code{Y}: a vector of length N; the estimated seasonal component. It is the Bayesian model averaging of all the individual sampled seasonal curve.
		 \item \code{SD}: [Standard Deviation] a vector of length N; the estimated standard deviation of the estimated seasonal component.
		 \item \code{CI}: [Standard Deviation] a matrix of dimension \code{N x 2}; the estimated credible interval of the estimated seasonal curve. One vector of the matrix is for the upper envelope and another for the lower envelope.
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
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{tetris}}, \code{\link[Rbeast]{geeLandsat}} 
}

 
\examples{

library(Rbeast)

\donttest{

######################################################################################
# Note that the BEAST algorithm is currently implemented to handle only regular time 
# series. 'beast.irreg' accepts irregular time series but internally it aggregates them
# into regular ones prior to applying the BEAST model. For the aggregation, both the 
# "time" and "deltat" args are needed to specify individual times of data points and the
# regular time interval desired. If there is  a cyclic componet, 'period' should also be given; 
# if not, a possible value is guessed via auto-correlation


######################################################################################
# 'ohio' is a data.frame on an irregular Landsat time series of reflectances & ndvi 
# (e.g., surface greenness) at an Ohio site. It has multiple columns of alternative date 
# formats, such as year, month, day, doy (date of year), rdate (R's date class), and
# time (fractional year)

 data(ohio)
 str(ohio)
 plot(ohio$rdate, ohio$ndvi,type='o') # ndvi is irregularly spaced and unordered in time
 
######################################################################################
# Below, 'time' is given as numeric values, which can be of any arbitray unit. Although
# here 1/12 can be interepreted as 1/12 year or 1 month, BEAST itself doesn't care about
# the time unit. So, the unit of 1/12 is irrelevant for BEAST. 'freq' or 'period' is missing
# and a guess of it is used.

 o=beast.irreg(ohio$ndvi, time=ohio$time,deltat=1/12) 
 plot(o)
 print(o)

######################################################################################
# Aggregrate the time series at a monthly interval (deltat=1/12) and explictly provide
# the 'freq' or 'period' arg

 o=beast.irreg(ohio$ndvi, time=ohio$time,deltat=1/12, period=1.0) 
#o=beast.irreg(ohio$ndvi, time=ohio$time,deltat=1/12, freq  =12) 

}

\dontrun{
######################################################################################
# Aggregate the time series at a half-monthly time interval, and the 'freq' becomes 24 
# while the period is still 1. That is, PERIOD (1.0)=deltat(1/24) X  freq (24)

 o=beast.irreg(ohio$ndvi, time=ohio$time,deltat=1/24, freq   = 24) 
#o=beast.irreg(ohio$ndvi, time=ohio$time,deltat=1/24, period = 1) 

######################################################################################
# 'time' is given as R's dates. The unit is YEAR. 1/12 refers to 1/12 year or 1 month

 o=beast.irreg(ohio$ndvi, time=ohio$rdate,deltat=1/12) 

######################################################################################
# 'time' is given as data strings. The unit is YEAR. 1/12 refers to 1/12 year or 1 month


 o=beast.irreg(ohio$ndvi, time=ohio$datestr1,deltat=1/12)  #"LT4-1984-03-27"  (YYYY-MM-DD)
 o=beast.irreg(ohio$ndvi, time=ohio$datestr2,deltat=1/12)  #"LT4-1984087ndvi" (YYYYDOY)
 o=beast.irreg(ohio$ndvi, time=ohio$datestr3,deltat=1/12)  #"1984,, 3/ 27"    (YYYY M D)
 


######################################################################################
# 'time' is given as data strings, with a format specifier

 

 TIME =list()
 TIME$datestr = ohio$datestr1
 TIME$strfmt  = "LT4-YYYY-MM-DD"   # "LT4-1984-03-27"
 o=beast.irreg(ohio$ndvi, time=TIME,deltat=1/12)  
 
 TIME =list()
 TIME$datestr = ohio$datestr2
 TIME$strfmt  = "LT4-YYYYDOYndvi"   # LT4-1984087ndvi
 o=beast.irreg(ohio$ndvi, time=TIME,deltat=1/12)   
 

######################################################################################
# 'time' is given as  a list object


 TIME = list()
 
 TIME$year  = ohio$Y
 TIME$month = ohio$M
 TIME$day   = ohio$D
 o=beast.irreg(ohio$ndvi, time=TIME,deltat=1/12)   
 
 TIME = list() 
 TIME$year  = ohio$Y
 TIME$doy   = ohio$doy 
 o=beast.irreg(ohio$ndvi, time=TIME, deltat=1/12)    
 

}

}



\references{
	\enumerate{
		\item Zhao, K., Wulder, M.A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick, B., Zhang, X. and Brown, M., 2019. Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment, 232, p.111181 (the beast algorithm paper).
		\item Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection. Remote Sensing of Environment, 132, pp.102-119 (the Bayesian MCMC scheme used in beast). 
		\item Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine. ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261(a beast application paper). 
	}
}

\keyword{misc}
