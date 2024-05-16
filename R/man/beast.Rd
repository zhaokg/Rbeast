\name{beast}
\docType{methods}
\alias{beast}
\alias{Beast}
\alias{BEAST}
\alias{rbeast}
\alias{Rbeast}
\alias{RBEAST}
\title{Bayesian changepoint detection detection and time series decomposition for trend, periodicity or seasonality, and abrupt changes}
\description{
   A Bayesian model averaging algorithm called BEAST to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations. BEAST is useful for changepoint detection (e.g., breakpoints or structural breaks), nonlinear trend analysis, time series decomposition, and time series segmentation.
   \cr \cr
   \if{html}{\figure{Nile.jpg}{options: width="95\%" alt="Figure: Nile.jpg"}}
   \if{latex}{\figure{Nile.pdf}{options: width=10cm}}
  }
\usage{
  beast(   
         y, 
         start           = 1,
         deltat          = 1, 
         season          = c("harmonic", "svd", "dummy", "none"), 
         period          = NULL,  		   
         scp.minmax      = c(0,10), sorder.minmax   = c(0,5),  	  
         tcp.minmax      = c(0,10), torder.minmax   = c(0,1), 	   
	     sseg.min        = NULL,    sseg.leftmargin = NULL,  sseg.rightmargin = NULL,
	     tseg.min        = NULL,    tseg.leftmargin = NULL,  tseg.rightmargin = NULL, 
         method          = c('bayes','bic', 'aic','aicc','hic'),		 
         detrend         = FALSE, 
         deseasonalize   = FALSE,
         mcmc.seed       = 0,      
         mcmc.burnin     = 200, 
         mcmc.chains     = 3,
         mcmc.thin       = 5,
         mcmc.samples    = 8000,
         ci              = FALSE,
         precValue       = 1.5,
         precPriorType   = c('componentwise','uniform','constant','orderwise'),
         print.options   = TRUE,
         print.progress  = TRUE,		 
         quiet           = FALSE,		 
         gui             = FALSE,
         ...
      )
}
\arguments{

  \item{y}{ a vector for an evenly-spaced regular time series. Missing values such as NA and NaN are allowed.
			\itemize{
				\item  If \code{y} is irregular or unordered in time (e.g., multiple years of daily data spanning across leap years: 365 points in some years, and 366 in others), use the \code{\link[Rbeast]{beast.irreg}} function instead.
				
				\item  If \code{y} is a matrix or 3D array consisting of multiple regular or irregular time series (e.g., stacked images), use \code{\link[Rbeast]{beast123}} instead.
				
				\item  If \code{y} is an object of class 'ts','xts', or 'zoo', its time attributes (i.e.,start, end, frequency) will be used to specify the next several args such as   \code{start},\code{detlta},\code{period}, and \code{season}: No need to provide them explicitly; even if provided, the values are ignored to honor the time attributes of \code{y}. For example, if \code{y} has a frequency = 1, \code{season = 'none'} is always assumed; if \code{y} has a frequency > 1 (i.e., with a periodic component) but \code{season='none'} is specified by the user, 'none' will be replaced by 'harmonic'.
            }  

			If a list of multiple time series is provided for \code{y}, the multivariate version of the BEAST algorithm will be invoked to decompose the multiple time series and detect common changepoints altogether. This feature is \code{experimental} only and under further development. Check \code{\link[Rbeast]{ohio}} for a working example.			
  }
  

  \item{start}{ numeric (default to 1.0) or \code{Date}; the time of the 1st datapoint of \code{y}. It can be specified as a scalar (e.g., 2021.0644),  a vector of three values in the order of Year, Month, and Day (e.g., \code{c(2021,1,24)} ), or a R's Date object (e.g., \code{as.Date('2021-1-24')} ).}
  
  \item{deltat}{ numeric (default to 1.0) or string; the time interval between consecutive data points. Its unit should be consistent with \code{start}. If \code{start} takes a numeric scalar, the unit is arbitrary and irrelevant to beast (e.g., 2021.3 can be of any unit: Year 2021.3, 2021.3 meters, 2021.3 degrees ...). If \code{start} is a vector of Year, Month, and Day or an R's Date, \code{deltat} has the unit of YEAR. For example, if \code{start=c(2021,1,24)} for a monthly time series, \code{start} is converted to a fractional year 2021+(24-0.5)/365=2021.0644 and \code{deltat=1/12} needs to be set in order to specify the monthly interval. Alternatively,  \code{deltat} can be provided as a string to specify whether its unit is day, month, or year. Examples include '7 days', '7d', '1/2 months', '1 mn', '1.0 year', and '1y'. }

  \item{season}{ characters (default to 'harmonic'); specify if \code{y} has a periodic component or not. Four strings are possible.
				\itemize{
					\item  \code{'none'}: \code{y} is trend-only; no periodic components are present in the time series. The args for the seasonal component (i.e.,\code{sorder.minmax}, \code{scp.minmax} and \code{sseg.max}) will be irrelevant and ignored. 
					
					\item  \code{'harmonic'}: \code{y} has a periodic/seasonal component. The term \code{season} is a misnomer, being used here to broadly refer to any periodic variations present in \code{y}. The periodicity is NOT a model parameter estimated by BEAST but a known constant given by the user through \code{freq}. By default, the periodic component is modeled as a harmonic curve--a combination of sins and cosines.
					
					\item  \code{'dummy'}: the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a non-parametric curve. The harmonic order arg \code{sorder.minmax} is irrelevant and is ignored.
					\item  \code{'svd'}: (experimental feature) the same as \code{'harmonic'} except that the periodic/seasonal component is modeled as a linear combination of function bases derived from a Single-value decomposition. The SVD-based basis functions are more parsimonious than the harmonic sin/cos bases in parameterizing the seasonal variations; therefore, more subtle changepoints are likely to be detected.
				}  
  
  }
  
  \item{period}{ numeric or string. Specify the period for the periodic/seasonal component in \code{y}. Needed only for data with a periodic/cyclic component (i.e., season=\code{'harmonic'} or \code{'dummy'}) and not used for trend-only data (i.e., \code{season='none'}).  The period of the cyclic component should have a unit consisent with the unit of \code{deltat}. It holds that \code{period=deltat*freq} where \code{freq} is the number of data samples per period. (Note that the \code{freq} argument in earlier versions becomes obsolete and now is replaced by \code{period}. \code{freq} is still supported but\code{period} takes precedence if both are provided.)  \code{period} or the number of data points per period is not a BEAST model parameter and it has to be specified by the user. But if \code{period} is missing, BEAST first attempts to guess its value via auto-correlation before fitting the model. If \code{period} <= 0, \code{season='none'} is assumed, and  the trend-only model is fitted without a seasonal/cyclic component. If needed, use a string to specify whether the unit of period is day, month, or year. Examples are '1.0 year', '12 months', '365d', '366 days'.}

  \item{scp.minmax}{ a vector of 2 integers (>=0); the min and max number of seasonal changepoints (scp) allowed in segmenting the seasonal component. \code{scp.minmax} is used only if \code{y} has a seasonal component (i.e., \code{season='harmonic' or 'dummy' }) and ignored for trend-only data. If the min and max changepoint numbers are equal, BEAST assumes a constant number of scp and won't infer the posterior probability of the number of changepoints, but it still estimates the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur). If both the min and max numbers are set to 0, no changepoints are allowed; then a global harmonic model is used to fit the seasonal component, but still, the most likely harmonic order will be inferred if sorder.minmax[1] is not equal to sorder.minmax[2].}
  
  \item{sorder.minmax}{ a vector of 2 integers (>=1); the min and max harmonic orders considered to fit the seasonal component. \code{sorder.minmax} is used only used if the time series has a seasonal component (i.e., \code{season='harmonic'}) and ignored for trend-only data or when \code{season='dummy'}. If the min and max orders are equal (\code{sorder.minmax[1]=sorder.minmax[2]}), BEAST assumes a constant harmonic order used and won't infer the posterior probability of harmonic orders. }
  
  %%%%%%%%%%%%%%%%%%%%%	
   
  \item{tcp.minmax}{ a vector of 2 integers (>=0); the min and max number of trend changepoints (tcp) allowed in segmenting the trend component. If the min and max changepoint numbers are equal, BEAST assumes a constant number of changepoints and won't infer the posterior probability of the number of changepoints for the trend, but it still estimates the occurrence probability of the changepoints over time (i.e., the most likely times at which these changepoints occur in the trend). If both the min and max numbers are set to 0, no changepoints are allowed; then a global polynomial trend is used to fit the trend component, but still, the most likely polynomial order will be inferred if torder.minmax[1] is not equal to torder.minmax[2].}

  \item{torder.minmax}{ a vector of 2 integers (>=0); the min and max orders of the polynomials considered to fit the trend component. The 0-th order corresponds to a constant term/a flat line and the 1st order is a line. If \code{torder.minmax[1]=torder.minmax[2]}, BEAST assumes a constant polynomial order used and won't infer the posterior probability of polynomial orders. }	
  
  %%%%%%%%%%%%%%%%%%%%%	
  
  \item{sseg.min}{ an integer (>0); the min segment length allowed between two neighboring season changepoints. That is, when fitting a piecewise harmonic seasonal model, two changepoints are not allowed to occur within a time window of length \code{sseg.min}. \code{sseg.min} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{sseg.min*deltat}. \code{sseg.min} defaults to NULL and its value will be given a default value in reference to \code{freq}. }	
  
  \item{sseg.leftmargin}{ an integer (>=0); the number of leftmost data points excluded for seasonal changepoint detection. That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the starting window/segment of length \code{sseg.leftmargin}. \code{sseg.leftmargin} must be an unitless integer--the number of time intervals/data points so that the time window in the original unit is \code{sseg.leftmargin*deltat}. If missing, \code{sseg.leftmargin} defaults to \code{sseg.min}.}
  
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
		} %   \itemize
  } %  \item{method}
     
  \item{detrend}{ logical; If \code{TRUE}, a global trend is first fitted and removed from the time series before running BEAST; after BEAST finishes, the global trend is added back to the BEAST result. }
  
  \item{deseasonalize}{ logical; If \code{TRUE}, a global seasonal model is first fitted and removed from the time series before running BEAST; after BEAST finishes, the global seasonal curve is added back to the BEAST result. \code{deseasonalize} is ignored if \code{season='none'} (i.e., trend-only data). }		      
  
  \item{mcmc.seed}{ integer (>=0); the seed for the random number generator used for Monte Carlo Markov Chain (mcmc).  If \code{mcmc.seed=0}, an arbitrary seed is picked and the fitting results vary across runs. If fixed to the same non-zero integer, the result can be re-produced for different runs. But the results from the same seed may still vary if run on different computers because the random generator library depends on CPU's instruction sets.}	 
 
  \item{mcmc.chains}{ integer (>0); the number of  MCMC chains.}
  
  \item{mcmc.thin}{ integer (>0); a factor to thin chains (e.g., if thinningFactor=5, samples will be taken every 3 iterations)}
  
  \item{mcmc.burnin}{ integer (>0); the number of burn-in samples discarded at the start of each chain}
  
  \item{mcmc.samples}{ integer (>=0); the number of samples collected per MCMC chain. The total number of iterations is \code{(burnin+samples*thin)*chains}.}

  \item{ci}{ boolean; If \code{TRUE}, credible intervals (i.e., out$season$CI or out$trend$CI) will be computed for the estimated seasonal and trend components. Computing CI is time-consuming, due to sorting, so set \code{ci} to FALSE if a symmetric credible interval (i.e., out$trend$SD and out$season$SD) suffices.}	
  
%%%%%%%%%%%%%%%%%%%%%	
  
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
  
  \item{gui}{ boolean. If \code{TRUE}, BEAST will be run with a GUI window to show an animation of the MCMC sampling in the model space step by step; as an experimental feature, "\code{gui=TRUE}" works only for Windows x64 systems not Windows 32 or Linux/Mac.}	
 	
  \item{\dots}{additional parameters. There are many more settings for the implementation but not made available in the beast() interface; please use the function \code{beast123()} instead }
 
}

\value{
  The output is an object of class "beast". It is a list, consisting of the following variables. Its structure is the same as the outputs from the other two alternative functions \code{\link[Rbeast]{beast.irreg}} and \code{\link[Rbeast]{beast123}}. In the explanations below, we assume the input \code{y} is a single time series of length \code{N}:

  \item{time}{a vector of size \code{1xN}: the times at the \code{N} sampled locations. By default, it is simply set to \code{1:N} if the input arguments \code{delta}, \code{start}, or \code{time} are missing.}
  \item{data}{a vector, matrix, or 3D array; this is a copy of the input \code{y} if extra$dumpInputData = TRUE. If extra$dumpInputData=FALSE, it is set to NULL. If the original input \code{y} is irregular (as in \code{\link[Rbeast]{beast.irreg}}), the copy here is the regular version aggregated from the original at the time interval specified by deltat (in \code{\link[Rbeast]{beast.irreg}} or metadata$deltaTime (in \code{\link[Rbeast]{beast123}}). }
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
		 \item \code{cpOccPr}: [ChangePoint OCCurence PRobability]. a vector of length N; it gives a probability distribution of having a changepoint in the trend at each point of time. Plotting \code{cpOccPr} will depict a continious curve of probability-of-being-changepoint. \emph{ \bold{ Of particular note, in the curve, a higher peak indicates a higher chance of being a changepoint only at that particular SINGLE point in time and does not necessarily mean a higher chance of observing a changepoint AROUND that time. For example, a window of cpOccPr values \code{c(0,0,0.5,0,0)} (i.e., the peak prob is 0.5 and the summed prob is 0.5) is less likely to be a changepoint compared to another window \code{c(0.1,0.2,0.21,0.2,0.1)} (i.e., the peak prob is 0.21 but the summed prob is 0.71}}).
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
		 \item \code{Y}: a vector of length N; the estimated seasonal component. It is the Bayesian model averaging of all the individual sampled signal.
		 \item \code{SD}: [Standard Deviation] a vector of length N; the estimated standard deviation of the estimated seasonal component.
		 \item \code{CI}: [Standard Deviation] a matrix of dimension \code{N x 2}; the estimated credible interval of the estimated seasonal signal. One vector of the matrix is for the upper envelope and another for the lower envelope.
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
  \code{\link[Rbeast]{beast}},    \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{tetris}},   \code{\link[Rbeast]{geeLandsat}} 
}

 
\examples{

library(Rbeast)
\donttest{

#------------------------------------Example 1----------------------------------------#
# 'googletrend_beach' is the monthly Google Trend popularity of searching for 'beach' 
# in the US from 2004 to 2022. Sudden changes in the time series coincide with known 
# extreme  weather events (e.g., 2006 North American Blizzard, 2011 US hottest summer 
# on record, Record warm January in 2016) or the covid19 outbreak.

 out <- beast(googletrend_beach)
 
 plot(out)                          
 plot(out, vars=c('t','slpsgn') ) # plot the trend and probability of slope sign only.
                                  # In the slpsgn panel, the upper red portion refers to
                                  # probability of trend slope being positive, the middle 
                                  # green to the prob of slope being zero, and the lower 
                                  # blue to the probability of slope being negative.
                                  # Run "?plot.beast" for details on the plot function.

#------------------------------------Example 2----------------------------------------#
# Yellowstone is a half-monthly satellite time series of 774 NDVI(vegetation greeness) 
# observations starting from July 1-15,1981(i.e., start=c(1981,7,7)) at a Yellowstone
# forest site. It has 24 data points per year (i.e., freq=24). Note that the beast 
# function  hanldes only evenly-spaced regular time series. Irregular data need to be
# first  aggegrated at a regular time interval of your choice--the aggregation 
# functionality is implemented in beast.irreg() and beast123().

 data(Yellowstone)
 plot( 1981.5+(0:773)/24, Yellowstone, type='l')  # A sudden drop in greenness in 1988 
                                                  # due to the 1988 Yellowstone Fire

# Yellowstone is not a object of class 'ts' but a pure vector without time attributes. 
# Below, no extra argument is supplied, so default values (i.e.,start=1, deltat=1) are 
# used and the time is 1:774. 'period' is missing and so is guessed via auto-correlation.
# Use of auto-correlation to compute the period of a cyclic time series is not always 
# reliable, so it is suggested to always supply 'period' directly, as in Example 2 and
# Example 3.

 o = beast(Yellowstone)   # By defualt, the times assumed to be 1:length(Yellowstone) 
                          # and a periodic component is assumed (season='harmonic')
 plot(o)
 
#o = beast(Yellowstone, quiet=TRUE)                        # print no warning messages
#o = beast(Yellowstone, quiet=TRUE, print.progress=FALSE)  # print nothing

#------------------------------------Example 3----------------------------------------#
# The time info such as start,delta,and period is explicitly provided. 'start' can be 
# given as (1) a fractional number, (2) a vector comprising year, month,& day, or (3)
# a R's Date. In (1), the unit of start and deltat does not necessarily refer to time and can 
# be arbitrary (e.g., a sequence of data observed at evenly-spaced distaces along a 
# transect or a elevation gradient)

 # (1) Unknown unit such that 1981.5137 can be interpreted arbitrarily
 o=beast(Yellowstone, start=1981.5137,            deltat=1/24,  period=1.0) 
 
 # Use a string to explictly specify a time unit so that times are intepreted as dates
 # o=beast(Yellowstone, start=1981.5137, deltat='1/24 year', period=1.0) # 1.0 = 1 yr
 # o=beast(Yellowstone, start=1981.5137, deltat='0.5 mon',   period=1.0) # 1.0 = 1 yr
 # o=beast(Yellowstone, start=1981.5137, deltat=1/24, period='1 yr')    # 1/24 = 1/24 yr
 # o=beast(Yellowstone, start=1981.5137, deltat=1/24, period='365 days')# 1/24 = 1/24 yr
  
 # (2) start is provided as YMD, the unit is year: deltat=1/24 year=0.5 month
 # o=beast(Yellowstone, start=c(1981,7,7),         deltat=1/24,  period=1.0)

 # (3) start is provided as Date, the unit is year: deltat=1/24 year=0.5 month
 #o=beast(Yellowstone, start=as.Date('1981-7-7'), deltat=1/24, period=1.0)  

 print(o)                           # o is a R LIST object with many fields
 str(o)                             # See a list of fields in o
 
 plot(o)                            # plot many variables 
 plot(o, vars=c('y','s','t') )      # plot the Y, seasonal, and trend components only
 plot(o, vars=c('s','scp','samp','t','tcp','tslp'))# Plot some selected variables in 
                                                   # 'o'. Type "?plot.beast" to see 
                                                   # more about vars
 plot(o, vars=c('s','t'),col=c('red','blue') )     # Specify colors of selected subplots

 plot(o$time, o$season$Y,type='l') # directly plot output: the fitted season
 plot(o$time, o$season$cpOccPr)    # directly plot output: season chgpt  prob
 plot(o$time, o$trend$Y,type='l')  # directly plot output: the fitted trend
 plot(o$time, o$trend$cpOccPr)     # directly plot output: trend chgpt occurrence prob
 plot(o$time, o$season$order)      # directly plot output: avg harmonic order
 
 plot(o, interactive=TRUE)         # manually choose which variables to plot

}

\donttest{

#------------------------------------Example 4----------------------------------------#
# Specify other arguments explicitly.  Default values are used for missing parameters.
# Note that beast(), beast.irreg(), and beast123() call the same internal C/C++ library,
# so in beast(), the input parameters will be converted to metadata, prior, mcmc, and 
# extra parameters as explained for the beast123() function. Or type 'View(beast)' to 
# check the parameter assignment in the code.
 
 
 # In R's terminology, the  number of datapoints per period is also called 'freq'. In this
 # version, the 'freq' argument is obsolete and replaced by 'period'.
 
 # period=deltat*number_of_datapoints_per_period = 1.0*24=24 because deltat is set to 1.0 by
 # default and this signal has 24 samples per period. 
 out = beast(Yellowstone, period=24.0, mcmc.samples=5000, tseg.min=20)
 
 # period=deltat*number_of_datapoints_per_period = 1/24*24=1.0.
 # out = beast(Yellowstone, deltat=1/24 period=1.0, mcmc.samples=5000, tseg.min=20)
  
 out = beast( 
     Yellowstone,            # Yellowstone: a pure numeric vector wo time info
     start   = 1981.51, 
     deltat  = 1/24,         
     period  = 1.0,           # period=delta*number_of_datapoints_per_period
     season  = 'harmonic',    # periodic compnt exisits,fitted as a harmonic curve 
     scp.minmax     = c(0,3), # min and max numbers of seasonal changpts allowed
     sorder.minmax  = c(1,5), # min and max harmonic orders allowed
     sseg.min       = 24,     # the min length of segments btw neighboring chnpts
	                          # '24' means 24 datapoints; the unit is datapoint.
     sseg.leftmargin= 40,     # no seasonal chgpts allowed in the starting 40 datapoints
     tcp.minmax     = c(0,10),# min and max numbers of changpts allowed in the trend
     torder.minmax  = c(0,1), #  min and maxx polynomial orders to fit trend
     tseg.min       = 24,     # the min length of segments btw neighboring trend chnpts
     tseg.leftmargin= 10,     # no trend chgpts allowed in the starting 10 datapoints
     deseasonalize  = TRUE, # remove the global seasonality before fitting the beast model
     detrend        = TRUE, # remove the global trend before fitting the beast model
     mcmc.seed      = 0,    # a seed for mcmc's random nummber generator; use a
                           # non-zero integer to reproduce results across runs
     mcmc.burnin    = 500, # number of initial iterations discarded
     mcmc.chains    = 2,    # number of chains
     mcmc.thin      = 3,    # include samples every 3 iterations
     mcmc.samples   = 6000  # number of samples taken per chain
                           # total iteration: (500+3*6000)*2		 
     )
 plot(out)
 plot(out, interactive=TRUE)
 }
 
#------------------------------------Example 5----------------------------------------#
# Run an interactive GUI to visualize how BEAST is samplinig from the possible model 
# spaces in terms of the numbers and timings of seasonal and trend changepoints.
# The GUI inferface allows changing the option parameters interactively. This GUI is 
# only available on Win x64 machines, not Mac or Linux.

\dontrun{
 beast(Yellowstone, period=24, gui=TRUE) 
}

#------------------------------------Example 6----------------------------------------#
# Apply beast to trend-only data. 'Nile' is the ANNUAL river flow of the river
# Nile at Aswan since 1871. It is a 'ts' object; its time attributes (start=1871, 
# end=1970,frequency=1) are used to replace the user-supplied start,deltat, and freq, 
# if any. 

\donttest{
 data(Nile)  
 plot(Nile)     
 attributes(Nile) # a ts object with time attributes (i.e., tsp=(start,end,freq)
 
 o = beast(Nile)  # start=1871, delta=1, and freq=1 taken from Nile itself
 plot(o)
 
 o = beast(Nile,             # the same as above. The user-supplied values (i.e., 2023,
           start=2023,       # 9999) are ignored bcz Nile carries its own time attributes.
           period=9999,      # Its frequency tag is 1 (i.e., trend-only), so season='none'
           season='harmonic' # is used instead of the supplied 'harmonic'
		   )
} 
 
#------------------------------------Example 7----------------------------------------#
# NileVec is  a pure data vector. The first run below is WRONG bcz NileVec was assumed
# to have a perodic component by default and beast gets a best estimate of freq=6 while 
# the true value is freq=1. To fit a trend-only model, season='none' has to be explicitly
# specified, as in the 2nd & 3rd funs.
\donttest{
 NileVec = as.vector(Nile) # NileVec is not a ts obj but a pure numeric data vector
 o       = beast(NileVec)  # WRONG WAY to call: No time attributes available to interpret 
                           # NileVec. By default, beast assumes season='harmonic', start=1,
                           # & deltat=1. 'freq' is missing and guessed to be 6 (WRONG).    
						   
 plot(o)                   # WRONG Results: The result has a suprious seasonal component 
							  
 o=beast(NileVec,season='none') # The correct way to call: Use season='none' for trend-only 
                                # analysis; the default time is the integer indices
                                # "1:length(NileVec)'. 
 print(o$time)							 
								
 o=beast(NileVec,               # Recommended way to call: The true time attributes are 
         start  = 1871,         # given explicitly through start and deltat (or freq if 
         deltat = 1,            # there is a  cyclic/seasonal cmponent). 
         season = 'none')  
 print(o$time)			 
 plot(o)
}


#------------------------------------Example 8----------------------------------------#
# beast can handle missing data. co2 is a monthly time series (i.e.,freq=12) starting
# from Jan 1959. We generate some missing values at random indices
 
\dontrun{

 data(co2)  
 attributes(co2)                          # A ts object with time attributes (i.e., tsp)
 badIdx      = sample( 1:length(co2), 50) # Get a set of random indices
 co2[badIdx] = NA                         # Insert some data gaps   

 out=beast(co2) # co2 is a ts object and its 'tsp' time attributes are used to get the
                # true time info. No need to specify 'start','deltat', & freq explicity.
				
 out=beast(co2,                  # The supplied time/period values will be ignored bcz
           start  = c(1959,1,15),# co2 is a ts object; the correct period = 1 will be 
           deltat = 1/12,        # used.
           period = 365)  
 print(out)
 plot(out)
} 



#------------------------------------Example 9----------------------------------------#
# Apply beast to time seris-like sequence data: the unit of sequences is not 
# necessarily time.
 
\donttest{
  data(CNAchrom11) # DNA copy number alterations in Chromesome 11 for cell line GM05296
                   # The data is orderd by genomic position (not time), and the values
                   # are the log2-based intensity ratio of copy numbers between the sample
                   # the reference. A value of zero means no gain or loss in copy number.
  o = beast(CNAchrom11,season='none') # season is a misnomer here bcz the data has nothing
                                      # to do with time. Regardless, we fit only a trend.
  plot(o)									  
 
} 


#------------------------------------Example 10---------------------------------------#
# Apply beast to time seris-like data: the unit of sequences is not necessarily time.
 
\donttest{
  # Age of Death of Successive Kings of England
  # If the data link is deprecated, install the time series data library instead,
  # which is available at https://pkg.yangzhuoranyang.com/tsdl/
  # install.packages("devtools")
  # devtools::install_github("FinYang/tsdl")
  # kings = tsdl::tsdl[[293]]
  
  kings = scan("http://robjhyndman.com/tsdldata/misc/kings.dat",skip=3)
  out   = beast(kings,season='none')
  plot(out) 
  
} 

 
#------------------------------------Example 11---------------------------------------#
# Another example from the tsdl data library
 

\donttest{
  # Number of monthly births in New York from Jan 1946 to Dec 1959
  # If the data link becomes invalid, install the time series data package instead
  # install.packages("devtools")
  # devtools::install_github("FinYang/tsdl")
  # kings = tsdl::tsdl[[534]]
  
  births = scan("http://robjhyndman.com/tsdldata/data/nybirths.dat") 
  out    = beast(births,start=c(1946,1,15), deltat=1/12 )  
  plot(out) # the result is wrong bcz the guessed freq via auto-correlation by beast 
            # is 2 rather than 12, so we recommend always specifying 'freq' explicitly
            # for those time series with a periodic component, as shown below.
  out    = beast(births,start=c(1946,1,15), deltat=1/12, freq  =12 )  
  out    = beast(births,start=c(1946,1,15), deltat=1/12, period=1.0 )  
  plot(out)  
  
}

#------------------------------------Example 12---------------------------------------#
#    Daily confirmed COVID-19 new cases and deaths across the globe
 
 \dontrun{
 data(covid19)
 plot(covid19$date, covid19$newcases, type='l')
 
 newcases = sqrt( covid19$newcases )  # Apply a square root-transformation
 
 # This ts varies periodically every 7 days. 7 days can't be precisely represented 
 # in  the unit of year bcz some years has 365 days and others has 366. BEAST can hanlde 
 # this in two ways.


 #(1) Use the date number as the time unit--the num of days lapsed since 1970-01-01. 
  
  datenum  = as.numeric(covid19$date) 
  o        = beast(newcases, start=min(datenum), deltat=1, period=7) 
  o$time   = as.Date(o$time, origin='1970-01-01') # Convert from integers to Date.
  plot(o)
 
 #(2) Use strings to explicitly specify deltat and period with a unit. 
 
  startdate = covid19$date[1]
  o         = beast(newcases, start=startdate, deltat='1day', period='7days') 
  plot(o)
 
 }
 
#------------------------------------Example 13---------------------------------------#
# The old API interface of beast is still made available but NOT recommended. It is 
# kept mainly to ensure the working of the sample code on Page 475 in the text
# Ecological Metods by Drs. Southwood and Henderson.

\dontrun{

  # The interface as shown here will be deprecated and NOT recommended.
  beast(Yellowstone, 24)  #24 is the freq: number of datapoints per period
  
   
  # Specify the model or MCMC parameters through opt as in Rbeast v0.2
  opt=list()             #Create an empty list to append individual model parameters
  opt$period=24          #Period of the cyclic component (i.e.,freq in the new version)
  opt$minSeasonOrder=2   #Min harmonic order allowed in fitting season component
  opt$maxSeasonOrder=8   #Max harmonic order allowed in fititing season component
  opt$minTrendOrder=0    #Min polynomial order allowed to fit trend (0 for constant)
  opt$maxTrendOrder=1    #Max polynomial order allowed to fit trend (1 for linear term)
  opt$minSepDist_Season=20#Min separation time btw neighboring season changepoints 
  opt$minSepDist_Trend=20 #Min separation time btw neighboring trend  changepoints
  opt$maxKnotNum_Season=4 #Max number of season changepoints allowed 
  opt$maxKnotNum_Trend=10 #Max number of trend changepoints allowed  
  opt$omittedValue=NA    #A customized value to indicate bad/missing values in the time 
                          #series, in additon to those NA or NaN values.
  					
  # The following parameters used to configure the reverisible-jump MCMC (RJMCC) sampler
  opt$chainNumber=2       #Number of parallel MCMC chains 
  opt$sample=1000         #Number of samples to be collected per chain
  opt$thinningFactor=3    #A factor to thin chains  
  opt$burnin=500          #Number of burn-in samples discarded at the start 
  opt$maxMoveStepSize=30  #For the move proposal, the max window allowed in jumping from 
                           #the current changepoint
  opt$resamplingSeasonOrderProb=0.2 #The probability of selecting a re-sampling proposal 
                                    #(e.g., resample seasonal harmonic order)
  opt$resamplingTrendOrderProb=0.2  #The probability of selecting a re-sampling proposal 
                                    #(e.g., resample trend polynomial order)
								   
  opt$seed=65654    #A seed for the random generator: If seed=0,random numbers differ
                    #for different BEAST runs. Setting seed to a chosen non-zero integer 
                    #will allow reproducing the same result for different BEAST runs.
 
  beast(Yellowstone, opt)  
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
