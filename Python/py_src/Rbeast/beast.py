from .             import Rbeast as cbeast              # Rbeast here is the Rbeast-xxx.pyd C extension module
from .cvt_to_numpy import force_convert_to_numpy_ndarray

def beast(Y,                     \
          start            = 1,
          deltat           = 1,
          season           = 'harmonic',  # 'harmonic','dummy','svd','none'
          period           = float('nan'),
          scp_minmax       = [0, 10],     # the min and max numbers of seasonal changepoints
          sorder_minmax    = [0, 5],
          sseg_minlength   = None,    # an integer
          sseg_leftmargin  = None,    # an integer
          sseg_rightmargin = None,    # an integer
          tcp_minmax       = [0, 10], # the min and max numbers of trend changepoints
          torder_minmax    = [0, 1 ],
          tseg_minlength   = None,    # an integer
          tseg_leftmargin  = None,    # an integer
          tseg_rightmargin = None,    # an integer
          s_complexfct     = 0.0,     # control the complexity of the fitted seasonal curve
          t_complexfct     = 0.0,     # control the complexity of the fitted trend curve
          method           = 'bayes', # 'bayes','bic','aic','aicc','hic', 'bic0.25','bic0.5'
          detrend          = False,
          deseasonalize    = False,
          mcmc_seed        = 0,
          mcmc_burbin      = 200,
          mcmc_chains      = 3,
          mcmc_thin        = 5,
          mcmc_samples     = 8000,
          precValue        = 1.5,
          precPriorType    = 'componentwise',  # componentwise','uniform','constant','orderwise'
          hasOutlier       = False,
          ocp_minmax       = [0, 10],	  # the min and max numbers of outlier-type changepoints	  
          print_param      = True,
          print_progress   = True,
          print_warning    = True,
          quiet            = False,
          gui              = False,
          dump_ci          = False,
          dump_mcmc        = False,		  
          **kwargs ):
      """
      
      
######################################################################################################
Bayesian changepoint detection and time series decomposition for regular or irregular time series data
    
The fitted model is:

  Y = trend + error             if data has no periodic/seasonal variation (i.e., season='none')
  Y = trend + seasonal + error  if data has periodic/seasonal variation 
  Y = trend + outlier  + error  if data is trend-only (no seasonal variation) but with potential outliers
  Y = trend + seasonal + outlier + error if data has periodic/seasonal variation and also has outliers

where trend is a piecewise linear or polynomial function with an unknown number of trend changepoints to be
nferred; seasonal is a piecewise periodic function with an unknown number of seasonal changepoints to be 
inferred; and the outlier component refers to potential anomalous spikes or dips at isolated data points and
is included only if hasOutlier=True (in beast or beast_irreg) or metadata.hasOutlierCmpnt=True (in beast123) 
######################################################################################################
  
--------------------------------------------------------------------------------------------------
*Quick Examples*:
--------------------------------------------------------------------------------------------------
import Rbeast as rb
    
nile,yr =  rb.load_example('nile')             # annual flow of the Nile river
rb.beast( nile, start=1871, season='none' )
    
beach, year = rb.load_example('googletrend')   # the google search trend for 'beach'
o = rb.beast(beach, start= 2004, deltat=1/12, period = 1.0)
rb.plot(o) 
    
--------------------------------------------------------------------------------------------------    
*Input arguments*:
--------------------------------------------------------------------------------------------------
Y:  a regular time series; it should be a numeric vector. For ireggular time series, use 
    'beast_irreg' or 'beast123' instead. For multiple time series or stacked time series 
    images such as satellite data, use 'beast123'.
 
... :  the remaining arguments are many paired keywords and values to specifiy time information 
    or parameters for the beast algorithm. Check the R version of BEAST for detailed explanations
    (https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf). Below is a brief description.

--------------------------------------------------------------------------------------------------  
*Possible Keywords*:
--------------------------------------------------------------------------------------------------   
start: 
         the start time of the regular time series. Possille values are
         (1) a numeric scalar (e.g., 2002.33); the unit is arbitrary 
         (2) a vector of two values (Year, Month) (e.g., [2002,4] ) 
         (3) a vector of three values (Year, Month, Day) (e.g., [2002,4,15] )
         (4) a date string (e.g., "2002-4-15", "2002/04/15", or "2002/4")
         (5) a datenum as a object variable ( e.g. rb.args( datenum=731320, origin='python'),
             rb.args(year=2002, month=4, day=15) )         
deltat: 
         a number or string; the time interval between consecutive datapoints. Its unit must be consistent with "start". 
         (1) If start is numeric, the unit is unkown and arbitrary and  will be ignored by BEAST (e.g., 2021.3 can be be of any unit, a time, a distance, etc.).
         (2) If start is a vector of [year, month, day] or a date string, deltat is interpreted in units of years. For example, for a monthly series with 
          start = [2021, 1, 24], start is internally converted to a fractional year 2021 + (24 - 0.5)/365 = 2021.0644 and deltat = 1/12} can be used to
          specify a 1-month interval (1/12 yr). Alternatively, use a string to specify deltat  (e.g., '1/12 year', '1.0 month', '30 days').
period:  
         a number or string to specify the period if peridodic/seasonal variations 
         are present in the data. If period is given a zero, negative value or 'none' 
         it suggests no seasonal/periodic component in the signal. (season='none'
         also suggests no periodic component).
         In earlier versions, 'freq' was used to specify the period and
         now deprecated in this version. If period is given as anumber, the unit of 'period', 
         if any, should be consistent with the unit of 'deltat'. If given as a string,
         the unit of period needs to be expicilty specified (e.g., '1 year', '12 mon', '365 days')
season: 
         a string specifier. Possible values - (1) 'none':  trend-only data with no 
         seasonality; (2)'harmonic': the seasonal/peridoic  component modelled via 
         harmonic curves; (3)'dummy': the seasonal component  modelled via a dummy 
         basis (i.e., pulse-like bases); (4)'svd': svd-derived  bases (experimental 
         feature)
scp_minmax: 
         a vector of two integers (e.g.,[0,5]); the min and max number of
         seasonal changepoints allowed
sorder_minmax: 
         a vector of two integers (e.g.,[1,3]); the min and max harmonic orders of
         seasonal changepoints (scp) allowed
sseg_minlength: 
         an integer; the min length of the segment for the seasonal component 
         i.e., the min distance between neighorbing changepoints)
sseg_leftmargin: 
         an integer;  the number of leftmost data points excluded for seasonal changepoint detection.
         That is,  no changepoints are allowed in the starting window/segment of length sseg_leftmargin. 
         sseg_leftmargin must be an unitless integer–the number of time intervals/data points so that the
         time window in the original unit is sseg_leftmargin*deltat. If missing, sseg_leftmargin defaults
         to the minimum segment length 'sseg_min'
sseg_rightmargin: 
         an integer;  the number of rightmost data points excluded for seasonal changepoint detection.
         That is,  no changepoints are allowed in the ending window/segment of length sseg_rightmargin. 
         sseg_rightmargin must be an unitless integer–the number of time intervals/data points so that the
         time window in the original unit is sseg_rightmargin*deltat. If missing, sseg_rightmargin defaults
         to the minimum segment length 'sseg_min'		 
tcp_minmax: 
         a vector of two integers (e.g.,[0,5]); the min and max numbers of
         trend changepoints (tcp) allowed
torder_minmax: 
         a vector of two integers (e.g.,[1,3]); the min and max orders of
         polynomials used to model the trend
tseg_minlength: 
         an integer; the min length of the segment for the trend component (i.e.,
         the min distance between neighorbing changepoints)
tseg_leftmargin: 
         an integer; the number of leftmost data points excluded for trend changepoint detection.
         That is,  no trend changepoints are allowed in the starting window/segment of length tseg_leftmargin. 
         tseg_leftmargin must be an unitless integer–the number of time intervals/data points so that the
         time window in the original unit is tseg_leftmargin*deltat. 
tseg_rightmargin: 
         an integer;  the number of rightmost data points excluded for trend changepoint detection.
         That is,  no trend changepoints are allowed in the ending window/segment of length tseg_rightmargin. 
         tseg_rightmargin must be an unitless integer–the number of time intervals/data points so that the
         time window in the original unit is tseg_rightmargin*deltat.
s_complexfct:
          Numeric (defaulted to 0.0); a hyperprior parameter--newly added in Version 0.1.24--controlling the complexity of
          the seasonal curve (i.e., the model dimension or the number of seasonal changepoints). A prior of the form 
          "p(k) ~ exp[lambda*(k+1)]" is placed on the number of seasonal changepoints k, where lambda is s_complexfct 
          (i.e., prior.seasonComplexityFactor in the beast123() function). Setting lambda = 0 (i.e., s_complexfct=0) yields 
          a non-informative prior "p(k) ~ 1.0" where all model dimensions are equally likely a priori. Users may tune s_complexfct
          (for beast and beast_irreg)  or seasonComplexityFactor (for beast123) in the range of [-20, 20]} or an
          even wider range: Negative values (e.g., lambda = -15.9) favor fewer changepoints (simpler seasonal curves), whereas 
          positive values (e.g., lambda = 5.76) favor more changepoints (more complex curves).
t_complexfct: 
          Numeric (defaulted to 0.0); analogous to s_complexfct, but for the trend component and the number of trend changepoints.
method: 
         a string specifying which method to formulat model posterior probability. Possible values are
         (1) 'bayes': the full Bayesian formulation (this is the default)  
         (2)'bic':  approximation of posterior probability using the Bayesian information criterion (bic)
         (3)'aic':  approximation of posterior probability using the Akaike information criterion (aic)
         (4)'aicc': approximation of posterior probability using the corrected Akaike information criterion (aicc)
         (5)'hic':  approximation of  posterior probability using the Hannan–Quinn information criterion  (hic)
         (6)'bic0.25':  approximation using the Bayesian information criterion adopted from Kim et al. (2016) <doi: 
              10.1016/j.jspi.2015.09.008>; bic0.25=n*ln(SSE)+0.25k*ln(n) with less complexity penelaty than the standard BIC.
         (7)'bic0.50': the same as above except that the penalty factor is 0.50.
         (8)'bic1.5':  the same as above except that the penalty factor is 1.5.
         (9)'bic2':    the same as above except that the penalty factor is 2.0.
deseasonalize: 
         boolean; if true, the input time series will be first de-seasonalized before applying
         beast by removing a global seasonal component
detrend: 
         boolean; if true, the input time series will be first de-trend before applying 
         beast by removing a global trend  
mcmc_seed: 
         a seed for the random number generator; set it to a non-zero integer to
         reproduce the results among different runs
mcmc_samples: 
         number of MCMC samples collected; the larger, the better
mcmc_thin: 
         a thinning factor for MCMC chains: take every 'mcmc.thin'-th sample
mcmc_burnin: 
         the number of initial samples of each chain to be discarded
mcmc_chains: 
         the number of MCMC chains; the larger, the better but with more computation. 
precValue:
         numeric (>0); the hyperparameter of the precision prior; the default value is 1.5. precValue
         is useful only when precPriorType='constant', as further explained below
precPriorType:
         a string taking one of 'constant', 'uniform',  'componentwise' (default), and 'orderwise'.
         (1) 'constant':  the precision parameter used to parameterize the model coefficients is fixed to
           a const specified by precValue. In other words, precValue is a user-defined hyperparameter 
           and the fitting result may be sensitive to the chosen values of precValue.
         (2) 'uniform':  the precision parameter used to parameterize the model coefficients is a random variable;
           its initial value is specified by precValue. In other words, precValue will be inferred by the MCMC,
           so the fitting result will be insensitive to the chose inital value of precValue.
         (3) 'componentwise': multiple precision parameters are used to parameterize the model coefficients for
           individual components (e.g., one for season and another for trend); their initial values is specified 
           by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result will be 
           insensitive to the choice in precValue.
         (4) 'orderwise'}: multiple precision parameters are used to parameterize the model coefficients not just for 
           individual components but also for individual orders of each component; their initial values is specified 
           by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result will be 
           insensitive to the choice in precValue. 
hasOutlier:
        boolean; if true, the model with an outlier component will be fitted (if season='none',
        Y=trend+outlier+error, or if season ~= 'none', Y=trend+season+outlier+error).		
ocp_minmax:
        a vector of 2 integers (>=0); the min and max numbers of outlier-type changepoints (ocp) allowed in the time series.
        Ocp refers to spikes or dips at isolated times that can't be modeled as trends or seasonal terms.
print_param: 
         boolean; if true, print the beast paramers. 		 
print_progress: 
         boolean; if true, print a progress bar
print_warning: 
         boolean; if true, print warning messages
quiet:
         boolean; if true, supress all the messages and printing		
dump_ci: 
         boolean; if true, credible intervals (i.e., out.season.CI or out.trend.CI) will be computed 
         for the estimated seasonal and trend components. Computing CI is time-consuming, due to sorting, 
         so set dump_ci=Flase if a symmetric credible interval (i.e., out.trend.SD and out.season.SD) suffices.	  		 
dump_mcmc:
         boolean; if true, dump the sampled models in the MCMC chains
gui: 
         boolean; if true, show a gui to demostrate the MCMC sampling; runs only 
         on Windows not Linux or MacOS
		 
######################################################################################################
The keywords for beast() are converted to 'metadata', 'prior','mcmc', and 'extra' options used 
in the beast123() interface. Examples are:
         deseasonalize <-> metadata.deseasonalize
         hasOutlier    <-> metadata.hasOutlierCmpnt		 
         scp_minmax[0] <-> prior.seasonMinOrder
         scp_minmax[1] <-> prior.seasonMaxOrder
         sseg_min      <-> prior.seasonMinSepDist
         mcmc_seed     <-> mcmc.seed
         tcp_minmax[0] <-> prior.trendMinKnotNumber
         tcp_minmax[1] <-> prior.trendMaxKnotNumber
         dump_ci       <-> extra.computeCredible
Experts should use the the beast123 function.
######################################################################################################		 
 
-------------------------------------------------------------------------------------------------- 
*Result/Output*: The output is a struct variable; example of the fields include
--------------------------------------------------------------------------------------------------
        marg_lik: marginal likilood; the larger, the better
        sig2    : variance  of error
        trend   : the trend component; a struct variable (say, T)
        season  : the season componet; a stuct variable  (say,S)
        The subfields of trend or season:
        .ncpPr        : the prob distribution for number of changepoints
        .ncp          : mean number of changepoints in trend or seasonality
        .ncp_meidan   : median number of changepoints
        .ncp_mode     : mode from ncpPr
        .ncp_pct90    : 90% percentile from ncpPr
        .cpOccPr      : changepoint occurrance probability over time
        .cp           : list of all possible changepoints (many are not sigficant)
        .cpPr         : occurrence probability of the changepoints in cp
        .cpAbruptChange: the sudden changes in trend or seasonlity at cp
        .cpCI         : confidence interval of the cps
        .Y            : the fitted trend or seasonality 
        .SD           : standard deviation of the fitted Y
        .CI           : Credible interval of the fittted Y
        .order   : the mean harmonic or polynomial orders estimated to fit the seasonal and trend     
        trend.slp     : slope of the trend 
        trend.slpSD   : standard dev of the estimated slope
        trend.slpSgnPosPr: time-varying probability of the slope being postive
        trend.slpSgnZeroPr: time-varying probability of the slope being 0
        season.amp     : amplitue of the estiamted seasonality overtime
        season.ampSD   : standard ev of the estiamated amplitude

-------------------------------------------------------------------------------------------------- 
More help:  
--------------------------------------------------------------------------------------------------
       This terse help doc sucks (I know); so far, the best details are still the
       R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
       Python doesn't allow a '.' in variable names, so Python's equivalent to R's 
       beast(Y,start=1987,tcp.minmax=c(0,5)) is beast(Y, start=1987,  tcp_minmax=[0, 5]).
       
--------------------------------------------------------------------------------------------------
Examples:
--------------------------------------------------------------------------------------------------
import Rbeast as rb

# Nile river annual streamflow: trend-only data
Nile, yr=rb.load_example('nile')      
o = rb.beast(Nile, start=1871, season='none') 
rb.print(o)
rb.plot(o)

# Explicitly specify deltat=1. beast knows nothing about the unit
# of 1871 and 1.0 (i.e., 1871 years, 1871 seconds, or 1871 meters?) 
o=rb.beast(Nile, start=1871, deltat=1.0, season='none') 
 
# start is given a Date 1871-1 (Year-Mon). Then, the time unit is
# fractional/decimal year. delta=1.0 means 1.0 year
o=rb.beast(Nile, start=[1871,1], deltat=1.0, season='none') 
 
# period=0 means a trend-only signal, which is equivalent to season='none'
o=rb.beast(Nile, start=1871, deltat=1, period=0) 
 
# Use a string to specify a unit for deltat or period (e.g., deltat='1 year')
# The time unit is also fractional year. 1871 means Year 1871
o=rb.beast(Nile, start=1871, deltat='1 year', period=0) 
 
# Use a string to specify a unit for delta or period (e.g., deltat='12 mo')
# The time unit is fractional year. 1871 means Year 1871
o=rb.beast(Nile, start=1871, deltat='12 mo', period=0) 
 
#Do not print the options 
o=rb.beast(Nile, start=1871, deltat=1.0, season='none',print_options=False)
 
# Show a gui window to demostrate the beast sampling process in
# real-time (for Windows only not Linux and MacOS)% 
rb.beast(Nile,season='none', gui=True, start=1871) 
 

## Monthly google search popularity of the keyword 'beach'
beach, time = rb.load_example('googletrend')   

# If not specfified, season='harmonic' is assumed by default, with a seasonal 
# component fitted.
o = rb.beast(beach, start=[2004,1], deltat=1/12)  # deltat = 1/12 yr =1 month
rb.print(o)
rb.plot(o)
rb.plot(o,ncpStat='median')
 
# delta  = 1/12: for dates, the default unit is year, so delta=1/12yr=1 month;       
# period = 1.0 means 1 year
o=rb.beast(beach, start='2004-1-15', deltat=1/12, period=1.0)  
 
# period='12 month': use a string to explicitly specify the unit      
o=rb.beast(beach, start=[2004,1], deltat=1/12,      period='12 month')
o=rb.beast(beach, start=[2004,1], deltat='1 month', period='365 days')


# "t_complexfct" or "s_complexfct" is a hyperparameter controlling the complexity of
#  the fitted trend or seasonal curve. Postive values for them tend to give more changepoints
#  and negative values give fewer changepoints. These two parameters are introduced
#  in Rbeast v1.24.0. In earlier values, they were not expoosed in the function
#  API and were simply assumed to be zeros (i.e., non-informative priors on the model dimension/the number
# of changepoints). Again, postive values favor more changepoints, and negative values farvor fewer changepoints
 
o1 = rb.beast(beach, start=[2004,1], deltat=1/12, t_complexfct = 2.7  )   # More changepoints than the non-informative prior
o2 = rb.beast(beach, start=[2004,1], deltat=1/12, t_complexfct = -3.1  )  # Fewer changepoints than the non-informative prior
rb.plot(o1)
rb.plot(o2)

## Monthly air co2 data since 1959: deltaTime=1/12 year
co2, time = rb.load_example('co2')     
o = rb.beast(co2, start=[1959,1,15],  deltat=1/12,  period=1.0)
rb.print(o)
rb.plot(o)
rb.plot(o,ncpStat='median')
 
## Daily covid-19 infection statistics 
covid=rb.load_example('covid19')    
Y       = covid.newcases
datestr = covid.date
import numpy as np
Y       = np.sqrt(Y)     # a sqrt-root transformation

% the min length of seasonal segments is set to 30 data points
o=rb.beast(Y, start=[2020,1,22],  deltat=1/365,  period=7/365, sseg_min=30)
rb.print(o)
rb.plot(o)
 
# Use a string to specify delta with a unit
o=rb.beast(Y, start=[2020,1,22], deltat='1.0 day',  period='7days', sseg_min=30)
 
# Use a string to specify delta with a unit
o=rb.beast(Y, start=[2020,01,22], deltat='1.0 day',  period='7days', sseg_min=30)
 
Contact info: To report bug or get help, do not hesitate to contact Kaiguang Zhao
at zhao.1423@osu.edu.
      """

    # For simplicity, just convert Y into a numpy ndarray 
      Y = force_convert_to_numpy_ndarray(Y)

    # ......Start of  'MetaData' ......
      metadata                  = lambda: None   ###Just get an empty object###
      metadata.isRegularOrdered = True
      metadata.season           = season
      metadata.startTime        = start
      metadata.deltaTime        = deltat
      #metadata.whichDimIsTime   = 1
      metadata.period           = period
      if (season != 'none' and period != period and deltat == deltat):
         if 'freq' in kwargs:
            freq            = kwargs['freq']
            metadata.period = deltat * freq
      metadata.missingValue     = float('nan')
      metadata.maxMissingRate   = 0.7500
      metadata.deseasonalize    = deseasonalize
      metadata.detrend          = detrend
      metadata.hasOutlierCmpnt  = hasOutlier
    #........End of  MetaData ........

    #......Start of  'prior' ......
      prior                      = lambda: None   ###Just get an empty object###
      prior.modelPriorType	 = 1
      if season !='none' or season == None:
            prior.seasonMinOrder   = sorder_minmax[0]
            prior.seasonMaxOrder   = sorder_minmax[1]
            prior.seasonMinKnotNum = scp_minmax[0]
            prior.seasonMaxKnotNum = scp_minmax[1]
            prior.seasonMinSepDist = sseg_minlength
            prior.seasonLeftMargin  = sseg_leftmargin
            prior.seasonRightMargin = sseg_rightmargin
      prior.trendMinOrder	  = torder_minmax[0]
      prior.trendMaxOrder	  = torder_minmax[1]
      prior.trendMinKnotNum   = tcp_minmax[0]
      prior.trendMaxKnotNum   = tcp_minmax[1]
      prior.trendMinSepDist   = tseg_minlength
      prior.trendLeftMargin   = tseg_leftmargin
      prior.trendRightMargin  = tseg_rightmargin  
      if hasOutlier:
            prior.outlierMinKnotNum = ocp_minmax[0]	  	  
            prior.outlierMaxKnotNum = ocp_minmax[1]	  
      prior.K_MAX            = 0
      prior.precValue        = precValue
      prior.precPriorType    = precPriorType
      prior.seasonModelPriorFactor = s_complexfct;  # new parameter introduced in Rbeast v1.24.0
      prior.trendComplexityFactor  = t_complexfct;  # new parameter introduced in Rbeast v1.24.0  
    #......End of  pripr ......
    
    #......Start of  'mcmc' ......
      mcmc = lambda: None   ###Just get an empty object###
      mcmc.seed                      = mcmc_seed
      mcmc.samples                   = mcmc_samples
      mcmc.thinningFactor            = mcmc_thin
      mcmc.burnin                    = mcmc_burbin
      mcmc.chainNumber               = mcmc_chains
      mcmc.maxMoveStepSize           = 6
      mcmc.trendResamplingOrderProb  = 0.1000
      mcmc.seasonResamplingOrderProb = 0.1700
      mcmc.credIntervalAlphaLevel    = 0.950
    #......End of  mcmc ......
    
    #......Start of  'extra' ......    
      extra = lambda: None   ###Just get an empty object###
      extra.dumpInputData        = True
      extra.whichOutputDimIsTime = 1
      if 'ci' in kwargs:
            dump_ci  = kwargs['ci']      
      extra.computeCredible      = dump_ci
      extra.fastCIComputation    = True
      extra.computeSeasonOrder   = True
      extra.computeTrendOrder    = True
      extra.computeSeasonChngpt  = True
      extra.computeTrendChngpt   = True
      extra.computeSeasonAmp     = False
      extra.computeTrendSlope    = True
      extra.tallyPosNegSeasonJump= False
      extra.tallyPosNegTrendJump = False
      extra.tallyIncDecTrendJump = False
      extra.printProgress        = print_progress
      extra.printParameter       = print_param
      extra.printWarning         = print_warning  
      extra.quiet                = quiet
      extra.dumpMCMCSamples      = dump_mcmc
      extra.consoleWidth         = 85
      extra.numThreadsPerCPU     = 2
      extra.numParThreads        = 0
      if 'cputype' in kwargs.keys():
            extra.cputype = kwargs.get('cputype')		   
    #......End of  extra ......


      # For developers only, debug the C code on local machine if do_local_mode is not None
      do_local_mode =  kwargs.get('local')    
      
      if gui:
        o = cbeast.Rbeast('beastv4demo',Y, metadata, prior, mcmc, extra)
      elif do_local_mode is not None:
        # - Debug the C code on local machine if do_local_mode is not None
        # - y:/testold/Rbeast.pyd is the compiled C extennsion on my local machine
        import importlib
        spec   = importlib.util.spec_from_file_location('Rbeast', 'y:/testold/Rbeast.pyd')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)         
        class xxx:
           pass
        module.setClassObjects(xxx) 
        o = module.Rbeast('beast_'+method,Y, metadata, prior, mcmc, extra) 
      else:
        o = cbeast.Rbeast('beast_'+method,Y, metadata, prior, mcmc, extra) 
     
      return (o)

