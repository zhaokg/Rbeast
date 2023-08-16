from .             import Rbeast as cb
from .cvt_to_numpy import force_convert_to_numpy

def beast(Y, \
          start          = 1,
          deltat         = 1,
          season         = 'harmonic',  # 'harmonic','dummy','svd','none'
          period         = float('nan'),
          scp_minmax     = [0, 10],
          sorder_minmax  = [0, 5],
          sseg_minlength = None,  # an integer
          tcp_minmax     = [0, 10],
          torder_minmax  = [0, 1],
          tseg_minlength = None,  # an integer
          detrend        = False,
          deseasonalize  = False,
          mcmc_seed      = 0,
          mcmc_burbin    = 200,
          mcmc_chains    = 3,
          mcmc_thin      = 5,
          mcmc_samples   = 8000,
          ci             = False,
          precValue      = 1.5,
          precPriorType  = 'componentwise',  # componentwise','uniform','constant','orderwise'
          print_options  = True,
          print_progress = True,
          quiet          = False,
          hasOutlier     = False,
          ocp_max        = 10,
          gui            = False,
          **kwargs
          ):
      """
################################################################################################
Bayesian changepoint detection and time series decomposition for a regular 1D signal or time series
    
Y = trend + error is fitted if the data has no periodic/seasonal variation (i.e., season='none')
Y = trend + seasonal + error is fitted if the data has  no periodic/seasonal variation 

 *Quick Examples*:
--------------------------------------------------------------------------------------------------
import Rbeast as rb
    
nile,yr =  rb.load_example('nile')             # annual flow of the Nile river
rb.beast( nile, start=1871, season='none' )
    
beach, year = rb.load_example('googletrend')         # the google search trend for 'beach'
o = rb.beast(beach, start= 2004, deltat=1/12, freq =12)
rb.plot(o) 
    
    
*Input arguments*:
--------------------------------------------------------------------------------------------------
Y:  a regular time series; it should be a numeric vector. For ireggular 
    time series, use 'beast_irreg' or 'beast123' instead. For multiple time 
    series or stacked time series images such as satellite data, use 'beast123'.
 
     ... :  the remaining arguments are many paired keywords and values to specifiy time information 
    or parameters for the beast algorithm. Check the R version of BEAST for detailed explanations
    (https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf). Below is a brief description.
 
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
         a number or string; the time interval between consecutive datapoints.. Use a string to specify the time unit
         (e.g., '1/12 year', '1.0 month', '30 days'). Possile units are year,
freq:  Deprecated. Replaced with 'period'. See below
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
         a string specifier. Possible values - 'none':  trend-only data with no 
         seasonality; 'harmonic': the seasonal/peridoic  component modelled via 
         harmonic curves; 'dummy': the seasonal component  modelled via a dummy 
         basis (i.e., pulse-like bases); 'svd': svd-derived  bases (experimental 
         feature)
scp_minmax: 
         a vector of two integers (e.g.,[0,5]); the min and max number of
         seasonal changepoints allowed
sorder_minmax: 
         a vector of two integers (e.g.,[1,3]); the min and max harmonic orders of
         seasonal changepoints (scp) allowed
sseg_min: 
         an integer; the min length of the segment for the seasonal component 
         i.e., the min distance between neighorbing changepoints)
tcp_minmax: 
         a vector of two integers (e.g.,[0,5]); the min and max number of
         trend changepoints (tcp) allowed
torder_minmax: 
         a vector of two integers (e.g.,[1,3]); the min and max orders of
         polynomials used to model the trend
tseg_min: 
         an integer; the min length of the segment for the trend component (i.e.,
         the min distance between neighorbing changepoints)%
deseasonalize: 
         boolean; if true, the input time series will be first
         de-seasonalized before applying beast by removing a global seasonal 
         component
detrend: 
         boolean; if true, the input time series will be first
         de-trend before applying beast by removing a global trend %
mcmc_seed: 
         a seed for the random number generator; set it to a non-zero
         integer to reproduce the results among different runs
mcmc_samples: 
         number of MCMC samples collected; the larger, the better
mcmc_thin: 
         a thinning factor for MCMC chains: take every 'mcmc.thin'-th sample
mcmc_burnin: 
         the number of initial samples of each chain to be discarded
mcmc_chains: 
         the number of MCMC chains%
print_progress: 
         boolean; if true, a progress bar is shown
print_options: 
         boolean; if true, print the beast paramers. 
quiet:
         boolean; if true, supress all the messages and printing
gui: 
        boolean; if true, show a gui to demostrate the MCMC sampling; runs only 
        on Windows not Linux or MacOS
 
The keywords for beast() are converted to 'metadata', 'prior','mcmc', and 'extra' options used 
         in the beast123() interface. Some examples are:
             deseasonalize <-> metadata.deseasonalize
             scp_minmax[0] <-> prior.seasonMinOrder
             scp_minmax[1] <-> prior.seasonMaxOrder
             sseg_min      <-> prior.seasonMinSepDist
             mcmc_seed     <-> mcmc.seed
             tcp_minmax[0] <-> prior.trendMinKnotNumber
             tcp_minmax[1] <-> prior.trendMaxKnotNumber
Experts should use the the beast123 function.
 
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
 
More help:  
--------------------------------------------------------------------------------------------------
       This terse help doc sucks (I know); so far, the best details are still the
       R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
       Python doesn't allow a '.' in variable names, so Python's equivalent to R's 
       beast(Y,start=1987,tcp.minmax=c(0,5)) is beast(Y, start=1987,  tcp_minmax=[0, 5]).
       
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
    
      Y = force_convert_to_numpy(Y)

    # ......Start of displaying 'MetaData' ......
      metadata = lambda: None  ###Just get an empty object###
                 
    #......Start of displaying 'MetaData' ......
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
    #........End of displaying MetaData ........

    #......Start of displaying 'prior' ......
      prior                   = lambda: None   ###Just get an empty object###
      prior.modelPriorType	  = 1
      if season !='none' or season == None:
            prior.seasonMinOrder   = sorder_minmax[0]
            prior.seasonMaxOrder   = sorder_minmax[1]
            prior.seasonMinKnotNum = scp_minmax[0]
            prior.seasonMaxKnotNum = scp_minmax[1]
            prior.seasonMinSepDist = sseg_minlength
      prior.trendMinOrder	  = torder_minmax[0]
      prior.trendMaxOrder	  = torder_minmax[1]
      prior.trendMinKnotNum  = tcp_minmax[0]
      prior.trendMaxKnotNum  = tcp_minmax[1]
      prior.trendMinSepDist  = tseg_minlength
      prior.K_MAX            = 500
      prior.precValue        = precValue
      prior.precPriorType    = precPriorType
    #......End of displaying pripr ......
    #......Start of displaying 'mcmc' ......
      mcmc = lambda: None   ###Just get an empty object###
      mcmc.seed                      =  mcmc_seed
      mcmc.samples                   = mcmc_samples
      mcmc.thinningFactor            = mcmc_thin
      mcmc.burnin                    = mcmc_burbin
      mcmc.chainNumber               = mcmc_chains
      mcmc.maxMoveStepSize           = 6
      mcmc.trendResamplingOrderProb  = 0.1000
      mcmc.seasonResamplingOrderProb = 0.1700
      mcmc.credIntervalAlphaLevel    = 0.950
    #......End of displaying mcmc ......
    #......Start of displaying 'extra' ......
      extra = lambda: None   ###Just get an empty object###
      extra.dumpInputData        = True
      extra.whichOutputDimIsTime = 1
      extra.computeCredible      = True
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
      extra.printProgressBar     = print_progress
      extra.printOptions         = print_options
      extra.quiet                 = quiet
      extra.consoleWidth         = 85
      extra.numThreadsPerCPU     = 2
      extra.numParThreads        = 0
    #......End of displaying extra ......
      if gui:
        o = cb.Rbeast('beastv4demo',Y, metadata, prior, mcmc, extra)
      else:
        #import importlib
        #spec   = importlib.util.spec_from_file_location('Rbeast', 'y:/testold/Rbeast.pyd')
        #module = importlib.util.module_from_spec(spec)
        #spec.loader.exec_module(module)         
        #class xxx:
        #   pass
        #module.setClassObjects(xxx)      
        o = cb.Rbeast('beastv4',Y, metadata, prior, mcmc, extra)      
      return (o)


