function out = beast123(Y, metadata, prior, mcmc, extra, method)
%  
%   Run 'help beast123' to see the following
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  <strong>USAGE: out = beast123(Y, metadata, prior, mcmc, extra, method) </strong>
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   ----------------------------------------------------------------------------------------
%   <strong>***Y***</strong>:  input time series being regular or inrregular
%   ----------------------------------------------------------------------------------------
%
%        It could be a numeric vector, 2D matrix (e.g., multiple time series of the same length), or 
%        3D array time series (e.g., stacked  satellite images). The interface for beast123 is similar to the
%        documentation of its R version at <a href="matlab:web('https://rdrr.io/cran/Rbeast/man/beast123.html')">rdrr.io/cran/Rbeast/man/beast123.html</a>.
%
%   ----------------------------------------------------------------------------------------
%   <strong>***metadata***</strong>:  a struct variable specifying metadata for the input Y 
%   ----------------------------------------------------------------------------------------
%
%   metadata.isRegularOrdered :
%         Deprecated. Now use 'startTime' and 'deltaTime' to specify times for evenly-spaced/regulardata, 
%         and use 'time' to  provide times of individual data points for irregular or regular data
%   metadata.startTime        :
%        numeric (default to 1.0); the time of the 1st datapoint. It can be specified as a scalar
%        (e.g., 2021.0644), a vector of three values in the order of Year, Month, and  Day (e.g., [2021,1,24] )
%   metadata.deltaTime        : 
%        numeric (default to 1.0) or string; the time interval between consecutive data points. 
%        Its unit should be consistent with startTime or start. If start takes a numeric scalar, the unit
%		 is arbitrary and irrelevant to beast (e.g., 2021.3 can be of any unit: Year 2021.3, 2021.3
%        meters, 2021.3 degrees, ...). If start is a vector of [Year, Month, and Day], deltat has 
%        the unit of YEAR. For example, if start=[2021,1,2] for a monthly time series, start is 
%        converted to a fractional year 2021+(24-0.5)/365=2021.0644 and deltat=1/12 needs to be
%        set in order to specify the monthly interval. Alternatively, deltat can be provided as
%        a string to specify whether its unit is day, month, or year. Examples include '7 days',
%        '7d', '1/2 months', '1 mn', '1.0 year', and '1yr'. 
%   metadata.time             : 
%       Individual times of data points.  It can be a vector of numbers or date strings; it can also be a 
%       list of vectors of year, months, and days. Possible formats include:
%        (1) A vector of numerical values (e.g., [1984.23, 1984.27, 1984.36, ...] ). The unit of the
%           times is irrelevant to BEAST as long as it is consistent with the unit used for specifying 
%           startTime (start), deltaTime (deltat), and period.
%        (2) A vector of char strings. Examples are:
%           ["1984-03-27", "1984-04-10", "1984-05-12"] % a vector of strings
%           {'1984-03-27', '1984-04-10', '1984-05-12'} % a cell array of characters
%           {"1984/03/27", "1984,04,10", "1984 05 12"} (the delimiters differ as long as the YMD order is consistent)
%           {"LT4-1984-03-27", "LT4-1984-04-10", "LT4-1984+05,12"}
%           {"LT4-1984087ndvi", "LT4-1984101ndvi", "LT4-1984133ndvi"}
%           ["1984,,abc 3/ 27", "1984,,ddxfdd 4/ 10" "ggd1984,, 5/ ttt 12"]
%          BEAST uses several heuristics to automatically parse the date strings without a format specifier but
%          may fail due to ambiguity (e.g., in "LC8-2020-09-20-1984", no way to tell if 2020 or 1984 is the year). 
%          To ensure correctness, use a struct object as explained below to provide a date format specifier.
%        (3) A struct object time=struct(datestr=..., strfmat='...') consisting of a vector of date strings 
%          (time.datestr) and a format specifier (time.strFmt). The string time.strFmt specifies how to parse dateStr.
%           Three formats are currently supported:
%            (3a). All the date strings have a fixed pattern in terms of the relative positions of Year, Month, and Day.
%              For example, to extract 2001/12/02 etc from 
%              time.dateStr = {'P23R34-2001.1202333xd', 'O93X94-2002.1108133fd', 'TP3R34-2009.0122333td'} 
%              use time.strFmt='P23R34-yyyy.mmdd333xd' where yyyy, mm, and dd are the specifiers and other positions are 
%              wildcards and can be filled with any other letters different from yyyy, mm and dd.
%            (3b). All the date strings have a fixed pattern in terms of the relative positions of year and doy. 
%               For example, to extract 2001/045(day of year) from 'P23R342001888045', use strFmt='123123yyyy888doy'
%               where yyyy and doy are the specifiers and other positions are wildcards and can be filled with any other
%               letters different from yyyy, and doy. 'doy' must be three digit in length.
%            (3c). All the date strings have a fixed pattern in terms of the separation characters between year, month,
%               and day. For example, to extract 2002/12/02 from '2002,12/02', ' 2002 , 12/2', '2002,12 /02 ', use
%               strFmt='Y,M/D' where the whitespaces are ignored. To get 2002/12/02 from '2–12, 2012 ', use
%               strmFmt='D–M,Y'.
%        (4) A struct object of vectors to specify individual dates of the time series. Use time.year,time.month,and 
%          time.day to specify the dates; or alternatively use time.year and time.doy where each value of the doy vector is
%          a number within 1 and 365/366. Each vector must have the same length as the time dimension of Y.
%
%       <strong>NOTE about startTime (start), deltaTime (dT), and time:</strong>
%        BEAST currently handles only regular time series; irregular inputs will be first aggregated  before applying BEAST.
%       (a)For irregular time series, metadata.time is mandantory; start and dT should also be provided to specify 
%          the desired start and dT of the aggregated time series. If start and dT are missing, some best guesses 
%          will be used.
%       (b) For regular time series, the times can be specified by either (1) start and dT or (2) time.
%          If start, dT, and time are all missing, the default times are the indices 1:length(Y). 
%          If start, dT, and time are all provided, start and dT are ignored if the regular interval inferred from time
%          is the same as the supplied dT, but if the inferred dT is not the same as the supplied dT, the input regular 
%          time series will be re-sampled/aggregated to another regular time series at the supplied dT.%                         
%   metadata.whichDimIsTime   : 
%       which dim of the 2D or 3D input refer to time
%   metadata.season           : 
%        string (default to 'harmonic'); specify if y has a periodic component or not. Four values are possible.
%        (1)'none'    : y is trend-only; no periodic components are present in the time series. The args for 
%          the seasonal component (i.e.,sorder.minmax, scp.minmax and sseg.max) will be irrelevant and ignored.
%        (2)'harmonic': y has a periodic/seasonal component. The term season is a misnomer, being used here
%          to broadly refer to any periodic variations present in y. The periodicity is NOT a model parameter
%          estimated by BEAST but a known constant given by the user through freq. By default, the periodic 
%          component is modeled as a harmonic curve–a combination of sins and cosines.
%        (3)'dummy'  : the same as 'harmonic' except that the periodic/seasonal component is modeled as a 
%          non-parametric curve. The harmonic order arg sorder.minmax is irrelevant and is ignored.
%        (4)'svd'    : (experimental feature) the same as 'harmonic' except that the periodic/seasonal component
%          is modeled as a linear combination of function bases derived from a Single-value decomposition. The
%          SVD-based basis functions are more parsimonious than the harmonic sin/cos bases in parameterizing 
%          the seasonal variations; therefore, more subtle changepoints are likely to be detected.    
%   metadata.period           : 
%        numeric or string. Specify the period for the periodic/seasonal component in y. Needed 
%        only for data with a periodic/cyclic component (i.e., season='harmonic', 'dummy', or 'svd')
%        and not used for trend-only data (i.e., season='none'). The period of the cyclic component 
%		 should have a unit consisent with the unit of deltat. It holds that period=deltat*freq where
%        freq is the number of data samples per period. period or the number of data points per period
%        is not a BEAST model parameter and it has to be specified by the user. But if period is 
%        missing, BEAST first attempts to guess its value via auto-correlation before fitting the 
%        model. If period <= 0, season='none' is assumed, and the trend-only model is fitted without
%        a seasonal/cyclic component. If needed, use a string to specify whether the unit of period 
%        is day, month, or year. Examples are '1.0 year', '12 months', '365d', '366 days'.
%        Note: in earlier versions, 'freq' was used to specify the period and now deprecated in this 
%        version. The unit of 'period', if any, should be consistent with the unit of 'deltat' or deltaTime.  
%   metadata.missingValue     : 
%        a value indicating bad/missing data
%   metadata.maxMissingRate   : 
%        the max missingness rate beyond which the  ts will be ignored
%   metadata.deseasonalize    : 
%        boolean; if true, the input ts is first de-seasonalize by removing a global seasonalcomponent, prior 
%        to applying BEAST; after BEAST is done, the the global seasonal component will be added back
%        model.
%   metadata.detrend          : 
%        boolean; if true, the input ts is first de-trended by removing a global trend component, prior to
%        applying BEAST; after BEAST is done, the the global trend component will be added back
%   metadata.hasOutlierCmpnt  : (experimental feature)
%        boolean; if true, the model with an outlier component  will be fitted:
%        (a) Y = Trend + error            if season ='none'  and hasOutlier=false 
%		 (b) Y = Trend + Season  + error) if season~='none'  and hasOutlier=false 
%        (c) Y = Trend + Outlier + error  if season = 'none' and hasOutlier=true 
%		 (d) Y = Trend + Season  + Outlier + error if season~='none' and hasOutlier=true 
%
%   ----------------------------------------------------------------------------------------
%   <strong>***prior***</strong>:  prior hyperparameter for the BEAST model
%   ----------------------------------------------------------------------------------------
%  
%   prior.seasonMinKnotNum : 
%        an integer; min number of seasonal changepints (scp) allowed
%   prior.seasonMaxKnotNum : 
%       an integer;  max number of seasonal changepints (scp)allowed
%       NOTE: seasonMinKnotNum and seasonMaxKnotNum are used only if y has a seasonal component (i.e., season ~= 'none')
%        and ignored for tend-only data (i.e., season == 'none' ). If the min and max changepoint numbers are equal, 
%        BEAST assumes a constant number of scp and won't infer the posterior probability of the number of
%        changepoints, but it still estimates the occurrence probability of the changepoints over time (i.e.,
%        the most likely times at which these changepoints occur). If both the min and max scp numbers are set to 0,
%        no changepoints are allowed; then a global harmonic model is used to fit the seasonal component, but 
%       still, the most likely harmonic order will be inferred if seasonMinOrder is not equal to seasonMaxOrder.
%   prior.seasonMinOrder  : 
%        an integer; min harmonic order considered for the seasonal compnt
%   prior.seasonMaxOrder  : 
%        an integer; max harmonic order considered for the seasonal compnt
%        NOTE:  If the min and max orders are equal (seasonMinOrder==seasonMaxOrder) , BEAST assumes 
%        a constant harmonic order used and won't infer the posterior probability of harmonic orders.
%   prior.seasonMinSepDist : 
%       an integer; the min segment length of a seaosnal segment (i.e, the min seperation distance
%       between neighorboring seasonal changepoints)
%   prior.seasonLeftMargin : 
%        an integer (>=0); the number of leftmost data points excluded for seasonal changepoint detection.
%        That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the 
%        starting window/segment of length seasonLeftMargin. seasonLeftMargin must be an unitless integer 
%        (the number of time intervals/data points) so that the time window in the original unit is 
%        seasonLeftMargin*deltat. If missing, seasonLeftMargindefaults to seasonMinSepDist.
%   prior.seasonRightMargin : 
%        an integer (>=0); the number of rightmost data points excluded for seasonal changepoint detection.
%        That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the 
%        ending window/segment of length seasonRightMargin. seasonRightMargin must be an unitless integer 
%        (the number of time intervals/data points) so that the time window in the original unit is 
%        seasonRightMargin*deltat. If missing, seasonRightMargin to seasonMinSepDist.%
%   prior.trendMinKnotNum  :
%         an integer (>=0);min number of trend changepints allowed
%   prior.trendMaxKnotNum  : 
%         an integer (>=0); max number of trend changepints allowed
%         NOTE: If the min and max changepoint numbers are equal, BEAST assumes a constant number of changepoints and 
%           won't infer the posterior probability of the number of changepoints for the trend, but it still estimates
%           the occurrence probability of the changepoints over time (i.e., the most likely times at which these 
%           changepoints occur in the trend). If both the min and max numbers are set to 0, no changepoints are 
%           allowed; then a global polynomial trend is used to fit the trend component, but still, the most likely 
%           polynomial order will be inferred if trendMinOrder is not equal to trendMaxOrder.
%   prior.trendMinOrder	   : 
%         an integer (>=0); min polyonimal order considered for the trend component
%   prior.trendMaxOrder	   : 
%         an integer (>=0); max polyonimal order considered for the trend  component
%         NOTE:  The 0-th order corresponds to a constant term/a flat line and the 1st order is a line. 
%          If trendMinOrder==trendMaxOrder, BEAST assumes a constant polynomial order used and won't infer
%          the posterior probability of polynomial orders.
%   prior.trendMinSepDist  : 
%         an integer; the min segment length of a trend segment(i.e, the min seperation distance between
%         neighorboring trend changepoints)%
%   prior.trendLeftMargin  : 
%        an integer (>=0); the number of leftmost data points excluded for trend changepoint detection.
%        That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the 
%        starting window/segment of length trendLeftMargin. trendLeftMargin must be an unitless integer 
%        (the number of time intervals/data points) so that the time window in the original unit is 
%        trendLeftMargin*deltat. If missing, trendLeftMargin to trendMinSepDist.
%   prior.trendRightMargin : 
%        an integer (>=0); the number of rightmost data points excluded for trend changepoint detection.
%        That is, when fitting a piecewise polynomial trend model, no changepoints are allowed in the 
%        ending window/segment of length trendRightMargin. trendRightMargin must be an unitless integer 
%        (the number of time intervals/data points) so that the time window in the original unit is 
%        trendRightMargin*deltat. If missing, trendRightMargin to trendMinSepDist.%
%   prior.outlierMaxKnotNum : 
%        integer; needed only if metadata.hasOutlierCmpnt==true to specify the maximum number of outliers (i.e., 
%        outlier-type changepoints) allowed in the time series%
%   prior.precValue         : 
%        numeric (>0); the hyperparameter of the precision prior; the default value is 1.5.
%        precValue is useful only when precPriorType='constant', as further explained below
%   prior.precPriorType    : 
% 	      string. It takes one of 'constant', 'uniform', 'componentwise' (the default), and 'orderwise'. 
%         Below are the differences between them.
%         (1) 'constant': the precision parameter used to parameterize the model coefficients is fixed to
%            a constant specified by precValue. In other words, precValue is a user-defined hyperparameter 
%            and the fitting result may be VERY sensitive to the chosen values of precValue.
%         (2) 'uniform': the precision parameter used to parameterize the model coefficients is a random variable;
%            its initial value is specified by precValue. In other words, precValue will be inferred by the MCMC,
%            so the fitting result will be insensitive to the initial choice of precValue.
%         (3) 'componentwise': multiple precision parameters are used to parameterize the model coefficients for 
%            individual components (e.g., one for season and another for trend); their initial values is specified 
%            by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result will be 
%           insensitive to the initial choice in precValue.
%         (4) 'orderwise': multiple precision parameters are used to parameterize the model coefficients not just 
%            for individual components but also for individual orders of each component; their initial values is 
%            specified by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result
%            will be insensitive to the initial choice in precValue.
%
%   ----------------------------------------------------------------------------------------
%   <strong>***mcmc***</strong> :   parameters to set up the MCMC sampler
%   ----------------------------------------------------------------------------------------
%
%   mcmc.seed             : 
%        the seed for the random number generator            
%   mcmc.samples          : 
%        the number of samples to be collected
%   mcmc.thinningFactor   : 
%         chain thinning factor; take every 'thinningFactor'-th sample from the chain
%   mcmc.burnin           : 
%        number of initial samples discarded
%   mcmc.chainNumber      : 
%         number of chains
%   mcmc.maxMoveStepSize  : 
%         a moving step
%   mcmc.trendResamplingOrderProb  : 
%         probability to re-sample trend order
%   mcmc.seasonResamplingOrderProb : 
%         probability to re-sample seasonal order
%   mcmc.credIntervalAlphaLevel    : 
%        signifiance level for computing credible interval
%
%   ----------------------------------------------------------------------------------------
%   <strong>***extra***</strong> :   parameters to control computaiton or outputs
%   ----------------------------------------------------------------------------------------
%
%   extra.dumpInputData       :
%         boolean; if true, dump the input time series (o.data)
%   extra.whichOutputDimIsTime: 
%         an integer (>=1); which dim of the output array refer to time
%   extra.computeCredible     : 
%         boolean; if true, compute credible intervals. Unless needed for assessing the curve fitting, it is better 
%         to disable it for changepoint analysis because computing CI (e.g., o.trend.CI and o.season.CI)is expensive.
%   extra.fastCIComputation   : 
%        boolean; if true, a fast version of algorithm used to compute CI. But it is still expensive.
%        Set extra.computeCredible=false if CI (e.g., o.trend.CI and o.season.CI) is of litte interest
%   extra.computeSeasonOrder  : 
%        boolean;if true, estimate harmonic orders needed to adequately fit the seasonal component
%       (o.season.order)
%   extra.computeTrendOrder   : 
%        boolean; if true, estimate polynomial orders needed to adequately fit the trend component
%        (o.trend.order)
%   extra.computeSeasonChngpt  : 
%        boolean; if true, dump the scp detected (o.season.cp, o.season.cpCI, and cpAbruptChange)
%   extra.computeTrendChngpt   : 
%        boolean; if true, dump the tcp detected (o.trend.cp, o.trend.cpCI, and cpAbruptChange)
%   extra.computeSeasonAmp     : 
%        boolean; if true, estimate the seasonal ammplitue (o.seasin.amp)
%   extra.computeTrendSlope    : 
%        boolean; if true, estimate time-varying slope and the associated probabilities of slopes being positive or negative
%   extra.tallyPosNegSeasonJump:
%   extra.tallyPosNegTrendJump :
%   extra.tallyIncDecTrendJump :
%   extra.printProgressBar     : 
%        boolean; if true, show a progress bar
%   extra.printOptions         : 
%        boolean; if true, print the BEAST parameters at the start
%   extra.quite                : 
%        boolean; if true, print the BEAST parameters at the start
%   extra.consoleWidth         : 
%       an integer; the console/terminal width
%   extra.numThreadsPerCPU     : 
%        an integer; when processing multipplspecify the numbers of concurrent threads per cpu core            
%   extra.numParThreads        : 
%        an integer; specify the number of total concurrent  threads
%        <strong>Note</strong>: When handling many time series, BEAST can use multiple concurrent threads. 
%        extra.numParThreads specifies the total number of concurrent threads. If numParThreads=0, the actual
%        number of threads will be extra.numThreadsPerCPU * cpuCoreNumber; that is, each CPU core will generate
%        a number 'numThreadsPerCPU' of threads. On Windows 64, BEAST is group-aware and will affine or distribute
%        the threads to all the NUMA node. But currently, up to 256 CPU cores are supported. 
%
%        Mulithreading is IMPORTANT for remote sensing/earth science applicaitons.  For stacked images, do not use
%        an external parralell library (e.g. forpar) to call beast() in a loop pixel by pixel. Rather, use beast123(),
%        which supports faster multithreading internally.
%
%   <strong>Note</strong>:
%       beast, beast_irreg, and beast123 calls the same library(i.e., Rbeast.mex); they 
%       are just wrappers to the core BEAST algorithm.In particular, BEAST currently handles 
%       only regular time series; irregular inputs will be first aggregated  before applying BEAST.
%
%   ----------------------------------------------------------------------------------------
%   <strong>***method***</strong>: specify which method is used to formulate model posterior probability.
%   ----------------------------------------------------------------------------------------
%        a string (default to 'bayes'); four values are possible.
%        (1) 'bayes': the full Bayesian formulation as described in Zhao et al. (2019).
%        (2) 'bic'  : approximation of posterior probability using the Bayesian information criterion (bic).
%        (3) 'aic'  : approximation of posterior probability using the Akaike information criterion (aic).
%        (4) 'aicc' : approximation of posterior probability using the corrected Akaike information criterion (aicc).
%        (5) 'hic'  : approximation of posterior probability using the Hannan-Quinn information criterion (hic)
%
%
%   <strong>More help</strong>:  
%      The terse doc sucks (I know); so far, the best details are still the
%      R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
%      Matlab doesn't support keyword-style args, so Matlab's equivalent to R's beast(freq=1) 
%      beast('freq',1).
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   <strong>Note about beast(), beast_irreg, and beast123() </strong>:
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   beast, beast_irreg, and beast123 calls the same library(i.e., Rbeast.mex); they are just wrappers to 
%   the core BEAST algorithm.In particular, BEAST currently handles only regular time series; irregular 
%   inputs will be first aggregated  before applying BEAST.
% 
%   The keywords for beast() or beast_irreg() are converted to 'metadata', 'prior','mcmc', and 'extra' options
%   used in the beast123() interface. Some examples of the mapping are:
%            start               <->  metadata.startTime
%            deltat              <->  metadata.deltaTime
%            time                <->  metadata.time
%            season              <->  metadata.season
%            period              <->  metadata.period
%            deseasonalize       <->  metadata.deseasonalize            
%            detrend             <->  metadata.detrend  
%            hasOutlier          <->  metadata.hasOutlierCmpnt 
%            scp.minmax(1)       <->  prior.seasonMinKnotNum
%            scp.minmax(2)       <->  prior.seasonMaxOrder
%            sorder.minmax(1)    <->  prior.seasonMinOrder
%            sorder.minmax(2)    <->  prior.seasonMaxOrder
%            sseg.min            <->  prior.seasonMinSepDist
%	         sseg.leftmargin     <->  prior.seasonLeftMargin  
%	         sseg.rightmargin    <->  prior.seasonRightMargin  
%            tcp.minmax(1)       <->  prior.trendMinKnotNum
%            tcp.minmax(2)       <->  prior.trendMaxOrder
%            torder.minmax(1)    <->  prior.trendMinOrder
%            torder.minmax(2)    <->  prior.trendMaxOrder
%            tseg.min            <->  prior.trendMinSepDist
%	         tseg.leftmargin     <->  prior.trendLeftMargin  
%	         tseg.rightmargin    <->  prior.trendRightMargin
%	         precValue           <->  prior.precValue
%	         precPriorType       <->  prior.precPriorType
%	         ocp.max             <->  prior.outlierMaxKnotNum  
%            mcmc.seed           <->  mcmc.seed           
%            mcmc.samples        <->  mcmc.samples    
%            mcmc.thinningFactor <->  mcmc.thinningFactor 
%            mcmc.burnin         <->  mcmc.burnin 
%            mcmc.chainNumber    <->  mcmc.chainNumber  
%            ci                  <->  extra.computeCredible
%            print.progress      <->  extra.printProgressBar
%            print.options       <->  extra.printOptions 
%            quiet               <->  extra.quiet
%
%   <strong>Experts should use the the beast123 function.</strong>
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%   <strong>Examples</strong> :   
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       load('Nile')                % Nile river annual streamflow: trend-only data
%       metadata           = [];    % an emepty object to stuff new attribute fields 
%       metadata.season    = 'none' % trend-only
%       metadata.startTime = 1871;
%       metadata.deltaTime = 1;%       
%       o = beast123(Nile,metadata) % Default values will be used if parameters are missing 
%       printbeast(o)
%       plotbeast(o)
%
%
%       load('ohioNDVI')     % irregular Landsat NDVI time series
%       metadata            = [];          % an emepty object to stuff new attribute fields
%       metadata.season     = 'harmonic';  % seasonality is present
%       metadata.time       =  ohio.time;
%       metadata.deltaTime  = 1/12;        % aggregate into a monthly ts
%       metadata.period     = 1.0;         % period= 1/12 x 12 (freq)%     
%       o=beast123(ohio.ndvi,metadata)     % Default values will be used if parameters are missing
%       printbeast(o)
%       plotbeast(o, 'ncpStat','median')       
%       
%       load('ohioNDVI')     % irregular Landsat NDVI time series
%       ohio.datestr1        % strings of times for ohio.ndvi
%       metadata=[];          % an emepty object to stuff new attribute fields
%       metadata.time         =  [];
%       metadata.time.dateStr = ohio.datestr1;
%       metadata.time.strFmt  = '????yyyy?mm?dd';
%       metadata.deltaTime    = 1/12;   % aggregate into a monthly ts
%       metadata.period       = 1.0;    % period= 1/12 x 12 (freq)
%       % Default values will be used if parameters are missing
%       o = beast123(ohio.ndvi,metadata,[],[], struct('dumpInputData', true)) 
%       plotbeast(o, 'ncpStat','median')   
%
%       ohio.datestr2      % strings of times for ohio.ndvi
%       metadata              =[];
%       metadata.time         = [];
%       metadata.time.dateStr = ohio.datestr2;
%       metadata.time.strFmt  = 'LC8-yyyydoyxdvi';
%       metadata.deltaTime    = 1/12;   % aggregate into a monthly ts
%       metadata.period       = 1.0;    % period= 1/12 x 12 (freq)
%       o = beast123(ohio.ndvi,metadata) 
%       plotbeast(o, 'ncpStat','median')  
%       %See https://rdrr.io/cran/Rbeast/man/beast123.html for more details
%       %about the accepted formats of date strings
%
%       load('simData.mat') 
%       % A toy example of 3 time series: I decides to be lazy here bcz the 3
%       % time series are essentially the same. We just assume they are
%       % different for illustring the use of beast123 to handle 2d matrix      
%       simData.Y % a 774×3 matrix: 774 is the time dimesnion      
%       metadata                = [];      
%       metadata.whichDimIsTime = 1;      % 774 is the ts length 
%       o = beast123(simData.Y ,metadata) % default values used for missing parameters
%       printbeast(o,  1)                 % print the result for the first ts
%       plotbeast(o, 'index', 1)          % plot the result for the first ts
%       printbeast(o,  2)                 % print the result for the 2nd ts
%       plotbeast(o, 'index', 2)          % plot the result for the 2nd ts      
%
%      
%       Ytran   = simData.Y' % a 3x774 matrix: 774 is the time dimesnion      
%       metadata                = [];      
%       metadata.whichDimIsTime = 2 ; % 774 is the ts length 
%       o = beast123(Ytran,metadata)  % default values used for missing parameters
%       printbeast(o,  2)             % print the result for the 2nd ts
%       plotbeast(o, 'index', 2)      % plot the result for the 2nd ts   
%
%       load('imageStack.mat') 
%       % A toy example of stacked time series images: unevely-spaced in time
%       NDVI3D       = imageStack.ndvi    % a 12x9x1066 3D cube
%       NDVIdatestr  = imageStack.datestr % 1066 is the time series length%
%       metadata                =[];      
%       metadata.time           =[];
%       metadata.time.dateStr   = NDVIdatestr
%       metadata.time.strFmt    = 'LT05_018032_20110726.yyyy-mm-dd';
%       metadata.deltaTime      =  1/12;  % aggregated at a monthly interval
%       metadata.period         = 1.0;    % the period is 1.0 (year)
%       extra                   =  [];
%       extra.dumpInputData     = true;   % get a copy of the aggregated input
%       extra.numThreadsPerCPU   = 2;     % 2 threads per CPU core
%       o = beast123(NDVI3D,metadata,[],[], extra) 
%       imagesc(o.sig2)
%       imagesc(o.trend.ncpPr(:,:,1:3))
%       printbeast(o,[2,4]) %print the result at row 2 and col 4     
%       plotbeast(o,'index',[2,4]) %plot the result at row 2 and col 4     
%
%       NDVI_fyear  = imageStack.fyear  % the time in the unit of fractional year
%       metadata                =[];      
%       metadata.time           = NDVI_fyear
%       metadata.deltaTime      = 1/12;  % aggregated at a monthly interval
%       metadata.period         = 1.0;    % the period is 1.0 (year)
%       extra                   =  [];
%       extra.numThreadsPerCPU   = 2;     % 2 threads per CPU core
%       o = beast123(NDVI3D,metadata,[],[], extra) 
%
%   <strong>Contact info</strong>: To report bug or get help, do not hesitate to contact Kaiguang Zhao
%   at <strong>zhao.1423@osu.edu</strong>.
%
%   See also beast_irreg, beast, printbeast, plotbeast, extractbeast

%Check the second argument -- the option parameter
    if (nargin<5)     extra    =[];   end 
    if (nargin<4)     mcmc     =[];   end
    if (nargin<3)     prior    =[];   end
    if (nargin<2)     metadata =[];   end
	
	if (nargin<6)     method   ='bayes';  end
    
    %if nargin>=2 && (isstring(metadata) || ischar(metadata))
    %       seasonType   = char(metadata);
    %       metadata     = struct('season',seasonType);
    %end    
    
   if (nargin==0)          
       error("The input should not be empty; at least a time series needs to be provided.");
   end  
%%
    out=Rbeast( strcat('beast_',method), Y, metadata, prior, mcmc, extra);
end
