function out = beast(y, varargin)
%  
%  Run 'help beast' to see the following
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   <strong>USAGE: out = beast(y, varargin)</strong>
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   <strong>y</strong>: a regular time series as a numeric vector. For irregular time series, use 'beast_irreg'
%      or 'beast123' instead. For multiple time series or stacked time series images such as 
%      satellite data, use 'beast123'.
%
%   <strong>varargin</strong>:  a series of paired keywords and values to specifiy time information or BEAST 
%      parameters. The keywords mimic the R version of beast <a href="matlab:web('https://rdrr.io/cran/Rbeast/man/beast.html')">rdrr.io/cran/Rbeast/man/beast.html</a>. 
%      Unlike R, older Matlab versions don't support keyword-style arguments, so the beast parameters
%      should be provided in the following forms:
%
%     <strong>beast( Nile, 'season','none' )</strong>           % Equivalent to beast(Nile, season='none') in R/Python
%     <strong>beast( Nile, 'start', 1871, 'deltat', 1, 'season','none' )</strong>  
%     <strong>beast( Yellowstone, 'start', [1981,7,7], 'tcp.minmax', [0,10], 'deltat', 1/24 )</strong> 
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   <strong>*Possible Keywords*</strong>:
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      
%   <strong>start</strong>: 
%        numeric (default to 1.0); the time of the 1st datapoint of y. It can be specified as a scalar 
%       (e.g., 2021.0644), a vector of three values in the order of Year, Month, and Day (e.g., [2021,1,24])
%   <strong>deltat</strong>: 
%        numeric (default to 1.0) or string; the time interval between consecutive data points. Its unit
%        should be consistent with start. If start takes a numeric scalar, the unit is arbitrary and 
%        irrelevant to beast (e.g., 2021.3 can be of any unit: Year 2021.3, 2021.3 meters, 2021.3 degrees,
%        ...). If start is a vector of Year, Month, and Day, deltat has the unit of YEAR. For example, 
%        if start=[2021,1,2] for a monthly time series, start is converted to a fractional year 
%        2021+(24-0.5)/365=2021.0644 and deltat=1/12 needs to be set in order to specify the monthly interval. 
%        Alternatively, deltat can be provided as a string to specify whether its unit is day, month, or year. 
%        Examples include '7 days', '7d', '1/2 months', '1 mn', '1.0 year', and '1yr'. 
%   <strong>freq</strong>:  Deprecated. Replaced with 'period'. See below
%   <strong>period</strong>:  
%        numeric or string. Specify the period for the periodic/seasonal component in y. Needed only for 
%        data with a periodic/cyclic component (i.e., season='harmonic', 'dummy', or 'svd') and not used 
%        for trend-only data (i.e., season='none'). The period of the cyclic component should have a unit
%        consisent with the unit of deltat. It holds that period=deltat*freq where freq is the number of 
%        data samples per period. period or the number of data points per period is not a BEAST model parameter
%        and it has to be specified by the user. But if period is missing, BEAST first attempts to guess its 
%        value via auto-correlation before fitting the model. If period <= 0, season=='none' is assumed, and 
%        the trend-only model is fitted without a seasonal/cyclic component. If needed, use a string to specify
%        whether the unit of period is day, month, or year. Examples are '1.0 year', '12 months', '365d', '366 days'.
%        Note: in earlier versions, 'freq' was used to specify the period and now deprecated in this version. The
%        unit of 'period', if any, should be consistent with the unit of 'deltat'.
%   <strong>season</strong>: 
%        string (default to 'harmonic'); specify wheather y has a periodic component. Four values are possible.
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
%   <strong>scp.minmax</strong>: 
%        a vector of two integers (e.g.,[0,5]); the min and max numbers of seasonal changepoints (scp) allowed.
%        scp.minmax is used only if y has a seasonal component (i.e., season ~= 'none' ) and ignored for 
%        trend-only data (i.e., season == 'none' ). If the min and max changepoint numbers are equal, BEAST 
%        assumes a constant number of scp and won't infer the posterior probability of the number of changepoints,
%        but it still estimates the occurrence probability of the changepoints over time (i.e.,the most likely times
%        at which these changepoints occur). If both the min and max numbers are set to 0, no changepoints are allowed;
%        then a global harmonic model is used to fit the seasonal component, but still, the most likely harmonic order
%        will be inferred if sorder.minmax[1] is not equal to sorder.minmax[2]. 
%   <strong>sorder.minmax</strong>: 
%        a vector of two integers (e.g.,[1,3]); the min and max harmonic orders allowed to fit the seasonal component
%        If the min and max orders are equal (sorder.minmax[1]=sorder.minmax[2]), BEAST assumes a constant 
%        harmonic order used and won't infer the posterior probability of harmonic orders.
%   <strong>sseg.min</strong>: 
%        an integer; the min length of the segments in the seasonal component (i.e., the min distance between 
%        neighboring seasonal changepoints)
%   <strong>sseg.leftmargin</strong>: 
%        an integer (>=0); the number of leftmost data points excluded for seasonal changepoint detection.
%        That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the starting 
%        window/segment of length sseg.leftmargin. sseg.leftmargin must be an unitless integer (the number of time
%        intervals/data points) so that the time window in the original unit is  sseg.leftmargin*deltat. 
%        If missing, sseg.leftmargin defaults to sseg.min.
%   <strong>sseg.rightmargin</strong>: 
%        an integer (>=0); the number of rightmost data points excluded for seasonal changepoint detection.
%        That is, when fitting a piecewise harmonic seasonal model, no changepoints are allowed in the ending 
%        window/segment of length sseg.rightmargin. sseg.rightmargin must be an unitless integer (the number of time 
%        intervals/data points) so that the time window in the original unit is sseg.rightmargin*deltat. 
%        If missing, sseg.rightmargin defaults to sseg.min.
%   <strong>tcp.minmax</strong>: 
%        a vector of two integers (e.g.,[0,5]); the min and max number of trend changepoints (tcp) allowed.
%       If the min and max changepoint numbers are equal, BEAST assumes a constant number of changepoints and 
%       won't infer the posterior probability of the number of changepoints for the trend, but it still estimates
%       the occurrence probability of the changepoints over time (i.e., the most likely times at which these 
%       changepoints occur in the trend). If both the min and max numbers are set to 0, no changepoints are 
%       allowed; then a global polynomial trend is used to fit the trend component, but still, the most likely 
%      polynomial order will be inferred if torder.minmax[1] is not equal to torder.minmax[2].
%   <strong>torder.minmax</strong>: 
%        a vector of two integers (e.g.,[1,3]); the min and max orders of polynomials used to model the trend.
%        The 0-th order corresponds to a constant term/a flat line and the 1st order is a line. 
%        If torder.minmax[1]=torder.minmax[2], BEAST assumes a constant polynomial order used and won't infer
%        the posterior probability of polynomial orders.
%   <strong>tseg.min</strong>: 
%        an integer; the min length of the segment for the trend component (i.e., the min distance between 
%        neighboring changepoints)
%   <strong>tseg.leftmargin</strong>: 
%        an integer (>=0); the number of leftmost data points excluded for trend changepoint detection.
%        That is, when fitting a piecewise trend model, no changepoints are allowed in the starting 
%        window/segment of length tseg.leftmargin. tseg.leftmargin must be an unitless integer (the number of
%        time intervals/data points) so that the time window in the original unit is tseg.leftmargin*deltat. 
%        If missing, tseg.leftmargin defaults to tseg.min.
%   <strong>tseg.rightmargin</strong>: 
%        an integer (>=0); the number of rightmost data points excluded for trend changepoint detection.
%        That is, when fitting a piecewise trend model, no changepoints are allowed in the ending 
%        window/segment of length tseg.rightmargin. tseg.rightmargin must be an unitless integer (the number of 
%        time intervals/data points) so that the time window in the original unit is tseg.rightmargin*deltat. 
%        If missing, tseg.rightmargin defaults to tseg.min.
%   <strong>method</strong>: 
%        a string (default to 'bayes'); specify which method is used to model the posterior probability.
%        (1) 'bayes': the full Bayesian formulation as described in Zhao et al. (2019).
%        (2) 'bic'  : approximation of posterior probability using the Bayesian information criterion (bic).
%        (3) 'aic'  : approximation of posterior probability using the Akaike information criterion (aic).
%        (4) 'aicc' : approximation of posterior probability using the corrected Akaike information criterion (aicc).
%        (5) 'hic'  : approximation of posterior probability using the Hannan-Quinn information criterion (hic)
%        (6) 'bic0.25':  approximation using the Bayesian information criterion adopted from Kim et al. (2016) <doi: 
%             10.1016/j.jspi.2015.09.008>; bic0.25=n*ln(SSE)+0.25k*ln(n) with less complexity penelaty than the standard BIC.
%        (7) 'bic0.50': the same as above except that the penalty factor is 0.50.
%        (8) 'bic1.5':  the same as above except that the penalty factor is 1.5.
%        (9) 'bic2':    the same as above except that the penalty factor is 2.0.
%   <strong>deseasonalize</strong>: 
%       boolean; if true, the input ts is first de-seasonalize by removing a global seasonalcomponent, prior 
%       to applying BEAST; after BEAST is done, the the global seasonal component will be added back.
%   <strong>detrend</strong>: 
%       boolean; if true, the input ts is first de-trended by removing a global trend component, prior to
%       applying BEAST; after BEAST is done, the the global trend component will be added back.
%   <strong>mcmc.seed</strong>: 
%        a seed for the random number generator; set it to a non-zero integer to reproduce the results 
%        among different runs
%   <strong>mcmc.samples</strong>: 
%         an integer; number of MCMC samples collected; the larger, the better
%   <strong>mcmc.thin</strong>: 
%        a thinning factor for MCMC chains: take every 'mcmc.thin'-th sample
%   <strong>mcmc.burnin</strong>: 
%        the number of initial samples of each chain to be discarded
%   <strong>mcmc.chains</strong>: 
%        the number of MCMC chains
%   <strong>precValue</strong>: 
%        numeric (>0); the hyperparameter of the precision prior; the default value is 1.5.
%        precValue is useful only when precPriorType='constant', as further explained below
%   <strong>precPriorType</strong>: 
% 	      a string. It takes one of 'constant', 'uniform', 'componentwise' (the default), and 'orderwise'. 
%         (1) 'constant'     : the precision parameter used to parameterize the model coefficients is fixed to
%            a constant specified by precValue. In other words, precValue is a user-defined hyperparameter 
%            and the fitting result may be VERY sensitive to the chosen values of precValue.
%         (2) 'uniform'      : the precision parameter used to parameterize the model coefficients is a random 
%            variable; its initial value is specified by precValue. In other words, precValue will be inferred by 
%            the MCMC, so the fitting result will be insensitive to the initial choice of precValue.
%         (3) 'componentwise': multiple precision parameters are used to parameterize the model coefficients for 
%            individual components (e.g., one for season and another for trend); their initial values is specified 
%            by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result will be 
%           insensitive to the initial choice in precValue.
%         (4) 'orderwise'    : multiple precision parameters are used to parameterize the model coefficients not just 
%            for individual components but also for individual orders of each component; their initial values is 
%            specified by precValue. In other words, precValue will be inferred by the MCMC, so the fitting result
%            will be insensitive to the initial choice in precValue.
%   <strong>hasOutlier</strong>: 
%        boolean; if true, the model with an outlier component  will be fitted:
%        (a) Y = Trend + error                     if season = 'none'  and hasOutlier=false 
%		 (b) Y = Trend + Season  + error)          if season~= 'none'  and hasOutlier=false 
%        (c) Y = Trend + Outlier + error           if season = 'none'  and hasOutlier=true 
%		 (d) Y = Trend + Season  + Outlier + error if season~= 'none'  and hasOutlier=true 
%       where the outlier component is extreme spikes or dips at isolated points of time.
%   <strong>ocp.minmax</strong>: 
%        a vector of 2 integers (>=0); the min and max numbers of outlier-type changepoints (ocp) allowed in the time series
%        trend component. Ocp refers to spikes or dips at isolated times that can't be modeled as trends or seasonal terms.
%   <strong>print.progress</strong>: 
%        boolean; if true, a progress bar is shown
%   <strong>print.param</strong>: 
%        boolean; if true, print the BEAST paramers. 
%   <strong>print.warning</strong>: 
%        boolean; if true, print warning messages. 
%   <strong>quiet</strong>: 
%        boolean. If TRUE, warning messages are suppressed and not printed.
%   <strong>dump.ci</strong>: 
%         boolean; if true, compute credible intervals. Unless needed for assessing the curve fitting, it is better 
%         to disable it for changepoint analysis because computing CI (e.g., o.trend.CI and o.season.CI)is expensive.
%   <strong>dump.mcmc</strong>: 
%        boolean. If TRUE, dump the sampled models in the MCMC chains.
%   <strong>gui</strong>: 
%       boolean; if true, show a gui to demostrate the MCMC sampling; runs only on Windows not Linux or MacOS
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   <strong>Note about beast(), beast_irreg, and beast123() </strong>:
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%            sseg.leftmargin     <->  prior.seasonLeftMargin  
%            seg.rightmargin    <->  prior.seasonRightMargin  
%            tcp.minmax(1)       <->  prior.trendMinKnotNum
%            tcp.minmax(2)       <->  prior.trendMaxOrder
%            torder.minmax(1)    <->  prior.trendMinOrder
%            torder.minmax(2)    <->  prior.trendMaxOrder
%            tseg.min            <->  prior.trendMinSepDist
%            tseg.leftmargin     <->  prior.trendLeftMargin  
%            tseg.rightmargin    <->  prior.trendRightMargin
%            precValue           <->  prior.precValue
%            precPriorType       <->  prior.precPriorType
%            ocp.minmax(2)       <->  prior.outlierMaxKnotNum  
%            mcmc.seed           <->  mcmc.seed           
%            mcmc.samples        <->  mcmc.samples    
%            mcmc.thinningFactor <->  mcmc.thinningFactor 
%            mcmc.burnin         <->  mcmc.burnin 
%            mcmc.chainNumber    <->  mcmc.chainNumber  
%            dump.ci             <->  extra.computeCredible
%            print.progress      <->  extra.printProgress
%            print.param         <->  extra.printParameter
%            quiet               <->  extra.quiet
%            dump.mcmc           <->  extra.dumpMCMCSamples
%
%   <strong>Experts should use the the beast123 function.</strong>
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  <strong>*Result/Output*</strong>: The output is a struct variable; example of the fields include
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    marg_lik: marginal likilood; the larger, the better
%    sig2    : variance  of error
%    trend   : the trend component; a struct variable (say, T)
%    season  : the season componet; a stuct variable  (say,S)
%    outlier : the outlier componet; a stuct variable  (say,O)
%    The subfields of trend or season:
%      .ncpPr        : the prob distribution for number of changepoints
%      .ncp          : mean number of changepoints in trend or seasonality
%      .ncp_meidan   : median number of changepoints
%      .ncp_mode     : mode from ncpPr
%      .ncp_pct90    : 90% percentile from ncpPr
%      .cpOccPr      : changepoint occurrance probability over time
%      .cp           : list of all possible changepoints (many are not signficant)
%      .cpPr         : occurrence probability of the changepoints in cp
%      .cpAbruptChange: the sudden changes in trend or seasonlity at cp
%      .cpCI         : confidence interval of the cps
%      .Y            : the fitted trend or seasonality 
%      .SD           : standard deviation of the fitted Y
%      .CI           : Credible interval of the fittted Y
%      .order        : the mean harmonic or polynomial orders estimated to fit the seasonal and trend     
%     out.trend.slp         : slope of the trend 
%     out.trend.slpSD       : standard dev of the estimated slope
%     out.trend.slpSgnPosPr : time-varying probability of the slope being postive
%     out.trend.slpSgnZeroPr: time-varying probability of the slope being 0
%     out.season.amp        : amplitue of the estiamted seasonality overtime
%     out.season.ampSD      : standard ev of the estiamated amplitude
%
%   <strong>More help</strong>:  
%      This terse help doc sucks (I know); so far, the best details are still the
%      R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
%      Matlab doesn't support keyword-style args, so Matlab's equivalent to R's 
%      beast(Y,<strong>start</strong>=1987,<strong>freq</strong>=1) is beast(Y,<strong>'start'</strong>, 1987, <strong>'freq'</strong>,1).
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%   <strong>Examples</strong>:
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       % Nile river annual streamflow: trend-only data
%       load('Nile.mat')              
%       o = beast(Nile, 'start', 1871, 'season','none') 
%       printbeast(o)
%       plotbeast(o)
%       
%       % Explicitly specify deltat=1. BEAST knows nothing about the unit
%       % of 1871 and 1.0 (i.e., 1871 years, 1871 seconds, or 1871 meters?) 
%       o = beast(Nile, 'start', 1871, 'deltat', 1.0, 'season','none') 
%
%       % start is given a date 1871-1 (Year-Mon). The time unit is
%       % then fractional/decimal year. delta=1.0 means 1.0 year
%       o = beast(Nile, 'start', [1871,1], 'deltat', 1.0, 'season','none') 
%
%       % period=0 means a trend-only signal, which is equivalent to season='none'
%       o = beast(Nile, 'start', 1871, 'deltat', 1, 'period', 0) 
%
%       % Use a string to specify a unit for deltat or period (e.g., deltat='1 year')
%       % The time unit is also fractional year. 1871 means Year 1871
%       o = beast(Nile, 'start', 1871, 'deltat', '1 year', 'period', 0) 
%
%       % Use a string to specify a unit for delta or period (e.g., deltat='12 mo')
%       % The time unit is fractional year. 1871 means Year 1871
%       o = beast(Nile, 'start', 1871, 'deltat', '12 mo', 'period', 0) 
%
%       % Fit the data with an extra outlier component:l Y=trend+outlier+error 
%       o = beast(Nile, 'start', 1871, 'deltat',1.0,'season','none','hasOutlier',true)
%       plotbeast(o)
%		
%       % Do not print the options 
%       o=beast(Nile, 'start', 1871, 'deltat',1.0,'season','none','print.options',false)
%
%       % Show a gui window to demostrate the BEAST sampling process in
%       % real-time (for Windows only not Linux and MacOS)% 
%       beast(Nile,'season','none', 'gui',true) 
%
%       %% Monthly google trend of the search word 'beach'
%       load('googletrend.mat')   
%       o = beast(beach, 'start', [2004,1],'deltat', 1/12) %deltat = 1/12 yr =1 month
%       printbeast(o)
%       plotbeast(o)
%       plotbeast(o,'ncpStat','median')
%
%       % delta  = 1/12: for dates, the default unit is year, so delta=1/12yr=1 month;       
%       % period = 1.0 means 1 year
%       o = beast(beach, 'start', [2004,1],'deltat', 1/12, 'period',1.0)  
%
%       % period='12 month': use a string to explicitly specify the unit              
%       o = beast(beach, 'start', [2004,1],'deltat', 1/12, 'period','12 month')
%       o = beast(beach, 'start', [2004,1],'deltat', '1 month', 'period','365 days')
%
%       %% Monthly air co2 data since 1959: deltaTime=1/12 year
%       load('co2.mat')     
%       o = beast(co2, 'start', [1959,1,15], 'deltat', 1/12, 'period',1.0)
%       printbeast(o)
%       plotbeast(o)
%       plotbeast(o,'ncpStat','median')
%
%      %% Daily covid-19 infection statistics 
%       load('covid19.mat')    
%       Y    = sqrt(double(covid19.newcases));
%       Date = covid19.datestr;
%       % the min length of seasonal segments is set to 30 data points
%       o    = beast(Y, 'start',[2020,01,22], 'deltat', 1/365, 'period', 7/365,'sseg.min',30)
%       printbeast(o)
%       plotbeast(o)
%       plotbeast(o,'ncpStat','median')
%
%       % Use a string to specify delta with a unit
%       o = beast(Y, 'start',[2020,01,22], 'deltat', '1.0 day',  'period', '7days','sseg.min',30)
%
%       % Convert and aggregate the daily data into a weekly time series (i.e., deltaT=7 days)
%       % then fit a trend-only model with no periodic component
%       o = beast(Y, 'time',Date, 'deltat', '7days',  'period', 'none')
%       %o= beast(Y, 'time',Date, 'deltat', '7days',  'period', 0 ) % equivalent to period='none'
%       plotbeast(o)
%
%   <strong>Contact info</strong>: To report bug or get help, do not hesitate to contact Kaiguang Zhao
%   at <strong>zhao.1423@osu.edu</strong>.
%
%   See also beast123, beast_irreg, printbeast, plotbeast, extractbeast


%% Check the second argument -- the option parameter
   n = length(varargin);
   if mod(n,2) ~= 0
        msg=[ "the optional arg list must be paired keywords and values; the number of extra args must be even.\n", ...
              " Examples:  beast(y) \n", ...
              "            beast(y,'season','none','start',1980) \n", ...
             ];
        msg = sprintf(strcat(msg{:}));
        error(msg);
   end
    
   KeyList   = varargin(1:2:n);
   KeyList   = cellfun(@char,KeyList,'UniformOutput',false); % convert strings--if any--to chars,
   KeyList   = lower(KeyList);
   ValList   = varargin(2:2:n);
      
%% get values from keys. The last arg is the default value if the key is missing from varagin/KeyList
  
   start           = GetValueByKey(KeyList, ValList, 'start',  []);
   deltat          = GetValueByKey(KeyList, ValList, 'deltat', []);
   time            = GetValueByKey(KeyList, ValList, 'time',   []);    
   period          = GetValueByKey(KeyList, ValList, 'period',  []); 
   nsamples_per_period  = GetValueByKey(KeyList, ValList, 'freq',  []); 
    
   season          = GetValueByKey(KeyList, ValList, 'season',        'harmonic'); 
   sorder_minmax   = GetValueByKey(KeyList, ValList, 'sorder.minmax', [1,5]); 
   scp_minmax      = GetValueByKey(KeyList, ValList, 'scp.minmax',    [0,10]); 
   sseg_min        = GetValueByKey(KeyList, ValList, 'sseg.min',      []); 
   sseg_leftmargin = GetValueByKey(KeyList, ValList, 'sseg.leftmargin',  []); 
   sseg_rightmargin= GetValueByKey(KeyList, ValList, 'sseg.rightmargin', []); 
   
   deseasonalize   = GetValueByKey(KeyList, ValList, 'deseasonalize', false); 
   detrend         = GetValueByKey(KeyList, ValList, 'detrend', false); 
   
   torder_minmax   = GetValueByKey(KeyList, ValList, 'torder.minmax', [0,1]); 
   tcp_minmax      = GetValueByKey(KeyList, ValList, 'tcp.minmax',    [0,10]); 
   tseg_min        = GetValueByKey(KeyList, ValList, 'tseg.min',      []);
   tseg_leftmargin = GetValueByKey(KeyList, ValList, 'tseg.leftmargin',  []); 
   tseg_rightmargin= GetValueByKey(KeyList, ValList, 'tseg.rightmargin', []); 

   precValue       = GetValueByKey(KeyList, ValList, 'precValue',       1.5); 
   precPriorType   = GetValueByKey(KeyList, ValList, 'precPriorType',   'componentwise');    
   
   hasOutlierCmpnt = GetValueByKey(KeyList, ValList, 'hasOutlier',        []); 
   ocp_minmax      = GetValueByKey(KeyList, ValList, 'ocp.minmax',        [0,10]); 
      
   mcmc_seed       = GetValueByKey(KeyList, ValList, 'mcmc.seed',     0);         
   mcmc_samples    = GetValueByKey(KeyList, ValList, 'mcmc.samples',  8000);
   mcmc_thin       = GetValueByKey(KeyList, ValList, 'mcmc.thin',     5); 
   mcmc_burnin     = GetValueByKey(KeyList, ValList, 'mcmc.burnin',   200);
   mcmc_chainNumber= GetValueByKey(KeyList, ValList, 'mcmc.chains',   3);  

   printProgress   = GetValueByKey(KeyList, ValList, 'print.progress', true);     
   printParameter   = GetValueByKey(KeyList, ValList, 'print.param',  true);    
   printWarning     = GetValueByKey(KeyList, ValList, 'print.warning',  true);
   quiet            = GetValueByKey(KeyList, ValList, 'quiet',          false);   
   gui              = GetValueByKey(KeyList, ValList, 'gui',            false); 
   ci               = GetValueByKey(KeyList, ValList, 'dump.ci',        false);      
   mcmc_dump        = GetValueByKey(KeyList, ValList, 'dump.mcmc',   false);     
   
   methods          = GetValueByKey(KeyList, ValList, 'method',        'bayes'); 
   

%% Convert the opt parameters to the individual option parameters (e.g., metadata, prior, mcmc, and extra)

   %......Start of displaying 'MetaData' ......
   metadata = [];
   metadata.isRegularOrdered = true;
   metadata.season           = season;
   metadata.time             = time;
   metadata.startTime        = start;
   metadata.deltaTime        = deltat;
   if isempty(period) && ~isempty(deltat) && ~isempty(nsamples_per_period) && ~strcmp(season, 'none')
       period=nsamples_per_period*deltat;
   end   
   metadata.period           = period;

 
   if strcmp(metadata.season, 'svd')
      % if isempty(freq)|| freq <= 1.1 || isnan(freq)
      %     error("When season=svd, freq must be specified and larger than 1.");
      % end
      % metadata.svdTerms = svdbasis(y, freq, deseasonalize);
   end
   metadata.missingValue     = NaN;
   metadata.maxMissingRate   = 0.75;
   metadata.deseasonalize    = deseasonalize;
   metadata.detrend          = detrend;
   metadata.hasOutlierCmpnt  = hasOutlierCmpnt;
%........End of displaying MetaData ........

%......Start of displaying 'prior' ......
   prior = [];
   prior.modelPriorType	  = 1;
   if ~strcmp(metadata.season, 'none')              
       prior.seasonMinOrder    = sorder_minmax(1);
       prior.seasonMaxOrder    = sorder_minmax(2);
       prior.seasonMinKnotNum  = scp_minmax(1);
       prior.seasonMaxKnotNum  = scp_minmax(2);   
       prior.seasonMinSepDist  = sseg_min;
       prior.seasonLeftMargin  = sseg_leftmargin;
       prior.seasonRightMargin = sseg_rightmargin;
   end   
   prior.trendMinOrder	  = torder_minmax(1);
   prior.trendMaxOrder	  = torder_minmax(2);
   prior.trendMinKnotNum  = tcp_minmax(1);
   prior.trendMaxKnotNum  = tcp_minmax(2);
   prior.trendMinSepDist  = tseg_min;
   prior.trendLeftMargin  = tseg_leftmargin;
   prior.trendRightMargin = tseg_rightmargin;

   if hasOutlierCmpnt
      prior.outlierMinKnotNum = ocp_minmax(1);   
      prior.outlierMaxKnotNum = ocp_minmax(2);
   end
        
   prior.precValue        = precValue;
   prior.precPriorType    = precPriorType;
%......End of displaying pripr ......

%......Start of displaying 'mcmc' ......
   mcmc = [];
   mcmc.seed                      = mcmc_seed;
   mcmc.samples                   = mcmc_samples;
   mcmc.thinningFactor            = mcmc_thin;
   mcmc.burnin                    = mcmc_burnin;
   mcmc.chainNumber               = mcmc_chainNumber;
   
   %mcmc.maxMoveStepSize           = 28
   mcmc.trendResamplingOrderProb  = 0.1000;
   mcmc.seasonResamplingOrderProb = 0.1700;
   mcmc.credIntervalAlphaLevel    = 0.950;
%......End of displaying mcmc ......

%......Start of displaying 'extra' ......
   extra = [];
   extra.dumpInputData        = true;
   extra.whichOutputDimIsTime = 1;
   extra.computeCredible      = ci;
   extra.fastCIComputation    = true;
   extra.computeSeasonOrder   = true;
   extra.computeTrendOrder    = true;
   extra.computeSeasonChngpt  = true;
   extra.computeTrendChngpt   = true;
   extra.computeSeasonAmp     = ~strcmp(metadata.season, 'svd');
   extra.computeTrendSlope    = true;
   extra.tallyPosNegSeasonJump= false;
   extra.tallyPosNegTrendJump = false;
   extra.tallyIncDecTrendJump = false;
   extra.printProgress        = printProgress;
   extra.printParameter       = printParameter;
   extra.printWarning         = printWarning;
   extra.quiet                = quiet;
   extra.consoleWidth         = 70;
   extra.numThreadsPerCPU     = 2;
   extra.numParThreads        = 0;
   extra.dumpMCMCSamples      = mcmc_dump;
%......End of displaying extra ......


 if (gui)
    out=Rbeast(' beastv4demo',            y, metadata, prior, mcmc, extra);
 else
    out=Rbeast( strcat('beast_',methods), y, metadata, prior, mcmc, extra);
 end
 
end


%% Functions to return a default value if the field is missing from opt
function value=GetValueByKey(KeyList, ValList, key,defaultValue)
   idx = find(strcmp(KeyList,lower(key)));
   if isempty(idx)
       value = defaultValue;
   else
       value = ValList{idx(1)};
   end
end
 