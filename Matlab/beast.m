function out = beast(y, varargin)
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run 'help beast' to see the following
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   USAGE: out=<strong>beast(y, ...) </strong>
%
%   <strong>y</strong>:  a regular time series; it should be a numeric vector. For ireggular 
%   time series, use 'beast_irreg' or 'beast123' instead. For multiple time 
%   series or stacked time series images such as satellite data, use 'beast123'.
%
%   <strong> ... </strong>:  a series of paired keywords and values to  specifiy time information 
%   or parameters for the BEAST algorithm. The keywords mimic the R version of beast
%   <a href="matlab:
%   web(https://rdrr.io/cran/Rbeast/man/beast.html')">rdrr.io/cran/Rbeast/man/beast.html</a>. Unlike R, Matlab doesn't support keyword-style 
%   arguments, so the beast parameters should be provided in the following forms:
%
%   <strong>beast( Nile, 'start', 1871, 'deltat', 1, 'season','none' )</strong>
%   <strong>beast( Yellowstone, 'start', [1981,7,7], 'tcp.minmax', [0,10], 'deltat', 1/24 )</strong> 
%
%   <strong>*Possible Keywords*</strong>:
%      
%   <strong>start</strong>: 
%        the start time of the regular time series
%   <strong>deltat</strong>: 
%        the time interval between consecutive datapoints (e.g., 1/12
%        for monthly time series if the time unit is year).
%   <strong>freq</strong>:  Deprecated. Replaced with 'period'. See below
%   <strong>period</strong>:  
%        a numeric value to specify the period if peridodic/seasonal variations 
%        are present in the data. If period is given a zero, negative value or 'none' 
%        it suggests no seasonal/periodic component in the signal. (season='none'
%        also suggests no periodic component).
%        Note: in earlier versions, 'freq' was used to specify the period and
%        now deprecated in this version. The unit of 'period', if any
%        should be consistent with the unit of 'deltat'..
%   <strong>season</strong>: 
%        a string specifier. Possible values - 'none':  trend-only data with no 
%        seasonality; 'harmonic': the seasonal/peridoic  component modelled via 
%        harmonic curves; 'dummy': the seasonal component  modelled via a dummy 
%        basis (i.e., pulse-like bases); 'svd': svd-derived  bases (experimental 
%        feature)
%   <strong>scp.minmax</strong>: 
%        a vector of two integers (e.g.,[0,5]); the min and max number of
%        seasonal changepoints allowed
%   <strong>sorder.minmax</strong>: 
%        a vector of two integers (e.g.,[1,3]); the min and max harmonic orders of
%        seasonal changepoints (scp) allowed
%   <strong>sseg.min</strong>: 
%        an integer; the min length of the segment for the seasonal component 
%        i.e., the min distance between neighorbing changepoints)
%   <strong>tcp.minmax</strong>: 
%        a vector of two integers (e.g.,[0,5]); the min and max number of
%        trend changepoints (tcp) allowed
%   <strong>torder.minmax</strong>: 
%        a vector of two integers (e.g.,[1,3]); the min and max orders of
%        polynomials used to model the trend
%   <strong>tseg.min</strong>: 
%        an integer; the min length of the segment for the trend component (i.e.,
%        the min distance between neighorbing changepoints)%
%   <strong>deseasonalize</strong>: 
%        boolean; if true, the input time series will be first
%        de-seasonalized before applying BEAST by removing a global seasonal 
%        component
%   <strong>detrend</strong>: 
%        boolean; if true, the input time series will be first
%        de-trend before applying BEAST by removing a global trend %
%   <strong>mcmc.seed</strong>: 
%        a seed for the random number generator; set it to a non-zero
%        integer to reproduce the results among different runs
%   <strong>mcmc.samples</strong>: 
%        number of MCMC samples collected; the larger, the better
%   <strong>mcmc.thin</strong>: 
%        a thinning factor for MCMC chains: take every 'mcmc.thin'-th sample
%   <strong>mcmc.burnin</strong>: 
%        the number of initial samples of each chain to be discarded
%   <strong>mcmc.chains</strong>: 
%        the number of MCMC chains%
%   <strong>print.progress</strong>: 
%        boolean; if true, a progress bar is shown
%   <strong>print.options</strong>: 
%        boolean; if true, print the BEAST paramers. 
%   <strong>gui</strong>: 
%       boolean; if true, show a gui to demostrate the MCMC sampling; runs only 
%       on Windows not Linux or MacOS
%
%   The keywords for beast() are converted to 'metadata', 'prior','mcmc', and 'extra' options used 
%        in the beast123() interface. Some examples are:
%            deseasonalize <-> metadata.deseasonalize
%            scp.minmax(1) <-> prior.seasonMinOrder
%            scp.minmax(2) <-> prior.seasonMaxOrder
%            sseg.min      <-> prior.seasonMinSepDist
%            mcmc.seed     <-> mcmc.seed
%            tcp.minmax(1) <-> prior.trendMinKnotNumber
%   <strong>Experts should use the the beast123 function.</strong>
%
%   <strong>*Result/Output*</strong>: The output is a struct variable; example of the fields include
%
%       marg_lik: marginal likilood; the larger, the better
%       sig2    : variance  of error
%       trend   : the trend component; a struct variable (say, T)
%       season  : the season componet; a stuct variable  (say,S)
%       The subfields of trend or season:
%       .ncpPr        : the prob distribution for number of changepoints
%       .ncp          : mean number of changepoints in trend or seasonality
%       .ncp_meidan   : median number of changepoints
%       .ncp_mode     : mode from ncpPr
%       .ncp_pct90    : 90% percentile from ncpPr
%       .cpOccPr      : changepoint occurrance probability over time
%       .cp           : list of all possible changepoints (many are not sigficant)
%       .cpPr         : occurrence probability of the changepoints in cp
%       .cpAbruptChange: the sudden changes in trend or seasonlity at cp
%       .cpCI         : confidence interval of the cps
%       .Y            : the fitted trend or seasonality 
%       .SD           : standard deviation of the fitted Y
%       .CI           : Credible interval of the fittted Y
%       .order   : the mean harmonic or polynomial orders estimated to fit the seasonal and trend     
%       trend.slp     : slope of the trend 
%       trend.slpSD   : standard dev of the estimated slope
%       trend.slpSgnPosPr: time-varying probability of the slope being postive
%       trend.slpSgnZeroPr: time-varying probability of the slope being 0
%       season.amp     : amplitue of the estiamted seasonality overtime
%       season.ampSD   : standard ev of the estiamated amplitude
%
%   <strong>More help</strong>:  
%      This terse help doc sucks (I know); so far, the best details are still the
%      R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
%      Matlab doesn't support keyword-style args, so Matlab's equivalent to R's 
%      beast(Y,<strong>start</strong>=1987,<strong>freq</strong>=1) is beast(Y,<strong>'start'</strong>, 1987, <strong>'freq'</strong>,1).
%      
%   <strong>Examples</strong>:
%       % Nile river annual streamflow: trend-only data
%       load('Nile.mat')              
%       o = beast(Nile, 'start', 1871, 'season','none') 
%       printbeast(o)
%       plotbeast(o)
%       
%       % Explicitly specify deltat=1. BEAST knows nothing about the unit
%       % of 1871 and 1.0 (i.e., 1871 years, 1871 seconds, or 1871 meters?) 
%       o=beast(Nile, 'start', 1871, 'deltat', 1.0, 'season','none') 
%
%       % start is given a date 1871-1 (Year-Mon). The time unit is
%       % then fractional/decimal year. delta=1.0 means 1.0 year
%       o=beast(Nile, 'start', [1871,1], 'deltat', 1.0, 'season','none') 
%
%       % period=0 means a trend-only signal, which is equivalent to season='none'
%       o=beast(Nile, 'start', 1871, 'deltat', 1, 'period', 0) 
%
%       % Use a string to specify a unit for deltat or period (e.g., deltat='1 year')
%       % The time unit is also fractional year. 1871 means Year 1871
%       o=beast(Nile, 'start', 1871, 'deltat', '1 year', 'period', 0) 
%
%       % Use a string to specify a unit for delta or period (e.g., deltat='12 mo')
%       % The time unit is fractional year. 1871 means Year 1871
%       o=beast(Nile, 'start', 1871, 'deltat', '12 mo', 'period', 0) 
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
%       o=beast(beach, 'start', [2004,1],'deltat', 1/12, 'period',1.0)  
%
%       % period='12 month': use a string to explicitly specify the unit              
%       o=beast(beach, 'start', [2004,1],'deltat', 1/12, 'period','12 month')
%       o=beast(beach, 'start', [2004,1],'deltat', '1 month', 'period','365 days')
%
%       %% Monthly air co2 data since 1959: deltaTime=1/12 year
%       load('co2.mat')     
%       o=beast(co2, 'start', [1959,1,15], 'deltat', 1/12, 'period',1.0)
%       printbeast(o)
%       plotbeast(o)
%       plotbeast(o,'ncpStat','median')
%
%      %% Daily covid-19 infection statistics 
%       load('covid19.mat')    
%       Y    = sqrt(double(covid19.newcases));
%       Date = covid19.datestr;
%       % the min length of seasonal segments is set to 30 data points
%       o=beast(Y, 'start',[2020,01,22], 'deltat', 1/365, 'period', 7/365,'sseg.min',30)
%       printbeast(o)
%       plotbeast(o)
%       plotbeast(o,'ncpStat','median')
%
%       % Use a string to specify delta with a unit
%       o=beast(Y, 'start',[2020,01,22], 'deltat', '1.0 day',  'period', '7days','sseg.min',30)
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
    n=length(varargin);
    if mod(n,2)~=0
        msg=[ "the optional arg list must be paired keywords and values; the number of extra args must be even.\n", ...
              " Examples:  beast(y) \n", ...
              "            beast(y,'season','none','start',1980) \n", ...
             ];
        msg=sprintf(strcat(msg{:}));
        error(msg);
    end
    
    KeyList   = varargin(1:2:n);
    KeyList   = cellfun(@char,KeyList,'UniformOutput',false); % convert strings--if any--to chars,
    ValList   = varargin(2:2:n);
   %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % get values from keys. The last arg is the default value if the key is
   % missing from varagin/KeyList
  
   start  = GetValueByKey(KeyList, ValList, 'start',  []);
   deltat = GetValueByKey(KeyList, ValList, 'deltat', []);
   num_samples_per_period  = GetValueByKey(KeyList, ValList, 'freq',  []); 
   period   = GetValueByKey(KeyList, ValList, 'period',  []); 
   time    = GetValueByKey(KeyList, ValList, 'time',  []); 
   
   season          = GetValueByKey(KeyList, ValList, 'season',[]); %'harmonic'
   sorder_minmax   = GetValueByKey(KeyList, ValList, 'sorder.minmax', [1,5]); 
   scp_minmax      = GetValueByKey(KeyList, ValList, 'scp.minmax',    [0,10]); 
   sseg_min        = GetValueByKey(KeyList, ValList, 'sseg.min',      []); 
   
   deseasonalize   = GetValueByKey(KeyList, ValList, 'deseasonalize', false); 
   detrend         = GetValueByKey(KeyList, ValList, 'detrend', false); 
   
   torder_minmax   = GetValueByKey(KeyList, ValList, 'torder.minmax', [0,1]); 
   tcp_minmax      = GetValueByKey(KeyList, ValList, 'tcp.minmax',    [0,10]); 
   tseg_min        = GetValueByKey(KeyList, ValList, 'tseg.min',     []); 
   
   ocp             = GetValueByKey(KeyList, ValList, 'ocp',  []); 
   hasOutlierCmpnt = ~isempty(ocp);
   
   mcmc_seed       =GetValueByKey(KeyList, ValList, 'mcmc.seed',  0);         
   mcmc_samples    =GetValueByKey(KeyList, ValList, 'mcmc.samples',  8000);
   mcmc_thin       =GetValueByKey(KeyList, ValList, 'mcmc.thin',  5); 
   mcmc_burnin     =GetValueByKey(KeyList, ValList, 'mcmc.burnin',  200);
   mcmc_chainNumber=GetValueByKey(KeyList, ValList, 'mcmc.chains',   3);  
   
   ci               =GetValueByKey(KeyList, ValList, 'ci',   false);   
   printProgressBar =GetValueByKey(KeyList, ValList, 'print.progress',  true);     
   printOptions     =GetValueByKey(KeyList, ValList, 'print.options',  true);      
   gui              = GetValueByKey(KeyList, ValList, 'gui',  false); 
%% Convert the opt parameters to the individual option parameters (e.g.,
%  metadata, prior, mcmc, and extra)

   %......Start of displaying 'MetaData' ......
   metadata = [];
   metadata.isRegularOrdered = true;
   metadata.season           = season;
   metadata.startTime        = start;
   metadata.deltaTime        = deltat;
   if isempty(period) && ~isempty(deltat) && ~isempty(num_samples_per_period) && ~strcmp(season, 'none')
       period=num_samples_per_period*deltat;
   end   
   metadata.period           = period;
   metadata.time             = time;
 
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
       prior.seasonMinOrder   = sorder_minmax(1);
       prior.seasonMaxOrder   = sorder_minmax(2);
       prior.seasonMinKnotNum = scp_minmax(1);
       prior.seasonMaxKnotNum = scp_minmax(2);   
       prior.seasonMinSepDist = sseg_min;
   end   
   prior.trendMinOrder	  = torder_minmax(1);
   prior.trendMaxOrder	  = torder_minmax(2);
   prior.trendMinKnotNum  = tcp_minmax(1);
   prior.trendMaxKnotNum  = tcp_minmax(2);
   prior.trendMinSepDist  = tseg_min;
   
   prior.outlierMaxKnotNum=ocp;
        
   prior.precValue        = 1.500000;
   prior.precPriorType    = 'componentwise';
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
   extra.printProgressBar     = printProgressBar;
   extra.printOptions         = printOptions;
   extra.consoleWidth         = 70;
   extra.numThreadsPerCPU     = 2;
   extra.numParThreads        = 0;
%......End of displaying extra ......
%%
 if (gui)
    out=Rbeast('beastv4demo',y,metadata, prior,mcmc, extra);
 else
    out=Rbeast('beastv4',y,metadata, prior,mcmc, extra);
 end
end


%% Functions to return a default value if the field is missing from opt
function value=GetValueByKey(KeyList, ValList, key,defaultValue)
   idx=find(strcmp(KeyList,key));
   if isempty(idx)
       value=defaultValue;
   else
       value=ValList{idx(1)};
   end
end

 



 