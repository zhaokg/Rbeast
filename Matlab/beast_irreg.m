function out = beast_irreg(y, varargin)
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run 'help beast_irreg' to see the following
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   USAGE: out=<strong>beast_irreg(y, ...) </strong>
%
%   <strong>y</strong>:  iregular time series; it should be a numeric vector. For reggular 
%   time series, use 'beast' or 'beast123' instead. For multiple time series 
%    or stacked time series images such as satellite data, use 'beast123'.
%
%   <strong> ... </strong>:  a series of paired keywords and values to  specifiy time information 
%   or parameters for the BEAST algorithm. The keywords mimic the R version of beast
%   <a href="matlab:
%   web(https://rdrr.io/cran/Rbeast/man/beast.irreg.html')">rdrr.io/cran/Rbeast/man/beast.irreg.html</a>. Unlike R, Matlab doesn't support keyword-style 
%   arguments, so the beast parameters should be provided in the following forms:
%
%   <strong>beast(Nile, 'start', 1871, 'deltat', 1, 'season','none')</strong>
%   <strong>beast(Yellowstone, 'start', [1981,7,7], 'tcp.minmax', [0,10], 'deltat', 1/12)</strong> 
%
%   <strong>Possible Keywords</strong>:
%      
%   <strong>time</strong>: 
%        this is the only parameter that differs from those of the beast() function.     
%        'time' specifies the times of individual points of the input 'y'
%   <strong>deltat</strong>: 
%        the time interval between consecutive datapoints (e.g., 1/12
%        for monthly time series if the time unit is year). 'deltat' is
%        used to aggregate the unevely input into regular time seires
%        spaced by delta.
%        Note: why 'deltat' is needed for irregular time series:
%        internally, BEAST is implemented to handle only evenly-spaced data
%        for the sake of faster computation. So BEAST will aggregate it
%        into regular data at the specified interval of 'deltat' before
%        fitting the BEAST model. In the aggregration, data gaps/missing
%        values are possible, which is of no issue because BEAST can handle
%        missing values.
%   <strong>freq</strong>:  Deprecated. Replaced with 'period'. See below
%   <strong>period</strong>:  
%        a numeric value to specify the period if peridodic/seasonal variations 
%        are present in the data. If period is given a zero, negative value or 'none' 
%        it suggests no seasonal/periodic component in the signal. (season='none'
%        also suggests no periodic component).
%        Note: in earlier versions, 'freq' was used to specify the period and
%        now deprecated in this version. The unit of 'period', if any
%        should be consistent with the unit of 'deltat'.
%   <strong>season</strong>: 
%        a string specifier. Possible values - 'none':  trend-only data with no 
%        seasonality; 'harmonic': the seasonal/peridoic  component modelled via 
%        harmonic curves; 'dummy': the seasonal component  modelled via a dummy 
%        basis (i.e., pulse-like bases); 'svd': svd-derived  bases (available only in 
%        R not Matlab)
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
%   <strong>Note</strong>:
%       beast, beast_irreg, and beast123 calls the same library
%       (Rbeast.mex); they are just wrappers to the core BEAST algorithm.
%       In particular, BEAST currently handles only regular time series;
%       irregular inputs will be aggregated first before applying BEAST.
%
%   <strong>Result/Output</strong>: The output is a struct variable; example of the fields include
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
%      
%   <strong>Examples</strong>:
%
%       load('ohioNDVI')   %an unevenly-spaced Landsat NDVI time series
%       ohio               %a struct variable consiting of multiple fields
%       y    = ohio.ndvi
%       time = ohio.time
%       plot(time,y,'o')
%       % time is in the unit of fractional/decimal year
%       o=beast_irreg(y, 'time',time,'deltat', 1/12, 'period',1.0)
%       printbeast(o)
%       plotbeast(o)
%
%       
%       o=beast_irreg(y, 'time',ohio.datestr1) %parse date strings automatically
%       o=beast_irreg(y, 'time',ohio.datestr2) %parse date strings automatically
%       o=beast_irreg(y, 'time',ohio.datestr3) %parse date strings automatically
%       o=beast_irreg(y, 'time',ohio.datestr4) %parse date strings automatically
% 
%       time=[]
%       time.year=ohio.Y
%       time.month=ohio.M
%       time.day =ohio.D
%       o=beast_irreg(y, 'time', time) % provide numeric vectors of year, month, and days
%
%       time=[]
%       time.datestr=ohio.datestr1
%       time.strfmt= '????yyyy?mm?dd'
%       o=beast_irreg(y, 'time', time) % specify the pattern of data strings
%       
%
%   <strong>Contact info</strong>: To report bug or get help, do not hesitate to contact Kaiguang Zhao
%   at <strong>zhao.1423@osu.edu</strong>.
%
%   See also beast123, beast, printbeast, plotbeast, extractbeast


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
  
   time     = GetValueByKey(KeyList, ValList, 'time',   []);   
   start     = GetValueByKey(KeyList, ValList, 'start', []);
   deltat   = GetValueByKey(KeyList, ValList, 'deltat', []);
   num_samples_per_period  = GetValueByKey(KeyList, ValList, 'freq',  []); 
   period   = GetValueByKey(KeyList, ValList, 'period',  []); 
   
   season          = GetValueByKey(KeyList, ValList, 'season',[]); %'harmonic'
   sorder_minmax   = GetValueByKey(KeyList, ValList, 'sorder.minmax', [1,5]); 
   scp_minmax      = GetValueByKey(KeyList, ValList, 'scp.minmax',    [0,10]); 
   sseg_min        = GetValueByKey(KeyList, ValList, 'sseg.min',      []); 
   
   deseasonalize   = GetValueByKey(KeyList, ValList, 'deseasonalize', false); 
   detrend         = GetValueByKey(KeyList, ValList, 'detrend', false); 
   
   torder_minmax   = GetValueByKey(KeyList, ValList, 'torder.minmax', [0,1]); 
   tcp_minmax      = GetValueByKey(KeyList, ValList, 'tcp.minmax',    [0,10]); 
   tseg_min        = GetValueByKey(KeyList, ValList, 'tseg.min',  []); 

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
   %metadata.isRegularOrdered = false;
   metadata.season           = season;
   if strcmp(metadata.season, 'svd')
    % warning("season=svd is supported only for regular time series in beast()! 'harmonic' is used instead.");
	% metadata.season  = 'harmonic';
   end
   metadata.time             = time;
   metadata.startTime        = start;
   metadata.deltaTime        = deltat;
   if isempty(period) && ~isempty(deltat) && ~isempty(num_samples_per_period) && ~strcmp(season, 'none')
       period=num_samples_per_period*deltat;
   end   
   metadata.period           = period;  
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
   extra.computeSeasonAmp     = true;
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

 