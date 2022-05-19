function out = beast123(Y, metadata, prior, mcmc, extra)
%  
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run 'help beast123' to see the following
%   USAGE: out=<strong>beast123(Y, metadata, prior, mcmc, extra) </strong>
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   <strong>Y</strong>:  input time series being regular or inrregular; it could be a numeric vector,
%    2D matrix (e.g., multiple time series of the same length), or 3D array time series (e.g.,
%    stacked  satellite images). The interface for beast123 is similar to the
%    documentation of its R version at <a href="matlab:
%   web(https://rdrr.io/cran/Rbeast/man/beast123.html')">rdrr.io/cran/Rbeast/man/beast123.html</a>.
%
%   <strong>metadata</strong>:  a struct variable specifying metadata for the input Y 
%
%   metadata.isRegularOrdered : if true, Y is regular ts; if false, Y is
%                               unevenly spaced time series
%   metadata.season           : a string specifier. 'none' - trend-only  data,
%                               'harmonic' - harmonic model for the seasonal component,     
%                               'dummy' - dummy model for the seasonal component     
%   metadata.period           : the period of the seasonal/periodic  comonent
%   metadata.startTime        : the start time of regular input
%   metadata.deltaTime        : the time interval between consecutive data
%                               points. For regular time seires, it
%                               determins the times of all the data points;
%                               for irregular time series, deltaTime
%                               specifies the interval at which the
%                               irregular input is aggregated/resampled
%   metadata.time             : for irregular inputs only, specify the
%                               individual times of data points. Both
%                               string and numeric vectors are allowed,
%                               check rdrr.io/cran/Rbeast/man/beast123.html
%                               for more details.
%   metadata.whichDimIsTime   : which dim of the 2D or 3D input refer to
%                               time
%   metadata.missingValue     : a value indicating bad/missing data
%   metadata.maxMissingRate   : the max missingness rate beyond which the
%                               ts will be ignored
%   metadata.deseasonalize    : if true, the input ts is first
%                              de-seasonalize by removing a global seasonal
%                              component, prior to applying BEAST
%   metadata.detrend          : if true, the input ts is first
%                              de-trended by removing a global trend 
%                              component, prior to applying BEAST
%
%
%   <strong>prior</strong>:  prior hyperparameter for the BEAST model
%  
%
%   prior.seasonMinOrder  : min harmonic order considered for the seasonal
%                         compnt
%   prior.seasonMaxOrder  : max harmonic order considered for the seasonal
%                         compnt
%   prior.seasonMinKnotNum : min number of seasonal changepints allowed
%   prior.seasonMaxKnotNum : max number of seasonal changepints allowed
%   prior.seasonMinSepDist : the min segment length of a seaosnal segment
%                           (i.e, the min seperation distance between
%                           neighorboring seasonal changepoints)
%   prior.trendMinOrder	   : min polyonimal order considered for the trend
%                         compnt
%   prior.trendMaxOrder	   : max polyonimal order considered for the trend
%                         compnt
%   prior.trendMinKnotNum  : min number of trend changepints allowed
%   prior.trendMaxKnotNum  : max number of trend changepints allowed
%   prior.trendMinSepDist  : the min segment length of a trend segment
%                           (i.e, the min seperation distance between
%                           neighorboring trend changepoints)
%
%   <strong>mcmc</strong> :   parameters to set up the MCMC sampler
%
%   mcmc.seed             : the seed for the random number generator            
%   mcmc.samples          : the number of samples to be collected
%   mcmc.thinningFactor   : chain thinning factor; take every
%                           'thinningFactor'-th sample from the chain
%   mcmc.burnin           : number of initial samples discarded
%   mcmc.chainNumber      : number of chains
%   mcmc.maxMoveStepSize  : a moving step
%   mcmc.trendResamplingOrderProb  : probability to re-sample trend order
%   mcmc.seasonResamplingOrderProb : probability to re-sample seasonal order
%   mcmc.credIntervalAlphaLevel    : signifiance level for computing
%                                    credible interval
%
%   <strong>extra</strong> :   parameters to control computaiton or outputs
%   extra.dumpInputData      : if true, dump the input time series (o.data)
%   extra.whichOutputDimIsTime: which dim of the output array refer to time
%   extra.computeCredible     : compute credible intervals (unless needed
%                               for assessing the curve fitting, it is
%                               better to disable it for changepoint
%                               analysis because computing CI is
%                               expensive.)(o.trend.CI and o.season.CI)
%   extra.fastCIComputation   : if true, a fast version of algorithm used to
%                              compute CI. But it is still expensive.
%                              Set computeCredible=false if CI is of litte
%                              interest
%   extra.computeSeasonOrder  : if true, estimate harmonic orders needed to
%                              adequately fit the seasonal component
%                              (o.season.order)
%   extra.computeTrendOrder   : if true, estimate polynomial orders needed to
%                              adequately fit the trend component
%                              (o.trend.order)
%   extra.computeSeasonChngpt  : if true, dump the scp detected
%                              (o.season.cp, o.season.cpCI, and cpAbruptChange)
%   extra.computeTrendChngpt   : if true, dump the tcp detected
%                              (o.trend.cp, o.trend.cpCI, and cpAbruptChange)
%   extra.computeSeasonAmp     : if true, estimate the seasonal ammplitue
%                              (o.seasin.amp)
%   extra.computeTrendSlope    : if true, estimate time-varying slope and
%                               the associated probabilities of slopes being positive or negative
%   extra.tallyPosNegSeasonJump:
%   extra.tallyPosNegTrendJump :
%   extra.tallyIncDecTrendJump :
%   extra.printProgressBar     : if true, show a progress bar
%   extra.printOptions         : if true, print the BEAST parameters at the
%                               start
%   extra.consoleWidth         : the console/terminal width
%   extra.numThreadsPerCPU     : IMPORTANT for remote sensing/earth science
%                               applicaitons: specify the numbers of
%                               concurrent threads per cpu core
%   extra.numParThreads        : specify the number of total concurrent
%                              threads
%
%
%   <strong>Note</strong>:
%       beast, beast_irreg, and beast123 calls the same library
%       (i.e., Rbeast.mex); they are just wrappers to the core BEAST algorithm.
%       In particular, BEAST currently handles only regular time series;
%       irregular inputs will be first aggregated  before applying BEAST.
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
%      The terse doc sucks (I know); so far, the best details are still the
%      R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
%      Matlab doesn't support keyword-style args, so Matlab's equivalent to R's beast(freq=1) 
%      beast('freq',1).
%           
%   <strong>Examples</strong> :   
%
%       load('Nile')   % Nile river annual streamflow: trend-only data
%       metadata=[];
%       metadata.isRegularOrdered=true;
%       metadata.season     = 'none'  % trend-only
%       metadata.startTime  =1871;
%       metadata.deltaTime =1;
%       % Default values will be used if parameters are missing
%       o=beast123(Nile,metadata) 
%       printbeast(o)
%       plotbeast(o)
%
%
%       load('ohioNDVI')   % irregular Landsat NDVI time series
%       metadata=[];
%       metadata.isRegularOrdered =false;
%       metadata.season     =  'harmonic'  % seasonality is present
%       metadata.time       =  ohio.time;
%       metadata.deltaTime  = 1/12;   % aggregate into a monthly ts
%       metadata.period     = 1.0;    % period= 1/12 x 12 (freq)
%       % Default values will be used if parameters are missing
%       o=beast123(ohio.ndvi,metadata,[],[], struct('dumpInputData', true)) 
%       printbeast(o)
%       plotbeast(o, 'ncpStat','median')       
%       
%       load('ohioNDVI')   % irregular Landsat NDVI time series
%       ohio.datestr1      % strings of times for ohio.ndvi
%       metadata=[];
%       metadata.isRegularOrdered =false;
%       metadata.time       =  [];
%       metadata.time.dateStr =ohio.datestr1;
%       metadata.time.strFmt  ='????yyyy?mm?dd';
%       metadata.deltaTime  = 1/12;   % aggregate into a monthly ts
%       metadata.period     = 1.0;    % period= 1/12 x 12 (freq)
%       % Default values will be used if parameters are missing
%       o=beast123(ohio.ndvi,metadata,[],[], struct('dumpInputData', true)) 
%       plotbeast(o, 'ncpStat','median')   
%
%       ohio.datestr2      % strings of times for ohio.ndvi
%       metadata=[];
%       metadata.isRegularOrdered =false;
%       metadata.time       =  [];
%       metadata.time.dateStr =ohio.datestr2;
%       metadata.time.strFmt  ='LC8-yyyydoyxdvi';
%       metadata.deltaTime  = 1/12;   % aggregate into a monthly ts
%       metadata.period     = 1.0;    % period= 1/12 x 12 (freq)
%       o=beast123(ohio.ndvi,metadata) 
%       plotbeast(o, 'ncpStat','median')  
%       %See https://rdrr.io/cran/Rbeast/man/beast123.html for more details
%       %about the accepted formats of date strings
%
%       load('simData.mat') 
%       % A toy example of 3 time series: I decides to be lazy here bcz the 3
%       % time series are essentially the same. We just assume they are
%       % different for illustring the use of beast123 to handle 2d matrix      
%       simData.Y % a 774×3 matrix: 774 is the time dimesnion      
%       metadata=[];      
%       metadata.whichDimIsTime =1 % 774 is the ts length 
%       o=beast123(simData.Y ,metadata) 
%       printbeast(o,  1)         %print the result for the first ts
%       plotbeast(o, 'index', 1)  %plot the result for the first ts
%       printbeast(o,  2)         %print the result for the 2nd ts
%       plotbeast(o, 'index', 2)  %plot the result for the 2nd ts      
%
%      
%       Ytran   =simData.Y' % a 3x774 matrix: 774 is the time dimesnion      
%       metadata=[];      
%       metadata.whichDimIsTime = 2 % 774 is the ts length 
%       o=beast123(Ytran ,metadata) 
%       printbeast(o,  2)         %print the result for the 2nd ts
%       plotbeast(o, 'index', 2)  %plot the result for the 2nd ts   
%
%      load('imageStack.mat') 
%       % A toy example of stacked time series images: unevely-spaced in time
%      NDVI3D=imageStack.ndvi    % a 12x9x1066 3D cube
%      TIME  =imageStack.datestr % 1066 is the time series length%
%      metadata=[];      
%      metadata.isRegularOrdered =false % irregular input
%      metadata.whichDimIsTime =3 % 1066 is the ts length 
%      metadata.time=[];
%      metadata.time.dateStr=TIME
%      metadata.time.strFmt='LT05_018032_20110726.yyyy-mm-dd';
%      metadata.deltaTime  =1/12; % aggregated at a monthly interval
%      metadata.period     =1.0;  % the period is 1.0 (year)
%      extra=[];
%      extra.dumpInputData   =true % get a copy of the aggregated input
%      extra.numThreadsPerCPU = 2; % 2 threads per CPU core
%      o=beast123(NDVI3D,metadata,[],[], extra) 
%      imagesc(o.sig2)
%      imagesc(o.trend.ncpPr(:,:,1:3))
%      printbeast(o,[2,4]) %print the result at row 2 and col 4     
%      plotbeast(o,'index',[2,4]) %plot the result at row 2 and col 4     
%
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
    
    if nargin>=2 && (isstring(metadata) || ischar(metadata))
           seasonType   = char(metadata);
           metadata     = struct('season',seasonType);
    end    
    
   if (nargin==0)          
       error("The input should not be empty; at least a time series needs to be provided.");
   end  
%%
    out=Rbeast('beastv4',Y,metadata, prior,mcmc, extra);
end
