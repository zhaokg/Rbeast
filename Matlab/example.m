%% Get a test dataset
% Download testDat.mat from gitHub and then specify the path where the file is saved.
% The folder used below is just a placeholder.
load('a:\testData.mat');

%% Example1: Take a look at the downloaded dataset
% "testData" contains two time series: simData and landsat

% First, let us look at simData. It is a struct variable consisting of
% ".Y":          A simulated time-series signal
% ".trueSeason": The true seasonal component used in the simulation
% ".trueTrend":  The true trend component used in the simulation

% Of particulate note, this is the simulated time series used for Example 1 
% in our RSE paper (Zhao et al., 2019).

Y     = simData.Y;
S_true = simData.trueSeason;
T_true = simData.trueTrend;

clf
subplot(3,1,1)
plot(Y);
hold on;
plot(S_true+T_true);

subplot(3,1,2);
plot(S_true);

subplot(3,1,3);
plot(T_true);
%% Set up the parameters needed for the BEAST algorithm
% Some of these parameters are the model specficiation parameters of BEAST
% (e.g., minSeasonOrder, maxSeasonOrder, minSetpDist_trend,
% minSepDist_Season); other parameters are just some input variables to
% control simulation behaviors or program outputs (e.g., samples,
% thinningFactor, seed, computeCredible).
% 
opt.period    = 24;  
opt.minSeasonOrder = 1;
opt.maxSeasonOrder = 6;
opt.minTrendOrder=0;
opt.maxTrendOrder=1;
opt.minSepDist_Trend  =  24;
opt.minSepDist_Season =  24;
opt.maxKnotNum_Trend =  8;
opt.maxKnotNum_Season = 8; 
opt.maxMoveStepSize   = 40;
opt.samples = 10000;
opt.thinningFactor = 1;
opt.burnin = 200;
opt.chainNumber=2;
opt.resamplingTrendOrderProb=0.2;
opt.resamplingSeasonOrderProb=0.17;
opt.omissionValue=-999;
opt.seed=100;
opt.computeCredible=0;
opt.computeSlopeSign=1;
opt.algorithm='beast';
opt.computeHarmonicOrder=1;
opt.computeTrendOrder=1;
opt.computeChangepoints=1;
%opt.timeDimensionIndex=3;
%% Run BEAST on "Y"
tic
out=beast_default(Y, opt);
toc
%% Plot the result
clf

%plot the detected seasonal component
subplot(6,1,1); plot(out.s)     

%plot the seasonal-changepoint probability (i.e., the probability of observinng a seasoanl changepoint over time)
subplot(6,1,2); plot(out.sProb) 

%Go back to Subplot#1 and plot the most likely locations of seasonal
%changepoints (scp)
%------------------------------------------------------------------------------------------
subplot(6,1,1);
% out.scp saves the most probabe scp locations. It is filled with NaNs if
% the detected scp number is less than the specified
%"opt.maxKnotNum_Season".
numOfScp = length( find(~isnan(out.scp) )) ;
hold on;
for index=1:numOfScp
    scpLoc = out.scp(index);
    plot([scpLoc,scpLoc],[-2,2] ,'r');
end
%------------------------------------------------------------------------------------------

%BEAST also estimates the harmonic orders needed to sufficiently to
%approximate the seasoanl component. The result is outputted to out.horder.
%"horder" varies with time.
subplot(6,1,3);plot(out.horder);

%plot the detected trend
subplot(6,1,4);plot(out.t); 
%------------------------------------------------------------------------------------------
hold on;
% out.tcp saves the most probabe locations of trend changepoints (tcp). It is filled with NaNs if
% the detected tcp number is less than the specified "opt.maxKnotNum_Trend".
numOfTcp = length( find(~isnan(out.tcp) )) ;
hold on;
for index=1:numOfTcp
    tcpLoc = out.tcp(index);
    plot([tcpLoc,tcpLoc],[-.5,.5] ,'r');
end
%------------------------------------------------------------------------------------------

%plot the tcp probability (i.e., the probability of detecting a tcp over time)
subplot(6,1,5);plot(out.tProb);  

%BEAST also estimates the trend order needed to sufficiently to
%approximate the seasoanl component. The result is outputted to out.torder.
%"torder" varies with time.
% By "trend order", the linear segment can be either a constant (zero-th
% order) or a sloped line (1st order). "torder" gives the mean order over
% the sampled models.
subplot(6,1,6);plot(out.torder);  



