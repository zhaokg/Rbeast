#from . import Rbeast as rb

beast.irreg <- function(  
                    y,    
                    time, 
					deltat     = NULL , 
					period     = NULL,                  
                    season     = c('harmonic','svd','dummy','none'),
					scp.minmax = c(0,10), sorder.minmax=c(0,5), 
					tcp.minmax = c(0,10), torder.minmax=c(0,1), 
					sseg.min   = NULL, sseg.leftmargin = NULL,  sseg.rightmargin = NULL, 
					tseg.min   = NULL, tseg.leftmargin = NULL,  tseg.rightmargin = NULL, 
					method         = c('bayes','bic', 'aic','aicc','hic'),
                    detrend        = FALSE,
                    deseasonalize  = FALSE,
                    mcmc.seed      = 0,  mcmc.burnin=200, mcmc.chains=3, mcmc.thin=5,mcmc.samples=8000,
					ci             = FALSE,
                    precValue      = 1.5,
                    precPriorType  = c('componentwise','uniform','constant','orderwise'),					
                    print.options  = TRUE,
                    print.progress = TRUE,					
					quiet          = FALSE,
                    gui            = FALSE,
					...)
{

  #time=df$date
  #deltat=1/24
  #freq   = 24
  #season ='harmonic';
  #scp.minmax=c(0,10); sorder.minmax=c(0,5); sseg.min=3  
  #tcp.minmax=c(0,10); torder.minmax=c(0,1); tseg.min=3 
  #detrend = FALSE; deseasonalize=FALSE
  #mcmc.seed=0;  mcmc.bunrin=200; mcmc.chains=3; mcmc.thin=5; mcmc.samples=8000
  #print.options=TRUE
  #print.progress  =TRUE
  #gui=FALSE
  
  # list is supported in this version for the multivariate cases
  #if ( !hasArg("y") || is.list(y) )  {  
  if ( !hasArg("y") )  {  
    stop(" The input 'y' is missing!")
    #stop("Something is wrong with the input 'y'. Make sure that y is a vector")	
    invisible(return(NULL))         
  }  
  
  if ( is.matrix(y) )  {
    dims=dim(y);
	if (dims[1]>1 && dims[2]>1) {	
		stop("If there are multiple time series to process (e.g., stacked images), pls use the beast123() function. Type ?beast123 for more information!")
		invisible(return(NULL))
	}	
	y=as.vector(y);
  }  
  yClass=class(y); 
  if ( sum(yClass=='ts')>0 || sum(yClass=='zoo')>0 || sum(yClass=='xts')>0 )  {    
	y=as.vector(y);
  }  
  if (length(y)==1) {
  	stop("Something is wrong with the input 'y'. Make sure that y is a vector")
	invisible(return(NULL))
  }

 
  season        = match.arg(season)
  precPriorType = match.arg(precPriorType)
 # tmplist=list(...)
#......Start of displaying 'MetaData' ......
   metadata = list()
   metadata$isRegularOrdered = FALSE
   metadata$season           = season   
   metadata$time             = time
   #metadata$startTime       = start
   metadata$deltaTime        = deltat
   if ( season != 'none'){
	metadata$period           = period;
   }   
   #metadata$whichDimIsTime   = 1
   metadata$deseasonalize     = deseasonalize
   metadata$detrend           = detrend
   metadata$missingValue      = NaN
   metadata$maxMissingRate    = 0.7500
   if ( hasArg("maxMissingRate") ) {
		metadata$maxMissingRate    = list(...)[['maxMissingRate']]
   }   
   if ( hasArg('hasOutlier') ) {   
        hasOutlier =list(...)[['hasOutlier']]
		metadata$hasOutlierCmpnt=as.logical(hasOutlier)		           
   }		     
   if ( is.null(period) && is.numeric(deltat) && season != 'none' && hasArg('freq')){
        freq                = list(...)[['freq']]  
		metadata$period     = deltat*freq		       
   }
		
#......End of displaying MetaData ......
   prior = list()
   prior$modelPriorType	  = 1   
   if (!(season=='none')){
    prior$seasonMinOrder   = sorder.minmax[1]
	prior$seasonMaxOrder   = sorder.minmax[2]
    prior$seasonMinKnotNum = scp.minmax[1]
    prior$seasonMaxKnotNum = scp.minmax[2]
    if (!is.null(sseg.min) && !is.na(sseg.min))   prior$seasonMinSepDist = sseg.min
	if (!is.null(sseg.leftmargin)  && !is.na(sseg.leftmargin))  prior$seasonLeftMargin = sseg.leftmargin
	if (!is.null(sseg.rightmargin) && !is.na(sseg.rightmargin)) prior$seasonRightMargin = sseg.rightmargin
   }   
   prior$trendMinOrder	  = torder.minmax[1]
   prior$trendMaxOrder	  = torder.minmax[2]
   prior$trendMinKnotNum  = tcp.minmax[1]
   prior$trendMaxKnotNum  = tcp.minmax[2]
   if (!is.null(tseg.min) && !is.na(tseg.min))                 prior$trendMinSepDist = tseg.min   
   if (!is.null(tseg.leftmargin)  && !is.na(tseg.leftmargin))  prior$trendLeftMargin  = tseg.leftmargin
   if (!is.null(tseg.rightmargin) && !is.na(tseg.rightmargin)) prior$trendRightMargin = tseg.rightmargin
   
   if ( hasArg('ocp') ) {   
		metadata$hasOutlierCmpnt = TRUE	
        prior$outlierMaxKnotNum	 =list(...)[['ocp']] 
   }  
   prior$K_MAX            = 0
   prior$precValue        = precValue
   prior$precPriorType    = precPriorType
   
 
#......End of displaying pripr ......

#......Start of displaying 'mcmc' ......
   mcmc = list()
   mcmc$seed            = mcmc.seed
   mcmc$samples         = mcmc.samples
   mcmc$thinningFactor  = mcmc.thin
   mcmc$burnin          = mcmc.burnin
   mcmc$chainNumber     = mcmc.chains
   mcmc$maxMoveStepSize = 0;
   mcmc$trendResamplingOrderProb  = 0.1000
   mcmc$seasonResamplingOrderProb = 0.1700
   mcmc$credIntervalAlphaLevel    = 0.950
#......End of displaying mcmc ......

#......Start of displaying 'extra' ......
   extra = list()
   extra$dumpInputData        = TRUE
   #extra$whichOutputDimIsTime = 1
   extra$computeCredible      = ci
   extra$fastCIComputation    = TRUE
   extra$computeSeasonOrder   = TRUE
   extra$computeTrendOrder    = TRUE
   extra$computeSeasonChngpt  = TRUE
   extra$computeTrendChngpt   = TRUE
   extra$computeSeasonAmp     = TRUE
   extra$computeTrendSlope    = TRUE
   extra$tallyPosNegSeasonJump= TRUE
   extra$tallyPosNegTrendJump = TRUE
   extra$tallyIncDecTrendJump = TRUE
   extra$printProgressBar     = print.progress
   extra$printOptions         = print.options
   extra$consoleWidth         = 0
   extra$quiet                 = quiet
   #extra$numThreadsPerCPU     = 2
   #extra$numParThreads        = 0
 
  if (gui && !base::interactive()) {
	warning('R is not running in the inteactive mode. Resetting gui to FALSE.');
	gui = FALSE
 }


 method = match.arg(method) 
 funstr = ifelse(gui,"beastv4demo", paste0("beast_",method)) 
 
  if ( hasArg("cputype") )  {
    cputype = list(...)[['cputype']]  
	cputype = switch(cputype, sse=1, avx2=2, avx512=3);	
	ANS     = .Call( BEASTV4_rexFunction, list(funstr,y,metadata,prior,mcmc,extra,cputype),   212345)   		   
 } else {
    if (hasArg("local")){ # run the local developer's version of Rbeast
		#ANS  = .Call( "rexFunction1", list(funstr,y,metadata,prior,mcmc,extra),   212345, PACKAGE="Rbeast.mexw64")  
	} else{
	    ANS  = .Call( BEASTV4_rexFunction, list(funstr,y,metadata,prior,mcmc,extra),   212345)   		   
	}	    		   
 }
 		   
 invisible(return(ANS))    
}