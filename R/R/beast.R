beast <- function(  y,                        
				    start      = 1, 
					deltat     = 1, 
					season     = c('harmonic','svd', 'dummy','none'),
					period     = NULL,                  
					scp.minmax = c(0,10), sorder.minmax=c(0,5), 
					tcp.minmax = c(0,10), torder.minmax=c(0,1), 
					sseg.min = NULL, sseg.leftmargin = NULL,  sseg.rightmargin = NULL, 
					tseg.min = NULL, tseg.leftmargin = NULL,  tseg.rightmargin = NULL, 
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
  
  # start  = 1; deltat = 1;  freq   = NA       
  # season ='harmonic';           
  # scp.minmax=c(0,10); sorder.minmax = c(0,5); sseg.min=3  
  # tcp.minmax=c(0,10); torder.minmax = c(0,1); tseg.min=3 
  # detrend   = FALSE;  deseasonalize = FALSE
  # mcmc.seed = 0;      mcmc.burin    = 200; mcmc.chains=3; mcmc.thin=5; mcmc.samples=8000
  
  season        = match.arg(season)
  precPriorType = match.arg(precPriorType)
  
  # list is supported in this version for the multivariate cases
  if ( !hasArg("y") )  {                # if ( !hasArg("y") || is.list(y) ) 
    stop(" The input 'y' is missing!") 	# stop("Something is wrong with the input 'y'. Make sure that y is a vector")
    invisible(return(NULL))         
  }  
  
  if ( is.matrix(y) )  {
    dims = dim(y)
	if (dims[1]>1 && dims[2]>1) {	
		stop("If there are multiple time series to process (e.g., stacked images), pls use the beast123() function. Type ?beast123 for more information!")
		invisible(return(NULL))
	}	
	y=as.vector(y);
  } 

  if ( length(y) == 1 ) {
  	stop("Something is wrong with the input 'y'. Make sure that y is a vector")
	invisible(return(NULL))
  }  
   
  yclass = class(y); 
  if (  sum( yclass=='ts' | yclass=='zoo' | yclass=='xts' ) > 0  ){      
     if ( sum(yclass=='ts') > 0 ) {
		  tsp     = attributes(y)$tsp
		  start   = tsp[1]
		  end     = tsp[2]
		  deltat  = (end-start)/(length(y)-1)
		  freq_ts = tsp[3] 
		  period  = freq_ts *deltat
		  if ( freq_ts == 1 && season != 'none'){
		    season = 'none'	;
            period = NULL;
		    if (!quiet) warning("The input is a object of class 'ts' with a frequency of 1: trend-only data with no periodic component; season='none' is used.");		   
		  } else if (freq_ts >1 && season=='none'){
		    season = 'harmonic'
			msg    = sprintf("The input is a object of class 'ts' with a frequency of %d (i.e., with a periodic component); season='harmonic' is used instead.", freq_ts)
			if (!quiet)  warning(msg);
		  }
	 } else {
	    if (!quiet) warning("The input is a object of class 'zoo' or 'xts'. Its time attributes are not ignored, and only the data vector is used.");  
    	y = as.vector(y);
	 }     
  }
 

  #################################################################################
  syscall = sys.call()
  if( length(syscall) == 3  ) {
       call3    = syscall[[3]]
       arg2     = eval(call3,envir = parent.frame())
	   argname2 = names(syscall)[[3]]
	   
	   if( is.null(argname2) ){argname2=''}
       if (is.numeric(arg2) && length(arg2)==1 && argname2=='') {  
		  if (!quiet) warning(sprintf('Switching to the old interface of Rbeast v0.2: beast(Y,freq=%d)\n', arg2));
		  freq_old=arg2;
          invisible( return( beast.old(y,freq_old) ) )
       }
     
       if (is.list(arg2)) {
        
          valname2 = deparse(syscall[[3]])
          if (  argname2 == '' && valname2 =='opt' )     {
              s=sprintf('Switching to the old interface of Rbeast v0.2: beast(Y,opt)\n'); 
              if (!quiet)warning(s);
			  opt=arg2;
			  invisible( return( beast.old(y,opt) )  )
          } else if (argname2 == 'option' )    {
              s  =sprintf('Switching to the old interace of Rbeast v0.2: beast(Y,option)\n');  
			  if (!quiet) warning(s);
			  opt=arg2;
			  invisible( return( beast.old(y,option=opt) )  )
          }
        } 
        
     }
	 
   # eval(substitute(alist(...)))
   # substitute(alist(...))
   # ...names()
   # ...length()
   # ...1   
   # tmplist = list(...)	
   
  #################################################################################  

 #......Start of displaying 'MetaData' ......
   metadata = list()
   metadata$isRegularOrdered = TRUE
   metadata$season           = season   
   metadata$startTime        = start
   metadata$deltaTime        = deltat
   metadata$period           = period;
 
   # This is the old interface that uses freq to specify the period
   if ( is.null(period) && is.numeric(deltat) && season != 'none' && hasArg('freq')){
        freq                = list(...)[['freq']]  
		metadata$period     = deltat*freq		       
   }
   
   #if ( season=='svd' ){       
       # freq       = metadata$period/deltat
	   #	notInteger = abs(as.integer(freq) - freq) > 0.00001
	   #	if (notInteger || is.null(freq) || freq<=1.1 || is.na(freq) ) {
	   #  	stop('When season=svd, freq (i.e., period/deltat) must be an integer larger than 1.')
	   #    invisible(return(NULL))
	   #}
      # metadata$svdTerms=svdbasis(y,freq,deseasonalize)
   #}
   
   metadata$deseasonalize     = deseasonalize
   metadata$detrend           = detrend
   #metadata$whichDimIsTime   = 1
   metadata$missingValue      = NaN
   metadata$maxMissingRate    = 0.75
   if ( hasArg('hasOutlier') ) {
        hasOutlier               = list(...)[['hasOutlier']]   
		metadata$hasOutlierCmpnt = as.logical(hasOutlier)		           
   }
   
#......End of displaying MetaData ......
   prior = list()
   prior$modelPriorType	  = 1   
   if (!(season=='none')){
    prior$seasonMinOrder   = sorder.minmax[1]
	prior$seasonMaxOrder   = sorder.minmax[2]
    prior$seasonMinKnotNum = scp.minmax[1]
    prior$seasonMaxKnotNum = scp.minmax[2]
	if (!is.null(sseg.min) && !is.na(sseg.min))                 prior$seasonMinSepDist = sseg.min
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
		metadata$hasOutlierCmpnt= TRUE	
        prior$outlierMaxKnotNum	= list(...)[['ocp']] 
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
   mcmc$maxMoveStepSize = 0    # if set to zero, a default value in reference to freq will be used
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
   extra$computeSeasonAmp     = season != 'svd' # Not implemented for SVD
   extra$computeTrendSlope    = TRUE
   extra$tallyPosNegSeasonJump= TRUE
   extra$tallyPosNegTrendJump = TRUE
   extra$tallyIncDecTrendJump = TRUE
   extra$printProgressBar     = print.progress
   extra$printOptions         = TRUE
   extra$consoleWidth         = 0
   extra$quiet                = quiet
   if ( hasArg('randcoeff') ) {    
        extra$useMeanOrRndBeta	= list(...)[['randcoeff']] 
   }  else {
        extra$useMeanOrRndBeta=TRUE
   }
      
   #extra$numThreadsPerCPU     = 2
   #extra$numParThreads        = 0
 
 
 if (gui && !base::interactive()) {
	warning('R is not running in the inteactive mode. Resetting gui to FALSE.');
	gui = FALSE
 }
 
 if (is.list(y)){
	gui = FALSE
 }
 
 method = match.arg(method) 
 funstr = ifelse(gui,"beastv4demo", paste0("beast_",method)) 
 
 if ( hasArg("cputype") )  {
    cputype = list(...)[['cputype']]
	cputype = switch(cputype, sse=1, avx2=2, avx512=3);
	ANS    = .Call( BEASTV4_rexFunction, list(funstr,y,metadata,prior,mcmc,extra,cputype),   212345)   		   
 } else {
    if (hasArg("local")){ 
	  # run the local developer's version of Rbeast
	    dyn.load('y:/testold/Rbeast.mexw64')
		ANS  = .Call( "rexFunction1",      list(funstr,y,metadata,prior,mcmc,extra),   212345, PACKAGE="Rbeast.mexw64")  
	} else{
	    ANS  = .Call( BEASTV4_rexFunction, list(funstr,y,metadata,prior,mcmc,extra),   212345)   		   
	}
 }
 invisible(return(ANS))    
}

# meanfilter <- function(x,n=5){filter(x,rep(1,n), sides=2)}

 
 
 
