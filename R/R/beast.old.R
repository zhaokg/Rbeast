beast.old <- function(y, option=list() ) {
  
  if (!hasArg("y") || is.list(y) || length(y)==1){
    warning("Something is wrong with the input 'Y'. Make sure the data par is a vector or matrix.")
    return(NULL)
  }
  
  if( is.numeric(option)&&(length(option)==1) )   {
          
    #......Start of displaying 'MetaData' ......
    metadata = list()
    metadata$isRegularOrdered = TRUE
    metadata$period           = option
  
    #......Start of displaying 'extra' ......
    extra = list()
    extra$dumpInputData        = TRUE
    extra$computeCredible      = TRUE
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
    extra$printProgressBar     = TRUE
    extra$printOptions         = TRUE
    extra$consoleWidth         = 0
    #extra$numThreadsPerCPU     = 2
    #extra$numParThreads        = 0
    
       
    #ANS    = .Call( "rexFunction1", list("beastv4",y,metadata,prior=NULL,mcmc=NULL,extra),   212345,PACKAGE="Rbeast.mexw64")  
    ANS    = .Call( BEASTV4_rexFunction, list("beast_bayes",y,metadata,prior=NULL,mcmc=NULL,extra),   212345)   		   
    invisible(return(ANS))    
  }
  else if (is.list(option)) {
    opt=option
    metadata = list()
    metadata$isRegularOrdered = TRUE
    metadata$season           = 'harnomic'   
    if (!is.null(opt$period))        metadata$period   = opt$period
    if (!is.null(opt$omittedValue))  metadata$missingValue   = opt$omittedValue
  
    #......End of displaying MetaData ......
    prior = list()
    if (!is.null(opt$minSeasonOrder))       prior$seasonMinOrder   = opt$minSeasonOrder
    if (!is.null(opt$maxSeasonOrder))       prior$seasonMaxOrder   = opt$maxSeasonOrder
    if (!is.null(opt$minTrendOrder))        prior$trendMinOrder   = opt$minTrendOrder
    if (!is.null(opt$maxTrendOrder))        prior$trendMaxOrder   = opt$maxTrendOrder   
    if (!is.null(opt$maxKnotNum_Season))    prior$seasonMaxKnotNum   = opt$maxKnotNum_Season
    if (!is.null(opt$maxKnotNum_Trend))     prior$trendMaxKnotNum   = opt$maxKnotNum_Trend  
    if (!is.null(opt$minSepDist_Season))    prior$seasonMinSepDist   = opt$minSepDist_Season
    if (!is.null(opt$minSepDist_Trend))     prior$trendMinSepDist   = opt$minSepDist_Trend  
    
 
    mcmc=list()
    if (!is.null(opt$seed))       mcmc$seed   = opt$seed
    if (!is.null(opt$chainNumber))       mcmc$chainNumber   = opt$chainNumber
    if (!is.null(opt$sample))       mcmc$samples   = opt$sample
    if (!is.null(opt$thinningFactor))       mcmc$thinningFactor   = opt$thinningFactor
    if (!is.null(opt$burnin))       mcmc$burnin   = opt$burnin
    if (!is.null(opt$maxMoveStepSize))       mcmc$maxMoveStepSize   = opt$maxMoveStepSize
    if (!is.null(opt$resamplingSeasonOrderProb)) mcmc$seasonResamplingOrderProb=opt$resamplingSeasonOrderProb
    if (!is.null(opt$resamplingTrendOrderProb)) mcmc$trendResamplingOrderProb= opt$resamplingTrendOrderProb
    #......Start of displaying 'mcmc' ......
 
     
    #......Start of displaying 'extra' ......
    extra = list()
    extra$dumpInputData        = TRUE
    #extra$whichOutputDimIsTime = 1
    extra$computeCredible      = TRUE
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
    extra$printProgressBar     = TRUE
    extra$printOptions         = TRUE
    extra$consoleWidth         = 0
    #extra$numThreadsPerCPU     = 2
    #extra$numParThreads        = 0
	
	#ANS    = .Call( "rexFunction1", list("beastv4",y,metadata,prior=NULL,mcmc=NULL,extra),   212345,PACKAGE="Rbeast.mexw64")  	   
    ANS    = .Call( BEASTV4_rexFunction, list("beast_bayes",y,metadata,prior=NULL,mcmc=NULL,extra),   212345)   		   
    invisible(return(ANS))  
  }
  else{
    stop("The input arguments are unrecongized.")
    invisible(return(NULL))  
  }
    
 
 

 

  
 
 
  
  
  
}