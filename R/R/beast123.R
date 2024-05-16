beast123 <- function(Y, metadata=list(), prior=list(), mcmc=list(), extra=list(),
                     season = c('harmonic','svd','dummy','none'),
					 method = c('bayes','bic', 'aic','aicc','hic'), ... )
{
  # tmplist = list(...)
  # if (!hasArg("Y") || is.list(Y) || length(Y)==1)  {
  
  if (!hasArg("Y"))  {
    warning("Something is wrong with the input 'Y'. Make sure the Y par is a vector, a matrix, or an 3D array.")
    invisible(return(NULL))
  }

 season = match.arg(season)
 
 if ( season=='none') {
      if( is.numeric(metadata) )     {
    			tmp      = metadata
    			metadata = list();
    			if (length(tmp)==1){
    			  metadata$startTime=tmp(1)	
    			} else if( length(tmp)==2 ){
    			  metadata$startTime=tmp(1)
    			  metadata$deltaTime=tmp(2)
    			}         
      } else if (is.list(metadata)) {
        # do nothing 
      }
      else{
        # this branch should be never touched!
        metatata=list()
      }
 } else   {
      if( is.numeric(metadata) )    {
        tmp       = metadata
        metadata  = list();
        if(length(tmp)==1){
          metadata$period=tmp(1)	
        } 		 
      } else if (is.list(metadata)) {
        # do nothing 
      }  else{
        # this branch should be never touched!
        metadata=list()
      }
      
}

  if (is.null(metadata$season))
    metadata$season=season;
  
 method = match.arg(method) 
 funstr =  paste0("beast_",method)
 
  if ( hasArg("cputype") )  {
    cputype = list(...)[['cputype']]  
    cputype = switch(cputype, sse=1, avx2=2, avx512=3);	
	ANS     = .Call( BEASTV4_rexFunction, list(funstr,Y,metadata,prior,mcmc,extra,cputype),   212345)   		   
 } else {
    if (hasArg("local")){ # run the local developer's version of Rbeast
		#ANS  = .xxxCall( "rexFunction1", list(funstr,Y,metadata,prior,mcmc,extra),   212345, PACKAGE="Rbeast.mexw64")  
	} else{
	    ANS  = .Call( BEASTV4_rexFunction, list(funstr,Y,metadata,prior,mcmc,extra),   212345)   		   
	}		   
 }
 
 invisible(return(ANS))
    
}

meanfilter <- function(x,n=5){filter(x,rep(1,n), sides=2)}

 
 
 
