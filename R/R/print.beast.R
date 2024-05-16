printf = function(...) { 
   cat(sprintf(...))
}
  
print.beast <- function (x,index=1,...) {

  
  if ( length(x$marg_lik)> 1 ){
  # More than one time series is present
	 rows=x$nrows
	 cols=x$ncols
	 nTS=rows*cols
	  
	 call  = match.call()
     oName = as.character(call$x)
  
	 if( length(index)==1){
		 if(index>nTS){	 
		      stop(sprintf("Invalid index: The input object conatains a total of %d time series.",nTS))
	      }
		  printf("\nResult for time series #%d (total number of time series in '%s': %d)\n\n", index, oName, nTS) 		 
     } 
	 else {
	     if(index[1]<1 || index[1]>rows || index[2]<1 || index[2]>cols){	 
		         stop(sprintf("Invalid indices of [%d, %d]: The valid row and col dimensions of the input object are %d and %d.",index[1],index[2],rows, cols))
	     }
		 printf("\nResult for time series #(%d,%d) (total number of time series in '%s': %d=%dx%d)\n\n", index[1],index[2],oName, nTS,rows, cols)		 
     }  
  }  
  

  

 
  
  .Call( BEASTV4_rexFunction, list("print",x,index),   212345)   	
  return(invisible(NULL))
  
  #The following R code is deprecated and has been ported in C. But they are kept here
  #as a reference and they will be never run.
  
  call  = match.call()
  oName = as.character(call$x)
 
  #https://stackoverflow.com/questions/50561768/r-get-argument-names-from-function-call
  #argName=formalArgs(deparse(substitute(x))[[1]])
  #https://stackoverflow.com/questions/4959463/getting-the-object-name-for-s3-print-method-failing
  #argName=deparse(substitute(x))
  
  #argName=strtrim(argName)
  if (is.null(x)) { return( invisible(NULL))  }
  
  xdim  = dim(x$R2)
  n     = xdim[1]
  m     = xdim[2]
  nTS   = n*m
  is3D  = (m>1)
 
  if ( is.null( attributes(x)$tsextract ) ) {
	   x = tsextract(x,index)
  } else {
       is3D = NULL
  } 
 
  if (!is.null(is3D)) {  
	  if (is3D & length(index)>1)
		 printf("\nResult for time series #(%d,%d) (total number of time series in '%s': %d)\n\n", index[1],index[2],oName, nTS)
	  else
		 printf("\nResult for time series #%d (total number of time series in '%s': %d)\n\n", index, oName, nTS) 
  }
     
  
  hasSeason   =!is.null(x$season)
  hasOutlier  =!is.null(x$outlier)
  hasHarmonic =!is.null(x$season$order)
  
  s1    = '                                                ';
  s2    = '************************************************';

  ###################################################################################
  ####            SEASONAL CHANGEPOINTS                                             $##
  ###################################################################################   	
 
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("+                      SEASONAL CHANGEPOINTS                      +\n")
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n")
  
  if (!hasSeason) {
    cat(" No seasonal/periodic component present (i.e., season='none')\n")
  } else {
    y     = x$season
    n     = length(y$ncpPr)
    maxPr = max(y$ncpPr)
    maxIx = which(y$ncpPr==maxPr)[1]
    n     = min(n,99)
    printf('An ascii plot of the probability distribution for number of chgpts\n');
    printf('------------------------------------------------------------------\n');
    for (i in 1:n){
	  p=max(1,ceiling( y$ncpPr[i]/maxPr*(nchar(s1)-1)) )
      s=s1; substr(s,1,p)=s2	  
      printf("Pr(ncp=%-2d)=%.3f|%s|\n",i-1, y$ncpPr[i],s)      
    }
	
    printf('----------Summary for number of seasonal changepoints-------------\n')
    printf('ncp_max    : %-4d | A parameter you set (e.g., maxSeasonKnotNum) |\n',length(y$ncpPr)-1 );
    printf('ncp_mode   : %-4d | Pr(ncp=%2d)=%3.2f; there is a %3.1f%% probability|\n',maxIx-1,min(maxIx-1,99),maxPr,maxPr*100)
    printf('	          | that the seasonal componet has %2d chngpt(s). |\n', min(maxIx-1,99) );
    printf('ncp_mean   : %-4.2f | Sum[ncp*Pr(ncp)]                             |\n',          y$ncp);
    printf('ncp_median : %-4.2f | Median number of changepoints                |\n',          y$ncp_median);
    printf('ncp_pct90  : %-4.2f | 90%% perctile for number of changepoints     |\n',         y$ncp_pct90);
    printf('------------------------------------------------------------------\n');  
    cat('\n')
  
    printf("List of probable seasonal changepoints: Please combine the ncp reported \nabove to determine which CPs below are meaningful.\n" )
    printf("---------------------------------------.\n")
    printf("scp#   |time (cp)        |prob(cpPr)   |\n")
    printf("-------|-----------------|-------------|\n")
	 ncp = sum( !is.na (y$cp) )
    for (i in 1:ncp){
     printf("%-7d|%-17.4f|%13.5f|\n", i,y$cp[i], y$cpPr[i] )  
    }
    printf("---------------------------------------'\n")   
   cat('\n')  
 
    
    
  }  
  cat('\n\n')
  
 
  
  ###################################################################################
  ####            TREND CHANGEPOINTS                                               $##
  ###################################################################################  
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("+                      TREND CHANGEPOINTS                         +\n")
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n")

  y=x$trend
  n=length(y$ncpPr)
  maxPr=max(y$ncpPr)
  maxIx=which(y$ncpPr==maxPr)[1]
  n=min(n,99)
  printf('An ascii plot of the probability distribution for number of chgpts\n');
  printf('------------------------------------------------------------------\n');
  for (i in 1:n){
	  p=max(1,ceiling( y$ncpPr[i]/maxPr*(nchar(s1)-1)) )
      s=s1; substr(s,1,p)=s2	  
      printf("Pr(ncp=%-2d)=%.3f|%s|\n",i-1, y$ncpPr[i],s)   
  }
  printf('----------Summary for number of trend changepoints----------------\n')
  printf('ncp_max    : %-4d | A parameter you set (e.g., maxTrendKnotNum)  |\n',length(y$ncpPr)-1 );
  printf('ncp_mode   : %-4d | Pr(ncp=%2d)=%3.2f; there is a %3.1f%% probability|\n',maxIx-1,min(maxIx-1,99),maxPr,maxPr*100)
  printf('	          | that the trend componet has %2d chngpt(s).    |\n', min(maxIx-1,99) );
  printf('ncp_mean   : %-4.2f | Sum[ncp*Pr(ncp)]                             |\n',          y$ncp);
  printf('ncp_median : %-4.2f | Median number of changepoints                |\n',          y$ncp_median);
  printf('ncp_pct90  : %-4.2f | 90%% perctile for number of changepoints      |\n',         y$ncp_pct90);
  printf('------------------------------------------------------------------\n');  
  cat('\n')
  
  printf("List of probable trend changepoints: Please combine the ncp reported \nabove to determine which CPs below are meaningful.\n" )
  printf("---------------------------------------.\n")
  printf("tcp#   |time (cp)        |prob(cpPr)   |\n")
  printf("-------|-----------------|-------------|\n")
  ncp = sum( !is.na (y$cp) )
  for (i in 1:ncp){
   printf("%-7d|%-17.4f|%13.5f|\n", i,y$cp[i], y$cpPr[i] )  
  }
  printf("---------------------------------------'\n")   
  cat('\n') 
  
   
  ###################################################################################
  ####            OUTLIER CHANGEPOINTS                                             $##
  ###################################################################################     
  if (hasOutlier) {

  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("+                      OUTLIER CHANGEPOINTS                       +\n")
  cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n")  
 
    y     = x$outlier
    n     = length(y$ncpPr)
    maxPr = max(y$ncpPr)
    maxIx = which(y$ncpPr==maxPr)[1]
    n     = min(n,99)
    printf('An ascii plot of the probability distribution for number of chgpts\n');
    printf('---------------------------------------------------------------\n');
    for (i in 1:n){
	  p=max(1,ceiling( y$ncpPr[i]/maxPr*(nchar(s1)-1)) )
      s=s1; substr(s,1,p)=s2	  
      printf("Pr(ncp=%-2d)=%.3f|%s|\n",i-1, y$ncpPr[i],s)      
    }
    printf('----------Summary for number of trend changepoints----------------\n')
    printf('ncp_max    : %-4d | A parameter you set (e.g., maxTrendKnotNum)  |\n',length(y$ncpPr)-1 );
    printf('ncp_mode   : %-4d | Pr(ncp=%2d)=%3.2f; there is a %3.1f%% probability|\n',maxIx-1,min(maxIx-1,99),maxPr,maxPr*100)
    printf('	          | that the trend componet has %2d chngpt(s).    |\n', min(maxIx-1,99) );
    printf('ncp_mean   : %-4.2f | Sum[ncp*Pr(ncp)]                             |\n',          y$ncp);
    printf('ncp_median : %-4.2f | Median number of changepoints                |\n',          y$ncp_median);
    printf('ncp_pct90  : %-4.2f | 90%% perctile for number of changepoints      |\n',         y$ncp_pct90);
    printf('------------------------------------------------------------------\n');  
    cat('\n')
  
    printf("List of probable trend changepoints: Please combine the ncp reported \nabove to determine which CPs below are meaningful.\n" )
    printf("---------------------------------------.\n")
    printf("tcp#   |time (cp)        |prob(cpPr)   |\n")
    printf("-------|-----------------|-------------|\n")
    ncp = sum( !is.na (y$cp) )
   for (i in 1:ncp){
    printf("%-7d|%-17.4f|%13.5f|\n", i,y$cp[i], y$cpPr[i] )  
   }
   printf("---------------------------------------'\n")   
      
 
  }   
  printf("Note: the beast output object '%s' is a LIST. Type 'str(%s)' to see all \nthe elements in it. Or use 'plot(%s)' or 'plot(%s,interactive=TRUE)' to \nplot the model output.\n",oName,oName,oName,oName)

  
}
 