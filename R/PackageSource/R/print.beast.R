print.beast<-function (x,...)
{
#https://stackoverflow.com/questions/50561768/r-get-argument-names-from-function-call
  #argName=formalArgs(deparse(substitute(x))[[1]])
#https://stackoverflow.com/questions/4959463/getting-the-object-name-for-s3-print-method-failing
 #argName=deparse(substitute(x))
 
 #argName=strtrim(argName)
  if (is.null(x))
  {
     return( invisible(NULL))
  }
  
  nTS= length(x$tN);
  cat(paste("\nDecomposition results for a total of", sprintf("%4d",nTS)," time series: \n"))
  if (attributes(x)$algorithm=='beastTrend')
  {
	  for (i in 1:nTS)
	  {
		s=paste("Time series #",as.character(i),": ",  
				sprintf("%4.1f",x$sN[i]) , "No seasonal component&",
				sprintf("%4.1f",x$tN[i]) , "trend changepoints \n"          );
		cat(s)
	}
  }
  else
  {
	  for (i in 1:nTS)
	  {
		s=paste("    Time series #",as.character(i),": ",  
				sprintf("%4.1f",x$sN[i]) , " seasonal changepoints &",
				sprintf("%4.1f",x$tN[i]) , " trend changepoints \n"          );
		cat(s)
	  }
  }
  
  cat("\n\n======\n")
  s=paste("The beast output variable, say 'y', is a LIST object. Type names(y)", 
  "to see a list of elements in y. The current 'y' contains the following elements:\n")
  cat(s)
  namesList=names(x);
  namesList[1]=paste("  ", namesList[1]);
  cat(namesList,sep="\n  ")
  cat("\nCheck individual elements to see the model outputs (e.g, Type 'y$t' to see the fitted trend). ")  
  cat("Or type 'plot(y)' to plot the model output.\n")  
  
}
