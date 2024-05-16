tsextract = function(x, index=1){
 
  if ( sum(class(x) == 'beast') == 0 ){
    stop("The input x must have a 'beast' class type (i.e., an output from the BEAST function).")
  }
  
  if ( length(x$marg_lik)==1 ){
  # only one time series is present
    return(invisible(x))   #invisible will return back
  }  
  
  rows=x$nrows
  cols=x$ncols
  nTS=rows*cols
  if( length(index)==1){
     if(index>nTS){	 
		stop(sprintf("Invalid index: The input object conatains a total of %d time series.",nTS))
	 }
  } else {
	if(index[1]<1 || index[1]>rows || index[2]<1 || index[2]>cols){	 
		stop(sprintf("Invalid indices of [%d, %d]: The valid row and col dimensions of the input object are %d and %d.",index[1],index[2],rows, cols))
	}
  }  
  
  o = .Call(BEASTV4_rexFunction,list('tsextract',x,index),12345);
  invisible(o)
}
 
 "[.beast" = function(x, index1=1, index2=NULL, ...){
 
        if (is.character(index1)){
		   return(x[index1])
		} else {
		   if (is.null(index2)) return(tsextract(x,index1))
		   else		        return(tsextract( x, c(index1,index2 ) ) )		
		}
		
 }