beast <- function(data, option=list(),demoGUI=FALSE,...)
{
  if (!hasArg("data") || is.list(data) || length(data)==1)
  {
    warning("Something is wrong with the input 'data'. Make sure the data par is a vector or matrix.")
    return(NULL)
  }

 
  OPTION=NULL;
  if (hasArg("option"))
  {

       if( is.numeric(option)&&(length(option)==1) )
       {
         OPTION=list();
         OPTION$period=option
       }
       else if (is.list(option)) {
         
         OPTION=option
       }
       else{
         OPTION=list()
       }
       
  }
  else
  { 
   
    if (hasArg("PERIOD"))
    {
       tmpList = list(...)
       OPTION  = list()
       OPTION$period=tmpList$period
      }
    else if (hasArg("period"))
    { 
      tmpList = list(...)
      OPTION=list()
      OPTION$period=tmpList$period
     }
    else
    {
		OPTION=list()
    }
     
  }
  
  #if(season)
  #{
  #res=.Call(SARAH_beastST_multipleChain_fast, data, OPTION)
  #}
  #else
  #{
  #res=.Call(SARAH_beastTrend_multipleChain_fast, data, OPTION)
  #}
  
  if (!hasArg("demoGUI"))
  {
  demoGUI=FALSE
  }
  
  if(!demoGUI)
  {

  OPTION$computeChangepoints=1
  ANS=.Call(SARAH_beast2, data, OPTION)
  #cat("\n\n==================================================================\n")
  #s=paste("The beast output variable (e.g., x) is a LIST object. Type names(x)", 
  #"to see a list of elements in x.  \n\nThe current 'x' contains the following #elements: ")
  #cat(s)
  #namesList=names(ANS);
  #cat(namesList)
  #cat(". ")
  #cat("\n\nCheck individual elements to see the model outputs (e.g,type x$t to see #the fitted trend) or plot(x) to draw the model decomposition result.")
 #cat("\n==================================================================")
	
  return(ANS) 
  }
 else
 {
  if(is.loaded("WinMainDemoST") || is.loaded("GUI_beast") )
  {
    OPTION$computeCredible=1
	OPTION$fastCIComputation=1  
    ANS=.Call("GUI_beast", data, OPTION) 
  }
  else
  {
	warning("'demoGUI=TRUE' only works if the current system is Windows x64.")
  }

 }
 
  
    
}

meanfilter <- function(x,n=5){filter(x,rep(1,n), sides=2)}

 
 
 
