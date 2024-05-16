get.names  = function(o){  
  TAG    = list()
  TAGSTR = c()
    
  i      = 1
  NAMES  = names(o)
  for (name in NAMES) {
    e = o[[name]]
    if (is.list(e)){
      subNAMES = names(e)
      for (subname in subNAMES){
        #e1        = e[[subname]]     
        TAG[[i]]   = c(name,subname)
        TAGSTR[i]  = paste(name,subname,sep='$')   
        i  = i+1		
      }
    } else {      
      TAG[[i]]  = name
      TAGSTR[i] = name  
      i  = i+1	  
    }	
	
  }
  
  len = nchar(TAGSTR)
  idx = which(len==max(len))[[1]]
  x   = list(tag=TAG, tagstr=TAGSTR, maxstr=TAGSTR[idx])  
  invisible(x)
}

