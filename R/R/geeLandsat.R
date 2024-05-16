geeLandsat = function(lon=NA,lat=NA,radius=100, stat='mean',timeout=700){
  
  #lon=-83.229320
  #lat=38.697866;
  #radius=100
  #stat='mean'
  
  oldtimeout = options('timeout')
  
  if(is.na(lon) | is.na(lat) | abs(lon) > 180 | abs(lat) >90)
    stop("lon and lat must be numeric with its typical range.")
  
  radius = abs(radius)
  if(radius>500)
    stop("radius must be less than 500 meters.")
  radius=as.integer(radius)
  
  stats=c('mean','min','max','median')
  if ( sum(stats==tolower(stat))==0 )
    stop("The given spatial aggeration method is invalid. Possile values include mean,min,max, and median.")
 
  
  addr=sprintf("https://zhaokg.pythonanywhere.com/landsat?lat=%f&lon=%f&r=%d&stat=%s", lat,lon,radius,stat)
  print(addr)
  options(timeout=timeout)
  
  #con = url(addr,'rb')
  tic=Sys.time()
  con=NULL
  while (is.null(con)  && as.numeric(Sys.time()-tic) < timeout ){
		tryCatch(   suppressWarnings( {con=url(addr,'rb')} ),
              error=function(e){con=NULL; print("Timeout: try again ...");print(e)} 
			);  
  }

  #open(con)
  #https://stackoverflow.com/questions/26584227/how-can-i-read-float-data-from-a-binary-file-using-r
  nElem=readBin(con,'double',1,size=4)
  
  if ( abs(nElem -as.integer(nElem)) > 0.00001)
    stop("The returned data from the seriver has an incorrect format.")
  
  x=readBin(con,'double',nElem*10,size=4)
  close(con)
  options(timeout=oldtimeout)
  
  if(length(x)!=nElem*10)
    stop("The returned data from the seriver has an incorrect format.")
  
  x=data.frame(matrix(x,nrow=nElem))
  colnames(x)=c('year','mon','day','sensor', 'blue','green','red','nir','swir1','swir2')
  s = x$sensor
  s[x$sensor==4]='LT4'
  s[x$sensor==5]='LT4'
  s[x$sensor==7]='LE7'
  s[x$sensor==8]='LC8'
  x$sensor = s
  x$ndvi=(x$nir-x$red)/(x$nir+x$red)
  x$date=as.Date(paste(x$year,x$mon,x$day,sep = '-'))
  return(x)
 
 
}