 
# When loaded, the namespace asNamespace('Rbeast') of the package is locked.
# The global varaiables in the package can't be modified using <<-. So
# Solution: 
#1.unlockBinding('dpi',envir=asNamespace('Rbeast'))
#2.Create a package environment
#https://stackoverflow.com/questions/34254716/how-to-define-hidden-global-variables-inside-r-packages
#https://stackoverflow.com/questions/34254716/how-to-define-hidden-global-variables-inside-r-packages
#https://stackoverflow.com/questions/28246952/global-variable-in-a-package-which-approach-is-more-recommended

pkgEnv     = new.env(parent=emptyenv())
pkgEnv$dpi = NA

#get('dpi', envir=pkgEnv)
#assign('dpi',22,envir=pkgEnv)


##################################################################
get.devfun   = function(dev) {
    dev <- if (exists(dev, asNamespace("grDevices")))  
      get(dev, asNamespace("grDevices"))
    else if (exists(dev, .GlobalEnv)) 
      get(dev, .GlobalEnv)
	else 
	  NA
}
devfun.exist = function (dev) {   
   is.function(get.devfun(dev)) 
}
new.interactiveWindow = function(width, height) {

  oldDev =  getOption('device')
  #identical(oldDev,get.devfun('windows')) ) 
  #get.devfun('RStudioGD')
  
  if ( is.function(oldDev) ) {
     
	  if (.Platform$OS.type == "windows")  {
	     options(device = get.devfun('windows')) 
	     dev.new(width=width, height=height)
	  }  else {
	     if ( devfun.exist('x11')) {
		   options(device = get.devfun('x11')) 
		   dev.new(width=width, height=height,type = "Xlib")  
		 } else if ( devfun.exist('quartz')) {
		    options(device = get.devfun('quartz')) 
		    dev.new(width=width, height=height)  		 
		 } else {
		   stop('Can\'t open an interactive window')
		 }	     
	  }   
		
  }  
  
   if ( is.character(oldDev) ) {
      #identical(oldDev,get.devfun('windows')) ) 
	  if (.Platform$OS.type == "windows")  {
	     #options(device = 'windows') 
		  # Changed to the function bcz some package like IRanges
		  # redefines "windows" and dev.new won't work if specified with
		  # the string "windows".
		 options(device = get.devfun('windows')) 
	     dev.new(width=width, height=height)
	  }  else {
	     if ( devfun.exist('x11')) {
		   options(device = 'x11') 
		   # X11 font can't be loaded on Ubuntu Linux for Xlib
		     if(base::capabilities('cairo')) {
			    dev.new(width=width, height=height,type = "cairo")  
			 } else {
			 	dev.new(width=width, height=height,type = "Xlib")  
			 }			
		 } else if ( devfun.exist('quartz')) {
		    options(device = 'quartz') 
		    dev.new(width=width, height=height)  		 
		 } else {
		   stop('Can\'t open an interactive window')
		 }	     
	  }   
		
  }   
  
  options(device=oldDev)

}
##################################################################
get.dpi = function(){
# set the global vaiabe dpi
  dpi = pkgEnv$dpi
  
  if( !is.na(dpi) ) {
	return (dpi)
  }

  #create a tempeporay window to cacluate the dot per inch
  #dev.new(  width=nx/2, height=ny/2)
  #win.graph(width=nx/3, height=(7+ny)/3)
  #x11(10,10)
  #dev.new(10, 10)
  {
    #sysName=Sys.info()['sysname']
    #if ( toupper(substr(sysName,1,3)) == 'LIN' )
    #  x11(width=650/dpi, height=650/dpi/nx*(8+ny))
    #else if ( toupper(substr(sysName,1,3)) == 'DAR')
    #  quartz(width=650/dpi, height=650/dpi/nx*(8+ny))
    #else
    #  windows(width=650/dpi, height=650/dpi/nx*(8+ny))
  }
  
   # R CRAN check problem:
   #Found an obsolete/platform-specific call in the following functions:
   #  'get.dpi' 'minesweeper' 'plot.interactive'
   #Found the platform-specific device:
   #  'x11'   dev.new() is the preferred way to open a new device, in the unlikely
   
   
   #if (.Platform$OS.type == "windows")        X11(width=2, height=2)
   #else                                       X11(width=2, height=2,type = "Xlib")  
  
   #Possible solution:
   #  https://stat.ethz.ch/pipermail/r-devel/2013-October/067812.html
   # https://github.com/statnet/ndtv/blob/master/ndtv/R/export.movie.R
   # https://github.com/cran/dynatopmodel/blob/0d7d074b022452f925c8de4271d4a843882ac84c/R/disp_util.r
   
  new.interactiveWindow(width=5,height=5)
  
  wsize = dev.size('px');
  dpi   =  wsize[1]/5
  pkgEnv$dpi = dpi
  dev.off()
  
  return (dpi)
}
 
####################################################################





#prepare the curIrcle pattern for the bomb icon
circleX = cos( seq(1,50)/50 *2*3.1415935  );
circleY = sin( seq(1,50)/50 *2*3.1415935 );

####################################################################
plotij.close <- function(nx,ny,I,J, w)  {
  
  I=I-1;
  J=J-1;
  r=0.1*w;
   
  for (i in I){
		for (j in J){
		    polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
					c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
					col='#BBBBBB',border=NA )
			polygon( c(i*w ,i*w,i*w+r, i*w+r, (i+1)*w-r,(i+1)*w ), 
					c(ny*w-j*w ,ny*w-(j+1)*w, ny*w-(j+1)*w +r ,ny*w-j*w -r,ny*w-j*w -r,ny*w-j*w ),
					col='#F8F8F8',border=NA )   
			polygon( c(i*w ,(i+1)*w,(i+1)*w, (i+1)*w-r, (i+1)*w-r,i*w+r ), 
					c(ny*w-(j+1)*w ,ny*w-(j+1)*w, ny*w-j*w  ,ny*w-j*w-r, ny*w-(j+1)*w+r,ny*w-(j+1)*w +r),
					col='#999999',border=NA )  
		}
	}  
}

plotij.open <- function(nx,ny,i,j, w)  {
  i=i-1;
  j=j-1;
  r=0.1*w;
  polygon( c(i*w  ,i*w ,(i+1)*w ,(i+1)*w  ), 
             c(ny*w-j*w ,ny*w-(j+1)*w ,ny*w-(j+1)*w  ,ny*w-j*w ),
             col='#BBBBBB',border='#777777' )    
}
plotij.bomb <- function(nx,ny,i,j, w)  {
  i=i-1;
  j=j-1;
  r=0.1*w;
  polygon( c(i*w  ,i*w ,(i+1)*w ,(i+1)*w  ), 
           c(ny*w-j*w ,ny*w-(j+1)*w ,ny*w-(j+1)*w  ,ny*w-j*w ),
           col='#BBBBBB',border='#777777' )    
    
  polygon( ((i+1)-0.5)*w  + circleX * w/4  , (ny-(j+1)+0.5)*w  + circleY * w/4, col='black')
  lines(  c( ((i+1)-0.5)*w,     ((i+1)-0.5)*w+w/3) ,         c( (ny-(j+1)+0.5)*w ,      (ny-(j+1)+0.5)*w + w/3 ))
  lines(  c( ((i+1)-0.5)*w+w/3, ((i+1)-0.5)*w+w/3 +w/12  ) , c( (ny-(j+1)+0.5)*w + w/3, (ny-(j+1)+0.5)*w + w/3 -w/12 ) )    
}

plotij.flag <- function(nx,ny,i,j, w)  {
  i=i-1;
  j=j-1;
  r=0.1*w;
    polygon(c(i*w +r ,i*w +r,(i+1)*w -r,(i+1)*w -r ), 
            c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
            col='#BBBBBB',border=NA )     
    polygon(c(i*w ,i*w,i*w+r, i*w+r, (i+1)*w-r,(i+1)*w ), 
            c(ny*w-j*w ,ny*w-(j+1)*w, ny*w-(j+1)*w +r ,ny*w-j*w -r,ny*w-j*w -r,ny*w-j*w ),
            col='#F8F8F8',border=NA )   
    polygon(c(i*w ,(i+1)*w,(i+1)*w, (i+1)*w-r, (i+1)*w-r,i*w+r ), 
            c(ny*w-(j+1)*w ,ny*w-(j+1)*w, ny*w-j*w  ,ny*w-j*w-r, ny*w-(j+1)*w+r,ny*w-(j+1)*w +r),
            col='#999999',border=NA )       
    
    xc = (i+0.5)*w - w/6;
    yc = (ny-(j+1)+0.5)*w;
    a  = w/2;
    b  = w/9;
    c  = w/4;
    polygon( c( xc  ,xc, xc+a, xc+a  ), c( yc+b+c, yc+b,yc+b,yc+b+c ), col='#FF0000',border=NA )   
    
    d=w/4;
    lines( c( xc,xc ) , c(yc+b+c,yc-d) )
    e  = w/3;
    yc = yc-d
    f  = w/7;
    xc = (i+0.5)*w 
    polygon( c( xc-e  ,xc -e, xc+e, xc+e  ),  c( yc, yc-f,yc-f,yc ),   col='#000000',border=NA )    
}
plotij.bombing <-function(nx,ny,i,j, w)  {
  i=i-1;
  j=j-1;
  r=0.1*w;
 
 polygon( c(i*w  ,i*w ,(i+1)*w ,(i+1)*w  ), 
          c(ny*w-j*w ,ny*w-(j+1)*w ,ny*w-(j+1)*w  ,ny*w-j*w ),
           col='red',border='#777777' )    
  polygon( ((i+1)-0.5)*w  + circleX * w/4  , (ny-(j+1)+0.5)*w  + circleY * w/4, col='black')
  lines(  c( ((i+1)-0.5)*w, ((i+1)-0.5)*w+w/3) , c( (ny-(j+1)+0.5)*w ,(ny-(j+1)+0.5)*w + w/3 ))
  lines(  c( ((i+1)-0.5)*w+w/3, ((i+1)-0.5)*w+w/3 +w/12 ) , c( (ny-(j+1)+0.5)*w + w/3,(ny-(j+1)+0.5)*w + w/3 -w/12 ) )
}

plotij.down <-function(nx,ny,i,j, w)  {  
	i=i-1;
	j=j-1;
	r=0.1*w;
 
    polygon( c(i*w  ,i*w ,(i+1)*w ,(i+1)*w  ), 
             c(ny*w-j*w ,ny*w-(j+1)*w ,ny*w-(j+1)*w  ,ny*w-j*w ),
             col='#BBBBBB',border='#AAAAAA' )      
}
####################################################################
plotbutton <- function (x,y,w1,w2,r,str, cex) {  
  polygon( c(x +r , x +r,    x+w1 -r, x+w1 -r ), 
           c(y - r, y-w2 +r, y-w2 +r ,y -r),
           col='#BBBBBB',border=NA )
  polygon( c(x ,x,    x+r,    x+r,  x+w1-r,    x+w1 ), 
           c(y ,y-w2, y-w2+r ,y-r,  y-r,       y  ),
          col='#F0F0F0',border=NA )   
  polygon(c(x ,   x+w1, x+w1, x+w1-r, x+w1-r,x+r ), 
          c(y-w2 ,y-w2, y  ,  y -r,   y-w2+r,y-w2 +r),
          col='#999999',border=NA )  
  text(x+w1/2-0.01,y-w2/2-0.01,str,font=1,cex=cex,col='#000000')
  text(x+w1/2,     y-w2/2,     str,font=1,cex=cex,col='#FF0000')
}

pplot.field <-function(nx,ny, yExtra, w){
 
    par(  oma=c(0,0,0,0) )      #the margin between the device oboarard and figure brder
    par(  mar=c(.0,.0,.0,.0) )  #the margin between the figure boarder and plot border 	
    #plot( c(0,nx*w),c(-yExtra,(ny+0)*w),type='n',ann=FALSE , axes=FALSE )
	
	#legend(0, -1,instructStr ,cex=1.2 ,text.col="blue", box.col="red", bg="yellow")   
    #text(0,   -4,            instructStr,font=1,adj=c(0,0),cex=1.1)    
	
	# the outter box
    r =.8
	x0=-r;
	x1=nx+r
	y0=-yExtra-r
	y1=ny+r	   	 
	plot(c(0,nx*w),c(-yExtra,(ny+0)*w),type='p',ann=FALSE,axes=FALSE,
			xaxs='i',yaxs='i',xlim=c(x0, x1),ylim=c(y0,y1) )

    polygon( c(x0, x0, x1, x1 ),  c(y1 ,y0, y0, y1 ),   col='#CCCCCC',border=NA )  
	
   # the field boondary
	r=0.15;
	x0=-r;
	y0=-r
	x1=nx+r
	y1=ny+r
	polygon( c(x0+r, x0 +r, x1 -r, x1 -r ),    c(y1 - r, y0 +r, y0 +r, y1 -r),  col='#BBBBBB',border=NA )    #middle
    polygon( c(x0, x0, x0+r, x0+r, x1-r, x1 ), c(y1 ,y0,  y0+r, y1 -r,y1 -r, y1 ),col='#999999',border=NA )   #upper
    polygon( c(x0 ,x1, x1,  x1-r,x1-r,x0+r ),  c(y0 ,y0, y1, y1-r, y0+r,y0+r),   col='#F8F8F8',border=NA )     #bottom

	plotij.close(nx,ny,1:nx,1:ny, w )  
			
   # The bottom box	
	r  = 0.15;
    x0 = -r
	x1 = nx+r
	y0 =-yExtra
	y1 =-0.5
    polygon( c(x0, x0,   x1,  x1 ),   c(y1 ,y0,   y0, y1 ),   col='#AAAAAA',border='#000000' ) 
	
}

# plot.count will be considerd as a class generic method
pplot.count <- function(nx,ny, yExtra, w, count) {
	count = ifelse(count>999,999,count);
	#bnum=paste('No. of bombs left: ', formatC(sum(img)-sum(right),format='d',digits=3),' '  ) ;
	#bnum= paste('Bombs left:', formatC(count,format='d',digits=3),sep=''  ) ;
	bnum=sprintf("Bombs left: %4d ",count)
	
	# The counter and the surrouding box
	sw  = strwidth( bnum);	
	sh  = strheight(bnum);	
	cex = min(nx*0.6/sw, (yExtra-0.5)/sh )
	cex = cex*.8	
	sh  = sh*cex;
	sw  = sw*cex;
	x0  = 0
	x1 = sw
	del =((yExtra-0.5)-sh)/2
	y0 = -yExtra+del
	y1 = -0.5 -del 
 
    polygon( c(x0, x0,   x1,  x1 ),   c(y1 ,y0,   y0, y1 ),         col='#BBBBBB',border='#CCCCCC' )  	
	text(   0, -0.5-(yExtra-0.5)/2, bnum, font=1, adj=c(0,0.5),cex=cex)   
	#legend(.5, -3.5,bnum,xjust=0, yjust=0.5, cex=1.2,text.col="blue", box.col="red",bg="yellow")   
}

pplot.buttons <-function(nx,ny, yExtra, w) {	
	# Buttons
	r  = 0.1
	w1 = min((nx-2*r)*0.25,5);
	w2 = 1.5
 
	sw   = strwidth('Restart');	 
	sh   = strheight('Restart');	 
	cex1 = min((w1-2*r)/sw, w2/sh) *0.8;
	plotbutton(nx-w1,      -0.5-(yExtra-0.5)/2+w2/2, w1, w2,r,'Quit',cex1)
	plotbutton(nx-w1-r-w1, -0.5-(yExtra-0.5)/2+w2/2, w1, w2,r,'Restart',cex1) 
}

minesweeper <- function(height=15, width=12, prob=0.1) {
  
  if (! base::interactive()) {
    cat("minesweeper() runs only in the interactive command mode.")
    invisible(return(NULL))
  }
  
  #print(parent.env( environment()))
  #print(parent.env(parent.env( environment())))
  #print(parent.env(a))
  
  height=ifelse(height<=0,12, height)
  width =ifelse(width<=0, 25, width)
  prob  =ifelse(prob<=0, 0.1, prob)
  prob  =ifelse(prob>1,  0.7, prob)
  ny  <- height
  nx  <- width  
  prob =  1-prob;
  
  
  #set the grid size & should be always kept 1 
  w   =   1;
  dpi <- get.dpi() 
  
  #use dpi to creat a window with a width of about 600 pixels  
  gridsize = 25
  yExtra   = 3;  
  numCex   = 1.0
  WIDTH    = gridsize*nx;
  HEIGHT   = gridsize*(ny+yExtra);
  
  if (HEIGHT> 850 || WIDTH>1500) {
	  gridsize = 20
	  WIDTH    = gridsize*nx;
	  HEIGHT   = gridsize*(ny+yExtra);
	  numCex   = 20/25
  }  
  if (HEIGHT>850) {
    ny        <<- floor(850/gridsize)
	HEIGHT    = gridsize*(ny+yExtra);
  }
  if (WIDTH>1500) {
    nx       <<- floor(1500/gridsize)
	WIDTH    = gridsize*nx;
  }
 
  #if (.Platform$OS.type == "windows")     x11(width=WIDTH/dpi , height=HEIGHT/dpi           )
  #else                                    x11(width=WIDTH/dpi , height=HEIGHT/dpi, type = "Xlib")
  new.interactiveWindow(width=WIDTH/dpi , height=HEIGHT/dpi)
   
  #Initialize variables to store status of the grid 
  img = matrix(as.numeric(runif(nx*ny) > prob ),nrow=ny);
  num = img*0;   # the number of bombs in 3x3 window
  #couting the number of bomb for each spot 
  for (i in 1:ny)
    for (j in 1:nx)
    { num[i,j]= sum( img[ max((i-1),1):min((i+1),ny)      ,  max((j-1),1):min((j+1),nx)     ]) }
  
  left         = img*0; # status indicating whetehr a spot is cleared or not 
  right        = img*0; # status indicating whether a spot is flagged or not 
  process0     = img*0; # status indicating whether a spot is cleared or not when wheel-clicking
  taskFailure  = 0;
  taskFinished = 0;
 
  curI      = 0;
  curJ      = 0;
  curButton = 0;
  
  clickRestart = FALSE;
  clickQuit    = FALSE
 
 ########################################################################
 # Configure plotting area 
 ##################################################################
 pplot.field(nx,ny, yExtra, w)
 count=sum(img)-sum(right)
 pplot.count(nx,ny, yExtra, w,  count)
 pplot.buttons(nx,ny,yExtra,w);
 
 #######################################################################
 numCol=c('blue','green','red','brown','gold','darkviolet','purple','pink' )
 plotij.number <- function (i,j,w ){
    idx=num[j,i]
    if (idx> 0)   {
      text( (i-0.5)*w ,(ny-j+0.5)*w, toString(idx),col=numCol[idx], font=2,cex=numCex )
    }
}
  
 ##################################################################
  mousedown <-function(button,x,y){
    
    x=grconvertX(x, "ndc", "user")
    y=grconvertY(y, "ndc", "user")
    #points(x, y)  
    
    if ( !(x < (nx*w) & x > 0 & y < (ny*w) & y > 0 ) )    {
	
		r  = 0.1
		w1 = min((nx-2*r)*0.25,5);
		w2 = 1.5
		
		x0=nx-w1-r-w1
		y0=-0.5-(yExtra-0.5)/2 +w2/2		
 
      #hit the restart button
      if ( x > (x0+r) & x < (x0+w1-r) & y < (y0-r) & y > (y0-w2+r)  )
      {           
        clickRestart  <<- TRUE
        xloc = x0
        yloc = y0            
        polygon(c(xloc + r, xloc +r,   xloc+w1 -r,xloc+w1 -r ), 
                c(yloc - r, yloc-w2 +r,yloc-w2 +r ,yloc -r),
                col='#AA9999',border=NA )
        text(xloc+w1/2-0.01,yloc-w2/2-0.01,'Restart',font=1,cex=.7,col='#000000')
        text(xloc+w1/2,     yloc-w2/2,      'Restart',font=1,cex=.7,col='#FF0000')
        
      }
      
      
 
	 x0=nx-w1
	 y0= -0.5-(yExtra-0.5)/2 +w2/2
      #hit the quit button
      if ( x > (x0+r) & x < (x0+w1-r) & y < (y0-r) & y > (y0-w2+r)  )
      {
        clickQuit<<-TRUE
        xloc=x0
        yloc=y0            
        polygon(c(xloc +r ,xloc +r,xloc+w1 -r,xloc+w1 -r ), 
                c(yloc - r,yloc-w2 +r,yloc-w2 +r ,yloc -r),
                col='#AA9999',border=NA )
        text(xloc+w1/2-0.01,yloc-w2/2-0.01,'Quit',font=1,cex=.7,col='#000000')
        text(xloc+w1/2,yloc-w2/2,'Quit',font=1,cex=.7,col='#FF0000')
        
      }
      

      invisible(return(NULL))
    }
    
    if (taskFailure ==1 | taskFinished==1)   {
      invisible(return(NULL))
    }
    
    
    curI      <<- ceiling(x/w);
    curJ      <<- ceiling((ny*w-y)/w);    
    curButton <<-button;
    if (left[curJ,curI] == 0 & button[1] !=1)     { 
      plotij.down(nx,ny,curI,curJ,w)
    }
    
    return(NULL)
  }
  
  
  #$#$$$$$$$$$$$$$$$$$$$$$$$$$$################
  mouseup <-function(button,x,y){
    
    x=grconvertX(x, "ndc", "user")
    y=grconvertY(y, "ndc", "user")
    button=curButton;
    #points(x, y)
    
	if (clickQuit)   {
      dev.off()
    }
   
    #if the Restart button has been clicked
    if (clickRestart)     {
      
      #reset the clickResart status variable
      clickRestart <<- FALSE
      
      #Initialize variables to store status of the grid 
      img <<- matrix(as.numeric(runif(nx*ny) > prob ),nrow=ny);
      num <<- img*0;   # the number of bombs in 3x3 window
      #couting the number of bomb for each spot 
      for (i in 1:ny)
        for (j in 1:nx)
        { num[i,j]<<- sum( img[ max((i-1),1):min((i+1),ny)      ,  max((j-1),1):min((j+1),nx)     ]) }
      
      left        <<- img*0;  # status indicating whetehr a spot is cleared or not 
      right       <<- img*0; # status indicating whether a spot is flagged or not 
      process0    <<- img*0; # status indicating whether a spot is cleared or not when wheel-clicking
      taskFailure <<- 0;
      taskFinished<<- 0;
      
      curI<<-0;
      curJ<<-0;
      curButton<<-0;      
  
     #Once a plot is created, you can add to the plot, but nothing can be removed. You need to redraw the plot without the legend.	 #https://stackoverflow.com/questions/7365464/remove-legend-in-r
     pplot.field(nx,ny, yExtra, w)
	 count=sum(img)-sum(right)
	 pplot.count(nx,ny, yExtra, w,  count)
	 pplot.buttons(nx,ny,yExtra,w);
  
      return(NULL)
    }
    
    #if click the region out of spot canvas
    if ( !(x < (nx*w) & x > 0 & y < (ny*w) & y > 0 ) )    {      
      invisible(return(NULL))
    }
    
    if (taskFailure ==1 | taskFinished==1) {
      invisible(return(NULL))
    }
    
    i= ceiling(x/w);
    j= ceiling((ny*w-y)/w);
    
    if (i!=curI |j !=curJ ){
      if (left[curJ,curI] ==0)  { 
        
        plotij.close(nx,ny, curI,curJ,w)
        if (right[curJ,curI] ==0)
           plotij.close(nx,ny, curI,curJ,w)
        else if (right[curJ,curI] ==1)
           plotij.flag(nx,ny, curI,curJ,w)
        
      }
      invisible(return(NULL))
    }
    
    if (button[1]==0) {
      #####################
      if (right[j,i]==1 ) {
        plotij.flag(nx,ny,i,j,w)
        invisible(return(NULL))
      }
      
      if ( left[j,i]==0 ) { # a never-touched grid
        if (img[j,i]==0)    # a safe try
        {
          plotij.open( nx,ny,i,j,w); 
          plotij.number(i,j,w)
          left[j,i] <<-1;
          
          if (num[j,i] ==0) 
          {process0[j,i] <<-1;
          hitzero(i,j)}
          
        }
        else # a bad try
        {
          idx    = which(img==1);
          rowlist= row(img)[idx];
          collist= col(img)[idx];
          for ( n in 1:length(idx) )
          {      
            #cat('aa ', rowlist[n],collist[n], '\n')
            plotij.bomb(nx,ny,collist[n],rowlist[n],w)
          }
          plotij.bombing(nx,ny,i,j,w)
          
          taskFailure <<- 1;
          
          text(nx/2,ny/2.,'Game over!',font=2,cex=3.5)
          text(nx/2*0.99,ny/2*0.99,'Game over!',font=2,cex=3.5,col='#FF0000')
        }
      }
      #####################    
    }
    else if (button[1]==2)
    {
      
      if (left[j,i]==0) # not revealed yet
      {
        if (right[j,i]==0) # not flagged yet
        {
          plotij.flag(nx,ny,i,j,w)
          right[j,i]<<-1;
          ###GAME Passed
          if( all(img==right))
          {
            taskFinished<<-1;
            text(nx/2,ny/2.,'Congratulations! \n You win!',font=2,cex=3.5)
            text(nx/2*0.99,ny/2*0.99,'Congratulations! \n You win!',font=2,cex=3.5,col='#FF0000')
          }
		
			count=sum(img)-sum(right)
			pplot.count(nx,ny, yExtra, w,  count)    
	
        }
        else # already flaged
        {          
          plotij.close(nx,ny,i,j,w)
          right[j,i]<<-0;
          
		  count=sum(img)-sum(right)
		  pplot.count(nx,ny, yExtra, w,  count)             
        }
        
        
      }
      
    }
    else if (button[1]==1)
    {
      
      if (left[j,i]==1) #  revealed yet
      {        
        hitflag(i,j)
      }
      
      
    }
    
    invisible(return(NULL))
  }
  
  
  ####################################
  hitzero <- function (i,j){    
    for(M in -1:1)
    {
      for (N in -1:1){
        
        curJ=(j+M);
        curI=(i+N);
        if (M ==0 & N==0) {next}
        if (curJ < 1 | curJ > ny | curI <1 | curI > nx) {next}
        
        if (right[curJ,curI] ==0)
        {
          plotij.open(nx,ny,curI,curJ,w) 
		  plotij.number(curI,curJ,w)	      
          left[curJ,curI]<<-1;       
        }
        
        if (num[curJ,curI]==0 &  process0[curJ,curI]==0)         
        {process0[curJ,curI]<<-1;
        hitzero(curI,curJ)}
        
      }
      
    }   
    
  }
  
  ###############################################
  hitflag <- function (i,j){    
    imgsub   = img[ max(j-1,1):min(j+1,ny) ,  max(i-1,1):min(i+1,nx)        ]
    rightsub = right[ max(j-1,1):min(j+1,ny) ,  max(i-1,1):min(i+1,nx)        ]
    
    if (all(imgsub==rightsub) & max(imgsub) == 1)     {
      
      for(M in -1:1)       {
        
        for (N in -1:1){
          
          curJ=(j+M);
          curI=(i+N);
          if (M ==0 & N==0) {next}
          if (curJ < 1 | curJ > ny | curI <1 | curI > nx) {next}
          
          if (right[curJ,curI] == 0)           {
            plotij.open(nx,ny,curI,curJ,w);
			plotij.number(curI,curJ,w)			
            left[curJ,curI] <<-1;       
          }
          
        }
        
      }
      
    }
  }
  
  instructStr ='MineSweeper in R:\n(1) LEFT-click to clear\n(2) RIGHT-click to flag \n(3) MIDDEL-click(wheel) a cleared&numbered spot to open neighor spots \n(4) Click Restart for a new game\n'
  getGraphicsEvent(instructStr, onMouseDown = mousedown, onMouseUp = mouseup);  
  invisible(NA)
}
