#get.devfun   = function(dev) {
#    dev <- if (exists(dev, asNamespace("grDevices")))  
#      get(dev, asNamespace("grDevices"))
#    else if (exists(dev, .GlobalEnv)) 
#      get(dev, .GlobalEnv)
#	else 
#	  NA
#}
#devfun.exist = function (dev) {   
#   is.function(get.devfun(dev)) 
#}
#new.interactiveWindow = function(width, height) {

#  oldDev =  getOption('device')
#  #identical(oldDev,get.devfun('windows')) ) 
#  #get.devfun('RStudioGD')
  
#  if ( is.function(oldDev) ) {
     
#	  if (.Platform$OS.type == "windows")  {
#	     options(device = get.devfun('windows')) 
#	     dev.new(width=width, height=height)
#	  }  else {
#	     if ( devfun.exist('x11')) {
#		   options(device = get.devfun('x11')) 
#		   dev.new(width=width, height=height,type = "Xlib")  
#		 } else if ( devfun.exist('quartz')) {
# 		    options(device = get.devfun('quartz')) 
#		    dev.new(width=width, height=height)  		 
#		 } else {
#		   stop('Can\'t open an interactive window')
#		 }	     
#	  }   
		
#  }  
  
#   if ( is.character(oldDev) ) {
#      #identical(oldDev,get.devfun('windows')) ) 
#	  if (.Platform$OS.type == "windows")  {
#	     options(device = 'windows') 
#	     dev.new(width=width, height=height)
#	  }  else {
#	     if ( devfun.exist('x11')) {
#		   options(device = 'x11') 
#		   dev.new(width=width, height=height,type = "Xlib")  
#		 } else if ( devfun.exist('quartz')) {
#		    options(device = 'quartz') 
#		    dev.new(width=width, height=height)  		 
#		 } else {
##		   stop('Can\'t open an interactive window')
#		 }	     
#	  }   
		
#  }   
  
#  options(device=oldDev)

#}
##################################################################

tetris <- function(height=25, width=14, speed=0.6) {

  if (! base::interactive()) {
    cat("teris() runs only in the interactive command mode.")
    invisible(return(NULL))
  }
  speed=max(speed,0.05)
  speed=min(speed,2)
   
  
  ####################################################################
  ####################################################################
  plotij.close <- function(I,J, w) {
    I=I-1;
    J=J-1;
    r=0.1*w;
    for (k in 1:length(I)){
      i=I[k]
      j=J[k]
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

  plotij.black <- function(i,j, w)  {
    i=i-1;
    j=j-1;
    r=0.1*w;
    polygon( c(i*w  ,i*w ,(i+1)*w ,(i+1)*w  ), 
             c(ny*w-j*w ,ny*w-(j+1)*w ,ny*w-(j+1)*w  ,ny*w-j*w ),
             col='#111111',border='#111111' )    
  }
  plotij.allgrid <-function(){
    plotij.rect(0,0,nx,ny);  
    for (i in 1:nx) {
      for (j in 1:ny){
        if (img[j,i] >0)  plotij.close(i,j,1)
      }
    }  
  }
  
  ####################################################################
  plotij.close <- function(I,J, w) {
    I=I-1;
    J=J-1;
    for (k in 1:length(I)){
      i=I[k]
      j=J[k]
      r=0.05*w;
      polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
               c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
               col='#2D3025',border=NA )
      r=0.15*w;
      polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
               c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
               col='#9EAD86',border=NA )
      r=0.26*w;
      polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
               c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
               col='#2D3025',border=NA )
    }  
  }


  plotij.black <- function(i,j, w) {
    
    i=i-1
    j=j-1
    r=0.05*w;
    polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
             c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
             col='#889373',border=NA )
    r=0.15*w;
    polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
             c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
             col='#9EAD86',border=NA )
    r=0.26*w;
    polygon( c(i*w +r ,     i*w +r,        (i+1)*w -r,       (i+1)*w -r ), 
             c(ny*w-j*w - r,ny*w-(j+1)*w +r,ny*w-(j+1)*w +r ,ny*w-j*w -r),
             col='#889373',border=NA )
  
  }
  
  plotij.allgrid <-function(){
    for (i in 1:nx) {
      for (j in 1:ny){
        if (img[j,i] >0)  plotij.close(i,j,1)
        else              plotij.black(i,j,1)
      }
    }  
  }
  ####################################################################
  plotij.rect <- function(x0,y0,x1,y1)  {
    rect(x0, y0,x1,y1,col='#111111',border='#777777' )    
  }
  
  pplot.field <-function(nx,ny, yExtra, w){
    
    par(  oma=c(0,0,0,0) )      #the margin between the device oboarard and figure brder
    par(  mar=c(.0,.0,.0,.0) )  #the margin between the figure boarder and plot border 	
    #plot( c(0,nx*w),c(-yExtra,(ny+0)*w),type='n',ann=FALSE , axes=FALSE )
    
    # the outter box
    r =.8
    x0=-r;
    x1=nx+r
    y0=-yExtra-r
    y1=ny+r	  
    #Style "i" (internal) just finds an axis with pretty labels that fits within the original data range.
    plot(c(0,nx*w),c(-yExtra,(ny+0)*w),type='p',ann=FALSE,axes=FALSE,
           xaxs='i',yaxs='i',xlim=c(x0, x1),ylim=c(y0,y1) )
    
    polygon( c(x0, x0, x1, x1 ),  c(y1 ,y0, y0, y1 ),   col='#CCCCCC',border=NA )  
    
    # the field boondary
    r=0.15;
    x0=-r;
    y0=-r
    x1=nx+r
    y1=ny+r
    polygon( c(x0+r, x0 +r, x1 -r, x1 -r ),    c(y1 - r, y0 +r, y0 +r, y1 -r),  col='#BBBBBB',border=NA )   #middle
    polygon( c(x0, x0, x0+r, x0+r, x1-r, x1 ), c(y1 ,y0,  y0+r, y1 -r,y1 -r, y1 ),col='#999999',border=NA ) #upper
    polygon( c(x0 ,x1, x1,  x1-r,x1-r,x0+r ),  c(y0 ,y0, y1, y1-r, y0+r,y0+r),   col='#F8F8F8',border=NA )  #bottom
    
    #plotij.close(1:nx,1:ny, w )  
    #plotij.rect(0,0,nx,ny);  
    
    # The bottom box	
    r  = 0.15;
    x0 = -r
    x1 = nx+r
    y0 =-yExtra
    y1 =-0.5
    polygon( c(x0, x0,   x1,  x1 ),   c(y1 ,y0,   y0, y1 ),   col='#AAAAAA',border='#000000' ) 
  }
  
  
 
  ####################################################################
  pplot.buttons <-function(nx,ny, yExtra, w) {	
  
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
    ####################################################################
    # Buttons
    r  = 0.1
    w1 = min((nx-2*r)*0.25,5);
    w2 = 1.5
    
    sw   = strwidth('Restart');	 
    sh   = strheight('Restart');	 
    cex1 = min((w1-2*r)/sw, w2/sh) *0.8;
    plotbutton(nx-w1,      -0.5-(yExtra-0.5)/2+w2/2, w1, w2,r,'Quit',       cex1)
    plotbutton(nx-w1-r-w1, -0.5-(yExtra-0.5)/2+w2/2, w1, w2,r, buttonLabel, cex1) 
    
  }
  
  pplot.count <- function(nx,ny, yExtra, w, count) {
    count = ifelse(count>999,999,count);
    #bnum=paste('No. of bombs left: ', formatC(sum(img)-sum(right),format='d',digits=3),' '  ) ;
    #bnum= paste('Bombs left:', formatC(count,format='d',digits=3),sep=''  ) ;
    bnum=sprintf("Score: %4d ",count)
    
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
  
  if ( !is.loaded("TetrisSetTimer")){
    stop("tetris runs only on Windows not Linux or Mac.")
	invisible(return(NULL))
  }
  
  #https://www.chegg.com/homework-help/questions-and-answers/python-questions-defaltposofshape-shape-def-poseifshifted-positions-shift-q27747438
  fig     =list()
  fig[[1]]=list( c(2,6,10,14), c(9,10,11,12),c(3,7,11,15),c(5,6,7,8)) #i
  fig[[2]]=list( c(5,6,9,10))  #o
  fig[[3]]=list( c(2,5,6,10), c(5,6,7,10),c(2,6,7,10),c(2,5,6,7))#t
  fig[[4]]=list( c(1,2,6,10), c(5,6,7,9), c(2,6,10,11),c(3,5,6,7))#j
  fig[[5]]=list( c(2,6,9,10), c(5,6,7,11),c(2,3,6,10),c(1,5,6,7))#l
  fig[[6]]=list( c(2,5,6,9),  c(5,6,10,11),c(3,6,7,10),c(1,2,6,7))#s
  fig[[7]]=list( c(1,5,6,10), c(6,7,9,10),c(2,6,7,11),c(2,3,5,6))#z
  
  
  #wall kick: https://tetris.fandom.com/wiki/SRS
  
  #JTSZ
  kick37=list()
  kick37[[1]]=matrix(c(0,0,-1,0,-1,1,0,-2,-1,-2),nrow=2)
  kick37[[2]]=matrix(c(0,0,1,0, 1,-1,0,-2,1, 2),nrow=2)
  kick37[[3]]=matrix(c(0,0,1,0, 1, 1,0,-2,1, -2),nrow=2)
  kick37[[4]]=matrix(c(0,0,-1,0, -1, -1,0,2, -1, 2),nrow=2)
  #I
  kick1=list()
  kick1[[1]]=matrix(c(0,0, -2,0, 1,0, -2,-1, 1,2),nrow=2)
  kick1[[2]]=matrix(c(0,0,-1,0, 2,0, -1,2,2, -1),nrow=2)
  kick1[[3]]=matrix(c(0,0,2,0, -1,0, 2,1,-1, -2),nrow=2)
  kick1[[4]]=matrix(c(0,0,1,0, -2,0, 1,-2, -2, 2),nrow=2)
  
  #Global Variable
  ##################################################################
  ny     = 0
  nx     = 0
  img    =0
  curType=1
  curRot =1
  curI   = round(nx/2);
  curJ   = 1;
  buttonLabel='Start'
  DONE    = 0
  SCORE  = 0
  isRemoveLineMode =FALSE
  ##################################################################
  
  height=ifelse(height<=0,12, height)
  width =ifelse(width<=0, 25, width)
  ny  = height
  nx  = width  
  
  #set the grid size & should be always kept 1 
  w   =   1;
  dpi = get.dpi() 
  
  #use dpi to creat a window with a width of about 600 pixels  
  gridsize = 35
  yExtra   = 3;  
  numCex   = 1.0
  WIDTH    = gridsize*nx;
  HEIGHT   = gridsize*(ny+yExtra);
  if (HEIGHT> 850 || WIDTH>1500) {
    gridsize = 30
    WIDTH    = gridsize*nx;
    HEIGHT   = gridsize*(ny+yExtra);
    numCex   = 30/35
  }  
  if (HEIGHT>850) {  ny <<- floor(850/gridsize);   HEIGHT   = gridsize*(ny+yExtra);}
  if (WIDTH>1500) {  nx <<- floor(1500/gridsize);  WIDTH    = gridsize*nx;         }
  new.interactiveWindow(width=WIDTH/dpi , height=HEIGHT/dpi)
  
  clickRestart = FALSE;
  clickQuit    = FALSE
  

  curEnv=environment();

  #Global Variable
  ##################################################################
  #ny     = 0
  #nx     = 0
  img    =- matrix(as.numeric(runif(nx*ny) > 100 ),nrow=ny);
  curType= 1
  curRot = 1
  curI   = round(nx/2);
  curJ   = 1;
  ##################################################################

  ########################################################################
  # Configure plotting area 
  ##################################################################
  {
    #plot.field(  nx,ny, yExtra, w)
    #plot.buttons(nx,ny,yExtra,w); 
    #plot.count(nx,ny,yExtra,w,SCORE);
  }
  replot.allfigs =function(){
    pplot.field(nx,ny, yExtra, w)
    plotij.allgrid()
    pplot.buttons(nx,ny,yExtra,w);
    pplot.count(nx,ny,yExtra,w,SCORE);
#    for (i in 1:nx) {
#      for (j in 1:ny){
#        if (img[j,i] >0)  plotij.close(i,j,1)
#      }
#    }
  }
  
  replot.allfigs()
  
  
  plotij.number <- function (i,j,w ){
    text( (i-0.5)*w ,(ny-j+0.5)*w, '*', font=2,cex=numCex ,col='red')
  }
  plot.star=function(){
    return(NULL)
    for (i in 1:nx)
      for (j in 1:ny){
        if (img[j,i] >0){
          plotij.number(i,j,w )
        }
      }
  }
  plot.fig=function(i,j, t,rot){
    idx=fig[[t]][[rot]]
    col=ceiling(idx/4)
    row=idx-(col-1)*4
    col=col-2;
    
    col=col+i
    row=row-1+j  #row=(ny+1)-(row+j)
    plotij.close(col,row, 1) 
    #cat('xxxxxxxxxx \n')
    #print(col)
    #print(row)
  }
  remove.fig=function(i,j, t,rot){
    idx=fig[[t]][[rot]]
    col=ceiling(idx/4)
    row=idx-(col-1)*4
    
    col=col-2+i
    row=row-1+j  #row=(ny+1)-(row+j)
    for (k in 1:length(col)){
      plotij.black(col[k],row[k], 1)    
    }
  }

  checkCollision<-function(i,j,rot=curRot){
    idx=fig[[curType]][[rot]]
    col=ceiling(idx/4)
    row=idx-(col-1)*4
    
    col=col-2+i
    row=row-1+j  #row=(ny+1)-(row+j)
    if( sum(col<1)>0 | sum(col>nx)>0){   
      return(1)   # collison with the vertical wall
     } 
    if( sum(row>ny)>0 )             {     
      return(2)  # collison with the top
    } 
    idx=row+(col-1)*ny
    if( sum(img[idx]) >0.001)       {      return(3)    } #collison with other figs
    return(FALSE)
  }
  
  
  freezefig=function(i,j){
    idx=fig[[curType]][[curRot]]
    col=ceiling(idx/4)
    row=idx-(col-1)*4
    col=col-2+i
    row=row-1+j  #row=(ny+1)-(row+j)
    
    idx=row+(col-1)*ny
    if (sum(idx>nx*ny)>0.01 ) { cat("error occured somewhere!")  }
    
    img[idx] <<- 1
  }
  #################################################
  isAnyFullLine = function(){
    rowsum = rowSums(img)
    return (  sum(rowsum == nx) >0.001 )
  }
  remove.oneline=function(ROW) {
    # First, find the bottom black line/top non-full-back line
    figTopLineIdx = which(rowSums(img)!=0)[1]
    # Move the upper part down by one line
    img[2:ROW,]<<-img[1:(ROW-1),]

    #replot the upper part
    for (j in figTopLineIdx:ROW ){
      for (i in 1:nx ) {
        if (img[j,i] >0) plotij.close(i,j,1)
        else             plotij.black(i,j,1)
      } 
    }
    
    SCORE<<-SCORE+1
  }
  remove.full.lines=function( )  {
    #function not used
    anyFullLine=FALSE
    for (i in 2:ny){
      line=img[i,];
      if (sum(line)==nx){ 
        remove.oneline(i);
        anyFullLine=TRUE
        #Sys.sleep(0.3)
      }
    }
    
    if (anyFullLine){
      replot.allfigs()
    }
    
  }
 
  gen.newfig=function(){
    #iotjlsz
    curType <<- ceiling(runif(1)*7)
    curRot  <<- ifelse(curType==2,1, ceiling(runif(1)*4))
    curI    <<- round(nx/2);
    curJ    <<- 1;

    if( checkCollision(curI,curJ)){
      # Game over
      .Call("TetrisSetTimer",0L,speed, curEnv )
      DONE <<-TRUE
      text(nx/2,ny/2.,'Game over!',font=2,cex=3.5)
      text(nx/2*0.99,ny/2*0.99,'Game over!',font=2,cex=3.5,col='#FF0000')
    }
  }
  
  timer = function(prob){
    
    if (DONE){return(NULL)}
    
    if (isRemoveLineMode ){
      idxFullLine = which(rowSums(img) == nx )
      remove.oneline(idxFullLine[1])
      if (length(idxFullLine)-1 ==0){
        # all full lines have been removed
        isRemoveLineMode<<-FALSE
        replot.allfigs()
      }
      return(1)
    }
    
    if( checkCollision(curI,curJ+1) ){
      freezefig(curI, curJ)
      gen.newfig()              # generate a new fig
      isRemoveLineMode <<- isAnyFullLine() # the generated  new fig won't be displayed until all full liines are cleared
      plot.star()
      return(1)
    } else {
      remove.fig(curI, curJ, curType,curRot)
      curJ <<- curJ+1
      plot.fig(  curI, curJ, curType,curRot)
      return(2)
    }
    
  }
  ##################################################################
  mousedown <-function(button,x,y){
    x=grconvertX(x, "ndc", "user")
    y=grconvertY(y, "ndc", "user")
    
    if ( !(x < (nx*w) & x > 0 & y < (ny*w) & y > 0 ) )    {
  
      r  = 0.1
      w1 = min((nx-2*r)*0.25,5);
      w2 = 1.5
      
      x0=nx-w1-r-w1
      y0=-0.5-(yExtra-0.5)/2 +w2/2		
      
      #Hit the restart button
      if ( x > (x0+r) & x < (x0+w1-r) & y < (y0-r) & y > (y0-w2+r)  )   {   
        buttonLabel <<- 'Restart' # before the first run, buttonLabel is 'start'
        .Call("TetrisSetTimer",1L,speed, curEnv )
        DONE  <<- FALSE
        SCORE <<- 0
        img   <<- matrix(as.numeric(runif(nx*ny) > 100 ),nrow=ny);
        #img[20:ny,]<<-1
        #img[, 2]<<-0
        curType <<- ceiling(runif(1)*7)
        curRot  <<- ifelse(curType==2,1, ceiling(runif(1)*4))
        isRemoveLineMode <<- FALSE
        replot.allfigs()
      }

      x0=nx-w1
      y0= -0.5-(yExtra-0.5)/2 +w2/2
      #hit the quit button
      if ( x > (x0+r) & x < (x0+w1-r) & y < (y0-r) & y > (y0-w2+r)  )     {
        .Call("TetrisSetTimer",0L,speed, curEnv )
        dev.off()
      }
      
    }
    #curI      <<- ceiling(x/w);
    #curJ      <<- ceiling((ny*w-y)/w);    
    #curButton <<-button;
    return(NULL)
  }
  #$#$$$$$$$$$$$$$$$$$$$$$$$$$$################
  mouseup <-function(button,x,y){
    x=grconvertX(x, "ndc", "user")
    y=grconvertY(y, "ndc", "user")
    invisible(return(NULL))
  }
  keyleft <- function(key) {
    if (DONE || isRemoveLineMode ){return(NULL)}
    
    if ( checkCollision(curI-1,curJ,curRot) ){
      plot.star()
      return(NULL)
    }
    remove.fig(curI, curJ,curType, curRot)
    curI <<- curI-1
    plot.fig(curI, curJ, curType,curRot);
  }
  keyright <- function(key) {
    if (DONE || isRemoveLineMode ){return(NULL)}
    if ( checkCollision(curI+1,curJ,curRot) ){
      plot.star()
      return(NULL)
    }
    remove.fig(curI, curJ,curType, curRot)
    curI<<-curI+1
    plot.fig(curI, curJ, curType,curRot);
  }
  keydwon <- function(key) {
    if (DONE || isRemoveLineMode ){return(NULL)}
    if ( checkCollision(curI,curJ+1,curRot) ){
      plot.star()
      return(NULL)
    }
    remove.fig(curI, curJ,curType, curRot)
    curJ<<-curJ+1
    plot.fig(curI, curJ, curType,curRot);
    NULL
  }  
  keyup <- function(key) {
    if (DONE || isRemoveLineMode ){return(NULL)}
    if (curType==2){  return(NULL)  }
    if (curType==1)   kik=kick1[[curRot]]
    else              kik=kick37[[curRot]]
  
    newRot=ifelse(curRot==4, 1, curRot+1)
    goodKick=-999
    for (i in 1:5){
      off=kik[ ,i]
      if (!checkCollision(curI+off[1],curJ+off[2],newRot)){
        goodKick=i
        break;
      }
    }
    
    if (goodKick>0){
      remove.fig(curI, curJ, curType, curRot)
      curI  <<- curI+kik[1,goodKick]
      curJ  <<- curJ+kik[2,goodKick]
      curRot<<- newRot
      plot.fig(curI, curJ, curType,curRot);
    }
    NULL
  }  
  
  keyspace <- function(key) {
    if (DONE || isRemoveLineMode ){return(NULL)}
    finalRow=curJ;
    for (j in (curJ+1):ny) {
      if ( checkCollision(curI,j,curRot) ){
        break;
      }
      finalRow=j;
    }
    remove.fig(curI, curJ,curType, curRot)
    curJ <<- finalRow
    plot.fig(curI, curJ, curType,curRot);
    
    freezefig(curI, curJ)
    #remove.full.lines()
    #generate a new fig
    gen.newfig()
    isRemoveLineMode <<- isAnyFullLine() # the generated  new fig won't be displayed until all full liines are cleared
    plot.star()
    
  }  
  
  instructStr ='Tetris in R:\n(1) LEFT arrow: move left\n(2) RIGHT arrow: move rigth \n(3) Up: rotate \n(4) Down: speed up\n(5) Space: sink to the bottom\n'
  cat(instructStr)
  getGraphicsEvent('Rtetris-Rbeast', onMouseDown = mousedown, onMouseUp = mouseup,onKeybd =NULL);  
  invisible(NA)
}

 
