#grid.ls(viewports = TRUE)
#current.vpTree()
#pushViewport(D)
#popViewport()

# change to pplot.interactive because the R checker complains about 
#  "Mismatches for apparent methods not registered". Basically, R 
# think plot.interactive is a class generic methods
pplot.interactive =function (o, index=1, ncpStat='mode') {
  
if ( length(o$marg_lik)> 1 ) {
    # more than time series is present
	o=tsextract(o,index)
}  
 
t  = o$time; 
y  = o$data

if (is.null(y)){
  if (is.null(o$season))     y=o$trend$Y
  else                       y=o$trend$Y+o$season$Y
}

t2t=c(t, rev(t))
  

CI = 0
Y  = 0
SD = 0

Prob1 =0;
Prob  =0; 
Order =0
Slp   =0
SlpSD =0
SlpSignPos=0
SlpSignZero=0

cp   =0        
cpCI =0        
ncp  =0    
ncp_mean   = 0      
ncp_mode   = 0      
ncp_median = 0      
ncp_pct90  = 0     
ncp_pct10  = 0    
ncpPr    = 0        
cpPr     = 0
cpChange = 0

":" = function(i,j){
  if(i>j)   return(NULL)
  else      return( seq(i,j))
}

get.Yts    = function(){
  if (is.null(o$season))     Yts=o$trend$Y
  else                       Yts=o$trend$Y+o$season$Y
}
get.T      = function(){

  Y  <<- o$trend$Y;
  SD <<- c( Y-o$trend$SD,  rev(Y+o$trend$SD) ); 
  if ( !is.null(o$trend$CI) )		CI <<-c(o$trend$CI[,1], rev(o$trend$CI[,2])) 
  else                          CI <<-SD

  
  Slp          <<- o$trend$slp
  SlpSD        <<- c(Slp-o$trend$slpSD,  rev(Slp+o$trend$slpSD));         
  SlpSignPos   <<- o$trend$slpSgnPosPr
  SlpSignPos   <<- o$trend$slpSgnZeroPr
  Order        <<- o$trend$order;
  
}
get.tcp    = function(){

  cmpnt = o$trend
  
  cp   <<- cmpnt$cp;         
  cpCI <<- cmpnt$cpCI;  
  ncp  <<- switch(ncpStat, mode=cmpnt$ncp_mode, median=cmpnt$ncp_median,mean=cmpnt$ncp,
                 pct90=cmpnt$ncp_pct90,pct10=cmpnt$ncp_pct10,max=sum(!is.nan(cp))  )
  ncp_mean   <<- cmpnt$ncp
  ncp_mode   <<- cmpnt$ncp_mode
  ncp_median <<- cmpnt$ncp_median
  ncp_pct90  <<- cmpnt$ncp_pct90
  ncp_pct10  <<- cmpnt$ncp_pct10
  
  ncpPr     <<- cmpnt$ncpPr
  cpPr      <<- cmpnt$cpPr
  cpChange  <<- cmpnt$cpAbruptChange
  
  Prob   <<- cmpnt$cpOccPr;    
  Prob1  <<- c(Prob,Prob-Prob)
}
get.pos_tcp= function(){
  cp  <<-o$trend$pos_cp;         
  cpCI <<-o$trend$pos_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$trend$pos_ncpPr
  cpPr <<-o$trend$pos_cpPr
  cpChange<<-o$trend$pos_cpAbruptChange
  
  Prob <<-o$trend$pos_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}
get.neg_tcp= function(){
  cp  <<-o$trend$neg_cp;         
  cpCI <<-o$trend$neg_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$trend$neg_ncpPr
  cpPr <<-o$trend$neg_cpPr
  cpChange<<-o$trend$neg_cpAbruptChange
  
  Prob <<-o$trend$neg_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}
get.inc_tcp = function(){
  cp  <<-o$trend$inc_cp;         
  cpCI <<-o$trend$inc_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$trend$inc_ncpPr
  cpPr <<-o$trend$inc_cpPr
  cpChange<<-o$trend$inc_cpAbruptChange
  
  Prob <<-o$trend$inc_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}

get.dec_tcp=function(){
  cp  <<-o$trend$dec_cp;         
  cpCI <<-o$trend$dec_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$trend$dec_ncpPr
  cpPr <<-o$trend$dec_cpPr
  cpChange<<-o$trend$dec_cpAbruptChange
  
  Prob <<-o$trend$dec_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}


########################################################
#    Functions for seasonal component
########################################################
Amp=0
AmpSD=0
get.S      = function(){

  Y  <<- o$season$Y;
  SD <<- c( Y-o$season$SD,  rev(Y+o$season$SD) ); 
  if ( !is.null(o$season$CI) )	  CI <<-c(o$season$CI[,1], rev(o$season$CI[,2])) 
  else                            CI <<-SD

  Amp    <<-o$season$amp
  AmpSD <<-  c(Amp-o$season$ampSD,  rev(Amp+o$season$ampSD));         
  Order<<-o$season$order;
  
}
get.scp    = function(){
  cmpnt = o$season
  
  cp   <<- cmpnt$cp;         
  cpCI <<- cmpnt$cpCI;  
  ncp  <<- switch(ncpStat, mode=cmpnt$ncp_mode, median=cmpnt$ncp_median,mean=cmpnt$ncp,
                 pct90=cmpnt$ncp_pct90,pct10=cmpnt$ncp_pct10,max=sum(!is.nan(cp))  )
  ncp_mean   <<- cmpnt$ncp
  ncp_mode   <<- cmpnt$ncp_mode
  ncp_median <<- cmpnt$ncp_median
  ncp_pct90  <<- cmpnt$ncp_pct90
  ncp_pct10  <<- cmpnt$ncp_pct10
  
  ncpPr     <<- cmpnt$ncpPr
  cpPr      <<- cmpnt$cpPr
  cpChange  <<- cmpnt$cpAbruptChange
  
  Prob   <<- cmpnt$cpOccPr;    
  Prob1  <<- c(Prob,Prob-Prob)   
  Prob1 <<-c(Prob,Prob-Prob)
}
get.pos_scp= function(){
  cp  <<-o$season$pos_cp;         
  cpCI <<-o$season$pos_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$season$pos_ncpPr
  cpPr <<-o$season$pos_cpPr
  cpChange<<-o$season$pos_cpAbruptChange
  
  Prob <<-o$season$pos_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}
get.neg_scp= function(){
  cp  <<-o$season$neg_cp;         
  cpCI <<-o$season$neg_cpCI;  
  ncp <<-sum(!is.nan(cp))
  ncpPr<<-o$season$neg_ncpPr
  cpPr <<-o$season$neg_cpPr
  cpChange<<-o$season$neg_cpAbruptChange
  
  Prob <<-o$season$neg_cpOccPr;    
  Prob1 <<-c(Prob,Prob-Prob)
}


RED   = function(alpha=1){gpar(col='red',   fill='red',    alpha=alpha)}
BLK   = function(alpha=1){gpar(col='black', fill='black',  alpha=alpha)}
RDASH = function(alpha=1){gpar(col='red',   fill='red', lty=2,   alpha=alpha)}
BDASH = function(alpha=1){gpar(col='black', fill='black', lty=2, alpha=alpha)}

Xs = 0
Ys = 0
H  = 0
H1 = 0
H2 = 0

GLIST   = gList()
VP      = vpList()
GTREE   = NULL

clear.dataplot = function(){
  
  tmp=grid.get("dataplot")
  if (!is.null(tmp)) {
    tryCatch( {grid.remove("dataplot")  },   error=function(e){cat("Error in removing dataplot!")}     )
  }
  #tryCatch( { upViewport(0);  seekViewport("top1");  popViewport();  },   error=function(e){}     )
  #tryCatch( { upViewport(0);  seekViewport("top2");  popViewport();  },   error=function(e){}     )
  
  GLIST   <<- gList()
  VP      <<- vpList()
  GTREE   <<- NULL
}



create.vp1=function(x=NA,y=NA, xscale=NULL, yscale=NULL){
  
  clear.dataplot()
  Xs <<-x0+ncol*w+0.08;
  Ys <<-.5
  W  =1-Xs-0.02;
  H  <<-y0-Ys
  H2 <<-H*0.
  H1 <<-H*1
  
  n = sum(!is.na(y))
  if (n==1){ yt=y[!is.na(y)];   y=c(yt,-yt)}
  if (n==0){ y=c(-1,1)}
  
  if (is.null(xscale)){
    vp=dataViewport(x,y,clip='off', name='topdata1')
    xscale=vp$xscale
    yscale=vp$yscale
  } 
  
  VP[[1]]<<-viewport(Xs,Ys+H2, width=W, height=H1,just=c(0,0),xscale=xscale,yscale=yscale, name='top1')
  
}
create.vp12 = function(x1,y1, x2, y2){
  
  clear.dataplot()
  
  Xs <<-x0+ncol*w+0.08;
  Ys <<-.5
  W  = 1-Xs-0.02;
  H  <<-y0-Ys
  H2 <<-H*0.3
  H1 <<-H*.7
  
  n1 = sum(!is.na(y1))
  if (n1==1){ yt=y1[!is.na(y1)];   y1=c(yt,-yt)}
  if (n1==0){ y1=c(-1,1)}
  
  n2= sum(!is.na(y2))
  if (n2==1){ yt=y2[!is.na(y2)];   y2=c(yt,-yt)}
  if (n2==0){ y2=c(-1,1)}
  
  vp=dataViewport(x1,y1,clip='off', name='topdata1')
  xscale=vp$xscale
  yscale=vp$yscale
  VP[[1]]<<-viewport(Xs,Ys+H2, width=W, height=H1,just=c(0,0),xscale=xscale,yscale=yscale, name='top1')
  
  vp=dataViewport(x2,y2,clip='off', name='topdata2')
  xscale=vp$xscale
  yscale=vp$yscale
  VP[[2]]<<-viewport(Xs,Ys, width=W, height=H2,just=c(0,0),xscale=xscale,yscale=yscale, name='top2')
  
} 
create.vp123 = function(x1,y1, x2, y2,x3,y3){
  
  clear.dataplot()
  
  Xs <<-x0+ncol*w+0.08;
  Ys <<-.5
  W  = 1-Xs-0.02;
  H  <<-y0-Ys
  H2 <<-H*0.3
  H1 <<-H*.7
  
  n1 = sum(!is.na(y1))
  if (n1==1){ yt=y1[!is.na(y1)];   y1=c(yt,-yt)}
  if (n1==0){ y1=c(-1,1)}
  
  n2= sum(!is.na(y2))
  if (n2==1){ yt=y2[!is.na(y2)];   y2=c(yt,-yt)}
  if (n2==0){ y2=c(-1,1)}
  
  vp=dataViewport(x1,y1,clip='off', name='topdata1')
  xscale=vp$xscale
  yscale=vp$yscale
  VP[[1]]<<-viewport(Xs,Ys+H2, width=W, height=H1,just=c(0,0),xscale=xscale,yscale=yscale, name='top1')
  
  vp=dataViewport(x2,y2,clip='off', name='topdata2')
  xscale=vp$xscale
  yscale=vp$yscale
  VP[[2]]<<-viewport(Xs,Ys, width=W, height=H2,just=c(0,0),xscale=xscale,yscale=yscale, name='top2')
  
  vp=dataViewport(x3,y3,clip='off')
  VP[[3]]<<-viewport(Xs+.04, Ys+H2+H1/2, width=H1/2.5, height=H1/2.5,just=c(0,0),
                     xscale=vp$xscale, yscale=vp$yscale,name='top12')
  
}


g.polygon =function(i,...){   g=polygonGrob(...,vp=VP[[i]]);     GLIST<<-gList(GLIST,g)  }
g.segments=function(i,...){   g=segmentsGrob(...,vp=VP[[i]]);     GLIST<<-gList(GLIST,g)  }
g.points  =function(i,...){   g=pointsGrob(...,vp=VP[[i]]);       GLIST<<-gList(GLIST,g)  }
g.lines   =function(i,...){   g=linesGrob(...,vp=VP[[i]]);        GLIST<<-gList(GLIST,g)  }
g.rect    =function(i,...){   g=rectGrob(...,vp=VP[[i]]);         GLIST<<-gList(GLIST,g)  }
g.text    =function(i,...){   g=textGrob(...,vp=VP[[i]]);                     GLIST<<-gList(GLIST,g)  } 
g.xaxis   =function(i,...){   g=xaxisGrob(...,vp=VP[[i]]); g=grid.force(g);   GLIST<<-gList(GLIST,g)  }
g.yaxis   =function(i,...){   g=yaxisGrob(...,vp=VP[[i]]); g=grid.force(g);   GLIST<<-gList(GLIST,g)  } 
g.draw    =function(){ GTREE<<-gTree(children = GLIST,name="dataplot");grid.draw(GTREE)}
g.del     =function(){}



YL=0
YU=0;
get.yscale=function(i){ YRNG=VP[[i]]$yscale; YL<<-YRNG[1];YU<<-YRNG[2] }

#go.plot1=function(){upViewport(0);downViewport('topdata1');  get.yscale() }
#go.plot2=function(){upViewport(0);downViewport('topdata2');  get.yscale() }

set.axis = function(i, xlabel,ylabel,xcol=BLK(), ycol=BLK()){
  
  xname=as.character(runif(1)*1000000);
  yname=as.character(runif(1)*1000000);
  
  g.rect(i)
  if (is.null(xlabel) ) {
    g.yaxis(i,name=yname,gp=ycol)
    #grid.edit(c(ylabel,"labels"),rot=90, just=0.5,gp=gpar(cex=0.8))
    t=GLIST[[length(GLIST)]];
    t$children$labels$rot =90
    t$children$labels$just=.5
    t$children$labels$gp  =gpar(cex=0.8)
    GLIST[[length(GLIST)]] <<-t
    
    g.text(i, ylabel,x=unit(-2.5,"line"), gp=ycol,rot=90)
  } else {
    g.yaxis(i,name=yname,gp=ycol)
    t=GLIST[[length(GLIST)]];
    t$children$labels$rot =90
    t$children$labels$just=.5
    t$children$labels$gp  =gpar(cex=0.8)
    GLIST[[length(GLIST)]] <<-t
    
    g.xaxis(i,name=xname,gp=xcol)
    g.text(i,ylabel,x=unit(-2.5,"line"), gp=ycol,rot=90)
    g.text(i,xlabel,y=unit(-2.5,"line"), gp=xcol) 
  }
}

pnull =function() { 
  create.vp1(t,y);
  
  nms=label$tag[[ClickedButton]]
  nms=paste(nms,collapse="$")
  str=sprintf("%s=NULL: no output available.", nms)
  g.text(1, str,  gp=BLK(),default.units='npc')
  g.rect(1)
  g.draw()
 # cat(ClickedButton,'\n')
}
ptime =function() { 
  create.vp1(t,y);
  g.points(1, t,  y, gp=BLK(),default.units='native',pch=1)
  g.lines( 1, t,  y, gp=BLK(),default.units='native')
  set.axis(1,'time','data',xcol=RED())
  g.draw()
}

pdata =function() {
  create.vp1(t,y);
  g.points(1, t,  y, gp=RED(),default.units='native',pch=1)
  g.lines(1,  t,   y, gp=RED(),default.units='native')
  set.axis(1, 'time','data',ycol=RED())
  g.draw()
}

pmarg_lik =function() {
  
  if (is.null(o$data)) {
    pdata()
    return(NULL)
  }
  
  Yts=get.Yts()
  err=y-Yts

  create.vp12(t,y, t, err);
  ################################################################################
  g.points(1, t,  y,    gp=BLK(),default.units='native',pch=1)
  g.lines(1,  t,  Yts,  gp=RED(),default.units='native')
  set.axis(1,NULL,'data')
  ################################################################################
  g.lines(2, t, err,   gp=BLK(), default.units='native')
  set.axis(2,'time','residual')
  ################################################################################
  g.draw()
 
 
}

pR2   = function() {
  
  if (is.null(o$data)) {
    pdata()
    return(NULL)
  }
  
  Yts=get.Yts()
  err=y-Yts
  create.vp123(t,y, t, err,y,y);
  #create.vp12(t,y, t, err);
  ################################################################################
  g.points(1, t,  y,    gp=BLK(),default.units='native',pch=1)
  g.lines(1,  t,  Yts,  gp=BLK(),default.units='native')
  set.axis(1,NULL,'data')
  ################################################################################
  g.lines(2, t, err,   gp=BLK(), default.units='native')
  set.axis(2,'time','residual')
  ################################################################################
  get.yscale(3)
  g.rect(3,   gp=gpar(fill='white'))
  g.points(3, Yts,        y,         gp=RED(),default.units='native',pch=1)
  g.lines(3,  c(YL,YU),   c(YL,YU),  gp=BLK(),default.units='native')
  
  set.axis(3,'data','fitted')
  
  g.yaxis(3,name='yaxis1',gp=RED())
  g.xaxis(3,name='xaxis1',gp=RED())
  #grid.force()
  #grid.edit(c("yaxis1","labels"), rot=90, just=0.5,gp=gpar(cex=0.5))
  #grid.edit(c("xaxis1","labels"), just=0.5,        gp=gpar(cex=0.5))
  #g.text(3,"data",  x=unit(-2.5,"line"),rot=90,   gp=gpar(cex=0.5))
  #g.text(3,"fitted",y=unit(-2.5,"line"),          gp=gpar(cex=0.5))
  g.draw()
}
pRMSE = function() {
  
  if (is.null(o$data)) {
    pdata()
    return(NULL)
  }
  
  
  Yts=get.Yts()
  err=y-Yts
  create.vp12(t,y, t, err);
  ################################################################################
  g.points(1, t,  y,    gp=BLK(),default.units='native',pch=1)
  g.lines(1,  t,  Yts,  gp=BLK(),default.units='native')
  set.axis(1,NULL,'data')
  ################################################################################
  g.lines( 2,t, err,   gp=RED(), default.units='native')
  set.axis(2,'time','residual')
  ################################################################################
  g.draw()
  
}
psig2 = function() {
   pRMSE()
}

ptrend.Y =function() {
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=RED(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.SD =function() {
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, SD, gp=RED(0.3),default.units='native')
  g.polygon(1,t2t, SD, gp=gpar(col='red',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
}
ptrend.CI =function() {
  
  get.T()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=RED(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='red',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
}

ptrend.slp =function() {
  
  get.T()
  
  create.vp12(t,CI, t, SlpSD);
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines( 2,t,   Slp, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','slope')
  ################################################################################
  g.draw()
  
}
ptrend.slpSD =function() {
  
  get.T()
  
  create.vp12(t,CI, t, SlpSD);
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,SlpSD, gp=RED(0.3),default.units='native')
  g.lines(2, t,   Slp, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','slope')
  ################################################################################
  g.draw()
  
}
ptrend.slpSgnPosPr =function() {
  
  get.T()
  
  create.vp12(t,CI, t, c(0, 1.01)) ;
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines(2, t,   t-t+0.5, gp=gpar(col='black',lty=2), default.units='native')
  g.text(2,'Prob=0.5',   x=t[1], y=0.5, just=c(0,0),   default.units='native')
  g.lines(2, t,   SlpSignPos, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Prob(slope>0)')
  ################################################################################
  g.draw()
  
}
ptrend.slpSgnZeroPr =function() {
  
  get.T()
  
  create.vp12(t,CI, t, c(0, 1.01)) ;
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines(2, t,   t-t+0.5, gp=gpar(col='black',lty=2), default.units='native')
  g.text(2,'Prob=0.5',   x=t[1], y=0.5, just=c(0,0),   default.units='native')
  g.lines(2, t,   SlpSignZero, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Prob(slope>0)')
  ################################################################################
  g.draw()
  
}
ptrend.order =function( ) {
  get.T()
  get.tcp()
  
  create.vp12(t,CI, t, c(0,Order));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines(2, t,  Order, gp=RED(), default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','order')
  ################################################################################
  g.draw()
}

#########################################################
#########################################################
ptrend.ncp =function() {
  
  get.tcp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_mean, ncp_mean), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
ptrend.ncp_median =function() {
  
  get.tcp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_median, ncp_median), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
ptrend.ncp_mode =function() {
  
  get.tcp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_mode, ncp_mode), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
ptrend.ncp_pct90 =function() {
  
  get.tcp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_pct90, ncp_pct90), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
ptrend.ncp_pct10 =function() {
  
  get.tcp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_pct10, ncp_pct10), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
ptrend.ncpPr =function() {
  
  
  get.tcp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

ptrend.cpOccPr =function( ) {
  
  get.T()
  get.tcp()
 
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
 
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

ptrend.cp =function( ) {
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.cpPr =function( ) {
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

ptrend.cpAbruptChange =function( ) {
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

ptrend.cpCI =function() {
  
  
  get.T()
  get.tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

#########################################################
#########################################################

#########################################################
#########################################################

ptrend.pos_ncp =function() {
  
  get.pos_tcp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

ptrend.pos_ncpPr =function() {
  
  
  get.pos_tcp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

ptrend.pos_cpOccPr =function( ) {
  
  get.T()
  get.pos_tcp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

ptrend.pos_cp =function( ) {
  
  get.T()
  get.pos_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.pos_cpPr =function( ) {
  
  get.T()
  get.pos_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

ptrend.pos_cpAbruptChange =function( ) {
  
  get.T()
  get.pos_tcp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

ptrend.pos_cpCI =function() {
  
  
  get.T()
  get.pos_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}






#########################################################
#########################################################

ptrend.neg_ncp =function() {
  
  get.neg_tcp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

ptrend.neg_ncpPr =function() {
  
  
  get.neg_tcp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

ptrend.neg_cpOccPr =function( ) {
  
  get.T()
  get.neg_tcp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

ptrend.neg_cp =function( ) {
  
  get.T()
  get.neg_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.neg_cpPr =function( ) {
  
  get.T()
  get.neg_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

ptrend.neg_cpAbruptChange =function( ) {
  
  get.T()
  get.neg_tcp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

ptrend.neg_cpCI =function() {
  
  
  get.T()
  get.neg_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

#########################################################
#########################################################

ptrend.inc_ncp =function() {
  
  get.inc_tcp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

ptrend.inc_ncpPr =function() {
  
  
  get.inc_tcp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

ptrend.inc_cpOccPr =function( ) {
  
  get.T()
  get.inc_tcp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

ptrend.inc_cp =function( ) {
  
  get.T()
  get.inc_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.inc_cpPr =function( ) {
  
  get.T()
  get.inc_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

ptrend.inc_cpAbruptChange =function( ) {
  
  get.T()
  get.inc_tcp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

ptrend.inc_cpCI =function() {
  
  
  get.T()
  get.inc_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}


#########################################################
#########################################################

ptrend.dec_ncp =function() {
  
  get.dec_tcp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

ptrend.dec_ncpPr =function() {
  
  
  get.dec_tcp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

ptrend.dec_cpOccPr =function( ) {
  
  get.T()
  get.dec_tcp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

ptrend.dec_cp =function( ) {
  
  get.T()
  get.dec_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

ptrend.dec_cpPr =function( ) {
  
  get.T()
  get.dec_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

ptrend.dec_cpAbruptChange =function( ) {
  
  get.T()
  get.dec_tcp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

ptrend.dec_cpCI =function() {
  
  
  get.T()
  get.dec_tcp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}


#######################################################################################
#####################################################################################
# SEASO COMPONENT
######################################################################################
###################################################################################

pseason.Y =function() {
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=RED(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

pseason.SD =function() {
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, SD, gp=RED(0.3),default.units='native')
  g.polygon(1,t2t, SD, gp=gpar(col='red',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
}
pseason.CI =function() {
  
  get.S()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=RED(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='red',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
}

pseason.amp =function() {
  
  get.S()
  
  create.vp12(t,CI, t, AmpSD);
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines( 2,t,   Amp, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','slope')
  ################################################################################
  g.draw()
  
}
pseason.ampSD =function() {
  get.S()
  
  create.vp12(t,CI, t, AmpSD);
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.polygon(1,t2t, CI, gp=gpar(col='black',fill=NA),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,AmpSD, gp=RED(0.3),default.units='native')
  g.lines(2, t,   Amp, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','slope')
  ################################################################################
  g.draw()
  
}
pseason.order =function( ) {
  get.S()
  get.scp()
  
  create.vp12(t,CI, t, c(0,Order));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.lines(2, t,  Order, gp=RED(), default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','order')
  ################################################################################
  g.draw()
}

#########################################################
#########################################################
pseason.ncp =function() {
  
  get.scp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_mean, ncp_mean), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
pseason.ncp_median =function() {
  
  get.scp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_median, ncp_median), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
pseason.ncp_mode =function() {
  
  get.scp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_mode, ncp_mode), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
pseason.ncp_pct90 =function() {
  
  get.scp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_pct90, ncp_pct90), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
pseason.ncp_pct10 =function() {
  
  get.scp()
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  g.lines(1,c(ncp_pct10, ncp_pct10), c(0, min( max(ncpPr)*1.5,1 ) ), gp=gpar(col='red'),default.units='native')
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}
pseason.ncpPr =function() {
  
  
  get.scp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

pseason.cpOccPr =function( ) {
  
  get.S()
  get.scp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

pseason.cp =function( ) {
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

pseason.cpPr =function( ) {
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

pseason.cpAbruptChange =function( ) {
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

pseason.cpCI =function() {
  
  
  get.S()
  get.scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    print(cpCI)
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

#########################################################
#########################################################

#########################################################
#########################################################

pseason.pos_ncp =function() {
  
  get.pos_scp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

pseason.pos_ncpPr =function() {
  
  
  get.pos_scp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

pseason.pos_cpOccPr =function( ) {
  
  get.S()
  get.pos_scp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

pseason.pos_cp =function( ) {
  
  get.S()
  get.pos_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

pseason.pos_cpPr =function( ) {
  
  get.S()
  get.pos_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

pseason.pos_cpAbruptChange =function( ) {
  
  get.S()
  get.pos_scp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

pseason.pos_cpCI =function() {
  
  
  get.S()
  get.pos_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}






#########################################################
#########################################################

pseason.neg_ncp =function() {
  
  get.neg_scp()  
  ncpMax=length(ncpPr)
  
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6, ncpPr[i],just=c(0,0),gp=gpar(col='black',fill='gray'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()
}

pseason.neg_ncpPr =function() {
  
  
  get.neg_scp()
  ncpMax=length(ncpPr)
  create.vp1( xscale=c(-0.5, ncpMax+0.5), yscale=c(0,    min( max(ncpPr)*1.5,1 )  )    )   
  for (i in 1:ncpMax){
    g.rect(1,(i-1)-0.3, 0,0.6 , ncpPr[i],just=c(0,0), gp=gpar(col='black',fill='red'),default.units='native')
  }
  set.axis(1,"Number of changepoints",'Probability')
  g.draw()  
  
}

pseason.neg_cpOccPr =function( ) {
  
  get.S()
  get.neg_scp()
  
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=RED(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=RED(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()  
}

pseason.neg_cp =function( ) {
  
  get.S()
  get.neg_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1,t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=RDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}

pseason.neg_cpPr =function( ) {
  
  get.S()
  get.neg_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines( 1, t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    #g.rect(cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=gpar(col='black',fill='#777777',alpha=0.2),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    str=sprintf("%.3f",cpPr[i])
    g.text(1, str,x=cp[i],y=YL, just=c(0,0), rot=90, default.units='native',gp=RED())
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines( 2,t,   Prob, gp=BLK(),   default.units='native')
  get.yscale(2)
  for (i in 1:ncp)  g.segments(2,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
  
  
}

pseason.neg_cpAbruptChange =function( ) {
  
  get.S()
  get.neg_scp()
  create.vp12(t,CI, t, cpChange);
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  I=1; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(1, NULL,'trend')
  ################################################################################
  g.lines(2, t,    t-t,  gp=BLK(),default.units='native')
  for (i in 1:ncp){
    g.segments(2,cp[i],cp[i],y0=0,y1=cpChange[i],default.units='native',gp=RED())
  }
  set.axis(2,'time','Abrupt change')
  ################################################################################
  g.draw()
  
}

pseason.neg_cpCI =function() {
  
  
  get.S()
  get.neg_scp()
  create.vp12(t,CI, t, c(0.2,Prob));
  ################################################################################
  g.polygon(1, t2t, CI, gp=BLK(0.3),default.units='native')
  g.lines(1,  t,    Y,  gp=BLK(),   default.units='native')
  get.yscale(1)
  for (i in 1:ncp){
    g.rect(1, cpCI[i,1], YL, cpCI[i,2]-cpCI[i,1], YU-YL,just=c(0,0), gp=RED(0.3),default.units='native')
    g.segments(1, cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
    
  }
  set.axis(1,NULL,'trend')
  ################################################################################
  g.polygon(2,t2t,Prob1,gp=BLK(0.1),default.units='native')
  g.lines(2, t,   Prob, gp=BLK(),   default.units='native')
  I=2; get.yscale(I)
  for (i in 1:ncp)  g.segments(I,cp[i],cp[i],y0=YL,y1=YU,default.units='native',gp=BDASH())
  set.axis(2,'time','Probability')
  ################################################################################
  g.draw()
}




splitString = function(text) {
  strings <- strsplit(text, " ")[[1]]
  if (length(strings) < 2)
    return(text)
  newstring <- strings[1]
  linewidth <- stringWidth(newstring)
  gapwidth <- stringWidth(" ")
  availwidth <- convertWidth(unit(1, "npc"),
                             "in", valueOnly=TRUE)
  for (i in 2:length(strings)) {
    width <- stringWidth(strings[i])
    if (convertWidth(linewidth + gapwidth + width,
                     "in", valueOnly=TRUE) < availwidth) {
      sep <- " "
      linewidth <- linewidth + gapwidth + width     }
    else {
      sep <- "\n"
      linewidth <- width     }
    newstring <- paste(newstring, strings[i], sep=sep)
  }
  newstring
}
 
create.button  = function(x0,y0,w,h,r,str){
  x1=x0+w;
  y1=y0+h;
  a=polygonGrob(  c(x0+r, x0 +r, x1 -r, x1 -r ),   
                  c(y1 - r, y0 +r, y0 +r, y1 -r), 
                  default.units = 'npc',
                  gp=gpar( col=NA,fill='#BBBBBB'),name='middle')    #middle
  b=polygonGrob(  c(x0, x0, x0+r, x0+r, x1-r, x1 ),
                  c(y1 ,y0,  y0+r, y1 -r,y1 -r, y1 ),
                  default.units = 'npc',
                  gp=gpar( col=NA,fill='#F8F8F8'),name='upper')   #upper
  c=polygonGrob(  c(x0 ,x1, x1,  x1-r,x1-r,x0+r ),
                  c(y0 ,y0, y1, y1-r, y0+r,y0+r),  
                  default.units = 'npc',
                  gp=gpar( col=NA,fill='#999999') ,name='bottom')  #bottom
  
  t1   = textGrob(str, x0+r, (y0+y1)/2, just=c(0.,0.5),name='t1')
  t2   = textGrob(str, x0+r, (y0+y1)/2, just=c(0.,0.5),name='t2',gp=gpar(col='red'))
  box  = polygonGrob(  c(x0, x0, x1, x1 ),   
                       c(y1, y0, y0 , y1 ), 
                       default.units = 'npc',
                       gp=gpar( col='#BBBBBB',fill='#666666'),name='box' )    #middle
  gTree(children = gList(box,a,b,c,t2,t1),name=str)
}


create.menugrob= function(nrow=38, x0=0.01, y0=0.96, h=0.025, label ){
  
  ncol=length(label$tagstr)/nrow
  ncol=ceiling(ncol)
  
  r=h/10
  wd=convertWidth(stringWidth(label$maxstr),'npc')
  w = as.numeric(wd)+2*r
  
  G=gList()
  N=length(label$tagstr)
  
  g=0;
  for (i in 1:ncol){
    s=(i-1)*nrow+1
    e=min(i*nrow,N)
    for (j in 1:(e-s+1)){
      g=g+1
      G[[g]]=create.button(x0+(i-1)*w,y0-h*j,w,h,r, label$tagstr[s+j-1])
    }
  }
  
  list(G=G,ncol=ncol,w=w)
}


plot.button <-function(idclicked=-1,label){
  
  if(idclicked==ClickedButton) {
    return(0)
  }
  if (ClickedButton>0) 
    grid.reorder(label$tagstr[ClickedButton],c('box','middle','upper','bottom','t2','t1'))
  grid.reorder(label$tagstr[idclicked],c("middle","upper","bottom","box","t1","t2"))
  
}

 
get.fun=function(){
  FUN1=list()
  for (i in 1:length(label$tagstr))  {
    nm = label$tag[[i]]
    if (length(nm)>=2)      { 
      if (is.null(  o[[nm[1]]][[nm[2]]] ) ){
        FUN1[[i]]=pnull
        next
      }
      f=sprintf("p%s.%s",nm[1],nm[2])  
    }
    else{
      if (is.null(  o[[nm[1]]] ) ){
        FUN1[[i]]=pnull
        next
      }
      f=sprintf("p%s",nm[1])  
    } 
    
    tryCatch({FUN1[[i]]=get(f);}, error=function(e){FUN1[[i]]=pnull } )
    
  }
  FUN1
}


mousedown <- function(button,x,y) {
  
  y1 = y0-nrow*h
  x1 = x0+ncol*w
  
  if(  !(x > x0 && x<x1 && y<y0 && y>y1)) {
    return(NULL)
  }
  
  
  col=ceiling((x-x0)/w)
  row=ceiling((y0-y)/h)
  i  = nrow*(col-1)+row
  if (i>length(label$tagstr))    return(NULL)
  if(ClickedButton==i)         return(NULL)
  
  plot.button(i,label)
  ClickedButton <<- i
 
  cat(label$tagstr[i],"\n") 
  f=FUN[[i]];
  
  if (is.function(f)){
    f()
  } else {
    cat('No plot function defined for this variable\n');
  }
  
  
}

keydown <- function(key) {
  
  if (key=='q' || key=='Q' ){
	dev.off()
	return(NULL)
  }
  
  
  if (key!='Down' && key!='Up' ){
    return(NULL)
  }  
  
  if(ClickedButton<=0){
    i=1
  } else {
    i=ClickedButton+(key=="Down")-(key=="Up")
    i=ifelse(i>length(label$tagstr),1, i)
    i=ifelse(i<1,length(label$tagstr) , i)
  }
  plot.button(i, label)
  ClickedButton<<-i
  
  cat(label$tagstr[i],"\n") 
  f=FUN[[i]];
  if (is.function(f)){
    #f()
  }  
  f()
  
  
  NULL
}


#if (.Platform$OS.type == "windows")        x11(width=2, height=2)
#else                                       x11(width=2, height=2,type = "Xlib")  
#wsize = dev.size('px');
#dpi   = wsize[1]/2
#dev.off()
##use dpi to creat a window with a width of about 600 pixels  
WIDTH=1200
HEIGHT=800
#if (.Platform$OS.type == "windows")     x11(width=WIDTH/dpi , height=HEIGHT/dpi           )
#else                                    x11(width=WIDTH/dpi , height=HEIGHT/dpi, type = "Xlib")

dpi = get.dpi() 
new.interactiveWindow(width=WIDTH/dpi , height=HEIGHT/dpi)

x0   <- 0.01
y0   <- 0.96
h    <- 0.025
nrow <- 38
ClickedButton <- -1;

label <- get.names(o)
FUN   <- get.fun()
DF    = create.menugrob(nrow=nrow,x0=x0, y0=y0, h=h,label=label)
ncol <- DF$ncol
w    <- DF$w
G    <- DF$G

grid.newpage()
grid.draw(G)
grid.reorder(label$tagstr[1],c("box","middle","upper","bottom","t2","t1"))

setGraphicsEventHandlers(prompt = "Hit q to quit", onMouseDown = mousedown,onKeybd = keydown )
getGraphicsEvent() 



}