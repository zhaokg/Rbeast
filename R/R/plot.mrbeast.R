
plot.mrbeast<-function(
  x, 
  index = 1,
  vars  = c('st','s','scp','sorder','t','tcp','torder','o','ocp','error'),  
  col   = NULL, 
  main  ="BEAST decomposition and changepoint detection",
  xlab  ='Time',
  ylab  = NULL,
  cex.main =1,
  cex.lab  =1,
  relative.heights= NULL,
  interactive     = FALSE,
  ncpStat   = c('median','mode','mean','pct90','max'),
  
  ... ) 
{
  
  ncpStat=match.arg(ncpStat)
  
  
  #vars = c('st','s','scp','sorder','t','tcp','torder','error')  
  #vars = c('st','s','scp','sorder','samp','t','tcp','torder', 'tslp','o','ocp','error'),  
  #col  = NULL 
  #main ="BEAST decomposition and changepoint detection"
  #xlab ='Time'
  #ylab = NULL
  #cex.main = 1.2
  #cex.lab  = 1
  #relative.heights=NULL
  #interactive=FALSE 
  
  vars     = tolower(vars)
  vars_log = vars=='st'|vars=='s'|vars=='t'|vars=='scp'|vars=='tcp'|vars=='sorder'|vars=='torder'|vars=='error'|vars=='o'|vars=='ocp'|vars=='samp'|vars=='tslp'|vars=='slpsgn';  
  vars     = vars[vars_log]
  
  if (length(ylab)== length(vars_log)){  
    ylab = ylab[vars_log]
  } else{
    ylab = vector(mode='character',length(vars))
    ylab[vars=='st']  ='Y'
    ylab[vars=='s']   ='season'
    ylab[vars=='t']   ='trend'
    ylab[vars=='o']   ='outlier'
    ylab[vars=='scp'] ="Pr(scp)"
    ylab[vars=='tcp'] ="Pr(tcp)"
    ylab[vars=='ocp'] ="Pr(ocp)"
    ylab[vars=='sorder']= expression("order"[S])
    ylab[vars=='torder']= expression("order"[T])
    ylab[vars=='samp']  = 'amplitude'
    ylab[vars=='tslp']  = 'slope'
    ylab[vars=='slpsgn'] = "slpSign"
    ylab[vars=='error']  = "error"
    
  }
  if (length(col)== length(vars_log)){  
    col = col[vars_log]
  } else{
    col = vector(mode='character',length(vars))
    col[vars=='st']    ='#111111'
    col[vars=='s']     ='red'
    col[vars=='scp']   ='red'
    col[vars=='sorder']='red'
    col[vars=='samp']   ='red'
    col[vars=='t']     ='green'
    col[vars=='tcp']   ='green'
    col[vars=='torder']='green'
    col[vars=='tslp']    ='green'
    col[vars=='slpsgn']  ='black'  #the border only
    col[vars=='o']     ='blue'
    col[vars=='ocp']   ='blue'
    col[vars=='error'] ='darkgray'
  }
  
  
  if ( length(x$marg_lik)> 1 ) {
    x=tsextract(x,index)   # more than 1 time series is present
  }  
  q  = length(x$trend$Y)
  
  heights=relative.heights
  if (length(heights)== length(vars_log)){  
    heights = heights[vars_log]
  } else{
    heights = (1:length(vars))
    heights =heights-heights+1
    heights[vars=='st']=.8/q
    heights[vars=='s']= .8/q
    heights[vars=='t']=.8/q
    heights[vars=='o']=.8/q
    heights[vars=='scp']=.2
    heights[vars=='tcp']=.2
    heights[vars=='ocp']=.2
    heights[vars=='sorder']=.2
    heights[vars=='torder']=.2
    heights[vars=='samp']   =.4/q
    heights[vars=='tslp']   =.4/q
    heights[vars=='slpsgn'] =.4/q
    heights[vars=='error']  =0.4/q
  }
  
  hasSeason   =!is.null(x$season)
  hasOutlier  =!is.null(x$outlier)
  hasHarmonic =!is.null(x$season$order)
  hasTOrder   =!is.null(x$trend$order)
  hasData     =!is.null(x$data)
  hasAmp      =!is.null(x$season$amp)
  hasSlp      =!is.null(x$trend$slp)
  
  
  idx =1:length(vars) > 0
  if(!hasAmp)      {  idx   = idx & !(vars=='samp')}
  if(!hasSlp)      {  idx   = idx & !(vars=='tslp'|vars=='slpsgn')}
  if(!hasSeason)   {  idx   = idx & !(vars=='st'|vars=='s'|vars=='sorder'|vars=='scp') }
  if(!hasHarmonic) {  idx   = idx & !(vars=='sorder') }
  if(!hasTOrder)   {  idx   = idx & !(vars=='torder') }
  if(!hasOutlier)  {  idx   =  idx &!(vars=='o'|vars=='ocp') }
  if(!hasData)     {  idx   =  idx &!(vars=='error') }
  
  col         = col[idx]
  vars        = vars[idx]
  ylab        = ylab[idx]
  heights     = heights[idx]
  

 
 
  nPlots=length(vars)
  
  newHeights=NULL;
  for (i in 1:nPlots){
     v=vars[i]
    if (v=='st' || v=='s'  || v== 't' ||  v == 'o' ||  v == 't' ||  v == 'tslp' ||  v == 'samp'  ||  v == 'error')  {
         newHeights=c(newHeights, rep( heights[i], q) );
    } else{
        newHeights=c(newHeights, heights[i]);
    }
  }
  heights=newHeights;
  
 
  
  newYlab=NULL;
  for (i in 1:nPlots){
     v=vars[i]
    if (v=='st' || v=='s'  || v== 't' ||  v == 'o' ||  v == 't' ||  v == 'tslp' ||  v == 'samp'  ||  v == 'error')  {
         newYlab=c(newYlab, rep( ylab[i], q) );
    } else{
        newYlab=c(newYlab, ylab[i]);
    }
  }
  ylab=newYlab;
  
  if(nPlots==0){
    stop("No valid variable names speciffied int the 'vars' argument. Possible names include 'st','t','s','sorder','torder','scp','tcp','samp','tslp','o', 'ocp', and 'error'. ");
  }
  
  #######################################################
  #  Functions and variables to load the outputs
  ########################################################
  
  
  
  t     = x$time; 
  data  = x$data
  t2t   = c(t, rev(t))
  N     = length(t)
  CI    = 0
  Y     = 0
  SD    = 0
  
  
  Prob2Prob  = 0;
  Prob   = 0; 
  Order  = 0
  Slp    =0
  SlpSD  =0
  SlpSignPos =0
  SlpSignZero=0
  
  cp   =0        
  cpCI =0        
  ncp  =0        
  ncpPr=0        
  cpPr =0
  cpChange=0
  HasChangePoint=FALSE;
  
  ":" = function(i,j){ if(i>j){return(NULL)} ;  seq(i,j) }
  
  Yts   = 0
  YtsSD = 0
  Yerr  = 0
  
  get.Yts  = function() {
    
    sig2  = diag(x$sig2)
    
    Yts   <<- list();
    YtsSD <<- list();
    Yerr  <<- list();
    for (i in 1:q) {
 
      Yts[[i]]   <<- x$trend$Y[[i]]
      SD2        = x$trend$SD[[i]]^2 +sig2[i]
      
      if (!is.null(x$season))  {
        Yts[[i]]  <<- Yts[[i]]+x$season$Y[[i]]
        SD2       =  SD2+x$season$SD[[i]]^2
      }
      if (!is.null(x$outlier))  {
        Yts[[i]]  <<- Yts[[i]]+x$outlier$Y[[i]]
        SD2   =  SD2+x$outlier$SD[[i]]^2
      }	
      SD = sqrt(SD2)
      YtsSD[[i]] <<- c(Yts[[i]]-SD, rev(Yts[[i]]+SD) )  
      Yerr[[i]]  <<- data[[i]]-Yts[[i]]
    }
    
  }
  
  get.T      = function(){
    
    Y <<- list()
    SD <<- list()
    CI <<- list()
    Slp <<- list()
    SlpSD <<- list()
    SlpSignPos <<- list()
    SlpSignZero  <<- list()
    for (i in 1:q){
        Y[[i]] <<-x$trend$Y[[i]];
        SD[[i]] <<-c(Y[[i]]-x$trend$SD[[i]],  rev(Y[[i]]+x$trend$SD[[i]])); 
        if ( !is.null(x$trend$CI) )		  CI[[i]] <<-c(x$trend$CI[[i]][,1], rev(x$trend$CI[[i]][,2])) 
        else                            CI[[i]] <<-SD[[i]]
        
        if (hasSlp){
          Slp[[i]]         <<- x$trend$slp[[i]]
          SlpSD[[i]]       <<- c(Slp[[i]]-x$trend$slpSD[[i]],  rev(Slp[[i]]+x$trend$slpSD[[i]]));         
          SlpSignPos[[i]]  <<- x$trend$slpSgnPosPr[[i]]
          SlpSignZero[[i]] <<- x$trend$slpSgnZeroPr[[i]]
        }
    }
    Order       <<- x$trend$order;
    HasChangePoint <<- !is.null(x$trend$cp);
  }
  
  get.tcp    = function(){
    cmpnt   =  x$trend
    cp      <<- cmpnt$cp;         
    cpCI    <<- cmpnt$cpCI;  
    ncp     <<- switch(ncpStat, mode=cmpnt$ncp_mode, median=cmpnt$ncp_median,mean=cmpnt$ncp,pct90=cmpnt$ncp_pct90,max=sum(!is.nan(cp)))
    ncp     <<- base::round(ncp)
    ncpPr   <<- cmpnt$ncpPr
    cpPr    <<- cmpnt$cpPr
    cpChange<<- cmpnt$cpAbruptChange
    HasChangePoint <<- !is.null(x$trend$cp);
    
    Prob        <<- cmpnt$cpOccPr;    
    Prob2Prob   <<- c(Prob,Prob-Prob)
  }
  
  
  ###########################################################
  Amp   = 0
  AmpSD = 0
  get.S = function(){    
    Y  <<- list()
    SD <<- list()
    CI <<- list()
    Amp <<- list()
    AmpSD <<- list()
    for (i in 1:q){
        Y[[i]]     <<-x$season$Y[[i]];
        SD[[i]]     <<-c(Y[[i]] -x$season$SD[[i]],  rev(Y[[i]] +x$season$SD[[i]])); 
        if ( !is.null(x$season$CI) )	CI[[i]] <<-c(x$season$CI[[i]][,1], rev(x$season$CI[[i]][,2]))  
        else                          CI[[i]] <<-SD[[i]]
        
        if (hasAmp){
          Amp[[i]]    <<-x$season$amp[[i]]
          AmpSD[[i]]  <<-c(Amp[[i]]-x$season$ampSD[[i]],  rev(Amp[[i]]+x$season$ampSD[[i]]));    
        }
    }
    Order  <<-x$season$order;
    HasChangePoint <<- !is.null(x$season$cp);
    
  }
  get.scp    = function(){
    cmpnt   =  x$season
    cp      <<-cmpnt$cp;         
    cpCI    <<-cmpnt$cpCI;  
    ncp     <<-switch(ncpStat, mode=cmpnt$ncp_mode, median=cmpnt$ncp_median,mean=cmpnt$ncp,pct90=cmpnt$ncp_pct90,max=sum(!is.nan(cp)))
    ncp     <<-base::round(ncp)
    ncpPr   <<-cmpnt$ncpPr
    cpPr    <<-cmpnt$cpPr
    cpChange<<-cmpnt$cpAbruptChange
    
    Prob         <<-x$season$cpOccPr;    
    Prob2Prob   <<-c(Prob,Prob-Prob)
    HasChangePoint <<- !is.null(x$season$cp);
  }
  get.O      = function(){   
    Y  <<-list()
    SD <<-list()
    
    for (i in 1:q){
        Y[[i]]         <<-x$outlier$Y[[i]]  ;
        SD[[i]]       <<-c(Y[[i]]  -x$outlier$SD[[i]]  ,  rev(Y[[i]]  +x$outlier$SD[[i]]  )); 
        if ( !is.null(x$outlier$CI) )	  CI[[i]]   <<-c(x$outlier$CI[[i]][,1], rev(x$outlier$CI[[i]][,2]))  
        else                              CI[[i]]   <<-SD[[i]]  
    }
  }
  get.ocp    = function(){
    #cp   <<-x$season$cp;         
    #cpCI <<-x$season$cpCI;  
    #ncp <<-sum(!is.nan(cp))
    #ncpPr<<-x$season$ncpPr
    #cpPr <<-x$season$cpPr
    Prob  <<-x$outlier$cpOccPr;    
    Prob2Prob <<-c(Prob,Prob-Prob)
  }
  
  
  YSIDE = 2;
  padj  = -1.2
  tcl   = 0.2
  plotId = 1;
  
  plot.st=function(col=c(0.2,0.2,0.2), ylabel){
    
    for (i in 1:q) {
      
      subplot(plotId,heights )
      plotId <<- plotId+1;
      
      if (hasData){
        Yall=c(YtsSD[[i]],data[[i]])
        plot(   c(t[1],t[N]), c(min(Yall,na.rm = TRUE),max(Yall,na.rm = TRUE)),  type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
      } else {
        plot(   t2t, YtsSD[[i]],  type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');  
      }
      
     ### polygon(t2t, YtsSD[[i]],  col  = rgb(.5,0.5,0.5,.5), border = NA);
      if (hasData){
        points( t, data[[i]], type = 'p', col='#777777');  
      }
      points( t, Yts[[i]],  type = 'l', lwd=2,col=rgb(col[1],col[2],col[3]));
      axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
      axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
      YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	  
      title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
      
    }
    
  }
  
  plot.y =function(col=c(0,1,0), ylabel){
    alpha=0.1
    for (i in 1:q) {
        
        subplot(plotId,heights )
        plotId <<- plotId+1;
        

        if (hasData && !hasSeason){
          Yall=c(CI[[i]],data[[i]])
          plot( c(t[1],t[N]), c(min(Yall,na.rm = TRUE),max(Yall,na.rm = TRUE)), type = 'n',ann=FALSE, xaxt='n', yaxt='n');
        } else {
          plot(   t2t, CI[[i]],   type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
        }
        
        if (hasData && !hasSeason){
          points(t,  data[[i]], type = 'p', col='#777777');
        }
        ###  polygon(t2t, CI[[i]],   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
        points( t,   Y[[i]],    type = 'l',col=rgb(col[1],col[2],col[3] ) );
        yext =par('usr')
        if (HasChangePoint) {
          for (i in 1 : ncp) {
            points(c(cp[i], cp[i]), yext[3:4],type = 'l', lty = 'dashed', col = 'grey10',lwd =1, bg = 'grey10');
          }    
        }
        axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
        axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
        YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
		title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
    }
  }
  
  plot.prob <- function( col=c(0,1,0), ylabel ){
    
    subplot(plotId,heights )
    plotId <<- plotId+1;
    
    alpha=0.2
    plot( c(t2t[1],t2t), c(0.22,Prob2Prob),type = 'n', ann=FALSE, xaxt='n', yaxt='n');
    polygon(t2t, Prob2Prob, col  = rgb(col[1],col[2],col[3],alpha), border = NA);
    points( t,   Prob,  col  = rgb(col[1],col[2],col[3])  ,       lwd = 1,type = 'l' );
    yext =par('usr')
    if (HasChangePoint) {
      for (i in 1 : ncp) {
        points(c(cp[i], cp[i]), yext[3:4],type = 'l', lty = 'dashed', col = 'grey10',lwd =1, bg = 'grey10');
      }    
    }	
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  
  plot.order <- function( col=c(0,1,0), ylabel ){    	
    subplot(plotId,heights )
    plotId <<- plotId+1;
    
    maxOrder=max(1.01,max(Order));	
    plot(   c(t[1],t[length(t)]),   c(0,maxOrder),type = 'n', ann=FALSE, xaxt='n', yaxt='n');
    points( t,               Order,  col  = rgb(col[1],col[2],col[3]) ,       lwd = 1,type = 'l' );
    yext =par('usr')
    if (HasChangePoint) {
      for (i in 1 : ncp) {
        points(c(cp[i], cp[i]), yext[3:4],type = 'l', lty = 'dashed', col = 'grey10',lwd =1, bg = 'grey10');
      }    
    }	
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  plot.amp =function(col=c(0,1,0), ylabel){
    alpha=0.5
    plot(   t, Amp,   type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
    #polygon(t2t, AmpSD,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
    points( t,   Amp,    type = 'l',col=rgb(col[1],col[2],col[3] ));
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  
  plot.slp =function(col=c(0,1,0), ylabel){
    alpha=0.5
    plot(   t2t, SlpSD,   type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
    polygon(t2t, SlpSD,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
    points( t,   Slp,    type = 'l',col=rgb(0,0,0,0.8 ));
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  plot.slpsgn = function(col=c(0,1,0), ylabel){
    alpha=0.5
    SlpSignNeg = 1-SlpSignPos-SlpSignZero
    plot( c(t[1],t[length(t)]), c(0,1),   type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
    #polygon(t2t, SlpSD,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
    y2y   = c(t-t, rev(SlpSignNeg));           polygon(t2t, y2y, col=rgb(0,0,1,alpha),border=NA)
    y2y   = c(SlpSignNeg,   rev(1-SlpSignPos)); polygon(t2t, y2y, col=rgb(0,1,0,alpha),border=NA) 
    y2y   = c(1-SlpSignPos, t-t+1);            polygon(t2t, y2y, col=rgb(1,0,0,alpha),border=NA) 
    points( t,  t-t+0.5 , type = 'l', lty=2, col=rgb(col[1],col[2],col[3] ));
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  plot.o =function(col=c(0,1,0), ylabel){
    
    alpha=0.5
    for (i in 1:q) {
      
      subplot(plotId,heights )
      plotId <<- plotId+1;
        
      plot(   t2t, CI[[i]],   type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
      #polygon(t2t, CI,   col  = rgb(col[1],col[2],col[3],alpha), border = NA);
      #points( t,   Y,    type = 'l',col='#333333');
      segments(t,t-t,t,Y[[i]], col= rgb(col[1],col[2],col[3]))
      axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
      axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
      YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	  title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
    }
  }
  
  plot.oprob <- function( col=c(0,1,0), ylabel ){
    subplot(plotId,heights )
    plotId <<- plotId+1;
    
    alpha=0.2
    plot( c(t2t[1],t2t), c(0.22,Prob2Prob),type = 'n', ann=FALSE, xaxt='n', yaxt='n');
    # polygon(t2t, Prob2Prob, col  = rgb(col[1],col[2],col[3],alpha), border = NA);
    # points( t,   Prob,  col  = rgb(col[1],col[2],col[3])  ,       lwd = 1,type = 'l' );
    segments(t,t-t,t,Prob, col= rgb(col[1],col[2],col[3]),lwd=1.5)
    axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
    axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
    YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
	title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
  }
  
  plot.error=function(col=c(0.2,0.2,0.2), ylabel){
  
    for (i in 1:q) { 
		subplot(plotId,heights )
		plotId <<- plotId+1;
	  
		plot(   t,   Yerr[[i]],  type = 'n',ann = FALSE, xaxt = 'n', yaxt = 'n');
		lines(t, t-t, col= rgb(col[1],col[2],col[3]))
		segments(t,t-t,t,Yerr[[i]], col= rgb(col[1],col[2],col[3]))
		axis(1,     labels = FALSE, padj = padj,                        tcl = tcl);
		axis(YSIDE, labels = TRUE,  padj = ifelse(YSIDE==2,-padj,padj), tcl = tcl)
		YSIDE<<-(YSIDE==2)*4+(YSIDE==4)*2
		title(ylab = ylab[plotId-1],cex=cex.lab, mgp = c(1.25, 1.25, 0));
	}
  }
  
  
  
  plot.new()
  opar= par(c('oma','mar','fig'))
  par(oma=c(0,0,0,0))
  par(mar=c(0,2,0,2))
  par(fig=c(0,1,0,1))
  
  
  lm= strheight('a',units='figure')*dev.size()[2]/dev.size()[1]*4.5
  rm= strheight('a',units='figure')*dev.size()[2]/dev.size()[1]*2.3
  bm= strheight('a',units='figure')*4.0
  tm= strheight('a',units='figure')*4.0
  
  #######################################################
  #  Create a subplot given the relative heights of the vertical plots
  ########################################################
  subplot <- function(i, heights){
    Hmiddle  =1-tm-bm
    h       = (heights/sum(heights))*Hmiddle
    hcum    = cumsum(h)
    plt     = c(lm, 1-rm, 1-tm-hcum[i], 1-tm-hcum[i]+h[i])
    par(plt=plt,new=TRUE)
  }
  
  
  for (i in 1:length(vars)){
    
    ytitle = ylab[i]
    var    = tolower(vars[i])
    clr    = as.numeric(col2rgb( col[i]  ))/255
    
    
    if (var=='st')    {  get.Yts();           plot.st(clr,    ytitle)  }
    if (var=='s' )    {  get.S(); get.scp();  plot.y(clr,     ytitle)  }
    if (var=='t' )    {  get.T(); get.tcp();  plot.y(clr,     ytitle)  }
    if (var=='scp')   {  get.S(); get.scp();  plot.prob(clr,  ytitle) }
    if (var=='tcp')   {  get.T(); get.tcp();  plot.prob(clr,  ytitle) }
    if (var=='sorder'){  get.S(); get.scp();  plot.order(clr, ytitle)  }
    if (var=='torder'){  get.T(); get.tcp();  plot.order(clr, ytitle)  }
    if (var=='samp')  {  get.S();             plot.amp(clr,   ytitle)  }
    if (var=='tslp')  {  get.T();             plot.slp(clr,   ytitle)  }    
    if (var=='slpsgn'){  get.T();             plot.slpsgn(clr, ytitle)  }   
    if (var=='o')     {  get.O();             plot.o(clr, ytitle)  }
    if (var=='ocp')   {  get.ocp();           plot.oprob(clr, ytitle)  }
    if (var=='error') {  get.Yts();           plot.error(clr, ytitle)  }
    
    #title(ylab = ytitle,cex=cex.lab, mgp = c(1.25, 1.25, 0));
    
    if(i==1){
      #mtext(main,line=.5,cex=cex.main,font=2)
    }
    
    if(i==nPlots){
      title(xlab = xlab ,mgp = c(1.2, 1.2, 0));  
      if (sum(class(t)=='Date') >0 ){
        axis.Date(1,  t, labels = TRUE, padj =padj, tcl =tcl,cex.axis=1);	  	
      }	 else{
        axis(  1,   labels = TRUE, padj =padj, tcl =tcl,cex.axis=1);
      }
      
    }
    
  }
  
  par(opar)
  #text(mean(time) / 2.0, 'NO SEASONAL COMPONET.IGNORE THIS POLOT');
  
}

#plot.mrbeast(o,vars=c('st','s','scp','t','tcp','o','ocp'))
