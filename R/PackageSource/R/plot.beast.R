plot.beast<-function (x, index, ...) {
 #https://stackoverflow.com/questions/39670646/r-github-package-w-devtools-warning-unknown-macro-item
 
  if(missing(index)) {index=1}
  indx=index;
  if (attributes(x)$algorithm=='beastTrend')
   {
   N=length(x$t)
   x$s=x$t-x$t
   x$sProb=x$tProb-x$tProb
       x$sSD=x$tSD-x$tSD
   }
   
  lm=0.1;  rm=0.02;  tm=0.1;  bm=0.1
  mm=0.08;  wd=(1-tm-bm-mm)/2
  wd1=wd*0.6;  wd2=wd*0.4
  
  N=length(x$s[,indx])
  
  #x11()
  #===========================================================
  y=x$s[,indx];  sD=x$sSD[,indx]; sD[is.nan(sD)]=0;   sD[is.na(sD)]=0
  xx=c(1:N, seq(N, 1, by=-1) );   yy=c(y-sD, rev(y+sD))
  
  scp=x$scp[,indx];
  scp=scp[!is.na(scp)]
  scpNum=length(scp)
  
  
  tcp=x$tcp[,indx];
  tcp=tcp[!is.na(tcp)]
  tcpNum=length(tcp)
  
  
  
  #--The seasonal compoent---#
  par(plt=c(lm, 1-rm,1-tm-wd1,1-tm)  )
  plot(xx,yy,type='n',ann=F, xaxt="n",yaxt="n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "skyblue")
  axis(1, labels=F, tcl=0.25)
  axis(2, labels=T, padj=1.5, tcl=0.25)
  title(ylab="Seasonal", mgp=c(1.2,0,0) )
  
  polygon(xx,yy, col="#B4B4B4", border="#848484")
  #par(new=TRUE)
  points(1:N,y, type='l', col='red',xlab=NULL,lwd=2,ann=FALSE, bg='red')
  
  if (scpNum > 0)
  {
    for (i in 1:scpNum)
    {
      cpt=scp[i]
      points(c(cpt,cpt),c(par("usr")[3],par("usr")[4])   , type='l', lty="dashed", col='grey10',xlab=NULL,lwd=2,ann=FALSE, bg='grey10')
    }
  }

  
  y=x$sProb[,indx]
  par(plt=c(lm, 1-rm,1-tm-wd,1-tm-wd1), new=TRUE)
  plot(1:N,y,type='n',ann=F, xaxt="n",yaxt="n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "skyblue")
  par(new=TRUE)
  plot(1:N,y,xlab="",ylab="", type='l',col='black',lwd=2, xaxt="n", yaxt="n"  )
  axis(1, labels=T, padj=-1.5, tcl=0.25)
  axis(2, labels=T, padj=1.5,  tcl=0.25)
  title(xlab="Time",ylab="Prob", mgp=c(1.2,1.2,0) )
  
  #=========================================================== 
  y=x$t[,indx];  sD=x$tSD[,indx]; sD[is.nan(sD)]=0;   sD[is.na(sD)]=0
  xx=c(1:N, seq(N, 1, by=-1) );   yy=c(y-sD, rev(y+sD))
  
  par(plt=c(lm, 1-rm,1-tm-wd-mm-wd1,1-tm-wd-mm) , new=TRUE)
  plot(xx,yy,type='n',ann=F, xaxt="n",yaxt="n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "darkolivegreen3")
  axis(1, labels=F, tcl=0.25)
  axis(2, labels=T, padj=1.5, tcl=0.25)
  title(ylab="Trend", mgp=c(1.2,0,0) )
  polygon(xx,yy, col="#B4B4B4", border="#848484")
  #par(new=TRUE)
  points(1:N,y, type='l',col='red',xlab=NULL,lwd=2,ann=FALSE, bg='red')
  
  if (tcpNum > 0)
  {
    for (i in 1:tcpNum)
    {
     cpt=tcp[i]
     points(c(cpt,cpt),c(par("usr")[3],par("usr")[4])   , type='l', lty="dashed", col='grey10',xlab=NULL,lwd=2,ann=FALSE, bg='grey10')
    }
  }
  
  
  
  y=x$tProb[,indx]
  par(plt=c(lm, 1-rm, bm,bm+wd2), new=TRUE)
  plot(1:N,y,type='n',ann=F, xaxt="n",yaxt="n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "darkolivegreen3")
  par(new=TRUE)
  plot(1:N,y,xlab="",ylab="", type='l',col='black',lwd=2, xaxt="n", yaxt="n"  )
  axis(1, labels=T, padj=-1.5, tcl=0.25)
  axis(2, labels=T, padj=1.5,  tcl=0.25)
  title(xlab="Time",ylab="Prob", mgp=c(1.2,1.2,0) )
  

  #if (requireNamespace("plotly"))
  {
  
   # require("")
   # d=data.frame(x=1:N, s=x$s[,1],t=x$t[,1], s1=x$sCI[,1,1], sProb=x$sProb[,1], tProb=x$tProb[,1])
   # ax <- list(   zeroline = FALSE,   showline = FALSE, showgrid=FALSE,  mirror = "ticks",
   #               gridcolor = toRGB("gray50"),   gridwidth = 2,   zerolinecolor = toRGB("red"),
   #               zerolinewidth = 1,   linecolor = toRGB("black"),   linewidth = 2)
    
   # p1 <- plot_ly(d, x = ~x, y = ~s, name = 'trace 0', type = 'scatter', mode = 'lines')  %>%
    #  add_trace(x=~x, y=~s1, type='scatter', mode='lines', fill='tonexty', fillcolor='rgb(19,10,10)')
    
    #p2 <- plot_ly(d, x = ~x, y = ~sProb, name = 'trace 0',  fill="tozerox", type = 'scatter', mode = 'lines') 
    #p3 <- plot_ly(d, x = ~x, y = ~t, name = 'trace 0', type = 'scatter', mode = 'lines') %>%
    #  layout(xaxis=ax, yaxis=ax, plot_bgcolor='rgb(29,229,229)')
    #p4 <- plot_ly(d, x = ~x, y = ~tProb, name = 'trace 0', type = 'scatter', fill="tozerox", mode = 'lines') %>%
    #  layout(xaxis=ax, yaxis=ax, plot_bgcolor='rgb(129,229,29)')
    
    #p<-subplot(p1,p2,p3,p4,shareX=TRUE, heights =c(0.25,0.2,0.24,0.2), nrows=4)
    #print(p)
  }
  return( invisible(NULL));
}
