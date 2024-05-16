
getseason_meanfilter = function(x, p , residual=FALSE) {

  n  = length(x)
  p1 = round(p/2)
  p2 = p-p1-1

  #############################################################
  # a mean filiting with a window size p
  trend=x
  for (i in (p1+1):(n-p2)) {
    trend[i]=mean(x[(i-p1):(i+p2)])
  }
 
  # fill the start
  t  = ((p1+1):p)/p
  y  = trend [(p1+1):p]
  X  = cbind(t-t+1,t,t^2)
  beta=structure( qr.solve(X, y) , dim=c(3,1))
  
  t  = (1:p1)/p
  X  = cbind(t-t+1,t,t^2)
  dim(beta)=c(3,1)
  trend[1:p1]=X %*% beta
 
  # fill the ending
  t=(n-p1-p2):(n-p2)
  t=(t-n+p2+p1)/p
  y=trend[(n-p2-p1):(n-p2)]
  X=cbind(t-t+1,t,t^2)
  beta=structure( qr.solve(X, y) , dim=c(3,1))
  
  t = (n-p2+1):n
  t = (t-n+p2+p1)/p
  X = cbind(t-t+1,t,t^2)
  trend[(n-p2+1):n]=X %*% beta
  #############################################################

  season = x- trend
  if (residual){
    seasonAvg=season
    for (i in 1:p){
      idx           =seq(i, n, p)
      seasonAvg[idx]=mean(season[idx])
    }
    SSS=season-seasonAvg  
  } else {
    SSS=season
  }
  
  return(SSS)
}



getseason_stl= function(x, p, residual=FALSE) {
  
  #res     = stl( ts(as.vector(x),frequency=p), 'per')
  #season  = as.vector(x)- as.vector( res$time.series[,2])
  res   =NA
  season=NA
  
  return(season)
}

getseason_polyfit= function(x, p, residual=FALSE) {
  
  n  = length(x)
  p1 = round(p/2)
  p2 = p-p1-1
 
  
  maxTrendOder=8
  t          =(1:n)
  XXX =matrix(NA, n, p1*2+1+maxTrendOder)
  for (i in 1:p1) {
    ttt=2*pi*i*t/p  
    XXX[,(i-1)*2+1]=sin(ttt)
    XXX[, i*2]     =cos(ttt)
  }
  
  t   = (1:n)/n
  XXX[,p1*2+1] =1
  y   = x
  
  bestAIC  =1e300
  bestOrder=0
  for ( order in 1:maxTrendOder ) {
    xdim       = p1*2+1+order
    XXX[, xdim]=t^order
    beta       = structure( qr.solve(XXX[,1:xdim], y) , dim=c(xdim,1))
    yfit       = XXX[,1:xdim] %*% beta
    res        = as.vector(y-yfit)
    newAIC     = n*log(sum(res*res))+2*(order+1)
   
    if (newAIC > bestAIC+2){
      break
    } else {
      if (newAIC < bestAIC){
        bestAIC   =newAIC
        bestOrder=order
        
      }
    }
  } 
  
  XXX     = XXX[, 1:(p1*2+bestOrder+1)]
  beta    = qr.solve(XXX, y) 
  beta    = structure( beta[(p1*2+1):(p1*2+bestOrder+1)] ,dim=c(bestOrder+1,1))
  trend  =XXX[,(p1*2+1):(p1*2+bestOrder+1)] %*% beta
  
  
  season = x-trend
  
  if (residual){
    seasonAvg=season
    for (i in 1:p){
      idx           =seq(i, n, p)
      seasonAvg[idx]=mean(season[idx])
    }
    cat('res..........\n')
    SSS=season-seasonAvg  
  } else {
    SSS=season
  }
  return(SSS)
}

getseason_polyfit_bad= function(x, p, goodIdx, residual=FALSE) {
  n  = length(x)
  p1 = round(p/2)
  p2 = p-p1-1
  
  ngood=sum(goodIdx)
  
  
  maxTrendOder=8
  t   =(1:n)
  XXX =matrix(NA, n, p1*2+1+maxTrendOder)
  for (i in 1:p1) {
    ttt=2*pi*i*t/p  
    XXX[,(i-1)*2+1]=sin(ttt)
    XXX[, i*2]     =cos(ttt)
  }
  
  t   = (1:n)/n
  XXX[,p1*2+1] =1
  y   = x
  
  bestAIC  =1e300
  bestOrder=0
  for ( order in 1:maxTrendOder ) {
    xdim       = p1*2+1+order
    XXX[, xdim]=t^order
    beta       = structure( qr.solve(XXX[goodIdx,1:xdim], y[goodIdx]) , dim=c(xdim,1))
    yfit      = XXX[goodIdx,1:xdim] %*% beta
    res       = as.vector(y[goodIdx]-yfit)
    newAIC    = ngood*log(sum(res*res))+2*(order+1)

    if (newAIC > bestAIC+2){
      break
    } else {
      if (newAIC < bestAIC){
        bestAIC   =newAIC
        bestOrder=order
      
      }
    }
  } 
  
  XXX     = XXX[, 1:(p1*2+bestOrder+1)]
  beta    = qr.solve(XXX[goodIdx,], y[goodIdx]) 
  xfull   = XXX %*% structure( beta ,dim=c(p1*2+bestOrder+1,1) )
  xfull[goodIdx] =x[goodIdx]
  
  beta   = structure( beta[(p1*2+1):(p1*2+bestOrder+1)] ,dim=c(bestOrder+1,1))
  trend  = XXX[,(p1*2+1):(p1*2+bestOrder+1)] %*% beta
  
 
  season = xfull-trend
  
  if (residual){
    seasonAvg=season
    for (i in 1:p){
      idx           =seq(i, n, p)
      seasonAvg[idx]=mean(season[idx])
    }
    SSS=season-seasonAvg  
  
  } else {
    SSS=season
  }
  return(SSS)
}


svdbasis = function(x, p , residual=FALSE) {
  
  x  = as.vector(x)
  n  = length(x)
  
  goodidx = !is.na(x)
  ngood   = sum(goodidx)
  
 
  
  if(n==ngood){
    SSS=getseason_polyfit(x,p,residual)  
  } else {
    SSS=getseason_polyfit_bad(x,p,goodidx,residual)  
  }
 

  M    =  floor(n/p)
  SSS  = structure ( SSS[1:(M*p)], dim=c(p,M))
  v    = svd(SSS,p )
  u    = v$u
  
  U     = matrix(NA, n,p)
  M1    =floor((n+(p-1))/p)
  for (i in 1:p){
    ui=u[,i]
    ui=(ui-mean(ui))/sd(ui)
    u1=rep(u[,i],M1)
    u1=u1[1:n]
    U[,i]=u1;
  }
  
  return(U)
}

 


