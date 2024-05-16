\name{simdata}
\alias{simdata}
\docType{data}
\title{Simulated time series to test BEAST}
\description{
  \code{simdata} is a 300 x 3 matrix, consisting three time series of length 300. Currently, the three time series are the same. It is used to illustrate BEAST can handle multiple time series at a single function call.
 of BEAST.  
}
\usage{
    data(simdata)
}
\source{
Rbeast v0.9.2
}
 
\examples{

 library(Rbeast)
 data(simdata)
 plot(simdata[,1],type='l')
 
 \dontrun{ 
 #out=beast123(simdata) # Error: whichDimIsTime has to be specified. See below
 out=beast123(simdata, metadata=list(whichDimIsTime=1))  
 
 plot(out,1)
 plot(out,2)
 plot(out,3)
 }
}
\keyword{misc}

\references{
\enumerate{
\item Zhao, K., Wulder, M.A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick, B., Zhang, X. and Brown, M., 2019. Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment, 232, p.111181 (the beast algorithm paper).
\item Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection. Remote Sensing of Environment, 132, pp.102-119 (the Bayesian MCMC scheme used in beast). 
\item Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine. ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261(a beast application paper). 
}

}
