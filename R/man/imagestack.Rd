\name{imagestack}
\alias{imagestack}
\docType{data}
\title{Decades of Landsat NDVI time series over a small area in Ohio}
\description{
  \code{imagestack} is a LIST containing Landsat-derived NDVI image chips at an Ohio site
}
\usage{
    data(imagestack)
}
\source{
 Landsat images courtesy of the U.S. Geological Survey
}
 
\examples{

\donttest{
 
 data(imagestack) 
 imagestack$datestr # A string vector containg the observation dates of individual ndvi images
 \dontrun{
 imagestack$ndvi    # NDVI images collected over the past several deccades
 }
 plot(imagestack$ndvi[3,4,],type='l') # Plot the raw data at a pixel


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
