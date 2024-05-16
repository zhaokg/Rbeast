\name{ohio}
\alias{ohio}
\alias{Ohio}
\docType{data}
\title{An irregular Landsat NDVI time series at an Ohio site}
\description{
  \code{ohio} is a data.frame object comprising decades of Landsat-observed surface reflectances and NDVI at an Ohio site
}
\usage{
    data(ohio)
}
\source{
  Landsat images courtesy of the U.S. Geological Survey
}
 
\examples{
 library(Rbeast)
 data(ohio) # Landsat surface references and NDVI at a single pixel observed over time
 str(ohio)
 
\dontrun{
 # ohio$ndvi is a single irregular time series
 y = ohio$ndvi
 o = beast.irreg(y, time=ohio$time,deltat=1/12) 
 plot(o)
 print(o)
 
 # ohio also contains irregular time series of individual spectral bands
 # Below, run the multivariate version of the BEAST algorithm to decompose
 # the 5 time series and detect common changepoints altogether
 
 y = list(ohio$blue, ohio$green, ohio$red, ohio$nir, ohio$swir1);
 o = beast.irreg(y, time=ohio$time,deltat=1/12, freq=12) 
 plot(o)
 print(o) 
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
