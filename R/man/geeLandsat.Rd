\name{geeLandsat}
\alias{geeLandsat}
\title{Landsat reflectance and NDVI time series from Google Earth Engine}
\description{
    Get Landsat reflectance and NDVI time series from Google Earth Engine given longitude and latitude 
}
\usage{
   geeLandsat(lon=NA, lat=NA, radius=100, stat='mean',timeout=700)
}
\arguments{

  \item{lon}{numeric within [-180,180]}
  \item{lat}{numeric within [-90,  90]}
  \item{radius}{a positive number ( <=500 meters ); the radius of a buffer around the given latitude and longitude for aggregation. If \code{radius=0}, the single pixel at the lat and lon will be retrieved}
  \item{stat}{character; if \code{radius>0}, used to specify the spatial aggregation method for pixels in the buffer. Possible values are 'mean','min','max', or 'median'.}
  \item{timeout}{integer; the seconds elapsed to wait for connection timeout. See the note for an explanation.}
}
  
\value{
    a data.frame object consisting of dates, sensor type, reflectances, and NDVI for the requested location. It contains only valid and clear-sky values as obtained by referring to the standard clouds flags.
   }
 
\note{
    As a poor man's scheme to interact with Google Earth Engine, \code{geeLandsat} should be used only for occasional retrieval of Landsat time series at a few sites, NOT for batch downloading for thousands of sites in a R loop. This procedure is provided to get example time series for testing BEAST. Behind the scene, this function calls to a free Python-based server using my own GEE credential. Normally it takes several seconds to retrieve one time series, but as a free cloud service, the Python server only offers 100 seconds of free CPU time per day, with throttling applied. So it may take up to a few mins to get a time series on your end. It may fail due to connection timeout; if so, give it a few tries. If you need to retrieve data for thousands or millions of sites, please contact the author.
   }
 
\seealso{
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{tetris}}
}
\examples{


 library(Rbeast)
\dontrun{ 
 df = geeLandsat(lon=-80.983877,lat= 40.476882) #if it fails, try a few more times before giving up
 print(df)
}

 
}


\references{
\enumerate{
\item Zhao, K., Wulder, M.A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick, B., Zhang, X. and Brown, M., 2019. Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm. Remote Sensing of Environment, 232, p.111181 (the beast algorithm paper).
\item Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection. Remote Sensing of Environment, 132, pp.102-119 (the Bayesian MCMC scheme used in beast). 
\item Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine. ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261(a beast application paper). 
}

}

\keyword{misc}
