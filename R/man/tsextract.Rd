\name{tsextract}
\alias{tsextract}
\title{Bayesian changepoint detection and time series decomposition}
\description{
		Extract the result of a single time series from an object of class beast
}
\usage{
    tsextract( x, index = 1 ) 
}
\arguments{

  \item{x}{a "beast" object returned  by \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}}, or \code{\link[Rbeast]{beast123}}. It may contain one or many time series.}
  
  \item{index}{an integer (default to 1 ) or a vector of two integers to specify the index of the time series to extract if \code{x} contains results for multiple time series. If \code{x} has 1 time series, \code{index} should be always 1. If \code{x} is returned by \code{\link[Rbeast]{beast123}} applied to a 2D input,\code{index} should be a single index. If \code{x} is from \code{\link[Rbeast]{beast123}} applied to 3D arrays of time series (e.g., stacked satellite images), \code{index} can be a linear index or two subscripts to specify the row and column of the desired pixel/grid.}
  
}
  
\value{
   A LIST object of the result for the chosen time series, which contains the same field as \code{x}.
   }
 
\note{
   Use this function only to manually and interactively examine individual times series. If the purpose is to loop through \code{x}, the use of direct indexing is much faster. For example, if \code{x} is a beast object for a 300x200x1000 3D array (row x col x time), use x$trend$Y[20,40,] to get the fitted trend at the pixel of row 20 and col 40.}
\seealso{
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{tetris}}, \code{\link[Rbeast]{geeLandsat}} 
} 
\examples{

 library(Rbeast)
 data(simdata)
 
 \donttest{
 # handle only the 1st ts
 out=beast(simdata[,1]) 
 }
 
 \dontrun{
 # handle all the ts
 out=beast123(simdata, metadata=list(whichDimIsTime=1))  
 
 plot(out,1)
 plot(out,2)
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

