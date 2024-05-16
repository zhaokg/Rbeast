\name{googletrend_beach}
\alias{googletrend_beach}
\alias{googletrend_beach}
\docType{data}
\title{A monthly Google Trend time series of the US search interest in the word "beach"}
\description{
  \code{googletrend_beach} is a ts object comprising monthly search interest in "beach" from the United States, as reported from Google Trends. Sudden changes in the search trend are attributed to extreme weather events or the covid19 outbreak
}
\usage{
    data(googletrend_beach)
}
\source{
  https://trends.google.com/trends/explore?date=all&geo=US&q=beach
}
 
\examples{
 library(Rbeast)
 data(googletrend_beach) # A monthly ts starting from Jan 2004
 \donttest{
 o = beast(googletrend_beach )  
 plot(o)
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
