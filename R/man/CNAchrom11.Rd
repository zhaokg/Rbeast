\name{CNAchrom11}
\alias{CNAchrom11}
\alias{CNAchrom11}
\docType{data}
\title{DNA copy number alteration data in array-based CGH data for Chromesome 11}
\description{
  \code{CNAchrom11} is a vector of the log2 intensity ratios for cell line GM03576 for Chromosome 11, obtained from Snijders et al. (2001).
}
\usage{
    data(CNAchrom11)
}
\source{
 Snijders et al. (2001), Assembly of microarrays for genome-wide measurement of DNA copy number, Nature Genetics, 29, 263-264 (http://www.nature.com/ng/journal/v29/n3/full/ng754.html).
}
 
\examples{

 library(Rbeast)
 data(CNAchrom11)  
 
 \donttest{
 o = beast(CNAchrom11, season='none') # no periodic component
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
