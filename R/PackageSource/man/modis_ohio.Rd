\name{modis_ohio}
\alias{modis_ohio}
\docType{data}
\title{14 years' MODIS EVI data at a pixel in Southern Ohio}
\description{
  \code{modis_ohio} is a vector comprising 14 years' MODIS EVI data at a pixel in Southern Ohio.
}
\usage{
    data(modis_ohio)
}
\source{
Rbeast v0.2.1
}
 
\examples{

\donttest{

library(Rbeast)
data(modis_ohio)

plot(modis_ohio,type='l')
result=beast(modis_ohio)
plot(result)

}
}
\keyword{misc}
