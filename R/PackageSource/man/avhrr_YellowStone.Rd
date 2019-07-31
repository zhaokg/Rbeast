\name{avhrr_YellowStone}
\alias{avhrr_YellowStone}
\docType{data}
\title{30 years' AVHRR NDVI data at a Yellostone site}
\description{
  \code{avhrr_YellowStone} is a vector comprising 30 years' AVHRR NDVI data at a Yellostone site
}
\usage{
    data(avhrr_YellowStone)
}
\source{
Rbeast v0.2.1
}
 
\examples{
library(Rbeast)
data(avhrr_YellowStone)
plot(avhrr_YellowStone,type='l')
\donttest{
result=beast(avhrr_YellowStone)
plot(result)
}
}
\keyword{misc}
