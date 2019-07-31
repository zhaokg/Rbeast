\name{simdata}
\alias{simdata}
\docType{data}
\title{Simulated time series to test BEAST}
\description{
  \code{simdata} is a 300 x 3 matrix, consiting three time series of length 300. Currently, the three time series are the same. It is
  used to illustrate BEAST can handle multiple time series at a singel function call.
 of BEAST.  
}
\usage{
    data(simdata)
}
\source{
Rbeast v0.2.1
}
 
\examples{

library(Rbeast)
data(simdata)
plot(simdata,type='l')
\donttest{
result=beast(simdata)
plot(result,1)
plot(result,2)
plot(result,3)
}
}
\keyword{misc}
