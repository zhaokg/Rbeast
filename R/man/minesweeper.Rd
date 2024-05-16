\name{minesweeper}
\docType{methods}
\alias{minesweeper}
 
\title{The Minesweeper game in R}
\description{
   A poor man's implementation of the minesweeper game in R. Yes, you are right: it has nothing to do with time series decomposition, changepoint detection, and time series segmentation. Its only remote connection to Rbeast is that this is a practice script I wrote to learn R graphics for implementing Rbeast.
   \cr \cr
   \if{html}{\figure{minesweeper.png}{options: width="35\%" alt="Figure: minesweeper.png"}}
   \if{latex}{\figure{minesweeper.pdf}{options: width=5cm}}
  }
\usage{
	minesweeper(height=15, width=12, prob=0.1) 
}
\arguments{
  
  \item{height}{integer; number of rows of the mine grid along the vertical direction.  }
  \item{width}{integer; number of columns of the mine grid along the horizontal direction.}
  \item{prob}{ numeric; a fraction between 0 and 1 to specify the probability of mine occurrence in the mine grid.}   
 
}

\value{
	Instructions:
		\itemize{
		 \item  LEFT-click to clear a spot.
		 \item  RIGHT-click to flag a spot.
		 \item  MIDDLE-click(wheel) a cleared and numbered spot to open neighbor spots, if flagged correctly.
		 \item  Click Restart for a new game
		}
 
}

\note{
    An interactive graphics window is needed to run this function correctly. So it won't run in RStudio's plot pane. The function will use the x11() or x11(type='Xlib') graphic device to open a pop-up window.
}
 
\seealso{
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{tetris}}, \code{\link[Rbeast]{geeLandsat}} 
}
\examples{
\donttest{

 library(Rbeast)
 
 \dontrun{
 minesweeper() 
 
 # A mine field of size 20x25 with rougly a 15% of the spots being mines.
 minesweeper(20,25,0.15) 
 }
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
