\name{tetris}
\docType{methods}
\alias{tetris}
\alias{Tetris}
\title{The Tetris game in R}
\description{
   A poor man's implementation of the Tetris game in R. Yes, you are right again: it has nothing to do with time series decomposition, changepoint detection, and time series segmentation. Its only remote connection to Rbeast is that this is a practice script I wrote to learn R graphics for implementing Rbeast.
   \cr \cr
   \if{html}{\figure{tetris.png}{options: width="30\%" alt="Figure: tetris.png"}}
   \if{latex}{\figure{tetris.pdf}{options: width=5cm}}
  }
\usage{
	tetris(height=25, width=14, speed=0.6) 
}
\arguments{
  
  \item{height}{integer; number of rows of the mine grid along the vertical direction.}
  \item{width}{integer; number of columns of the mine grid along the horizontal direction.}
  \item{speed}{ numeric; a time interval between 0.05 and 2 seconds, specifying how fast the tetriminos moves down. The smaller, the faster.}   
 
}

\value{
	Instructions:
		\itemize{
		 \item  Left arrow to move left.
		 \item  Right arrow to move right.
		 \item  Up arrow to rotate.
		 \item  Down arrow to speed up.
		 \item  Space key to sink to the bottom.
		}
 
}

\note{
    This function works only under the Windows OS not Linux or Mac. An interactive graphics window is needed to run this function correctly. So it won't run in RStudio's plot pane. The function will use the x11() or x11(type='Xlib') graphic device to open a pop-up window.
}
 
\seealso{
  \code{\link[Rbeast]{beast}}, \code{\link[Rbeast]{beast.irreg}},  \code{\link[Rbeast]{beast123}}, \code{\link[Rbeast]{minesweeper}},  \code{\link[Rbeast]{geeLandsat}} 
}
\examples{
\donttest{

 library(Rbeast)
 
 \dontrun{
 tetris() 
 
 # A  field of size 20x25 with blocks moving down every 0.1 sec.
 tetris(20,25,0.1)  
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
