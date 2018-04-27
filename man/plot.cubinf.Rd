\name{plot.cubinf}
\alias{plot.cubinf}

\title{Plot diagnostics for an "cubinf" object
}
\description{
      Diagnostic plots for elements of class "cubinf". Six plots are currently available:
      Residuals vs Fitted Values, Sqrt of abs(Residuals) vs Fitted Values, 
      Response vs Fitted Values" and  QQline of Residuals. 
      r-f spread plot is not available in R .and Cook's distances are not available for 
      cubinf object.
}
\usage{
\method{plot}{cubinf}(x, residuals = NULL, smooths = FALSE, rugplot = FALSE, id.n = 0, 
  ask = TRUE, num=0,  ...)
}
\arguments{
  \item{x}{
An object of class "cubinf", usually, a result of a call to \code{\link{glm}} where 
method="cubinf"
}
  \item{residuals}{The residuals to be used in the plots if not null
}
  \item{smooths}{logical indicating if a smoother should be added to most plots
}
  \item{rugplot}{logical indicating if a "rug" representation of the data shuold be added 
to the plot
}
  \item{id.n}{number of points to be labelled in each plot, starting with the most extreme.
}

  \item{ask}{If ask=TRUE, plot.lm.robust operates in interactive mode
}
  \item{\dots}{Optional arguments for \code{\link[graphics]{par}}
}
  \item{num}{Integer between 0 and 6. If num>0, plot the num-th choice in the previous list
  in batch mode. 
}
}

\references{
Kuensch, H.R., Stefanski L.A., Carroll R.J. (1989).
Conditionally unbiased bounded-influence estimation in general regression models,
with application to generalized linear models. 
Journal of the American Statistical Association, 84, 460-466.

Marazzi, A. (1993).
Algorithms, Routines, and S-functions for robust Statistics.
Chapman and Hall, New York.
}

\seealso{
\code{\link{glm.fit}}, \code{\link[graphics]{plot.default}}
}

\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 lVol <-log(Vol); lRate <- log(Rate)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
 summary(z.cub)
\dontrun{
 plot(z.cub, smooth=TRUE, ask=TRUE)
}
}

