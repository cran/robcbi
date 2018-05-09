\name{plot.cubinf}
\alias{plot.cubinf}

\title{Diagnostic plots for an object of class "cubinf"
}
\description{
      Six plots are available:
      Residuals vs Fitted Values, Sqrt of abs(Residuals) vs Fitted Values, 
      Response vs Fitted Values" and  QQline of Residuals. 
      r-f spread plot is not available and Cook's distances are not available for objects of class "cubinf".
}
\usage{
\method{plot}{cubinf}(x, residuals = NULL, smooths = FALSE, rugplot = FALSE, id.n = 0, 
  ask = TRUE, num=0,  ...)
}
\arguments{
  \item{x}{
An object of class "cubinf"
}
  \item{residuals}{The residuals to be used in the plots if not null.
}
  \item{smooths}{Logical indicating if a smoother should be added to most plots.
}
  \item{rugplot}{Logical indicating if a "rug" representation of the data should be added 
to the plot.
}
  \item{id.n}{Number of points to be labelled in each plot, starting with the most extreme.
}

  \item{ask}{If ask=TRUE, the function operates in interactive mode.
}
  \item{\dots}{Optional arguments for \code{\link[graphics]{par}.}
}
  \item{num}{Integer between 0 and 6. If num>0, plot the num-th choice in the previous list
  in batch mode. 
}
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

