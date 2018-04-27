\name{plot.fits.compare}
\alias{plot.fits.compare}

\title{
Plot linear model fits comparison 
}
\description{
Plots the results of a call to fits.compare(). 
Plotting the "fits.compare" object results in a sequence of graphical
displays. These displays are designed to be of use in comparing two sets of 
parameter estimates in linear models.
}

\usage{
\method{plot}{fits.compare}(x, xplots = FALSE, ..., ask = TRUE)
}


\arguments{
  \item{x}{
an object inheriting from class "fits.compare", the result of a call to 
\code{fits.compare}. 
}

  \item{xplots}{
  If TRUE, the graphics are displayed.
}

  \item{\dots}{
 further arguments passed to or from \code{plot} method.
}
  \item{ask}{
Graphical parameter, if TRUE (and the R session is interactive) the user is asked for 
input, before a new figure is drawn
}
}
\details{
It is not recommended to compare objects with different structure.
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
The model fitting function \code{\link{glm}}, \code{\link{cubinf}}
}

\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 lVol <-log(Vol); lRate <- log(Rate)
 z.glm <- glm(Resp~lVol+lRate,family=binomial)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
 comp <- fits.compare(z.glm,z.cub)
 comp
\dontrun{
 plot(comp)
}
}

\keyword{stats }

