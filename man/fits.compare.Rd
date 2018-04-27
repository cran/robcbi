\name{fits.compare}
\alias{fits.compare}
\alias{print.fits.compare}
\title{
Generalized linear model fits comparison 
}
\description{
The fits.compare() function accepts a sequence of objects of class
"glm", "cubinf", or "aov" (with optional names), and creates a class "fits.compare"
object. The print.fits.compare() function prints summaries of each
of the input objects in a manner suitable for comparing the
input models. 
}
\usage{
\method{print}{fits.compare}(x, digits = max(3, .Options$digits - 3), ...)

fits.compare(...)
}

\arguments{
  \item{x}{
an object inheriting from class "fits.compare", the result of a call to 
\code{fits.compare}. 
}
  \item{digits}{
minimal number of \emph{significant digits}.
}

  \item{\dots}{
In \code{fits.compare()}, \code{\dots} is a sequence of objects of class "lm", "lm.robust", 
or "aov". 
Otherwise \code{\dots} is further arguments passed to or from \code{print}  method.
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
}
\keyword{stats }

