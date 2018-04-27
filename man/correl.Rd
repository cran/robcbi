\name{correl}
\alias{correl}
\alias{covar}
\alias{Rank}
\alias{rscale}
\alias{weights}
\title{
Generic fonctions for objects of class "glm" or "cubinf" 
}
\description{
Compute respectively the correlation matrix or the covariance matrix of the parameter 
estimates,  the rank, the scale estimate and the weights
}
\usage{
correl(object, tl = 1e-10)

covar(object)
 
Rank(object)

rscale(object)

weights(object)
}
\arguments{
  \item{object}{
an object inheriting from class "glm" or "gubinf", the result of a call to \code{glm}
}
  \item{tl}{
The tolerance for a scale factor (denominator) close to zero.
}

}

\details{
The generic accessor functions \code{coef}, \code{residuals}, \code{fitted}, 
\code{formula}, \code{deviance},\code{rscale}, \code{covar}, \code{correl}, 
\code{weights}, \code{Rank} can be used to extract various useful features 
of the value returned by \code{glm}. They are not all available for objects of class 
\code{cubinf}.
\code{scale.estimate} in a previous version of the package is replaced by \code{rscale}.
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
The model fitting function \code{\link{glm}}
}

\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 lVol <-log(Vol); lRate <- log(Rate)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
 correl(z.cub)
 covar(z.cub)
 Rank(z.cub)
 rscale(z.cub)
 weights(z.cub)
}

\keyword{stats }

