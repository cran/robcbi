\name{correl}
\alias{correl}
\alias{covar}
\alias{Rank}
\alias{rscale}
\alias{weights}
\title{
Generic functions for objects of classes "glm" and "cubinf" 
}
\description{
Correlation and covariance matrix of the parameter estimates,  
rank, scale estimate, and weights.
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
An object inheriting from class "glm" or "cubinf".
}
  \item{tl}{
Tolerance for a scale factor (denominator) close to zero.
}

}

\details{
The generic functions \code{coef}, \code{residuals}, \code{fitted}, 
\code{formula}, \code{deviance}, \code{rscale}, \code{covar}, \code{correl}, \code{weights}, \code{Rank} 
can be used to extract elements from an object returned by \code{glm}. 
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

