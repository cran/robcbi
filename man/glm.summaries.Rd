\name{glm.summaries}
\alias{covar.glm}
\alias{Rank.glm}
\alias{rscale.glm}
\alias{weights.glm}

\title{
Accessor functions for objects the class "glm"
}
\description{
Covariance matrix of the coefficient estimates, rank, scale estimate and the weights 
for class "glm" objects. All these functions are methods.

}

\usage{
\method{covar}{glm}(object)
 
\method{Rank}{glm}(object)
 
\method{rscale}{glm}(object)
 
\method{weights}{glm}(object)
} 

\arguments{
  \item{object}{
An object inheriting from class "glm".}

}

\details{
The generic accessor functions \code{coef}, \code{residuals}, \code{fitted}, \code{formula}, 
\code{deviance}, \code{rscale}, \code{covar}, \code{correl}, \code{weights} 
and \code{Rank} can be used to extract elements from an object returned by \code{glm}. 
}

\seealso{
The model fitting function \code{\link{glm}}
}

\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 lVol <-log(Vol); lRate <- log(Rate)
 z.glm <- glm(Resp~lVol+lRate,family=binomial)
 summary(z.glm)
 covar(z.glm)
 Rank(z.glm)
 rscale(z.glm)
 weights(z.glm)
}

\keyword{stats }

