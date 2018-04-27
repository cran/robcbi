\name{glm.summaries}
\alias{covar.glm}
\alias{Rank.glm}
\alias{rscale.glm}
\alias{weights.glm}

\title{
New accessor functions for linear model fits of glm class
}
\description{
All these functions are methods, respectively the covariance matrix of the 
parameter estimates, the rank, the scale estimate and the computed weights 
for class "glm" objects. 

}

\usage{
\method{covar}{glm}(object)
 
\method{Rank}{glm}(object)
 
\method{rscale}{glm}(object)
 
\method{weights}{glm}(object)
} 

\arguments{
  \item{object}{
an object inheriting from class glm, usually the result of a call to \code{glm}}

}

\details{
The generic accessor functions \code{coef}, \code{residuals}, \code{fitted}, \code{formula}, 
\code{deviance}, \code{rscale}, \code{covar}, \code{correl}, \code{weights} 
and \code{Rank} can be used to extract various useful features of the value returned by 
 \code{glm}. 
}


\references{
Marazzi, A. (1993). \emph{Algorithms, Routines, and S-functions for robust Statistics.}
Chapman and Hall, New York.

Kuensch, H.R., Stefanski L.A., Carroll R.J. (1989).
Conditionally unbiased bounded-influence estimation in general regression models,
with application to generalized linear models. 
Journal of the American Statistical Association, 84, 460-466.
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

