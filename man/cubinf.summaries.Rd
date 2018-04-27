\name{cubinf.summaries}
\alias{residuals.cubinf}
\alias{summary.cubinf}
\alias{covar.cubinf}
\alias{deviance.cubinf}
\alias{family.cubinf}
\alias{Rank.cubinf}
\alias{rscale.cubinf}
\alias{weights.cubinf}

\title{
Accessing linear model fits of cubinf class
}
\description{
All these functions are methods, respectively the summary, the covariance matrix of the 
parameter estimates,the deviance, the family (model of function used by glm), the residuals, 
the rank, the scale estimate and the computed weights for class "cubinf" objects. 

}

\usage{
\method{residuals}{cubinf}(object, type = c("deviance", "pearson", "response"), ...)

\method{summary}{cubinf}(object, ...)
 
\method{covar}{cubinf}(object)
 
\method{deviance}{cubinf}(object, ...)

\method{family}{cubinf}(object, ...)
 
\method{Rank}{cubinf}(object)
 
\method{rscale}{cubinf}(object)
 
\method{weights}{cubinf}(object)
} 

\arguments{
  \item{object}{
an object inheriting from class cubinf, usually the result of a call to \code{glm} were 
\code{method="cubinf"}. 
}
 \item{type}{A character string for the residuals type}
 \item{\dots}{Optional arguments. For summary  can be correlation=TRUE}
}


\details{
The generic accessor functions \code{coef}, \code{residuals}, \code{fitted}, \code{formula}, 
\code{deviance}, \code{rscale}, \code{r.squared}, \code{covar}, \code{correl}, \code{weights} 
and \code{Rank} can be used to extract various useful features of the value returned by 
 \code{glm}. 
The old function  
The class "lm" functions \code{effects}, \code{alias}, \code{add1}, \code{drop1}, code{kappa},
\code{proj}, \code{step}, \code{influence}, \code{anova} and \code{sensitivity} are not
implemented to objects of class "cubinf".
}

\value{
For summary.cubinf a list is returned with the following components: 
\item{terms }{The terms object used in fitting this model. }
\item{coefficients }{A matrix with three columns, containing the coefficients, their 
 standard errors and the corresponding t-statistics.}
\item{dispersion }{Dispersion coefficient } 
\item{df }{The number of degrees of freedom for the model and for the residuals. }
\item{deviance.resid }{The deviance residuals}
\item{family }{The family function used: binomial() or poisson()}
\item{cov.unscaled }{Unscaled covariance matrix; i.e, a matrix such that multiplying it by 
 an estimate of the error variance produces an estimated covariance matrix for the 
 coefficients. }
\item{correlation }{Correlation coefficient matrix for the coefficients in the model.}
\item{deviance }{The vector of deviance. }
\item{null.deviance }{The vector of deviance for the model with only an intercept. }
\item{iter }{The number of iterations required to compute the coefficient estimates. }
\item{nas }{A logical vector whose i-th component is equal to TRUE if the i-th coefficient 
is NA. }
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
The model fitting function \code{\link{glm}}, \code{\link{cubinf}}
}

\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 lVol <-log(Vol); lRate <- log(Rate)
 z.glm <- glm(Resp~lVol+lRate,family=binomial)
 summary(z.glm)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
 summary(z.cub)
 weights(z.cub)
 covar(z.cub)
 deviance(z.cub)
 Rank(z.cub)
 residuals(z.cub)
 rscale(z.cub)
}

\keyword{robust }
\keyword{stats }

