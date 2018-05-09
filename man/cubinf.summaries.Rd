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
Functions required by the corresponding access functions 
}
\description{Auxiliary functions for residuals(), summary(), covar(), deviance(), family(), Rank(), rscale(), weights().
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
 \item{object}{An object inheriting from class "cubinf".} 
 \item{type}{A character string for the residual type.}
 \item{\dots}{Optional arguments. For summary, it can be correlation=TRUE.}
}


\details{
The generic functions \code{coef}, \code{residuals}, \code{fitted}, \code{formula}, 
\code{deviance}, \code{rscale}, \code{r.squared}, \code{covar}, \code{correl}, \code{weights} 
and \code{Rank} can be used to extract elements from an object of class "cubinf" returned by \code{glm}. 
The class "lm" functions \code{effects}, \code{alias}, \code{add1}, \code{drop1}, code{kappa},
\code{proj}, \code{step}, \code{influence}, \code{anova} and \code{sensitivity} are not
implemented to objects of class "cubinf".
}

\value{
summary.cubinf returns a list with the following components: 
\item{call }{The model formula used in glm. }
\item{terms }{Terms object used in fitting the model. }
\item{coefficients }{A matrix with three columns, containing the coefficients, their 
 standard errors and the corresponding t-statistics.}
\item{dispersion }{Dispersion coefficient } 
\item{df }{Degrees of freedom for model and residuals. }
\item{deviance.resid }{Deviance residuals}
\item{family }{The family function used: binomial or poisson}
\item{cov.unscaled }{Unscaled covariance matrix of coefficient estimates.}
\item{correlation }{Correlation matrix of coefficient estimates.}
\item{deviance }{Deviance.}
\item{null.deviance }{Null deviance. }
\item{iter }{Number of iterations of the main algorithm. }
\item{nas }{A logical vector whose i-th component is TRUE if the i-th coefficient is NA. }
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

