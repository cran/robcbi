\name{robust.print.summaries}
\alias{print.cubinf}
\alias{print.cubinf.i}
\alias{print.summary.cubinf}
\alias{print.glm.i}

\title{
Print methods for objects of class "cubinf", "cubinf.i", "summary.cubinf" or "glm.i"
}

\description{
Printing linear model fits provided by glm or with method="cubinf"
}
\usage{
\method{print}{cubinf}(x, ai = FALSE, ci = FALSE, A.mat = FALSE, ...)

\method{print}{summary.cubinf}(x, ...)

\method{print}{glm.i}(x, ...)
}

\arguments{
  \item{x}{
  An object result of a call to \code{summary.cubinf} (first usage),
  to \code{glm} with method="cubinf" (second usage),
  to \code{rscale.cubinf} or to \code{summary.cubinf} or to \code{weights.cubinf}
  or to \code{covar.cubinf} (third usage) and respectively to \code{rscale.glm} or 
  or to \code{covar.glm} or to \code{weights.glm}.
}
  \item{ai}{Vector with components a_i=ufact/|Ax_i| (where x_i^T 
    denotes the ith row of the model matrix).
}
  \item{ci}{Vector of the final bias corrections.
}
  \item{A.mat}{The final value of the matrix A.
}
  \item{...}{
 Further optional arguments according to the print method. Implicit argument in all these 
 functions is digits = max(3, .Options$digits - 3).
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
The model fitting function \code{\link{glm}}, \code{\link{cubinf}}
}
\examples{
library(robcbi)
 
## Dobson (1990) Page 93: Randomized Controlled Trial :
 counts <- c(18,17,15,20,10,20,25,13,12)
 outcome <- gl(3,1,9)
 treatment <- gl(3,3)
 print(d.AD <- data.frame(treatment, outcome, counts))
 zD93 <- glm(counts ~ outcome + treatment, family=poisson,method="cubinf",ufact=3.2)
 summary(zD93)
 print(zD93)
 wi <- weights(zD93)
 print(wi)
}

\keyword{stats }

