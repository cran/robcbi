\name{cubinf.missing}
\alias{anova}
\alias{add1}
\alias{drop1}
\alias{step}

\title{
Missing methods for object of class "cubinf".
}
\description{
These functions are not implemented for object of class "cubinf" 
}
\usage{

anova(object, ...)

add1(object, ...)
 
drop1(object, ...)
 
step(object, ...)
 
}

\arguments{
  \item{object}{
an object inheriting from class "glm" or "cubinf", the result of a call to \code{glm}}
  \item{\dots}{
Optional argument according to the method}
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

\keyword{robust }

