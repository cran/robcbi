\name{predict.cubinf}
\alias{predict.cubinf}

\title{
Prediction methods for objects of class "cubinf"
}
\description{
Predictions provided by a model fit when method is "cubinf".}
\usage{
\method{predict}{cubinf}(object, newdata, type = c("link", "response", "terms"), 
   se.fit = FALSE, terms = labels(object$terms), ...)
}

\arguments{
  \item{object}{An object of class "cubinf" for which predictions are desired.
}
  \item{newdata}{Specify the explanatory variables to used.
}
  \item{type}{The prediction type.
}
  \item{se.fit}{Logical to specify if standard errors are returned or not.
}
  \item{terms}{The terms in newdata.
}
  \item{\dots}{Additional arguments affecting the predictions produced.
}
}

\value{
The value returned depends on type.
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
 \code{\link{predict.glm}}
}
\examples{
 library(robcbi)
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
 df <- data.frame(lVol = log(Vol), lRate = log(Rate), Resp = Resp)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,data=df,method="cubinf",ufact=3.2)
 set.seed(123)
 rVol  <- runif(20,0.4,3.7); rRate <- runif(20,0.3,3.75)
 newdat <- data.frame(lVol=log(rVol),lRate=log(rRate))
 predict(z.cub, newdat, type="response")
}

\keyword{stats }
\keyword{robust } 
