\name{QQline}
\alias{QQline}

\title{
Add a theoretical QQ-line in a plot
}
\description{
Adds a QQ-line for the values in x in the current plot.
}
\usage{
QQline(x, ...)
}
\arguments{
  \item{x}{The sample for QQ-line
}
  \item{\dots}{Graphical parameters
}
}

\value{
The intercept and the slope of the QQ-line are returned
}
\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. 
Wadsworth & Brooks/Cole. 
}

\keyword{stats }

\examples{
library(robcbi)
data(Finney)
Vol  <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
lVol <-log(Vol); lRate <- log(Rate)
z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
x    <- residuals(z.cub, type="deviance")
\dontrun{
qqnorm(x, ylab = "Deviance Residuals")
QQline(x, lty = 2)
}
}
