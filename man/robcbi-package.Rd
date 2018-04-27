\name{robcbi-package}
\alias{robcbi-package}
\docType{package}
\title{
Robust fit for discrete generalized linear model
}
\description{
Conditionally unbiased bounded influence estimates as described in Kuensch et al.  
(1989) in three special cases of the generalized linear model: Bernoulli, Binomial,  
and Poisson distributed responses. 
}
\details{
\tabular{ll}{
Package: \tab cubinf\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-07-04\cr
License: \tab  GPL (>= 2)\cr
}
}

\author{
A. Marazzi  <Alfio.Marazzi@chuv.ch>

Maintainer: A. Randriamiharisoa <Alex.Randriamiharisoa@chuv.ch>
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

\examples{
 library(robcbi)
# First example
 data(Finney)
 Vol <- Finney$Vol; Rate <- Finney$Rate; Resp <- Finney$Resp
\dontrun{
 plot(Vol,Rate,type="n")
 points(Vol[Resp==0],Rate[Resp==0],pch=5, cex=1.2)
 points(Vol[Resp==1],Rate[Resp==1],pch=16,cex=1.2)
}
 lVol <-log(Vol); lRate <- log(Rate)
 z.glm <- glm(Resp~lVol+lRate,family=binomial)
 summary(z.glm)
 z.cub <- glm(Resp~lVol+lRate,family=binomial,method="cubinf", ufact=3.2)
 summary(z.cub)
 weights(z.cub)
\dontrun{
 plot(z.cub, smooth=TRUE, ask=TRUE)
}
 comp <- fits.compare(z.glm,z.cub)
 comp
\dontrun{
 plot(comp)
}
# Second example
 data(Breslow)
\dontrun{ 
 help(Breslow)
}
 y  = Breslow$sumY
 x1 = Breslow$Age10
 x2 = Breslow$Base4
 x3 = rep(0,length(y))
 x3[Breslow$Trt=="progabide"] = 1
 x4 = x2*x3
 CBA = glm(y~x1+x2+x3+x4,family=poisson,method=cubinf,ufact=3.2)
\dontrun{
 plot(CBA,num=5)
}
 weights(CBA)
#
# compute the $R_n^2$ statistic (Section 2.5) to compare CBA 
# with a reduced model with three variables:
#
 CBA.red = update(CBA, .~.-x3-x4)
 np = 5       # number of parameters of the full model
 nq = 3       # number of parameters of the reduced model
 CVR = covar(CBA)
 CFF = coef(CBA)
 K22 = CVR[(nq+1):np,(nq+1):np]
 cff = as.matrix(CFF[(nq+1):np])
 Rn2 = t(cff)\%*\%solve(K22)\%*\%cff
 Rn2
}

\keyword{ package }

