\name{Finney}
\alias{Finney}
\docType{data}
\title{
Finney data from 'Annals of Eugenics' 1947
}
\description{
 Finney data over 39 observations on occurence or not of vaso-constriction.
 The data were obtained in a study of the effect of the rate and volume 
 of air inspired on a transcient vaso-constriction in the skin of the digits. 
 The R function, plotFdat, for plotting the variables is included in the list. 
}
\usage{
data(Finney)
}

\format{
 A list with the following components
  \describe{
    \item{Resp}{Occurence (Resp=1) or not (Resp=0) of vaso-constriction}
    \item{Vol}{Volume of air inspired}
    \item{Rate}{Observed Rate}
    \item{plotFdat}{R function for plotting the response and explanatory variables}
}}

\source{
Finney (1947), The truncated binomial distribution, Annals of Eugenics, 14: 319-328.
}

\references{
Kuensch, H.R., Stefanski L.A., Carroll R.J. (1989).
Conditionally unbiased bounded-influence estimation in general regression models,
with application to generalized linear models. 
Journal of the American Statistical Association, 84, 460-466.
}

\examples{
library(robcbi)
data(Finney)
Vol   <- Finney$Vol; Rate <- Finney$Rate; Resp  <- Finney$Resp
lVol  <- log(Vol);  lRate <- log(Rate)
}
\keyword{datasets}
