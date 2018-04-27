\name{Breslow}
\alias{Breslow}
\docType{data}
\title{
Breslow Data
}
\description{
 Patients suffering from simple or complex partial seizures were randomized to receive 
 either the antiepileptic drug progabide or a placebo. At each of four successive
 postrandomization clinic visits, the number of seizures occuring over the previous 
 two weeks was reported
}
\usage{
data(Breslow)
}

\format{
 A data frame with 59 observations and the following 4 variables 
  \describe{
    \item{Trt}{the treatment: a factor with levels "placebo" and "progabide"}
    \item{sumY}{aninteger value, the sum of seizures during the 1st, 2nd, 3rd and 4th two week periods}
    \item{Age10}{Age divided by 10}
    \item{Base4}{the eight-week baseline seizure count divided by 4}
}}


\references{
Breslow, N. E., and Clayton, D. G. (1993), "Approximate Inference in Generalized Linear Mixed Models," 
Journal of the American Statistical Association, Vol. 88, No. 421, pp. 9-25. 

Thrall, P. F., and Vail, S. C. (1990), "Some Covariance Models for Longitudinal Count Data With 
Overdispersion," Biometrics, Vol. 46, pp. 657-671. 
}

\examples{
library(robcbi)
data(Breslow)
y  <- Breslow$sumY
x1 <- Breslow$Age10
x2 <- Breslow$Base4
x3 <- rep(0,length(y))
x3[Breslow$Trt=="progabide"] <- 1
}
\keyword{datasets}
