\name{cubinf.control}
\alias{cubinf.control}

\title{
Control parameters for the function cubinf
}
\description{
Allows the user to set parameters affecting the estimation of 
the discrete GLMs implemented in cubinf. Most control parameters 
are parameters of the ROBETH subroutine GYMAIN (Marazzi, 1993). 
}
\usage{
cubinf.control(tlo = 0.001, tua = 1e-06, mxx = 30, mxt = 10, mxf = 10, ntm = 0, gma = 1, 
               iug = 1, ipo = 1, ilg = 2, icn = 1, icv = 1, ufact = 0, cpar = 1.5,
               null.dev=TRUE, ...)
}

\arguments{
  \item{tlo}{Relative precision for the convergence criterion of the main algorithm (GYMAIN) 
    called by cubinf. The relative precision for the convergence criterion in the lower level 
    steps (theta-step, A-step and c-step) is '10*tlo'. 
}
  \item{tua}{Tolerance used for the determination of the pseudo-rank. 
}
  \item{mxx}{Maximum number of cycles for the main algorithm. 
}
  \item{mxt}{Maximum number of iterations for the theta-step. 
}
  \item{mxf}{Maximum number of iterations for the A-step. 
}
  \item{ntm}{Parameter to control iteration monitoring. When the number of iterations in the 
     theta-step reaches a multiple of 'ntm', the current parameter values as well as the 
     corresponding value of the objective function are printed. 
}
  \item{gma}{Relaxation factor for the theta-step. 
}
  \item{iug}{Parameter for the choice of the u-function in the A-step. 
    See Marazzi, 1993, for details. 
}
  \item{ipo}{Parameter for the choice of the steplength algorithm in the theta-step.
    If 'ipo=1', a quadratic comparison function is minimized.
    If 'ipo=2', the Goldstein-Armijo step length algorithm is used. 
}
  \item{ilg}{Parameter for the choice of the algorithm in the c-step.
    If 'ilg=1', the H-algorithm is used.
    If 'ilg=2', the W-algorithm is used.
}
  \item{icn}{Parameter for the choice of the convergence criterion for the theta-step and the main algorithm.
     If 'icn=1', convergence is assumed when the change in each coefficient is less than 
     the tolerance ('10*tlo') times an estimate of the coefficient variance. 
     See Marazzi (1993, p. 281), for the other options ('icn=2' and 'icn=3'). 
}
  \item{icv}{Parameter for the choice of the convergence criterion for the A-step. 
    If 'icv=1', convergence is assumed when the norm of the difference between two consecutive 
    values of A is less than the tolerance (10*tol).
    See Marazzi (1993, p.288 and p. 301), for another option ('icv=2'). 
}
  \item{ufact}{The tuning constant b is set equal to ufact*sqrt(p), where p is the dimension of 
    the observation vectors. The default value of b is 1.1*sqrt(p); this value is used when 'ufact=0' on input. 
}
  \item{cpar}{Parameter used in determining an initial value of theta (standard Mallows estimate, see Marazzi, 1993, p281). 
}
  \item{null.dev}{If 'null.dev=TRUE', the null deviance is computed. The null deviance is the deviance 
   of the model with no predictors. 
}
 \item{\dots }{Further named control arguments as singular.ok or qr.out used in the case where the
   x matrix is singular} 
}

\value{
List of control parameters. 
}
\references{
Marazzi, A. (1993).
Algorithms, Routines, and S-functions for robust Statistics.
Chapman and Hall, New York.
}

\seealso{
\code{\link{cubinf}}
}
\examples{
#To compute the classical estimates using cubinf, set:
control <- cubinf.control(ufact=300)
}

