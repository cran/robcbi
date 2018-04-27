\name{cubinf}
\alias{cubinf}

\title{
M-estimates for discrete Generalized Linear Models
}
\description{
Conditionally unbiased bounded influence estimates as described in Kuensch et al. 
(1989) in three special cases of the generalized linear model: Bernoulli, Binomial,
and Poisson distributed responses. Produces an object of class "cubinf" which is a 
robust generalized linear fit of the data. 
}
\usage{
cubinf(x, y, weights = NULL,  start=NULL, etastart=NULL, mustart=NULL,
       offset = NULL, family = binomial(), control = cubinf.control(...),  
       intercept = FALSE, ...)
}

\arguments{
  \item{x}{Vector or matrix of explanatory variable(s).
           Columns represent variables and rows represent observations.
}
  \item{y}{Vector of observations on the response variable. In the case of Binomial responses, 
           y is a two column matrix:
           the 1st column contains the number of successes, 
           the 2nd column the number of failures.
           The Bernoulli case, is treated as a special Binomial case. However, the response y 
           is a categorical variable (not a matrix with two colums) with two levels.
}
  \item{weights}{Optional weigths for weighted regression. Component values must be 
        non negative integers. 
}
  \item{start}{Starting values for the parameters in the linear predictor. Not used but
        required for compatibility with the glm R function. 
}
  \item{etastart}{Starting values for the linear predictor. Not used but required for 
        compatibility with the glm R function. 
}
  \item{mustart}{Starting values for the vector of means. Not used but required for 
        compatibility with the glm R function.
}

  \item{offset}{Optional offset added to the linear predictor. 
}
  \item{family}{A family object. Only two options are available for cubinf: 
        'family=binomial()' and 'family=poisson()'. 
}

  \item{control}{A list of control parameters for the numerical algorithms. 
   See cubinf.control() for the possible control parameters and their default settings
}
  \item{intercept}{Logical flag: if TRUE, an intercept term is added to the model. 
}
\item{\dots }{Further named control arguments as singular.ok or qr.out used in the case 
 where the x matrix is singular} 
}

\value{
A list with the following components: 
  \item{coefficients}{Coefficient estimates. }
  \item{residuals}{Working residuals. }
  \item{rsdev}{Deviance residuals. }
  \item{fitted.values }{Fitted values. }
  \item{cov}{Estimated covariance matrix of the coefficients. }
  \item{rank }{The rank of the model matrix. }
  \item{df.residuals}{The degrees of freedom in the residuals. }
  \item{ci}{Vector (c) of the final bias corrections. }
  \item{A}{The final value of the matrix A. }
  \item{ai}{Vector with components a_i=ufact/|Ax_i| (where x_i^T 
    denotes the ith row of the model matrix) }
  \item{converged}{A logical value. FALSE if the maximum number of iterations was 
   reached.  }
  \item{control}{The control parameters used.  }
  \item{prior.weights}{The input vector w (when some of its components are different 
   from 1). }
  \item{family}{The family object used in the call to cubinf
    'ics=1' for the Bernoulli case.
    'ics=2' for the Binomial case.
    'ics=3' for the Poisson case. }
  \item{linear.predictors}{Components of the linear predictor (the model matrix 
   multiplied by the coefficient vector). }
  \item{iter}{The number of iterations required in the main algorithm.  }
  \item{y }{The coded value of the response. }
  \item{gradient }{Vector of the final unscaled negative gradient of the objective 
   function. }
  \item{inv.hessian }{Vector of the final inverse of the Hessian matrix in compact 
   storage mode. }
}
\details{
Initial values of the coefficients (theta), the matrix A and the bias correction c are 
computed using the ROBETH subroutine GINTAC (Marazzi, 1993). Then an initial covariance 
matrix (for the convergence criterion) is computed by means of the ROBETH subroutines 
GFEDCA and KTASKW. 
Finally, the main algorithm (subroutine GYMAIN) alternates between improving values of
- theta, for fixed A and c (theta-step, subroutine GYTSTP),
- c, for fixed theta and A (c-step, subroutine GYCSTP),
- A, for fixed theta and c (A-step, subroutine GYASTP).

For the different available options see the function cubinf.control. 
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
glm(..., method="cubinf"), \code{\link{cubinf.control}}
}
\examples{
  library(robcbi)
  y <- c(5,10,15,20,30,40,60,80,100)
  x <- matrix(
	c(0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1),
    nrow=9,byrow=FALSE)
  z <- cubinf(x,y, family=poisson, control=list(ufact=3.2), intercept=TRUE)
  z$iter
  z$coeff
  z <- cubinf(x,y, family=poisson, control=list(ufact=30), intercept=TRUE)
  z$iter
  z$coeff
}

\keyword{ stats }
\keyword{ robust }
