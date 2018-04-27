print.cubinf <-
function(x, ai=FALSE, ci=FALSE, A.mat=FALSE, ...)
{
        digits <- max(3, .Options$digits-3)
        zl <- list(...)
        if (length(zl) >0) {
          pos <- grep("digit", names(zl), fixed="TRUE")
          if (length(pos)>0) {
            digits <- as.numeric(zl[pos])
            zl <- zl[-pos]
          }
        }
#   cat("\nCall:  \n")
#    str(x$ccall)
        coef <- x$coef
        if(any(nas <- is.na(coef))) {
                if(is.null(names(coef))) names(coef) <- paste("b", 1:length(
                                coef), sep = "")        
                cat("\nCoefficients: (", sum(nas),
                        " not defined because of singularities)\n", sep = "")
        }
        else cat("\nCoefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)

        rank <- x$rank
        if(is.null(rank))
                rank <- sum(!nas)
        nobs <- length(x$residuals)
        rdf <- x$df.resid
        if(is.null(rdf))
                rdf <- nobs - rank
        cat("\nDegrees of Freedom:", nobs-1, "Total (i.e. Null);", rdf, "Residual\n")
        cat("Residual Deviance:", format(x$deviance), "\n")
        if (A.mat) {
          A <- x$A
          if(!is.null(A)) {
                  p <- dim(A)[2]
                  if(p > 1) {
                          cat("\nLower-triangular A matrix :\n")
                          ll <- row(A) >= col(A)
                          A  <- format(round(A, digits))
                          A[!ll] <- ""
                          print(A[drop = FALSE], digits=digits, quote = FALSE)
                  }
          }
        }
        if (ai) {
             cat("\nContants a_1,..., a_n :\n")
             print(x$ai, digits=digits, quote=FALSE)
        }
        if (ci) {
             cat("\nContants c_1,..., c_n :\n")
             print(x$ci, digits=digits, quote=FALSE)
        }
        invisible(x)
}
