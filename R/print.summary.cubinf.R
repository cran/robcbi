print.summary.cubinf <-
function(x, ...) # digits=max(3,.Options$digit-3, prefix = "")
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
        nas <- x$nas
        coef <- x$coef
        correl <- x$correl
        if(any(nas)) {
                nc <- length(nas)
                cnames <- names(nas)
                coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
                coef1[!nas,  ] <- coef
                coef <- coef1
                if(!is.null(correl)) {
                        correl1 <- matrix(NA, nc, nc, dimnames = list(cnames,
                                cnames))
                        correl1[!nas, !nas] <- correl
                        correl <- correl1
                }
        }
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

        dresid <- x$deviance.resid
        n <- length(dresid)
        rdf <- x$df[2]
        if(rdf > 5) {
                cat("Deviance Residuals:\n")
                rq <- quantile(as.vector(dresid), na.rm=TRUE)
                names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
                print(rq, digits = digits)
        }
        else if(rdf > 0) {
                cat("Deviance Residuals:\n")
                print(dresid, digits = digits)
        }
        if(any(nas))
                cat("\nCoefficients: (", sum(nas),
                        " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        print(coef, digits = digits)
        cat(paste("\n(Dispersion Parameter for", names(x$dispersion),
                "family taken to be", format(round(x$dispersion, digits)),
                ")\n"))
        int <- attr(x$terms, "intercept")
        if(is.null(int))
                int <- 1
        cat("\n    Null Deviance:", format(round(x$null.deviance, digits)),
                "on", n - int, "degrees of freedom\n")
        cat("\nResidual Deviance:", format(round(x$deviance, digits)), "on",
                round(rdf, digits), "degrees of freedom\n")
        cat("\nNumber of Iterations:", format(trunc(x$iter)),
                "\n")
        if(!is.null(correl)) {
                p <- dim(correl)[2]
                if(p > 1) {
                        cat("\nCorrelation of Coefficients:\n")
                        ll <- lower.tri(correl)
                        correl[ll] <- format(round(correl[ll], digits))
                        correl[!ll] <- ""
                        print(correl[-1,  - p, drop = FALSE], quote = FALSE, digits =
                                digits)
                 }
        }
        invisible(NULL)
}
