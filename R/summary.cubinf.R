summary.cubinf <-
function(object, ...)
{
        coef <- object$coef
        resid <- object$residuals
        dresid <- residuals(object, "deviance")
        famname <- object$family
        family <- get(famname)
        family <- family()
        n <- length(resid)
        p <- object$rank
        if(is.null(p))
                p <- sum(!is.na(coef))
        if(!p) {
                warning("This model has zero rank --- no summary is provided")
                return(object)
        }
        nsingular <- length(coef) - p
        rdf <- object$df.resid
        if(is.null(rdf))
                rdf <- n - p
#       famname <- object$family#["name"]
        dispersion <- 1
        names(dispersion) <- famname

        covun <- object$cov
        var <- diag(covun)
        nas <- is.na(coef)
        cnames <- names(coef[!nas])
        coef <- matrix(rep(coef[!nas], 3), ncol = 3)
        dimnames(coef) <- list(cnames, c("Value", "Std. Error", "t value"))
        coef[, 2] <- sqrt(var)
        coef[, 3] <- coef[, 1]/coef[, 2]
        zl <- list(...)
        correlation <- TRUE
        if (length(zl) > 0) {
          pos <- grep("cor",names(zl), fixed=TRUE)
          if (length(pos)>0) {correlation <- as.logical(zl[pos])
                              zl <- zl[-pos] }
        }
        if(correlation) {
                cor <- covun
                for (i in 1:nrow(cor)) {
                  if (var[i]<1.e-10)
                    {str <- paste("Variance number",i,"smaller than 1.e-10",
                                "(set to 1.e-10)")
                     print(str)}
                  cor[i,1:i] <- cor[i,1:i]/sqrt(var[i]*var[1:i])
                  cor[1:i,i] <- cor[i,1:i]
                }
                dimnames(cor) <- list(cnames, cnames)
        }
        else cor <- NULL
        dimnames(covun) <- list(cnames, cnames)
        ocall <- object$call
#       if(!is.null(form <- object$formula)) {
#               if(is.null(ocall$formula))
#                       ocall <- match.call(get("glm"), ocall)
#               ocall$formula <- form
#       }
        ans <-  list(call = ocall, terms = object$terms, coefficients = coef, residuals=resid,
                dispersion = dispersion, df = c(p, rdf), deviance.resid = dresid, 
                family=family, cov.unscaled = covun, correlation = cor, deviance =
                deviance(object), null.deviance = object$null.deviance, iter =
                object$iter, nas=nas)
       class(ans) <-  c("summary.cubinf", "summary.glm")
       return(ans)
}
