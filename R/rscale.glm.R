rscale.glm <-
function(object)
 {
  famname <- object$family["name"]
  if(is.null(famname)) famname <- "Gaussian"
  z <- 1
  if (famname!="Poisson" && famname !="Binomial")
        { resid <- object$residuals
          n <- length(resid)
          p <- object$rank
          rdf <- object$df.resid
          if(is.null(rdf)) rdf <- n - p
          z <- sum(resid^2)/rdf}
 names(z) <- famname
 class(z) <- "glm.i"
 z}
