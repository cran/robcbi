print.fits.compare <-
function(x, digits = max(3, .Options$digits - 3), ...)
{
	names <- rep(" ", length(x))
	cat("\nCalls: \nName\n")
	nmod <- length(x)
	modsum <- vector("list", nmod)
        lgm <- rep(0,nmod)
	for(i in 1:nmod) {
		names[i] <- x[[i]]$name
                x[[i]]$lgm <- 0
                cl <- class(x[[i]])[1] 
		if(cl == "aov")
	           modsum[[i]] <- summary.lm(x[[i]])
		else if(cl == "glm" )
		   {modsum[[i]] <- summary.lm(x[[i]])
                    x[[i]]$lgm <- 1}
                else if(cl == "cubinf")
		   {modsum[[i]] <- summary.cubinf(x[[i]])
                    x[[i]]$lgm <- 2}
		else modsum[[i]] <- summary(x[[i]])
	}
   allgm <- sapply(x,"[[","lgm")
	names <- format(names)
	for(i in 1:nmod) {
		cat(names[i], "   ")
		dput(modsum[[i]]$call)	#                cat("\n")
	}
#        ii <<- 0
	rsid  <- function(x)
        {
#         ii <<- ii+1 
          if (x$lgm>=1) 
	   z <- residuals(x, type="response")
          else 
           z <- residuals(x)  
         z 
       	}
	rq <- t(sapply(lapply(x, rsid), quantile,na.rm=TRUE))
	dimnames(rq)[[1]] <- names
	rq <- cbind(rq, rq[, 1:4, drop = FALSE])
	dimnames(rq)[[2]] <- c("Min", "1Q", "Median", "3Q", "Max", "Nobs", 
		"Resid df", "Model Parameters", "Est. Parameters")
	aname <- function(x)
        {
         dimnames(x$coefficients)[1] 
       	}
        unames <- unique(unlist(sapply(modsum, aname)))
	lnames <- length(unames)
        Ncl <- 4
        if (any(allgm!=0)) Ncl <- 3
        dn2 <- dimnames(modsum[[1]]$coefficients)[[2]]
	coeff <- matrix(NA, ncol = Ncl, nrow = nmod * lnames)
	dimnames(coeff) <- list(as.list(paste(rep(unames, rep(nmod, length(
		unames))), names, sep = ":")), dn2[1:Ncl]) 
	for(i in 1:nmod) {
		j <- match(dimnames(modsum[[i]]$coefficients)[[1]], unames) - 1
                z <- modsum[[i]]$coefficients 
		coeff[j * nmod + i, 1:Ncl] <- z[,1:Ncl]
		rq[i, 7] <- modsum[[i]]$df[2]
		rq[i, 8] <- modsum[[i]]$df[1]
                if (allgm[i]==0) {
		  rq[i, 6] <- length(modsum[[i]]$residuals)
		  rq[i, 9] <- modsum[[i]]$df[3]}
                else {
                 rq[i,6] <- length(residuals(x[[i]])) 
                 z  <- x[[i]]$coef
                 nas <- is.na(z)
                 rq[i,9] <- length(z[!nas])}
	}
	cat("\n\nResidual Statistics\n")
	print(rq[, 1:5, drop = FALSE], digits = digits, ...)
	cat("\n\nNumber of Parameter in each Model\n")
	print(rq[, 6:9, drop = FALSE], digits = digits, ...)
	cat("\n\nCoefficients:\n")
	print(format(round(coeff, digits = digits)), quote = FALSE, ...)
	cat("\n\nResidual Scale Estimates:\n")
          for(i in 1:nmod) {
            if(allgm[i] >= 1) { 
                z <- rscale(x[[i]])
		cat(names[i], ":", format(signif(z, digits)),"\n")
            }
	    else 
		cat(names[i], ":", format(signif(modsum[[i]]$sigma, digits)), 
			"on", rq[i, 7], "degrees of freedom\n")
          }
          ii <- 0
	  for(i in 1:nmod) {
            if (allgm[i] == 0) {
             ii <- ii+1 
             if (ii==1) 
       cat("\nProportion of variation in response(s) explained by model(s):\n")
		cat(names[i], ":", format(signif(modsum[[i]]$r.squared, digits)
			), "\n")}
	  }
          ii <- 0
          for(i in 1:nmod) {
            if (allgm[i] == 0) {
              ii <- ii+1
              if (ii == 1) cat("\nF-statistics (NA if not available):\n")
              if (is.null(modsum[[i]]$fstatistic[1]))
                cat(names[i],": NA \n")
              else 
		cat(names[i], ":", format(signif(modsum[[1]]$fstatistic[1], 
			digits)), "on", modsum[[i]]$fstatistic[2], "and", 
			modsum[[i]]$fstatistic[3], "degrees of freedom, the ", 
			"p-value is", format(signif(1 - pf(modsum[[i]]$
			fstatistic[1], modsum[[i]]$fstatistic[2], modsum[[i]]$
			fstatistic[3]), digits)), "\n")
           }
	}
	cat("\nCorrelation of Coefficients:\n")
	for(i in 1:nmod) {
                if (allgm[i]>=1)
                 cor <- correl(x[[i]])
                else 
                 cor <- modsum[[i]]$correlation
		if(!is.null(cor)) {
			cat("\nModel = ", names[i], "\n")
			p <- dim(cor)[2]
			if(p > 1) {
				ll <- lower.tri(cor)
				cor[ll] <- format(round(cor[ll], digits))
				cor[!ll] <- ""
				print(cor[-1,  - p, drop = FALSE], quote = FALSE, 
				  digits = digits, ...)
			}
		}
	}
	invisible(modsum)
}

