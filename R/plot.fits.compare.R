plot.fits.compare <-
function(x, xplots = FALSE, ..., ask = TRUE)
{
#
#
#   Save and set par() values
#
	lgx <- length(x)
	oldpar <- par(no.readonly=TRUE)
	on.exit(par(oldpar))
	oldmar <- c(3.6,3.6,2.1,1.1) #par()$mar
	newmar <- oldmar
#	newmar[3] <- oldmar[3] + lgx + 2
	oldmfcol <- par()$mfcol
	par(mfcol = c(1, 1))
	par(mar = newmar)
	oldask <- par()$ask
 	par(ask = ask)

#
# frac is the fraction to eliminate in the vertical direction (for
#              the legend)
# frax is the horizontal fraction to eliminate
#
	frac <- (lgx + 2)/16
	frax <- 2/16
	xnames <- rep("", lgx)
	if(xplots) {
		aname <- function(x)
		{
			names(x$coefficients)
		}
		unames <- unique(unlist(sapply(x, aname)))
		k <- match("(Intercept)", unames, nomatch = 0)
		if(k)
			unames <- unames[ - k]
		lnames <- length(unames)
		xmin <- rep(Inf, lnames)
		xmax <- rep( - Inf, lnames)
	}
#
#
# set up the lists we will need below
#
	form <- vector("list", lgx)
	f <- vector("list", lgx)
	r <- vector("list", lgx)
	den <- vector("list", lgx)
	y <- vector("list", lgx)
	q <- vector("list", lgx)
	yname <- vector("list", lgx)
	fname <- vector("list", lgx)
	xscal <- rep(0, lgx)
	if(xplots) {
		indep <- vector("list", lgx)
	}
	denmax <- 0
	denmin <- 1e+20
#
#
# Gather the data into the lists
#
	for(i in 1:lgx) {
                xx <- class(x[[i]])[1]
                if (xx=="glm" || xx=="cubinf") 
                {tmp <-  x[[i]]$dispersion
                 if (!is.null(tmp)) xscal[i] <- tmp
                } else {
                 tmp <- x[[i]]$scale
    	         if (!is.null(tmp)) xscal[i] <- tmp}
                rsi <- residuals(x[[i]])
                n.na <- sum(is.na(rsi))
                if (n.na != 0) rsi <- rsi[!is.na(rsi)]
                if(is.null(tmp))
		xscal[i] <- sum(rsi^2/(x[[i]]$df.residual-n.na))
		xnames[i] <- x[[i]]$name
		n <- length(rsi) + 0.5
		form[[i]] <- formula(x[[i]])
		f[[i]] <- predict(x[[i]])
		r[[i]] <- rsi
		den[[i]] <- density(rsi)
		denmax <- max(denmax, den[[i]]$y)
		denmin <- min(denmin, den[[i]]$y)
		q[[i]] <- qnorm(ppoints(length(r[[i]])))
		y[[i]] <- f[[i]] + r[[i]]
		yname[[i]] <- deparse(form[[i]][[2]])
 		fname[[i]] <- paste("Fitted :", deparse(form[[i]][[3]]),
                    collapse = " ")
		if(xplots) {
			xx <- model.matrix(x[[i]])
			indep[[i]] <- xx
			for(j in 1:lnames) {
				k <- match(unames[j], dimnames(xx)[[2]], 
				  nomatch = 0)
				if(k) {
				  xmax[j] <- max(xmax[j], xx[, k])
				  xmin[j] <- min(xmin[j], xx[, k])
				}
			}
		}
	}
#
#
# Define some functions for sapply
#
	minabs <- function(x)
	{
		min(abs(x), na.rm = TRUE)
	}
	maxabs <- function(x)
	{
		max(abs(x), na.rm = TRUE)
	}
	choose <- function(i)
	{
		if(i == 1) 16
		else i
	}
#
#
# Get the minimums and maximums
#
	fmax <- max(sapply(f, max, na.rm = TRUE))
	fmin <- min(sapply(f, min, na.rm = TRUE))
	ymax <- max(sapply(y, max, na.rm = TRUE))
	ymin <- min(sapply(y, min, na.rm = TRUE))
	rmin <- min(sapply(r, min, na.rm = TRUE))
	rmax <- max(sapply(r, max, na.rm = TRUE))
	qmin <- min(sapply(q, min, na.rm = TRUE))
	qmax <- max(sapply(q, max, na.rm = TRUE))#
#
# Histograms and qq plots
#
	par(mar = oldmar)
	par(mfcol = c(1, 1))
#	fraq <- 1-frac+0.05
#	am <- matrix(c(0,1,fraq,1,0,1,0,fraq),byrow=TRUE,ncol=4)
#	split.screen(am)
	split.screen(c(1,2))
	split.screen(c(lgx,1),1)
	rrmin <- rmin - 0.07*(rmax-rmin)
	rrmax <- rmax + 0.07*(rmax-rmin)
	zz <- hist(unlist(list(r[[1]],c(rrmin,rrmax))), plot=FALSE, nclass = 10)
	breaks <- zz$breaks
	iscreen <- 3
	for(i in 1:lgx) {
		screen(iscreen)
        par(mar=c(3.6,3.6,2.1,1.1), cex=0.7, tck=-0.01, mgp=c(1.5,0.4,0)) 
		iscreen <- iscreen + 1
 	    hist(r[[i]], probability = TRUE, xlab = "residual", 
 	        ylab = "density", main = xnames[i],  
#            ylim = c( denmin, denmax), 
            density=-1,breaks = breaks)
	}

	split.screen(c(lgx,1),2)
	par(xpd = FALSE)
	for(i in 1:lgx) {
#		if(i != 1) par(new = TRUE)
        screen(iscreen)
        par(mar=c(3.1,3.1,2.1,1.1), cex=0.7, tck=-0.01, mgp=c(1.5,0.4,0)) 
		qqnorm(r[[i]], ylab = "Residual", main = xnames[i], xlim = c(qmin, 
			qmax), ylim = c(rmin, rmax), pch=choose(i) )
		iscreen <- iscreen + 1
	}
	close.screen(all.screens=TRUE)
	par(mfrow = c(1, 1))
	par(mar=newmar, ask=ask)
	par(usr = c(0, 1, 0, 1))
#
#
# Fitted versus observed
#
#    xl <- 0.05-frax; yl <- 1+frac
#	if (xl<0) xl <- 0.05; if(yl>1) yl <- 0.995
#	par(mar = newmar)
#	par(mfcol = c(1, 1))
    par(cex=0.7, tck=-0.01, mgp=c(1.5,0.4,0))  
	plot(f[[1]], y[[1]], xlab = "Fitted Response", ylab = 
		"Observed Response", pch = 16, xlim = c(fmin, fmax), ylim = c(
		ymin, ymax), ...)
	legend("topleft", legend = xnames, lty =
		1:lgx, pch=c(16,2:lgx))
	abline(0, 1, lty = 2)
	for(i in 2:lgx) {
		points(f[[i]], y[[i]], pch = i, ...)
	}
#
# Fitted versus residuals
#
#    par(mar=c(3.1,3.1,2.1,1.1), cex=0.6, tck=-0.02, mgp=c(1.5,0.3,0)) 
#	legend(fmin - frax * (fmax - fmin), ymax + frac * (ymax - ymin), legend
#		 = xnames, pch = c(16, 2:lgx))
	plot(f[[1]], r[[1]], xlab = "Fitted Response", ylab = "Residual", pch
		 = 16, xlim = c(fmin, fmax), ylim = c(rmin, rmax),...)
	legend("topleft", legend = xnames, pch = c(16, 2:lgx))
	for(i in 2:lgx) {
		points(f[[i]], r[[i]], pch = i, ...)
	}		 
#
#
# x plots
#
	par(mfrow = c(1, 1))
	par(mar = newmar)
	if(xplots) {
		for(j in 1:length(unames)) {
			plt <- TRUE
			xin <- 1:lgx
			for(i in 1:lgx) {
				k <- match(unames[j], dimnames(indep[[i]])[[2]],
				  nomatch = 0)
				if(k) {
				  if(plt) {
				    plt <- FALSE    
                    par(mar=c(3.6,3.6,2.1,1.1), cex=0.7, tck=-0.02, mgp=c(1.5,0.4,0)) 
				    plot(indep[[i]][, k], r[[i]], xlab = unames[
				      j], ylab = "Residual", pch = choose(i), 
				      xlim = c(xmin[j], xmax[j]), ylim = c(rmin,
				      rmax), ...)
				  }
				  else points(indep[[i]][, k], r[[i]], pch = 
				      choose(i), ...)
				}
				else xin[i] <-  - xin[i]
			}
#			legend(xmin[j] - frax * (xmax[j] - xmin[j]), rmax + 
#				frac * (rmax - rmin), legend = xnames[xin], 
#				pch = c(16, 2:length(xin)))
			legend("topleft", legend = xnames[xin], 
				pch = c(16, 2:length(xin)))

		}
	}
	invisible()
}
