plot.cubinf <-
function(x, residuals = NULL, smooths = FALSE, rugplot = FALSE, id.n = 0, 
	ask = TRUE, num=0, ...)
{
        glm.obj <- x
	Residuals <- resid(glm.obj, type = "deviance")
	if(!is.null(residuals)) {
		if(length(residuals) == 1 && residuals)
			residuals <- Residuals
		else Residuals <- residuals
	}
	fits <- fitted(glm.obj)
	preds <- predict(glm.obj)
	response <- glm.obj$y
	form <- formula(glm.obj)
	response.name <- deparse(form[[2]])
	model <- deparse(form[[3]])
	add.ons <- function(x, y, smooths = TRUE, rugplot = TRUE, id.n = 3)
	{
		if(smooths) {
			prediction <- loess.smooth(x, y, span = 1, degree
				 = 1)
			lines(prediction)
		}
		if(rugplot) {
			jx <- jitter(x[!is.na(x)])
			xlim <- range(jx)
			rug(jx)
		}
		if(id.n) {
# Identify id.n greatest y-values (in absolute value)
			n <- length(y)
			oy <- order(abs(y))
			which <- oy[(n - id.n + 1):n]
			text(x[which], y[which], as.character(which), adj
				 = -0.05)
		}
	}
	choices <- c("All", "Residuals vs Fitted Values", 
		"Sqrt of abs(Residuals) vs Predictions", 
		"Response vs Fitted Values", 
		"Normal QQplot of Pearson Residuals",
                "Residuals vs Weights")
	choices <- substring(choices, 1, 40)	#truncate long names
	tmenu <- paste("plot:", choices)
	fit.lab <- paste("Fitted :", model, sep = " ")
	pred.lab <- paste("Predicted :", model, sep = " ") 
	pick <- 2
	ask.now <- ask
	while(pick <= (length(tmenu) + 2) ) {
            if (num>0) {pick <- num+1}
            else {
		if(ask.now)
			pick <- menu(tmenu, title = 
				"\nMake a plot selection (or 0 to exit):"
				) + 1}
		switch(pick,
			invisible(return(glm.obj)),
			{
# Plot all choices one by one
				ask.now <- FALSE
			}
			,
			{
# Residuals vs Fitted Values
				x <- fits
				y <- Residuals
				plot(x, y, xlab = fit.lab, ylab = if(
				  is.null(residuals)) 
				    "Deviance Residuals" else deparse(
				    substitute(residuals)), ...)
				abline(h = 0, lty = 2)
				add.ons(x, y, smooths = smooths, rugplot
				   = rugplot, id.n = id.n)
			}
			,
			{
# Sqrt of abs(Residuals) vs Predicted Values
				x <- preds
				y <- sqrt(abs(Residuals))
				plot(x, y, xlab = pred.lab, ylab = "Sqrt(|residuals|)", ...)
				add.ons(x, y, smooths = smooths, rugplot
				   = rugplot, id.n = id.n)
			}
			,
			{
# Response vs Fitted Values
				x <- fits
				y <- response
				xylims <- range(x, y)
				plot(x, y, xlab = fit.lab, ylab = 
				  response.name, xlim = xylims, ylim = 
				  xylims, ...)
				abline(0, 1, lty = 2)
				add.ons(x, y, smooths = smooths, rugplot
				   = rugplot, id.n = FALSE)
			}
			,
			{
# Normal QQplot of Residuals
				x <- resid(glm.obj, type = "pearson")
				qqnorm(x, ylab = "Pearson Residuals")
				QQline(x, lty = 2)
			}
                        ,
                        {
# Residuals vs Weights 
                             x <- weights(glm.obj)
                             y <- resid(glm.obj, type="deviance") 
                             plot(x,y, ylab="Residuals", xlab="Weights", ...)
			     add.ons(x, y, smooths=FALSE, rugplot=FALSE, id.n=id.n) 
                        }
			)
                if(num>0) return()
		if(!ask.now)
			pick <- pick + 1
		if(pick > (length(tmenu) + 2))  
			{ask.now <- ask
                         pick <- 1}
                     
	}
	invisible(glm.obj)
}
