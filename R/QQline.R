QQline <-
function(x, ...)
{
	data.quartiles <- quantile(x, c(0.25, 0.75))
	norm.quartiles <- qnorm(c(0.25, 0.75))
	b <- (data.quartiles[2] - data.quartiles[1])/(norm.quartiles[2] - 
		norm.quartiles[1])
	a <- data.quartiles[1] - norm.quartiles[1] * b
	abline(a, b, ...)
	ans <- as.vector(c(a, b))
	names(ans) <- c("intercept", "slope")
	invisible(ans)
}
