covar.cubinf <-
function(object)
{
 z   <- object$cov
 class(z) <- "cubinf.i"
 z}
