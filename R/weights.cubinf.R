weights.cubinf <-
function(object)
{
 res <- object$y -object$ci - object$fitted.values  
 z <- pmin(1,object$ai/abs(res)) 
 class(z) <- "cubinf.i"
 z}
