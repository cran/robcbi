weights.glm <-
function(object)
{
 z <- object$weights
 class(z) <- "glm.i"
 z}
