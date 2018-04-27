covar.glm <-
function(object)
{
 R <- object$R
 p <- object$rank
 if(p < max(dim(R))) R <- R[1:p, 1:p]
 rinv <- diag(p) 
 rinv <- backsolve(R, rinv)
 z <- rinv %*% t(rinv) 
 dn <- dimnames(R) 
 dimnames(z) <- dn
 class(z) <- "glm.i"
 z}
