correl <-
function(object, tl=1.e-10)
 {
  z    <- covar(object)
  cl   <- class(z) 
  zvar <- diag(z)
  for (i in 1:nrow(z))
  {
       if (zvar[i]<tl)
       {str <- paste("Variance number",i,"smaller than",tl,
                     "(set to",tl,")")
        print(str)}
        z[i,1:i] <- z[i,1:i]/sqrt(zvar[i]*zvar[1:i])
        z[1:i,i] <- z[i,1:i]
  }
  class(z) <- cl
  z}
