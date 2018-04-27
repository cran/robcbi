residuals.cubinf <-
function(object, type = c("deviance", "pearson", "response"), ...)
{
        type <- match.arg(type)
        y    <- object$y
        mu   <- object$fitted
        ics  <- object$ics
        if (ics==2) { ni <- object$ni; y <- y/ni} 
        resp <- y - mu
        family <- object$family
        if (is.character(family)) family <- get(family)
        family <- family()
        variance <- family$variance
        dev.resids <- family$dev.resids
        switch(type,
                pearson =  resp/sqrt(variance(mu)),
                deviance = {
                   w <- object$prior.w
                   if(is.null(w))  w <- rep(1, length(mu))
#                  dev <- object$control$dev 
#                  if (dev=="S+") 
#                    new.dev <- sum(dev.resids(y, mu, w))
#                  else {
                     new.dev <- object$rsdev 
#                  } 
                   new.dev
                }
                ,
                response = resp)
}
