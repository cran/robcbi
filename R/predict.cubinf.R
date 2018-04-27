predict.cubinf <-
function(object, newdata=NULL, type = c("link", "response", "terms"), se.fit = FALSE,
        terms = labels(object$terms), ...)
{
        type <- match.arg(type)
        family <- object$family
        if (is.character(family)) family <- get(family,mode="function",envir=parent.frame())
        if (is.function(family)) family <- family()
        if(!se.fit) {
#No standard errors
                if(missing(newdata)) {
				   pred <- switch(type, link = object$linear.predictors,
                           response = object$fitted,
                           terms = predict.lm(object, se.fit = se.fit, 
                                   scale = 1, type = "terms", terms = terms))
                }
                else {
				   pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                           type = ifelse(type == "link", "response", type), 
                           terms = terms, na.action = na.pass)
				   switch(type, response = {pred <- family$linkinv(pred)},
                                link = {
                                  type <- "response"
                                  NextMethod("predict")
                                }
                                ,
                                NextMethod("predict"))}
        }
        else {
                pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                        type = ifelse(type == "link", "response", type), 
                        terms = terms, na.action = na.pass)
				fit <- pred$fit
                se.fit <- pred$se.fit
                switch(type, response = {
#                      pred <- NextMethod("predict")
                       pred$fit <- family$linkinv(fit)
                       pred$se.fit <- se.fit/abs(family$deriv(pred$fit))
                        }
                        ,
                        link = {
                                type <- "response"
                                NextMethod("predict")
                        }
                        ,
                        NextMethod("predict"))
        }
	pred
}
