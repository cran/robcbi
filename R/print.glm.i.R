print.glm.i <-
function(x, ...)
 {
#       x <- object
        atrass <- attr(x,"assign")
        if (!is.null(atrass)) attr(x,"assign") <- NULL
        atrsng <- attr(x,"singular")
        if (!is.null(atrsng)) attr(x,"singular") <- NULL
        class(x) <- NULL
        NextMethod(...)
        if (!is.null(atrass)) cat(attr(x,"assign"),"\n")
        if (!is.null(atrsng)) cat(attr(x,"singular"),"\n")
        invisible(x)
 }
