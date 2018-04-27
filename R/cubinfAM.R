cubinf.control <- function(tlo = 0.001, tua = 1.e-06, mxx = 30, mxt = 10, mxf = 10 ,
         ntm = 0, gma = 1, iug = 1, ipo = 1, ilg = 2, icn = 1, icv = 1, 
         ufact=0, cpar=1.5, null.dev = TRUE, ...)
         {list(tlo = tlo, tua = tua, mxx = mxx, mxt = mxt, mxf = mxf, 
          ntm = ntm, gma = gma, iug=iug, ipo=ipo, ilg=ilg, icn = icn,  
          icv = icv, ufact=ufact, cpar=cpar, null.dev=null.dev, ...)}

 
cubinf <- function(x, y, weights = NULL, start=NULL, etastart=NULL, mustart=NULL, 
         offset = NULL, family = binomial(), control = cubinf.control(...), intercept=FALSE, ...){
#
        x   <- as.matrix(x)
        n   <- nobs <- nrow(x)
        if (is.null(weights)) weights <- rep.int(1,n)
        if(any(weights != round(weights))) stop("Weights must be integer number")
        if(is.character(family)) family <- get(family,mode="function",envir=parent.frame())
        if(is.function(family)) family <- family()
        nn  <- family$family #["name"]
        if (nn == "gaussian") stop("Use lm(formula, ..., method=\"robust\") for the gaussian case")
        wi  <- weights
        ics <- 0
        if (nn == "binomial") {
           ics <- 2 
           if (is.matrix(y)) { 
             if (dim(y)[2]>2) stop("Only binomial response matrices (2 columns)")
             ni <- as.vector(y %*% c(1,1))
             y  <- y[,1]
             ly <- length(y)} 
           else {
             ics <- 1
             if (is.factor(y)) y <- as.numeric(y != levels(y)[1])
             else y <- as.vector(y)
             ly <- length(y) 
             ni <- rep(1,ly)
           } 
        }
        if (nn == "poisson") { 
           if (any(y < 0)) stop("Negative values not allowed for the Poisson family")
           ics <- 3; ly <- length(y)
           ni  <- rep(1,ly)}
        if (ics == 0) stop(paste(nn,": family not implemented for method='cubinf'", sep=" "))
        eta <- ci <- ai <- rsdev <- y 
        yy  <- y; nni <- ni 
        dn  <- dimnames(x)
        xn  <- dn[[2]]
        yn  <- dn[[1]]
#
# Data preprocessing for weights >= 0
        if (intercept & any(x[,1]!=1)) x <- cbind(1,x)
        p   <- ncol(x)
        EMPTY <- p==0
        ncov <- p*(p+1)/2        

        if (is.null(offset) || length(offset)==1) offset <- rep(0,ly) 
        zero <- wi == 0
        if (any(zero)) {
           pos <- !zero
           x0 <- x[zero, , drop=FALSE]
           y0 <- y[zero]
           x  <- x[pos, , drop=FALSE]
           y  <- y[pos]
           ni <- ni[pos]
           wi <- wi[pos]
           offset <- offset[pos]
        }
        nw <- length(wi) 
        ind <- rep.int(1:nw,wi)
        x  <- x[ind, , drop=FALSE]
        y  <- y[ind]
        offset <- offset[ind]
        ni <- ni[ind] 
#
# Initializations
        qrx <- qr(x,tol=1e-7,LAPACK=FALSE)[c("qr", "rank", "pivot", "qraux")]
        qrx$tol <- 1e-7
        rank <- qrx$rank
        piv <- 1:rank
        control <- do.call("cubinf.control", control)
        tlo <- control$tlo
        tua <- control$tua
        mxx <- control$mxx
        mxt <- control$mxt
        mxf <- control$mxf
        ntm <- control$ntm
        gma <- control$gma
        iug <- control$iug
        ipo <- control$ipo
        ilg <- control$ilg
        icn <- control$icn
        icv <- control$icv
        null.dev <- control$null.dev
        tmp <- control$singular.ok
        if (!is.null(tmp)) singular.ok <- tmp else singular.ok <- FALSE
        tmp <- control$qr.out
        if (!is.null(tmp)) qr.out <- tmp  else qr.out <- FALSE
        if (rank < p) {
                if (!singular.ok) stop(paste("x is singular, rank(x)=", rank))
                else {piv <- qrx$pivot[1:rank]
                      x <- x[, piv, drop=FALSE]
                      if (any(zero)) x0 <- x0[,piv,drop=FALSE]
                      xn <- dimnames(x)[[1]] }
        }
#       old <- comval()
        ufact <- control$ufact
        if (ufact==0) ufact <- 1.1
        upar <- ufact*sqrt(rank)
        cpar <- control$cpar
#       dev  <- control$dev
 
# Deviance for the model reduced to the constant term.

        if(null.dev) {Null.dev <- cubinf.null(x, y, ni, offset, ics, family, control)
                      if (nrow(x)==1 & intercept) return(list(deviance=Null.dev)) } 
        else Null.dev <- NULL 
#
#  Initial theta, A (A0) and c (c0)
#       tmp <- tempfile("cub"); zf <- file(tmp, open="wt")
        sink("tmpzzz", type="output") #Redirection of message 460 in RYWALG
        if (ufact >= 20) cpar <- 20*cpar
        z      <- gintac(x, y, ni, offset, icase=ics, tolt=10*tlo,tola=10*tlo, b=upar, c=cpar)
        theta0 <- z$theta[1:rank]; A0 <- z$a;  c0 <- z$ci
#
#  Initial cut off points a_i (wa)
        wa   <- upar/pmax(1.e-3,z$dist)
#
#  Initial covariance matrix of coefficient estimates
        vtheta <- as.vector(x%*%theta0)
#       z    <- gfedca(vtheta, c0, wa, ni, offset, ics)
#       zc   <- ktaskw(x, z$D, z$E, f=1/n, f1=1, iainv=0); covi <- zc$cov
        z    <- gfedcaAM(vtheta, c0, wa, ni, offset, ics)
        zc   <- covarAM(x, z$D, z$E,intercept=FALSE)$Cov
        covi <- zc[row(zc) <= col(zc)]
        iii <- cumsum(1:p)
        adiag <- round(covi[iii],4)
        jjj <- (iii)[adiag==0]
        if (length(jjj)>0) { covi[jjj] <- 1; cat("Info: initial cov re-defined\n",covi[1:ncov],"\n")} 
        if (icn != 1) {
           zc <- mchl(covi, rank)
           zi <- minv(zc$a, rank)
           covi <- mtt1(zi$r, rank)$b}
#
#  Final theta, A, c (ci) and a(wa)
        zf   <- gymain(x, y, ni, covi, A0, theta0, offset, b=upar, gam=gma, 
                tau=tua, icase=ics, iugl=iug, iopt=ipo, ialg=ilg, icnvt=icn, 
                icnva=icv, maxit=mxx, maxtt=mxt, maxta=mxf, maxtc=mxt, 
                nitmnt=ntm, nitmna=ntm, tol=tlo, tolt=10*tlo, tola=10*tlo, tolc=10*tlo)
        sink()
        nit  <- zf$nit; converged <- TRUE
        if (mxx > 1 && nit == mxx) {cat("\nWarning: Maximum number of iterations [mxx=", mxx,"] reached.\n")
                                   converged <- FALSE}
        coefs <- zf$theta
#
#  Deviance
#
        zd   <- glmdev(y, ni, zf$ci, zf$wa, zf$vtheta, offset, icase = ics)
#
#  Final covariance matrix of coeff. estimates (modified march 2018)
#       sink("tmpzzz", type="output") #Redirection of message 450 in KTASKW
#       z    <- gfedca(zf$vtheta, zf$ci, zf$wa, ni, offset, ics)
#       zc   <- ktaskw(x, z$D, z$E, f=1/n, f1=1, iainv=0)
        z    <- gfedcaAM(zf$vtheta, zf$ci, zf$wa, ni, offset, ics)
        zc   <- covarAM(x, z$D, z$E,intercept=FALSE)$Cov
#       sink()
        A    <- matrix(0, nrow=rank, ncol=rank)
        cov  <- zc[1:rank,1:rank]
        i2 <- 0
        for(i in 1:rank) {
                i1 <- i2 + 1
                i2 <- i1 + i - 1
                A[i, 1:i]  <- zf$a[i1:i2]
#               cov[i,1:i] <- zc$cov[i1:i2]
#               cov[1:i,i] <- zc$cov[i1:i2]
        }
        xn <- dimnames(x)[[2]]
        xn <- xn[piv]
        attributes(coefs) <- NULL
        attributes(A) <- NULL
        attributes(cov) <- NULL
        attr(A,"dim") <- c(rank,rank)
        attr(cov,"dim") <- c(rank,rank)
        names(coefs) <- xn
        dimnames(A)   <- list(xn,xn)
        dimnames(cov) <- list(xn,xn)
        asgn <- attr(x, "assign")
        ai <- zf$wa; ci <- zf$ci; rsdev <- zd$li-zd$sc
#       zl <- lrfctd(ics,y,ci,vtheta,offset,ai,ni,1,1,1)
#       Li <- zl$f0; li <- zl$f1; lip <- zl$f2 # rs <- zf$rs
#
# Compute eta, mu and residuals. 
        dni  <- c(ind[1],diff(ind)) 
        lll  <- dni!=0 
        iii  <- cumsum(dni[lll])
        lll  <- as.vector(1*lll)
        jjj  <- (1:n)*lll
        eta[iii] <- zf$vtheta[jjj]
        if(any(offset!=0)) offset[iii] <- offset[jjj]
        ci[iii] <- zf$ci[jjj]
        ai[iii] <- zf$wa[jjj] 
#       Li[iii] <- Li[jjj]
#       li[iii] <- li[jjj]
#       lip[iii] <- lip[jjj]
#       rs[iii] <- rs[jjj]
        rsdev[iii] <- zd$li[jjj]-zd$sc[jjj]
        dni  <- nni 
        ni   <- rep(1,length(eta))
        ni[iii] <- dni[jjj] 
        if (any(zero)) {
          eta[zero] <- as.vector(x0 %*% coefs)
          ci[zero]  <- 0
          ai[zero]  <- 0
#         Li[zero]  <- 0
#         li[zero]  <- 0
#         lip[zero] <- 0
#         rs[zero]  <- 0
          offset[zero] <- 0
          rsdev[zero] <- 0}   
        mu   <- family$linkinv(eta+offset)
        names(eta) <- yn
        if(rank < p) {
                coefs[piv[ - seq(rank)]] <- NA
                pasgn <- asgn
                newpos <- match(1:p, piv)
                names(newpos) <- xn
                for(j in names(asgn)) {
                        aj <- asgn[[j]]
                        aj <- aj[ok <- (nj <- newpos[aj]) <= rank]
                        if(length(aj)) {
                                asgn[[j]] <- aj
                                pasgn[[j]] <- nj[ok]
                        }
                        else asgn[[j]] <- pasgn[[j]] <- NULL
                }
                cnames <- xn[piv]
        }

        new.dev <- zd$dev 
#       new.dev <- sum(family$dev.resids(yy/ni, mu, w = rep(1.0, n))) 
        resp <- yy/nni - mu
        if (any(zero)) {resp[zero] <- 0}
        names(ai) <- yn
        df.residual <- ly - rank - sum(weights==0)
        fit <- list(coefficients = coefs, fitted.values=mu,
#                   effects = effects, weights=weights,
                    ci=ci, rank = rank, assign = asgn,
                    df.residual = df.residual, control=control)
        if(rank < p) { if(df.residual > 0) fit$assign.residual <- (rank + 1):n}
        if(qr.out) fit$qr <- qrx
        fit$A   <- A
        fit$ai  <- ai
        fit$cov <- cov
        fit$class <- "cubinf"
        fit$converged <- converged 
        if (!all(weights==1)) fit$prior.weights <- weights
        if (!all(ni==1)) fit$ni  <- ni 
        rsdev <- sign(y-mu)*sqrt(2*abs(rsdev)*weights) 
        attributes(zf$grad) <- NULL
        attributes(zf$hessnv) <- NULL
        residuals <- yy/nni - ci - mu
        attributes(residuals) <- NULL

# Restore common values for ROBETH
#       dfcomn(ipsi = old$ipsi, c = old$c, d = old$d, beta = old$bta)
        c(fit, list(family = family$family, ics=ics, linear.predictors = eta,
          deviance = new.dev, null.deviance=Null.dev, iter = nit, qr=qrx,
          y = yy/nni, contrasts = attr(x,"contrasts"), rsdev = rsdev,
          gradient = zf$grad, inv.hessian = zf$hessnv, residuals = residuals))
}

cubinf.null <-
function(x, y, ni, offset, ics, family, control)
{
        control <- do.call("cubinf.control", control)
        tlo <- control$tlo
        tua <- control$tua
        mxx <- control$mxx
        mxt <- control$mxt
        mxf <- control$mxf
        ntm <- control$ntm
        gma <- control$gma
        iug <- control$iug
        ipo <- control$ipo
        ilg <- control$ilg
        icn <- control$icn
        icv <- control$icv
        rank <- 1
        ufact <- control$ufact
        if (ufact==0) ufact <- 1.1
        upar <- ufact
        cpar <- control$cpar
        if (ufact >= 20) cpar <- 20*cpar
        ly   <- length(y)
        w    <- rep(1,ly)
        ai   <- ci <- rep(0,ly)
#       intl <- attr(x, "term.labels")
#       int  <- if(is.null(intl)) FALSE else as.logical(match(intl[1], c("(Int.)",
#                       "(Intercept)"), FALSE))
        linkinv <- family$linkinv
        dev.resids <- family$dev.resids
#       if(!int) {eta <- rep(0,ly); mu <- linkinv(eta+offset) 
#                 cval <- 0.5; if (ics==3) cval <- 1
#                 ci <- rep(cval,ly); ai <- rep(9999.,ly) } else {
           X    <- matrix(rep(1,ly),ncol=1)
           if (ufact >= 20) cpar <- 20*cpar
           sink("tmpzzz", type="output")#Redirection of message 460 in RYWALG
           z    <- gintac(X, y, ni, offset, icase = ics, tolt=10*tlo,
                   tola=10*tlo, b = upar, c = cpar)
#           sink(); unlink(tmp)     
           t0   <- z$theta[1]; A0 <- z$a; c0 <- z$ci
           wa   <- upar/pmax(1.e-3,z$dist)
           vtheta <- rep(t0,ly)
#          z    <- gfedca(vtheta, c0, wa, ni, offset, ics)
           z    <- gfedcaAM(vtheta, c0, wa, ni, offset, ics)
           zc   <- ktaskw(X, z$D, z$E, f=1/ly, f1=1, iainv=0)
           covi <- zc$cov
           if (icn != 1) covi <- 1/covi
           zf   <- gymain(X, y, ni, covi, A0, t0, offset, b=upar, gam=gma, 
                   tau=tua, icase=ics, iugl=iug, iopt=ipo, ialg=ilg, icnvt=icn,
                   icnva=icv, maxit=mxx, maxtt=mxt, maxta=mxf, maxtc=mxt, 
                   nitmnt=ntm, nitmna=ntm, tol=tlo, tolt=10*tlo, tola=10*tlo, 
                   tolc=10*tlo)
           sink()
           ai   <- zf$wa
           ci   <- zf$ci
           eta  <- zf$vtheta
#          mu   <- linkinv(eta+offset)
#       }
          zd   <- glmdev(y,ni,ci,ai,eta,offset=offset,icase=ics)
  zd$dev
}


gfedcaAM <- function(vtheta, c0, wa, ni, offset,ics, precision=0) {
psi <- function(x,c) { max(-c,min(x,c) ) }
vtheta1 <- vtheta+offset	
n <- length(vtheta)
E  <- D <- rep(0,n)
for ( i in 1:n) {
sumE <- sumD <- 0; j <- 0; termE <-  termD <- 100
if (ics==1 | ics==2) {probi   <- exp(vtheta1[i])/(1+exp(vtheta1[i]))
                      lambdai <- ni[i]*probi }
if (ics==3)           {lambdai <- exp(vtheta1[i])}
while (max(termE,termD) > precision) {
if (ics==1 | ics==2) {lpij  <- dbinom(j,ni[i], probi, log = TRUE) }
if (ics==3)          {lpij  <- dpois(j,lambdai, log = TRUE) }
tmpsi   <- psi( j-c0[i]-lambdai, c=wa[i])
termE   <- log( tmpsi^2 ) +  lpij
termE   <- exp(termE)
sumE    <- sumE + termE
tmpsi   <- tmpsi*(j-lambdai)
if (tmpsi>0) {
 termD   <- log(tmpsi) +  lpij
 termD   <- exp(termD)
 sumD    <- sumD + termD
} else {
 termD   <- tmpsi*exp(lpij)
 sumD    <- sumD + termD
 termD   <- abs(termD)
} 
j <- j+1
}
E[i] <- sumE; D[i] <- sumD}
list(D=D,E=E) }

covarAM <- function(X,D,E,intercept=TRUE,tol=sqrt(.Machine$double.eps)) {
XI    = X
if (intercept) XI = cbind(1,X)
S1   <- t(XI)%*%(D*XI)
S2   <- t(XI)%*%(E*XI)
#SI   = solve(S1)
Xsvd <- svd(S1)
Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
       SI <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
       SI <- array(0, dim(S1)[2L:1L])
    else 
      SI <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
         t(Xsvd$u[, Positive, drop = FALSE]))
Cov  <- SI%*%S2%*%SI
list(Cov=Cov) }


