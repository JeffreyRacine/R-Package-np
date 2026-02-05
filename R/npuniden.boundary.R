npuniden.boundary <- function(X=NULL,
                              Y=NULL,
                              h=NULL,
                              a=min(X),
                              b=max(X),
                              bwmethod=c("cv.ls","cv.ml"),
                              cv=c("grid-hybrid","numeric"),
                              grid=NULL,
                              kertype=c("gaussian1","gaussian2",
                                  "beta1","beta2",
                                  "fb","fbl","fbu",
                                  "rigaussian","gamma"),
                              nmulti=5,
                              proper=FALSE) {
    kertype <- match.arg(kertype)
    cv <- match.arg(cv)
    bwmethod <- match.arg(bwmethod)
    if(!is.null(grid) && any(grid<=0)) stop(" the grid vector must contain positive values")
    if(is.null(X)) stop("you must pass a vector X")
    if(kertype=="gamma" || kertype=="rigaussian") b <- Inf
    if(kertype=="fbl") b <- Inf
    if(kertype=="fbu") a <- -Inf
    if(a>=b) stop("a must be less than b")
    if(any(X<a)) stop("X must be >= a")
    if(any(X>b)) stop("X must be <= b")
    if(!is.null(Y) && any(Y<a)) stop("Y must be >= a")
    if(!is.null(Y) && any(Y>b)) stop("Y must be <= b")
    if(is.null(Y)) Y <- X
    if(!is.null(h) && h <= 0) stop("bandwidth h must be positive")
    if(nmulti < 1) stop("number of multistarts nmulti must be positive")
    if(!is.logical(proper)) stop("proper must be either TRUE or FALSE")
    if(kertype=="gaussian2" && (!is.finite(a) || !is.finite(b))) stop("finite bounds are required for kertype gaussian2")
    h.opt <- NULL
    if(kertype=="gaussian1") {
        ## Gaussian reweighted boundary kernel function (bias of O(h))
        kernel <- function(x,X,h,a=0,b=1) {
            dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h)))
        }
    } else if(kertype=="gaussian2") {
        ## Gaussian reweighted second-order boundary kernel function
        ## (bias of O(h^2)). Instability surfaces for extremely large
        ## bandwidths relative to range of the data, so we shrink to
        ## the uniform when h exceeds 10,000 times the range (b-a)
        kernel <- function(x,X,h,a=0,b=1) {
            z <- (x-X)/h
            z.a <- (a-x)/h
            z.b <- (b-x)/h
            pnorm.zb.m.pnorm.za <- (pnorm(z.b)-pnorm(z.a))
            mu.1 <- (dnorm(z.a)-dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            mu.2 <- 1+(z.a*dnorm(z.a)-z.b*dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            mu.3 <- ((z.a**2+2)*dnorm(z.a)-(z.b**2+2)*dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            aa <- mu.3/(mu.3-mu.1*mu.2)
            bb <- -mu.1/(mu.3-mu.1*mu.2)
            if((b-a)/h > 1e-04) {
                (aa+bb*z**2)*dnorm(z)/(h*pnorm.zb.m.pnorm.za)
            } else {
                rep(1/(b-a),length(X))
            }
        }
    } else if(kertype=="beta1") {
        ## Chen (1999), Beta 1 kernel function (bias of O(h), function
        ## of f' and f'', no division by h), need to rescale to
        ## integrate to 1 on [a,b]
        kernel <- function(x,X,h,a=0,b=1) {
            X <- (X-a)/(b-a)
            x <- (x-a)/(b-a)
            dbeta(X,x/h+1,(1-x)/h+1)/(b-a)
        }
    } else if(kertype=="beta2") {
        ## Chen (1999), Beta 2 kernel function (bias of O(h), function
        ## of f'' only, no division by h), need to rescale to
        ## integrate to 1 on [a,b]
        rho <- function(x,h) {2*h**2+2.5-sqrt(4*h**4+6*h**2+2.25-x**2-x/h)}
        kernel <- function(x,X,h,a=0,b=1) {
            X <- (X-a)/(b-a)
            x <- (x-a)/(b-a)
            if(x < 2*h && h < (b-a)) {
                dbeta(X,rho(x,h),(1-x)/h)/(b-a)
            } else if((2*h <= x && x <= 1-2*h) || h >= (b-a)) {
                dbeta(X,x/h,(1-x)/h)/(b-a)
            } else if(x > 1-2*h && h < (b-a)) {
                dbeta(X,x/h,rho(1-x,h))/(b-a)
            }
        }
    } else if(kertype=="gamma") {
        ## Gamma kernel function for x in [a,Inf]
        kernel <- function(x,X,h,a=0,b=1) {
            ## No division by h, rescale to lie in [0,Inf], b is a
            ## dummy, not used but needed to avoid warning about
            ## function kernel having different named arguments
            X <- X-a
            x <- x-a
            dgamma(X,x/h+1,1/h)
        }
    } else if(kertype=="rigaussian") {
        ## Reverse inverse Gaussian for x in [a,Inf]
        kernel <- function(x,X,h,a=0,b=1) {
            ## No division by h, rescale to lie in [0,Inf], b is a
            ## dummy, not used but needed to avoid warning about
            ## function kernel having different named arguments
            X <- X - a
            x <- x - a
            x.res <- sqrt(x**2+h*x)
            k <- exp(-x.res/(2*h)*(X/x.res+x.res/X-2))/sqrt(2*pi*h*X)
            k[is.nan(k)] <- 0
            k
        }
    } else if(kertype=="fb") {
        ## Floating boundary kernel (Scott (1992), Page 46), left and
        ## right bound, truncated biweight in interior
        kernel <- function(x,X,h,a=0,b=1) {
            t <- (X-x)/h
            if(x < a+h && h < (b-a)) {
                c <- (a-x)/h
                ifelse(c <= t & t <= 2+c,.75*(c+1-1.25*(1+2*c)*(t-c)^2)*(t-(c+2))^2,0)/h
            } else if((a+h <= x && x <= b-h) || h >= (b-a)) {
                z.a <- (a-x)/h
                z.b <- (b-x)/h  
                rw <- (3*(z.b^5-z.a^5)-10*(z.b^3-z.a^3)+15*(z.b-z.a))/16
                rw[rw>1] <- 1
                ifelse(abs(t)<1,(15/16)*(1-t**2)**2/(h*rw),0)
            } else if(x > b-h && h < (b-a)) {
                c <- (b-x)/h
                ifelse(c-2 <= t & t <= c,.75*(1-c+1.25*(-1+2*c)*(t-c)^2)*(t-(c-2))^2,0)/h
            }
        }
    } else if(kertype=="fbl") {
        ## Floating boundary kernel (Scott (1992), Page 46), left bound
        kernel <- function(x,X,h,a=0,b=1) {
            t <- (X-x)/h
            if(x < a+h) {
                c <- (a-x)/h
                ifelse(c <= t & t <= 2+c,.75*(c+1-1.25*(1+2*c)*(t-c)^2)*(t-(c+2))^2,0)/h
            } else {
                ifelse(abs(t)<1,(15/16)*(1-t**2)**2/h,0)
            }
        }
    } else if(kertype=="fbu") {
        kernel <- function(x,X,h,a=0,b=1) {
            ## Floating boundary kernel (Scott (1992), Page 46), right bound
            t <- (X-x)/h
            if(x <= b-h) {
                ifelse(abs(t)<1,(15/16)*(1-t**2)**2/h,0)
            } else {
                c <- (b-x)/h
                ifelse(c-2 <= t & t <= c,.75*(1-c+1.25*(-1+2*c)*(t-c)^2)*(t-(c-2))^2,0)/h
            }
        }
    }
    int.kernel.squared <- function(X,h,a=a,b=b) {
        ## Use numeric integration to compute Kappa, the integral of
        ## the square of the kernel function needed for the asymptotic
        ## standard error of the density estimate seq(a,b) will barf
        ## on -Inf or Inf, trap these cases and use extendrange
        if(is.finite(a) && is.finite(b)) X.seq <- seq(a,b,length=1000)
        if(is.finite(a) && !is.finite(b)) X.seq <- seq(a,extendrange(X,f=10)[2],length=1000)
        if(!is.finite(a) && is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],b,length=1000)
        if(!is.finite(a) && !is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],extendrange(X,f=10)[2],length=1000)
        sapply(1:length(X),function(i){integrate.trapezoidal(X.seq,h*kernel(X[i],X.seq,h,a,b)**2)[length(X.seq)]})
    }
    fhat <- function(X,Y,h,a=0,b=1,proper=FALSE) {
        f <- sapply(1:length(Y),function(i){mean(kernel(Y[i],X,h,a,b))})
        if(proper) {
            if(is.finite(a) && is.finite(b)) X.seq <- seq(a,b,length=1000)
            if(is.finite(a) && !is.finite(b)) X.seq <- seq(a,extendrange(X,f=10)[2],length=1000)
            if(!is.finite(a) && is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],b,length=1000)
            if(!is.finite(a) && !is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],extendrange(X,f=10)[2],length=1000)
            f.seq <- sapply(1:length(X.seq),function(i){mean(kernel(X.seq[i],X,h,a,b))})
            if(any(f.seq<0)) {
                f <- f - min(f.seq)
                f.seq <- f.seq - min(f.seq)
            }
            int.f.seq <- integrate.trapezoidal(X.seq,f.seq)[length(X.seq)]
            f <- f/int.f.seq
        }
        return(f)
    }
    Fhat <- function(Y,f,a,b,proper=FALSE) {
        ## Numerical integration of f, check for aberrant values, if
        ## on range of data ensure F\in[0,1], if not make sure value
        ## is proper (negative boundary kernel functions can cause
        ## unwanted artifacts)
        f[is.na(f)] <- 0
        F <- integrate.trapezoidal(Y,f)
        if(proper) {
            if(min(Y)==a && max(Y)==b) {
                F <- (F-min(F))/(max(F)-min(F))
            } else {
                if(min(F)<0) F <- F+min(F)
                if(max(F)>1) F <- F/max(F)
            }
        }
        F
    }
    fhat.loo <- function(X,h,a=0,b=1) {
        sapply(1:length(X),function(i){mean(kernel(X[i],X[-i],h,a,b))})
    }
    if(bwmethod=="cv.ml") {
        ## Likelihood cross-validation function (maximizing)
        fnscale <- list(fnscale = -1)
        cv.function <- function(h,X,a=0,b=1) {
            f.loo <- fhat.loo(X,h,a,b)
            return(sum(log(ifelse(f.loo > 0 & is.finite(f.loo), f.loo, .Machine$double.xmin))))
        }
    } else {
        ## Least-squares cross-validation function (minimizing)
        fnscale <- list(fnscale = 1) 
        cv.function <- function(h,X,a=0,b=1) {
            cv.ls <- (integrate.trapezoidal(X,fhat(X,X,h,a,b)**2)[order(X)])[length(X)]-2*mean(fhat.loo(X,h,a,b))
            ifelse(is.finite(cv.ls),cv.ls,sqrt(sqrt(.Machine$double.xmax)))
        }
    }
    ## Grid search and then numeric optimization search (no
    ## multistarting, but sound starting point always used for
    ## subsequent refinement by optim)
    if(is.null(h) && cv == "grid-hybrid") {
        ## First establish a sound starting value using grid search,
        ## then use that starting value for numeric search
        if(is.null(grid)) {
            rob.spread <- c(sd(X),IQR(X)/1.349)
            rob.spread <- min(rob.spread[rob.spread>0])
            constant <- rob.spread*length(X)**(-0.2)
            h.vec <- c(seq(0.25,1.75,length=10),2^(1:25))*constant
            cv.vec <- sapply(1:length(h.vec),function(i){cv.function(h.vec[i],X,a,b)})
            foo <- optim(h.vec[ifelse(bwmethod=="cv.ml",which.max(cv.vec),which.min(cv.vec))],
                         cv.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=ifelse(kertype=="beta2",(b-a)/4,Inf),
                         control = fnscale,
                         X=X,
                         a=a,
                         b=b)
            h.opt <- foo$par
            cv.opt <- foo$value
        } else {
            cv.vec <- sapply(1:length(grid),function(i){cv.function(grid[i],X,a,b)})
            foo <- optim(grid[ifelse(bwmethod=="cv.ml",which.max(cv.vec),which.min(cv.vec))],
                         cv.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=ifelse(kertype=="beta2",(b-a)/4,Inf),
                         control = fnscale,
                         X=X,
                         a=a,
                         b=b)
            h.opt <- foo$par
            cv.opt <- foo$value
        }
    }
    if(is.null(h.opt)) {
        ## Manual inputted bandwidth
        f <- fhat(X,Y,h,a,b,proper=proper)
        ## Numerical integration via the trapezoidal rule
        F <- Fhat(Y,f,a,b,proper=proper)
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(Y,h,a,b)/(h*length(f)))),
                    sd.F=sqrt(abs(F*(1-F)/length(F))),
                    h=h))
    } else {
        ## Search bandwidth
        f <- fhat(X,Y,h.opt,a,b,proper=proper)
        ## Numerical integration via the trapezoidal rule
        F <- Fhat(Y,f,a,b,proper=proper)
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(Y,h.opt,a,b)/(h.opt*length(f)))),
                    sd.F=sqrt(abs(F*(1-F)/length(F))),
                    h=h.opt,
                    nmulti=nmulti,
                    cv.opt=cv.opt))
    }
}
