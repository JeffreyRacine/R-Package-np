npuniden.boundary <- function(X=NULL,
                              h=NULL,
                              a=0,
                              b=1,
                              kertype=c("gaussian1","gaussian2","beta1","beta2","gamma"),
                              cv=c("grid-hybrid","numeric"),
                              grid=NULL,
                              nmulti=5) {
    kertype <- match.arg(kertype)
    cv <- match.arg(cv)
    if(!is.null(grid) && any(grid<=0)) stop(" the grid vector must contain positive values")
    if(is.null(X)) stop("you must pass a vector X")
    if(a>=b) stop("a must be less than b")
    if(any(X<a)) stop("X must be >= a")
    if(any(X>b)) stop("X must be <= b")
    if(!is.null(h) && h <= 0) stop("bandwidth h must be positive")
    if(nmulti < 1) stop("number of multistarts nmulti must be positive")
    if(kertype=="gaussian2" && !is.finite(a) && !is.finite(b)) stop("at least one finite bound needed with kertype gaussian2")
    h.opt <- NULL
    if(kertype=="gaussian1") {
        ## Gaussian reweighted boundary kernel function (bias of O(h))
        kernel <- function(x,X,h,a=0,b=1) {
            dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h)))
        }
        kernel.int <- function(x,X,h,a=0,b=1) {
            z.a <- (a-x)/h
            z.b <- (b-x)/h            
            z <- pmax(pmin((x-X)/h,(b-x)/h),(a-x)/h)
            (pnorm(z)-pnorm(z.a))/(pnorm(z.b)-pnorm(z.a))
        }
    } else if(kertype=="gaussian2") {
        ## Gaussian reweighted second-order boundary kernel function
        ## (bias of O(h^2))
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
            (aa+bb*z**2)*dnorm(z)/(h*pnorm.zb.m.pnorm.za)
        }
        kernel.int <- function(x,X,h,a=0,b=1) {
            z.a <- (a-x)/h
            z.b <- (b-x)/h
            z <- pmax(pmin((x-X)/h,(b-x)/h),(a-x)/h)
            pnorm.zb.m.pnorm.za <- (pnorm(z.b)-pnorm(z.a))
            mu.1 <- (dnorm(z.a)-dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            mu.2 <- 1+(z.a*dnorm(z.a)-z.b*dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            mu.3 <- ((z.a**2+2)*dnorm(z.a)-(z.b**2+2)*dnorm(z.b))/(pnorm.zb.m.pnorm.za)
            aa <- mu.3/(mu.3-mu.1*mu.2)
            bb <- -mu.1/(mu.3-mu.1*mu.2)
            ((pnorm(z)-pnorm(z.a))*(aa+bb)+(dnorm(z.a)*z.a-dnorm(z)*z)*bb)/(pnorm.zb.m.pnorm.za)
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
        kernel.int <- function(x,X,h,a=0,b=1) {
            X <- (X-a)/(b-a)
            x <- (x-a)/(b-a)
            1-pbeta(X,x/h+1,(1-x)/h+1)/(b-a)
        }
    } else if(kertype=="beta2") {
        ## Chen (1999), Beta 2 kernel function (bias of O(h), function
        ## of f'' only, no division by h), need to rescale to
        ## integrate to 1 on [a,b]
        rho <- function(x,h) {2*h**2+2.5-sqrt(4*h**4+6*h**2+2.25-x**2-x/h)}
        kernel <- function(x,X,h,a=0,b=1) {
            X <- (X-a)/(b-a)
            x <- (x-a)/(b-a)
            if(x < 2*h) {
                dbeta(X,rho(x,h),(1-x)/h)/(b-a)
            } else if(2*h <= x & x <= 1-2*h) {
                dbeta(X,x/h,(1-x)/h)/(b-a)
            } else if(x > 1-2*h) {
                dbeta(X,x/h,rho(1-x,h))/(b-a)
            }
        }
        kernel.int <- function(x,X,h,a=0,b=1) {
            X <- (X-a)/(b-a)
            x <- (x-a)/(b-a)
            if(x < 2*h) {
                1-pbeta(X,rho(x,h),(1-x)/h)/(b-a)
            } else if(2*h <= x & x <= 1-2*h) {
                1-pbeta(X,x/h,(1-x)/h)/(b-a)
            } else if(x > 1-2*h) {
                1-pbeta(X,x/h,rho(1-x,h))/(b-a)
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
        kernel.int <- function(x,X,h,a=0,b=1) {
            ## No division by h, rescale to lie in [0,Inf], b is a
            ## dummy, not used but needed to avoid warning about
            ## function kernel having different named arguments
            X <- X-a
            x <- x-a
            1-pgamma(X,x/h+1,1/h)
        }
    }
    int.kernel.squared <- function(X,h,a=1,b=1) {
        ## Use numeric integration to compute Kappa, the integral of
        ## the square of the kernel function needed for the asymptotic
        ## standard error of the density estimate
        ## seq(a,b) will barf on -Inf or Inf, trap these cases and use extendrange
        if(is.finite(a) && is.finite(b)) X.seq <- seq(a,b,length=100)
        if(is.finite(a) && !is.finite(b)) X.seq <- seq(a,extendrange(X,f=10)[2],length=1000)
        if(!is.finite(a) && is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],b,length=1000)
        if(!is.finite(a) && !is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],extendrange(X,f=10)[2],length=1000)
        if(kertype=="gaussian1" || kertype=="gaussian2" || kertype=="epanechnikov") {
            ## Gaussian and Epanechnikov divide by h in kernel, need
            ## to adjust for integration by multiplying kernel squared
            ## by h (integral of kernel squared wrt X is h times
            ## integral of kernel squared wrt (x-X)/h, so kernel
            ## squared has 1/h^2 term, adjust by multiplying kernel
            ## squared by h)
            sapply(1:length(X),function(i){integrate.trapezoidal(X.seq,h*kernel(X[i],X.seq,h,a,b)**2)[length(X.seq)]})
        } else {
            sapply(1:length(X),function(i){integrate.trapezoidal(X.seq,h*kernel(X[i],X.seq,h,a,b)**2)[length(X.seq)]})
        }
    }
    fhat <- function(X,h,a=0,b=1) {
        sapply(1:length(X),function(i){mean(kernel(X[i],X,h,a,b))})
    }
    fhat.loo <- function(X,h,a=0,b=1) {
        sapply(1:length(X),function(i){mean(kernel(X[i],X[-i],h,a,b))})
    }
    Fhat <- function(X,h,a=0,b=1) {
        sapply(1:length(X),function(i){mean(kernel.int(X[i],X,h,a,b))})
    }
    ## For search, if gaussian2 use gaussian1
    if(kertype=="gaussian2") {
        ## Gaussian reweighted boundary kernel function (bias of O(h))
        kernel <- function(x,X,h,a=0,b=1) {
            dnorm((x-X)/h)/(h*(pnorm((b-x)/h)-pnorm((a-x)/h)))
        }
    }
    ## Likelihood cross-validation function
    cv.ml.function <- function(h,X,a=0,b=1) {
        f <- fhat.loo(X,h,a,b)
        return(sum(log(ifelse(f > 0 & is.finite(f), f, .Machine$double.xmin))))
    }
    ## Numeric optimization bandwidth search with multistarting
    if(is.null(h) && cv == "numeric") {
        h.vec <- numeric(nmulti)
        cv.vec <- numeric(nmulti)
        rob.spread <- c(sd(X),IQR(X)/1.349)
        rob.spread <- min(rob.spread[rob.spread>0])
        h.init <- runif(nmulti,0.5,1.5)*rob.spread*length(X)**(-0.2)
        ## We are maximizing, so set the pre-optimization maximum for
        ## the objective function to -Inf
        cv.opt <- -Inf
        ## Run the optimizer nmulti times, each time starting from
        ## different random initial values
        for(i in 1:nmulti) {
            foo <- optim(h.init[i],
                         cv.ml.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=Inf,
                         control = list(fnscale = -1),
                         X=X,
                         a=a,
                         b=b)
            ## If we improve upon the previous maximum, save lambda
            ## and reset the maximum (cv.opt)
            h.vec[i] <- foo$par
            cv.vec[i] <- foo$value
            if(foo$value > cv.opt) {
                cv.opt <- foo$value
                h.opt <- foo$par
            }
        }
    }
    ## Grid search and then numeric optimization search (no
    ## multistarting)
    if(is.null(h) && cv == "grid-hybrid") {
        ## First establish a sound starting value using grid search,
        ## then use that starting value for numeric search
        if(is.null(grid)) {
            rob.spread <- c(sd(X),IQR(X)/1.349)
            rob.spread <- min(rob.spread[rob.spread>0])
            constant <- rob.spread*length(X)**(-0.2)
            h.vec <- c(seq(0.25,1.75,length=10),2^(1:25))*constant
            cv.vec <- sapply(1:length(h.vec),function(i){cv.ml.function(h.vec[i],X,a,b)})
            foo <- optim(h.vec[which.max(cv.vec)],
                         cv.ml.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=Inf,
                         control = list(fnscale = -1),
                         X=X,
                         a=a,
                         b=b)
            h.opt <- foo$par
            cv.opt <- foo$value
        } else {
            cv.vec <- sapply(1:length(grid),function(i){cv.ml.function(grid[i],X,a,b)})
            foo <- optim(grid[which.max(cv.vec)],
                         cv.ml.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=Inf,
                         control = list(fnscale = -1),
                         X=X,
                         a=a,
                         b=b)
            h.opt <- foo$par
            cv.opt <- foo$value
        }
    }

    if(kertype=="gaussian2") {
        ## Gaussian reweighted second-order boundary kernel function
        ## (bias of O(h^2))
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
            (aa+bb*z**2)*dnorm(z)/(h*pnorm.zb.m.pnorm.za)
        }
    }

    if(is.null(h.opt)) {
        ## Manual inputted bandwidth
        f <- fhat(X,h,a,b)
        F <- Fhat(X,h,a,b)
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(X,h,a,b)/(h*length(f)))),
                    sd.F=sqrt(abs(F*(1-F)/length(F))),
                    h=h))
    } else {
        f <- fhat(X,h.opt,a,b)
        F <- Fhat(X,h.opt,a,b)
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(X,h.opt,a,b)/(h.opt*length(f)))),
                    sd.F=sqrt(abs(F*(1-F)/length(F))),
                    h=h.opt,
                    nmulti=nmulti,
                    cv.opt=cv.opt))
    }
}
