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
                              nmulti=1,
                              proper=FALSE) {
    kertype <- match.arg(kertype)
    cv <- match.arg(cv)
    bwmethod <- match.arg(bwmethod)
    if(!is.null(grid) && anyNA(grid)) stop("grid must not contain missing values")
    if(!is.null(grid) && any(grid<=0, na.rm = TRUE)) stop(" the grid vector must contain positive values")
    if(is.null(X)) stop("you must pass a vector X")
    if(anyNA(X)) stop("X must not contain missing values")
    if(kertype=="gamma" || kertype=="rigaussian") b <- Inf
    if(kertype=="fbl") b <- Inf
    if(kertype=="fbu") a <- -Inf
    if(a>=b) stop("a must be less than b")
    if(any(X<a, na.rm = TRUE)) stop("X must be >= a")
    if(any(X>b, na.rm = TRUE)) stop("X must be <= b")
    if(!is.null(Y) && anyNA(Y)) stop("Y must not contain missing values")
    if(!is.null(Y) && any(Y<a, na.rm = TRUE)) stop("Y must be >= a")
    if(!is.null(Y) && any(Y>b, na.rm = TRUE)) stop("Y must be <= b")
    if(is.null(Y)) Y <- X
    if(!is.null(h) && h <= 0) stop("bandwidth h must be positive")
    if(nmulti < 1) stop("number of multistarts nmulti must be positive")
    if(!is.logical(proper)) stop("proper must be either TRUE or FALSE")
    if(kertype=="gaussian2" && (!is.finite(a) || !is.finite(b))) stop("finite bounds are required for kertype gaussian2")
    beta.kernel <- kertype %in% c("beta1", "beta2")
    if(beta.kernel && (!is.finite(a) || !is.finite(b)))
        stop("finite bounds are required for beta kernels")
    if(kertype=="beta2") {
        if(!is.null(h) && (!is.finite(h) || h > 1/4))
            stop("beta2 bandwidth h must satisfy 0 < h <= 1/4 on normalized support")
        if(!is.null(grid) && any(!is.finite(grid) | grid > 1/4))
            stop("beta2 grid values must satisfy 0 < h <= 1/4 on normalized support")
    }
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
            if(x < 2*h) {
                dbeta(X,rho(x,h),(1-x)/h)/(b-a)
            } else if(2*h <= x && x <= 1-2*h) {
                dbeta(X,x/h,(1-x)/h)/(b-a)
            } else if(x > 1-2*h) {
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
                out <- numeric(length(t))
                mask <- (c <= t) & (t <= 2 + c)
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- .75 * (c + 1 - 1.25 * (1 + 2 * c) * (tm - c)^2) * (tm - (c + 2))^2
                }
                out / h
            } else if((a+h <= x && x <= b-h) || h >= (b-a)) {
                z.a <- (a-x)/h
                z.b <- (b-x)/h  
                rw <- (3*(z.b^5-z.a^5)-10*(z.b^3-z.a^3)+15*(z.b-z.a))/16
                rw[rw>1] <- 1
                out <- numeric(length(t))
                mask <- abs(t) < 1
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- (15 / 16) * (1 - tm^2)^2 / (h * rw)
                }
                out
            } else if(x > b-h && h < (b-a)) {
                c <- (b-x)/h
                out <- numeric(length(t))
                mask <- (c - 2 <= t) & (t <= c)
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- .75 * (1 - c + 1.25 * (-1 + 2 * c) * (tm - c)^2) * (tm - (c - 2))^2
                }
                out / h
            }
        }
    } else if(kertype=="fbl") {
        ## Floating boundary kernel (Scott (1992), Page 46), left bound
        kernel <- function(x,X,h,a=0,b=1) {
            t <- (X-x)/h
            if(x < a+h) {
                c <- (a-x)/h
                out <- numeric(length(t))
                mask <- (c <= t) & (t <= 2 + c)
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- .75 * (c + 1 - 1.25 * (1 + 2 * c) * (tm - c)^2) * (tm - (c + 2))^2
                }
                out / h
            } else {
                out <- numeric(length(t))
                mask <- abs(t) < 1
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- (15 / 16) * (1 - tm^2)^2 / h
                }
                out
            }
        }
    } else if(kertype=="fbu") {
        kernel <- function(x,X,h,a=0,b=1) {
            ## Floating boundary kernel (Scott (1992), Page 46), right bound
            t <- (X-x)/h
            if(x <= b-h) {
                out <- numeric(length(t))
                mask <- abs(t) < 1
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- (15 / 16) * (1 - tm^2)^2 / h
                }
                out
            } else {
                c <- (b-x)/h
                out <- numeric(length(t))
                mask <- (c - 2 <= t) & (t <= c)
                if(any(mask)) {
                    tm <- t[mask]
                    out[mask] <- .75 * (1 - c + 1.25 * (-1 + 2 * c) * (tm - c)^2) * (tm - (c - 2))^2
                }
                out / h
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
        geometry <- .np_quadrature_prepare(X.seq)
        sapply(seq_along(X), function(i){.np_quadrature_total(X.seq,
            h*kernel(X[i],X.seq,h,a,b)**2, geometry)})
    }
    fhat <- function(X,Y,h,a=0,b=1) {
        sapply(seq_along(Y), function(i){mean(kernel(Y[i],X,h,a,b))})
    }
    final.fit <- function(h) {
        density <- function(y) fhat(X,y,h,a,b)
        f <- density(Y)
        # Only the ordinary, unbounded Gaussian kernel has this analytic CDF.
        if(kertype=="gaussian1" && !is.finite(a) && !is.finite(b)) {
            F <- vapply(Y,function(y) mean(pnorm((y-X)/h)),0.0)
            return(list(f=f,F=F))
        }
        branch.points <- switch(kertype,
            beta2=c(a+2*h*(b-a),b-2*h*(b-a)),
            fb=c(a+h,b-h),fbl=a+h,fbu=b-h,numeric())
        integration.scale <- if(beta.kernel) h*(b-a) else h
        integration.widths <- rep(integration.scale,length(X))
        if(beta.kernel) {
            normalized.X <- (X-a)/(b-a)
            integration.widths <- (b-a)*sqrt(h*(normalized.X*(1-normalized.X)+h))
        } else if(kertype=="gamma") {
            integration.widths <- sqrt(h*(X-a)+h*h)
        } else if(kertype=="rigaussian") {
            integration.widths[X>a] <- sqrt(h*(X[X>a]-a))
            integration.scale <- min(integration.widths)
        }
        integrand <- if(proper) function(y) pmax(density(y),0) else density
        integral <- .np_density_integral(integrand,Y,a,b,
            X,integration.scale,branch.points,positive=proper,widths=integration.widths)
        F <- integral$F
        if(proper) {
            if(!is.finite(integral$total) || integral$total<=0)
                stop("proper density requires finite positive whole-support mass")
            f <- pmax(f,0)/integral$total
            F <- F/integral$total
        }
        list(f=f,F=F)
    }
    fhat.loo <- function(X,h,a=0,b=1) {
        n <- length(X)
        if (n <= 1L) return(rep(NA_real_, n))
        sapply(seq_along(X), function(i){
            kv <- kernel(X[i], X, h, a, b)
            (sum(kv) - kv[i])/(n - 1L)
        })
    }
    if(bwmethod=="cv.ml") {
        ## Likelihood cross-validation function (maximizing)
        fnscale <- list(fnscale = -1)
        cv.function <- function(h,X,a=0,b=1) {
            f.loo <- fhat.loo(X,h,a,b)
            good <- (f.loo > 0) & is.finite(f.loo)
            f.safe <- f.loo
            f.safe[!good] <- .Machine$double.xmin
            return(sum(log(f.safe)))
        }
    } else {
        ## Least-squares cross-validation function (minimizing)
        fnscale <- list(fnscale = 1) 
        cv.geometry <- if(is.null(h)) .np_quadrature_prepare(X) else NULL
        cv.function <- function(h,X,a=0,b=1) {
            cv.ls <- .np_quadrature_total(X,fhat(X,X,h,a,b)**2,cv.geometry)-2*mean(fhat.loo(X,h,a,b))
            if (is.finite(cv.ls)) cv.ls else sqrt(sqrt(.Machine$double.xmax))
        }
    }
    cv.cache <- .np_objective_exact_cache_new(npObjectiveCacheEnabled())
    cv.function.uncached <- cv.function
    cv.function <- function(h,X,a=0,b=1) {
        cache.hit <- .np_objective_exact_cache_get(cv.cache, h)
        if (isTRUE(cache.hit$hit))
            return(cache.hit$value)
        value <- cv.function.uncached(h,X,a,b)
        .np_objective_exact_cache_put(cv.cache, cache.hit$token, value)
        value
    }
    ## Grid search and then numeric optimization search (no
    ## multistarting, but sound starting point always used for
    ## subsequent refinement by optim)
    if(is.null(h) && cv == "grid-hybrid") {
        ## First establish a sound starting value using grid search,
        ## then use that starting value for numeric search
        if(is.null(grid)) {
            # Beta h is dimensionless, as are the kernel's shape parameters.
            search.X <- if(beta.kernel) (X-a)/(b-a) else X
            rob.spread <- c(sd(search.X),IQR(search.X)/1.349)
            rob.spread <- min(rob.spread[rob.spread>0])
            constant <- rob.spread*length(X)**(-0.2)
            h.vec <- c(seq(0.25,1.75,length=10),2^(1:25))*constant
            if(kertype=="beta2")
                h.vec <- unique(c(h.vec[is.finite(h.vec) &
                    h.vec >= sqrt(.Machine$double.eps) & h.vec <= 1/4], 1/4))
            cv.vec <- sapply(seq_along(h.vec), function(i){cv.function(h.vec[i],X,a,b)})
            start.idx <- if (bwmethod=="cv.ml") which.max(cv.vec) else which.min(cv.vec)
            upper.bound <- if (kertype=="beta2") 1/4 else Inf
            foo <- optim(h.vec[start.idx],
                         cv.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=upper.bound,
                         control = fnscale,
                         X=X,
                         a=a,
                         b=b)
            h.opt <- foo$par
            cv.opt <- foo$value
        } else {
            cv.vec <- sapply(seq_along(grid), function(i){cv.function(grid[i],X,a,b)})
            start.idx <- if (bwmethod=="cv.ml") which.max(cv.vec) else which.min(cv.vec)
            upper.bound <- if (kertype=="beta2") 1/4 else Inf
            foo <- optim(grid[start.idx],
                         cv.function,
                         method="L-BFGS-B",
                         lower=sqrt(.Machine$double.eps),
                         upper=upper.bound,
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
        fit <- final.fit(h)
        f <- fit$f
        F <- fit$F
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(Y,h,a,b)/(h*length(X)))),
                    sd.F=sqrt(abs(F*(1-F)/length(X))),
                    h=h))
    } else {
        ## Search bandwidth
        fit <- final.fit(h.opt)
        f <- fit$f
        F <- fit$F
        return(list(f=f,
                    F=F,
                    sd.f=sqrt(abs(f*int.kernel.squared(Y,h.opt,a,b)/(h.opt*length(X)))),
                    sd.F=sqrt(abs(F*(1-F)/length(X))),
                    h=h.opt,
                    nmulti=nmulti,
                    cv.opt=cv.opt))
    }
}
