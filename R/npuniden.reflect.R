npuniden.reflect <- function(X = NULL,
                             h = NULL,
                             a = 0,
                             b = 1,
                             ...) {
    ## Data reflection - note that the bandwidth should be based on
    ## the original sample of n observations (SCOTT (1992), Page 149)
    if(is.null(X)) stop("you must pass a vector X")
    if(!is.finite(a) && !is.finite(b)) stop("one of either a or b must be finite")
    if(a>=b) stop("a must be less than b")
    if(any(X<a)) stop("X must be >= a")
    if(any(X>b)) stop("X must be <= b")
    if(!is.null(h) && h <= 0) stop("bandwidth h must be positive")
    if(is.null(h)) {
        bw <- npudensbw(~X,...)
        hh <- bw$bw
    } else {
        hh <- h
    }
    if(is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*a,-X+2*b)
        f.reflect <- npudens(tdat=X.reflect,edat=X,bws=hh,...)
        f <- 3*fitted(f.reflect)
        std <- 3*se(f.reflect)
    } else if(is.finite(a) && !is.finite(b)) {
        X.reflect <- c(X,-X+2*a)
        f.reflect <- npudens(tdat=X.reflect,edat=X,bws=hh,...)
        f <- 2*fitted(f.reflect)
        std <- 2*se(f.reflect)
    } else if(!is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*b)
        f.reflect <- npudens(tdat=X.reflect,edat=X,bws=hh,...)
        f <- 2*fitted(f.reflect)
        std <- 2*se(f.reflect)
    }
    cdf <- integrate.trapezoidal(X,f)
    if(is.null(h)) {
        return(list(f=f,F=cdf,sd.f=std,sd.F=sqrt(cdf*(1-cdf)/length(cdf)),h=hh,nmulti=length(bw$fval.history),cv.opt=bw$fval))
    } else {
        return(list(f=f,F=cdf,sd.f=std,sd.F=sqrt(cdf*(1-cdf)/length(cdf)),h=hh))
    }
}