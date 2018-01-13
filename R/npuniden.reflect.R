npuniden.reflect <- function(X = NULL,
                             Y = NULL,
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
    if(!is.null(Y) && any(Y<a)) stop("Y must be >= a")
    if(!is.null(Y) && any(Y>b)) stop("Y must be <= b")
    if(is.null(Y)) Y <- X
    if(!is.null(h) && h <= 0) stop("bandwidth h must be positive")
    if(is.null(h)) {
        bw <- npudensbw(~X,...)
        hh <- bw$bw
    } else {
        hh <- h
    }
    if(is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*a,-X+2*b)
        f.reflect <- npudens(tdat=X.reflect,edat=Y,bws=hh,...)
        f <- 3*fitted(f.reflect)
        cdf <- 3*fitted(npudist(tdat=X.reflect,edat=Y,bws=hh,...))-1
        std <- 3*se(f.reflect)
    } else if(is.finite(a) && !is.finite(b)) {
        X.reflect <- c(X,-X+2*a)
        f.reflect <- npudens(tdat=X.reflect,edat=Y,bws=hh,...)
        f <- 2*fitted(f.reflect)
        cdf <- 2*fitted(npudist(tdat=X.reflect,edat=Y,bws=hh,...))-1      
        std <- 2*se(f.reflect)
    } else if(!is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*b)
        f.reflect <- npudens(tdat=X.reflect,edat=Y,bws=hh,...)
        f <- 2*fitted(f.reflect)
        cdf <- 2*fitted(npudist(tdat=X.reflect,edat=Y,bws=hh,...))      
        std <- 2*se(f.reflect)
    }
    if(is.null(h)) {
        return(list(f=f,
                    F=cdf,
                    sd.f=std,
                    sd.F=sqrt(abs(cdf*(1-cdf)/length(cdf))),
                    h=hh,
                    nmulti=length(bw$fval.history),
                    cv.opt=bw$fval))
    } else {
        return(list(f=f,
                    F=cdf,
                    sd.f=std,
                    sd.F=sqrt(abs(cdf*(1-cdf)/length(cdf))),
                    h=hh))
    }
}
