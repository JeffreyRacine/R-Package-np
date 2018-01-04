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
    integrate.trapezoidal <- function(x,y) {
        n <- length(x)
        rank.x <- rank(x)
        order.x <- order(x)
        y <- y[order.x]
        x <- x[order.x]
        int.vec <- numeric(length(x))
        int.vec[1] <- 0
        int.vec[2:n] <- cumsum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2)
        return((int.vec[rank.x])[n])
    }
    if(is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*a,-X+2*b)
        f.reflect <- npudens(~X.reflect,bws=hh,...)
        f <- 3*predict(f.reflect,newdata=data.frame(X.reflect=X))
    } else if(is.finite(a) && !is.finite(b)) {
        X.reflect <- c(X,-X+2*a)
        f.reflect <- npudens(~X.reflect,bws=hh,...)
        f <- 2*predict(f.reflect,newdata=data.frame(X.reflect=X))
    } else if(!is.finite(a) && is.finite(b)) {
        X.reflect <- c(X,-X+2*b)
        f.reflect <- npudens(~X.reflect,bws=hh,...)
        f <- 2*predict(f.reflect,newdata=data.frame(X.reflect=X))
    }
    if(is.null(h)) {
        return(list(f=f,h=hh,nmulti=length(bw$fval.history),cv.opt=bw$fval,int=integrate.trapezoidal(X,f)))
    } else {
        return(list(f=f,h=hh,int=integrate.trapezoidal(X,f)))
    }
}
