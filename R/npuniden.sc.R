npuniden.sc <- function(X=NULL,
                        Y=NULL,
                        h=NULL,
                        a=0,
                        b=1,
                        constraint=c("mono.incr",
                            "mono.decr",
                            "concave",
                            "convex",
                            "log-concave",
                            "log-convex")) {

    ## Gaussian kernel function, derivatives up to order two

    kernel <- function(x,X,h,a=0,b=1,deriv=0) {

        if(deriv < 0 || deriv > 2) stop("deriv must be one of 0, 1, or 2")
        
        ## Note: the different sign/position of x in z, z.a, and z.b
        ## changes the sign of the derivative from 1/h to -1/h
        
        z <- (x-X)/h
        z.a <- (a-x)/h
        z.b <- (b-x)/h
        
        K.z <- dnorm(z)
        K.z.d1 <- -z*K.z/h
        K.z.d2 <- (z^2-1)*K.z/h^2    
        
        K.z.a <- dnorm(z.a)
        K.z.b <- dnorm(z.b)
        
        G.z.a <- pnorm(z.a)
        G.z.b <- pnorm(z.b)
        
        G.z.b.m.G.z.a <- G.z.b-G.z.a
        
        G.z.a.d1 <- -K.z.a/h
        G.z.b.d1 <- -K.z.b/h
        
        G.z.b.d1.m.G.z.a.d1 <- G.z.b.d1-G.z.a.d1
        
        G.z.a.d2 <- -z.a*K.z.a/h^2
        G.z.b.d2 <- -z.b*K.z.b/h^2
        
        G.z.b.d2.m.G.z.a.d2 <- G.z.b.d2-G.z.a.d2
        
        ## If the above change due to the use of different base kernel
        ## functions, the formulae below remain unaffected
        
        if(deriv==0) {
            K <- K.z/G.z.b.m.G.z.a
            return(K/h)
        } else if(deriv==1) {
            K <- K.z.d1/G.z.b.m.G.z.a - K.z*G.z.b.d1.m.G.z.a.d1/G.z.b.m.G.z.a^2
            return(K/h)        
        } else if(deriv==2) {
            K <- K.z.d2/G.z.b.m.G.z.a -
                (2*K.z.d1*G.z.b.d1.m.G.z.a.d1+K.z*G.z.b.d2.m.G.z.a.d2)/G.z.b.m.G.z.a^2 +
                    2*K.z*G.z.b.d1.m.G.z.a.d1^2/G.z.b.m.G.z.a^3
            return(K/h)
        }
    }
    
    W.kernel <- function(x,X,h,a=0,b=1,deriv=0) {
        sapply(1:length(x),function(i){kernel(x[i],X,h,a,b,deriv)})
    }
    
    constraint <- match.arg(constraint)
    
    if(is.null(X)) stop("you must pass a vector X")
    if(a>=b) stop("a must be less than b")
    if(any(X<a)) stop("X must be >= a")
    if(any(X>b)) stop("X must be <= b")
    if(!is.null(Y) && any(Y<a)) stop("Y must be >= a")
    if(!is.null(Y) && any(Y>b)) stop("Y must be <= b")
    if(is.null(Y)) Y <- X
    if(is.null(h)) stop("you must provide a bandwidth")
    if(h <= 0) stop("bandwidth h must be positive")

    ## Always require the unrestricted weight matrix

    A <- W.kernel(Y,X,h,a=a,b=b,deriv=0)
    f <- colMeans(A)

    ## First or second order derivative needed

    deriv <- 1
    if(constraint=="concave" || constraint=="convex" || constraint=="log-concave" || constraint=="log-convex") deriv <- 2

    ## Second derivative uses z.a and z.b times other arguments, which
    ## must be finite
    
    if(deriv==2 && a==-Inf) a <- extendrange(X,f=100)[1]
    if(deriv==2 && b==Inf) b <- extendrange(X,f=100)[2]

    A.deriv <- W.kernel(Y,X,h,a=a,b=b,deriv=deriv)
    if(constraint=="log-concave" || constraint=="log-convex") {
        f.deriv <- (colMeans(A.deriv)*f - colMeans(W.kernel(Y,X,h,a=a,b=b,deriv=1))^2)/(f^2)
    } else {
        f.deriv <- colMeans(A.deriv)
    }
    
    sign.deriv <- 1
    if(constraint=="mono.decr" || constraint=="concave" || constraint=="log-concave") sign.deriv <- -1

    n <- length(Y)
    n.train <- length(X)
    
    output.QP <- NULL
    constant <- 1
    attempts <- 0
    while((is.null(output.QP) || any(is.nan(output.QP$solutions))) && attempts < 25) {
        output.QP <- tryCatch(solve.QP(Dmat=diag(n.train)/constant,
                                       dvec=rep(0,n.train),
                                       Amat=cbind(rep(1,n.train),A,sign.deriv*A.deriv),
                                       bvec=c(0,-f,-sign.deriv*f.deriv),
                                       meq=1),
                              error = function(e) NULL)
        constant <- constant*10
        attempts <- attempts+1
    }

    if(is.null(output.QP) || any(is.nan(output.QP$solutions))) stop(" solve.QP was unable to find a solution ")
    if(constraint=="log-concave" || constraint=="log-convex") {
        f.sc <- as.numeric(exp(log(f)+t(A)%*%output.QP$solution))
        ## Constrained derivatives, i.e., derivatives of log(f), not f
        ## - if you wish derivatives of f uncomment the following
        ## line.
        ## f.deriv <- colMeans(A.deriv)
        f.sc.deriv <- f.deriv+t(A.deriv)%*%output.QP$solution
    } else {
        f.sc <- as.numeric(f+t(A)%*%output.QP$solution)
        f.sc.deriv <- as.numeric(f.deriv+t(A.deriv)%*%output.QP$solution)
    }

    corr.factor <- integrate.trapezoidal(Y,f.sc)[length(Y)]/integrate.trapezoidal(Y,f)[length(Y)]
    #f.sc.deriv <- f.sc.deriv/corr.factor    
    f.sc <- f.sc/corr.factor
    
    return(list(f=f,f.sc=f.sc,f.deriv=f.deriv,f.sc.deriv=f.sc.deriv))

}
