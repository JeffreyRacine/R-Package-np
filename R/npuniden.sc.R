npuniden.sc <- function(X=NULL,
                        Y=NULL,
                        h=NULL,
                        a=0,
                        b=1,
                        lb=NULL,
                        ub=NULL,
                        extend.range=0,
                        num.grid=0,
                        function.distance=TRUE,
                        integral.equal=FALSE,
                        constraint=c("density",
                            "mono.incr",
                            "mono.decr",
                            "concave",
                            "convex",
                            "log-concave",
                            "log-convex")) {

    constraint <- match.arg(constraint)

    ## Gaussian1 kernel function, derivatives up to order two

    kernel <- function(x,X,h,a=0,b=1,deriv=0) {

        if(deriv < 0 || deriv > 2) stop("deriv must be one of 0, 1, or 2")
        
        ## Note: the different sign/position of x in z, z.a, and z.b
        ## changes the sign of the derivative from 1/h to -1/h
        
        z <- (x-X)/h
        z.a <- (a-x)/h
        z.b <- (b-x)/h
        
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
        
        K.z <- dnorm(z)
        K.z.d1 <- -z*K.z/h
        K.z.d2 <- (z^2-1)*K.z/h^2

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
    
    int.kernel.squared <- function(X,h,a=a,b=b) {
        ## Use numeric integration to compute Kappa, the integral of
        ## the square of the kernel function needed for the asymptotic
        ## standard error of the density estimate seq(a,b) will barf
        ## on -Inf or Inf, trap these cases and use extendrange
        if(is.finite(a) && is.finite(b)) X.seq <- seq(a,b,length=1000)
        if(is.finite(a) && !is.finite(b)) X.seq <- seq(a,extendrange(X,f=10)[2],length=1000)
        if(!is.finite(a) && is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],b,length=1000)
        if(!is.finite(a) && !is.finite(b)) X.seq <- seq(extendrange(X,f=10)[1],extendrange(X,f=10)[2],length=1000)
        sapply(1:length(X),function(i){integrate.trapezoidal(X.seq,h*kernel(X[i],X.seq,h,a,b,deriv=0)^2)[length(X.seq)]})
    }

    W.kernel <- function(x,X,h,a=0,b=1,deriv=0) {
        sapply(1:length(x),function(i){kernel(x[i],X,h,a,b,deriv)})
    }
    
    if(is.null(X)) stop("you must pass a vector X")
    if(a>=b) stop("a must be less than b")
    if(any(X<a)) stop("X must be >= a")
    if(any(X>b)) stop("X must be <= b")
    if(!is.null(Y) && any(Y<a)) stop("Y must be >= a")
    if(!is.null(Y) && any(Y>b)) stop("Y must be <= b")
    if(is.null(h)) stop("you must provide a bandwidth")
    if(h <= 0) stop("bandwidth h must be positive")
    if(num.grid < 0) stop("num.grid must be a non-negative integer")
    if(constraint=="density" && is.null(lb) && is.null(ub)) stop("you must provide lower and/or upper bounds when constraining the density")
    if(!is.null(lb) && any(lb<0)) stop("lower bound must be non-negative")
    if(!is.null(ub) && any(ub<0)) stop("upper bound must be non-negative")
    if(!is.null(lb) && !is.null(ub) && any(ub<lb)) stop("upper bound must be greater than or equal to lower bound")

    ## First elements are X if Y=NULL or Y, rest are to make sure
    ## constraints are imposed on the bulk of the support and a bit
    ## outside governed by f=extend.range

    if(num.grid==0) {
        grid <- NULL
    } else {
        grid <- seq(max(c(a,extendrange(X,f=extend.range)[1])),min(c(b,extendrange(X,f=extend.range)[2])),,num.grid)
    }

    X.grid <- c(Y,X,grid)
    
    ## We always need the unconstrained weight matrix and estimate

    A <- W.kernel(X.grid,X,h,a=a,b=b,deriv=0)
    f <- colMeans(A)

    ## First or second order derivative needed, place them both in *.deriv

    deriv <- 0
    if(constraint=="mono.decr" || constraint=="mono.incr") deriv <- 1
    if(constraint=="concave" || constraint=="convex" || constraint=="log-concave" || constraint=="log-convex") deriv <- 2

    ## Second derivative uses z.a and z.b times other arguments, which
    ## must be finite, so set to very small and large values depending
    ## on the range of the training data
    
    if(deriv==2 && a==-Inf) a <- extendrange(X,f=100)[1]
    if(deriv==2 && b==Inf) b <- extendrange(X,f=100)[2]

    ## The derivative for the log-density estimate differs from that
    ## for the raw density, and also requires the first derivative

    A.deriv <- W.kernel(X.grid,X,h,a=a,b=b,deriv=deriv)
    if(constraint=="log-concave" || constraint=="log-convex") {
        f.deriv <- (colMeans(A.deriv)*f - colMeans(W.kernel(X.grid,X,h,a=a,b=b,deriv=1))^2)/(f^2)
    } else {
        f.deriv <- colMeans(A.deriv)
    }

    ## The sign of the constraint matrix Amat and bvec switch
    ## depending on whether the constraint is >= or <=, so automate
    ## this
    
    sign.deriv <- 1
    if(constraint=="mono.decr" || constraint=="concave" || constraint=="log-concave") sign.deriv <- -1

    ## Length of the grid vector and training data vector
    
    n.grid <- length(X.grid)
    n.train <- length(X)

    ## Solve the quadratic program

    solve.QP.flag <- TRUE
    output.QP <- NULL
    constant <- c(1,1/10,10,1/100,100,1/1000,1000,1/10000,10000,1/100000,100000)
    attempts.max <- length(constant)
    attempts <- 1
    while((is.null(output.QP) || any(is.na(output.QP$solution))) && attempts <= attempts.max) {
        if(attempts==attempts.max && !function.distance) {
            warning("solve.QP was unable to find a solution with function.distance=FALSE, restarting with function.distance=TRUE", immediate. = TRUE)
            attempts <- 1
            function.distance <- TRUE
        }
        if(function.distance) {
            ## Non-identity forcing matrix minimizes the squared
            ## function difference distance
            Dmat <- (A%*%t(A)+sqrt(.Machine$double.eps)*diag(n.train))/constant[attempts]
        } else {
            ## Identity forcing matrix minimizes the squared weight
            ## distance
            Dmat <- diag(n.train)/constant[attempts]
        }
        ## The unconstrained weight vector contains zeros
        dvec <- rep(0,n.train)
        if(constraint=="log-concave" || constraint=="log-convex") {
            ## The log-concave/convex estimate will be non-negative,
            ## no need to impose this constraint (will be slack)
            Amat <- cbind(rep(1,n.train),sign.deriv*A.deriv)
            bvec <- c(0,-sign.deriv*f.deriv)
        } else if(constraint!="density") {
            ## For other constraints we need to enforce non-negativity
            ## of the estimate - in either case, the unconstrained
            ## weights are zero hence the "sum to zero" constraint
            Amat <- cbind(rep(1,n.train),A,sign.deriv*A.deriv)
            bvec <- c(0,-f,-sign.deriv*f.deriv)
        } else {
            ## Constrain the density
            Amat <- cbind(rep(1,n.train),A,-A)
            bvec <- c(0,lb-f,f-ub)
        }
        output.QP <- tryCatch(solve.QP(Dmat=Dmat,
                                       dvec=dvec,
                                       Amat=Amat,
                                       bvec=bvec,
                                       meq=1),
                              error = function(e) NULL)
        ## Sometimes rescaling Dmat can overcome deficiencies in
        ## solve.QP
        attempts <- attempts+1
    }

    ## If solve.QP cannot find a solution issue an immediate warning
    ## but return the unconstrained vector

    if(is.null(output.QP) || any(is.na(output.QP$solution))) {
        warning("solve.QP was unable to find a solution, unconstrained estimate returned", immediate. = TRUE)
        output.QP$solution <- rep(0,n.train)
        solve.QP.flag <- FALSE
    }
    
    if(constraint=="log-concave" || constraint=="log-convex") {
        f.sc <- as.numeric(exp(log(f)+t(A)%*%output.QP$solution))
        ## Constrained derivatives, i.e., derivatives of log(f), not f
        ## - if you wish derivatives of f uncomment the following
        ## line
        ## f.deriv <- colMeans(A.deriv)
        f.sc.deriv <- f.deriv+t(A.deriv)%*%output.QP$solution
    } else {
        f.sc <- as.numeric(f+t(A)%*%output.QP$solution)
        f.sc.deriv <- as.numeric(f.deriv+t(A.deriv)%*%output.QP$solution)
    }

    ## If the option integral.equal=TRUE is set, adjust the
    ## constrained estimate to have the same integral as the
    ## unconstrained estimate

    corr.factor <- 1
    int.f.sc <- integrate.trapezoidal(X.grid,f.sc)[n.grid]
    int.f <- integrate.trapezoidal(X.grid,f)[n.grid]
    if(integral.equal) corr.factor <- int.f.sc/int.f
    f.sc <- f.sc/corr.factor

    ## The constraints are imposed on X.grid, but we only want the
    ## constrained estimate at the sample points X (or evaluation
    ## points Y if Y is specified)

    if(is.null(Y)) {
        index <- 1:n.train
    } else {
        index <- 1:length(Y)
    }

    se.f <- sqrt(abs(f*int.kernel.squared(X.grid,h,a,b)/(h*length(f))))[index]
    se.f.sc <- sqrt(abs(f.sc*int.kernel.squared(X.grid,h,a,b)/(h*length(f.sc))))[index]    

    F <- integrate.trapezoidal(X.grid[index],f[index])
    F.sc <- integrate.trapezoidal(X.grid[index],f.sc[index])
    f <- f[index]
    f.sc <- f.sc[index]
    f.deriv <- f.deriv[index]
    f.sc.deriv <- f.sc.deriv[index]

    ## Return a list
    
    return(list(f=f,
                f.sc=f.sc,
                se.f=se.f,
                se.f.sc=se.f.sc,
                f.deriv=f.deriv,
                f.sc.deriv=f.sc.deriv,
                F=F,
                F.sc=F.sc,
                se.F=sqrt(abs(F*(1-F)/length(F))),
                se.F.sc=sqrt(abs(F.sc*(1-F.sc)/length(F.sc))),                
                f.integral=int.f,
                f.sc.integral=int.f.sc,
                solve.QP=solve.QP.flag,
                attempts=attempts))

}
