## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## ... optional arguments for glpreg/glpcv()

## This function returns a list with the following elements:

## phi: the IV estimator of phi(z) corresponding to the estimated
## deriviative phihat(z)
## phi.prime: the IV derivative estimator
## num.iterations: number of iterations taken by Landweber-Fridman
## norm.stop: the stopping rule for each Landweber-Fridman iteration

## First, a series of functions for local polynomial kernel regression
## Functions for generalized local polynomial regression

## This function returns the weight matrix for a local polynomial

## supports mixed data types. It presumes that Y is in column 1. Basic
## error checking is undertaken. j.reg= strips off weights for mean
## (1), partials up to order p, and cross-partials. All partials and
## cross partials are wrt continuous regressors, and cross-partials
## require k > 1 and p > 1. Shrinking towards the local constant mean,
## first, and second partial derivatives is implemented for regions
## where the local polynomial estimator is ill-conditioned (sparse
## data, small h etc.).

## This function will compute the integral using the trapezoidal rule
## and the cumsum function as we need to compute this in a
## computationally efficient manner.

integrate.trapezoidal <- function(x,y) {
  n <- length(x)
  rank.x <- rank(x)
  order.x <- order(x)
  y <- y[order.x]
  x <- x[order.x]
  int.vec <- numeric(length(x))
  int.vec[1] <- 0
  int.vec[2:n] <- cumsum((x[2:n] - x[2:n-1]) * (y[2:n] + y[2:n-1]) / 2)
  return(int.vec[rank.x])
}

npregivderiv <- function(y,
                         z,
                         w,
                         x=NULL,
                         zeval=NULL,
                         weval=NULL,
                         xeval=NULL,
                         p=1,
                         nmulti=1,
                         random.seed=42,
                         optim.maxattempts = 10,
                         optim.method=c("Nelder-Mead", "BFGS", "CG"),
                         optim.reltol=sqrt(.Machine$double.eps),
                         optim.abstol=.Machine$double.eps,
                         optim.maxit=500,
                         iterate.max=1000,
                         iterate.tol=1.0e-04,
                         constant=0.5,
                         start.phi.zero=FALSE,
                         smooth.while.iterating=TRUE,
                         stop.on.increase=TRUE,
                         smooth.residuals=TRUE,
                         ...) {

  ## First internal to this function we adopt the identical code in
  ## npregiv. This code is internal as the code in crs() will
  ## supersede it (i.e. npglpreg()).

  require(np)

  Kmat.lp <- function(j.reg=0,
                      mydata.train=NULL,
                      mydata.eval=NULL,
                      bws=NULL,
                      p=0,
                      shrink=TRUE,
                      warning.immediate=TRUE) {

    ## Basic error checking...

    if(is.null(mydata.train)) stop("You must provide training data")
    if(is.null(mydata.eval)) mydata.eval <- mydata.train
    if(is.null(bws)) stop("You must provide a bandwidth object")

    n.train=nrow(mydata.train)
    n.eval=nrow(mydata.eval)

    X.train <- as.data.frame(mydata.train)
    X.eval <- as.data.frame(mydata.eval)

    ## Check whether it appears that training and evaluation data are
    ## conformable...

    if(ncol(X.train)!=ncol(X.eval))
      stop("Error: training and evaluation data have unequal number of columns\n")

    X.col.numeric <- sapply(1:ncol(X.train),function(i){is.numeric(X.train[,i])})

    ## k represents the number of numeric regressors, this will return
    ## zero if there are none

    k <- ncol(as.data.frame(X.train[,X.col.numeric]))

    ## Determine the number of cross-product terms. It depends on the
    ## number of numeric regressors and on the order of the
    ## polynomial. The cross products can max out when p > k, p must be
    ## greater than one.

    num.cp <- 0

    if(k > 0) {
      X.train.numeric <- as.data.frame(X.train[,X.col.numeric])
      X.eval.numeric <- as.data.frame(X.eval[,X.col.numeric])
      if(p>=2 && k>=2) for(i in 2:min(k,p)) num.cp <- num.cp + ncol(combn(1:k,i))
    }

    if(j.reg<0||j.reg>(p*k+num.cp))
      stop(paste("Error: j.reg (integer) is invalid\n[min = ", 0, ", max = ",  p*k+num.cp, "]\n",sep=""))

    if(p < 0)
      stop(paste("Error: p (order of polynomial) must be a non-negative integer\np is (", p, ")\n",sep=""))

    Kmat <- matrix(NA,nrow=n.eval,ncol=n.train)

    iota <- rep(1,n.train)

    ## If there are no continuous regressors the local polynomial
    ## estimator collapses to the local constant estimator.

    if(k==0) {
      W <- as.matrix(iota)
    } else {
      ## First column will always be one if we create W this way
      W <- matrix(1,n.train,p*k+1+num.cp)
    }

    for(j in 1:n.eval) {

      if(k > 0 && p > 0) {
        ## Check that we are not conducting local constant
        ## estimation. If we are not (k>0) then create X_{ij}-x_{j},
        ## i=1,..,n.train for the numeric regressors

        X.eval.numeric.diff <- X.eval.numeric[rep(j,n.train),]

        ## If there is only one numeric regressor there are no
        ## crossproducts, but there are always the polynomial terms up
        ## to order p
        l.beg <- 1
        for(l in 1:p) {
          ## Compute the direct polynomial terms. Note that the first
          ## column of W is iota.
          W[,((l-1)*k+2):(l*k+1)] <- as.matrix(X.train.numeric**l)-
            as.matrix(X.eval.numeric.diff**l)
          ## Compute the cross-product terms - this requires that k > 1

          if(k >= 2 && l >= 2 && l <= k) {
            l.end <- l.beg + ncol(as.matrix(combn(1:k,l))) -1

            ## Initialize the cross-product terms if l <= k
            ## ncol(combn(1:k,l))) columns taken which shrinks as l
            ## increases...

            W[,(p*k+1+l.beg):(p*k+1+l.end)] <- as.matrix(X.train.numeric[,combn(1:k,l)[1,]])
            W.tmp <- as.matrix(X.eval.numeric.diff[,combn(1:k,l)[1,]])
            ## Generate cross-products
            for(jj in 2:l) {
              W[,(p*k+1+l.beg):(p*k+1+l.end)] <- W[,(p*k+1+l.beg):(p*k+1+l.end)]*
                as.matrix(X.train.numeric[,combn(1:k,l)[jj,]])
              W.tmp <- W.tmp*as.matrix(X.eval.numeric.diff[,combn(1:k,l)[jj,]])
            }

            W[,(p*k+1+l.beg):(p*k+1+l.end)] <- W[,(p*k+1+l.beg):(p*k+1+l.end)]-W.tmp
            l.beg <- l.end + 1
          }
        }
      }

      if(p == 0) {

        Wmat.sum <- npksum(exdat=X.eval[j,],
                           txdat=X.train,
                           bws=bws,
                           ukertype="liracine",
                           okertype="liracine")$ksum

      } else {

        Wmat.sum <- npksum(exdat=X.eval[j,],
                           txdat=X.train,
                           tydat=W,
                           weights=W,
                           bws=bws,
                           ukertype="liracine",
                           okertype="liracine")$ksum[,,1]

      }

      ## Same in both cases as we need to unroll hence incorporate W
      ## below. We switch roles of tx and ex to get vector instead of
      ## sum. K is a vector of length n.train.

      K <- npksum(txdat=X.eval[j,],
                  exdat=X.train,
                  bws=bws,
                  ukertype="liracine",
                  okertype="liracine")$ksum

      ## p == 0

      Wmat.sum <- as.matrix(Wmat.sum)

      nc <- ncol(Wmat.sum)

      ## No singularity problems...

      if(tryCatch(Wmat.sum.inv <- as.matrix(solve(Wmat.sum)),
                  error = function(e){
                    return(matrix(FALSE,nc,nc))
                  })[1,1]!=FALSE) {

        Kmat[j,] <- sapply(1:n.train,
                            function(i){(Wmat.sum.inv %*% W[i,]*K[i])[(j.reg+1)]})

      } else {

        if(shrink==FALSE) {

          ## If we do not explicitly engage ridging then we do not fail
          ## and terminate, rather, we return NA when Wmat.sum is
          ## singular

          Kmat[j,] <- NA

        } else {

          ## Ridging

          epsilon <- 1/n.train
          ridge <- 0

          while(tryCatch(as.matrix(solve(Wmat.sum+diag(rep(ridge,nc)))),
                         error = function(e){
                           return(matrix(FALSE,nc,nc))
                         })[1,1]==FALSE) {
            ridge <- ridge + epsilon
          }

          Wmat.sum.inv <- as.matrix(solve(Wmat.sum+diag(rep(ridge,nc))))

          ## Add for debugging...

          warning(paste("Ridging obs. ", j, ", ridge = ", signif(ridge,6),sep=""),immediate.=warning.immediate,call.=!warning.immediate)

          ## Now we will create matrices for ridging the mean, first,
          ## and second derivatives (the latter only when p>=2). For
          ## column 1, ridging for the mean. For columns 2:k+1 the first
          ## partial, etc.

          if(p ==1 ) {

            Ridge.mat <- matrix(0,nrow=n.train,ncol=(k+1))

          } else {

            Ridge.mat <- matrix(0,nrow=n.train,ncol=(2*k+1))

          }

          ## Now for the mean - here we explicitly use y*K

          M.vector <- y*K
          sK <- max(.Machine$double.eps,sum(K))

          g.hat.weights <- M.vector/sK

          Ridge.mat[,1] <- ridge*g.hat.weights

          ## First partials... note that this presumes a second order
          ## Gaussian kernel.

          for(m in 1:k) {

            h <- max(.Machine$double.eps,bws[X.col.numeric][m])
            Z <- (as.data.frame(X.train[,X.col.numeric])[,m]-as.data.frame(X.eval[j,X.col.numeric])[,m])/h

            M.1.vector <- M.vector*Z/h
            K.1.vector <- K*Z/h

            sK1divsK <- sum(K.1.vector)/sK

            fp.hat.weights <- M.1.vector/sK-g.hat.weights*sK1divsK

            Ridge.mat[,m+1] <- ridge*fp.hat.weights

            ## Second derivatives only when p >= 2

            if(p >= 2) {

              K.2.vector <- K*(Z^2-1)/h^2
              sK2divsK <- sum(K.2.vector)/sK
              M.2.vector <- M.vector*(Z^2-1)/(h^2)

              sp.hat.weights <- M.2.vector/sK -
                M.1.vector/sK*sK1divsK -
                  fp.hat.weights*sK1divsK -
                    g.hat.weights*(sK2divsK-sK1divsK^2)

              Ridge.mat[,k+m+1] <- ridge*sp.hat.weights

            }

          }

          Kmat[j,] <- sapply(1:n.train,
                              function(i){
                                WK <- W[i,]*K[i]
                                ## Mean
                                WK[1] <- WK[1] + Ridge.mat[i,1]
                                ## First partials
                                WK[2:(k+1)] <- WK[2:(k+1)] + Ridge.mat[i,2:(k+1)]
                                ## Second partials if p >= 2
                                if(p >= 2) WK[(k+2):(2*k+1)] <- WK[(k+2):(2*k+1)] + Ridge.mat[(k+2):(2*k+1)]
                                WKinvWK <- Wmat.sum.inv %*% WK
                                return(WKinvWK[(j.reg+1)])
                              })

        }

      }

    }

    return(Kmat)

  }

  ## No Zero Denominator, used in C code for kernel estimation...

  NZD <- function(a) {
    sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.xmin,a[i]) else max(.Machine$double.xmin,a[i])})
  }

  mypoly <- function(X,degree) {

    if(missing(X)) stop(" X required")
    if(missing(degree)) stop(" degree required")
    if(degree < 1) stop("degree must be a positive integer")

    P <- NULL
    for(i in 1:degree) P <- cbind(P,X**i)

    return(as.matrix(P))

  }

  ## W.glp is a modified version of the polym() function (stats). The
  ## function accepts a vector of degrees and provides a generalized
  ## polynomial with varying polynomial order.

  W.glp <- function(xdat = NULL,
                    degree = NULL) {

    if(is.null(xdat)) stop("Error: You must provide data")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

    xdat <- as.data.frame(xdat)

    xdat.col.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})

    k <- ncol(as.data.frame(xdat[,xdat.col.numeric]))

    if(k > 0) {
      xdat.numeric <- as.data.frame(xdat[,xdat.col.numeric])
    }

    if(length(degree) != ncol(xdat.numeric)) stop(" degree vector and number of numeric predictors incompatible")

    if(all(degree == 0) | k == 0) {

      ## Local constant OR no continuous variables

      return(matrix(1,nrow=nrow(xdat.numeric),ncol=1))

    } else {

      degree.list <- list()
      for(i in 1:k) degree.list[[i]] <- 0:degree[i]
      z <- do.call("expand.grid", degree.list, k)
      s <- rowSums(z)
      ind <- (s > 0) & (s <= max(degree))
      z <- z[ind, ,drop=FALSE]
      if(!all(degree==max(degree))) {
        for(j in 1:length(degree)) {
          d <- degree[j]
          if((d < max(degree)) & (d > 0)) {
            s <- rowSums(z)
            d <- (s > d) & (z[,j,drop=FALSE]==matrix(d,nrow(z),1,byrow=TRUE))
            z <- z[!d, ]
          }
        }
      }
      res <- rep.int(1,nrow(xdat.numeric))
      if(degree[1] > 0) res <- cbind(1, mypoly(xdat.numeric[,1], degree[1]))[, 1 + z[, 1]]
      if(k > 1) for (i in 2:k) if(degree[i] > 0) res <- res * cbind(1, mypoly(xdat.numeric[,i], degree[i]))[, 1 + z[, i]]
      res <- as.matrix(res)
      colnames(res) <- apply(z, 1L, function(x) paste(x, collapse = "."))
      return(as.matrix(cbind(1,res)))

    }

  }

  glpreg <- function(tydat=NULL,
                     txdat=NULL,
                     eydat=NULL,
                     exdat=NULL,
                     bws=NULL,
                     degree=NULL,
                     leave.one.out=FALSE,
                     ...) {

    ## Don't think this error checking is robust

    if(is.null(tydat)) stop("Error: You must provide y data")
    if(is.null(txdat)) stop("Error: You must provide X data")
    if(is.null(bws)) stop("Error: You must provide a bandwidth object")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

    miss.ex = missing(exdat)
    miss.ey = missing(eydat)

    if (miss.ex){
      exdat <- txdat
    }

    txdat <- as.data.frame(txdat)
    exdat <- as.data.frame(exdat)

    maxPenalty <- sqrt(.Machine$double.xmax)

    n.train <- nrow(txdat)
    n.eval <- nrow(exdat)

    ## Check whether it appears that training and evaluation data are
    ## conformable

    if(ncol(txdat)!=ncol(exdat))
      stop("Error: training and evaluation data have unequal number of columns\n")

    if(all(degree == 0)) {

      ## Local constant using only one call to npksum

      if(leave.one.out == TRUE) {

        ## exdat not supported with leave.one.out, but this is only used
        ## for cross-validation hence no exdat

        tww <- npksum(txdat = txdat,
                      weights = as.matrix(data.frame(1,tydat)),
                      tydat = rep(1,length(tydat)),
                      bws = bws,
                      bandwidth.divide = TRUE,
                      leave.one.out = leave.one.out,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

      } else {

        tww <- npksum(txdat = txdat,
                      exdat = exdat,
                      weights = as.matrix(data.frame(1,tydat)),
                      tydat = rep(1,length(tydat)),
                      bws = bws,
                      bandwidth.divide = TRUE,
                      leave.one.out = leave.one.out,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum
      }

      ## Note that as bandwidth approaches zero the local constant
      ## estimator undersmooths and approaches each sample realization,
      ## so use the convention that when the sum of the kernel weights
      ## equals 0, return y. This is unique to this code.

      mhat <- tww[2,]/NZD(tww[1,])

      return(list(mean = mhat))

    } else {

      W <- W.glp(txdat,degree)
      W.eval <- W.glp(exdat,degree)

      ## Local polynomial via smooth coefficient formulation and one
      ## call to npksum

      if(leave.one.out == TRUE) {

        ## exdat not supported with leave.one.out, but this is only used
        ## for cross-validation hence no exdat

        tww <- npksum(txdat = txdat,
                      tydat = as.matrix(cbind(tydat,W)),
                      weights = W,
                      bws = bws,
                      bandwidth.divide = TRUE,
                      leave.one.out = leave.one.out,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

      } else {

        tww <- npksum(txdat = txdat,
                      exdat = exdat,
                      tydat = as.matrix(cbind(tydat,W)),
                      weights = W,
                      bws = bws,
                      bandwidth.divide = TRUE,
                      leave.one.out = leave.one.out,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

      }

      tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[1,,]
      tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n.eval))[-1,,]

      coef.mat <- matrix(maxPenalty,ncol(W),n.eval)
      epsilon <- 1.0/n.eval
      ridge <- double(n.eval)
      doridge <- !logical(n.eval)

      nc <- ncol(tww[,,1])

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                       tyw[,i]+c(ridge.val,rep(0,nc-1))),
                 error = function(e){
                   ridge[i] <<- ridge[i]+epsilon
                   doridge[i] <<- TRUE
                   return(rep(maxPenalty,nc))
                 })
      }

      while(any(doridge)){
        iloo <- (1:n.eval)[doridge]
        coef.mat[,iloo] <- sapply(iloo, ridger)
      }

      mhat <- sapply(1:n.eval, function(i) {
        W.eval[i,, drop = FALSE] %*% coef.mat[,i]
      })

      return(list(mean = mhat,grad = t(coef.mat[-1,])))

    }

  }

  minimand.cv.ls <- function(bws=NULL,
                             ydat=NULL,
                             xdat=NULL,
                             degree=NULL,
                             W=NULL,
                             ...) {

    ## Don't think this error checking is robust

    if(is.null(ydat)) stop("Error: You must provide y data")
    if(is.null(xdat)) stop("Error: You must provide X data")
    if(is.null(W)) stop("Error: You must provide a weighting matrix W")
    if(is.null(bws)) stop("Error: You must provide a bandwidth object")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

    xdat <- as.data.frame(xdat)

    n <- length(ydat)

    maxPenalty <- sqrt(.Machine$double.xmax)

    if(any(bws<=0)) {

      return(maxPenalty)

    } else {

      if(all(degree == 0)) {

        ## Local constant via one call to npksum

        tww <- npksum(txdat = xdat,
                      weights = as.matrix(data.frame(1,ydat)),
                      tydat = rep(1,n),
                      bws = bws,
                      leave.one.out = TRUE,
                      bandwidth.divide = TRUE,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

        mean.loo <- tww[2,]/NZD(tww[1,])

        if (!any(mean.loo == maxPenalty)){
          fv <- mean((ydat-mean.loo)^2)
        } else {
          fv <- maxPenalty
        }

        return(ifelse(is.finite(fv),fv,maxPenalty))

      } else {

        ## Generalized local polynomial via smooth coefficient
        ## formulation and one call to npksum

        tww <- npksum(txdat = xdat,
                      tydat = as.matrix(cbind(ydat,W)),
                      weights = W,
                      bws = bws,
                      leave.one.out = TRUE,
                      bandwidth.divide = TRUE,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

        tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
        tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

        mean.loo <- rep(maxPenalty,n)
        epsilon <- 1.0/n
        ridge <- double(n)
        doridge <- !logical(n)

        nc <- ncol(tww[,,1])

        ## Test for singularity of the generalized local polynomial
        ## estimator, shrink the mean towards the local constant mean.

        ridger <- function(i) {
          doridge[i] <<- FALSE
          ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
          W[i,, drop = FALSE] %*% tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                  tyw[,i]+c(ridge.val,rep(0,nc-1))),
                  error = function(e){
                    ridge[i] <<- ridge[i]+epsilon
                    doridge[i] <<- TRUE
                    return(rep(maxPenalty,nc))
                  })
        }

        while(any(doridge)){
          iloo <- (1:n)[doridge]
          mean.loo[iloo] <- sapply(iloo, ridger)
        }

        if (!any(mean.loo == maxPenalty)){
          fv <- mean((ydat-mean.loo)^2)
        } else {
          fv <- maxPenalty
        }

        return(ifelse(is.finite(fv),fv,maxPenalty))

      }

    }

  }

  minimand.cv.aic <- function(bws=NULL,
                              ydat=NULL,
                              xdat=NULL,
                              degree=NULL,
                              W=NULL,
                              ...) {

    ## Don't think this error checking is robust

    if(is.null(ydat)) stop("Error: You must provide y data")
    if(is.null(xdat)) stop("Error: You must provide X data")
    if(!all(degree==0)) if(is.null(W)) stop("Error: You must provide a weighting matrix W")
    if(is.null(bws)) stop("Error: You must provide a bandwidth object")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))

    xdat <- as.data.frame(xdat)

    n <- length(ydat)

    maxPenalty <- sqrt(.Machine$double.xmax)

    if(any(bws<=0)) {

      return(maxPenalty)

    } else {

      ## This computes the kernel function when i=j (i.e., K(0))

      kernel.i.eq.j <- npksum(txdat = xdat[1,],
                              weights = as.matrix(data.frame(1,ydat)[1,]),
                              tydat = 1,
                              bws = bws,
                              bandwidth.divide = TRUE,
                              ukertype="liracine",
                              okertype="liracine",
                              ...)$ksum[1,1]

      if(all(degree == 0)) {

        ## Local constant via one call to npksum

        tww <- npksum(txdat = xdat,
                      weights = as.matrix(data.frame(1,ydat)),
                      tydat = rep(1,n),
                      bws = bws,
                      bandwidth.divide = TRUE,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

        ghat <- tww[2,]/NZD(tww[1,])

        trH <- kernel.i.eq.j*sum(1/NZD(tww[1,]))

        aic.penalty <- (1+trH/n)/(1-(trH+2)/n)

        if (!any(ghat == maxPenalty) & (aic.penalty > 0)){
          fv <- log(mean((ydat-ghat)^2)) + aic.penalty
        } else {
          fv <- maxPenalty
        }

        return(ifelse(is.finite(fv),fv,maxPenalty))

      } else {

        ## Generalized local polynomial via smooth coefficient
        ## formulation and one call to npksum

        tww <- npksum(txdat = xdat,
                      tydat = as.matrix(cbind(ydat,W)),
                      weights = W,
                      bws = bws,
                      bandwidth.divide = TRUE,
                      ukertype="liracine",
                      okertype="liracine",
                      ...)$ksum

        tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
        tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

        ghat <- rep(maxPenalty,n)
        epsilon <- 1.0/n
        ridge <- double(n)
        doridge <- !logical(n)

        nc <- ncol(tww[,,1])

        ## Test for singularity of the generalized local polynomial
        ## estimator, shrink the mean towards the local constant mean.

        ridger <- function(i) {
          doridge[i] <<- FALSE
          ridge.val <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
          W[i,, drop = FALSE] %*% tryCatch(solve(tww[,,i]+diag(rep(ridge[i],nc)),
                  tyw[,i]+c(ridge.val,rep(0,nc-1))),
                  error = function(e){
                    ridge[i] <<- ridge[i]+epsilon
                    doridge[i] <<- TRUE
                    return(rep(maxPenalty,nc))
                  })
        }

        while(any(doridge)){
          ii <- (1:n)[doridge]
          ghat[ii] <- sapply(ii, ridger)
        }

        trH <- kernel.i.eq.j*sum(sapply(1:n,function(i){
          W[i,, drop = FALSE] %*% solve(tww[,,i]+diag(rep(ridge[i],nc))) %*% t(W[i,, drop = FALSE])
        }))

        if (!any(ghat == maxPenalty)){
          fv <- log(mean((ydat-ghat)^2)) + (1+trH/n)/(1-(trH+2)/n)
        } else {
          fv <- maxPenalty
        }

        return(ifelse(is.finite(fv),fv,maxPenalty))

      }

    }

  }

  glpcv <- function(ydat=NULL,
                    xdat=NULL,
                    degree=NULL,
                    bwmethod=c("cv.ls","cv.aic"),
                    nmulti=1,
                    random.seed=42,
                    optim.maxattempts = 10,
                    optim.method=c("Nelder-Mead", "BFGS", "CG"),
                    optim.reltol=sqrt(.Machine$double.eps),
                    optim.abstol=.Machine$double.eps,
                    optim.maxit=500,
                    debug=FALSE,
                    ...) {

    ## Save seed prior to setting

    if(exists(".Random.seed", .GlobalEnv)) {
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
    } else {
      exists.seed = FALSE
    }

    set.seed(random.seed)

    if(debug) system("rm optim.debug bandwidth.out optim.out")

    ## Don't think this error checking is robust

    if(is.null(ydat)) stop("Error: You must provide y data")
    if(is.null(xdat)) stop("Error: You must provide X data")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))
    if(!is.null(nmulti) && nmulti < 1) stop(paste("Error: nmulti must be a positive integer (minimum 1)\nnmulti is (", nmulti, ")\n",sep=""))

    bwmethod = match.arg(bwmethod)

    optim.method <- match.arg(optim.method)
    optim.control <- list(abstol = optim.abstol,
                          reltol = optim.reltol,
                          maxit = optim.maxit)

    maxPenalty <- sqrt(.Machine$double.xmax)

    xdat <- as.data.frame(xdat)

    num.bw <- ncol(xdat)

    if(is.null(nmulti)) nmulti <- min(5,num.bw)

    ## Which variables are categorical, which are discrete...

    xdat.numeric <- sapply(1:ncol(xdat),function(i){is.numeric(xdat[,i])})

    ## First initialize initial search values of the vector of
    ## bandwidths to lie in [0,1]

    if(debug) write(c("cv",paste(rep("x",num.bw),seq(1:num.bw),sep="")),file="optim.debug",ncolumns=(num.bw+1))

    ## Pass in the local polynomial weight matrix rather than
    ## recomputing with each iteration.

    W <- W.glp(xdat,degree)

    sum.lscv <- function(bw.gamma,...) {

      ## Note - we set the kernel for unordered and ordered regressors
      ## to the liracine kernel (0<=lambda<=1) and test for proper
      ## bounds in sum.lscv.

      if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
        lscv <- minimand.cv.ls(bws=bw.gamma,ydat=ydat,xdat=xdat,...)
      } else {
        lscv <- maxPenalty
      }

      if(debug) write(c(lscv,bw.gamma),file="optim.debug",ncolumns=(num.bw+1),append=TRUE)
      return(lscv)
    }

    sum.aicc <- function(bw.gamma,...) {

      ## Note - we set the kernel for unordered and ordered regressors
      ## to the liracine kernel (0<=lambda<=1) and test for proper
      ## bounds in sum.lscv.

      if(all(bw.gamma>=0)&&all(bw.gamma[!xdat.numeric]<=1)) {
        aicc <- minimand.cv.aic(bws=bw.gamma,ydat=ydat,xdat=xdat,...)
      } else {
        aicc <- maxPenalty
      }

      if(debug) write(c(aicc,bw.gamma),file="optim.debug",ncolumns=(num.bw+1),append=TRUE)
      return(aicc)
    }

    ## Multistarting

    fv.vec <- numeric(nmulti)

    ## Pass in the W matrix rather than recomputing it each time

    for(iMulti in 1:nmulti) {

      num.numeric <- ncol(as.data.frame(xdat[,xdat.numeric]))

      ## First initialize to values for factors (`liracine' kernel)

      init.search.vals <- runif(ncol(xdat),0,1)

      for(i in 1:ncol(xdat)) {
        if(xdat.numeric[i]==TRUE) {
          init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
        }
      }

      ## Initialize `best' values prior to search

      if(iMulti == 1) {
        fv <- maxPenalty
        numimp <- 0
        bw.opt <- init.search.vals
        best <- 1
      }

      if(bwmethod == "cv.ls" ) {

        suppressWarnings(optim.return <- optim(init.search.vals,
                                               fn=sum.lscv,
                                               method=optim.method,
                                               control=optim.control,
                                               degree=degree,
                                               W=W,
                                               ...))

        attempts <- 0
        while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
          init.search.vals <- runif(ncol(xdat),0,1)
          if(xdat.numeric[i]==TRUE) {
            init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
          }
          attempts <- attempts + 1
          optim.control$abstol <- optim.control$abstol * 10.0
          optim.control$reltol <- optim.control$reltol * 10.0
  #        optim.control <- lapply(optim.control, '*', 10.0) ## Perhaps do not want to keep increasing maxit??? Jan 31 2011
          suppressWarnings(optim.return <- optim(init.search.vals,
                                                 fn=sum.lscv,
                                                 method=optim.method,
                                                 control=optim.control,
                                                 degree=degree,
                                                 W=W,
                                                 ...))
        }

      } else {

        suppressWarnings(optim.return <- optim(init.search.vals,
                                               fn=sum.aicc,
                                               method=optim.method,
                                               control=optim.control,
                                               degree=degree,
                                               W=W,
                                               ...))

        attempts <- 0
        while((optim.return$convergence != 0) && (attempts <= optim.maxattempts)) {
          init.search.vals <- runif(ncol(xdat),0,1)
          if(xdat.numeric[i]==TRUE) {
            init.search.vals[i] <- runif(1,.5,1.5)*(IQR(xdat[,i])/1.349)*nrow(xdat)^{-1/(4+num.numeric)}
          }
          attempts <- attempts + 1
          optim.control$abstol <- optim.control$abstol * 10.0
          optim.control$reltol <- optim.control$reltol * 10.0
  #        optim.control <- lapply(optim.control, '*', 10.0) ## Perhaps do not want to keep increasing maxit??? Jan 31 2011
          suppressWarnings(optim.return <- optim(init.search.vals,
                                                 fn = sum.aicc,
                                                 method=optim.method,
                                                 control = optim.control,
                                                 W=W,
                                                 ...))
        }
      }

      if(optim.return$convergence != 0) warning(" optim failed to converge")

      fv.vec[iMulti] <- optim.return$value

      if(optim.return$value < fv) {
        bw.opt <- optim.return$par
        fv <- optim.return$value
        numimp <- numimp + 1
        best <- iMulti
        if(debug) {
          if(iMulti==1) {
            write(cbind(iMulti,t(bw.opt)),"bandwidth.out",ncolumns=(1+length(bw.opt)))
            write(cbind(iMulti,fv),"optim.out",ncolumns=2)
          } else {
            write(cbind(iMulti,t(bw.opt)),"bandwidth.out",ncolumns=(1+length(bw.opt)),append=TRUE)
            write(cbind(iMulti,fv),"optim.out",ncolumns=2,append=TRUE)
          }
        }
      }

    }

    ## Restore seed

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

    return(list(bw=bw.opt,fv=fv,numimp=numimp,best=best,fv.vec=fv.vec))

  }

  ## All code above was for generalized polynomial kernel regression

  console <- newLineConsole()

  ## Basic error checking

  if(!is.logical(start.phi.zero)) stop("start.phi.zero must be logical (TRUE/FALSE)")
  if(!is.logical(smooth.while.iterating)) stop("smooth.while.iterating must be logical (TRUE/FALSE)")
  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")  
  if(!is.logical(smooth.residuals)) stop("smooth.residuals must be logical (TRUE/FALSE)")  

  optim.method <- match.arg(optim.method)

  if(p < 0) stop("The order of the local polynomial must be a positive integer")
  if(nmulti < 1) stop("The number of multistarts must be a positive integer")
  if(optim.maxattempts < 1) stop("The maximum number of optim attempts must be a positive integer")
  if(optim.reltol <= 0) stop("optim.reltol must be positive")
  if(optim.abstol <= 0) stop("optim.abstol must be positive")
  if(optim.maxit <= 0) stop("optim.maxit must be a positive integer")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(iterate.tol <= 0) stop("iterate.tol must be positive")
  if(constant <= 0 || constant >= 1) stop("constant must lie in the range (0,1)")

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(!is.null(x) && NROW(y) != NROW(x)) stop("y and x have differing numbers of rows")

  ## Check for evaluation data

  if(is.null(zeval)) zeval <- z
  if(is.null(weval)) weval <- w
  if(is.null(weval)) xeval <- x

  ## Need to determine how many x, w, z are numeric

  z <- data.frame(z)
  w <- data.frame(w)
  if(!is.null(x)) z <- data.frame(z,x)
  if(!is.null(xeval)) zeval <- data.frame(zeval,xeval)

  z.numeric <- sapply(1:NCOL(z),function(i){is.numeric(z[,i])})
  num.z.numeric <- NCOL(as.data.frame(z[,z.numeric]))

  w.numeric <- sapply(1:NCOL(w),function(i){is.numeric(w[,i])})
  num.w.numeric <- NCOL(as.data.frame(w[,w.numeric]))

  ## Landweber-Fridman

  ## We begin the iteration computing phi.prime.0

  ## Note - here I am only treating the univariate case, so let's
  ## throw a stop with warning for now...

  if(NCOL(z) > 1) stop(" This version supports univariate z only (beta after all)")

  ## For all results we need the density function for Z and the
  ## survivor function for Z (1-CDF of Z)

  console <- printClear(console)
  console <- printPop(console)
  if(is.null(x)) {
    console <- printPush(paste("Computing optimal smoothing for f(z) and S(z) for iteration 1...",sep=""),console)
  } else {
    console <- printPush(paste("Computing optimal smoothing for f(z) and S(z) for iteration 1...",sep=""),console)
  }

  ## Let's compute the bandwidth object for the unconditional
  ## density for the moment. Use the normal-reference rule for speed
  ## considerations.

  bw <- npudensbw(dat=z,bwmethod="normal-reference")
  model.fz <- npudens(tdat=z,bws=bw)
  f.z <- predict(model.fz,newdata=zeval)
  model.Sz <- npudist(tdat=z,bws=bw)
  S.z <- 1-predict(model.Sz,newdata=zeval)

  ## Potential alternative starting rule (consistent with
  ## npregiv). Here we start with E(Y|Z) rather than zero

  if(!start.phi.zero) {

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(y|z) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing for E(y|z,x) for iteration 1...",sep=""),console)
    }

    h <- glpcv(ydat=y,
               xdat=z,
               degree=rep(p, num.z.numeric),
               nmulti=nmulti,
               random.seed=random.seed,
               optim.maxattempts=optim.maxattempts,
               optim.method=optim.method,
               optim.reltol=optim.reltol,
               optim.abstol=optim.abstol,
               optim.maxit=optim.maxit,
               ...)

    if(p == 0) {

      phi.prime <- gradients(npreg(tydat=y,
                                   txdat=z,
                                   exdat=zeval,
                                   bws=h$bw,
                                   gradients=TRUE))[,1]

    } else {

      grad.object <- glpreg(tydat=y,
                            txdat=z,
                            exdat=zeval,
                            bws=h$bw,
                            degree=rep(p, num.z.numeric),
                            ...)$grad

      ## Not sure why this object switches rows and columns, but too
      ## some time to track down (waste!).

      if(p == 1) {
        phi.prime <- grad.object[1,]
      } else {
        phi.prime <- grad.object[,1]
      }

    }

    ## Step 1 - begin iteration - for this we require \varphi_0. To
    ## compute \varphi_{0,i}, we require \mu_{0,i}. For j=0 (first
    ## term in the series), \mu_{0,i} is Y_i.

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1...",sep=""),console)

    ## For stopping rule...
    
    hyw <- glpcv(ydat=y,
                 xdat=w,
                 degree=rep(p, num.w.numeric),
                 nmulti=nmulti,
                 random.seed=random.seed,
                 optim.maxattempts=optim.maxattempts,
                 optim.method=optim.method,
                 optim.reltol=optim.reltol,
                 optim.abstol=optim.abstol,
                 optim.maxit=optim.maxit,
                 ...)
    
    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    exdat=weval,
                    bws=hyw$bw,
                    degree=rep(p, num.w.numeric),
                    ...)$mean

  } else {
    
    ## Step 1 - begin iteration - for this we require \varphi_0. To
    ## compute \varphi_{0,i}, we require \mu_{0,i}. For j=0 (first
    ## term in the series), \mu_{0,i} is Y_i.
    
    mu <- y
    
    ## We also require the mean of \miu_{0,i} shortly...
    
    mean.mu <- mean(mu)
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1...",sep=""),console)

    ## Next, we regress require \mu_{0,i} W (for first iteration mu is y)
    
    hyw <- glpcv(ydat=y,
                 xdat=w,
                 degree=rep(p, num.w.numeric),
                 nmulti=nmulti,
                 random.seed=random.seed,
                 optim.maxattempts=optim.maxattempts,
                 optim.method=optim.method,
                 optim.reltol=optim.reltol,
                 optim.abstol=optim.abstol,
                 optim.maxit=optim.maxit,
                 ...)
    
    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    exdat=weval,
                    bws=hyw$bw,
                    degree=rep(p, num.w.numeric),
                    ...)$mean
    
    ## We require the mean of the fitted values
    
    mean.predicted.E.mu.w <- mean(E.y.w)
    
    ## We need the mean of the fitted values for this (we readily
    ## compute the CDF not the survivor, so anything that is weighted
    ## by the survivor kernel can be expressed as the mean of that
    ## being weighted minus the weighting using the CDF kernel).
    ## Next, we need the weighted sum of the survivor kernel where the
    ## weights are E[\mu_{0,i}|W]. We can write this as the mean of the
    ## \mu_{0,i} minus the weighted sum using the CDF kernel, i.e. if
    ## K is a CDF kernel, then n^{-1}\sum_j \bar K() \mu_{0,i} =
    ## n^{-1}\sum_j (1- K()) \mu_{0,i} = n^{-1}\sum_j\mu_{0,i}-
    ## n^{-1}\sum_j K() \mu_{0,i}
    
    ## Now we compute T^* applied to E.y.w, and this is phi.prime.0 for
    ## j=0.
    
    ## CDF weighted sum (but we need survivor weighted sum...)
    
    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(E.y.w),
                                   operator="integral",
                                   bws=bw$bw)$ksum/length(y)
    
    survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

    phi.prime <- (survivor.weighted.average - S.z*mean.mu)/f.z

  }  
  
  ## Now we can compute phi.0 by integrating phi.prime.0 up to each
  ## sample realization (here we use the trapezoidal rule)

  ## NOTE - this presumes univariate z case... in general this would
  ## be a continuous variable's index

  phi <- integrate.trapezoidal(z[,1],phi.prime)

  ## In the definition of phi we have the integral minus the mean of
  ## the integral with respect to z, so subtract the mean here

  phi <- phi - mean(phi) + mean(y)

  console <- printClear(console)
  console <- printPop(console)
  console <- printPush(paste("Computing optimal smoothing for E(phi|w) (stopping rule) for iteration 1...",sep=""),console)

  ## For the stopping rule, we require E.phi.w

  h.E.phi.w <- glpcv(ydat=phi,
                     xdat=w,
                     degree=rep(p, num.w.numeric),
                     nmulti=nmulti,
                     random.seed=random.seed,
                     optim.maxattempts=optim.maxattempts,
                     optim.method=optim.method,
                     optim.reltol=optim.reltol,
                     optim.abstol=optim.abstol,
                     optim.maxit=optim.maxit,
                     ...)

  E.phi.w <- glpreg(tydat=phi,
                    txdat=w,
                    eydat=phi,
                    exdat=weval,
                    bws=h.E.phi.w$bw,
                    degree=rep(p, num.w.numeric),
                    ...)$mean

  norm.stop <- numeric()

  norm.stop[1] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

#  plot(w[,1],E.y.w,col="red")
#  points(w[,1],E.phi.w,col="red")  

  ## Now we compute mu.0 (a residual of sorts)

  mu <- y - phi

  ## Now we repeat this entire process using mu = y = phi.0 rather than y

  mean.mu <- mean(mu)

  if(smooth.residuals) {
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing for E(mu|w) (stopping rule) for iteration 1...",sep=""),console)

    ## Additional smoothing on top of the stopping rule required, but
    ## we have computed the stopping rule so reuse the bandwidth
    ## vector to be passed below. Here we compute the bandwidth
    ## optimal for the regression of mu on w.
    
    h.E.mu.w <- glpcv(ydat=mu,
                       xdat=w,
                       degree=rep(p, num.w.numeric),
                       nmulti=nmulti,
                       random.seed=random.seed,
                       optim.maxattempts=optim.maxattempts,
                       optim.method=optim.method,
                       optim.reltol=optim.reltol,
                       optim.abstol=optim.abstol,
                       optim.maxit=optim.maxit,
                       ...)
    
    ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w

    predicted.E.mu.w <- glpreg(tydat=mu,
                               txdat=w,
                               eydat=mu,
                               exdat=weval,
                               bws=h.E.mu.w$bw,
                               degree=rep(p, num.w.numeric),
                               ...)$mean

  } else {

    predicted.E.mu.w <- E.y.w - E.phi.w

  }

  ## We again require the mean of the fitted values

  mean.predicted.E.mu.w <- mean(predicted.E.mu.w)

  ## Now we compute T^* applied to mu

  cdf.weighted.average <- npksum(txdat=z,
                                 exdat=zeval,
                                 tydat=as.matrix(predicted.E.mu.w),
                                 operator="integral",
                                 bws=bw$bw)$ksum/length(y)

  survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

  T.star.mu <- (survivor.weighted.average-S.z*mean.mu)/f.z

  ## Now we update phi.prime.0, this provides phi.prime.1, and now
  ## we can iterate until convergence... note we replace phi.prime.0
  ## with phi.prime.1 (i.e. overwrite phi.prime)

  phi.prime <- phi.prime + constant*T.star.mu

  ## This we iterate...

  for(j in 2:iterate.max) {

    ## Save previous in case stop norm increases

    phi.j.m.1 <- phi
    phi.prime.j.m.1 <- phi.prime

    if(smooth.while.iterating) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing optimal smoothing for E(phi|w) for iteration ", j,"...",sep=""),console)
    } else {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing E(phi|w) for iteration ", j,"...",sep=""),console)
    }

    ## NOTE - this presumes univariate z case... in general this would
    ## be a continuous variable's index

    phi <- integrate.trapezoidal(z[,1],phi.prime)

    ## In the definition of phi we have the integral minus the mean of
    ## the integral with respect to z, so subtract the mean here

    phi <- phi - mean(phi) + mean(y)

    ## For the stopping rule, we require E.phi.w

    if(smooth.while.iterating) {

        h.E.phi.w <- glpcv(ydat=phi,
                           xdat=w,
                           degree=rep(p, num.w.numeric),
                           nmulti=nmulti,
                           random.seed=random.seed,
                           optim.maxattempts=optim.maxattempts,
                           optim.method=optim.method,
                           optim.reltol=optim.reltol,
                           optim.abstol=optim.abstol,
                           optim.maxit=optim.maxit,
                           ...)

    }

    E.phi.w <- glpreg(tydat=phi,
                      txdat=w,
                      eydat=phi,
                      exdat=weval,
                      bws=h.E.phi.w$bw,
                      degree=rep(p, num.w.numeric),
                      ...)$mean

    norm.stop[j] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

#    plot(w[,1],E.y.w,col="red")
#    points(w[,1],E.phi.w,col="red")  

    ## Now we compute mu.0 (a residual of sorts)

    mu <- y - phi

    ## Now we repeat this entire process using mu = y = phi.0 rather than y

    mean.mu <- mean(mu)

    if(smooth.residuals) {
    
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing optimal smoothing for E(mu|w) for iteration ", j,"...",sep=""),console)
      
      ## Additional smoothing on top of the stopping rule required, but
      ## we have computed the stopping rule so reuse the bandwidth
      ## vector to be passed below. Here we compute the bandwidth
      ## optimal for the regression of mu on w.
      
      h.E.mu.w <- glpcv(ydat=mu,
                        xdat=w,
                        degree=rep(p, num.w.numeric),
                        nmulti=nmulti,
                        random.seed=random.seed,
                        optim.maxattempts=optim.maxattempts,
                        optim.method=optim.method,
                        optim.reltol=optim.reltol,
                        optim.abstol=optim.abstol,
                        optim.maxit=optim.maxit,
                        ...)
      
      ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w
      
      predicted.E.mu.w <- glpreg(tydat=mu,
                                 txdat=w,
                                 eydat=mu,
                                 exdat=weval,
                                 bws=h.E.mu.w$bw,
                                 degree=rep(p, num.w.numeric),
                                 ...)$mean
      
    } else {
      
      predicted.E.mu.w <- E.y.w - E.phi.w
      
    }
    
    mean.predicted.E.mu.w <- mean(predicted.E.mu.w)

    ## Now we compute T^* applied to mu

    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(predicted.E.mu.w),
                                   operator="integral",
                                   bws=bw$bw)$ksum/length(y)

    survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

    T.star.mu <- (survivor.weighted.average-S.z*mean.mu)/f.z

    ## Now we update, this provides phi.prime.1, and now we can
    ## iterate until convergence...

    phi.prime <- phi.prime + constant*T.star.mu

    ## If stopping rule criterion increases or we are below stopping
    ## tolerance then break

    if(norm.stop[j] < iterate.tol) break()
    if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
      phi <- phi.j.m.1 
      phi.prime <- phi.prime.j.m.1
      break()
    }

  }

  console <- printClear(console)
  console <- printPop(console)

  if(j == iterate.max) warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")

  return(list(phi=phi,phi.prime=phi.prime,num.iterations=j,norm.stop=norm.stop))

}

