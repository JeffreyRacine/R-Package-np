## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## xeval: optional evaluation data for the exogenous predictors

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## p: order of the local polynomial kernel estimator (p=0 is local
##    constant, p=1 local linear etc.)

## This function returns a list with at least the following elements:

## phi: the IV estimator of phi(z)
## convergence: a character string indicating whether/why iteration terminated

npregiv <- function(y,
                    z,
                    w,
                    x=NULL,
                    zeval=NULL,
                    xeval=NULL,
                    alpha=NULL,
                    alpha.iter=NULL,
                    alpha.max=1.0e-01,
                    alpha.min=1.0e-10,
                    alpha.tol=.Machine$double.eps^0.25,
                    bw=NULL,
                    constant=0.5,
                    iterate.diff.tol=1.0e-08,
                    iterate.max=1000,
                    iterate.Tikhonov=TRUE,
                    iterate.Tikhonov.num=1,
                    method=c("Landweber-Fridman","Tikhonov"),
                    nmulti=NULL,
                    optim.abstol=.Machine$double.eps,
                    optim.maxattempts = 10,
                    optim.maxit=500,
                    optim.method=c("Nelder-Mead", "BFGS", "CG"),
                    optim.reltol=sqrt(.Machine$double.eps),
                    p=1,
                    penalize.iteration=TRUE,
                    random.seed=42,
                    return.weights.phi=FALSE,
                    return.weights.phi.deriv.1=FALSE,
                    return.weights.phi.deriv.2=FALSE,
                    smooth.residuals=TRUE,
                    start.from=c("Eyz","EEywz"),
                    starting.values=NULL,
                    stop.on.increase=TRUE,
                    ...) {

  ptm.start <- proc.time()
  cl <- match.call()

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr> to reproduce illustrations in
  ## the following papers:

  ## A) Econometrica (2011), Volume 79, pp. 1541-1565
  ## "Nonparametric Instrumental Regression"
  ## S. Darolles, Y. Fan, J.P. Florens, E. Renault

  ## B) Econometrics Journal (2010), volume 13, pp. S1-S27. doi:
  ## 10.1111/j.1368-423X.2010.00314.x

  ## "The practice of non-parametric estimation by solving inverse
  ## problems: the example of transformation models"

  ## FREDERIQUE FEVE AND JEAN-PIERRE FLORENS
  ## IDEI and Toulouse School of Economics, Universite de Toulouse
  ## Capitole 21 alle de de Brienne, 31000 Toulouse, France. E-mails:
  ## feve@cict.fr, florens@cict.fr

  ## It was modified by Jeffrey S. Racine <racinej@mcmaster.ca> and all
  ## errors remain my responsibility. I am indebted to Samuele and the
  ## Toulouse School of Economics for their generous hospitality.

  ## First we require two functions, the first that conducts Regularized
  ## Tikhonov Regression' (aka Ridge Regression)

  ## This function conducts regularized Tikhonov regression which
  ## corresponds to (3.9) in Feve & Florens (2010).

  ## This function accepts as arguments

  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)

  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011, forthcoming Econometrica)
  ## it maps W (the instrument) onto the space of the endogenous
  ## regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  tikh <- function(alpha,CZ,CY,Cr.r,cholesky=FALSE){
      if(cholesky) {
          return(chol2inv(chol(diag(alpha, nrow(CY)) + CY%*%CZ)) %*% Cr.r)
      } else {
          return(solve(diag(alpha, nrow(CY)) + CY%*%CZ) %*% Cr.r)
      }
  }

  ## Samuele indicates alternate form for estimator (visit to SUNY
  ## Stony Brook Feb 24 2015) that can be used with evaluation
  ## data. There is no need to carry around two versions of the same
  ## function, so with some thought we could jettison the above and
  ## use this throughout.

  tikh.eval <- function(alpha,CZ,CY,CY.eval,r,cholesky=FALSE){
    if(cholesky) {
        return(CY.eval%*%chol2inv(chol(diag(alpha, nrow(CY)) + CZ%*%CY)) %*% r)
    } else {
        return(CY.eval%*%solve(diag(alpha, nrow(CY)) + CZ%*%CY) %*% r)
    }
  }

  ## This function applies the iterated Tikhonov approach which
  ## corresponds to (3.10) in Feve & Florens (2010).

  ## This function accepts as arguments

  ## alpha: penalty
  ## CZ:    row-normalized kernel weights for the `independent' variable
  ## CY:    row-normalized kernel weights for the `dependent' variable
  ## Cr:    row-normalized kernel weights for the `instrument/endogenous' variable (see NOTE below)
  ## r:     vector of conditional expectations (z can be E(Z|z) - see NOTE below)

  ## NOTE: for Cr, in the transformation model case treated in Feve &
  ## Florens (2010) this maps Z onto the Y space. In the IV case
  ## (Darrolles, Fan, Florens & Renault (2011, forthcoming Econometrica)
  ## it maps W (the instrument) onto the space of the endogenous
  ## regressor Z.

  ## NOTE: for r, in the transformation model it will be equivalent to
  ## the vector of exogenous covariates, and in the endogenous case r is
  ## the conditional mean of y given the instrument W.

  ## This function returns TBA (need better error checking!)

  ## phi:   the vector of estimated values for the unknown function at the evaluation points

  ## SSalpha: (scalar) value of the sum of square residuals criterion
  ## which is a function of alpha (see (3.10) of Feve & Florens (2010)

  ittik <- function(alpha,CZ,CY,Cr.r,r,cholesky=FALSE) {
      if(cholesky) {
          invmat <- chol2inv(chol(diag(alpha, nrow(CY)) + CY%*%CZ))
      } else {
          invmat <- solve(diag(alpha, nrow(CY)) + CY%*%CZ)
      }
      invmat.Cr.r <- invmat %*% Cr.r
      phi <- invmat.Cr.r + alpha * invmat %*% invmat.Cr.r
      return(sum((CZ%*%phi - r)^2)/NZD(alpha))
  }

  ## This function returns the weight matrix for a local polynomial,
  ## and was rewritten 14/1/15 in Toulouse while visiting JP. It
  ## supports mixed data types. Basic error checking is
  ## undertaken. deriv = 0, strips off weights for mean, = p partials
  ## up to order 2. No cross-partials in this one. Basically useful
  ## for univariate case when deriv > 0 though could be refined - the
  ## old function was slower but had more capability (that basically
  ## went unused).

  ## Update - from ?npksum, "The option permutation.operator= can be
  ## used to `mix and match' operator strings to create a `hybrid'
  ## kernel, in addition to the kernel sum with no operators applied,
  ## one for each continuous dimension in the data. For example, for a
  ## two-dimensional data frame of numeric datatypes,
  ## permutation.operator=c("derivative") will return the usual kernel
  ## sum as if operator = c("normal","normal") in the ksum member, and
  ## in the p.ksum member, it will return kernel sums for operator =
  ## c("derivative","normal"), and operator =
  ## c("normal","derivative"). This makes the computation of gradients
  ## much easier."

  ## So, the upshot is that I could, for the multivariate case, add
  ## the derivative stuff.

  Kmat.lp <- function(deriv=0,
                      mydata.train=NULL,
                      mydata.eval=NULL,
                      bws=NULL,
                      p=0,
                      shrink=TRUE,
                      warning.immediate=TRUE,
                      ...) {

      ## 14/1/15, Toulouse - note that the weights herein ** DO NOT **
      ## shrink towards the lc estimator (neither for the function nor
      ## derivatives), unlike the function returned in
      ## glpreg(). However, they all appear to agree with the previous
      ## Kmat.lp with ** also ** did not shrink towards the lc
      ## estimator. This is noticeably faster, which ought to render
      ## Tikhonov faster as well.

      ## Basic error checking...

      if(is.null(mydata.train)) stop("You must provide training data")
      if(is.null(mydata.eval)) mydata.eval <- mydata.train
      if(is.null(bws)) stop("You must provide bandwidths")

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

      if(k > 0) {
          X.train.numeric <- as.data.frame(X.train[,X.col.numeric])
          X.eval.numeric <- as.data.frame(X.eval[,X.col.numeric])
      }

      if(deriv<0||deriv>2)
          stop(paste("Error: deriv= (integer) is invalid\n[min = ", 0, ", max = ",  p, "]\n",sep=""))

      if(length(p) != 1) {
          p <- p[1]
      }

      if(p < 0)
          stop(paste("Error: p (order of polynomial) must be a non-negative integer\np is (", p, ")\n",sep=""))

      K.x <- npksum(txdat=X.train,
                    exdat=X.eval,
                    bws=bws,
                    return.kernel.weights=TRUE,
                    ...)$kw

      if(p==0) {

          ## No shrinking necessary for local constant estimator

          if(deriv==0) {

              Kmat <- t(K.x)/NZD(rowSums(t(K.x)))

          } else if(deriv==1) {

              ## Note this is not general XXX Feb 25 2015, for
              ## univariate z only

              K.x.deriv <- npksum(txdat=X.train,
                                  exdat=X.eval,
                                  bws=bws,
                                  return.kernel.weights=TRUE,
                                  operator="derivative",
                                  ...)$kw/NZD(bws)

              rSk <- NZD(rowSums(t(K.x)))

              Kmat <- t(K.x.deriv)/NZD(rSk)-t(K.x)/NZD(rSk)*(rowSums(t(K.x.deriv))/NZD(rSk))

          }

      }

      if(p > 0) {

          ## Re-use this matrix, shrinking occurs here

          W.z <- W.glp(xdat=X.train.numeric,
                             degree=rep(p,NCOL(X.train.numeric)))

          if(is.null(mydata.eval)) {
              ## Guess we could avoid copy with conditional statement below using either W.z or W.z.eval
              W.z.eval <- W.z
          } else {
              W.z.eval <- W.glp(xdat=X.train.numeric,
                                      exdat=as.data.frame(X.eval.numeric),
                                      degree=rep(p,NCOL(X.train.numeric)))
          }

          nc <- ncol(W.z)

          WzkWz.inv <- list()

          for(i in 1:ncol(K.x)) {

              if(tryCatch(WzkWz.inv[[i]] <- as.matrix(chol2inv(chol(t(W.z)%*%(K.x[,i]*W.z)))),
                          error = function(e){
                              return(matrix(FALSE,nc,nc))
                          })[1,1]!=FALSE) {

              } else {

                  if(shrink==FALSE) {

                      ## If we do not explicitly engage ridging then we do not fail
                      ## and terminate, rather, we return NA when Wmat.sum is
                      ## singular

                      Kmat <- NA

                  } else {

                      ## Ridging

                      epsilon <- 1/n.train
                      ridge <- 0

                      while(tryCatch(as.matrix(chol2inv(chol((chol(t(W.z)%*%(K.x[,i]*W.z)+diag(rep(ridge,nc))))))),
                                     error = function(e){
                                         return(matrix(FALSE,nc,nc))
                                     })[1,1]==FALSE) {
                          ridge <- ridge + epsilon
                      }

                      WzkWz.inv[[i]] <- as.matrix(chol2inv(chol(t(W.z)%*%(K.x[,i]*W.z)+diag(rep(ridge,nc)))))

                      warning(paste("Ridging obs. ", i, ", ridge = ", signif(ridge,6),sep=""),
                              immediate.=warning.immediate,
                              call.=!warning.immediate)

                  }

              }

          }
      }

      if(p==1) {

          if(deriv==0) {
              Kmat <- matrix(NA, n.eval, n.train)
              for(i in 1:n.eval) {
                  Kmat[i,] <- W.z.eval[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]
              }
          }
          if(deriv==1) {
              W.z.deriv.1 <- W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 1)

              Kmat <- matrix(NA, n.eval, n.train)
              for(i in 1:n.eval) {
                  Kmat[i,] <- W.z.deriv.1[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]
              }
          }

      }

      if(p >= 2) {

          Kmat <- matrix(NA, n.eval, n.train)
          if(deriv==0) {
              for(i in 1:n.eval) {
                  Kmat[i,] <- W.z.eval[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]
              }
          }


          if(deriv==1) {

              W.z.deriv.1 <- W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 1)

              for(i in 1:n.eval) {
                  Kmat[i,] <- W.z.deriv.1[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]
              }

          }

          if(deriv==2) {

              W.z.deriv.2 <- W.glp(xdat=X.train.numeric,
                                         exdat=as.matrix(X.eval.numeric),
                                         degree=rep(p,NCOL(X.train.numeric)),
                                         gradient.vec = 2)

              for(i in 1:n.eval) {
                  Kmat[i,] <- W.z.deriv.2[i,,drop=FALSE]%*%WzkWz.inv[[i]]%*%t(W.z)*K.x[,i]
              }

          }

      }

      return(Kmat)
  }

  glpreg <- function(tydat=NULL,
                     txdat=NULL,
                     exdat=NULL,
                     bws=NULL,
                     degree=NULL,
                     leave.one.out=FALSE,
                     deriv=1,
                     ...) {

    ## Don't think this error checking is robust

    if(is.null(tydat)) stop("Error: You must provide y data")
    if(is.null(txdat)) stop("Error: You must provide X data")
    if(is.null(bws)) stop("Error: You must provide a bandwidth object")
    if(is.null(degree) | any(degree < 0)) stop(paste("Error: degree vector must contain non-negative integers\ndegree is (", degree, ")\n",sep=""))
    if(p>0 && (any(deriv < 0) || any(deriv > degree)))
      stop("deriv must lie between 0 and degree")

    miss.ex = missing(exdat)

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

      grad <- gradients(npreg(tydat=tydat,
                              txdat=txdat,
                              exdat = exdat,
                              bws = bws,
                              ukertype="liracine",
                              okertype="liracine",
                              gradients=TRUE,
                              ...))

      return(list(mean = mhat,
                  grad = grad))

    } else {

      W <- W.glp(xdat=txdat,
                       degree=degree)

      W.eval <- W.glp(xdat=txdat,
                            exdat=exdat,
                            degree=degree)

      W.eval.deriv <- W.glp(xdat=txdat,
                                  exdat=exdat,
                                  degree=degree,
                                  gradient.vec=rep(deriv,NCOL(txdat)))

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
      ridge.lc <- double(n.eval)
      doridge <- !logical(n.eval)

      nc <- ncol(tww[,,1])
      I.nc <- diag(nc)

      ## Test for singularity of the generalized local polynomial
      ## estimator, shrink the mean towards the local constant mean.

      ridger <- function(i) {
        doridge[i] <<- FALSE
        ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
        tryCatch(chol2inv(chol(tww[,,i] + ridge[i]*I.nc))%*%tyw[,i],
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
        (1-ridge[i])*W.eval[i,, drop = FALSE] %*% coef.mat[,i] + ridge.lc[i]
      })

      grad <- sapply(1:n.eval, function(i) {W.eval.deriv[i,-1, drop = FALSE] %*% coef.mat[-1,i]})

      return(list(mean = mhat,
                  grad = grad))

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

        if (!any(is.nan(mean.loo)) && !any(mean.loo == maxPenalty)){
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
        ridge.lc <- double(n)
        doridge <- !logical(n)

        nc <- ncol(tww[,,1])

        ## Test for singularity of the generalized local polynomial
        ## estimator, shrink the mean towards the local constant mean.

        ridger <- function(i) {
          doridge[i] <<- FALSE
          ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
          W[i,, drop = FALSE] %*% tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
                  error = function(e){
                    ridge[i] <<- ridge[i]+epsilon
                    doridge[i] <<- TRUE
                    return(rep(maxPenalty,nc))
                  })
        }

        while(any(doridge)){
          iloo <- (1:n)[doridge]
          mean.loo[iloo] <- (1-ridge[iloo])*sapply(iloo, ridger) + ridge.lc[iloo]
        }

        if (!any(is.nan(mean.loo)) && !any(mean.loo == maxPenalty)){
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
        ridge.lc <- double(n)
        doridge <- !logical(n)

        nc <- ncol(tww[,,1])

        ## Test for singularity of the generalized local polynomial
        ## estimator, shrink the mean towards the local constant mean.

        ridger <- function(i) {
          doridge[i] <<- FALSE
          ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
          W[i,, drop = FALSE] %*% tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
                  error = function(e){
                    ridge[i] <<- ridge[i]+epsilon
                    doridge[i] <<- TRUE
                    return(rep(maxPenalty,nc))
                  })
        }

        while(any(doridge)){
          ii <- (1:n)[doridge]
          ghat[ii] <- (1-ridge[ii])*sapply(ii, ridger) + ridge.lc[ii]
        }

        trH <- kernel.i.eq.j*sum(sapply(1:n,function(i){
          (1-ridge[i])*W[i,, drop = FALSE] %*% chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc)))) %*% t(W[i,, drop = FALSE]) + ridge[i]/NZD(tww[,,i][1,1])
        }))

        aic.penalty <- (1+trH/n)/(1-(trH+2)/n)

        if (!any(ghat == maxPenalty) & (aic.penalty > 0)){
          fv <- log(mean((ydat-ghat)^2)) + aic.penalty
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
                    nmulti=nmulti,
                    random.seed=42,
                    optim.maxattempts = 10,
                    optim.method=c("Nelder-Mead", "BFGS", "CG"),
                    optim.reltol=sqrt(.Machine$double.eps),
                    optim.abstol=.Machine$double.eps,
                    optim.maxit=500,
                    debug=FALSE,
                    bw.init=NULL,
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

    W <- W.glp(xdat=xdat,
                     degree=degree)

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

      if(iMulti == 1 && !is.null(bw.init)) {
        init.search.vals <- bw.init
      } else {

        ## First initialize to values for factors (`liracine' kernel)

        init.search.vals <- runif(ncol(xdat),0,1)

        for(i in 1:ncol(xdat)) {
          if(xdat.numeric[i]==TRUE) {
            init.search.vals[i] <- runif(1,.5,1.5)*EssDee(xdat[,i])*nrow(xdat)^{-1/(4+num.numeric)}
          }
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
            init.search.vals[i] <- runif(1,.5,1.5)*EssDee(xdat[,i])*nrow(xdat)^{-1/(4+num.numeric)}
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
            init.search.vals[i] <- runif(1,.5,1.5)*EssDee(xdat[,i])*nrow(xdat)^{-1/(4+num.numeric)}
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

    return(list(bw=bw.opt,
                fv=fv,
                numimp=numimp,
                best=best,
                fv.vec=fv.vec))

  }

  ## Here is where the function `npregiv' really begins:

  console <- newLineConsole()

  ## Basic error checking

  if(!is.logical(penalize.iteration)) stop("penalize.iteration must be logical (TRUE/FALSE)")
  if(!is.logical(smooth.residuals)) stop("smooth.residuals must be logical (TRUE/FALSE)")
  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")
  if(!is.logical(iterate.Tikhonov)) stop("iterate.Tikhonov must be logical (TRUE/FALSE)")
  if(iterate.Tikhonov.num < 1) stop("iterate.Tikhonov.num must be a positive integer")  

  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(iterate.diff.tol < 0) stop("iterate.diff.tol must be non-negative")
  if(constant <= 0 || constant >=1) stop("constant must lie in (0,1)")
  if(p < 0) stop("p must be a non-negative integer")

  if(!is.null(alpha) && alpha <= 0) stop("alpha must be positive")
  if(!is.null(alpha.iter) && alpha.iter <= 0) stop("alpha.iter must be positive")

  if(return.weights.phi.deriv.1 && !return.weights.phi) stop("must use return.weights.phi=TRUE when using return.weights.phi.deriv.1=TRUE")
  if(return.weights.phi.deriv.2 && !return.weights.phi) stop("must use return.weights.phi=TRUE when using return.weights.phi.deriv.2=TRUE")
  if(return.weights.phi.deriv.2 && p<2) stop("must use p >= 2 when using return.weights.phi.deriv.2=TRUE")

  start.from <- match.arg(start.from)
  method <- match.arg(method)

  nmulti.loop <- if(!is.null(nmulti)) nmulti else 1
  nmulti <- if(!is.null(nmulti)) nmulti else 5

  ## Need to determine how many x, w, z are numeric

  z <- data.frame(z)
  w <- data.frame(w)
  if(!is.null(x)) {
      z <- data.frame(z,x)
      ## JP points out that, with exogenous predictors, they must be
      ## part of both z and the instruments. The line below was added
      ## 20/1/15 in Toulouse.
      w <- data.frame(w,x)
      ## Obviously, if you have exogenous variables that are only in
      ## the instrument set, you can trivially accommodate this
      ## (append to w before invoking the function - added to man
      ## page)
      if(!is.null(zeval)&&!is.null(xeval)) zeval <- data.frame(zeval,xeval)
  }

  z.numeric <- sapply(1:NCOL(z),function(i){is.numeric(z[,i])})
  num.z.numeric <- NCOL(as.data.frame(z[,z.numeric]))

  w.numeric <- sapply(1:NCOL(w),function(i){is.numeric(w[,i])})
  num.w.numeric <- NCOL(as.data.frame(w[,w.numeric]))

  if(method=="Tikhonov") {

    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the
    ## first step we use E(y|w) as a proxy for phi(z) and use
    ## bandwidths for y on w).

    ## Convergence value returned for Landweber-Fridman but value
    ## required

    ## convergence <- NULL

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(bw)) {
        console <- printPush("Computing bandwidths and E(y|w)...",console)
    } else {
        console <- printPush("Computing E(y|w) using supplied bandwidths...",console)
    }

    if(is.null(bw)) {
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

        bw.E.y.w <- hyw$bw
    } else {
        bw.E.y.w <- bw$bw.E.y.w
    }

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix and E(y|w)...", console)

    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    bws=bw.E.y.w,
                    degree=rep(p, num.w.numeric),
                    ...)$mean

    KYW <- Kmat.lp(mydata.train=data.frame(w),
                   bws=bw.E.y.w,
                   p=rep(p, num.w.numeric),
                   ...)

    ## We conduct local polynomial kernel regression of E(y|w) on z

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(bw)) {
        console <- printPush("Computing bandwidths for E(E(y|w)|z)...", console)
    } else {
        console <- printPush("Computing E(E(y|w)|z) using supplied bandwidths...", console)
    }

    if(is.null(bw)) {
        hywz <- glpcv(ydat=E.y.w,
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

        bw.E.E.y.w.z <- hywz$bw
    } else {
        bw.E.E.y.w.z <- bw$bw.E.E.y.w.z
    }

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix and E(E(y|w)|z)...", console)

    E.E.y.w.z <- glpreg(tydat=E.y.w,
                        txdat=z,
                        bws=bw.E.E.y.w.z,
                        degree=rep(p, num.z.numeric),
                        ...)$mean

    KYWZ <- Kmat.lp(mydata.train=data.frame(z),
                    bws=bw.E.E.y.w.z,
                    p=rep(p, num.z.numeric),
                    ...)

    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov function) to
    ## determine the optimal alpha for the non-iterated scheme. Note
    ## that the function `optimize' accepts bounds on the search (in
    ## this case alpha.min to alpha.max))

    ## E(r|z)=E(E(phi(z)|w)|z)
    ## \phi^\alpha = (\alpha I+CzCw)^{-1}Cr x r

    if(!is.null(bw)) alpha <- bw$alpha

    if(is.null(alpha)&&is.null(bw)) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("Numerically solving for alpha...", console)
      alpha <- optimize(ittik, c(alpha.min, alpha.max), tol = alpha.tol, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    }

    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha.

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing initial phi(z) estimate...", console)
    phi <- as.vector(tikh(alpha, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z))
    phi.mat <- phi
    
    phi.eval.mat <- NULL
    if(!is.null(zeval)) {
        ## If there is evaluation data, KPHIWZ and KPHIWZ.eval will
        ## differ...
        
        KPHIWZ.eval <- Kmat.lp(mydata.train=data.frame(z),
                               mydata.eval=data.frame(z=zeval),
                               bws=bw.E.E.y.w.z,
                               p=rep(p, num.z.numeric),
                               ...)
        
        phi.eval <- as.vector(tikh.eval(alpha, CZ = KYW, CY = KYWZ, CY.eval = KPHIWZ.eval, r = E.y.w))
        phi.eval.mat <- cbind(phi.eval.mat,phi.eval)
    }    

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(bw)) {
        console <- printPush("Computing bandwidths for E(phi(z)|w)...", console)
    } else {
        console <- printPush("Computing E(phi(z)|w) using supplied bandwidths...", console)
    }

    bw.E.phi.w <- NULL
    bw.E.E.phi.w.z <- NULL

    for(i in 1:iterate.Tikhonov.num) {

      if(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num) {
          console <- printClear(console)
          console <- printPop(console)
          console <- printPush(paste("Iteration ",i," of ",iterate.Tikhonov.num,sep=""), console)
      }

      if(is.null(bw)) {
          hphiw <- glpcv(ydat=phi, ## 23/1/15 phi is sample
                         xdat=w,
                         degree=rep(p, num.w.numeric),
                         nmulti=nmulti.loop,
                         random.seed=random.seed,
                         optim.maxattempts=optim.maxattempts,
                         optim.method=optim.method,
                         optim.reltol=optim.reltol,
                         optim.abstol=optim.abstol,
                         optim.maxit=optim.maxit,
                         bw.init=bw.E.phi.w,
                         ...)

          bw.E.phi.w <- hphiw$bw
      } else {
          bw.E.phi.w <- bw$bw.E.phi.w
      }

      if(!(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num)) {
          console <- printClear(console)
          console <- printPop(console)
          console <- printPush("Computing weight matrix for E(phi(z)|w)...", console)
      }

      E.phi.w <- glpreg(tydat=phi,
                        txdat=w,
                        bws=bw.E.phi.w,
                        degree=rep(p, num.w.numeric),
                        ...)$mean

      KPHIW <- Kmat.lp(mydata.train=data.frame(w),
                       bws=bw.E.phi.w,
                       p=rep(p, num.w.numeric),
                       ...)

      if(!(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num)) {
          console <- printClear(console)
          console <- printPop(console)
          if(is.null(bw)) {
              console <- printPush("Computing bandwidths for E(E(phi(z)|w)|z)...", console)
          } else {
              console <- printPush("Computing E(E(phi(z)|w)|z) using supplied bandwidths...", console)
          }
      }

      if(is.null(bw)) {
          hphiwz <- glpcv(ydat=E.phi.w,
                          xdat=z,
                          degree=rep(p, num.z.numeric),
                          nmulti=nmulti.loop,
                          random.seed=random.seed,
                          optim.maxattempts=optim.maxattempts,
                          optim.method=optim.method,
                          optim.reltol=optim.reltol,
                          optim.abstol=optim.abstol,
                          optim.maxit=optim.maxit,
                          bw.init=bw.E.E.phi.w.z,
                          ...)

          bw.E.E.phi.w.z <- hphiwz$bw
      } else {
          bw.E.E.phi.w.z <- bw$bw.E.E.phi.w.z
      }

      E.E.phi.w.z <- glpreg(tydat=E.y.w,
                            txdat=z,
                            bws=bw.E.E.phi.w.z,
                            degree=rep(p, num.z.numeric),
                            ...)$mean

      if(!(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num)) {
          console <- printClear(console)
          console <- printPop(console)
          console <- printPush("Computing weight matrix for E(E(phi(z)|w)|z)...", console)
      }
      
      KPHIW <- Kmat.lp(mydata.train=data.frame(w),
                       bws=bw.E.phi.w,
                       p=rep(p, num.w.numeric),
                       ...)

      KPHIWZ <- Kmat.lp(mydata.train=data.frame(z),
                        bws=bw.E.E.phi.w.z,
                        p=rep(p, num.z.numeric),
                        ...)

      ## Next, we minimize the function ittik to obtain the optimal value
      ## of alpha (here we use the iterated Tikhonov approach) to
      ## determine the optimal alpha for the non-iterated scheme.

      if(!is.null(bw)) alpha.iter <- bw$alpha.iter

      if(!iterate.Tikhonov) {
          alpha.iter <- alpha
      } else {

          if(is.null(alpha.iter)&&is.null(bw)) {
              if(!(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num)) {
                  console <- printClear(console)
                  console <- printPop(console)
                  console <- printPush(paste("Iterating and recomputing the numerical solution for alpha (iteration ",i," of ",iterate.Tikhonov.num,")",sep=""), console)
              }
              alpha.iter <- optimize(ittik, c(alpha.min, alpha.max), tol = alpha.tol, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.phi.w.z, r = E.y.w)$minimum
          }
      }

      ## Finally, we conduct regularized Tikhonov regression using this
      ## optimal alpha and the updated bandwidths.
      
      if(!(iterate.Tikhonov.num > 1 && i < iterate.Tikhonov.num)) {
          console <- printClear(console)
          console <- printPop(console)
          console <- printPush("Computing final phi(z) estimate...", console)
      }

      phi <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ, r = E.y.w))
      phi.mat <- cbind(phi.mat,phi)    

      H <- NULL
      if(return.weights.phi) {
          H <- KPHIWZ%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
      }
      
      ## First derivative
      
      KPHIWZ.deriv.1 <- Kmat.lp(deriv=1,
                                mydata.train=data.frame(z),
                                bws=bw.E.E.phi.w.z,
                                p=rep(p, num.z.numeric),
                                ...)
      
      phi.deriv.1 <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ.deriv.1, r = E.y.w))
      
      H.deriv.1 <- NULL
      if(return.weights.phi.deriv.1) {
          H.deriv.1 <- KPHIWZ.deriv.1%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
      }
      
      ## Second derivative
      
      phi.deriv.2 <- NULL
      H.deriv.2 <- NULL
      
      if(p >= 2) {
          
          KPHIWZ.deriv.2 <- Kmat.lp(deriv=2,
                                    mydata.train=data.frame(z),
                                    bws=bw.E.E.phi.w.z,
                                    p=rep(p, num.z.numeric),
                                    ...)
          
          phi.deriv.2 <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ.deriv.2, r = E.y.w))
          
          
          if(return.weights.phi.deriv.2) {
              H.deriv.2 <- KPHIWZ.deriv.2%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
          }
          
      }
      
      ## If evaluation data are provided...
      
      phi.eval <- NULL
      phi.deriv.eval.1 <- NULL
      phi.deriv.eval.2 <- NULL
      H.eval <- NULL
      H.deriv.eval.1 <- NULL
      H.deriv.eval.2 <- NULL
      
      if(!is.null(zeval)) {
          ## If there is evaluation data, KPHIWZ and KPHIWZ.eval will
          ## differ...
          
          KPHIWZ.eval <- Kmat.lp(mydata.train=data.frame(z),
                                 mydata.eval=data.frame(z=zeval),
                                 bws=bw.E.E.phi.w.z,
                                 p=rep(p, num.z.numeric),
                                 ...)
          
          phi.eval <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ.eval, r = E.y.w))
          phi.eval.mat <- cbind(phi.eval.mat,phi.eval)
          
          if(return.weights.phi) {
              H.eval <- KPHIWZ.eval%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
          }
          
          KPHIWZ.eval.deriv.1 <- Kmat.lp(deriv=1,
                                         mydata.train=data.frame(z),
                                         mydata.eval=data.frame(z=zeval),
                                         bws=bw.E.E.phi.w.z,
                                         p=rep(p, num.z.numeric),
                                         ...)
          
          phi.deriv.eval.1 <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ.eval.deriv.1, r = E.y.w))
          
          if(return.weights.phi.deriv.1) {
              H.deriv.eval.1 <- KPHIWZ.eval.deriv.1%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
          }
          
          if(p >= 2) {
              
              KPHIWZ.eval.deriv.2 <- Kmat.lp(deriv=2,
                                             mydata.train=data.frame(z),
                                             mydata.eval=data.frame(z=zeval),
                                             bws=bw.E.E.phi.w.z,
                                             p=rep(p, num.z.numeric),
                                             ...)
              
              phi.deriv.eval.2 <- as.vector(tikh.eval(alpha.iter, CZ = KPHIW, CY = KPHIWZ, CY.eval = KPHIWZ.eval.deriv.2, r = E.y.w))
              
              if(return.weights.phi.deriv.2) {
                  H.deriv.eval.2 <- KPHIWZ.eval.deriv.2%*%solve(alpha.iter*diag(nrow(KPHIWZ)) + KPHIW%*%KPHIWZ)%*%KYW
              }
              
          }
      }
      
    }
    
    console <- printClear(console)
    console <- printPop(console)

    if((alpha.iter-alpha.min)/NZD(alpha.min) < 0.01) warning(paste("Tikhonov parameter alpha (",formatC(alpha.iter,digits=4,format="f"),") is close to the search minimum (",alpha.min,")",sep=""))
    if((alpha.max-alpha.iter)/NZD(alpha.max) < 0.01) warning(paste("Tikhonov parameter alpha (",formatC(alpha.iter,digits=4,format="f"),") is close to the search maximum (",alpha.max,")",sep=""))

    ret <- list(phi=phi,
                phi.eval=phi.eval,
                phi.mat=phi.mat,
                phi.eval.mat=phi.eval.mat,
                phi.deriv.1=as.matrix(phi.deriv.1),
                phi.deriv.eval.1=if(!is.null(phi.deriv.eval.1)){as.matrix(phi.deriv.eval.1)}else{NULL},
                phi.deriv.2=if(!is.null(phi.deriv.2)){as.matrix(phi.deriv.2)}else{NULL},
                phi.deriv.eval.2=if(!is.null(phi.deriv.eval.2)){as.matrix(phi.deriv.eval.2)}else{NULL},
                phi.weights=H,
                phi.deriv.1.weights=H.deriv.1,
                phi.deriv.2.weights=H.deriv.2,
                phi.eval.weights=H.eval,
                phi.deriv.eval.1.weights=H.deriv.eval.1,
                phi.deriv.eval.2.weights=H.deriv.eval.2,
                alpha=alpha,
                alpha.iter=alpha.iter,
                bw.E.y.w=bw.E.y.w,
                bw.E.E.y.w.z=bw.E.E.y.w.z,
                bw.E.phi.w=bw.E.phi.w,
                bw.E.E.phi.w.z=bw.E.E.phi.w.z,
                call=cl,
                y=y,
                z=z,
                w=w,
                x=x,
                zeval=zeval,
                xeval=xeval,
                p=p,
                nmulti=nmulti,
                method=method,
                ptm=proc.time() - ptm.start)
    class(ret) <- "npregiv"
    return(ret)

  } else {

    ## Landweber-Fridman

    ## For the stopping rule

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(bw)) {
        console <- printPush(paste("Computing bandwidths and E(y|w) for stopping rule...",sep=""),console)
    } else {
        console <- printPush(paste("Computing E(y|w) for stopping rule using supplied bandwidths...",sep=""),console)
    }

    norm.stop <- numeric()

    if(is.null(bw)) {
        h <- glpcv(ydat=y,
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
        bw.E.y.w <- h$bw
    } else {
        bw.E.y.w <- bw$bw.E.y.w
    }

    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    bws=bw.E.y.w,
                    degree=rep(p, num.w.numeric),
                    ...)$mean

    if(return.weights.phi) {

        if(p<0) stop("glp return weights not supported")
        if(NCOL(z) > 1) stop("dimension of z must be one for currently supported return weights")

        T.mat.r <- Kmat.lp(mydata.train=data.frame(w),
                           bws=bw.E.y.w,
                           p=rep(p, num.w.numeric),
                           ...)

    }

    ## We begin the iteration computing phi.0 and phi.1 directly, then
    ## iterate.

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(bw)) {
        console <- printPush(paste("Computing bandwidths and E(y|z) for iteration 0...",sep=""),console)
    } else {
        console <- printPush(paste("Computing E(y|z) for iteration 0 using supplied bandwidths...",sep=""),console)
    }

    if(is.null(starting.values)) {

      if(is.null(bw)) {
          h <- glpcv(ydat=if(start.from=="Eyz") y else E.y.w,
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
          bw.E.y.z <- h$bw
      } else {
          bw.E.y.z <- bw$bw.E.y.z
      }

      g <- glpreg(tydat=if(start.from=="Eyz") y else E.y.w,
                  txdat=z,
                  bws=bw.E.y.z,
                  degree=rep(p, num.z.numeric),
                  ...)

      phi.0 <- g$mean
      phi.0.deriv.1 <- g$grad
      if(p >= 2) {
          phi.0.deriv.2 <- glpreg(tydat=if(start.from=="Eyz") y else E.y.w,
                                  txdat=z,
                                  bws=bw.E.y.z,
                                  degree=rep(p, num.z.numeric),
                                  deriv=2,
                                  ...)$grad
      }

      if(!is.null(zeval)) {
          g <- glpreg(tydat=if(start.from=="Eyz") y else E.y.w,
                      txdat=z,
                      exdat=zeval,
                      bws=bw.E.y.z,
                      degree=rep(p, num.z.numeric),
                      ...)

          phi.eval.0 <- g$mean
          phi.eval.0.deriv.1 <- g$grad
          if(p >= 2) {
              phi.eval.0.deriv.2 <- glpreg(tydat=if(start.from=="Eyz") y else E.y.w,
                                           txdat=z,
                                           exdat=zeval,
                                           bws=bw.E.y.z,
                                           degree=rep(p, num.z.numeric),
                                           deriv=2,
                                           ...)$grad
          }
      } else {
          phi.eval.0 <- NULL
          phi.eval.0.deriv.1 <- NULL
          phi.eval.0.deriv.2 <- NULL
      }

      if(return.weights.phi) {

          H <- Kmat.lp(mydata.train=data.frame(z),
                       bws=bw.E.y.z,
                       p=rep(p, num.z.numeric),
                       ...)

          if(!is.null(zeval)) H.eval <- Kmat.lp(mydata.train=data.frame(z),
                                                mydata.eval=data.frame(z=zeval),
                                                bws=bw.E.y.z,
                                                p=rep(p, num.z.numeric),
                                                ...)

          if(p==0 || p==1) {

              if(return.weights.phi.deriv.1) {

                  H.deriv.1 <- Kmat.lp(deriv=1,
                                       mydata.train=data.frame(z),
                                       bws=bw.E.y.z,
                                       p=rep(p, num.z.numeric),
                                       ...)

                  if(!is.null(zeval)) H.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                                mydata.train=data.frame(z),
                                                                mydata.eval=data.frame(z=zeval),
                                                                bws=bw.E.y.z,
                                                                p=rep(p, num.z.numeric),
                                                                ...)
              }

          } else {

              if(return.weights.phi.deriv.1) {

                  H.deriv.1 <- Kmat.lp(deriv=1,
                                       mydata.train=data.frame(z),
                                       bws=bw.E.y.z,
                                       p=rep(p, num.z.numeric),
                                       ...)
                  if(!is.null(zeval)) {

                      H.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                mydata.train=data.frame(z),
                                                mydata.eval=data.frame(z=zeval),
                                                bws=bw.E.y.z,
                                                p=rep(p, num.z.numeric),
                                                ...)

                  }

              }

              if(return.weights.phi.deriv.2) {

                  H.deriv.2 <- Kmat.lp(deriv=2,
                                       mydata.train=data.frame(z),
                                       bws=bw.E.y.z,
                                       p=rep(p, num.z.numeric),
                                       ...)

                  if(!is.null(zeval)) {

                      H.deriv.eval.2 <- Kmat.lp(deriv=2,
                                                mydata.train=data.frame(z),
                                                mydata.eval=data.frame(z=zeval),
                                                bws=bw.E.y.z,
                                                p=rep(p, num.z.numeric),
                                                ...)

                  }

              }

          }

      }

    } else {

      ## Starting values input by user

      phi.0 <- starting.values

      if(return.weights.phi)  H <- NULL
      bw.E.y.z <- NULL

    }

    starting.values.phi <- phi.0

    console <- printClear(console)
    console <- printPop(console)
    if(smooth.residuals) {
        if(is.null(bw)) {
            console <- printPush(paste("Computing bandwidths and E[y-phi(z)|w] for iteration 1...",sep=""),console)
        } else {
            console <- printPush(paste("Computing E[y-phi(z)|w] for iteration 1 using supplied bandwidths...",sep=""),console)
        }
    } else {
        if(is.null(bw)) {
            console <- printPush(paste("Computing bandwidths and E[phi(z)|w] for iteration 1...",sep=""),console)
        } else {
            console <- printPush(paste("Computing E[phi(z)|w] for iteration 1 using supplied bandwidths...",sep=""),console)
        }
    }

    if(smooth.residuals) {

      resid <- y - phi.0

      if(is.null(bw)) {
          h <- glpcv(ydat=resid,
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
          bw.resid.w <- h$bw
      } else {
          bw.resid.w <- bw$bw.resid.w[1,]
      }

      resid.fitted <- glpreg(tydat=resid,
                             txdat=w,
                             bws=bw.resid.w,
                             degree=rep(p, num.w.numeric),
                             ...)$mean

    } else {

      if(is.null(bw)) {
          h <- glpcv(ydat=phi.0,
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
          bw.resid.w <- h$bw
      } else {
          bw.resid.w <- bw$bw.resid.w[1,]
      }

      resid.fitted <- E.y.w - glpreg(tydat=phi.0,
                                     txdat=w,
                                     bws=bw.resid.w,
                                     degree=rep(p, num.w.numeric),
                                     ...)$mean

    }

    if(return.weights.phi) {

        T.mat <- Kmat.lp(mydata.train=data.frame(w),
                         bws=bw.resid.w,
                         p=rep(p, num.z.numeric),
                         ...)

    }

    norm.stop[1] <- sum(resid.fitted^2)/NZD_pos(sum(E.y.w^2))

    console <- printClear(console)
    console <- printPop(console)
    if(smooth.residuals) {
        if(is.null(bw)) {
            console <- printPush(paste("Computing bandwidths and E[E(y-phi(z)|w)|z] for iteration 1...",sep=""),console)
        } else {
            console <- printPush(paste("Computing E[E(y-phi(z)|w)|z] for iteration 1 using supplied bandwidths...",sep=""),console)
        }
    } else {
        if(is.null(bw)) {
            console <- printPush(paste("Computing bandwidths and E[E(y|w) - E(phi(z)|w)|z] for iteration 1...",sep=""),console)
        } else {
            console <- printPush(paste("Computing E[E(y|w) - E(phi(z)|w)|z] for iteration 1 using supplied bandwidths...",sep=""),console)
        }
    }

    if(is.null(bw)) {
        h <- glpcv(ydat=resid.fitted,
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
        bw.resid.fitted.w.z <- h$bw
    } else {
        bw.resid.fitted.w.z <- bw$bw.resid.fitted.w.z[1,]
    }

    g <- glpreg(tydat=resid.fitted,
                txdat=z,
                bws=bw.resid.fitted.w.z,
                degree=rep(p, num.z.numeric),
                ...)

    phi.deriv.1.list <- vector(mode="list", length=iterate.max)
    phi.deriv.eval.1.list <- vector(mode="list", length=iterate.max)

    phi <- phi.0 + constant*g$mean
    phi.deriv.1 <- phi.0.deriv.1 + constant*g$grad

    phi.mat <- matrix(NA, length(phi), iterate.max)
    phi.mat[,1] <- phi
    phi.deriv.1.list[[1]] <- phi.deriv.1

    phi.deriv.2.list <- vector(mode="list", length=iterate.max)
    phi.deriv.eval.2.list <- vector(mode="list", length=iterate.max)
    phi.deriv.2 <- NULL
    phi.deriv.eval.2 <- NULL

    if(p >= 2) {
        phi.deriv.2 <- phi.0.deriv.2 + constant*glpreg(tydat=resid.fitted,
                                                       txdat=z,
                                                       bws=bw.resid.fitted.w.z,
                                                       degree=rep(p, num.z.numeric),
                                                       deriv=2,
                                                       ...)$grad
        phi.deriv.2.list[[1]] <- phi.deriv.2
    }

    if(!is.null(zeval)) {
        g <- glpreg(tydat=resid.fitted,
                    txdat=z,
                    exdat=zeval,
                    bws=bw.resid.fitted.w.z,
                    degree=rep(p, num.z.numeric),
                    ...)
        phi.eval <- phi.eval.0 + constant*g$mean
        phi.eval.mat <- matrix(NA, length(phi.eval), iterate.max)
        phi.eval.mat[,1] <- phi.eval

        phi.deriv.eval.1 <- phi.eval.0.deriv.1 + constant*g$grad
        phi.deriv.eval.1.list[[1]] <- phi.deriv.eval.1

        if(p >= 2) {
            phi.deriv.eval.2 <- phi.eval.0.deriv.2 + constant*glpreg(tydat=resid.fitted,
                                                                     txdat=z,
                                                                     bws=bw.resid.fitted.w.z,
                                                                     exdat=zeval,
                                                                     degree=rep(p, num.z.numeric),
                                                                     deriv=2,
                                                                     ...)$grad
            phi.deriv.eval.2.list[[1]] <- phi.deriv.eval.2
        }



    } else {
        phi.eval <- NULL
        phi.eval.mat <- NULL
        phi.deriv.eval.1 <- NULL
        phi.deriv.eval.1.list <- NULL
    }

    ## Need these list even when no weights for return

    phi.weights.list <- vector(mode="list", length=iterate.max)
    phi.deriv.1.weights.list <- vector(mode="list", length=iterate.max)
    phi.deriv.2.weights.list <- vector(mode="list", length=iterate.max)

    phi.eval.weights.list <- vector(mode="list", length=iterate.max)
    phi.deriv.eval.1.weights.list <- vector(mode="list", length=iterate.max)
    phi.deriv.eval.2.weights.list <- vector(mode="list", length=iterate.max)

    ## Pre-allocate bandwidth matrices
    bw.resid.w.mat <- matrix(NA, iterate.max, length(bw.resid.w))
    bw.resid.w.mat[1,] <- bw.resid.w
    bw.resid.fitted.w.z.mat <- matrix(NA, iterate.max, length(bw.resid.fitted.w.z))
    bw.resid.fitted.w.z.mat[1,] <- bw.resid.fitted.w.z

    if(return.weights.phi) {

        T.mat.adjoint <- Kmat.lp(mydata.train=data.frame(z),
                                 bws=bw.resid.fitted.w.z,
                                 p=rep(p, num.z.numeric),
                                 ...)

        if(!is.null(zeval)) T.mat.adjoint.eval <- Kmat.lp(mydata.train=data.frame(z),
                                                          mydata.eval=data.frame(z=zeval),
                                                          bws=bw.resid.fitted.w.z,
                                                          p=rep(p, num.z.numeric),
                                                          ...)
        if(p==0 || p==1) {

            if(return.weights.phi.deriv.1) {

                T.mat.adjoint.deriv.1 <- Kmat.lp(deriv=1,
                                                 mydata.train=data.frame(z),
                                                 bws=bw.resid.fitted.w.z,
                                                 p=rep(p, num.z.numeric),
                                                 ...)

                if(!is.null(zeval)) T.mat.adjoint.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                                          mydata.train=data.frame(z),
                                                                          mydata.eval=data.frame(z=zeval),
                                                                          bws=bw.resid.fitted.w.z,
                                                                          p=rep(p, num.z.numeric),
                                                                          ...)
            }

        } else {

            if(return.weights.phi.deriv.1) {

                T.mat.adjoint.deriv.1 <- Kmat.lp(deriv=1,
                                                 mydata.train=data.frame(z),
                                                 bws=bw.resid.fitted.w.z,
                                                 p=rep(p, num.z.numeric),
                                                 ...)

                if(!is.null(zeval)) {

                    T.mat.adjoint.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                          mydata.train=data.frame(z),
                                                          mydata.eval=data.frame(z=zeval),
                                                          bws=bw.resid.fitted.w.z,
                                                          p=rep(p, num.z.numeric),
                                                          ...)

                }

            }

            if(return.weights.phi.deriv.2) {

                T.mat.adjoint.deriv.2 <- Kmat.lp(deriv=2,
                                                 mydata.train=data.frame(z),
                                                 bws=bw.resid.fitted.w.z,
                                                 p=rep(p, num.z.numeric),
                                                 ...)

                if(!is.null(zeval)) {

                    T.mat.adjoint.deriv.eval.2 <- Kmat.lp(deriv=2,
                                                          mydata.train=data.frame(z),
                                                          mydata.eval=data.frame(z=zeval),
                                                          bws=bw.resid.fitted.w.z,
                                                          p=rep(p, num.z.numeric),
                                                          ...)

                }

            }

        }

        if(smooth.residuals) {
            if(!is.null(zeval)) {
                H.eval <- H.eval + constant*T.mat.adjoint.eval%*%(T.mat-T.mat%*%H)
                if(return.weights.phi.deriv.1) H.deriv.eval.1 <- H.deriv.eval.1 + constant*T.mat.adjoint.deriv.eval.1%*%(T.mat-T.mat%*%H)
                if(p>1 && return.weights.phi.deriv.2) H.deriv.eval.2 <- H.deriv.eval.2 + constant*T.mat.adjoint.deriv.eval.2%*%(T.mat-T.mat%*%H)
            }
            if(return.weights.phi.deriv.1) H.deriv.1 <- H.deriv.1 + constant*T.mat.adjoint.deriv.1%*%(T.mat-T.mat%*%H)
            if(p>1  && return.weights.phi.deriv.2) H.deriv.2 <- H.deriv.2 + constant*T.mat.adjoint.deriv.2%*%(T.mat-T.mat%*%H)
            H <- H + constant*T.mat.adjoint%*%(T.mat-T.mat%*%H)
        } else {
            if(!is.null(zeval)) {
                H.eval <- H.eval + constant*T.mat.adjoint.eval%*%(T.mat.r-T.mat%*%H)
                if(return.weights.phi.deriv.1) H.deriv.eval.1 <- H.deriv.eval.1 + constant*T.mat.adjoint.deriv.eval.1%*%(T.mat.r-T.mat%*%H)
                if(p>1 && return.weights.phi.deriv.2) H.deriv.eval.2 <- H.deriv.eval.2 + constant*T.mat.adjoint.deriv.eval.2%*%(T.mat.r-T.mat%*%H)
            }
            if(return.weights.phi.deriv.1) H.deriv.1 <- H.deriv.1 + constant*T.mat.adjoint.deriv.1%*%(T.mat.r-T.mat%*%H)
            if(p>1 && return.weights.phi.deriv.2) H.deriv.2 <- H.deriv.2 + constant*T.mat.adjoint.deriv.2%*%(T.mat.r-T.mat%*%H)
            H <- H + constant*T.mat.adjoint%*%(T.mat.r-T.mat%*%H)
        }

        phi.weights.list[[1]] <- H
        if(return.weights.phi.deriv.1) phi.deriv.1.weights.list[[1]] <- H.deriv.1
        if(p>1 && return.weights.phi.deriv.2) phi.deriv.2.weights.list[[1]] <- H.deriv.2
        if(!is.null(zeval)) {
            phi.eval.weights.list[[1]] <- H.eval
            if(return.weights.phi.deriv.1) phi.deriv.eval.1.weights.list[[1]] <- H.deriv.eval.1
            if(p>1 && return.weights.phi.deriv.2) phi.deriv.eval.2.weights.list[[1]] <- H.deriv.eval.2
        }

    }

    if(!is.null(bw)) iterate.max <- bw$norm.index

    ## In what follows we rbind() bandwidths to return and are careful
    ## about which ones are used when fed in, so we use h$bw below
    ## (but above all are named).

    for(j in 2:iterate.max) {

      console <- printClear(console)
      console <- printPop(console)

      if(smooth.residuals) {
          if(is.null(bw)) {
              console <- printPush(paste("Computing bandwidths and E[y-phi(z)|w] for iteration ", j,"...",sep=""),console)
          } else {
              console <- printPush(paste("Computing E[y-phi(z)|w] for iteration ", j," using supplied bandwidths...",sep=""),console)
          }
      } else {
          if(is.null(bw)) {
              console <- printPush(paste("Computing bandwidths and E[phi(z)|w] for iteration ", j,"...",sep=""),console)
          } else {
              console <- printPush(paste("Computing E[phi(z)|w] for iteration ", j," using supplied bandwidths...",sep=""),console)
          }
      }

      if(smooth.residuals) {

        resid <- y - phi

        if(is.null(bw)) {
            h <- glpcv(ydat=resid,
                       xdat=w,
                       degree=rep(p, num.w.numeric),
                       nmulti=nmulti.loop,
                       random.seed=random.seed,
                       optim.maxattempts=optim.maxattempts,
                       optim.method=optim.method,
                       optim.reltol=optim.reltol,
                       optim.abstol=optim.abstol,
                       optim.maxit=optim.maxit,
                       bw.init=bw.resid.w.mat[j-1,],
                       ...)
        } else {
            h <- NULL
            h$bw <- bw$bw.resid.w[j,]
        }

        resid.fitted <- glpreg(tydat=resid,
                               txdat=w,
                               bws=h$bw,
                               degree=rep(p, num.w.numeric),
                               ...)$mean


      } else {

        if(is.null(bw)) {
            h <- glpcv(ydat=phi,
                       xdat=w,
                       degree=rep(p, num.w.numeric),
                       nmulti=nmulti.loop,
                       random.seed=random.seed,
                       optim.maxattempts=optim.maxattempts,
                       optim.method=optim.method,
                       optim.reltol=optim.reltol,
                       optim.abstol=optim.abstol,
                       optim.maxit=optim.maxit,
                       bw.init=bw.resid.w.mat[j-1,],
                       ...)
        } else {
            h <- NULL
            h$bw <- bw$bw.resid.w[j,]
        }

        resid.fitted <- E.y.w - glpreg(tydat=phi,
                                       txdat=w,
                                       bws=h$bw,
                                       degree=rep(p, num.w.numeric),
                                       ...)$mean

      }

      bw.resid.w.mat[j,] <- h$bw

      if(return.weights.phi) {

          T.mat <- Kmat.lp(mydata.train=data.frame(w),
                           bws=h$bw,
                           p=rep(p, num.w.numeric),
                           ...)

      }

      norm.stop[j] <- ifelse(penalize.iteration,j*sum(resid.fitted^2)/NZD_pos(sum(E.y.w^2)),sum(resid.fitted^2)/NZD_pos(sum(E.y.w^2)))

      console <- printClear(console)
      console <- printPop(console)
      if(smooth.residuals) {
          if(is.null(bw)) {
              console <- printPush(paste("Computing bandwidths and E[E(y-phi(z)|w)|z] for iteration ", j,"...",sep=""),console)
          } else {
              console <- printPush(paste("Computing E[E(y-phi(z)|w)|z] for iteration ", j," using supplied bandwidths...",sep=""),console)
          }
      } else {
          if(is.null(bw)) {
              console <- printPush(paste("Computing bandwidths and E[E(y|z)-E(phi(z)|w)|z] for iteration ", j,"...",sep=""),console)
          } else {
              console <- printPush(paste("Computing E[E(y|z)-E(phi(z)|w)|z] for iteration ", j," using supplied bandwidths...",sep=""),console)
          }
      }

      if(is.null(bw)) {
          h <- glpcv(ydat=resid.fitted,
                     xdat=z,
                     degree=rep(p, num.z.numeric),
                     nmulti=nmulti.loop,
                     random.seed=random.seed,
                     optim.maxattempts=optim.maxattempts,
                     optim.method=optim.method,
                     optim.reltol=optim.reltol,
                     optim.abstol=optim.abstol,
                     optim.maxit=optim.maxit,
                     bw.init=bw.resid.fitted.w.z.mat[j-1,],
                     ...)
      } else {
          h$bw <- bw$bw.resid.fitted.w.z[j,]
      }

      g <- glpreg(tydat=resid.fitted,
                  txdat=z,
                  bws=h$bw,
                  degree=rep(p, num.z.numeric),
                  ...)

      phi <- phi + constant*g$mean
      phi.mat[,j] <- phi
      phi.deriv.1 <- phi.deriv.1 + constant*g$grad
      phi.deriv.1.list[[j]] <- phi.deriv.1

      if(p >= 2) {
          phi.deriv.2 <- phi.deriv.2 + constant*glpreg(tydat=resid.fitted,
                                                       txdat=z,
                                                       bws=h$bw,
                                                       degree=rep(p, num.z.numeric),
                                                       deriv=2,
                                                       ...)$grad
          phi.deriv.2.list[[j]] <- phi.deriv.2
      }

      if(!is.null(zeval)) {
          g <- glpreg(tydat=resid.fitted,
                      txdat=z,
                      exdat=zeval,
                      bws=h$bw,
                      degree=rep(p, num.z.numeric),
                      ...)
          phi.eval <- phi.eval + constant*g$mean
          phi.eval.mat[,j] <- phi.eval

          phi.deriv.eval.1 <- phi.deriv.eval.1 + constant*g$grad
          phi.deriv.eval.1.list[[j]] <- phi.deriv.eval.1
          if(p >= 2) {
              phi.deriv.eval.2 <- phi.deriv.eval.2 + constant*glpreg(tydat=resid.fitted,
                                                                     txdat=z,
                                                                     exdat=zeval,
                                                                     bws=h$bw,
                                                                     degree=rep(p, num.z.numeric),
                                                                     deriv=2,
                                                                     ...)$grad
              phi.deriv.eval.2.list[[j]] <- phi.deriv.eval.2
          }

      }

      bw.resid.fitted.w.z.mat[j,] <- h$bw

      if(return.weights.phi) {

          T.mat.adjoint <- Kmat.lp(mydata.train=data.frame(z),
                                   bws=h$bw,
                                   p=rep(p, num.z.numeric),
                                   ...)

          if(!is.null(zeval)) T.mat.adjoint.eval <- Kmat.lp(mydata.train=data.frame(z),
                                                            mydata.eval=data.frame(z=zeval),
                                                            bws=h$bw,
                                                            p=rep(p, num.z.numeric),
                                                            ...)
          if(p==0 || p==1) {

              if(return.weights.phi.deriv.1) {

                  T.mat.adjoint.deriv.1 <- Kmat.lp(deriv=1,
                                                   mydata.train=data.frame(z),
                                                   bws=h$bw,
                                                   p=rep(p, num.z.numeric),
                                                   ...)

                  if(!is.null(zeval)) T.mat.adjoint.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                                            mydata.train=data.frame(z),
                                                                            mydata.eval=data.frame(z=zeval),
                                                                            bws=h$bw,
                                                                            p=rep(p, num.z.numeric),
                                                                            ...)
              }

          } else {

              if(return.weights.phi.deriv.1) {

                  T.mat.adjoint.deriv.1 <- Kmat.lp(deriv=1,
                                                   mydata.train=data.frame(z),
                                                   bws=h$bw,
                                                   p=rep(p, num.z.numeric),
                                                   ...)

                  if(!is.null(zeval)) {

                      T.mat.adjoint.deriv.eval.1 <- Kmat.lp(deriv=1,
                                                            mydata.train=data.frame(z),
                                                            mydata.eval=data.frame(z=zeval),
                                                            bws=h$bw,
                                                            p=rep(p, num.z.numeric),
                                                            ...)

                  }

              }

              if(return.weights.phi.deriv.2) {

                  T.mat.adjoint.deriv.2 <- Kmat.lp(deriv=2,
                                                   mydata.train=data.frame(z),
                                                   bws=h$bw,
                                                   p=rep(p, num.z.numeric),
                                                   ...)

                  if(!is.null(zeval)) {

                      T.mat.adjoint.deriv.eval.2 <- Kmat.lp(deriv=2,
                                                            mydata.train=data.frame(z),
                                                            mydata.eval=data.frame(z=zeval),
                                                            bws=h$bw,
                                                            p=rep(p, num.z.numeric),
                                                            ...)

                  }

              }

          }

          if(smooth.residuals) {
              if(!is.null(zeval)) {
                  H.eval <- H.eval + constant*T.mat.adjoint.eval%*%(T.mat-T.mat%*%H)
                  if(return.weights.phi.deriv.1) H.deriv.eval.1 <- H.deriv.eval.1 + constant*T.mat.adjoint.deriv.eval.1%*%(T.mat-T.mat%*%H)
                  if(p>1 && return.weights.phi.deriv.2) H.deriv.eval.2 <- H.deriv.eval.2 + constant*T.mat.adjoint.deriv.eval.2%*%(T.mat-T.mat%*%H)
              }
              if(return.weights.phi.deriv.1) H.deriv.1 <- H.deriv.1 + constant*T.mat.adjoint.deriv.1%*%(T.mat-T.mat%*%H)
              if(p>1  && return.weights.phi.deriv.2) H.deriv.2 <- H.deriv.2 + constant*T.mat.adjoint.deriv.2%*%(T.mat-T.mat%*%H)
              H <- H + constant*T.mat.adjoint%*%(T.mat-T.mat%*%H)
          } else {
              if(!is.null(zeval)) {
                  H.eval <- H.eval + constant*T.mat.adjoint.eval%*%(T.mat.r-T.mat%*%H)
                  if(return.weights.phi.deriv.1) H.deriv.eval.1 <- H.deriv.eval.1 + constant*T.mat.adjoint.deriv.eval.1%*%(T.mat.r-T.mat%*%H)
                  if(p>1 && return.weights.phi.deriv.2) H.deriv.eval.2 <- H.deriv.eval.2 + constant*T.mat.adjoint.deriv.eval.2%*%(T.mat.r-T.mat%*%H)
              }
              if(return.weights.phi.deriv.1) H.deriv.1 <- H.deriv.1 + constant*T.mat.adjoint.deriv.1%*%(T.mat.r-T.mat%*%H)
              if(p>1 && return.weights.phi.deriv.2) H.deriv.2 <- H.deriv.2 + constant*T.mat.adjoint.deriv.2%*%(T.mat.r-T.mat%*%H)
              H <- H + constant*T.mat.adjoint%*%(T.mat.r-T.mat%*%H)
          }

          phi.weights.list[[j]] <- H
          if(return.weights.phi.deriv.1) phi.deriv.1.weights.list[[j]] <- H.deriv.1
          if(p>1 && return.weights.phi.deriv.2) phi.deriv.2.weights.list[[j]] <- H.deriv.2
          if(!is.null(zeval)) {
              phi.eval.weights.list[[j]] <- H.eval
              if(return.weights.phi.deriv.1) phi.deriv.eval.1.weights.list[[j]] <- H.deriv.eval.1
              if(p>1 && return.weights.phi.deriv.2) phi.deriv.eval.2.weights.list[[j]] <- H.deriv.eval.2
          }

      }

      console <- printClear(console)
      console <- printPop(console)
      if(is.null(bw)) console <- printPush(paste("Computing stopping rule for iteration ", j,"...",sep=""),console)

      ## The number of iterations in LF is asymptotically equivalent
      ## to 1/alpha (where alpha is the regularization parameter in
      ## Tikhonov).  Plus the criterion function we use is increasing
      ## for very small number of iterations. So we need a threshold
      ## after which we can pretty much confidently say that the
      ## stopping criterion is decreasing.  In Darolles et al. (2011)
      ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so
      ## called qualification of your regularization method. Take the
      ## worst case in which beta = 0 and then the number of
      ## iterations is ~ N^0.5.

      if(is.null(bw))  {

          if(j > round(sqrt(nrow(z))) && !is.monotone.increasing(norm.stop)) {

              ## If stopping rule criterion increases or we are below stopping
              ## tolerance then break

              if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
                  convergence <- "STOP_ON_INCREASE"
                  break()
              }
              if(abs(norm.stop[j-1]-norm.stop[j]) < iterate.diff.tol) {
                  convergence <- "ITERATE_DIFF_TOL"
                  break()
              }

          }

          convergence <- "ITERATE_MAX"

      }

    }

    ## Trim matrices and lists to the actual number of iterations performed
    if(j < iterate.max) {
        phi.mat <- phi.mat[, 1:j, drop=FALSE]
        if(!is.null(phi.eval.mat)) phi.eval.mat <- phi.eval.mat[, 1:j, drop=FALSE]
        norm.stop <- norm.stop[1:j]
        phi.deriv.1.list <- phi.deriv.1.list[1:j]
        phi.deriv.2.list <- phi.deriv.2.list[1:j]
        if(!is.null(phi.deriv.eval.1.list)) phi.deriv.eval.1.list <- phi.deriv.eval.1.list[1:j]
        if(!is.null(phi.deriv.eval.2.list)) phi.deriv.eval.2.list <- phi.deriv.eval.2.list[1:j]
        phi.weights.list <- phi.weights.list[1:j]
        phi.deriv.1.weights.list <- phi.deriv.1.weights.list[1:j]
        phi.deriv.2.weights.list <- phi.deriv.2.weights.list[1:j]
        phi.eval.weights.list <- phi.eval.weights.list[1:j]
        phi.deriv.eval.1.weights.list <- phi.deriv.eval.1.weights.list[1:j]
        phi.deriv.eval.2.weights.list <- phi.deriv.eval.2.weights.list[1:j]
        bw.resid.w.mat <- bw.resid.w.mat[1:j, , drop=FALSE]
        bw.resid.fitted.w.z.mat <- bw.resid.fitted.w.z.mat[1:j, , drop=FALSE]
    }

    bw.resid.w <- bw.resid.w.mat
    bw.resid.fitted.w.z <- bw.resid.fitted.w.z.mat

    ## Extract minimum, and check for monotone increasing function and
    ## then decreasing (and potentially increasing thereafter) portion
    ## of the stopping function, ignore the initial increasing portion,
    ## and take the min from where the initial inflection point occurs
    ## to the length of norm.stop

    phi.weights <- NULL
    phi.deriv.1.weights <- NULL
    phi.deriv.2.weights <- NULL
    
    phi.eval.weights <- NULL
    phi.deriv.eval.1.weights <- NULL
    phi.deriv.eval.2.weights <- NULL

    if(is.null(bw))  {

      norm.value <- norm.stop/(1:length(norm.stop))

      if(which.min(norm.stop) == 1 && is.monotone.increasing(norm.stop)) {
          warning("Stopping rule increases monotonically (consult model$norm.stop):\nThis could be the result of an inspired initial value (unlikely)\nNote: we suggest manually choosing phi.0 and restarting (e.g. instead set `starting.values' to E[E(Y|w)|z])")
          convergence <- "FAILURE_MONOTONE_INCREASING"
          phi <- starting.values.phi
          j <- 1
          while(norm.value[j+1] > norm.value[j]) j <- j + 1
          j <- j-1 + which.min(norm.value[j:length(norm.value)])
          phi <- phi.mat[,j]
          if(p>0) phi.deriv.1 <- phi.deriv.1.list[[j]]
          if(p>=2) phi.deriv.2 <- phi.deriv.2.list[[j]]
          if(return.weights.phi) phi.weights <- phi.weights.list[[j]]
          if(return.weights.phi.deriv.1) phi.deriv.1.weights <- phi.deriv.1.weights.list[[j]]
          if(p>=2 && return.weights.phi.deriv.2) phi.deriv.2.weights <- phi.deriv.2.weights.list[[j]]
          if(!is.null(zeval)) {
              phi.eval <- phi.eval.mat[,j]
              if(p>0) phi.deriv.eval.1 <- phi.deriv.eval.1.list[[j]]
              if(p>=2) phi.deriv.eval.2 <- phi.deriv.eval.2.list[[j]]
              if(return.weights.phi) phi.eval.weights <- phi.eval.weights.list[[j]]
              if(return.weights.phi.deriv.1) phi.deriv.eval.1.weights <- phi.deriv.eval.1.weights.list[[j]]
              if(p>=2 && return.weights.phi.deriv.2) phi.deriv.eval.2.weights <- phi.deriv.eval.2.weights.list[[j]]    
          }
      } else {
          ## Ignore the initial increasing portion, take the min to the
          ## right of where the initial inflection point occurs
          j <- 1
          while(norm.stop[j+1] > norm.stop[j]) j <- j + 1
          j <- j-1 + which.min(norm.stop[j:length(norm.stop)])
          phi <- phi.mat[,j]
          if(p>0) phi.deriv.1 <- phi.deriv.1.list[[j]]
          if(p>=2) phi.deriv.2 <- phi.deriv.2.list[[j]]
          if(return.weights.phi) phi.weights <- phi.weights.list[[j]]
          if(return.weights.phi.deriv.1) phi.deriv.1.weights <- phi.deriv.1.weights.list[[j]]
          if(p>=2 && return.weights.phi.deriv.2) phi.deriv.2.weights <- phi.deriv.2.weights.list[[j]]
          if(!is.null(zeval)) {
              phi.eval <- phi.eval.mat[,j]
              if(p>0) phi.deriv.eval.1 <- phi.deriv.eval.1.list[[j]]
              if(p>=2) phi.deriv.eval.2 <- phi.deriv.eval.2.list[[j]]
              if(return.weights.phi) phi.eval.weights <- phi.eval.weights.list[[j]]
              if(return.weights.phi.deriv.1) phi.deriv.eval.1.weights <- phi.deriv.eval.1.weights.list[[j]]
              if(p>=2 && return.weights.phi.deriv.2) phi.deriv.eval.2.weights <- phi.deriv.eval.2.weights.list[[j]]    
          }
      }
      if(j == iterate.max) warning("iterate.max reached: increase iterate.max or inspect norm.stop vector")      
    } else {
        ## bw passed in, set j to norm.index, push out weights etc.
        j <- bw$norm.index
        phi <- phi.mat[,j]
        if(p>0) phi.deriv.1 <- phi.deriv.1.list[[j]]
        if(p>=2) phi.deriv.2 <- phi.deriv.2.list[[j]]
        if(return.weights.phi) phi.weights <- phi.weights.list[[j]]
        if(return.weights.phi.deriv.1) phi.deriv.1.weights <- phi.deriv.1.weights.list[[j]]
        if(p>=2 && return.weights.phi.deriv.2) phi.deriv.2.weights <- phi.deriv.2.weights.list[[j]]
        if(!is.null(zeval)) {
            phi.eval <- phi.eval.mat[,j]
            if(p>0) phi.deriv.eval.1 <- phi.deriv.eval.1.list[[j]]
            if(p>=2) phi.deriv.eval.2 <- phi.deriv.eval.2.list[[j]]
            if(return.weights.phi) phi.eval.weights <- phi.eval.weights.list[[j]]
            if(return.weights.phi.deriv.1) phi.deriv.eval.1.weights <- phi.deriv.eval.1.weights.list[[j]]
            if(p>=2 && return.weights.phi.deriv.2) phi.deriv.eval.2.weights <- phi.deriv.eval.2.weights.list[[j]]    
        }
        norm.value <- NULL
        norm.stop <- NULL
        convergence <- NULL
    }
    
    console <- printClear(console)
    console <- printPop(console)

    ret <- list(phi=phi,
                phi.mat=phi.mat,
                phi.deriv.1=as.matrix(phi.deriv.1),
                phi.deriv.2=if(!is.null(phi.deriv.2)){as.matrix(phi.deriv.2)}else{NULL},
                phi.weights=phi.weights,
                phi.deriv.1.weights=phi.deriv.1.weights,
                phi.deriv.2.weights=phi.deriv.2.weights,
                phi.eval=phi.eval,
                phi.eval.mat=phi.eval.mat,
                phi.deriv.eval.1=if(!is.null(phi.deriv.eval.1)){as.matrix(phi.deriv.eval.1)}else{NULL},
                phi.deriv.eval.2=if(!is.null(phi.deriv.eval.2)){as.matrix(phi.deriv.eval.2)}else{NULL},
                phi.eval.weights=phi.eval.weights,
                phi.deriv.eval.1.weights=phi.deriv.eval.1.weights,
                phi.deriv.eval.2.weights=phi.deriv.eval.2.weights,
                norm.index=j,
                norm.stop=norm.stop,
                norm.value=norm.value,
                convergence=convergence,
                starting.values.phi=starting.values.phi,
                return.weights.phi=return.weights.phi,
                bw.E.y.w=bw.E.y.w,
                bw.E.y.z=bw.E.y.z,
                bw.resid.w=as.matrix(bw.resid.w),
                bw.resid.fitted.w.z=as.matrix(bw.resid.fitted.w.z),
                call=cl,
                y=y,
                z=z,
                w=w,
                x=x,
                zeval=zeval,
                xeval=xeval,
                p=p,
                nmulti=nmulti,
                method=method,
                ptm=proc.time() - ptm.start)
    class(ret) <- "npregiv"
    return(ret)

  }

}

print.npregiv <- function(x, ...) {
  summary.npregiv(x, ...)
  invisible(x)
}

summary.npregiv <- function(object, ...) {
  format_bw <- function(bw, label, names = NULL) {
    if(is.null(bw)) return()
    if(is.matrix(bw)) bw <- bw[nrow(bw), , drop = TRUE]
    bw <- as.numeric(bw)
    if(!is.null(names) && length(names) == length(bw)) {
      vals <- paste(paste(names, formatC(bw, digits=8, format="g"), sep=": "), collapse=", ")
    } else {
      vals <- paste(formatC(bw, digits=8, format="g"), collapse=", ")
    }
    cat(paste("\n", label, " ", vals, sep=""))
  }

  cat("Call:\n")
  print(object$call)

  if(is.null(object$alpha))
    cat("\nNonparametric Instrumental Kernel Regression\n",sep="")
  else
    cat("\nNonparametric Instrumental Kernel Regression (Tikhonov)\n",sep="")

  cat(paste("\nNumber of continuous endogenous predictors: ",format(NCOL(object$z)),sep=""),sep="")
  cat(paste("\nNumber of continuous instruments: ",format(NCOL(object$w)),sep=""),sep="")
  if(!is.null(object$x)) cat(paste("\nNumber of continuous exogenous predictors: ",format(NCOL(object$x)),sep=""),sep="")

  cat(paste("\nLocal polynomial order (p): ", format(object$p), sep=""))
  cat(paste("\nTraining observations: ", format(NROW(object$y)), sep=""))

  if(!is.null(object$alpha)) {
    cat(paste("\n\nRegularization method: Tikhonov",sep=""))
    cat(paste("\nTikhonov parameter (alpha): ", format(object$alpha,digits=8), sep=""))
    if(!is.null(object$alpha.iter)) cat(paste("\nIterated Tikhonov parameter (alpha.iter): ", format(object$alpha.iter,digits=8), sep=""))
  } else {
    cat(paste("\n\nRegularization method: Landweber-Fridman",sep=""))
    cat(paste("\nNumber of iterations: ", format(object$norm.index), sep=""))
    cat(paste("\nStopping rule value: ", format(object$norm.stop[length(object$norm.stop)],digits=8), sep=""))
  }

  w.names <- if(!is.null(object$w)) colnames(object$w) else NULL
  z.names <- if(!is.null(object$z)) colnames(object$z) else NULL

  if(is.null(object$alpha)) {
    format_bw(object$bw.E.y.w, "Bandwidth for E(y|w):", w.names)
    format_bw(object$bw.E.y.z, "Bandwidth for E(y|z):", z.names)
    format_bw(object$bw.resid.w, "Bandwidth for E(y-phi(z)|w):", w.names)
    format_bw(object$bw.resid.fitted.w.z, "Bandwidth for E(E(y-phi(z)|w)|z):", z.names)
  } else {
    format_bw(object$bw.E.y.w, "Bandwidth for E(y|w):", w.names)
    format_bw(object$bw.E.E.y.w.z, "Bandwidth for E(E(y|w)|z):", z.names)
    format_bw(object$bw.E.phi.w, "Bandwidth for E(phi(z)|w):", w.names)
    format_bw(object$bw.E.E.phi.w.z, "Bandwidth for E(E(phi(z)|w)|z):", z.names)
  }

  cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))
  cat("\n\n")
}

plot.npregiv <- function(x,
                         plot.data = FALSE,
                         deriv = FALSE,
                         ...) {

  object <- x

  ## We only support univariate endogenous predictor z
  if(NCOL(object$z) > 1) stop(" only univariate z is supported")

  z <- object$z[,1]
  y <- object$y
  phi <- object$phi

  zname <- names(object$z)[1]
  yname <- "y" ## Default

  if(deriv) {
    phi.prime <- object$phi.deriv.1[,1]

    plot(z[order(z)], phi.prime[order(z)],
         type="l",
         xlab=zname,
         ylab=paste("d", yname, "/d", zname, sep=""),
         ...)

  } else {

    if(plot.data) {
      plot(z, y,
           xlab=zname,
           ylab=yname,
           type="p",
           col="lightgrey",
           ...)
      lines(z[order(z)], phi[order(z)],
            lwd=2,
            ...)
    } else {
      plot(z[order(z)], phi[order(z)],
           type="l",
           xlab=zname,
           ylab=yname,
           lwd=2,
           ...)
    }
  }

}
