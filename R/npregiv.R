## $Id: npregiv.R,v 1.8 2011/05/30 19:17:35 jracine Exp jracine $

## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instrument
## xeval: optional evaluation data for the exogenous predictors

## alpha.min: minimum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## alpha.max: maximum value when conducting 1-dimensional search for
##            optimal Tikhonov regularization parameter alpha

## p: order of the local polynomial kernel estimator (p=0 is local
##    constant, p=1 local linear etc.)

## This function returns a list with the following elements:

## phihat: the IV estimator of phi(y)
## alpha:  the Tikhonov regularization parameter

npregiv <- function(y,
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
                    alpha=NULL,
                    alpha.min=1.0e-10,
                    alpha.max=1.0e-01,
                    alpha.tol=.Machine$double.eps^0.25,
                    iterate.max=100,
                    iterate.tol=1.0e-04,
                    constant=0.5,
                    method=c("Landweber-Fridman","Tikhonov"),
                    stop.on.increase=TRUE,
                    ...) {

  ## This function was constructed initially by Samuele Centorrino
  ## <samuele.centorrino@univ-tlse1.fr> to reproduce illustrations in
  ## the following papers:
  
  ## A) Econometrica (forthcoming, article date February 25 2011)
  
  ## "Nonparametric Instrumental Regression"
  ## S. Darolles, Y. Fan, J.P. Florens, E. Renault
  
  ## B) Econometrics Journal (2010), volume 13, pp. S1–S27. doi:
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
  
  tikh <- function(alpha,CZ,CY,Cr.r){
    return(solve(alpha*diag(length(Cr.r)) + CY%*%CZ) %*% Cr.r)
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
  
  ittik <- function(alpha,CZ,CY,Cr.r,r) {
    invmat <- solve(alpha*diag(length(Cr.r)) + CY%*%CZ)
    phi <- invmat %*% Cr.r + alpha * invmat %*% invmat %*% Cr.r        
    return(sum((CZ%*%phi - r)^2)/alpha)    
  }

  ## $Id: npregiv.R,v 1.8 2011/05/30 19:17:35 jracine Exp jracine $

  ## This function returns the weight matrix for a local polynomial

  ## supports mixed data types. It presumes that Y is in column 1. Basic
  ## error checking is undertaken. j.reg= strips off weights for mean
  ## (1), partials up to order p, and cross-partials. All partials and
  ## cross partials are wrt continuous regressors, and cross-partials
  ## require k > 1 and p > 1. Shrinking towards the local constant mean,
  ## first, and second partial derivatives is implemented for regions
  ## where the local polynomial estimator is ill-conditioned (sparse
  ## data, small h etc.).
  
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
    
  ## $Id: npregiv.R,v 1.8 2011/05/30 19:17:35 jracine Exp jracine $
  
  ## Functions for generalized local polynomial regression
  
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
                    nmulti=nmulti,
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

  ## Here is where the function `npregiv' really begins:

  console <- newLineConsole()
    
  ## Basic error checking

  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")  
  
  if(missing(y)) stop("You must provide y")
  if(missing(z)) stop("You must provide z")
  if(missing(w)) stop("You must provide w")
  if(NCOL(y) > 1) stop("y must be univariate")
  if(NROW(y) != NROW(z) || NROW(y) != NROW(w)) stop("y, z, and w have differing numbers of rows")
  if(iterate.max < 2) stop("iterate.max must be at least 2")
  if(p < 0) stop("p must be a non-negative integer")

  if(!is.null(alpha) && alpha <= 0) stop("alpha must be positive")

  method <- match.arg(method)
  
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

  if(method=="Tikhonov") {
  
    ## Now y=phi(z) + u, hence E(y|w)=E(phi(z)|w) so we need two
    ## bandwidths, one for y on w and one for phi(z) on w (in the
    ## first step we use E(y|w) as a proxy for phi(z) and use
    ## bandwidths for y on w).
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing bandwidths for E(y|w)...", console)

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

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix and E(y|w)...", console)  

    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    exdat=weval,
                    bws=hyw$bw,
                    degree=rep(p, num.w.numeric),
                    ...)$mean

    KYW <- Kmat.lp(mydata.train=data.frame(w), mydata.eval=data.frame(w=weval), bws=hyw$bw, p=rep(p, num.w.numeric))
    
    ## We conduct local polynomial kernel regression of E(y|w) on z
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing bandwidths for E(E(y|w)|z)...", console)

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

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix and E(E(y|w)|z)...", console)  

    E.E.y.w.z <- glpreg(tydat=E.y.w,
                        txdat=z,
                        eydat=E.y.w,
                        exdat=zeval,
                        bws=hywz$bw,
                        degree=rep(p, num.z.numeric),
                        ...)$mean

    KYWZ <- Kmat.lp(mydata.train=data.frame(z), mydata.eval=data.frame(z=zeval), bws=hywz$bw, p=rep(p, num.z.numeric))
    
    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov function) to
    ## determine the optimal alpha for the non-iterated scheme. Note
    ## that the function `optimize' accepts bounds on the search (in
    ## this case alpha.min to alpha.max))
    
    ## E(r|z)=E(E(phi(z)|w)|z)
    ## \phi^\alpha = (\alpha I+CzCw)^{-1}Cr x r

    if(is.null(alpha)) {
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
    phihat <- as.vector(tikh(alpha, CZ = KYW, CY = KYWZ, Cr.r = E.E.y.w.z))
    
    ## KYWZ and KYWS no longer used, save memory
    
    rm(KYWZ, KYW)
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing bandwidths for E(phi(z)|w)...", console)

    hphiw <- glpcv(ydat=phihat,
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

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix for E(phi(z)|w)...", console)

    E.phihat.w <- glpreg(tydat=phihat,
                         txdat=w,
                         eydat=phihat,
                         exdat=weval,
                         bws=hphiw$bw,
                         degree=rep(p, num.w.numeric),
                         ...)$mean

    KPHIW <- Kmat.lp(mydata.train=data.frame(w), mydata.eval=data.frame(w=weval), bws=hphiw$bw, p=rep(p, num.w.numeric))
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing bandwidths for E(E(phi(z)|w)|z)...", console)

    hphiwz <- glpcv(ydat=E.phihat.w,
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

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing weight matrix for E(E(phi(z)|w)|z)...", console)

    KPHIWZ <- Kmat.lp(mydata.train=data.frame(z), mydata.eval=data.frame(z=zeval), bws=hphiwz$bw, p=rep(p, num.z.numeric))
    
    ## Next, we minimize the function ittik to obtain the optimal value
    ## of alpha (here we use the iterated Tikhonov approach) to
    ## determine the optimal alpha for the non-iterated scheme.
    
    if(is.null(alpha)) {
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush("Iterating and recomputing the numerical solution for alpha...", console)
      alpha <- optimize(ittik, c(alpha.min, alpha.max), tol = alpha.tol, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z, r = E.y.w)$minimum
    }
    
    ## Finally, we conduct regularized Tikhonov regression using this
    ## optimal alpha and the updated bandwidths.
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush("Computing final phi(z) estimate...", console)
    phihat <- as.vector(tikh(alpha, CZ = KPHIW, CY = KPHIWZ, Cr.r = E.E.y.w.z))
    
    console <- printClear(console)
    console <- printPop(console)

    if((alpha-alpha.min)/alpha.min < 0.01) warning(paste("Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search minimum (",alpha.min,")",sep=""))
    if((alpha.max-alpha)/alpha.max < 0.01) warning(paste("Tikhonov parameter alpha (",formatC(alpha,digits=4,format="f"),") is close to the search maximum (",alpha.max,")",sep=""))
    
    return(list(phihat=phihat, alpha=alpha))
    
  } else {

    ## Landweber-Fridman

    ## We begin the iteration computing phi.0 and phi.1 directly, then
    ## interate.
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing bandwidths and E(y|z) for iteration 0...",sep=""),console)

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

    phi.0 <- glpreg(tydat=y,
                    txdat=z,
                    exdat=zeval,
                    bws=h$bw,
                    degree=rep(p, num.z.numeric),
                    ...)$mean
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing bandwidths and E(y-phi(z)|w) for iteration 1...",sep=""),console)

    resid <- y - phi.0

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

    resid.fitted <- glpreg(tydat=resid,
                           txdat=w,
                           eydat=resid,
                           exdat=weval,
                           bws=h$bw,
                           degree=rep(p, num.w.numeric),
                           ...)$mean

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing bandwidths and E(E(y-phi(z)|w)|z) for iteration 1...",sep=""),console)

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

    phi.j.m.1 <- phi.0 + constant*glpreg(tydat=resid.fitted,
                                         txdat=z,
                                         eydat=resid.fitted,
                                         exdat=zeval,
                                         bws=h$bw,
                                         degree=rep(p, num.z.numeric),
                                         ...)$mean
    
    ## For the stopping rule
    
    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing bandwidths and E(y|w) for stopping rule...",sep=""),console)

    norm.stop <- numeric()

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

    E.y.w <- glpreg(tydat=y,
                    txdat=w,
                    exdat=weval,
                    bws=h$bw,
                    degree=rep(p, num.w.numeric),
                    ...)$mean
    phihat <- phi.j.m.1

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing bandwidths and E(phi(z)|w) for stopping rule...",sep=""),console)

    h.E.phi.w <- glpcv(ydat=phi.j.m.1,
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

    E.phi.w <- glpreg(tydat=phi.j.m.1,
                      txdat=w,
                      eydat=phi.j.m.1,
                      exdat=weval,
                      bws=h.E.phi.w$bw,
                      degree=rep(p, num.w.numeric),
                      ...)$mean

    norm.stop[1] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

    for(j in 2:iterate.max) {

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing bandwidths and E(y-phi(z)|w) for iteration ", j,"...",sep=""),console)

      resid <- y - phi.j.m.1
      h.resid.w <- glpcv(ydat=resid,
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

      resid.fitted <- glpreg(tydat=resid,
                             txdat=w,
                             eydat=resid,
                             exdat=weval,
                             bws=h.resid.w$bw,
                             degree=rep(p, num.w.numeric),
                             ...)$mean

      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing bandwidths and E(E(y-phi(z)|w)|z) for iteration ", j,"...",sep=""),console)

      h.resid.fitted.z <- glpcv(ydat=resid.fitted,
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

      phi.j <- phi.j.m.1 + constant*glpreg(tydat=resid.fitted,
                                           txdat=z,
                                           eydat=resid.fitted,
                                           exdat=zeval,
                                           bws=h.resid.fitted.z$bw,
                                           degree=rep(p, num.z.numeric),
                                           ...)$mean
      console <- printClear(console)
      console <- printPop(console)
      console <- printPush(paste("Computing stopping rule for iteration ", j,"...",sep=""),console)

      ## For the stopping rule (use same smoothing as original)

      E.phi.w <- glpreg(tydat=phi.j,
                        txdat=w,
                        eydat=phi.j,
                        exdat=weval,
                        bws=h.E.phi.w$bw,
                        degree=rep(p, num.w.numeric),
                        ...)$mean

      norm.stop[j] <- mean(((E.y.w-E.phi.w)/E.y.w)^2)

      ## If stopping rule criterion increases or we are below stopping
      ## tolerance then break

      if(norm.stop[j] < iterate.tol) break()
      if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
        phi.j <- phi.j.m.1 
        break()
      }
      
      phi.j.m.1 <- phi.j

    }

    console <- printClear(console)
    console <- printPop(console)

    if(j == iterate.max) warning("iterate.max reached: increase iterate.max or inspect norm.stop vector")

    return(list(phihat=phi.j, num.iterations=j, norm.stop=norm.stop))

  }
  
}
