## This functions accepts the following arguments:

## y: univariate outcome
## z: endogenous predictors
## w: instruments
## x: exogenous predictors

## zeval: optional evaluation data for the endogenous predictors
## weval: optional evaluation data for the instruments
## xeval: optional evaluation data for the exogenous predictors

## ... optional arguments for npreg()

## This function returns a list with the following elements:

## phi: the IV estimator of phi(z) corresponding to the estimated
## derivative phihat(z)
## phi.prime: the IV derivative estimator
## phi.mat: the matrix with colums phi_1, phi_2 etc. over all iterations
## phi.prime.mat: the matrix with colums phi'_1, phi'_2 etc. over all iterations
## num.iterations: number of iterations taken by Landweber-Fridman
## norm.stop: the stopping rule for each Landweber-Fridman iteration
## convergence: a character string indicating whether/why iteration terminated

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

## This function will compute the cumulative integral at each sample
## realization using the trapezoidal rule and the cumsum function as
## we need to compute this in a computationally efficient manner.

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
                         random.seed=42,
                         iterate.max=1000,
                         iterate.break=TRUE,
                         constant=0.5,
                         start.from=c("Eyz","EEywz"),
                         starting.values=NULL,
                         stop.on.increase=TRUE,
                         smooth.residuals=TRUE,
                         ...) {

  console <- newLineConsole()

  ## Basic error checking

  if(!is.logical(stop.on.increase)) stop("stop.on.increase must be logical (TRUE/FALSE)")
  if(!is.logical(smooth.residuals)) stop("smooth.residuals must be logical (TRUE/FALSE)")

  start.from <- match.arg(start.from)

  if(iterate.max < 2) stop("iterate.max must be at least 2")
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

  if(!is.null(starting.values) && (NROW(starting.values) != NROW(zeval))) stop(paste("starting.values must be of length",NROW(zeval)))

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
  ## considerations, same smoothing for PDF and CDF.

  bw <- npudensbw(dat=z, bwmethod="normal-reference")
  model.fz <- npudens(tdat=z, bws=bw$bw)
  f.z <- predict(model.fz, newdata=zeval)
  model.Sz <- npudist(tdat=z, bws=bw$bw)
  S.z <- 1-predict(model.Sz, newdata=zeval)

  console <- printClear(console)
  console <- printPop(console)
  console <- printPush(paste("Computing optimal smoothing for E(y|w) (stopping rule) for iteration 1...",sep=""),console)

  ## For stopping rule...

  E.y.w <- npreg(tydat=y,
                 txdat=w,
                 exdat=weval,
                 ...)$mean

  ## Potential alternative starting rule (consistent with
  ## npregiv). Here we start with E(Y|Z) rather than zero

  if(is.null(starting.values)) {

    console <- printClear(console)
    console <- printPop(console)
    if(is.null(x)) {
      console <- printPush(paste("Computing optimal smoothing for E(y|z) for iteration 1...",sep=""),console)
    } else {
      console <- printPush(paste("Computing optimal smoothing for E(y|z,x) for iteration 1...",sep=""),console)
    }

    phi.prime <- npreg(tydat=if(start.from=="Eyz") y else E.y.w,
                       txdat=z,
                       exdat=zeval,
                       gradients=TRUE,
                       ...)$grad[,1]

  } else {

    phi.prime <- starting.values

  }

  ## Now we can compute phi.0 by integrating phi.prime.0 up to each
  ## sample realization (here we use the trapezoidal rule)

  ## NOTE - this presumes univariate z case... in general this would
  ## be a continuous variable's index

  phi <- integrate.trapezoidal(z[,1],phi.prime)

  starting.values.phi <- phi
  starting.values.phi.prime <- phi.prime

  ## In the definition of phi we have the integral minus the mean of
  ## the integral with respect to z, so subtract the mean here

  phi <- phi - mean(phi) + mean(y)

  console <- printClear(console)
  console <- printPop(console)
  console <- printPush(paste("Computing optimal smoothing for E(phi|w) (stopping rule) for iteration 1...",sep=""),console)

  norm.stop <- numeric()

  ## Now we compute mu.0 (a residual of sorts)

  mu <- y - phi

  ## Now we repeat this entire process using mu = y = phi.0 rather than y

  if(smooth.residuals) {

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing for E(mu|w) (stopping rule) for iteration 1...",sep=""),console)

    ## Additional smoothing on top of the stopping rule required, but
    ## we have computed the stopping rule so reuse the bandwidth
    ## vector to be passed below. Here we compute the bandwidth
    ## optimal for the regression of mu on w.

    ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w

    predicted.E.mu.w <- npreg(tydat=mu,
                              txdat=w,
                              exdat=weval,
                              ...)$mean

  } else {

    E.phi.w <- npreg(tydat=phi,
                     txdat=w,
                     exdat=weval,
                     ...)$mean

    predicted.E.mu.w <- E.y.w - E.phi.w

  }


  norm.stop[1] <- sum(predicted.E.mu.w^2)/NZD(sum(E.y.w^2))
  
  ## We again require the mean of the fitted values

  mean.predicted.E.mu.w <- mean(predicted.E.mu.w)

  ## Now we compute T^* applied to mu

  cdf.weighted.average <- npksum(txdat=z,
                                 exdat=zeval,
                                 tydat=as.matrix(predicted.E.mu.w),
                                 operator="integral",
                                 ukertype="liracine",
                                 okertype="liracine",
                                 bws=bw$bw,
                                 ...)$ksum/length(y)

  survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

  T.star.mu <- (survivor.weighted.average-S.z*mean.predicted.E.mu.w)/NZD(f.z)

  ## Now we update phi.prime.0, this provides phi.prime.1, and now
  ## we can iterate until convergence... note we replace phi.prime.0
  ## with phi.prime.1 (i.e. overwrite phi.prime)

  phi.prime <- phi.prime + constant*T.star.mu
  phi.mat <- phi
  phi.prime.mat <- phi.prime

  ## This we iterate...

  for(j in 2:iterate.max) {

    console <- printClear(console)
    console <- printPop(console)
    console <- printPush(paste("Computing optimal smoothing for E(phi|w) for iteration ", j,"...",sep=""),console)

    ## NOTE - this presumes univariate z case... in general this would
    ## be a continuous variable's index

    phi <- integrate.trapezoidal(z[,1],phi.prime)

    ## In the definition of phi we have the integral minus the mean of
    ## the integral with respect to z, so subtract the mean here

    phi <- phi - mean(phi) + mean(y)

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

      ## Next, we regress require \mu_{0,i} W using bws optimal for phi on w

      predicted.E.mu.w <- npreg(tydat=mu,
                                txdat=w,
                                eydat=mu,
                                exdat=weval,
                                ...)$mean

    } else {

      E.phi.w <- npreg(tydat=phi,
                       txdat=w,
                       eydat=phi,
                       exdat=weval,
                       ...)$mean

      predicted.E.mu.w <- E.y.w - E.phi.w

    }

    norm.stop[j] <- j*sum(predicted.E.mu.w^2)/NZD(sum(E.y.w^2))
    
    mean.predicted.E.mu.w <- mean(predicted.E.mu.w)

    ## Now we compute T^* applied to mu

    cdf.weighted.average <- npksum(txdat=z,
                                   exdat=zeval,
                                   tydat=as.matrix(predicted.E.mu.w),
                                   operator="integral",
                                   ukertype="liracine",
                                   okertype="liracine",
                                   bws=bw$bw,
                                   ...)$ksum/length(y)

    survivor.weighted.average <- mean.predicted.E.mu.w - cdf.weighted.average

    T.star.mu <- (survivor.weighted.average-S.z*mean.mu)/NZD(f.z)

    ## Now we update, this provides phi.prime.1, and now we can
    ## iterate until convergence...

    phi.prime <- phi.prime + constant*T.star.mu
    phi.mat <- cbind(phi.mat,phi)
    phi.prime.mat <- cbind(phi.prime.mat,phi.prime)

    ## The number of iterations in LF is asymptotically equivalent to
    ## 1/alpha (where alpha is the regularization parameter in
    ## Tikhonov).  Plus the criterion function we use is increasing
    ## for very small number of iterations. So we need a threshold
    ## after which we can pretty much confidently say that the
    ## stopping criterion is decreasing.  In Darolles et al. (2011)
    ## \alpha ~ O(N^(-1/(min(beta,2)+2)), where beta is the so called
    ## qualification of your regularization method. Take the worst
    ## case in which beta = 0 and then the number of iterations is ~
    ## N^0.5. Note that derivative estimation seems to require more
    ## iterations hence the heuristic sqrt(N)

    if(j > round(sqrt(nrow(z))) && !is.monotone.increasing(norm.stop)) {

      ## If stopping rule criterion increases or we are below stopping
      ## tolerance then break

      if(stop.on.increase && norm.stop[j] > norm.stop[j-1]) {
        convergence <- "STOP_ON_INCREASE"
        if(iterate.break) break()
      }

    }

    convergence <- "ITERATE_MAX"

  }

  ## Extract minimum, and check for monotone increasing function and
  ## issue warning in that case. Otherwise allow for an increasing
  ## then decreasing (and potentially increasing thereafter) portion
  ## of the stopping function, ignore the initial increasing portion,
  ## and take the min from where the initial inflection point occurs
  ## to the length of norm.stop.

  if(is.monotone.increasing(norm.stop)) {
    warning("Stopping rule increases monotonically (consult model$norm.stop):\nThis could be the result of an inspired initial value (unlikely)\nNote: we suggest manually choosing phi.0 and restarting (e.g., instead set `start.from' to EEywz or provide a vector of starting values")
    convergence <- "FAILURE_MONOTONE_INCREASING"
    j <- length(norm.stop)
    phi <- phi.mat[,1]
    phi.prime <- phi.prime.mat[,1]
  } else {
    ## Ignore the initial increasing portion, take the min to the
    ## right of where the initial inflection point occurs.
    j <- 1
    ## Climb the initial hill...
    while(norm.stop[j+1] >= norm.stop[j] & j < length(norm.stop)) j <- j + 1
    ## Descend into the first valley
    while(norm.stop[j+1] < norm.stop[j] & j < length(norm.stop)) j <- j + 1
    ## When you start to climb again, stop, previous location was min
    phi <- phi.mat[,j-1]
    phi.prime <- phi.prime.mat[,j-1]
  }
  
  console <- printClear(console)
  console <- printPop(console)

  if(j == iterate.max) warning(" iterate.max reached: increase iterate.max or inspect norm.stop vector")

  return(list(phi=phi,
              phi.prime=phi.prime,
              phi.mat=phi.mat,
              phi.prime.mat=phi.prime.mat,
              num.iterations=j,
              norm.stop=norm.stop,
              convergence=convergence,
              starting.values.phi=starting.values.phi,
              starting.values.phi.prime=starting.values.phi.prime))

}

