## Function that implements the entropy metric test for serial
## dependence described in "A Dependence Metric For Possibly Nonlinear
## Processes By C. W. Granger, E. Maasoumi and J. Racine Journal of
## Time Series Analysis (2004), Vol. 25, No. 5, 649-669.

.np_sdeptest_fixed_density <- function(data, bw) {
  as.numeric(npksum(
    txdat = data,
    bws = bw,
    bandwidth.divide = TRUE
  )$ksum / NROW(data))
}

npsdeptest <- function(data = NULL,
                       lag.num = 1,
                       method=c("integration","summation"),
                       bootstrap = TRUE,
                       boot.num = 399,
                       random.seed = 42) {
  
  ## Trap fatal errors

  if(is.data.frame(data)) stop(" you must enter a data vector (not data frame)")
  if(is.null(data)) stop(" you must enter a data vector")
  if((lag.num < 1) || (lag.num > length(data))) stop(" lag.num must be a positive integer less than the number of observations")
  if(ncol(data.frame(data)) != 1) stop(" data must have one dimension only")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  method <- match.arg(method)

  ## Save seed prior to setting

  seed.state <- .np_seed_enter(random.seed)
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)
  .np_progress_note("Computing bandwidths")


  ## If the variable is a time series convert to type numeric

  if(is.ts(data)) data <- as.numeric(data)

  ## Omit NAs

  data <- na.omit(data)

  ## Define the metric entropy function

  Srho.bivar <- function(x.dat = NULL,
                         y.dat = NULL,
                         bw.x = NULL,
                         bw.y = NULL,
                         bw.joint = NULL,
                         method=c("integration","summation")) {
    
    if(is.null(bw.x)||is.null(bw.y)||is.null(bw.joint)) stop(" you must provide numeric bandwidths for f(x), f(y) and f(x,y)")
    
    method <- match.arg(method)
    
    if(method=="summation") {
      
      ## Summation version:
      ## \sum_i(1-sqrt(f(x_i)f(y_i)/f(x_i,y_i)))^2
      ## We could code this by hand per below, however, there is no
      ## performance penalty imposed by sum() so no need.

      ## Issue of common support points for evaluation
      
      f.x <- .np_sdeptest_fixed_density(x.dat, bw.x)
      f.y <- .np_sdeptest_fixed_density(y.dat, bw.y)
      f.xy <- .np_sdeptest_fixed_density(cbind(x.dat, y.dat), bw.joint)
      summand <- f.x*f.y/f.xy

      ## In summation version we divide densities which can lead to
      ## numerical instability issues not present in the integrand
      ## version (which uses differences instead). We check for this
      ## case, remove offending points, and warn. This traps -Inf,
      ## Inf, and NaN.

      if(!all(is.finite(summand))) {
        .np_warning(" non-finite value in summation-based statistic: integration recommended")
        summand <- summand[is.finite(summand)]
      }

      return(0.5*mean((1-sqrt(summand))**2))
      
    } else {
      
      ## Integration version:
      ## \int\int (sqrt(f(x,y))-sqrt(f(x))sqrt(f(y)))^2 dx dy
      
      lo.default <- c(min(x.dat)-10.0*IQR(x.dat),min(y.dat)-10.0*IQR(y.dat))
      up.default <- c(max(x.dat)+10.0*IQR(x.dat),max(y.dat)+10.0*IQR(y.dat))
      
      return(.np_entropy_bivariate_integral(
        x.dat = x.dat,
        y.dat = y.dat,
        bw.x = bw.x,
        bw.y = bw.y,
        bw.joint = bw.joint,
        lower = lo.default,
        upper = up.default
      ))
      
    }

  } ## end of Srho.bivar function

  ## Compute the metric entropy for lags 1 through lag.num

  Srho.vec <- numeric()
  ## `Portmanteau' cumulant of all lags
  Srho.cumulant.vec <- numeric()
  
  bw.y <- numeric()
  bw.y.lag <- numeric()
  bw.joint.y <- numeric()
  bw.joint.y.lag <- numeric()      
  
  ## Save the bandwidths for resampling exercise...
  lag.progress <- .np_progress_begin("Constructing metric entropy by lag", total = lag.num, surface = "lag")
  
  for (k in seq_len(lag.num)) {
    ## Create y and y.lag
    
    tmp <- ts.intersect(as.ts(data),lag(as.ts(data),k))
    y <- as.numeric(tmp[,1])
    y.lag <- as.numeric(tmp[,2])
    rm(tmp)

    ## Compute and save bandwidths (save for bootstrapping if
    ## requested)

    bw.y[k] <- .np_progress_with_legacy_suppressed(npudensbw(~y))$bw
    bw.y.lag[k] <- .np_progress_with_legacy_suppressed(npudensbw(~y.lag))$bw
    bw.joint <- .np_progress_with_legacy_suppressed(npudensbw(~y+y.lag))$bw
    bw.joint.y[k] <- bw.joint[1]
    bw.joint.y.lag[k] <- bw.joint[2]
    
    Srho.vec[k] <- Srho.bivar(y,y.lag,bw.y[k],bw.y.lag[k],bw.joint,method=method)
    lag.progress <- .np_progress_step(lag.progress, done = k, detail = paste("lag", k))

  }

  lag.progress <- .np_progress_end(lag.progress, detail = paste("lag", lag.num))

  for (k in seq_len(lag.num)) Srho.cumulant.vec[k] <- sum(Srho.vec[seq_len(k)])

  ## Bootstrap if requested - null is independence so simple iid
  ## index resampling under replacement is sufficient

  if(bootstrap) {

    Srho.vec.boot <- numeric()
    ## `Portmanteau' cumulant of all lags
    Srho.cumulant.vec.boot <- numeric()

    ## Matrix for resamples

    Srho.bootstrap.mat <- matrix(NA,boot.num,(lag.num))
		Srho.cumulant.bootstrap.mat <- matrix(NA,boot.num,(lag.num))
    progress <- .np_progress_begin("Bootstrap replications", total = boot.num, surface = "bootstrap")

    for (b in seq_len(boot.num)) {

      ## Resample under the null

      resampled.ts <- as.ts(data[sample.int(length(data), replace = TRUE)])

      for (k in seq_len(lag.num)) {
        tmp <- ts.intersect(resampled.ts, lag(resampled.ts, k))
        y <- as.numeric(tmp[,1])
        y.lag <- as.numeric(tmp[,2])
        Srho.vec.boot[k] <- Srho.bivar(y, y.lag, bw.y[k], bw.y.lag[k], c(bw.joint.y[k], bw.joint.y.lag[k]), method = method)
        Srho.bootstrap.mat[b,k] <- Srho.vec.boot[k]
        ## `Portmanteau' cumulant of all lags
        Srho.cumulant.vec.boot[k] <- sum(Srho.vec.boot[seq_len(k)])
        Srho.cumulant.bootstrap.mat[b,k] <- Srho.cumulant.vec.boot[k]
      }

      progress <- .np_progress_step(progress, done = b)

    }

    progress <- .np_progress_end(progress)

    ## Compute P-values

    P.vec <- numeric()
    P.cumulant.vec <- numeric()
    
    for (k in seq_len(lag.num)) {
      P.vec[k] <- mean(Srho.bootstrap.mat[,k] > Srho.vec[k])
      P.cumulant.vec[k] <- mean(Srho.cumulant.bootstrap.mat[,k] > Srho.cumulant.vec[k])
    }

  }

  ## Restore seed
  
  .np_seed_exit(seed.state, remove_if_absent = TRUE)
  
  if(bootstrap) {
    
    sdeptest(Srho = Srho.vec,
             Srho.cumulant = Srho.cumulant.vec,
             Srho.bootstrap.mat = Srho.bootstrap.mat,
             Srho.cumulant.bootstrap.mat = Srho.cumulant.bootstrap.mat,
             P = P.vec,
             P.cumulant = P.cumulant.vec,
             bootstrap = bootstrap,
             boot.num = boot.num,
             lag.num = lag.num,
             bw.y = bw.y,
             bw.y.lag = bw.y.lag,
             bw.joint = cbind(bw.joint.y,bw.joint.y.lag))
    
  } else {
    
    sdeptest(Srho = Srho.vec,
             Srho.cumulant = Srho.cumulant.vec,
             Srho.bootstrap.mat = NULL,
             Srho.cumulant.bootstrap.mat = NULL,
             P = NULL,
             P.cumulant = NULL,
             bootstrap = bootstrap,
             boot.num = NULL,
             lag.num = lag.num,
             bw.y = bw.y,
             bw.y.lag = bw.y.lag,
             bw.joint = cbind(bw.joint.y,bw.joint.y.lag))
    
  }

}
