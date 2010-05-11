## Function that implements the test for asymmetry applicable to both
## real-values and categorical univariate data described in Maasoumi,
## E. and J.S. Racine (2009), "A Robust Entropy-Based Test of
## Asymmetry for Discrete and Continuous Processes," Econometric
## Reviews, Volume 28, pp 246-261.

npsymtest <- function(data = NULL,
                      method = c("integration","summation"),
                      boot.num = 399,
                      bw = NULL,
                      boot.method = c("iid", "geom"),
                      random.seed = 42,
                      ...) {

  if(is.data.frame(data)) stop(" you must enter a data vector (not data frame)")
  if(is.null(data)) stop(" you must enter a data vector")
  if(ncol(data.frame(data)) != 1) stop(" data must have one dimension only")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")

  boot.method <- match.arg(boot.method)
  method <- match.arg(method)

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  console <- newLineConsole()
  console <- printPush(paste(sep="", "Working..."), console = console)
  console <- printPop(console)

  ## If of type ts convert to numeric to handle time series data

  if(is.ts(data)) data <- as.numeric(data)

  ## Remove NAs

  data <- na.omit(data)

  ## Note - this function accepts numeric and factors/ordered,
  ## however, factors must be integers.

  if(is.numeric(data)) {
    ## Numerical stability - integrate does best centered around zero,
    ## so standardize numeric data - this affects nothing else.
    data <- (data-mean(data))/sd(data)
    data.rotate <- -(data-mean(data))+mean(data)
    if(is.null(bw)) bw <- bw.SJ(data)
  } else {
    if(is.ordered(data)) {
      ## Rotate around the median, ordered case.
      data.levels <- levels(data)
      tmp <- as.numeric(data.matrix(data))
      location <- median(sort(unique(tmp)))
      data.rotate <- ordered(-(tmp-location)+location,levels=data.levels)
    } else {
      ## unordered case.
      data.levels <- levels(data)
      tmp <- as.numeric(data.matrix(data))
      location <- median(sort(unique(tmp)))
      data.rotate <- factor(-(tmp-location)+location,levels=data.levels)
    }
    ## Optimal bandwidth for factor (unordered) used for both
    ## (plug-in).
    if(is.null(bw)) {
      n <- length(data)
      c <- length(unique(data))
      xeval <- unique(data)
      p <- fitted(npudens(tdat=data,edat=xeval,bws=0,...))
      sum.Lambda3 <- c/(c-1)*sum(p*(1-p))
      sum.Lambda2.minus.Lambda1.sq <- c^2/((c-1)^2)*sum(p*(1-p))
      n.sum.Lambda1.sq <- n*sum(((1-c*p)/(c-1))^2)
      bw <- sum.Lambda3/(sum.Lambda2.minus.Lambda1.sq+n.sum.Lambda1.sq)
    }
  }
	
  ##  Impose the null and generate a data of resampled statistics

  if(is.numeric(data)) {
    data.null <- c(data,data.rotate)
  } else {
    if(is.ordered(data)) {
      data.null <- ordered(c(as.character(data),as.character(data.rotate)))
    } else {
      data.null <- factor(c(as.character(data),as.character(data.rotate)))
    }
  }

  ## Hellinger distance (Srho) between the densities for the data and
  ## the rotated data. This function accepts numeric and
  ## factor/ordered.
  
  Srho.sym <- function(data,
                       data.rotate,
                       bw,
                       method=c("integration","summation")) {
    
    if(ncol(data.frame(data)) != 1)  stop(" data must have one dimension only") 

    if(is.numeric(data)) {
      if(method=="summation") {
        ## Summation version uses
        ##  0.5\sum_i(1-sqrt(f.data.rotate/f.data)^2). Note that here
        ##  we must evaluate on common sample, so we use `data' for
        ##  the evaluation points. We could also combine both data and
        ##  data.rotate for the evaluation points if we chose (would
        ##  this improve power?).
        f.data <- fitted(npudens(tdat=data,edat=data,bws=bw,...))
        f.data.rotate <- fitted(npudens(tdat=data.rotate,edat=data,bws=bw,...))
        ## In summation version we divide densities which can lead to
        ## numerical instability issues not present in the integrand
        ## version (which uses differences instead). We check for this
        ## case, remove offending points, and warn. This traps -Inf,
        ## Inf, and NaN.
        summand <- f.data.rotate/f.data
        if(!all(is.finite(summand))) {
          warning(" non-finite value in summation-based statistic: integration recommended")
          summand <- summand[is.finite(summand)]
        }
        return(0.5*mean((1-sqrt(summand))**2))
      } else {
        ## Integration version
        h <- function(x,data,data.rotate) {
          f.data <- fitted(npudens(tdat=data,edat=x,bws=bw,...))
          f.data.rotate <- fitted(npudens(tdat=data.rotate,edat=x,bws=bw,...))
          return(0.5*(sqrt(f.data)-sqrt(f.data.rotate))**2)
        }
        return.integrate <- integrate(h,-Inf,Inf,subdivisions=1e+05,stop.on.error=FALSE,data=data,data.rotate=data.rotate)
        if(return.integrate$message != "OK") warning(return.integrate$message)
        return(return.integrate$value)
      }
    } else {
      xeval <- unique(data)
      p.data <- fitted(npudens(tdat=data,edat=xeval,bws=bw,...))
      p.data.rotate <- fitted(npudens(tdat=data.rotate,edat=xeval,bws=bw,...))   
      ## Sum the information function over the unique probabilities and
      ## return.
      return(sum(0.5*(sqrt(p.data)-sqrt(p.data.rotate))**2))
    }
  }
  
  ## Compute the test statistic

  test.stat <- Srho.sym(data,data.rotate,bw,method=method)

  ## Function to be fed to boot - accepts data that gets
  ## permuted/rearranged to define resampled data. The sole difference
  ## between boot.fun.boot and boot.fun.tsboot is the order of
  ## arguments.

  B.counter <- 0

  ## Function to be fed to tsboot - accepts a vector of integers
  ## corresponding to all observations in the sample (1,2,...) that
  ## get permuted/rearranged to define resampled data.

	boot.fun <- function(ii,data.null,bw) {
    console <<- printClear(console)
    console <<- printPush(paste(sep="", "Bootstrap replication ",
                                    B.counter, "/", boot.num, "..."), console = console)
    null.sample1 <- data.null[ii]
    if(is.numeric(data.null)) {
      null.sample2 <- -(null.sample1-mean(null.sample1))+mean(null.sample1)
    } else {
      if(is.ordered(data)) {
        ## Rotate around the median, ordered case
        null.sample1.levels <- levels(null.sample1)
        tmp <- as.numeric(data.matrix(null.sample1))
        location <- median(sort(unique(tmp)))
        null.sample2 <- ordered(-(tmp-location)+location,levels=data.levels)
      } else {
        ## unordered case
        null.sample1.levels <- levels(null.sample1)
        tmp <- as.numeric(data.matrix(null.sample1))
        location <- median(sort(unique(tmp)))
        null.sample2 <- factor(-(tmp-location)+location,levels=data.levels)
      }
    }
    B.counter <<- B.counter + 1
    return(Srho.sym(null.sample1,null.sample2,bw,method=method))
	}

  ## Need to bootstrap integers for data.null to accommodate both
  ## numeric and factor/ordered.

  if(boot.method == "iid")  {

    ## data.null is of length 2*n - if we use boot() we automatically
    ## get resamples of length 2*n - however, we want resamples from
    ## data.null of length n. Solution here is to use the block
    ## bootstrap with a block length of 1 which is presumed to
    ## generate an iid bootstrap.

    resampled.stat <- tsboot(tseries = 1:(2*length(data)),
                             statistic = boot.fun,
                             R = boot.num,
                             n.sim = length(data),
                             l = 1,
                             sim = "fixed",
                             data.null = data.null,
                             bw = bw)$t

  } else {

    if(is.numeric(data)) {
      ## Optimal block length is based upon `data' (not data.null).
      boot.blocklen <- b.star(data,round=TRUE)[1,1]
    } else {
      boot.blocklen <- b.star(as.numeric(data.matrix(data)),round=TRUE)[1,1]
    }

    resampled.stat <- tsboot(tseries = 1:(2*length(data)),
                             statistic = boot.fun,
                             R = boot.num,
                             n.sim = length(data),
                             l = boot.blocklen,
                             sim = "geom",
                             data.null = data.null,
                             bw = bw)$t

  }

  console <- printClear(console)
  console <- printPop(console)  

  p.value <- mean(ifelse(resampled.stat > test.stat, 1, 0))

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
  symtest(Srho = test.stat,
          Srho.bootstrap = resampled.stat,
          P = p.value,
          boot.num = boot.num,
          data.rotate = data.rotate,
          bw = bw)

}
