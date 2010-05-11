## Function that implements the univariate test for equality of
## densities described in Maasoumi, E.~and J.S.~Racine (2002),
## "Entropy and Predictability of Stock Market Returns," Journal of
## Econometrics, March, Volume 107, Issue 2, pp 291-312.

npunitest <- function(data.x = NULL,
                      data.y = NULL,
                      method = c("integration","summation"),
                      bootstrap = TRUE,
                      boot.num = 399,
                      bw.x = NULL,
                      bw.y = NULL,                         
                      random.seed = 42,
                      ...) {

  if(is.null(data.x) || is.null(data.y)) stop(" you must enter data vectors for x and y")
  if(is.data.frame(data.x) || is.data.frame(data.y)) stop(" you must enter data vectors (not data frames)")
  if(class(data.x) != class(data.y)) stop(" data vectors must be of same data type")
  if((ncol(data.frame(data.x)) != 1) ||( ncol(data.frame(data.y)) != 1)) stop(" data vectors must have one dimension only")
  if(boot.num < 9) stop(" number of bootstrap replications must be >= 9")
  if(is.numeric(data.x) && (max(data.x) < min(data.y) || max(data.y) < min(data.x))) warning("non-overlapping empirical distributions (see `Details' in ?npunidist)")

  method <- match.arg(method)

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  console <- newLineConsole()
  console <- printPush("Computing bandwidths...", console = console)

  ## If of type ts convert to numeric to handle time series data

  if(is.ts(data.x)) data.x <- as.numeric(data.x)
  if(is.ts(data.y)) data.y <- as.numeric(data.y)

  ## Remove NAs

  data.x <- na.omit(data.x)
  data.y <- na.omit(data.y)

  ## Note - this function accepts numeric and factors/ordered,
  ## however, factors must be integers.

  if(is.numeric(data.x)) {
    if(is.null(bw.x)) bw.x <- bw.SJ(data.x)
    if(is.null(bw.y)) bw.y <- bw.SJ(data.y)
  } else {
    ## Optimal bandwidth for factor (unordered) used for both
    ## (plug-in).
    if(is.null(bw.x)) {
      n <- length(data.x)
      c <- length(unique(data.x))
      xeval <- unique(data.x)
      p <- fitted(npudens(tdat=data.x,edat=xeval,bws=0,...))
      sum.Lambda3 <- c/(c-1)*sum(p*(1-p))
      sum.Lambda2.minus.Lambda1.sq <- c^2/((c-1)^2)*sum(p*(1-p))
      n.sum.Lambda1.sq <- n*sum(((1-c*p)/(c-1))^2)
      bw.x <- sum.Lambda3/(sum.Lambda2.minus.Lambda1.sq+n.sum.Lambda1.sq)
    }
    if(is.null(bw.y)) {
      n <- length(data.y)
      c <- length(unique(data.y))
      yeval <- unique(data.y)
      p <- fitted(npudens(tdat=data.y,edat=xeval,bws=0,...))
      sum.Lambda3 <- c/(c-1)*sum(p*(1-p))
      sum.Lambda2.minus.Lambda1.sq <- c^2/((c-1)^2)*sum(p*(1-p))
      n.sum.Lambda1.sq <- n*sum(((1-c*p)/(c-1))^2)
      bw.y <- sum.Lambda3/(sum.Lambda2.minus.Lambda1.sq+n.sum.Lambda1.sq)
    }
  }
	
  ## Hellinger distance (Srho) between the densities. This function
  ## accepts numeric and factor/ordered.
  
  Srho.univar <- function(data.x,
                          data.y,
                          bw.x,
                          bw.y,
                          method=c("integration","summation")) {

    if((ncol(data.frame(data.x)) != 1) || (ncol(data.frame(data.y)) != 1))  stop(" data vectors must have one dimension only") 
    if(is.numeric(data.x)) {
      if(method=="summation") {
        ## Summation - we use positive support data (data.x) when
        ## evaluating. This is proper and avoids division by
        ## zero. Note the difference between integration and summation
        ## in this case when the variables do not share common support.
        f.data.x <- fitted(npudens(tdat=data.x,edat=data.x,bws=bw.x,...))
        f.data.y <- fitted(npudens(tdat=data.y,edat=data.x,bws=bw.y,...))
        summand <- f.data.y/f.data.x
        ## In summation version we divide densities which can lead to
        ## numerical instability issues not present in the integrand
        ## version (which uses differences instead). We check for this
        ## case, remove offending points, and warn. This traps -Inf,
        ## Inf, and NaN.
        if(!all(is.finite(summand))) {
          warning(" non-finite value in summation-based statistic: integration recommended")
          summand <- summand[is.finite(summand)]
        }
        Srho <- 0.5*mean((1-sqrt(summand))**2)
        if(Srho < 0 || Srho > 1) warning(" numerical instability in summation-based statistic: integration recommended")
        return(Srho)
      } else {
        ## Integration
        h <- function(x,data.x,data.y) {
          f.data.x <- fitted(npudens(tdat=data.x,edat=x,bws=bw.x,...))
          f.data.y <- fitted(npudens(tdat=data.y,edat=x,bws=bw.y,...))
          return(0.5*(sqrt(f.data.x)-sqrt(f.data.y))**2)
        }
        return.integrate <- integrate(h,-Inf,Inf,subdivisions=1e+05,stop.on.error=FALSE,data.x=data.x,data.y=data.y)      
        if(return.integrate$message != "OK") warning(return.integrate$message)
        return(return.integrate$value)
      }
    } else {
      xeval <- unique(data.x)
      p.data.x <- fitted(npudens(tdat=data.x,edat=xeval,bws=bw.x,...))
      p.data.y <- fitted(npudens(tdat=data.y,edat=xeval,bws=bw.y,...))   
      ## Sum the information function over the unique probabilities and
      ## return.
      return(sum(0.5*(sqrt(p.data.x)-sqrt(p.data.y))**2))
    }
  }

  
  ## Compute the test statistic

  console <- printClear(console)
  console <- printPush("Computing test statistic...", console = console)

  test.stat <- Srho.univar(data.x,data.y,bw.x,bw.y,method=method)

  console <- printClear(console)
  console <- printPop(console)  

  if(bootstrap) {
    
    if(method == "summation") {
      
      ## For summation since we evaluate Srho at the support data.x
      ## (not data.y) we resample from data.x only
      
      data.null <- data.x
      
    } else {
      
      ##  For integration impose the null by resampling from pooled
      ##  data.x and data.y
      
      if(is.numeric(data.x)) {
        data.null <- c(data.x,data.y)
      } else {
        if(is.ordered(data.x)) {
          data.null <- ordered(c(as.character(data.x),as.character(data.y)))
        } else {
          data.null <- factor(c(as.character(data.x),as.character(data.y)))
        }
      }
      
    }
    
    resampled.stat <- numeric(boot.num)

    for(b in 1:boot.num) {
      
      console <- printClear(console)
      console <- printPush(paste(sep="", "Bootstrap replication ",
                                  b, "/", boot.num, "..."), console = console)

      ## Need to think this through... is the null one density? If so
      ## resample from that density for both x and y?
      
      ## Conduct simple iid bootstrap resamples
      
      data.null.x <- data.null[sample(1:length(data.null),length(data.x),replace=TRUE)]
      data.null.y <- data.null[sample(1:length(data.null),length(data.y),replace=TRUE)]

      resampled.stat[b] <- Srho.univar(data.null.x,data.null.y,bw.x,bw.y,method=method)
      
    }
    
    p.value <- mean(ifelse(resampled.stat > test.stat, 1, 0))
    
  }
  
  console <- printClear(console)
  console <- printPop(console)  

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  if(bootstrap) {
    unitest(Srho = test.stat,
               Srho.bootstrap = resampled.stat,
               P = p.value,
               bootstrap = bootstrap,
               boot.num = boot.num,
               bw.x = bw.x,
               bw.y = bw.y)
  } else {
    unitest(Srho = test.stat,
               Srho.bootstrap = NULL,
               P = NULL,
               bootstrap = bootstrap,
               boot.num = NULL,
               bw.x = bw.x,
               bw.y = bw.y)
  }

}
