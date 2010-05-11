## Function that implements the entropy metric test for serial
## dependence described in "A Dependence Metric For Possibly Nonlinear
## Processes By C. W. Granger, E. Maasoumi and J. Racine Journal of
## Time Series Analysis (2004), Vol. 25, No. 5, 649-669.

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

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }
  
  set.seed(random.seed)

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
      
      f.x <- fitted(npudens(tdat=x.dat,bws=bw.x))
      f.y <- fitted(npudens(tdat=y.dat,bws=bw.y))
      f.xy <- fitted(npudens(tdat=cbind(x.dat,y.dat),bws=bw.joint))
      summand <- f.x*f.y/f.xy

      ## In summation version we divide densities which can lead to
      ## numerical instability issues not present in the integrand
      ## version (which uses differences instead). We check for this
      ## case, remove offending points, and warn. This traps -Inf,
      ## Inf, and NaN.

      if(!all(is.finite(summand))) {
        warning(" non-finite value in summation-based statistic: integration recommended")
        summand <- summand[is.finite(summand)]
      }

      return(0.5*mean((1-sqrt(summand))**2))
      
    } else {
      
      ## Integration version:
      ## \int\int (sqrt(f(x,y))-sqrt(f(x))sqrt(f(y)))^2 dx dy
      
      lo.default <- c(min(x.dat)-10.0*IQR(x.dat),min(y.dat)-10.0*IQR(y.dat))
      up.default <- c(max(x.dat)+10.0*IQR(x.dat),max(y.dat)+10.0*IQR(y.dat))
      
      Srho.integrand <- function(xy) {
        
        ## We code this up by hand using a second order gaussian
        ## kernel. The multidimensional integration routines for some
        ## reason crawl to almost a halt if we call the functions in the
        ## np package. Good to have noticed this at this juncture, bad
        ## to not have general code though could branch at this point?
        
        ##      f.x <- fitted(npudens(tdat=x.dat,edat=xy[1],bws=bw.x))
        ##      f.y <- fitted(npudens(tdat=y.dat,edat=xy[2],bws=bw.y))
        ##      f.xy <- fitted(npudens(tdat=cbind(x.dat,y.dat),edat=cbind(xy[1],xy[2]),bws=bw.joint))


        f.x <- mean(dnorm((xy[1]-x.dat)/bw.x))/bw.x
        f.y <- mean(dnorm((xy[2]-y.dat)/bw.y))/bw.y
        f.xy <- mean(dnorm((xy[1]-x.dat)/bw.joint[1])*dnorm((xy[2]-y.dat)/bw.joint[2])
                    / (bw.joint[1]*bw.joint[2]))
        
        return((sqrt(f.xy)-sqrt(f.x)*sqrt(f.y))**2)
        
      }
      
      return(0.5*adaptIntegrate(Srho.integrand,
                                lowerLimit=lo.default,
                                upperLimit=up.default)$integral)
      
    }

  } ## end of Srho.bivar function

  ## Compute the metric entropy for lags 1 through lag.num

  console <- newLineConsole()

  Srho.vec <- numeric()
  ## `Portmanteau' cumulant of all lags
  Srho.cumulant.vec <- numeric()
  
  bw.y <- numeric()
  bw.y.lag <- numeric()
  bw.joint.y <- numeric()
  bw.joint.y.lag <- numeric()      
  
  ## Save the bandwidths for resampling exercise...
  
  for(k in 1:lag.num) {
    
    console <- printClear(console)
    console <- printPop(console)  
    
    ## Create y and y.lag
    
    tmp <- ts.intersect(as.ts(data),lag(as.ts(data),k))
    y <- as.numeric(tmp[,1])
    y.lag <- as.numeric(tmp[,2])
    rm(tmp)

    ## Compute and save bandwidths (save for bootstrapping if
    ## requested)

    bw.y[k] <- npudensbw(~y)$bw
    bw.y.lag[k] <- npudensbw(~y.lag)$bw
    bw.joint <- npudensbw(~y+y.lag)$bw
    bw.joint.y[k] <- bw.joint[1]
    bw.joint.y.lag[k] <- bw.joint[2]

    console <- printClear(console)
    console <- printPush(paste(sep="", "Constructing metric entropy at lag ", k, "/", lag.num, "..."), console = console)
    
    Srho.vec[k] <- Srho.bivar(y,y.lag,bw.y[k],bw.y.lag[k],bw.joint,method=method)

  }

  for(k in 1:lag.num) Srho.cumulant.vec[k] <- sum(Srho.vec[1:k])

  ## Bootstrap if requested - null is independence so simple iid
  ## bootstrap (sample(x,replace=TRUE)) will work

  if(bootstrap) {

    Srho.vec.boot <- numeric()
    ## `Portmanteau' cumulant of all lags
    Srho.cumulant.vec.boot <- numeric()

    ## Matrix for resamples

    Srho.bootstrap.mat <- matrix(NA,boot.num,(lag.num))
		Srho.cumulant.bootstrap.mat <- matrix(NA,boot.num,(lag.num))

    for(b in 1:boot.num) {

      console <- printClear(console)
      console <- printPush(paste(sep="", "Bootstrap replication ",
                                 b, "/", boot.num, "..."), console)

      ## Resample under the null

      tmp <- as.ts(sample(data,replace=TRUE))
      tmp <- ts.intersect(tmp,lag(tmp,k))
      y <- as.numeric(tmp[,1])
      y.lag <- as.numeric(tmp[,2])
      rm(tmp)

      for(k in 1:lag.num) {
        Srho.vec.boot[k] <- Srho.bivar(y,y.lag,bw.y[k],bw.y.lag[k],c(bw.joint.y[k],bw.joint.y.lag[k]),method=method)
        Srho.bootstrap.mat[b,k] <- Srho.vec.boot[k]
				## `Portmanteau' cumulant of all lags
				Srho.cumulant.vec.boot[k] <- sum(Srho.vec.boot[1:k])
				Srho.cumulant.bootstrap.mat[b,k] <- Srho.cumulant.vec.boot[k]          
      }

    }

    ## Compute P-values

    P.vec <- numeric()
    P.cumulant.vec <- numeric()
    
    for(k in 1:lag.num) {
      P.vec[k] <- mean(ifelse(Srho.bootstrap.mat[,k]>Srho.vec[k],1,0))
      P.cumulant.vec[k] <- mean(ifelse(Srho.cumulant.bootstrap.mat[,k]>Srho.cumulant.vec[k],1,0))
    }

  }

  console <- printClear(console)
  console <- printPop(console)
  
  ## Restore seed
  
  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
  
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
