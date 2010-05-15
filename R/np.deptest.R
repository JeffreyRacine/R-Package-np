## Function that implements the bivariate dependence metric described
## in Maasoumi, E.~and J.S.~Racine (2002), "Entropy and Predictability
## of Stock Market Returns," Journal of Econometrics, March, Volume
## 107, Issue 2, pp 291-312.

npdeptest <- function(data.x = NULL,
                      data.y = NULL,
                      method=c("integration","summation"),
                      bootstrap = TRUE,
                      boot.num = 399,
                      random.seed = 42) {
  
  ## Trap fatal errors

  if(is.data.frame(data.x)||is.data.frame(data.y)) stop(" you must enter two data vectors (and not data frames)")
  if(is.factor(data.x)||is.factor(data.y)) stop(" does not support factors")
  if(is.null(data.x)||is.null(data.y)) stop(" you must enter x and y data vectors")
  if(ncol(data.frame(data.x)) != 1) stop(" data must have one dimension only")
  if(length(data.x)!=length(data.y)) stop(" data vectors must be of equal length")
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

  if(is.ts(data.x)) data.x <- as.numeric(data.x)
  if(is.ts(data.y)) data.y <- as.numeric(data.y)

  ## Remove any NAs from paired data

  tmp <- na.omit(data.frame(data.x,data.y))
  data.x <- tmp$data.x
  data.y <- tmp$data.y
  rm(tmp)

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
      
      ## Summation version: \sum_i(1-sqrt(f(x_i)f(y_i)/f(x_i,y_i)))^2
      ## We could code this by hand per below, however, there is no
      ## performance penalty imposed by sum() so no need.
      
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
  
  console <- newLineConsole()
  console <- printClear(console)
  console <- printPop(console)  
  
  ## Compute and save bandwidths (save for bootstrapping if requested)

  bw.data.x <- npudensbw(~data.x)$bw
  bw.data.y <- npudensbw(~data.y)$bw
  bw.joint <- npudensbw(~data.x+data.y)$bw
  
  console <- printClear(console)
  console <- printPush(paste(sep="", "Constructing metric entropy..."), console = console)
  
  Srho.vec <- Srho.bivar(data.x,data.y,bw.data.x,bw.data.y,bw.joint,method=method)

  ## Bootstrap if requested - null is independence so simple iid
  ## bootstrap (sample(x,replace=TRUE)) will work

  if(bootstrap) {

    Srho.vec.boot <- numeric()

    for(b in 1:boot.num) {

      console <- printClear(console)
      console <- printPush(paste(sep="", "Bootstrap replication ",
                                 b, "/", boot.num, "..."), console)
      
      ## Break systematic relationship between x and y (null)
      
      data.x.boot <- sample(data.x,replace=TRUE)
      
      Srho.vec.boot[b] <- Srho.bivar(data.x.boot,data.y,bw.data.x,bw.data.y,bw.joint,method=method)

    }

    ## Compute P-values

    P <- mean(ifelse(Srho.vec.boot>Srho.vec,1,0))

  }

  console <- printClear(console)
  console <- printPop(console)

  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  if(bootstrap) {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = Srho.vec.boot,
            P = P,
            bootstrap = bootstrap,
            boot.num = boot.num,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  } else {

    deptest(Srho = Srho.vec,
            Srho.bootstrap.vec = NULL,
            P = NULL,
            bootstrap = bootstrap,
            boot.num = NULL,
            bw.data.x = bw.data.x,
            bw.data.y = bw.data.y,
            bw.joint = bw.joint)

  }

}
