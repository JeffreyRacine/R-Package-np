npquantile <- function(bws=NULL,
                       dat=NULL,
                       tau=c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),
                       num.eval=5000) {

  ## Note - options for npudist are passed in/taken from the bandwidth
  ## object.

  ## Some basic error checking.

  if(is.null(bws)) stop("must provide bw object")
  if(is.null(dat)) stop("must provide data")
  if(class(dat) != "numeric") stop("dat must be numeric and univariate")

  if(any(tau<=0)|any(tau>=1)) stop("tau must lie in the open interval (0,1)")
  if(length(bws$xnames)>1) stop("bw object must be univariate")
  if(class(bws)!="dbandwidth") stop("bw object must be a npudistbw() object")
  if(num.eval < 100) stop("num.eval must be >= 100")

  ## Create grid from which quasi-inverse is extracted - extend the
  ## range of dat for evaluation grid, also add empirical quantiles to
  ## grid. For finer grain grid increase length=5000 below.

  dat.er <- extendrange(dat,f=1)
  dat.eval <- sort(c(seq(dat.er[1],dat.er[2],length=num.eval),
                     quantile(dat,tau)))

  F <- fitted(npudist(tdat=dat,
                      edat=dat.eval,
                      bws=bws))

  ## Now compute the quasi-inverse from the estimated F for the
  ## evaluation points. If tau is input and any value lies beyond the
  ## CDF values for the evaluation points, reset them to the min/max
  ## CDF values for the evaluation data (otherwise the quantiles are
  ## undefined).

  dat.tau <- numeric(length(tau))
  
  for(i in 1:length(tau)) {
    tau[tau<min(F)] <- min(F)
    tau[tau>max(F)] <- max(F)        
    dat.tau[i] <-  min(dat.eval[F>=tau[i]])
  }

  return(dat.tau)

}
