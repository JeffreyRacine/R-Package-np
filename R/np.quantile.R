npquantile <- function(x=NULL,
                       tau=c(0.01,0.05,0.25,0.50,0.75,0.95,0.99),
                       num.eval=10000,
                       bws=NULL,
                       f=1,
                       ...) {

  ## Some basic error checking.

  if(is.null(x)) stop("must provide data")
  if(class(x) != "numeric") stop("x must be numeric and univariate")

  if(any(tau<0)|any(tau>1)) stop("tau must lie in the closed interval [0,1]")
  if(length(bws$xnames)>1) stop("bw object must be univariate")
  if(num.eval < 100) stop("num.eval must be >= 100")

  if(is.null(bws)) bws <- npudistbw(~x,...)
  if(class(bws)!="dbandwidth") stop("bw object must be a npudistbw() object")

  ## Create grid from which quasi-inverse is extracted - extend the
  ## range of x for evaluation grid, also add empirical quantiles to
  ## grid.

  x.er <- extendrange(x,f=f)
  x.eval <- na.omit(sort(c(seq(x.er[1],x.er[2],length=num.eval),
                           quantile(x,tau,na.rm=TRUE))))

  F <- fitted(npudist(tdat=x,
                      edat=x.eval,
                      bws=bws))

  ## Now compute the quasi-inverse from the estimated F for the
  ## evaluation points. If tau is input and any value lies beyond the
  ## CDF values for the evaluation points, reset them to the min/max
  ## CDF values for the evaluation data (otherwise the quantiles are
  ## undefined).

  x.tau <- numeric(length(tau))
  
  for(i in 1:length(tau)) {
    tau[tau<min(F)] <- min(F)
    tau[tau>max(F)] <- max(F)        
    x.tau[i] <-  min(x.eval[F>=tau[i]])
  }

  names(x.tau) <- paste(tau*100,"%",sep="")

  return(x.tau)

}
