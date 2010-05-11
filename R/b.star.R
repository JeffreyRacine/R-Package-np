## $Id: ppw.R,v 1.47 2008/12/12 14:52:17 jracine Exp jracine $

## Original code in Matlab by A. Patton, R translation and
## modifications by C. Parmeter and J. Racine.
##
## We are grateful to Andrew Patton and Dimitris Politis for their
## assistance and feedback. Kindly report features, deficiencies, and
## improvements to racinej@mcmaster.ca.
##
## The citation is A. Patton, D.N. Politis, and H. White (2008,
## forthcoming), "CORRECTION TO `Automatic Block-Length Selection for
## the Dependent Bootstrap' by D.N. Politis and H. White". This is
## based on the article by Politis, D.N., and H. White (2004),
## "Automatic block-length selection for the dependent bootstrap."
## Econometric Reviews, vol. 23.
##
## INPUTS:  data, an n x k matrix.
##
## OUTPUTS: b.star, a 2 x k vector of optimal bootstrap block lengths
## for the stationary bootstrap and circular bootstrap (BstarSB,
## BstarCB).

## The function lam() is used to construct a "flat-top" lag window for
## spectral estimation based on Politis, D.N. and J.P. Romano (1995),
## "Bias-Corrected Nonparametric Spectral Estimation", Journal of Time
## Series Analysis, vol. 16, No. 1.

lam <- function(s){
  return((abs(s)>=0)*(abs(s)<0.5)+2*(1-abs(s))*(abs(s)>=0.5)*(abs(s)<=1))
}

## The function b.star() returns the optimal bootstrap block
## lengths. Note that an example for usage appears at the bottom of
## this file. If you use this function as input into a routine such as
## tsboot() in the boot library (Angelo Canty and Brian Ripley
## (2008). boot: Bootstrap R (S-Plus) Functions. R package version
## 1.2-34.) you ought to use the option round=TRUE.

b.star <- function(data,
                   Kn = NULL,
                   mmax= NULL,
                   Bmax = NULL,
                   c = NULL,
                   round = FALSE){

  ## Convert the data object to a data frame to handle both vectors
  ## and matrices.

  data <- data.frame(data)
  n <- nrow(data)
  k <- ncol(data)

  ## Set Defaults. Note that in footnote c, page 59, for Kn Politis
  ## and White (2004) use max(5,log10(n)). Since this must be an
  ## integer we use ceiling(log10(n)).

  if (is.null(Kn)) Kn <- max(5,ceiling(log10(n)))
  if (is.null(mmax)) mmax <- ceiling(sqrt(n))+Kn
  if (is.null(Bmax)) Bmax <- ceiling(min(3*sqrt(n),n/3))
  if (is.null(c)) c <- qnorm(0.975)

  ## Create two vectors of length k in which we store results.

  BstarSB <- numeric(length=k)
  BstarCB <- numeric(length=k)

  ## Now we loop through each variable in data (i.e., column,
  ## data[,i]).

  for(i in 1:k) {

    ## We first obtain the autocorrelations rho(1),...,rho(mmax) (we
    ## need to drop the first autocorrelation as it is rho(0), hence
    ## acf[-1]). This is the default in acf [type="correlation"]. Note
    ## that Patton uses sample correlations after dropping the first
    ## mmax observations, while we instead use the acf to obtain
    ## rho(k).

    rho.k <- acf(data[,i],
                 lag.max = mmax,
                 type = "correlation",
                 plot = FALSE)$acf[-1]

    ## Next we compute mhat. The use of c*sqrt(log10(n)/n) for
    ## critical values is given in footnote c of Politis and White
    ## (2004, page 59), and the approach for determining mhat is
    ## described in footnote c.

    rho.k.crit <- c*sqrt(log10(n)/n)

    ## Compute the number of insignificant runs following each rho(k),
    ## k=1,...,mmax.
    
    num.insignificant <- sapply(1:(mmax-Kn+1),
                                function(j){
                                  sum((abs(rho.k) < rho.k.crit)[j:(j+Kn-1)])
                                })

    ## If there are any values of rho(k) for which the Kn proceeding
    ## values of rho(k+j), j=1,...,Kn are all insignificant, take the
    ## smallest rho(k) such that this holds (see footnote c for
    ## further details).

    if(any(num.insignificant==Kn)) {
      mhat <- which(num.insignificant==Kn)[1]
    } else {

      ## If no runs of length Kn are insignificant, take the smallest
      ## value of rho(k) that is significant.

      if(any(abs(rho.k) > rho.k.crit)) {
        
        lag.sig <- which(abs(rho.k) > rho.k.crit)
        k.sig <- length(lag.sig)
        
        if(k.sig == 1) {
          
          ## When only one lag is significant, mhat is the sole
          ## significant rho(k).
          
          mhat <- lag.sig
          
        } else {

          ## If there are more than one significant lags but no runs
          ## of length Kn, take the largest value of rho(k) that is
          ## significant.
          
          mhat <- max(lag.sig)
          
        }
        
      } else {
        
        ## When there are no significant lags, mhat must be the
        ## smallest positive integer (footnote c), hence mhat is set
        ## to one.
        
        mhat <- 1
        
      }

    }

    ## Compute M (mhat is at least one).

    M <- ifelse(2*mhat > mmax, mmax, 2*mhat)

    ## We compute BstarSB and BstarCB using the formulas in the above
    ## references. Now we require the autocovariance R(k) (hence
    ## type="covariance" in the acf call). Note that Patton uses
    ## sample covariances after dropping the first mmax observations,
    ## while we instead use the acf with type="covariance" to obtain
    ## R(k). Note also that we require R(0) hence we do not drop it as
    ## we did for rho(k) via acf(...)$acf[-1].
      
    kk <- seq(-M,M)

    R.k <- ccf(data[,i], data[,i],
               lag.max = M,
               type = "covariance",
               plot = FALSE)$acf
    
    Ghat <- sum(lam(kk/M)*abs(kk)*R.k)
    DCBhat <- 4/3*sum(lam(kk/M)*R.k)^2
    DSBhat <- 2*sum(lam(kk/M)*R.k)^2
    BstarSB[i] <- ((2*Ghat^2)/DSBhat)^(1/3)*n^(1/3)
    BstarCB[i] <- ((2*(Ghat^2)/DCBhat)^(1/3))*(n^(1/3))
      
  }

  ## The user can choose whether they want rounded values returned or
  ## not. BstarCB is rounded up, BstarSB simply rounded but both must
  ## be positive integers.

  if(round == FALSE) {
    
    BstarSB <- ifelse(BstarSB > Bmax, Bmax, BstarSB)
    BstarCB <- ifelse(BstarCB > Bmax, Bmax, BstarCB)

  } else {

    BstarSB <- ifelse(BstarSB > Bmax, Bmax, ifelse(BstarSB < 1, 1, round(BstarSB)))
    BstarCB <- ifelse(BstarCB > Bmax, Bmax, ifelse(BstarCB < 1, 1, ceiling(BstarCB)))

  }
  
  return(cbind(BstarSB,BstarCB)) 
  
}

## Here is a simple example with an n x 2 matrix containing n=10^5
## observations, where column 2 of x is more persistent than column
## 1. This requires that you first install the forecast library (i.e.,
## install.packages("forecast")).
##
##  library(forecast)
##  set.seed(123)
##  x <- cbind(arima.sim(n = 100000, list(ar = c(.5,.0), ma = c(0,0)),sd = 1),
##             arima.sim(n = 100000, list(ar = c(.5,.4), ma = c(0,0)),sd = 1))
##  b.star(x)  
##  b.star(x,round=TRUE)
##>   b.star(x)  
##       BstarSB   BstarCB
##[1,]  50.39272  57.68526
##[2,] 251.62894 288.04323
##>   b.star(x,round=TRUE)
##     BstarSB BstarCB
##[1,]      50      58
##[2,]     252     289
