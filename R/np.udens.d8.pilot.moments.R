# Plot-only complete covariance under the unchanged-bandwidth positive pilot.
.np_udens_d8_pilot_prepare <- function(bws, x, max.bytes = 256 * 1024^2) {
  state <- .np_udens_d8_prepare(bws, x, max.bytes)
  state$construction <- "density-d8-positive-pilot-v2"
  state
}

.np_udens_d8_pilot_coefficients <- function(q, power) {
  a <- 1 + q/2
  original <- c(a*(a+1)*(a+2)/6, -(a+1)*(a+2)/2, (a+2)/2, -1/6)
  polynomial <- original
  if (power == 2L) {
    polynomial <- numeric(7L)
    for (i in 1:4) for (j in 1:4)
      polynomial[i+j-1L] <- polynomial[i+j-1L] + original[i]*original[j]
  }
  out <- numeric(length(polynomial)); c <- 1/(1+power)
  for (k in 0:(length(polynomial)-1L)) for (j in 0:k)
    out[j+1L] <- out[j+1L] + polynomial[k+1L]*choose(k,j)*
      exp(lgamma(q/2+k)-lgamma(q/2+j))*c^(k-j)/(1+power)^(2*j)
  out
}

.np_udens_d8_polynomial_gaussian <- function(t, coefficient, rate, log.normal) {
  degree <- length(coefficient)-1L
  value <- matrix(0,nrow(t),ncol(t))
  ordinary <- is.finite(t) & t <= 100
  u <- t[ordinary]; p <- rep(coefficient[degree+1L],length(u))
  for (k in rev(seq_len(degree))) p <- coefficient[k]+u*p
  value[ordinary] <- p*exp(log.normal-rate*u)
  far <- is.finite(t) & !ordinary
  u <- t[far]; p <- rep(coefficient[1L],length(u))
  for (k in 2:(degree+1L)) p <- coefficient[k]+p/u
  value[far] <- sign(p)*exp(log.normal-rate*u+degree*log(u)+log(abs(p)))
  if (any(!is.finite(value))) stop("nonfinite D8 pilot moment",call.=FALSE)
  value
}

.np_udens_d8_pilot_moments <- function(state, grid, tile = 32L, donors = 256L) {
  m <- nrow(grid); n <- state$n; q <- state$q
  center <- original <- mu <- second <- numeric(m)
  a <- 1+q/2; co <- c(a*(a+1)*(a+2)/6,-(a+1)*(a+2)/2,(a+2)/2,-1/6)
  c1 <- .np_udens_d8_pilot_coefficients(q,1L)
  c2 <- .np_udens_d8_pilot_coefficients(q,2L)
  log.normal <- -q*log(2*pi)/2-sum(log(state$h))
  for (first in seq.int(1L,m,by=tile)) {
    ii <- seq.int(first,min(m,first+tile-1L))
    for (donor in seq.int(1L,n,by=donors)) {
      jj <- seq.int(donor,min(n,donor+donors-1L));t <- matrix(0,length(ii),length(jj))
      for (j in seq_len(q)) t <- t+(outer(grid[[j]][ii],state$x[[j]][jj],"-")/state$h[j])^2/2
      original[ii] <- original[ii]+rowSums(exp(log.normal-t))
      center[ii] <- center[ii]+rowSums(.np_udens_d8_polynomial_gaussian(t,co,1,log.normal))
      mu[ii] <- mu[ii]+rowSums(.np_udens_d8_polynomial_gaussian(t,c1,.5,log.normal-q*log(2)/2))
      second[ii] <- second[ii]+rowSums(.np_udens_d8_polynomial_gaussian(t,c2,2/3,2*log.normal-q*log(3)/2))
    }
    .np_plot_progress_notify(stage="Constructing complete D8 pilot moments",done=max(ii),total=m,unit="targets")
  }
  mu <- mu/n; variance <- (second/n-mu^2)/n
  if (any(!is.finite(variance)) || any(variance<=0))
    stop("D8 pilot covariance requires positive finite scales",call.=FALSE)
  list(center=center/n,original.center=original/n,mean=mu,se=sqrt(variance))
}

.np_udens_d8_pilot_column <- function(state, grid, mean, column, tile = 32L, donors = 256L) {
  m <- nrow(grid); n <- state$n; q <- state$q; total <- numeric(m)
  log.normal <- -q*log(2*pi)-q*log(3)/2-2*sum(log(state$h))
  for (first in seq.int(1L,m,by=tile)) {
    ii <- seq.int(first,min(m,first+tile-1L))
    for (donor in seq.int(1L,n,by=donors)) {
      jj <- seq.int(donor,min(n,donor+donors-1L))
      A <- B <- D <- exponent <- matrix(0,length(ii),length(jj))
      for (j in seq_len(q)) {
        left <- outer(grid[[j]][ii],state$x[[j]][jj],"-")/state$h[j]
        right <- matrix(rep((grid[[j]][column]-state$x[[j]][jj])/state$h[j],each=length(ii)),length(ii))
        u <- (2*left-right)/3;v <- (2*right-left)/3
        A <- A+u^2;B <- B+v^2;D <- D+u*v
        exponent <- exponent+((left-right)^2+left^2+right^2)/6
      }
      polynomial <- .np_udens_d8_cross_polynomial(q,A,B,D)
      product <- polynomial*exp(log.normal-exponent)
      if (any(!is.finite(product))) stop("nonfinite complete D8 covariance product",call.=FALSE)
      total[ii] <- total[ii]+rowSums(product)
    }
  }
  out <- (total/n-mean*mean[column])/n
  if (any(!is.finite(out))) stop("nonfinite complete D8 covariance column",call.=FALSE)
  out
}

# Complete radial Gaussian-product polynomial; fixed pilot variance H.
# Coefficients are generated exact Gaussian moments, not fitted constants.
.np_udens_d8_cross_polynomial <- function(q, A, B, D) {
  if (q == 1L) return(A*(A*(A*(B*(B*((1.0/2304.0)*B - 1.0/144.0) + 17.0/576.0) - 43.0/1296.0) + B*(B*(-1.0/144.0*B + (1.0/192.0)*D + 137.0/1152.0) - 1.0/18.0*D - 29.0/54.0) + (17.0/144.0)*D + 1631.0/2592.0) + B*(B*((17.0/576.0)*B - 1.0/18.0*D - 29.0/54.0) + D*((1.0/72.0)*D + 5.0/8.0) + 6517.0/2592.0) + D*(-2.0/27.0*D - 112.0/81.0) - 35.0/12.0) + B*(B*(-43.0/1296.0*B + (17.0/144.0)*D + 1631.0/2592.0) + D*(-2.0/27.0*D - 112.0/81.0) - 35.0/12.0) + D*(D*((1.0/162.0)*D + 133.0/324.0) + 665.0/216.0) + 11725.0/3888.0)
  if (q == 2L) return(A*(A*(A*(B*(B*((1.0/2304.0)*B - 1.0/128.0) + 11.0/288.0) - 65.0/1296.0) + B*(B*(-1.0/128.0*B + (1.0/192.0)*D + 43.0/288.0) - 1.0/16.0*D - 37.0/48.0) + (11.0/72.0)*D + 19.0/18.0) + B*(B*((11.0/288.0)*B - 1.0/16.0*D - 37.0/48.0) + D*((1.0/72.0)*D + 85.0/108.0) + 223.0/54.0) + D*(-1.0/12.0*D - 2) - 461.0/81.0) + B*(B*(-65.0/1296.0*B + (11.0/72.0)*D + 19.0/18.0) + D*(-1.0/12.0*D - 2) - 461.0/81.0) + D*(D*((1.0/162.0)*D + 14.0/27.0) + 418.0/81.0) + 5396.0/729.0)
  if (q == 3L) return(A*(A*(A*(B*(B*((1.0/2304.0)*B - 5.0/576.0) + 83.0/1728.0) - 47.0/648.0) + B*(B*(-5.0/576.0*B + (1.0/192.0)*D + 211.0/1152.0) - 5.0/72.0*D - 115.0/108.0) + (83.0/432.0)*D + 1447.0/864.0) + B*(B*((83.0/1728.0)*B - 5.0/72.0*D - 115.0/108.0) + D*((1.0/72.0)*D + 209.0/216.0) + 1849.0/288.0) + D*(-5.0/54.0*D - 25.0/9.0) - 6643.0/648.0) + B*(B*(-47.0/648.0*B + (83.0/432.0)*D + 1447.0/864.0) + D*(-5.0/54.0*D - 25.0/9.0) - 6643.0/648.0) + D*(D*((1.0/162.0)*D + 23.0/36.0) + 1757.0/216.0) + 61915.0/3888.0)
  stop("invalid D8 pilot dimension", call. = FALSE)
}
