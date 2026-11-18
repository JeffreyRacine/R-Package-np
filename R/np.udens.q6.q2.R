# Bivariate Q6 plot inference. Common empirical donors preserve dependence;
# the normal guide has two mean and three covariance influence directions.
.np_udens_q6_q2_spec <- function(bws) {
  if (!inherits(bws, "bandwidth") ||
      !identical(bws[["type", exact=TRUE]], "fixed") ||
      !identical(bws[["ckerbound", exact=TRUE]], "none") ||
      !isTRUE(bws[["ckertype", exact=TRUE]] %in% c("uniform", "epanechnikov")) ||
      !isTRUE(bws[["ckerorder", exact=TRUE]] == 2L) ||
      !isTRUE(bws[["ndim", exact=TRUE]] == 2L) ||
      !isTRUE(bws[["ncon", exact=TRUE]] == 2L))
    stop("q=2 symmetric smoothed-q6 requires two continuous, whole-line, fixed-bandwidth order-two uniform or Epanechnikov coordinates", call.=FALSE)
  invisible(NULL)
}

.np_udens_q6_public_spec <- function(bws, calibration, band=NULL) {
  if (isTRUE(bws[["ndim", exact=TRUE]] == 2L) &&
      isTRUE(calibration %in% c("symmetric-sd", "asymmetric-sd"))) {
    .np_udens_q6_q2_spec(bws)
    if (identical(calibration, "asymmetric-sd") &&
        !isTRUE(band %in% c("pointwise", "simultaneous")))
      stop("q=2 asymmetric smoothed-q6 requires pointwise or simultaneous bands", call.=FALSE)
  } else .np_udens_q6_spec(bws)
}

.np_udens_q6_q2_geometry <- function(x, h, grid, kernel) {
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 2L || nrow(x) < 10L ||
      any(!is.finite(x)) || !is.numeric(h) || length(h) != 2L ||
      any(!is.finite(h) | h <= 0) || !is.matrix(grid) ||
      !is.numeric(grid) || ncol(grid) != 2L || !nrow(grid) ||
      any(!is.finite(grid)) || !isTRUE(kernel %in% c("uniform", "epanechnikov")))
    stop("invalid q=2 Q6 geometry", call.=FALSE)
  S <- crossprod(sweep(x, 2L, colMeans(x))) / nrow(x)
  if (any(!is.finite(S)) || !is.finite(det(S)) || det(S) <= 0)
    stop("q=2 Q6 requires a positive definite bivariate guide", call.=FALSE)
  invisible(NULL)
}

.np_udens_q6_q2_prepare <- function(bws, x, seed) {
  .np_udens_q6_q2_spec(bws)
  if (!is.data.frame(x) || ncol(x) != 2L || any(!nzchar(names(x))) ||
      anyDuplicated(names(x)) || !all(vapply(x, function(z)
        is.numeric(z) && !is.object(z), TRUE)))
    stop("q=2 Q6 requires two ordinary numeric continuous coordinates", call.=FALSE)
  h <- as.double(.np_physical_bandwidth(bws))
  .np_udens_q6_q2_geometry(as.matrix(x), h, as.matrix(x[1L,,drop=FALSE]), bws$ckertype)
  npValidateNonNegativeInteger(seed, "random.seed")
  list(x=x, h=h, n=nrow(x), q=2L, kernel=bws$ckertype,
    plot.random.seed=seed, plot.law="smoothed-q6", plot.calibration="symmetric-sd",
    construction="density-compact-q6-q2-symmetric-v1")
}

.np_udens_q6_q2_active <- function(state)
  isTRUE(state[["construction", exact=TRUE]] %in%
    c("density-compact-q6-q2-symmetric-v1", "density-compact-q6-q2-asymmetric-v1"))

.np_udens_q6_q2_budget <- function(n, m, axes, B, chunk=32L) {
  # Account for paired axis caches, task seeds, samples, working target blocks,
  # all returned roots in a block, and both refinement levels.
  if (any(!is.finite(c(n,m,axes,B,chunk))) || any(c(n,m,axes,B,chunk) < 1) ||
      any(c(n,m,axes,B,chunk) != floor(c(n,m,axes,B,chunk))))
    stop("invalid q=2 Q6 workspace request", call.=FALSE)
  n <- as.double(n); m <- as.double(m); B <- as.double(B)
  base <- 8 * (12*n*sum(axes) + 20*n*chunk + 40*m +
                10*626*ceiling(B/chunk))
  per.target <- 8 * (24*n + 32*chunk + 8*B)
  block <- min(m, 64, floor((192*1024^2-base)/per.target))
  if (!is.finite(block) || block < 1 || n*chunk > .Machine$integer.max)
    stop("q=2 Q6 inference exceeds its 192 MiB workspace bound", call.=FALSE)
  list(block=as.integer(block), bytes=base+block*per.target)
}

.np_udens_q6_q2_moments <- function(x, h) {
  mu <- colMeans(x); z <- sweep(x, 2L, mu)
  noise <- function(j, k) if (k==0L) 1 else if (k==2L) h[j]^2 else if (k==4L) 6*h[j]^4 else 0
  moment <- function(i, j) {
    ans <- 0
    for (a in 0:i) for (b in 0:j)
      ans <- ans + choose(i,a)*choose(j,b)*mean(z[,1L]^(i-a)*z[,2L]^(j-b))*noise(1L,a)*noise(2L,b)
    ans
  }
  S <- crossprod(z)/nrow(x) + diag(h^2)
  powers <- rbind(c(1L,0L),c(0L,1L),c(2L,0L),c(1L,1L),c(0L,2L))
  means <- c(0,0,S[1L,1L],S[1L,2L],S[2L,2L]); V <- matrix(0,5L,5L)
  for (i in 1:5) for (j in 1:5)
    V[i,j] <- moment(sum(powers[c(i,j),1L]),sum(powers[c(i,j),2L])) - means[i]*means[j]
  list(mu=mu, S=S, V=V)
}

.np_udens_q6_q2_normal <- function(grid, mu, S, score=FALSE) {
  if (any(!is.finite(S)) || !is.finite(det(S)) || det(S) <= 0)
    stop("q=2 Q6 guide covariance is not positive definite", call.=FALSE)
  z <- sweep(grid, 2L, mu); inv <- solve(S); v <- z %*% inv
  density <- exp(-rowSums(z*v)/2)/(2*pi*sqrt(det(S)))
  if (!score) return(density)
  cbind(v[,1L],v[,2L],(v[,1L]^2-inv[1L,1L])/2,
        v[,1L]*v[,2L]-inv[1L,2L],(v[,2L]^2-inv[2L,2L])/2)*density
}

.np_udens_q6_q2_rule <- function(kernel, order) {
  rule <- .np_udens_q6_rule(order)
  span <- if (kernel=="uniform") 1 else sqrt(5)
  u <- span*rule$x; w <- span*rule$w
  list(u=u, w=w, L=.np_udens_q6_kernel(u,kernel))
}

.np_udens_q6_q2_axis <- function(target, x, h, mu, rule) {
  out <- lapply(seq_len(4L),function(i) matrix(0,length(target),length(x)))
  for (a in seq_along(target)) {
    v <- target[a]+h*rule$u
    qp <- dt(outer(v,x,"-")/(h*sqrt(2/3)),6)/(h*sqrt(2/3))
    for (j in 0:2) out[[j+1L]][a,] <- drop((rule$w*rule$L*(v-mu)^j) %*% qp)
    out[[4L]][a,] <- drop((rule$w*rule$L^2/h) %*% qp)
  }
  names(out) <- c("M0","M1","M2","J")
  out
}

.np_udens_q6_q2_convolution <- function(grid, mu, S, h, rule, score=FALSE) {
  out <- if (score) matrix(0,nrow(grid),5L) else numeric(nrow(grid))
  w <- rule$w*rule$L
  for (i in seq_along(w)) for (j in seq_along(w))
    out <- out + w[i]*w[j]*.np_udens_q6_q2_normal(
      sweep(grid,2L,h*c(rule$u[i],rule$u[j]),"+"),mu,S,score)
  out
}

.np_udens_q6_q2_level <- function(x, h, grid, kernel, order, block) {
  rule <- .np_udens_q6_q2_rule(kernel,order); moment <- .np_udens_q6_q2_moments(x,h)
  axes <- lapply(1:2,function(j) unique(grid[,j]))
  index <- lapply(1:2,function(j) match(grid[,j],axes[[j]]))
  cache <- lapply(1:2,function(j) .np_udens_q6_q2_axis(axes[[j]],x[,j],h[j],moment$mu[j],rule))
  out <- list(variance=numeric(nrow(grid)),center=numeric(nrow(grid)),
              functional=numeric(nrow(grid)),ma=numeric(nrow(grid)),
              X=matrix(0,nrow(grid),5L),D=matrix(0,nrow(grid),5L))
  observed <- crossprod(sweep(x,2L,moment$mu))/nrow(x)
  for (first in seq.int(1L,nrow(grid),by=block)) {
    ix <- first:min(first+block-1L,nrow(grid)); target <- grid[ix,,drop=FALSE]
    A <- lapply(cache[[1L]],function(z) z[index[[1L]][ix],,drop=FALSE])
    C <- lapply(cache[[2L]],function(z) z[index[[2L]][ix],,drop=FALSE])
    ma <- rowMeans(A$M0*C$M0); va <- rowMeans(A$J*C$J)-ma^2
    X <- cbind(rowMeans(A$M1*C$M0),rowMeans(A$M0*C$M1),
      rowMeans(A$M2*C$M0)-moment$S[1L,1L]*ma,
      rowMeans(A$M1*C$M1)-moment$S[1L,2L]*ma,
      rowMeans(A$M0*C$M2)-moment$S[2L,2L]*ma)
    D <- .np_udens_q6_q2_normal(target,moment$mu,moment$S,TRUE) -
      .np_udens_q6_q2_convolution(target,moment$mu,moment$S,h,rule,TRUE)
    out$variance[ix] <- (va+2*rowSums(X*D)+rowSums((D%*%moment$V)*D))/nrow(x)
    raw <- rowMeans(.np_udens_q6_kernel(outer(target[,1L],x[,1L],"-")/h[1L],kernel)*
                    .np_udens_q6_kernel(outer(target[,2L],x[,2L],"-")/h[2L],kernel))/prod(h)
    out$center[ix] <- raw - .np_udens_q6_q2_convolution(target,moment$mu,observed,h,rule) +
      .np_udens_q6_q2_normal(target,moment$mu,observed)
    out$functional[ix] <- ma - .np_udens_q6_q2_convolution(target,moment$mu,moment$S,h,rule) +
      .np_udens_q6_q2_normal(target,moment$mu,moment$S)
    out$ma[ix] <- ma; out$X[ix,] <- X; out$D[ix,] <- D
    .np_plot_progress_notify(stage="Constructing q=2 Q6 reference",
      done=max(ix),total=nrow(grid),unit="targets")
  }
  out
}

.np_udens_q6_q2_reference <- function(x,h,grid,kernel,B) {
  .np_udens_q6_q2_geometry(x,h,grid,kernel)
  # Sort unique pairs, preserving a separate mapping to every supplied row.
  unique.grid <- unique(grid)
  unique.grid <- unique.grid[order(unique.grid[,1L],unique.grid[,2L]),,drop=FALSE]
  # R's coordinate equality treats signed zero as the same target. Pair
  # integer coordinate IDs rather than formatting floating-point bit patterns.
  key <- function(z) paste(match(z[,1L],unique.grid[,1L]),
                           match(z[,2L],unique.grid[,2L]),sep=":")
  index <- match(key(grid),key(unique.grid))
  axes <- vapply(1:2,function(j) length(unique(unique.grid[,j])),0L)
  budget <- .np_udens_q6_q2_budget(nrow(x),nrow(unique.grid),axes,B)
  low <- .np_udens_q6_q2_level(x,h,unique.grid,kernel,24L,budget$block)
  high <- .np_udens_q6_q2_level(x,h,unique.grid,kernel,32L,budget$block)
  checks <- vapply(names(high),function(field)
    .np_udens_q6_check(low[[field]],high[[field]],paste("q=2 original",field)),0.0)
  if (any(!is.finite(high$variance)) || any(high$variance<=0))
    stop("q=2 Q6 complete marginal variance is not positive",call.=FALSE)
  list(x=x,h=h,kernel=kernel,grid=unique.grid,index=index,center=high$center,
    reference.center=high$functional,se=sqrt(high$variance),
    diagnostics=checks,workspace=budget)
}

.np_udens_q6_q2_curve <- function(Y,grid,h,kernel,order) {
  a <- Y[[1L]]; b <- Y[[2L]]
  if (!is.matrix(a) || !identical(dim(a),dim(b)))
    stop("invalid paired q=2 Q6 sample",call.=FALSE)
  mu1 <- colMeans(a); mu2 <- colMeans(b)
  z1 <- sweep(a,2L,mu1); z2 <- sweep(b,2L,mu2)
  s11 <- colMeans(z1*z1); s12 <- colMeans(z1*z2); s22 <- colMeans(z2*z2)
  determinant <- s11*s22-s12*s12
  if (any(!is.finite(determinant)) || any(determinant<=0))
    stop("q=2 Q6 conditional guide covariance is not positive definite",call.=FALSE)
  t1 <- matrix(grid[,1L],nrow(grid),ncol(a)); t2 <- matrix(grid[,2L],nrow(grid),ncol(a))
  density <- function(t1,t2) {
    v1 <- sweep(t1,2L,mu1); v2 <- sweep(t2,2L,mu2)
    numerator <- sweep(v1*v1,2L,s22,"*")-2*sweep(v1*v2,2L,s12,"*")+sweep(v2*v2,2L,s11,"*")
    sweep(exp(-sweep(numerator,2L,2*determinant,"/")),2L,2*pi*sqrt(determinant),"/")
  }
  rule <- .np_udens_q6_q2_rule(kernel,order); w <- rule$w*rule$L
  convolution <- matrix(0,nrow(grid),ncol(a))
  for (i in seq_along(w)) for (j in seq_along(w))
    convolution <- convolution+w[i]*w[j]*density(t1+h[1L]*rule$u[i],t2+h[2L]*rule$u[j])
  raw <- matrix(0,nrow(grid),ncol(a))
  for (j in seq_len(ncol(a)))
    raw[,j] <- rowMeans(.np_udens_q6_kernel(outer(grid[,1L],a[,j],"-")/h[1L],kernel)*
                       .np_udens_q6_kernel(outer(grid[,2L],b[,j],"-")/h[2L],kernel))/prod(h)
  raw-convolution+density(t1,t2)
}
