# Plot-only order-six compact correction with complete normal-guide influence.
# The fitted order-two Epanechnikov density and its selected bandwidth are intact.
.np_udens_m4_request <- function(bws, control, center, method, band,
                                calibration = NULL) {
  inherits(bws, "bandwidth") && is.null(control) && is.null(calibration) &&
    identical(center, "bias-corrected") && identical(method, "asymptotic") &&
    identical(band, "simultaneous") &&
    identical(bws[["ckertype", exact = TRUE]], "epanechnikov") &&
    identical(bws[["type", exact = TRUE]], "fixed") &&
    identical(bws[["ckerbound", exact = TRUE]], "none") &&
    isTRUE(bws[["ckerorder", exact = TRUE]] == 2L) &&
    isTRUE(bws[["ndim", exact = TRUE]] == 1L) &&
    isTRUE(bws[["ncon", exact = TRUE]] == 1L)
}

.np_udens_m4_active <- function(state)
  identical(state[["construction", exact = TRUE]], "density-m4-normal-guide-v1")

.np_udens_m4_prepare <- function(bws, x, seed) {
  if (!.np_udens_m4_request(bws, NULL, "bias-corrected", "asymptotic", "simultaneous") ||
      !is.data.frame(x) || ncol(x) != 1L || nrow(x) < 10L ||
      !nzchar(names(x)[1L]) || !is.numeric(x[[1L]]) || any(!is.finite(x[[1L]])))
    stop("M4 corrected bands require at least ten finite continuous observations", call. = FALSE)
  h <- as.double(.np_physical_bandwidth(bws))
  if (length(h) != 1L || !is.finite(h) || h <= 0)
    stop("M4 corrected bands require a positive fixed bandwidth", call. = FALSE)
  npValidateNonNegativeInteger(seed, "random.seed")
  list(x = x, h = h, n = nrow(x), q = 1L, plot.random.seed = seed,
       plot.law = "gaussian", construction = "density-m4-normal-guide-v1")
}

.np_udens_m4_budget <- function(n, m, block = 128L, B = 0) {
  # Charge both refinement matrices, factor/eigen copies and product batches.
  bytes <- 8 * (12 * as.double(m)^2 + 8 * as.double(m) * block +
    4 * as.double(n) + 8 * as.double(block)^2 + as.double(B))
  if (!is.finite(bytes) || bytes > 192 * 1024^2)
    stop("M4 inference workspace bound", call. = FALSE)
  bytes
}

.np_udens_m4_kernel <- function(r) {
  span <- sqrt(5)
  ifelse(abs(r) < span,
    3 / (4 * span) * (1 - r*r/5) * (2.734375 - 3.28125*r*r + .721875*r^4), 0)
}

.np_udens_m4_moments <- function(x, h, grid, order = 48L,
                                block = 128L, progress = function(...) NULL) {
  n <- length(x); m <- length(grid)
  if (!is.numeric(x) || !is.numeric(grid) || n < 10L || m < 1L ||
      length(h) != 1L || !is.finite(h) || h <= 0 ||
      any(!is.finite(c(x, grid))) || !order %in% c(24L, 48L) ||
      length(block) != 1L || !is.finite(block) || block < 1 || block != floor(block))
    stop("invalid M4 geometry", call. = FALSE)
  forecast <- .np_udens_m4_budget(n, m, block)
  mu <- mean(x); v <- mean((x-mu)^2); vq <- v + h*h
  if (!is.finite(v) || v <= 0 || !is.finite(vq) || vq <= 0)
    stop("invalid M4 normal guide", call. = FALSE)
  j <- seq_len(order-1L); J <- matrix(0, order, order)
  J[cbind(j,j+1L)] <- j/sqrt(4*j*j-1); J <- J + t(J)
  rule <- eigen(J, symmetric = TRUE); ix <- order(rule$values)
  nodes <- rule$values[ix]; weights <- 2*rule$vectors[1L,ix]^2
  cuts <- sort(unique(c(grid-sqrt(5)*h, grid+sqrt(5)*h)))
  if (any(!is.finite(cuts)) || length(cuts) < 2L)
    stop("unrepresentable M4 support partition", call. = FALSE)
  ma <- original.guide <- numeric(m)
  CAA <- matrix(0,m,m); X <- Dint <- matrix(0,m,2L)
  node.count <- 0L
  for (i in seq_len(length(cuts)-1L)) {
    left <- cuts[i]; right <- cuts[i+1L]; mid <- (left+right)/2
    if (!any(abs(mid-grid) < sqrt(5)*h)) next
    u <- mid + (right-left)/2*nodes; w <- (right-left)/2*weights
    A <- .np_udens_m4_kernel(outer(grid,u,function(t,z)(z-t)/h))/h
    Q <- numeric(length(u))
    for (first in seq.int(1L,n,by=block)) {
      donors <- first:min(n,first+block-1L)
      Q <- Q + rowSums(dnorm(outer(u,x[donors],"-")/h))/h/n
    }
    p <- w*Q; ma <- ma + drop(A %*% p)
    CAA <- CAA + tcrossprod(sweep(A,2,sqrt(p),"*"))
    IF <- cbind(u-mu, ((u-mu)^2/vq-1)/2)
    X <- X + A %*% (IF*p)
    pq <- dnorm(u,mu,sqrt(vq))
    Dint <- Dint + A %*% (cbind((u-mu)/vq,(u-mu)^2/vq-1)*(w*pq))
    original.guide <- original.guide + drop(A %*% (w*dnorm(u,mu,sqrt(v))))
    node.count <- node.count + length(u)
    progress("complete covariance",i,length(cuts)-1L)
  }
  CAA <- CAA - tcrossprod(ma)
  gq <- dnorm(grid,mu,sqrt(vq))
  D <- cbind(gq*(grid-mu)/vq,gq*((grid-mu)^2/vq-1)) - Dint
  m3 <- mean((x-mu)^3); m4 <- mean((x-mu)^4)+6*h*h*v+3*h^4
  V <- matrix(c(vq,m3/(2*vq),m3/(2*vq),(m4/vq^2-1)/4),2L)
  C <- (CAA + X %*% t(D) + D %*% t(X) + D %*% V %*% t(D))/n
  observed <- numeric(m)
  for (first in seq.int(1L,n,by=block)) {
    donors <- first:min(n,first+block-1L)
    observed <- observed + rowSums(.np_udens_m4_kernel(
      outer(grid,x[donors],function(t,z)(z-t)/h))/h)/n
  }
  center <- observed-original.guide+dnorm(grid,mu,sqrt(v))
  if (any(!is.finite(c(center,C))) || any(diag(C) <= 0))
    stop("invalid M4 moments", call. = FALSE)
  list(C=C, center=center, nodes=node.count, mu=mu, v=v, vq=vq,
       ma=ma, D=D, V=V, X=X, forecast.bytes=forecast)
}

.np_udens_m4_inference <- function(state, grid, alpha, B, progress) {
  original <- grid[[1L]]; targets <- sort(unique(original))
  restore <- match(original,targets); m <- length(targets)
  B <- .np_plot_gaussian_integer(B,"density draws",.Machine$integer.max)
  .np_udens_m4_budget(state$n,m,B=B)
  .np_plot_gaussian_draws(B,alpha,m)
  low <- .np_udens_m4_moments(state$x[[1L]],state$h,targets,24L,progress=progress)
  high <- .np_udens_m4_moments(state$x[[1L]],state$h,targets,48L,progress=progress)
  covariance.error <- max(abs(low$C-high$C)); covariance.limit <- 1e-12+1e-9*max(abs(high$C))
  center.error <- max(abs(low$center-high$center))
  if (covariance.error > covariance.limit || center.error > 1e-9)
    stop("M4 quadrature refinement failed", call. = FALSE)
  rm(low)
  eig <- eigen(high$C,symmetric=TRUE)
  if (max(eig$values) <= 0 || min(eig$values) < -1e-10*max(eig$values))
    stop("M4 covariance is not positive semidefinite", call. = FALSE)
  # The frozen method retains every nonnegative eigenvalue, including near-null ones.
  factor <- sweep(eig$vectors,2,sqrt(pmax(eig$values,0)),"*")
  if (max(abs(tcrossprod(factor)-high$C)) > covariance.limit)
    stop("M4 covariance factor reconstruction failed", call. = FALSE)
  maxima <- numeric(B)
  for (first in seq.int(1L,B,by=128L)) {
    ids <- first:min(B,first+127L)
    z <- matrix(rnorm(m*length(ids)),m)
    values <- factor %*% z
    if (any(!is.finite(values))) stop("nonfinite M4 Gaussian process",call.=FALSE)
    maxima[ids] <- apply(abs(values),2L,max)
    progress("Gaussian maximum",max(ids),B)
  }
  critical <- unname(quantile(maxima,1-alpha,type=7L))
  center <- high$center[restore]; scale <- sqrt(diag(high$C))[restore]
  limits <- array(c(center-critical,center+critical),c(length(original),2L,1L),
    dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
  bands <- list(simultaneous=limits)
  list(center=center,se=scale,bands=bands,intervals=bands,absolute.intervals=bands,
       maximum=maxima,critical=critical,covariance=NULL,
       checks=list(covariance.error=covariance.error,covariance.limit=covariance.limit,
                   center.error=center.error,minimum.eigenvalue=min(eig$values)),
       settings=list(construction=state$construction,law="gaussian",h=state$h,n=state$n,q=1L,
         alpha=alpha,B=B,grid.size=m,family="unique supplied finite grid",
         family.rows=restore,available.bands="simultaneous",primary="absolute",
         quantile.type=7L,maximum.statistic="unstudentized absolute Gaussian maximum",
         covariance="complete corrected-curve influence under equal-weight Gaussian mixture pilot",
         correction="order-six compact operator with fitted normal guide",
         calibration.status="explicit simultaneous-only M4 density plot inference"))
}

.np_udens_m4_public_errors <- function(state, grid, alpha, B, band) {
  if (!.np_udens_m4_active(state) || !identical(band,"simultaneous") ||
      !is.data.frame(grid) || ncol(grid) != 1L || !nrow(grid) ||
      !identical(names(grid),names(state$x)) || !is.numeric(grid[[1L]]) ||
      any(!is.finite(grid[[1L]])) || !is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("invalid M4 simultaneous-band request",call.=FALSE)
  progress <- function(stage,done,total).np_plot_progress_notify(
    stage=paste("Compact density M4:",stage),done=done,total=total,unit="work")
  inference <- .np_plot_corrected_activity_run("Constructing compact density simultaneous intervals",
    .np_with_seed(state$plot.random.seed,
      .np_udens_m4_inference(state,grid,alpha,B,progress)))
  list(boot.err=cbind(inference$critical,inference$critical,inference$center),
       boot.all.err=NULL,inference=inference,bxp=list())
}
