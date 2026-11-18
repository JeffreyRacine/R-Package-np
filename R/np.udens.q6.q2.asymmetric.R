# Directional variability is centered; the calibrated corrected root is not.
.np_udens_q6_q2_asymmetric_active <- function(state)
  identical(state[["construction", exact=TRUE]], "density-compact-q6-q2-asymmetric-v1")

.np_udens_q6_q2_asymmetric_prepare <- function(bws,x,seed) {
  state <- .np_udens_q6_q2_prepare(bws,x,seed)
  state$plot.calibration <- "asymmetric-sd"
  state$construction <- "density-compact-q6-q2-asymmetric-v1"
  state
}

.np_udens_q6_q2_directional <- function(root,alpha) {
  if (!is.matrix(root) || !is.numeric(root) || !nrow(root) || ncol(root)<2L ||
      any(!is.finite(root)) || !is.numeric(alpha) || length(alpha)!=1L ||
      !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("invalid q=2 Q6 directional roots",call.=FALSE)
  u <- root-rowMeans(root)
  plus <- sqrt(2*rowMeans(pmax(u,0)^2))
  minus <- sqrt(2*rowMeans(pmax(-u,0)^2))
  if (any(!is.finite(plus)) || any(!is.finite(minus)) || any(plus<=0 | minus<=0))
    stop("q=2 Q6 directional variability must be finite and positive",call.=FALSE)
  z <- pmax(pmax(root,0)/plus,pmax(-root,0)/minus)
  if (any(!is.finite(z))) stop("nonfinite q=2 Q6 directional pivot",call.=FALSE)
  list(plus=plus,minus=minus,
    pointwise=apply(z,1L,quantile,probs=1-alpha,type=7,names=FALSE),
    maximum=apply(z,2L,max))
}

.np_udens_q6_q2_asymmetric_reduce <- function(reference,tasks,alpha) {
  B <- sum(vapply(tasks,function(z) z$bsz,0L)); m <- nrow(reference$grid)
  plus <- minus <- critical <- numeric(m); maximum <- numeric(B); error <- 0
  for (first in seq.int(1L,m,by=reference$workspace$block)) {
    ix <- first:min(first+reference$workspace$block-1L,m)
    part <- reference; part$grid <- reference$grid[ix,,drop=FALSE]
    for (field in c("center","reference.center","se")) part[[field]] <- reference[[field]][ix]
    part$directional <- TRUE
    draws <- .np_udens_q6_q2_dispatch(tasks,part)
    if (!is.matrix(draws) || !identical(dim(draws),c(B,length(ix)+1L)) ||
        any(!is.finite(draws)) || any(draws[,ncol(draws)]<0))
      stop("incomplete q=2 Q6 directional replicate collection",call.=FALSE)
    z <- .np_udens_q6_q2_directional(t(draws[,seq_along(ix),drop=FALSE]),alpha)
    plus[ix] <- z$plus; minus[ix] <- z$minus; critical[ix] <- z$pointwise
    maximum <- pmax(maximum,z$maximum)
    error <- max(error,draws[,ncol(draws)])
    .np_plot_progress_notify(stage="Calibrating directional q=2 Q6 bands",
      done=max(ix),total=m,unit="targets")
  }
  ix <- reference$index
  list(center=reference$center[ix],reference.se=reference$se[ix],
    plus=plus[ix],minus=minus[ix],pointwise=critical[ix],
    simultaneous=unname(quantile(maximum,1-alpha,type=7)),
    maximum=maximum,refinement=error)
}

.np_udens_q6_q2_asymmetric_public <- function(state,grid,alpha,B,band) {
  if (!.np_udens_q6_q2_asymmetric_active(state) ||
      !isTRUE(band %in% c("pointwise","simultaneous")) ||
      !is.data.frame(grid) || !identical(names(grid),names(state$x)) ||
      !all(vapply(grid,function(z) is.numeric(z) && !is.object(z),TRUE)) ||
      !is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("invalid q=2 asymmetric Q6 pointwise/simultaneous request",call.=FALSE)
  B <- .np_plot_gaussian_integer(B,"q=2 asymmetric Q6 B",.Machine$integer.max)
  if (B<2L) stop("q=2 asymmetric Q6 requires at least two conditional draws",call.=FALSE)
  x <- as.matrix(state$x); target <- as.matrix(grid)
  .np_udens_q6_q2_geometry(x,state$h,target,state$kernel)
  .np_udens_q6_q2_budget(nrow(x),nrow(unique(target)),
    vapply(1:2,function(j) length(unique(target[,j])),0L),B)
  inference <- .np_plot_corrected_activity_run("Constructing directional q=2 Q6 density intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_q2_reference(x,state$h,target,state$kernel,B)
      tasks <- .np_udens_q6_q2_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),state$n,B)
      z <- .np_udens_q6_q2_asymmetric_reduce(reference,tasks,alpha)
      make <- function(k) array(c(z$center-k*z$plus,z$center+k*z$minus),
        c(nrow(grid),2L,1L),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      critical <- z[c("pointwise","simultaneous")]
      bands <- lapply(critical[band],make); m <- nrow(reference$grid)
      resolution <- .np_band_resolution(B,rep("density",m),alpha,rep(TRUE,m),band)
      resolution$calibration <- "smoothed-q6 asymmetric-sd"
      resolution$quantile.type <- 7L
      list(center=z$center,se=sqrt((z$plus^2+z$minus^2)/2),
        directional.plus=z$plus,directional.minus=z$minus,reference.se=z$reference.se,
        critical=z$pointwise,bands=bands,intervals=bands,absolute.intervals=bands,
        covariance=NULL,resolution=resolution,maximum=z$maximum,companion.critical=critical,
        checks=c(reference$diagnostics,conditional=z$refinement),
        settings=list(construction=state$construction,law="smoothed-q6",h=state$h,
          n=state$n,q=2L,alpha=alpha,B=B,grid.size=m,family.rows=reference$index,
          family="unique supplied finite grid",available.bands=band,primary="absolute",
          quantile.type=7L,calibration="asymmetric-sd",workspace.bytes=reference$workspace$bytes,
          reference="paired empirical mixture plus independent h_j*sqrt(2/3)*Student6 coordinate noise",
          root.center="full corrected functional T_h(Q6)",
          scale="centered conditional positive/negative semideviations times sqrt(2)",
          marginal.statistic="uncentered complete corrected root divided by its directional scale",
          correction="product order-six compact operator with refitted bivariate normal guide",
          calibration.status="explicit directional Q6; conservative empirical-grid q=2 confirmation"))
    }))
  .np_band_resolution_warn(inference$resolution)
  k <- inference$companion.critical[[band]]
  list(boot.err=cbind(k*inference$directional.plus,k*inference$directional.minus,inference$center),
       boot.all.err=NULL,inference=inference,bxp=list())
}

.np_udens_q6_q2_asymmetric_subset <- function(raw,rows) {
  raw$boot.err <- raw$boot.err[rows,,drop=FALSE]
  for (field in c("center","se","critical","directional.plus","directional.minus","reference.se"))
    raw$inference[[field]] <- raw$inference[[field]][rows]
  for (field in c("bands","intervals","absolute.intervals"))
    raw$inference[[field]] <- lapply(raw$inference[[field]],function(z) z[rows,,,drop=FALSE])
  raw$inference$companion.critical$pointwise <- raw$inference$companion.critical$pointwise[rows]
  raw$inference$settings$family.rows <- raw$inference$settings$family.rows[rows]
  raw$inference$settings$plot.rows <- rows
  raw$inference$settings$family <- "union of all supplied plot slices"
  raw
}
