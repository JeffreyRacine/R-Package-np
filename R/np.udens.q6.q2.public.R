.np_udens_q6_q2_subset <- function(raw,rows) {
  raw$boot.err <- raw$boot.err[rows,,drop=FALSE]
  if (!is.null(raw$boot.all.err))
    raw$boot.all.err <- lapply(raw$boot.all.err,function(z) z[rows,,drop=FALSE])
  for (field in c("center","se","critical")) raw$inference[[field]] <- raw$inference[[field]][rows]
  for (field in c("bands","intervals","absolute.intervals"))
    raw$inference[[field]] <- lapply(raw$inference[[field]],function(z) z[rows,,,drop=FALSE])
  raw$inference$raw.bonferroni <- raw$inference$raw.bonferroni[rows,,,drop=FALSE]
  for (field in c("pointwise","bonferroni","bonferroni.raw"))
    raw$inference$companion.critical[[field]] <- raw$inference$companion.critical[[field]][rows]
  raw$inference$settings$family.rows <- raw$inference$settings$family.rows[rows]
  raw$inference$settings$plot.rows <- rows
  raw$inference$settings$family <- "union of all supplied plot slices"
  raw
}
# One paired process supplies all four symmetric band summaries. Targets are
# blocked only for storage; each block replays the identical conditional samples.
.np_udens_q6_q2_reduce <- function(reference,tasks,alpha) {
  B <- sum(vapply(tasks,function(z) z$bsz,0L)); m <- nrow(reference$grid)
  critical <- bonferroni <- numeric(m); maximum <- numeric(B); error <- 0
  for (first in seq.int(1L,m,by=reference$workspace$block)) {
    ix <- first:min(first+reference$workspace$block-1L,m)
    part <- reference; part$grid <- reference$grid[ix,,drop=FALSE]
    for (field in c("center","reference.center","se")) part[[field]] <- reference[[field]][ix]
    draws <- .np_udens_q6_q2_dispatch(tasks,part)
    if (!is.matrix(draws) || !identical(dim(draws),c(B,length(ix)+1L)) ||
        any(!is.finite(draws)) || any(draws<0))
      stop("incomplete q=2 Q6 replicate collection",call.=FALSE)
    for (j in seq_along(ix)) {
      critical[ix[j]] <- unname(quantile(draws[,j],1-alpha,type=7))
      bonferroni[ix[j]] <- unname(quantile(draws[,j],1-alpha/m,type=7))
    }
    maximum <- pmax(maximum,apply(draws[,seq_along(ix),drop=FALSE],1L,max))
    error <- max(error,draws[,ncol(draws)])
    .np_plot_progress_notify(stage="Calibrating paired q=2 Q6 bands",
      done=max(ix),total=m,unit="targets")
  }
  simultaneous <- unname(quantile(maximum,1-alpha,type=7)); ix <- reference$index
  list(center=reference$center[ix],se=reference$se[ix],pointwise=critical[ix],
    simultaneous=simultaneous,bonferroni.raw=bonferroni[ix],
    bonferroni=pmax(bonferroni[ix],simultaneous),maximum=maximum,refinement=error)
}

.np_udens_q6_q2_public <- function(state,grid,alpha,B,band) {
  if (!.np_udens_q6_q2_active(state) || !isTRUE(band %in% c("pointwise","simultaneous","bonferroni","all")) ||
      !is.data.frame(grid) || !identical(names(grid),names(state$x)) ||
      !all(vapply(grid,function(z) is.numeric(z) && !is.object(z),TRUE)) ||
      !is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5)
    stop("invalid q=2 Q6 band request",call.=FALSE)
  B <- .np_plot_gaussian_integer(B,"q=2 Q6 B",.Machine$integer.max)
  x <- as.matrix(state$x); target <- as.matrix(grid)
  .np_udens_q6_q2_geometry(x,state$h,target,state$kernel)
  .np_udens_q6_q2_budget(nrow(x),nrow(unique(target)),
    vapply(1:2,function(j) length(unique(target[,j])),0L),B)
  inference <- .np_plot_corrected_activity_run("Constructing paired q=2 Q6 density intervals",
    .np_with_seed(state$plot.random.seed, {
      reference <- .np_udens_q6_q2_reference(x,state$h,target,state$kernel,B)
      tasks <- .np_udens_q6_q2_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),state$n,B)
      z <- .np_udens_q6_q2_reduce(reference,tasks,alpha)
      make <- function(k) {
        width <- k*z$se
        array(c(z$center-width,z$center+width),c(nrow(grid),2L,1L),
          dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      }
      requested <- if (band=="all") c("pointwise","simultaneous","bonferroni") else band
      critical <- z[c("pointwise","simultaneous","bonferroni","bonferroni.raw")]
      bands <- lapply(critical[requested],make)
      m <- nrow(reference$grid)
      resolution <- .np_band_resolution(B,rep("density",m),alpha,rep(TRUE,m),requested)
      resolution$calibration <- "smoothed-q6 symmetric-sd"
      resolution$quantile.type <- 7L
      resolution$interpolated.index <- 1+(B-1)*(1-alpha/ifelse(resolution$band=="bonferroni",m,1L))
      list(center=z$center,se=z$se,critical=z$pointwise,bands=bands,intervals=bands,
        absolute.intervals=bands,covariance=NULL,resolution=resolution,
        maximum=z$maximum,raw.bonferroni=make(z$bonferroni.raw),
        companion.critical=critical,checks=c(reference$diagnostics,conditional=z$refinement),
        settings=list(construction=state$construction,law="smoothed-q6",h=state$h,
          n=state$n,q=2L,alpha=alpha,B=B,grid.size=m,family.rows=reference$index,
          family="unique supplied finite grid",available.bands=requested,primary="absolute",
          quantile.type=7L,calibration="symmetric-sd",workspace.bytes=reference$workspace$bytes,
          bonferroni.construction="type-7 marginal absolute quantile at alpha/m with simultaneous hull; raw endpoints retained",
          marginal.statistic="absolute complete corrected bootstrap root divided by fixed original scale",
          reference="paired empirical mixture plus independent h_j*sqrt(2/3)*Student6 coordinate noise",
          root.center="full corrected functional T_h(Q6)",
          covariance="complete five-direction bivariate guide influence and cross terms under Q6",
          correction="product order-six compact operator with refitted bivariate normal guide",
          calibration.status="explicit symmetric Q6 density plot inference; conservative q=2 confirmation"))
    }))
  .np_band_resolution_warn(inference$resolution)
  widths <- lapply(inference$companion.critical[names(inference$bands)],
    function(k) { width <- k*inference$se; cbind(width,width) })
  list(boot.err=cbind(widths[[if(band=="all")"pointwise" else band]],inference$center),
       boot.all.err=if(band=="all") widths else NULL,inference=inference,bxp=list())
}
