# Explicit positive beta-binomial smoothing law for corrected density plots.
# It retains the complete boundary operator, original h and fixed target grid.
.np_udens_bb_budget <- function(n,m,B=1L,chunk=16L) {
  d <- as.double(c(n,m,B,chunk))
  if(any(!is.finite(d)) || any(d!=floor(d)) || n<10 || m<1 || B<1 || chunk<1 ||
     any(d>.Machine$integer.max) || n*as.double(chunk)>1e6)
    stop("invalid beta-binomial workspace dimensions",call.=FALSE)
  bytes <- .np_udens_u5_budget(n,m,block=8L)+
    8*(6*as.double(B)*m+12*as.double(n)*chunk+8000*ceiling(B/chunk)+192000*6+30*m)
  if(!is.finite(bytes) || bytes>192*1024^2)
    stop("beta-binomial inference exceeds its 192 MiB workspace bound",call.=FALSE)
  bytes
}

.np_udens_bb_curve <- function(plan,x,lx=log(x),lq=log1p(-x),fit=NULL) {
  if(is.null(fit)) fit <- .np_udens_bb_fit(x,lx,lq)
  ab <- fit$theta+1
  low <- .np_udens_bb_guide(plan,ab,80,32L)
  high <- .np_udens_bb_guide(plan,ab,100,48L)
  refinement <- .np_udens_u5_check(low$value,high$value,"beta-binomial complete guide expectation") /
    (1e-12+1e-9*max(abs(high$value)))
  observed <- .Call(C_np_plot_udens_bb_guide,plan$grid,plan$h,
                   as.double(lx),as.double(lq),rep(1/length(x),length(x)),plan$operator)
  center <- observed-high$value+dbeta(plan$grid,ab[1],ab[2])
  if(any(!is.finite(center))) stop("nonfinite beta-binomial corrected curve",call.=FALSE)
  list(center=center,guide=ab,refinement=refinement)
}

.np_udens_bb_reference <- function(x,h,grid,progress=function(...)NULL) {
  if(!is.numeric(x) || is.object(x) || !is.null(dim(x)) || any(!is.finite(x)) ||
     any(x<=0|x>=1)) stop("beta-binomial observations must be strictly interior numeric values",call.=FALSE)
  .np_udens_bb_budget(length(x),length(grid))
  plan <- .np_udens_bb_guide_plan(sort(unique(grid)),h);targets <- plan$grid
  fit <- .np_udens_beta_q1_fit(x);ab <- fit$theta+1
  if(max(abs(digamma(ab)-digamma(sum(ab))-c(mean(log(x)),mean(log1p(-x)))))>=2e-8)
    stop("beta-binomial original guide score exceeds its frozen tolerance",call.=FALSE)
  n <- length(x);k <- ceiling(2*n^(2/5));j <- 0:k
  A <- ab[1]+j;BB <- ab[2]+k-j
  weights <- rowMeans(vapply(x,function(t) dbinom(j,k,t),numeric(k+1L)))
  if(any(!is.finite(weights)) || any(weights<0)) stop("invalid beta-binomial mixture weights",call.=FALSE)
  .np_udens_bb_check(sum(weights),1,"reference mass")
  log.moments <- .np_udens_bb_log_moments(ab,k,weights)
  logQ <- attr(log.moments,"mean",exact=TRUE)
  fq <- .np_udens_u5_fitmom(logQ,ab)
  z <- list(k=k,a=ab[1],b=ab[2],aQ=fq$ab[1],bQ=fq$ab[2],h=h,bb.log.moments=log.moments)
  original <- .np_udens_bb_curve(plan,x,fit=fit)
  reference.center <- se <- numeric(length(targets));checks <- list()
  for(first in seq.int(1L,length(targets),by=8L)) {
    ix <- first:min(first+7L,length(targets));g <- targets[ix]
    low <- .np_udens_bb_covariance(z,x,g,64L)
    high <- .np_udens_bb_covariance(z,x,g,96L)
    check <- vapply(c("C","center","EA","B","component_means"),function(key)
      .np_udens_bb_check(low[[key]],high[[key]],paste("reference",key)),0.)
    if(any(!is.finite(diag(high$C))) || any(diag(high$C)<=0))
      stop("beta-binomial complete variance must be positive",call.=FALSE)
    reference.center[ix] <- high$center;se[ix] <- sqrt(diag(high$C))
    checks[[length(checks)+1L]] <- c(check,high$positive.checks$refinement,
      influence.mean=max(abs(high$psi_mean)),relative.correction=max(high$positive.checks$correction))
    progress("complete original variance",max(ix),length(targets))
  }
  list(x=x,n=n,h=h,grid=targets,index=match(grid,targets),plan=plan,center=original$center,
    reference.center=reference.center,se=se,weights=weights,A=A,BB=BB,
    original.guide=ab,reference.guide=fq$ab,degree=k,checks=checks,
    original.refinement=original$refinement)
}

.np_udens_bb_active <- function(state)
  identical(state[["construction",exact=TRUE]],"density-beta-binomial-symmetric-v1")

.np_udens_bb_spec <- function(bws) {
  if(!inherits(bws,"bandwidth") || !identical(bws[["type",exact=TRUE]],"fixed") ||
     !identical(bws[["ckerbound",exact=TRUE]],"fixed") ||
     !identical(bws[["ckertype",exact=TRUE]],"epanechnikov") ||
     !isTRUE(bws[["ckerorder",exact=TRUE]]==4L) ||
     !isTRUE(bws[["ndim",exact=TRUE]]==1L) || !isTRUE(bws[["ncon",exact=TRUE]]==1L))
    stop("smoothed-beta-binomial requires q=1 continuous fixed-bandwidth order-four Epanechnikov density with finite known fixed bounds",call.=FALSE)
  lower <- bws[["ckerlb",exact=TRUE]];upper <- bws[["ckerub",exact=TRUE]]
  if(!is.numeric(lower) || !is.numeric(upper) || length(lower)!=1L || length(upper)!=1L ||
     any(!is.finite(c(lower,upper))) || upper<=lower)
    stop("smoothed-beta-binomial requires finite ordered fixed bounds",call.=FALSE)
  invisible(NULL)
}

.np_udens_bb_prepare <- function(bws,x,seed) {
  .np_udens_bb_spec(bws)
  # Reuse only the existing validated unit-support geometry, not its owner.
  out <- .np_udens_u5_prepare(bws,x,seed)
  out$normal.pmzsd <- NULL;out$plot.law <- "smoothed-beta-binomial"
  out$plot.calibration <- "symmetric-sd";out$construction <- "density-beta-binomial-symmetric-v1"
  .np_udens_bb_budget(out$n,1L)
  out
}

.np_udens_bb_reduce <- function(reference,draws,alpha,band) {
  B <- nrow(draws);m <- length(reference$grid)
  if(!is.matrix(draws) || !identical(dim(draws),c(B,m+1L)) ||
     !is.numeric(draws) || B<1L || any(!is.finite(draws)) || any(draws[,m+1L]<0) ||
     !is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5 ||
     !isTRUE(band %in% c("pointwise","simultaneous","bonferroni","all")))
    stop("incomplete beta-binomial root collection",call.=FALSE)
  roots <- abs(draws[,seq_len(m),drop=FALSE])
  pointwise <- vapply(seq_len(m),function(j)unname(quantile(roots[,j],1-alpha,type=7)),0.)
  maximum <- apply(roots,1L,max);simultaneous <- unname(quantile(maximum,1-alpha,type=7))
  bonferroni.raw <- vapply(seq_len(m),function(j)unname(quantile(roots[,j],1-alpha/m,type=7)),0.)
  # This explicit companion hull follows the existing SD-band convention.
  bonferroni <- pmax(bonferroni.raw,simultaneous)
  list(pointwise=pointwise,simultaneous=simultaneous,bonferroni=bonferroni,
       bonferroni.raw=bonferroni.raw,maximum=maximum,refinement=max(draws[,m+1L]))
}

.np_udens_bb_public <- function(state,grid,alpha,B,band) {
  if(!.np_udens_bb_active(state) || !is.data.frame(grid) || ncol(grid)!=1L ||
     !nrow(grid) || !identical(names(grid),names(state$x)) ||
     !is.numeric(grid[[1L]]) || any(!is.finite(grid[[1L]])) ||
     !is.numeric(alpha) || length(alpha)!=1L || !is.finite(alpha) || alpha<=0 || alpha>=.5 ||
     !isTRUE(band %in% c("pointwise","simultaneous","bonferroni","all")))
    stop("invalid beta-binomial band request",call.=FALSE)
  target <- (grid[[1L]]-state$lower)/state$width
  if(any(!is.finite(target)) || any(target<=0|target>=1))
    stop("beta-binomial targets must be strictly interior; support endpoints are excluded",call.=FALSE)
  B <- .np_plot_gaussian_integer(B,"beta-binomial B",.Machine$integer.max)
  workspace <- .np_udens_bb_budget(state$n,length(unique(target)),B)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste("Beta-binomial density:",stage),done=done,total=total,unit="work")
  inference <- .np_plot_corrected_activity_run("Constructing beta-binomial density intervals",
    .np_with_seed(state$plot.random.seed,{
      reference <- .np_udens_bb_reference(state$normalized,state$h,target,progress)
      tasks <- .np_udens_bb_plan(get(".Random.seed",.GlobalEnv,inherits=FALSE),reference,B)
      draws <- .np_udens_bb_dispatch(tasks,reference)
      z <- .np_udens_bb_reduce(reference,draws,alpha,band)
      ix <- reference$index;center <- reference$center[ix]/state$width;se <- reference$se[ix]/state$width
      make <- function(critical) {
        width <- if(length(critical)==1L) critical*se else critical[ix]*se
        array(c(center-width,center+width),c(nrow(grid),2L,1L),
              dimnames=list(NULL,c("lower","upper"),as.character(alpha)))
      }
      requested <- if(band=="all")c("pointwise","simultaneous","bonferroni") else band
      bands <- setNames(lapply(requested,function(k)make(z[[k]])),requested)
      resolution <- .np_band_resolution(B,rep("density",length(reference$grid)),alpha,
                                        rep(TRUE,length(reference$grid)),requested)
      resolution$calibration <- "smoothed-beta-binomial symmetric-sd";resolution$quantile.type <- 7L
      list(center=center,se=se,bands=bands,intervals=bands,absolute.intervals=bands,
        critical=list(pointwise=z$pointwise[ix],simultaneous=z$simultaneous,bonferroni=z$bonferroni[ix]),
        raw.bonferroni=make(z$bonferroni.raw),maximum=z$maximum,resolution=resolution,covariance=NULL,
        checks=list(original=reference$checks,guide=reference$original.refinement,
          conditional=z$refinement,workspace=workspace),
        reference=list(original.guide=reference$original.guide,variance.guide=reference$reference.guide,
          degree=reference$degree,weights=reference$weights),
        settings=list(construction=state$construction,law="smoothed-beta-binomial",calibration="symmetric-sd",
          kernel="epanechnikov",kernel.order=4L,n=state$n,q=1L,B=B,alpha=alpha,
          bandwidth=state$h*state$width,bounds=c(state$lower,state$lower+state$width),
          grid.size=length(reference$grid),family.rows=ix,family="unique supplied interior targets",
          available.bands=requested,primary="absolute",quantile.type=7L,
          root.center="complete corrected functional of the positive beta-binomial law",
          scale="fixed original complete influence variance, including guide and cross terms",
          bonferroni.construction="type-7 marginal absolute quantile at alpha/m with simultaneous hull; raw endpoints retained",
          correction="degree-five support-reproducing operator with refitted beta guide"))
    }))
  .np_band_resolution_warn(inference$resolution)
  widths <- lapply(inference$bands,function(b)cbind(inference$center-b[,1L,1L],b[,2L,1L]-inference$center))
  list(boot.err=cbind(widths[[if(band=="all")"pointwise" else band]],inference$center),
       boot.all.err=if(band=="all")widths else NULL,inference=inference,bxp=list())
}

.np_udens_bb_options <- function(dots, allow = FALSE) {
  value <- function(public, internal) {
    z <- dots[[public, exact = TRUE]]
    if (is.null(z)) dots[[internal, exact = TRUE]] else z
  }
  law <- value("bootstrap", "plot.errors.boot.method")
  if (!identical(law, "smoothed-beta-binomial")) return(invisible(FALSE))
  if (!isTRUE(allow)) stop("smoothed-beta-binomial is available only for supported unconditional density plots", call. = FALSE)
  ctrl <- dots[["boot.control", exact = TRUE]]
  calibration <- if (!is.null(ctrl)) .np_plot_boot_calibration(ctrl) else
    dots[["plot.errors.boot.calibration", exact = TRUE]]
  if (!identical(value("errors", "plot.errors.method"), "bootstrap") ||
      !identical(value("center", "plot.errors.center"), "bias-corrected") ||
      !isTRUE(value("band", "plot.errors.type") %in% c("pointwise", "simultaneous", "bonferroni", "all")) ||
      !identical(calibration, "symmetric-sd"))
    stop("smoothed-beta-binomial requires explicit bootstrap errors, bias-corrected center, pointwise, simultaneous, bonferroni or all band and symmetric-sd calibration", call. = FALSE)
  if (!is.null(dots[["bias.control", exact = TRUE]]) ||
      isTRUE(dots[["gradients", exact = TRUE]]) || isTRUE(dots[["gradient", exact = TRUE]]) ||
      "gaussian.draws" %in% names(dots) ||
      !is.null(dots[["plot.errors.boot.blocklen", exact = TRUE]]) ||
      (!is.null(ctrl) && (!is.null(ctrl$blocklen) ||
        !identical(ctrl$nonfixed, "exact") || !identical(ctrl$wild, "mammen"))) ||
      (!is.null(dots[["plot.errors.boot.nonfixed", exact = TRUE]]) &&
        !identical(dots[["plot.errors.boot.nonfixed", exact = TRUE]], "exact")) ||
      (!is.null(dots[["plot.errors.boot.wild", exact = TRUE]]) &&
        !identical(dots[["plot.errors.boot.wild", exact = TRUE]], "mammen")))
    stop("smoothed-beta-binomial does not accept custom bias, gradient, Gaussian-draw or dependent/nonfixed/wild controls", call. = FALSE)
  B <- value("B", "plot.errors.boot.num")
  if (!is.null(B) && (!is.numeric(B) || length(B) != 1L || !is.finite(B) ||
      B < 1 || B != floor(B) || B > .Machine$integer.max))
    stop("beta-binomial B must be a positive representable integer", call. = FALSE)
  invisible(TRUE)
}
