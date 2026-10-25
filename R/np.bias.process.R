# Calibration is estimator-neutral. A family adapter supplies complete rows,
# centers and scales. One observation draw is shared across the whole family.
.np_bias_draws <- function(n, size, law) {
  switch(law,
    gaussian=matrix(rnorm(n*size),n),
    mammen={p<-(sqrt(5)+1)/(2*sqrt(5));
      matrix(ifelse(runif(n*size)<p,(1-sqrt(5))/2,(1+sqrt(5))/2),n)},
    rademacher=matrix(ifelse(runif(n*size)<.5,-1,1),n),
    stop("unknown bias multiplier law",call.=FALSE))
}

# The same bounded innovation/chunk contract is used by the MPI transport.
.np_bias_process_innovations <- function(n, first, last, B, law, blocklen) {
  if (law == "inid") return(rmultinom(last-first+1L,n,rep(1/n,n))-1)
  if (law %in% c("fixed","geom"))
    return(.np_block_counts_drawer(n,B,blocklen,law)(first,last)-1)
  .np_bias_draws(n,last-first+1L,law)
}

.np_bias_process_chunk <- function(F, first, last, B, law, blocklen) {
  F %*% .np_bias_process_innovations(ncol(F),first,last,B,law,blocklen)
}

.np_bias_process_draws <- function(F, B, law, blocklen) {
  result <- matrix(0,nrow(F),B)
  for (first in seq.int(1L,B,by=32L)) {
    last <- min(B,first+31L)
    result[,first:last] <- .np_bias_process_chunk(F,first,last,B,law,blocklen)
  }
  result
}

.np_bias_block_selection <- function(residual, method, blocklen=NULL) {
  column <- if(method=="geom") "BstarSB" else "BstarCB"
  if(!is.null(blocklen)) {
    if(!is.numeric(blocklen)||length(blocklen)!=1L||!is.finite(blocklen)||blocklen<1||
       blocklen>length(residual)||blocklen!=floor(blocklen))stop("invalid explicit block length",call.=FALSE)
    return(list(length=as.integer(blocklen),source="override",column=column,method=method))
  }
  if(any(!is.finite(residual))||var(residual)<=0)stop("block selection requires nonconstant finite residuals",call.=FALSE)
  selected <- b.star(matrix(residual,ncol=1L),round=TRUE)
  value <- unname(selected[1L,column])
  if(!is.finite(value)||value<1||value>length(residual)||value!=floor(value))
    stop("b.star returned an invalid block length",call.=FALSE)
  list(length=as.integer(value),source="raw P0 residuals",column=column,
    method=method,selector="b.star",both=selected)
}

# Retain exact marginal order statistics and one draw maximum per family.
# The full grid-by-draw matrix is never retained by inference.
.np_bias_limits_prepare <- function(groups, alpha, gaussian, B) {
  index <- split(seq_along(groups),groups)
  list(index=index,alpha=alpha,gaussian=gaussian,
    maximum=matrix(0,length(index),B),
    point=matrix(0,length(groups),length(alpha)),
    bonf=matrix(0,length(groups),length(alpha)))
}

.np_bias_limits_update <- function(state, absolute, rows) {
  if(nrow(absolute)!=length(rows) || ncol(absolute)!=ncol(state$maximum) ||
     any(!is.finite(absolute)) || any(absolute<0))
    stop("invalid corrected process draws",call.=FALSE)
  for(g in seq_along(state$index)) {
    local <- which(rows %in% state$index[[g]])
    if(!length(local))next
    maximum <- state$maximum[g,]
    for(j in local)maximum <- pmax(maximum,absolute[j,])
    state$maximum[g,] <- maximum
    if(!state$gaussian)for(j in local) {
      q <- unname(quantile(absolute[j,],c(1-state$alpha,1-state$alpha/length(state$index[[g]])),type=1L))
      state$point[rows[j],] <- q[seq_along(state$alpha)]
      state$bonf[rows[j],] <- q[length(state$alpha)+seq_along(state$alpha)]
    }
  }
  state
}

.np_bias_limits_finish <- function(state) {
  alpha <- state$alpha;index <- state$index
  limits <- setNames(lapply(c("pointwise","simultaneous","bonferroni"),function(b)
    array(0,c(nrow(state$point),2L,length(alpha)),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))),
    c("pointwise","simultaneous","bonferroni"))
  raw.max <- matrix(0,length(alpha),length(index),dimnames=list(as.character(alpha),names(index)))
  for(g in seq_along(index)) {
    ii <- index[[g]];m <- length(ii)
    raw.max[,g] <- unname(quantile(state$maximum[g,],1-alpha,type=1L))
    if(state$gaussian) {
      point <- matrix(rep(qnorm(alpha/2,lower.tail=FALSE),each=m),m,length(alpha))
      bonf <- matrix(rep(qnorm(alpha/(2*m),lower.tail=FALSE),each=m),m,length(alpha))
      simultaneous <- pmin(bonf,pmax(point,matrix(rep(raw.max[,g],each=m),m,length(alpha))))
    } else {
      point <- state$point[ii,,drop=FALSE];bonf <- state$bonf[ii,,drop=FALSE]
      simultaneous <- matrix(rep(raw.max[,g],each=m),m,length(alpha))
      if(any(point>simultaneous))stop("maximum marginal-order invariant failed",call.=FALSE)
      bonf <- pmax(bonf,simultaneous)
    }
    for(band in names(limits)) {
      critical <- switch(band,pointwise=point,simultaneous=simultaneous,bonferroni=bonf)
      limits[[band]][ii,1L,] <- -critical
      limits[[band]][ii,2L,] <- critical
    }
  }
  list(limits=limits,raw_max_critical=raw.max)
}

.np_bias_limits <- function(absolute,groups,alpha,gaussian=FALSE) {
  state <- .np_bias_limits_prepare(groups,alpha,gaussian,ncol(absolute))
  .np_bias_limits_finish(.np_bias_limits_update(state,absolute,seq_len(nrow(absolute))))
}

# Replaying a draw stream for bounded target batches does not consume a new
# independent stream. Advance the caller once, by one complete draw sequence.
.np_bias_stream <- function(state, code) {
  previous <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
  on.exit(assign(".Random.seed",previous,.GlobalEnv))
  assign(".Random.seed",state,.GlobalEnv)
  value <- force(code)
  list(value=value,end=get(".Random.seed",.GlobalEnv,inherits=FALSE))
}

# Partition both axes without changing the target-major family order or draws.
.np_bias_target_blocks <- function(n.grid, n.targets, anchors = 25L, rows = 500L) {
  if (any(!is.finite(c(n.grid,n.targets,anchors,rows))) ||
      any(c(n.grid,n.targets,anchors,rows) < 1) ||
      any(c(n.grid,n.targets,anchors,rows) != floor(c(n.grid,n.targets,anchors,rows))) ||
      n.grid > .Machine$integer.max / n.targets)
    stop("invalid or unrepresentable bias target grid", call.=FALSE)
  count <- min(n.targets, rows)
  grid.block <- min(anchors, rows %/% count)
  out <- list()
  for (first in seq.int(1L,n.grid,by=grid.block)) {
    jj <- seq.int(first,min(n.grid,first+grid.block-1L))
    for (target in seq.int(1L,n.targets,by=count)) {
      tt <- seq.int(target,min(n.targets,target+count-1L))
      out[[length(out)+1L]] <- list(grid=jj,
        rows=as.vector(outer(jj,(as.double(tt)-1)*n.grid,"+")))
    }
  }
  out
}

.np_bias_inference <- function(state, grid, targets="mean", alpha=.05,
                               B=9999L, law="gaussian", blocklen=NULL, fit=NULL,
                               studentize=TRUE) {
  if (isTRUE(state$categorical.only) && !law %in% c("gaussian","mammen","inid","rademacher","fixed","geom"))
    stop("categorical-only bias correction does not support this calibration law",call.=FALSE)
  higher <- !is.null(state$gradient.order) && any(vapply(targets,function(t)
    startsWith(t,"d") && .np_bias_derivative_order(state,as.integer(substring(t,2L)))>1L,logical(1)))
  if (higher && !law %in% c("gaussian","mammen","inid","rademacher","fixed","geom"))
    stop("higher-order bias correction does not support this calibration law",call.=FALSE)
  if (isTRUE(state$pooled) &&
      !law %in% c("gaussian","mammen","inid"))
    stop("pooled mixed-data bias correction currently supports Gaussian, Mammen or inid calibration only",call.=FALSE)
  if(!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<2||B!=floor(B)||B>.Machine$integer.max)
    stop("bias calibration draws must be an integer >= 2",call.=FALSE)
  if(!length(alpha)||any(!is.finite(alpha))||any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))
    stop("alpha must contain distinct values in (0, 0.5)",call.=FALSE)
  if(!law%in%c("gaussian","mammen","rademacher","inid","fixed","geom"))stop("unsupported corrected inference law",call.=FALSE)
  if(!identical(names(grid),names(state$x))||!nrow(grid)||anyNA(grid))stop("invalid bias evaluation grid",call.=FALSE)
  if(!is.character(targets)||!length(targets)||anyNA(targets)||anyDuplicated(targets))
    stop("invalid or duplicate bias target family",call.=FALSE)
  blocks <- .np_bias_target_blocks(nrow(grid),length(targets))
  map <- expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
  M <- nrow(map);B <- as.integer(B)
  center <- original <- se <- score.se <- numeric(M)
  student <- isTRUE(studentize) && law %in% c("gaussian","mammen")
  calibration.state <- .np_bias_limits_prepare(map$target,alpha,law=="gaussian" && !student,B)
  prepared <- NULL
  regularization <- adaptation <- list()
  adapted.rows <- if (!is.null(state$adaptation))
    .np_bias_adaptation_capture(state,fit,grid,targets) else NULL
  score <- law%in%c("inid","fixed","geom")
  selection <- if(law%in%c("fixed","geom")) .np_bias_block_selection(state$residual_raw,law,blocklen) else NULL
  if(!exists(".Random.seed",.GlobalEnv,inherits=FALSE))runif(1L)
  rng <- get(".Random.seed",.GlobalEnv,inherits=FALSE);rng.kind <- RNGkind();end <- rng
  for(batch in blocks) {
    jj <- batch$grid;rr <- batch$rows
    local.map <- map[rr,,drop=FALSE];local.map$row <- match(local.map$row,jj)
    value <- if (is.null(adapted.rows)) .np_bias_regression_rows(state,grid[jj,,drop=FALSE],local.map) else
      adapted.rows(jj,local.map)
    if (!is.null(value$adaptation)) adaptation[[length(adaptation)+1L]] <-
      list(rows=rr,grid=jj,details=value$adaptation)
    diagnostic <- value$main.regularization
    if(law!="gaussian" && any(diagnostic$ridge>0))
      stop("bootstrap calibration for regularized generalized-NN main rows remains unqualified; Gaussian inference is available with the qualifications documented in np_bias_control",call.=FALSE)
    diagnostic$row <- jj[diagnostic$row]
    regularization[[length(regularization)+1L]] <- diagnostic
    center[rr] <- value$center;original[rr] <- value$original_center;se[rr] <- value$se
    divisor <- ifelse(value$se==0,1,value$se)
    if (student) {
      if (is.null(prepared)) prepared <- .np_bias_student_prepare(state)
      part <- .np_bias_stream(rng,
        .np_bias_student_draws(prepared,value$W,value$se,B,law))
    } else {
      F <- if(score)sweep(value$W,2L,state$residual,"*") else value$factor
      if(score)F <- F-rowMeans(F)
      part <- .np_bias_stream(rng,
        .np_bias_process_draws(F,B,law,if(is.null(selection))NULL else selection$length))
    }
    .np_plot_progress_notify(stage="Constructing bias-corrected intervals",
      done=max(jj),total=nrow(grid),unit="grid")
    end <- part$end
    if(score) {
      s <- apply(part$value,1L,sd)
      if(any(!is.finite(s))||any(s<=0 & value$se>0))stop("nonpositive score-bootstrap scale",call.=FALSE)
      score.se[rr] <- s
      divisor <- ifelse(s==0,1,s)
      if(law!="inid")se[rr] <- s
    }
    calibration.state <- .np_bias_limits_update(calibration.state,
      if(student)part$value else abs(part$value/divisor),rr)
  }
  assign(".Random.seed",end,.GlobalEnv)
  calibrated <- .np_bias_limits_finish(calibration.state)
  bands <- lapply(calibrated$limits,function(x) {
    for(k in seq_along(alpha))for(side in 1:2)x[,side,k] <- center+se*x[,side,k]
    x
  })
  list(construction="iterated",version=if(student)2L else 1L,estimate=original,center=center,
    correction=original-center,se=se,bands=bands,map=map,grid=grid,
    settings=state$settings,law=law,block.selection=selection,draws=B,alpha=alpha,
    main.regularization=do.call(rbind,regularization),
    adaptation=if(length(adaptation))adaptation else NULL,
    covariance=list(operator="complete corrected observation rows",
      variance=if(law%in%c("fixed","geom"))"empirical residual-score process variance" else
        if(isTRUE(state$categorical.only))"within-profile residual variance" else "smoothed P0 residual variance",
      studentization=if(isTRUE(state$categorical.only) && student)"replicate-specific within-profile residual SD" else
        if(score)"empirical residual-score process SD" else
        if(isTRUE(state$categorical.only))"within-profile residual SD" else
        if(student)"replicate-specific smoothed P0 residual SD" else "smoothed complete-row SD",
      calibration=if(student)"residual-pilot studentized multiplier process" else
        if(score)"residual score process" else "fixed-factor multiplier process"),
    studentization=if(student)list(construction=if(isTRUE(state$categorical.only))
      "finite-profile-studentized-v1" else "residual-pilot-v1",
      draw.block=256L,reuse.block=1024L,training.block=16L) else NULL,
    calibration=calibrated,score.se=if(score)score.se else NULL,
    rng=list(kind=rng.kind,start=rng,end=end,convention="common observation draws across all target blocks"),
    scope="finite grid within each scalar-target family",owner="bounded corrected observation factors")
}

.np_bias_plot_errors <- function(state, grid, target, method, bootstrap, wild,
                                 alpha, B, gaussian.draws, blocklen, band, fit=NULL) {
  law <- if(method=="asymptotic")"gaussian" else if(bootstrap=="wild")wild else bootstrap
  inference <- .np_plot_activity_run("Constructing bias-corrected intervals",
    .np_bias_inference(state,grid,target,alpha,
      if(method=="asymptotic")gaussian.draws else B,law,blocklen,fit,
      studentize=band!="pmzsd"))
  all <- lapply(inference$bands,function(x)cbind(inference$center-x[,1,1],x[,2,1]-inference$center))
  width <- if(band=="pmzsd")cbind(inference$se,inference$se)*qnorm(alpha/2,lower.tail=FALSE) else
    all[[if(band=="all")"pointwise" else band]]
  if(is.null(width))stop("unsupported bias interval type",call.=FALSE)
  list(boot.err=cbind(width,inference$center),boot.all.err=if(band=="all")all else NULL,
    inference=inference,bxp=list())
}
