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
  n <- ncol(F)
  .npRmpi_bootstrap_fanout_enabled(n=n,B=B,chunk.size=32L,what="corrected process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,32L,with.seeds=FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
    invisible(.np_bias_process_innovations(n,task$start,task$start+task$bsz-1L,B,law,blocklen)))
  worker <- function(task,F,B,law,blocklen) {
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_bias_process_chunk(F,task$start,task$start+task$bsz-1L,B,law,blocklen))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(F),
    what="corrected process",profile.where="mpi.applyLB:bias-process",
    progress.label="Calibrating corrected process",F=F,B=B,law=law,blocklen=blocklen))
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

.np_bias_limits <- function(absolute,groups,alpha,gaussian=FALSE) {
  if(any(!is.finite(absolute))||any(absolute<0))stop("invalid corrected process draws",call.=FALSE)
  index <- split(seq_len(nrow(absolute)),groups)
  limits <- setNames(lapply(c("pointwise","simultaneous","bonferroni"),function(b)
    array(0,c(nrow(absolute),2L,length(alpha)),dimnames=list(NULL,c("lower","upper"),as.character(alpha)))),
    c("pointwise","simultaneous","bonferroni"))
  raw.max <- matrix(0,length(alpha),length(index),dimnames=list(as.character(alpha),names(index)))
  for(g in seq_along(index)) {
    ii <- index[[g]];m <- length(ii)
    maximum <- apply(absolute[ii,,drop=FALSE],2L,max)
    raw.max[,g] <- unname(quantile(maximum,1-alpha,type=1L))
    for(j in ii) {
      if(gaussian) {
        point <- qnorm(alpha/2,lower.tail=FALSE)
        bonf <- qnorm(alpha/(2*m),lower.tail=FALSE)
        simultaneous <- pmin(bonf,pmax(point,raw.max[,g]))
      } else {
        q <- unname(quantile(absolute[j,],c(1-alpha,1-alpha/m),type=1L))
        point <- q[seq_along(alpha)];bonf <- q[length(alpha)+seq_along(alpha)]
        simultaneous <- raw.max[,g]
        if(any(point>simultaneous))stop("maximum marginal-order invariant failed",call.=FALSE)
        bonf <- pmax(bonf,simultaneous)
      }
      for(band in names(limits)) {
        critical <- switch(band,pointwise=point,simultaneous=simultaneous,bonferroni=bonf)
        limits[[band]][j,1L,] <- -critical
        limits[[band]][j,2L,] <- critical
      }
    }
  }
  list(limits=limits,raw_max_critical=raw.max)
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

.np_bias_inference <- function(state, grid, targets="mean", alpha=.05,
                               B=9999L, law="gaussian", blocklen=NULL) {
  if(!is.numeric(B)||length(B)!=1L||!is.finite(B)||B<2||B!=floor(B)||B>.Machine$integer.max)
    stop("bias calibration draws must be an integer >= 2",call.=FALSE)
  if(!length(alpha)||any(!is.finite(alpha))||any(alpha<=0|alpha>=.5)||anyDuplicated(alpha))
    stop("alpha must contain distinct values in (0, 0.5)",call.=FALSE)
  if(!law%in%c("gaussian","mammen","rademacher","inid","fixed","geom"))stop("unsupported corrected inference law",call.=FALSE)
  if(!identical(names(grid),names(state$x))||!nrow(grid)||anyNA(grid))stop("invalid bias evaluation grid",call.=FALSE)
  map <- expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE)
  if(anyDuplicated(targets))stop("duplicate bias target family",call.=FALSE)
  M <- nrow(map);B <- as.integer(B)
  center <- original <- se <- score.se <- numeric(M)
  draws <- matrix(0,M,B)
  score <- law%in%c("inid","fixed","geom")
  selection <- if(law%in%c("fixed","geom")) .np_bias_block_selection(state$residual_raw,law,blocklen) else NULL
  if(!exists(".Random.seed",.GlobalEnv,inherits=FALSE))runif(1L)
  rng <- get(".Random.seed",.GlobalEnv,inherits=FALSE);rng.kind <- RNGkind();end <- rng
  for(first in seq.int(1L,nrow(grid),by=25L)) {
    jj <- seq.int(first,min(nrow(grid),first+24L));rr <- which(map$row%in%jj)
    local.map <- map[rr,,drop=FALSE];local.map$row <- match(local.map$row,jj)
    value <- .np_bias_regression_rows(state,grid[jj,,drop=FALSE],local.map)
    center[rr] <- value$center;original[rr] <- value$original_center;se[rr] <- value$se
    F <- if(score)sweep(value$W,2L,state$residual,"*") else value$factor
    if(score)F <- F-rowMeans(F)
    divisor <- ifelse(value$se==0,1,value$se)
    part <- .np_bias_stream(rng,
      .np_bias_process_draws(F,B,law,if(is.null(selection))NULL else selection$length))
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
    draws[rr,] <- abs(part$value/divisor)
  }
  assign(".Random.seed",end,.GlobalEnv)
  calibrated <- .np_bias_limits(draws,map$target,alpha,law=="gaussian")
  bands <- lapply(calibrated$limits,function(x) {
    for(k in seq_along(alpha))for(side in 1:2)x[,side,k] <- center+se*x[,side,k]
    x
  })
  list(construction="iterated",version=1L,estimate=original,center=center,
    correction=original-center,se=se,bands=bands,map=map,grid=grid,
    settings=state$settings,law=law,block.selection=selection,draws=B,alpha=alpha,
    calibration=calibrated,score.se=if(score)score.se else NULL,
    rng=list(kind=rng.kind,start=rng,end=end,convention="common observation draws across all target blocks"),
    scope="finite grid within each scalar-target family",owner="bounded corrected observation factors")
}

.np_bias_plot_errors <- function(state, grid, target, method, bootstrap, wild,
                                 alpha, B, gaussian.draws, blocklen, band) {
  law <- if(method=="asymptotic")"gaussian" else if(bootstrap=="wild")wild else bootstrap
  inference <- .np_plot_activity_run("Constructing bias-corrected intervals",
    .np_bias_inference(state,grid,target,alpha,
      if(method=="asymptotic")gaussian.draws else B,law,blocklen))
  all <- lapply(inference$bands,function(x)cbind(inference$center-x[,1,1],x[,2,1]-inference$center))
  width <- if(band=="pmzsd")cbind(inference$se,inference$se)*qnorm(alpha/2,lower.tail=FALSE) else
    all[[if(band=="all")"pointwise" else band]]
  if(is.null(width))stop("unsupported bias interval type",call.=FALSE)
  list(boot.err=cbind(width,inference$center),boot.all.err=if(band=="all")all else NULL,
    inference=inference,bxp=list())
}
