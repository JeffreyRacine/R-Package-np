# Rank-independent paired donor and coordinate-noise streams.
.np_udens_q6_q2_plan <- function(start,n,B,chunk=32L) {
  .np_udens_q6_q2_budget(n,1L,c(1L,1L),B,chunk)
  .np_bias_replay_check("inid",NULL,TRUE,seed=start)
  .np_bias_replay_check("gaussian",NULL,TRUE,seed=start)
  n <- as.integer(n); B <- as.integer(B); chunk <- as.integer(chunk)
  starts <- seq.int(1L,B,by=chunk)
  stream <- .np_bias_stream(start, {
    tasks <- vector("list",length(starts))
    for (i in seq_along(starts)) {
      size <- min(chunk,B-starts[i]+1L)
      state <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      invisible(sample.int(n,n*size,replace=TRUE))
      tasks[[i]] <- list(start=starts[i],bsz=size,seed=NA_integer_,rng_state=state,
        donor=state,donor.end=get(".Random.seed",.GlobalEnv,inherits=FALSE),
        noise=vector("list",2L),noise.end=vector("list",2L))
    }
    for (j in 1:2) for (i in seq_along(tasks)) {
      tasks[[i]]$noise[[j]] <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      invisible(rnorm(n*tasks[[i]]$bsz))
      tasks[[i]]$noise.end[[j]] <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      .np_plot_progress_notify(stage="Preparing paired q=2 Q6 draws",
        done=(j-1L)*length(tasks)+i,total=2L*length(tasks),unit="blocks")
    }
    tasks
  })
  attr(stream$value,"rng_final_state") <- stream$end
  stream$value
}

.np_udens_q6_q2_samples <- function(task,x,h) {
  n <- nrow(x)
  donor <- .np_bias_stream(task$donor,sample.int(n,n*task$bsz,replace=TRUE))
  if (!identical(donor$end,task$donor.end))
    stop("q=2 Q6 donor replay consumed an inconsistent stream",call.=FALSE)
  lapply(1:2,function(j) {
    noise <- .np_bias_stream(task$noise[[j]],rnorm(n*task$bsz))
    if (!identical(noise$end,task$noise.end[[j]]))
      stop("q=2 Q6 noise replay consumed an inconsistent stream",call.=FALSE)
    probability <- pnorm(noise$value)
    if (any(!is.finite(probability)) || any(probability<=0 | probability>=1))
      stop("Q6 normal transform reached an unrepresentable probability",call.=FALSE)
    value <- x[donor$value,j]+h[j]*sqrt(2/3)*qt(probability,df=6)
    if (any(!is.finite(value))) stop("nonfinite q=2 Q6 conditional sample",call.=FALSE)
    matrix(value,n,task$bsz)
  })
}

.np_udens_q6_q2_chunk <- function(task,reference) {
  Y <- .np_udens_q6_q2_samples(task,reference$x,reference$h)
  low <- .np_udens_q6_q2_curve(Y,reference$grid,reference$h,reference$kernel,24L)
  high <- .np_udens_q6_q2_curve(Y,reference$grid,reference$h,reference$kernel,32L)
  error <- .np_udens_q6_check(low,high,"q=2 conditional curve")
  if (isTRUE(reference[["directional", exact=TRUE]]))
    return(cbind(t(sweep(high,1L,reference$reference.center)),
      refinement=rep(error,task$bsz)))
  roots <- abs(sweep(sweep(high,1L,reference$reference.center),1L,reference$se,"/"))
  cbind(t(roots),refinement=rep(error,task$bsz))
}

.np_udens_q6_q2_dispatch <- function(tasks,reference) {
  .npRmpi_bootstrap_fanout_enabled(n=nrow(reference$x),
    B=sum(vapply(tasks,function(z) z$bsz,0L)),chunk.size=32L,
    what="q=2 Q6 density process")
  worker <- function(task,reference) .np_udens_q6_q2_chunk(task,reference)
  .npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(reference$grid)+1L,
    what="q=2 Q6 density process",profile.where="mpi.applyLB:q6-density-q2",
    progress.label="Constructing paired q=2 Q6 bands",
    required.bindings=list(.np_udens_q6_q2_chunk=.np_udens_q6_q2_chunk,
      .np_udens_q6_q2_samples=.np_udens_q6_q2_samples,
      .np_udens_q6_q2_curve=.np_udens_q6_q2_curve,
      .np_udens_q6_q2_rule=.np_udens_q6_q2_rule,
      .np_udens_q6_rule=.np_udens_q6_rule,
      .np_udens_q6_kernel=.np_udens_q6_kernel,
      .np_udens_q6_check=.np_udens_q6_check),reference=reference)
}
