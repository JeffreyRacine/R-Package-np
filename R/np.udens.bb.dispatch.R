# All mixture labels precede all beta variates, as in the frozen research law.
# Plans retain primitive stream states, never a B-by-n bank of observations.
.np_udens_bb_plan <- function(start,reference,B,chunk=16L) {
  n <- reference$n;.np_udens_bb_budget(n,length(reference$grid),B,chunk)
  .np_bias_replay_check("inid",NULL,TRUE,seed=start)
  starts <- seq.int(1L,B,by=chunk)
  stream <- .np_bias_stream(start,{
    tasks <- vector("list",length(starts))
    for(i in seq_along(starts)) {
      size <- min(chunk,B-starts[i]+1L)
      before <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      invisible(sample.int(length(reference$weights),n*size,replace=TRUE,prob=reference$weights))
      tasks[[i]] <- list(start=starts[i],bsz=size,seed=NA_integer_,rng_state=before,
        label=before,label.end=get(".Random.seed",.GlobalEnv,inherits=FALSE))
    }
    for(i in seq_along(tasks)) {
      label <- .np_bias_stream(tasks[[i]]$label,
        sample.int(length(reference$weights),n*tasks[[i]]$bsz,replace=TRUE,prob=reference$weights))
      if(!identical(label$end,tasks[[i]]$label.end)) stop("inconsistent beta-binomial label replay",call.=FALSE)
      tasks[[i]]$beta <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      invisible(.np_udens_bb_draw(reference$A[label$value],reference$BB[label$value]))
      tasks[[i]]$beta.end <- get(".Random.seed",.GlobalEnv,inherits=FALSE)
      .np_plot_progress_notify(stage="Preparing beta-binomial draws",done=i,total=length(tasks),unit="blocks")
    }
    tasks
  })
  attr(stream$value,"rng_final_state") <- stream$end
  stream$value
}

.np_udens_bb_samples <- function(task,reference) {
  label <- .np_bias_stream(task$label,
    sample.int(length(reference$weights),reference$n*task$bsz,replace=TRUE,prob=reference$weights))
  draw <- .np_bias_stream(task$beta,
    .np_udens_bb_draw(reference$A[label$value],reference$BB[label$value]))
  if(!identical(label$end,task$label.end) || !identical(draw$end,task$beta.end))
    stop("inconsistent beta-binomial primitive replay",call.=FALSE)
  draw$value
}

.np_udens_bb_chunk <- function(task,reference) {
  sample <- .np_udens_bb_samples(task,reference);m <- length(reference$grid)
  out <- matrix(0,task$bsz,m+1L)
  for(b in seq_len(task$bsz)) {
    ix <- seq.int(1L+(b-1L)*reference$n,length.out=reference$n)
    curve <- .np_udens_bb_curve(reference$plan,sample[ix,1L],sample[ix,2L],sample[ix,3L])
    out[b,seq_len(m)] <- (curve$center-reference$reference.center)/reference$se
    out[b,m+1L] <- curve$refinement
  }
  out
}

.np_udens_bb_dispatch <- function(tasks,reference) {
  B <- sum(vapply(tasks,function(t)t$bsz,0L));m <- length(reference$grid)
  out <- matrix(0,B,m+1L)
  for(task in tasks) {
    ix <- seq.int(task$start,length.out=task$bsz)
    out[ix,] <- .np_udens_bb_chunk(task,reference)
    .np_plot_progress_notify(stage="Constructing beta-binomial bands",done=max(ix),total=B,unit="draws")
  }
  out
}
