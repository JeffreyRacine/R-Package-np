# Exact conditional score variance for the literal circular block law.
# Target rows are bounded by inference; storage is rows by n, never n by n.
.np_bias_score_variance <- function(F, blocklen, law) {
  F <- F-rowMeans(F)
  n <- ncol(F)
  if (law!="fixed") stop("unknown score-studentized block law",call.=FALSE)
  if (blocklen==1L) return(rowSums(F^2))
  if (blocklen==n) return(numeric(nrow(F)))
  variance <- function(width) {
    if (!width) return(numeric(nrow(F)))
    S <- matrix(0,nrow(F),n)
    for (k in seq_len(width))
      S <- S+F[,((seq_len(n)+k-2L)%%n)+1L,drop=FALSE]
    rowMeans(S^2)
  }
  (n%/%blocklen)*variance(blocklen)+variance(n%%blocklen)
}

.np_bias_score_student_chunk <- function(F, active, observed, first, last,
                                          B, law, blocklen, studentize) {
  n <- ncol(F)
  out <- matrix(0,nrow(F),last-first+1L)
  if (blocklen==1L) {
    counts <- .np_bias_process_innovations(n,first,last,B,law,blocklen)+1
    numerator <- F%*%counts
    variance <- if(studentize) (F^2)%*%counts-numerator^2/n else
      matrix(observed,nrow(F),ncol(counts))
    if (any(!is.finite(variance)) || any(variance[active,,drop=FALSE]<=0))
      stop("nonpositive replicate block-one score variance",call.=FALSE)
    out[active,] <- abs(numerator[active,,drop=FALSE]/sqrt(variance[active,,drop=FALSE]))
  } else {
    draw <- .np_block_replicate_drawer(n,n,blocklen,law,TRUE)
    for (b in seq_len(ncol(out))) {
      Z <- F[,draw(),drop=FALSE]
      numerator <- rowSums(Z)
      variance <- if(studentize) .np_bias_score_variance(Z,blocklen,law) else observed
      if (any(!is.finite(variance)) || any(variance[active]<=0))
        stop("nonpositive replicate dependent score variance",call.=FALSE)
      out[active,b] <- abs(numerator[active]/sqrt(variance[active]))
    }
  }
  out
}

.np_bias_score_student_draws <- function(F, B, law, blocklen, active, studentize=TRUE) {
  F <- F-rowMeans(F)
  observed <- .np_bias_score_variance(F,blocklen,law)
  if (length(active)!=nrow(F) || anyNA(active) || any(!is.finite(observed)) ||
      any(observed[active]<=0) || any(observed[!active]!=0))
    stop("nonpositive or inconsistent selected-block score variance",call.=FALSE)
  absolute <- .np_bias_score_student_dispatch(F,active,observed,B,law,blocklen,studentize)
  list(absolute=absolute,se=sqrt(observed))
}

.np_bias_score_student_dispatch <- function(F, active, observed, B, law, blocklen, studentize) {
  n <- ncol(F)
  .npRmpi_bootstrap_fanout_enabled(n=n,B=B,chunk.size=32L,what="studentized score process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,32L,with.seeds=FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
    invisible(.np_bias_process_innovations(n,task$start,task$start+task$bsz-1L,B,law,blocklen)))
  worker <- function(task,factors,active,observed,B,law,blocklen,studentize) {
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_bias_score_student_chunk(factors,active,observed,task$start,
      task$start+task$bsz-1L,B,law,blocklen,studentize))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(F),
    what="studentized score process",profile.where="mpi.applyLB:bias-score-student",
    progress.label="Calibrating studentized score process",progress.unit="draws",
    factors=F,active=active,observed=observed,B=B,law=law,blocklen=blocklen,
    studentize=studentize))
}
