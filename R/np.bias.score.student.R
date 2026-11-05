# Exact conditional score variance for the literal circular block law.
# Target rows are bounded by inference; storage is rows by n, never n by n.
.np_bias_score_variance <- function(F, blocklen, law) {
  F <- F-rowMeans(F)
  n <- ncol(F)
  if (!law %in% c("fixed","geom")) stop("unknown score-studentized block law",call.=FALSE)
  if (blocklen==1L) return(rowSums(F^2))
  if (law=="geom") {
    # Finite-n circular autocovariances; the end-of-sequence weights matter.
    z <- stats::mvfft(t(F))
    ac <- Re(stats::mvfft(Conj(z)*z,inverse=TRUE))/n
    k <- seq_len(n-1L)
    w <- c(n,2*(n-k)*(1-1/blocklen)^k)/n
    return(drop(crossprod(w,ac)))
  }
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
                                          B, law, blocklen, studentize, signed=FALSE) {
  n <- ncol(F)
  out <- matrix(0,nrow(F),last-first+1L)
  if (blocklen==1L) {
    counts <- .np_bias_process_innovations(n,first,last,B,law,blocklen)+1
    numerator <- F%*%counts
    variance <- if(studentize) (F^2)%*%counts-numerator^2/n else
      matrix(observed,nrow(F),ncol(counts))
    if (any(!is.finite(variance)) || any(variance[active,,drop=FALSE]<=0))
      stop("nonpositive replicate block-one score variance",call.=FALSE)
    out[active,] <- if(signed) numerator[active,,drop=FALSE]/sqrt(variance[active,,drop=FALSE]) else
      abs(numerator[active,,drop=FALSE]/sqrt(variance[active,,drop=FALSE]))
  } else {
    draw <- .np_block_replicate_drawer(n,n,blocklen,law,TRUE)
    for (b in seq_len(ncol(out))) {
      Z <- F[,draw(),drop=FALSE]
      numerator <- rowSums(Z)
      variance <- if(studentize) .np_bias_score_variance(Z,blocklen,law) else observed
      if (any(!is.finite(variance)) || any(variance[active]<=0))
        stop("nonpositive replicate dependent score variance",call.=FALSE)
      out[active,b] <- if(signed) numerator[active]/sqrt(variance[active]) else
        abs(numerator[active]/sqrt(variance[active]))
    }
  }
  out
}

.np_bias_score_student_draws <- function(F, B, law, blocklen, active, studentize=TRUE, signed=FALSE) {
  F <- F-rowMeans(F)
  observed <- .np_bias_score_variance(F,blocklen,law)
  if (length(active)!=nrow(F) || anyNA(active) || any(!is.finite(observed)) ||
      any(observed[active]<=0) || any(observed[!active]!=0))
    stop("nonpositive or inconsistent selected-block score variance",call.=FALSE)
  value <- .np_bias_score_student_dispatch(F,active,observed,B,law,blocklen,studentize,signed=signed)
  if(signed) list(signed=value,se=sqrt(observed)) else list(absolute=value,se=sqrt(observed))
}

.np_bias_score_student_dispatch <- function(F, active, observed, B, law, blocklen, studentize, signed=FALSE) {
  n <- ncol(F)
  .npRmpi_bootstrap_fanout_enabled(n=n,B=B,chunk.size=32L,what="studentized score process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B,32L,with.seeds=FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks,function(task)
    invisible(.np_bias_process_innovations(n,task$start,task$start+task$bsz-1L,B,law,blocklen)),
    normal.draws=identical(law,"geom") && blocklen>1L)
  worker <- function(task,factors,active,observed,B,law,blocklen,studentize,signed) {
    .npRmpi_bootstrap_use_task_rng(task)
    t(.np_bias_score_student_chunk(factors,active,observed,task$start,
      task$start+task$bsz-1L,B,law,blocklen,studentize,signed=signed))
  }
  t(.npRmpi_bootstrap_run_fanout(tasks,worker,ncol.out=nrow(F),
    what="studentized score process",profile.where="mpi.applyLB:bias-score-student",
    progress.label="Calibrating studentized score process",progress.unit="draws",
    factors=F,active=active,observed=observed,B=B,law=law,blocklen=blocklen,
    studentize=studentize,signed=signed))
}
