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
  out <- matrix(0,nrow(F),B)
  for (first in seq.int(1L,B,by=32L)) {
    last <- min(B,first+31L)
    out[,first:last] <- .np_bias_score_student_chunk(F,active,observed,first,last,
      B,law,blocklen,studentize)
  }
  out
}
