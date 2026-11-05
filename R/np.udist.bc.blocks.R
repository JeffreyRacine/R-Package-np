# Literal fixed-length and stationary blocks; first-series selection is an unchanged contract.
.np_udist_bc_block_request <- function(x, law, blocklen = NULL) {
  law <- match.arg(law, c("fixed", "geom"))
  automatic <- is.null(blocklen)
  selected <- if (automatic) b.star(x, round = TRUE) else NULL
  if (automatic) blocklen <- unname(selected[1L, 1L])
  if (!is.numeric(blocklen) || length(blocklen) != 1L || !is.finite(blocklen) ||
      blocklen < 1 || blocklen > nrow(x) || blocklen != floor(blocklen))
    stop("invalid CDF block length", call. = FALSE)
  list(length = as.integer(blocklen), law = law, automatic = automatic,
    source = if (automatic) "b.star(xdat, round=TRUE)[1,1]" else "explicit override",
    series = names(x)[1L], column = "BstarSB", selected = selected,
    caveat = "The first input series may not represent dependence elsewhere in multivariate data.")
}

.np_udist_bc_block_chunk <- function(F, raw, first, last, B, law, blocklen) {
  n <- ncol(F)
  if (blocklen == 1L) {
    N <- .np_bias_process_innovations(n, first, last, B, law, blocklen) + 1
    return(.np_udens_bc_inid_pivot(F, raw, N, accurate = TRUE))
  }
  draw <- .np_block_replicate_drawer(n, n, blocklen, law, TRUE)
  out <- matrix(0, nrow(F), last-first+1L); degenerate <- integer(nrow(F))
  for (b in seq_len(ncol(out))) {
    index <- draw(); values <- F[, index, drop = FALSE]
    numerator <- rowSums(values)
    variance <- .np_bias_score_variance(values, blocklen, law)
    if (any(!is.finite(variance)) || any(variance < 0))
      stop("invalid CDF dependent replicate variance", call. = FALSE)
    zero <- which(variance == 0)
    if (length(zero)) for (j in zero) {
      same <- all(raw[j, index] == raw[j, index[1L]])
      if (!same || numerator[j] == 0)
        stop("undefined or uncertified CDF block pivot", call. = FALSE)
      degenerate[j] <- degenerate[j] + 1L
    }
    out[, b] <- abs(numerator/sqrt(variance))
    if (anyNA(out[, b])) stop("undefined CDF block pivot", call. = FALSE)
  }
  list(absolute = out, degeneracy = degenerate, cancellation = integer(nrow(F)))
}

.np_udist_bc_block_dispatch <- function(F, raw, B, law, blocklen) {
  n <- ncol(F)
  .npRmpi_bootstrap_fanout_enabled(n = n, B = B, chunk.size = 32L, what = "CDF block process")
  tasks <- .npRmpi_bootstrap_chunk_tasks(B, 32L, with.seeds = FALSE)
  tasks <- .npRmpi_bootstrap_attach_rng_stream(tasks, function(task) {
    if (blocklen == 1L) invisible(.np_bias_process_innovations(n, task$start,
      task$start+task$bsz-1L, B, law, blocklen)) else {
      draw <- .np_block_replicate_drawer(n, n, blocklen, law, TRUE)
      for (b in seq_len(task$bsz)) draw()
    }
  }, normal.draws = identical(law, "geom") && blocklen > 1L)
  worker <- function(task, F, raw, B, law, blocklen) {
    .npRmpi_bootstrap_use_task_rng(task)
    z <- .np_udist_bc_block_chunk(F, raw, task$start, task$start+task$bsz-1L, B, law, blocklen)
    out <- t(z$absolute); attr(out, "cdf.pivot") <- z[c("degeneracy", "cancellation")]
    out
  }
  reducer <- function(out, parts, tasks) {
    d <- lapply(parts, attr, which = "cdf.pivot", exact = TRUE)
    if (any(vapply(d, is.null, logical(1L)))) stop("CDF worker omitted diagnostics", call. = FALSE)
    attr(out, "cdf.pivot") <- list(degeneracy = Reduce(`+`, lapply(d, `[[`, "degeneracy")),
      cancellation = Reduce(`+`, lapply(d, `[[`, "cancellation")))
    out
  }
  out <- .npRmpi_bootstrap_run_fanout(tasks, worker, ncol.out = nrow(F),
    what = "CDF block process", profile.where = "mpi.applyLB:cdf-block",
    progress.label = "Calibrating CDF block process", progress.unit = "draws",
    required.bindings = list(.np_udist_bc_block_chunk = .np_udist_bc_block_chunk),
    metadata.reducer = reducer, F = F, raw = raw, B = B, law = law, blocklen = blocklen)
  d <- attr(out, "cdf.pivot", exact = TRUE)
  list(absolute = t(out), degeneracy = d$degeneracy, cancellation = d$cancellation)
}
