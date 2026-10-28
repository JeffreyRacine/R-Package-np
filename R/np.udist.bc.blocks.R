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
  out <- matrix(0, nrow(F), B); degenerate <- cancellation <- integer(nrow(F))
  for (first in seq.int(1L, B, by = 32L)) {
    last <- min(B, first + 31L)
    z <- .np_udist_bc_block_chunk(F, raw, first, last, B, law, blocklen)
    out[, first:last] <- z$absolute
    degenerate <- degenerate + z$degeneracy; cancellation <- cancellation + z$cancellation
    .np_plot_progress_notify(stage = "Calibrating CDF block process", done = last, total = B, unit = "draws")
  }
  list(absolute = out, degeneracy = degenerate, cancellation = cancellation)
}
