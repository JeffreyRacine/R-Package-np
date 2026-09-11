# Requested uniform-density uncertainty only: full marked moving faces.
.np_ann_face_window <- function(z, at) {
  n <- length(z)
  rank <- findInterval(at, z)
  m <- min(floor((n - 1)^(2/3)), floor(min(rank - 1L, n - rank)/2))
  if (m < 1L) return(NULL)
  lo <- rank - m
  hi <- rank + m
  width <- z[hi] - z[lo]
  if (!is.finite(width) || width <= 0) return(NULL)
  density <- ((hi - lo)/n)/width
  if (!is.finite(density) || density <= 0) return(NULL)
  list(lo = lo, hi = hi, width = width, density = density)
}

.np_ann_uniform_prepare <- function(kbw, data) {
  con <- which(kbw[["icon", exact = TRUE]])
  other <- lapply(con, function(j) {
    keep <- setdiff(seq_len(ncol(data)), j)
    if (!length(keep)) return(list(keep = keep, bw = NULL))
    # The full bandwidth already uses the canonical density-family mapping.
    # Omit the irrelevant order argument for a uniform kernel.
    bw <- kbandwidth.numeric(bw = kbw[["bw", exact = TRUE]][keep],
      bwtype = "adaptive_nn", ckertype = "uniform", ckerbound = "none",
      ukertype = kbw[["ukertype", exact = TRUE]],
      okertype = kbw[["okertype", exact = TRUE]],
      nobs = nrow(data), xdati = untangle(data[keep]), xnames = names(data)[keep])
    list(keep = keep, bw = bw)
  })
  list(other = other, k = kbw[["bw", exact = TRUE]][con])
}

.np_ann_uniform_faces <- function(state, data, evaluation, e) {
  n <- nrow(data)
  m <- nrow(evaluation)
  p <- ncol(state$x)
  cut <- coef <- matrix(NA_real_, m, 2L * p)
  for (j in seq_len(p)) {
    z <- state$x[state$order[[j]], j]
    other <- state$uniform$other[[j]]
    marks <- if (!length(other$keep)) matrix(1, n, m) else
      npksum(bws = other$bw, txdat = data[other$keep],
        exdat = evaluation[other$keep], operator = "normal",
        bandwidth.divide = TRUE, return.kernel.weights = TRUE,
        .np.internal.bandwidth.divide.weights = TRUE)[["kw", exact = TRUE]]
    xq <- e[, j]
    for (q in seq_len(m)) {
      Fq <- findInterval(xq[q], z)/n
      prob <- Fq + c(-1, 1) * state$uniform$k[j]/(n - 1)
      for (h in 1:2) {
        column <- 2L * (j - 1L) + h
        # Strictly absent support faces remain at infinity locally; they
        # contribute zero. Equality at the probability endpoint is nonregular.
        if (prob[h] < 0 || prob[h] > 1) {
          cut[q, column] <- xq[q]
          coef[q, column] <- 0
          next
        }
        if (prob[h] == 0 || prob[h] == 1) next
        # Type7 sample quantile on an already sorted vector, without a new sort.
        index <- 1 + (n - 1) * prob[h]
        lo <- floor(index)
        s <- z[lo] + (index - lo) * (z[lo + 1L] - z[lo])
        center <- (xq[q] + s)/2
        radius <- abs(xq[q] - s)/2
        fs <- .np_ann_face_window(z, s)
        gt <- .np_ann_face_window(z, center)
        if (is.null(fs) || is.null(gt) || !is.finite(radius) || radius <= 0) next
        ids <- state$order[[j]][seq.int(gt$lo, gt$hi)]
        mark <- marks[ids, q]
        mark[c(1L, length(mark))] <- mark[c(1L, length(mark))]/2
        marked <- (sum(mark)/n)/gt$width
        value <- (if (h == 1L) -1 else 1) * (0.25 * (marked/fs$density))/radius
        if (!is.finite(value)) next
        cut[q, column] <- s
        coef[q, column] <- value
      }
    }
  }
  list(cut = cut, coefficient = coef)
}
