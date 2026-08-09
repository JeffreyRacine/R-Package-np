.np_entropy_points_per_bandwidth <- 4L
.np_entropy_tail_bandwidths <- 8
.np_entropy_workspace_bytes <- 16 * 1024^2

.np_entropy_count_chunk_size <- function(support.length,
                                         bytes.per.support,
                                         max.chunk) {
  support.length <- as.double(support.length)
  bytes.per.support <- as.double(bytes.per.support)
  max.chunk <- as.integer(max.chunk)
  if (length(support.length) != 1L || !is.finite(support.length) ||
      support.length < 1 || length(bytes.per.support) != 1L ||
      !is.finite(bytes.per.support) || bytes.per.support <= 0 ||
      length(max.chunk) != 1L || is.na(max.chunk) || max.chunk < 1L) {
    stop("invalid entropy count-chunk dimensions")
  }
  as.integer(max(
    1L,
    min(max.chunk, floor(
      .np_entropy_workspace_bytes / (bytes.per.support * support.length)
    ))
  ))
}

.np_entropy_uses_default_fixed_gaussian <- function(dots) {
  ## Only recognize explicit options that are exactly equivalent to the
  ## default estimator.  Every other density option remains on the general
  ## path; in particular, nearest-neighbour bandwidths must be refitted on the
  ## duplicate-preserving bootstrap sample rather than compressed as weights.
  if (!length(dots))
    return(TRUE)
  dot.names <- names(dots)
  allowed <- c("bwtype", "ckertype", "ckerorder")
  if (is.null(dot.names) || anyNA(dot.names) || any(!nzchar(dot.names)) ||
      anyDuplicated(dot.names) || any(!dot.names %in% allowed))
    return(FALSE)

  if (!is.null(dots$bwtype) &&
      !identical(dots$bwtype, "fixed"))
    return(FALSE)
  if (!is.null(dots$ckertype) &&
      !identical(dots$ckertype, "gaussian"))
    return(FALSE)
  if (!is.null(dots$ckerorder) &&
      (length(dots$ckerorder) != 1L || is.na(dots$ckerorder) ||
       !is.numeric(dots$ckerorder) || dots$ckerorder != 2))
    return(FALSE)
  TRUE
}

.np_entropy_trapezoid_weights <- function(grid) {
  n <- length(grid)
  if (n < 4L || is.unsorted(grid) || anyNA(grid) ||
      !all(is.finite(grid))) {
    stop("entropy quadrature requires at least four ordered finite points")
  }

  dx <- diff(grid)
  weights <- numeric(n)
  weights[1L] <- dx[1L] / 2
  weights[n] <- dx[n - 1L] / 2
  weights[2L:(n - 1L)] <- (dx[-length(dx)] + dx[-1L]) / 2

  ## Match the endpoint correction used by integrate.trapezoidal().
  correction <- dx[1L] / 12
  weights[1L] <- weights[1L] - correction
  weights[2L] <- weights[2L] + correction
  weights[n - 1L] <- weights[n - 1L] + correction
  weights[n] <- weights[n] - correction
  weights
}

.np_entropy_quadrature_axis <- function(
    lower,
    upper,
    bandwidth,
    points.per.bandwidth = .np_entropy_points_per_bandwidth) {
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("entropy integration bandwidths must be finite and positive")
  }
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("entropy integration bounds must be finite and increasing")
  }

  intervals <- ceiling(
    (upper - lower) * points.per.bandwidth / bandwidth
  )
  if (!is.finite(intervals) || intervals > .Machine$integer.max - 1L) {
    stop("entropy quadrature grid is too large for the supplied bandwidth")
  }
  intervals <- max(3L, as.integer(intervals))
  grid <- seq(lower, upper, length.out = intervals + 1L)
  list(grid = grid, weights = .np_entropy_trapezoid_weights(grid))
}

.np_entropy_univariate_domain <- function(data, bandwidths) {
  bandwidths <- as.double(bandwidths)
  if (!length(bandwidths) || anyNA(bandwidths) ||
      any(!is.finite(bandwidths)) || any(bandwidths <= 0)) {
    stop("entropy integration bandwidths must be finite and positive")
  }
  support <- range(data)
  tail <- .np_entropy_tail_bandwidths * max(bandwidths)
  c(support[1L] - tail, support[2L] + tail)
}

.np_entropy_univariate_density <- function(grid, data, bandwidth) {
  as.numeric(npksum(
    txdat = data,
    exdat = grid,
    bws = bandwidth,
    bandwidth.divide = TRUE
  )$ksum / length(data))
}

.np_entropy_univariate_gaussian_summation <- function(data.x,
                                                       data.y,
                                                       bw.x,
                                                       bw.y) {
  density.x <- .np_entropy_univariate_density(data.x, data.x, bw.x)
  density.y <- .np_entropy_univariate_density(data.x, data.y, bw.y)
  summand <- density.y / density.x
  if (!all(is.finite(summand))) {
    .np_warning(
      " non-finite value in summation-based statistic: integration recommended"
    )
    summand <- summand[is.finite(summand)]
  }
  value <- 0.5 * mean((1 - sqrt(summand))^2)
  if (value < 0 || value > 1) {
    .np_warning(
      " numerical instability in summation-based statistic: integration recommended"
    )
  }
  value
}

.np_entropy_univariate_gaussian_integral <- function(data.x,
                                                      data.y,
                                                      bw.x,
                                                      bw.y) {
  domain <- .np_entropy_univariate_domain(c(data.x, data.y), c(bw.x, bw.y))
  axis <- .np_entropy_quadrature_axis(
    domain[1L], domain[2L], min(bw.x, bw.y)
  )
  density.x <- .np_entropy_univariate_density(axis$grid, data.x, bw.x)
  density.y <- .np_entropy_univariate_density(axis$grid, data.y, bw.y)
  unname(crossprod(
    axis$weights,
    0.5 * (sqrt(density.x) - sqrt(density.y))^2
  )[1L])
}

.np_entropy_count_densities <- function(grid,
                                        support,
                                        counts,
                                        bandwidth,
                                        sample.size) {
  value <- npksum(
    txdat = support,
    tydat = counts,
    exdat = grid,
    bws = bandwidth,
    bandwidth.divide = TRUE
  )$ksum
  matrix(value, nrow = length(grid), ncol = ncol(counts)) / sample.size
}

.np_entropy_univariate_iid_bootstrap <- function(data.null,
                                                 size.x,
                                                 size.y,
                                                 bw.x,
                                                 bw.y,
                                                 boot.num,
                                                 progress) {
  domain <- .np_entropy_univariate_domain(data.null, c(bw.x, bw.y))
  axis <- .np_entropy_quadrature_axis(
    domain[1L], domain[2L], min(bw.x, bw.y)
  )
  support.length <- length(data.null)
  bytes.per.replication <- 16 * support.length + 64 * length(axis$grid)
  chunk.size <- as.integer(max(
    1L,
    min(16L, floor(.np_entropy_workspace_bytes / bytes.per.replication))
  ))
  values <- numeric(boot.num)

  for (start in seq.int(1L, boot.num, by = chunk.size)) {
    index <- start:min(boot.num, start + chunk.size - 1L)
    counts.x <- matrix(0, nrow = support.length, ncol = length(index))
    counts.y <- matrix(0, nrow = support.length, ncol = length(index))

    for (column in seq_along(index)) {
      progress <- .np_progress_step(progress, done = index[column])
      counts.x[, column] <- tabulate(
        sample.int(support.length, size.x, replace = TRUE),
        nbins = support.length
      )
      counts.y[, column] <- tabulate(
        sample.int(support.length, size.y, replace = TRUE),
        nbins = support.length
      )
    }

    density.x <- .np_entropy_count_densities(
      axis$grid, data.null, counts.x, bw.x, size.x
    )
    density.y <- .np_entropy_count_densities(
      axis$grid, data.null, counts.y, bw.y, size.y
    )
    values[index] <- as.numeric(crossprod(
      axis$weights,
      0.5 * (sqrt(density.x) - sqrt(density.y))^2
    ))
  }

  list(values = values, progress = progress)
}

.np_entropy_univariate_iid_summation_bootstrap <- function(data.null,
                                                           size.x,
                                                           size.y,
                                                           bw.x,
                                                           bw.y,
                                                           boot.num,
                                                           progress) {
  support.length <- length(data.null)
  chunk.size <- .np_entropy_count_chunk_size(
    support.length, bytes.per.support = 80, max.chunk = 16L
  )
  values <- numeric(boot.num)

  for (start in seq.int(1L, boot.num, by = chunk.size)) {
    index <- start:min(boot.num, start + chunk.size - 1L)
    counts.x <- matrix(0, nrow = support.length, ncol = length(index))
    counts.y <- matrix(0, nrow = support.length, ncol = length(index))

    for (column in seq_along(index)) {
      progress <- .np_progress_step(progress, done = index[column])
      counts.x[, column] <- tabulate(
        sample.int(support.length, size.x, replace = TRUE),
        nbins = support.length
      )
      counts.y[, column] <- tabulate(
        sample.int(support.length, size.y, replace = TRUE),
        nbins = support.length
      )
    }

    values[index] <- .Call(
      "C_np_entropy_univariate_summation_counts",
      as.double(data.null),
      counts.x,
      counts.y,
      as.double(c(bw.x, bw.y)),
      PACKAGE = "np"
    )
  }

  list(values = values, progress = progress)
}

.np_entropy_symmetric_gaussian_integral <- function(data, bandwidth) {
  location <- mean(data)
  half.width <- max(abs(data - location)) +
    .np_entropy_tail_bandwidths * bandwidth
  axis <- .np_entropy_quadrature_axis(
    location - half.width,
    location + half.width,
    bandwidth
  )
  density <- .np_entropy_univariate_density(axis$grid, data, bandwidth)
  unname(crossprod(
    axis$weights,
    0.5 * (sqrt(density) - sqrt(rev(density)))^2
  )[1L])
}

.np_entropy_symmetric_gaussian_summation <- function(data, bandwidth) {
  location <- mean(data)
  density <- .np_entropy_univariate_density(data, data, bandwidth)
  density.rotate <- .np_entropy_univariate_density(
    2 * location - data, data, bandwidth
  )
  summand <- density.rotate / density
  if (!all(is.finite(summand))) {
    .np_warning(
      " non-finite value in summation-based statistic: integration recommended"
    )
    summand <- summand[is.finite(summand)]
  }
  0.5 * mean((1 - sqrt(summand))^2)
}

.np_entropy_tsboot_counts_drawer <- function(n,
                                              B,
                                              blocklen,
                                              sim = c("fixed", "geom"),
                                              n.sim,
                                              endcorr = TRUE) {
  sim <- match.arg(sim)
  n <- as.integer(n)
  B <- as.integer(B)
  blocklen <- as.integer(blocklen)
  n.sim <- as.integer(n.sim)
  if (n < 1L || B < 1L || n.sim < 1L || length(blocklen) != 1L ||
      is.na(blocklen) || blocklen < 1L || blocklen > n) {
    stop("invalid entropy block-bootstrap dimensions")
  }

  ts.array <- utils::getFromNamespace("ts.array", "boot")
  make.ends <- utils::getFromNamespace("make.ends", "boot")
  draws <- ts.array(
    n = n,
    n.sim = n.sim,
    R = B,
    l = blocklen,
    sim = sim,
    endcorr = if (identical(sim, "geom")) TRUE else isTRUE(endcorr)
  )
  starts <- as.matrix(draws$starts)
  lengths <- draws$lengths

  function(start, stopi) {
    start <- as.integer(start)
    stopi <- as.integer(stopi)
    if (start < 1L || stopi < start || stopi > B)
      stop("invalid entropy block-bootstrap chunk bounds")

    replications <- start:stopi
    counts <- matrix(0, nrow = n, ncol = length(replications))
    for (column in seq_along(replications)) {
      replication <- replications[column]
      if (identical(sim, "fixed") && identical(blocklen, 1L)) {
        indices <- as.integer(starts[replication, seq_len(n.sim)])
      } else {
        ends <- if (identical(sim, "geom")) {
          cbind(starts[replication, ], lengths[replication, ])
        } else {
          cbind(starts[replication, ], lengths)
        }
        indices <- apply(ends, 1L, make.ends, n)
        indices <- if (is.list(indices)) {
          as.integer(unlist(indices)[seq_len(n.sim)])
        } else {
          as.integer(indices)[seq_len(n.sim)]
        }
      }
      counts[, column] <- tabulate(indices, nbins = n)
    }
    counts
  }
}

.np_entropy_symmetric_summation_bootstrap <- function(data.null,
                                                       sample.size,
                                                       bandwidth,
                                                       boot.num,
                                                       blocklen,
                                                       sim,
                                                       progress) {
  support.length <- length(data.null)
  chunk.size <- .np_entropy_count_chunk_size(
    support.length, bytes.per.support = 8, max.chunk = 64L
  )
  draw.counts <- .np_entropy_tsboot_counts_drawer(
    n = support.length,
    B = boot.num,
    blocklen = blocklen,
    sim = sim,
    n.sim = sample.size
  )
  values <- numeric(boot.num)

  for (start in seq.int(1L, boot.num, by = chunk.size)) {
    index <- start:min(boot.num, start + chunk.size - 1L)
    counts <- draw.counts(index[1L], index[length(index)])
    values[index] <- .Call(
      "C_np_entropy_symmetric_summation_counts",
      as.double(data.null),
      counts,
      as.double(bandwidth),
      PACKAGE = "np"
    )
    for (done in index)
      progress <- .np_progress_step(progress, done = done)
  }

  list(values = values, progress = progress)
}

.np_entropy_bivariate_gaussian_summation <- function(x.dat,
                                                      y.dat,
                                                      bw.x,
                                                      bw.y,
                                                      bw.joint) {
  value <- .Call(
    "C_np_entropy_bivariate_summation",
    as.double(x.dat),
    as.double(y.dat),
    as.double(c(bw.x, bw.y, bw.joint)),
    PACKAGE = "np"
  )
  if (!is.finite(value)) {
    .np_warning(
      " non-finite value in summation-based statistic: integration recommended"
    )
  }
  value
}

.np_entropy_bivariate_gaussian_summation_xindex <- function(x.dat,
                                                             y.dat,
                                                             index,
                                                             bw.x,
                                                             bw.y,
                                                             bw.joint) {
  value <- .Call(
    "C_np_entropy_bivariate_summation_xindex",
    as.double(x.dat),
    as.double(y.dat),
    index,
    as.double(c(bw.x, bw.y, bw.joint)),
    PACKAGE = "np"
  )
  if (any(!is.finite(value))) {
    .np_warning(
      " non-finite value in summation-based statistic: integration recommended"
    )
  }
  value
}

.np_entropy_bivariate_domain <- function(x.dat,
                                         y.dat,
                                         bw.x,
                                         bw.y,
                                         bw.joint) {
  tail.x <- .np_entropy_tail_bandwidths * max(bw.x, bw.joint[1L])
  tail.y <- .np_entropy_tail_bandwidths * max(bw.y, bw.joint[2L])
  list(
    lower = c(min(x.dat) - tail.x, min(y.dat) - tail.y),
    upper = c(max(x.dat) + tail.x, max(y.dat) + tail.y)
  )
}

.np_entropy_bivariate_integral <- function(x.dat,
                                           y.dat,
                                           bw.x,
                                           bw.y,
                                           bw.joint,
                                           lower,
                                           upper) {
  x.dat <- as.double(x.dat)
  y.dat <- as.double(y.dat)
  bandwidths <- as.double(c(bw.x, bw.y, bw.joint))
  if (length(x.dat) < 1L || length(y.dat) != length(x.dat)) {
    stop("entropy integration data must have equal positive lengths")
  }
  if (length(bandwidths) != 4L || anyNA(bandwidths) ||
      any(!is.finite(bandwidths)) || any(bandwidths <= 0)) {
    stop("entropy integration bandwidths must be finite and positive")
  }

  axis.x <- .np_entropy_quadrature_axis(
    lower[1L], upper[1L], min(bandwidths[c(1L, 3L)]),
    points.per.bandwidth = 2L
  )
  axis.y <- .np_entropy_quadrature_axis(
    lower[2L], upper[2L], min(bandwidths[c(2L, 4L)]),
    points.per.bandwidth = 2L
  )
  max.points <- as.integer(max(
    1L, floor(.np_entropy_workspace_bytes / 64)
  ))
  x.block.size <- as.integer(min(length(axis.x$grid), max.points))
  total <- 0
  compensation <- 0

  for (x.start in seq.int(1L, length(axis.x$grid), by = x.block.size)) {
    x.index <- x.start:min(
      length(axis.x$grid), x.start + x.block.size - 1L
    )
    y.block.size <- as.integer(max(
      1L, floor(max.points / length(x.index))
    ))
    for (y.start in seq.int(1L, length(axis.y$grid), by = y.block.size)) {
      y.index <- y.start:min(
        length(axis.y$grid), y.start + y.block.size - 1L
      )
      points <- rbind(
        rep(axis.x$grid[x.index], times = length(y.index)),
        rep(axis.y$grid[y.index], each = length(x.index))
      )
      weights <-
        rep(axis.x$weights[x.index], times = length(y.index)) *
        rep(axis.y$weights[y.index], each = length(x.index))
      values <- .Call(
        "C_np_entropy_gaussian_integrand",
        points,
        x.dat,
        y.dat,
        bandwidths,
        PACKAGE = "np"
      )
      block.value <- sum(weights * values)
      adjusted <- block.value - compensation
      next.total <- total + adjusted
      compensation <- (next.total - total) - adjusted
      total <- next.total
    }
  }

  0.5 * total
}
