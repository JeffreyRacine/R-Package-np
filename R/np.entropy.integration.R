.np_entropy_points_per_bandwidth <- 4L
.np_entropy_tail_bandwidths <- 8
.np_entropy_workspace_bytes <- 16 * 1024^2

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
  bytes.per.replication <- 80 * support.length
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
      data.null, data.null, counts.x, bw.x, size.x
    )
    density.y <- .np_entropy_count_densities(
      data.null, data.null, counts.y, bw.y, size.y
    )
    summand <- density.y / density.x
    finite <- is.finite(summand)

    for (column in seq_along(index)) {
      if (any(counts.x[, column] > 0 & !finite[, column])) {
        .np_warning(
          " non-finite value in summation-based statistic: integration recommended"
        )
      }
    }

    term <- (1 - sqrt(summand))^2
    term[!finite] <- 0
    denominator <- colSums(counts.x * finite)
    values[index] <- 0.5 * colSums(counts.x * term) / denominator

    unstable <- which(values[index] < 0 | values[index] > 1)
    for (unused in unstable) {
      .np_warning(
        " numerical instability in summation-based statistic: integration recommended"
      )
    }
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

.np_entropy_fixed_density <- function(data, bandwidth) {
  as.numeric(npksum(
    txdat = data,
    bws = bandwidth,
    bandwidth.divide = TRUE
  )$ksum / NROW(data))
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
