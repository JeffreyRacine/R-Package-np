# Recompute each normalized conditional CDF under the same multinomial counts.
.np_conditional_pairs_block <- function(task, plan) {
  ratio <- function(numerator, weight) {
    total <- crossprod(task$value, weight)
    if (any(!is.finite(total)) || any(total <= 0))
      stop("resampled conditional CDF has zero or nonfinite predictor mass", call. = FALSE)
    crossprod(task$value, numerator) / total[, plan$ix, drop = FALSE]
  }
  sweep(2 * ratio(plan$main, plan$w) - ratio(plan$wide, plan$jw),
        2L, plan$center, "-")
}

# Exact extrema of each signed donor contribution over the response coordinate.
.np_conditional_pairs_range <- function(z) {
  a <- 2 * z$w
  b <- z$jw
  low <- pmin(0 * a, a - b)
  high <- pmax(0 * a, a - b)
  both <- which(a > 0 & b > 0)
  squared <- 4 * (log(a[both]) - log(b[both]) + log(2) / 2)
  index <- both[squared > 0]
  roots <- sqrt(squared[squared > 0])
  for (v in list(-roots, roots)) {
    value <- a[index] * pnorm(v) - b[index] * pnorm(v / sqrt(2))
    low[index] <- pmin(low[index], value)
    high[index] <- pmax(high[index], value)
  }
  span <- colSums(high - low)
  if (any(!is.finite(span)) || any(span <= 0))
    stop("undefined signed conditional CDF contribution range", call. = FALSE)
  q <- sweep(high - low, 2L, span, "/")
  list(lower = colSums(low)[z$ix], span = span[z$ix],
       q = q[, z$ix, drop = FALSE], qmax = apply(q, 2L, max)[z$ix])
}

.np_conditional_pairs_score <- function(x, y, grid, h, B, alpha, original, bytes) {
  n <- length(x)
  m <- nrow(grid)
  bytes <- max(bytes, 8 * (10 * n * m + 2 * n * B + 8 * m * B + 12 * n * 128))
  if (!is.finite(bytes) || bytes > 256 * 1024^2)
    stop("conditional pairs numerical workspace exceeds the 256 MiB budget", call. = FALSE)
  z <- .np_conditional_jack_context(.np_conditional_den_context(x, y, grid, h))
  delta <- outer(y, z$gy, function(a, b) b - a)
  numerator <- function(weight, width) {
    weight[, z$ix, drop = FALSE] * pnorm(delta / width)[, z$iy, drop = FALSE]
  }
  main <- numerator(z$w, h[2L])
  wide <- numerator(z$jw, z$H[2L])
  center <- 2 * colSums(main) - colSums(wide)
  score <- .np_conditional_range_interval(center, .np_conditional_pairs_range(z),
                                          "sharp-score", alpha)
  prob <- rep(1 / n, n)
  tasks <- lapply(seq.int(1L, B, by = 128L), function(first) {
    size <- min(128L, B - first + 1L)
    list(start = first, bsz = size, value = stats::rmultinom(size, n, prob))
  })
  attr(tasks, "rng_final_state") <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  plan <- list(ix = z$ix, w = z$w, jw = z$jw, main = main, wide = wide, center = center)
  roots <- .np_conditional_dispatch(tasks, plan, "pairs-score", m)
  probs <- c(alpha / (2*m), alpha/2, 1-alpha/2, 1-alpha/(2*m))
  cuts <- vapply(seq_len(m), function(i)
    quantile(roots[i, ], probs, type = 1L, names = FALSE), numeric(4L))
  point <- cbind(lower = pmin(center - cuts[3L, ], score[, 1L]),
                 upper = pmax(center - cuts[2L, ], score[, 2L]))
  raw <- unname(quantile(apply(abs(roots), 2L, max), 1-alpha, type = 1L))
  radius <- max(raw, center - point[, 1L], point[, 2L] - center)
  joint <- cbind(lower = center - radius, upper = center + radius)
  bonf <- cbind(lower = pmin(center - cuts[4L, ], joint[, 1L]),
                upper = pmax(center - cuts[1L, ], joint[, 2L]))
  list(center = center, original = original, se = apply(roots, 1L, sd),
       bands = list(pointwise = point, simultaneous = joint, bonferroni = bonf),
       critical = radius, critical.raw = raw, critical.units = "CDF",
       alpha = alpha, B = B, physical.bandwidth = h, grid = grid,
       law = "inid", correction = "pairs-score", calibration = "percentile",
       primary = "fixed-X sharp-score envelope and basic inid pairs companion",
       conditioning = "bandwidths fixed; pairs resample X and Y jointly; marginal score conditions on X",
       endpoint.construction = "jackknife center; hull of sharp-score and signed basic pairs limits; one joint radius",
       scope = "experimental LC conditional CDF, fixed Gaussian2, requested finite family",
       resolution.limited = B*alpha/(2*m) < 1, estimated.numerical.bytes = bytes)
}
