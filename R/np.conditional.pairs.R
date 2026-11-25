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

.np_conditional_pairs_prepare <- function(x, y, grid, h, alpha) {
  z <- .np_conditional_jack_context(.np_conditional_den_context(x, y, grid, h, pairs.score = TRUE))
  delta <- outer(y, z$gy, function(a, b) b - a)
  numerator <- function(weight, width) {
    weight[, z$ix, drop = FALSE] * pnorm(delta / width)[, z$iy, drop = FALSE]
  }
  main <- numerator(z$w, h[2L])
  wide <- numerator(z$jw, z$H[2L])
  center <- 2 * colSums(main) - colSums(wide)
  score <- .np_conditional_range_interval(center, .np_conditional_pairs_range(z),
                                          "sharp-score", alpha)
  plan <- list(ix = z$ix, w = z$w, jw = z$jw, main = main, wide = wide, center = center)
  list(plan = plan, score = score)
}

.np_conditional_pairs_allocation <- function(n, m, B) {
  limit <- 256 * 1024^2
  # Four complete integer count banks allow for the coordinator list and its
  # transport copies. Coercion and worker temporaries use 128-draw blocks.
  fixed <- 2 * n * B + 12 * n * 128 + 256 * n + 16 * m + 8 * B
  per.target <- 24 * n + 6 * B
  tile <- min(m, 250L, floor((limit / 8 - fixed) / per.target))
  center.bytes <- 8 * (24 * n * m + 256 * n + 16 * m)
  # Bound repeated count transmissions as well as simultaneously live storage.
  if (tile < 1 || m > 12 * tile || center.bytes > limit) return(NULL)
  list(B = B, estimated.bytes = max(center.bytes, 8 * (fixed + tile * per.target)),
       tile.size = as.integer(tile))
}

.np_conditional_pairs_support <- function(x, targets, h, tasks = NULL) {
  for (width in c(h, sqrt(2) * h)) {
    for (target in targets) {
      w <- dnorm((target - x) / width)
      mass <- sum(w)
      if (!is.finite(mass) || mass <= 0) stop("undefined predictor mass", call. = FALSE)
      if (is.null(tasks)) next
      w <- w / mass
      # If every observation contributes a normal positive number, every
      # multinomial resample has positive mass, regardless of multiplicities.
      if (all(w >= .Machine$double.xmin)) next
      for (task in tasks) {
        total <- crossprod(task$value, w)
        if (any(!is.finite(total) | total <= 0))
          stop("resampled conditional CDF has zero or nonfinite predictor mass", call. = FALSE)
      }
    }
  }
  invisible(NULL)
}

.np_conditional_pairs_score <- function(x, y, grid, h, B, alpha, original, bytes, tile.size = NULL) {
  n <- length(x)
  m <- nrow(grid)
  tiled <- !is.null(tile.size)
  if (!tiled) {
    bytes <- max(bytes, 8 * (10 * n * m + 2 * n * B + 8 * m * B + 12 * n * 128))
    if (!is.finite(bytes) || bytes > 256 * 1024^2)
      stop("conditional pairs numerical workspace exceeds the 256 MiB budget", call. = FALSE)
    tile.size <- m
  } else {
    allocation <- .np_conditional_pairs_allocation(n, m, B)
    if (length(tile.size) != 1L || !is.numeric(tile.size) || !is.finite(tile.size) ||
        tile.size < 1 || tile.size != floor(tile.size) || B > 9999L ||
        is.null(allocation) || tile.size > allocation$tile.size || m > 12 * tile.size)
      stop("invalid or over-budget conditional pairs tile size", call. = FALSE)
    targets <- unique(grid$x)
    .np_conditional_pairs_support(x, targets, h[1L])
  }

  probs <- c(alpha / (2*m), alpha/2, 1-alpha/2, 1-alpha/(2*m))
  center <- se <- numeric(m)
  point <- matrix(NA_real_, m, 2L, dimnames = list(NULL, c("lower", "upper")))
  cuts <- matrix(NA_real_, 4L, m)
  draw.max <- numeric(B)
  for (first in seq.int(1L, m, by = tile.size)) {
    idx <- seq.int(first, min(m, first + tile.size - 1L))
    tile <- .np_conditional_pairs_prepare(x, y, grid[idx, , drop = FALSE], h, alpha)
    if (first == 1L) {
      # Every target sees exactly the same coordinator-generated pairs draws.
      prob <- rep(1 / n, n)
      tasks <- lapply(seq.int(1L, B, by = 128L), function(first) {
        size <- min(128L, B - first + 1L)
        list(start = first, bsz = size, value = stats::rmultinom(size, n, prob))
      })
      attr(tasks, "rng_final_state") <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    if (tiled) {
      if (first == 1L) .np_conditional_pairs_support(x, targets, h[1L], tasks)
      tile$plan$tile.label <- sprintf("Conditional inference, target tile %d of %d",
                                     (first - 1L) %/% tile.size + 1L, ceiling(m / tile.size))
    }
    roots <- .np_conditional_dispatch(tasks, tile$plan, "pairs-score", length(idx))
    center[idx] <- tile$plan$center
    if (tiled) {
      for (i in seq_along(idx)) {
        v <- roots[i, ]
        se[idx[i]] <- sd(v)
        cuts[, idx[i]] <- quantile(v, probs, type = 1L, names = FALSE)
        draw.max <- pmax(draw.max, abs(v))
      }
    } else {
      se[idx] <- apply(roots, 1L, sd)
      cuts[, idx] <- vapply(seq_along(idx), function(i)
        quantile(roots[i, ], probs, type = 1L, names = FALSE), numeric(4L))
      draw.max <- pmax(draw.max, apply(abs(roots), 2L, max))
    }
    point[idx, 1L] <- pmin(center[idx] - cuts[3L, idx], tile$score[, 1L])
    point[idx, 2L] <- pmax(center[idx] - cuts[2L, idx], tile$score[, 2L])
    rm(roots, tile)
  }
  raw <- unname(quantile(draw.max, 1-alpha, type = 1L))
  radius <- max(raw, center - point[, 1L], point[, 2L] - center)
  joint <- cbind(lower = center - radius, upper = center + radius)
  bonf <- cbind(lower = pmin(center - cuts[4L, ], joint[, 1L]),
                upper = pmax(center - cuts[1L, ], joint[, 2L]))
  list(center = center, original = original, se = se,
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
