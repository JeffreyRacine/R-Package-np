# Recompute each normalized conditional CDF under the same multinomial counts.
.np_conditional_pairs_block <- function(task, plan) {
  denominator <- function(weight) {
    total <- crossprod(task$value, weight)
    if (any(!is.finite(total)) || any(total <= 0))
      stop("resampled conditional CDF has zero or nonfinite predictor mass", call. = FALSE)
    total[, plan$ix, drop = FALSE]
  }
  stable <- plan[["stable"]]
  if (is.null(stable)) {
    # Keep only one denominator live at a time on the ordinary path.
    ratio <- function(numerator, weight) {
      mass <- denominator(weight)
      crossprod(task$value, numerator) / mass
    }
    return(sweep(2 * ratio(plan$main, plan$w) - ratio(plan$wide, plan$jw),
                 2L, plan$center, "-"))
  }
  main.mass <- denominator(plan$w)
  wide.mass <- denominator(plan$jw)
  roots <- function(data) {
    sweep(2 * (crossprod(task$value, data$main) / main.mass) -
          crossprod(task$value, data$wide) / wide.mass, 2L, data$center, "-")
  }
  cbind(roots(plan), roots(stable))
}

# Integrated unit-variance second-order Epanechnikov kernel.
.np_conditional_epan_cdf <- function(u) {
  t <- u / sqrt(5)
  t[] <- pmin(1, pmax(-1, t))
  tail <- .25 * (1 - abs(t))^2 * (2 + abs(t))
  ifelse(t < 0, tail, 1 - tail)
}

# Exact extrema of each signed donor contribution over the response coordinate.
.np_conditional_pairs_range <- function(z, response.kernel = "gaussian") {
  a <- 2 * z$w
  b <- z$jw
  low <- pmin(0 * a, a - b)
  high <- pmax(0 * a, a - b)
  if (identical(response.kernel, "uniform")) {
    if (!isTRUE(min(a, b) >= 0))
      stop("uniform donor ranges require nonnegative conditioning weights", call. = FALSE)
    # A signed difference of two uniform distribution functions is linear
    # between its knots. The outer knots are already in low and high.
    value <- a - b * (.5 + .5 / sqrt(2))
    reflected <- a - b - value
    low <- pmin(low, value, reflected)
    high <- pmax(high, value, reflected)
  } else if (identical(response.kernel, "epanechnikov")) {
    if (!isTRUE(min(a, b) >= 0))
      stop("Epanechnikov donor ranges require nonnegative conditioning weights", call. = FALSE)
    # On the common support the derivative is linear in the squared argument.
    # From the inner to the outer knot it is nonpositive. Only the positive
    # crossing and its reflection can enlarge the range beyond its endpoints.
    scale <- sqrt(2)
    index <- which(b > 0 & a > b / scale)
    numerator <- a[index] - b[index] / scale
    denominator <- a[index] - b[index] / scale^3
    root <- sqrt(5 * numerator / denominator)
    value <- a[index] * .np_conditional_epan_cdf(root) -
      b[index] * .np_conditional_epan_cdf(root / scale)
    reflected <- (a[index] - b[index]) - value
    low[index] <- pmin(low[index], value, reflected)
    high[index] <- pmax(high[index], value, reflected)
  } else {
    both <- which(a > 0 & b > 0)
    squared <- 4 * (log(a[both]) - log(b[both]) + log(2) / 2)
    index <- both[squared > 0]
    roots <- sqrt(squared[squared > 0])
    for (v in list(-roots, roots)) {
      value <- a[index] * pnorm(v) - b[index] * pnorm(v / sqrt(2))
      low[index] <- pmin(low[index], value)
      high[index] <- pmax(high[index], value)
    }
  }
  span <- colSums(high - low)
  if (any(!is.finite(span)) || any(span <= 0))
    stop("undefined signed conditional CDF contribution range", call. = FALSE)
  q <- sweep(high - low, 2L, span, "/")
  list(lower = colSums(low)[z$ix], span = span[z$ix],
       q = q[, z$ix, drop = FALSE], qmax = apply(q, 2L, max)[z$ix])
}

# Evaluate each conditioning row at the fitted and widened physical bandwidth.
.np_conditional_beta_conditioning <- function(x, targets, h, bounds) {
  if (length(bounds) != 2L || any(!is.finite(bounds)) || bounds[1L] >= bounds[2L] ||
      any(x < bounds[1L] | x > bounds[2L]) || any(targets < bounds[1L] | targets > bounds[2L]))
    stop("beta conditioning requires data and targets inside known finite support", call. = FALSE)
  targets <- unique(targets)
  train <- data.frame(x = x)
  evaluation <- data.frame(x = targets)
  geometry <- untangle(train)
  weight <- function(width) {
    bw <- kbandwidth.numeric(width, ckertype = "beta", ckerorder = 2L,
                             ckerbound = "fixed", ckerlb = bounds[1L], ckerub = bounds[2L],
                             nobs = length(x), xdati = geometry, xnames = "x")
    values <- .np_conditional_literal_density_matrix(bw, train, evaluation, "normal")
    mass <- rowSums(values)
    if (any(!is.finite(mass)) || any(mass <= 0) || any(values < 0))
      stop("undefined beta conditioning mass", call. = FALSE)
    t(values / mass)
  }
  list(targets = targets, main = weight(h), wide = weight(sqrt(2) * h))
}

.np_conditional_pairs_context <- function(x, y, grid, h, conditioning) {
  if (is.null(conditioning))
    return(.np_conditional_jack_context(.np_conditional_den_context(x, y, grid, h, pairs.score = TRUE)))
  gx <- unique(grid$x)
  group <- match(gx, conditioning$targets)
  list(ix = match(grid$x, gx), w = conditioning$main[, group, drop = FALSE],
       jw = conditioning$wide[, group, drop = FALSE])
}

# The empirical event is shared by all three response kernels. Process each
# conditioning coordinate once after all pairs tiles have been reduced.
.np_conditional_local_wilson <- function(x, y, grid, bandwidth, alpha, conditioning = NULL) {
  ord <- order(y)
  end <- findInterval(grid$y, y[ord]) + 1L
  level <- qnorm(1 - alpha / 2)
  limits <- matrix(0, nrow(grid), 2L)
  groups <- split(seq_len(nrow(grid)), match(grid$x, unique(grid$x)))
  for (ids in groups) {
    target <- grid$x[ids[1L]]
    weight <- if (is.null(conditioning)) dnorm((target - x[ord]) / bandwidth)
              else conditioning$main[ord, match(target, conditioning$targets)]
    cumulative <- c(0, cumsum(weight))
    total <- cumulative[length(cumulative)]
    probability <- cumulative[end[ids]] / total
    scale <- level^2 * sum((weight / total)^2)
    half <- sqrt(scale * (probability * (1 - probability) + scale / 4))
    middle <- probability + scale / 2
    limits[ids, ] <- cbind(middle - half, middle + half) / (1 + scale)
  }
  limits
}

# Build the two response tables once. Workers need only the usual donor plan.
.np_conditional_beta_pairs_provider <- function(response, targets, bandwidth, bounds, conditioning = NULL) {
  width <- diff(bounds)
  if (length(bounds) != 2L || any(!is.finite(bounds)) || !is.finite(width) || width <= 0 ||
      any(response < bounds[1L] | response > bounds[2L]) ||
      any(targets < bounds[1L] | targets > bounds[2L]))
    stop("beta pairs inference requires data and targets inside known finite support", call. = FALSE)
  targets <- sort(unique(targets))
  donor <- (response - bounds[1L]) / width
  donor.upper <- (bounds[2L] - response) / width
  left <- (targets - bounds[1L]) / width
  right <- (bounds[2L] - targets) / width
  edge <- seq(0, 1, length.out = 65L)
  primitive <- function(s) {
    precision <- (width / s)^2
    if (!is.finite(precision))
      stop("beta pairs bandwidth exceeds representable precision", call. = FALSE)
    shape1 <- 1 + precision * c(donor, edge)
    shape2 <- 1 + precision * c(donor.upper, rev(edge))
    values <- vapply(seq_along(targets), function(j) {
      if (left[j] <= .5) pbeta(left[j], shape1, shape2)
      else pbeta(right[j], shape2, shape1, lower.tail = FALSE)
    }, numeric(length(shape1)))
    if (any(!is.finite(values)))
      stop("beta pairs response integrals must be finite", call. = FALSE)
    values
  }
  fitted <- primitive(bandwidth)
  wider <- primitive(sqrt(2) * bandwidth)
  obs <- seq_along(response)
  edges <- length(response) + seq_along(edge)
  bins <- seq_len(length(edge) - 1L)
  function(x, y, grid, h, alpha, response.kernel) {
    z <- .np_conditional_pairs_context(x, y, grid, h, conditioning)
    iy <- match(grid$y, targets)
    main <- z$w[, z$ix, drop = FALSE] * fitted[obs, iy, drop = FALSE]
    wide <- z$jw[, z$ix, drop = FALSE] * wider[obs, iy, drop = FALSE]
    center <- 2 * colSums(main) - colSums(wide)
    score <- matrix(0, length(center), 2L)
    for (j in unique(iy)) {
      idx <- which(iy == j)
      if (targets[j] == bounds[1L] || targets[j] == bounds[2L]) {
        score[idx, ] <- as.numeric(targets[j] == bounds[2L])
        next
      }
      a <- 2 * z$w[, z$ix[idx], drop = FALSE]
      b <- z$jw[, z$ix[idx], drop = FALSE]
      first <- fitted[edges, j]
      second <- wider[edges, j]
      # Each primitive decreases in its donor coordinate. Opposite ends of
      # each bin therefore enclose the signed contribution at every donor.
      lower <- upper <- a
      lower[] <- Inf
      upper[] <- -Inf
      for (i in bins) {
        lower <- pmin(lower, a * first[i + 1L] - b * second[i])
        upper <- pmax(upper, a * first[i] - b * second[i + 1L])
      }
      allowance <- 2e-13 * (abs(lower) + abs(upper) + a + b)
      lower <- lower - allowance
      ranges <- upper + allowance - lower
      span <- colSums(ranges)
      if (any(!is.finite(span) | span <= 0))
        stop("undefined beta signed donor range", call. = FALSE)
      q <- sweep(ranges, 2L, span, "/")
      score[idx, ] <- .np_conditional_range_interval(center[idx], list(
        lower = colSums(lower), span = span, q = q, qmax = apply(q, 2L, max)),
        "sharp-score", alpha)
    }
    list(plan = list(ix = z$ix, w = z$w, jw = z$jw,
                     main = main, wide = wide, center = center), score = score)
  }
}

.np_conditional_pairs_prepare <- function(x, y, grid, h, alpha, response.kernel = "gaussian") {
  z <- .np_conditional_jack_context(.np_conditional_den_context(x, y, grid, h, pairs.score = TRUE))
  delta <- outer(y, z$gy, function(a, b) b - a)
  primitive <- switch(response.kernel, gaussian = pnorm, epanechnikov = .np_conditional_epan_cdf,
    uniform = function(u) stats::punif(u, -1, 1),
    stop("unsupported response kernel for conditional CDF pairs-score", call. = FALSE))
  numerator <- function(weight, width) {
    weight[, z$ix, drop = FALSE] * primitive(delta / width)[, z$iy, drop = FALSE]
  }
  main <- numerator(z$w, h[2L])
  wide <- numerator(z$jw, z$H[2L])
  center <- 2 * colSums(main) - colSums(wide)
  score <- .np_conditional_range_interval(center, .np_conditional_pairs_range(z, response.kernel),
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

.np_conditional_pairs_support <- function(x, targets, h, tasks = NULL, conditioning = NULL) {
  for (scale in seq_len(2L)) {
    width <- if (scale == 1L) h else sqrt(2) * h
    for (target in targets) {
      if (is.null(conditioning)) w <- dnorm((target - x) / width)
      else {
        weights <- if (scale == 1L) conditioning$main else conditioning$wide
        w <- weights[, match(target, conditioning$targets)]
      }
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

.np_conditional_pairs_score <- function(x, y, grid, h, B, alpha, original, bytes, tile.size = NULL, response.kernel = "gaussian", pointwise = NULL, response.bounds = NULL, response.order = 2L, conditioning.bounds = NULL, profile.weights = NULL, simultaneous = NULL) {
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
    if (is.null(conditioning.bounds) && is.null(profile.weights)) .np_conditional_pairs_support(x, targets, h[1L])
  }

  # The explicit standardized maximum needs a second numerical root per
  # target and draw. Both use the same counts, denominator and worker task.
  if (identical(simultaneous, "rms-envelope") && !is.null(profile.weights)) {
    cap <- 256 * 1024^2
    nx <- length(unique(grid$x))
    ny <- length(unique(grid$y))
    reserved <- 64 * B + 32 * m + 8 * (12 * n * nx + 8 * n * ny + 24 * n + 32 * 10000)
    if (identical(pointwise, "local-wilson")) reserved <- reserved + 128 * (n + m)
    fixed <- 2 * as.double(n) * B + 12 * n * 128 + 256 * n + 16 * m + 8 * B
    # Per target, four simultaneous copies of the two-channel draw output
    # bound collection and transposition; donor setup is independent of B.
    per.target <- 48 * n + 8 * B
    size <- min(tile.size, floor(((cap - reserved) / 8 - fixed) / per.target))
    if (size < 1L || m > 12 * size)
      stop("RMS envelope exceeds the 256 MiB numerical workspace or transport limit", call. = FALSE)
    tile.size <- as.integer(size)
    tiled <- tiled || size < m
    bytes <- max(bytes, 8 * (fixed + size * per.target))
  }
  if (identical(pointwise, "local-wilson")) {
    # Reserve sorting, cumulative weights and the second set of endpoint arrays.
    bytes <- bytes + 128 * (as.double(n) + m)
    if (bytes > 256 * 1024^2)
      stop("local Wilson numerical workspace exceeds the 256 MiB budget", call. = FALSE)
  }
  rms.envelope <- identical(simultaneous, "rms-envelope")
  if (!is.null(simultaneous) && (!rms.envelope || is.null(profile.weights)))
    stop("unsupported pairs-score simultaneous calculation", call. = FALSE)
  if (rms.envelope) {
    bytes <- bytes + 64 * as.double(B) + 32 * m
    if (bytes > 256 * 1024^2)
      stop("RMS envelope exceeds the 256 MiB numerical workspace", call. = FALSE)
    draw.rms.max <- numeric(B)
    root.rms <- numeric(m)
  }
  conditioning <- NULL
  if (!is.null(conditioning.bounds) || !is.null(profile.weights)) {
    targets <- unique(grid$x)
    bytes <- bytes + 8 * (12 * as.double(n) * length(targets) + 12 * n)
    if (bytes > 256 * 1024^2)
      stop(if (is.null(profile.weights)) "beta conditioning numerical workspace exceeds the 256 MiB budget" else "conditioning profiles exceed the 256 MiB budget", call. = FALSE)
    conditioning <- if (is.null(profile.weights)) .np_conditional_beta_conditioning(x, targets, h[1L], conditioning.bounds) else profile.weights()
  }
  if (response.order > 2L || (!is.null(conditioning) && identical(response.kernel, "gaussian"))) {
    targets.y <- unique(grid$y)
    # Two native response tables, their setup and the unit-scale support grid.
    bytes <- bytes + 8 * (8 * as.double(n) * length(targets.y) + 12 * n + 32 * 10000)
    if (bytes > 256 * 1024^2)
      stop("higher-order pairs numerical workspace exceeds the 256 MiB budget", call. = FALSE)
    prepare <- .np_conditional_high_pairs_provider(y, targets.y, h[2L], response.kernel, response.order, conditioning, rms = rms.envelope)
  } else if (identical(response.kernel, "beta")) {
    targets.y <- unique(grid$y)
    # Two response tables plus their construction, workspace and outward bounds.
    bytes <- bytes + 8 * (8 * (as.double(n) + 65) * length(targets.y) + 12 * n)
    if (bytes > 256 * 1024^2)
      stop("beta pairs numerical workspace exceeds the 256 MiB budget", call. = FALSE)
    prepare <- .np_conditional_beta_pairs_provider(y, targets.y, h[2L], response.bounds, conditioning)
  } else if (is.null(conditioning)) {
    prepare <- .np_conditional_pairs_prepare
  } else {
    stop("no qualified response range for the supplied conditioning weights", call. = FALSE)
  }
  probs <- c(alpha / (2*m), alpha/2, 1-alpha/2, 1-alpha/(2*m))
  if (!is.null(profile.weights)) original <- numeric(m)
  center <- se <- numeric(m)
  point <- matrix(NA_real_, m, 2L, dimnames = list(NULL, c("lower", "upper")))
  cuts <- matrix(NA_real_, 4L, m)
  draw.max <- numeric(B)
  for (first in seq.int(1L, m, by = tile.size)) {
    idx <- seq.int(first, min(m, first + tile.size - 1L))
    tile <- prepare(x, y, grid[idx, , drop = FALSE], h, alpha, response.kernel)
    if (first == 1L) {
      # Every target sees exactly the same coordinator-generated pairs draws.
      prob <- rep(1 / n, n)
      tasks <- lapply(seq.int(1L, B, by = 128L), function(first) {
        size <- min(128L, B - first + 1L)
        list(start = first, bsz = size, value = stats::rmultinom(size, n, prob))
      })
      attr(tasks, "rng_final_state") <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    if (first == 1L && (tiled || !is.null(conditioning)))
      .np_conditional_pairs_support(x, targets, h[1L], tasks, conditioning)
    if (tiled) {
      tile$plan$tile.label <- sprintf("Conditional inference, target tile %d of %d",
                                     (first - 1L) %/% tile.size + 1L, ceiling(m / tile.size))
    }
    roots <- .np_conditional_dispatch(tasks, tile$plan, "pairs-score",
                                      length(idx) * if (rms.envelope) 2L else 1L)
    if (!is.null(profile.weights)) original[idx] <- colSums(tile$plan$main)
    center[idx] <- tile$plan$center
    if (tiled || rms.envelope) {
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
    if (rms.envelope) {
      for (i in seq_along(idx)) {
        absolute <- abs(roots[length(idx) + i, ])
        peak <- max(absolute)
        if (peak == 0) next
        # The bootstrap roots are bounded. Relative scaling also protects tiny roots.
        relative <- absolute / peak
        scale <- sqrt(mean(relative^2))
        root.rms[idx[i]] <- peak * scale
        if (root.rms[idx[i]] == 0)
          stop("positive bootstrap RMS is below numerical resolution", call. = FALSE)
        draw.rms.max <- pmax(draw.rms.max, relative / scale)
      }
    }
    rm(roots, tile)
  }
  raw <- unname(quantile(draw.max, 1-alpha, type = 1L))
  radius <- max(raw, center - point[, 1L], point[, 2L] - center)
  joint <- cbind(lower = center - radius, upper = center + radius)
  bonf <- cbind(lower = pmin(center - cuts[4L, ], joint[, 1L]),
                upper = pmax(center - cuts[1L, ], joint[, 2L]))
  ans <- list(center = center, original = original, se = se,
       bands = list(pointwise = point, simultaneous = joint, bonferroni = bonf),
       critical = radius, critical.raw = raw, critical.units = "CDF",
       alpha = alpha, B = B, physical.bandwidth = h, grid = grid,
       law = "inid", correction = "pairs-score", calibration = "percentile",
       primary = "fixed-X sharp-score envelope and basic inid pairs companion",
       conditioning = "bandwidths fixed; pairs resample X and Y jointly; marginal score conditions on X",
       endpoint.construction = "jackknife center; hull of sharp-score and signed basic pairs limits; one joint radius",
       scope = if (!is.null(profile.weights))
         "experimental LC conditional CDF, two fixed-bandwidth predictors, Gaussian2 response; requested panel"
       else if (!is.null(conditioning))
         sprintf("experimental LC conditional CDF, fixed known-support beta2 X and %s2 Y, requested finite family", response.kernel)
       else if (response.order > 2L)
         sprintf("experimental LC conditional CDF, fixed Gaussian2 X and %s order %d Y, requested finite family", response.kernel, response.order)
       else if (identical(response.kernel, "beta"))
         "experimental LC conditional CDF, fixed Gaussian2 X and known-support beta2 Y, requested finite family"
       else if (identical(response.kernel, "uniform"))
         "experimental LC conditional CDF, fixed Gaussian2 X and uniform Y, requested finite family"
       else if (identical(response.kernel, "epanechnikov"))
         "experimental LC conditional CDF, fixed Gaussian2 X and Epanechnikov2 Y, requested finite family"
       else "experimental LC conditional CDF, fixed Gaussian2, requested finite family",
       resolution.limited = B*alpha/(2*m) < 1, estimated.numerical.bytes = bytes)
  if (identical(response.kernel, "beta"))
    ans$range.construction <- "conservative enclosure from monotone beta response intervals"
  else if (response.order > 2L || !is.null(conditioning))
    ans$range.construction <- "conservative chord enclosure of signed donor extrema"
  if (identical(pointwise, "local-wilson")) {
    local <- .np_conditional_local_wilson(x, y, grid, h[1L], alpha, conditioning)
    local[, 1L] <- pmin(center, center - cuts[3L, ], local[, 1L])
    local[, 2L] <- pmax(center, center - cuts[2L, ], local[, 2L])
    extended <- local[, 1L] < joint[, 1L] | local[, 2L] > joint[, 2L]
    joint[, 1L] <- pmin(joint[, 1L], local[, 1L])
    joint[, 2L] <- pmax(joint[, 2L], local[, 2L])
    bonf[, 1L] <- pmin(bonf[, 1L], joint[, 1L])
    bonf[, 2L] <- pmax(bonf[, 2L], joint[, 2L])
    ans$bands <- list(pointwise = local, simultaneous = joint, bonferroni = bonf)
    ans$pointwise.calibration <- "local-wilson"
    ans$nested.coordinates <- which(extended)
    ans$primary <- "positive local Wilson and basic pairs pointwise interval"
    ans$conditioning <- "fixed bandwidths; joint pairs resampling; Wilson uses local homogeneity"
    ans$endpoint.construction <- "local Wilson and basic pairs pointwise hull; original sharp-score joint radius with outward nesting"
  }
  if (rms.envelope) {
    multiplier <- unname(quantile(draw.rms.max, 1-alpha, type = 1L))
    width <- multiplier * root.rms
    joint[, 1L] <- pmin(joint[, 1L], center - width)
    joint[, 2L] <- pmax(joint[, 2L], center + width)
    bonf[, 1L] <- pmin(bonf[, 1L], joint[, 1L])
    bonf[, 2L] <- pmax(bonf[, 2L], joint[, 2L])
    ans$bands$simultaneous <- joint
    ans$bands$bonferroni <- bonf
    ans$settings <- list(critical.components = list(raw.radius = radius,
      rms.multiplier = multiplier, root.rms = root.rms))
    ans$endpoint.construction <- paste(ans$endpoint.construction,
      "joint hull also includes the maximum of RMS-scaled pairs deviations", sep = ". ")
  }
  if (!is.null(ans$range.construction))
    ans$endpoint.construction <- paste(ans$endpoint.construction, ans$range.construction, sep = ". ")
  ans
}
