# Evaluate the fitted fixed-bandwidth integrated response kernel at one scale.
.np_conditional_high_integral <- function(y, target, width, kernel, order) {
  data <- data.frame(y = y)
  bw <- kbandwidth.numeric(width, ckertype = kernel, ckerorder = order,
                           nobs = length(y), xdati = untangle(data), xnames = "y")
  t(.np_conditional_literal_density_matrix(bw, data, data.frame(y = target), "integral"))
}

# Convex upper and concave lower support of (1-t)*F(u) - t*F(u/s).
# Chords enclose both support functions between adjacent weight ratios.
.np_conditional_high_envelope <- function(kernel, order, scale) {
  step <- 1 / 512
  extent <- if (kernel == "gaussian") 16 * scale else (sqrt(5) + step) * scale
  knot <- seq.int(0, ceiling(extent / step)) * step
  if (kernel == "epanechnikov")
    knot <- sort(unique(c(knot, sqrt(5), sqrt(5) * scale)))
  response <- function(h) .np_conditional_high_integral(c(0, 1), knot, h, kernel, order)[1L, ]
  main <- response(1)
  wide <- response(scale)
  # On each smooth piece the magnitude of the second derivative is below
  # 512 for all six integrated kernels. The extra 2e-6 covers the rounded
  # native compact support coefficients at the knots and the Gaussian tail.
  allowance <- 512 * step^2 / 8 + 2e-6
  nb <- 256L
  bounds <- vapply((0:nb) / nb, function(t) {
    v <- (1 - t) * main - t * wide
    terminal <- 1 - 2 * t
    r <- range(0, terminal, v)
    c(min(r[1L], terminal - r[2L]) - allowance,
      max(r[2L], terminal - r[1L]) + allowance)
  }, numeric(2L))
  function(a, b) {
    total <- a + b
    pos <- nb * (b / total)
    pos[total == 0] <- 0
    i <- pmin(as.integer(pos), nb - 1L)
    frac <- pos - i
    value <- function(j) total * (bounds[j, i + 1L] * (1 - frac) + bounds[j, i + 2L] * frac)
    list(lower = value(1L), upper = value(2L))
  }
}

.np_conditional_high_pairs_provider <- function(y, targets, h, kernel, order, conditioning = NULL, rms = FALSE) {
  scale <- 2^(1 / order)
  targets <- sort(unique(targets))
  main.y <- .np_conditional_high_integral(y, targets, h, kernel, order)
  wide.y <- .np_conditional_high_integral(y, targets, h * scale, kernel, order)
  envelope <- .np_conditional_high_envelope(kernel, order, scale)
  function(x, y, grid, h, alpha, response.kernel) {
    z <- .np_conditional_pairs_context(x, y, grid, h, conditioning)
    iy <- match(grid$y, targets)
    main <- z$w[, z$ix, drop = FALSE] * main.y[, iy, drop = FALSE]
    wide <- z$jw[, z$ix, drop = FALSE] * wide.y[, iy, drop = FALSE]
    center <- 2 * colSums(main) - colSums(wide)
    bounds <- envelope(2 * z$w, z$jw)
    ranges <- bounds$upper - bounds$lower
    span <- colSums(ranges)
    if (any(!is.finite(span) | span <= 0))
      stop("undefined signed higher-order response range", call. = FALSE)
    q <- sweep(ranges, 2L, span, "/")
    score <- .np_conditional_range_interval(center, list(lower = colSums(bounds$lower)[z$ix],
      span = span[z$ix], q = q[, z$ix, drop = FALSE], qmax = apply(q, 2L, max)[z$ix]),
      "sharp-score", alpha)
    plan <- list(ix = z$ix, w = z$w, jw = z$jw,
                 main = main, wide = wide, center = center)
    if (rms) {
      # Working around an observed donor avoids forming tiny deviations by
      # subtracting conditional probabilities which have rounded to one.
      centered <- function(weights, response) {
        anchor <- max.col(t(weights), ties.method = "first")
        shift <- response[cbind(anchor[z$ix], iy)]
        W <- weights[, z$ix, drop = FALSE]
        deviation <- sweep(response[, iy, drop = FALSE], 2L, shift, "-")
        total <- colSums(W)
        mean <- colSums(W * deviation) / total
        numerator <- W * sweep(deviation, 2L, mean, "-")
        list(numerator = numerator, original = colSums(numerator) / total)
      }
      a <- centered(z$w, main.y)
      b <- centered(z$jw, wide.y)
      plan$stable <- list(main = a$numerator, wide = b$numerator,
                          center = 2 * a$original - b$original)
    }
    list(plan = plan, score = score)
  }
}
