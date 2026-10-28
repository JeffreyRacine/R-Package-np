# Fixed Gaussian2 numeric CDF contributions have this certified signed range.
# The score interval targets their expectation at fixed bandwidth and target;
# smoothing bias and data-selected targets require separate coverage evidence.
.np_udist_marginal_range <- function(q) {
  t <- sqrt(1 + 2/q)
  a <- pnorm(-t) - (q/2) * t * dnorm(t)
  c(lower = a, upper = 1 - a)
}

.np_udist_marginal_score <- function(center, n, q, alpha) {
  bounds <- .np_udist_marginal_range(q)
  span <- bounds[2L] - bounds[1L]
  p <- (center - bounds[1L])/span
  if (any(!is.finite(p)) || any(p < 0 | p > 1))
    stop("CDF center outside certified contribution range", call. = FALSE)
  z <- qnorm(alpha/2, lower.tail = FALSE)
  denominator <- 1 + z*z/n
  midpoint <- (p + z*z/(2*n))/denominator
  radius <- z * sqrt(p*(1-p)/n + z*z/(4*n*n))/denominator
  cbind(lower = bounds[1L] + span*(midpoint-radius),
        upper = bounds[1L] + span*(midpoint+radius))
}

# Preserve the calibrated process and RNG. Only the deterministic marginal
# envelope and the minimum uniform outward nesting envelope change endpoints.
.np_udist_marginal_envelope <- function(answer, widths, n, q, alpha) {
  center <- answer$center
  m <- length(center)
  answer$reference.bands <- answer$intervals
  answer$reference.critical <- answer$critical
  expansion <- numeric(length(alpha))
  for (k in seq_along(alpha)) {
    score <- .np_udist_marginal_score(center, n, q, alpha[k])
    point <- cbind(pmin(center-widths$pointwise[, k], score[, 1L]),
                   pmax(center+widths$pointwise[, k], score[, 2L]))
    joint <- max(widths$simultaneous[, k], center-point[, 1L],
                 point[, 2L]-center)
    simultaneous <- cbind(center-joint, center+joint)
    bonferroni <- cbind(pmin(center-widths$bonferroni[, k], simultaneous[, 1L]),
                        pmax(center+widths$bonferroni[, k], simultaneous[, 2L]))
    if (m == 1L) simultaneous <- bonferroni <- point
    expansion[k] <- max(center-simultaneous[, 1L], simultaneous[, 2L]-center) -
      max(widths$simultaneous[, k])
    endpoints <- list(pointwise = point, simultaneous = simultaneous,
                      bonferroni = bonferroni)
    for (band in names(answer$intervals)) {
      answer$intervals[[band]][, , k] <- endpoints[[band]]
      answer$critical[[band]][, k] <- pmax(center-endpoints[[band]][, 1L],
                                           endpoints[[band]][, 2L]-center)
    }
  }
  answer$settings$marginal.calibration <- "signed-contribution score envelope"
  answer$settings$marginal.range <- .np_udist_marginal_range(q)
  answer$settings$marginal.target <- "fixed-bandwidth contribution expectation"
  answer$settings$nesting.expansion <- expansion
  answer
}
