adaptive_distribution_epanechnikov2_cdf <- function(u) {
  edge <- sqrt(5)
  ifelse(u <= -edge, 0,
         ifelse(u >= edge, 1,
                0.5 + 3 / (4 * edge) * (u - u^3 / 15)))
}

adaptive_distribution_fold_radius <- function(value, held_out, donor, k) {
  keep <- setdiff(seq_along(value), c(held_out, donor))
  distances <- sort(abs(value[keep] - value[donor]), method = "radix")
  lookup_k <- min(k, length(distances))
  distances[[lookup_k]] * k / lookup_k
}

adaptive_distribution_literal_objective <- function(value, k, kernel,
                                                       ordered = NULL,
                                                       lambda = NULL) {
  value <- as.matrix(value)
  n <- nrow(value)
  objective <- 0
  for (held_out in seq_len(n)) {
    donors <- setdiff(seq_len(n), held_out)
    radius <- vapply(seq_len(ncol(value)), function(coordinate) {
      vapply(donors, function(donor) {
        adaptive_distribution_fold_radius(
          value[, coordinate], held_out, donor, k[[coordinate]])
      }, numeric(1L))
    }, numeric(length(donors)))
    if (ncol(value) == 1L)
      radius <- matrix(radius, ncol = 1L)
    for (evaluation in setdiff(seq_len(n), held_out)) {
      weight <- vapply(seq_along(donors), function(index) {
        donor <- donors[[index]]
        u <- (value[evaluation, ] - value[donor, ]) / radius[index, ]
        continuous <- if (kernel == "gaussian") {
          prod(stats::pnorm(u))
        } else {
          prod(adaptive_distribution_epanechnikov2_cdf(u))
        }
        categorical <- if (is.null(ordered)) {
          1
        } else {
          evaluation_level <- as.integer(ordered[[evaluation]])
          donor_level <- as.integer(ordered[[donor]])
          distance <- abs(evaluation_level - donor_level)
          if (evaluation_level == donor_level) {
            1 - 0.5 * lambda
          } else if (evaluation_level < donor_level) {
            0.5 * lambda^distance
          } else {
            1 - lambda^distance
          }
        }
        continuous * categorical
      }, numeric(1L))
      fitted <- mean(weight)
      indicator <- all(value[held_out, ] <= value[evaluation, ]) &&
        (is.null(ordered) ||
         as.integer(ordered[[held_out]]) <= as.integer(ordered[[evaluation]]))
      objective <- objective + (as.double(indicator) - fitted)^2
    }
  }
  objective / (n * (n - 1L))
}

adaptive_distribution_package_objective <- function(dat, k, kernel,
                                                      tree, accelerate,
                                                      lambda = NULL) {
  bandwidth <- if (is.null(lambda)) k else c(k, lambda)
  state <- npudistbw(
    dat = dat, bws = bandwidth, bwmethod = "cv.cdf",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = kernel, okertype = "wangvanryzin",
    bandwidth.compute = FALSE)
  old_options <- options(np.tree = tree,
                         np.macMseries.accelerate = accelerate,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)
  as.numeric(np:::npudistbw.dbandwidth(
    dat = dat, bws = state, bandwidth.compute = TRUE,
    nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
    eval.only = TRUE, do.full.integral = TRUE)$fval)
}

test_that("adaptive empirical CDF-CV uses pair-local donor radii", {
  fixtures <- list(
    p1 = list(
      value = matrix(c(-2.2, -1.1, -0.35, 0.15, 0.8, 1.55, 2.7),
                     ncol = 1L),
      k = 2L),
    p2 = list(
      value = cbind(
        c(-1.9, -1.1, -0.45, 0.05, 0.55, 1.05, 1.8, 2.65),
        c(0.75, -0.35, 1.2, 0.1, -0.85, 0.5, -0.55, 0.25)),
      k = c(2L, 3L)),
    p3 = list(
      value = cbind(
        c(-2, -1.3, -0.65, -0.1, 0.35, 0.9, 1.45, 2.1, 2.8),
        c(0.9, -0.4, 1.1, 0.2, -0.8, 0.55, -0.6, 0.35, -0.05),
        c(-0.7, 0.6, -0.2, 1.25, 0.1, -1, 0.85, -0.45, 0.4)),
      k = c(2L, 3L, 2L)))

  for (fixture in fixtures) {
    dat <- as.data.frame(fixture$value)
    names(dat) <- paste0("x", seq_len(ncol(dat)))
    for (kernel in c("gaussian", "epanechnikov")) {
      expected <- adaptive_distribution_literal_objective(
        fixture$value, fixture$k, kernel)
      observed <- unlist(lapply(c(FALSE, TRUE), function(tree) {
        vapply(c(FALSE, TRUE), function(accelerate) {
          adaptive_distribution_package_objective(
            dat, fixture$k, kernel, tree, accelerate)
        }, numeric(1L))
      }))
      expect_equal(observed, rep(expected, length(observed)),
                   tolerance = 2e-10)
    }
  }
})

test_that("adaptive empirical CDF-CV composes categorical factors", {
  value <- c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1, 5, 6.2)
  group <- ordered(c(1, 2, 1, 3, 2, 3, 1, 2, 3, 1), levels = 1:3)
  dat <- data.frame(value = value, group = group)
  k <- 3L
  lambda <- 0.2
  expected <- adaptive_distribution_literal_objective(
    matrix(value, ncol = 1L), k, "gaussian", group, lambda)

  expect_equal(
    adaptive_distribution_package_objective(
      dat, k, "gaussian", tree = TRUE, accelerate = FALSE,
      lambda = lambda),
    expected, tolerance = 2e-10)
})

test_that("adaptive empirical CDF-CV honors saturation and replay", {
  old_options <- options(np.tree = FALSE,
                         np.extendednn = FALSE,
                         np.macMseries.accelerate = FALSE,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)
  dat <- data.frame(
    x = c(-2.8, -1.5, -0.7, -0.1, 0.4, 1.2, 2.7, 4.1))
  k <- nrow(dat) - 2L
  state <- npudistbw(
    dat = dat, bws = k, bwmethod = "cv.cdf",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = "gaussian", bandwidth.compute = FALSE)
  replay <- unserialize(serialize(state, NULL, version = 3L))
  evaluate <- function(bw) {
    as.numeric(np:::npudistbw.dbandwidth(
      dat = dat, bws = bw, bandwidth.compute = TRUE,
      nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
      eval.only = TRUE, do.full.integral = TRUE,
      invalid.penalty = "dbmax")$fval)
  }
  expected <- adaptive_distribution_literal_objective(
    as.matrix(dat), k, "gaussian")

  expect_equal(evaluate(replay), expected, tolerance = 2e-10)
  expect_identical(evaluate(replay), evaluate(state))
  saturated <- state
  saturated$bw[] <- nrow(dat) - 1L
  options(np.extendednn = TRUE)
  saturated_expected <- adaptive_distribution_literal_objective(
    as.matrix(dat), saturated$bw, "gaussian")
  expect_equal(evaluate(saturated), saturated_expected, tolerance = 2e-10)
  options(np.extendednn = FALSE)
  expect_identical(evaluate(saturated), .Machine$double.xmax)

  zero_radius <- npudistbw(
    dat = data.frame(x = c(0, 0, 0, 0, 1, 2)),
    bws = 2L, bwmethod = "cv.cdf", bwtype = "adaptive_nn",
    bwscaling = FALSE, ckertype = "gaussian",
    bandwidth.compute = FALSE)
  expect_identical(
    as.numeric(np:::npudistbw.dbandwidth(
      dat = data.frame(x = c(0, 0, 0, 0, 1, 2)),
      bws = zero_radius, bandwidth.compute = TRUE,
      nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
      eval.only = TRUE, do.full.integral = TRUE,
      invalid.penalty = "dbmax")$fval),
    .Machine$double.xmax)
})
