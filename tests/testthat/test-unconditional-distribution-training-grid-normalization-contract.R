udist_training_grid_objective <- function(dat, bw, tree, accelerate) {
  old_options <- options(np.tree = tree,
                         np.macMseries.accelerate = accelerate,
                         np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)
  evaluated <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")(
    npRmpi:::npudistbw.dbandwidth(
      dat = dat, bws = bw, bandwidth.compute = TRUE,
      nmulti = 1L, powell.remin = FALSE, bwsolver = "powell",
      eval.only = TRUE, do.full.integral = TRUE)
  )
  as.numeric(evaluated$fval)
}

udist_literal_training_grid <- function(value, cdf_weight) {
  value <- as.matrix(value)
  n <- nrow(value)
  objective <- 0
  for (held_out in seq_len(n)) {
    donors <- setdiff(seq_len(n), held_out)
    for (evaluation in setdiff(seq_len(n), held_out)) {
      fitted <- mean(vapply(donors, function(donor) {
        product <- 1
        for (coordinate in seq_len(ncol(value)))
          product <- product * cdf_weight(
            value[evaluation, coordinate], value[donor, coordinate],
            coordinate, held_out)
        product
      }, numeric(1L)))
      indicator <- all(value[held_out, ] <= value[evaluation, ])
      objective <- objective + (as.double(indicator) - fitted)^2
    }
  }
  objective / (n * (n - 1))
}

test_that("unconditional continuous CDF-CV averages the admitted ordered pairs", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  fixtures <- list(
    p1 = matrix(c(-2.1, -1.2, -0.3, 0.2, 0.9, 1.7, 2.8), ncol = 1L),
    p2 = cbind(seq(-1.8, 2.2, length.out = 11L),
               sin(seq_len(11L) * 0.67) + seq_len(11L) / 19))
  bandwidths <- list(p1 = 0.7, p2 = c(0.65, 0.55))

  for (fixture_name in names(fixtures)) {
    value <- fixtures[[fixture_name]]
    dat <- as.data.frame(value)
    names(dat) <- paste0("x", seq_len(ncol(dat)))
    bandwidth <- bandwidths[[fixture_name]]
    for (kernel_name in c("gaussian", "epanechnikov")) {
      kernel_cdf <- if (kernel_name == "gaussian") {
        function(evaluation, donor, coordinate, held_out) {
          stats::pnorm((evaluation - donor) / bandwidth[coordinate])
        }
      } else {
        function(evaluation, donor, coordinate, held_out) {
          u <- (evaluation - donor) / bandwidth[coordinate]
          edge <- sqrt(5)
          if (u <= -edge) return(0)
          if (u >= edge) return(1)
          0.5 + 3 / (4 * edge) * (u - u^3 / 15)
        }
      }
      expected <- udist_literal_training_grid(value, kernel_cdf)
      bw <- npudistbw(
        dat = dat, bws = bandwidth, bwmethod = "cv.cdf",
        bwtype = "fixed", bwscaling = FALSE,
        ckertype = kernel_name, ckerorder = 2L,
        bandwidth.compute = FALSE)
      for (tree in c(FALSE, TRUE))
        for (accelerate in c(FALSE, TRUE))
          expect_equal(
            udist_training_grid_objective(dat, bw, tree, accelerate),
            expected, tolerance = 1e-9)
    }
  }
})

test_that("extended adaptive CDF-CV uses two-exclusion radii and the empirical pair finalizer", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  old_options <- options(np.extendednn = TRUE)
  on.exit(options(old_options), add = TRUE)

  value <- cbind(
    seq(-1.7, 2.2, length.out = 9L),
    sin(seq_len(9L) * 0.71) + seq_len(9L) / 19
  )
  dat <- as.data.frame(value)
  names(dat) <- c("x1", "x2")
  extended_k <- c(nrow(value) + 2L, nrow(value) + 3L)
  expected <- udist_literal_training_grid(
    value,
    function(evaluation, donor, coordinate, held_out) {
      donor_index <- match(donor, value[, coordinate])
      retained <- setdiff(seq_len(nrow(value)), c(held_out, donor_index))
      bandwidth <-
        max(abs(value[retained, coordinate] - donor)) *
        extended_k[[coordinate]] / (nrow(value) - 2L)
      stats::pnorm(
        (evaluation - donor) / bandwidth
      )
    }
  )
  state <- npudistbw(
    dat = dat, bws = extended_k, bwmethod = "cv.cdf",
    bwtype = "adaptive_nn", bwscaling = FALSE,
    ckertype = "gaussian", bandwidth.compute = FALSE
  )

  expect_equal(
    udist_training_grid_objective(dat, state, FALSE, FALSE),
    expected, tolerance = 1e-10
  )
})

test_that("beta CDF-CV uses the same empirical training-grid normalization", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  value <- c(0.04, 0.12, 0.27, 0.43, 0.61, 0.78, 0.93)
  dat <- data.frame(x = value)
  bandwidth <- 0.19
  beta_cdf <- function(evaluation, donor, coordinate, held_out) {
    concentration <- 1 / bandwidth^2
    stats::pbeta(evaluation,
                 1 + donor * concentration,
                 1 + (1 - donor) * concentration)
  }
  expected <- udist_literal_training_grid(matrix(value, ncol = 1L), beta_cdf)
  bw <- npudistbw(
    dat = dat, bws = bandwidth, bwmethod = "cv.cdf",
    bwtype = "fixed", bwscaling = FALSE,
    ckertype = "beta", ckerorder = 2L,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1,
    bandwidth.compute = FALSE)

  expect_equal(udist_training_grid_objective(dat, bw, FALSE, FALSE),
               expected, tolerance = 1e-9)
})

test_that("ordered-profile and dense CDF-CV share the empirical finalizer", {
  old_npRmpi_local_mode <- getOption("npRmpi.local.regression.mode", FALSE)
  options(npRmpi.local.regression.mode = TRUE)
  on.exit(options(npRmpi.local.regression.mode = old_npRmpi_local_mode),
          add = TRUE)
  value <- ordered(c(1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4),
                   levels = 1:4)
  numeric_value <- as.numeric(value)
  dat <- data.frame(x = value)
  lambda <- 0.31
  ordered_cdf <- function(evaluation, donor, coordinate, held_out) {
    if (evaluation == donor) return(1 - 0.5 * lambda)
    distance <- abs(evaluation - donor)
    if (evaluation < donor) 0.5 * lambda^distance else 1 - lambda^distance
  }
  expected <- udist_literal_training_grid(
    matrix(numeric_value, ncol = 1L), ordered_cdf)
  bw <- npudistbw(
    dat = dat, bws = lambda, bwmethod = "cv.cdf",
    bwtype = "fixed", bwscaling = FALSE,
    okertype = "wangvanryzin", bandwidth.compute = FALSE)

  dense <- udist_training_grid_objective(dat, bw, FALSE, FALSE)
  profile <- udist_training_grid_objective(dat, bw, TRUE, FALSE)
  expect_equal(dense, expected, tolerance = 1e-12)
  expect_equal(profile, expected, tolerance = 1e-12)
  expect_identical(profile, dense)
})
