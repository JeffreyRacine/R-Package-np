# These route fixtures deliberately use small explicit draw counts. Precision
# guidance itself is tested below; other warnings remain visible.
local({
  low_budget <- function(w) {
    if (grepl("below the recommended Monte Carlo precision count", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  }
  test_that <- function(desc, code) {
    testthat::test_that(desc, withCallingHandlers(force(code), warning = low_budget))
  }
test_that("Gaussian plot calibration supports singularity and canonical draw blocks", {
  C <- matrix(c(1, .6, .6, 1), 2L)
  a <- .np_with_seed(412L, .np_plot_gaussian_calibrate(C, draws = 199L, block = 1L))
  b <- .np_with_seed(412L, .np_plot_gaussian_calibrate(C, draws = 199L, block = 37L))
  expect_equal(a, b, tolerance = 1e-14)
  expect_identical(.np_plot_gaussian_calibrate(matrix(0, 2, 2))$critical, 0)
  singular <- .np_with_seed(412L,
    .np_plot_gaussian_calibrate(matrix(1, 2, 2), draws = 199L))
  expect_equal(singular$rank, 1L)
  expect_error(.np_plot_gaussian_calibrate(matrix(c(1, 2, 2, 1), 2)), "indefinite")
  expect_error(.np_plot_gaussian_calibrate(C, draws = 1.5), "integer")
  expect_error(.np_plot_gaussian_calibrate(C, max.bytes = 1), "budget")
  expect_identical(formals(.np_plot_gaussian_calibrate)$draws, 9999L)
  expect_identical(.np_plot_gaussian_draws(99999L, .05, 2500L), 99999L)
  expect_warning(expect_identical(.np_plot_gaussian_draws(99L, .05, 50L), 99L),
                 "7299.*retaining")
  expect_identical(.np_plot_gaussian_draws(1000001L, .05, 50L), 1000001L)
  A <- matrix(c(1, 2, 0, 1, 0, 0), 2L)
  a <- .np_with_seed(412L, .np_plot_gaussian_factor_calibrate(A, draws = 199L, block = 1L))
  b <- .np_with_seed(412L, .np_plot_gaussian_factor_calibrate(A, draws = 199L, block = 37L))
  expect_equal(a, b, tolerance = 1e-14)
})

test_that("regression joint bands preserve marginal intervals and mean payloads", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- x$x + sin(seq_len(24L)) / 4
  bw <- npregbw(xdat = x, ydat = y, bws = .4, regtype = "ll",
                bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise",
                neval = 7L, output = "data")[[1L]]
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 399L,
              neval = 7L, output = "data")[[1L]]
  expect_identical(all$mean, point$mean)
  expect_identical(all$merr, point$merr)
  expect_equal(all$band.errors$pointwise, abs(point$merr), ignore_attr = TRUE)
  expect_true(all(is.finite(all$band.errors$simultaneous)))
  expect_identical(all$band.calibration$scope, "finite grid within panel")
  grad <- plot(bw, errors = "asymptotic", band = "all", gradients = TRUE,
               neval = 7L, gaussian.draws = 399L, output = "data")[[1L]]
  expect_true(all(is.finite(grad$band.errors$simultaneous)))
  expect_equal(dim(grad$band.errors$simultaneous), c(7L, 2L))
  expect_error(plot(bw, gaussian.draws = 399L, output = "data"), "requires asymptotic")
  expect_error(plot(bw, errors = "asymptotic", band = "all", B = 199L),
               "bootstrap controls")
})


test_that("single-index joint bands condition on fitted coefficients", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L), z = sin(seq_len(24L)))
  y <- x$x + .2 * x$z + cos(seq_len(24L)) / 4
  bw <- npindexbw(xdat = x, ydat = y, bws = c(1, .2, .4),
                  bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise",
                neval = 7L, output = "data")[[1L]]
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 399L,
              neval = 7L, output = "data")[[1L]]
  expect_identical(all$mean, point$mean)
  expect_identical(all$merr, point$merr)
  expect_true(all(is.finite(all$band.errors$simultaneous)))
  expect_match(all$band.calibration$conditioning, "fitted index coefficients")
  grad <- plot(bw, errors = "asymptotic", band = "all", gradients = TRUE,
               neval = 7L, gaussian.draws = 399L, output = "data")[[1L]]
  expect_equal(dim(grad$band.errors$simultaneous), c(7L, 2L, 2L))
  expect_true(all(is.finite(grad$band.errors$simultaneous)))
  expect_equal(grad$band.errors$simultaneous[, , 2L],
               .2 * grad$band.errors$simultaneous[, , 1L])
})


test_that("smooth-coefficient mean bands retain the accepted marginal errors", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 30L), w = sin(seq_len(30L)))
  z <- data.frame(z = cos(seq_len(30L)))
  y <- x$x + (1 + z$z) * x$w + sin(2 * seq_len(30L)) / 4
  bw <- npscoefbw(xdat = x, ydat = y, zdat = z, bws = .6,
                  bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise", perspective = FALSE,
                neval = 7L, output = "data")
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 399L,
              perspective = FALSE, neval = 7L, output = "data")
  for (i in seq_along(all)) {
    expect_identical(all[[i]]$mean, point[[i]]$mean)
    expect_identical(all[[i]]$merr, point[[i]]$merr)
    expect_true(all(is.finite(all[[i]]$band.errors$simultaneous)))
  }
  grad <- plot(bw, errors = "asymptotic", band = "all", gradients = TRUE,
               perspective = FALSE, neval = 7L, gaussian.draws = 399L, output = "data")
  expect_equal(length(grad), 3L)
  for (panel in grad) {
    expect_equal(dim(panel$band.errors$simultaneous), c(7L, 2L))
    expect_true(all(is.finite(panel$band.errors$simultaneous)))
  }
})

for (family in c("npreg", "npindex", "npscoef")) {
  test_that(paste(family, "joint bands normalize legacy values before admission"), {
    skip_on_cran()
    skip_if_not(spawn_mpi_slaves())
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
    old <- options(np.messages = FALSE)
    on.exit(options(old), add = TRUE)
    x <- data.frame(x = seq(-1, 1, length.out = 36L), w = sin(seq_len(36L)))
    z <- data.frame(z = cos(seq_len(36L)))
    y <- x$x + (1 + z$z) * x$w + sin(2 * seq_len(36L)) / 4
    bw <- switch(family,
      npreg = npregbw(xdat = x["x"], ydat = y, bws = .5,
                      bandwidth.compute = FALSE),
      npindex = npindexbw(xdat = x, ydat = y, bws = c(1, .2, .5),
                          bandwidth.compute = FALSE),
      npscoef = npscoefbw(xdat = x, ydat = y, zdat = z, bws = .6,
                          bandwidth.compute = FALSE)
    )
    args <- list(x = bw, neval = 7L, output = "data")
    if (family == "npscoef") args$perspective <- FALSE
    for (band in c("all", "simultaneous")) {
      full <- do.call(plot, c(args, list(errors = "asymptotic", band = band,
                                        gaussian.draws = 399L)))
      legacy <- do.call(plot, c(args, list(plot.errors.method = "asy",
                                          plot.errors.type = substr(band, 1L, 1L),
                                          gaussian.draws = 399L)))
      for (i in seq_along(full)) {
        expect_identical(legacy[[i]]$band.calibration, full[[i]]$band.calibration)
        expect_identical(legacy[[i]]$merr, full[[i]]$merr)
        expect_identical(legacy[[i]]$band.errors, full[[i]]$band.errors)
      }
    }
    full <- do.call(plot, c(args, list(errors = "asymptotic", band = "all")))
    legacy <- do.call(plot, c(args, list(plot.errors.method = "asy",
                                        plot.errors.type = "a")))
    for (i in seq_along(full))
      expect_identical(legacy[[i]]$band.calibration, full[[i]]$band.calibration)
    expect_error(do.call(plot, c(args, list(errors = "asy", band = "all"))),
                 "errors")
    expect_error(do.call(plot, c(args, list(errors = "asymptotic", band = "a"))),
                 "band")
    expect_error(do.call(plot, c(args, list(plot.errors.method = "asy",
                                           plot.errors.type = "a",
                                           gaussian.draws = 0))), "integer")
  })
}

test_that("Gaussian mean surfaces calibrate every evaluation point", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 30L), w = sin(seq_len(30L)))
  y <- x$x + x$w + cos(seq_len(30L)) / 4
  models <- list(
    npregbw(xdat = x, ydat = y, bws = c(.8, .8), bandwidth.compute = FALSE),
    npscoefbw(xdat = x["x"], zdat = x["w"], ydat = y, bws = .8,
               bandwidth.compute = FALSE)
  )
  for (bw in models) {
    out <- plot(bw, errors = "asymptotic", band = "all", neval = 5L,
                gaussian.draws = 999L, output = "data")[[1L]]
    expect_identical(out$band.calibration$grid.size, 25L)
    expect_equal(dim(out$band.errors$simultaneous), c(25L, 2L))
    expect_true(all(is.finite(out$band.errors$simultaneous)))
  }
})

test_that("legacy error request resolution preserves scalar-choice defaults", {
  expect_identical(.np_plot_normalize_error_request("asy", "s"),
                   list(method = "asymptotic", type = "simultaneous"))
  for (empty in list(NULL, character(), NA_character_))
    expect_identical(.np_plot_normalize_error_request(empty, empty),
                     list(method = "none", type = "simultaneous"))
  expect_identical(.np_plot_normalize_error_request(c("none", "asymptotic"),
                                                   c("pmzsd", "all")),
                   list(method = "none", type = "pmzsd"))
  expect_error(.np_plot_normalize_error_request("invalid", "all"), "arg")
  expect_error(.np_plot_normalize_error_request("asymptotic", "p"), "arg")
})

test_that("NN joint normalization retains training radii", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- x$x^2 + sin(seq_len(24L)) / 4
  for (kind in c("generalized_nn", "adaptive_nn")) {
    bw <- npregbw(xdat = x, ydat = y, bws = 14L, bwtype = kind,
                  regtype = "ll", bandwidth.compute = FALSE)
    S <- npreghat(bw, txdat = x)
    E <- diag(nrow(x)) - S
    expected <- drop(E %*% y) / sqrt(rowSums(E^2))
    a <- .np_plot_regression_joint_covariance(bw, x, y, x[1:5, , drop = FALSE])
    expect_equal(a$residual, expected, tolerance = 1e-10)
    if (kind == "generalized_nn") {
      external <- diag(nrow(x)) - npreghat(bw, txdat = x, exdat = x)
      wrong <- drop(external %*% y) / sqrt(rowSums(external^2))
      expect_gt(max(abs(expected - wrong)), 1e-6)
    }
    out <- plot(bw, errors = "asymptotic", band = "all", gradients = TRUE,
                gaussian.draws = 199L, neval = 5L, output = "data")[[1L]]
    expect_true(all(is.finite(out$band.errors$simultaneous)))
  }
})

# Unique contract: practical reconstruction tolerance is tied to uncertainty,
# with response-unit invariance and material/nonfinite failures kept explicit.
test_that("partially linear reconstruction uses the plotted uncertainty scale", {
  check <- .np_plot_plreg_check_reconstruction
  for (scale in c(1e-12, 1, 1e12)) {
    expect_identical(check(scale * c(0, 2), scale * c(1e-5, 2.005),
                           scale * c(.1, 1)), c(FALSE, FALSE))
    expect_identical(check(scale * c(0, 2), scale * c(1e-5, 2.02),
                           scale * c(.1, 1)), c(FALSE, TRUE))
  }
  expect_identical(check(c(0, 2), c(0, 2), c(0, 0)), c(FALSE, FALSE))
  expect_false(check(2, 2 + .Machine$double.eps * 4, 0))
  expect_true(check(0, 1e-8, 0))
  expect_error(check(0, 0, -1), "invalid fitted values or standard errors")
  expect_error(check(0, Inf, 1), "invalid fitted values or standard errors")
  expect_error(check(c(0, 1), 0, 1), "invalid fitted values or standard errors")
})

test_that("partially linear categorical influence uses the fitted mean owner", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)
  z <- data.frame(z = factor(rep(letters[1:4], 9L)))
  x <- data.frame(x = sin(seq_len(36L)) + as.integer(z$z))
  y <- .5 * x$x + cos(seq_len(36L))
  bw <- npplregbw(xdat = x, ydat = y, zdat = z,
                  bws = matrix(c(.35, .7425), 2L),
                  regtype = "lc", bandwidth.compute = FALSE)
  ex <- x[1:4, , drop = FALSE]; ez <- z[1:4, , drop = FALSE]
  influence <- .np_plot_plreg_influence(bw, x, y, z, ex, ez, factor = TRUE)
  fit <- npplreg(bw, txdat = x, tydat = y, tzdat = z, exdat = ex, ezdat = ez)
  expect_equal(influence$mean, fit$mean, tolerance = 1e-12)
})

test_that("partially linear plots include shared-response covariance", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  z <- data.frame(z = seq(-1, 1, length.out = 24L))
  x <- data.frame(x = z$z + sin(seq_len(24L)))
  y <- .5 * x$x + z$z^2 + cos(seq_len(24L)) / 4
  bw <- npplregbw(xdat = x, ydat = y, zdat = z, bws = matrix(c(.6, .6), 2L),
                  regtype = "ll", bandwidth.compute = FALSE)
  ex <- x[1:5, , drop = FALSE]; ez <- z[1:5, , drop = FALSE]
  a <- .np_plot_plreg_influence(bw, x, y, z, ex, ez, factor = TRUE)
  shifted <- .np_plot_plreg_influence(bw, x, y + 5 * x$x, z, ex, ez, factor = TRUE)
  expect_equal(a$se, shifted$se, tolerance = 1e-10)
  expect_equal(shifted$mean - a$mean, 5 * ex$x, tolerance = 1e-10)
  H <- npplreghat(bw, txdat = x, tzdat = z, exdat = ex, ezdat = ez, output = "matrix")
  S <- npplreghat(bw, txdat = x, tzdat = z, output = "matrix")
  E <- diag(nrow(x)) - S
  factor <- sweep(H, 2L, drop(E %*% y) / sqrt(rowSums(E^2)), "*")
  expect_equal(a$influence, unclass(factor), ignore_attr = TRUE, tolerance = 1e-10)
  out <- plot(bw, errors = "asymptotic", band = "all", perspective = FALSE,
              gaussian.draws = 199L, neval = 5L, output = "data")
  point <- plot(bw, errors = "asymptotic", band = "pointwise", perspective = FALSE,
                neval = 5L, output = "data")
  for (j in seq_along(out)) {
    expect_identical(out[[j]]$mean, point[[j]]$mean)
    expect_identical(out[[j]]$merr, point[[j]]$merr)
    expect_true(all(is.finite(out[[j]]$band.errors$simultaneous)))
    expect_match(out[[j]]$band.calibration$conditioning, "includes coefficient estimation")
  }
})

# Unique contract: categorical endpoint covariance and degenerate contrast
# rows. Broad kernel/bandwidth/family matrices stay in the external harness.
test_that("categorical joint bands retain endpoint covariance and level semantics", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(u = factor(rep(c("c", "a", "b"), 12L),
                            levels = c("c", "a", "b")),
                  o = ordered(rep(c("low", "mid", "high"), each = 12L),
                              levels = c("low", "mid", "high")))
  y <- as.integer(x$u) / 3 + as.integer(x$o) / 4 + sin(seq_len(36L)) / 5
  bw <- npregbw(xdat = x, ydat = y, bws = c(.2, .3),
                regtype = "lc", bandwidth.compute = FALSE)
  ex <- x[c(1L, 2L, 14L, 28L), , drop = FALSE]
  for (j in seq_len(2L)) {
    frames <- npCategoricalFirstDifferenceFrames(ex, j, "test")
    D <- npreghat(bw, txdat = x, exdat = frames$upper, output = "matrix") -
      npreghat(bw, txdat = x, exdat = frames$lower, output = "matrix")
    joint <- .np_plot_regression_joint_covariance(bw, x, y, ex,
                                                  categorical.index = j)
    expect_equal(joint$estimate, drop(D %*% y), ignore_attr = TRUE)
    expect_equal(joint$covariance, tcrossprod(sweep(D, 2L, joint$residual, "*")),
                 ignore_attr = TRUE)
  }
  levels <- plot(bw, errors = "asymptotic", band = "all", output = "data",
                 perspective = FALSE, gaussian.draws = 199L)
  grad <- plot(bw, errors = "asymptotic", band = "all", gradients = TRUE,
               output = "data", perspective = FALSE, gaussian.draws = 199L)
  point <- plot(bw, errors = "asymptotic", band = "pointwise", gradients = TRUE,
                output = "data", perspective = FALSE)
  for (j in seq_len(2L)) {
    expect_true(all(is.finite(levels[[j]]$band.errors$simultaneous)))
    expect_true(all(is.finite(grad[[j]]$band.errors$simultaneous)))
    expect_identical(grad[[j]]$grad, point[[j]]$grad)
    expect_identical(grad[[j]]$gerr, point[[j]]$gerr)
  }
  expect_equal(grad[[1L]]$band.errors$simultaneous[1L, ], c(0, 0))
  expect_equal(grad[[2L]]$band.errors$simultaneous[1L, ],
               grad[[2L]]$band.errors$simultaneous[2L, ])
})

})

test_that("Gaussian calibration has the independent maximum probability", {
  B <- 9999L
  m <- 4L
  result <- .np_with_seed(71L, .np_plot_gaussian_calibrate(diag(m), draws = B))
  # For independent standard normals P(max |Z_j| <= c) = (2 Phi(c)-1)^m.
  # The empirical type-1 quantile has an exact beta law on that CDF scale.
  probability <- (2 * pnorm(result$critical) - 1)^m
  k <- ceiling(.95 * B)
  limits <- qbeta(c(.0005, .9995), k, B + 1L - k)
  expect_gt(probability, limits[1L])
  expect_lt(probability, limits[2L])
  expect_lt((2 * pnorm(qnorm(.975)) - 1)^m, limits[1L])
  expect_lt((2 * pnorm(result$critical / 2) - 1)^m, limits[1L])
})

test_that("joint unavailability preserves only the declared compatibility boundary", {
  expect_warning(expect_null(.np_plot_joint_try(
    .np_plot_joint_unavailable("identified limitation"))), "NA simultaneous")
  expect_error(.np_plot_joint_try(stop("unexpected failure")), "unexpected failure")
  expect_error(.np_plot_joint_try(.np_plot_gaussian_calibrate(
    matrix(c(1, 2, 2, 1), 2))), "indefinite")
  expect_warning(expect_identical(.np_plot_gaussian_draws(7298L, .05, 2500L), 7298L), "7299")
  expect_warning(.np_plot_gaussian_draws(9999L, .01, 50L), "38031")
  expect_warning(.np_plot_gaussian_draws(7299L, .05, 2500L), NA)
})

# Protect NA row retention and the finite display range after joint rejection.
test_that("unavailable joint limits retain rows and render fitted centers", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  errors <- matrix(NA_real_, 3L, 2L)
  kept <- .np_plot_joint_interval_rows(errors, c(1, NA, 3), TRUE, "simultaneous")
  expect_identical(dim(kept), c(2L, 2L))
  expect_true(all(is.na(kept)))
  expect_identical(.np_plot_joint_panel_range(c(1, NA, 3), errors), c(1, 3))
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- x$x + sin(seq_len(24L)) / 4
  bw <- npregbw(xdat = x, ydat = y, bws = .5, regtype = "ll",
                bandwidth.compute = FALSE)
  for (gradient in c(FALSE, TRUE)) {
    expect_warning(out <- plot(bw, errors = "asymptotic", band = "simultaneous",
      gradients = gradient, neval = 513L, output = "data")[[1L]], "memory/grid budget")
    limits <- if (gradient) out$gerr else out$merr
    expect_identical(dim(limits), c(513L, 2L))
    expect_true(all(is.na(limits)))
  }
})

local({
test_that("Gaussian Monte Carlo critical values retain analytical nesting", {
  low <- qnorm(.025, lower.tail = FALSE)
  high <- qnorm(.025 / 2, lower.tail = FALSE)
  expect_identical(.np_plot_gaussian_critical(0, .05, 2), low)
  expect_identical(.np_plot_gaussian_critical(100, .05, 2), high)
  expect_identical(.np_plot_gaussian_critical(2, .05, 2), 2)
  withr::local_seed(4L)
  out <- .np_plot_gaussian_factor_calibrate(matrix(1, 2, 1))
  state <- .Random.seed
  set.seed(4L)
  invisible(rnorm(9999L))
  expect_identical(.Random.seed, state)
  expect_lt(out$critical.raw, low)
  expect_identical(out$critical, low)
  expect_identical(.np_plot_gaussian_factor_calibrate(matrix(0, 2, 1))$critical, 0)
})

test_that("PL influence retains the fitted scalar ridge decision", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.tree = FALSE, np.macMseries.accelerate = FALSE,
                           np.messages = FALSE))
  data(wage1, package = "npRmpi", envir = environment())
  z <- wage1["tenure"]; y <- wage1$lwage
  bw <- npregbw(xdat = z, ydat = y, bws = 0.42347409342680487,
    regtype = "lp", degree = 2L, bernstein.basis = FALSE,
    bandwidth.compute = FALSE)
  ez <- data.frame(tenure = 42.204081632653065)
  before <- .np_regression_direct(bw, z, y, ez)$mean
  operator <- .np_plot_plreg_operator(bw, z, y, ez)
  H <- operator(1L)
  withr::local_seed(1003L)
  perturbations <- matrix(sample(c(-1, 1), length(y) * 19L, TRUE), length(y))
  direct <- vapply(seq_len(ncol(perturbations)), function(j)
    .np_plot_with_local_compiled_eval(.np_regression_direct(bw, z, y + perturbations[, j], ez))$mean - before,
    numeric(1L))
  reference.sd <- sd(direct)
  expect_lt(abs(drop(H %*% y) - before), .01 * reference.sd)
  expect_lt(sd(direct - drop(H %*% perturbations)), .01 * reference.sd)
  after <- .np_plot_with_local_compiled_eval(.np_regression_direct(bw, z, y, ez, fit.projection = TRUE))
  expect_identical(after$mean, before)
})
})


test_that("PL plot projections retain the canonical degree-zero owner", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  withr::local_seed(661L)
  n <- 40L
  z <- data.frame(z = runif(n)); y <- sin(z$z) + rnorm(n)
  for (type in c("fixed", "adaptive_nn", "generalized_nn")) {
    bw <- npregbw(xdat = z, ydat = y, bws = if (type == "fixed") .3 else 25L,
      bwtype = type, regtype = "lp", degree = 0L, bandwidth.compute = FALSE)
    operator <- .np_plot_plreg_operator(bw, z, y)
    H <- do.call(rbind, lapply(split(seq_len(n), ceiling(seq_len(n) / 7)), operator))
    fitted <- .np_regression_direct(bw, z, y)$mean
    expect_equal(drop(H %*% y), fitted, tolerance = 1e-10)
  }
})
