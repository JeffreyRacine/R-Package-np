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
  expect_identical(formals(.np_plot_gaussian_calibrate)$draws, 99999L)
  expect_identical(.np_plot_gaussian_draws(99999L, .05, 2500L), 99999L)
  expect_warning(expect_identical(.np_plot_gaussian_draws(99L, .05, 50L), 99L),
                 "1999.*retaining")
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
