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
})

test_that("regression joint bands preserve marginal intervals and mean payloads", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- x$x + sin(seq_len(24L)) / 4
  bw <- npregbw(xdat = x, ydat = y, bws = .4, regtype = "ll",
                bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise",
                neval = 7L, output = "data")[[1L]]
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 199L,
              neval = 7L, output = "data")[[1L]]
  expect_identical(all$mean, point$mean)
  expect_identical(all$merr, point$merr)
  expect_equal(all$band.errors$pointwise, abs(point$merr), ignore_attr = TRUE)
  expect_true(all(is.finite(all$band.errors$simultaneous)))
  expect_identical(all$band.calibration$scope, "finite grid within panel")
  expect_error(plot(bw, gaussian.draws = 199L, output = "data"), "requires asymptotic")
  expect_error(plot(bw, errors = "asymptotic", band = "all", B = 199L),
               "bootstrap controls")
})


test_that("single-index joint bands condition on fitted coefficients", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L), z = sin(seq_len(24L)))
  y <- x$x + .2 * x$z + cos(seq_len(24L)) / 4
  bw <- npindexbw(xdat = x, ydat = y, bws = c(1, .2, .4),
                  bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise",
                neval = 7L, output = "data")[[1L]]
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 199L,
              neval = 7L, output = "data")[[1L]]
  expect_identical(all$mean, point$mean)
  expect_identical(all$merr, point$merr)
  expect_true(all(is.finite(all$band.errors$simultaneous)))
  expect_match(all$band.calibration$conditioning, "fitted index coefficients")
})


test_that("smooth-coefficient mean bands retain the accepted marginal errors", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 30L), w = sin(seq_len(30L)))
  z <- data.frame(z = cos(seq_len(30L)))
  y <- x$x + (1 + z$z) * x$w + sin(2 * seq_len(30L)) / 4
  bw <- npscoefbw(xdat = x, ydat = y, zdat = z, bws = .6,
                  bandwidth.compute = FALSE)
  point <- plot(bw, errors = "asymptotic", band = "pointwise", perspective = FALSE,
                neval = 7L, output = "data")
  all <- plot(bw, errors = "asymptotic", band = "all", gaussian.draws = 199L,
              perspective = FALSE, neval = 7L, output = "data")
  for (i in seq_along(all)) {
    expect_identical(all[[i]]$mean, point[[i]]$mean)
    expect_identical(all[[i]]$merr, point[[i]]$merr)
    expect_true(all(is.finite(all[[i]]$band.errors$simultaneous)))
  }
})
