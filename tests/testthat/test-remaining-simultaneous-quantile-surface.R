test_that("explicit scalar quantile surfaces preserve their complete grid", {
  set.seed(413L)
  n <- 36L
  x <- data.frame(first = seq(-1, 1, length.out = n),
                  second = ordered(rep(c("10", "20"), n/2L), levels = c("10", "20")))
  y <- .4 * x$first + rnorm(n, sd = .55)
  bw <- npcdistbw(xdat = x, ydat = y, bws = c(.42, .65, .15),
                  bandwidth.compute = FALSE)
  fit <- npqreg(bw, txdat = x, tydat = y, tau = .4, se = TRUE)
  seed <- .Random.seed
  out <- plot(fit, perspective = TRUE, errors = "asymptotic", band = "all",
              gaussian_draws = 9999L, neval = 3L, xtrim = .15,
              output = "data", view = "fixed", layout = FALSE, data_overlay = FALSE)
  expect_length(out$panels, 1L)
  expect_identical(names(plot_coordinates(out, role = "X")), names(x))
  expect_true(is.ordered(plot_coordinates(out, role = "X")$second))
  expect_identical(as.character(plot_coordinates(out, role = "X")$second), rep(c("10", "20"), each = 3L))
  reference <- npqreg(bw, txdat = x, tydat = y, exdat = plot_coordinates(out, role = "X"), tau = .4, se = TRUE)
  expect_equal(as.numeric(plot_target(out, band = "pointwise")$estimate), as.numeric(fitted(reference)), tolerance = 2e-11)
  expect_identical(out$inference[[1]]$grid_size, 6L)
  expect_true(all(is.finite(plot_band(out, "simultaneous"))))
  expect_equal(as.numeric(plot_band(out, "pointwise")),
               as.numeric(cbind(fitted(reference) - qnorm(.975) * se(reference), fitted(reference) + qnorm(.975) * se(reference))),
               tolerance = 2e-11)
  expect_identical(.Random.seed, seed)
  expect_length(plot(fit, errors = "none", neval = 3L, output = "data",
                     view = "fixed", data_overlay = FALSE)$panels, 2L)
  expect_length(plot(fit, tau = c(.7, .3), perspective = TRUE, errors = "none",
                     neval = 3L, output = "data", view = "fixed", data_overlay = FALSE)$panels, 2L)

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
  plot.new()
  old <- par(c("cex", "mfrow"))
  tail <- suppressWarnings(plot(fit, tau = 1e-10, perspective = TRUE,
    errors = "asymptotic", band = "all", gaussian_draws = 9999L, neval = 3L,
    output = "plot-data", view = "fixed", layout = FALSE, data_overlay = FALSE))
  expect_true(all(is.finite(plot_target(tail, band = "pointwise")$estimate)))
  expect_true(all(is.na(plot_band(tail, "pointwise"))))
  expect_identical(par(c("cex", "mfrow")), old)
  expect_false(par("new"))
})
