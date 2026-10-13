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
              gaussian.draws = 9999L, neval = 3L, xtrim = .15,
              output = "data", view = "fixed", layout = FALSE, data.overlay = FALSE)
  expect_length(out, 1L)
  expect_identical(names(out[[1L]]$xeval), names(x))
  expect_true(is.ordered(out[[1L]]$xeval$second))
  expect_identical(as.character(out[[1L]]$xeval$second), rep(c("10", "20"), each = 3L))
  reference <- npqreg(bw, txdat = x, tydat = y, exdat = out[[1L]]$xeval, tau = .4, se = TRUE)
  expect_equal(as.numeric(out[[1L]]$quantile), as.numeric(fitted(reference)), tolerance = 2e-11)
  expect_identical(out[[1L]]$band.calibration$grid.size, 6L)
  expect_true(all(is.finite(out[[1L]]$band.errors$simultaneous)))
  expect_equal(as.numeric(out[[1L]]$band.errors$pointwise),
               as.numeric(cbind(qnorm(.975) * se(reference), qnorm(.975) * se(reference))),
               tolerance = 2e-11)
  expect_identical(.Random.seed, seed)
  expect_length(plot(fit, errors = "none", neval = 3L, output = "data",
                     view = "fixed", data.overlay = FALSE), 2L)
  expect_length(plot(fit, tau = c(.7, .3), perspective = TRUE, errors = "none",
                     neval = 3L, output = "data", view = "fixed", data.overlay = FALSE), 2L)

  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({grDevices::dev.off(); unlink(file)}, add = TRUE)
  plot.new()
  old <- par(c("cex", "mfrow"))
  tail <- suppressWarnings(plot(fit, tau = 1e-10, perspective = TRUE,
    errors = "asymptotic", band = "all", gaussian.draws = 9999L, neval = 3L,
    output = "plot-data", view = "fixed", layout = FALSE, data.overlay = FALSE))
  expect_true(all(is.finite(tail[[1L]]$quantile)))
  expect_true(all(is.na(tail[[1L]]$quanterr)))
  expect_identical(par(c("cex", "mfrow")), old)
  expect_false(par("new"))
})
