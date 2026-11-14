test_that("constant automatic surface limits expand without changing ordinary ranges", {
  expand <- getFromNamespace(".np_plot_expand_surface_range", "npRmpi")
  for (value in c(0, -1, 1, -.Machine$double.xmax, .Machine$double.xmax)) {
    limits <- expand(c(value, value))
    expect_true(all(is.finite(limits)))
    expect_lt(limits[1L], limits[2L])
    expect_true(limits[1L] <= value && limits[2L] >= value)
  }
  expect_identical(expand(c(1, 2)), c(1, 2))
  expect_identical(expand(c(NA_real_, NA_real_)), c(NA_real_, NA_real_))
  expect_identical(expand(c(0, 0)), c(-.04, .04))
})

test_that("constant conditional surfaces and replay preserve numerical data", {
  skip_on_cran()
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(31)
  x <- c(runif(30, 0, .1), runif(30, .9, 1))
  d <- data.frame(x = x, y = x + rnorm(60, sd = .2))
  for (kernel in c("epanechnikov", "gaussian")) {
    fit <- npcdens(y ~ x, data = d, bws = c(.0005, .04),
      bandwidth.compute = FALSE, cxkertype = kernel, cykertype = "uniform")
    before <- fit
    for (errors in c("none", "asymptotic")) {
      args <- list(x = fit, view = "fixed", neval = 6L, ytrim = .03,
        errors = errors, output = "data")
      if (errors != "none") args$band <- "pointwise"
      seed <- .Random.seed
      data <- suppressWarnings(do.call(plot, args))
      estimate <- data$panels[[1L]]$values$estimate
      expect_length(unique(estimate[is.finite(estimate)]), 1L)
      args$output <- "plot-data"
      drawn <- suppressWarnings(do.call(plot, args))
      expect_identical(drawn$panels, data$panels)
      expect_no_error(suppressWarnings(plot(data)))
      args$zlim <- c(0, 1)
      explicit <- suppressWarnings(do.call(plot, args))
      expect_identical(explicit$panels[[1L]]$values, data$panels[[1L]]$values)
      expect_identical(explicit$panels[[1L]]$intervals, data$panels[[1L]]$intervals)
      expect_identical(explicit$panels[[1L]]$display$style$zlim, c(0, 1))
      args$zlim <- c(0, 0)
      expect_error(suppressWarnings(do.call(plot, args)), "invalid 'z' limits")
      expect_identical(.Random.seed, seed)
    }
    expect_identical(fit, before)
  }
})
