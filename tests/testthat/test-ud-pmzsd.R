test_that("normal corrected bands use alpha and retain provider evidence", {
  normal <- getFromNamespace(".np_plot_bias_normal", "npRmpi")
  x <- list(center = c(-.2, .5, 1.1), se = c(.1, 0, .2),
    bands = list(pointwise = "retained"), settings = list(law = "inid"))
  for (alpha in c(.05, .1)) {
    z <- normal(x, alpha)
    expect_identical(z$bands$pointwise, x$bands$pointwise)
    expect_lt(max(abs(z$bands$pmzsd[, 1, 1] - (x$center - qnorm(1-alpha/2)*x$se))), 1e-14)
    expect_lt(max(abs(z$bands$pmzsd[, 2, 1] - (x$center + qnorm(1-alpha/2)*x$se))), 1e-14)
    expect_identical(z$settings$display.band, "pmzsd")
    expect_false(z$normal.band$calibrated.envelope)
  }
  expect_error(normal(x, .5), "alpha")
  x$se[1] <- -1
  expect_error(normal(x, .05), "nonnegative corrected scales")
})

test_that("omitted corrected density and CDF bands retain pmzsd semantics", {
  skip_if_not(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2, 2, length.out = 25))
  for (family in c("udens", "udist")) {
    bw <- get(paste0("np", family, "bw"))(dat=x, bws=.6, bandwidth.compute=FALSE)
    for (method in c("asymptotic", "bootstrap")) {
      args <- list(x=bw, xdat=x, errors=method, center="bias-corrected",
        neval=3L, output="data", alpha=.1)
      if (method == "bootstrap") args$B <- 49L
      set.seed(628028)
      entry <- .Random.seed
      cap <- capture_plot_provider(do.call(plot, args))
      a <- cap$provider$d1
      flat <- as.data.frame(cap$result)
      expect_identical(entry, .Random.seed)
      args$band <- "pmzsd"
      b <- capture_plot_provider(do.call(plot, args))$provider$d1
      expect_identical(a$inference, b$inference)
      expect_equal(flat$lower,
        a$inference$bands$pmzsd[,1L,1L], tolerance=1e-14)
      expect_equal(flat$upper,
        a$inference$bands$pmzsd[,2L,1L], tolerance=1e-14)
    }
  }
})
