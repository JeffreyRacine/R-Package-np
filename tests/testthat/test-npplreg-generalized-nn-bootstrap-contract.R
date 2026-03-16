test_that("npplreg generalized-nn cv.aic bootstrap representative session contract holds", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260311)

  options(np.messages = FALSE)

  n <- 24L
  tx <- data.frame(x1 = runif(n))
  tz <- data.frame(z1 = runif(n))
  y <- 0.7 * tx$x1 + sin(2 * pi * tz$z1) + rnorm(n, sd = 0.06)

  bws.nn <- matrix(c(2, 9), nrow = 2L)
  cases <- list(
    list(regtype = "lp", boot_method = "inid")
  )

  for (case in cases) {
    regtype <- case$regtype
    bw.args <- list(
      xdat = tx,
      ydat = y,
      zdat = tz,
      regtype = regtype,
      bwmethod = "cv.aic",
      bwtype = "generalized_nn",
      bws = bws.nn,
      bandwidth.compute = FALSE
    )

    if (identical(regtype, "lp")) {
      bw.args$basis <- "glp"
      bw.args$degree <- 1L
      bw.args$bernstein.basis <- FALSE
    }

    bw <- do.call(npplregbw, bw.args)
    out <- plot(
      bw,
      xdat = tx,
      ydat = y,
      zdat = tz,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = case$boot_method,
      plot.errors.boot.num = 3L,
      plot.errors.alpha = 0.1,
      perspective = FALSE
    )

    expect_true(is.list(out), info = sprintf("regtype=%s boot=%s", regtype, case$boot_method))
    expect_true(length(out) > 0L, info = sprintf("regtype=%s boot=%s", regtype, case$boot_method))
    expect_true(all(c("plr1", "plr2") %in% names(out)),
                info = sprintf("regtype=%s boot=%s", regtype, case$boot_method))
  }
})
