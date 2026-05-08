test_that("npplreg generalized-nn cv.aic bootstrap plot contract holds across regtypes and methods", {
  set.seed(20260311)

  n <- 60L
  tx <- data.frame(x1 = runif(n))
  tz <- data.frame(z1 = runif(n))
  y <- 0.7 * tx$x1 + sin(2 * pi * tz$z1) + rnorm(n, sd = 0.06)

  boot_methods <- c("inid", "fixed", "geom", "wild")
  regtypes <- c("lc", "ll", "lp")
  bws.nn <- matrix(c(2, 9), nrow = 2L)

  for (regtype in regtypes) {
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

    for (boot_method in boot_methods) {
      out <- plot(
        bw,
        xdat = tx,
        ydat = y,
        zdat = tz,
        output = "data",
        errors = "bootstrap",
        bootstrap = boot_method,
        B = 7L,
        alpha = 0.1,
        perspective = FALSE
      )

      expect_true(is.list(out), info = sprintf("regtype=%s boot=%s", regtype, boot_method))
      expect_true(length(out) > 0L, info = sprintf("regtype=%s boot=%s", regtype, boot_method))
      expect_true(all(c("plr1", "plr2") %in% names(out)),
                  info = sprintf("regtype=%s boot=%s", regtype, boot_method))
    }
  }
})
