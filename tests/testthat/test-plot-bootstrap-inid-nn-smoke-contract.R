test_that("nearest-neighbor plot helpers run for regression and density/distribution families", {
  skip_if_not_installed("np")

  set.seed(329)
  n <- 45
  x <- data.frame(x = runif(n))
  y <- sin(2 * pi * x$x) + rnorm(n, sd = 0.08)
  yframe <- data.frame(y = y)

  run_plot <- function(bw, ...) {
    suppressWarnings(plot(
      bw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5,
      ...
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    rbw <- npregbw(xdat = x, ydat = y, regtype = "ll", nmulti = 1, bwtype = bt)
    expect_type(run_plot(rbw, xdat = x, ydat = y), "list")
    expect_type(run_plot(rbw, xdat = x, ydat = y, gradients = TRUE), "list")

    ubw <- npudensbw(dat = x, nmulti = 1, bwtype = bt)
    expect_type(run_plot(ubw), "list")

    dbw <- npudistbw(dat = x, nmulti = 1, bwtype = bt)
    expect_type(run_plot(dbw), "list")

    cbw <- npcdensbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
    expect_type(run_plot(cbw, xdat = x, ydat = yframe, view = "fixed"), "list")

    cdbw <- npcdistbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
    expect_type(run_plot(cdbw, xdat = x, ydat = yframe, view = "fixed"), "list")
  }
})

test_that("nearest-neighbor frozen bootstrap plots run across regression and unsupervised families", {
  skip_if_not_installed("np")

  set.seed(32908)
  n <- 45
  x <- data.frame(x = runif(n))
  y <- sin(2 * pi * x$x) + rnorm(n, sd = 0.08)
  yframe <- data.frame(y = y)

  run_plot <- function(bw, ...) {
    suppressWarnings(plot(
      bw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.nonfixed = "frozen",
      plot.errors.boot.num = 5,
      ...
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    rbw <- npregbw(xdat = x, ydat = y, regtype = "ll", nmulti = 1, bwtype = bt)

    for (boot.method in c("inid", "fixed", "geom")) {
      expect_type(
        run_plot(rbw, xdat = x, ydat = y, plot.errors.boot.method = boot.method),
        "list"
      )

      ubw <- npudensbw(dat = x, nmulti = 1, bwtype = bt)
      expect_type(run_plot(ubw, plot.errors.boot.method = boot.method), "list")

      dbw <- npudistbw(dat = x, nmulti = 1, bwtype = bt)
      expect_type(run_plot(dbw, plot.errors.boot.method = boot.method), "list")

      cbw <- npcdensbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
      expect_type(
        run_plot(cbw, xdat = x, ydat = yframe, view = "fixed", plot.errors.boot.method = boot.method),
        "list"
      )

      cdbw <- npcdistbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
      expect_type(
        run_plot(cdbw, xdat = x, ydat = yframe, view = "fixed", plot.errors.boot.method = boot.method),
        "list"
      )
    }
  }
})
