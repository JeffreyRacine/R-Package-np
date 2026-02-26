test_that("plot contract: regression and conditional estimators return data payloads", {
  skip_if_not_installed("np")

  set.seed(101)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- rnorm(n)

  rbw <- npregbw(y ~ x + z, nmulti = 1)
  rout <- suppressWarnings(
    plot(
      rbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(rout, "list")
  expect_true(length(rout) > 0)
  expect_true(all(vapply(rout, inherits, logical(1), "npregression")))
  expect_true(all(vapply(rout, function(xi) !is.null(xi$mean), logical(1))))
  expect_true(all(vapply(rout, function(xi) !is.null(xi$merr), logical(1))))
})

test_that("plot contract: npcdens/npcdist support bootstrap all in data mode", {
  skip_if_not_installed("np")

  set.seed(102)
  n <- 80
  x <- runif(n)
  y <- rnorm(n)

  cbw <- npcdensbw(y ~ x, nmulti = 1, bwmethod = "cv.ls")
  cout <- suppressWarnings(
    plot(
      cbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(cout, "list")
  expect_true(length(cout) > 0)
  expect_true(all(vapply(cout, inherits, logical(1), "condensity")))
  expect_true(all(vapply(cout, function(xi) !is.null(xi$condens), logical(1))))
  expect_true(all(vapply(cout, function(xi) !is.null(xi$conderr), logical(1))))

  dbw <- npcdistbw(y ~ x, nmulti = 1, bwmethod = "cv.ls")
  dout <- suppressWarnings(
    plot(
      dbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.num = 9,
      plot.errors.type = "all",
      view = "fixed"
    )
  )

  expect_type(dout, "list")
  expect_true(length(dout) > 0)
  expect_true(all(vapply(dout, inherits, logical(1), "condistribution")))
  expect_true(all(vapply(dout, function(xi) !is.null(xi$condist), logical(1))))
  expect_true(all(vapply(dout, function(xi) !is.null(xi$conderr), logical(1))))
})

test_that("plot contract: plot.errors.alpha is enforced in (0,0.5)", {
  skip_if_not_installed("np")

  set.seed(103)
  n <- 60
  x <- runif(n)
  y <- rnorm(n)
  bw <- npregbw(y ~ x, nmulti = 1)

  expect_error(
    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.alpha = 0.5
      )
    ),
    "must lie in \\(0,0\\.5\\)"
  )

  expect_error(
    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.alpha = 0
      )
    ),
    "must lie in \\(0,0\\.5\\)"
  )
})

test_that("plot contract: default condensity plot call stays scalar-safe", {
  skip_if_not_installed("np")

  set.seed(104)
  n <- 80
  x <- rnorm(n)
  y <- rnorm(n)

  bw <- npcdensbw(y ~ x, bws = c(1.0, 1.0), bandwidth.compute = FALSE)
  fit <- npcdens(bws = bw)

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_error(plot(fit, perspective = FALSE), NA)
})

test_that("plot contract: npscoef supports coef=TRUE plot path", {
  skip_if_not_installed("np")

  set.seed(105)
  n <- 60
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.15))

  fit <- npscoef(y ~ x | z, regtype = "ll", betas = TRUE)
  out <- suppressWarnings(
    plot(
      fit,
      coef = TRUE,
      coef.index = 1,
      perspective = FALSE,
      neval = 20,
      plot.behavior = "plot-data",
      plot.errors.method = "none"
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_true(all(vapply(out, inherits, logical(1), "smoothcoefficient")))
  expect_true(all(vapply(out, function(xi) !is.null(xi$mean), logical(1))))
})

test_that("plot contract: npplreg supports coef=TRUE plot-data payload", {
  skip_if_not_installed("np")

  set.seed(106)
  n <- 80
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- 1 + 0.7 * x + sin(z) + rnorm(n, sd = 0.15)
  xdat <- data.frame(x = x)
  zdat <- data.frame(z = z)

  fit <- npplreg(y ~ x | z, regtype = "ll")
  out <- suppressWarnings(
    plot(
      fit,
      xdat = xdat,
      ydat = y,
      zdat = zdat,
      coef = TRUE,
      plot.behavior = "plot-data",
      plot.errors.method = "none"
    )
  )

  expect_type(out, "list")
  expect_true(all(c("coefficients", "coefficient.stderr", "fit") %in% names(out)))
  expect_true(is.numeric(out$coefficients))
  expect_true(length(out$coefficients) >= 1L)
})

test_that("plot contract: bootstrap defaults are wild for regression-class and inid for unsupervised engines", {
  skip_if_not_installed("np")

  reg.engines <- c(
    ".np_plot_rbandwidth_engine",
    ".np_plot_scbandwidth_engine",
    ".np_plot_plbandwidth_engine",
    ".np_plot_sibandwidth_engine"
  )
  unsup.engines <- c(
    ".np_plot_bandwidth_engine",
    ".np_plot_dbandwidth_engine",
    ".np_plot_conbandwidth_engine",
    ".np_plot_condbandwidth_engine"
  )

  for (nm in reg.engines) {
    fn <- getFromNamespace(nm, "np")
    defaults <- eval(formals(fn)$plot.errors.boot.method)
    expect_identical(defaults[1L], "wild")
  }

  for (nm in unsup.engines) {
    fn <- getFromNamespace(nm, "np")
    defaults <- eval(formals(fn)$plot.errors.boot.method)
    expect_identical(defaults[1L], "inid")
  }
})
