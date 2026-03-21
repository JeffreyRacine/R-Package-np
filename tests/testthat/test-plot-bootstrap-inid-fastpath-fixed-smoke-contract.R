test_that("npscoef plot bootstrap inid supports ll/lp basis variants", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(3233)
  n <- 70
  x <- runif(n)
  z <- runif(n)
  y <- (0.6 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.08)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      zdat = tz,
      ydat = y,
      bws = 0.22,
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npscoefbw, bw.args)
    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        zdat = tz,
        perspective = FALSE,
        plot.behavior = "data",
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0, info = cfg$label)
  }
})

test_that("density/distribution plot bootstrap rejects wild selector", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(328)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n)

  ubw <- npudensbw(dat = data.frame(x = x), bws = 0.5, bandwidth.compute = FALSE)
  expect_error(
    suppressWarnings(plot(
      ubw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )),
    "not supported for unconditional density/distribution estimators"
  )

  cbw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.6, 0.6),
    bandwidth.compute = FALSE
  )
  expect_error(
    suppressWarnings(plot(
      cbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )),
    "not supported for conditional density/distribution estimators"
  )
})

test_that("plot contract defaults nonfixed bootstrap mode to exact", {
  engines <- c(
    ".np_plot_rbandwidth_engine",
    ".np_plot_scbandwidth_engine",
    ".np_plot_plbandwidth_engine",
    ".np_plot_sibandwidth_engine",
    ".np_plot_bandwidth_engine",
    ".np_plot_dbandwidth_engine",
    ".np_plot_conbandwidth_engine",
    ".np_plot_condbandwidth_engine"
  )

  for (nm in engines) {
    fn <- getFromNamespace(nm, "npRmpi")
    defaults <- eval(formals(fn)$plot.errors.boot.nonfixed)
    expect_identical(defaults[1L], "exact")
  }
})
