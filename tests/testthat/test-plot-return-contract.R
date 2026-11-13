library(npRmpi)

test_that("plot return contract: 3D plot-data matches data mode for regression and conditional estimators", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  with_plot_device <- function(expr) {
    pdf(file = tempfile(fileext = ".pdf"))
    on.exit(dev.off(), add = TRUE)
    force(expr)
  }

  set.seed(107)
  n <- 80
  x <- rnorm(n)
  z <- rnorm(n)
  y <- x - 0.5 * z + rnorm(n, sd = 0.2)

  rfit <- npreg(
    bws = npregbw(
      xdat = data.frame(x = x, z = z),
      ydat = y,
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE
    ),
    txdat = data.frame(x = x, z = z),
    tydat = y
  )
  rdata <- plot(rfit, output = "data", view = "fixed")
  rplotdata <- with_plot_device(plot(rfit, output = "plot-data", view = "fixed"))

  expect_s3_class(rplotdata, "np_plot_data")
  expect_equal(rplotdata$panels, rdata$panels)
  expect_equal(rplotdata$targets, rdata$targets)
  expect_equal(rplotdata$inference, rdata$inference)
  expect_equal(rplotdata$coverage_families, rdata$coverage_families)

  cfit <- npcdens(
    bws = npcdensbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein.basis = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  undercoverage <- "the supplied y-grid may under-cover response support"
  cdata <- npRmpi_expect_warning_messages(
    plot(cfit, output = "data", view = "fixed"),
    undercoverage
  )
  cplotdata <- npRmpi_expect_warning_messages(
    with_plot_device(plot(cfit, output = "plot-data", view = "fixed")),
    undercoverage
  )

  expect_s3_class(cplotdata, "np_plot_data")
  expect_equal(cplotdata$panels, cdata$panels)
  expect_equal(cplotdata$targets, cdata$targets)
  expect_equal(cplotdata$inference, cdata$inference)
  expect_equal(cplotdata$coverage_families, cdata$coverage_families)

  dfit <- npcdist(
    bws = npcdistbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein.basis = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  ddata <- npRmpi_expect_warning_messages(
    plot(dfit, output = "data", view = "fixed"),
    character()
  )
  dplotdata <- npRmpi_expect_warning_messages(
    with_plot_device(plot(dfit, output = "plot-data", view = "fixed")),
    character()
  )

  expect_s3_class(dplotdata, "np_plot_data")
  expect_equal(dplotdata$panels, ddata$panels)
  expect_equal(dplotdata$targets, ddata$targets)
  expect_equal(dplotdata$inference, ddata$inference)
  expect_equal(dplotdata$coverage_families, ddata$coverage_families)
})

test_that("plot return contract: bounded conditional bootstrap data stays available", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(109)
  n <- 60
  x <- runif(n)
  y <- runif(n)

  cfit <- npcdens(y ~ x, cykerbound = "range")
  cout <- suppressWarnings(plot(
    cfit,
    output = "data",
    view = "fixed",
    data.overlay = FALSE,
    errors = "bootstrap",
    bootstrap = "inid",
    B = 5L,
    band = "pointwise"
  ))

  expect_s3_class(cout, "np_plot_data")
  expect_gt(nrow(plot_target(cout, band = "pointwise")), 0L)
  expect_true(all(is.finite(plot_target(cout)$estimate)))
  expect_true(all(is.finite(plot_band(cout, "pointwise"))))

  dfit <- npcdist(y ~ x, cykerbound = "range")
  dout <- suppressWarnings(plot(
    dfit,
    output = "data",
    view = "fixed",
    data.overlay = FALSE,
    errors = "bootstrap",
    bootstrap = "inid",
    B = 5L,
    band = "pointwise"
  ))

  expect_s3_class(dout, "np_plot_data")
  expect_gt(nrow(plot_target(dout, band = "pointwise")), 0L)
  expect_true(all(is.finite(plot_target(dout)$estimate)))
  expect_true(all(is.finite(plot_band(dout, "pointwise"))))
})

test_that("plot return contract: remaining public plot families return plot-data payloads", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  with_plot_device <- function(expr) {
    pdf(file = tempfile(fileext = ".pdf"))
    on.exit(dev.off(), add = TRUE)
    suppressWarnings(force(expr))
  }

  expect_plot_modes_match <- function(object, ...) {
    data.out <- suppressWarnings(plot(object, output = "data", ...))
    plot.out <- with_plot_device(plot(object, output = "plot-data", ...))

    expect_s3_class(plot.out, "np_plot_data")
    expect_s3_class(data.out, "np_plot_data")
    expect_gt(length(plot.out$panels), 0L)
    expect_equal(plot.out$panels, data.out$panels)
    expect_equal(plot.out$targets, data.out$targets)
    expect_equal(plot.out$inference, data.out$inference)
    expect_equal(plot.out$coverage_families, data.out$coverage_families)
  }

  set.seed(108)
  n <- 60
  x <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.5 * x2 + 1.2 * z + rnorm(n, sd = 0.05)

  bw.ud <- npudensbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.ud, perspective = FALSE)

  fit.ud <- npudens(bws = bw.ud, tdat = data.frame(x = x))
  expect_plot_modes_match(fit.ud, perspective = FALSE)

  bw.ui <- npudistbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.ui, perspective = FALSE)

  fit.ui <- npudist(bws = bw.ui, tdat = data.frame(x = x))
  expect_plot_modes_match(fit.ui, perspective = FALSE)

  bw.reg <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.reg, perspective = FALSE)
  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.reg.nn <- do.call(npregbw, list(
      xdat = data.frame(x = x),
      ydat = y,
      bws = 7,
      bwtype = bt,
      bandwidth.compute = FALSE
    ))
    fit.reg.nn <- npreg(bws = bw.reg.nn, txdat = data.frame(x = x), tydat = y)
    expect_plot_modes_match(
      fit.reg.nn,

      perspective = FALSE,
      errors = "none"
    )
  }
  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(
    bw.cd,

    perspective = TRUE,
    view = "fixed"
  )

  bw.cf <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(
    bw.cf,

    perspective = TRUE,
    view = "fixed"
  )

  fit.qr <- npqreg(
    bws = bw.cf,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    tau = 0.5
  )
  expect_plot_modes_match(
    fit.qr,

    perspective = TRUE,
    view = "fixed"
  )

  bw.si <- npindexbw(
    xdat = data.frame(x = x, x2 = x2),
    ydat = y,
    bws = c(0.25, 0.25, 1),
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.si)

  fit.si <- npindex(bws = bw.si, txdat = data.frame(x = x, x2 = x2), tydat = y)
  expect_plot_modes_match(fit.si)

  bw.sc <- npscoefbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.sc, perspective = FALSE)

  fit.sc <- npscoef(
    bws = bw.sc,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    se = FALSE
  )
  expect_plot_modes_match(fit.sc, perspective = FALSE)

  fit.pl <- npplreg(
    bws = npplregbw(
      xdat = data.frame(z = z),
      zdat = data.frame(x = x),
      ydat = y,
      bws = matrix(c(0.25, 0.25), nrow = 2),
      bandwidth.compute = FALSE
    ),
    txdat = data.frame(z = z),
    tzdat = data.frame(x = x),
    tydat = y
  )
  expect_plot_modes_match(fit.pl, perspective = FALSE)
})

test_that("plot return contract: npRmpi npscoef fitted perspective path preserves semantic z eval name", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(108)
  n <- 120
  x <- runif(n)
  z <- runif(n, min = -2, max = 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.15))

  fit <- npscoef(
    y ~ x | z,
    nomad = TRUE,
    nmulti = 1L
  )

  npRmpi.ns <- asNamespace("npRmpi")
  cap <- new.env(parent = emptyenv())
  cap$tzdat.names <- NULL
  cap$ezdat.names <- NULL

  trace(
    what = ".np_scoef_fit_internal",
    where = npRmpi.ns,
    tracer = bquote({
      if (!missing(tzdat) && !missing(ezdat) && !is.null(tzdat) && !is.null(ezdat)) {
        assign("tzdat.names", colnames(as.data.frame(tzdat)), envir = .(cap))
        assign("ezdat.names", colnames(as.data.frame(ezdat)), envir = .(cap))
      }
    }),
    print = FALSE
  )
  on.exit(try(untrace(".np_scoef_fit_internal", where = npRmpi.ns), silent = TRUE), add = TRUE)

  out <- suppressWarnings(plot(
    fit,
    view = "fixed",
    renderer = "base",
    output = "data",
    data.overlay = FALSE,
    errors = "none"
  ))

  expect_type(out, "list")
  expect_identical(cap$tzdat.names, "z")
  expect_identical(cap$ezdat.names, "z")
})

test_that("plot return contract: rgl plot-data returns the usual data payload", {
  old.opts <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.opts), add = TRUE)
  skip_if_not_installed("rgl")
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(110)
  n <- 50
  x <- rnorm(n)
  z <- rnorm(n)
  y <- x - 0.5 * z + rnorm(n, sd = 0.1)

  rfit <- npreg(
    bws = npregbw(
      xdat = data.frame(x = x, z = z),
      ydat = y,
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE
    ),
    txdat = data.frame(x = x, z = z),
    tydat = y
  )

  rdata <- npRmpi_expect_warning_messages(
    plot(
      rfit,
      output = "data",
      renderer = "rgl",
      view = "fixed",
      data.overlay = FALSE
    ),
    character()
  )

  rplotdata <- npRmpi_expect_warning_messages(
    plot(
      rfit,
      output = "plot-data",
      renderer = "rgl",
      view = "fixed",
      data.overlay = FALSE
    ),
    character()
  )

  expect_s3_class(rplotdata, "np_plot_data")
  expect_equal(rplotdata$panels, rdata$panels)
  expect_equal(rplotdata$targets, rdata$targets)
  expect_equal(rplotdata$inference, rdata$inference)
  expect_equal(rplotdata$coverage_families, rdata$coverage_families)

  cfit <- npcdens(
    bws = npcdensbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein.basis = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  undercoverage <- "the supplied y-grid may under-cover response support"
  cdata <- npRmpi_expect_warning_messages(
    plot(
      cfit,
      output = "data",
      renderer = "rgl",
      view = "fixed"
    ),
    undercoverage
  )

  cplotdata <- npRmpi_expect_warning_messages(
    plot(
      cfit,
      output = "plot-data",
      renderer = "rgl",
      view = "fixed"
    ),
    undercoverage
  )

  expect_s3_class(cplotdata, "np_plot_data")
  expect_equal(cplotdata$panels, cdata$panels)
  expect_equal(cplotdata$targets, cdata$targets)
  expect_equal(cplotdata$inference, cdata$inference)
  expect_equal(cplotdata$coverage_families, cdata$coverage_families)
})
