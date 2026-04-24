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
    suppressWarnings(force(expr))
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
  rdata <- plot(rfit, plot.behavior = "data", view = "fixed")
  rplotdata <- with_plot_device(plot(rfit, plot.behavior = "plot-data", view = "fixed"))

  expect_type(rplotdata, "list")
  expect_named(rplotdata, names(rdata))
  expect_s3_class(rplotdata$r1, "npregression")
  expect_equal(rplotdata$r1$mean, rdata$r1$mean)
  expect_equal(rplotdata$r1$eval, rdata$r1$eval)

  cfit <- npcdens(
    bws = npcdensbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  cdata <- plot(cfit, plot.behavior = "data", view = "fixed")
  cplotdata <- with_plot_device(plot(cfit, plot.behavior = "plot-data", view = "fixed"))

  expect_type(cplotdata, "list")
  expect_named(cplotdata, names(cdata))
  expect_s3_class(cplotdata$cd1, "condensity")
  expect_equal(cplotdata$cd1$condens, cdata$cd1$condens)
  expect_equal(cplotdata$cd1$condens.raw, cdata$cd1$condens.raw)
  expect_equal(cplotdata$cd1$xeval, cdata$cd1$xeval)
  expect_equal(cplotdata$cd1$yeval, cdata$cd1$yeval)

  dfit <- npcdist(
    bws = npcdistbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )
  ddata <- suppressWarnings(plot(dfit, plot.behavior = "data", view = "fixed"))
  dplotdata <- with_plot_device(plot(dfit, plot.behavior = "plot-data", view = "fixed"))

  expect_type(dplotdata, "list")
  expect_named(dplotdata, names(ddata))
  expect_s3_class(dplotdata$cd1, "condistribution")
  expect_equal(dplotdata$cd1$condist, ddata$cd1$condist)
  expect_equal(dplotdata$cd1$condist.raw, ddata$cd1$condist.raw)
  expect_equal(dplotdata$cd1$xeval, ddata$cd1$xeval)
  expect_equal(dplotdata$cd1$yeval, ddata$cd1$yeval)
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
    plot.behavior = "data",
    view = "fixed",
    plot.data.overlay = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 5L,
    plot.errors.type = "pointwise"
  ))

  expect_s3_class(cout$cd1, "condensity")
  expect_true(all(is.finite(cout$cd1$condens)))
  expect_true(all(is.finite(cout$cd1$conderr)))

  dfit <- npcdist(y ~ x, cykerbound = "range")
  dout <- suppressWarnings(plot(
    dfit,
    plot.behavior = "data",
    view = "fixed",
    plot.data.overlay = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 5L,
    plot.errors.type = "pointwise"
  ))

  expect_s3_class(dout$cd1, "condistribution")
  expect_true(all(is.finite(dout$cd1$condist)))
  expect_true(all(is.finite(dout$cd1$conderr)))
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

  expect_plot_modes_match <- function(object, fields, ...) {
    data.out <- suppressWarnings(plot(object, plot.behavior = "data", ...))
    plot.out <- with_plot_device(plot(object, plot.behavior = "plot-data", ...))

    expect_type(plot.out, "list")
    expect_named(plot.out, names(data.out))

    for (nm in names(data.out)) {
      expect_s3_class(plot.out[[nm]], class(data.out[[nm]])[1L])
      for (field in fields) {
        if (field %in% names(data.out[[nm]]) && field %in% names(plot.out[[nm]])) {
          expect_equal(plot.out[[nm]][[field]], data.out[[nm]][[field]])
        }
      }
    }
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
  expect_plot_modes_match(bw.ud, fields = c("eval", "dens"), perspective = FALSE)

  fit.ud <- npudens(bws = bw.ud, tdat = data.frame(x = x))
  expect_plot_modes_match(fit.ud, fields = c("eval", "dens"), perspective = FALSE)

  bw.ui <- npudistbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.ui, fields = c("eval", "dist"), perspective = FALSE)

  fit.ui <- npudist(bws = bw.ui, tdat = data.frame(x = x))
  expect_plot_modes_match(fit.ui, fields = c("eval", "dist"), perspective = FALSE)

  bw.reg <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.reg, fields = c("eval", "mean", "merr"), perspective = FALSE)
  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.reg.nn <- npregbw(
      xdat = data.frame(x = x),
      ydat = y,
      bws = 7,
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    fit.reg.nn <- npreg(bws = bw.reg.nn, txdat = data.frame(x = x), tydat = y)
    expect_plot_modes_match(
      fit.reg.nn,
      fields = c("eval", "mean", "merr"),
      perspective = FALSE,
      plot.errors.method = "none"
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
    fields = c("condens", "xeval", "yeval", "conderr"),
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
    fields = c("condist", "xeval", "yeval", "conderr"),
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
    fields = c("quantile", "xeval", "quanterr"),
    perspective = TRUE,
    view = "fixed"
  )

  bw.si <- npindexbw(
    xdat = data.frame(x = x, x2 = x2),
    ydat = y,
    bws = c(0.25, 0.25, 1),
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.si, fields = c("index", "mean", "grad"))

  fit.si <- npindex(bws = bw.si, txdat = data.frame(x = x, x2 = x2), tydat = y)
  expect_plot_modes_match(fit.si, fields = c("index", "mean", "grad"))

  bw.sc <- npscoefbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.sc, fields = c("eval", "mean", "grad"), perspective = FALSE)

  fit.sc <- npscoef(
    bws = bw.sc,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    tydat = y,
    errors = FALSE
  )
  expect_plot_modes_match(fit.sc, fields = c("eval", "mean", "grad"), perspective = FALSE)

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
  expect_plot_modes_match(fit.pl, fields = c("mean", "merr"), perspective = FALSE)
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
    plot.behavior = "data",
    plot.data.overlay = FALSE,
    plot.errors.method = "none"
  ))

  expect_type(out, "list")
  expect_identical(cap$tzdat.names, "z")
  expect_identical(cap$ezdat.names, "z")
})

test_that("plot return contract: rgl plot-data returns the usual data payload", {
  skip_if_not_installed("rgl")
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  old.opts <- options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
  on.exit(options(old.opts), add = TRUE)

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

  rdata <- suppressWarnings(plot(
    rfit,
    plot.behavior = "data",
    renderer = "rgl",
    view = "fixed",
    plot.data.overlay = FALSE
  ))

  rplotdata <- suppressWarnings(plot(
    rfit,
    plot.behavior = "plot-data",
    renderer = "rgl",
    view = "fixed",
    plot.data.overlay = FALSE
  ))

  expect_type(rplotdata, "list")
  expect_named(rplotdata, names(rdata))
  expect_s3_class(rplotdata$r1, "npregression")
  expect_equal(rplotdata$r1$mean, rdata$r1$mean)
  expect_equal(rplotdata$r1$eval, rdata$r1$eval)

  cfit <- npcdens(
    bws = npcdensbw(
      xdat = data.frame(x = x),
      ydat = data.frame(y = y),
      bws = c(0.55, 0.55),
      bandwidth.compute = FALSE,
      regtype = "lp",
      degree = 4L,
      bernstein = TRUE
    ),
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  cdata <- suppressWarnings(plot(
    cfit,
    plot.behavior = "data",
    renderer = "rgl",
    view = "fixed"
  ))

  cplotdata <- suppressWarnings(plot(
    cfit,
    plot.behavior = "plot-data",
    renderer = "rgl",
    view = "fixed"
  ))

  expect_type(cplotdata, "list")
  expect_named(cplotdata, names(cdata))
  expect_s3_class(cplotdata$cd1, "condensity")
  expect_equal(cplotdata$cd1$condens, cdata$cd1$condens)
  expect_equal(cplotdata$cd1$xeval, cdata$cd1$xeval)
  expect_equal(cplotdata$cd1$yeval, cdata$cd1$yeval)
})
