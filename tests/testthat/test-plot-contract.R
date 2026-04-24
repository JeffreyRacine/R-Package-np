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

test_that("plot contract: 3D plot-data matches data mode for regression and conditional estimators", {
  skip_if_not_installed("np")

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
  expect_warning(plot(rfit, plot.behavior = "data", view = "fixed"), NA)
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
  expect_warning(plot(cfit, plot.behavior = "data", view = "fixed"), NA)
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

test_that("plot contract: remaining public plot families return plot-data payloads", {
  skip_if_not_installed("np")

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

  fit.qr.alt <- npqreg(
    bws = bw.cf,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    tau = 0.4
  )
  bw.qr.alt.fixed <- plot(
    bw.cf,
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed",
    quantreg = TRUE,
    tau = 0.4
  )
  fit.qr.alt.fixed <- plot(
    fit.qr.alt,
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed"
  )
  bw.qr.alt.slice <- plot(
    bw.cf,
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    plot.behavior = "data",
    perspective = FALSE,
    quantreg = TRUE,
    tau = 0.4
  )
  fit.qr.alt.slice <- plot(
    fit.qr.alt,
    plot.behavior = "data",
    perspective = FALSE
  )

  expect_true(all(vapply(fit.qr.alt.fixed, inherits, logical(1), "qregression")))
  expect_true(all(vapply(fit.qr.alt.slice, inherits, logical(1), "qregression")))
  expect_true(all(vapply(fit.qr.alt.fixed, function(xi) identical(xi$tau, 0.4), logical(1))))
  expect_true(all(vapply(fit.qr.alt.slice, function(xi) identical(xi$tau, 0.4), logical(1))))
  expect_equal(fit.qr.alt.fixed[[1L]]$quantile, bw.qr.alt.fixed[[1L]]$quantile)
  expect_equal(fit.qr.alt.fixed[[1L]]$xeval, bw.qr.alt.fixed[[1L]]$xeval)
  expect_equal(fit.qr.alt.fixed[[1L]]$quanterr, bw.qr.alt.fixed[[1L]]$quanterr)
  expect_equal(fit.qr.alt.slice[[1L]]$quantile, bw.qr.alt.slice[[1L]]$quantile)
  expect_equal(fit.qr.alt.slice[[1L]]$xeval, bw.qr.alt.slice[[1L]]$xeval)
  expect_equal(fit.qr.alt.slice[[1L]]$quanterr, bw.qr.alt.slice[[1L]]$quanterr)

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
    tydat = y
  )
  expect_plot_modes_match(fit.sc, fields = c("eval", "mean", "grad"), perspective = FALSE)

  bw.pl <- npplregbw(
    xdat = data.frame(z = z),
    zdat = data.frame(x = x),
    ydat = y,
    bws = matrix(c(0.25, 0.25), nrow = 2),
    bandwidth.compute = FALSE
  )
  expect_plot_modes_match(bw.pl, fields = c("mean", "merr"), perspective = FALSE)

  fit.pl <- npplreg(bws = bw.pl, txdat = data.frame(z = z), tzdat = data.frame(x = x), tydat = y)
  expect_plot_modes_match(fit.pl, fields = c("mean", "merr"), perspective = FALSE)
})

test_that("plot contract: rgl plot-data returns the usual data payload", {
  skip_if_not_installed("np")
  skip_if_not_installed("rgl")

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

test_that("plot contract: scbandwidth perspective bootstrap uses z evaluation column", {
  skip_if_not_installed("np")

  set.seed(107)
  n <- 50
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.05)
  xdat <- data.frame(x = x)
  zdat <- data.frame(z = z)

  bw <- npscoefbw(xdat = xdat, ydat = y, zdat = zdat, regtype = "lc", nmulti = 1)
  neval <- 12L

  np.ns <- asNamespace("np")
  cap <- new.env(parent = emptyenv())
  cap$exdat <- NULL
  cap$ezdat <- NULL

  tq <- getFromNamespace("trim.quantiles", "np")
  x1.eval <- seq(tq(xdat[, 1L], 0.0)[1L], tq(xdat[, 1L], 0.0)[2L], length.out = neval)
  x2.eval <- seq(tq(zdat[, 1L], 0.0)[1L], tq(zdat[, 1L], 0.0)[2L], length.out = neval)
  x.eval <- expand.grid(x1.eval, x2.eval)
  colnames(x.eval) <- c("x", "z")

  trace(
    what = "compute.bootstrap.errors.scbandwidth",
    where = np.ns,
    tracer = bquote({
      assign("exdat", exdat, envir = .(cap))
      assign("ezdat", ezdat, envir = .(cap))
    }),
    print = FALSE
  )
  on.exit(try(untrace("compute.bootstrap.errors.scbandwidth", where = np.ns), silent = TRUE), add = TRUE)

  suppressWarnings(
    plot(
      bw,
      xdat = xdat,
      ydat = y,
      zdat = zdat,
      perspective = TRUE,
      view = "fixed",
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 9,
      plot.errors.type = "pointwise",
      neval = neval
    )
  )

  expect_equal(cap$exdat, x.eval[, 1, drop = FALSE])
  expect_equal(cap$ezdat, x.eval[, 2, drop = FALSE])
})

test_that("plot contract: npscoef fitted perspective path preserves semantic z eval name", {
  skip_if_not_installed("np")

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

  np.ns <- asNamespace("np")
  cap <- new.env(parent = emptyenv())
  cap$tzdat.names <- NULL
  cap$ezdat.names <- NULL

  trace(
    what = ".np_scoef_fit_internal",
    where = np.ns,
    tracer = bquote({
      if (!missing(tzdat) && !missing(ezdat) && !is.null(tzdat) && !is.null(ezdat)) {
        assign("tzdat.names", colnames(as.data.frame(tzdat)), envir = .(cap))
        assign("ezdat.names", colnames(as.data.frame(ezdat)), envir = .(cap))
      }
    }),
    print = FALSE
  )
  on.exit(try(untrace(".np_scoef_fit_internal", where = np.ns), silent = TRUE), add = TRUE)

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

  for (nm in c(reg.engines, unsup.engines)) {
    fn <- getFromNamespace(nm, "np")
    defaults <- eval(formals(fn)$plot.errors.boot.nonfixed)
    expect_identical(defaults[1L], "exact")
  }
})
