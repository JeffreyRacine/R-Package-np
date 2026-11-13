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
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      B = 9,
      band = "all",
      view = "fixed"
    )
  )

  expect_s3_class(rout, "np_plot_data")
  expect_true(length(rout$panels) > 0)
  expect_true(all(rout$targets$quantity == "mean"))
  expect_true(all(vapply(rout$panels, function(p) nrow(p$values) > 0L &&
    setequal(unique(p$intervals$band), c("pointwise", "simultaneous", "bonferroni")), logical(1L))))
})

test_that("plot contract: legend defaults are role-sized and user-overridable", {
  legend_args <- getFromNamespace(".np_plot_legend_args", "np")

  defaults <- list(x = "topright", legend = "A", cex = 0.8, bty = "n")
  out <- legend_args(defaults, legend = TRUE)
  expect_equal(out$cex, 0.8)
  expect_equal(out$bty, "n")

  overridden <- legend_args(
    defaults,
    legend = list(x = "bottomleft", cex = 0.55, bty = "o")
  )
  expect_equal(overridden$x, "bottomleft")
  expect_equal(overridden$cex, 0.55)
  expect_equal(overridden$bty, "o")
  expect_equal(overridden$legend, "A")
})

test_that("plot contract: package legends are user-controllable on band-all plots", {
  set.seed(20260511)
  n <- 36L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  rbw <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.35,
    bandwidth.compute = FALSE
  )

  old.dev <- grDevices::dev.cur()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit({
    grDevices::dev.off()
    if (old.dev > 1L)
      grDevices::dev.set(old.dev)
  }, add = TRUE)

  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = FALSE
  )))
  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = "bottomright"
  )))
  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = list(x = "bottomleft", bty = "o", cex = 0.55)
  )))
  expect_error(
    suppressWarnings(plot(
      rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
      errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
      neval = 6L, legend = 1
    )),
    "legend must be TRUE/FALSE"
  )
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
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      B = 9,
      band = "all",
      view = "fixed"
    )
  )

  expect_s3_class(cout, "np_plot_data")
  expect_true(length(cout$panels) > 0)
  expect_true(all(cout$targets$quantity == "conditional_density"))
  expect_true(all(vapply(cout$panels, function(p) nrow(p$values) > 0L &&
    setequal(unique(p$intervals$band), c("pointwise", "simultaneous", "bonferroni")), logical(1L))))

  dbw <- npcdistbw(y ~ x, nmulti = 1, bwmethod = "cv.ls")
  dout <- suppressWarnings(
    plot(
      dbw,
      output = "data",
      perspective = FALSE,
      errors = "bootstrap",
      B = 9,
      band = "all",
      view = "fixed"
    )
  )

  expect_s3_class(dout, "np_plot_data")
  expect_true(length(dout$panels) > 0)
  expect_true(all(dout$targets$quantity == "conditional_cdf"))
  expect_true(all(vapply(dout$panels, function(p) nrow(p$values) > 0L &&
    setequal(unique(p$intervals$band), c("pointwise", "simultaneous", "bonferroni")), logical(1L))))

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
  rdata <- plot(rfit, output = "data", view = "fixed")
  expect_warning(plot(rfit, output = "data", view = "fixed"), NA)
  rplotdata <- with_plot_device(plot(rfit, output = "plot-data", view = "fixed"))

  expect_type(rplotdata, "list")
  expect_named(rplotdata, names(rdata))
  expect_s3_class(rplotdata, "np_plot_data")
  expect_identical(rplotdata, rdata)

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
  cdata <- suppressWarnings(plot(cfit, output = "data", view = "fixed"))
  expect_warning(plot(cfit, output = "data", view = "fixed"), "proper=TRUE projected")
  cplotdata <- with_plot_device(plot(cfit, output = "plot-data", view = "fixed"))

  expect_type(cplotdata, "list")
  expect_named(cplotdata, names(cdata))
  expect_s3_class(cplotdata, "np_plot_data")
  expect_identical(cplotdata, cdata)

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
  ddata <- suppressWarnings(plot(dfit, output = "data", view = "fixed"))
  dplotdata <- with_plot_device(plot(dfit, output = "plot-data", view = "fixed"))

  expect_type(dplotdata, "list")
  expect_named(dplotdata, names(ddata))
  expect_s3_class(dplotdata, "np_plot_data")
  expect_identical(dplotdata, ddata)
})

test_that("plot contract: remaining public plot families return plot-data payloads", {
  skip_if_not_installed("np")

  with_plot_device <- function(expr) {
    pdf(file = tempfile(fileext = ".pdf"))
    on.exit(dev.off(), add = TRUE)
    suppressWarnings(force(expr))
  }

  expect_plot_modes_match <- function(object, fields, ...) {
    data.out <- suppressWarnings(plot(object, output = "data", ...))
    plot.out <- with_plot_device(plot(object, output = "plot-data", ...))

    expect_type(plot.out, "list")
    expect_named(plot.out, names(data.out))

    expect_s3_class(plot.out, "np_plot_data")
    expect_identical(plot.out, data.out)
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
    output = "data",
    perspective = TRUE,
    view = "fixed",
    quantreg = TRUE,
    tau = 0.4
  )
  fit.qr.alt.fixed <- plot(
    fit.qr.alt,
    output = "data",
    perspective = TRUE,
    view = "fixed"
  )
  bw.qr.alt.slice <- plot(
    bw.cf,
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    output = "data",
    perspective = FALSE,
    quantreg = TRUE,
    tau = 0.4
  )
  fit.qr.alt.slice <- plot(
    fit.qr.alt,
    output = "data",
    perspective = FALSE
  )

  for (pair in list(list(fit.qr.alt.fixed, bw.qr.alt.fixed),
                    list(fit.qr.alt.slice, bw.qr.alt.slice))) {
    expect_s3_class(pair[[1L]], "np_plot_data")
    expect_true(all(pair[[1L]]$targets$tau == .4))
    # Fitted-quantile and bandwidth-quantreg routes have distinct incumbent
    # overlay defaults. Their numerical tables must agree exactly.
    expect_equal(lapply(pair[[1L]]$panels, function(p) p[c("grid", "variables", "values", "intervals")]),
                 lapply(pair[[2L]]$panels, function(p) p[c("grid", "variables", "values", "intervals")]))
    expect_equal(pair[[1L]]$targets, pair[[2L]]$targets)
  }

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
    output = "data",
    renderer = "rgl",
    view = "fixed",
    data.overlay = FALSE
  ))

  rplotdata <- suppressWarnings(plot(
    rfit,
    output = "plot-data",
    renderer = "rgl",
    view = "fixed",
    data.overlay = FALSE
  ))

  expect_type(rplotdata, "list")
  expect_named(rplotdata, names(rdata))
  expect_s3_class(rplotdata, "np_plot_data")
  expect_identical(rplotdata, rdata)

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

  cdata <- suppressWarnings(plot(
    cfit,
    output = "data",
    renderer = "rgl",
    view = "fixed"
  ))

  cplotdata <- suppressWarnings(plot(
    cfit,
    output = "plot-data",
    renderer = "rgl",
    view = "fixed"
  ))

  expect_type(cplotdata, "list")
  expect_named(cplotdata, names(cdata))
  expect_s3_class(cplotdata, "np_plot_data")
  expect_identical(cplotdata, cdata)
})

test_that("plot contract: alpha is enforced in (0,0.5)", {
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
        output = "data",
        perspective = FALSE,
        errors = "bootstrap",
        alpha = 0.5
      )
    ),
    "must lie in \\(0, 0\\.5\\)"
  )

  expect_error(
    suppressWarnings(
      plot(
        bw,
        output = "data",
        perspective = FALSE,
        errors = "bootstrap",
        alpha = 0
      )
    ),
    "must lie in \\(0, 0\\.5\\)"
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
  withr::local_pdf(NULL)

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
      output = "plot-data",
      errors = "none"
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_s3_class(out, "np_plot_data")
  expect_true(all(out$targets$quantity == "coefficient"))
  expect_true(all(vapply(out$panels, function(p) nrow(p$values) > 0L, logical(1L))))
})

test_that("plot contract: npplreg supports coef=TRUE plot-data payload", {
  skip_if_not_installed("np")
  withr::local_pdf(NULL)

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
      output = "plot-data",
      errors = "none"
    )
  )

  expect_type(out, "list")
  expect_s3_class(out, "np_plot_data")
  d <- as.data.frame(out)
  expect_true(is.numeric(d$estimate))
  expect_equal(d$estimate, as.vector(coef(fit)))
  expect_gte(nrow(d), 1L)
  expect_false("fit" %in% names(out))
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
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 9,
      band = "pointwise",
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
    output = "data",
    data.overlay = FALSE,
    errors = "none"
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
