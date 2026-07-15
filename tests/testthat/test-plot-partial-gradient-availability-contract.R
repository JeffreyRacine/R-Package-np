np_partial_gradient_local <- function(expr) {
  expr <- substitute(expr)
  if ("npRmpi" %in% loadedNamespaces()) {
    local <- get(".npRmpi_with_local_regression", envir = asNamespace("npRmpi"))
    return(eval(substitute(LOCAL(EXPR), list(LOCAL = local, EXPR = expr)),
                envir = parent.frame()))
  }
  eval(expr, envir = parent.frame())
}

test_that("mixed-degree gradient plots retain empty and categorical panels", {
  old <- options(np.messages = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  set.seed(20260715)
  n <- 42L
  x <- data.frame(
    x1 = runif(n, -1, 1),
    f = factor(rep(c("a", "b"), length.out = n)),
    x2 = runif(n, -1, 1)
  )
  y <- x$x1^2 + x$x2 + 0.3 * (x$f == "b") + rnorm(n, sd = 0.03)
  bw <- np_partial_gradient_local(npregbw(
    xdat = x, ydat = y, bws = c(0.4, 0.5, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(1L, 2L),
    basis = "glp"
  ))

  expect_warning(
    out <- np_partial_gradient_local(plot(
      bw, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = c(2L, 1L), errors = "none", output = "data",
      neval = 6L
    )),
    "x1.*requested order 2.*degree 1"
  )
  expect_length(out, 3L)
  expect_length(out[[1L]]$grad, 6L)
  expect_true(all(is.na(out[[1L]]$grad)))
  expect_true(any(is.finite(out[[2L]]$grad)))
  expect_true(any(is.finite(out[[3L]]$grad)))

  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    np_partial_gradient_local(suppressWarnings(plot(
      bw, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = c(2L, 1L), errors = "none", neval = 6L,
      common.scale = FALSE
    )))
    np_partial_gradient_local(suppressWarnings(plot(
      bw, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = c(2L, 1L), errors = "none", neval = 6L,
      common.scale = TRUE
    )))
  })
})

test_that("an all-unavailable gradient plot fails before RNG or device mutation", {
  old <- options(np.messages = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  set.seed(20260716)
  n <- 36L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- x$x1 + x$x2 + rnorm(n, sd = 0.04)
  bw <- np_partial_gradient_local(npregbw(
    xdat = x, ydat = y, bws = c(0.4, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(1L, 1L),
    basis = "glp"
  ))
  seed.before <- .Random.seed
  device.before <- grDevices::dev.cur()

  expect_error(
    np_partial_gradient_local(plot(
      bw, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = c(2L, 2L), errors = "none", output = "data"
    )),
    paste0(
      "no requested component is available.*",
      "x1.*requested order 2.*degree 1.*",
      "x2.*requested order 2.*degree 1.*",
      "Lower gradient.order.*refit with sufficient polynomial degree"
    )
  )
  expect_identical(.Random.seed, seed.before)
  expect_identical(grDevices::dev.cur(), device.before)
})

test_that("derived local-polynomial families preserve partial availability", {
  old <- options(np.messages = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  set.seed(20260717)
  n <- 40L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- x$x1^2 + x$x2 + rnorm(n, sd = 0.04)

  bw.lsq <- np_partial_gradient_local(nplsqregbw(
    xdat = x, ydat = y, bws = c(0.4, 0.4), bandwidth.compute = FALSE,
    scale = rep(0.4, n), delta = 0.5, tau = c(0.25, 0.75),
    regtype = "lp", degree = c(2L, 1L), basis = "glp"
  ))
  lsq <- np_partial_gradient_local(suppressWarnings(nplsqreg(
    bw.lsq, gradients = TRUE, gradient.order = c(2L, 2L)
  )))
  lsq.plot <- np_partial_gradient_local(suppressWarnings(plot(
    lsq, gradients = TRUE, gradient_order = c(2L, 2L),
    errors = "none", output = "data", neval = 5L,
    perspective = FALSE
  )))
  expect_true(any(is.finite(lsq.plot[[1L]]$quantgrad[, 1L, ])))
  expect_true(all(is.na(lsq.plot[[2L]]$quantgrad[, 2L, ])))

  bw.q <- np_partial_gradient_local(npcdistbw(
    xdat = x, ydat = data.frame(y = y), bws = c(0.45, 0.45, 0.45),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(0L, 2L),
    basis = "glp"
  ))
  qfit <- np_partial_gradient_local(suppressWarnings(npqreg(
    bw.q, txdat = x, tydat = y, tau = c(0.25, 0.75), gradients = TRUE
  )))
  qplot <- np_partial_gradient_local(suppressWarnings(plot(
    qfit, gradients = TRUE, errors = "asymptotic", output = "data",
    neval = 5L, perspective = FALSE
  )))
  expect_true(all(is.na(qplot[[1L]]$quantgrad[, 1L, ])))
  expect_true(all(is.na(qplot[[1L]]$quantgerr[, 1L, ])))
  expect_true(any(is.finite(qplot[[2L]]$quantgrad[, 2L, ])))
})

test_that("conditional and conmode plots retain categorical effects beside NA derivatives", {
  old <- options(np.messages = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  set.seed(20260718)
  n <- 42L
  x <- data.frame(
    x1 = runif(n, -1, 1),
    f = factor(rep(c("a", "b"), length.out = n)),
    x2 = runif(n, -1, 1)
  )
  y <- x$x1^2 + x$x2 + 0.35 * (x$f == "b") + rnorm(n, sd = 0.04)
  y.df <- data.frame(y = y)

  check.conditional <- function(bw) {
    out <- np_partial_gradient_local(suppressWarnings(plot(
      bw, xdat = x, ydat = y.df, gradients = TRUE,
      gradient.order = c(2L, 1L), errors = "none", output = "data",
      neval = 5L, perspective = FALSE
    )))
    for (slice in out) {
      expect_true(all(is.na(slice$congrad[, 1L])))
      expect_true(any(is.finite(slice$congrad[, 2L])))
      expect_true(any(is.finite(slice$congrad[, 3L])))
    }
  }
  check.conditional(np_partial_gradient_local(npcdensbw(
    xdat = x, ydat = y.df, bws = c(0.45, 0.5, 0.45, 0.45),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(1L, 2L),
    basis = "glp"
  )))
  check.conditional(np_partial_gradient_local(npcdistbw(
    xdat = x, ydat = y.df, bws = c(0.45, 0.5, 0.45, 0.45),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(1L, 2L),
    basis = "glp"
  )))

  yc <- factor(ifelse(y > median(y), "high", "low"),
               levels = c("low", "high"))
  bw.cm <- np_partial_gradient_local(npcdensbw(
    xdat = x, ydat = data.frame(y = yc),
    bws = c(0.05, 0.45, 0.45, 0.05), bandwidth.compute = FALSE,
    regtype = "lp", degree = c(0L, 2L), basis = "glp"
  ))
  cm <- np_partial_gradient_local(suppressWarnings(npconmode(
    bw.cm, txdat = x, tydat = yc, probabilities = TRUE,
    gradients = TRUE, gradient.level = "low"
  )))
  cm.plot <- np_partial_gradient_local(suppressWarnings(plot(
    cm, gradients = TRUE, view = "fixed", neval = 5L, output = "data"
  )))
  expect_true(all(is.na(cm.plot[[1L]]$effect)))
  expect_true(any(abs(cm.plot[[2L]]$effect) > 1e-10))
  expect_true(any(is.finite(cm.plot[[3L]]$effect)))
})

test_that("partial availability spans Bernstein, NN, and interleaved factors", {
  old <- options(np.messages = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  set.seed(20260720)
  n <- 40L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- 1 + x$x1 + x$x2^2 + rnorm(n, sd = 0.02)

  bw.bern <- np_partial_gradient_local(npregbw(
    xdat = x, ydat = y, regtype = "lp", degree = c(1L, 2L),
    basis = "tensor", bernstein.basis = TRUE, bwtype = "fixed",
    bws = c(0.42, 0.42), bandwidth.compute = FALSE
  ))
  bern <- np_partial_gradient_local(suppressWarnings(plot(
    bw.bern, xdat = x, ydat = y, gradients = TRUE,
    gradient.order = c(2L, 1L), errors = "none", output = "data",
    neval = 6L, perspective = FALSE
  )))
  expect_true(all(is.na(bern[[1L]]$grad)))
  expect_true(any(is.finite(bern[[2L]]$grad)))

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    bw.nn <- np_partial_gradient_local(npregbw(
      xdat = x, ydat = y, regtype = "lp", degree = c(1L, 2L),
      basis = "glp", bwtype = bwtype, bws = c(9, 9),
      bandwidth.compute = FALSE
    ))
    out.nn <- np_partial_gradient_local(suppressWarnings(plot(
      bw.nn, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = c(2L, 1L), errors = "none", output = "data",
      neval = 6L, perspective = FALSE
    )))
    expect_true(all(is.na(out.nn[[1L]]$grad)), info = bwtype)
    expect_true(any(is.finite(out.nn[[2L]]$grad)), info = bwtype)
  }

  mixed <- data.frame(
    x1 = runif(n, -1, 1),
    u = factor(rep(c("a", "b"), length.out = n)),
    o = ordered(rep(c("low", "mid", "high", "mid"), length.out = n),
                levels = c("low", "mid", "high")),
    x2 = runif(n, -1, 1)
  )
  y.mixed <- mixed$x1 + mixed$x2^2 +
    0.25 * (mixed$u == "b") + 0.15 * as.integer(mixed$o) +
    rnorm(n, sd = 0.02)
  bw.mixed <- np_partial_gradient_local(npregbw(
    xdat = mixed, ydat = y.mixed, regtype = "lp",
    degree = c(1L, 2L), basis = "glp", bwtype = "fixed",
    bws = c(0.45, 0.5, 0.5, 0.45), bandwidth.compute = FALSE
  ))
  mixed.out <- np_partial_gradient_local(suppressWarnings(plot(
    bw.mixed, xdat = mixed, ydat = y.mixed, gradients = TRUE,
    gradient.order = c(2L, 1L), errors = "none", output = "data",
    neval = 6L, perspective = FALSE
  )))
  expect_named(mixed.out, c("rg1", "rg2", "rg3", "rg4"))
  expect_true(all(is.na(mixed.out[[1L]]$grad)))
  expect_true(any(is.finite(mixed.out[[2L]]$grad)))
  expect_true(any(is.finite(mixed.out[[3L]]$grad)))
  expect_true(any(is.finite(mixed.out[[4L]]$grad)))

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plot.data <- np_partial_gradient_local(suppressWarnings(plot(
    bw.mixed, xdat = mixed, ydat = y.mixed, gradients = TRUE,
    gradient.order = c(2L, 1L), errors = "none", output = "plot-data",
    neval = 6L, perspective = FALSE, common.scale = FALSE,
    xlab = "custom x", ylab = "custom y", ylim = c(-3, 3)
  )))
  expect_named(plot.data, names(mixed.out))
  expect_true(all(is.na(plot.data[[1L]]$grad)))
})

test_that("common-scale all-band rendering retains unavailable terminal panels", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  if ("npRmpi" %in% loadedNamespaces())
    old <- c(old, options(npRmpi.autodispatch = FALSE))
  on.exit(options(old), add = TRUE)

  capture_warnings <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(
      force(expr),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(value = value, warnings = warnings)
  }

  expect_availability_only <- function(result, where,
                                       allow.simultaneous = FALSE) {
    availability <- grepl("x1.*requested order 2.*degree 1",
                          result$warnings)
    simultaneous <- grepl("asymptotic simultaneous confidence bands",
                          result$warnings)
    allowed <- availability |
      (isTRUE(allow.simultaneous) & simultaneous)
    expect_true(any(availability), info = where)
    expect_true(all(allowed), info = where)
    expect_true(!any(grepl("no non-missing arguments to (min|max)",
                           result$warnings)),
                info = where)
  }

  set.seed(20260721)
  n <- 60L
  x <- data.frame(
    group = factor(rep(c("a", "b"), length.out = n)),
    x2 = runif(n, -1, 1),
    x1 = runif(n, -1, 1)
  )
  y <- 0.4 * (x$group == "b") + x$x1 + x$x2^2 +
    rnorm(n, sd = 0.05)
  y.df <- data.frame(y = y)
  xc <- x[c("x2", "x1")]

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  bw.reg <- np_partial_gradient_local(npregbw(
    xdat = x, ydat = y, bws = c(0.5, 0.4, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    basis = "glp"
  ))
  if (!"npRmpi" %in% loadedNamespaces()) {
    reg.boot <- capture_warnings(np_partial_gradient_local(plot(
      bw.reg, xdat = x, ydat = y, gradients = TRUE,
      gradient.order = 2L, errors = "bootstrap", band = "all",
      common.scale = TRUE, B = 79L, neval = 2L
    )))
    expect_availability_only(reg.boot, "regression bootstrap common scale")
  }

  bw.reg.asym <- np_partial_gradient_local(npregbw(
    xdat = xc, ydat = y, bws = c(0.4, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    basis = "glp"
  ))
  reg.asym <- capture_warnings(np_partial_gradient_local(plot(
    bw.reg.asym, xdat = xc, ydat = y, gradients = TRUE,
    gradient.order = 2L, errors = "asymptotic", band = "all",
    common.scale = TRUE, neval = 2L
  )))
  expect_availability_only(reg.asym, "regression asymptotic common scale",
                           allow.simultaneous = TRUE)

  reg.panel <- capture_warnings(np_partial_gradient_local(plot(
    bw.reg.asym, xdat = xc, ydat = y, gradients = TRUE,
    gradient.order = 2L, errors = "asymptotic", band = "all",
    common.scale = FALSE, neval = 2L, ylim = c(-3, 3)
  )))
  expect_availability_only(reg.panel, "regression panel scale",
                           allow.simultaneous = TRUE)

  bw.cd <- np_partial_gradient_local(npcdensbw(
    xdat = xc, ydat = y.df, bws = c(0.4, 0.4, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    basis = "glp"
  ))
  cd <- capture_warnings(np_partial_gradient_local(plot(
    bw.cd, xdat = xc, ydat = y.df, gradients = TRUE,
    gradient.order = 2L, errors = "asymptotic", band = "all",
    common.scale = TRUE, neval = 2L, perspective = FALSE
  )))
  expect_availability_only(cd, "conditional density",
                           allow.simultaneous = TRUE)

  bw.dst <- np_partial_gradient_local(npcdistbw(
    xdat = xc, ydat = y.df, bws = c(0.4, 0.4, 0.4),
    bandwidth.compute = FALSE, regtype = "lp", degree = c(2L, 1L),
    basis = "glp"
  ))
  dst <- capture_warnings(np_partial_gradient_local(plot(
    bw.dst, xdat = xc, ydat = y.df, gradients = TRUE,
    gradient.order = 2L, errors = "asymptotic", band = "all",
    common.scale = TRUE, neval = 2L, perspective = FALSE
  )))
  expect_availability_only(dst, "conditional distribution",
                           allow.simultaneous = TRUE)
})
