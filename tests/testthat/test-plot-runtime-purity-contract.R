library(npRmpi)

with_public_trace_counter <- function(fun, expr) {
  ns <- asNamespace("npRmpi")
  ctr <- new.env(parent = emptyenv())
  ctr$n <- 0L

  trace(
    what = fun,
    where = ns,
    tracer = bquote(assign("n", get("n", envir = .(ctr)) + 1L, envir = .(ctr))),
    print = FALSE
  )
  on.exit(try(untrace(fun, where = ns), silent = TRUE), add = TRUE)

  force(expr)
  ctr$n
}

test_that("session-route object-fed plot paths avoid public estimator re-entry", {
  npRmpi.init(nslaves = 1, quiet = TRUE)
  on.exit(npRmpi.quit(), add = TRUE)

  set.seed(9410)
  n.reg <- 48
  xdat.reg <- data.frame(x = runif(n.reg))
  ydat.reg <- sin(2 * pi * xdat.reg$x) + rnorm(n.reg, sd = 0.05)

  bw.reg <- npregbw(xdat = xdat.reg, ydat = ydat.reg, regtype = "lc", nmulti = 1)
  fit.reg <- npreg(txdat = xdat.reg, tydat = ydat.reg, bws = bw.reg)

  calls <- with_public_trace_counter("npreg", {
    out.bw <- plot(bw.reg, xdat = xdat.reg, ydat = ydat.reg, plot.behavior = "data")
    out.fit <- plot(fit.reg, xdat = xdat.reg, ydat = ydat.reg, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)

  set.seed(9411)
  n <- 42
  dat <- data.frame(x = runif(n))

  dens.bw <- npudensbw(dat = dat, bws = 0.18, bandwidth.compute = FALSE, bwtype = "fixed")
  dist.bw <- npudistbw(dat = dat, bws = 0.18, bandwidth.compute = FALSE, bwtype = "fixed")
  dens.fit <- npudens(tdat = dat, bws = dens.bw)
  dist.fit <- npudist(tdat = dat, bws = dist.bw)

  dens.calls <- with_public_trace_counter("npudens", {
    out.bw <- plot(dens.bw, xdat = dat, plot.behavior = "data")
    out.fit <- plot(dens.fit, xdat = dat, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  dist.calls <- with_public_trace_counter("npudist", {
    out.bw <- plot(dist.bw, xdat = dat, plot.behavior = "data")
    out.fit <- plot(dist.fit, xdat = dat, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(dens.calls, 0L)
  expect_identical(dist.calls, 0L)

  set.seed(9412)
  n <- 20
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- 1.5 * xdat$x + cos(2 * pi * zdat$z) + rnorm(n, sd = 0.05)

  bw <- npplregbw(xdat = xdat, ydat = ydat, zdat = zdat, regtype = "lc", nmulti = 1)
  fit <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat, bws = bw)

  calls <- with_public_trace_counter("npplreg", {
    out.bw <- plot(
      bw,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      coef = TRUE,
      plot.behavior = "data",
      perspective = FALSE
    )
    out.fit <- plot(
      fit,
      xdat = xdat,
      ydat = ydat,
      zdat = zdat,
      coef = TRUE,
      plot.behavior = "data",
      perspective = FALSE
    )
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)

  set.seed(9413)
  n <- 52
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  ydat <- sin(xdat$x1 + xdat$x2) + rnorm(n, sd = 0.05)

  bw <- npindexbw(xdat = xdat, ydat = ydat, method = "ichimura", nmulti = 1)
  fit <- npindex(txdat = xdat, tydat = ydat, bws = bw)

  calls <- with_public_trace_counter("npindex", {
    out.bw <- plot(bw, xdat = xdat, ydat = ydat, plot.behavior = "data")
    out.fit <- plot(fit, xdat = xdat, ydat = ydat, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)
})
