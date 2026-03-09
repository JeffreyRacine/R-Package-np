with_public_trace_counter <- function(pkg, fun, expr) {
  ns <- asNamespace(pkg)
  ctr <- new.env(parent = emptyenv())
  ctr$n <- 0L

  trace(
    what = fun,
    where = ns,
    tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
    print = FALSE
  )
  on.exit(try(untrace(fun, where = ns), silent = TRUE), add = TRUE)

  force(expr)
  ctr$n
}

test_that("object-fed regression plots avoid public npreg re-entry", {
  set.seed(9210)
  n <- 48
  xdat <- data.frame(x = runif(n))
  ydat <- sin(2 * pi * xdat$x) + rnorm(n, sd = 0.05)

  bw <- npregbw(xdat = xdat, ydat = ydat, regtype = "lc", nmulti = 1)
  fit <- npreg(txdat = xdat, tydat = ydat, bws = bw)

  calls <- with_public_trace_counter("np", "npreg", {
    out.bw <- plot(bw, xdat = xdat, ydat = ydat, plot.behavior = "data")
    out.fit <- plot(fit, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)
})

test_that("object-fed unconditional density/distribution plots avoid public re-entry", {
  set.seed(9211)
  n <- 42
  dat <- data.frame(x = runif(n))

  dens.bw <- npudensbw(dat = dat, bws = 0.18, bandwidth.compute = FALSE, bwtype = "fixed")
  dist.bw <- npudistbw(dat = dat, bws = 0.18, bandwidth.compute = FALSE, bwtype = "fixed")
  dens.fit <- npudens(tdat = dat, bws = dens.bw)
  dist.fit <- npudist(tdat = dat, bws = dist.bw)

  dens.calls <- with_public_trace_counter("np", "npudens", {
    out.bw <- plot(dens.bw, xdat = dat, plot.behavior = "data")
    out.fit <- plot(dens.fit, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  dist.calls <- with_public_trace_counter("np", "npudist", {
    out.bw <- plot(dist.bw, xdat = dat, plot.behavior = "data")
    out.fit <- plot(dist.fit, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(dens.calls, 0L)
  expect_identical(dist.calls, 0L)
})

test_that("object-fed partially linear plots avoid public npplreg re-entry", {
  set.seed(9212)
  n <- 54
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- 1.5 * xdat$x + cos(2 * pi * zdat$z) + rnorm(n, sd = 0.05)

  bw <- npplregbw(xdat = xdat, ydat = ydat, zdat = zdat, regtype = "lc", nmulti = 1)
  fit <- npplreg(txdat = xdat, tydat = ydat, tzdat = zdat, bws = bw)

  calls <- with_public_trace_counter("np", "npplreg", {
    out.bw <- plot(bw, xdat = xdat, ydat = ydat, zdat = zdat, plot.behavior = "data")
    out.fit <- plot(fit, plot.behavior = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)
})
