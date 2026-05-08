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
    out.bw <- plot(bw, xdat = xdat, ydat = ydat, output = "data")
    out.fit <- plot(fit, output = "data")
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
    out.bw <- plot(dens.bw, xdat = dat, output = "data")
    out.fit <- plot(dens.fit, output = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  dist.calls <- with_public_trace_counter("np", "npudist", {
    out.bw <- plot(dist.bw, xdat = dat, output = "data")
    out.fit <- plot(dist.fit, output = "data")
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
    out.bw <- plot(bw, xdat = xdat, ydat = ydat, zdat = zdat, output = "data")
    out.fit <- plot(fit, output = "data")
    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
  })

  expect_identical(calls, 0L)
})

test_that("object-fed conditional density/distribution plots avoid public re-entry", {
  set.seed(9213)
  n <- 46
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = rnorm(n))

  dens.bw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.22, 0.22), bandwidth.compute = FALSE)
  dist.bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.22, 0.22), bandwidth.compute = FALSE)
  dens.fit <- npcdens(txdat = xdat, tydat = ydat, bws = dens.bw)
  dist.fit <- npcdist(txdat = xdat, tydat = ydat, bws = dist.bw)

  dens.calls <- with_public_trace_counter("np", "npcdens", {
    out.grid.bw <- plot(dens.bw, xdat = xdat, ydat = ydat, output = "data", view = "fixed")
    out.grid.fit <- plot(dens.fit, output = "data", view = "fixed")
    out.slice.bw <- plot(dens.bw, xdat = xdat, ydat = ydat, output = "data", perspective = FALSE, errors = "asymptotic")
    out.slice.fit <- plot(dens.fit, output = "data", perspective = FALSE, errors = "asymptotic")
    expect_type(out.grid.bw, "list")
    expect_type(out.grid.fit, "list")
    expect_type(out.slice.bw, "list")
    expect_type(out.slice.fit, "list")
  })

  dist.calls <- with_public_trace_counter("np", "npcdist", {
    out.grid.bw <- plot(dist.bw, xdat = xdat, ydat = ydat, output = "data", view = "fixed")
    out.grid.fit <- plot(dist.fit, output = "data", view = "fixed")
    out.slice.bw <- plot(dist.bw, xdat = xdat, ydat = ydat, output = "data", perspective = FALSE, errors = "asymptotic")
    out.slice.fit <- plot(dist.fit, output = "data", perspective = FALSE, errors = "asymptotic")
    expect_type(out.grid.bw, "list")
    expect_type(out.grid.fit, "list")
    expect_type(out.slice.bw, "list")
    expect_type(out.slice.fit, "list")
  })

  expect_identical(dens.calls, 0L)
  expect_identical(dist.calls, 0L)
})

test_that("object-fed quantile plots avoid public npqreg re-entry", {
  set.seed(9214)
  n <- 46
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = rnorm(n))

  bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.22, 0.22), bandwidth.compute = FALSE)
  fit <- npqreg(txdat = xdat, tydat = ydat, bws = bw, tau = 0.4)

  calls <- with_public_trace_counter("np", "npqreg", {
    out.grid.bw <- plot(bw, xdat = xdat, ydat = ydat, output = "data", view = "fixed", quantreg = TRUE, tau = 0.4)
    out.grid.fit <- plot(fit, output = "data", view = "fixed")
    out.slice.bw <- plot(bw, xdat = xdat, ydat = ydat, output = "data", perspective = FALSE, quantreg = TRUE, tau = 0.4)
    out.slice.fit <- plot(fit, output = "data", perspective = FALSE)
    expect_type(out.grid.bw, "list")
    expect_type(out.grid.fit, "list")
    expect_type(out.slice.bw, "list")
    expect_type(out.slice.fit, "list")
  })

  expect_identical(calls, 0L)
})

test_that("smooth coefficient coef/asymptotic plots avoid public npscoef re-entry", {
  set.seed(9215)
  n <- 52
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  ydat <- sin(2 * pi * zdat$z) + xdat$x * (1 + zdat$z) + rnorm(n, sd = 0.05)

  bw <- npscoefbw(xdat = xdat, ydat = ydat, zdat = zdat, regtype = "ll", nmulti = 1)
  fit <- npscoef(bws = bw, txdat = xdat, tydat = ydat, tzdat = zdat, errors = TRUE, betas = TRUE)

  calls <- with_public_trace_counter("np", "npscoef", {
    out.bw.coef <- suppressWarnings(
      plot(
        bw,
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        coef = TRUE,
        perspective = FALSE,
        output = "data",
        errors = "none"
      )
    )
    out.fit.coef <- suppressWarnings(
      plot(
        fit,
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        coef = TRUE,
        perspective = FALSE,
        output = "data",
        errors = "none"
      )
    )
    out.bw.asym <- suppressWarnings(
      plot(
        bw,
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        coef = FALSE,
        perspective = FALSE,
        output = "data",
        errors = "asymptotic"
      )
    )
    out.fit.asym <- suppressWarnings(
      plot(
        fit,
        xdat = xdat,
        ydat = ydat,
        zdat = zdat,
        coef = FALSE,
        perspective = FALSE,
        output = "data",
        errors = "asymptotic"
      )
    )
    expect_type(out.bw.coef, "list")
    expect_type(out.fit.coef, "list")
    expect_type(out.bw.asym, "list")
    expect_type(out.fit.asym, "list")
  })

  expect_identical(calls, 0L)
})
