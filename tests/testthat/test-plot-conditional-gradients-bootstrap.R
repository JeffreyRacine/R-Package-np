# Numerical-provider proof: public plot schema is checked by the interception
# helper and dedicated schema/geometry tests; no provider fields are public.
test_that("conditional density/distribution gradient bootstrap inid works for bandwidth objects", {
  skip_if_not_installed("np")

  library(np)

  set.seed(20260225)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  bw.cdist <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))

  for (obj in list(bw.cd, bw.cdist)) {
    out <- suppressWarnings(
      capture_plot_provider(plot(
        obj,
        xdat = xdat,
        ydat = ydat,
        output = "data",
        perspective = FALSE,
        gradients = TRUE,
        errors = "bootstrap",
        bootstrap = "inid",
        B = 5L
      ))$provider
    )
    expect_type(out, "list")
    expect_true(length(out[[1L]]$bxp) > 0L)
  }
})

test_that("conditional density/distribution gradient bootstrap inid works for fitted objects", {
  skip_if_not_installed("np")

  library(np)

  set.seed(20260226)
  n <- 24L
  x1 <- factor(sample(c("a", "b"), n, replace = TRUE))
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- suppressWarnings(npcdensbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  bw.cdist <- suppressWarnings(npcdistbw(xdat = xdat, ydat = ydat, nmulti = 1L))
  fit.cd <- npcdens(bws = bw.cd)
  fit.cdist <- npcdist(bws = bw.cdist)

  for (obj in list(fit.cd, fit.cdist)) {
    out <- suppressWarnings(
      capture_plot_provider(plot(
        obj,
        output = "data",
        perspective = FALSE,
        gradients = TRUE,
        errors = "bootstrap",
        bootstrap = "inid",
        B = 5L
      ))$provider
    )
    expect_type(out, "list")
    expect_true(length(out[[1L]]$bxp) > 0L)
  }
})
