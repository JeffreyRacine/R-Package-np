test_that("conditional density/distribution gradient bootstrap inid fails fast for bandwidth objects", {
  skip_if_not_installed("np")

  set.seed(20260225)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.45, 0.45, 0.45),
    bandwidth.compute = FALSE
  )
  bw.cdist <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.45, 0.45, 0.45),
    bandwidth.compute = FALSE
  )

  expect_error(
    suppressWarnings(
      plot(
        bw.cd,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 5
      )
    ),
    "inid conditional helper unavailable",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        bw.cdist,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 5
      )
    ),
    "inid conditional helper unavailable",
    fixed = TRUE
  )
})

test_that("conditional density/distribution gradient bootstrap inid fails fast for fitted objects", {
  skip_if_not_installed("np")

  set.seed(20260226)
  n <- 35L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.5, 0.5, 0.5),
    bandwidth.compute = FALSE
  )
  bw.cdist <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.5, 0.5, 0.5),
    bandwidth.compute = FALSE
  )

  fit.cd <- npcdens(bws = bw.cd)
  fit.cdist <- npcdist(bws = bw.cdist)

  expect_error(
    suppressWarnings(
      plot(
        fit.cd,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 5
      )
    ),
    "inid conditional helper unavailable",
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        fit.cdist,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 5
      )
    ),
    "inid conditional helper unavailable",
    fixed = TRUE
  )
})
