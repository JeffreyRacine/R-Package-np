test_that("lp bernstein bootstrap plot errors fail fast", {
  suppressPackageStartupMessages(library(npRmpi))
  npRmpi.init(nslaves = 1)
  on.exit(try(npRmpi.quit(force = TRUE), silent = TRUE), add = TRUE)

  set.seed(4201)
  n <- 120
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)

  fit <- npreg(y ~ x, regtype = "lp", degree = 3L, bernstein = TRUE)

  expect_error(
    suppressWarnings(
      plot(
        fit,
        plot.behavior = "data",
        neval = 40,
        gradients = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 9L,
        plot.errors.type = "pointwise"
      )
    ),
    "bootstrap plot errors for regtype='lp' with bernstein.basis=TRUE are unsupported",
    fixed = TRUE
  )
})
