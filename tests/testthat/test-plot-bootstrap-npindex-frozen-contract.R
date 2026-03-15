library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npindex nonfixed frozen bootstrap supports mean and gradient slices", {
  set.seed(42)
  n <- 80L
  x <- runif(n)
  z <- rnorm(n)
  y <- x + rnorm(n)

  fit <- quiet_eval(
    npindex(
      y ~ x + z,
      regtype = "ll",
      bwtype = "adaptive_nn",
      nmulti = 1
    )
  )

  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    try(grDevices::dev.off(), silent = TRUE)
    unlink(tf)
  }, add = TRUE)

  expect_no_error(
    capture.output(plot(
      fit,
      neval = 20L,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.nonfixed = "frozen",
      plot.errors.boot.num = 41L,
      plot.errors.type = "pointwise",
      gradients = TRUE
    ))
  )
})
