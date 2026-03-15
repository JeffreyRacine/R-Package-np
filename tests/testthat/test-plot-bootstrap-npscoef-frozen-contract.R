library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npscoef nonfixed frozen bootstrap supports generalized and adaptive routes", {
  run_case <- function(bwtype, boot.method) {
    set.seed(4201)
    n <- 80L
    x <- runif(n)
    z <- runif(n)
    y <- x * (1 + z) + rnorm(n, sd = 0.05)

    fit <- quiet_eval(
      npscoef(
        y ~ x | z,
        regtype = "ll",
        bwtype = bwtype,
        nmulti = 1,
        betas = TRUE
      )
    )

    tf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tf)
    on.exit({
      try(grDevices::dev.off(), silent = TRUE)
      unlink(tf)
    }, add = TRUE)

    args <- list(
      x = fit,
      neval = 20L,
      coef = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = boot.method,
      plot.errors.boot.nonfixed = "frozen",
      plot.errors.boot.num = 41L,
      plot.errors.type = "pointwise"
    )
    if (identical(boot.method, "geom"))
      args$plot.errors.boot.blocklen <- 4L

    expect_no_error(capture.output(do.call(plot, args)))
  }

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    run_case(bwtype = bwtype, boot.method = "inid")
    run_case(bwtype = bwtype, boot.method = "geom")
  }
})
