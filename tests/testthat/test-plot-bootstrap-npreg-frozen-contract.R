library(np)

quiet_eval <- function(expr) {
  value <- NULL
  capture.output(value <- force(expr))
  value
}

test_that("npreg nonfixed frozen plot bootstrap supports continuous gradients", {
  run_case <- function(bwtype, boot.method) {
    set.seed(42)
    n <- 80L
    x <- runif(n)
    y <- x + rnorm(n)

    g <- quiet_eval(
      npreg(
        y ~ x,
        regtype = "ll",
        degree = c(1),
        bernstein = TRUE,
        bwtype = bwtype,
        nmulti = 1
      )
    )

    tf <- tempfile(fileext = ".pdf")
    grDevices::pdf(tf)
    on.exit({
      grDevices::dev.off()
      unlink(tf)
    }, add = TRUE)

    args <- list(
      x = g,
      neval = 20L,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = boot.method,
      plot.errors.boot.nonfixed = "frozen",
      plot.errors.boot.num = 41L,
      plot.errors.type = "pointwise",
      gradients = TRUE
    )
    if (identical(boot.method, "geom"))
      args$plot.errors.boot.blocklen <- 4L

    expect_no_error(do.call(plot, args))
  }

  for (bwtype in c("generalized_nn", "adaptive_nn")) {
    run_case(bwtype = bwtype, boot.method = "inid")
    run_case(bwtype = bwtype, boot.method = "geom")
  }
})
