test_that("npindex plot bootstrap inid supports nearest-neighbor bwtypes", {
  skip_if_not_installed("np")

  set.seed(3291)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  run_plot <- function(bw) {
    suppressWarnings(plot(
      bw,
      xdat = tx,
      ydat = y,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.lc <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt,
      regtype = "lc"
    )
    out.lc <- run_plot(bw.lc)
    expect_type(out.lc, "list")
    expect_true(length(out.lc) > 0)

    bw.ll <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt,
      regtype = "ll"
    )
    out.ll <- run_plot(bw.ll)
    expect_type(out.ll, "list")
    expect_true(length(out.ll) > 0)
  }
})
