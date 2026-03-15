test_that("npcdens plot-data centers stay on predict and bootstrap offsets remain ordered", {
  if (!spawn_mpi_slaves(1L)) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  check_mode <- function(bwtype, nonfixed_mode) {
    set.seed(42)
    n <- 100L
    x <- runif(n, -1, 1)
    y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))

    bw <- npcdensbw(y ~ x, nmulti = 1L, bwtype = bwtype)
    fit <- npcdens(bws = bw)

    out <- suppressWarnings(plot(
      fit,
      plot.behavior = "data",
      view = "fixed",
      perspective = FALSE,
      neval = 40L,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 39L,
      plot.errors.type = "pointwise",
      plot.errors.boot.nonfixed = nonfixed_mode
    ))

    expect_true(length(out) > 0L)
    for (obj in out) {
      expect_s3_class(obj, "condensity")
      pred <- predict(fit, exdat = obj$xeval, eydat = obj$yeval)
      lower <- obj$condens + obj$conderr[, 1L]
      upper <- obj$condens + obj$conderr[, 2L]

      expect_equal(obj$condens, as.numeric(pred), tolerance = 1e-8)
      expect_true(all(lower <= upper + 1e-10, na.rm = TRUE))
    }

    out
  }

  fixed.exact <- check_mode("fixed", "exact")
  generalized.exact <- check_mode("generalized_nn", "exact")
  generalized.frozen <- check_mode("generalized_nn", "frozen")
  adaptive.exact <- check_mode("adaptive_nn", "exact")
  adaptive.frozen <- check_mode("adaptive_nn", "frozen")

  expect_equal(
    lapply(generalized.exact, `[[`, "condens"),
    lapply(generalized.frozen, `[[`, "condens"),
    tolerance = 1e-12
  )
  expect_equal(
    lapply(adaptive.exact, `[[`, "condens"),
    lapply(adaptive.frozen, `[[`, "condens"),
    tolerance = 1e-12
  )
  expect_true(all(vapply(fixed.exact, function(obj) is.matrix(obj$conderr), logical(1))))
})
