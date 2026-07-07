skip_slow_npreg_adaptive_lc_exdat_search <- function() {
  skip_if_not(
    identical(Sys.getenv("NP_RUN_SLOW_NPREG_ADAPTIVE_LC_EXDAT_MPI"), "true"),
    paste(
      "set NP_RUN_SLOW_NPREG_ADAPTIVE_LC_EXDAT_MPI=true to run the slow",
      "computed adaptive-nn lc exdat bandwidth-search smoke"
    )
  )
}

check_adaptive_lc_exdat_mean <- function(bw, tx, y, ex) {
  fit.ex <- npreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    exdat = ex,
    gradients = FALSE,
    warn.glp.gradient = FALSE
  )
  hat.apply <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = y,
    output = "apply"
  )
  hat.matrix <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    output = "matrix"
  )

  expect_equal(as.vector(hat.apply), as.vector(fit.ex$mean), tolerance = 1e-8)
  expect_equal(as.vector(hat.matrix %*% y), as.vector(fit.ex$mean), tolerance = 1e-8)
}

test_that("adaptive-nn lc manual bandwidth exdat mean matches npreg and npreghat on boundary grid", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260315)
  n <- 24L
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 8L))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    bws = 5,
    bandwidth.compute = FALSE
  )

  check_adaptive_lc_exdat_mean(bw, tx, y, ex)
})

test_that("adaptive-nn lc computed bandwidth exdat mean matches npreg and npreghat on boundary grid", {
  skip_slow_npreg_adaptive_lc_exdat_search()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260315)
  n <- 100L
  x <- runif(n, -1, 1)
  y <- x^2 + rnorm(n, sd = 0.25 * stats::sd(x))
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 25L))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lc",
    bwtype = "adaptive_nn",
    nmulti = 3L
  )

  check_adaptive_lc_exdat_mean(bw, tx, y, ex)
})
