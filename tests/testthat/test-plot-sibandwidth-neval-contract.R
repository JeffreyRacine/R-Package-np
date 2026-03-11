test_that("sibandwidth plot respects bounded neval in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  n <- 32
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(tx$x1 + 2 * tx$x2) + rnorm(n, sd = 0.05)
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.25),
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )

  out.mean <- suppressWarnings(plot(
    bw,
    xdat = tx,
    ydat = y,
    neval = 13L,
    plot.behavior = "data",
    perspective = FALSE,
    gradients = FALSE
  ))[[1]]

  expect_equal(length(out.mean$index), 13L)
  expect_equal(length(out.mean$mean), 13L)
  expect_false(out.mean$trainiseval)

  out.grad <- suppressWarnings(plot(
    bw,
    xdat = tx,
    ydat = y,
    neval = 13L,
    plot.behavior = "data",
    perspective = FALSE,
    gradients = TRUE
  ))[[1]]

  expect_equal(length(out.grad$index), 13L)
  expect_equal(nrow(out.grad$grad), 13L)
  expect_equal(dim(out.grad$glerr), c(13L, ncol(tx)))
  expect_equal(dim(out.grad$gherr), c(13L, ncol(tx)))
  expect_false(out.grad$trainiseval)
  expect_true(all(is.finite(out.grad$grad)))
})

test_that("sibandwidth bootstrap helpers honor bounded neval in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  n <- 30
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- cos(tx$x1 - tx$x2) + rnorm(n, sd = 0.06)

  for (bt in c("fixed", "adaptive_nn")) {
    h <- if (identical(bt, "fixed")) 0.25 else 5L
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, h),
      bandwidth.compute = FALSE,
      regtype = "ll",
      bwtype = bt
    )

    for (boot.method in c("wild", "inid", "fixed", "geom")) {
      out <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = 11L,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = boot.method,
        plot.errors.boot.num = 9L
      ))[[1]]

      expect_equal(length(out$mean), 11L, info = paste(bt, boot.method, "mean length"))
      expect_equal(dim(out$merr), c(11L, 2L), info = paste(bt, boot.method, "merr shape"))
      expect_false(out$trainiseval, info = paste(bt, boot.method, "trainiseval"))
      expect_true(all(is.finite(out$mean)), info = paste(bt, boot.method, "mean finite"))
      expect_true(all(is.finite(out$merr)), info = paste(bt, boot.method, "merr finite"))
    }
  }
})
