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
    output = "data",
    perspective = FALSE,
    gradients = FALSE
  ))

  expect_s3_class(out.mean, "np_plot_data")
  expect_equal(nrow(out.mean$panels[[1L]]$grid), 13L)
  expect_equal(nrow(plot_target(out.mean)), 13L)
  expect_false(nrow(out.mean$panels[[1L]]$grid) == n)

  out.grad <- suppressWarnings(plot(
    bw,
    xdat = tx,
    ydat = y,
    neval = 13L,
    output = "data",
    perspective = FALSE,
    gradients = TRUE
  ))

  expect_equal(nrow(out.grad$panels[[1L]]$grid), 13L)
  expect_equal(nrow(out.grad$targets), ncol(tx))
  expect_setequal(out.grad$targets$variable, names(tx))
  for (target in seq_len(ncol(tx))) {
    d <- plot_target(out.grad, target = target)
    expect_equal(nrow(d), 13L)
    expect_true(all(is.finite(d$estimate)))
    expect_false(any(c("lower", "upper") %in% names(d)))
    expect_identical(out.grad$request$errors, "none")
    expect_equal(nrow(out.grad$panels[[1L]]$intervals), 0L)
  }
})

test_that("sibandwidth bootstrap helpers honor bounded neval in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260311)
  n <- 30
  tx <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- cos(tx$x1 - tx$x2) + rnorm(n, sd = 0.06)

  for (bt in c("fixed", "adaptive_nn")) {
    # Bootstrap rows are discrete: keep this success fixture inside the
    # literal ANN-radius domain; zero-radius errors have dedicated contracts.
    h <- if (identical(bt, "fixed")) 0.25 else n - 1L
    bw <- do.call(npindexbw, list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, h),
      bandwidth.compute = FALSE,
      regtype = "ll",
      bwtype = bt
    ))

    for (boot.method in c("wild", "inid", "fixed", "geom")) {
      out <- suppressWarnings(plot(
        bw,
        xdat = tx,
        ydat = y,
        neval = 11L,
        output = "data",
        perspective = FALSE,
        gradients = FALSE,
        errors = "bootstrap",
        bootstrap = boot.method,
        B = 9L
      ))

      d <- plot_target(out, band = "pmzsd")
      expect_equal(nrow(d), 11L, info = paste(bt, boot.method, "mean length"))
      expect_equal(dim(as.matrix(d[c("lower", "upper")])), c(11L, 2L))
      expect_false(nrow(out$panels[[1L]]$grid) == n)
      expect_true(all(is.finite(d$estimate)), info = paste(bt, boot.method, "mean finite"))
      expect_true(all(is.finite(as.matrix(d[c("lower", "upper")]))),
                  info = paste(bt, boot.method, "endpoints finite"))
    }
  }
})
