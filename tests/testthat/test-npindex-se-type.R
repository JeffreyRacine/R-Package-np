test_that("single-index uncertainty choices preserve lean refitting", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(73)
  x <- data.frame(x1 = runif(48, -1, 1), x2 = runif(48, -1, 1))
  y <- sin(x$x1 + .4*x$x2) + .1*cos(seq_len(48))
  dat <- cbind(x, y = y)
  b <- npindexbw(y ~ x1 + x2, data = dat, bws = c(1, .4, .8),
                 bandwidth.compute = FALSE)
  initial <- npindex(bws = b, se = FALSE)
  expect_error(npindex(bws = b, se.type = "unknown"), "arg")
  testthat::local_mocked_bindings(
    npindexbw = function(...) stop("unexpected repeated bandwidth selection"),
    .package = "npRmpi"
  )
  default <- npindex(bws = initial$bws)
  expect_true(all(is.finite(vcov(default))))
  expect_true(all(is.finite(se(default))))
  expect_false(default$gradients)
  with.grad <- npindex(bws = initial$bws, gradients = TRUE)
  expect_identical(vcov(default), vcov(with.grad))
  # npreg's asymptotic owner and the lean kernel-sum owner can differ by ulps.
  expect_equal(fitted(default), fitted(initial), tolerance = 1e-14)
  expect_error(gradients(initial), "initial$bws", fixed = TRUE)
  expect_error(se(initial), "initial$bws", fixed = TRUE)
  expect_error(vcov(identity(initial)), "Replace 'object'", fixed = TRUE)
  rng <- .Random.seed
  light <- npindex(bws = b, se = FALSE, gradients = TRUE)
  expect_identical(.Random.seed, rng)
  expect_null(light$merr)
  expect_null(light$gerr)
  expect_null(light$betavcov)
  expect_true(all(is.finite(gradients(light))))
  boot <- npindex(bws = b, se.type = "bootstrap", gradients = TRUE, B = 3L)
  expect_identical(boot$se.type, "bootstrap")
  expect_true(all(is.finite(boot$mean.gerr)))
  rng <- .Random.seed
  pred.asym <- predict(boot, newdata = x[1:5, ], se.fit = TRUE,
                       se.type = "asymptotic")
  expect_identical(.Random.seed, rng)
  expect_true(all(is.finite(pred.asym$se.fit)))
  pred.boot <- predict(boot, newdata = x[1:5, ], se.fit = TRUE, B = 3L)
  expect_false(identical(.Random.seed, rng))
  expect_equal(pred.boot$fit, pred.asym$fit, tolerance = 1e-14)
})

test_that("adaptive-NN mean and bootstrap use the same donor normalization", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(619L)
  x <- data.frame(x1 = runif(48, -1, 1), x2 = runif(48, -1, 1),
                  x3 = runif(48, -1, 1))
  index <- x$x1 + .6*x$x2 - .25*x$x3
  y <- sin(index) + .1*cos(seq_len(48))
  b <- npindexbw(xdat = x, ydat = y, bws = c(1, .6, -.25, 32),
                 bwtype = "adaptive_nn", bandwidth.compute = FALSE)
  mean.only <- npindex(bws = b, txdat = x, tydat = y, se = FALSE)
  reference <- npreg(txdat = data.frame(index = index), tydat = y,
                      bws = 32, bwtype = "adaptive_nn", se = FALSE)
  expect_equal(as.vector(fitted(mean.only)), as.vector(fitted(reference)), tolerance = 1e-10)
  H <- npindexhat(bws = b, txdat = x, exdat = x, output = "matrix")
  expect_equal(as.vector(H %*% y), as.vector(fitted(reference)), tolerance = 1e-10)
  expect_equal(as.vector(rowSums(H)), rep(1, nrow(x)), tolerance = 1e-14)
  plot.mean <- .np_plot_singleindex_hat_apply_index(b, data.frame(index = index),
                                                    data.frame(index = index), y)
  expect_equal(as.vector(plot.mean), as.vector(fitted(reference)), tolerance = 1e-10)
  requested.type <- "bootstrap"
  set.seed(42)
  boot.mean <- npindex(bws = b, txdat = x, tydat = y, se.type = requested.type, B = 3L)
  rng <- .Random.seed
  set.seed(42)
  boot.grad <- npindex(bws = b, txdat = x, tydat = y, se.type = "bootstrap",
                       gradients = TRUE, B = 3L)
  expect_identical(.Random.seed, rng)
  expect_equal(as.vector(fitted(boot.mean)), as.vector(fitted(boot.grad)), tolerance = 1e-10)
  expect_equal(as.vector(se(boot.mean)), as.vector(se(boot.grad)), tolerance = 1e-10)
})

test_that("one normalized index coefficient needs no information inversion", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = seq(-1, 1, length.out = 32))
  y <- sin(x$x1) + .1*cos(seq_len(32))
  expect_warning(b <- npindexbw(xdat = x, ydat = y, bws = c(1, .8),
                                bandwidth.compute = FALSE),
                 "xdat has one dimension", fixed = TRUE)
  testthat::local_mocked_bindings(
    .np_index_covariance_inverse = function(...) stop("unrequested inverse"),
    .package = "npRmpi"
  )
  for (grad in c(FALSE, TRUE)) {
    fit <- npindex(bws = b, txdat = x, tydat = y, gradients = grad)
    expect_identical(unname(vcov(fit)), matrix(0, 1, 1))
    expect_true(all(is.finite(se(fit))))
  }
  x$x2 <- cos(seq_len(32))
  b <- npindexbw(xdat = x, ydat = y, bws = c(1, .2, .8),
                 bandwidth.compute = FALSE)
  fit <- npindex(bws = b, txdat = x, tydat = y, se = FALSE, gradients = TRUE)
  expect_null(fit$betavcov)
})
