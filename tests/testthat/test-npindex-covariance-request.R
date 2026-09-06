test_that("single-index coefficient covariance does not request bootstrap SDs", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260906)
  n <- 40L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
                  x3 = runif(n, -1, 1))
  index <- x$x1 + .6 * x$x2 - .25 * x$x3
  for (method in c("ichimura", "kleinspady")) {
    y <- if (method == "ichimura") sin(index) + .1 * cos(seq_len(n)) else
      as.double(index + sin(seq_len(n) * 2.3) > 0)
    bw <- npindexbw(xdat = x, ydat = y, method = method,
                    bws = c(1, .6, -.25, .8), bandwidth.compute = FALSE)
    for (external in c(FALSE, TRUE)) {
      args <- list(bws = bw, txdat = x, tydat = y, gradients = TRUE, B = 3L)
      if (external) args$exdat <- x[seq_len(7L), , drop = FALSE]
      rng <- .Random.seed
      fit <- do.call(npindex, args)
      expect_identical(.Random.seed, rng)
      expect_true(all(is.finite(vcov(fit))))
      expect_identical(unname(vcov(fit)[1L, ]), rep(0, 3L))
      expect_identical(unname(vcov(fit)[, 1L]), rep(0, 3L))
      expect_true(all(is.finite(se(fit))))
      expect_true(all(is.finite(gradients(fit, se = TRUE))))
      expect_identical(fit$se.type, "asymptotic")
      expect_null(fit$mean.gerr)
      full <- do.call(npindex, c(args, list(se = TRUE, se.type = "bootstrap")))
      expect_identical(vcov(fit), vcov(full))
      expect_identical(fitted(fit), fitted(full))
      expect_identical(gradients(fit), gradients(full))
      expect_true(all(is.finite(se(full))))
      expect_true(all(is.finite(gradients(full, se = TRUE))))
    }
    plain <- npindex(bws = bw, txdat = x, tydat = y, se = FALSE)
    expect_null(plain$betavcov)
    expect_error(vcov(plain), "npindex(bws = plain$bws, gradients = FALSE, se = TRUE)", fixed = TRUE)
    expect_error(gradients(plain), "npindex(bws = plain$bws, gradients = TRUE, se = FALSE)", fixed = TRUE)
    default <- npindex(bws = bw, txdat = x, tydat = y)
    expect_true(all(is.finite(vcov(default))))
    expect_identical(default$gradients, FALSE)
    expect_true(all(is.finite(se(default))))
  }
})

test_that("single-index refit hints never evaluate object expressions", {
  expect_match(.np_index_refit_hint(quote(my_fit), se = TRUE), "my_fit$bws", fixed = TRUE)
  expect_match(.np_index_refit_hint(quote(models[[1]]), se = TRUE), "Replace 'object'", fixed = TRUE)
  expect_match(.np_index_refit_hint(quote(`my fit`), gradients = TRUE), "`my fit`$bws", fixed = TRUE)
  expect_error(.np_index_covariance_inverse(matrix(0, 1, 1)),
               "free-coefficient information matrix is singular", fixed = TRUE)
  expect_identical(.np_index_covariance_inverse(diag(2)), diag(2))
})
