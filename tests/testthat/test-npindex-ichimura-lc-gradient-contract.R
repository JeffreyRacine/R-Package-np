test_that("npRmpi npindex ichimura lc gradients preserve the scalar-index fit contract", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 1000L
  x1 <- runif(n, min = -1, max = 1)
  x2 <- runif(n, min = -1, max = 1)
  y <- x1 - x2 + rnorm(n)
  dat <- data.frame(x1 = x1, x2 = x2, y = y)

  bw <- npindexbw(formula = y ~ x1 + x2, data = dat)
  fit <- npindex(bws = bw, gradients = TRUE)

  index.df <- data.frame(index = as.vector(as.matrix(dat[c("x1", "x2")]) %*% bw$beta))
  oracle <- do.call(npreg, list(
    txdat = index.df,
    tydat = dat$y,
    exdat = index.df,
    bws = bw$bw,
    bwtype = bw$type,
    ckertype = bw$ckertype,
    ckerorder = bw$ckerorder,
    regtype = bw$regtype,
    gradients = TRUE,
    warn.glp.gradient = FALSE
  ))

  expect_true(all(is.finite(fit$grad)))
  expect_true(all(is.finite(fit$betavcov)))
  expect_equal(as.numeric(fit$mean), as.numeric(oracle$mean), tolerance = 0)
  expect_equal(
    unname(as.matrix(fit$grad)),
    unname(as.matrix(oracle$grad) %*% t(as.vector(bw$beta))),
    tolerance = 0
  )
})
