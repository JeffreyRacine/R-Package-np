test_that("Klein-Spady validation preserves numeric binary values without coercion", {
  check <- getFromNamespace(".npindex_check_binary_response", "npRmpi")
  expect_null(check(c(0, 1), "test", require.both = TRUE))
  expect_null(check(c(0L, 1L), "test", require.both = TRUE))
  expect_null(check(c(1, 1), "test"))
  expect_null(check(c(0, NA, 1), "test"))
  expect_error(check(c(1, 1), "test", require.both = TRUE), "both 0 and 1")
  expect_error(check(c(0, NA, 1), "test", require.both = TRUE), "both 0 and 1")
  for (y in list(factor(c(0, 1)), ordered(c(0, 1)), c("0", "1"),
                 c(FALSE, TRUE), c(0+0i, 1+0i), c(0, .5, 1), c(0, 2),
                 c(0, Inf), c(0, -Inf), matrix(c(0, 1), 2L))) {
    expect_error(check(y, "test"), "requires a numeric", fixed = TRUE)
  }
})

test_that("Klein-Spady rejects invalid response inputs before search and held fitting", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260906)
  n <- 40L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1),
                  x3 = runif(n, -1, 1))
  index <- x$x1 + .6*x$x2 - .25*x$x3
  y <- as.double(index + sin(seq_len(n)*2.3) > 0)
  b <- npindexbw(xdat = x, ydat = y, method = "kleinspady",
                 bws = c(1, .6, -.25, .8), bandwidth.compute = FALSE)
  for (bad in list(factor(y), ordered(y), as.character(y), y == 1,
                  y + 1, replace(y, 1L, .5), replace(y, 1L, Inf))) {
    expect_error(npindexbw(xdat = x, ydat = bad, method = "kleinspady"),
                 "requires a numeric", fixed = TRUE)
    expect_error(npindexbw(xdat = x, ydat = bad, bws = b),
                 "requires a numeric", fixed = TRUE)
    expect_error(npindex(bws = b, txdat = x, tydat = bad, se = FALSE),
                 "requires a numeric", fixed = TRUE)
  }
  for (bad in list(factor(y[1:7]), rep(.5, 7))) {
    expect_error(npindex(bws = b, txdat = x, tydat = y,
                         exdat = x[1:7, ], eydat = bad, se = FALSE),
                 "requires a numeric", fixed = TRUE)
  }
  fit <- npindex(bws = b, txdat = x, tydat = y, gradients = TRUE)
  expect_true(all(is.finite(vcov(fit))))
  integer.fit <- npindex(bws = b, txdat = x, tydat = as.integer(y),
                         gradients = TRUE)
  expect_identical(fitted(fit), fitted(integer.fit))
  expect_identical(vcov(fit), vcov(integer.fit))
  single.class <- npindex(bws = b, txdat = x, tydat = y,
    exdat = x[1:7, ], eydat = rep(1, 7), se = FALSE)
  expect_true(all(is.finite(fitted(single.class))))
  dat <- data.frame(y = factor(y), x)
  expect_error(npindex(y ~ x1+x2+x3, data = dat, method = "kleinspady"),
               "requires a numeric", fixed = TRUE)
  x$x3 <- factor(rep(c("a", "b"), length.out = n))
  b.i <- npindexbw(xdat = x, ydat = sin(index), bws = c(1,.6,-.25,.8),
                   bandwidth.compute = FALSE)
  expect_true(all(is.finite(fitted(npindex(bws = b.i, txdat = x,
    tydat = sin(index), se = FALSE)))))
})
