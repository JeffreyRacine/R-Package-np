# Cheap admission proof here; optimized/MPI route combinations are external sentinels.
test_that("binary bandwidth admission uses the joint complete sample", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(906)
  d <- data.frame(x = rnorm(84), z = runif(84), y = rep(0:1, 42))
  d$x[3] <- NA_real_
  d$y[7] <- NA_real_
  for (categorical in c(FALSE, TRUE)) {
    x <- d[c("x", "z")]
    if (categorical) x$g <- factor(rep(letters[1:3], 28))
    start <- if (categorical) c(1, .2, .3, -.1, .8) else c(1, .2, .8)
    b <- npindexbw(xdat = x, ydat = d$y, bws = start,
                  method = "kleinspady", bandwidth.compute = FALSE)
    expect_equal(b$nobs, 82L)
    expect_equal(as.integer(b$rows.omit), c(3L, 7L))
    fresh <- npindexbw(xdat = x[complete.cases(x, d$y), , drop = FALSE],
                      ydat = d$y[complete.cases(x, d$y)], bws = start,
                      method = "kleinspady", bandwidth.compute = FALSE)
    expect_identical(b$beta, fresh$beta)
    expect_identical(b$bw, fresh$bw)
    refit <- npindexbw(xdat = x, ydat = d$y, bws = b, bandwidth.compute = FALSE)
    expect_equal(refit$nobs, 82L)
    lost <- rep(0, nrow(x)); lost[3] <- 1
    expect_error(npindexbw(xdat = x, ydat = lost, bws = start,
                          method = "kleinspady", bandwidth.compute = FALSE),
                 "both 0 and 1")
    expect_error(npindexbw(xdat = x, ydat = lost, bws = b,
                          bandwidth.compute = FALSE), "both 0 and 1")
  }
})

test_that("binary domain validation remains strict before numeric conversion", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(907)
  x <- data.frame(x = rnorm(80), z = runif(80))
  y <- rep(0:1, 40)
  for (invalid in list(factor(y), as.logical(y), as.complex(y), replace(y, 2, 2)))
    expect_error(npindexbw(xdat = x, ydat = invalid, bws = c(1, .2, .8),
                          method = "kleinspady", bandwidth.compute = FALSE),
                 "numeric.*0/1")
  expect_error(npindexbw(xdat = x, ydat = rep(NA_real_, 80), bws = c(1, .2, .8),
                        method = "kleinspady", bandwidth.compute = FALSE),
               "no rows without NAs")
  b <- npindexbw(xdat = x, ydat = y, bws = c(1, .2, .8),
                method = "kleinspady", bandwidth.compute = FALSE)
  f <- npindex(b, exdat = x[1:20, ], eydat = rep(0, 20), se = FALSE)
  expect_true(all(is.finite(fitted(f))))
  expect_equal(f$diagnostics.nobs, 20L)
})
