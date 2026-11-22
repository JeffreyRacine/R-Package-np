test_that("partial-linearity test is exported and both result methods are registered", {
  ns <- asNamespace("np")
  expect_identical(getExportedValue("np", "nppltest"), ns$nppltest)
  for (generic in c("print", "summary"))
    expect_identical(getS3method(generic, "pltest", envir = ns),
                     get(paste0(generic, ".pltest"), envir = ns, inherits = FALSE))
  expect_error(nppltest(NULL, B = TRUE), "B must")
  expect_error(nppltest(NULL), "object must")
})

test_that("mixed training operator uses the fitted categorical weights", {
  set.seed(481)
  n <- 30L
  x <- data.frame(x = rnorm(n), w = factor(rep(1:2, length.out = n)))
  z <- data.frame(group = factor(rep(1:2, each = n/2)), v = rnorm(n))
  y <- rnorm(n)
  responses <- cbind(y, cos(y), seq_len(n)/n)
  for (reg in c("lc", "ll")) {
    bw <- npplregbw(xdat = x, ydat = y, zdat = z, bandwidth.compute = FALSE,
                    regtype = reg, bws = rbind(c(.496, .8), c(.2, .55), c(.496, .6)))
    H <- .pltest_operator(bw, x, y, z)$H
    oracle <- apply(responses, 2L, function(v)
      as.numeric(fitted(npplreg(bws = bw, txdat = x, tydat = v, tzdat = z))))
    expect_equal(unname(H %*% responses), unname(oracle), tolerance = 1e-10)
    expect_equal(.pltest_residual_norms(H), rowSums((diag(n) - H)^2))
  }
})
