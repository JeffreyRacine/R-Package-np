test_that("density geometry does not survive a search data owner", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = TRUE,
                 np.largeh.rel.tol = 0.001)
  on.exit(options(old), add = TRUE)
  # Nonconstant narrow support primes the all-large cache without invoking
  # the unrelated constant-data warning contract.
  inputs <- list((1:8) * 1e-6, rep(0:1, each = 4L), (1:8) / 8)
  for (x in rep(inputs, 2L)) {
    dat <- data.frame(x = x)
    bw <- npudensbw(dat = dat, bws = 0.5, bwscaling = FALSE,
                   bandwidth.compute = FALSE)
    raw <- np:::npudensbw.bandwidth(dat = dat, bws = bw,
      eval.only = TRUE, invalid.penalty = "dbmax")
    weight <- dnorm(outer(x, x, "-") / 0.5) / 0.5
    if (diff(range(x)) / 0.5 <= sqrt(-2 * log(1 - 0.001)))
      weight[] <- dnorm(0) / 0.5
    diag(weight) <- 0
    expect_equal(as.double(raw$fval), sum(log(rowSums(weight) / 7)),
                 tolerance = 2e-12)
  }
})
