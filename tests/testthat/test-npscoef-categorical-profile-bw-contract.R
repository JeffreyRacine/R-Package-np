test_that("npscoef all-categorical coefficient profiles preserve bandwidth CV", {
  old.tree <- getOption("np.tree")
  old.compress <- getOption("np.categorical.compress")
  old.messages <- getOption("np.messages")
  on.exit({
    options(np.tree = old.tree)
    options(np.categorical.compress = old.compress)
    options(np.messages = old.messages)
  }, add = TRUE)

  set.seed(1417)
  n <- 192L
  z1 <- factor(rbinom(n, 1L, 0.45))
  z2 <- ordered(sample(0:3, n, replace = TRUE))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 1 + 0.5 * x1 - 0.25 * x2 +
    as.numeric(z1) + as.numeric(z2) + rnorm(n, sd = 0.2)
  dat <- data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)

  options(np.messages = FALSE, np.tree = FALSE,
          np.categorical.compress = FALSE)
  dense <- npscoefbw(y ~ x1 + x2 | z1 + z2,
                     data = dat, regtype = "lc", nmulti = 1)

  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npscoefbw(y ~ x1 + x2 | z1 + z2,
                       data = dat, regtype = "lc", nmulti = 1)

  expect_true(is.finite(profile$fval))
  expect_true(all(is.finite(profile$bw)))
  expect_true(abs(profile$fval - dense$fval) < 1e-3)
})
