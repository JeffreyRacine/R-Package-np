test_that("single-index stochastic terms and subsets use the prepared values", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(929L)
  d <- data.frame(x = runif(32), z = rnorm(32), y = rnorm(32))
  jittered <- function(x) x + rnorm(length(x), sd = .01)
  f <- y ~ jittered(x) + z
  set.seed(930L)
  x <- data.frame(jittered = jittered(d$x), z = d$z)
  names(x) <- attr(terms(f), "term.labels")
  expected.rng <- .Random.seed
  h <- c(1, .2, .8)
  for (route in c("constructor", "one-call")) {
    set.seed(930L)
    if (route == "constructor") {
      bw <- npindexbw(f, data = d, bws = h, bandwidth.compute = FALSE)
      expect_identical(.Random.seed, expected.rng)
      # This is a separate native-data fit, not an independent formula refit.
      fit <- npindex(bw, txdat = x, tydat = d$y, se = FALSE)
    } else {
      fit <- npindex(f, data = d, bws = h, se = FALSE)
      expect_identical(.Random.seed, expected.rng)
      bw <- fit$bws
    }
    oracle <- npindex(bw, txdat = x, tydat = d$y, se = FALSE)
    expect_identical(fitted(fit), fitted(oracle))
  }
  d$y[7L] <- NA_real_
  bw <- npindexbw(y ~ x + z, data = d, subset = x > .2,
    na.action = na.omit, bws = h, bandwidth.compute = FALSE)
  mf <- model.frame(y ~ x + z, data = d, subset = x > .2, na.action = na.omit)
  fit <- npindex(bw, se = FALSE)
  oracle <- npindex(bw, txdat = mf[c("x", "z")], tydat = mf$y, se = FALSE)
  expect_identical(fitted(fit), fitted(oracle))
  expect_identical(fit$ntrain, nrow(mf))
  expect_identical(bw$rows.omit, as.vector(attr(mf, "na.action")))
})
