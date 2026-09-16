test_that("manual quantile widths do not launch a default search", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(903L)
  d <- data.frame(y = rnorm(24), x = runif(24))
  h <- c(.7, .6)
  for (named in c(FALSE, TRUE)) for (compute in list(NULL, FALSE)) {
    args <- if (named) list(formula = y ~ x, data = d, bws = h) else
      list(y ~ x, data = d, bws = h)
    if (!is.null(compute)) args$bandwidth.compute <- compute
    rng <- .Random.seed
    value <- do.call(npqreg, args)
    expect_identical(.Random.seed, rng)
    expect_identical(unname(value$bws$ybw), h[1L])
    expect_identical(unname(value$bws$xbw), h[2L])
    expect_true(is.na(value$bws$fval))
    expect_true(is.na(value$bws$num.feval))
    oracle <- npqreg(txdat = d["x"], tydat = d$y, bws = h)
    expect_identical(value$quantile, oracle$quantile)
  }
})

test_that("automatic quantile searches retain constructor results and RNG", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(904L)
  d <- data.frame(y = rnorm(24), x = runif(24))
  for (named in c(FALSE, TRUE)) {
    args <- if (named) list(formula = y ~ x, data = d, nmulti = 1L) else
      list(y ~ x, data = d, nmulti = 1L)
    set.seed(905L)
    value <- do.call(npqreg, args)
    rng <- .Random.seed
    set.seed(905L)
    bw <- npcdistbw(y ~ x, data = d, nmulti = 1L)
    expect_identical(.Random.seed, rng)
    for (field in c("xbw", "ybw", "fval", "num.feval"))
      expect_identical(value$bws[[field]], bw[[field]])
    expect_true(is.finite(value$bws$fval))
    expect_gt(value$bws$num.feval, 0)
    expect_identical(value$quantile, npqreg(bws = bw)$quantile)
  }
})

test_that("existing explicit search-control contracts are retained", {
  ns <- environment(npqreg)
  if (exists(".npRmpi_autodispatch_materialize_call", ns, inherits = FALSE)) {
    old <- options(np.messages = FALSE)
    on.exit(options(old), add = TRUE)
    set.seed(906L)
    d <- data.frame(y = rnorm(24), x = runif(24))
    h <- c(.7, .6)
    expected.error <- tryCatch(
      npcdistbw(y ~ x, data = d, bws = h, bandwidth.compute = NA),
      error = identity)
    actual.error <- tryCatch(
      npqreg(y ~ x, data = d, bws = h, bandwidth.compute = NA),
      error = identity)
    expect_s3_class(actual.error, "error")
    expect_identical(conditionMessage(actual.error), conditionMessage(expected.error))
    set.seed(907L)
    value <- npqreg(y ~ x, data = d, bws = h, bandwidth.compute = TRUE, nmulti = 1L)
    rng <- .Random.seed
    set.seed(907L)
    bw <- npcdistbw(y ~ x, data = d, bws = h, bandwidth.compute = TRUE, nmulti = 1L)
    expect_identical(.Random.seed, rng)
    for (field in c("xbw", "ybw", "fval", "num.feval"))
      expect_identical(value$bws[[field]], bw[[field]])
    expect_true(is.finite(value$bws$fval))
    expect_gt(value$bws$num.feval, 0)
    expect_identical(value$quantile, npqreg(bws = bw)$quantile)
  } else {
    # Serial's existing estimator-level contract treats supplied widths as
    # fixed even if a constructor search control is also present.
    old <- options(np.messages = FALSE)
    on.exit(options(old), add = TRUE)
    d <- data.frame(y = sin(seq_len(24)), x = seq_len(24) / 24)
    fit <- npqreg(y ~ x, data = d, bws = c(.7, .6), bandwidth.compute = TRUE)
    expect_identical(unname(fit$bws$ybw), .7)
    expect_identical(unname(fit$bws$xbw), .6)
    expect_true(is.na(fit$bws$fval))
  }
})
