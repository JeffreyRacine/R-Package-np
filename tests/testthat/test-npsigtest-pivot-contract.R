test_that("npsigtest has one logical pivot policy for all tested predictors", {
  plan <- getFromNamespace(".np_npsig_pivot_plan", "np")
  x <- data.frame(u = factor(c("a", "b", "a")),
    o = ordered(c("low", "high", "low")), x = c(.1, .4, .9))
  for (joint in c(FALSE, TRUE)) for (pivot in c(FALSE, TRUE)) {
    p <- plan(pivot, x, 1:3, joint)
    expect_identical(p$requested, pivot)
    expect_identical(p$effective, c(u = pivot, o = pivot, x = pivot))
  }
  for (bad in list(NULL, NA, 1, "TRUE", c(TRUE, FALSE)))
    expect_error(plan(bad, x, 1:3, FALSE), "pivot")
})

test_that("npsigtest statistic retains its denominator and rejects undefined SEs", {
  statistic <- getFromNamespace(".np_npsig_statistic", "np")
  fit <- list(grad = matrix(c(0, 2, 3, 4), 2L,
                            dimnames = list(NULL, c("u", "x"))),
              gerr = matrix(c(0, 1, 1.5, 2), 2L))
  structural <- matrix(c(TRUE, FALSE, FALSE, FALSE), 2L)
  expect_identical(statistic(fit, 1:2, FALSE), mean(fit$grad^2))
  expect_identical(statistic(fit, 1:2, TRUE, structural), mean(c(0, 2, 2, 2)^2))
  expect_identical(statistic(fit, 2L, TRUE), mean((fit$grad[, 2L]/fit$gerr[, 2L])^2))
  expect_error(statistic(fit, 1L, TRUE), "zero standard error.*'u'.*row 1")
  for (bad in c(0, -1, NA_real_, Inf)) {
    fit$gerr[2L, 1L] <- bad
    expect_error(statistic(fit, 1:2, TRUE, structural, "bootstrap replication 7"),
                 "bootstrap replication 7.*standard error.*'u'.*row 2")
  }
  fit$gerr <- NULL
  expect_error(statistic(fit, 1:2, TRUE, structural), "standard errors are unavailable")
  fit$grad[2L, 1L] <- Inf
  expect_error(statistic(fit, 1:2, FALSE), "non-finite gradient estimates")
})

test_that("omitted pivot is TRUE for unordered, ordered, continuous and joint tests", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(270827)
  n <- 48L
  x <- data.frame(u = factor(rep(letters[1:3], length.out = n)),
    o = ordered(rep(1:4, each = 12L)), x = seq(-1, 1, length.out = n))
  y <- .4 * (x$u == "b") + x$x^2 + rnorm(n, sd = .2)
  bw <- npregbw(xdat = x, ydat = y, bws = c(.3, .3, .6),
    bandwidth.compute = FALSE, regtype = "ll")
  for (joint in c(FALSE, TRUE)) {
    implicit <- npsigtest(bw, B = 9L, joint = joint, random.seed = 81)
    explicit <- npsigtest(bw, B = 9L, joint = joint, pivot = TRUE, random.seed = 81)
    expect_identical(implicit$In, explicit$In)
    expect_identical(implicit$In.bootstrap, explicit$In.bootstrap)
    expect_identical(implicit$pivot.effective, c(u = TRUE, o = TRUE, x = TRUE))
    expect_output(print(implicit), "Pivot = TRUE", fixed = TRUE)
  }
  expect_error(npsigtest(bw, B = 9L, pivot = NULL), "pivot")
  raw <- npsigtest(bw, B = 9L, pivot = FALSE)
  expect_identical(raw$pivot.effective, c(u = FALSE, o = FALSE, x = FALSE))
  legacy <- raw
  legacy$pivot <- NULL
  expect_output(print(legacy), "automatic -> FALSE", fixed = TRUE)
})
