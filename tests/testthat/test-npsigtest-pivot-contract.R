test_that("npsigtest pivot policy is method-aware", {
  plan <- getFromNamespace(".np_npsig_pivot_plan", "np")
  xdat <- data.frame(
    category = factor(c("a", "b", "a")),
    ordered = ordered(c("low", "high", "low")),
    continuous = c(0.1, 0.4, 0.9)
  )

  individual <- plan(NULL, xdat, 1:3, joint = FALSE)
  expect_null(individual$requested)
  expect_identical(
    individual$effective,
    c(category = FALSE, ordered = FALSE, continuous = TRUE)
  )
  expect_identical(
    plan(NULL, xdat, c(1L, 3L), joint = TRUE)$effective,
    c(category = FALSE, continuous = FALSE)
  )
  expect_identical(
    plan(FALSE, xdat, c(1L, 3L), joint = FALSE)$effective,
    c(category = FALSE, continuous = FALSE)
  )
  expect_identical(
    plan(TRUE, xdat, 3L, joint = FALSE)$effective,
    c(continuous = TRUE)
  )
  expect_error(
    plan(TRUE, xdat, c(1L, 3L), joint = FALSE),
    "published categorical test is unstandardized",
    fixed = TRUE
  )
})

test_that("npsigtest statistic helper implements literal defined arithmetic", {
  statistic <- getFromNamespace(".np_npsig_statistic", "np")
  fit <- list(
    grad = matrix(c(1, 2, 3, 4), nrow = 2L),
    gerr = matrix(c(0.5, 1, 1.5, 2), nrow = 2L)
  )

  expect_identical(statistic(fit, 1L, FALSE), mean(fit$grad[, 1L]^2))
  expect_identical(
    statistic(fit, 2L, TRUE),
    mean((fit$grad[, 2L] / fit$gerr[, 2L])^2)
  )

  fit$gerr[1L, 2L] <- 0
  expect_error(statistic(fit, 2L, TRUE), "non-positive or non-finite")
  fit$gerr[1L, 2L] <- NA_real_
  expect_error(statistic(fit, 2L, TRUE), "non-positive or non-finite")
  fit$gerr <- NULL
  expect_error(statistic(fit, 2L, TRUE), "standard errors are unavailable")
  fit$grad[1L, 2L] <- Inf
  expect_error(statistic(fit, 2L, FALSE), "non-finite gradient estimates")
})

test_that("npsigtest public pivot modes agree with explicit counterparts", {
  old_options <- options(np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(270827)
  n <- 45L
  z <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- seq(-1, 1, length.out = n)
  y <- 0.4 * (z == "b") + x^2 + rnorm(n, sd = 0.15)
  bw <- npregbw(
    xdat = data.frame(z, x),
    ydat = y,
    bws = c(0.35, 0.3),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  categorical.auto <- npsigtest(bw, boot.num = 9, index = 1L, random.seed = 81)
  categorical.raw <- npsigtest(
    bw, boot.num = 9, index = 1L, pivot = FALSE, random.seed = 81
  )
  expect_identical(categorical.auto$In, categorical.raw$In)
  expect_identical(categorical.auto$In.bootstrap, categorical.raw$In.bootstrap)
  expect_identical(categorical.auto$pivot.effective, c(z = FALSE))

  continuous.auto <- npsigtest(bw, boot.num = 9, index = 2L, random.seed = 82)
  continuous.pivot <- npsigtest(
    bw, boot.num = 9, index = 2L, pivot = TRUE, random.seed = 82
  )
  expect_identical(continuous.auto$In, continuous.pivot$In)
  expect_identical(continuous.auto$In.bootstrap, continuous.pivot$In.bootstrap)
  expect_identical(continuous.auto$pivot.effective, c(x = TRUE))

  joint.auto <- npsigtest(
    bw, boot.num = 9, index = 1:2, joint = TRUE, random.seed = 83
  )
  joint.raw <- npsigtest(
    bw, boot.num = 9, index = 1:2, joint = TRUE,
    pivot = FALSE, random.seed = 83
  )
  expect_identical(joint.auto$In, joint.raw$In)
  expect_identical(joint.auto$In.bootstrap, joint.raw$In.bootstrap)
  expect_identical(joint.auto$pivot.effective, c(z = FALSE, x = FALSE))

  expect_error(
    npsigtest(bw, boot.num = 9, index = 1L, pivot = TRUE),
    "published categorical test is unstandardized",
    fixed = TRUE
  )
  expect_output(print(joint.auto), "automatic -> FALSE", fixed = TRUE)

  legacy <- continuous.pivot
  legacy$pivot.effective <- NULL
  legacy.output <- capture.output(print(legacy))
  expect_true(any(grepl("Pivot = TRUE", legacy.output, fixed = TRUE)))
  expect_false(any(grepl("automatic", legacy.output, fixed = TRUE)))
})
