test_that("npsigtest has one logical pivot policy for all tested predictors", {
  plan <- getFromNamespace(".np_npsig_pivot_plan", "npRmpi")
  xdat <- data.frame(
    category = factor(c("a", "b", "a")),
    ordered = ordered(c("low", "high", "low")),
    continuous = c(0.1, 0.4, 0.9)
  )

  for (joint in c(FALSE, TRUE)) for (pivot in c(FALSE, TRUE)) {
    policy <- plan(pivot, xdat, 1:3, joint)
    expect_identical(policy$requested, pivot)
    expect_identical(policy$effective,
      c(category = pivot, ordered = pivot, continuous = pivot))
  }
  for (bad in list(NULL, NA, 1, "TRUE", c(TRUE, FALSE)))
    expect_error(plan(bad, xdat, 1:3, FALSE), "pivot")
})

test_that("npsigtest statistic helper implements literal defined arithmetic", {
  statistic <- getFromNamespace(".np_npsig_statistic", "npRmpi")
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
  expect_error(statistic(fit, 2L, TRUE), "zero standard error.*'2'.*row 1")
  fit$gerr[1L, 2L] <- NA_real_
  expect_error(statistic(fit, 2L, TRUE), "non-finite standard error.*'2'.*row 1")
  fit$gerr <- NULL
  expect_error(statistic(fit, 2L, TRUE), "standard errors are unavailable")
  fit$grad[1L, 2L] <- Inf
  expect_error(statistic(fit, 2L, FALSE), "non-finite gradient estimates")
})

test_that("npsigtest public pivot modes agree under MPI autodispatch", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_options <- options(np.messages = FALSE, npRmpi.autodispatch = TRUE)
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

  categorical.auto <- npsigtest(bw, B = 9, index = 1L, random.seed = 81)
  categorical.pivot <- npsigtest(
    bw, B = 9, index = 1L, pivot = TRUE, random.seed = 81
  )
  categorical.raw <- npsigtest(
    bw, B = 9, index = 1L, pivot = FALSE, random.seed = 81
  )
  expect_identical(categorical.auto$In, categorical.pivot$In)
  expect_identical(categorical.auto$In.bootstrap, categorical.pivot$In.bootstrap)
  expect_identical(categorical.auto$pivot.effective, c(z = TRUE))
  raw.fit <- npreg(bws = bw, txdat = data.frame(z, x), tydat = y,
                   gradients = TRUE, se = FALSE)
  expect_identical(categorical.raw$In, mean(raw.fit$grad[, 1L]^2))
  expect_identical(categorical.raw$pivot.effective, c(z = FALSE))

  continuous.auto <- npsigtest(bw, B = 9, index = 2L, random.seed = 82)
  continuous.pivot <- npsigtest(
    bw, B = 9, index = 2L, pivot = TRUE, random.seed = 82
  )
  expect_identical(continuous.auto$In, continuous.pivot$In)
  expect_identical(continuous.auto$In.bootstrap, continuous.pivot$In.bootstrap)
  expect_identical(continuous.auto$pivot.effective, c(x = TRUE))

  joint.auto <- npsigtest(
    bw, B = 9, index = 1:2, joint = TRUE, random.seed = 83
  )
  joint.pivot <- npsigtest(
    bw, B = 9, index = 1:2, joint = TRUE,
    pivot = TRUE, random.seed = 83
  )
  expect_identical(joint.auto$In, joint.pivot$In)
  expect_identical(joint.auto$In.bootstrap, joint.pivot$In.bootstrap)
  expect_identical(joint.auto$pivot.effective, c(z = TRUE, x = TRUE))
  joint.raw <- npsigtest(
    bw, B = 9, index = 1:2, joint = TRUE,
    pivot = FALSE, random.seed = 83
  )
  expect_identical(joint.raw$In, mean(raw.fit$grad^2))

  expect_error(
    npsigtest(bw, B = 9, index = 1L, pivot = NULL),
    "pivot"
  )
  expect_output(print(joint.auto), "Pivot = TRUE", fixed = TRUE)

  legacy <- continuous.pivot
  legacy$pivot.effective <- NULL
  legacy.output <- capture.output(print(legacy))
  expect_true(any(grepl("Pivot = TRUE", legacy.output, fixed = TRUE)))
  expect_false(any(grepl("automatic", legacy.output, fixed = TRUE)))
})
