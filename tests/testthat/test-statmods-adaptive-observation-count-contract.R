suppressPackageStartupMessages(library(np))

test_that("adaptive duplicate mass promotes zero-radius neighborhoods to the nearest positive radius", {
  x <- data.frame(x = c(0, 0, 0, 1))
  y <- c(0, 1, 2, 3)

  bw2 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "adaptive_nn",
    bws = 2,
    bandwidth.compute = FALSE
  )
  bw3 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "adaptive_nn",
    bws = 3,
    bandwidth.compute = FALSE
  )

  fit2 <- npreg(bws = bw2, exdat = x)
  fit3 <- npreg(bws = bw3, exdat = x)

  expect_equal(fit2$mean, fit3$mean, tolerance = 0)
})

test_that("adaptive observation-count bandwidths handle tied supports that formerly failed", {
  x <- data.frame(x = c(0, 1, 1, 1, 2))
  y <- c(0, 1, 2, 3, 4)

  bw <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "adaptive_nn",
    bws = 3,
    bandwidth.compute = FALSE
  )

  fit <- npreg(bws = bw, exdat = x)

  expect_length(fit$mean, nrow(x))
  expect_true(all(is.finite(fit$mean)))
})

test_that("generalized duplicate mass promotes zero-radius neighborhoods to the nearest positive radius", {
  x <- data.frame(x = c(0, 0, 0, 1))
  y <- c(0, 1, 2, 3)

  bw2 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 2,
    bandwidth.compute = FALSE
  )
  bw3 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 3,
    bandwidth.compute = FALSE
  )

  fit2 <- npreg(bws = bw2, exdat = x)
  fit3 <- npreg(bws = bw3, exdat = x)

  expect_equal(fit2$mean, fit3$mean, tolerance = 0)
})

test_that("generalized observation-count bandwidths handle low-support ties and out-of-support evaluation", {
  x <- data.frame(x = c(0, 0, 0, 1, 1, 2))
  y <- c(0, 1, 2, 3, 4, 5)
  ex <- data.frame(x = c(0, 0.5, 1, 1.5, 2))

  bw <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 3,
    bandwidth.compute = FALSE
  )

  fit <- npreg(bws = bw, exdat = ex)

  expect_length(fit$mean, nrow(ex))
  expect_true(all(is.finite(fit$mean)))
})
