suppressPackageStartupMessages(library(np))

test_that("adaptive duplicate mass rejects a literal zero-radius neighborhood", {
  x <- data.frame(x = c(0, 0, 0, 1))
  y <- c(0, 1, 2, 3)

  bw2 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "adaptive_nn",
    bws = 2,
    bandwidth.compute = FALSE
  )
  expect_error(
    npreg(bws = bw2, exdat = x),
    "zero literal radius",
    fixed = TRUE
  )
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

test_that("generalized duplicate mass reports a literal zero-radius query", {
  x <- data.frame(x = c(0, 0, 0, 1))
  y <- c(0, 1, 2, 3)

  bw2 <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 2,
    bandwidth.compute = FALSE
  )
  expect_error(
    npreg(bws = bw2, exdat = x),
    "zero literal radius",
    fixed = TRUE
  )
})

test_that("generalized external queries reject zero-radius tied supports", {
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

  expect_error(
    npreg(bws = bw, exdat = ex),
    "zero literal radius",
    fixed = TRUE
  )
})
