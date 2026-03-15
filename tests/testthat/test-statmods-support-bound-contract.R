suppressPackageStartupMessages(library(np))

test_that("nonfixed search uses observation-count bounds for adaptive and generalized", {
  set.seed(123)

  x <- sample(rep(1:8, each = 10), 60, replace = FALSE)
  y <- x + rnorm(60)
  n_kmax <- length(x) - 1L

  bw.adap <- npregbw(ydat = y, xdat = data.frame(x = x), bwtype = "adaptive_nn")
  bw.gen <- npregbw(ydat = y, xdat = data.frame(x = x), bwtype = "generalized_nn")

  expect_gte(unname(bw.adap$bw), 1)
  expect_lte(unname(bw.adap$bw), n_kmax)
  expect_gte(unname(bw.gen$bw), 1)
  expect_lte(unname(bw.gen$bw), n_kmax)
})

test_that("explicit k may exceed empirical support size when observation count allows it", {
  x <- c(0, 0, 0, 1)
  y <- c(0, 1, 2, 3)

  bw.adap <- npregbw(
    ydat = y,
    xdat = data.frame(x = x),
    bwtype = "adaptive_nn",
    bws = length(x) - 1L,
    bandwidth.compute = FALSE
  )
  bw.gen <- npregbw(
    ydat = y,
    xdat = data.frame(x = x),
    bwtype = "generalized_nn",
    bws = length(x) - 1L,
    bandwidth.compute = FALSE
  )

  expect_no_error(
    npreg(bws = bw.adap, exdat = data.frame(x = x))
  )

  expect_no_error(
    npreg(bws = bw.gen, exdat = data.frame(x = x)),
  )
})

test_that("nonfixed regression bandwidth search fails fast for constant regressors", {
  expect_error(
    suppressWarnings(
      npregbw(ydat = rnorm(50), xdat = data.frame(x = rep(1, 50)), bwtype = "adaptive_nn")
    ),
    "at least two distinct continuous regressor values per dimension"
  )
})

test_that("nonfixed conditional bandwidth search fails fast for constant response-side support", {
  x <- data.frame(x = rep(1:5, 8))
  y <- rep(1, 40)

  expect_error(
    suppressWarnings(
      npcdensbw(ydat = y, xdat = x, bwtype = "adaptive_nn")
    ),
    "at least two distinct continuous variable values per dimension"
  )
})
