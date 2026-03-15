suppressPackageStartupMessages(library(np))

test_that("nonfixed regression bandwidth search respects empirical support bounds", {
  set.seed(123)

  x <- sample(rep(1:8, each = 10), 60, replace = FALSE)
  y <- x + rnorm(60)
  support_kmax <- length(unique(x)) - 1L

  bw.adap <- npregbw(ydat = y, xdat = data.frame(x = x), bwtype = "adaptive_nn")
  bw.gen <- npregbw(ydat = y, xdat = data.frame(x = x), bwtype = "generalized_nn")

  expect_gte(unname(bw.adap$bw), 1)
  expect_lte(unname(bw.adap$bw), support_kmax)
  expect_gte(unname(bw.gen$bw), 1)
  expect_lte(unname(bw.gen$bw), support_kmax)
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
