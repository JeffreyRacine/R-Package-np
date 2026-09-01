suppressPackageStartupMessages(library(np))

test_that("adaptive NN direct estimators reject literal zero radii", {
  x <- data.frame(x = c(0, 0, 0, 1, 2))
  y <- c(0, 1, 2, 3, 4)
  response <- data.frame(y = c(-2, -1, 0, 1, 2))

  regression_bw <- npregbw(
    xdat = x, ydat = y, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  )
  density_bw <- npudensbw(
    dat = x, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  )
  distribution_bw <- npudistbw(
    dat = x, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  )
  conditional_density_bw <- npcdensbw(
    xdat = x, ydat = response, bwtype = "adaptive_nn", bws = c(4, 2),
    bandwidth.compute = FALSE
  )
  conditional_distribution_bw <- npcdistbw(
    xdat = x, ydat = response, bwtype = "adaptive_nn", bws = c(4, 2),
    bandwidth.compute = FALSE
  )

  expect_error(npreg(bws = regression_bw, exdat = x),
               "zero literal radius", fixed = TRUE)
  expect_error(npudens(bws = density_bw, tdat = x),
               "zero literal radius", fixed = TRUE)
  expect_error(npudist(bws = distribution_bw, tdat = x),
               "zero literal radius", fixed = TRUE)
  expect_error(
    npcdens(bws = conditional_density_bw, txdat = x, tydat = response,
            exdat = x, eydat = response),
    "zero literal radius", fixed = TRUE
  )
  expect_error(
    npcdist(bws = conditional_distribution_bw, txdat = x, tydat = response,
            exdat = x, eydat = response),
    "zero literal radius", fixed = TRUE
  )
})

test_that("beta and kernel-sum direct owners retain literal zero status", {
  x <- data.frame(x = c(0.1, 0.1, 0.1, 0.6, 0.9))
  y <- seq_len(nrow(x))

  for (mode in c("generalized_nn", "adaptive_nn")) {
    beta_bw <- npregbw(
      xdat = x, ydat = y, bwtype = mode, bws = 2,
      bandwidth.compute = FALSE, ckertype = "beta",
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    expect_error(npreg(bws = beta_bw, exdat = x),
                 "zero literal radius", fixed = TRUE)
    expect_error(
      npksum(txdat = x, exdat = x, bws = 2, bwtype = mode,
             ckertype = "gaussian"),
      "zero literal radius", fixed = TRUE
    )
    expect_error(
      npksum(txdat = x, exdat = x, bws = 2, bwtype = mode,
             ckertype = "beta", ckerbound = "fixed",
             ckerlb = 0, ckerub = 1),
      "zero literal radius", fixed = TRUE
    )
  }
})

test_that("positive represented subnormal NN radii remain valid geometry", {
  positive_subnormal <- .Machine$double.xmin / 2
  train <- c(0, positive_subnormal, 1)

  expect_gt(positive_subnormal, 0)
  expect_lt(positive_subnormal, .Machine$double.xmin)
  expect_true(.Call(
    "C_np_regression_k1_geometry_validate",
    as.double(train), as.double(train), TRUE,
    PACKAGE = "np"
  ))
})
