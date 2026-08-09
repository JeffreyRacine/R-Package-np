test_that("fused bivariate summation matches direct npudens densities", {
  helper <- getFromNamespace(
    ".np_entropy_bivariate_gaussian_summation", "np"
  )
  x <- c(-2.5, -1, -1, -0.2, 0.3, 1.8, 4.5)
  y <- c(3.2, 0.1, 0.1, -1.7, 0.8, 2.4, 8.0)
  bw.x <- 0.73
  bw.y <- 1.21
  bw.joint <- c(0.81, 1.08)

  density.x <- fitted(npudens(tdat = x, bws = bw.x))
  density.y <- fitted(npudens(tdat = y, bws = bw.y))
  density.joint <- fitted(npudens(
    tdat = cbind(x, y), bws = bw.joint
  ))
  reference <- 0.5 * mean(
    (1 - sqrt(density.x * density.y / density.joint))^2
  )

  expect_equal(
    helper(x, y, bw.x, bw.y, bw.joint),
    reference,
    tolerance = 128 * .Machine$double.eps
  )
})

test_that("npsdeptest summation route retains its public result shape", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(1927)
  data <- c(rexp(22) - 1, -2, -2, 7)
  out <- npsdeptest(
    data = data,
    lag.num = 2,
    method = "summation",
    bootstrap = TRUE,
    boot.num = 9,
    random.seed = 4817
  )

  expect_s3_class(out, "sdeptest")
  expect_identical(dim(out$Srho.bootstrap.mat), c(9L, 2L))
  expect_identical(dim(out$Srho.cumulant.bootstrap.mat), c(9L, 2L))
  expect_length(out$P, 2L)
  expect_length(out$P.cumulant, 2L)
})
