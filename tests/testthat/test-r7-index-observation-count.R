test_that("index observation count accepts vector and matrix indices", {
  bws <- list(beta = c(1, .5), bw = .6)
  for (index in list(1:5, matrix(1:5, ncol = 1L))) {
    obj <- singleindex(bws, index = index, mean = 1:5, ntrain = 20L,
                       se = TRUE, merr = rep(.1, 5))
    expect_identical(obj$nobs, 5L)
    expect_identical(predict(obj, se.fit = TRUE)$df, 5L)
  }
})
