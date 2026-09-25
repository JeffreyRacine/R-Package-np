test_that("quantile omission metadata separates training and output domains", {
  make <- getFromNamespace("qregression", "npRmpi")
  b <- list(xbw = .3, ybw = .3)
  for (training in c(FALSE, TRUE)) {
    fit <- make(b, data.frame(x = 1:5), .5, 1:5, ntrain = 22L,
                trainiseval = training, train.rows.omit = c(3L, 7L),
                eval.rows.omit = if (training) integer(0) else 2L)
    expect_identical(fit$train.rows.omit, c(3L, 7L))
    expect_identical(fit$eval.rows.omit, if (training) integer(0) else 2L)
    expect_identical(fit$rows.omit, if (training) c(3L, 7L) else 2L)
    expect_identical(fit$nobs.omit, if (training) 2L else 1L)
    expect_identical(fitted(fit), 1:5)
  }
})
