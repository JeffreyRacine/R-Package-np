test_that("single-index bandwidth reentry retains untransformed native training", {
  withr::local_options(npRmpi.local.regression.mode = TRUE,
                       npRmpi.autodispatch.disable = TRUE, np.messages = FALSE)
  x <- data.frame(a = seq(.05, .95, length.out = 36), b = sin(seq_len(36)))
  y <- cos(x$a) + x$b
  x$a[3] <- NA_real_
  y[7] <- NA_real_
  initial <- npindexbw(xdat = x, ydat = y, bws = c(1, .25, .35),
                       method = "ichimura", bandwidth.compute = FALSE)
  value <- npindexbw(xdat = x, ydat = y, bws = initial, bandwidth.compute = FALSE)
  expect_identical(value[[".np.native.training"]], list(xdat = x, ydat = y))
  expect_identical(value$beta, initial$beta)
  expect_identical(value$bw, initial$bw)
})
