test_that("LSQ training responses obey the numeric-vector contract", {
  x <- data.frame(x = seq_len(12L) / 12)
  y <- sin(x$x)
  invalid <- list(as.character(y), y > .5, factor(y), ordered(y),
                  as.complex(y), as.list(y))
  for (value in invalid)
    expect_error(.nplsqreg_prepare_train_data(x, value),
                 "numeric (integer or double) response vector", fixed = TRUE)

  integer.y <- seq_len(nrow(x))
  expect_identical(.nplsqreg_prepare_train_data(x, integer.y)$ydat,
                   as.double(integer.y))
  expect_identical(.nplsqreg_prepare_train_data(x, y)$ydat, y)
  y[c(3L, 7L)] <- NA_real_
  prepared <- .nplsqreg_prepare_train_data(x, y, seq_len(nrow(x)))
  expect_identical(prepared$ydat, y[-c(3L, 7L)])
  expect_identical(prepared$scale, seq_len(nrow(x))[-c(3L, 7L)])
  expect_identical(as.integer(prepared$omit), c(3L, 7L))
})

test_that("public LSQ native and formula entries reject coerced responses", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 24L
  x <- seq(.05, .95, length.out = n)
  for (y in list(as.character(sin(x)), x > .5)) {
    expect_error(nplsqregbw(xdat = data.frame(x), ydat = y, bws = .2,
      scale = rep(1, n), delta = .5, regtype = "ll",
      bandwidth.compute = FALSE), "numeric (integer or double)", fixed = TRUE)
    expect_error(nplsqregbw(bws = y ~ x, data = data.frame(x, y), bw = .2,
      scale = rep(1, n), delta = .5, regtype = "ll",
      bandwidth.compute = FALSE), "numeric (integer or double)", fixed = TRUE)
  }
})
