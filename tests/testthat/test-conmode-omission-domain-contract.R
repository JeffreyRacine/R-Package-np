test_that("conditional mode records active and training omission domains separately", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24))
  y <- factor(rep(c("a", "b"), 12))
  x$x[c(3, 7)] <- NA_real_
  e <- data.frame(x = c(-.6, NA, .1, .4, .7))
  b <- npcdensbw(xdat = x, ydat = y, bws = c(.2, .5),
                 bandwidth.compute = FALSE)
  fit <- npconmode(b, txdat = x, tydat = y, exdat = e,
                   probabilities = TRUE)
  expect_identical(as.integer(fit$rows.omit), 2L)
  expect_identical(as.integer(fit$omit), 2L)
  expect_identical(fit$nobs.omit, 1L)
  expect_identical(as.integer(fit$train.rows.omit), c(3L, 7L))
  expect_identical(fit$train.nobs.omit, 2L)
  expect_identical(as.integer(fit$eval.rows.omit), 2L)
  expect_identical(fit$eval.nobs.omit, 1L)
  expect_true(is.na(fit$conmode[2]))
  expect_identical(length(fit$conmode), nrow(e))
  clean <- npconmode(b, txdat = x[-c(3, 7), , drop = FALSE],
                    tydat = y[-c(3, 7)], exdat = e[-2, , drop = FALSE],
                    probabilities = TRUE)
  expect_identical(fit$conmode[-2], clean$conmode)
  expect_equal(unname(fit$probabilities[-2, , drop = FALSE]),
               unname(clean$probabilities), tolerance = 0)
  training <- npconmode(b, txdat = x, tydat = y)
  expect_identical(as.integer(training$rows.omit), c(3L, 7L))
  expect_identical(training$train.rows.omit, training$rows.omit)
  expect_null(training$eval.nobs.omit)
  expect_null(training$eval.rows.omit)
})

test_that("formula and native newdata keep their own omission coordinates", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 24),
                  y = factor(rep(c("a", "b"), 12)))
  d$x[c(3, 7)] <- NA_real_
  e <- data.frame(x = c(-.6, NA, .1, .4, .7))
  b <- npcdensbw(y ~ x, data = d, bws = c(.2, .5),
                 bandwidth.compute = FALSE, na.action = na.exclude)
  for (native in c(FALSE, TRUE)) {
    fit <- if (native) npconmode(b, exdat = e, probabilities = TRUE) else
      npconmode(b, newdata = e, probabilities = TRUE)
    expect_identical(as.integer(fit$rows.omit), 2L)
    expect_identical(as.integer(fit$omit), 2L)
    expect_identical(fit$nobs.omit, 1L)
    expect_identical(as.integer(fit$train.rows.omit), c(3L, 7L))
    expect_identical(fit$train.nobs.omit, 2L)
    expect_identical(as.integer(fit$eval.rows.omit), 2L)
    expect_identical(fit$eval.nobs.omit, 1L)
    expect_true(is.na(fit$conmode[2]))
    expect_identical(length(fit$conmode), nrow(e))
    saved <- unserialize(serialize(fit, NULL))
    expect_identical(predict(saved, newdata = e, type = "prob"),
                     predict(fit, newdata = e, type = "prob"))
  }
  training <- npconmode(b)
  expect_identical(as.integer(training$rows.omit), c(3L, 7L))
  expect_identical(training$train.nobs.omit, 2L)
  expect_null(training$eval.nobs.omit)
})
