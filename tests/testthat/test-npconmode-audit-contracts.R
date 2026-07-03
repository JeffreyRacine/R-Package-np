.npconmode_audit_data <- function(n = 36L) {
  x <- seq(-1, 1, length.out = n)
  p <- plogis(1.4 * x)
  data.frame(
    x = x,
    y = factor(ifelse(seq_len(n) / n <= p, "b", "a"), levels = c("a", "b"))
  )
}

.npconmode_audit_bw <- function(d) {
  npcdensbw(
    xdat = d["x"],
    ydat = d["y"],
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
}

.npconmode_expect_missing_row <- function(fit, row, n) {
  expect_equal(length(fit$conmode), n)
  expect_true(is.na(fit$conmode[row]))
  expect_true(is.na(fit$condens[row]))
  expect_true(is.na(fit$conderr[row]))
  if (!is.null(fit$probabilities))
    expect_true(all(is.na(fit$probabilities[row, ])))
  if (!is.null(fit$probability.errors))
    expect_true(all(is.na(fit$probability.errors[row, ])))
  if (!is.null(fit$probability.gradients))
    expect_true(all(is.na(fit$probability.gradients[row, ])))
  expect_identical(fit$nobs, n)
}

test_that("npconmode modal selector preserves row length and first positive ties", {
  select <- getFromNamespace(".npConmodeSelect", "np")
  selected.factor <- getFromNamespace(".npConmodeSelectedFactor", "np")

  pmat <- matrix(c(
    NaN, NaN,
    0.0, 0.0,
    -0.3, -0.1,
    0.4, 0.4,
    0.2, 0.6
  ), ncol = 2L, byrow = TRUE)
  perr <- matrix(seq_len(length(pmat)) / 100, nrow = nrow(pmat))

  out <- select(pmat, perr)

  expect_identical(out$indices, c(0L, 0L, 0L, 1L, 2L))
  expect_true(all(is.na(out$condens[1:3])))
  expect_true(all(is.na(out$conderr[1:3])))
  expect_equal(out$condens[4:5], c(0.4, 0.6))

  classes <- selected.factor(out$indices, factor(c("a", "b"), levels = c("a", "b")))
  expect_equal(length(classes), nrow(pmat))
  expect_true(all(is.na(classes[1:3])))
  expect_identical(as.character(classes[4:5]), c("a", "b"))
})

test_that("npconmode rejects conditional-distribution bandwidth objects clearly", {
  d <- .npconmode_audit_data()
  cdf.bw <- npcdistbw(
    xdat = d["x"],
    ydat = data.frame(y = as.numeric(d$y == "b")),
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE
  )

  expect_error(
    npconmode(bws = cdf.bw, txdat = d["x"], tydat = d["y"]),
    "expected conditional density bandwidths instead of conditional distribution bandwidths",
    fixed = TRUE
  )
})

test_that("npconmode direct route pads omitted training rows", {
  d <- .npconmode_audit_data()
  bw <- .npconmode_audit_bw(d)
  tx <- d["x"]
  tx$x[4L] <- NA_real_

  fit <- npconmode(
    bws = bw,
    txdat = tx,
    tydat = d["y"],
    probabilities = TRUE,
    gradients = TRUE,
    level = "b"
  )

  .npconmode_expect_missing_row(fit, row = 4L, n = nrow(d))
  expect_identical(fit$rows.omit, 4L)
  expect_identical(fit$nobs.omit, 1L)
  expect_true(is.null(fit$eval.rows.omit))
  expect_equal(nrow(predict(fit, type = "prob")), nrow(d))
  expect_equal(NROW(gradients(fit)), nrow(d))
})

test_that("npconmode direct route pads omitted evaluation rows separately", {
  d <- .npconmode_audit_data()
  bw <- .npconmode_audit_bw(d)
  ex <- data.frame(x = c(-0.75, NA_real_, 0.75))

  fit <- npconmode(
    bws = bw,
    txdat = d["x"],
    tydat = d["y"],
    exdat = ex,
    probabilities = TRUE,
    gradients = TRUE,
    level = "b"
  )

  .npconmode_expect_missing_row(fit, row = 2L, n = nrow(ex))
  expect_identical(fit$rows.omit, NA_integer_)
  expect_identical(fit$nobs.omit, 0L)
  expect_identical(fit$eval.rows.omit, 2L)
  expect_identical(fit$eval.nobs.omit, 1L)
  expect_equal(nrow(predict(fit, type = "prob")), nrow(ex))
  expect_equal(NROW(gradients(fit)), nrow(ex))
})

test_that("npconmode formula route separates training and evaluation omissions", {
  d <- .npconmode_audit_data()
  d.train <- d
  d.train$x[5L] <- NA_real_

  fit.train <- npconmode(
    y ~ x,
    data = d.train,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE,
    regtype = "ll",
    probabilities = TRUE,
    gradients = TRUE,
    level = "b"
  )

  .npconmode_expect_missing_row(fit.train, row = 5L, n = nrow(d.train))
  expect_identical(as.integer(fit.train$rows.omit), 5L)
  expect_identical(fit.train$nobs.omit, 1L)

  nd <- data.frame(x = c(-0.6, NA_real_, 0.6))
  fit.eval <- npconmode(
    y ~ x,
    data = d,
    newdata = nd,
    bws = c(0.45, 0.45),
    bandwidth.compute = FALSE,
    regtype = "ll",
    probabilities = TRUE,
    gradients = TRUE,
    level = "b"
  )

  .npconmode_expect_missing_row(fit.eval, row = 2L, n = nrow(nd))
  expect_identical(fit.eval$rows.omit, NA_integer_)
  expect_identical(fit.eval$nobs.omit, 0L)
  expect_identical(as.integer(fit.eval$eval.rows.omit), 2L)
  expect_identical(fit.eval$eval.nobs.omit, 1L)
})
