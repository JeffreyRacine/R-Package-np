test_that("named bws formula dispatch matches positional npindex route", {
  x1 <- seq(-0.95, 0.95, length.out = 24)
  x2 <- cos(seq(0, pi, length.out = 24))
  dat <- data.frame(
    y = 1 + 0.7 * x1 - 0.4 * x2 + sin(x1 + 0.5 * x2),
    x1 = x1,
    x2 = x2
  )
  ex <- data.frame(
    x1 = seq(-0.8, 0.8, length.out = 8),
    x2 = seq(0.8, -0.8, length.out = 8)
  )

  fit.pos <- npindex(
    y ~ x1 + x2,
    data = dat,
    newdata = ex,
    method = "ichimura"
  )
  fit.named <- npindex(
    bws = y ~ x1 + x2,
    data = dat,
    newdata = ex,
    method = "ichimura"
  )

  expect_equal(as.numeric(fit.named$mean), as.numeric(fit.pos$mean), tolerance = 0)
  expect_equal(as.numeric(fit.named$merr), as.numeric(fit.pos$merr), tolerance = 0)
})

test_that("formula bws reentry with explicit native data preserves the formula RHS", {
  dat <- data.frame(
    y = seq_len(8),
    x1 = seq(0.1, 0.8, length.out = 8),
    x2 = seq(0.8, 0.1, length.out = 8),
    noise = rev(seq_len(8))
  )
  captured <- new.env(parent = emptyenv())

  fake_bw <- function(xdat, ydat, ...) {
    captured$bw_xnames <- names(xdat)
    captured$bw_ydat <- ydat
    structure(list(xnames = names(xdat), beta = rep.int(1, ncol(xdat))),
              class = "sibandwidth")
  }
  fake_fit <- function(bws, txdat, tydat, ...) {
    dots <- list(...)
    captured$fit_xnames <- names(txdat)
    captured$fit_ydat <- tydat
    captured$fit_exnames <- names(dots$exdat)
    list(bws = bws, mean = numeric(0), merr = numeric(0),
         gradients = FALSE, residuals = FALSE)
  }

  testthat::with_mocked_bindings({
    getFromNamespace("npindex.default", "np")(
      bws = y ~ x1 + x2,
      txdat = dat[c("x1", "x2", "noise")],
      tydat = dat$y,
      exdat = dat[1:3, c("x1", "x2", "noise")]
    )
  },
    npindexbw = fake_bw,
    npindex = fake_fit,
    .package = "np"
  )

  expect_identical(captured$bw_xnames, c("x1", "x2"))
  expect_identical(captured$fit_xnames, c("x1", "x2"))
  expect_identical(captured$fit_exnames, c("x1", "x2"))
  expect_identical(unname(captured$bw_ydat), dat$y)
  expect_identical(unname(captured$fit_ydat), dat$y)
})

test_that("formula bws reentry treats explicit tydat as materialized response data", {
  dat <- data.frame(
    y = seq(2, 9, length.out = 8),
    x1 = seq(0.1, 0.8, length.out = 8),
    x2 = seq(0.8, 0.1, length.out = 8),
    noise = rev(seq_len(8))
  )
  captured <- new.env(parent = emptyenv())

  fake_bw <- function(xdat, ydat, ...) {
    captured$bw_xnames <- names(xdat)
    captured$bw_ydat <- ydat
    structure(list(xnames = names(xdat), beta = rep.int(1, ncol(xdat))),
              class = "sibandwidth")
  }
  fake_fit <- function(bws, txdat, tydat, ...) {
    dots <- list(...)
    captured$fit_xnames <- names(txdat)
    captured$fit_exnames <- names(dots$exdat)
    captured$fit_ydat <- tydat
    list(bws = bws, mean = numeric(0), merr = numeric(0),
         gradients = FALSE, residuals = FALSE)
  }

  testthat::with_mocked_bindings({
    getFromNamespace("npindex.default", "np")(
      bws = log(y) ~ .,
      txdat = dat[c("x1", "x2", "noise")],
      tydat = log(dat$y),
      exdat = dat[1:3, c("x1", "x2", "noise")]
    )
  },
    npindexbw = fake_bw,
    npindex = fake_fit,
    .package = "np"
  )

  expect_identical(captured$bw_xnames, c("x1", "x2", "noise"))
  expect_identical(captured$fit_xnames, c("x1", "x2", "noise"))
  expect_identical(captured$fit_exnames, c("x1", "x2", "noise"))
  expect_identical(unname(captured$bw_ydat), log(dat$y))
  expect_identical(unname(captured$fit_ydat), log(dat$y))
})
