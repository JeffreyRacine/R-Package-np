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
    getFromNamespace("npindex.default", "npRmpi")(
      bws = y ~ x1 + x2,
      txdat = dat[c("x1", "x2", "noise")],
      tydat = dat$y,
      exdat = dat[1:3, c("x1", "x2", "noise")]
    )
  },
    .npRmpi_require_active_slave_pool = function(...) invisible(TRUE),
    npindexbw = fake_bw,
    npindex = fake_fit,
    .package = "npRmpi"
  )

  expect_identical(captured$bw_xnames, c("x1", "x2"))
  expect_identical(captured$fit_xnames, c("x1", "x2"))
  expect_identical(captured$fit_exnames, c("x1", "x2"))
  expect_identical(unname(captured$bw_ydat), dat$y)
  expect_identical(unname(captured$fit_ydat), dat$y)
})
