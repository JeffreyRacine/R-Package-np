test_that("public single-index external rows do not depend on inference requests", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 24L))
  y <- 1 + sin(x$x)
  bw <- suppressWarnings(npindexbw(xdat = x, ydat = y, bws = c(1, .5),
    ckertype = "epanechnikov", bandwidth.compute = FALSE))
  ex <- data.frame(x = c(0, .25, 10))
  capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(force(expr), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
  }
  off <- capture(npindex(bws = bw, txdat = x, tydat = y, exdat = ex, se = FALSE))
  on <- capture(npindex(bws = bw, txdat = x, tydat = y, exdat = ex))
  grad <- capture(npindex(bws = bw, txdat = x, tydat = y, exdat = ex,
    se = FALSE, gradients = TRUE))
  expect_equal(off$value$mean, on$value$mean, tolerance = 1e-12)
  expect_equal(off$value$mean, grad$value$mean, tolerance = 1e-12)
  expect_true(is.na(off$value$mean[3L]))
  for (result in list(off, on, grad))
    expect_length(result$warnings, 1L)
  set.seed(18)
  boot <- capture(npindex(bws = bw, txdat = x, tydat = y, exdat = ex,
    se = TRUE, se.type = "bootstrap", B = 1L))
  expect_true(is.na(boot$value$merr[3L]))
  expect_length(boot$warnings, 1L)
})

test_that("index normalization is lazy and does not disguise signed cancellation", {
  tww <- array(1e-300, c(2L, 2L, 2L))
  tww[1L, 2L, ] <- c(2e-300, 3e-300)
  expect_identical(.np_index_normalized_mean(tww, stop("unused"),
    stop("unused"), stop("unused"), TRUE), c(2, 3))
  zero <- array(0, c(2L, 2L, 1L))
  bw <- list(bw = 1, type = "fixed", ckertype = "gaussian",
             ckerorder = 4L, ckerbound = "none")
  x <- data.frame(index = c(-1, 1))
  testthat::local_mocked_bindings(
    .np_index_kernel_sum = function(...) list(kw = matrix(c(1, -1), ncol = 1L)),
    .package = "np")
  expect_error(.np_index_normalized_mean(zero, x, x[1, , drop = FALSE], bw, TRUE),
    "required ratio is undefined", fixed = TRUE)
  expect_error(.np_index_normalized_mean(zero, x, x[1, , drop = FALSE], bw),
    "required ratio is undefined", fixed = TRUE)
  healthy <- list(mean = c(1, 2))
  expect_identical(.np_index_fit_rows(healthy, stop("unused"), stop("unused"),
    stop("unused"), TRUE), healthy)
})

test_that("index bootstrap draw support does not overwrite the original point fit", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(seq(-.5, .5, length.out = 11L), 4))
  y <- 1 + x$x + .2*sin(seq_len(nrow(x)))
  set.seed(18)
  plan <- boot::boot(data.frame(x, y), function(data, indices) 0, R = 4L)
  draws <- boot::boot.array(plan, indices = TRUE)
  expected.rng <- .Random.seed
  expect_true(any(rowSums(draws == nrow(x)) == 0L))
  for (reg in c("lc", "ll", "lp")) for (grad in c(FALSE, TRUE)) {
    bw <- suppressWarnings(do.call(npindexbw, c(list(xdat = x, ydat = y,
      bws = c(1, .8), regtype = reg, bandwidth.compute = FALSE,
      ckertype = "epanechnikov"), if (reg == "lp") list(degree = 2L) else list())))
    point <- npindex(bws = bw, txdat = x, tydat = y, se = FALSE, gradients = grad)
    set.seed(18)
    warnings <- character()
    fit <- withCallingHandlers(npindex(bws = bw, txdat = x, tydat = y,
      se = TRUE, se.type = "bootstrap", B = 4L, gradients = grad),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    expect_identical(fit$mean, point$mean)
    expect_true(all(is.finite(fit$mean)))
    expect_true(is.na(fit$merr[12L]))
    expect_true(all(is.finite(fit$merr[1:11])))
    expect_length(warnings, 1L)
    expect_identical(.Random.seed, expected.rng)
    if (grad) expect_true(all(is.na(fit$gerr[12L, ])))
  }
})
