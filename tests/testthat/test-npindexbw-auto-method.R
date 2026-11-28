test_that("automatic index searches use the expanded index width", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2104)
  x <- as.data.frame(matrix(rnorm(80L * 8L), 80L))
  y <- sin(x[[1L]] + .4 * x[[2L]]) + rnorm(80L, sd = .25)
  original <- optim
  calls <- character()
  testthat::local_mocked_bindings(optim = function(par, fn, ..., method) {
    calls <<- c(calls, method)
    original(par, fn, ..., method = method)
  })
  check <- function(xx, response, args, expected) {
    opts <- c(list(xdat = xx, ydat = response, nmulti = 1L), args)
    calls <<- character()
    auto <- do.call(npindexbw, opts)
    expect_true(length(calls) > 0L && all(calls == expected))
    opts$optim.method <- expected
    calls <<- character()
    reference <- do.call(npindexbw, opts)
    fields <- c("pomethod", "beta", "bw", "fval", "fval.vector", "num.feval")
    expect_identical(auto[fields], reference[fields])
    expect_identical(auto$pomethod, expected)
    expect_true(length(calls) > 0L && all(calls == expected))
  }
  check(x[1:7], y, list(), "Nelder-Mead")
  check(x, as.numeric(y > 0), list(method = "kleinspady"), "BFGS")
  check(x, y, list(optim.threshold = 9L), "Nelder-Mead")
  check(x[1:2], y, list(optim.method = "CG", optim.threshold = 1L), "CG")
  f <- data.frame(a = x[[1L]], g = factor(rep(letters[1:4], 20L)))
  check(f, y, list(optim.threshold = 4L), "BFGS")
})

test_that("invalid index optimizer thresholds have an informative error", {
  for (value in list(NULL, NA_integer_, 0L, -1L, 1.5, Inf, "8", TRUE,
                     c(7L, 8L), matrix(8L), .Machine$integer.max + 1))
    expect_error(npindexbw(xdat = data.frame(x = 1:8), ydat = 1:8,
                           bandwidth.compute = FALSE, optim.threshold = value),
                 "optim.threshold.*positive integer")
})
