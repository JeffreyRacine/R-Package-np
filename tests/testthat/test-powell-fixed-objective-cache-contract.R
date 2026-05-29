fixed_cache_regression_data <- function(kind, n = 80L) {
  set.seed(101)
  if (identical(kind, "continuous")) {
    dat <- data.frame(x1 = runif(n), x2 = runif(n))
    dat$y <- sin(dat$x1) + dat$x2 + rnorm(n, sd = 0.2)
    return(dat)
  }
  if (identical(kind, "mixed")) {
    dat <- data.frame(x1 = runif(n), xo = ordered(sample(1:3, n, TRUE)))
    dat$y <- dat$x1 + as.integer(dat$xo) / 3 + rnorm(n, sd = 0.2)
    return(dat)
  }
  dat <- data.frame(
    xo = ordered(sample(1:4, n, TRUE)),
    xu = factor(sample(letters[1:3], n, TRUE))
  )
  dat$y <- as.integer(dat$xo) + as.integer(dat$xu) + rnorm(n, sd = 0.2)
  dat
}

fixed_cache_formula <- function(kind) {
  if (identical(kind, "continuous"))
    return(y ~ x1 + x2)
  if (identical(kind, "mixed"))
    return(y ~ x1 + xo)
  y ~ xo + xu
}

run_fixed_cache_bw <- function(kind, cache) {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.objective.cache = cache)
  on.exit(options(old), add = TRUE)
  dat <- fixed_cache_regression_data(kind)
  np::npregbw(
    fixed_cache_formula(kind),
    data = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lc",
    nmulti = 1
  )
}

expect_fixed_cache_parity <- function(kind) {
  cached <- run_fixed_cache_bw(kind, TRUE)
  uncached <- run_fixed_cache_bw(kind, FALSE)

  expect_equal(cached$bw, uncached$bw, tolerance = 0)
  expect_equal(cached$fval, uncached$fval, tolerance = 1e-12)
  expect_equal(cached$num.feval, uncached$num.feval)

  expect_gt(as.numeric(cached$num.feval.fast[1L]), 0)
  expect_equal(as.numeric(uncached$num.feval.fast[1L]), 0)
  expect_equal(unname(cached$nn.cache[["objective.enabled"]]), 1)
  expect_gt(unname(cached$nn.cache[["objective.hits"]]), 0)
  expect_equal(
    as.numeric(cached$num.feval.fast[1L]),
    unname(cached$nn.cache[["objective.hits"]])
  )
  expect_equal(unname(uncached$nn.cache[["objective.enabled"]]), 0)
  expect_equal(unname(uncached$nn.cache[["objective.hits"]]), 0)
}

test_that("np.objective.cache preserves fixed continuous Powell regression parity", {
  expect_fixed_cache_parity("continuous")
})

test_that("np.objective.cache preserves fixed mixed Powell regression parity", {
  expect_fixed_cache_parity("mixed")
})

test_that("np.objective.cache preserves fixed categorical Powell regression parity", {
  expect_fixed_cache_parity("categorical")
})
