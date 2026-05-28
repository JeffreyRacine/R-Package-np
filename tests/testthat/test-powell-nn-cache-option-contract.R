test_that("np.objective.cache controls continuous NN Powell caching", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.objective.cache = TRUE)
  on.exit(options(old), add = TRUE)

  run_bw <- function(cache) {
    set.seed(42)
    dat <- data.frame(x1 = runif(80), x2 = runif(80))
    dat$y <- dat$x1 + dat$x2 + rnorm(80)
    options(np.objective.cache = cache)
    np::npregbw(
      y ~ x1 + x2,
      data = dat,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      regtype = "lc",
      nmulti = 1
    )
  }

  cached <- run_bw(TRUE)
  uncached <- run_bw(FALSE)

  expect_equal(cached$bw, uncached$bw)
  expect_equal(cached$fval, uncached$fval, tolerance = 1e-12)

  expect_equal(unname(cached$nn.cache[["enabled"]]), 1)
  expect_gt(unname(cached$nn.cache[["hits"]]), 0)
  expect_gte(cached$num.feval.fast, unname(cached$nn.cache[["hits"]]))

  expect_equal(unname(uncached$nn.cache[["enabled"]]), 0)
  expect_equal(unname(uncached$nn.cache[["hits"]]), 0)
})

test_that("np.objective.cache default ignores legacy environment off switch", {
  old <- Sys.getenv("NP_NN_POWELL_CACHE_INSTRUMENT", unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("NP_NN_POWELL_CACHE_INSTRUMENT")
    } else {
      Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = old)
    }
  }, add = TRUE)

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "off")
  expect_true(np:::npObjectiveCacheEnabled())

  Sys.setenv(NP_NN_POWELL_CACHE_INSTRUMENT = "on")
  expect_true(np:::npObjectiveCacheEnabled())
})
