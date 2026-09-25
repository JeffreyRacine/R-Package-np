test_that("NN terminal context leaves successful and unrelated calls unchanged", {
  context <- getFromNamespace(".np_with_nn_radius_context", "npRmpi")
  expect_identical(context(42, stop("names forced"), bws = stop("bws forced"),
                           ntrain = stop("ntrain forced")), 42)
  expect_error(context(stop("unrelated"), NULL, bws = stop("bws forced")),
               "^unrelated$")
})

test_that("NN failed folds explain the existing effective-sample limit", {
  context <- getFromNamespace(".np_with_nn_radius_context", "npRmpi")
  old <- options(np.extendednn = FALSE); on.exit(options(old), add = TRUE)
  b <- list(type = "generalized_nn", icon = TRUE, bw = 29, nobs = 30L)
  original <- simpleError("C_np_regression_lp_apply_conditional: LP hat helper failed")
  e <- tryCatch(context(stop(original), "x", bws = b, ntrain = 30L,
                         leave.one.out = TRUE), error = identity)
  expect_identical(class(e), class(original))
  expect_match(conditionMessage(e), "effective training size 29", fixed = TRUE)
  expect_match(conditionMessage(e), "bandwidth exceeds n-1", fixed = TRUE)
  unchanged <- tryCatch(context(stop(original), "x", bws = b, ntrain = 30L),
                        error = identity)
  expect_identical(unchanged, original)
})
