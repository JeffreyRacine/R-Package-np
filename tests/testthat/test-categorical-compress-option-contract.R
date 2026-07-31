test_that("categorical compression accepts only an explicit logical scalar", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  options(np.categorical.compress = TRUE)
  expect_true(np:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))
  expect_false(np:::npUseCategoricalCompress(ncon = 1L, ncat = 1L))
  expect_false(np:::npUseCategoricalCompress(ncon = 0L, ncat = 0L))

  options(np.categorical.compress = FALSE)
  expect_false(np:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))

  invalid <- list(NA, logical(0L), c(TRUE, FALSE), 1, "TRUE")
  for (value in invalid) {
    options(np.categorical.compress = value)
    expect_error(
      np:::npUseCategoricalCompress(ncon = 0L, ncat = 1L),
      "option 'np.categorical.compress' must be a single non-missing logical value",
      fixed = TRUE
    )
    expect_false(np:::npUseCategoricalCompress(ncon = 1L, ncat = 1L))
    expect_false(np:::npUseCategoricalCompress(ncon = 0L, ncat = 0L))
  }
})

test_that("categorical compression defaults to enabled when unset", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  options(np.categorical.compress = NULL)
  expect_true(np:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))
})

test_that("public categorical routes reject an invalid compression option", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  dat <- factor(rep(letters[1:3], each = 8L))
  options(np.categorical.compress = TRUE)
  bws <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE
  )

  options(np.categorical.compress = "TRUE")
  expect_error(
    npudens(bws = bws, tdat = dat),
    "option 'np.categorical.compress' must be a single non-missing logical value",
    fixed = TRUE
  )
})
