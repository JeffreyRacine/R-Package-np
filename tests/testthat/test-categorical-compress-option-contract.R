test_that("categorical compression accepts only an explicit logical scalar", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  options(np.categorical.compress = TRUE)
  expect_true(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))
  expect_false(npRmpi:::npUseCategoricalCompress(ncon = 1L, ncat = 1L))
  expect_false(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 0L))

  options(np.categorical.compress = FALSE)
  expect_false(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))

  invalid <- list(NA, logical(0L), c(TRUE, FALSE), 1, "TRUE")
  for (value in invalid) {
    options(np.categorical.compress = value)
    expect_error(
      npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 1L),
      "option 'np.categorical.compress' must be a single non-missing logical value",
      fixed = TRUE
    )
    expect_false(npRmpi:::npUseCategoricalCompress(ncon = 1L, ncat = 1L))
    expect_false(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 0L))
  }
})

test_that("categorical compression defaults to enabled when unset", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  options(np.categorical.compress = NULL)
  expect_true(npRmpi:::npUseCategoricalCompress(ncon = 0L, ncat = 1L))
})

test_that("categorical route resolution rejects an invalid compression option", {
  old <- getOption("np.categorical.compress")
  on.exit(options(np.categorical.compress = old), add = TRUE)

  options(np.categorical.compress = "TRUE")
  expect_error(
    npRmpi:::npDoTreeOrCategoricalCompress(ncon = 0L, ncat = 1L),
    "option 'np.categorical.compress' must be a single non-missing logical value",
    fixed = TRUE
  )
})
