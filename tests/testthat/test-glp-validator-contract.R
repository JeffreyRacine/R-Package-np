test_that("GLP validators handle non-numeric and missing inputs deterministically", {
  expect_null(npRmpi:::npValidateGlpDegree("lc", "a", 1L))
  expect_identical(npRmpi:::npValidateGlpDegree("lp", NULL, 0L), integer(0))

  expect_error(
    npRmpi:::npValidateGlpDegree("lp", "a", 1L),
    "finite non-negative integers"
  )
  expect_error(
    npRmpi:::npValidateGlpDegree("lp", c(1, NA_real_), 2L),
    "finite non-negative integers"
  )

  expect_null(npRmpi:::npValidateGlpGradientOrder("lc", "a", 1L))
  expect_identical(npRmpi:::npValidateGlpGradientOrder("lp", NULL, 0L), integer(0))

  expect_error(
    npRmpi:::npValidateGlpGradientOrder("lp", "a", 1L),
    "finite positive integers"
  )
  expect_error(
    npRmpi:::npValidateGlpGradientOrder("lp", c(1, NA_real_), 2L),
    "finite positive integers"
  )
})
