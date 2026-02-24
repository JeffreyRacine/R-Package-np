test_that("GLP validators handle non-numeric and missing inputs deterministically", {
  expect_null(np:::npValidateGlpDegree("lc", "a", 1L))
  expect_identical(np:::npValidateGlpDegree("lp", NULL, 0L), integer(0))

  expect_error(
    np:::npValidateGlpDegree("lp", "a", 1L),
    "finite non-negative integers"
  )
  expect_error(
    np:::npValidateGlpDegree("lp", c(1, NA_real_), 2L),
    "finite non-negative integers"
  )

  expect_null(np:::npValidateGlpGradientOrder("lc", "a", 1L))
  expect_identical(np:::npValidateGlpGradientOrder("lp", NULL, 0L), integer(0))

  expect_error(
    np:::npValidateGlpGradientOrder("lp", "a", 1L),
    "finite positive integers"
  )
  expect_error(
    np:::npValidateGlpGradientOrder("lp", c(1, NA_real_), 2L),
    "finite positive integers"
  )
})
