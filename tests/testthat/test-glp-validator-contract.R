dim_basis <- getFromNamespace("dim_basis", "npRmpi")
dimBS <- getFromNamespace("dimBS", "npRmpi")

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
  expect_warning(
    expect_identical(npRmpi:::npValidateGlpDegree("lp", 30L, 1L), 30L),
    "unusually large polynomial degree values above 25"
  )
  expect_error(
    npRmpi:::npValidateGlpDegree("lp", 101L, 1L),
    "\\[0,100\\]"
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

test_that("dim_basis validates integer-like vectors before coercion", {
  expect_error(dim_basis(degree = "a"), "degree must contain finite non-negative integers")
  expect_error(dim_basis(degree = 1L, segments = "a"), "segments must contain finite positive integers")
  expect_error(
    dim_basis(degree = integer(0), include = c("a"), categories = c(0L)),
    "include must contain finite non-negative integers"
  )
  expect_error(
    dim_basis(degree = integer(0), include = c(0L), categories = c("a")),
    "categories must contain finite non-negative integers"
  )
})

test_that("dimBS compatibility wrapper matches dim_basis", {
  expect_identical(
    dimBS(basis = "tensor", degree = c(2L, 3L), segments = c(1L, 2L)),
    dim_basis(basis = "tensor", degree = c(2L, 3L), segments = c(1L, 2L))
  )
  expect_identical(
    dimBS(basis = "glp",
          kernel = FALSE,
          degree = c(2L, 1L),
          segments = c(1L, 1L),
          include = c(1L, 0L, 1L),
          categories = c(3L, 2L, 4L)),
    dim_basis(basis = "glp",
              kernel = FALSE,
              degree = c(2L, 1L),
              segments = c(1L, 1L),
              include = c(1L, 0L, 1L),
              categories = c(3L, 2L, 4L))
  )
})
