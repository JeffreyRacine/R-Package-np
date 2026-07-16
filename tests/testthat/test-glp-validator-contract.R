dim_basis <- getFromNamespace("dim_basis", "npRmpi")
dimBS <- getFromNamespace("dimBS", "npRmpi")

test_that("GLP validators handle non-numeric and missing inputs deterministically", {
  expect_null(npRmpi:::npValidateGlpDegree("lc", "a", 1L))
  expect_error(
    npRmpi:::npValidateGlpDegree("lp", NULL, 0L),
    "degree must be 0 when regtype='lp' has no continuous predictors"
  )
  expect_identical(npRmpi:::npValidateGlpDegree("lp", 0L, 0L), integer(0))
  expect_identical(npRmpi:::npValidateGlpDegree("lp", integer(0), 0L), integer(0))

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

test_that("dim_basis rejects capacity overflow before native integer arithmetic", {
  imax <- .Machine$integer.max

  expect_identical(
    dim_basis(basis = "additive",
              kernel = TRUE,
              degree = 1L,
              segments = 1L,
              include = 2L,
              categories = as.integer(1073741824)),
    1
  )
  expect_error(
    dim_basis(basis = "additive", degree = imax - 1L, segments = 2L),
    "basis dimension exceeds supported capacity"
  )
  expect_error(
    dim_basis(basis = "additive", degree = imax, segments = 2L),
    "basis dimension exceeds supported capacity"
  )
  expect_error(
    dim_basis(basis = "additive",
              kernel = FALSE,
              degree = integer(0),
              include = 2L,
              categories = as.integer(1073741824)),
    "basis dimension exceeds supported capacity"
  )
})

test_that("dim_basis rejects dimensions beyond supported capacity", {
  expect_error(
    dim_basis(basis = "tensor",
              degree = c(46340L, 46340L),
              segments = c(1L, 1L)),
    "basis dimension exceeds supported capacity"
  )
  expect_error(
    dim_basis(basis = "glp",
              degree = c(65535L, 65535L),
              segments = c(1L, 1L)),
    "LP basis dimension exceeds supported term capacity \\(100000\\)"
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
