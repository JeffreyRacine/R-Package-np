test_that("mirrored nmulti helper returns the new default cap", {
  expect_identical(npRmpi:::npDefaultNmulti(1L), 1L)
  expect_identical(npRmpi:::npDefaultNmulti(2L), 2L)
  expect_identical(npRmpi:::npDefaultNmulti(6L), 2L)
})
