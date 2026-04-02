test_that("npseed accepts valid integer-representable scalars", {
  expect_silent(npseed(42))
  expect_silent(npseed(-42))
  expect_silent(npseed(42.0))
})

test_that("npseed rejects invalid scalar inputs before reaching C", {
  expect_error(npseed(NA_integer_), "'seed' must be a single finite numeric value")
  expect_error(npseed(Inf), "'seed' must be a single finite numeric value")
  expect_error(npseed(NaN), "'seed' must be a single finite numeric value")
  expect_error(npseed(numeric()), "'seed' must be a single finite numeric value")
  expect_error(npseed(c(1, 2)), "'seed' must be a single finite numeric value")
  expect_error(npseed(1.5), "'seed' must be representable as a non-negative integer after abs\\(\\)")
  expect_error(npseed(.Machine$integer.max + 1), "'seed' must be representable as a non-negative integer after abs\\(\\)")
  expect_error(npseed("1"), "'seed' must be a single finite numeric value")
})

test_that("C_np_set_seed rejects direct misuse", {
  expect_error(.Call("C_np_set_seed", NA_integer_, PACKAGE = "np"), "seed must be finite")
  expect_error(.Call("C_np_set_seed", Inf, PACKAGE = "np"), "seed must be finite")
  expect_error(.Call("C_np_set_seed", 1.5, PACKAGE = "np"), "seed must be representable as a non-negative integer after abs\\(\\)")
  expect_error(.Call("C_np_set_seed", c(1L, 2L), PACKAGE = "np"), "seed must have length 1")
  expect_error(.Call("C_np_set_seed", TRUE, PACKAGE = "np"), "seed must be numeric")
})
