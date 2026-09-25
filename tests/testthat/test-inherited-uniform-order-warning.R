test_that("retained uniform orders preserve metadata without a fresh advisory", {
  ns <- asNamespace("npRmpi")
  convert <- get("kbandwidth.default", ns)
  info <- get("untangle", ns)(data.frame(x = seq(.1, .9, length.out = 10)))
  b <- list(bandwidth = list(x = .3), type = "fixed", ckertype = "uniform",
            ckerorder = 4L, ckerbound = "none", ckerlb = -Inf, ckerub = Inf,
            ukertype = "aitchisonaitken", okertype = "liracine",
            nobs = 10L, xdati = info, ydati = NULL, xnames = "x", ynames = NULL)
  expect_warning(value <- convert(b), NA)
  expect_identical(value$ckerorder, 4L)
  expect_identical(value$bw, c(x = .3))
  expect_warning(testthat::with_mocked_bindings(convert(b),
    kbandwidth.numeric = function(...) { warning("unrelated warning"); 1 },
    .package = "npRmpi"), "unrelated warning", fixed = TRUE)
  # A duplicate explicit order still fails as before; conversion is not a retry.
  expect_error(convert(b, ckerorder = 2L), "matched by multiple actual arguments")
})
