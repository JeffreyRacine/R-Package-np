test_that("fixed uniform iterated correction remains unavailable at all degrees", {
  b <- list(type="fixed", ckertype="uniform", ckerorder=2L, ckerbound="none",
            ncon=1L, regtype="lc", basis="glp", degree=0L)
  for (bound in c("none", "fixed")) for (reg in c("lc", "ll", "lp")) {
    b$ckerbound <- bound
    b$regtype <- reg
    for (degree in list(0L, 1L, 2L, c(0L,1L))) {
      b$degree <- degree
      for (cal in c("absolute","separate-tails","tail-envelope","basic-envelope","basic-local-envelope"))
        expect_identical(.np_bias_regression_kernel_class(b,cal), "unsupported")
    }
  }
})

test_that("fixed rows with saved precision failures preserve live explicit scope", {
  b <- list(type="fixed", regtype="lp", basis="glp", degree=2L, ncon=1L,
            ckertype="epanechnikov", ckerorder=2L, ckerbound="none",
            bernstein.basis=FALSE)
  resolve <- function(b, calibration, gradients=FALSE)
    .np_bias_regression_calibration(b, calibration, "bootstrap", "simultaneous",
                                    "inid", "mammen", gradients)
  for (g in c(FALSE, TRUE)) {
    expect_error(resolve(b, "default", g), "explicit basic-local-envelope")
    expect_identical(resolve(b, "basic-local-envelope", g), "basic-local-envelope")
  }
  b$bernstein.basis <- TRUE
  expect_identical(resolve(b, "default"), "basic-local-envelope")
  b$ckertype <- "beta"
  b$ckerbound <- "fixed"
  expect_identical(resolve(b, "default"), "basic-local-envelope")
  expect_error(resolve(b, "absolute", TRUE), "explicit basic-local-envelope")
  expect_identical(resolve(b, "basic-local-envelope", TRUE), "basic-local-envelope")
  b$degree <- 0L
  expect_identical(resolve(b, "default", TRUE), "basic-local-envelope")
})
