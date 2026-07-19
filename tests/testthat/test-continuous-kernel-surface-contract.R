test_that("continuous kernel choices expose the supported kernel set", {
  legacy_supported <- c("gaussian", "epanechnikov", "uniform")
  beta_supported <- c(legacy_supported, "beta")
  constructors <- list(
    plbandwidth = "ckertype",
    sibandwidth = "ckertype",
    scbandwidth = "ckertype"
  )

  for (constructor in names(constructors)) {
    fun <- getFromNamespace(constructor, "npRmpi")
    defaults <- formals(fun)
    for (argument in constructors[[constructor]]) {
      expect_identical(
        eval(defaults[[argument]], envir = baseenv()),
        legacy_supported,
        info = paste(constructor, argument)
      )
    }
  }

  expect_identical(
    eval(formals(getFromNamespace("dbandwidth", "np"))$ckertype,
         envir = baseenv()),
    beta_supported
  )
  expect_identical(
    eval(formals(getFromNamespace("rbandwidth", "np"))$ckertype,
         envir = baseenv()),
    beta_supported
  )
  expect_identical(
    eval(formals(getFromNamespace("kbandwidth.numeric", "np"))$ckertype,
         envir = baseenv()),
    beta_supported
  )
  expect_identical(
    eval(formals(getFromNamespace("bandwidth", "np"))$ckertype,
         envir = baseenv()),
    beta_supported
  )
  for (constructor in c("conbandwidth", "condbandwidth")) {
    defaults <- formals(getFromNamespace(constructor, "np"))
    for (argument in c("cxkertype", "cykertype")) {
      expect_identical(
        eval(defaults[[argument]], envir = baseenv()),
        beta_supported,
        info = paste(constructor, argument)
      )
    }
  }
})

test_that("continuous kernel native codes preserve the reserved slot", {
  expect_identical(npRmpi:::CKER_GAUSS, 0)
  expect_identical(npRmpi:::CKER_EPAN, 4)
  expect_identical(npRmpi:::CKER_UNI, 8)
  expect_identical(npRmpi:::CKER_RESERVED, 9)
  expect_length(npRmpi:::int.kernels, 10L)
  expect_true(is.na(npRmpi:::int.kernels[[10L]]))
})
