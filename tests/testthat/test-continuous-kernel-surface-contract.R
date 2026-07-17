test_that("continuous kernel choices expose the supported kernel set", {
  supported <- c("gaussian", "epanechnikov", "uniform")
  constructors <- list(
    bandwidth = "ckertype",
    dbandwidth = "ckertype",
    rbandwidth = "ckertype",
    plbandwidth = "ckertype",
    sibandwidth = "ckertype",
    scbandwidth = "ckertype",
    conbandwidth = c("cxkertype", "cykertype"),
    condbandwidth = c("cxkertype", "cykertype"),
    kbandwidth.numeric = "ckertype"
  )

  for (constructor in names(constructors)) {
    fun <- getFromNamespace(constructor, "npRmpi")
    defaults <- formals(fun)
    for (argument in constructors[[constructor]]) {
      expect_identical(
        eval(defaults[[argument]], envir = baseenv()),
        supported,
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
