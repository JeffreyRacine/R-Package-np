test_that("partial-linearity test is exported and both result methods are registered", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi", "nppltest"), ns$nppltest)
  for (generic in c("print", "summary"))
    expect_identical(getS3method(generic, "pltest", envir = ns),
                     get(paste0(generic, ".pltest"), envir = ns, inherits = FALSE))
  expect_error(nppltest(NULL, B = TRUE), "B must")
  expect_error(nppltest(NULL), "object must")
})
