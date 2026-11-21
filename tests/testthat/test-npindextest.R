test_that("index test is exported with registered methods and validates arguments", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi","npindextest"), ns$npindextest)
  for (generic in c("print","summary"))
    expect_identical(getS3method(generic,"indextest",envir=ns),
                     get(paste0(generic,".indextest"),ns,inherits=FALSE))
  expect_error(npindextest(NULL,B=2L),"at least 9")
  expect_error(npindextest(NULL),"object must")
  expect_error(npindextest(NULL,bw.method="cv"),"arg")
})
