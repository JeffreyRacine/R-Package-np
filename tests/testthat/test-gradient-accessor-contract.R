test_that("gradient accessors fail fast when gradients are unavailable", {
  npreg.obj <- structure(
    list(
      grad = NA,
      gerr = NA,
      bws = list(regtype = "ll", ncon = 1L, icon = TRUE, degree = 1L)
    ),
    class = "npregression"
  )
  expect_error(gradients(npreg.obj), "gradients are not available")
  expect_error(gradients(npreg.obj, errors = TRUE), "gradient standard errors are not available")

  si.obj <- structure(list(grad = NA, gerr = NA), class = "singleindex")
  expect_error(gradients(si.obj), "gradients are not available")
  expect_error(gradients(si.obj, errors = TRUE), "gradient standard errors are not available")

  cd.obj <- structure(list(congrad = NA, congerr = NA), class = "condensity")
  expect_error(gradients(cd.obj), "gradients are not available")
  expect_error(gradients(cd.obj, errors = TRUE), "gradient standard errors are not available")

  cdf.obj <- structure(list(congrad = NA, congerr = NA), class = "condistribution")
  expect_error(gradients(cdf.obj), "gradients are not available")
  expect_error(gradients(cdf.obj, errors = TRUE), "gradient standard errors are not available")
})
