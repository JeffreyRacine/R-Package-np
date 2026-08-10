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
  expect_error(gradients(npreg.obj, se = TRUE), "gradient standard errors were not computed")

  si.obj <- structure(list(grad = NA, gerr = NA), class = "singleindex")
  expect_error(gradients(si.obj), "gradients are not available")
  expect_error(gradients(si.obj, se = TRUE), "gradient standard errors were not computed")

  cd.obj <- structure(list(congrad = NA, congerr = NA), class = "condensity")
  expect_error(gradients(cd.obj), "gradients are not available")
  expect_error(gradients(cd.obj, se = TRUE), "gradient standard errors were not computed")

  cdf.obj <- structure(list(congrad = NA, congerr = NA), class = "condistribution")
  expect_error(gradients(cdf.obj), "gradients are not available")
  expect_error(gradients(cdf.obj, se = TRUE), "gradient standard errors were not computed")
})
