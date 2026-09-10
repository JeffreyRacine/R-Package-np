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

test_that("missing optional outputs give safe family-specific no-search refits", {
  specs <- list(
    npregression = c("npreg", "grad", "gerr"),
    condensity = c("npcdens", "congrad", "congerr"),
    condistribution = c("npcdist", "congrad", "congerr"),
    lsqregression = c("nplsqreg", "quantgrad", "quantgerr"),
    qregression = c("npqreg", "quantgrad", "quantgerr"))
  for (cls in names(specs)) {
    spec <- specs[[cls]]
    fields <- setNames(list(NA, NA), spec[2:3])
    model <- structure(fields, class = cls)
    expect_error(gradients(model), paste0(
      "Refit without repeating bandwidth search: ", spec[[1]],
      "(bws = model$bws, gradients = TRUE)."), fixed = TRUE)
    expect_error(gradients(model, se = TRUE), paste0(
      "Refit without repeating bandwidth search: ", spec[[1]],
      "(bws = model$bws, gradients = TRUE, se = TRUE)."), fixed = TRUE)
  }
  for (spec in list(c("npregression", "npreg"), c("smoothcoefficient", "npscoef"),
                   c("lsqregression", "nplsqreg"), c("condensity", "npcdens"),
                   c("condistribution", "npcdist"))) {
    model <- structure(list(se = FALSE), class = spec[[1]])
    expect_error(se(model), paste0(
      "standard errors were not computed. Refit without repeating bandwidth search: ",
      spec[[2]], "(bws = model$bws, se = TRUE)."), fixed = TRUE)
  }
  model <- structure(list(), class = "conmode")
  expect_error(gradients(model),
    "npconmode(bws = model$bws, gradients = TRUE).", fixed = TRUE)
  expect_error(gradients(model, se = TRUE),
    "gradient standard errors are not available for conmode objects", fixed = TRUE)
  expect_identical(
    .np_se_refit_hint(quote(`my fit`), "npreg", "gradients = TRUE"),
    "Refit without repeating bandwidth search: npreg(bws = `my fit`$bws, gradients = TRUE).")
  expect_match(.np_se_refit_hint(quote(models[[1]]), "npreg"),
               "Replace 'object' with your fitted model.", fixed = TRUE)
  expect_identical(.np_index_refit_hint(quote(model), gradients = TRUE, se = TRUE),
    "Refit without repeating bandwidth search: npindex(bws = model$bws, gradients = TRUE, se = TRUE).")
})
