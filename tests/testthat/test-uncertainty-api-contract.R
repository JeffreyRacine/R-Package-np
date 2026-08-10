test_that("public uncertainty controls use one canonical naming contract", {
  fit.methods <- list(
    getS3method("npreg", "formula"),
    getS3method("npreg", "default"),
    getS3method("npreg", "rbandwidth"),
    getS3method("npscoef", "formula"),
    getS3method("npscoef", "default"),
    getS3method("npscoef", "scbandwidth"),
    getS3method("npindex", "formula"),
    getS3method("npindex", "default"),
    getS3method("npindex", "sibandwidth"),
    getS3method("nplsqreg", "formula"),
    getS3method("nplsqreg", "default"),
    getS3method("nplsqreg", "lsqregressionbandwidth")
  )
  for (fun in fit.methods) {
    expect_true("se" %in% names(formals(fun)))
    expect_identical(formals(fun)$se, FALSE)
    expect_false("errors" %in% names(formals(fun)))
  }

  gradient.classes <- c(
    "npregression", "condensity", "condistribution", "lsqregression",
    "qregression", "conmode", "singleindex"
  )
  for (class in gradient.classes) {
    fun <- getS3method("gradients", class)
    expect_true("se" %in% names(formals(fun)))
    expect_identical(formals(fun)$se, FALSE)
    expect_false("errors" %in% names(formals(fun)))
  }

  predict.reg <- getS3method("predict", "npregression")
  expect_true("se.fit" %in% names(formals(predict.reg)))
  expect_identical(formals(predict.reg)$se.fit, FALSE)
  expect_true("errors" %in% names(formals(getS3method("plot", "npcopula"))))
})

test_that("legacy Boolean errors spelling fails instead of translating", {
  expect_error(npreg(errors = TRUE), "use 'se' instead", fixed = TRUE)
  expect_error(npscoef(errors = TRUE), "use 'se' instead", fixed = TRUE)
  expect_error(npindex(errors = TRUE), "use 'se' instead", fixed = TRUE)
  expect_error(nplsqreg(errors = TRUE), "use 'se' instead", fixed = TRUE)

  gradient.classes <- c(
    "npregression", "condensity", "condistribution", "lsqregression",
    "qregression", "conmode", "singleindex"
  )
  for (class in gradient.classes) {
    object <- structure(list(), class = class)
    expect_error(
      gradients(object, errors = TRUE),
      "use 'se' instead",
      fixed = TRUE,
      info = class
    )
  }

  expect_error(
    gradients(structure(list(), class = "qregression"), errors = FALSE),
    "use 'se' instead",
    fixed = TRUE
  )
  expect_error(
    coef(structure(list(), class = "plregression"), errors = TRUE),
    "use 'se' instead",
    fixed = TRUE
  )
})

test_that("valid uncertainty extraction remains controlled only by se", {
  qreg <- structure(
    list(quantgrad = 1, quantgerr = 2),
    class = "qregression"
  )
  expect_identical(gradients(qreg, se = FALSE), 1)
  expect_identical(gradients(qreg, se = TRUE), 2)

  plreg <- structure(
    list(xcoef = 3, xcoeferr = 4),
    class = "plregression"
  )
  expect_identical(coef(plreg, se = FALSE), 3)
  expect_identical(coef(plreg, se = TRUE), 4)
})

test_that("missing uncertainty state fails helpfully", {
  reg <- structure(list(se = FALSE, merr = NULL), class = "npregression")
  expect_error(se(reg), "se=TRUE", fixed = TRUE)

  sc <- structure(list(se = FALSE, merr = NA_real_), class = "smoothcoefficient")
  expect_error(se(sc), "se=TRUE", fixed = TRUE)

  si <- structure(list(merr = NULL), class = "singleindex")
  expect_error(se(si), "se=TRUE", fixed = TRUE)

  pl <- structure(list(xcoef = 1, xcoeferr = NULL, xcoefvcov = NULL),
                  class = "plregression")
  expect_error(coef(pl, se = TRUE), "not computed", fixed = TRUE)
  expect_error(vcov(pl), "not computed", fixed = TRUE)
})
