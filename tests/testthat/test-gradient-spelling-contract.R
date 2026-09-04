test_that("retired gradient spelling fails without forcing its value", {
  pkg <- environmentName(environment(npreg))
  for (route in c("npreg", "npcdens", "npcdist", "npindex", "nplsqreg")) {
    fun <- getFromNamespace(route, pkg)
    expect_error(fun(gradient_order = stop("retired value forced")),
                 "unused argument 'gradient_order'; did you mean 'gradient.order'?",
                 fixed = TRUE, info = route)
    expect_error(fun(gradient.order = 1L, gradient_order = NULL),
                 "unused argument 'gradient_order'", fixed = TRUE, info = route)
  }
})

test_that("gradient accessors reject only the retired spelling at entry", {
  pkg <- environmentName(environment(npreg))
  supported <- c("npregression", "condensity", "condistribution", "singleindex",
                 "lsqregression")
  unsupported <- c("qregression", "conmode", "npregiv", "npregivderiv")
  for (cls in c(supported, unsupported)) {
    fun <- getFromNamespace(paste0("gradients.", cls), pkg)
    for (value in list(NULL, 1L, quote(stop("retired value forced")))) {
      call <- as.call(c(list(fun, x = quote(stop("object forced"))),
                        list(gradient_order = value)))
      err <- tryCatch(eval(call), error = conditionMessage)
      expect_match(err, "unused argument 'gradient_order'", fixed = TRUE, info = cls)
      expect_identical(grepl("did you mean", err, fixed = TRUE), cls %in% supported)
    }
    expect_error(fun(stop("object forced"), gradient.order = 1L,
                     gradient_order = 1L), "unused argument 'gradient_order'",
                 fixed = TRUE, info = cls)
  }
})

test_that("plot spelling checks preserve capabilities and canonical values", {
  pkg <- environmentName(environment(npreg))
  validate <- getFromNamespace(".np_plot_validate_public_dots", pkg)
  engine <- getFromNamespace(".np_plot_engine_for_bws", pkg)
  normalize <- getFromNamespace(".np_plot_normalize_public_dots", pkg)
  supported <- c("rbandwidth", "conbandwidth", "condbandwidth", "sibandwidth")
  for (cls in c(supported, "scbandwidth", "plbandwidth", "bandwidth", "dbandwidth")) {
    bws <- structure(list(), class = cls)
    err <- tryCatch(validate(alist(gradient_order = stop("value forced")),
                             method = engine(bws), bws = bws), error = conditionMessage)
    expect_match(err, "unused argument 'gradient_order'", fixed = TRUE, info = cls)
    expect_identical(grepl("did you mean", err, fixed = TRUE), cls %in% supported)
    if (cls %in% supported)
      expect_invisible(validate(alist(gradient.order = 1L), method = engine(bws)))
  }
  for (value in list(NULL, 1.5, Inf, NA_real_, c(2L, 1L)))
    expect_identical(normalize(list(gradient.order = value)), list(gradient.order = value))
  err <- tryCatch(validate(alist(gradient_order = 1L),
                           method = engine(structure(list(), class = "condbandwidth")),
                           context = "plot.qregression"), error = conditionMessage)
  expect_false(grepl("did you mean", err, fixed = TRUE))
})

test_that("unsupported quantile orders and single-index limits are unchanged", {
  pkg <- environmentName(environment(npreg))
  reject <- getFromNamespace(".npqreg_reject_gradient_order_dots", pkg)
  for (name in c("gradient.order", "gradient_order")) {
    err <- tryCatch(reject(setNames(list(1L), name)), error = conditionMessage)
    expect_match(err, name, fixed = TRUE)
    expect_false(grepl("did you mean", err, fixed = TRUE))
  }
  check <- getFromNamespace(".np_singleindex_reject_higher_gradient_order", pkg)
  expect_invisible(check(list(gradient.order = 1L)))
  expect_error(check(list(gradient.order = 2L)), "only first-order", fixed = TRUE)
})
