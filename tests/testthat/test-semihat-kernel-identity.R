test_that("prepared regression kernels obey the caller policy", {
  
  x <- data.frame(x = seq(.05, .95, length.out = 40L))
  y <- sin(x$x * 4)
  source <- list(regtype="lp", basis="glp", degree=1L,
    bernstein.basis=TRUE, method="cv.ls", type="fixed",
    ckertype="gaussian", ckerorder=2L, ckerbound="fixed",
    ckerlb=0, ckerub=1, ynames="y")
  make <- getFromNamespace(".np_semihat_make_regbw_state", "np")
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    source$ckertype <- kernel
    out <- make(source, x, y, .25)
    expect_identical(out$ckertype, kernel)
    expect_equal(out$ckerorder, 2)
    expect_identical(out$ckerbound, "fixed")
    expect_equal(out$ckerlb, 0)
    expect_equal(out$ckerub, 1)
  }
  for (kernel in list(NULL, character())) {
    source$ckertype <- kernel
    expect_identical(make(source, x, y, .25)$ckertype, "gaussian")
  }
  source$ckertype <- "beta"
  expect_error(make(source, x, y, .25), "kernel.*not supported.*preparation")
  out <- make(source, x, y, .25, kernel.choices = c("gaussian", "beta"))
  expect_identical(out$ckertype, "beta")
  source$ckertype <- "not-a-kernel"
  expect_error(make(source, x, y, .25), "kernel.*not supported.*preparation")
  source$ckertype <- "gaussian"
  for (method in c("ichimura", "kleinspady", "manual")) {
    source$method <- method
    expect_identical(make(source, x, y, .25)$method, "cv.ls")
  }
})

test_that("bias preparation preserves main and auxiliary kernel identities", {
  
  x <- data.frame(x = seq(.05, .95, length.out = 40L))
  source <- list(regtype="lp", basis="glp", degree=1L,
    bernstein.basis=TRUE, method="cv.ls", type="fixed",
    ckertype="beta", ckerorder=2, ckerbound="fixed",
    ckerlb=0, ckerub=1, ynames="y")
  make <- getFromNamespace(".np_bias_regbw", "np")
  main <- make(source, x, .25)
  expect_identical(main$ckertype, "beta")
  expect_equal(main$ckerorder, 2)
  expect_identical(main$ckerbound, "fixed")
  auxiliary <- make(source, x, .25, degree=2L, bernstein=TRUE)
  expect_identical(auxiliary$ckertype, "gaussian")
  expect_equal(auxiliary$ckerorder, 2)
  expect_identical(auxiliary$ckerbound, "none")
  expect_null(auxiliary$call)
  for (kernel in c("epanechnikov", "uniform")) {
    retained <- source
    retained$ckertype <- kernel
    expect_identical(make(retained, x, .25)$ckertype, kernel)
  }
  original <- getFromNamespace(".np_semihat_make_regbw_state", "np")
  for (field in c("ckertype", "ckerorder", "ckerbound")) {
    testthat::with_mocked_bindings({
      expect_error(make(source, x, .25), "kernel identity")
    }, .np_semihat_make_regbw_state = function(...) {
      out <- original(...)
      out[[field]] <- switch(field, ckertype="gaussian", ckerorder=4L,
                             ckerbound="range")
      out
    }, .package="np")
  }
  source$ckerorder <- 2L
  testthat::with_mocked_bindings({
    expect_identical(make(source, x, .25)$ckerorder, 2)
  }, .np_semihat_make_regbw_state = function(...) {
    out <- original(...)
    out$ckerorder <- as.double(out$ckerorder)
    out
  }, .package="np")
})

test_that("semiparametric constructors retain beta refusals", {
  
  x <- data.frame(x = seq(.05, .95, length.out = 40L))
  z <- data.frame(z = rev(x$x))
  y <- sin(x$x * 4)
  expect_error(npindexbw(xdat=data.frame(x=x$x, x2=sin(x$x)), ydat=y, bws=c(1, .2, .25),
    bandwidth.compute=FALSE, ckertype="beta"), "arg.*should be one of")
  expect_error(npscoefbw(xdat=x, zdat=z, ydat=y, bws=.25,
    bandwidth.compute=FALSE, ckertype="beta"), "arg.*should be one of")
})

test_that("single-index preparation refuses explicit invalid kernel metadata", {
  x <- data.frame(index = seq(.05, .95, length.out = 40L))
  y <- sin(4 * x$index)
  bw <- npindexbw(xdat=x, ydat=y, bws=c(1, .25), regtype="ll",
                  bandwidth.compute=FALSE)
  make <- getFromNamespace(".np_indexhat_kbw", "np")
  baseline <- make(bw, x)
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    current <- bw
    current[["ckertype"]] <- kernel
    expect_identical(make(current, x)[["ckertype"]], kernel)
  }
  for (kernel in list(NULL, character())) {
    current <- bw
    current[["ckertype"]] <- kernel
    expect_equal(make(current, x), baseline)
  }
  for (kernel in list("not-a-kernel", "beta", NA_character_, FALSE)) {
    current <- bw
    current[["ckertype"]] <- kernel
    expect_error(make(current, x), "kernel.*not supported.*single-index preparation")
  }
  # An unrelated prefix must not shadow a missing exact kernel field.
  current <- bw
  current[["ckertype"]] <- NULL
  current[["ckertype.annotation"]] <- "not-a-kernel"
  expect_equal(make(current, x), baseline)
})
