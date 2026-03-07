test_that("npudensbw formula and default paths preserve current metadata shape", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(7)
  dat <- data.frame(x = runif(12))

  bw_formula <- np::npudensbw(
    ~ x,
    data = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )
  bw_default <- np::npudensbw(
    dat = dat["x"],
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )

  expect_s3_class(bw_formula, "bandwidth")
  expect_s3_class(bw_default, "bandwidth")

  expect_true(all(c("call", "formula", "terms", "nobs.omit") %in% names(bw_formula)))
  expect_false("rows.omit" %in% names(bw_formula))
  expect_false("formula" %in% names(bw_default))
  expect_false("terms" %in% names(bw_default))
  expect_true(all(c("call", "rows.omit", "nobs.omit") %in% names(bw_default)))

  expect_identical(as.character(bw_formula$call[[1L]]), "npudensbw.formula")
  expect_identical(as.character(bw_default$call[[1L]]), "npudensbw.NULL")
  expect_identical(attr(bw_formula$terms, "term.labels"), "x")

  expect_identical(bw_formula$nobs.omit, 0L)
  expect_identical(bw_default$rows.omit, NA)
  expect_identical(bw_default$nobs.omit, 0)

  expect_identical(bw_formula$varnames$x, "x")
  expect_identical(bw_default$varnames$x, "x")
})

test_that("npudensbw formula path preserves omitted-row metadata", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(x = c(0.1, 0.2, NA, 0.4, 0.5, 0.6))

  bw <- np::npudensbw(
    ~ x,
    data = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )

  expect_identical(bw$rows.omit, 3L)
  expect_identical(bw$nobs.omit, 1L)
  expect_identical(bw$nobs, 5L)
})

test_that("npudensbw formula path rejects responses with current message", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(y = rnorm(8), x = runif(8))

  expect_error(
    np::npudensbw(y ~ x, data = dat, bwtype = "fixed", nmulti = 1),
    "invalid density formula"
  )
})

test_that("npudensbw rejects mismatched bandwidth object types with current message", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  numeric_bw <- np::npudensbw(
    dat = data.frame(x = runif(8)),
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    bws = 0.4
  )

  expect_error(
    np::npudensbw(
      dat = data.frame(x = factor(letters[1:8])),
      bws = numeric_bw,
      bandwidth.compute = FALSE
    ),
    "supplied bandwidths do not match 'dat' in type"
  )
})
