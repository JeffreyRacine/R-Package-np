test_that("npregbw formula and default paths preserve current metadata shape", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  dat <- data.frame(
    y = rnorm(12),
    x = runif(12)
  )

  bw_formula <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )
  bw_default <- np::npregbw(
    xdat = dat["x"],
    ydat = dat$y,
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )

  expect_s3_class(bw_formula, "rbandwidth")
  expect_s3_class(bw_default, "rbandwidth")

  expect_true(all(c("call", "formula", "terms", "nobs.omit") %in% names(bw_formula)))
  expect_false("rows.omit" %in% names(bw_formula))
  expect_false("formula" %in% names(bw_default))
  expect_false("terms" %in% names(bw_default))
  expect_true(all(c("call", "rows.omit", "nobs.omit") %in% names(bw_default)))

  expect_identical(as.character(bw_formula$call[[1L]]), "npregbw.formula")
  expect_identical(as.character(bw_default$call[[1L]]), "npregbw.NULL")
  expect_identical(attr(bw_formula$terms, "term.labels"), "x")

  expect_identical(bw_formula$nobs.omit, 0L)
  expect_identical(bw_default$rows.omit, NA)
  expect_identical(bw_default$nobs.omit, 0)

  expect_identical(bw_formula$varnames$x, "x")
  expect_identical(bw_formula$varnames$y, "y")
  expect_identical(bw_default$varnames$x, "x")
  expect_identical(bw_default$varnames$y, "dat$y")
})

test_that("npregbw formula path preserves omitted-row metadata", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(
    y = c(1, 2, NA, 4, 5, 6),
    x = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )

  expect_identical(bw$rows.omit, 3L)
  expect_identical(bw$nobs.omit, 1L)
  expect_identical(bw$nobs, 5L)
})

test_that("npregbw rejects mismatched bandwidth object types with current message", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(1)
  y <- rnorm(8)
  numeric_bw <- np::npregbw(
    xdat = data.frame(x = runif(8)),
    ydat = y,
    regtype = "lc",
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    bws = 0.4
  )

  expect_error(
    np::npregbw(
      xdat = data.frame(x = factor(letters[1:8])),
      ydat = y,
      bws = numeric_bw,
      bandwidth.compute = FALSE
    ),
    "supplied bandwidths do not match 'xdat' in type"
  )
})
