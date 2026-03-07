test_that("npcdensbw formula and default paths preserve current metadata shape", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(11)
  dat <- data.frame(
    y = rnorm(12),
    x = runif(12)
  )

  bw_formula <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )
  bw_default <- np::npcdensbw(
    xdat = dat["x"],
    ydat = dat["y"],
    regtype = "lc",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    nmulti = 1
  )

  expect_s3_class(bw_formula, "conbandwidth")
  expect_s3_class(bw_default, "conbandwidth")

  expect_true(all(c("call", "formula", "terms", "variableNames", "nobs.omit") %in% names(bw_formula)))
  expect_false("rows.omit" %in% names(bw_formula))
  expect_false("formula" %in% names(bw_default))
  expect_false("terms" %in% names(bw_default))
  expect_false("variableNames" %in% names(bw_default))
  expect_true(all(c("call", "rows.omit", "nobs.omit") %in% names(bw_default)))

  expect_identical(as.character(bw_formula$call[[1L]]), "npcdensbw.formula")
  expect_identical(as.character(bw_default$call[[1L]]), "npcdensbw.NULL")
  expect_identical(attr(bw_formula$terms, "term.labels"), c("y", "x"))
  expect_identical(bw_formula$variableNames$response, "y")
  expect_identical(bw_formula$variableNames$terms, "x")

  expect_identical(bw_formula$nobs.omit, 0L)
  expect_identical(bw_default$rows.omit, NA)
  expect_identical(bw_default$nobs.omit, 0)

  expect_identical(bw_formula$varnames$x, "x")
  expect_identical(bw_formula$varnames$y, "y")
  expect_identical(bw_default$varnames$x, "x")
  expect_identical(bw_default$varnames$y, "y")
})

test_that("npcdensbw formula path preserves omitted-row metadata", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(
    y = c(1, 2, NA, 4, 5, 6),
    x = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  )

  bw <- np::npcdensbw(
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

test_that("npcdensbw rejects mismatched x bandwidth object types with current message", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  y <- rnorm(8)
  x_numeric <- data.frame(x = runif(8))
  x_factor <- data.frame(x = factor(letters[1:8]))

  bw_numeric <- np::npcdensbw(
    xdat = x_numeric,
    ydat = data.frame(y = y),
    regtype = "lc",
    bwtype = "fixed",
    bandwidth.compute = FALSE,
    bws = c(0.4, 0.4)
  )

  expect_error(
    np::npcdensbw(
      xdat = x_factor,
      ydat = data.frame(y = y),
      bws = bw_numeric,
      bandwidth.compute = FALSE
    ),
    "supplied bandwidths do not match 'xdat' in type"
  )
})
