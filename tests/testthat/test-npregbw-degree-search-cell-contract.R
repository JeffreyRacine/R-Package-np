test_that("npregbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(28)))
  dat$y <- dat$x^2 + rnorm(nrow(dat), sd = 0.05)

  bw0 <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expected_fval <- min(bw0$fval, bw1$fval)

  expect_s3_class(auto, "rbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, expected_fval + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))
  expect_identical(nrow(auto$degree.search$trace), auto$degree.search$n.unique)
  expect_identical(auto$degree.search$n.cached, auto$degree.search$n.visits - auto$degree.search$n.unique)

  manual <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npregbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(26),
    x2 = runif(26)
  )
  dat$y <- dat$x1 + dat$x2^2 + rnorm(nrow(dat), sd = 0.05)

  exhaustive <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  coordinate <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = TRUE,
    degree.restarts = 1L,
    degree.max.cycles = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(coordinate$degree.search$mode, "coordinate")
  expect_true(isTRUE(coordinate$degree.search$completed))
  expect_true(isTRUE(coordinate$degree.search$certified))
  expect_equal(as.integer(coordinate$degree), as.integer(exhaustive$degree))
  expect_equal(coordinate$fval, exhaustive$fval, tolerance = 1e-10)
  expect_lte(coordinate$degree.search$best.fval, coordinate$degree.search$baseline.fval + 1e-10)
  expect_identical(nrow(coordinate$degree.search$trace), coordinate$degree.search$n.unique)
  expect_identical(coordinate$degree.search$n.cached, coordinate$degree.search$n.visits - coordinate$degree.search$n.unique)
})

test_that("npregbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(y = rnorm(24), x = runif(24))

  expect_error(
    np::npregbw(
      y ~ x,
      data = dat,
      regtype = "lc",
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "rbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_false(isTRUE(bw$bernstein.basis))
})

test_that("npreg forwards automatic LP degree search through npregbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(24))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  fit <- local({
    suppressPackageStartupMessages(library(np))
    npreg(
      y ~ x,
      data = dat,
      regtype = "lp",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    )
  })

  expect_s3_class(fit, "npregression")
  expect_s3_class(fit$bws, "rbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})
