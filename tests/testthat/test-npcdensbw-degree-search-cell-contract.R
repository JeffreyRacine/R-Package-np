test_that("npcdensbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(24)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw0 <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- np::npcdensbw(
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

  expect_s3_class(auto, "conbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_gte(auto$fval, max(bw0$fval, bw1$fval) - 1e-10)
  expect_gte(auto$degree.search$best.fval, auto$degree.search$baseline.fval - 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))
  expect_identical(nrow(auto$degree.search$trace), auto$degree.search$n.unique)
  expect_identical(auto$degree.search$n.cached, auto$degree.search$n.visits - auto$degree.search$n.unique)

  manual <- np::npcdensbw(
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

test_that("npcdensbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(
    x1 = runif(22),
    x2 = runif(22)
  )
  dat$y <- dat$x1 + dat$x2^2 + rnorm(nrow(dat), sd = 0.08)

  exhaustive <- np::npcdensbw(
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
  coordinate <- np::npcdensbw(
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
})

test_that("npcdensbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(y = rnorm(20), x = runif(20))

  expect_error(
    np::npcdensbw(
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

  bw <- np::npcdensbw(
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

  expect_s3_class(bw, "conbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)
})

test_that("npcdens forwards automatic LP degree search through npcdensbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(20))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  fit <- local({
    suppressPackageStartupMessages(library(np))
    npcdens(
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

  expect_false(is.null(fit$bws))
  expect_s3_class(fit$bws, "conbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})
