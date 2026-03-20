test_that("npcdensbw bounded cv.ls nomad search prefers valid positive-score solutions", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  dat <- data.frame(x = runif(350L))
  dat$y <- rbeta(nrow(dat), 1, 1)

  bw <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.min = 0L,
    degree.max = 3L,
    degree.verify = FALSE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 3L,
    cxkerbound = "range",
    cykerbound = "range"
  )

  expect_true(is.finite(bw$fval))
  expect_gt(bw$fval, 0)
  expect_equal(
    bw$fval,
    np:::.npcdensbw_eval_only(data.frame(x = dat$x), data.frame(y = dat$y), bw)$objective,
    tolerance = 1e-12
  )
})
