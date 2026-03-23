test_that("default nmulti cap follows min(2, p) on public bandwidth routes", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  dat <- data.frame(
    y = rnorm(14),
    x1 = runif(14),
    x2 = runif(14)
  )

  bw_reg_uni <- np::npregbw(
    y ~ x1,
    data = dat,
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.aic"
  )
  bw_reg_multi <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.aic"
  )
  bw_cd_multi <- np::npcdensbw(
    xdat = dat[c("x1", "x2")],
    ydat = dat["y"],
    regtype = "lc",
    bwtype = "fixed",
    bwmethod = "cv.ls"
  )

  expect_identical(length(bw_reg_uni$fval.history), 1L)
  expect_identical(length(bw_reg_multi$fval.history), 2L)
  expect_identical(length(bw_cd_multi$fval.history), 2L)
})
