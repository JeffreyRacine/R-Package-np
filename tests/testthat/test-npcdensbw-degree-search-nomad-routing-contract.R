test_that("npcdensbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.08)

  bw <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$best.fval, bw$degree.search$baseline.fval - 1e-8)
})
