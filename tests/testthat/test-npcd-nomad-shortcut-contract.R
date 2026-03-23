test_that("npcdens nomad shortcut matches the explicit density preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  dat <- data.frame(x = sort(runif(14)), y = sort(runif(14)))

  bw_short <- np::npcdensbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )

  bw_long <- np::npcdensbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  fit_short <- local({
    npcdensbw <- np::npcdensbw
    np::npcdens(
      y ~ x,
      data = dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    )
  })

  fit_long <- local({
    npcdensbw <- np::npcdensbw
    np::npcdens(
      y ~ x,
      data = dat,
      regtype = "lp",
      search.engine = "nomad+powell",
      degree.select = "coordinate",
      bernstein.basis = TRUE,
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      nmulti = 1L
    )
  })

  expect_lt(max(abs(unname(bw_short$xbw) - unname(bw_long$xbw))), 5e-4)
  expect_lt(max(abs(unname(bw_short$ybw) - unname(bw_long$ybw))), 5e-4)
  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_true(isTRUE(bw_short$bernstein.basis))
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 5e-4)
  expect_true(is.list(fit_short$bws$nomad.shortcut))
})

test_that("npcdist nomad shortcut matches the explicit distribution preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  dat <- data.frame(x = sort(runif(14)), y = sort(runif(14)))

  bw_short <- np::npcdistbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )

  bw_long <- np::npcdistbw(
    y ~ x,
    data = dat,
    regtype = "lp",
    search.engine = "nomad+powell",
    degree.select = "coordinate",
    bernstein.basis = TRUE,
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = FALSE,
    bwtype = "fixed",
    nmulti = 1L
  )

  fit_short <- local({
    npcdistbw <- np::npcdistbw
    np::npcdist(
      y ~ x,
      data = dat,
      nomad = TRUE,
      degree.max = 1L,
      nmulti = 1L
    )
  })

  fit_long <- local({
    npcdistbw <- np::npcdistbw
    np::npcdist(
      y ~ x,
      data = dat,
      regtype = "lp",
      search.engine = "nomad+powell",
      degree.select = "coordinate",
      bernstein.basis = TRUE,
      degree.min = 0L,
      degree.max = 1L,
      degree.verify = FALSE,
      bwtype = "fixed",
      nmulti = 1L
    )
  })

  expect_equal(unname(bw_short$xbw), unname(bw_long$xbw), tolerance = 1e-8)
  expect_equal(unname(bw_short$ybw), unname(bw_long$ybw), tolerance = 1e-8)
  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_true(isTRUE(bw_short$bernstein.basis))
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 1e-8)
  expect_true(is.list(fit_short$bws$nomad.shortcut))
})

test_that("npcdens and npcdist accept nomad.nmulti through the NOMAD shortcut route", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  dat <- data.frame(x = sort(runif(14)), y = sort(runif(14)))

  bw_cd <- np::npcdensbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )
  bw_cdf <- np::npcdistbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )

  expect_true(is.list(bw_cd$nomad.shortcut))
  expect_true(is.list(bw_cdf$nomad.shortcut))
  expect_true(isTRUE(bw_cd$bernstein.basis))
  expect_true(isTRUE(bw_cdf$bernstein.basis))
})

test_that("npcdens and npcdist nomad shortcuts fail fast on incompatible settings", {
  x <- data.frame(x = seq(0, 1, length.out = 8))
  y <- data.frame(y = seq(0, 1, length.out = 8))

  expect_error(
    np::npcdensbw(xdat = x, ydat = y, nomad = TRUE, regtype = "ll"),
    "nomad=TRUE requires regtype='lp'"
  )
  expect_error(
    np::npcdensbw(xdat = x, ydat = y, nomad = TRUE, bernstein.basis = FALSE),
    "requires bernstein.basis=TRUE"
  )
  expect_error(
    np::npcdistbw(xdat = x, ydat = y, nomad = TRUE, bwtype = "adaptive_nn"),
    "requires bwtype='fixed'"
  )
  expect_error(
    np::npcdistbw(xdat = x, ydat = y, nomad = TRUE, search.engine = "cell"),
    "requires search.engine='nomad' or 'nomad\\+powell'"
  )
  expect_error(
    np::npcdensbw(xdat = x, ydat = y, nomad.nmulti = 1L),
    "nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad\\+powell'"
  )
})
