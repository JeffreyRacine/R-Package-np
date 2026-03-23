np_prepare_nomad_shortcut <- getFromNamespace(".np_prepare_nomad_shortcut", "np")
np_bw_call_uses_nomad_degree_search <- getFromNamespace(".np_bw_call_uses_nomad_degree_search", "np")

test_that(".np_prepare_nomad_shortcut fills only missing preset arguments", {
  out <- np_prepare_nomad_shortcut(
    nomad = TRUE,
    call_names = c("xdat", "ydat", "nomad", "degree.max"),
    preset = list(
      regtype = "lp",
      search.engine = "nomad+powell",
      degree.select = "coordinate",
      bernstein.basis = TRUE,
      degree.min = 0L,
      degree.max = 10L,
      degree.verify = FALSE,
      bwtype = "fixed"
    ),
    values = list(
      regtype = NULL,
      search.engine = NULL,
      degree.select = NULL,
      bernstein.basis = NULL,
      degree.min = NULL,
      degree.max = 1L,
      degree.verify = NULL,
      bwtype = NULL
    ),
    where = "npregbw"
  )

  expect_true(out$enabled)
  expect_identical(out$values$degree.max, 1L)
  expect_identical(out$values$regtype, "lp")
  expect_identical(out$metadata$user.supplied, "degree.max")
  expect_true(all(c("regtype", "search.engine", "degree.select", "bernstein.basis",
                    "degree.min", "degree.verify", "bwtype") %in%
                    out$metadata$auto.filled))
})

test_that("nomad shortcut marks automatic bandwidth calls for progress routing", {
  out <- local({
    x <- rnorm(10)
    y <- x + rnorm(10)
    np_bw_call_uses_nomad_degree_search(
      quote(npreg(y ~ x, nomad = TRUE, degree.max = 1L, nmulti = 1L)),
      caller_env = environment()
    )
  })

  expect_true(out)
})

test_that("npregbw nomad shortcut matches the explicit regression preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw_short <- np::npregbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )

  bw_long <- np::npregbw(
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

  expect_s3_class(bw_short, "rbandwidth")
  expect_equal(unname(bw_short$bw), unname(bw_long$bw), tolerance = 1e-8)
  expect_identical(as.integer(bw_short$degree), as.integer(bw_long$degree))
  expect_identical(bw_short$degree.search$mode, bw_long$degree.search$mode)
  expect_true(is.list(bw_short$nomad.shortcut))
  expect_identical(bw_short$nomad.shortcut$user.supplied, "degree.max")
  expect_true(all(c("regtype", "search.engine", "degree.select", "bernstein.basis",
                    "degree.min", "degree.verify", "bwtype") %in%
                    bw_short$nomad.shortcut$auto.filled))
})

test_that("npreg and npregbw accept nomad.nmulti through the NOMAD shortcut route", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- np::npregbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )
  fit <- np::npreg(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L,
    nomad.nmulti = 2L
  )

  expect_s3_class(bw, "rbandwidth")
  expect_true(is.list(fit$bws$nomad.shortcut))
  expect_identical(as.integer(fit$bws$degree), as.integer(bw$degree))
})

test_that("npreg nomad shortcut matches the explicit regression preset", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  dat <- data.frame(x = sort(runif(14)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  fit_short <- np::npreg(
    y ~ x,
    data = dat,
    nomad = TRUE,
    degree.max = 1L,
    nmulti = 1L
  )

  fit_long <- np::npreg(
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

  expect_equal(fitted(fit_short), fitted(fit_long), tolerance = 1e-8)
  expect_identical(as.integer(fit_short$bws$degree), as.integer(fit_long$bws$degree))
  expect_true(is.list(fit_short$bws$nomad.shortcut))
})

test_that("npregbw nomad shortcut fails fast on incompatible explicit settings", {
  x <- data.frame(x = seq(0, 1, length.out = 8))
  y <- sin(2 * pi * x$x)

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad = TRUE, regtype = "ll"),
    "nomad=TRUE requires regtype='lp'"
  )

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad = TRUE, bwtype = "adaptive_nn"),
    "nomad=TRUE currently requires bwtype='fixed'"
  )

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad = TRUE, degree = 1L),
    "does not support an explicit degree"
  )

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad = TRUE, bernstein.basis = FALSE),
    "requires bernstein.basis=TRUE"
  )

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad = TRUE, search.engine = "cell"),
    "requires search.engine='nomad' or 'nomad\\+powell'"
  )

  expect_error(
    np::npregbw(xdat = x, ydat = y, nomad.nmulti = 1L),
    "nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad\\+powell'"
  )

  expect_error(
    np::npregbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      degree.select = "coordinate",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      nomad.nmulti = 1L
    ),
    "nomad.nmulti is only supported when regtype='lp', automatic degree search is active, and search.engine is 'nomad' or 'nomad\\+powell'"
  )
})
