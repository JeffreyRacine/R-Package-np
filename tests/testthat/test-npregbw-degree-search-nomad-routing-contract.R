test_that("npregbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  bw <- np::npregbw(
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
})

test_that("npregbw NOMAD degree search fails fast when crs is unavailable", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  dat <- data.frame(x = runif(16), y = rnorm(16))

  expect_error(
    with_np_degree_bindings(
      list(.npregbw_nomad_native_require_crs = function()
        stop("crs missing", call. = FALSE)),
      np::npregbw(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    ),
    "crs missing"
  )
})

test_that("npregbw NOMAD degree search fails fast on native solver errors after baseline", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260519)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- dat$x + rnorm(nrow(dat), sd = 0.05)

  expect_error(
    with_np_degree_bindings(
      list(npNomadNativeSearchRegression = function(...)
        stop("forced native NOMAD failure", call. = FALSE)),
      np::npregbw(
        y ~ x,
        data = dat,
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        bwmethod = "cv.ls",
        nmulti = 1L
      )
    ),
    "forced native NOMAD failure"
  )
})

test_that("npregbw forwards nomad.nmulti to inner native NOMAD multistarts", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  inner.nmulti.seen <- integer(0)
  native.search <- getFromNamespace("npNomadNativeSearchRegression", "np")

  with_np_degree_bindings(
    list(npNomadNativeSearchRegression = function(...) {
      args <- list(...)
      inner.nmulti.seen <<- c(inner.nmulti.seen,
                              as.integer(args$inner.start.count[1L]))
      do.call(native.search, args)
    }),
    np::npregbw(
      y ~ x,
      data = dat,
      regtype = "lp",
      degree.select = "coordinate",
      search.engine = "nomad",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L,
      nomad.nmulti = 2L
    )
  )

  expect_true(any(inner.nmulti.seen == 2L))
})
