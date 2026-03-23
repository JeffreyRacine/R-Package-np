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
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
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

test_that("npregbw forwards nomad.nmulti to inner crs::snomadr multistarts", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260323)
  dat <- data.frame(x = sort(runif(16)))
  dat$y <- sin(2 * pi * dat$x) + rnorm(nrow(dat), sd = 0.05)

  assign("np_nomad_inner_nmulti_seen", integer(0), envir = .GlobalEnv)
  on.exit(rm(list = "np_nomad_inner_nmulti_seen", envir = .GlobalEnv), add = TRUE)
  trace(
    what = crs::snomadr,
    tracer = quote({
      .GlobalEnv[["np_nomad_inner_nmulti_seen"]] <- c(
        .GlobalEnv[["np_nomad_inner_nmulti_seen"]],
        as.integer(nmulti)
      )
    }),
    print = FALSE
  )
  on.exit(untrace(crs::snomadr), add = TRUE)

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

  expect_identical(get("np_nomad_inner_nmulti_seen", envir = .GlobalEnv), 2L)
})
