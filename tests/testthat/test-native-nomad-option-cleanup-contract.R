native_cache_off_opts <- list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 1)
native_cache_unknown_opts <- list(EVAL_USE_CACHE = NA, MAX_BB_EVAL = 1)

test_that("native NOMAD option builders reject cache-off solves before solver entry", {
  builders <- list(
    np:::.np_nomad_native_option_vectors,
    np:::.npudensbw_nomad_native_option_vectors,
    np:::.npcdensbw_nomad_shadow_native_option_vectors,
    np:::.npcdistbw_nomad_native_option_vectors,
    np:::.npregbw_nomad_native_option_vectors
  )

  for (builder in builders) {
    expect_error(
      builder(native_cache_off_opts),
      "requires EVAL_USE_CACHE = TRUE",
      fixed = TRUE
    )
    expect_error(
      builder(native_cache_unknown_opts),
      "requires EVAL_USE_CACHE = TRUE",
      fixed = TRUE
    )
    expect_no_error(builder(list(EVAL_USE_CACHE = TRUE, MAX_BB_EVAL = 1)))
  }
})

test_that("native R callback failure does not poison the next package solve", {
  bad <- np:::.np_nomad_native_r_callback_search(
    eval.f = function(x) stop("intentional package R callback failure"),
    x0 = 0,
    bbin = 0L,
    lb = -1,
    ub = 1,
    random.seed = 42L,
    option.names = c("MAX_BB_EVAL", "NB_THREADS_PARALLEL_EVAL"),
    option.values = c("1", "1")
  )
  expect_equal(bad$value$status, "error")
  expect_match(bad$value$message, "callback", ignore.case = TRUE)

  good <- np:::.np_nomad_native_r_callback_search(
    eval.f = function(x) sum((x - 0.25)^2),
    x0 = 0,
    bbin = 0L,
    lb = -1,
    ub = 1,
    random.seed = 42L,
    option.names = c("MAX_BB_EVAL", "NB_THREADS_PARALLEL_EVAL"),
    option.values = c("10", "1")
  )
  expect_equal(good$value$status, "ok")
  expect_equal(good$value$native_status, 0L)
  expect_equal(good$value$result_status, 0L)
  expect_true(is.finite(good$value$objective))
})

test_that("native npreg solver option failure does not poison the next serial solve", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.developer.native.nomad.diagnostics = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(111)
  n <- 40L
  dat <- data.frame(x = runif(n))
  dat$y <- sin(dat$x) + rnorm(n, sd = 0.2)

  expect_error(
    np::npregbw(
      y ~ x,
      data = dat,
      nomad = TRUE,
      nmulti = 1L,
      degree.max = 2L,
      nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)
    ),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )

  fit <- np::npregbw(
    y ~ x,
    data = dat,
    nomad = TRUE,
    nmulti = 1L,
    degree.max = 2L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )

  expect_true(is.finite(fit$fval))
  expect_false(is.null(attr(fit, "native.nomad.diagnostics")))
})

test_that("native bandwidth constructors reject cache-off options before solver work", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.developer.native.nomad.diagnostics = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(913)
  n <- 30L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = rnorm(n))

  expect_error(
    np::npudensbw(dat = x, bwsolver = "mads", nmulti = 1L,
                  nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )
  expect_true(is.finite(np::npudensbw(
    dat = x, bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )$fval))

  expect_error(
    np::npudistbw(dat = x, bwsolver = "mads", nmulti = 1L,
                  nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )
  expect_true(is.finite(np::npudistbw(
    dat = x, bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )$fval))

  expect_error(
    np::npcdensbw(ydat = y, xdat = x, bwsolver = "mads", nmulti = 1L,
                  nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )
  expect_true(is.finite(np::npcdensbw(
    ydat = y, xdat = x, bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )$fval))

  expect_error(
    np::npcdistbw(ydat = y, xdat = x, bwsolver = "mads", nmulti = 1L,
                  nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )
  expect_true(is.finite(np::npcdistbw(
    ydat = y, xdat = x, bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )$fval))
})

test_that("conditional-density native shadow cleanup survives degree search before option rejection", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE
  )
  on.exit(options(old), add = TRUE)

  set.seed(42)
  n <- 60L
  x <- runif(n)
  y <- rbeta(n, 1, 1)
  # npcdens() rebuilds its bandwidth call in the caller frame.
  npcdensbw <- np::npcdensbw

  fit <- np::npcdens(
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    cxkerbound = "range",
    cykerbound = "range",
    regtype = "lp",
    degree.min = 0L,
    degree.max = 3L,
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.verify = FALSE,
    bwtype = "fixed",
    bwmethod = "cv.ml",
    nmulti = 2L
  )
  expect_true(is.finite(fit$bws$fval))

  x2 <- data.frame(x = runif(30L))
  y2 <- data.frame(y = rnorm(30L))
  expect_error(
    np::npcdensbw(ydat = y2, xdat = x2, bwsolver = "mads", nmulti = 1L,
                  nomad.opts = list(EVAL_USE_CACHE = FALSE, MAX_BB_EVAL = 5)),
    "requires EVAL_USE_CACHE = TRUE",
    fixed = TRUE
  )
  expect_true(is.finite(np::npcdensbw(
    ydat = y2, xdat = x2, bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 8)
  )$fval))
})
