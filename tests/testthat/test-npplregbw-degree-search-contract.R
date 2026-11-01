test_that("degree-search metadata preserves structured child degrees", {
  children <- list(yzbw = c(1, 2), x1 = c(0, 1))
  metadata <- .np_degree_search_metadata(
    list(
      method = "nomad",
      direction = "min",
      completed = TRUE,
      baseline = list(degree = children, objective = 2),
      best = list(degree = children, objective = 1)
    ),
    default_direction = "min"
  )

  expect_identical(metadata$baseline.degree, lapply(children, as.integer))
  expect_identical(metadata$best.degree, lapply(children, as.integer))
  expect_identical(names(metadata$best.degree), names(children))
})

test_that("npplreg child searches use the canonical native regression owner", {
  ns <- asNamespace("np")
  child.search <- get(".npplregbw_child_specific_nomad_search", ns,
                      inherits = FALSE)
  child.args <- get(".npplregbw_child_nomad_call_args", ns, inherits = FALSE)
  child.search.text <- paste(deparse(body(child.search)), collapse = "\n")
  child.args.text <- paste(deparse(body(child.args)), collapse = "\n")

  expect_match(child.search.text, ".npplregbw_child_nomad_call_args",
               fixed = TRUE)
  expect_match(child.search.text, "do.call(.npregbw_auxiliary, child.args)", fixed = TRUE)
  expect_match(child.args.text, "nomad = TRUE", fixed = TRUE)
  expect_match(child.args.text, "search.engine = degree.search$engine",
               fixed = TRUE)
})

test_that("categorical child search executes the canonical constructor and native NOMAD", {
  skip_if_not_installed("crs")

  ns <- asNamespace("np")
  state <- new.env(parent = emptyenv()); state$private.factor <- 0L; state$native <- 0L
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE,
                           np.test.child.owner = state))
  trace(".npregbw_default_impl", where = ns, print = FALSE, tracer = quote({
    if (!.np.context$public && is.factor(ydat)) {
      witness <- getOption("np.test.child.owner")
      witness$private.factor <- witness$private.factor + 1L
    }
  }))
  withr::defer(untrace(".npregbw_default_impl", where = ns))
  trace("npNomadNativeSearchRegression", where = ns, print = FALSE, tracer = quote({
    .Call <- function(.NAME, ...) {
      if (identical(.NAME, "C_np_regression_nomad_native_search")) {
        witness <- getOption("np.test.child.owner")
        witness$native <- witness$native + 1L
      }
      base::.Call(.NAME, ...)
    }
  }))
  withr::defer(untrace("npNomadNativeSearchRegression", where = ns))
  n <- 24L
  z <- data.frame(z = seq(.05, .95, length.out = n))
  x <- data.frame(x = factor(rep(c("30", "10", "50"), 8),
                            levels = c("50", "10", "30", "70")))
  b <- npplregbw(xdat = x, ydat = sin(z$z), zdat = z,
    regtype = "lp", degree.select = "coordinate", degree.min = 0L,
    degree.max = 1L, search.engine = "nomad", nmulti = 1L, itmax = 2L,
    nomad.opts = list(MAX_BB_EVAL = 8L))
  expect_gt(state$private.factor, 0L)
  expect_gte(state$native, 2L)
  expect_s3_class(b, "plbandwidth")
  expect_error(npreg(b$bw[[2L]]), "must be a numeric")
  expect_error(npregbw(xdat = z, ydat = x$x, bws = .3,
                       bandwidth.compute = FALSE), "must be a numeric")
})

test_that("npplregbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 28
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  bw0 <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  bw1 <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  auto <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(auto, "plbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npplregbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z1 = runif(n),
    z2 = runif(n)
  )
  y <- 1 + 0.5 * xdat$x + sin(2 * pi * zdat$z1) + zdat$z2^2 + rnorm(n, sd = 0.08)

  exhaustive <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )
  coordinate <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
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

test_that("npplregbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  y <- 1 + xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  expect_error(
    npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lc",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    npplregbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      bandwidth.compute = FALSE,
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bws = matrix(0.2, nrow = 2L, ncol = 1L)
    ),
    "bandwidth.compute=TRUE"
  )

  bw <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
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

  expect_s3_class(bw, "plbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)
})

test_that("npplreg forwards automatic LP degree search through npplregbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  dat <- data.frame(
    x = rnorm(n),
    z = runif(n)
  )
  dat$y <- 1 + 0.75 * dat$x + sin(2 * pi * dat$z) + rnorm(n, sd = 0.08)

  fit <- npplreg(
    y ~ x | z,
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

  expect_s3_class(fit, "plregression")
  expect_s3_class(fit$bws, "plbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npplregbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  auto <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(auto, "plbandwidth")
  expect_false(is.null(auto$degree.search))
  expect_identical(auto$degree.search$mode, "nomad+powell")
  expect_equal(length(auto$degree), 1L)
  expect_true(is.finite(auto$fval))
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-8)
  expect_lte(auto$fval, auto$degree.search$baseline.fval + 1e-8)
  expect_equal(auto$fval, auto$degree.search$best.fval, tolerance = 1e-8)
})

test_that("npplreg explicit plbandwidth route preserves NOMAD child payload names", {
  skip_if_not_installed("crs")

  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260322)
  n <- 40L
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd = 0.1)

  bw <- npplregbw(
    xdat = data.frame(x1 = x1),
    zdat = data.frame(x2 = x2),
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad+powell",
    degree.max = 2L,
    nmulti = 1L
  )

  expect_identical(names(bw$bw), c("yzbw", "x1"))

  fit <- npplreg(
    bws = bw,
    txdat = data.frame(x1 = x1),
    tzdat = data.frame(x2 = x2),
    tydat = y
  )

  expect_s3_class(fit, "plregression")
  expect_equal(nrow(fit$evalx), n)
})
