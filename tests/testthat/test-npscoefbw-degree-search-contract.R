with_np_degree_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("np")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked) {
      unlockBinding(name, ns)
    }
    assign(name, bindings[[name]], envir = ns)
    if (was_locked) {
      lockBinding(name, ns)
    }
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked) {
        unlockBinding(name, ns)
      }
      assign(name, old[[name]], envir = ns)
      if (was_locked) {
        lockBinding(name, ns)
      }
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

test_that("npscoefbw exhaustive degree search matches manual profile minimum", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 28
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- (1 + zdat$z^2) * xdat$x + rnorm(n, sd = 0.08)

  bw0 <- npscoefbw(
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
  bw1 <- npscoefbw(
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
  auto <- npscoefbw(
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

  expect_s3_class(auto, "scbandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npscoefbw(
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

test_that("npscoefbw coordinate search can be exhaustively certified on a small grid", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(
    z1 = runif(n),
    z2 = runif(n)
  )
  y <- (1 + zdat$z1 + zdat$z2^2) * xdat$x + rnorm(n, sd = 0.08)

  exhaustive <- npscoefbw(
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
  coordinate <- npscoefbw(
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

test_that("npscoefbw automatic degree search enforces pilot guardrails", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  y <- (1 + zdat$z) * xdat$x + rnorm(n, sd = 0.08)

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
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

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      bandwidth.compute = FALSE,
      degree.select = "exhaustive",
      degree.min = 0L,
      degree.max = 1L,
      bws = 0.2
    ),
    "bandwidth.compute=TRUE"
  )

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 4L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "scbandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)

  expect_error(
    npscoefbw(
      xdat = xdat,
      zdat = zdat,
      ydat = y,
      regtype = "lp",
      degree.select = "coordinate",
      degree.min = 0L,
      degree.max = 1L,
      cv.iterate = TRUE,
      bwtype = "fixed",
      bwmethod = "cv.ls",
      nmulti = 1L
    ),
    "cv.iterate=FALSE"
  )
})

test_that("npscoef forwards automatic LP degree search through npscoefbw", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  dat <- data.frame(
    x = runif(n),
    z = runif(n)
  )
  dat$y <- (1 + dat$z^2) * dat$x + rnorm(n, sd = 0.08)

  fit <- npscoef(
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

  expect_s3_class(fit, "smoothcoefficient")
  expect_s3_class(fit$bws, "scbandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npscoefbw automatic degree search defaults to NOMAD plus Powell", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- (1 + zdat$z^2) * xdat$x + rnorm(n, sd = 0.08)

  bw <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "scbandwidth")
  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(is.finite(bw$nomad.time))
})

test_that("npscoefbw NOMAD degree search fails fast when crs is unavailable", {
  old_opts <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- (1 + zdat$z) * xdat$x + rnorm(n, sd = 0.08)

  expect_error(
    with_np_degree_bindings(
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
      npscoefbw(
        xdat = xdat,
        zdat = zdat,
        ydat = y,
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
