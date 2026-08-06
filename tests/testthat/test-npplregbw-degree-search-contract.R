skip_slow_npplregbw_degree_search <- function() {
  skip_on_cran()
}

with_nprmpi_npplreg_degree_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("npRmpi")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked)
      unlockBinding(name, ns)
    assign(name, bindings[[name]], envir = ns)
    if (was_locked)
      lockBinding(name, ns)
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked)
        unlockBinding(name, ns)
      assign(name, old[[name]], envir = ns)
      if (was_locked)
        lockBinding(name, ns)
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

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
  ns <- asNamespace("npRmpi")
  child.search <- get(".npplregbw_child_specific_nomad_search", ns,
                      inherits = FALSE)
  child.args <- get(".npplregbw_child_nomad_call_args", ns, inherits = FALSE)
  child.search.text <- paste(deparse(body(child.search)), collapse = "\n")
  child.args.text <- paste(deparse(body(child.args)), collapse = "\n")

  expect_match(child.search.text, ".npplregbw_child_nomad_call_args",
               fixed = TRUE)
  expect_match(child.search.text, "do.call(npregbw, child.args)", fixed = TRUE)
  expect_match(child.args.text, "nomad = TRUE", fixed = TRUE)
  expect_match(child.args.text, "search.engine = degree.search$engine",
               fixed = TRUE)
})

test_that("npplregbw exhaustive degree search matches manual profile minimum", {
  skip_slow_npplregbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
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
  skip_slow_npplregbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
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
  skip_slow_npplregbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
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
  skip_slow_npplregbw_degree_search()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
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

test_that("npplregbw NOMAD degree search backend improves over the baseline", {
  skip_slow_npplregbw_degree_search()
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  bw <- npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    bwmethod = "cv.ls",
    nmulti = 1L
  )

  expect_s3_class(bw, "plbandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
})

test_that("npplregbw automatic degree search defaults to NOMAD plus Powell", {
  skip_slow_npplregbw_degree_search()
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + 0.75 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.08)

  bw <- npplregbw(
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

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))
  expect_equal(as.double(bw$total.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
  expect_equal(as.double(bw$degree.search$optim.time),
               as.double(bw$nomad.time + bw$powell.time),
               tolerance = 1e-8)
})

test_that("npplregbw NOMAD route completes under an active pool", {
  skip_slow_npplregbw_degree_search()
  skip_if_not_installed("crs")

  close_mpi_slaves(force = TRUE)

  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(42)
  n <- 24L
  xdat <- data.frame(
    x1 = rnorm(n),
    x2 = factor(rbinom(n, 1, 0.5))
  )
  zdat <- data.frame(
    z1 = factor(rbinom(n, 1, 0.5)),
    z2 = rnorm(n)
  )
  y <- 1 + xdat$x1 + as.numeric(xdat$x2) + as.numeric(zdat$z1) + sin(zdat$z2) + rnorm(n, sd = 0.2)

  bw <- npplregbw(
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

  expect_identical(bw$degree.search$mode, "nomad+powell")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_true(is.finite(bw$nomad.time))
  expect_true(is.finite(bw$powell.time))

  expect_true(is.finite(bw$fval))
  expect_true(is.finite(bw$degree.search$best.fval))
})

test_that("npplreg explicit plbandwidth route preserves NOMAD child payload names", {
  skip_slow_npplregbw_degree_search()
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
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
