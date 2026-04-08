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

test_that("npplreg collective child evaluator matches the serial oracle at a fixed point", {
  close_mpi_slaves(force = TRUE)

  trace.file <- tempfile("npplreg-collective-trace-", fileext = ".log")
  old.trace.opt <- getOption("npRmpi.transport.trace.file")
  on.exit({
    options(npRmpi.transport.trace.file = old.trace.opt)
    unlink(trace.file)
  }, add = TRUE)

  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  ns <- asNamespace("npRmpi")
  options(npRmpi.transport.trace.file = trace.file)
  get(".npRmpi_bcast_cmd_expr", envir = ns, inherits = FALSE)(
    substitute(options(npRmpi.transport.trace.file = TF), list(TF = trace.file)),
    comm = 1L,
    caller.execute = TRUE
  )

  set.seed(20260408)
  n <- 24L
  xdat <- data.frame(x1 = runif(n), x2 = rnorm(n))
  zdat <- data.frame(z = sort(runif(n)))
  y <- 1 + xdat$x1 - 0.5 * xdat$x2 + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.05)

  reg.args <- list(
    regtype = "lp",
    basis = "glp",
    degree = 1L,
    bernstein.basis = TRUE,
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed"
  )
  outer.args <- reg.args[setdiff(names(reg.args), "bandwidth.compute")]
  bws <- matrix(0, nrow = 1L + ncol(xdat), ncol = ncol(zdat))
  degree.search <- list(
    engine = "nomad+powell",
    start.degree = 1L,
    candidates = list(c(0L, 1L)),
    lower = 0L,
    upper = 1L,
    basis = "glp",
    nobs = as.integer(n),
    start.user = FALSE,
    bernstein.basis = TRUE
  )

  state <- get(".npplregbw_nomad_prepare_state", envir = ns, inherits = FALSE)(
    xdat = xdat,
    ydat = y,
    zdat = zdat,
    bws = bws,
    reg.args = reg.args,
    outer.args = outer.args,
    degree.search = degree.search
  )

  bw.matrix <- matrix(c(0.7, 0.8, 0.9), nrow = length(state$child.templates), ncol = ncol(zdat))

  serial <- get(".npplregbw_eval_child_payload_serial", envir = ns, inherits = FALSE)(
    zdat = zdat,
    reg.args = reg.args,
    degree.search = degree.search,
    child.responses = state$child.responses,
    child.templates = state$child.templates,
    child.setup = state$child.setup,
    bw.matrix = bw.matrix,
    degree = 1L,
    penalty.multiplier = 10
  )

  collective.expr <- substitute(
    get(".npplregbw_eval_child_payload_collective", envir = asNamespace("npRmpi"), inherits = FALSE)(
      ZDAT,
      REGARGS,
      DEGREESEARCH,
      CHILDRESPONSES,
      CHILDTEMPLATES,
      CHILDSETUP,
      BWMATRIX,
      DEGREE,
      PENALTY,
      COMM,
      EVALID
    ),
    list(
      ZDAT = zdat,
      REGARGS = reg.args,
      DEGREESEARCH = degree.search,
      CHILDRESPONSES = state$child.responses,
      CHILDTEMPLATES = state$child.templates,
      CHILDSETUP = state$child.setup,
      BWMATRIX = bw.matrix,
      DEGREE = 1L,
      PENALTY = 10,
      COMM = 1L,
      EVALID = 1L
    )
  )

  collective <- get(".npRmpi_bcast_cmd_expr", envir = ns, inherits = FALSE)(
    collective.expr,
    comm = 1L,
    caller.execute = TRUE
  )

  expect_equal(collective$objective, serial$objective, tolerance = 1e-10)
  expect_equal(collective$num.feval, serial$num.feval, tolerance = 1e-10)
  expect_equal(collective$num.feval.fast, serial$num.feval.fast, tolerance = 1e-10)

  trace.lines <- readLines(trace.file, warn = FALSE)
  start.lines <- grep("role=npplreg.nomad\tevent=collective.eval.start", trace.lines, value = TRUE)
  done.lines <- grep("role=npplreg.nomad\tevent=collective.eval.done", trace.lines, value = TRUE)

  expect_true(length(start.lines) >= 2L)
  expect_true(length(done.lines) >= 2L)
  expect_true(any(grepl("rank=0", start.lines, fixed = TRUE)))
  expect_true(any(grepl("rank=1", start.lines, fixed = TRUE)))
  expect_true(any(grepl("child_indices=1($|\t)", start.lines)))
  expect_true(any(grepl("child_indices=2,3($|\t)", start.lines)))
})

test_that("npplregbw exhaustive degree search matches manual profile minimum", {
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
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
})

test_that("npplregbw automatic degree search defaults to NOMAD plus Powell", {
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

test_that("npplregbw NOMAD route uses the active pool for child objective evaluation", {
  skip_if_not_installed("crs")

  close_mpi_slaves(force = TRUE)

  trace.file <- tempfile("npplreg-nomad-route-trace-", fileext = ".log")
  old.trace.opt <- getOption("npRmpi.transport.trace.file")
  on.exit({
    options(npRmpi.transport.trace.file = old.trace.opt)
    unlink(trace.file)
  }, add = TRUE)

  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  ns <- asNamespace("npRmpi")
  options(npRmpi.transport.trace.file = trace.file)
  get(".npRmpi_bcast_cmd_expr", envir = ns, inherits = FALSE)(
    substitute(options(npRmpi.transport.trace.file = TF), list(TF = trace.file)),
    comm = 1L,
    caller.execute = TRUE
  )

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

  trace.lines <- readLines(trace.file, warn = FALSE)
  start.lines <- grep("role=npplreg.nomad\tevent=collective.eval.start", trace.lines, value = TRUE)
  done.lines <- grep("role=npplreg.nomad\tevent=collective.eval.done", trace.lines, value = TRUE)

  expect_true(length(start.lines) >= 2L)
  expect_true(length(done.lines) >= 2L)
  expect_true(any(grepl("rank=0", start.lines, fixed = TRUE)))
  expect_true(any(grepl("rank=1", start.lines, fixed = TRUE)))
  expect_true(any(grepl("child_indices=1($|\t)", start.lines)))
  expect_true(any(grepl("child_indices=2,3($|\t)", start.lines)))
})

test_that("npplreg explicit plbandwidth route preserves NOMAD child payload names", {
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

test_that("npplregbw NOMAD degree search fails fast when crs is unavailable", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 20
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = runif(n))
  y <- 1 + xdat$x + zdat$z + rnorm(n, sd = 0.08)

  expect_error(
    with_nprmpi_npplreg_degree_bindings(
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
      npplregbw(
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
