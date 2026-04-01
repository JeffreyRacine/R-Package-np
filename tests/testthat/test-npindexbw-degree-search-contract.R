with_nprmpi_npindex_degree_bindings <- function(bindings, code) {
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

capture_nprmpi_npindex_degree_messages_only <- function(expr) {
  messages <- character()
  withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  messages
}

nprmpi_npindex_degree_progress_time_values <- function(values) {
  force(values)
  i <- 0L
  function() {
    i <<- min(i + 1L, length(values))
    values[[i]]
  }
}

test_that("npindexbw exhaustive degree search matches manual Ichimura profile minimum", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 30
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- xdat$x1 + 0.75 * xdat$x2
  y <- sin(index) + 0.25 * index^2 + rnorm(n, sd = 0.05)
  start.bws <- c(1, 0.75, 0.35)

  bw0 <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    nmulti = 1L
  )
  bw1 <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    nmulti = 1L
  )
  auto <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_s3_class(auto, "sibandwidth")
  expect_true(isTRUE(auto$bernstein.basis))
  expect_identical(auto$degree.search$mode, "exhaustive")
  expect_true(isTRUE(auto$degree.search$completed))
  expect_true(isTRUE(auto$degree.search$certified))
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
  expect_lte(auto$degree.search$best.fval, auto$degree.search$baseline.fval + 1e-10)
  expect_true(all(c("degree", "fval", "status", "cached") %in% names(auto$degree.search$trace)))

  manual <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree = 1L,
    bwtype = "fixed",
    nmulti = 1L
  )
  expect_null(manual$degree.search)
})

test_that("npindexbw coordinate search can be exhaustively certified on a small grid", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 28
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- xdat$x1 + 0.5 * xdat$x2
  y <- cos(index) + 0.2 * index^2 + rnorm(n, sd = 0.05)
  start.bws <- c(1, 0.5, 0.3)

  exhaustive <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    nmulti = 1L
  )
  coordinate <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    degree.verify = TRUE,
    degree.restarts = 1L,
    degree.max.cycles = 4L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_identical(coordinate$degree.search$mode, "coordinate")
  expect_true(isTRUE(coordinate$degree.search$completed))
  expect_true(isTRUE(coordinate$degree.search$certified))
  expect_equal(as.integer(coordinate$degree), as.integer(exhaustive$degree))
  expect_equal(coordinate$fval, exhaustive$fval, tolerance = 1e-10)
})

test_that("npindexbw automatic degree search honors Klein-Spady objective direction", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 32
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- xdat$x1 - 0.6 * xdat$x2
  p <- plogis(1.25 * index - 0.4 * index^2)
  y <- rbinom(n, size = 1L, prob = p)
  start.bws <- c(1, -0.6, 0.35)

  bw0 <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "kleinspady",
    regtype = "lp",
    degree = 0L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    nmulti = 1L
  )
  bw1 <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "kleinspady",
    regtype = "lp",
    degree = 1L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    nmulti = 1L
  )
  auto <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "kleinspady",
    regtype = "lp",
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_identical(auto$method, "kleinspady")
  expect_lte(auto$fval, min(bw0$fval, bw1$fval) + 1e-10)
})

test_that("npindexbw automatic degree search enforces pilot guardrails", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  y <- xdat$x1 + xdat$x2 + rnorm(n, sd = 0.05)

  expect_error(
    npindexbw(
      xdat = xdat,
      ydat = y,
      bws = c(1, 1, 0.3),
      method = "ichimura",
      regtype = "lc",
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      nmulti = 1L
    ),
    "automatic degree search currently requires regtype='lp'"
  )

  expect_error(
    npindexbw(
      xdat = xdat,
      ydat = y,
      bws = c(1, 1, 0.3),
      method = "ichimura",
      regtype = "lp",
      bandwidth.compute = FALSE,
      degree.select = "exhaustive",
      search.engine = "cell",
      degree.min = 0L,
      degree.max = 1L,
      bwtype = "fixed",
      nmulti = 1L
    ),
    "bandwidth.compute=TRUE"
  )

  bw <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = c(1, 1, 0.3),
    method = "ichimura",
    regtype = "lp",
    bernstein.basis = FALSE,
    degree.select = "exhaustive",
    search.engine = "cell",
    degree.min = 0L,
    degree.max = 4L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_s3_class(bw, "sibandwidth")
  expect_false(isTRUE(bw$bernstein.basis))
  expect_lte(max(as.integer(bw$degree)), 4L)
})

test_that("npindexbw NOMAD degree search backend improves over the baseline", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- xdat$x1 + 0.6 * xdat$x2
  y <- sin(index) + 0.2 * index^2 + rnorm(n, sd = 0.05)
  start.bws <- c(1, 0.6, 0.3)

  bw <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree.select = "coordinate",
    search.engine = "nomad",
    degree.min = 0L,
    degree.max = 2L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_s3_class(bw, "sibandwidth")
  expect_identical(bw$degree.search$mode, "nomad")
  expect_true(isTRUE(bw$degree.search$completed))
  expect_gte(bw$degree.search$n.unique, 1L)
  expect_lte(bw$degree.search$best.fval, bw$degree.search$baseline.fval + 1e-10)
})

test_that("npindexbw automatic degree search defaults to NOMAD plus Powell", {
  skip_if_not_installed("crs")
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 26
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- xdat$x1 + 0.6 * xdat$x2
  y <- sin(index) + 0.2 * index^2 + rnorm(n, sd = 0.05)
  start.bws <- c(1, 0.6, 0.3)

  bw <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = start.bws,
    method = "ichimura",
    regtype = "lp",
    degree.select = "coordinate",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
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

test_that("npindexbw NOMAD degree search fails fast when crs is unavailable", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 20
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  y <- xdat$x1 + xdat$x2 + rnorm(n, sd = 0.05)

  expect_error(
    with_nprmpi_npindex_degree_bindings(
      list(.np_nomad_require_crs = function() stop("crs missing", call. = FALSE)),
      npindexbw(
        xdat = xdat,
        ydat = y,
        bws = c(1, 1, 0.3),
        method = "ichimura",
        regtype = "lp",
        degree.select = "coordinate",
        search.engine = "nomad",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        nmulti = 1L
      )
    ),
    "crs missing"
  )
})

test_that("npindex forwards automatic LP degree search through npindexbw", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = FALSE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 24
  dat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  index <- dat$x1 + 0.5 * dat$x2
  dat$y <- sin(index) + rnorm(n, sd = 0.05)

  fit <- npindex(
    y ~ x1 + x2,
    data = dat,
    method = "ichimura",
    regtype = "lp",
    degree.select = "exhaustive",
    degree.min = 0L,
    degree.max = 1L,
    bwtype = "fixed",
    nmulti = 1L
  )

  expect_s3_class(fit, "singleindex")
  expect_s3_class(fit$bws, "sibandwidth")
  expect_false(is.null(fit$bws$degree.search))
  expect_identical(fit$bws$degree.search$mode, "exhaustive")
})

test_that("npindexbw automatic degree search emits staged progress output", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  old_opts <- options(np.messages = TRUE, np.tree = FALSE, npRmpi.autodispatch = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260319)
  n <- 20
  xdat <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1)
  )
  y <- xdat$x1 + 0.5 * xdat$x2 + rnorm(n, sd = 0.05)

  msgs <- with_nprmpi_npindex_degree_bindings(
    list(
      .np_progress_is_interactive = function() TRUE,
      .np_progress_renderer_for_surface = function(surface, capability) "legacy",
      .np_progress_now = nprmpi_npindex_degree_progress_time_values(seq(0, 20, by = 0.5))
    ),
    capture_nprmpi_npindex_degree_messages_only(
      get("npindexbw", envir = asNamespace("npRmpi"), inherits = FALSE)(
        xdat = xdat,
        ydat = y,
        bws = c(1, 0.5, 0.3),
        method = "ichimura",
        regtype = "lp",
        degree.select = "exhaustive",
        degree.min = 0L,
        degree.max = 1L,
        bwtype = "fixed",
        nmulti = 1L
      )
    )
  )

  expect_true(any(grepl("Automatic polynomial degree search baseline \\(0\\)", msgs)))
  expect_true(any(grepl("Selecting degree and bandwidth", msgs, fixed = TRUE)))
  expect_true(any(grepl("exhaustive", msgs)))
  expect_true(any(grepl("best (", msgs, fixed = TRUE)))
})
