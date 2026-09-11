test_that("conditional gradient bootstrap uses fanout with an active pool", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_conditional_gradient", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  called$centers <- 0L

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 2L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .np_plot_conditional_eval = function(...) {
      called$centers <- called$centers + 1L
      list(congrad = matrix(c(0.25, 0.75), ncol = 1L))
    },
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size, with.seeds = TRUE) {
      expect_false(with.seeds)
      list(
        list(start = 1L, bsz = 2L, seed = 11L),
        list(start = 3L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            master_local_chunk = FALSE,
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "conditional-gradient")
      expect_identical(ncol.out, 2L)
      expect_true(master_local_chunk)
      expect_false("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_conditional_gradient_local = function(xdat,
                                                            ydat,
                                                            exdat,
                                                            eydat,
                                                            bws,
                                                            B,
                                                            cdf,
                                                            gradient.index,
                                                            gradient.order = 1L,
                                                            counts = NULL,
                                                            counts.drawer = NULL,
                                                            progress.label = NULL,
                                                            center = NULL) {
      expect_identical(center, c(0.25, 0.75))
      list(
        t = matrix(rep(seq_len(B), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(0.25, 0.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = data.frame(y = 1:3),
    exdat = data.frame(x = 1:2),
    eydat = data.frame(y = 1:2),
    bws = list(xdati = list(icon = TRUE)),
    B = 3L,
    cdf = FALSE,
    gradient.index = 1L
  )

  expect_true(called$fanout)
  expect_identical(called$centers, 1L)
  expect_identical(out$t0, c(0.25, 0.75))
  expect_identical(dim(out$t), c(3L, 2L))
})

test_that("quantile gradient bootstrap uses fanout with fixed counts", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_quantile_gradient", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  called$centers <- 0L
  counts <- matrix(c(1, 0, 2, 0, 2, 1), nrow = 3L, ncol = 2L)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 1L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .np_plot_quantile_eval = function(...) {
      called$centers <- called$centers + 1L
      list(quantgrad = matrix(c(1.25, 1.75), ncol = 1L))
    },
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size, with.seeds = TRUE) {
      expect_false(with.seeds)
      list(
        list(start = 1L, bsz = 1L, seed = 11L),
        list(start = 2L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "quantile-gradient")
      expect_identical(ncol.out, 2L)
      expect_true("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_quantile_gradient_local = function(xdat,
                                                         ydat,
                                                         exdat,
                                                         bws,
                                                         B,
                                                         tau,
                                                         gradient.index,
                                                         counts = NULL,
                                                         counts.drawer = NULL,
                                                         progress.label = NULL,
                                                         center = NULL) {
      expect_identical(center, c(1.25, 1.75))
      list(
        t = matrix(rep(colSums(counts), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(1.25, 1.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = 1:3,
    exdat = data.frame(x = 1:2),
    bws = list(),
    B = 2L,
    tau = 0.5,
    gradient.index = 1L,
    counts = counts
  )

  expect_true(called$fanout)
  expect_identical(called$centers, 1L)
  expect_identical(out$t0, c(1.25, 1.75))
  expect_identical(dim(out$t), c(2L, 2L))
})

test_that("quantile level bootstrap uses fanout with fixed counts", {
  boot <- getFromNamespace(".npRmpi_inid_boot_from_quantile_level", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE
  called$centers <- 0L
  counts <- matrix(c(1, 0, 2, 0, 2, 1), nrow = 3L, ncol = 2L)

  local_mocked_bindings(
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 1L,
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .np_plot_quantile_eval = function(...) {
      called$centers <- called$centers + 1L
      list(quantile = c(2.25, 2.75))
    },
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size, with.seeds = TRUE) {
      expect_false(with.seeds)
      list(
        list(start = 1L, bsz = 1L, seed = 11L),
        list(start = 2L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_run_fanout = function(tasks,
                                            worker,
                                            ncol.out,
                                            what = "bootstrap",
                                            required.bindings = list(),
                                            ...) {
      called$fanout <- TRUE
      expect_identical(what, "quantile-level")
      expect_identical(ncol.out, 2L)
      expect_true("counts.mat" %in% names(required.bindings))
      do.call(rbind, lapply(tasks, worker))
    },
    .np_inid_boot_from_quantile_level_local = function(xdat,
                                                       ydat,
                                                       exdat,
                                                       bws,
                                                       B,
                                                       tau,
                                                       counts = NULL,
                                                       counts.drawer = NULL,
                                                       progress.label = NULL,
                                                       center = NULL) {
      expect_identical(center, c(2.25, 2.75))
      list(
        t = matrix(rep(colSums(counts), each = 2L), nrow = B, ncol = 2L, byrow = TRUE),
        t0 = c(2.25, 2.75)
      )
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- boot(
    xdat = data.frame(x = 1:3),
    ydat = 1:3,
    exdat = data.frame(x = 1:2),
    bws = list(),
    B = 2L,
    tau = 0.5,
    counts = counts
  )

  expect_true(called$fanout)
  expect_identical(called$centers, 1L)
  expect_identical(out$t0, c(2.25, 2.75))
  expect_identical(dim(out$t), c(2L, 2L))
})

test_that("private bootstrap centers preserve values and reject malformed shape", {
  center <- getFromNamespace(".np_plot_bootstrap_supplied_center", "npRmpi")
  value <- structure(c(NA_real_, 2), names = c("a", "b"), private.note = "retained")
  expect_identical(center(value, 2L, "test"), value)
  for (bad in list(matrix(1:2, nrow = 1L), array(1:2, c(2L, 1L, 1L)),
                   1, list(1, 2), c("1", "2"), c(1i, 2i)))
    expect_error(center(bad, 2L, "test"), "invalid private bootstrap center")
})

test_that("local bootstrap centers are reused without changing replicate work", {
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  fit <- function(y, nout) {
    calls$n <- calls$n + 1L
    rep.int(sum(y), nout)
  }
  local_mocked_bindings(
    .np_plot_conditional_eval = function(ydat, exdat, ...) {
      list(congrad = matrix(fit(ydat[[1L]], nrow(exdat)), ncol = 1L))
    },
    .np_plot_quantile_eval = function(tydat, exdat, tau, gradients, ...) {
      value <- fit(tydat, nrow(exdat) * length(tau))
      if (gradients) list(quantgrad = array(value, c(nrow(exdat), 1L, length(tau))))
      else list(quantile = matrix(value, nrow(exdat), length(tau)))
    },
    .package = "npRmpi"
  )
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  x <- data.frame(x = 1:3)
  counts <- matrix(c(1, 1, 1, 2, 0, 1), nrow = 3L)
  cases <- list(
    list(fun = ".np_inid_boot_from_conditional_gradient_local",
      args = list(xdat = x, ydat = x, exdat = x[1:2, , drop = FALSE],
        eydat = x[1:2, , drop = FALSE], bws = list(xdati = list(icon = TRUE)),
        B = 2L, cdf = FALSE, gradient.index = 1L, counts = counts)),
    list(fun = ".np_inid_boot_from_quantile_level_local",
      args = list(xdat = x, ydat = 1:3, exdat = x[1:2, , drop = FALSE],
        bws = list(), B = 2L, tau = c(.25, .75), counts = counts)),
    list(fun = ".np_inid_boot_from_quantile_gradient_local",
      args = list(xdat = x, ydat = 1:3, exdat = x[1:2, , drop = FALSE],
        bws = list(), B = 2L, tau = c(.25, .75), gradient.index = 1L, counts = counts)))
  for (case in cases) {
    f <- getFromNamespace(case$fun, "npRmpi")
    calls$n <- 0L
    original <- do.call(f, case$args)
    expect_identical(calls$n, 3L)
    calls$n <- 0L
    reused <- do.call(f, c(case$args, list(center = original$t0)))
    expect_identical(calls$n, 2L)
    expect_identical(reused, original)
    expect_error(do.call(f, c(case$args, list(center = matrix(original$t0)))),
                 "invalid private bootstrap center")
    bad.B <- case$args
    bad.B$B <- 0L
    expect_error(do.call(f, c(bad.B, list(center = original$t0))), "invalid")
  }
})
