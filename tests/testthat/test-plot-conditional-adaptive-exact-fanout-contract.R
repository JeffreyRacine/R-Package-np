test_that("conditional adaptive-NN exact bootstrap partitions replications across MPI fanout", {
  core <- getFromNamespace(".np_inid_boot_from_ksum_conditional_exact", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$fanout <- FALSE

  xdat <- data.frame(x = c(1, 2, 3, 4))
  ydat <- data.frame(y = c(2, 4, 6, 8))
  exdat <- data.frame(x = c(1.5, 2.5))
  eydat <- data.frame(y = c(3, 5))
  bws <- list(type = "adaptive_nn")
  counts <- matrix(
    c(
      1, 0, 2,
      0, 1, 0,
      1, 1, 1,
      2, 2, 1
    ),
    nrow = 4L,
    ncol = 3L
  )

  local_mocked_bindings(
    .np_con_inid_ksum_eligible = function(bws) TRUE,
    .npRmpi_has_active_slave_pool = function(comm = 1L) TRUE,
    .npRmpi_bootstrap_tune_chunk_size = function(B, chunk.size, comm = 1L, include.master = TRUE) 2L,
    .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) {
      list(
        list(start = 1L, bsz = 2L, seed = 11L),
        list(start = 3L, bsz = 1L, seed = 12L)
      )
    },
    .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
    .npRmpi_bootstrap_run_fanout = function(tasks, worker, ncol.out, what = "bootstrap", ...) {
      called$fanout <- TRUE
      expect_identical(what, "inid-ksum-conditional-adaptive-exact-counts")
      expect_identical(ncol.out, 2L)
      do.call(rbind, lapply(tasks, worker))
    },
    .np_con_make_kbandwidth_x = function(...) list(ok = TRUE),
    .np_con_make_kbandwidth_xy = function(...) list(ok = TRUE),
    .np_conditional_exact_fit_or_stop = function(fit.expr, ...) fit.expr(),
    .np_ksum_conditional_eval_exact = function(xdat, ydat, ...) {
      c(sum(xdat$x), sum(ydat$y))
    },
    .package = "npRmpi"
  )
  withr::local_options(npRmpi.mpi.initialized = TRUE)

  out <- core(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    bws = bws,
    B = ncol(counts),
    cdf = FALSE,
    counts = counts,
    progress.label = "unit"
  )

  expect_true(called$fanout)
  expect_identical(out$t0, c(10, 20))
  expect_identical(dim(out$t), c(3L, 2L))
})

