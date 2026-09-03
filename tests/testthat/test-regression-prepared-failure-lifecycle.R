test_that("terminal prepared NN degree failure does not kill workers", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- as.data.frame(lapply(1:2, function(j) c(rep(0, 32L),
    sin((1:16) * sqrt(j + 1)) + (1:16) / (16 * (j + 1)))))
  y <- sin((1:48) / 3) + (1:48) / 90
  for (type in c("generalized_nn", "adaptive_nn")) {
    expect_error(do.call(npregbw, list(xdat = x, ydat = y, bws = c(7, 7),
      bwtype = type, regtype = "lp", nomad = TRUE, search.engine = "nomad",
      degree.select = "coordinate", degree.min = 0L, degree.max = 2L,
      degree.start = c(1L, 1L), degree.verify = FALSE, nmulti = 1L,
      itmax = 20L, nomad.opts = list(MAX_BB_EVAL = 30L))),
      "native npreg NOMAD degree-search route did not return a raw-valid solution")
    ranks <- mpi.bcast.cmd(npRmpi:::mpi.allgather.Robj(mpi.comm.rank(1L)),
      caller.execute = TRUE)
    expect_identical(as.integer(ranks), seq.int(0L, mpi.comm.size(1L) - 1L))
    cleared <- mpi.bcast.cmd(npRmpi:::mpi.allgather.Robj(tryCatch(
      npRmpi:::npRmpiPreparedObjectiveEvalRegressionRaw(c(44, 44), c(1L, 1L)),
      error = function(e) conditionMessage(e))), caller.execute = TRUE)
    expect_identical(as.character(cleared), rep("prepared npreg objective state is not active",
      mpi.comm.size(1L)))
    valid <- do.call(npregbw, list(xdat = x, ydat = y, bws = c(44, 44),
      bwtype = type, regtype = "lp", nomad = TRUE, search.engine = "nomad",
      degree.select = "coordinate", degree.min = 0L, degree.max = 2L,
      degree.start = c(1L, 1L), degree.verify = FALSE, nmulti = 1L, itmax = 20L,
      nomad.opts = list(MAX_BB_EVAL = 30L)))
    expect_true(.np_nn_raw_objective_valid(valid$fval))
  }
})
