.native_beta_self_map_call <- function(bandwidth.mode,
                                       train.is.eval = 1L,
                                       num.eval = 4L) {
  constant <- function(name) get(name, envir = asNamespace("npRmpi"),
                                 inherits = FALSE)
  options <- list(
    num_obs_train = 3L,
    num_obs_eval = num.eval,
    num_uno = 0L,
    num_ord = 0L,
    num_con = 1L,
    int_LARGE_SF = constant("SF_ARB"),
    BANDWIDTH_reg_extern = bandwidth.mode,
    int_MINIMIZE_IO = constant("IO_MIN_TRUE"),
    kerneval = constant("CKER_COORDINATE"),
    ukerneval = constant("UKER_AIT"),
    okerneval = constant("OKER_WANG"),
    miss.ex = train.is.eval,
    leave.one.out = 0L,
    bandwidth.divide = 0L,
    mcv.numRow = 0L,
    wncol = 0L,
    yncol = 0L,
    int_do_tree = constant("DO_TREE_NO"),
    return.kernel.weights = 0L,
    permutation.operator = constant("OP_NOOP"),
    compute.score = 0L,
    compute.ocg = 0L,
    suppress.parallel = 0L,
    continuous.kernel.family = constant("CKER_FAMILY_BETA"),
    continuous.kernel.order = 2L,
    divide.returned.kernel.weights = 0L
  )

  .Call(
    "C_np_kernelsum",
    numeric(), numeric(), c(0.1, 0.5, 0.9), numeric(), numeric(),
    numeric(), numeric(), numeric(),
    if (bandwidth.mode == constant("BW_FIXED")) 0.2 else 2,
    numeric(), 0, 0L, as.integer(unlist(options, use.names = FALSE)), 1,
    as.integer(num.eval), 0L, 0L, 0, 1,
    PACKAGE = "npRmpi"
  )
}

test_that("beta native self maps reject inconsistent layouts", {
  expect_error(
    .native_beta_self_map_call(getFromNamespace("BW_FIXED", "npRmpi")),
    "invalid beta kernel-sum dimensions", fixed = TRUE
  )
  expect_error(
    .native_beta_self_map_call(getFromNamespace("BW_GEN_NN", "npRmpi")),
    "invalid beta kernel-sum dimensions", fixed = TRUE
  )
  expect_error(
    .native_beta_self_map_call(
      getFromNamespace("BW_GEN_NN", "npRmpi"),
      train.is.eval = 2L,
      num.eval = 3L
    ),
    "invalid beta kernel-sum dimensions", fixed = TRUE
  )
})

test_that("ordinary beta self maps retain all bandwidth topologies", {
  had.pool <- .mpi_pool_active()
  if (!spawn_mpi_slaves(1L))
    skip("MPI slave pool is unavailable")
  if (!had.pool)
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  x <- data.frame(x = c(0.02, 0.14, 0.37, 0.63, 0.86, 0.98))
  for (mode in c("fixed", "generalized_nn", "adaptive_nn")) {
    fit <- npksum(
      bws = if (identical(mode, "fixed")) 0.16 else 3,
      txdat = x,
      bwtype = mode,
      ckertype = "beta",
      ckerorder = 8,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      leave.one.out = TRUE,
      return.kernel.weights = TRUE
    )

    expect_true(all(is.finite(fit$ksum)))
    expect_true(all(diag(fit$kw) == 0))
  }
})
