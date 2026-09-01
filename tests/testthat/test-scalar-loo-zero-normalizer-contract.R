.scalar_loo_mpi_observe <- function(value) {
  expr <- substitute(local({
    x.local <- data.frame(x = seq(0, 1, length.out = 8L))
    y.local <- sin(seq_len(8L))
    binary.local <- rep(c(0, 1), 4L)
    make.bw <- function(response) npregbw(
      xdat = x.local, ydat = response, bws = VALUE,
      bandwidth.compute = FALSE, bwtype = "fixed", bwscaling = FALSE,
      ckertype = "epanechnikov", ckerorder = 2L,
      bwmethod = "cv.ls", regtype = "lc"
    )
    reg.core <- get(".npregbw_call_fixed_degree_core",
                    envir = asNamespace("npRmpi"), inherits = FALSE)
    check.core <- get(".nplsqreg_call_fixed_degree_core",
                      envir = asNamespace("npRmpi"), inherits = FALSE)
    cvls.bw <- make.bw(y.local)
    cvks.bw <- make.bw(binary.local)
    run.reg <- function(response, bandwidth, objective, penalty) reg.core(
      xdat = x.local, ydat = response, bws = bandwidth,
      invalid.penalty = penalty, penalty.multiplier = 10,
      eval.only = TRUE, objective = objective
    )
    run.check <- function(penalty) check.core(
      xdat = x.local, ydat = y.local, scale = rep(1, 8L), tau = 0.5,
      bws = cvls.bw, delta = 0.5,
      delta.bounds = c(1e-4, 1 - 1e-4),
      opt.args = list(invalid.penalty = penalty, penalty.multiplier = 10),
      bandwidth.compute = FALSE
    )
    run <- function(penalty) list(
      CVLS = run.reg(y.local, cvls.bw, "ls", penalty),
      CVCHECK = run.check(penalty),
      CVKS = run.reg(binary.local, cvks.bw, "ks", penalty)
    )
    list(raw = run("dbmax"), mapped = run("baseline"))
  }), list(VALUE = as.double(value)))
  getFromNamespace(".npRmpi_bcast_cmd_expr", "npRmpi")(
    expr, comm = 1L, caller.execute = TRUE
  )
}

test_that("npRmpi scalar normalized LOO rejects an exact zero normalizer", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  observed <- .scalar_loo_mpi_observe(1e-8)
  for (mode in names(observed$raw)) {
    expect_identical(as.double(observed$raw[[mode]]$objective),
                     .Machine$double.xmax, info = mode)
    expect_identical(as.double(observed$raw[[mode]]$invalid.history[[1L]]),
                     1, info = mode)
    expect_true(is.finite(observed$mapped[[mode]]$objective), info = mode)
    expect_lt(observed$mapped[[mode]]$objective, .Machine$double.xmax)
  }
})

test_that("npRmpi scalar normalized LOO preserves supported rows", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  observed <- .scalar_loo_mpi_observe(0.35)
  for (mode in names(observed$raw)) {
    expect_identical(as.double(observed$raw[[mode]]$invalid.history[[1L]]),
                     0, info = mode)
    expect_identical(as.double(observed$mapped[[mode]]$objective),
                     as.double(observed$raw[[mode]]$objective), info = mode)
  }
})
