test_that("conditional-density degree search recovers only automatic invalid NN starts", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE, np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = c(rep(0, 32), sin((1:16)*sqrt(2))+(1:16)/32),
    x2 = c(rep(0, 32), sin((1:16)*sqrt(3))+(1:16)/48))
  y <- data.frame(y = sin((1:48)/3)+(1:48)/90)
  for (method in c("cv.ml", "cv.ls")) for (type in c("generalized_nn", "adaptive_nn")) {
    controls <- list(xdat = x, ydat = y, bwtype = type, bwmethod = method,
      regtype = "lp", nomad = TRUE, search.engine = "nomad",
      degree.select = "coordinate", degree.min = 0L, degree.max = 2L,
      degree.start = c(1L,1L), degree.verify = FALSE,
      nmulti = 1L, itmax = 20L, powell.remin = FALSE,
      nomad.opts = list(MAX_BB_EVAL = 30L))
    bw <- do.call(npcdensbw, controls)
    restarts <- bw$nomad.restart.results
    expect_length(restarts, 2L)
    expect_true(isTRUE(restarts[[2L]]$recovery))
    expect_identical(restarts[[2L]]$recovery_witness$evaluations, 3L)
    expect_equal(tail(restarts[[2L]]$start, 2L), c(1,1), tolerance = 0)
    expect_identical(bw$nomad.best.restart, 2L)
    expect_true(restarts[[1L]]$native$compiled_callback_calls > 0L)
    expect_equal(bw$num.feval, 4 + sum(vapply(restarts,
      function(z) z$native$total_num.feval, numeric(1L))), tolerance = 0)
    raw <- .npcdensbw_eval_only(xdat = x, ydat = y, bws = bw, invalid.penalty = "dbmax")
    expect_true(.np_nn_raw_objective_valid(raw$objective))
    expect_equal(bw$fval, raw$objective, tolerance = 2e-12)
    expect_true(all(bw$degree %in% 0:2))
    controls$bws <- rep(7, 3)
    controls$nomad.opts <- list(MAX_BB_EVAL = 1L)
    expect_error(do.call(npcdensbw, controls), "did not return a raw-valid solution")
    controls$bws <- rep(if (type == "adaptive_nn") 46 else 47, 3)
    explicit <- do.call(npcdensbw, controls)
    expect_false(any(vapply(explicit$nomad.restart.results,
      function(z) isTRUE(z$recovery), logical(1L))))
    expect_true(.np_nn_raw_objective_valid(explicit$fval))
  }
})

test_that("prepared conditional-density degree searches transport extended NN bounds", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  observed <- new.env(parent=emptyenv())
  original <- .npRmpi_bcast_cmd_expr
  testthat::local_mocked_bindings(
    .npRmpi_bcast_cmd_expr=function(expr,...) {
      if(is.call(expr) && is.call(expr[[1L]]) &&
         identical(expr[[1L]][[1L]],as.name("get")) &&
         identical(expr[[1L]][[2L]],"npRmpiPreparedObjectiveSearchConditionalDensity")) {
        observed$setup <- expr[[3L]]
        observed$upper <- expr[[9L]]
      }
      original(expr,...)
    }, .package="npRmpi")
  n <- 24L
  x <- data.frame(x=sin(seq_len(n)*sqrt(2))+seq_len(n)/n,
                  u=factor(rep(c("a","b"),n/2)))
  y <- cos(seq_len(n)/3)+x$x
  for(type in c("generalized_nn","adaptive_nn")) for(extended in c(FALSE,TRUE)) {
    options(np.extendednn=extended)
    observed$setup <- NULL
    bw <- npcdensbw(xdat=x,ydat=y,bwtype=type,regtype="lp",nomad=TRUE,
      search.engine="nomad",nmulti=1L,nomad.nmulti=1L,degree.min=0L,
      degree.max=2L,degree.start=1L,degree.verify=FALSE,random.seed=42L,
      nomad.opts=list(MAX_BB_EVAL=4L))
    expect_type(observed$setup,"list")
    expect_equal(observed$setup$cont_extendednn_upper,
      observed$upper[seq_along(observed$setup$cont_flat)],tolerance=0)
    expect_true(is.finite(bw$fval))
  }
})
