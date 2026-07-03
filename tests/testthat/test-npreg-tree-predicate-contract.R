test_that("npreg bandwidth tree predicate rejects non-LC nearest-neighbor CV routes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(3301)
  n <- 28L
  x <- data.frame(x1 = seq(-1, 1, length.out = n), x2 = runif(n))
  y <- sin(x$x1) + 0.25 * x$x2 + rnorm(n, sd = 0.03)

  bw.lp.nn <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = rep(1L, ncol(x)),
    bwmethod = "cv.aic",
    bwtype = "generalized_nn",
    bws = rep(7, ncol(x)),
    bandwidth.compute = FALSE
  )
  bw.lc.nn <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwmethod = "cv.aic",
    bwtype = "generalized_nn",
    bws = rep(7, ncol(x)),
    bandwidth.compute = FALSE
  )
  bw.lp.fixed <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = rep(1L, ncol(x)),
    bwmethod = "cv.aic",
    bwtype = "fixed",
    bws = rep(0.45, ncol(x)),
    bandwidth.compute = FALSE
  )

  old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  expect_identical(
    npRmpi:::.npregbw_tree_code(
      bw.lp.nn,
      ncon = bw.lp.nn$ncon,
      ncat = bw.lp.nn$nuno + bw.lp.nn$nord
    ),
    npRmpi:::DO_TREE_NO
  )
  expect_identical(
    npRmpi:::.npregbw_tree_code(
      bw.lc.nn,
      ncon = bw.lc.nn$ncon,
      ncat = bw.lc.nn$nuno + bw.lc.nn$nord
    ),
    npRmpi:::DO_TREE_YES
  )
  expect_identical(
    npRmpi:::.npregbw_tree_code(
      bw.lp.fixed,
      ncon = bw.lp.fixed$ncon,
      ncat = bw.lp.fixed$nuno + bw.lp.fixed$nord
    ),
    npRmpi:::DO_TREE_YES
  )

  tree.obj <- npRmpi:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective
  options(np.tree = FALSE)
  serial.obj <- npRmpi:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective

  expect_equal(tree.obj, serial.obj, tolerance = 1e-12)
})
