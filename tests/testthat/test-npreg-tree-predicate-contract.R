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
      ncat = bw.lp.nn$nuno + bw.lp.nn$nord,
      regtype.engine = bw.lp.nn[["regtype.engine"]]
    ),
    npRmpi:::DO_TREE_NO
  )
  expect_identical(
    npRmpi:::.npregbw_tree_code(
      bw.lc.nn,
      ncon = bw.lc.nn$ncon,
      ncat = bw.lc.nn$nuno + bw.lc.nn$nord,
      regtype.engine = bw.lc.nn[["regtype.engine"]]
    ),
    npRmpi:::DO_TREE_YES
  )
  expect_identical(
    npRmpi:::.npregbw_tree_code(
      bw.lp.fixed,
      ncon = bw.lp.fixed$ncon,
      ncat = bw.lp.fixed$nuno + bw.lp.fixed$nord,
      regtype.engine = bw.lp.fixed[["regtype.engine"]]
    ),
    npRmpi:::DO_TREE_YES
  )

  tree.obj <- npRmpi:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective
  options(np.tree = FALSE)
  serial.obj <- npRmpi:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective

  expect_equal(tree.obj, serial.obj, tolerance = 1e-12)
})

test_that("wide fixed LP AUTO retains the profitable MPI tree owner", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(3302)
  n <- 48L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x$x1) + x$x2 + rnorm(n, sd = 0.03)
  bw <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = c(2L, 2L),
    ckertype = "epanechnikov",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE
  )
  route <- function() {
    npRmpi:::.npregbw_tree_code(
      bw,
      ncon = bw$ncon,
      ncat = bw$nuno + bw$nord,
      regtype.engine = bw[["regtype.engine"]]
    )
  }

  old_opts <- options(np.tree = "auto", np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)
  expect_identical(route(), npRmpi:::DO_TREE_YES)
  auto.obj <- npRmpi:::.npregbw_eval_only(x, y, bw)$objective

  options(np.tree = TRUE)
  expect_identical(route(), npRmpi:::DO_TREE_YES)
  tree.obj <- npRmpi:::.npregbw_eval_only(x, y, bw)$objective
  expect_identical(auto.obj, tree.obj)

  options(np.tree = FALSE)
  expect_identical(route(), npRmpi:::DO_TREE_NO)
  dense.obj <- npRmpi:::.npregbw_eval_only(x, y, bw)$objective
  expect_equal(tree.obj, dense.obj, tolerance = 1e-10)
})
