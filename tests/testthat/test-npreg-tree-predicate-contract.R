test_that("npreg bandwidth tree predicate rejects non-LC nearest-neighbor CV routes", {
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
    np:::.npregbw_tree_code(
      bw.lp.nn,
      ncon = bw.lp.nn$ncon,
      ncat = bw.lp.nn$nuno + bw.lp.nn$nord,
      regtype.engine = bw.lp.nn[["regtype.engine"]],
      basis.engine = bw.lp.nn[["basis"]],
      degree.engine = bw.lp.nn[["degree"]]
    ),
    np:::DO_TREE_NO
  )
  expect_identical(
    np:::.npregbw_tree_code(
      bw.lc.nn,
      ncon = bw.lc.nn$ncon,
      ncat = bw.lc.nn$nuno + bw.lc.nn$nord,
      regtype.engine = bw.lc.nn[["regtype.engine"]],
      basis.engine = bw.lc.nn[["basis"]],
      degree.engine = bw.lc.nn[["degree"]]
    ),
    np:::DO_TREE_YES
  )
  expect_identical(
    np:::.npregbw_tree_code(
      bw.lp.fixed,
      ncon = bw.lp.fixed$ncon,
      ncat = bw.lp.fixed$nuno + bw.lp.fixed$nord,
      regtype.engine = bw.lp.fixed[["regtype.engine"]],
      basis.engine = bw.lp.fixed[["basis"]],
      degree.engine = bw.lp.fixed[["degree"]]
    ),
    np:::DO_TREE_YES
  )

  tree.obj <- np:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective
  options(np.tree = FALSE)
  serial.obj <- np:::.npregbw_eval_only(x, y, bw.lp.fixed)$objective

  expect_equal(tree.obj, serial.obj, tolerance = 1e-12)
})

test_that("wide fixed LP AUTO is conservative without overriding explicit tree requests", {
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
    np:::.npregbw_tree_code(
      bw,
      ncon = bw$ncon,
      ncat = bw$nuno + bw$nord,
      regtype.engine = bw[["regtype.engine"]],
      basis.engine = bw[["basis"]],
      degree.engine = bw[["degree"]]
    )
  }

  old_opts <- options(np.tree = "auto", np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)
  expect_identical(route(), np:::DO_TREE_NO)

  options(np.tree = TRUE)
  expect_identical(route(), np:::DO_TREE_YES)
  tree.obj <- np:::.npregbw_eval_only(x, y, bw)$objective

  options(np.tree = FALSE)
  expect_identical(route(), np:::DO_TREE_NO)
  dense.obj <- np:::.npregbw_eval_only(x, y, bw)$objective
  expect_equal(tree.obj, dense.obj, tolerance = 1e-10)
})
