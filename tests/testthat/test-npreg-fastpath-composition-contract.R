library(npRmpi)

test_that("npregbw fixed LP tree predicate follows bandwidth-selection policy", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)
  old_dispatch <- options(npRmpi.autodispatch = FALSE)
  on.exit(options(old_dispatch), add = TRUE)

  set.seed(20260518)
  n <- 96L
  x <- as.data.frame(matrix(runif(n * 2L), nrow = n))
  names(x) <- c("x1", "x2")
  y <- rowSums(sin(2 * pi * as.matrix(x))) + rnorm(n, sd = 0.1)

  bw.lp <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lp",
    degree = rep.int(3L, ncol(x)),
    bernstein.basis = TRUE,
    bwmethod = "cv.ls",
    bws = rep(0.35, ncol(x)),
    bandwidth.compute = FALSE
  )
  bw.lc <- npregbw(
    xdat = x,
    ydat = y,
    regtype = "lc",
    bwmethod = "cv.ls",
    bws = rep(0.35, ncol(x)),
    bandwidth.compute = FALSE
  )

  old_opts <- options(np.tree = TRUE, np.categorical.compress = FALSE)
  on.exit(options(old_opts), add = TRUE)

  expect_equal(
    npRmpi:::.npregbw_tree_code(bw.lp, ncon = bw.lp$ncon, ncat = bw.lp$nuno + bw.lp$nord),
    npRmpi:::DO_TREE_YES
  )
  expect_equal(
    npRmpi:::.npregbw_tree_code(bw.lc, ncon = bw.lc$ncon, ncat = bw.lc$nuno + bw.lc$nord),
    npRmpi:::DO_TREE_YES
  )

  options(np.tree = FALSE)
  expect_equal(
    npRmpi:::.npregbw_tree_code(bw.lp, ncon = bw.lp$ncon, ncat = bw.lp$nuno + bw.lp$nord),
    npRmpi:::DO_TREE_NO
  )
})
