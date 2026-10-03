# Evaluate through the actual rank-owned route when this test has an MPI pool.
cgnn_multivariate_evaluate <- function(X, Y, b, invalid.penalty) {
  pkg <- getNamespaceName(environment(npcdensbw))
  if (pkg != "npRmpi" || .mpi_suite_local_mode_owned())
    return(get(".npcdensbw_eval_only", asNamespace(pkg))(X, Y, b,
      invalid.penalty = invalid.penalty))
  command <- substitute(local({
    old <- options(np.tree = TREE)
    on.exit(options(old), add = TRUE)
    get(".npcdensbw_eval_only", asNamespace("npRmpi"))(
      X, Y, B, invalid.penalty = P, force.local = FALSE)
  }), list(X = X, Y = Y, B = b, P = invalid.penalty,
           TREE = getOption("np.tree")))
  get(".npRmpi_bcast_cmd_expr", asNamespace(pkg))(
    command, comm = 1L, caller.execute = TRUE)
}

test_that("multivariate compact-X GNN CVLS rejects insufficient deleted support", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  pkg <- getNamespaceName(environment(npcdensbw))
  evaluate <- cgnn_multivariate_evaluate
  set.seed(2701)
  x1 <- rnorm(60); x2 <- rnorm(60)
  Y <- data.frame(y = sin(x1) + .5*x2 + rnorm(60, sd = .3))
  for (tree in c(FALSE, TRUE)) for (bernstein in c(FALSE, TRUE)) {
    options(np.tree = tree)
    for (scale in c(1, 4, .25)) {
      X <- data.frame(x1 = scale*x1, x2 = scale*x2)
      b <- npcdensbw(xdat = X, ydat = Y, bws = c(17, 16, 11),
        bwtype = "generalized_nn", bwmethod = "cv.ls", regtype = "lp",
        degree = c(2L, 2L), bernstein.basis = bernstein,
        cxkertype = "uniform", bandwidth.compute = FALSE)
      # Independent strict-support census: fewer than six donors bounds rank
      # below the total-degree quadratic basis width, in every deleted fold.
      counts <- vapply(seq_len(60), function(i) {
        keep <- setdiff(seq_len(60), i)
        a <- abs(X$x1[keep] - X$x1[i]); b <- abs(X$x2[keep] - X$x2[i])
        sum(a < sort(a)[16L] & b < sort(b)[11L])
      }, 0L)
      expect_true(all(counts < 6L))
      expect_identical(evaluate(X, Y, b, invalid.penalty = "dbmax")$objective,
                       -.Machine$double.xmax)
      valid <- npcdensbw(xdat = X, ydat = Y, bws = c(17, 30, 30),
        bwtype = "generalized_nn", bwmethod = "cv.ls", regtype = "lp",
        degree = c(2L, 2L), bernstein.basis = bernstein,
        cxkertype = "uniform", bandwidth.compute = FALSE)
      value <- evaluate(X, Y, valid, invalid.penalty = "dbmax")$objective
      expect_equal(value, 0.368663925437, tolerance = 1e-10)
    }
  }
})

test_that("multivariate support counts complete basis rows instead of observations", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  evaluate <- cgnn_multivariate_evaluate
  # Four distinct locations repeated twelve times: observation count is ample,
  # but even the full support has rank at most four for six quadratic terms.
  X <- data.frame(x1 = rep(c(-1, -1, 1, 1), 12),
                  x2 = rep(c(-1, 1, -1, 1), 12))
  Y <- data.frame(y = sin(seq_len(48)))
  for (tree in c(FALSE, TRUE)) for (kernel in c("uniform", "epanechnikov")) {
    options(np.tree = tree)
    b <- npcdensbw(xdat = X, ydat = Y, bws = c(16, 48, 48),
      bwtype = "generalized_nn", bwmethod = "cv.ls", regtype = "lp",
      degree = c(2L, 2L), cxkertype = kernel, bandwidth.compute = FALSE)
    expect_identical(evaluate(X, Y, b, invalid.penalty = "dbmax")$objective,
                     -.Machine$double.xmax)
  }
})

test_that("zero-degree coordinates cannot inflate the structural rank bound", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  evaluate <- cgnn_multivariate_evaluate
  X <- data.frame(x1 = rep(c(-1, 1), 24), x2 = seq_len(48)/48)
  Y <- data.frame(y = cos(seq_len(48)))
  for (tree in c(FALSE, TRUE)) {
    options(np.tree = tree)
    b <- npcdensbw(xdat = X, ydat = Y, bws = c(16, 48, 48),
      bwtype = "generalized_nn", bwmethod = "cv.ls", regtype = "lp",
      degree = c(2L, 0L), cxkertype = "uniform", bandwidth.compute = FALSE)
    expect_identical(evaluate(X, Y, b, invalid.penalty = "dbmax")$objective,
                     -.Machine$double.xmax)
  }
})
