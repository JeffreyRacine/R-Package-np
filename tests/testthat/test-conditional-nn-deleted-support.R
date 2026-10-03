# Keep this proof on the actual objective owner, including pooled MPI runs.
conditional_nn_support_value <- function(X, Y, b, cdf = FALSE) {
  pkg <- getNamespaceName(environment(npcdensbw))
  name <- if (cdf) ".npcdistbw_eval_only" else ".npcdensbw_eval_only"
  if (pkg != "npRmpi" || .mpi_suite_local_mode_owned())
    return(get(name, asNamespace(pkg))(xdat = X, ydat = Y, bws = b, invalid.penalty = "dbmax")$objective)
  command <- substitute(local({
    old <- options(np.tree = TREE)
    on.exit(options(old), add = TRUE)
    get(NAME, asNamespace("npRmpi"))(xdat = X, ydat = Y, bws = B,
      invalid.penalty = "dbmax", force.local = FALSE)$objective
  }), list(X = X, Y = Y, B = b, NAME = name, TREE = getOption("np.tree")))
  get(".npRmpi_bcast_cmd_expr", asNamespace(pkg))(
    command, comm = 1L, caller.execute = TRUE)
}

test_that("compact NN CVML and CDF reject insufficient deleted designs in all units", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2701)
  x1 <- rnorm(60); x2 <- rnorm(60)
  Y <- data.frame(y = sin(x1) + .5*x2 + rnorm(60, sd = .3))
  for (tree in c(FALSE, TRUE)) for (scale in c(1, 4, .25)) {
    options(np.tree = tree)
    X <- data.frame(x1 = scale*x1, x2 = scale*x2)
    # Independent necessary-support bound for the six quadratic terms.
    counts <- vapply(seq_len(60), function(i) {
      keep <- setdiff(seq_len(60), i)
      a <- abs(X$x1[keep] - X$x1[i]); b <- abs(X$x2[keep] - X$x2[i])
      sum(a < sort(a)[15L] & b < sort(b)[15L])
    }, 0L)
    expect_true(any(counts < 6L))
    for (cdf in c(FALSE, TRUE)) for (k in c(15L, 30L)) {
      make <- if (cdf) npcdistbw else npcdensbw
      b <- make(xdat = X, ydat = Y, bws = c(15, k, k),
        bwtype = "generalized_nn", regtype = "lp", degree = c(2L, 2L),
        cxkertype = "uniform", bandwidth.compute = FALSE)
      value <- conditional_nn_support_value(X, Y, b, cdf)
      if (k == 15L) expect_identical(value, if (cdf) .Machine$double.xmax else -.Machine$double.xmax)
      else {
        expect_true(is.finite(value) && abs(value) < 1e100)
        if (!cdf) expect_equal(value, -26.8959180024696, tolerance = 1e-10)
      }
    }
  }
})

test_that("adaptive and bounded conditional owners retain deleted-support admission", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2701)
  X <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  Y <- data.frame(y = sin(X$x1) + .5*X$x2 + rnorm(60, sd = .3))
  set.seed(1)
  B <- data.frame(x1 = runif(80), x2 = runif(80))
  BY <- data.frame(y = sin(3*B$x1) + B$x2 + rnorm(80, sd = .3))
  for (tree in c(FALSE, TRUE)) for (cdf in c(FALSE, TRUE)) {
    options(np.tree = tree)
    make <- if (cdf) npcdistbw else npcdensbw
    b <- make(xdat = X, ydat = Y, bws = c(17, 15, 15),
      bwtype = "adaptive_nn", regtype = "lp", degree = c(2L, 2L),
      cxkertype = "uniform", bandwidth.compute = FALSE)
    expect_identical(conditional_nn_support_value(X, Y, b, cdf),
                     if (cdf) .Machine$double.xmax else -.Machine$double.xmax)
    b <- make(xdat = B, ydat = BY, bws = c(20, 20, 20),
      bwtype = "generalized_nn", regtype = "lp", degree = c(2L, 2L),
      cxkertype = "uniform", cxkerbound = "fixed",
      cxkerlb = c(0, 0), cxkerub = c(1, 1), bandwidth.compute = FALSE)
    expect_identical(conditional_nn_support_value(B, BY, b, cdf),
                     if (cdf) .Machine$double.xmax else -.Machine$double.xmax)
  }
})

test_that("categorical predictors do not bypass the compact NN support check", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(1)
  x <- rnorm(100); f <- factor(sample(letters[1:3], 100, TRUE))
  X <- data.frame(x = x, f = f)
  Y <- data.frame(y = sin(x) + as.integer(f)/2 + rnorm(100, sd = .3))
  for (tree in c(FALSE, TRUE)) for (cdf in c(FALSE, TRUE)) {
    options(np.tree = tree)
    make <- if (cdf) npcdistbw else npcdensbw
    b <- make(xdat = X, ydat = Y, bws = c(15, 2, .5),
      bwtype = "generalized_nn", regtype = "lp", degree = 2L,
      cxkertype = "epanechnikov", bandwidth.compute = FALSE)
    expect_identical(conditional_nn_support_value(X, Y, b, cdf),
                     if (cdf) .Machine$double.xmax else -.Machine$double.xmax)
  }
})
