# These are deliberately rank-local numerical contracts. The separate
# guarded-count contract below exercises actual collective ownership.
# CF195: deterministic column-order witnesses, not fitted-value transcripts.
test_that("conditional fixed LP objectives are stable under predictor permutation", {
  old <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  for (n in c(40L, 100L)) {
    set.seed(if (n == 40L) 2401L else 24101L)
    x1 <- rt(n, 2); x2 <- rexp(n)
    X <- data.frame(x1 = x1, x2 = x2)
    Y <- data.frame(y = 10 * x1 + 2 * x2 + rnorm(n, sd = 3))
    for (family in c("dens_ls", "dens_ml", "dist_ls")) {
      h2 <- if (family == "dens_ls") .4 else .3162829
      values <- vapply(list(1:2, 2:1), function(perm) {
        xd <- X[, perm, drop = FALSE]
        args <- list(xdat = xd, ydat = Y,
          bws = c(8.3936718, c(.92641032, h2)[perm]), bwscaling = FALSE,
          bandwidth.compute = FALSE, regtype = "ll")
        if (family != "dist_ls")
          args$bwmethod <- if (family == "dens_ls") "cv.ls" else "cv.ml"
        constructor <- if (family == "dist_ls") npcdistbw else npcdensbw
        bw <- do.call(constructor, args)
        evaluator <- getFromNamespace(if (family == "dist_ls")
          ".npcdistbw_eval_only" else ".npcdensbw_eval_only", "npRmpi")
        evaluator(xdat = xd, ydat = Y, bws = bw,
                  invalid.penalty = "dbmax", force.local = TRUE)$objective
      }, 0.)
      expect_true(all(is.finite(values)), info = paste(n, family))
      expect_true(abs(diff(values)) / (1 + abs(values[1])) < 1e-10,
                info = paste(n, family))
    }
  }
})

test_that("global and ordinary conditional LP deletion agree across degrees", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  for (degree in c(1L, 2L, 5L)) {
    set.seed(301910L + degree)
    X <- data.frame(x = runif(100, -.8, .8))
    Y <- data.frame(y = sin(X$x) + rnorm(100, sd = .3))
    bw <- npcdensbw(xdat = X, ydat = Y, bws = c(.4, 2),
      bwscaling = FALSE, bandwidth.compute = FALSE, regtype = "lp",
      degree = degree, basis = "glp", bernstein.basis = FALSE,
      cxkertype = "uniform", cykertype = "gaussian", bwmethod = "cv.ml")
    values <- vapply(c(FALSE, TRUE), function(global) {
      options(np.largeh = global)
      getFromNamespace(".npcdensbw_eval_only", "npRmpi")(
        X, Y, bw, invalid.penalty = "dbmax", force.local = TRUE)$objective
    }, 0.)
    expect_true(all(is.finite(values)), info = paste("degree", degree))
    expect_true(abs(diff(values)) / (1 + abs(values[1])) < 1e-10,
              info = paste("degree", degree))
  }
})
