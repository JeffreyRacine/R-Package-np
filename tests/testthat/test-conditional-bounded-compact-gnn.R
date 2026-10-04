test_that("bounded positive compact GNN CVLS retains the whole-response criterion", {
  old <- options(np.messages = FALSE, np.extendednn = TRUE)
  on.exit(options(old), add = TRUE)
  set.seed(2701)
  x <- data.frame(x = runif(30))
  y <- data.frame(y = sin(3*x$x) + rnorm(30, sd = .3))
  evaluate <- function(args)
    get(".npcdensbw_eval_only", asNamespace("npRmpi"))(
      x, y, do.call(npcdensbw, args), invalid.penalty = "dbmax")$objective
  for (kernel in c("uniform", "epanechnikov"))
    for (reg in c("lc", "ll", "lp")) for (tree in c(FALSE, TRUE)) {
      options(np.tree = tree)
      args <- list(xdat = x, ydat = y, bws = c(8, 20),
                   bwmethod = "cv.ls", bwtype = "generalized_nn",
                   regtype = reg, cxkertype = kernel, bandwidth.compute = FALSE)
      if (reg == "lp") args <- c(args, list(degree = 2L, bernstein.basis = FALSE))
      reference <- evaluate(args)
      expect_true(is.finite(reference) && abs(reference) < 1e100)
      # X-row normalization cancels; response support and deleted radii do not change.
      for (bounds in list(list(cxkerbound = "range"),
                          list(cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1),
                          list(cxkerbound = "fixed", cxkerlb = -1e6, cxkerub = 1e6)))
        expect_lte(abs(evaluate(c(args, bounds)) - reference), 1e-9)
    }
})
