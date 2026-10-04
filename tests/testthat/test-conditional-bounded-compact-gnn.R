test_that("bounded positive compact GNN CVLS retains the whole-response criterion", {
  old <- options(np.messages = FALSE, np.extendednn = TRUE,
                 np.tree = getOption("np.tree"))
  on.exit(options(old), add = TRUE)
  set.seed(2701)
  x <- data.frame(x = runif(30))
  y <- data.frame(y = sin(3*x$x) + rnorm(30, sd = .3))
  evaluate <- function(args)
    get(".npcdensbw_eval_only", asNamespace("np"))(
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


test_that("bounded higher-order Epanechnikov GNN CVLS uses deleted whole-line integrals", {
  old <- options(np.messages = FALSE, np.extendednn = TRUE,
                 np.tree = getOption("np.tree"))
  on.exit(options(old), add = TRUE)
  set.seed(2701)
  x <- data.frame(x = runif(60))
  y <- data.frame(y = sin(3*x$x) + rnorm(60, sd = .3))
  # Independently integrated deleted-sample mixtures over the whole real line,
  # splitting at each NN-radius change. These are criterion values, not fits
  # obtained from the package or from the former bounded objective.
  reference <- c(`2` = .715905433065280, `4` = .765391073036826,
                 `6` = .707738424778271, `8` = .601157485736127)
  ns <- asNamespace(getNamespaceName(environment(npcdensbw)))
  for (order in c(2L, 4L, 6L, 8L)) for (tree in list(FALSE, TRUE, "auto")) {
    options(np.tree = tree)
    args <- list(xdat = x, ydat = y, bws = c(12, 18),
      bwmethod = "cv.ls", bwtype = "generalized_nn", regtype = "lc",
      cxkertype = "epanechnikov", cxkerorder = order, bandwidth.compute = FALSE)
    for (bounds in list(list(cxkerbound = "range"),
                        list(cxkerbound = "fixed", cxkerlb = -1e6, cxkerub = 1e6))) {
      b <- do.call(npcdensbw, c(args, bounds))
      value <- get(".npcdensbw_eval_only", ns)(x, y, b,
                                               invalid.penalty = "dbmax")$objective
      expect_true(is.finite(value) && abs(value) < 1e100)
      expect_lte(abs(value - reference[as.character(order)]), 1e-9)
    }
  }
})
