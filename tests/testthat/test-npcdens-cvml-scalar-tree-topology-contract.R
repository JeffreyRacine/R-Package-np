test_that("scalar conditional CVML is invariant across canonical tree topologies", {
  old_options <- options(np.messages = FALSE, np.acceleration = FALSE)
  on.exit(options(old_options), add = TRUE)

  make_objective <- function(x, y, kernel, tree, engine) {
    options(np.tree = tree)
    p <- ncol(x)
    args <- list(
      xdat = x,
      ydat = y,
      # Public conditional bandwidth vectors are Y first, followed by X.
      bws = c(0.18, rep.int(0.31, p)),
      bandwidth.compute = FALSE,
      bwscaling = FALSE,
      bwtype = "fixed",
      bwmethod = "cv.ml",
      cxkertype = kernel
    )
    if (identical(engine, "lc")) {
      args$regtype <- "lc"
    } else {
      args$regtype <- "lp"
      args$basis <- "glp"
      args$degree <- rep.int(0L, p)
      args$bernstein.basis <- identical(engine, "bernstein")
    }
    bw <- do.call(npcdensbw, args)
    np:::.npcdensbw_eval_only(x, y, bw)$objective
  }

  for (p in 1:4) {
    set.seed(2026081100L + p)
    n <- 73L
    x <- as.data.frame(replicate(p, runif(n, -1, 1), simplify = FALSE))
    names(x) <- paste0("x", seq_len(p))
    score <- rowSums(x[, seq_len(min(p, 2L)), drop = FALSE])
    y <- data.frame(
      y = factor(ifelse(score + rnorm(n, sd = 0.4) > 0, "b", "a"))
    )

    for (kernel in c("epanechnikov", "uniform")) {
      values <- vapply(
        c("lc", "raw", "bernstein"),
        function(engine) {
          c(
            dense = make_objective(x, y, kernel, FALSE, engine),
            tree = make_objective(x, y, kernel, TRUE, engine),
            auto = make_objective(x, y, kernel, "auto", engine)
          )
        },
        numeric(3L)
      )
      scale <- max(1.0, abs(values))

      expect_identical(values[, "lc"], values[, "raw"])
      expect_identical(values[, "raw"], values[, "bernstein"])
      expect_identical(unname(values["tree", ]), unname(values["auto", ]))
      expect_lt(
        max(abs(values["dense", ] - values["tree", ])),
        512.0 * .Machine$double.eps * scale,
        label = paste("p", p, kernel)
      )
    }
  }
})
