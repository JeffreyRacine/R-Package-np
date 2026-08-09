test_that("conditional tree objectives preserve observation identities", {
  old.options <- options(np.messages = FALSE, np.largeh = FALSE,
                         np.categorical.compress = FALSE)
  on.exit(options(old.options), add = TRUE)

  npcdensbw <- getFromNamespace("npcdensbw", "np")
  npcdistbw <- getFromNamespace("npcdistbw", "np")
  density.eval <- getFromNamespace(".npcdensbw_eval_only", "np")
  distribution.eval <- getFromNamespace(".npcdistbw_eval_only", "np")

  set.seed(202608093L)
  n <- 36L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = 1 + 2 * x$x + rnorm(n, sd = 0.3))
  cases <- expand.grid(
    family = c("density-ls", "density-ml", "distribution-ls"),
    bwtype = c("fixed", "generalized_nn", "adaptive_nn"),
    basis = c("lc", "raw1", "bernstein1"),
    stringsAsFactors = FALSE
  )

  evaluate <- function(spec, use.tree) {
    options(np.tree = use.tree)
    args <- list(
      xdat = x,
      ydat = y,
      bws = if (identical(spec$bwtype, "fixed")) c(0.24, 0.22) else c(14, 15),
      bwtype = spec$bwtype,
      bwmethod = if (identical(spec$family, "density-ml")) "cv.ml" else "cv.ls",
      regtype = if (identical(spec$basis, "lc")) "lc" else "lp",
      bandwidth.compute = FALSE
    )
    if (!identical(spec$basis, "lc")) {
      args$degree <- 1L
      args$bernstein.basis <- identical(spec$basis, "bernstein1")
    }
    if (startsWith(spec$family, "density")) {
      bw <- do.call(npcdensbw, args)
      density.eval(x, y, bw)$objective
    } else {
      bw <- do.call(npcdistbw, args)
      distribution.eval(x, y, bws = bw, do.full.integral = TRUE)$objective
    }
  }

  for (i in seq_len(nrow(cases))) {
    tree.off <- evaluate(cases[i, ], FALSE)
    tree.on <- evaluate(cases[i, ], TRUE)
    expect_true(is.finite(tree.off))
    expect_equal(tree.on, tree.off, tolerance = 1e-10,
                 info = paste(cases[i, ], collapse = "/"))
  }
})
