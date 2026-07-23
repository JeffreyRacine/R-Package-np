test_that("Epanechnikov tree derivative-weight buffers match dense rows", {
  old <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(2026072201L)
  n <- 97L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(x$x1) + cos(x$x2)
  xe <- x[seq.int(2L, n, by = 2L), , drop = FALSE]

  for (spec in list(
    list(type = "fixed", order = 4L, bws = c(0.19, 0.31)),
    list(type = "generalized_nn", order = 6L, bws = c(19L, 27L))
  )) {
    options(np.tree = FALSE)
    bw <- npregbw(
      xdat = x, ydat = y, regtype = "lc", bwtype = spec$type,
      ckertype = "epanechnikov", ckerorder = spec$order,
      bws = spec$bws, bandwidth.compute = FALSE
    )

    invoke <- function(tree) {
      options(np.tree = tree)
      np:::npksum.default(
        bws = bw, txdat = x, exdat = xe,
        bandwidth.divide = TRUE,
        return.kernel.weights = TRUE,
        return.derivative.kernel.weights = TRUE,
        permutation.operator = "derivative"
      )
    }

    dense <- invoke(FALSE)
    tree <- invoke(TRUE)
    automatic <- invoke("auto")
    repeated <- invoke(TRUE)

    expect_identical(tree$p.kw, dense$p.kw)
    expect_identical(automatic$p.kw, dense$p.kw)
    expect_identical(repeated$p.kw, tree$p.kw)
    expect_true(all(is.finite(tree$p.kw)))

    for (dimension in seq_len(ncol(x))) {
      operator <- rep.int("normal", ncol(x))
      operator[[dimension]] <- "derivative"
      options(np.tree = FALSE)
      direct <- np:::npksum.default(
        bws = bw, txdat = x, exdat = xe,
        bandwidth.divide = TRUE,
        return.kernel.weights = TRUE,
        operator = operator
      )
      expect_identical(tree$p.kw[, , dimension], direct$kw)
    }

    for (name in c("ksum", "kw", "p.ksum"))
      expect_equal(tree[[name]], dense[[name]], tolerance = 1e-13)
  }
})
