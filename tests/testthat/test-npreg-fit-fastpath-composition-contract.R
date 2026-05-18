library(np)

test_that("npreg non-LC fit stream does not route through legacy tree rows", {
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
    np:::.npreg_fit_tree_code(bw.lp, ncon = bw.lp$ncon, ncat = bw.lp$nuno + bw.lp$nord),
    np:::DO_TREE_NO
  )
  expect_equal(
    np:::.npreg_fit_tree_code(bw.lc, ncon = bw.lc$ncon, ncat = bw.lc$nuno + bw.lc$nord),
    np:::DO_TREE_YES
  )

  fit.off <- npreg(bws = bw.lp, txdat = x, tydat = y)
  options(np.tree = FALSE)
  fit.base <- npreg(bws = bw.lp, txdat = x, tydat = y)
  expect_equal(fitted(fit.off), fitted(fit.base), tolerance = 1e-12)

  hat.off <- npreghat(bws = bw.lp, txdat = x, y = y, output = "apply")
  options(np.tree = TRUE)
  hat.on <- npreghat(bws = bw.lp, txdat = x, y = y, output = "apply")
  expect_equal(hat.on, hat.off, tolerance = 1e-12)
})
