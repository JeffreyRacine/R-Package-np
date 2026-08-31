wide_external_rhs <- function(n, width) {
  outer(seq_len(n), seq_len(width), function(i, j) {
    sin(i * (j + 1) / 17) + cos((i + 2) * j / 23) + i * j / 10000
  })
}

wide_external_bw <- function(x, y, type, bws) {
  npregbw(
    xdat = x, ydat = y, bws = bws,
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    bwtype = type, bwscaling = FALSE,
    regtype = "lp", degree = c(2L, 1L), degree.select = "manual",
    basis = "glp", bernstein.basis = FALSE,
    ckertype = "gaussian", ckerorder = 2L
  )
}

wide_external_expect_apply <- function(bw, tx, ex, width, tol = 1e-10) {
  H <- unclass(npreghat(
    bws = bw, txdat = tx, exdat = ex, output = "matrix"
  ))
  rhs <- wide_external_rhs(nrow(tx), width)
  observed <- npreghat(
    bws = bw, txdat = tx, exdat = ex, y = rhs, output = "apply"
  )

  expect_equal(unname(observed), unname(H %*% rhs), tolerance = tol)
}

wide_external_require_pool <- function() {
  if ("npRmpi" %in% loadedNamespaces()) {
    skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  }
  invisible(TRUE)
}

test_that("wide fixed GLP apply respects distinct external coordinates", {
  wide_external_require_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  n <- 48L
  x1 <- seq(-0.94, 0.94, length.out = n)
  x2 <- sin(seq(0.1, 2.8, length.out = n))
  tx <- data.frame(x1 = x1, x2 = x2)
  y <- sin(pi * x1) + 0.4 * x2 + seq_len(n) / 10000
  ex <- data.frame(
    x1 = rev(x1) * 0.91 + 0.015,
    x2 = rev(x2) * 0.87 - 0.025
  )
  bw <- wide_external_bw(tx, y, "fixed", c(0.46, 0.52))

  for (width in c(15L, 16L, 17L, 24L))
    wide_external_expect_apply(bw, tx, ex, width)
})

test_that("wide GLP ownership controls preserve copied, unequal-row, and GNN grids", {
  wide_external_require_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  n <- 48L
  x1 <- seq(-0.94, 0.94, length.out = n)
  x2 <- sin(seq(0.1, 2.8, length.out = n))
  tx <- data.frame(x1 = x1, x2 = x2)
  y <- sin(pi * x1) + 0.4 * x2 + seq_len(n) / 10000
  ex <- data.frame(
    x1 = rev(x1) * 0.91 + 0.015,
    x2 = rev(x2) * 0.87 - 0.025
  )
  copied <- data.frame(x1 = tx$x1 + 0, x2 = tx$x2 + 0)
  bw.fixed <- wide_external_bw(tx, y, "fixed", c(0.46, 0.52))
  bw.gnn <- wide_external_bw(tx, y, "generalized_nn", c(15L, 15L))

  wide_external_expect_apply(bw.fixed, tx, copied, 17L)
  wide_external_expect_apply(bw.fixed, tx, ex[-c(1L, 5L, 11L), ], 17L)
  wide_external_expect_apply(bw.gnn, tx, ex, 17L)
})

test_that("wide fixed GLP apply respects mixed external coordinates", {
  wide_external_require_pool()
  if ("npRmpi" %in% loadedNamespaces())
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  n <- 48L
  x1 <- seq(-0.94, 0.94, length.out = n)
  x2 <- sin(seq(0.1, 2.8, length.out = n))
  u <- factor(rep(c("a", "b", "c"), length.out = n))
  tx <- data.frame(x1 = x1, x2 = x2, u = u)
  y <- sin(pi * x1) + 0.4 * x2 + 0.15 * (u == "b")
  ex <- data.frame(
    x1 = rev(x1) * 0.91 + 0.015,
    x2 = rev(x2) * 0.87 - 0.025,
    u = factor(rep(c("c", "a", "b"), length.out = n), levels = levels(u))
  )
  bw <- npregbw(
    xdat = tx, ydat = y, bws = c(0.46, 0.52, 0.32),
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    bwtype = "fixed", bwscaling = FALSE,
    regtype = "lp", degree = c(2L, 1L), degree.select = "manual",
    basis = "glp", bernstein.basis = FALSE,
    ckertype = "gaussian", ckerorder = 2L
  )

  wide_external_expect_apply(bw, tx, ex, 17L)
})
