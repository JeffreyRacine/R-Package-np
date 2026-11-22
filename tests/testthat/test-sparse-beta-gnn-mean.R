test_that("fully pooled first-degree beta GNN means reproduce the weighted-design row", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  v <- readRDS(test_path("fixtures", "sparse-beta-gnn.rds"))
  # Seven boundary queries were checked against an 80-decimal reference.
  # The upper endpoint has the largest error with normal-equation arithmetic.
  for (bernstein in c(FALSE, TRUE)) {
    bw <- npregbw(xdat=v$x, ydat=v$y, bws=c(5, .5, .5), bandwidth.compute=FALSE,
                  regtype="lp", degree=1L, basis="glp", bernstein.basis=bernstein,
                  bwtype="generalized_nn", ckertype="beta", ckerorder=2L,
                  ckerbound="fixed", ckerlb=v$bounds[1], ckerub=v$bounds[2])
    row <- suppressWarnings(npreghat(bw, txdat=v$x, exdat=v$query, output="matrix"))
    fit <- suppressWarnings(npreg(bw, txdat=v$x, tydat=v$y, exdat=v$query, se=FALSE))
    expect_lt(max(abs(row - v$pooled_row)), 1e-8)
    expect_lt(abs(fit$mean - sum(v$pooled_row * v$y)), 1e-8)
  }
})
