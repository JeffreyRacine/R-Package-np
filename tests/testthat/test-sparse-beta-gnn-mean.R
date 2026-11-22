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

test_that("the quadratic beta GNN prepared mean retains the accepted row", {
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  v <- readRDS(test_path("fixtures", "sparse-beta-gnn.rds"))
  q <- v$quadratic
  x <- data.frame(x=(v$x$x-q$offset)/q$scale)
  at <- data.frame(x=(v$query$x-q$offset)/q$scale)
  bounds <- (v$bounds-q$offset)/q$scale
  for (bernstein in c(FALSE, TRUE)) {
    bw <- npregbw(xdat=x, ydat=v$y, bws=12, bandwidth.compute=FALSE,
                  regtype="lp", degree=2L, basis="glp", bernstein.basis=bernstein,
                  bwtype="generalized_nn", ckertype="beta", ckerorder=2L,
                  ckerbound="fixed", ckerlb=bounds[1], ckerub=bounds[2])
    spec <- list(d=1L, x=x, xcon=x, numeric.cols=1L, factor.cols=integer(), bws=bw)
    row <- suppressWarnings(.np_bias_main_rows(spec,at,0L))
    expect_identical(as.numeric(attr(row,"ridge.used")),0)
    expect_lt(max(abs(row-q$row)),1e-8)
  }
})
