test_that("saved epanechnikov ANN degree 3 fit and influence rows agree with the physical reference", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  v <- readRDS(test_path("fixtures", "r112-epanechnikov-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=53, bandwidth.compute=FALSE,
    regtype="lp", degree=3, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="epanechnikov", ckerorder=2)
  ex <- data.frame(x=v$eval)
  fit <- npreg(bw, txdat=x, tydat=v$y, exdat=ex, gradients=TRUE)
  H <- npreghat(bw, txdat=x, exdat=ex, s=0L, output="matrix")
  expect_lt(max(abs(as.double(fit$mean) - c(0.015300608505679068, 1.4389911138076479, 1.0108235641460825, 1.0025595903465754))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(697, 322, 120, 663))]
  expect_lt(max(abs(selected - c(0.1308964830337649, 0.028018602316903757, 0.08279279859518802, 0.27869249124428075))), 1e-8)
  H <- npreghat(bw, txdat=x, exdat=ex, s=1L, output="matrix")
  expect_lt(max(abs(as.double(fit$grad) - c(3.386560262251083, 5.5970346733446315, 1.780645990339114, -8.323603571237216))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(697, 558, 663, 663))]
  expect_lt(max(abs(selected - c(-13.618573116817286, 1.0701358612562375, 4.015630223082284, 37.22794910257105))), 1e-8)
})
