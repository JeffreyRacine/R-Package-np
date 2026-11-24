test_that("saved beta ANN degree 3 fit and influence rows agree with the physical reference", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  v <- readRDS(test_path("fixtures", "r112-cubic-beta-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=24, bandwidth.compute=FALSE,
    regtype="lp", degree=3, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="beta", ckerorder=2, ckerbound="fixed", ckerlb=-1, ckerub=1)
  ex <- data.frame(x=v$eval)
  fit <- npreg(bw, txdat=x, tydat=v$y, exdat=ex, gradients=TRUE)
  H <- npreghat(bw, txdat=x, exdat=ex, s=0L, output="matrix")
  expect_lt(max(abs(as.double(fit$mean) - c(-1.0504363387140758, 0.014800797480466872, 0.9991391777867784, 1.173516519593604))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(376, 218, 158, 190))]
  expect_lt(max(abs(selected - c(0.6492081715902706, 0.10261923851955485, 0.4786988559736247, 0.9826464447489853))), 1e-8)
  H <- npreghat(bw, txdat=x, exdat=ex, s=1L, output="matrix")
  expect_lt(max(abs(as.double(fit$grad) - c(69.53132571879681, -0.6275535883799216, 0.5133448696758263, 0.3100959079815286))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(376, 189, 390, 190))]
  expect_lt(max(abs(selected - c(-324.2844520796394, 4.229724219920694, -37.94126912085631, 1.2076162972133282))), 1e-8)
})
