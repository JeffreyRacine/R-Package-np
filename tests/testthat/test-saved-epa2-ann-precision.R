test_that("saved epanechnikov ANN degree 2 fit and influence rows agree with the physical reference", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  v <- readRDS(test_path("fixtures", "r112-epanechnikov-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=43, bandwidth.compute=FALSE,
    regtype="lp", degree=2, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="epanechnikov", ckerorder=2)
  ex <- data.frame(x=v$eval)
  fit <- npreg(bw, txdat=x, tydat=v$y, exdat=ex, gradients=TRUE)
  H <- npreghat(bw, txdat=x, exdat=ex, s=0L, output="matrix")
  expect_lt(max(abs(as.double(fit$mean) - c(0.015524963358829556, 1.442687400583536, 1.0011980661488433, 1.0210729834711072))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(697, 625, 101, 663))]
  expect_lt(max(abs(selected - c(0.09418090270370622, 0.030077682578600266, 0.061260747687858594, 0.21837981645090018))), 1e-8)
  H <- npreghat(bw, txdat=x, exdat=ex, s=1L, output="matrix")
  expect_lt(max(abs(as.double(fit$grad) - c(2.767750792984437, 4.171057132659753, 1.695496568438289, -1.3702585269610335))), 1e-8)
  selected <- H[cbind(seq_len(nrow(ex)), c(697, 275, 663, 663))]
  expect_lt(max(abs(selected - c(-6.190706973416302, -0.5375440089967234, 4.972864837911807, 17.027756623370184))), 1e-8)
})
