test_that("cubic Gaussian generalized-NN correction follows the fitted mean and slope", {
  if (!.mpi_pool_active() && !.mpi_suite_local_mode_owned()) {
    if (!spawn_mpi_slaves(n=1L)) skip("MPI workers unavailable")
    on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  }
  cases <- readRDS(testthat::test_path("fixtures","gaussian-cubic-gnn-rows.rds"))
  for (saved in cases) {
    x <- data.frame(x=saved$x)
    y <- sin(x$x)
    g <- data.frame(x=saved$at)
    b <- npregbw(xdat=x, ydat=y, bws=saved$k, bandwidth.compute=FALSE,
                 bwtype="generalized_nn", ckertype="gaussian", regtype="lp",
                 degree=3L, basis="glp", bernstein.basis=TRUE)
    state <- .np_bias_regression_prepare(b,x,y,gradient.order=1L,
                                        law="mammen",calibration="absolute")
    rows <- .np_bias_regression_rows(state,g,data.frame(row=1L,target=if(saved$nu) "d1" else "mean"))
    z <- x$x-saved$at
    h <- sort(abs(z),partial=saved$k)[saved$k]
    u <- z/h
    sw <- exp(-u^2/4)
    design <- cbind(1,u,u^2,u^3)
    physical <- qr.solve(design*sw,diag(sw))[saved$nu+1L,]/h^saved$nu
    expect_equal(rows$main.regularization$ridge,0)
    expect_lt(max(abs(rows$original-physical)),1e-8)
    expect_lt(abs(rows$original_center-sum(physical*y)),1e-8)
  }
})
