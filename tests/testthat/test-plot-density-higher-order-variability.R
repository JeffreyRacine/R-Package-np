# Pure admission checks only; independent kernel and public/MPI proofs are external.
test_that("signed fixed density kernels have a density-only joint admission", {
  b <- list(ckertype="gaussian", ckerorder=4L, type="fixed", ndim=1L,
            ncon=1L, ckerbound="none", ckerlb=-Inf, ckerub=Inf)
  for (kernel in c("gaussian","epanechnikov")) for (order in c(4L,6L,8L)) {
    z <- b; z$ckertype <- kernel; z$ckerorder <- order
    expect_true(.np_plot_density_joint_supported(z,FALSE))
    expect_false(.np_plot_density_joint_supported(z,FALSE,TRUE))
  }
  expect_false(.np_plot_density_higher_order_fixed(b,TRUE))
  for (q in 1:2) {
    z <- b; z$ndim <- z$ncon <- q
    z$ckerbound <- "fixed"; z$ckerlb <- rep(0,q); z$ckerub <- rep(1,q)
    expect_true(.np_plot_density_joint_supported(z,FALSE))
    z$ckerub[1L] <- Inf
    expect_false(.np_plot_density_joint_supported(z,FALSE))
  }
  for (change in list(list(type="generalized_nn"),list(type="adaptive_nn"),
                     list(ckerbound="range"),list(ndim=3L,ncon=3L),
                     list(ndim=2L,ncon=1L),list(ckerorder=3L),
                     list(ckertype="beta"),list(ckertype="uniform")))
    expect_false(.np_plot_density_joint_supported(modifyList(b,change),FALSE))
})


test_that("higher-order beta joint admission requires finite known fixed support", {
  for (q in 1:2) for (order in c(4L,6L,8L)) {
    b <- list(ckertype="beta", ckerorder=order, type="fixed", ndim=q,
              ncon=q, ckerbound="fixed", ckerlb=rep(-2,q), ckerub=rep(3,q))
    expect_true(.np_plot_density_joint_supported(b,FALSE))
    expect_false(.np_plot_density_higher_order_fixed(b,TRUE))
    expect_false(.np_plot_density_joint_supported(b,FALSE,TRUE))
    for (change in list(list(ckerbound="none"),list(ckerbound="range"),
                       list(ckerub=rep(Inf,q)),list(ckerub=rep(-2,q)),
                       list(type="generalized_nn"),list(type="adaptive_nn")))
      expect_false(.np_plot_density_joint_supported(modifyList(b,change),FALSE))
  }
})
