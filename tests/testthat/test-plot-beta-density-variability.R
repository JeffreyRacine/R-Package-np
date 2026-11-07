test_that("fixed beta density admission cannot leak into CDF or unsupported classes", {
  b <- list(ckertype="beta", type="fixed", ckerbound="fixed", ndim=1L, ncon=1L, ckerorder=2L)
  expect_true(.np_plot_density_joint_supported(b,FALSE))
  expect_false(.np_plot_density_joint_supported(b,FALSE,TRUE))
  for(change in list(list(type="adaptive_nn"),list(type="generalized_nn"),list(ckerbound="range"),list(ckerbound="none"),list(ndim=3L,ncon=3L),list(ndim=2L,ncon=1L),list(ckerorder=4L))) {
    z<-modifyList(b,change);expect_false(.np_plot_density_joint_supported(z,FALSE))
  }
  b$ndim<-b$ncon<-2L;expect_true(.np_plot_density_joint_supported(b,FALSE))
  expect_error(.np_plot_density_rescale(matrix(0,2,4),c(1,1)),"Unidentified correlation")
  expect_equal(.np_plot_density_rescale(matrix(0,2,4),c(0,0))$factor,matrix(0,2,4))
})
