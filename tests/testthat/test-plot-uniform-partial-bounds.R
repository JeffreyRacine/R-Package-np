test_that("uniform joint admission accepts only valid known bound geometries", {
  for (q in 1:2) for (lower in c(-Inf,0)) for (upper in c(1,Inf)) {
    b <- list(ckertype="uniform",ckerorder=2L,type="fixed",ndim=q,ncon=q,
              ckerbound="fixed",ckerlb=rep(lower,q),ckerub=rep(upper,q))
    expect_true(.np_plot_density_uniform_fixed(b,FALSE))
    expect_true(.np_plot_density_joint_supported(b,FALSE))
    expect_false(.np_plot_density_uniform_fixed(b,TRUE))
    expect_false(.np_plot_density_uniform_fixed(b,FALSE,TRUE))
    for (invalid in list(c(NA,1),c(NaN,1),c(0,NA),c(0,NaN),c(1,0),c(1,1),c(Inf,Inf),c(-Inf,-Inf),c(Inf,1),c(0,-Inf))) {
      z <- b; z$ckerlb[1L] <- invalid[1L]; z$ckerub[1L] <- invalid[2L]
      expect_false(.np_plot_density_uniform_fixed(z,FALSE))
    }
    for (change in list(list(type="generalized_nn"),list(ckerbound="range"),
                        list(ndim=3L,ncon=3L),list(ncon=q-1L),list(ckerorder=4L)))
      expect_false(.np_plot_density_uniform_fixed(modifyList(b,change),FALSE))
  }
})
