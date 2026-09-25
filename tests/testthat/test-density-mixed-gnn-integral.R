test_that("mixed GNN CVLS uses whole-support integration and deleted radii", {
  old <- options(np.messages=FALSE, np.tree=FALSE, np.extendednn=FALSE)
  on.exit(options(old), add=TRUE)
  dat <- data.frame(x=c(-1.31,-.89,-.41,-.17,.24,.52,1.06,1.57),
                    u=factor(rep(c("a","b"),4),levels=c("a","b","unused")))
  # Independent whole-line Gaussian quadrature, explicit AA product sum over
  # all declared categories, and literal deleted order statistics (k=2,6).
  expected <- c(-0.014690165539861397, 0.051809375654240458)
  for (j in 1:2) {
    b <- npudensbw(dat=dat,bws=c(c(2,6)[j],.3),bwtype="generalized_nn",
                    bwmethod="cv.ls",bandwidth.compute=FALSE)
    got <- npudensbw(dat=dat,bws=b,eval.only=TRUE,nmulti=1L,bwsolver="powell")
    expect_equal(unname(got$fval),expected[j],tolerance=2e-8)
  }
})
