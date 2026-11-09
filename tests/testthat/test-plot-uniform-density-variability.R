# Pure default-lane checks. Numerical grids and MPI proofs are external.
test_that("fixed uniform density admission is narrow and preserves ANN", {
  b <- list(ckertype="uniform", ckerorder=2L, type="fixed", ndim=1L,
            ncon=1L, ckerbound="none", ckerlb=-Inf, ckerub=Inf)
  expect_true(.np_plot_density_joint_supported(b,FALSE))
  expect_false(.np_plot_density_joint_supported(b,FALSE,TRUE))
  expect_false(.np_plot_density_uniform_fixed(b,TRUE))
  for (q in 1:2) {
    z <- b; z$ndim <- z$ncon <- q
    z$ckerbound <- "fixed"; z$ckerlb <- rep(0,q); z$ckerub <- rep(1,q)
    expect_true(.np_plot_density_uniform_fixed(z,FALSE))
    z$ckerub[1L] <- Inf
    expect_true(.np_plot_density_uniform_fixed(z,FALSE))
  }
  for (change in list(list(type="generalized_nn"),list(ckerbound="range"),
                     list(ndim=3L,ncon=3L),list(ndim=2L,ncon=1L),
                     list(ckerorder=4L)))
    expect_false(.np_plot_density_joint_supported(modifyList(b,change),FALSE))
  b$type <- "adaptive_nn"
  expect_true(.np_plot_density_joint_supported(b,FALSE))
  expect_true(.np_plot_density_joint_supported(b,FALSE,TRUE))
})

test_that("empirical correlation rescaling retains incumbent marginal variance", {
  # Independent overlapping/disjoint indicator windows, not a package kernel.
  K <- rbind(c(1,1,0,0),c(0,1,1,0),c(0,0,1,1),rep(0,4))
  A <- sweep(K,1,rowMeans(K))/4
  s <- c(.3,.4,.5,0)
  F <- .np_plot_density_rescale(A,s)$factor
  expected <- matrix(0,4,4)
  for (i in 1:3) for (j in 1:3)
    expected[i,j] <- s[i]*s[j]*sum(A[i,]*A[j,])/
      sqrt(sum(A[i,]^2)*sum(A[j,]^2))
  expect_equal(tcrossprod(F),expected,tolerance=1e-14)
  expect_equal(diag(tcrossprod(F)),s*s,tolerance=1e-14)
  expect_error(.np_plot_density_rescale(A,c(.3,.4,.5,.1)),"Unidentified correlation")
})


test_that("constant rounding residues do not invent uniform correlation", {
  testthat::local_mocked_bindings(.np_plot_density_unconditional =
    function(...) list(A=matrix(-6.938894e-18,1L,20L),estimate=1,
                      target="constant-window fixture"),
    .package = getNamespaceName(environment(.np_plot_density_joint)))
  b <- list(ckertype="uniform", ckerorder=2L, type="fixed", ndim=1L,
            ncon=1L, ckerbound="fixed", ckerlb=0, ckerub=1)
  x <- data.frame(x=seq(.025,.975,length.out=20L))
  expect_error(.np_plot_density_joint(b,x,data.frame(x=.5),.1),
               "unidentified", class="np_plot_joint_unavailable")
  expect_identical(.np_plot_density_joint(b,x,data.frame(x=.5),0)$critical,0)
})
