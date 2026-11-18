test_that("q2 Q6 admission is separate from q1 and asymmetric ownership", {
  b <- structure(list(type="fixed",ckerbound="none",ckertype="uniform",
    ckerorder=2L,ndim=2L,ncon=2L),class="bandwidth")
  expect_silent(.np_udens_q6_public_spec(b,"symmetric-sd"))
  expect_error(.np_udens_q6_public_spec(b,"asymmetric-sd"),"pointwise or simultaneous")
  for (field in c("type","ckerbound","ckertype","ckerorder","ndim","ncon")) {
    z <- b; z[[field]] <- switch(field,type="adaptive_nn",ckerbound="fixed",
      ckertype="gaussian",ckerorder=4L,ndim=3L,ncon=1L)
    expect_error(.np_udens_q6_q2_spec(z),"q=2 symmetric")
  }
})

test_that("q2 Q6 moment law retains common empirical donors", {
  x <- cbind(c(-2,-1,0,1,2,3),c(1,0,2,4,3,6)); h <- c(.4,.6)
  z <- .np_udens_q6_q2_moments(x,h); e <- sweep(x,2,colMeans(x))
  expect_equal(z$S,crossprod(e)/nrow(x)+diag(h*h),tolerance=1e-14)
  expect_equal(z$V[1L,2L],mean(e[,1L]*e[,2L]),tolerance=1e-14)
  expect_equal(z$V[3L,3L],
    mean(e[,1L]^4)+6*h[1L]^2*mean(e[,1L]^2)+6*h[1L]^4-z$S[1L,1L]^2,
    tolerance=1e-14)
  expect_gt(abs(z$V[1L,2L]),.1)
})

test_that("q2 Q6 primitives reject degenerate geometry and bounded storage", {
  x <- cbind(seq(-1,1,length.out=12),cos(seq_len(12)))
  expect_silent(.np_udens_q6_q2_geometry(x,c(.3,.4),x[1:2,,drop=FALSE],"uniform"))
  expect_error(.np_udens_q6_q2_geometry(cbind(x[,1],x[,1]),c(.3,.4),x,"uniform"),
    "positive definite")
  expect_error(.np_udens_q6_q2_budget(1e6,2500,c(50,50),999),"192 MiB")
  z <- .np_udens_q6_q2_budget(1000,2500,c(50,50),999)
  expect_gt(z$block,0L); expect_lte(z$bytes,192*1024^2)
})

test_that("q2 Q6 preserves every duplicate and signed-zero target row", {
  x <- cbind(seq(-1,1,length.out=12),cos(seq_len(12)))
  grid <- rbind(c(-0,.2),c(0,.2),c(.3,-.4),c(-0,.2))
  z <- .np_udens_q6_q2_reference(x,c(.4,.5),grid,"uniform",19L)
  expect_identical(z$index,c(1L,1L,2L,1L))
  expect_true(all(is.finite(z$center[z$index])))
  expect_identical(z$center[z$index][1L],z$center[z$index][2L])
})
