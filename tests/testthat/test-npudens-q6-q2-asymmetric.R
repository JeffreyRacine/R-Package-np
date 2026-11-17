test_that("q2 directional calibration preserves signed errors and reflection", {
  root <- rbind(c(-7,-1,1,2,3,9),c(8,9,10,11,13,17),c(-10,-7,-6,-5,-4,-3))
  z <- .np_udens_q6_q2_directional(root,.1)
  mirror <- .np_udens_q6_q2_directional(-root,.1)
  expect_equal(z$plus,mirror$minus,tolerance=1e-14)
  expect_equal(z$minus,mirror$plus,tolerance=1e-14)
  expect_equal(z$pointwise,mirror$pointwise,tolerance=1e-14)
  expect_equal(z$maximum,mirror$maximum,tolerance=1e-14)
  expect_equal(z$plus^2+z$minus^2,2*apply(root,1,var)*5/6,tolerance=1e-14)
  # All-positive and all-negative roots remain valid: only variability is centered.
  expect_gt(min(z$plus),0);expect_gt(min(z$minus),0)
  centered <- .np_udens_q6_q2_directional(root-rowMeans(root),.1)
  expect_equal(z$plus,centered$plus,tolerance=1e-14)
  expect_false(isTRUE(all.equal(z$pointwise,centered$pointwise)))
  symmetric <- rbind(c(-3,-2,-1,1,2,3),c(-6,-4,-2,2,4,6))
  a <- .np_udens_q6_q2_directional(symmetric,.1)
  se <- sqrt(rowMeans(symmetric^2))
  expect_equal(a$plus,se,tolerance=1e-14);expect_equal(a$minus,se,tolerance=1e-14)
  expect_equal(a$maximum,apply(abs(symmetric)/se,2,max),tolerance=1e-14)
  expect_error(.np_udens_q6_q2_directional(matrix(1,2,5),.05),"variability")
  bad<-root;bad[1,1]<-Inf;expect_error(.np_udens_q6_q2_directional(bad,.05),"invalid")
  expect_error(.np_udens_q6_q2_directional(root[,1,drop=FALSE],.05),"invalid")
  expect_error(.np_udens_q6_q2_directional(root,0),"invalid")
})

test_that("q2 asymmetric selection is narrow and both band restrictions are explicit", {
  b<-structure(list(bw=c(.4,.5),scaling=FALSE,type="fixed",ckerbound="none",
    ckertype="uniform",ckerorder=2L,ndim=2L,ncon=2L),class="bandwidth")
  for(band in c("pointwise","simultaneous"))
    expect_silent(.np_udens_q6_public_spec(b,"asymmetric-sd",band))
  for(band in c("bonferroni","all","pmzsd"))
    expect_error(.np_udens_q6_public_spec(b,"asymmetric-sd",band),"pointwise or simultaneous")
  for(band in c("pointwise","simultaneous","bonferroni","all"))
    expect_silent(.np_udens_q6_public_spec(b,"symmetric-sd",band))
  x<-data.frame(x=seq(-1,1,length.out=12),y=cos(seq_len(12)))
  a<-.np_udens_q6_q2_asymmetric_prepare(b,x,42L);s<-.np_udens_q6_q2_prepare(b,x,42L)
  expect_identical(a$h,s$h);expect_identical(a$x,s$x)
  expect_true(.np_udens_q6_q2_active(a));expect_true(.np_udens_q6_q2_active(s))
  expect_true(.np_udens_q6_q2_asymmetric_active(a));expect_false(.np_udens_q6_q2_asymmetric_active(s))
  expect_false(.np_udens_q6_q2_active(NULL))
  expect_error(.np_udens_q6_q2_asymmetric_public(a,x,.05,19L,"all"),"invalid")
  expect_error(.np_udens_q6_q2_asymmetric_public(a,x,.05,1L,"pointwise"),"two conditional")
})

test_that("signed q2 worker roots preserve exact original symmetric arithmetic and RNG", {
  x<-cbind(seq(-1,1,length.out=12),cos(seq_len(12)))
  for(kernel in c("uniform","epanechnikov")) {
    grid<-rbind(c(-.4,.3),c(.2,-.3),c(.6,.7))
    r<-.np_udens_q6_q2_reference(x,c(.4,.5),grid,kernel,17L)
    tasks<-.np_with_seed(42L,.np_udens_q6_q2_plan(.Random.seed,nrow(x),17L))
    before<-if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))get(".Random.seed",.GlobalEnv)else NULL
    old<-.np_udens_q6_q2_chunk(tasks[[1]],r)
    rr<-r;rr$directional<-TRUE
    signed<-.np_udens_q6_q2_chunk(tasks[[1]],rr)
    Y<-.np_udens_q6_q2_samples(tasks[[1]],x,c(.4,.5))
    expected<-.np_udens_q6_q2_curve(Y,r$grid,r$h,kernel,32L)-r$reference.center
    expect_identical(unname(signed[,1:3]),unname(t(expected)))
    expect_identical(unname(old[,1:3]),unname(t(abs(expected/r$se))))
    expect_identical(colnames(signed),colnames(old))
    expect_identical(old[,4],signed[,4])
    after<-if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))get(".Random.seed",.GlobalEnv)else NULL
    expect_identical(before,after)
  }
})
