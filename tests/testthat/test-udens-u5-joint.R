test_that('U5 joint factors retain dependent directions and covariance', {
 D<-rbind(sin(1:30),cos(1:30),sin(1:30)+cos(1:30))
 F<-.np_udens_u5_gram_factor(D)
 expect_identical(dim(F),c(3L,3L));expect_equal(tcrossprod(F),tcrossprod(D),tolerance=1e-12)
 expect_error(.np_udens_u5_gram_factor(matrix(NA_real_,3,30)),'invalid')
 expect_error(.np_udens_u5_gram_factor(matrix(1,30,3)),'invalid')
 for(n in c(100L,400L,1000L))expect_lt(.np_udens_u5_joint_budget(n,50L,9999L),192*1024^2)
 expect_error(.np_udens_u5_joint_budget(1000L,10000L,9999L),'workspace bound')
})
test_that('higher-order U5 reproduces all defining polynomial moments', {
 for(order in c(4L,6L,8L))for(t in c(.02,.5,.98)) {
  h<-.08;op<-.np_udens_u5_boundary('epanechnikov',max(-sqrt(5),-t/h),min(sqrt(5),(1-t)/h),kernel.order=order)
  for(k in 0:(order+1L))expect_equal(integrate(function(u)op$L(u)*(h*u)^k,op$l,op$r,abs.tol=1e-12)$value,as.numeric(k==0L),tolerance=1e-10)
 }
 expect_error(.np_udens_u5_boundary('epanechnikov',-1,1,kernel.order=2),'order')
})
