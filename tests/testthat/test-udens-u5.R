test_that('U5 admission changes only fixed higher-order normal corrected density', {
 b<-structure(list(type='fixed',ckerbound='fixed',ckertype='epanechnikov',ckerorder=4L,ndim=1L,ncon=1L),class='bandwidth')
 for(order in c(4L,6L,8L)){z<-b;z$ckerorder<-order;expect_true(.np_udens_u5_request(z,NULL,'bias-corrected','asymptotic','pmzsd'))}
 cases<-list(type=c('generalized_nn','adaptive_nn'),ckerbound=c('none','range'),ckertype=c('uniform','gaussian','beta'),ckerorder=c(2L,3L,10L),ndim=c(0L,2L,3L),ncon=c(0L,2L))
 for(k in names(cases))for(v in cases[[k]]){z<-b;z[[k]]<-v;expect_false(.np_udens_u5_request(z,NULL,'bias-corrected','asymptotic','pmzsd'))}
 for(band in c('pointwise','bonferroni','simultaneous','all'))expect_false(.np_udens_u5_request(b,NULL,'bias-corrected','asymptotic',band))
 expect_false(.np_udens_u5_request(b,NULL,'estimate','asymptotic','pmzsd'))
 expect_false(.np_udens_u5_request(b,NULL,'bias-corrected','bootstrap','pmzsd'))
 expect_false(.np_udens_u5_request(b,np_udens_bias_control(),'bias-corrected','asymptotic','pmzsd'))
 expect_false(.np_udens_u5_request(b,NULL,'bias-corrected','asymptotic','pmzsd','symmetric-sd'))
 expect_false(.np_udens_u5_active(NULL))
})
test_that('U5 retains support reproduction and fails invalid geometry', {
 h<-.08
 for(t in c(.02,.5,.98)) {
  op<-.np_udens_u5_boundary('epanechnikov',max(-sqrt(5),-t/h),min(sqrt(5),(1-t)/h))
  for(k in 0:5){v<-integrate(function(u)op$L(u)*(h*u)^k,op$l,op$r,abs.tol=1e-12)$value;expect_equal(v,as.numeric(k==0),tolerance=1e-10)}
 }
 expect_error(.np_udens_u5_budget(9,50),'workspace')
 expect_error(.np_udens_u5_budget(1e6,50),'workspace bound')
 expect_error(.np_udens_u5_moments(c(0,(1:20)/21),h,.5),'interior')
 expect_error(.np_udens_u5_moments((1:20)/21,h,c(0,.5)),'interior')
 expect_error(.np_udens_u5_moments((1:20)/21,0,.5),'positive bandwidth')
 expect_error(.np_udens_u5_check(1,1.01,'negative control'),'negative control')
 expect_error(.np_udens_u5_check(Inf,1,'negative control'),'negative control')
})
