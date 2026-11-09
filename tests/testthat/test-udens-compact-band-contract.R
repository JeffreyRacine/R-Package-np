test_that('compact cubic operators reproduce independent truncated moments', {
  for(kernel in c('uniform','epanechnikov'))for(h in c(.02,.3,.99999,1,1.00001,4)) {
    target<-c(.001,.2,.5,.8,.999);p<-.np_udens_compact_operator_plan(target,h,kernel)
    for(j in seq_along(target))for(k in 0:3) {
      one<-.np_udens_compact_operator_plan(target[j],h,kernel)
      z<-integrate(function(u)drop(.np_udens_compact_operator(one,u))*u^k,
        max(0,target[j]-p$span*h),min(1,target[j]+p$span*h),
        abs.tol=1e-11,rel.tol=1e-10)$value
      expect_lte(abs(z-target[j]^k),1e-12+1e-10*abs(target[j]^k))
    }
    u<-c(.0001,.11,.38,.61,.89,.9999)
    reflect<-.np_udens_compact_operator_plan(1-target,h,kernel)
    expect_equal(.np_udens_compact_operator(p,u),.np_udens_compact_operator(reflect,1-u),tolerance=1e-9)
  }
  expect_error(.np_udens_compact_operator_plan(.5,0,'uniform'),'geometry')
  expect_error(.np_udens_compact_operator_plan(.5,.1,'gaussian'),'geometry')
})

test_that('compact admission leaves incumbent owners and explicit controls unchanged', {
  b<-structure(list(ckertype='uniform',ckerorder=2L,ckerbound='fixed',ndim=1L,ncon=1L,type='fixed'),class='bandwidth')
  request<-function(bws=b,control=NULL,center='bias-corrected',method='asymptotic',calibration=NULL)
    .np_udens_compact_request(bws,control,center,method,calibration)
  expect_true(request());b$ckertype<-'epanechnikov';expect_true(request(b))
  for(kernel in c('beta','gaussian')) {z<-b;z$ckertype<-kernel;expect_false(request(z))}
  for(method in c('none','bootstrap'))expect_false(request(method=method))
  expect_false(request(center='estimate'));expect_false(request(calibration='legacy'))
  for(construction in c('legacy','positive-covariance','joint-scale'))expect_false(request(control=np_udens_bias_control(construction)))
  z<-b;z$ndim<-z$ncon<-2L;expect_false(request(z))
  z<-b;z$ckerorder<-4L;expect_false(request(z))
  z<-b;z$type<-'generalized_nn';expect_false(request(z))
  z<-b;z$ckerbound<-'none';expect_false(request(z))
  z<-b;z$ckerbound<-'range';expect_true(request(z))
  for(cl in c('rbandwidth','dbandwidth','conbandwidth')){z<-b;class(z)<-cl;expect_false(request(z))}
})

test_that('compact providers retain explicit bands and interior support', {
  x<-seq(.05,.95,length.out=11L)
  for(band in c('pmzsd','missing'))expect_error(.np_udens_compact_hidden(x,.1,.5,'uniform',.05,band),'request')
  expect_error(.np_udens_compact_hidden(x[-(1:2)],.1,.5,'uniform',.05,'pointwise'),'request')
  expect_error(.np_udens_compact_hidden(x,.1,0,'uniform',.05,'pointwise'),'request')
  expect_error(.np_udens_compact_hidden(x,.1,1,'uniform',.05,'pointwise'),'request')
  expect_error(.np_udens_compact_hidden(x,.1,.5,'uniform',0,'pointwise'),'request')
  for(n in 10:11) {
    p<-.np_udens_compact_prepare(x[seq_len(n)],.2,c(.2,.5,.8),'uniform',TRUE)
    expect_true(all(p$pilot_variance>0));expect_true(all(p$empirical_variance>0))
    expect_equal(diag(p$empirical),p$empirical_variance,tolerance=1e-12)
    expect_lt(max(abs(p$guide_residual_mean)),1e-9)
  }
})
