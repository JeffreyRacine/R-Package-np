test_that('fixed-bound Gaussian plot operator reproduces cubic polynomials', {
  gl<-.np_udens_beta_rule(c(1,1),128L)
  x<-c(.001,.2,.5,.8,.999)
  for(h in c(.05,.4,1,20)) {
    a<-1+.np_udens_fgd_operator(x,gl[,1],h)
    for(k in 0:3)expect_equal(drop(a%*%(gl[,2]*gl[,1]^k)),x^k,tolerance=1e-9)
    expect_equal(a,1+.np_udens_fgd_operator(1-x,1-gl[,1],h),tolerance=1e-9)
  }
  expect_error(.np_udens_fgd_operator(x,gl[,1],0),'geometry')
  expect_error(.np_udens_fgd_operator(x,NA_real_,.1),'geometry')
})

test_that('FG-D envelopes complete signed intervals rather than maximum scales', {
  # Both signed limits positive: a larger scale raises the lower endpoint.
  pilot<-list(lower=c(3,4),upper=c(5,6))
  empirical<-list(lower=c(2,3),upper=c(3,4))
  z<-.np_udens_fgd_envelope(pilot,empirical)
  expect_identical(z$lower,c(2,3));expect_identical(z$upper,c(5,6))
  expect_false(identical(z$lower,pilot$lower))
  expect_error(.np_udens_fgd_envelope(list(lower=2,upper=1),list(lower=0,upper=1)),'constituent')
})

test_that('FG-D rejects excessive owned workspace before allocation', {
  expect_error(.np_udens_fgd_disk(1000,1e6,2,9999,TRUE,TRUE),'budget')
  expect_error(.np_udens_fgd_hidden(matrix(.5,10,1),.1,matrix(0,1,1),.05,'asymptotic','pointwise'),'request')
})

test_that('FG-D selection preserves explicit and legacy routes', {
  bws<-structure(list(ckertype='gaussian',ckerbound='fixed'),class='bandwidth')
  request<-function(b=bws,control=NULL,center='bias-corrected',method='asymptotic',calibration=NULL)
    .np_udens_fgd_request(b,control,center,method,calibration)
  expect_true(request());expect_true(request(method='bootstrap',calibration='symmetric-sd'))
  expect_true(request(method='bootstrap',calibration='asymmetric-sd'))
  expect_false(request(method='bootstrap',calibration='legacy'))
  expect_false(request(control=np_udens_bias_control('legacy')))
  expect_false(request(control=np_udens_bias_control()))
  expect_false(request(center='estimate'));expect_false(request(method='none'))
  for(kernel in c('beta','epanechnikov','uniform')){b<-bws;b$ckertype<-kernel;expect_false(request(b))}
  b<-bws;b$ckerbound<-'none';expect_false(request(b))
  b$ckerbound<-'range';expect_true(request(b)) # admitted to the early fixed-bound refusal
  for(class in c('rbandwidth','dbandwidth','conbandwidth')){b<-bws;class(b)<-class;expect_false(request(b))}
})
