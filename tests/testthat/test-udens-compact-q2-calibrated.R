test_that('compact q2 calibrated admission preserves neighboring methods', {
  b<-structure(list(ckertype='uniform',type='fixed',ckerbound='fixed',ckerorder=2L,ndim=2L,ncon=2L),class='bandwidth')
  request<-function(b,method='asymptotic',band='pointwise',calibration=NULL,explicit=FALSE,control=NULL,center='bias-corrected')
    .np_udens_compact_q2_request(b,control,center,method,band,calibration,explicit)
  for(kernel in c('uniform','epanechnikov')) {
    z<-modifyList(b,list(ckertype=kernel))
    expect_true(request(z));expect_true(request(z,band='bonferroni'))
    # Normal keeps its separate owner; Gaussian-factor S/all remains unqualified.
    for(band in c('pmzsd','simultaneous','all'))expect_false(request(z,band=band))
    for(cal in c('symmetric-sd','asymmetric-sd'))for(band in c('pointwise','simultaneous','bonferroni','all')) {
      expect_false(request(z,'bootstrap',band,cal,TRUE))
      expect_false(request(z,'bootstrap',band,cal,FALSE))
    }
    expect_true(.np_udens_compact_normal_q2_request(z,NULL,'bias-corrected','asymptotic','pmzsd',NULL))
  }
  for(change in list(list(ckertype='gaussian'),list(ckertype='beta'),list(ndim=1L,ncon=1L),list(ndim=3L,ncon=3L),list(ncon=1L),list(type='adaptive_nn'),list(type='generalized_nn'),list(ckerbound='none'),list(ckerorder=4L)))
    expect_false(request(modifyList(b,change)))
  expect_false(request(b,center='estimate'));expect_false(request(b,control=np_udens_bias_control()))
  for(cal in list(NULL,'legacy','pmzsd','smoothed-percentile'))expect_false(request(b,'bootstrap','all',cal,TRUE))
  expect_false(.np_udens_compact_q2_active(NULL));expect_false(.np_udens_compact_q2_active(list(construction='density-fixed-compact-cubic-v1')))
})

test_that('compact q2 common-measure merge retains between-tile means', {
  a<-list(mean=c(1,2,3,4,5,6),variance=rep(0,6),cross=matrix(0,6,5),logmass=log(.2))
  b<-list(mean=c(4,3,2,1,0,-1),variance=rep(0,6),cross=matrix(0,6,5),logmass=log(.8))
  z<-.np_udens_compact_q2_raw_merge(a,b,1L);mu<-.2*a$mean+.8*b$mean
  expect_equal(z$mean,mu,tolerance=1e-14)
  expect_equal(z$variance,.2*(a$mean-mu)^2+.8*(b$mean-mu)^2,tolerance=1e-14)
  expected<-.2*outer(a$mean-mu,(a$mean-mu)[2:6])+.8*outer(b$mean-mu,(b$mean-mu)[2:6])
  expect_equal(z$cross,expected,tolerance=1e-14)
  expect_gt(max(abs(z$cross)),1);expect_equal(exp(z$logmass),1,tolerance=1e-14)
})

test_that('compact q2 admission rejects invalid geometry and excessive new work', {
  x<-cbind((1:20)/21,(20:1)/21);g<-matrix(c(.2,.8,.4,.6),2)
  expect_invisible(.np_udens_compact_q2_validate(x,c(.1,.2),g,'uniform'))
  for(boundary in c(0,1)){bad<-g;bad[1,1]<-boundary;expect_error(.np_udens_compact_q2_validate(x,c(.1,.2),bad,'uniform'),'request')}
  expect_error(.np_udens_compact_q2_validate(x,c(0,.2),g,'uniform'),'request')
  expect_error(.np_udens_compact_q2_work_forecast(20,9,c(8193,1),c(3,3)),'work bound')
  expect_error(.np_udens_compact_q2_work_forecast(20,9,c(8192,8192),c(3,3)),'work bound')
  expect_error(.np_udens_compact_q2_work_forecast(20,9,c(100,100),c(3,3),held=192*1024^2),'workspace')
})

# Retained-sample roots, independent signed tails, full grids and real MPI proof
# live in the bounded campaign evidence, not repeated in ordinary package checks.
