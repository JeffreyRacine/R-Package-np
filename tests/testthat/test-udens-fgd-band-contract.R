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
