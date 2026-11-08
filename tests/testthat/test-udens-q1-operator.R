test_that('q1 operator rules preserve the frozen guide and fixed budgets', {
  theta<-c(.043018331478746313,.23782979033327017)
  g<-c(.00267337365792877,.20589632514546927,.998465835946877)
  p<-.np_udens_q1_operator_plan(theta,.003439661669438287,g)
  expect_identical(p,.np_udens_q1_operator_plan(theta,.003439661669438287,c(rev(g),g)))
  for(outer in c(128L,256L)){
    z<-.np_udens_q1_operator_rule(theta,p,outer)
    ex<-.np_udens_bounded_normal_analytic(theta)
    expect_true(all(z$prob>0));expect_lte(z$nodes,8192)
    expect_equal(z$A,ex$A,tolerance=1e-12)
    expect_equal(unname(z$mean),ex$mean,tolerance=1e-12)
    expect_equal(unname(z$I),ex$I,tolerance=1e-11)
  }
  expect_error(.np_udens_q1_operator_plan(theta,1e-8,g),'node budget')
  expect_error(.np_udens_q1_operator_rule(theta,p,64L),'refinement')
  expect_error(.np_udens_q1_operator_plan(theta,.01,c(0,.5)),'geometry')
  # No tolerance can be bypassed by returning a previously accepted result.
  expect_error(.np_udens_beta_check(1,1.001,'operator mean',1e-10),'refinement failed')
})

test_that('q1 half-bandwidth first moments resolve the signed corrected score', {
  # Target/guide from the saved failure; deterministic compact fixture data.
  # Independent adaptive logit/angular integrals give mu=1.148821725796.
  x<-seq(.03,.97,length.out=1000)
  theta<-c(.043018331478746313,.23782979033327017)
  h<-.003439661669438287;g<-.20589632514546927
  a<-.np_udens_beta_q1_moments(x,h,g,theta,64L,128L,TRUE)
  b<-.np_udens_beta_q1_moments(x,h,g,theta,128L,256L,TRUE)
  expect_equal(unname(a$mu),1.148821725796,tolerance=1e-11)
  for(k in c('mu','B','center','pilot_variance','guide_variance','empirical_variance'))
    expect_silent(.np_udens_beta_check(a[[k]],b[[k]],k,if(k=='mu')1e-10 else 1e-9))
  expect_equal(diag(a$covariance),a$pilot_variance,tolerance=1e-14)
  th<-rbind(theta,theta+c(.02,-.02))
  c<-.np_udens_beta_q1_common(x,h,g,th,64L,128L)
  d<-.np_udens_beta_q1_common(x,h,g,th,128L,256L)
  expect_silent(.np_udens_beta_check(c$mu,d$mu,'refitted moments',1e-10))
  expect_silent(.np_udens_beta_check(c$mu[1,],a$mu,'original common row',1e-10))
})
