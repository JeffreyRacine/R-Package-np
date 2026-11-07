test_that("joint region certifies the returned normalizer and actual stationary choice", {
  ns<-asNamespace("np");cert<-get(".np_udens_beta_region_certificate",environment(get(".np_udens_beta_region_q2",ns)))
  p<-readRDS(test_path("fixtures","region-certificate-retained.rds"))
  a<-cert(p$runs,p$index,p$target,p$transform,p$sign)
  expect_gt(a$optimization_refinement,1e-10)
  expect_lt(a$refinement,1e-10);expect_lt(a$kkt,1e-4);expect_lt(a$ranking_gap,1e-10)
  expect_equal(a$value,p$oracle$value,tolerance=1e-10)
  bad<-cert;e<-new.env(parent=ns);e$.np_udens_beta_joint_normalizer<-function(...)get(".np_udens_beta_joint_normalizer",ns)(...)+1e-6;environment(bad)<-e
  expect_error(bad(p$runs,p$index,p$target,p$transform,p$sign),"independent-normalizer gate")
  expect_error(cert(p$runs,p$index,p$target+c(.25,-.11,.07,0,0),p$transform,p$sign),"stationarity")
  # A stationary maximum cannot pass the minimizing-start certificate.
  theta<-rep(0,5);target<-get(".np_udens_beta_joint_rule",ns)(theta,32L)$mean
  eta2<-c(-2,0,0,0,0);theta2<-c(exp(eta2[1:4])-1,eta2[5])
  runs<-list(list(eta=rep(0,5),theta=theta,direction=c(1,0,0,0,0),value=0),
             list(eta=eta2,theta=theta2,direction=c(-1,0,0,0,0),value=0))
  expect_error(cert(runs,1L,target,diag(5),1),"refined-ranking gate")
})
