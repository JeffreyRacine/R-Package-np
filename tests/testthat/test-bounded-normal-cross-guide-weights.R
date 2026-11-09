test_that("cross-guide weights retain pre-underflow log probabilities", {
  ns<-asNamespace('np');rule<-get(".np_udens_q1_operator_rule",ns)
  plan<-list(edges=seq(-40,40,by=2));tq<-c(3,38);th<-c(4,2.5)
  for(order in c(128L,256L)){
    q<-rule(tq,plan,order,keep.log=TRUE);plain<-rule(tq,plan,order)
    expect_identical(q[names(plain)],plain)
    expect_equal(exp(q$log.prob),q$prob,tolerance=1e-15)
    term<-drop(q$T%*%(th-tq))+q$A-lbeta(th[1]+1,th[2]+1)
    expect_true(any(is.nan(q$prob*exp(term))))
    weights<-exp(q$log.prob+term)
    expect_true(all(is.finite(weights)))
    expect_equal(sum(weights),1,tolerance=1e-10)
  }
})
