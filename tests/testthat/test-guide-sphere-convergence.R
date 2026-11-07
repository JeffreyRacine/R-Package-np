test_that("sphere roundoff acceptance recognizes the existing convergence threshold", {
  f<-get('.np_udens_beta_sphere_accept',asNamespace('npRmpi'))
  p<-readRDS(test_path('fixtures','sphere-convergence-retained.rds'))
  expect_false(do.call(f,p[1:5]))
  expect_true(do.call(f,p[1:6]))
  expect_lt(p$independent_kkt,p$gtol)
  for(tol in c(1e-10,1e-7,1e-4)) {
    for(relative in c(.99,1,1.01)) {
      expect_identical(f(0,.Machine$double.eps,0,relative*tol,1.1*tol,tol),relative<=1)
    }
    expect_false(f(0,.Machine$double.eps,0,tol,tol,tol))
    expect_false(f(0,100*.Machine$double.eps,0,.9*tol,1.1*tol,tol))
    expect_false(f(0,.Machine$double.eps,-100*.Machine$double.eps,.9*tol,1.1*tol,tol))
    expect_false(f(0,.Machine$double.eps,0,1.01*tol,1.1*tol,tol))
  }
  old<-function(value,trial,predicted,newnorm,oldnorm) {
    if(trial<=value+1e-4*predicted)return(TRUE)
    roundoff<-32*.Machine$double.eps*(1+abs(value))
    abs(predicted)<=roundoff && abs(trial-value)<=roundoff && newnorm<=.5*oldnorm
  }
  g<-expand.grid(value=c(-10,0,10),delta=c(-1e-5,0,1e-15,1e-4),predicted=c(-1,-1e-15,0),newnorm=c(0,1e-8,.1),oldnorm=c(1e-8,.01,1))
  for(i in seq_len(nrow(g))) {
    x<-as.list(g[i,]);a<-list(x$value,x$value+x$delta,x$predicted,x$newnorm,x$oldnorm)
    expect_identical(do.call(f,a),do.call(old,a))
  }
})
