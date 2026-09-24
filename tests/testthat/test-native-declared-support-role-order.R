test_that("declared support retains distinct roles and factor encodings", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L),"MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(xu=factor(rep(c("z","a","b","a"),2),levels=c("z","a","unused","b")),
    xo=ordered(rep(c(-2.5,.5,3.5,.5),2),levels=c(-2.5,-.5,.5,3.5)))
  y<-data.frame(yu=factor(rep(c("yes","no"),4),levels=c("unused0","yes","unused2","no","unused4")),
    yo=ordered(rep(c(1.5,4.5,2.5,4.5),2),levels=seq(.5,6.5)))
  aa<-function(z,l,m)ifelse(outer(as.integer(z),as.integer(z),"=="),1-l,l/(m-1))
  rly<-function(z,l) {
    v<-as.numeric(as.character(z));s<-as.numeric(levels(z))
    l^abs(outer(v,v,"-"))/rowSums(l^abs(outer(v,s,"-")))
  }
  W<-aa(x$xu,.25,4)*rly(x$xo,.4)
  K<-aa(y$yu,.3,5)*rly(y$yo,.2)
  deleted<-W;diag(deleted)<-0
  truth<-sum(log(colSums(deleted*K)/colSums(deleted)))
  a<-list(xdat=x,ydat=y,bws=c(.3,.2,.25,.4),bwmethod="cv.ml",regtype="lc",
    oxkertype="racineliyan",oykertype="racineliyan",bandwidth.compute=FALSE)
  b<-do.call(npcdensbw,a)
  expected<-list(as.double(1:5),seq(.5,6.5),as.double(1:4),c(-2.5,-.5,.5,3.5))
  expect_identical(unname(.np_native_categorical_support(b)),expected)
  expect_equal(.npcdensbw_eval_only(x,y,b)$objective,truth,tolerance=2e-10)
  fit<-fitted(npcdens(b,txdat=x,tydat=y,exdat=x,eydat=y))
  expect_equal(as.numeric(fit),colSums(W*K)/colSums(W),tolerance=2e-10)
  f<-npcdensbw(yu+yo~xu+xo,data=cbind(y,x),bws=c(.3,.2,.25,.4),
    bwmethod="cv.ml",regtype="lc",oxkertype="racineliyan",
    oykertype="racineliyan",bandwidth.compute=FALSE)
  expect_identical(unname(.np_native_categorical_support(f)),expected)
  expect_equal(.npcdensbw_eval_only(x,y,f)$objective,truth,tolerance=2e-10)
  x$xu<-factor(x$xu,levels=c("b","unused","z","a"));a$xdat<-x
  reordered<-do.call(npcdensbw,a)
  expect_equal(.npcdensbw_eval_only(x,y,reordered)$objective,truth,tolerance=2e-10)
})
