r21_scoef_adjoint <- function(W, E, kw, ridge) {
  H<-matrix(0,nrow(E),nrow(W))
  for(i in seq_len(nrow(E))){
    a<-crossprod(W,W*kw[,i]);v<-solve(t(a+diag(ridge,ncol(W))),E[i,])
    if(ridge>0)v[1L]<-v[1L]+ridge*v[1L]/a[1L,1L]
    H[i,]<-kw[,i]*drop(W%*%v)
  }
  H
}

test_that("smooth-coefficient level hats retain the categorical ridge correction", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(836);n<-30L
  z<-data.frame(z=ordered(rep(0:2,each=10),levels=0:2))
  x<-data.frame(x=rep(c(1,3,-1),each=10));y<-rnorm(n)+x$x;yy<-cbind(y,y+1)
  W<-cbind(1,x$x);kw<-outer(as.integer(z$z),as.integer(z$z),`==`)*1
  for(regtype in c("lc","lp"))for(compress in c(FALSE,TRUE)) {
    options(np.categorical.compress=compress)
    a<-list(xdat=x,zdat=z,ydat=y,bws=0,regtype=regtype,bandwidth.compute=FALSE)
    if(regtype=="lp")a$degree<-0L
    b<-do.call(npscoefbw,a)
    H<-npscoefhat(b,txdat=x,tzdat=z,output="matrix")
    expected<-r21_scoef_adjoint(W,W,kw,1/n)
    fit<-npscoef(b,txdat=x,tzdat=z,tydat=y,iterate=FALSE,se=FALSE)
    expect_equal(H,expected,tolerance=3e-12)
    expect_equal(as.vector(H%*%y),fitted(fit),tolerance=3e-12)
    expect_equal(rowSums(H),rep(1,n),tolerance=3e-12)
    expect_equal(npscoefhat(b,txdat=x,tzdat=z,y=y,output="apply"),drop(H%*%y),tolerance=3e-12)
    expect_equal(npscoefhat(b,txdat=x,tzdat=z,y=yy,output="apply"),H%*%yy,tolerance=3e-12,ignore_attr=TRUE)
    expect_equal(npscoefhat(b,txdat=x,tzdat=z,y=y,output="constraint"),t(H)*y,tolerance=0)
  }
})

test_that("smooth-coefficient profile and LOO ridge units agree", {
  old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  n<-18L;z<-data.frame(z=factor(rep(letters[1:3],each=6)))
  x<-data.frame(x=sin(seq_len(n)));y<-x$x+cos(seq_len(n)); W<-cbind(1,x$x)
  for(lambda in c(0,.35)){
    b<-npscoefbw(xdat=x,zdat=z,ydat=y,bws=lambda,regtype="lc",ukertype="liracine",bandwidth.compute=FALSE)
    kw<-ifelse(outer(as.integer(z$z),as.integer(z$z),`==`),1,lambda)/(1+2*lambda)
    for(loo in c(FALSE,TRUE))for(compress in c(FALSE,TRUE)){
      options(np.categorical.compress=compress)
      k<-kw;if(loo)diag(k)<-0
      expected<-r21_scoef_adjoint(W,W,k,.2)
      H<-npscoefhat(b,txdat=x,tzdat=z,output="matrix",ridge=.2,leave.one.out=loo)
      expect_equal(H,expected,tolerance=3e-12)
      expect_equal(npscoefhat(b,txdat=x,tzdat=z,y=y,output="apply",ridge=.2,leave.one.out=loo),
        drop(expected%*%y),tolerance=3e-12)
    }
  }
})

test_that("smooth-coefficient tensor and explicit ridge hats use the same adjoint", {
  old<-options(np.messages=FALSE,np.categorical.compress=FALSE);on.exit(options(old),add=TRUE)
  ns<-asNamespace(getNamespaceName(environment(npscoef)))
  state<-get(".npscoef_lp_state",ns);tensor<-get(".npscoef_row_tensor_design",ns)
  kernel<-get(".np_kernel_weights_direct",ns)
  set.seed(837);n<-19L
  x<-data.frame(x=rnorm(n));z<-data.frame(z=seq(.05,.95,length.out=n));y<-sin(4*z$z)+x$x
  for(degree in 0:2) for(bwtype in c("fixed","generalized_nn","adaptive_nn")) {
    b<-npscoefbw(xdat=x,zdat=z,ydat=y,bws=if(bwtype=="fixed").3 else 10,
      regtype="lp",degree=degree,bwtype=bwtype,bandwidth.compute=FALSE)
    ex<-x[3:5,,drop=FALSE];ez<-z[3:5,,drop=FALSE]
    st<-state(b,z,ez)
    W<-if(degree==0)cbind(1,x$x)else tensor(cbind(1,x$x),st$W.train)
    E<-if(degree==0)cbind(1,ex$x)else tensor(cbind(1,ex$x),st$W.eval)
    kw<-kernel(bws=if(degree==0)b else st$rbw,txdat=z,exdat=ez,bandwidth.divide=TRUE)
    for(rho in c(0,.2)){
      expected<-r21_scoef_adjoint(W,E,kw,rho)
      H<-npscoefhat(b,txdat=x,tzdat=z,exdat=ex,ezdat=ez,output="matrix",ridge=rho)
      expect_equal(H,expected,tolerance=3e-12)
      expect_equal(npscoefhat(b,txdat=x,tzdat=z,exdat=ex,ezdat=ez,y=cbind(y,y+1),output="apply",ridge=rho),
        expected%*%cbind(y,y+1),tolerance=3e-12,ignore_attr=TRUE)
    }
  }
})
