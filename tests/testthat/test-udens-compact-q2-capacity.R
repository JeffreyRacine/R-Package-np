test_that("compact conditional covariance retains the complete quadrature law", {
  # Literal weighted sums are independent of the native conditional QR owner.
  for(nu in c(4L,5L,6L,11L)) for(law in 0:1) {
    nv<-7L; nx<-2L;ny<-2L
    u<-seq(.05,.95,length.out=nu);v<-seq(.08,.92,length.out=nv)
    a<-outer(1:2,u,function(i,x)sin(3*i*x));b<-outer(1:2,v,function(i,x)cos(4*i*x))
    pa<-array(0,c(nx,nu,24));pb<-array(0,c(ny,nv,24))
    for(k in 1:24) {
      pa[,,k]<-if(k<=8)a else if((k-1)%%4==0)1 else 0
      pb[,,k]<-if(k<=4)b else if(k<=8)1 else if((k-1)%%4==0)-.25*b else 0
    }
    ix<-list(c(2L,1L,2L,1L,1L),c(1L,2L,2L,1L,1L));m<-length(ix[[1]])
    axes<-list(list(lv=log(u),l1v=log1p(-u),prob=rep(1/nu,nu)),
               list(lv=log(v),l1v=log1p(-v),prob=rep(1/nv,nv)))
    ratio<-outer(u,v,function(x,y).8*sin(5*x*y));K<-matrix(sin(1:(m*5))/5,m,5)
    prep<-list(cache=list(parts=list(pa,pb),index=ix),axes=axes,theta=c(0,0,0,0,.4),
      logz=.2,x=matrix(0,20,2),ratio=ratio,Iinv=diag(5))
    id<-as.matrix(expand.grid(seq_len(nu),seq_len(nv)));ui<-id[,1];vi<-id[,2]
    H<-a[ix[[1]],ui,drop=FALSE]*b[ix[[2]],vi,drop=FALSE]+a[ix[[1]],ui,drop=FALSE]+b[ix[[2]],vi,drop=FALSE]
    T<-cbind(log(u[ui]),log1p(-u[ui]),log(v[vi]),log1p(-v[vi]),(2*u[ui]-1)*(2*v[vi]-1))
    lw<-.4*T[,5]-.2+law*ratio[id];w<-exp(lw-max(lw));w<-w/sum(w)
    Z<-H-K%*%t(T);Z<-Z-drop(Z%*%w);expected<-tcrossprod(sweep(Z,2,sqrt(w),'*'))/20
    actual<-.np_udens_compact_q2_covariance_native(prep,K,law)
    expect_equal(actual,expected,tolerance=1e-12)
    expect_equal(actual,.np_udens_beta_covariance_native(prep,K,law),tolerance=1e-12)
    expect_identical(actual,t(actual))
    raw<-.np_udens_compact_q2_moments_native(prep$cache,axes,.4,.2,ratio,FALSE)[[law+1L]]
    expect_equal(raw$mean,drop(rbind(H,t(T))%*%w),tolerance=1e-12)
  }
  prep$cache$parts[[1]][1,1,24]<-1
  expect_error(.np_udens_compact_q2_covariance_native(prep,K,1L),"feature encoding")
  expect_error(.np_udens_compact_q2_moments_native(prep$cache,axes,.4,.2,ratio,FALSE),"feature encoding")
})
