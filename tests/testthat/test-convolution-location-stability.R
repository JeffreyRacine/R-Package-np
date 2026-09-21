convolution_kernel_oracle <- function(z, kernel, order) {
  if(kernel=="gaussian")
    return(dnorm(z)*switch(as.character(order),"2"=1,"4"=1.5-.5*z^2,
      "6"=1.875-1.25*z^2+.125*z^4,"8"=(105-105*z^2+21*z^4-z^6)/48))
  z2<-z*z
  value<-switch(as.character(order),
    "2"=.33541019662496845446-.067082039324993690892*z2,
    "4"=.008385254916*(-15+7*z2)*(-5+z2),
    "6"=.33541019662496845446*(2.734375+z2*(-3.28125+.721875*z2))*(1-.2*z2),
    "8"=.33541019662496845446*(3.5888671875+z2*(-7.8955078125+z2*(4.1056640625-.5865234375*z2)))*(1-.2*z2))
  ifelse(abs(z)<sqrt(5),value,0)
}

test_that("convolution pairs agree with centered integration after translation", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  z<-c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1);n<-length(z)
  for(type in c("fixed","generalized_nn","adaptive_nn"))
   for(kernel in c("gaussian","epanechnikov"))
    for(order in if(kernel=="gaussian")8L else c(2L,4L,6L,8L)) {
      h<-if(type=="fixed")rep(.31,n)else vapply(z,function(x)sort(abs(z-x))[4L],0)
      # ANN donor radii exclude the donor itself; raw external-query radii
      # include an exact match. These are deliberately different contracts.
      ht<-if(type=="adaptive_nn")vapply(z,function(x)sort(abs(z-x))[5L],0)else h
      C<-outer(seq_len(n),seq_len(n),Vectorize(function(i,j) {
        delta<-z[i]-z[j]
        lower<-if(kernel=="gaussian")-Inf else max(-sqrt(5)*ht[i],-delta-sqrt(5)*h[j])
        upper<-if(kernel=="gaussian")Inf else min(sqrt(5)*ht[i],-delta+sqrt(5)*h[j])
        if(lower>=upper)return(0)
        integrate(function(t)convolution_kernel_oracle(t/ht[i],kernel,order)*
          convolution_kernel_oracle((t+delta)/h[j],kernel,order),lower,upper,
          rel.tol=1e-10,abs.tol=1e-12)$value
      }))
      for(tree in c(FALSE,TRUE))for(shift in c(0,10,100))for(scale in c(1,1e-8,1e8)) {
        options(np.tree=tree)
        x<-data.frame(z=(z+shift)*scale)
        a<-list(txdat=x,exdat=x,bws=if(type=="fixed").31*scale else 4,
                bwtype=type,ckertype=kernel,ckerorder=order,operator="convolution",
                return.kernel.weights=TRUE)
        got<-do.call(npksum,a)
        expect_equal(got$kw/scale,C,tolerance=2e-10)
        expect_equal(as.numeric(got$ksum)/scale,colSums(C),tolerance=2e-10)
      }
    }
})

test_that("density CVLS is the independent overlap minus delete-one cross term", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  z<-c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1);h<-.31;n<-length(z)
  for(kernel in c("gaussian","epanechnikov"))for(order in c(2L,4L,6L,8L)) {
    C<-outer(z,z,Vectorize(function(a,b) {
      lo<-if(kernel=="gaussian")-Inf else max(-sqrt(5)*h,b-a-sqrt(5)*h)
      hi<-if(kernel=="gaussian")Inf else min(sqrt(5)*h,b-a+sqrt(5)*h)
      if(lo>=hi)return(0)
      integrate(function(t)convolution_kernel_oracle(t/h,kernel,order)*
        convolution_kernel_oracle((t+a-b)/h,kernel,order)/h^2,lo,hi,
        rel.tol=1e-10,abs.tol=1e-11)$value
    }))
    N<-convolution_kernel_oracle(outer(z,z,"-")/h,kernel,order)/h;diag(N)<-0
    expected<-mean(C)-2*sum(N)/(n*(n-1))
    for(shift in c(0,100)) {
      x<-data.frame(z=z+shift)
      bw<-npudensbw(dat=x,bws=h,ckertype=kernel,ckerorder=order,
                     bwmethod="cv.ls",bandwidth.compute=FALSE)
      got<-npudensbw(dat=x,bws=bw,eval.only=TRUE)
      expect_equal(-as.numeric(got$fval),expected,tolerance=2e-9)
    }
  }
})

test_that("conditional CVLS convolution remains translation invariant", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(88921);n<-19L;x<-data.frame(x=runif(n,-1,1));y<-sin(x$x)+rnorm(n,sd=.3)
  for(kernel in c("gaussian","epanechnikov"))for(order in c(2L,4L,6L,8L)) {
    a<-list(xdat=x,ydat=y,bws=c(.38,.45),bwmethod="cv.ls",bandwidth.compute=FALSE,
            cxkertype=kernel,cykertype=kernel,cxkerorder=order,cykerorder=order)
    bw<-do.call(npcdensbw,a)
    base<-.npcdensbw_eval_only(bws=bw,xdat=x,ydat=data.frame(y))$objective
    shifted<-.npcdensbw_eval_only(bws=bw,xdat=data.frame(x=x$x+50),
                                  ydat=data.frame(y=y+100))$objective
    expect_equal(shifted,base,tolerance=2e-9)
  }
})
