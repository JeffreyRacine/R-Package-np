uniform_overlap_oracle <- function(a, b, ha, hb)
  pmax(0, pmin(a+ha,b+hb)-pmax(a-ha,b-hb))/4

test_that("uniform convolution keeps bandwidth ownership at its caller", {
  old<-options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  z<-c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1);n<-length(z)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    he<-if(type=="fixed")rep(.31,n)else vapply(z,function(v)sort(abs(z-v))[4],0)
    ht<-if(type=="adaptive_nn")vapply(z,function(v)sort(abs(z-v))[5],0)else he
    C<-outer(seq_len(n),seq_len(n),Vectorize(function(i,j)
      uniform_overlap_oracle(z[i],z[j],ht[i],he[j])))
    for(tree in c(FALSE,TRUE))for(shift in c(0,100))for(scale in c(.2,1,3)) {
      options(np.tree=tree);x<-data.frame(z=(z+shift)*scale)
      for(divide in c(FALSE,TRUE)) {
        got<-npksum(txdat=x,exdat=x,bws=if(type=="fixed").31*scale else 4,
          bwtype=type,ckertype="uniform",operator="convolution",
          bandwidth.divide=divide,return.kernel.weights=TRUE)
        expected<-if(divide)C/outer(ht,he)/scale else C*scale
        expect_equal(as.numeric(got$ksum),colSums(expected),tolerance=2e-12)
        if(!divide)expect_equal(got$kw,C*scale,tolerance=2e-12)
      }
    }
  }
  x<-data.frame(a=z,b=rev(z));h<-c(.29,.71)
  C1<-outer(z,z,Vectorize(function(a,b)uniform_overlap_oracle(a,b,h[1],h[1])))
  C2<-outer(rev(z),rev(z),Vectorize(function(a,b)uniform_overlap_oracle(a,b,h[2],h[2])))
  got<-npksum(txdat=x,exdat=x,bws=h,ckertype="uniform",
    operator=c("convolution","convolution"),bandwidth.divide=TRUE)
  expect_equal(as.numeric(got$ksum),colSums(C1*C2)/prod(h)^2,tolerance=2e-12)
})

test_that("uniform density CVLS agrees with independent support overlaps", {
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  z<-c(-.8,-.45,-.2,.03,.22,.47,.68,.91,1.1);n<-length(z)
  for(h in c(.2,.7,1,2)) {
    C<-outer(z,z,Vectorize(function(a,b)uniform_overlap_oracle(a,b,h,h)))/h^2
    N<-(abs(outer(z,z,"-"))<h)/2/h;diag(N)<-0
    expected<-mean(C)-2*sum(N)/(n*(n-1))
    bw<-npudensbw(dat=data.frame(z),bws=h,ckertype="uniform",bwmethod="cv.ls",bandwidth.compute=FALSE)
    actual<--npudensbw(dat=data.frame(z),bws=bw,eval.only=TRUE)$fval
    expect_equal(as.numeric(actual),expected,tolerance=2e-13)
  }
})

test_that("uniform conditional CVLS uses one response normalization", {
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(2481);n<-17L;x<-runif(n,-1,1);y<-sin(x)+rnorm(n,sd=.3)
  hx<-.52;hy<-.39
  wx<-(abs(outer(x,x,"-"))<hx)/2/hx;diag(wx)<-0
  wy<-(abs(outer(y,y,"-"))<hy)/2/hy
  C<-outer(y,y,Vectorize(function(a,b)uniform_overlap_oracle(a,b,hy,hy)))/hy^2
  expected<-mean(vapply(seq_len(n),function(j) {
    w<-wx[,j]/sum(wx[,j]);sum(w*(C%*%w))-2*sum(w*wy[,j])
  },0))
  for(scale in c(.3,1,4)) {
    xx<-data.frame(x);yy<-data.frame(y=y*scale)
    bw<-npcdensbw(xdat=xx,ydat=yy,bws=c(hy*scale,hx),regtype="lc",
      cxkertype="uniform",cykertype="uniform",bwmethod="cv.ls",bandwidth.compute=FALSE)
    actual<--.npcdensbw_eval_only(xx,yy,bw)$objective
    expect_equal(as.numeric(actual),expected/scale,tolerance=2e-12)
  }
})
