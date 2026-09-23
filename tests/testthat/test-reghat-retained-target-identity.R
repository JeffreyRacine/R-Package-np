test_that("hat recomputation preserves retained target identity", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  for(bwtype in c("fixed","generalized_nn","adaptive_nn"))
    for(kernel in c("gaussian","beta")) for(duplicates in c(FALSE,TRUE)) {
      x<-data.frame(x=if(duplicates)rep(seq(.08,.92,length.out=12),each=2) else
        seq(.08,.92,length.out=24)^1.2)
      y<-sin(8*x$x)+cos(2*x$x)
      a<-list(xdat=x,ydat=y,bws=if(bwtype=="fixed").2 else 9L,
        bwtype=bwtype,ckertype=kernel,regtype="ll",bandwidth.compute=FALSE)
      if(kernel=="beta")a<-c(a,list(ckerbound="fixed",ckerlb=0,ckerub=1))
      b<-do.call(npregbw,a)
      h<-npreghat(b,txdat=x)
      same<-predict(h,deriv=0L)
      expect_identical(attr(same,"trainiseval"),TRUE)
      expect_equal(as.numeric(same),as.numeric(h),tolerance=0)
      expect_identical(predict(h),h)
      ref<-npreghat(b,txdat=x,s=1L)
      actual<-predict(h,s=1L)
      expect_identical(attr(actual,"trainiseval"),TRUE)
      expect_equal(as.numeric(actual),as.numeric(ref),tolerance=0)
      expect_equal(as.numeric(predict(h,s=1L,output="apply",y=y)),
        as.numeric(ref%*%y),tolerance=2e-10)
      ext<-npreghat(b,txdat=x,exdat=x,s=1L)
      given<-predict(h,newdata=x,s=1L)
      expect_identical(attr(given,"trainiseval"),FALSE)
      expect_equal(as.numeric(given),as.numeric(ext),tolerance=0)
      explicit<-predict(h,exdat=x,s=1L)
      expect_equal(as.numeric(explicit),as.numeric(ext),tolerance=0)
      replay<-predict(ext,deriv=1L)
      expect_identical(attr(replay,"trainiseval"),FALSE)
      expect_equal(as.numeric(replay),as.numeric(ext),tolerance=0)
      legacy<-ext;attr(legacy,"trainiseval")<-NULL
      expect_equal(as.numeric(predict(legacy,deriv=1L)),as.numeric(ext),tolerance=0)
    }
})

test_that("hat LOO replay omits only implicit evaluation coordinates", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(.1,.9,length.out=24));y<-sin(x$x)
  for(bwtype in c("fixed","generalized_nn","adaptive_nn")) {
    b<-npregbw(xdat=x,ydat=y,bws=if(bwtype=="fixed").2 else 9L,
      bwtype=bwtype,regtype="ll",bandwidth.compute=FALSE)
    h<-npreghat(b,txdat=x,leave.one.out=TRUE)
    actual<-predict(h,deriv=0L)
    expect_identical(attr(actual,"trainiseval"),TRUE)
    expect_equal(as.numeric(actual),as.numeric(h),tolerance=0)
    expect_error(predict(h,newdata=x,deriv=0L),"leave.one.out")
  }
})
