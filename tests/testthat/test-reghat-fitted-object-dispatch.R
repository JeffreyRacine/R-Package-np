test_that("fitted regression objects retain the registered hat input owner", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d<-data.frame(x=seq(-1,1,length.out=24),y=sin(seq(-1,1,length.out=24)))
  for(route in c("formula","native")) for(type in c("lc","ll"))
    for(bwtype in c("fixed","generalized_nn","adaptive_nn")) {
      aa<-list(bws=if(bwtype=="fixed").4 else 6,bwtype=bwtype,regtype=type,
        bandwidth.compute=FALSE)
      if(route=="formula") {aa$formula<-y~x;aa$data<-d}
      else {aa$xdat<-d["x"];aa$ydat<-d$y}
      b<-do.call(npregbw,aa)
      model<-npreg(b)
      ref<-npreghat(b,txdat=d["x"],y=d$y)
      for(output in c("matrix","apply","constraint")) {
        calls<-0L
        once<-function(){calls<<-calls+1L;model}
        actual<-npreghat(once(),output=output)
        oracle<-npreghat(b,txdat=d["x"],y=d$y,output=output)
        expect_identical(calls,1L)
        expect_identical(dim(actual),dim(oracle))
        expect_equal(as.numeric(actual),as.numeric(oracle),tolerance=0)
      }
      if(type=="lc" && bwtype=="fixed") {
        w<-outer(d$x,d$x,function(a,z)dnorm((a-z)/.4))
        literal<-t(sweep(w,2,colSums(w),"/"))
        expect_equal(as.numeric(ref),as.numeric(literal),tolerance=2e-12)
      }
      changed<-d$y+1
      expect_equal(as.numeric(npreghat(model,txdat=d["x"],y=changed,output="apply")),
        as.numeric(npreghat(b,txdat=d["x"],y=changed,output="apply")),tolerance=0)
      inherited<-model;class(inherited)<-c("r21model",class(model))
      expect_equal(as.numeric(npreghat(inherited)),as.numeric(ref),tolerance=0)
    }
})
