test_that("evaluation outcomes cannot expand or reorder fitted class support", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  for(ordered in c(FALSE,TRUE)) for(route in c("native","formula"))
    for(proper in c(FALSE,TRUE)) {
      lev<-if(ordered)c("0","2","5") else c("a","b","c")
      extra<-if(ordered)"9" else "d"
      d<-data.frame(y=factor(rep(lev,12),levels=lev,ordered=ordered),
                    x=seq(-1,1,length.out=36))
      b<-npcdensbw(y~x,data=d,bws=c(.25,.5),bandwidth.compute=FALSE,
                    oykertype="racineliyan")
      e<-data.frame(x=c(-.7,-.1,.35,.8))
      call<-list(bws=b,probabilities=TRUE,gradients=TRUE,se=TRUE,
                 level=lev[2],proper=proper)
      if(route=="native") {
        call$txdat<-d["x"];call$tydat<-d["y"];call$exdat<-e
      } else call$newdata<-e
      ref<-do.call(npconmode,call)
      expect_identical(colnames(ref$probabilities),lev)
      if(!ordered) {
        # Independent LC Gaussian / Aitchison-Aitken conditional mass.
        wx<-outer(d$x,e$x,function(x,z)dnorm((z-x)/.5))
        wy<-outer(as.character(d$y),lev,function(y,z)ifelse(y==z,.75,.125))
        oracle<-crossprod(wx,wy)/colSums(wx)
        expect_equal(unname(ref$probabilities),unname(oracle),tolerance=2e-12)
      }
      truth<-c(lev,lev[1])
      variants<-list(
        factor(truth,levels=c(lev,extra),ordered=ordered),
        factor(truth,levels=if(ordered)c(lev,extra) else c(extra,rev(lev)),ordered=ordered),
        factor(c(lev,extra),levels=c(lev,extra),ordered=ordered))
      for(y in variants) {
        args<-call
        if(route=="native")args$eydat<-data.frame(y=y)
        else args$newdata<-data.frame(e,y=y)
        fit<-do.call(npconmode,args)
        expect_identical(colnames(fit$probabilities),lev)
        expect_equal(fit$probabilities,ref$probabilities,tolerance=0)
        expect_equal(fit$probability.errors,ref$probability.errors,tolerance=0)
        expect_equal(fit$probability.gradients,ref$probability.gradients,tolerance=0)
        expect_equal(fitted(fit),fitted(ref),tolerance=0)
        expect_equal(sum(fit$confusion.matrix),nrow(e))
        expect_equal(fit$CCR.overall,mean(as.character(predict(fit,type="class"))==as.character(y)))
        expect_equal(predict(fit,type="prob"),ref$probabilities,tolerance=0)
      }
      expect_error(do.call(npconmode,c(call[names(call)!="level"],list(level=extra))),
                   "'level' must identify one response level",fixed=TRUE)
    }
})
