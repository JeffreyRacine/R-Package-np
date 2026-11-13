test_that("first-derivative fits and extractors reject unsupported explicit orders", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  for(family in c("npreg","npcdens","npcdist"))
    for(type in c("lc","ll")) for(mixed in c(FALSE,TRUE)) {
      d<-data.frame(x=seq(-1,1,length.out=24),y=sin(seq(-1,1,length.out=24)))
      if(mixed)d$f<-factor(rep(c("a","b"),12))
      bw.args<-list(formula=if(mixed)y~x+f else y~x,data=d,
        bws=if(family=="npreg")c(.5,if(mixed).2) else c(.4,.5,if(mixed).2),
        regtype=type,bandwidth.compute=FALSE)
      b<-do.call(get(paste0(family,"bw")),bw.args)
      a<-list(bws=b,gradients=TRUE,se=TRUE)
      fit<-do.call(get(family),a)
      partial<-family!="npreg" && type=="ll" && mixed
      for(uncertainty in c(FALSE,TRUE)) {
        value<-if(family=="npreg") {
          if(uncertainty)fit$gerr else fit$grad
        } else if(uncertainty)fit$congerr else fit$congrad
        expect_identical(gradients(fit,se=uncertainty),value)
        expect_identical(gradients(fit,se=uncertainty,gradient.order=1L),value)
        if(partial) {
          expect_warning(masked<-gradients(fit,se=uncertainty,gradient.order=2L),
            "requested order 2.*degree 1")
          expect_true(all(is.na(masked[,1L])))
          expect_identical(masked[,2L],value[,2L])
        } else expect_error(gradients(fit,se=uncertainty,gradient.order=2L),
          "only first derivatives|only for regtype|no available derivative")
        expect_error(gradients(fit,se=uncertainty,gradient.order=0L),
          "finite positive integers")
      }
      if(partial) {
        expect_warning(partial.fit<-do.call(get(family),c(a,list(gradient.order=2L))),
          "requested order 2.*degree 1")
        expect_true(all(is.na(partial.fit$congrad[,1L])))
        expect_identical(partial.fit$congrad[,2L],fit$congrad[,2L])
        expect_warning(predict(fit,gradients=TRUE,gradient.order=2L),
          "requested order 2.*degree 1")
        expect_warning(plot(fit,gradients=TRUE,gradient.order=2L,
          errors="none",output="data",neval=4L),"requested order 2.*degree 1")
      } else {
        expect_error(do.call(get(family),c(a,list(gradient.order=2L))),
        "only first derivatives|only for regtype|no available derivative")
      expect_error(predict(fit,gradients=TRUE,gradient.order=2L),
        "only first derivatives|only for regtype|no available derivative")
      expect_error(plot(fit,gradients=TRUE,gradient.order=2L,
        errors="none",output="data",neval=4L),
        "only first derivatives|only for regtype|no requested component")
      }
    }
})

test_that("categorical first differences keep their independent order contract", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d<-data.frame(y=sin(seq_len(24)),f=factor(rep(c("a","b"),12)))
  for(family in c("npreg","npcdens","npcdist")) {
    b<-do.call(get(paste0(family,"bw")),list(formula=y~f,data=d,
      bws=if(family=="npreg").2 else c(.4,.2),
      regtype="lc",bandwidth.compute=FALSE))
    fit<-do.call(get(family),list(bws=b,gradients=TRUE,se=TRUE))
    for(uncertainty in c(FALSE,TRUE)) {
      expect_identical(gradients(fit,se=uncertainty,gradient.order=1L),
        gradients(fit,se=uncertainty))
      expect_error(gradients(fit,se=uncertainty,gradient.order=2L),"first differences")
    }
  }
})

test_that("LP extraction retains higher orders and partial availability without refitting", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d<-data.frame(x=seq(-1,1,length.out=32),z=sin(seq_len(32)*2))
  d$y<-d$x^2+d$z+.1*cos(seq_len(32))
  for(family in c("npreg","npcdens","npcdist")) {
    b<-do.call(get(paste0(family,"bw")),list(formula=y~x+z,data=d,
      bws=if(family=="npreg")c(.6,.7) else c(.5,.6,.7),
      regtype="lp",degree=c(2L,1L),bandwidth.compute=FALSE))
    fit<-do.call(get(family),list(bws=b,gradients=TRUE,se=TRUE,gradient.order=c(2L,1L)))
    for(uncertainty in c(FALSE,TRUE)) {
      value<-gradients(fit,se=uncertainty)
      expect_identical(gradients(fit,se=uncertainty,gradient.order=c(2L,1L)),value)
      expect_warning(masked<-gradients(fit,se=uncertainty,gradient.order=c(2L,2L)),
        "requested order 2.*degree 1")
      expect_identical(masked[,1L],value[,1L])
      expect_true(all(is.na(masked[,2L])))
      expect_error(gradients(fit,se=uncertainty,gradient.order=c(1L,1L)),
        "differs from the derivative order stored")
    }
  }
})
