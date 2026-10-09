test_that("formula subset provenance survives returned and serialized bandwidths", {
  withr::local_options(np.messages=FALSE, np.tree=FALSE)
  set.seed(3886);d<-data.frame(x=runif(60),z=runif(60));d$y<-sin(d$x)+d$z+rnorm(60,sd=.1)
  specs<-list(npreg=list(npregbw,y~x,.3),npudens=list(npudensbw,~x,.3),
    npudist=list(npudistbw,~x,.3),npcdens=list(npcdensbw,y~x,c(.3,.3)),
    npcdist=list(npcdistbw,y~x,c(.3,.3)),npindex=list(npindexbw,y~x+z,c(1,.3,.3)),
    npplreg=list(npplregbw,y~x|z,matrix(.3,2,1)),npscoef=list(npscoefbw,y~x|z,.3))
  for(family in names(specs)) {
    spec<-specs[[family]];ctor<-spec[[1]];fo<-spec[[2]];widths<-spec[[3]];fit<-get(family)
    wrap<-function(...)ctor(fo,data=d,bws=widths,bandwidth.compute=FALSE,...)
    bw<-(function(){cutoff<-.4;wrap(subset=x>cutoff)})()
    expect_identical(bw$call$subset,quote(x>cutoff),info=family)
    selected<-d[d$x>.4,];controls<-if(family=="npindex")list(se=FALSE)else list()
    # The oracle uses exactly the saved selected data via explicit native roles.
    native<-if(family %in% c("npudens","npudist"))list(tdat=selected["x"]) else
      if(family %in% c("npplreg","npscoef"))list(txdat=selected["x"],tydat=selected$y,tzdat=selected["z"]) else
      if(family=="npindex")list(txdat=selected[c("x","z")],tydat=selected$y) else
      if(family %in% c("npcdens","npcdist"))list(txdat=selected["x"],tydat=selected["y"]) else
      list(txdat=selected["x"],tydat=selected$y)
    oracle<-do.call(fit,c(list(bws=bw),native,controls))
    for(serialized in c(FALSE,TRUE)) {
      b<-if(serialized)unserialize(serialize(bw,NULL))else bw
      cutoff<-.95 # unrelated binding must never steal ownership
      got<-do.call(fit,c(list(bws=b,data=d),controls))
      expect_equal(fitted(got),fitted(oracle),tolerance=1e-10,info=family)
      replacement<-d;replacement$x<-1-d$x
      selected2<-replacement[replacement$x>.4,]
      want<-do.call(ctor,list(fo,data=selected2,bws=widths,bandwidth.compute=FALSE))
      expected<-do.call(fit,c(list(bws=want),controls))
      got2<-do.call(fit,c(list(bws=b,data=replacement),controls))
      expect_equal(fitted(got2),fitted(expected),tolerance=1e-10,info=family)
      # Default reuse must not evaluate the subset again.
      expect_equal(fitted(do.call(fit,c(list(bws=b),controls))),fitted(oracle),tolerance=1e-10)
    }
  }
})

test_that("lexical dots preserve the original subset mask and scope", {
  withr::local_options(np.messages=FALSE)
  set.seed(1);d<-data.frame(x=runif(90),z=runif(90),g=rep(1:2,each=45));d$y<-sin(d$x)+d$z
  fit<-function(dd,...)npreg(y~x,data=dd,bws=.3,bandwidth.compute=FALSE,...)
  grouped<-function(dat,...)lapply(split(dat,dat$g),function(dd)fit(dd,...))
  mapped<-function(dat,...)Map(function(dd)fit(dd,...),split(dat,dat$g))
  sized<-function(dat,...)vapply(split(dat,dat$g),function(dd)length(fitted(fit(dd,...))),0L)
  expected<-as.integer(tapply(d$x>.4,d$g,sum))
  for(fn in list(grouped,mapped))expect_identical(unname(vapply(fn(d,subset=x>.4),function(f)length(fitted(f)),0L)),expected)
  expect_identical(unname(sized(d,subset=x>.4)),expected)
  local.fit<-function(...)local(fit(d,...))
  expect_length(fitted(local.fit(subset=x>.4)),sum(expected))
  probe<-new.env();probe$n<-0L
  bw<-(function(...){npregbw(y~x,data=d,bws=.3,bandwidth.compute=FALSE,...) })(subset={probe$n<-probe$n+1L;x>.4})
  expect_identical(probe$n,1L)
  f<-npreg(bw);expect_identical(probe$n,1L)
  expect_length(fitted(f),sum(expected))
  npreg(bw,data=d);expect_identical(probe$n,2L)
  expect_error(local.fit(subset=stop("owned-subset-error")),"owned-subset-error")
})

test_that("one-call index stores replayable subset provenance", {
  withr::local_options(np.messages=FALSE)
  set.seed(1);d<-data.frame(x=runif(90),z=runif(90));d$y<-sin(d$x)+d$z
  wrap<-function(...)npindex(y~x+z,data=d,bws=c(1,.3,.3),bandwidth.compute=FALSE,se=FALSE,...)
  f<-wrap(subset=x>.4)
  expect_identical(f$bws$call$subset,quote(x>.4))
  expect_equal(fitted(npindex(f$bws,data=d,se=FALSE)),fitted(f),tolerance=0)
  b<-unserialize(serialize(f$bws,NULL));expect_equal(fitted(npindex(b,data=d,se=FALSE)),fitted(f),tolerance=0)
})
