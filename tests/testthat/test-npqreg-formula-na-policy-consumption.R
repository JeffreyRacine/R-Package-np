test_that("quantile formula replacement consumes exactly one resolved NA policy", {
  old<-options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  d<-data.frame(x=seq(-1,1,length.out=18),y=sin(seq(-1,1,length.out=18)))
  b<-npcdistbw(y~x,data=d,bws=c(.4,.5),bandwidth.compute=FALSE)
  changed<-d
  changed$y<-changed$y+.2*cos(seq_len(18))
  changed$x[c(3L,13L)]<-NA_real_
  keep<-complete.cases(changed)
  forwarded<-function(...)npqreg(...)
  for(tau in list(.5,c(.3,.7)))for(kind in c("omit","exclude","null")) {
    visits<-0L
    policy<-switch(kind,omit=function(x){visits<<-visits+1L;na.omit(x)},
      exclude=function(x){visits<<-visits+1L;na.exclude(x)},null=NULL)
    # NULL is tested on complete rows: the numerical owner independently
    # requires complete data and does not promise model.frame's NA semantics.
    input<-if(kind=="null")changed[keep,] else changed
    fit<-do.call(forwarded,list(bws=b,data=input,na.action=policy,tau=tau))
    expect_equal(visits,if(kind=="null")0L else 1L)
    oracle<-npqreg(b,txdat=changed[keep,"x",drop=FALSE],
      tydat=changed[keep,"y",drop=FALSE],tau=tau)
    value<-fitted(fit)
    # npqreg's existing output owner restores omitted evaluation positions
    # for both policies; this repair does not change that row contract.
    if(kind!="null") {
      expect_true(all(is.na(if(is.null(dim(value)))value[!keep] else value[!keep,,drop=FALSE])))
      value<-if(is.null(dim(value)))value[keep] else value[keep,,drop=FALSE]
    }
    expect_identical(dim(value),dim(fitted(oracle)))
    expect_equal(as.numeric(value),as.numeric(fitted(oracle)),tolerance=0)
    expect_equal(fit$ntrain,sum(keep))
    expect_equal(visits,if(kind=="null")0L else 1L)
    nd<-data.frame(x=c(-.6,.1,.7))
    expect_equal(unname(predict(fit,newdata=nd,tau=tau)),
      unname(fitted(npqreg(b,txdat=changed[keep,"x",drop=FALSE],
        tydat=changed[keep,"y",drop=FALSE],exdat=nd,tau=tau))),tolerance=0)
    expect_equal(visits,if(kind=="null")0L else 1L)
  }
  rejected<-0L
  reject<-function(x){rejected<<-rejected+1L;stop("policy rejection sentinel")}
  expect_error(npqreg(b,data=changed,na.action=reject),"policy rejection sentinel")
  expect_equal(rejected,1L)
  expect_error(npqreg(b,data=d,na.action=na.omit,not.a.control=1),"unused argument")
  expect_error(npqreg(b,txdat=d["x"],tydat=d["y"],na.action=na.omit),"unused argument")
})
