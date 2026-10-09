# NA01: evaluation policy is independent of retained training and scoring rows.
local({
  old <- options(np.messages = FALSE, na.action = 'na.omit')
  on.exit(options(old))
  set.seed(731)
  d <- data.frame(x = rnorm(160), z = rnorm(160), g = factor(rep(letters[1:3], length.out = 160)))
  d$y <- sin(d$x + .35*d$z) + rnorm(160, sd = .35)
  d$yb <- rbinom(160, 1, plogis(d$x + .35*d$z))
  make <- function(method='ichimura', type='fixed') {
    dd <- d; dd$y <- if (method == 'ichimura') d$y else d$yb
    b <- npindexbw(y ~ x + z, data = dd, bws = c(1, .35, if(type=='fixed') .7 else 25),
                   bandwidth.compute = FALSE, method = method, bwtype = type)
    npindex(b, se = FALSE)
  }
  fields <- c('mean','merr','index','grad','gerr','mean.grad','betavcov','nobs',
              'ntrain','omit','rows.omit','nobs.omit','MSE','CCR.overall',
              'diagnostics.sample','diagnostics.nobs','training.MSE')
  same <- function(a,b) for(n in fields) expect_equal(a[[n]],b[[n]],info=n)

  test_that('explicit evaluation action wins without dropping scoring-only NAs', {
    for(method in c('ichimura','kleinspady')) for(global in c('na.omit','na.exclude')) {
      options(na.action=global);f<-make(method);e<-d[1:15,];e$y<-if(method=='ichimura')e$y else e$yb
      e$x[c(2,7)]<-NA;e$y[c(3,9)]<-NA
      ref<-predict(f,newdata=e[-c(2,7),],y.eval=FALSE)
      omitted<-predict(f,newdata=e,na.action=na.omit)
      excluded<-predict(f,newdata=e,na.action=na.exclude)
      expect_equal(omitted,ref);expect_length(excluded,15L)
      expect_true(all(is.na(excluded[c(2,7)])));expect_equal(excluded[-c(2,7)],ref)
      expect_true(all(is.finite(excluded[c(3,9)])))
      expect_error(predict(f,newdata=e,na.action=na.fail),'missing values')
      fit<-npindex(f$bws,newdata=e,se=FALSE,na.action=na.exclude)
      expect_equal(fit$diagnostics.nobs,11L);expect_equal(fit$ntrain,160L)
      # An existing global-policy call is an exact object-convention oracle.
      options(na.action='na.exclude');same(fit,npindex(f$bws,newdata=e,se=FALSE))
    }
  })
  test_that('omitted policies and native evaluation precedence remain unchanged', {
    f<-make();e<-d[1:15,];e$x[2]<-NA
    for(action in c('na.omit','na.exclude')){
      options(na.action=action);a<-predict(f,newdata=e);b<-predict(f,newdata=e,na.action=get(action,asNamespace('stats')))
      expect_equal(a,b)
    }
    options(na.action='na.omit');x<-data.frame(x=c(-.3,.4),z=c(.2,.5))
    expect_equal(predict(f,newdata=e,exdat=x,na.action=na.fail),predict(f,exdat=x))
    for(policy in list(NULL,na.pass)) {
      actual<-npindex(f$bws,newdata=e,se=FALSE,na.action=policy)
      expect_equal(actual$rows.omit,2L);expect_equal(actual$nobs.omit,1L)
      expect_equal(fitted(actual),predict(f,newdata=e[-2,],na.action=na.fail))
      expect_length(fitted(actual),14L)
    }
  })
  test_that('policy function values are realized once and training stays retained', {
    options(na.action='na.omit');f<-make();e<-d[1:15,];e$x[2]<-NA
    state<-new.env();state$force<-state$apply<-0L
    policy<-function(x){state$apply<-state$apply+1L;stats::na.exclude(x)}
    p<-predict(f,newdata=e,na.action={state$force<-state$force+1L;policy})
    expect_length(p,15L);expect_equal(state$force,1L);expect_equal(state$apply,1L)
    expect_equal(p,predict(f,newdata=e,na.action='na.exclude'))
    expect_equal(f$bws$nobs,160L)
  })
  test_that('combined replacement training and response do not consume evaluation action', {
    options(na.action='na.omit');f<-make();tr<-d[1:110,];tr$x[4]<-NA;ev<-d[121:135,];ev$x[2]<-NA
    z<-npindex(f$bws,data=tr,newdata=ev,se=FALSE,na.action=na.exclude)
    expect_length(fitted(z),15L);expect_true(is.na(fitted(z)[2]));expect_equal(z$ntrain,109L)
    expect_error(npindex(f$bws,data=tr,tydat=tr$y,newdata=ev,se=FALSE),
                 "partial response replacement cannot be combined")
    yy<-d$y+1;z2<-npindex(f$bws,tydat=yy,newdata=ev,se=FALSE,na.action=na.exclude)
    reference<-npindex(f$bws,newdata=ev,se=FALSE,na.action=na.exclude)
    expect_length(fitted(z2),15L);expect_equal(z2$ntrain,160L);expect_true(is.na(fitted(z2)[2]))
    expect_equal(fitted(z2)[-2],fitted(reference)[-2]+1,tolerance=1e-10)
  })
  test_that('inference and gradients retain training provenance and omission maps', {
    options(na.action='na.omit');f<-make();e<-d[1:15,];e$x[2]<-NA;e$y[3]<-NA
    z<-npindex(f$bws,newdata=e,na.action=na.exclude,se=TRUE,gradients=TRUE)
    expect_length(se(z),15L);expect_true(is.na(se(z)[2]));expect_equal(NROW(gradients(z)),15L)
    expect_true(all(is.na(gradients(z)[2,])));options(na.action='na.exclude')
    ref<-npindex(f$bws,newdata=e,se=TRUE,gradients=TRUE);same(z,ref)
    a<-predict(f,newdata=e,na.action=na.exclude,se.fit=TRUE);b<-predict(f,newdata=e,se.fit=TRUE)
    expect_equal(a,b);expect_true(is.finite(a$residual.scale))
    expect_output(summary(z),'Evaluation');expect_output(print(z),'Single Index')
  })
  test_that('transformed factor designs and serialization preserve evaluation rows', {
    options(na.action='na.omit');b<-npindexbw(y~I(x^2)+g,data=d,bws=c(1,.25,-.15,.7),bandwidth.compute=FALSE)
    f<-unserialize(serialize(npindex(b,se=FALSE),NULL));e<-d[1:15,];e$x[2]<-NA;e$g<-factor(e$g,levels=rev(levels(d$g)))
    p<-predict(f,newdata=e,na.action=na.exclude);expect_length(p,15L);expect_true(is.na(p[2]))
    expect_equal(p[-2],predict(f,newdata=e[-2,]))
  })
  test_that('NN bandwidth types retain the same valid-row estimator', {
    options(na.action='na.omit');e<-d[1:15,];e$x[2]<-NA
    for(type in c('generalized_nn','adaptive_nn')){
      f<-make(type=type);p<-predict(f,newdata=e,na.action=na.exclude)
      expect_length(p,15L);expect_equal(p[-2],predict(f,newdata=e[-2,]))
    }
  })
  test_that('time-series alignment is preserved before evaluation policy', {
    options(na.action='na.omit');yt<-ts(cumsum(rnorm(100)),start=1980,frequency=4);xt<-ts(rnorm(100),start=1980,frequency=4)
    dd<-data.frame(yt=I(yt),xt=I(xt));b<-npindexbw(diff(yt)~lag(yt,-1)+xt,data=dd,bws=c(1,.3,.9),bandwidth.compute=FALSE)
    ye<-ts(cumsum(rnorm(30)),start=2010,frequency=4);xe<-ts(rnorm(30),start=2010,frequency=4);xe[5]<-NA
    ev<-data.frame(yt=I(ye),xt=I(xe));z<-npindex(b,newdata=ev,se=FALSE,na.action=na.exclude)
    options(na.action='na.exclude');ref<-npindex(b,newdata=ev,se=FALSE);same(z,ref)
    expect_length(fitted(z),29L);expect_equal(sum(is.na(fitted(z))),1L)
  })
})
