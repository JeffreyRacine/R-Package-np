test_that("IV saved-data replay applies named frame overrides once", {
  old <- options(np.messages=FALSE, plot.par.mfrow=FALSE)
  on.exit(options(old), add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add=TRUE)
  frames <- list()
  original.window <- graphics::plot.window
  testthat::local_mocked_bindings(plot.window=function(xlim,ylim,...) {
    frames[[length(frames)+1L]] <<- list(xlim=xlim,ylim=ylim)
    original.window(xlim=xlim,ylim=ylim,...)
  }, .package="graphics")
  iv <- structure(list(z=data.frame(z=c(-2,0,2)), y=c(-10,0,10),
    phi=c(-1.5,0,1.5),phi.deriv.1=matrix(c(.25,.5,.75),ncol=1L),
    zeval=data.frame(z=c(-.5,.5)),phi.eval=c(-1,1),
    phi.deriv.eval.1=matrix(c(.4,.6),ncol=1L)),class="npregiv")
  ivd <- structure(list(z=iv$z,y=iv$y,phi=iv$phi,phi.prime=c(.25,.5,.75),
    zeval=iv$zeval,phi.prime.eval=c(.4,.6)),class="npregivderiv")
  for(fit in list(iv,ivd)) for(derivative in c(FALSE,TRUE)) {
    for(overlay in c(FALSE,TRUE)) {
      out <- plot(fit,gradients=derivative,data.overlay=overlay,
        data.rug=TRUE,xlim=c(-1,2),ylim=c(-2,3),output="data")
      saved <- serialize(out,NULL)
      seed <- get0(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
      frames <- list()
      err <- tryCatch({plot(out,xlim=c(-.5,1.5),ylim=c(-1,2),
        xlab="override X",ylab="override Y");NULL},error=identity)
      expect_null(err)
      if(is.null(err)) {
        expect_length(frames,1L)
        expect_identical(frames[[1L]],list(xlim=c(-.5,1.5),ylim=c(-1,2)))
      }
      expect_identical(serialize(out,NULL),saved)
      expect_identical(get0(".Random.seed",envir=.GlobalEnv,inherits=FALSE),seed)
    }
  }
})

test_that("IV rendering frame precedence covers its three owners", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  frames <- list(); labels <- list()
  original.window <- graphics::plot.window
  original.title <- graphics::title
  testthat::local_mocked_bindings(plot.window=function(xlim,ylim,...) {
    frames[[length(frames)+1L]] <<- list(xlim=xlim,ylim=ylim)
    original.window(xlim=xlim,ylim=ylim,...)
  },title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...) {
    labels[[length(labels)+1L]] <<- list(xlab=xlab,ylab=ylab)
    original.title(main=main,sub=sub,xlab=xlab,ylab=ylab,...)
  },.package="graphics")
  for(overlay in c(FALSE,TRUE)) for(continuous in c(FALSE,TRUE)) {
    spec <- list(curve.args=list(x=c(1,2,3),y=c(1,3,2),type="l",
      xlab="saved X",ylab="saved Y",xlim=c(0,4),ylim=c(0,4)),
      frame.args=list(xlim=c(.5,3.5),ylim=c(-1,5),xlab="frame X",ylab="frame Y"),
      overlay=overlay,continuous=continuous,overlay.x=c(1,2,3),
      overlay.y=c(0,2,1),overlay.args=list(col="grey"),rug=FALSE)
    copy <- serialize(spec,NULL);frames <- list();labels <- list()
    err <- tryCatch({.np_iv_plot_render(spec);NULL},error=identity)
    expect_null(err)
    if(is.null(err)) {
      expect_identical(frames[[1L]]$xlim,c(.5,3.5))
      expect_identical(frames[[1L]]$ylim,c(-1,5))
      expect_identical(labels[[1L]]$xlab,"frame X")
      expect_identical(labels[[1L]]$ylab,"frame Y")
    }
    expect_identical(serialize(spec,NULL),copy)
  }
})


test_that("IV replay preserves explicit NULL frame overrides", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  iv <- structure(list(z=data.frame(z=c(-2,0,2)),y=c(-10,0,10),
    phi=c(-1.5,0,1.5),zeval=data.frame(z=c(-.5,.5)),phi.eval=c(-1,1)),class="npregiv")
  frames <- list(); original.window <- graphics::plot.window
  testthat::local_mocked_bindings(plot.window=function(xlim,ylim,...) {
    frames[[length(frames)+1L]] <<- list(xlim=xlim,ylim=ylim)
    original.window(xlim=xlim,ylim=ylim,...)
  },.package="graphics")
  out <- plot(iv,data.overlay=FALSE,xlim=c(-9,9),ylim=c(-8,8),output="data")
  copy <- serialize(out,NULL)
  frames <- list(); plot(out,xlim=NULL,ylim=NULL)
  expect_identical(frames[[1L]],list(xlim=c(-.5,.5),ylim=c(-1,1)))
  expect_identical(serialize(out,NULL),copy)
})

test_that("IV replay overrides match fresh primitive ownership and automatic ranges", {
  old<-options(np.messages=FALSE,plot.par.mfrow=FALSE);on.exit(options(old),add=TRUE)
state<-new.env(parent=emptyenv());state$events<-list();record<-function(kind,args)state$events[[length(state$events)+1L]]<-list(kind=kind,args=args)
orig.window<-graphics::plot.window;orig.title<-graphics::title;orig.xy<-graphics::plot.xy;orig.axis<-graphics::axis
local_mocked_bindings(plot.window=function(xlim,ylim,...){record('window',list(xlim=xlim,ylim=ylim,...));orig.window(xlim=xlim,ylim=ylim,...)},title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...){record('title',list(main=main,sub=sub,xlab=xlab,ylab=ylab,...));orig.title(main=main,sub=sub,xlab=xlab,ylab=ylab,...)},plot.xy=function(xy,type,...){record('xy',list(xy=xy,type=type,...));orig.xy(xy,type,...)},axis=function(side,at=NULL,labels=TRUE,...){record('axis',list(side=side,at=at,labels=labels,...));orig.axis(side=side,at=at,labels=labels,...)},.package='graphics')
grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
observe<-function(fun,args){state$events<-list();warnings<-character();error<-NULL;seed<-get0('.Random.seed',.GlobalEnv,inherits=FALSE);withCallingHandlers(tryCatch(do.call(fun,args),error=function(e)error<<-conditionMessage(e)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')});stopifnot(identical(seed,get0('.Random.seed',.GlobalEnv,inherits=FALSE)));list(draws=state$events,warnings=warnings,error=error)}
fixtures<-list(continuous=structure(list(z=data.frame(z=c(-2,0,2)),y=c(-10,0,10),phi=c(-1.5,0,1.5),phi.deriv.1=matrix(c(.25,.5,.75),ncol=1),zeval=data.frame(z=c(-.5,.5)),phi.eval=c(-1,1),phi.deriv.eval.1=matrix(c(.4,.6),ncol=1)),class='npregiv'))
f<-fixtures[[1]];f$z<-data.frame(z=factor(c('a','b','c')));f$zeval<-data.frame(z=factor(c('a','c'),levels=c('a','b','c')));f$phi.eval<-c(-20,30);fixtures$factor<-f
 ivd<-fixtures$continuous;ivd$phi.prime<-c(.25,.5,.75);ivd$phi.prime.eval<-c(.4,.6);class(ivd)<-'npregivderiv';fixtures$derivative_family<-ivd
 positive<-fixtures$continuous;positive$z[[1]]<-c(1,2,5);positive$zeval[[1]]<-c(2,3);positive$y<-c(1,2,10);positive$phi<-c(1,2,3);positive$phi.eval<-c(1.5,2.5);fixtures$positive<-positive
variants<-list(null=list(xlim=NULL,ylim=NULL),ylim=list(ylim=NULL),limits=list(xlim=c(-1,1),ylim=c(-2,2)),labels=list(xlab='new X',ylab='new Y',main='new'),style=list(col='red',lwd=5,type='b'),points=list(pch=3,cex=1.5,bg='green'),clear=list(col=NULL,lty=NULL,lwd=NULL,pch=NULL,cex=NULL,bg=NULL),none=list())
results<-list();i<-0L;for(n in names(fixtures))for(deriv in c(FALSE,TRUE))for(overlay in c(FALSE,TRUE))for(rug in c(FALSE,TRUE))for(v in names(variants)){

 fit<-fixtures[[n]];base<-list(gradients=deriv,data.overlay=overlay,data.rug=rug,xlim=c(-3,3),ylim=c(-40,40),xlab='saved X',ylab='saved Y',col='blue',lty=2,lwd=3)
 saved<-do.call(plot,c(list(x=fit,output='data'),base));copy<-serialize(saved,NULL);delta<-variants[[v]];fresh<-observe(plot,c(list(x=fit),utils::modifyList(base,delta,keep.null=TRUE)));cand<-observe(plot,c(list(x=saved),delta));stopifnot(identical(copy,serialize(saved,NULL)))
 i<-i+1L;results[[i]]<-list(id=paste(n,deriv,overlay,rug,v),fresh=fresh,after=cand,equal=identical(fresh,cand)); expect_identical(cand,fresh,info=paste(n,deriv,overlay,rug,v))
}

fit<-fixtures$positive
out<-plot(fit,data.overlay=TRUE,data.rug=TRUE,xlim=c(.5,6),ylim=c(.5,12),output="data")
expect_identical(observe(plot,list(x=out,xlim=NULL,ylim=NULL,log="xy")),
                 observe(plot,list(x=fit,data.overlay=TRUE,data.rug=TRUE,log="xy")))
# Stray and duplicated replay names retain their prior modifyList policy.
plain<-plot(fixtures$continuous,data.overlay=FALSE,output="data")
expect_identical(observe(plot,list(plain,NULL,"xy")),observe(plot,list(x=plain)))
expect_identical(observe(plot,c(list(x=plain),list(main="first",main="second"))),
                 observe(plot,list(x=plain,main="first")))
# These forms discriminate log filtering and NULL ownership from the ordinary grid.
for(n in c("continuous","factor","derivative_family"))
  for(deriv in c(FALSE,TRUE)) for(overlay in c(FALSE,TRUE))
    for(argument in c("type","xlab","ylab")) {
      fit<-fixtures[[n]]
      base<-list(gradients=deriv,data.overlay=overlay,data.rug=TRUE,
        type="b",xlab="retained X",ylab="retained Y")
      saved<-do.call(plot,c(list(x=fit,output="data"),base))
      copy<-serialize(saved,NULL);delta<-setNames(list(NULL),argument)
      expect_identical(observe(plot,c(list(x=saved),delta)),
        observe(plot,c(list(x=fit),utils::modifyList(base,delta,keep.null=TRUE))),
        info=paste(n,deriv,overlay,argument,"NULL"))
      expect_identical(serialize(saved,NULL),copy)
    }
mixed<-fixtures$continuous
for(overlay in c(FALSE,TRUE)) {
  saved<-plot(mixed,data.overlay=overlay,data.rug=TRUE,output="data")
  copy<-serialize(saved,NULL)
  expect_identical(observe(plot,list(x=saved,log="y",ylim=NULL)),
    observe(plot,list(x=mixed,data.overlay=overlay,data.rug=TRUE,log="y")),
    info=paste("mixed sign logarithmic range",overlay))
  expect_identical(serialize(saved,NULL),copy)
}
calls<-0L
plot(plain,xlim={calls<-calls+1L;c(-1,1)})
expect_identical(calls,1L)

})
