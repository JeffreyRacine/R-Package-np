test_that("literal labels reach conditional, quantile, categorical and IV renderers", {
  if (identical(environmentName(environment(.np_plot_literal_labels)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE,plot.par.mfrow=FALSE,rgl.useNULL=TRUE)
  withr::local_seed(751)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  d<-data.frame(x=seq(.05,.95,length.out=30),z=cos(seq(.1,5,length.out=30)),
                g=factor(rep(c("lo","hi"),15)))
  d$y<-sin(3*d$x)+.2*d$z
  fits<-list(
    density=npcdens(y~x,data=d,bws=c(.3,.3),bandwidth.compute=FALSE),
    distribution=npcdist(y~x,data=d,bws=c(.3,.3),bandwidth.compute=FALSE),
    quantile=npqreg(y~x,data=d,bws=c(.3,.3),bandwidth.compute=FALSE,tau=c(.25,.5)),
    lsq=nplsqreg(nplsqregbw(xdat=d["x"],ydat=d$y,bws=.3,delta=.1,
       scale=rep(1,nrow(d)),bandwidth.compute=FALSE),se=FALSE),
    categorical=npreg(y~g,data=d,bws=.2,bandwidth.compute=FALSE),
    iv=structure(list(z=data.frame(z=c(-2,0,2)),y=c(-10,0,10),
       phi=c(-1.5,0,1.5),zeval=data.frame(z=c(-.5,.5)),phi.eval=c(-1,1)),class="npregiv"))
  label<-quote(alpha^2)
  seen<-list();original<-graphics::title
  local_mocked_bindings(title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...) {
    seen[[length(seen)+1L]]<<-list(main=main,xlab=xlab,ylab=ylab)
    original(main=main,sub=sub,xlab=xlab,ylab=ylab,...)
  },.package="graphics")
  for(name in names(fits)) {
    args<-if(name=="iv")list(data.overlay=TRUE) else list(neval=4,perspective=FALSE)
    plain<-do.call(plot,c(list(fits[[name]],output="data"),args))
    seed<-.Random.seed;seen<-list()
    out<-do.call(plot,.np_plot_literal_labels(c(list(fits[[name]],output="plot-data",main=label,xlab=label,ylab=label),args)))
    expect_true(length(seen)>0,info=name)
    expect_true(any(vapply(seen,function(x)identical(x$main,label),logical(1))),info=name)
    expect_identical(lapply(out$panels,`[[`,"values"),lapply(plain$panels,`[[`,"values"))
    expect_identical(lapply(out$panels,`[[`,"intervals"),lapply(plain$panels,`[[`,"intervals"))
    copy<-serialize(out,NULL);expect_no_error(plot(unserialize(copy)))
    expect_identical(serialize(out,NULL),copy);expect_identical(.Random.seed,seed)
  }
})

test_that("surface handoffs carry literal labels to base and RGL owners", {
  if (identical(environmentName(environment(.np_plot_literal_labels)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  skip_if_not_installed("rgl")
  withr::local_options(np.messages=FALSE,rgl.useNULL=TRUE)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  on.exit(rgl::close3d(),add=TRUE)
  d<-data.frame(x=seq(.05,.95,length.out=30),z=cos(seq(.1,5,length.out=30)));d$y<-sin(3*d$x)+.2*d$z
  fit<-npreg(y~x+z,data=d,bws=c(.3,.4),bandwidth.compute=FALSE)
  label<-quote(alpha^2);seen<-list()
  original.base<-graphics::persp;original.rgl<-rgl::persp3d
  local_mocked_bindings(persp=function(...){a<-list(...);seen[[length(seen)+1L]]<<-a[c("main","xlab","ylab","zlab")];do.call(original.base,.np_plot_literal_labels(a))},.package="graphics")
  local_mocked_bindings(persp3d=function(...){a<-list(...);seen[[length(seen)+1L]]<<-a[c("main","xlab","ylab","zlab")];do.call(original.rgl,.np_plot_literal_labels(a))},.package="rgl")
  for(renderer in c("base","rgl")) {
    seen<-list()
    expect_no_error(out<-plot(fit,neval=4,renderer=renderer,main=label,xlab=label,ylab=label,zlab=label,output="plot-data"))
    expect_true(length(seen)>0)
    expected<-if(renderer=="rgl")as.expression(label) else label
    expect_true(all(vapply(unlist(seen,recursive=FALSE),function(x)is.null(x)||identical(x,expected),logical(1))))
    expect_true(any(vapply(seen,function(x)identical(x$main,expected),logical(1))))
    expect_no_error(plot(out,renderer=renderer))
  }
})

test_that("adjacent graphical owners retain language labels", {
  if (identical(environmentName(environment(.np_plot_literal_labels)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE,plot.par.mfrow=FALSE)
  withr::local_seed(752)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  d<-data.frame(x=seq(.05,.95,length.out=30),z=cos(seq(.1,5,length.out=30)),
                g=factor(rep(c("lo","hi"),15)))
  d$y<-sin(3*d$x)+.2*d$z
  label<-quote(alpha^2);symbol<-as.name("beta");seen<-list()
  original<-graphics::title
  local_mocked_bindings(title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...) {
    seen[[length(seen)+1L]]<<-list(main=main,sub=sub,xlab=xlab,ylab=ylab)
    original(main=main,sub=sub,xlab=xlab,ylab=ylab,...)
  },.package="graphics")
  index<-npindex(y~x+z,data=d,bws=c(1,.5,.3),bandwidth.compute=FALSE)
  cm<-npconmode(g~x,data=d,bws=c(.2,.3),bandwidth.compute=FALSE,probabilities=TRUE)
  cop<-npcopula(~x+z,data=d,bws=c(.3,.4),neval=4)
  run<-list(index=function()plot(index,neval=4,perspective=FALSE,xlab=label,main=label),
    conmode_all=function()expect_warning(plot(cm,neval=4,view="fixed",errors="asymptotic",band="all",
      gaussian.draws=99L,xlab=label,main=label), "gaussian.draws=99"),
    factor=function()plotFactor(factor(c("a","b")),c(.2,.4),xlab=label,main=label))
  for(view in c("contour","image","empirical","all"))
    run[[paste0("copula_",view)]]<-local({v<-view;function()plot(cop,view=v,xlab=label,ylab=symbol)})
  for(name in names(run)) {
    seen<-list();seed<-.Random.seed
    expect_no_error(run[[name]]())
    expect_true(any(vapply(seen,function(x)identical(x$xlab,label),logical(1))),info=name)
    expect_identical(.Random.seed,seed,info=name)
  }
  reg<-npreg(y~x,data=d,bws=.3,bandwidth.compute=FALSE)
  seen<-list();out<-plot(reg,neval=4,sub=label,output="plot-data")
  expect_true(any(vapply(seen,function(x)identical(x$sub,label),logical(1))))
  seen<-list();plot(out)
  expect_true(any(vapply(seen,function(x)identical(x$sub,label),logical(1))))
})

test_that("RGL public surface translates symbol labels only at its boundary", {
  if (identical(environmentName(environment(.np_plot_literal_labels)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  skip_if_not_installed("rgl")
  withr::local_options(np.messages=FALSE,rgl.useNULL=TRUE)
  on.exit(rgl::close3d(),add=TRUE)
  d<-data.frame(x=seq(.05,.95,length.out=30),z=cos(seq(.1,5,length.out=30)))
  d$y<-sin(3*d$x)+.2*d$z
  fit<-npreg(y~x+z,data=d,bws=c(.3,.4),bandwidth.compute=FALSE)
  symbol<-as.name("beta");seen<-list();original<-rgl::persp3d
  local_mocked_bindings(persp3d=function(...) {
    a<-list(...);seen[[length(seen)+1L]]<<-a$main
    do.call(original,.np_plot_literal_labels(a))
  },.package="rgl")
  out<-plot(fit,neval=4,renderer="rgl",main=symbol,output="plot-data")
  expect_true(length(seen)>0L)
  expect_true(all(vapply(seen,identical,logical(1),as.expression(symbol))))
  bytes<-serialize(out,NULL);seen<-list();plot(out,renderer="rgl")
  expect_true(length(seen)>0L)
  expect_true(all(vapply(seen,identical,logical(1),as.expression(symbol))))
  expect_identical(serialize(out,NULL),bytes)
})
