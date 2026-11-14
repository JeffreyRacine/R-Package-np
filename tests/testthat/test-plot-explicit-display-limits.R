test_that("line display limits reach graphics without changing plot statistics", {
  old <- options(np.messages=FALSE, plot.par.mfrow=FALSE)
  on.exit(options(old), add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add=TRUE)
  received <- list()
  original.window <- graphics::plot.window
  testthat::local_mocked_bindings(plot.window=function(xlim, ylim, ...) {
    received[[length(received)+1L]] <<- list(xlim=xlim,ylim=ylim)
    original.window(xlim=xlim,ylim=ylim,...)
  }, .package="graphics")
  d <- data.frame(x=seq(.05,.95,length.out=30L),
                  z=cos(seq(.1,5,length.out=30L)))
  d$y <- sin(3*d$x)+.2*d$z
  fits <- list(npudens(~x,data=d,bws=.3,bandwidth.compute=FALSE),
               npudist(~x,data=d,bws=.3,bandwidth.compute=FALSE),
               npreg(y~x,data=d,bws=.3,bandwidth.compute=FALSE))
  for(fit in fits) {
    before <- plot(fit,neval=5L,output="data")
    received <- list()
    after <- plot(fit,neval=5L,output="plot-data",xlim=c(.8,.2),ylim=c(-.2,2))
    expect_length(received,1L)
    expect_identical(received[[1L]],list(xlim=c(.8,.2),ylim=c(-.2,2)))
    expect_identical(before$panels[[1L]][c("grid","values","intervals")],
                     after$panels[[1L]][c("grid","values","intervals")])
    received <- list()
    plot(after,xlim=c(.7,.3),ylim=c(-.1,1.5))
    expect_identical(received[[1L]],list(xlim=c(.7,.3),ylim=c(-.1,1.5)))
  }
})

test_that("surface limits reach the frame and foreground and replay", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  received <- list(); original.persp <- getS3method("persp", "default", envir=asNamespace("graphics"))
  testthat::local_mocked_bindings(persp.default=function(x,y,z,xlim,ylim,...) {
    received[[length(received)+1L]] <<- list(xlim=xlim,ylim=ylim)
    original.persp(x=x,y=y,z=z,xlim=xlim,ylim=ylim,...)
  }, .package="graphics")
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  fit <- npudens(~x+z,data=d,bws=c(.3,.4),bandwidth.compute=FALSE)
  out <- plot(fit,neval=4L,view="fixed",output="plot-data",xlim=c(-.1,1.1),ylim=c(-1.1,1.1))
  expect_true(length(received)>=2L)
  for(args in received) expect_identical(args,list(xlim=c(-.1,1.1),ylim=c(-1.1,1.1)))
  received <- list(); plot(out,xlim=c(-.2,1.2),ylim=c(-1.2,1.2))
  expect_true(length(received)>=2L)
  for(args in received) expect_identical(args,list(xlim=c(-.2,1.2),ylim=c(-1.2,1.2)))
})

test_that("limit forwarding does not erase automatic ranges", {
  base <- list(xlim=c(-1,1),ylim=c(0,1),zlim=c(0,2))
  expect_identical(.np_plot_apply_limits(base),base)
  expect_identical(.np_plot_apply_limits(base,xlim=c(.8,.2)),
    list(xlim=c(.8,.2),ylim=c(0,1),zlim=c(0,2)))
})


test_that("copula replay retains an explicitly requested height limit", {

  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  received <- list()
  original.persp <- getS3method("persp", "default", envir=asNamespace("graphics"))
  testthat::local_mocked_bindings(persp.default=function(x,y,z,xlim,ylim,zlim,...) {
    received[[length(received)+1L]] <<- list(xlim=xlim,ylim=ylim,zlim=zlim)
    original.persp(x=x,y=y,z=z,xlim=xlim,ylim=ylim,zlim=zlim,...)
  }, .package="graphics")
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  fit <- npcopula(~x+z,data=d,bws=c(.3,.4))
  out <- plot(fit,view="fixed",output="data",xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.2,1.2))
  expect_identical(out$panels[[1L]]$display$style,
                   list(xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.2,1.2)))
  saved <- serialize(out,NULL)
  plot(out)
  expect_true(length(received)>=2L)
  for(args in received)expect_identical(args,list(xlim=c(-.1,1.1),ylim=c(-.1,1.1),zlim=c(-.2,1.2)))
  expect_identical(serialize(out,NULL),saved)
})


test_that("copula retained height limits follow the selected drawing", {

  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  fit <- npcopula(~x+z,data=d,bws=c(.3,.4),neval=5L)
  for(view in c("contour","image","empirical","all")) {
    out <- plot(fit,view=view,zlim=c(.2,.6),output="data")
    for(p in out$panels) expect_null(p$display$style$zlim)
    expect_warning(plot(out),NA)
  }
  received <- list(); original.persp <- getS3method("persp","default",envir=asNamespace("graphics"))
  testthat::local_mocked_bindings(persp.default=function(x,y,z,...) {
    received[[length(received)+1L]] <<- z
    original.persp(x=x,y=y,z=z,...)
  },.package="graphics")
  out <- plot(fit,view="fixed",zlim=c(.2,.6),output="plot-data")
  fresh <- received[[length(received)]]; original <- serialize(out,NULL)
  received <- list(); expect_warning(plot(out),NA)
  expect_identical(received[[length(received)]],fresh)
  expect_true(all(fresh>=.2 & fresh<=.6))
  expect_identical(serialize(out,NULL),original)
})

test_that("explicit Z-role replay limits take precedence over generic xlim", {

  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=5+cos(seq(.1,5,length.out=30L)))
  d$y <- sin(3*d$x)+.3*d$z
  fits <- list(npplreg(y~x|z,data=d,bws=matrix(c(.4,.4),2,1),bandwidth.compute=FALSE),
               npscoef(y~x|z,data=d,bws=.4,bandwidth.compute=FALSE))
  seen <- list(); original.window <- graphics::plot.window
  testthat::local_mocked_bindings(plot.window=function(xlim,ylim,...) {
    seen[[length(seen)+1L]] <<- xlim
    original.window(xlim=xlim,ylim=ylim,...)
  },.package="graphics")
  for(fit in fits) {
    out <- plot(fit,neval=4L,perspective=FALSE,output="data",xlim=c(-1,2),zlim=c(3,7))
    original <- serialize(out,NULL)
    seen <- list(); plot(out,xlim=c(-2,3),zlim=c(4,6))
    expect_identical(seen,list(c(-2,3),c(4,6)))
    seen <- list(); plot(out,xlim=c(-2,3))
    expect_identical(seen,list(c(-2,3),c(-2,3)))
    seen <- list(); plot(out)
    expect_identical(seen,list(c(-1,2),c(3,7)))
    seen <- list(); plot(out,xlim=c(-2,3),zlim=NULL)
    expect_identical(seen[[1L]],c(-2,3))
    expect_false(identical(seen[[2L]],c(-2,3)))
    expect_identical(serialize(out,NULL),original)
  }
})

test_that("multi-tau factor panels honor the same xlim as continuous panels", {

  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  d <- data.frame(x=seq(.05,.95,length.out=30L),g=factor(rep(letters[1:3],10)))
  d$y <- sin(3*d$x)+.2*as.integer(d$g)
  fit <- npqreg(y~x+g,data=d,bws=c(.4,.3,.2),tau=c(.25,.75))
  seen <- list(); original.window <- graphics::plot.window
  testthat::local_mocked_bindings(plot.window=function(xlim,ylim,...) {
    seen[[length(seen)+1L]] <<- xlim
    original.window(xlim=xlim,ylim=ylim,...)
  },.package="graphics")
  out <- plot(fit,neval=4L,xlim=c(-1,4),output="plot-data")
  expect_length(seen,2L)
  for(lim in seen) expect_identical(lim,c(-1,4))
  expect_s3_class(out,"np_plot_data")
})

test_that("the backdrop grid uses the actual perspective frame limits", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  seen <- NULL
  testthat::local_mocked_bindings(.np_plot_draw_box_grid_persp=function(xlim,ylim,zlim,...) {
    seen <<- list(xlim=xlim,ylim=ylim,zlim=zlim)
    invisible(NULL)
  },.package="np")
  axis <- c(0,.5,1)
  args <- list(x=axis,y=axis,z=outer(axis,axis,`+`),xlim=c(-1,2),ylim=c(-2,3),zlim=c(-1,3))
  .np_plot_render_surface_base_frame(args,c(0,1),c(0,1),c(0,2))
  expect_identical(seen,list(xlim=c(-1,2),ylim=c(-2,3),zlim=c(-1,3)))
})

test_that("rgl rejects unrenderable limits without leaking its device", {
  old <- options(rgl.useNULL=TRUE,browser=function(...) invisible(NULL),viewer=function(...) invisible(NULL))
  on.exit(options(old),add=TRUE)
  skip_if_not_installed("rgl")
  before <- rgl::rgl.dev.list()
  axis <- c(0,.5,1)
  args <- list(x=axis,y=axis,z=outer(axis,axis,`+`),xlab="x",ylab="y",zlab="z",main="",theta=0,phi=20)
  for(name in c("xlim","ylim","zlim")) {
    a <- args; a[[name]] <- c(0,Inf)
    expect_error(do.call(.np_plot_render_surface_rgl,a),paste0(name," must not contain infinite values"))
    expect_identical(rgl::rgl.dev.list(),before)
  }
  a <- args; a$xlim <- c(0,1e300)
  expect_error(do.call(.np_plot_render_surface_rgl,a),"finite coordinate range supported by rgl")
  expect_identical(rgl::rgl.dev.list(),before)
  a$xlim <- c(-1e8,1e8)
  expect_warning(do.call(.np_plot_render_surface_rgl,a),NA)
  expect_identical(rgl::rgl.dev.list(),before)
})


test_that("rgl preserves constant and missing automatic height ranges", {
  skip_if_not_installed("rgl")
  old <- options(np.messages=FALSE,rgl.useNULL=TRUE,
    browser=function(...) invisible(NULL),viewer=function(...) invisible(NULL))
  on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  dev <- rgl::rgl.dev.list()
  axis <- c(0,.5,1)
  args <- list(x=axis,y=axis,z=outer(axis,axis,`+`),xlab="x",ylab="y",zlab="z",main="",theta=0,phi=20)
  limits <- list(NULL,c(0,NA_real_),c(.5,.5),c(2,0),1,c(0,1,2))
  for(limit in limits) {
    a <- args; a["zlim"] <- list(limit)
    expect_no_error(do.call(.np_plot_render_surface_rgl,a))
    expect_identical(rgl::rgl.dev.list(),dev)
  }
  a <- args; a$z[2L,2L] <- NA_real_; a$zlim <- c(NA_real_,NA_real_)
  expect_no_error(do.call(.np_plot_render_surface_rgl,a))
  expect_identical(rgl::rgl.dev.list(),dev)
  d <- data.frame(x=seq(.05,.95,length.out=25L),w=cos(seq(.1,5,length.out=25L)),y=0)
  fit <- npreg(y~x+w,data=d,bws=c(.4,.6),bandwidth.compute=FALSE)
  fit.before <- fit
  data <- plot(fit,neval=4L,view="fixed",output="data",data.overlay=FALSE)
  expect_no_error(drawn <- plot(fit,neval=4L,view="fixed",output="plot-data",renderer="rgl",data.overlay=FALSE))
  expect_identical(data$panels[[1L]][c("grid","values","intervals")],
    drawn$panels[[1L]][c("grid","values","intervals")])
  expect_no_error(plot(data,renderer="rgl"))
  expect_identical(fit,fit.before)
  expect_identical(rgl::rgl.dev.list(),dev)
})
