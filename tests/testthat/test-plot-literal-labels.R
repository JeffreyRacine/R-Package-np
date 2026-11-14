test_that("graphical label transport retains values and other argument semantics", {
  labels <- list("plain", expression(alpha^2), quote(alpha^2), as.name("alpha"), NULL)
  for (label in labels) for (name in c("main","sub","xlab","ylab","zlab")) {
    args <- list(); args[name] <- list(label)
    out <- do.call(function(...)list(...), .np_plot_literal_labels(args))
    expect_identical(out, args)
  }
  for(label in list(quote(alpha^2),as.name("alpha"))) {
    out<-do.call(function(...)list(...),.np_plot_literal_labels(list(main=label),expression.labels=TRUE))
    expect_identical(out$main,as.expression(label))
  }
  count <- 0L; args <- list(main=quote({stop("must not evaluate")}),
                            col=quote({count <- count+1L; "red"}))
  out <- do.call(function(...)list(...), .np_plot_literal_labels(args))
  expect_identical(out$main,args$main);expect_identical(out$col,"red");expect_identical(count,1L)
  expect_identical(.np_plot_conmode_panel_title("Title","x"),"Title: x")
  expect_identical(.np_plot_conmode_panel_title(quote(alpha^2),"x"),quote(alpha^2))
})

test_that("conmode honors literal labels in fresh panels and saved replay", {
  if (identical(environmentName(environment(.np_plot_literal_labels)),"npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE,plot.par.mfrow=FALSE)
  withr::local_seed(75)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)),
                  g=factor(rep(c("lo","hi"),15L)))
  fits <- list(npconmode(g~x,data=d,bws=c(.2,.3),bandwidth.compute=FALSE,probabilities=TRUE),
               npconmode(g~x+z,data=d,bws=c(.2,.3,.4),bandwidth.compute=FALSE,probabilities=TRUE))
  seen <- list(); original <- graphics::title
  testthat::local_mocked_bindings(title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...) {
    seen[[length(seen)+1L]] <<- list(main=main,xlab=xlab,ylab=ylab)
    original(main=main,sub=sub,xlab=xlab,ylab=ylab,...)
  }, .package="graphics")
  for(fit in fits)for(label in list("Title",expression(alpha^2),quote(alpha^2),as.name("alpha"))) {
    seed <- .Random.seed; seen <- list(); count <- 0L
    out <- plot(fit,view="fixed",perspective=FALSE,neval=4L,output="plot-data",
                 main={count<-count+1L;label},xlab=label,ylab=label)
    fresh <- seen
    expect_identical(count,1L)
    expect_length(fresh,length(out$panels))
    for(i in seq_along(fresh)) {
      expect_identical(fresh[[i]]$xlab,label);expect_identical(fresh[[i]]$ylab,label)
      if(is.language(label))expect_identical(fresh[[i]]$main,label)
    }
    bytes <- serialize(out,NULL);seen <- list();plot(unserialize(bytes))
    expect_identical(seen,fresh);expect_identical(serialize(out,NULL),bytes)
    expect_identical(.Random.seed,seed)
  }
})

test_that("shared plot dispatch and replay retain call labels", {
  if (identical(environmentName(environment(.np_plot_literal_labels)),"npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE,plot.par.mfrow=FALSE)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  d<-data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  d$y<-sin(3*d$x)+.2*d$z
  fits<-list(npreg(y~x,data=d,bws=.3,bandwidth.compute=FALSE),
    npudens(~x,data=d,bws=.3,bandwidth.compute=FALSE),
    npudist(~x,data=d,bws=.3,bandwidth.compute=FALSE),
    npplreg(y~x|z,data=d,bws=matrix(c(.4,.4),2,1),bandwidth.compute=FALSE),
    npscoef(y~x|z,data=d,bws=.4,bandwidth.compute=FALSE))
  label<-quote(alpha^2)
  for(fit in fits) {
    plain<-plot(fit,neval=4L,perspective=FALSE,output="data")
    expect_no_error(out<-plot(fit,neval=4L,perspective=FALSE,output="plot-data",main=label,xlab=label,ylab=label))
    expect_no_error(plot(out,main=label,xlab=label,ylab=label))
    expect_identical(lapply(out$panels,`[[`,"values"),lapply(plain$panels,`[[`,"values"))
    expect_identical(lapply(out$panels,`[[`,"intervals"),lapply(plain$panels,`[[`,"intervals"))
  }
})
