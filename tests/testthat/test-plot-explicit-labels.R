test_that("PL and SC respect X and Z labels in fresh and saved panels", {
  skip_if_not(spawn_mpi_slaves(), "MPI workers unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages=FALSE,plot.par.mfrow=FALSE)
  on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL);on.exit(grDevices::dev.off(),add=TRUE)
  frames <- list(); original.title <- graphics::title
  testthat::local_mocked_bindings(title=function(main=NULL,sub=NULL,xlab=NULL,ylab=NULL,...) {
    frames[[length(frames)+1L]] <<- list(xlab=xlab,ylab=ylab,main=main)
    original.title(main=main,sub=sub,xlab=xlab,ylab=ylab,...)
  },.package="graphics")
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  d$y <- sin(3*d$x)+.2*d$z
  fits <- list(npplreg(y~x|z,data=d,bws=matrix(c(.4,.4),2L,1L),bandwidth.compute=FALSE),
    npscoef(y~x|z,data=d,bws=.4,bandwidth.compute=FALSE))
  for(fit in fits) for(common in c(TRUE, FALSE)) for(labels in list(list(x="X label",z="Z label"),
    list(x="",z=""),list(x=expression(X^2),z=expression(Z^2)))) {
    frames <- list()
    out <- plot(fit,neval=4L,perspective=FALSE,output="plot-data",common.scale=common,
      xlab=labels$x,zlab=labels$z,ylab="response",main="title")
    fresh <- frames;copy <- serialize(out,NULL)
    expected <- lapply(out$panels,function(p) {
      column <- if(!is.null(p$display$columns)) as.integer(p$display$columns[1L]) else match(p$display$coordinates[1L],names(p$grid))
      list(xlab=if(p$variables$role[column]=="Z")labels$z else labels$x,
        ylab="response",main="title")
    })
    expect_identical(unname(fresh),unname(expected))
    frames <- list();plot(out)
    expect_identical(unname(frames),unname(expected))
    expect_identical(serialize(out,NULL),copy)
    frames <- list(); plot(out,xlab="generic")
    generic <- lapply(expected,function(v) {v$xlab <- "generic";v})
    expect_identical(unname(frames),unname(generic))
    frames <- list(); plot(out,xlab="generic",zlab="specific")
    specific <- generic
    for(j in seq_along(out$panels)) {
      p <- out$panels[[j]]
      column <- if(!is.null(p$display$columns)) as.integer(p$display$columns[1L]) else match(p$display$coordinates[1L],names(p$grid))
      if(p$variables$role[column]=="Z") specific[[j]]$xlab <- "specific"
    }
    expect_identical(unname(frames),unname(specific))
    frames <- list(); plot(out,zlab=NULL)
    cleared <- expected
    for(j in seq_along(out$panels)) {
      p <- out$panels[[j]]
      column <- if(!is.null(p$display$columns)) as.integer(p$display$columns[1L]) else match(p$display$coordinates[1L],names(p$grid))
      if(p$variables$role[column]=="Z") cleared[[j]]$xlab <- p$variables$label[column]
    }
    expect_identical(unname(frames),unname(cleared))
    expect_identical(serialize(out,NULL),copy)
  }
})
