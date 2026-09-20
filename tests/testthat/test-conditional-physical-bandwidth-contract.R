test_that("conditional operators preserve physical bandwidth representation", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920322)
  n <- 40L
  x <- data.frame(x=runif(n),g=factor(rep(letters[1:2],length.out=n)))
  y <- data.frame(y=sin(4*x$x)+rnorm(n,sd=.4))
  for (family in c("npcdens","npcdist")) for (engine in c("lc","ll","lp")) {
    args <- list(xdat=x,ydat=y,bws=c(1,1,.4),bwscaling=TRUE,
                 bandwidth.compute=FALSE,regtype=engine)
    if (engine=="lp") args$degree <- 2L
    scaled <- do.call(get(paste0(family,"bw")),args)
    args$bwscaling <- FALSE
    args$bws <- c(scaled$bandwidth$y,scaled$bandwidth$x)
    physical <- do.call(get(paste0(family,"bw")),args)
    for (index in list(seq_len(n),c(1:25,1:15))) {
      a <- do.call(get(family),list(bws=scaled,txdat=x[index,],tydat=y[index,,drop=FALSE],
                                   exdat=x[1:5,],eydat=y[1:5,,drop=FALSE],se=TRUE,gradients=TRUE))
      b <- do.call(get(family),list(bws=physical,txdat=x[index,],tydat=y[index,,drop=FALSE],
                                   exdat=x[1:5,],eydat=y[1:5,,drop=FALSE],se=TRUE,gradients=TRUE))
      expect_equal(fitted(a),fitted(b),tolerance=1e-10)
      expect_equal(se(a),se(b),tolerance=1e-10)
      expect_equal(gradients(a),gradients(b),tolerance=1e-10)
      aa <- .np_conditional_eval_selected(scaled,x[index,],y[index,,drop=FALSE],
        x[1:5,],y[1:5,,drop=FALSE],cdf=family=="npcdist",gradients=TRUE)
      bb <- .np_conditional_eval_selected(physical,x[index,],y[index,,drop=FALSE],
        x[1:5,],y[1:5,,drop=FALSE],cdf=family=="npcdist",gradients=TRUE)
      field <- if (family=="npcdist") "condist" else "condens"
      expect_equal(aa[[field]],bb[[field]],tolerance=1e-10)
      expect_equal(aa$congrad,bb$congrad,tolerance=1e-10)
      expect_equal(aa$conderr,bb$conderr,tolerance=1e-10)
    }
    expect_equal(.npcdhat_make_xbw(scaled,x)$bw,physical$xbw,tolerance=1e-14)
    expect_equal(.npcdhat_make_ybw(scaled,y)$bw,physical$ybw,tolerance=1e-14)
    expect_equal(.np_con_make_kbandwidth_x(scaled,x)$bw,physical$xbw,tolerance=1e-14)
    expect_equal(.np_con_make_kbandwidth_xy(scaled,x,y)$bw,
                 c(physical$xbw,physical$ybw),tolerance=1e-14)
    expect_equal(.np_conditional_count_side_state(x,x[1:5,],scaled,"x")$bandwidth,
                 physical$xbw[physical$ixcon],tolerance=1e-14)
  }
})

test_that("bootstrap pilots keep scale factors separate from physical bandwidths", {
  withr::local_options(list(np.messages=FALSE))
  i <- seq_len(36L)
  x <- data.frame(x=sin(i*.73),g=factor(rep(letters[1:3],12L)))
  y <- cos(i*.27)+.2*x$x
  for(family in c("npreg","npcdens","npcdist","npscoef","npplreg")) {
    args <- list(xdat=x,ydat=y,bws=c(1,.4),bwscaling=TRUE,
                 bandwidth.compute=FALSE)
    if(family %in% c("npcdens","npcdist")) args$bws <- c(1,1,.4)
    if(family %in% c("npscoef","npplreg")) {
      args$zdat <- x
      args$xdat <- x["x"]
    }
    if(family=="npplreg") args$bws <- matrix(rep(c(1,.4),each=2L),2L)
    constructor <- get(paste0(family,"bw"))
    scaled <- do.call(constructor,args)
    args$bwscaling <- FALSE
    args$bws <- if(family %in% c("npcdens","npcdist"))
      c(scaled$bandwidth$y,scaled$bandwidth$x) else if(family=="npplreg")
      do.call(rbind,lapply(scaled$bw,function(z) z$bandwidth$x)) else
      unlist(scaled$bandwidth,use.names=FALSE)
    physical <- do.call(constructor,args)
    pilot <- switch(family,npreg=.np_plot_oversmooth_regression_bws,
      npcdens=function(z) .np_plot_oversmooth_conditional_bws(z,FALSE),
      npcdist=function(z) .np_plot_oversmooth_conditional_bws(z,TRUE),
      npscoef=.np_plot_oversmooth_scbandwidth_bws,
      npplreg=.np_plot_oversmooth_plbandwidth_bws)
    a <- pilot(scaled); b <- pilot(physical)
    expect_equal(a$bandwidth,b$bandwidth,tolerance=1e-14,info=family)
    expect_equal(a$sfactor,b$sfactor,tolerance=1e-14,info=family)
    if(family=="npplreg") for(child in a$bw)
      expect_equal(child$bw,child$sfactor$x,tolerance=1e-14)
    else if(family %in% c("npcdens","npcdist"))
      expect_equal(a$xbw,a$sfactor$x,tolerance=1e-14)
    else expect_equal(a$bw,unlist(a$sfactor,use.names=FALSE),tolerance=1e-14)
    fit.args <- list(txdat=args$xdat,tydat=y,se=FALSE)
    if(family %in% c("npscoef","npplreg")) fit.args$tzdat <- x
    if(family=="npscoef") fit.args$iterate <- FALSE
    expect_equal(fitted(do.call(get(family),c(list(bws=a),fit.args))),
                 fitted(do.call(get(family),c(list(bws=b),fit.args))),tolerance=1e-11)
    if(family %in% c("npcdens","npcdist")) {
      boot.args <- list(xdat=x,ydat=data.frame(y=y),exdat=x[1:3,],
        eydat=data.frame(y=rep(.1,3)),cdf=family=="npcdist",
        plot.errors.boot.method="inid",plot.errors.boot.blocklen=NULL,
        plot.errors.boot.num=2L,progress.label=NULL,gradient.index=2L)
      set.seed(920332)
      aa <- do.call(.np_plot_conditional_pilot_boot,c(list(bws=scaled),boot.args))
      seed <- .Random.seed
      set.seed(920332)
      bb <- do.call(.np_plot_conditional_pilot_boot,c(list(bws=physical),boot.args))
      expect_identical(.Random.seed,seed)
      for(field in c("t0","t","center"))
        expect_equal(aa[[field]],bb[[field]],tolerance=1e-11)
    }
  }
})
