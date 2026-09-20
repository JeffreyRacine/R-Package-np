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
