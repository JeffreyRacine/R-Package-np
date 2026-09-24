for (family in c("density", "regression", "conditional")) {
  test_that(paste(family, "automatic recovery reaches the deleted-sample endpoint"), {
    skip_on_cran()
    skip_if_not_installed("crs", minimum_version = "0.15.46")

    old <- options(np.messages=FALSE, np.tree=FALSE, np.extendednn=FALSE)
    on.exit(options(old), add=TRUE)
    x <- data.frame(x=c(rep(0,22),1:2))
    y <- sin(seq_len(24)/3)+seq_len(24)/90
    yd <- data.frame(y=c(rep(.2,22),.55,.9))
    n <- nrow(x)
    args <- switch(family, density=list(dat=x,bwmethod="cv.ml"),
      regression=list(xdat=x,ydat=y,regtype="lc",bwmethod="cv.ls"),
      conditional=list(xdat=x,ydat=yd,regtype="lc",bwmethod="cv.ml"))
    construct <- function(a) switch(family,
      density=do.call(npudensbw,a), regression=do.call(npregbw,a),
      conditional=do.call(npcdensbw,a))
    raw <- function(b) switch(family,
      density=npudensbw.bandwidth(dat=x,bws=b,eval.only=TRUE,nmulti=1L,
        invalid.penalty="dbmax")$fval,
      regression=.npregbw_eval_only(x,y,b,invalid.penalty="dbmax")$objective,
      conditional=.npcdensbw_eval_only(x,yd,b,invalid.penalty="dbmax")$objective)
    # Independent deleted-query Gaussian weights; no package kernel helper.
    wx <- vapply(seq_len(n), function(i) {
      h <- sort(abs(x$x[-i]-x$x[i]))[n-2L]
      z <- dnorm((x$x-x$x[i])/h)/h; z[i] <- 0; z
    }, numeric(n))
    expected <- switch(family,
      density=sum(log(colSums(wx)/(n-1L))),
      regression=mean((y-colSums(wx*y)/colSums(wx))^2),
      conditional={
        wy <- vapply(seq_len(n), function(i) {
          h <- sort(abs(yd$y[-i]-yd$y[i]))[n-2L]
          z <- dnorm((yd$y-yd$y[i])/h)/h; z[i] <- 0; z
        }, numeric(n))
        sum(log(colSums(wx*wy)/colSums(wx)))
      })
    for (type in c("generalized_nn", "adaptive_nn")) {
      args$bwtype <- type
      b <- construct(c(args,list(bws=rep(n-2L,if(family=="conditional")2L else 1L),
                                 bandwidth.compute=FALSE)))
      expect_true(is.finite(raw(b)) && abs(raw(b)) < .Machine$double.xmax)
      if(type=="generalized_nn")
        expect_equal(as.numeric(raw(b)),expected,tolerance=2e-12)
      b <- construct(c(args,list(bws=rep(n-1L,if(family=="conditional")2L else 1L),
                                 bandwidth.compute=FALSE)))
      expect_identical(abs(as.numeric(raw(b))),.Machine$double.xmax)
      for(solver in c("powell", "mads", "mads+powell")) {
        set.seed(42)
        fitted.bw <- construct(c(args,list(bwsolver=solver,nmulti=1L,itmax=20L,
          powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=30L))))
        expect_equal(as.numeric(if(family=="conditional")
                       c(fitted.bw[["ybw"]],fitted.bw[["xbw"]]) else fitted.bw[["bw"]]),
                     rep(n-2L,if(family=="conditional")2L else 1L),tolerance=0)
        expect_true(is.finite(fitted.bw$fval) &&
                    abs(fitted.bw$fval) < .Machine$double.xmax)
        expect_equal(as.numeric(fitted.bw$fval),as.numeric(raw(fitted.bw)),
                     tolerance=2e-12)
      }
      expect_error(construct(c(args,list(bws=rep(3,if(family=="conditional")2L else 1L),
        bwsolver="mads",nmulti=1L,itmax=20L,nomad.opts=list(MAX_BB_EVAL=1L)))),
        "did not return a raw-valid solution")
    }
  })
}

test_that("explicit extended NN counts bypass ordinary recovery", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version="0.15.46")
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(rep(0,22),1:2));y <- sin(seq_len(24)/3)+seq_len(24)/90
  yd <- data.frame(y=c(rep(.2,22),.55,.9))
  for(type in c("generalized_nn","adaptive_nn")) {
    controls <- list(bwtype=type,bwsolver="mads",nmulti=1L,itmax=20L,
      nomad.opts=list(MAX_BB_EVAL=1L))
    cases <- list(
      do.call(npudensbw,c(list(dat=x,bws=26,bwmethod="cv.ml"),controls)),
      do.call(npregbw,c(list(xdat=x,ydat=y,bws=26,regtype="lc"),controls)),
      do.call(npcdensbw,c(list(xdat=x,ydat=yd,bws=c(26,26),regtype="lc",
                              bwmethod="cv.ml"),controls)))
    for(i in seq_along(cases)) {
      b <- cases[[i]]
      counts <- if(i==3L)c(b[["ybw"]],b[["xbw"]]) else b[["bw"]]
      expect_identical(as.numeric(counts),rep(26,if(i==3L)2L else 1L))
      expect_true(is.finite(b$fval) && abs(b$fval)<.Machine$double.xmax)
      expect_false(any(vapply(b$nomad.restart.results,function(z)isTRUE(z$recovery),TRUE)))
    }
  }
})

test_that("CVAIC recovery retains the full-sample endpoint for both NN types", {
  skip_on_cran()
  skip_if_not_installed("crs", minimum_version = "0.15.46")

  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(rep(0,23),1)); y <- sin(seq_len(24)/3)+seq_len(24)/90
  for(type in c("generalized_nn","adaptive_nn"))
  for(solver in c("powell","mads","mads+powell")) {
    set.seed(42)
    b <- npregbw(xdat=x,ydat=y,bwtype=type,bwmethod="cv.aic",regtype="lc",
      bwsolver=solver,nmulti=1L,itmax=20L,powell.remin=FALSE,
      nomad.opts=list(MAX_BB_EVAL=30L))
    expect_identical(as.numeric(b$bw),23)
    expect_true(is.finite(b$fval) && b$fval < .Machine$double.xmax)
    expect_equal(as.numeric(b$fval),
      as.numeric(.npregbw_eval_only(x,y,b,invalid.penalty="dbmax")$objective),
      tolerance=2e-12)
  }
})
