a_capture_empty <- function(expr) {
  warning.state <- new.env(hash = FALSE, parent = emptyenv())
  warning.state$messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warning.state$messages <- c(warning.state$messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warning.state$messages[grepl("all computed kernel weights", warning.state$messages)])
}

test_that("LL and parameterized LP publish only undefined external components", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- expand.grid(x1 = seq(-1, 1, length.out = 8), x2 = seq(-1, 1, length.out = 8))
  y <- sin(x$x1) + x$x2^2 + .1*cos(seq_len(64))
  ex <- data.frame(x1 = c(10, 0, -.25, -10), x2 = c(0, 0, .25, 0))
  for (spec in list(list(regtype = "ll"),
                    list(regtype = "lp", degree = c(3L, 3L)),
                    list(regtype = "lp", degree = c(1L, 3L), bernstein.basis = TRUE),
                    list(regtype = "lp", degree = c(2L, 3L), basis = "additive"),
                    list(regtype = "lp", degree = c(1L, 3L), basis = "tensor"))) {
    bw <- do.call(npregbw, c(list(xdat = x, ydat = y, bws = c(.55,.55),
      bandwidth.compute = FALSE, ckertype = "epanechnikov"), spec))
    for (tree in list(FALSE, TRUE, "auto")) {
      options(np.tree = tree)
      got <- a_capture_empty(npreg(bw, txdat=x, tydat=y, exdat=ex, se=TRUE, gradients=TRUE))
      expect_length(got$warnings, 1L)
      expect_match(got$warnings, "2 external evaluation row")
      expect_match(got$warnings, "\\(1, 4\\)")
      for (field in c("mean","merr")) {
        expect_true(all(is.na(got$value[[field]][c(1,4)])))
        expect_true(all(is.finite(got$value[[field]][2:3])))
      }
      for (field in c("grad","gerr")) expect_true(all(is.na(got$value[[field]][c(1,4),])))
      supported <- suppressWarnings(npreg(bw,txdat=x,tydat=y,exdat=ex[2:3,],se=TRUE,gradients=TRUE))
      expect_equal(got$value$mean[2:3],supported$mean,tolerance=1e-12)
      expect_equal(got$value$merr[2:3],supported$merr,tolerance=1e-12)
      expect_equal(got$value$grad[2:3,,drop=FALSE],supported$grad,tolerance=1e-12)
      expect_equal(got$value$gerr[2:3,,drop=FALSE],supported$gerr,tolerance=1e-12)
      empty <- a_capture_empty(npreg(bw,txdat=x,tydat=y,exdat=ex[c(1,4),],se=TRUE,gradients=TRUE))
      expect_true(all(is.na(empty$value$mean)))
      expect_length(empty$warnings,1L)
      expect_null(attr(got$value,".np.empty.rows",exact=TRUE))
      pred <- a_capture_empty(predict(supported,newdata=ex,se.fit=TRUE))
      expect_equal(as.vector(pred$value$fit),got$value$mean,tolerance=1e-12)
      expect_equal(as.vector(pred$value$se.fit),got$value$merr,tolerance=1e-12)
    }
  }
})

test_that("categorical and nearest-neighbor frames preserve valid components", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=64),f=factor(rep(c("a","b"),32),levels=c("a","b","c")))
  y <- sin(x$x)
  ex <- data.frame(x=c(0,0),f=factor(c("a","c"),levels=levels(x$f)))
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    b <- npregbw(xdat=x,ydat=y,bws=c(if(type=="fixed") .4 else 24,0),
      bwtype=type,bandwidth.compute=FALSE,regtype="lp",degree=3L)
    result <- a_capture_empty(npreg(b,txdat=x,tydat=y,exdat=ex,se=TRUE,gradients=TRUE))
    expect_length(result$warnings,1L)
    expect_true(is.finite(result$value$mean[1L]))
    expect_true(is.na(result$value$mean[2L]))
    expect_true(all(is.na(result$value$grad[2L,])))
    expect_error(.npreg_complete(bws=b,txdat=x,tydat=y,exdat=ex,se=TRUE),"LP solve")
  }
  for (ordered in c(FALSE,TRUE)) {
    x$f <- factor(rep("b",64),levels=c("a","b"),ordered=ordered)
    ex <- x[32,,drop=FALSE]
    b <- npregbw(xdat=x,ydat=y,bws=c(.4,0),bandwidth.compute=FALSE,
                  regtype="lp",degree=3L,ckertype="epanechnikov")
    result <- a_capture_empty(npreg(b,txdat=x,tydat=y,exdat=ex,se=TRUE,gradients=TRUE))
    expect_length(result$warnings,1L)
    expect_true(is.finite(result$value$mean))
    expect_true(is.finite(result$value$merr))
    expect_true(is.finite(result$value$grad[1,1]))
    expect_true(is.finite(result$value$gerr[1,1]))
    expect_true(is.na(result$value$grad[1,2]))
    expect_true(is.na(result$value$gerr[1,2]))
    expect_error(npreg(b,txdat=x,tydat=y,gradients=TRUE),"LP solve")
  }
})

test_that("strict internals, nonfinite systems, radii and row metadata stay intentional", {
  old <- options(np.messages=FALSE,np.tree=FALSE,warn=0);on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(-1,1,length.out=64));y<-sin(x$x);ex<-data.frame(x=c(0,10))
  b<-npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE,regtype="lp",degree=3L,ckertype="epanechnikov")
  expect_error(.npreg_complete(bws=b,txdat=x,tydat=y,exdat=ex),"LP solve")
  expect_error(.np_regression_direct(bws=b,txdat=x,tydat=y,exdat=ex),"LP solve")
  y[32]<-Inf
  expect_error(npreg(b,txdat=x,tydat=y,exdat=ex),"LP solve")
  expect_error(npreg(b,txdat=x,tydat=y,exdat=ex,se=TRUE),"LP solve|HC0")
  y<-sin(x$x)
  metadata<-a_capture_empty(npreg(b,txdat=x,tydat=y,exdat=data.frame(x=c(0,NA,10))))
  expect_match(metadata$warnings,"\\(3\\)")
  expect_equal(as.integer(metadata$value$eval.rows.omit),2L)
  options(warn=2)
  expect_error(npreg(b,txdat=x,tydat=y,exdat=ex),"all computed kernel weights")
  options(warn=0)
  expect_true(is.finite(npreg(b,txdat=x,tydat=y,exdat=ex[1,,drop=FALSE])$mean))
  tied<-data.frame(x=c(rep(0,32),seq_len(32)))
  nn<-npregbw(xdat=tied,ydat=y,bws=2,bandwidth.compute=FALSE,regtype="lp",degree=3L,bwtype="generalized_nn")
  expect_error(npreg(nn,txdat=tied,tydat=y,exdat=data.frame(x=0)),"radius|radii")
})

test_that("beta and computed-zero Gaussian rows use facts without a numerical fallback", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  x<-data.frame(x=seq(.05,.95,length.out=64),f=factor(rep(c("a","b"),32),levels=c("a","b","c")))
  y<-sin(x$x);ex<-data.frame(x=c(.5,.5),f=factor(c("a","c"),levels=levels(x$f)))
  b<-npregbw(xdat=x,ydat=y,bws=c(.1,0),bandwidth.compute=FALSE,regtype="lp",degree=3L,
    ckertype="beta",ckerbound="range",ckerlb=0,ckerub=1)
  r<-a_capture_empty(npreg(b,txdat=x,tydat=y,exdat=ex,gradients=TRUE,se=TRUE))
  expect_length(r$warnings,1L)
  expect_true(is.finite(r$value$mean[1]))
  expect_true(is.na(r$value$mean[2]))
  expect_error(.npreg_complete(bws=b,txdat=x,tydat=y,exdat=ex),"LP solve")
  x<-x["x"]
  g<-npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE,regtype="lp",degree=3L)
  r<-a_capture_empty(npreg(g,txdat=x,tydat=y,exdat=data.frame(x=1e4),gradients=TRUE,se=TRUE))
  expect_true(is.na(r$value$mean))
  expect_length(r$warnings,1L)
  expect_false(grepl("underflow|outside.*support|ties",r$warnings))
})

test_that("composite callers publish once and required bootstrap errors remain terminal", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  set.seed(642)
  x<-data.frame(x=runif(64,-1,1));z<-data.frame(z=seq(-1,1,length.out=64))
  y<-sin(z$z)+.4*x$x+.1*cos(seq_len(64))
  ex<-data.frame(x=c(.1,.2));ez<-data.frame(z=c(0,10))
  bi<-suppressWarnings(npindexbw(xdat=z,ydat=y,bws=c(1,.4),bandwidth.compute=FALSE,
    regtype="lp",degree=3L,ckertype="epanechnikov"))
  index<-a_capture_empty(npindex(bi,txdat=z,tydat=y,exdat=ez,se=TRUE,gradients=TRUE))
  expect_length(index$warnings,1L)
  expect_true(is.na(index$value$mean[2L]))
  # A later required replicate must not publish the initial permissive fit.
  warning.state <- new.env(hash = FALSE, parent = emptyenv())
  warning.state$messages <- character()
  err<-withCallingHandlers(tryCatch(npindex(bi,txdat=z,tydat=y,exdat=ez,
    se=TRUE,gradients=TRUE,se.type="bootstrap",B=3L),error=identity),
    warning=function(w){warning.state$messages<-c(warning.state$messages,conditionMessage(w));invokeRestart("muffleWarning")})
  expect_s3_class(err,"error")
  expect_false(any(grepl("all computed kernel weights",warning.state$messages)))
  bp<-npplregbw(xdat=x,zdat=z,ydat=y,bws=matrix(.4,2,1),bandwidth.compute=FALSE,
    regtype="lp",degree=3L,ckertype="epanechnikov")
  p<-a_capture_empty(npplreg(bp,txdat=x,tydat=y,tzdat=z,exdat=ex,ezdat=ez))
  expect_length(p$warnings,1L)
  expect_true(is.na(p$value$mean[2L]))
  plotted<-a_capture_empty(.np_plot_plreg_asymptotic_fit(bp,x,y,z,ex,ez))
  expect_length(plotted$warnings,1L)
  expect_true(is.na(plotted$value$merr[2L]))
  bl<-nplsqregbw(xdat=z,ydat=y,bws=.4,bandwidth.compute=FALSE,scale=rep(1,64),
    delta=.5,tau=c(.25,.75),regtype="lp",degree=3L,ckertype="epanechnikov")
  q<-a_capture_empty(nplsqreg(bl,txdat=z,tydat=y,exdat=ez,se=TRUE,gradients=TRUE))
  expect_length(q$warnings,1L)
  expect_true(all(is.na(q$value$quantile[2L,])))
  expect_null(attr(q$value,".np.empty.rows",exact=TRUE))
  expect_true(all(vapply(q$value$tau.fits,function(f)is.null(attr(f,".np.empty.rows",exact=TRUE)),logical(1))))
})

test_that("formula frames retain meaningful original row identifiers", {
  old<-options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old),add=TRUE)
  dat<-data.frame(x=seq(-1,1,length.out=64));dat$y<-sin(dat$x)
  b<-npregbw(y~x,data=dat,bws=.3,bandwidth.compute=FALSE,regtype="lp",degree=3L,ckertype="epanechnikov")
  r<-a_capture_empty(npreg(b,newdata=data.frame(x=c(0,NA,10)),na.action=na.exclude))
  expect_length(r$warnings,1L)
  expect_match(r$warnings,"\\(3\\)")
  expect_equal(length(r$value$mean),3L)
  expect_true(all(is.na(r$value$mean[2:3])))
})
