test_that("automatic outcomes cannot block valid single-index predictions", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(731)
  d <- data.frame(x = rnorm(120), z = rnorm(120))
  d$y <- rbinom(120, 1, plogis(d$x + .5*d$z))
  b <- npindexbw(y ~ x + z, data = d, bws = c(1, .5, .4),
                 method = "kleinspady", bandwidth.compute = FALSE)
  f <- npindex(b, se = FALSE)
  e <- d[1:30, ]
  reference <- predict(f, newdata = e[c("x", "z")])
  for (v in list(NA, NA_real_, "unknown", e$y + 1L)) {
    nd <- e; nd$y <- v
    fit <- npindex(b, newdata = nd, se = FALSE)
    expect_identical(fitted(fit), reference)
    expect_identical(predict(f, newdata = nd), reference)
    expect_identical(fit$diagnostics.sample, "unavailable")
    expect_equal(fit$diagnostics.nobs, 0L)
    if (!is.double(v)) expect_error(npindex(b, newdata = nd, y.eval = TRUE, se = FALSE), "response")
  }
  expect_error(npindex(b, exdat = e[c("x", "z")], eydat = e$y+1L, se = FALSE), "response")
  bad <- e; bad$x <- "wrong"
  expect_error(predict(f, newdata = bad))
  nd <- e; nd$y[c(2, 5)] <- NA_real_
  scored <- npindex(b, newdata = nd, se = FALSE)
  expect_identical(fitted(scored), reference)
  expect_identical(scored$diagnostics.nobs, 28L)
  expect_equal(scored$CCR.overall, mean(round(reference[-c(2,5)]) == e$y[-c(2,5)]))
  expect_identical(predict(f, newdata = bad, exdat = e[c("x", "z")]), reference)
  native <- npindexbw(xdat = d[c("x", "z")], ydat = d$y,
    bws = c(1,.5,.4), method = "kleinspady", bandwidth.compute = FALSE)
  native$ynames <- "y"
  nd <- e; nd$y <- NA
  expect_identical(fitted(npindex(native, newdata = nd, se = FALSE)), reference)
  expect_error(npindex(native, newdata = nd, y.eval = TRUE, se = FALSE), "response")
})

test_that("differenced time-series outcomes align without choosing prediction rows", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(2); yts <- ts(cumsum(rnorm(60)), frequency = 4)
  b <- npindexbw(diff(yts) ~ lag(yts,-1) + lag(yts,-2),
    bws = c(1,.2,.8), bandwidth.compute = FALSE)
  yn <- ts(cumsum(rnorm(20)), frequency = 4, start = c(16,1))
  e <- data.frame(yts = yn)
  fit <- npindex(b, newdata = e, se = FALSE)
  optout <- npindex(b, newdata = e, y.eval = FALSE, se = FALSE)
  expect_identical(fitted(fit), fitted(optout))
  expect_length(fitted(fit), 19L)
  expect_identical(fit$diagnostics.nobs, 18L)
  expect_equal(fit$MSE, mean((as.numeric(diff(yn))[-1L] - fitted(fit)[1:18])^2))
  expect_identical(fitted(npindex(b, newdata = e, y.eval = TRUE, se = FALSE)), fitted(fit))
})

test_that("prediction residual scale retains training MSE provenance and legacy units", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(75); d <- data.frame(x=rnorm(120), z=rnorm(120))
  d$y <- sin(d$x+.5*d$z)+rnorm(120,sd=.2)
  b <- npindexbw(y~x+z,data=d,bws=c(1,.5,.4),bandwidth.compute=FALSE)
  f <- npindex(b, se=TRUE); e <- d[1:30, ]; e$y <- e$y+5
  for (nd in list(e, e[c("x","z")])) {
    p <- predict(f,newdata=nd,se.fit=TRUE)
    expect_equal(p$residual.scale,f$MSE)
  }
  external <- npindex(b,newdata=e,se=TRUE)
  expect_gt(external$MSE, f$MSE)
  expect_equal(predict(external,se.fit=TRUE)$residual.scale,f$MSE)
  expect_equal(predict(external,newdata=e,se.fit=TRUE)$residual.scale,f$MSE)
  expect_equal(f$training.MSE,mean((d$y-fitted(f))^2))
  replacement <- d; replacement$y <- .5*d$y
  refit <- npindex(b,data=replacement,se=TRUE)
  expect_equal(predict(f,data=replacement,newdata=e,se.fit=TRUE)$residual.scale,
               refit$MSE)
})

test_that("single-index plots preserve non-syntactic raw design names", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(91); x <- data.frame(rnorm(90),factor(rep(letters[1:3],30)),rnorm(90))
  names(x) <- c("my x", "grp-1", "w 2")
  y <- sin(x[[1]]) + .2*x[[3]]+rnorm(90,sd=.1)
  b <- npindexbw(xdat=x,ydat=y,bws=c(1,.2,.3,.4,.5),bandwidth.compute=FALSE)
  f <- npindex(b,se=FALSE)
  for (obj in list(b,f)) for (grad in c(FALSE,TRUE)) {
    p <- plot(obj,plot.behavior="data",plot.errors.method="none",neval=5,gradients=grad)
    expect_true(is.list(p))
    expect_gt(length(p),0L)
  }
})

test_that("fresh stored semiparametric bandwidths use the common complete sample", {
  old <- options(np.messages=FALSE,np.extendednn=FALSE); on.exit(options(old),add=TRUE)
  set.seed(94); x <- data.frame(x=rnorm(90),z=rnorm(90)); y<-rnorm(90)
  x$x[2]<-NA; x$z[5]<-NA; y[8]<-NA
  b <- npindexbw(xdat=x,ydat=y,bws=c(1,.5,.4),bandwidth.compute=FALSE)
  expect_identical(b$nobs,87L); expect_equal(b$rows.omit,c(2L,5L,8L)); expect_equal(b$nobs.omit,3L)
  p <- npplregbw(xdat=x["x"],zdat=x["z"],ydat=y,bws=matrix(.4,2,1),bandwidth.compute=FALSE)
  expect_identical(p$nobs,87L); expect_equal(p$rows.omit,c(2L,5L,8L))
  expect_true(all(vapply(p$bw,function(b)b$nobs==87L,logical(1))))
  expect_length(fitted(npindex(b,se=FALSE)),87L)
  expect_length(fitted(npplreg(p)),87L)
  keep <- complete.cases(x,y)
  reference <- npplregbw(xdat=x[keep,"x",drop=FALSE],zdat=x[keep,"z",drop=FALSE],
    ydat=y[keep],bws=matrix(.4,2,1),bandwidth.compute=FALSE)
  expect_identical(fitted(npplreg(p)),fitted(npplreg(reference)))
  options(np.extendednn=FALSE)
  expect_error(npindexbw(xdat=x,ydat=y,bws=c(1,.5,88),bwtype="generalized_nn",bandwidth.compute=FALSE), "exceeds n-1")
})

test_that("density and distribution NN refits validate against current complete rows", {
  old <- options(np.messages=FALSE,np.extendednn=FALSE); on.exit(options(old),add=TRUE)
  set.seed(95); x<-data.frame(x=rnorm(90)); y<-data.frame(y=rnorm(90)); xn<-x;xn$x[1:10]<-NA
  for (type in c("generalized_nn","adaptive_nn"))
    for (family in c("npudensbw","npudistbw","npcdensbw","npcdistbw")) {
      fn<-get(family); conditional<-family %in% c("npcdensbw","npcdistbw")
      args<-if(conditional) list(xdat=x,ydat=y,bws=c(85,85)) else list(dat=x,bws=85)
      args$bwtype<-type;args$bandwidth.compute<-FALSE
      b<-do.call(fn,args);args$bws<-b;args$bwtype<-NULL
      if(conditional)args$xdat<-xn else args$dat<-xn
      expect_error(do.call(fn,args),"exceeds n-1")
      options(np.extendednn=TRUE)
      expect_equal(do.call(fn,args)$nobs,80L)
      options(np.extendednn=FALSE)
    }
})


test_that("automatic factor outcomes cannot prevent Ichimura prediction", {
  old <- options(np.messages=FALSE); on.exit(options(old),add=TRUE)
  set.seed(5101); d <- data.frame(x=rnorm(90),z=rnorm(90)); d$y <- rpois(90,2)
  e <- d[1:20,]; e$y <- factor(e$y)
  for (factor.train in c(FALSE,TRUE)) {
    tr <- d
    if (factor.train) {tr$y <- factor(tr$y); e$y <- factor(rep("99",20))}
    b <- npindexbw(y ~ x + z, data=tr, bws=c(1,.5,.8), bandwidth.compute=FALSE)
    f <- npindex(b,se=FALSE)
    expected <- predict(f,newdata=e,y.eval=FALSE)
    expect_identical(predict(f,newdata=e),expected)
    expect_identical(fitted(npindex(b,newdata=e,se=FALSE)),expected)
    expect_error(predict(f,newdata=e,y.eval=TRUE), "factor|level")
    expect_error(npindex(b,exdat=e[c("x","z")],eydat=e$y,se=FALSE), "factor|level")
  }
})


test_that("factor dispatch retains only public design and training state", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(33252)
  d <- data.frame(x = rnorm(90), g = factor(rep(letters[1:3], 30)))
  d$y <- sin(d$x) + as.integer(d$g) / 4
  make <- function() {
    unrelated <- numeric(2e5)
    npindex(txdat = d[c("x", "g")], tydat = d$y,
      bws = c(1, .2, .3, .8), bandwidth.compute = FALSE, se = FALSE)
  }
  fit <- make()
  expect_null(fit$call)
  expect_lt(length(serialize(fit, NULL)), 8 * 2e5)
  expect_null(attr(fit$bws$.np.native.training$xdat, ".np.index.prepared"))
  expect_null(attr(fit$bws$call$xdat, ".np.index.prepared"))
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(predict(restored, newdata = d[1:8, ]),
                   predict(fit, newdata = d[1:8, ]))
})
