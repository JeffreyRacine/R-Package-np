test_that("derivative selectors preserve names and explicit replay arguments", {
  x <- data.frame(x = seq(-1, 1, length.out = 18), z = sin(seq_len(18)))
  y <- 2*x$x - 3*x$z
  bw <- npregbw(xdat = x, ydat = y, bws = c(.7, .8),
                regtype = "ll", bandwidth.compute = FALSE)
  e <- x[3:5, ]
  h <- npreghat(bw, txdat = x, exdat = e)
  expect_equal(as.numeric(predict(h, deriv = c(z = 1), y = y, output = "apply")),
               rep(-3, 3), tolerance = 1e-10)
  expect_equal(as.numeric(npreghat(bw, txdat = x, exdat = e,
                                  s = c(z = 1), y = y, output = "apply")),
               rep(-3, 3), tolerance = 1e-10)
  for (s in list(.9, -1, Inf, NA_real_, c(bogus = 1), c(x = 1, x = 0)))
    expect_error(npreghat(bw, txdat = x, s = s))
  expect_equal(as.numeric(predict(h, exdat = x[8:9, ], y = y, output = "apply")),
               y[8:9], tolerance = 1e-10)
  expect_equal(as.numeric(npreg(bw, tydat = 2*y)$mean),
               as.numeric(npreg(bw, txdat = x, tydat = 2*y)$mean))
  expect_equal(as.numeric(npreghat(bw, y = 2*y, output = "apply")),
               as.numeric(npreghat(bw, txdat = x, y = 2*y, output = "apply")))
  cb <- npcdensbw(xdat = x, ydat = data.frame(y), bws = c(.7,.8,.6),
                  regtype = "ll", bandwidth.compute = FALSE)
  expect_equal(as.numeric(npcdenshat(cb, txdat = x, tydat = data.frame(y),
    exdat = e, eydat = data.frame(y = y[3:5]), deriv = c(z = 1))),
    as.numeric(npcdenshat(cb, txdat = x, tydat = data.frame(y),
    exdat = e, eydat = data.frame(y = y[3:5]), s = c(0, 1))))
})

test_that("hat replay retains supported explicit ridge and matrix RHS shape", {
  x <- data.frame(x = seq(-1, 1, length.out = 18), z = sin(seq_len(18)))
  y <- x$x^2 + x$z^2
  bw <- npregbw(xdat = x, ydat = y, bws = c(.7,.8), regtype = "lp",
                degree = c(2,2), bandwidth.compute = FALSE)
  h <- npreghat(bw, txdat = x, exdat = x[3:5, ], s = c(1,1), ridge = .1)
  expect_equal(as.numeric(predict(h, newdata = x[3:5, ])), as.numeric(h))
  expect_equal(attr(predict(h, newdata = x[3:5, ]), "ridge.used"),
               attr(h, "ridge.used"))
  for (reg in c("lc", "ll", "lp")) {
    args <- list(xdat = x, ydat = y, bws = c(.7,1,.5), regtype = reg,
                 bandwidth.compute = FALSE)
    if (reg == "lp") args$degree <- 2
    b <- do.call(npindexbw, args)
    Y <- cbind(y, 2*y + 1)
    for (s in 0:1) {
      a <- npindexhat(b, txdat = x, exdat = x[3,,drop=FALSE],
                       y = Y, s = s, output = "apply")
      expect_identical(dim(a), c(1L,2L))
      for (j in 1:2) expect_equal(as.numeric(a[,j]), as.numeric(npindexhat(
        b, txdat = x, exdat = x[3,,drop=FALSE], y = Y[,j], s = s, output = "apply")))
    }
    expect_error(npindexhat(b, txdat = x, s = .9))
  }
})

test_that("copula replay retains inversion controls, order and sample cardinality", {
  x <- data.frame(x = seq(-1,1,length.out=20), z = sin(seq_len(20)))
  for (target in c("density", "distribution")) {
    bw <- if (target == "density") npudensbw(dat=x,bws=c(.5,.6),bandwidth.compute=FALSE) else
      npudistbw(dat=x,bws=c(.5,.6),bandwidth.compute=FALSE)
    u <- matrix(c(.25,.5,.5,.75), ncol=2)
    fit <- npcopula(bw, data=x, u=u, n.quasi.inv=20, er.quasi.inv=.2, se=TRUE)
    replay <- predict(fit, newdata=u, se.fit=TRUE)
    expect_equal(as.numeric(replay$fit), as.numeric(fitted(fit)), tolerance=1e-12)
    expect_equal(as.numeric(replay$se.fit), as.numeric(se(fit)), tolerance=1e-12)
    swapped <- npcopula(bw, data=x[2:1], u=u, n.quasi.inv=20, er.quasi.inv=.2, se=TRUE)
    expect_equal(as.numeric(fitted(swapped)), as.numeric(fitted(fit)))
    for (n in c(10L,15L,25L)) {
      d <- x[rep(seq_len(20),length.out=n), ]
      a <- npcopula(bw, data=d, se=TRUE)
      expect_identical(a$ntrain, n)
      expect_equal(length(fitted(a)), n)
      expect_equal(nrow(a$eval), n)
    }
  }
})

test_that("conditional mode replay retains requested probability projection", {
  x <- data.frame(x=seq(-1,1,length.out=20))
  y <- data.frame(y=factor(rep(0:1,each=10)))
  bw <- npcdensbw(xdat=x,ydat=y,bws=c(.4,.1),regtype="ll",bandwidth.compute=FALSE)
  e <- data.frame(x=c(-1.2,1.2))
  for (proper in c(FALSE, TRUE)) {
    fit <- npconmode(bw,txdat=x,tydat=y,exdat=e,proper=proper,probabilities=TRUE)
    a <- predict(fit,newdata=e,type="prob")
    b <- predict(fit,newdata=e,type="prob",proper=proper)
    expect_equal(a,b)
  }
})

test_that("explicit replacement responses override only retained defaults", {
  set.seed(631)
  d <- data.frame(x=runif(24),w=rnorm(24),z=runif(24),y=rnorm(24))
  for (family in c("npreg","npcdens","npcdist","npqreg","npindex","npplreg","npscoef","npreghat")) {
    conditional <- family %in% c("npcdens","npcdist","npqreg")
    semi <- family %in% c("npplreg","npscoef")
    root <- switch(family,npqreg="npcdist",npreghat="npreg",family)
    ctor <- get(paste0(root,"bw"),mode="function")
    fun <- get(family,mode="function")
    bandwidth <- if (conditional) rep(.5,3) else switch(family,
      npindex=c(.5,1,.5),npplreg=matrix(.5,3,1),npscoef=.5,c(.5,.5))
    for (missing.response in c(FALSE, TRUE)) {
      args <- c(list(xdat=d[c("x","w")],ydat=if(conditional)d["y"]else d$y),
                if(semi)list(zdat=d["z"])else list())
      b <- do.call(ctor,c(args,list(bws=bandwidth,bandwidth.compute=FALSE)))
      replacement <- 2*d$y+1
      if (missing.response) replacement[3] <- NA_real_
      user <- if (family=="npreghat")list(y=replacement,output="apply") else
        list(tydat=if(conditional)data.frame(y=replacement)else replacement,se=FALSE)
      a <- do.call(fun,c(list(bws=b),user))
      ref <- do.call(fun,c(list(bws=b,txdat=d[c("x","w")]),
        if(semi)list(tzdat=d["z"])else list(),user))
      value <- function(v) if(is.numeric(v))as.numeric(v)else as.numeric(fitted(v))
      expect_equal(value(a),value(ref),tolerance=1e-10,info=paste(family,missing.response))
    }
  }
})
