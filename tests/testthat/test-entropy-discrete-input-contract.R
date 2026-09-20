entropy_input_payload <- function(x) {
  x$timing.profile <- NULL
  x
}

test_that("entropy plug-in bandwidths belong to their complete declared sample", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  plugin <- function(x) {
    c <- nlevels(x); n <- length(x); p <- as.numeric(table(x))/n
    (c/(c-1)*sum(p*(1-p)))/
      (c^2/(c-1)^2*sum(p*(1-p)) + n*sum(((1-c*p)/(c-1))^2))
  }
  for (ord in c(FALSE, TRUE)) {
    x <- factor(rep(1:2,c(7,5)),levels=1:4,ordered=ord)
    y <- factor(rep(1:3,c(2,4,9)),levels=1:4,ordered=ord)
    for (bx in list(NULL,.2)) for (by in list(NULL,.3)) {
      a <- npunitest(x,y,bw.x=bx,bw.y=by,bootstrap=FALSE)
      expect_equal(a$bw.x, if(is.null(bx)) plugin(x) else bx, tolerance=1e-14)
      expect_equal(a$bw.y, if(is.null(by)) plugin(y) else by, tolerance=1e-14)
    }
    # The same selection owner is used by symmetry, without altering reflection.
    z <- factor(rep(1:3,c(7,5,3)),levels=1:4,ordered=ord)
    expect_equal(npsymtest(z,B=9)$bw,plugin(z),tolerance=1e-14)
  }
  constant <- factor(rep(1,12),levels=1:4)
  expect_error(npunitest(constant,factor(rep(1:4,3)),bootstrap=FALSE),
               "at least two distinct factor levels")
  expect_error(npsymtest(constant,B=9),"at least two distinct factor levels")
})

test_that("discrete entropy statistics and every draw use the common domain", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  for (ord in c(FALSE,TRUE)) {
    for (kernel in if(ord) c("wangvanryzin","liracine") else c("aitchisonaitken","liracine")) {
      x <- factor(rep(1:2,c(7,5)),levels=1:4,ordered=ord)
      y <- factor(rep(1:3,c(2,4,9)),levels=1:4,ordered=ord)
      grid <- factor(1:4,levels=1:4,ordered=ord)
      control <- if(ord) list(okertype=kernel) else list(ukertype=kernel)
      density <- function(z,h) fitted(do.call(npudens,c(list(tdat=z,edat=grid,bws=h),control)))
      metric <- function(a,b,ha=.15,hb=.25) .5*sum((sqrt(density(a,ha))-sqrt(density(b,hb)))^2)
      for (method in c("summation","integration")) {
        set.seed(715); saved <- .Random.seed
        actual <- do.call(npunitest,c(list(data.x=x,data.y=y,bw.x=.15,bw.y=.25,
          method=method,B=9,random.seed=13),control))
        expect_identical(.Random.seed,saved)
        expect_equal(actual$Srho,metric(x,y),tolerance=1e-13)
        swapped <- do.call(npunitest,c(list(data.x=y,data.y=x,bw.x=.25,bw.y=.15,
          method=method,bootstrap=FALSE),control))
        expect_equal(actual$Srho,swapped$Srho,tolerance=1e-13)
        set.seed(13)
        pool <- if(method=="summation") x else
          factor(c(as.character(x),as.character(y)),levels=levels(x),ordered=ord)
        draws <- replicate(9, {
          a <- pool[sample.int(length(pool),length(x),replace=TRUE)]
          b <- pool[sample.int(length(pool),length(y),replace=TRUE)]
          metric(a,b)
        })
        expect_equal(as.numeric(actual$Srho.bootstrap),draws,tolerance=1e-13)
        expect_equal(actual$P,mean(draws>actual$Srho),tolerance=1e-13)
        if(!ord && kernel=="aitchisonaitken") {
          f <- as.numeric(table(x))/length(x); g <- as.numeric(table(y))/length(y)
          p <- (1-.15)*f+.15/(nlevels(x)-1)*(1-f)
          q <- (1-.25)*g+.25/(nlevels(y)-1)*(1-g)
          expect_equal(actual$Srho,.5*sum((sqrt(p)-sqrt(q))^2),tolerance=1e-13)
        }
      }
    }
  }
})

test_that("categorical symmetry retains its reflection on the full declared grid", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  for (ord in c(FALSE,TRUE)) for(bt in c("iid","geom")) {
    z <- factor(rep(1:3,c(9,7,4)),levels=1:4,ordered=ord)
    grid <- factor(1:4,levels=1:4,ordered=ord)
    rotate <- function(z) {
      vals <- as.numeric(as.character(z))
      factor(2*median(sort(unique(vals)))-vals,levels=levels(z),ordered=ord)
    }
    metric <- function(z) {
      a <- fitted(npudens(tdat=z,edat=grid,bws=.2))
      b <- fitted(npudens(tdat=rotate(z),edat=grid,bws=.2))
      .5*sum((sqrt(a)-sqrt(b))^2)
    }
    pool <- factor(c(as.character(z),as.character(rotate(z))),
                   levels=levels(z),ordered=ord)
    set.seed(17)
    block <- if(bt=="iid") 1 else b.star(as.numeric(data.matrix(z)),round=TRUE)[1,1]
    reference <- boot::tsboot(seq_len(length(pool)),function(ii) metric(pool[ii]),
      R=9,n.sim=length(z),l=block,sim=if(bt=="iid") "fixed" else "geom")
    actual <- npsymtest(z,bw=.2,B=9,random.seed=17,boot.method=bt)
    expect_equal(actual$Srho,metric(z),tolerance=1e-13)
    expect_equal(as.numeric(actual$Srho.bootstrap),as.numeric(reference$t),tolerance=1e-13)
  }
})

test_that("numeric equality omits independent missing rows before support checks", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  x <- c(-1,-.3,.1,.5,1); y <- c(-.5,0,.4,1.2)
  for(method in c("integration","summation")) {
    set.seed(93); seed <- .Random.seed
    a <- npunitest(c(NA,x,NaN),c(y,NA),bw.x=.4,bw.y=.6,B=9,method=method)
    b <- npunitest(x,y,bw.x=.4,bw.y=.6,B=9,method=method)
    expect_identical(entropy_input_payload(a),entropy_input_payload(b))
    expect_identical(.Random.seed,seed)
    expect_identical(
      entropy_input_payload(npsymtest(c(NA,x,NaN),bw=.4,B=9,method=method)),
      entropy_input_payload(npsymtest(x,bw=.4,B=9,method=method)))
  }
  expect_error(npunitest(c(NA_real_,NaN),y,bw.x=.4,bw.y=.6),"non-missing")
  expect_error(npunitest(x,numeric(),bw.x=.4,bw.y=.6),"non-missing")
})

test_that("common factor support preserves labels and rejects ambiguous order", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  x <- factor(rep(c("a","b"),c(5,3)),levels=c("a","b"))
  y <- factor(rep(c("b","c"),c(3,6)),levels=c("b","c"))
  a <- npunitest(x,y,bw.x=.2,bw.y=.3,B=9)
  b <- npunitest(factor(x,levels=c("a","b","c")),
                factor(y,levels=c("a","b","c")),bw.x=.2,bw.y=.3,B=9)
  expect_identical(entropy_input_payload(a),entropy_input_payload(b))
  expect_error(npunitest(ordered(x,levels=c("a","b")),
    ordered(y,levels=c("b","c")),bw.x=.2,bw.y=.3),
    "unambiguous common level order")
  nx <- ordered(c(10,20,10),levels=c(10,20))
  ny <- ordered(c(20,30,30),levels=c(20,30))
  expect_equal(npunitest(nx,ny,bw.x=.2,bw.y=.3,bootstrap=FALSE)$Srho,
    npunitest(ordered(nx,levels=c(10,20,30)),ordered(ny,levels=c(10,20,30)),
      bw.x=.2,bw.y=.3,bootstrap=FALSE)$Srho,tolerance=1e-13)
})

test_that("categorical reflection never silently discards an observation", {
  old <- options(np.messages=FALSE); on.exit(options(old))
  for(ord in c(FALSE,TRUE)) {
    z <- factor(rep(c(1,2,4),c(7,5,3)),levels=1:4,ordered=ord)
    set.seed(123); seed <- .Random.seed
    expect_error(npsymtest(z,bw=.2,B=9), "reflection leaves the declared support")
    expect_identical(.Random.seed,seed)
    # Widening the declared numeric domain retains the reflection's zero.
    expanded <- factor(z,levels=0:4,ordered=ord)
    helper <- if (exists(".np_entropy_reflect_factor", inherits=TRUE))
      get(".np_entropy_reflect_factor", inherits=TRUE) else
      getFromNamespace(".np_entropy_reflect_factor", "npRmpi")
    reflected <- helper(expanded)
    expect_false(anyNA(reflected))
    expect_equal(as.numeric(as.character(reflected)),
                 c(rep(3,7),rep(2,5),rep(0,3)))
  }
})
