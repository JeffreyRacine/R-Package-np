test_that("smooth coefficient wild gradients reuse canonical projections", {
  pkg <- getNamespaceName(environment(npscoef))
  internal <- function(name) getFromNamespace(name, pkg)
  withr::local_options(np.messages = FALSE, np.plot.wild.hat.block.bytes = 8*64*2)
  withr::local_preserve_seed()
  set.seed(189918)
  n <- 64L
  x <- data.frame(x1 = runif(n, -.8, .8), x2 = rnorm(n))
  z <- data.frame(z = runif(n, -.8, .8), u = factor(rep(1:2, n/2)))
  y <- x$x1 * (1+z$z) + x$x2 + rnorm(n, sd=.2)
  idx <- c(3, 14, 26, 40, 57)
  for (rt in c("lc", "ll", "lp")) for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npscoefbw(xdat=x, ydat=y, zdat=z, bws=c(if(type=="fixed") .8 else 45, .25),
      bandwidth.compute=FALSE, regtype=rt, degree=if(rt=="lp") 2L else NULL, bwtype=type)
    # The incumbent wild owner uses explicit evaluation at the training values;
    # GNN fitted-row identity is a different pilot and must not be substituted.
    pilot <- fitted(npscoef(bw, txdat=x, tydat=y, tzdat=z, exdat=x, ezdat=z, iterate=FALSE))
    for (wild in c("rademacher", "mammen")) {
      set.seed(1981)
      u <- matrix(runif(n*7L), n, 7L)
      draws <- if(wild=="rademacher") ifelse(u<=.5, -1, 1) else
        ifelse(u <= (sqrt(5)+1)/(2*sqrt(5)), (1-sqrt(5))/2, (1+sqrt(5))/2)
      rng <- .Random.seed
      expected <- t(vapply(seq_len(7L), function(i)
        npscoef(bw, txdat=x, tydat=pilot+(y-pilot)*draws[,i], tzdat=z,
          exdat=x[idx,,drop=FALSE], ezdat=z[idx,,drop=FALSE], iterate=FALSE)$grad[,2L],
        numeric(length(idx))))
      tracker <- new.env(parent=emptyenv()); tracker$n <- 0L
      fit <- internal(".np_scoef_fit_internal")
      local({
        testthat::local_mocked_bindings(.np_scoef_fit_internal=function(...) {
          tracker$n <- tracker$n+1L
          if(tracker$n>2L) stop("unexpected per-response fit")
          fit(...)
        }, .package=pkg)
        set.seed(1981)
        actual <- internal(".np_wild_boot_from_scoef_exact")(txdat=x, ydat=y, tzdat=z,
          exdat=x[idx,,drop=FALSE], ezdat=z[idx,,drop=FALSE], bws=bw,
          B=7L, wild=wild, target="grad", gradient.index=2L)
        expect_equal(actual$t, expected, tolerance=1e-10)
        expect_identical(.Random.seed, rng)
        expect_identical(tracker$n, 2L)
      })
    }
  }
})

test_that("coefficient projection retains normalized positive ridge and empty rows", {
  pkg <- getNamespaceName(environment(npscoef))
  internal <- function(name) getFromNamespace(name, pkg)
  withr::local_options(np.messages=FALSE)
  withr::local_preserve_seed()
  set.seed(189716)
  n <- 96L
  x <- data.frame(x=runif(n,-.8,.8), z=runif(n,-.8,.8))
  y <- x$x + rnorm(n)
  idx <- c(5,20,50,80)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bw <- npscoefbw(xdat=x, ydat=y, bws=rep(if(type=="fixed") .7 else 60,2),
      bandwidth.compute=FALSE, regtype="lp",degree=c(2L,2L),bwtype=type)
    state <- internal(".np_scoef_fit_internal")(bw,txdat=x,tydat=y,
      exdat=x[idx,,drop=FALSE],iterate=FALSE,.np_coefficient_projection=1L)
    expect_true(any(state$ridge>0))
    H <- internal(".npscoef_coefficient_projection_block")(state,1L,length(idx))
    for(response in list(y,2*y+7,rnorm(n))) {
      ref <- npscoef(bw,txdat=x,tydat=response,exdat=x[idx,,drop=FALSE],iterate=FALSE)$grad[,1L]
      expect_equal(as.vector(H%*%response),ref,tolerance=1e-10)
    }
  }
  z <- data.frame(z=seq(-1,1,length.out=n))
  ez <- data.frame(z=c(0,10))
  bw <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=.5,bandwidth.compute=FALSE,
    regtype="lc",ckertype="epanechnikov")
  fit.args <- list(bws=bw,txdat=x,tydat=y,tzdat=z,exdat=x[1:2,,drop=FALSE],
    ezdat=ez,iterate=FALSE,.np_allow_undefined=TRUE)
  state <- do.call(internal(".np_scoef_fit_internal"),c(fit.args,list(.np_coefficient_projection=1L)))
  H <- internal(".npscoef_coefficient_projection_block")(state,1L,2L)
  expect_true(all(is.na(H[2L,])))
  expect_warning(ref <- do.call(internal(".np_scoef_fit_internal"),fit.args)$grad[,1L],
    "local fit is undefined")
  expect_equal(as.vector(H%*%y),ref,tolerance=1e-10)
})

test_that("categorical coefficient projections retain their moment kernel units", {
  pkg <- getNamespaceName(environment(npscoef))
  internal <- function(name) getFromNamespace(name, pkg)
  withr::local_options(np.messages=FALSE, np.plot.wild.hat.block.bytes=8*60*2)
  withr::local_preserve_seed()
  set.seed(6019)
  n <- 60L
  x <- data.frame(x=runif(n,-1,1),v=rnorm(n))
  y <- x$x + x$v + rnorm(n,sd=.3)
  for (compress in c(TRUE,FALSE)) for (uk in c("liracine","aitchisonaitken"))
    for (ok in c("liracine","wangvanryzin","racineliyan")) {
      withr::local_options(np.categorical.compress=compress)
      z <- data.frame(u=factor(rep(1:3,20)),o=ordered(rep(0:3,15)))
      e <- c(2L,5L,8L,14L)
      bw <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=c(.3,.2),
        bandwidth.compute=FALSE,ukertype=uk,okertype=ok)
      state <- internal(".np_scoef_fit_internal")(bw,txdat=x,tydat=y,tzdat=z,
        exdat=x[e,,drop=FALSE],ezdat=z[e,,drop=FALSE],iterate=FALSE,
        .np_coefficient_projection=1L)
      expect_identical(state$profile.weights,compress)
      if (compress) {
        expect_length(state$profile$train.id,n)
        expect_identical(nrow(state$profile$train.profile.codes),12L)
      } else expect_null(state$profile)
      H <- internal(".npscoef_coefficient_projection_block")(state,1L,length(e))
      for (response in list(y, 2*y+3)) {
        ref <- npscoef(bw,txdat=x,tydat=response,tzdat=z,exdat=x[e,,drop=FALSE],
                       ezdat=z[e,,drop=FALSE],iterate=FALSE)$grad[,1L]
        expect_equal(drop(H%*%response),ref,tolerance=1e-11)
      }
      # Independent fixed-design weighted normal equations, with no hat/helper.
      codes <- data.frame(u=as.integer(z$u),o=as.integer(z$o)-1L)
      D <- cbind(1,as.matrix(x))
      for (i in seq_along(e)) {
        wu <- if(uk=="liracine") ifelse(codes$u==codes$u[e[i]],1,.3) else
          ifelse(codes$u==codes$u[e[i]],.7,.15)
        d <- abs(codes$o-codes$o[e[i]])
        wo <- switch(ok,liracine=.2^d,
          wangvanryzin=ifelse(d==0,.8,.4*.2^d),
          racineliyan=.2^d/vapply(codes$o,function(a)sum(.2^abs(a-0:3)),numeric(1)))
        w <- wu*wo
        ref <- solve(crossprod(D,D*w),crossprod(D,w*y))[2L]
        expect_equal(drop(H[i,,drop=FALSE]%*%y),as.double(ref),tolerance=1e-11)
      }
      pilot <- fitted(npscoef(bw,txdat=x,tydat=y,tzdat=z,
                               exdat=x,ezdat=z,iterate=FALSE))
      set.seed(414)
      draws <- matrix(ifelse(runif(n*5L)<=.5,-1,1),n,5L)
      rng <- .Random.seed
      expected <- t(vapply(seq_len(5L),function(b)
        npscoef(bw,txdat=x,tydat=pilot+(y-pilot)*draws[,b],tzdat=z,
          exdat=x[e,,drop=FALSE],ezdat=z[e,,drop=FALSE],iterate=FALSE)$grad[,1L],
        numeric(length(e))))
      set.seed(414)
      got <- internal(".np_wild_boot_from_scoef_exact")(x,y,z,x[e,,drop=FALSE],
        z[e,,drop=FALSE],bw,B=5L,target="grad",gradient.index=1L)
      expect_equal(got$t,expected,tolerance=1e-11)
      expect_identical(.Random.seed,rng)
    }
})
