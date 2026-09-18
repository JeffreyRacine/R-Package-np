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
