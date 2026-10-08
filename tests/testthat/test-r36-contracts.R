test_that("smooth coefficient refits retain the z=x convention they used", {
  set.seed(3601); n <- 90L
  x <- data.frame(x=runif(n)); z <- data.frame(z=runif(n))
  y <- sin(4*x$x)+rnorm(n,sd=.2)
  make <- function() {
    xl <- x; zl <- z; yl <- y
    npscoefbw(xdat=xl,ydat=yl,zdat=zl,bws=.3,bandwidth.compute=FALSE)
  }
  b <- make(); X <- data.frame(x=runif(n)); Y <- cos(3*X$x)+rnorm(n,sd=.2)
  W <- cbind(1,X$x)
  oracle <- function(E,Z=E) vapply(seq_len(nrow(E)),function(i) {
    w <- dnorm((X$x-Z[[1]][i])/.3)
    sum(c(1,E$x[i])*solve(crossprod(W,W*w),crossprod(W,w*Y)))
  },0.0)
  E <- X[1:11,,drop=FALSE]
  progress.entry <- new.env(parent = emptyenv()); progress.entry$n <- 0L
  pooled <- npRmpi:::.npRmpi_autodispatch_active()
  if (pooled) {
    original.fit <- npRmpi:::.np_scoef_fit_internal
    testthat::local_mocked_bindings(
      .np_progress_is_interactive = function() TRUE,
      .np_scoef_fit_internal = function(...) {
        progress.entry$n <- progress.entry$n + 1L
        original.fit(...)
      }, .package = "npRmpi")
  }
  for (messages in c(FALSE, TRUE)) for (positional in c(FALSE,TRUE)) {
    withr::local_options(np.messages=messages)
    entered.before <- progress.entry$n
    f <- if(positional) npscoef(b,X,Y) else npscoef(b,txdat=X,tydat=Y)
    if (pooled && messages) expect_gt(progress.entry$n, entered.before)
    expect_equal(fitted(f),oracle(X),tolerance=1e-11)
    expect_equal(fitted(npscoef(b,txdat=X,tydat=Y,tzdat=NULL,newdata=E)),oracle(E),tolerance=1e-11)
    expect_null(f$bws$zdati); expect_null(f$bws$znames); expect_null(f$bws$varnames$z)
    expect_named(f$bws$dati, c("x", "y", "z"))
    expect_named(f$bws$varnames, c("x", "y", "z"))
    expect_null(f$bws$dati$z)
    for (field in c("sfactor", "bandwidth", "sumNum", "klist"))
      expect_named(f$bws[[field]], "x")
    expect_warning(capture.output(summary(f$bws)), NA)
    expect_equal(fitted(npscoef(f$bws)),oracle(X),tolerance=1e-11)
    expect_equal(predict(f,exdat=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(f,exdat=E,ezdat=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(f,newdata=E),oracle(E),tolerance=1e-11)
    expect_equal(predict(unserialize(serialize(f,NULL)),newdata=E),oracle(E),tolerance=1e-11)
  }
  # Explicit separate smoothing data retains its own role.
  f <- npscoef(b,txdat=X,tydat=Y,tzdat=z)
  expect_identical(f$bws$znames,b$znames)
  expect_equal(predict(f,exdat=X,ezdat=z),fitted(f),tolerance=1e-11)
  # A formula with a separate z continues to retain all its formula roles.
  d <- data.frame(x=x$x,z=z$z,y=y)
  bf <- npscoefbw(y~x|z,data=d,bws=.3,bandwidth.compute=FALSE)
  ff <- npscoef(bf)
  expect_equal(predict(ff,newdata=d),fitted(ff),tolerance=1e-11)
})

test_that("forwarded LSQ subsets retain their data mask and evaluate once", {
  set.seed(3620);n<-80L;d<-data.frame(x=runif(n),z=runif(n));d$y<-sin(4*d$x)+rnorm(n,sd=.2)
  b<-npregbw(y~x,data=d,bws=.2,bandwidth.compute=FALSE)
  wrapper<-function(...)nplsqreg(b,data=d,...,delta=.5,bandwidth.compute=FALSE)
  outer<-function(...)wrapper(...)
  for (fun in list(wrapper,outer)) {
    counter<-new.env();counter$n<-0L
    f<-fun(subset={counter$n<-counter$n+1L;x>.4},scale=rep(.3,n))
    keep<-d$x>.4
    oracle<-vapply(d$x[keep],function(e){w<-dnorm((d$x[keep]-e)/.2);sum(w*d$y[keep])/sum(w)},0.0)
    expect_identical(counter$n,1L);expect_equal(fitted(f),oracle,tolerance=1e-12)
  }
})

test_that("fit selectors use the first realization of forwarded controls", {
  # Pooled entry's value binder masks a selector-replay defect. This
  # discriminating test belongs to the local lane, which is run separately.
  skip_if(npRmpi:::.npRmpi_autodispatch_active(), "local selector-replay contract")
  set.seed(3621);n<-80L;x<-data.frame(x=runif(n));z<-data.frame(z=runif(n));y<-sin(4*x$x)+rnorm(n,sd=.2)
  for (tree in c(FALSE,TRUE)) {
    withr::local_options(np.tree=tree)
    counter<-new.env();counter$n<-0L
    f<-npreg(txdat=x,tydat=y,nmulti={counter$n<-counter$n+1L;1L},itmax=3L)
    expect_identical(counter$n,1L)
    counter$n<-0L
    f<-npreg(txdat=x,tydat=y,nomad={counter$n<-counter$n+1L;FALSE},nmulti=1L,itmax=3L)
    expect_identical(counter$n,1L)
  }
  counter<-new.env();counter$n<-0L
  f<-npscoef(txdat=x,tydat=y,tzdat=z,nmulti={counter$n<-counter$n+1L;1L})
  expect_identical(counter$n,1L)
})

test_that("native entry values and response labels survive dispatch", {
  withr::local_options(np.messages=FALSE)
  set.seed(3622);n<-80L;X<-data.frame(x=runif(n),z=rnorm(n));Y<-sin(4*X$x)+.2*X$z+rnorm(n,sd=.2)
  counter<-new.env();counter$n<-0L;s_m<-seq_len(40L);c_m<-'native call label'
  ref<-npregbw(xdat=X,ydat=Y,bws=c(.2,.3),bandwidth.compute=FALSE)
  b<-npregbw(xdat=X,ydat=Y,bws=c(.2,.3),bandwidth.compute=FALSE,subset={counter$n<-counter$n+1L;s_m})
  expect_identical(counter$n,1L);expect_identical(b$bw,ref$bw)
  # Native subset remains an accepted dot, not a formula row-selection request.
  set.seed(3623);a<-npregbw(xdat=X,ydat=Y,nmulti=1L,itmax=3L,subset=s_m)
  set.seed(3623);b<-npregbw(xdat=X,ydat=Y,nmulti=1L,itmax=3L)
  expect_identical(a$bw,b$bw)
  set.seed(3624);a<-npudensbw(dat=X,nmulti=1L,call=c_m)
  set.seed(3624);b<-npudensbw(dat=X,nmulti=1L)
  expect_identical(a$bw,b$bw)
  make<-function(){response<-Y;npindexbw(xdat=X,ydat=response,bws=c(1,.2,.3),bandwidth.compute=FALSE)}
  b<-make();expect_identical(b$ynames,'response')
  evaluation.y <- rev(Y) + .7
  E<-data.frame(X,response=evaluation.y)
  for(obj in list(b,unserialize(serialize(b,NULL)))) {
    f<-npindex(obj,newdata=E,se=FALSE)
    expect_identical(f$diagnostics.sample,'evaluation')
    expect_equal(f$MSE,mean((evaluation.y-fitted(f))^2),tolerance=1e-13)
  }
})

# The default lane adds two small ownership witnesses, not a kernel/search
# matrix. Existing R36 data-only subset controls above cover the adjacent route.
test_that("forwarded subsets retain caller-local bindings", {
  set.seed(3707)
  d <- data.frame(x = runif(90), grp = rep(1:3, c(45, 30, 15)))
  d$y <- sin(4*d$x) + rnorm(90, sd = .2)
  i <- 3L
  wrapper <- function(...) npreg(y ~ x, data = d, bws = .2,
                                 bandwidth.compute = FALSE, ...)
  expect_identical(vapply(1:2, function(i)
    length(fitted(wrapper(subset = d$grp == i))), 0L), c(45L, 30L))
  b <- npregbw(y ~ x, data = d, bws = .2, bandwidth.compute = FALSE)
  keep <- d$x > .9
  analysis <- function(fun) { keep <- d$x > .4; fun(subset = keep) }
  lsq <- function(...) nplsqreg(b, data = d, ..., delta = .5,
                               bandwidth.compute = FALSE)
  expect_length(fitted(analysis(function(...) lsq(...))), sum(d$x > .4))
  d$z <- seq_len(nrow(d))/nrow(d)
  index <- function(...) npindex(y ~ x + z, data = d, bws = c(1, .5, .2),
                                 bandwidth.compute = FALSE, se = FALSE, ...)
  expect_identical(vapply(1:2, function(i)
    length(fitted(index(subset = d$grp == i))), 0L), c(45L, 30L))
})

test_that("native one-call response labels do not contain data values", {
  set.seed(3708)
  X <- data.frame(x = runif(50)); response <- sin(X$x)
  fit <- npreg(txdat = X, tydat = response, bws = .2,
               bandwidth.compute = FALSE)
  expect_identical(fit$bws$ynames, "response")
  expect_identical(fit$bws$call$ydat, quote(response))
})

# Bounded transport witness: fail at the public seam before a missing master
# binding can become a collective hang. The external repair proof also runs
# the actual fixed/search tests and same-pool recovery at one/three workers.
test_that("native conditional-moment call labels are values before transport", {
  skip_on_cran()
  skip_if_not(npRmpi:::.npRmpi_autodispatch_active(), "requires the session pool")
  testthat::local_mocked_bindings(.npRmpi_distributed_call_impl =
    function(mc, ...) list(observed = mc), .package = "npRmpi")
  x <- data.frame(x = seq_len(20)); y <- sin(x$x)
  model <- lm(y ~ x, data = data.frame(x, y), x = TRUE, y = TRUE)
  count <- new.env(parent = emptyenv()); count$n <- 0L
  label <- "master-only label"
  out <- npcmstest(model = model, xdat = x, ydat = y, B = 9L,
    call = { count$n <- count$n + 1L; label })
  expect_identical(count$n, 1L)
  expect_identical(out$observed$call, label)
})

test_that("pooled unconditional vector names support named newdata", {
  skip_on_cran()
  skip_if_not(npRmpi:::.npRmpi_autodispatch_active(), "requires the session pool")
  income <- seq(.05, .95, length.out = 30)
  for (family in c("npudens", "npudist")) {
    b <- do.call(paste0(family, "bw"),
      list(dat = quote(income), bws = .2, bandwidth.compute = FALSE),
      envir = environment())
    expect_identical(b$xnames, "income")
    fit <- do.call(family, list(bws = b,
      newdata = data.frame(income = c(.2, .5, .8))))
    expect_length(fitted(fit), 3L)
  }
})
