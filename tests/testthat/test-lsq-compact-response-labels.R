test_that("LSQ response transport preserves executable-option isolation", {
  y <- seq_len(20L)
  expect_identical(.nplsqreg_bw_response_arg(y, list(nmulti=1L)), quote(ydat))
  for (option in list(quote({ydat <- 0; 1}), as.name("ydat"),
                      compiler::compile(quote({ydat <- 0; 1}))))
    expect_identical(.nplsqreg_bw_response_arg(y, list(control=option)), y)
  expect_identical(.nplsqreg_bw_response_arg(y, list(control=list(quote(ydat)))), quote(ydat))
  for (options in list(list(first=quote(ydat),last=1L),
                       list(first=1L,last=quote(ydat))))
    expect_identical(.nplsqreg_bw_response_arg(y, options), y)
  expect_identical(.nplsqreg_bw_response_arg(quote(d$y), list(nmulti=1L)), quote(d$y))
  # Exercise the real argument expressions at all three producers, not a copy.
  found <- list()
  walk <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("do.call")) &&
        identical(x[[2L]], as.name("npregbw"))) found[[length(found)+1L]] <<- x
    if (is.call(x) || is.pairlist(x)) for (i in seq_along(x))
      if (!identical(x[[i]], quote(expr=))) walk(x[[i]])
  }
  walk(body(nplsqregbw.default)); expect_length(found, 3L)
  for (call in found) {
    e <- new.env(parent=environment(nplsqregbw.default))
    e$xdat <- data.frame(x=seq_along(y)); e$ydat <- y; e$bws <- .3
    e$reg.dots <- e$fixed.reg.dots <- list(control=quote({ydat <- rev(ydat); 1}))
    e$npregbw <- function(xdat, ydat, ...) { list(...); ydat }
    expect_identical(eval(call,e), y)
    e$ydat <- y
    e$reg.dots <- e$fixed.reg.dots <- list(nmulti=1L)
    e$npregbw <- function(xdat, ydat, ...) list(arg=substitute(ydat),value=ydat)
    actual <- eval(call,e)
    expect_identical(actual$arg, quote(ydat))
    expect_identical(actual$value, y)
  }
  for (bw in list(quote({ydat <- rev(ydat); .3}),
                  compiler::compile(quote({ydat <- rev(ydat); .3})))) {
    e <- new.env(parent=environment(nplsqregbw.default))
    e$xdat <- data.frame(x=seq_along(y)); e$ydat <- y; e$bws <- bw
    e$reg.dots <- list(nmulti=1L)
    e$npregbw <- function(xdat, ydat, ...) { list(...); ydat }
    expect_identical(eval(found[[2L]],e), y)
  }
})

test_that("LSQ uses compact native names and retains formula names and data", {
  if (identical(environmentName(environment(.nplsqreg_bw_response_arg)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE)
  withr::local_seed(75)
  d <- data.frame(x=seq(.02,.98,length.out=40)); d$y <- sin(3*d$x)
  seed <- .Random.seed
  b <- nplsqregbw(xdat=d["x"], ydat=d$y, bws=.3, delta=.1,
                 scale=rep(1,40), bandwidth.compute=FALSE)
  expect_identical(b$reg.bws$ynames, "ydat")
  expect_identical(b$ydat, d$y)
  f <- nplsqreg(b, se=FALSE)
  expect_identical(fitted(f), fitted(nplsqreg(f$bws,se=FALSE)))
  expect_true(all(is.finite(predict(f,newdata=data.frame(x=c(.2,.7))))))
  count <- 0L
  named <- nplsqregbw(bws=y~x,data={count <- count+1L; d},bw=.3,delta=.1,
                      scale=rep(1,40),bandwidth.compute=FALSE)
  expect_identical(count, 1L)
  expect_identical(named$reg.bws$ynames, "y")
  expect_identical(named$ydat,d$y)
  expect_identical(.Random.seed,seed)
})


test_that("missing-start and singleton LSQ producers keep compact response names", {
  if (identical(environmentName(environment(.nplsqreg_bw_response_arg)), "npRmpi"))
    skip_if_not(spawn_mpi_slaves())
  withr::local_options(np.messages=FALSE)
  withr::local_seed(7501)
  d <- data.frame(x=seq(.01,.99,length.out=40)); d$y <- sin(3*d$x)
  a <- list(xdat=d["x"],ydat=d$y,delta=.1,scale=rep(1,40),
            bandwidth.compute=TRUE,nmulti=1L,itmax=5L,optim.control=list(maxit=2L))
  missing.start <- do.call(nplsqregbw,a)
  expect_identical(missing.start$reg.bws$ynames,"ydat")
  a$bws <- .3; a$nomad <- TRUE; a$degree.min <- a$degree.max <- 1L
  a$degree.restarts <- 0L
  singleton <- do.call(nplsqregbw,a)
  expect_identical(singleton$reg.bws$ynames,"ydat")
})
