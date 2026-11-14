test_that("missing-bandwidth regression avoids literal response labels", {
  withr::local_options(np.messages = FALSE)

  d <- data.frame(x = seq(.02,.98,length.out=40L))
  d$y <- sin(3*d$x)+.1*cos(seq_len(nrow(d)))
  observed <- numeric()
  original <- rbandwidth
  testthat::local_mocked_bindings(rbandwidth=function(ynames,...) {
    observed <<- c(observed,sum(nchar(ynames)))
    original(ynames=ynames,...)
  }, .package="np")
  for(type in c("lc","ll","lp")) {
    observed <- numeric(); nx <- ny <- 0L
    rx <- function(){nx <<- nx+1L;d["x"]}
    ry <- function(){ny <<- ny+1L;d$y}
    set.seed(74040)
    args <- list(xdat=quote(rx()),ydat=quote(ry()),regtype=type,nmulti=1L,itmax=5L)
    if(type=="lp")args <- c(args,list(degree=2L,bernstein.basis=TRUE))
    fit <- do.call(npregbw,args)
    expect_identical(nx,1L);expect_identical(ny,1L)
    expect_identical(fit$ynames,"ry()")
    expect_gt(length(observed),0L)
    expect_lt(max(observed),100)
  }
})


test_that("missing-bandwidth response forcing retains condition order", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(x = seq(.05, .95, length.out = 20L))
  events <- character()
  response <- function() {
    events <<- c(events, "response")
    stop("response-boundary")
  }
  for (extra in list(list(nomad.opts = "invalid"), list(remin = TRUE))) {
    events <- character()
    warnings <- character()
    expect_error(withCallingHandlers(
      do.call(npregbw, c(list(xdat = x, ydat = quote(response())), extra)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), "response-boundary")
    expect_identical(events, "response")
    expect_length(warnings, 0L)
  }
})


test_that("NULL regression preserves historical language response evaluation", {
  withr::local_options(np.messages = FALSE, np.tree = FALSE)
  x <- data.frame(x = seq(.05, .95, length.out = 20L))
  expression <- quote(sin(seq(.15, 2.85, length.out = 20L)))
  y <- eval(expression)
  set.seed(74141)
  reference <- npregbw(xdat = x, ydat = y, nmulti = 1L, itmax = 3L)
  seed <- .Random.seed
  set.seed(74141)
  actual <- npregbw(xdat = x, ydat = expression, nmulti = 1L, itmax = 3L)
  expect_identical(actual$bw, reference$bw)
  expect_identical(.Random.seed, seed)
  expect_identical(actual$ynames, "expression")
})

test_that("NULL response handoff preserves values beside executable options", {
  implementation <- .npregbw_null_impl
  calls <- as.list(body(implementation))[-1L]
  selected <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name("tbw")) && is.call(x[[3L]]) &&
    identical(x[[3L]][[1L]], as.name("do.call")), calls)
  expect_length(selected, 1L)
  handoff <- selected[[1L]][[3L]]
  y <- sin(seq(.1, 3, length.out = 20L))
  execute <- function(response, options, public = TRUE) {
    frame <- new.env(parent = environment())
    frame$xdat <- data.frame(x = seq_along(y))
    frame$ydat <- response
    frame$Y <- y
    frame$bws <- 0
    frame$dots <- options
    frame$.np.context <- list(public = public)
    consume <- function(xdat, ydat, bws, scale.factor.search.lower = NULL, nmulti = NULL, ...) {
      label <- substitute(ydat)
      force(scale.factor.search.lower)
      force(nmulti)
      list(value = ydat, label = label)
    }
    frame$npregbw.default <- frame$.npregbw_auxiliary <- consume
    eval(handoff, frame)
  }
  for (public in c(TRUE, FALSE)) {
    plain <- execute(y, list(scale.factor.search.lower = .1), public)
    expect_identical(plain$value, y)
    expect_identical(plain$label, quote(ydat))
    for (option in list(quote({ydat <- rev(ydat); .1}),
                        quote({rm(ydat); .1}),
                        compiler::compile(quote({ydat <- rev(ydat); .1})))) {
      actual <- execute(y, list(scale.factor.search.lower = option), public)
      expect_identical(actual$value, y)
      expect_identical(actual$label, y)
    }
    other <- execute(y, list(nmulti = quote({ydat <- rev(ydat); 1L})), public)
    expect_identical(other$value, y)
    expect_identical(other$label, y)
    nested <- execute(y, list(nomad.opts = list(unused = quote(ydat))), public)
    expect_identical(nested$value, y)
    expect_identical(nested$label, quote(ydat))
    for (response in list(quote(Y), quote(identity(Y)), compiler::compile(quote(Y))))
      expect_identical(execute(response, list(), public)$value, y)
  }
})

test_that("quoted search options cannot replace the resolved regression response", {
  withr::local_options(np.messages = FALSE, np.tree = FALSE)
  x <- data.frame(x = seq(.02, .98, length.out = 40L))
  y <- sin(3*x$x) + .1*cos(seq_len(nrow(x)))
  set.seed(74241)
  reference <- npregbw(xdat = x, ydat = y, scale.factor.search.lower = .1,
                       nmulti = 1L, itmax = 5L)
  after <- .Random.seed
  for (option in list(quote({ydat <- rev(ydat); .1}),
                      quote({rm(ydat); .1}),
                      compiler::compile(quote({ydat <- rev(ydat); .1})))) {
    set.seed(74241)
    actual <- npregbw(xdat = x, ydat = y, scale.factor.search.lower = option,
                     nmulti = 1L, itmax = 5L)
    expect_identical(actual$bw, reference$bw)
    expect_identical(actual$fval, reference$fval)
    expect_identical(.Random.seed, after)
  }
})
