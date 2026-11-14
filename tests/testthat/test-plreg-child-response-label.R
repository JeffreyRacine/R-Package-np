test_that("partially linear children avoid literal response labels", {
  withr::local_options(np.messages = FALSE)
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  d <- data.frame(z = seq(.04,.96,length.out=32L),
                  x = cos(seq(.1,5,length.out=32L)))
  d$x2 <- sin(seq(.2,8,length.out=nrow(d)))
  d$y <- sin(3*d$z) + .2*d$x - .1*d$x2
  observed <- numeric()
  original <- rbandwidth
  testthat::local_mocked_bindings(rbandwidth=function(ynames,...) {
    observed <<- c(observed,sum(nchar(ynames)))
    original(ynames=ynames,...)
  }, .package="npRmpi")
  nx <- ny <- nz <- 0L
  rx <- function(){nx <<- nx+1L;d[c("x", "x2")]}
  ry <- function(){ny <<- ny+1L;d$y}
  rz <- function(){nz <<- nz+1L;d["z"]}
  b <- npplregbw(xdat=rx(),ydat=ry(),zdat=rz(),
                bws=matrix(.3,3L,1L),bandwidth.compute=FALSE)
  expect_identical(c(nx,ny,nz),rep(1L,3L))
  expect_identical(b$ynames,"ydat")
  expect_gt(length(observed),0L)
  expect_lt(max(observed),100)
  for(i in seq_along(b$bw)) {
    samples <- b$bw[[i]][[".np.native.training",exact=TRUE]]
    expect_identical(samples$xdat,d["z"])
    expect_identical(samples$ydat,if(i==1L)d$y else d[[c("x", "x2")[i-1L]]])
    expect_identical(b$bw[[i]]$call[["ydat"]], samples$ydat)
    reference <- fitted(npreg(b$bw[[i]]))
    child <- b$bw[[i]]
    child[[".np.native.training"]] <- NULL
    expect_identical(fitted(npreg(child)), reference)
    expect_length(reference,nrow(d))
  }
})


test_that("PL call retention is limited to the admitted builders", {
  x <- data.frame(x = 1:3)
  y <- c(4, 5, 6)
  children <- list(list(call = quote(npregbw(xdat = x, ydat = y))))
  ordinary <- .npplregbw_retain_child_training(children, x, y, x)
  literal <- .npplregbw_retain_child_training(children, x, y, x,
                                            literal.response.call = TRUE)
  expect_identical(ordinary[[1L]]$call, children[[1L]]$call)
  expect_identical(literal[[1L]]$call[["ydat"]], y)
  absent <- list(list(), list(call = quote(npregbw(bws = bw))))
  retained <- .npplregbw_retain_child_training(absent, x, y, x,
                                             literal.response.call = TRUE)
  expect_null(retained[[1L]]$call)
  expect_identical(retained[[2L]]$call, absent[[2L]]$call)
})

test_that("PL NOMAD children retain the actual independent response columns", {
  withr::local_options(np.messages=FALSE)
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  d <- data.frame(z=seq(.04,.96,length.out=32L),
                  x=cos(seq(.1,5,length.out=32L)),
                  x2=sin(seq(.2,8,length.out=32L)))
  d$y <- sin(3*d$z)+.2*d$x-.1*d$x2
  expected <- list(d$y,d$x,d$x2)
  observed <- list()
  original <- .npregbw_auxiliary
  testthat::local_mocked_bindings(.npregbw_auxiliary=function(...) {
    args <- list(...)
    observed[[length(observed)+1L]] <<- args$ydat
    original(...)
  }, .package="npRmpi")
  b <- npplregbw(xdat=d[c("x","x2")],ydat=d$y,zdat=d["z"],
    bws=matrix(.3,3L,1L),nomad=TRUE,degree.min=0L,degree.max=1L,
    search.engine="nomad",nmulti=1L,itmax=3L,nomad.opts=list(MAX_BB_EVAL=8L))
  expect_identical(observed,expected)
  expect_length(b$bw,3L)
  for(i in seq_along(expected)) {
    expect_identical(b$bw[[i]]$call[["ydat"]],expected[[i]])
    expect_identical(b$bw[[i]][[".np.native.training",exact=TRUE]]$ydat,expected[[i]])
  }
})


test_that("PL fixed builder hands off the independent response columns", {
  withr::local_options(np.messages = FALSE)
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  d <- data.frame(z = seq(.04, .96, length.out = 32L),
                  x = cos(seq(.1, 5, length.out = 32L)),
                  x2 = sin(seq(.2, 8, length.out = 32L)))
  d$y <- sin(3 * d$z) + .2 * d$x - .1 * d$x2
  expected <- list(d$y, d$x, d$x2)
  observed <- list()
  original <- .npregbw_auxiliary
  testthat::local_mocked_bindings(.npregbw_auxiliary = function(...) {
    args <- list(...)
    observed[[length(observed) + 1L]] <<- args$ydat
    original(...)
  }, .package = "npRmpi")
  npplregbw(xdat = d[c("x", "x2")], ydat = d$y, zdat = d["z"],
            bws = matrix(.3, 3L, 1L), bandwidth.compute = FALSE)
  expect_identical(observed, expected)
})

test_that("PL NOMAD producer labels stay compact without a handoff wrapper", {
  withr::local_options(np.messages = FALSE)
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  d <- data.frame(z = seq(.04, .96, length.out = 32L),
                  x = cos(seq(.1, 5, length.out = 32L)),
                  x2 = sin(seq(.2, 8, length.out = 32L)))
  d$y <- sin(3 * d$z) + .2 * d$x - .1 * d$x2
  observed <- numeric()
  original <- rbandwidth
  testthat::local_mocked_bindings(rbandwidth = function(ynames, ...) {
    observed <<- c(observed, sum(nchar(ynames)))
    original(ynames = ynames, ...)
  }, .package = "npRmpi")
  npplregbw(xdat = d[c("x", "x2")], ydat = d$y, zdat = d["z"],
            bws = matrix(.3, 3L, 1L), nomad = TRUE,
            degree.min = 0L, degree.max = 1L, search.engine = "nomad",
            nmulti = 1L, itmax = 3L, nomad.opts = list(MAX_BB_EVAL = 8L))
  expect_gt(length(observed), 0L)
  expect_lt(max(observed), 100)
})


test_that("PL call reconstruction preserves already-literal response arguments", {
  x <- data.frame(x = 1:3)
  y <- c(4, 5, 6)
  retained <- rev(y)
  child <- list(call = as.call(list(as.name("npregbw"), ydat = y)))
  actual <- .npplregbw_retain_child_training(list(child), x, retained, x,
                                            literal.response.call = TRUE)
  expect_identical(actual[[1L]]$call[["ydat"]], y)
})

test_that("PL NOMAD handoff retains response values beside executable options", {
  withr::local_options(np.messages = FALSE)
  if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  set.seed(7404)
  x <- round(runif(60L), 4)
  z <- round(runif(60L), 4)
  d <- data.frame(x = x, z = z, x2 = round(cos(5*x) + z, 4))
  d$y <- round(sin(3*x) + .5*z + rnorm(60L, sd=.2), 4)
  args <- list(xdat = d[c("x", "x2")], ydat = d$y, zdat = d["z"],
    nomad = TRUE, nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 8L))
  set.seed(74041)
  ordinary <- do.call(npplregbw, args)
  final.rng <- .Random.seed
  for (option in list(
      list(ckertype = quote({child.responses[[1L]]$values <-
        rev(child.responses[[1L]]$values); "gaussian"})),
      list(ckertype = compiler::compile(quote({child.responses[[1L]]$values <-
        rev(child.responses[[1L]]$values); "gaussian"}))),
      list(nomad.opts = quote({child.responses[[1L]]$values <-
        rev(child.responses[[1L]]$values); list(MAX_BB_EVAL = 8L)})))) {
    set.seed(74041)
    call.args <- args
    call.args[names(option)] <- option
    actual <- do.call(npplregbw, call.args, quote = TRUE)
    expect_identical(lapply(actual$bw, function(x) x$bw),
                     lapply(ordinary$bw, function(x) x$bw))
    expect_identical(lapply(actual$bw, function(x) x$fval),
                     lapply(ordinary$bw, function(x) x$fval))
    expect_identical(.Random.seed, final.rng)
  }
})


test_that("PL fixed handoff preserves literal calls beside executable options", {
  withr::local_options(np.messages = FALSE)
  if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  d <- data.frame(z = seq(.04, .96, length.out = 32L),
                  x = cos(seq(.1, 5, length.out = 32L)),
                  x2 = sin(seq(.2, 8, length.out = 32L)))
  d$y <- sin(3 * d$z) + .2 * d$x - .1 * d$x2
  args <- list(xdat = d[c("x", "x2")], ydat = d$y, zdat = d["z"],
    bws = matrix(.3, 3L, 1L), bandwidth.compute = FALSE)
  for (option in list(
      list(ckertype = quote({ydat <- rev(ydat); "gaussian"})),
      list(ckertype = compiler::compile(quote({ydat <- rev(ydat); "gaussian"}))),
      list(bwmethod = quote({ydat <- rev(ydat); "cv.ls"})))) {
    result <- do.call(npplregbw, c(args, option), quote = TRUE)
    expect_identical(result$bw[[1L]]$call[["ydat"]], d$y)
  }
  result <- do.call(npplregbw, c(args, list(ckertype = quote({
    xdat <- xdat[nrow(xdat):1L, , drop = FALSE]; "gaussian"}))), quote = TRUE)
  expect_identical(lapply(result$bw, function(x) x$call[["ydat"]]),
                   list(yzbw = d$y, rev(d$x), d$x2))
})
