test_that("private regression dispatch retains its declared package leaf", {
  make <- getFromNamespace(".npregbw_owned_dispatch_call", "npRmpi")
  generic <- getFromNamespace(".npRmpi_autodispatch_as_generic_call", "npRmpi")
  for (owner in c("npregbw.NULL", "npregbw.default", "npregbw.rbandwidth")) {
    mc <- quote(npregbw.rbandwidth(xdat = stop("forced"), ydat = stop("forced"),
                                 bws = stop("forced"), degree = 2L))
    mc[[1L]] <- as.name(owner)
    context <- list(expanded = mc, public = FALSE)
    z <- make(context, owner)
    expect_identical(z[[1L]], as.name(owner))
    expect_identical(as.list(z)[-1L], as.list(mc)[-1L])
    expect_true(attr(z, ".npRmpi.bound.call.owner", exact = TRUE))
    context$public <- TRUE
    expect_identical(make(context, owner), generic("npregbw", mc))
    matched <- quote(npregbw.NULL(xdat = xdat, ydat = ydat, ... = pairlist(degree = 2L)))
    expect_identical(make(context, owner, matched), generic("npregbw", matched))
  }
  bad <- list(expanded = quote(foreign(xdat = stop("forced"))), public = FALSE)
  expect_error(make(bad, "npregbw.rbandwidth"), "private regression bandwidth dispatch owner")
  expect_error(make(bad, "foreign"), "private regression bandwidth dispatch owner")
  npregbw <- function(...) stop("ambient owner must not be used")
  expect_true(attr(make(list(expanded = quote(npregbw.default(xdat = x)), public = FALSE),
                        "npregbw.default"), ".npRmpi.bound.call.owner", exact = TRUE))
})


test_that("private regression fit dispatch retains its declared leaf", {
  make <- getFromNamespace(".npreg_owned_dispatch_call", "npRmpi")
  call <- quote(npreg.rbandwidth(bws = stop("forced"), txdat = stop("forced"),
                                tydat = stop("forced")))
  context <- list(expanded = call, public = FALSE)
  bound <- make(context)
  expect_identical(bound[[1L]], call[[1L]])
  expect_identical(as.list(bound)[-1L], as.list(call)[-1L])
  expect_true(attr(bound, ".npRmpi.bound.call.owner", exact = TRUE))
  context$public <- TRUE
  expect_identical(make(context), call)
  expect_error(make(list(expanded = quote(foreign()), public = FALSE)),
               "invalid private regression fit dispatch owner", fixed = TRUE)
})

test_that("private fit call site preserves ownership before transport", {
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  scope <- new.env(parent = ns)
  for (name in c(".npreg_auxiliary", ".npreg_fit_impl", ".np_regression_forward")) {
    fn <- get(name, ns); environment(fn) <- scope; scope[[name]] <- fn
  }
  scope$.npRmpi_autodispatch_call <- function(mc, caller_env, comm = 1L, owner.name = NULL)
    stop(structure(list(message = "captured fit dispatch", call = NULL,
                        invocation = mc, owner = owner.name),
                   class = c("fit_dispatch_capture", "error", "condition")))
  x <- data.frame(x = seq(.05, .95, length.out = 24))
  y <- sin(x$x)
  b <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE)
  got <- tryCatch(scope$.npreg_auxiliary(b, txdat = x, tydat = y),
                  fit_dispatch_capture = identity)
  expect_s3_class(got, "fit_dispatch_capture")
  expect_identical(got$owner, "npreg.rbandwidth")
  expect_identical(got$invocation[[1L]], quote(npreg.rbandwidth))
  expect_true(attr(got$invocation, ".npRmpi.bound.call.owner", exact = TRUE))
  expect_true(got$invocation$.np.defer.empty.rows)
})

test_that("namespace-only public partially linear fitting retains its child owner", {
  skip_on_cran()
  # Isolate package masking and coordinator progress from the suite's pool.
  result <- npRmpi_run_isolated_contract(c(
    "loadNamespace('npRmpi')",
    "npRmpi::npRmpi.init(nslaves=1L, quiet=TRUE)",
    "on.exit(npRmpi::npRmpi.quit(force=TRUE), add=TRUE)",
    "d <- data.frame(x=seq(.03,.97,length.out=31), z=sin(seq_len(31)))",
    "d$y <- sin(3*d$x)+.5*d$z+cos(seq_len(31))*.1",
    "options(np.messages=FALSE)",
    "b <- npRmpi::npplregbw(y~z|x,data=d,bws=matrix(.3,2,1),bandwidth.compute=FALSE)",
    "reference <- fitted(npRmpi::npplreg(bws=b))",
    "testthat::local_mocked_bindings(.np_progress_enabled=function(...) TRUE, .package='npRmpi')",
    "got <- fitted(npRmpi::npplreg(bws=b))",
    "stopifnot(isTRUE(all.equal(got,reference,tolerance=0)))",
    "cat('PRIVATE_FIT_PUBLIC_OWNER_PASS\n')"
  ), marker="PRIVATE_FIT_PUBLIC_OWNER_PASS", timeout=90L)
  if (is.null(result)) skip("installed subprocess library unavailable")
  expect_identical(result$status, 0L, info=paste(result$output, collapse="\n"))
  expect_true(result$witnessed)
})
