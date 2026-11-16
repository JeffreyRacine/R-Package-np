test_that("direct conditional NULL constructors publish reusable final calls", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  resolve <- get(".npRmpi_lease_resolve", ns)
  i <- seq_len(72L); x <- data.frame(x = sin(i * .7)); y <- cos(i * .3)
  for (family in c("npcdens", "npcdist")) {
    b <- get(paste0(family, "bw"), ns)(xdat = x, ydat = y,
      nmulti = 1L, itmax = 1L, powell.remin = FALSE)
    expect_identical(b$call$xdat, quote(x))
    expect_identical(b$call$ydat, quote(y))
    expect_false(is.null(resolve(b, family, "bws")))
    f <- get(family, ns)(bws = b, se = FALSE)
    expect_true(all(is.finite(fitted(f))))
    expect_false(is.null(resolve(b, family, "bws")))
  }
})

test_that("conditional selectors publish their final public call", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  native <- getFromNamespace(".npRmpi_bw_native_selector_plan", "npRmpi")
  finish <- getFromNamespace(".np_bw_selector_formula_call", "npRmpi")
  resolve <- getFromNamespace(".npRmpi_lease_resolve", "npRmpi")
  seen <- list(native = list(), null = list(), valid = character())
  testthat::local_mocked_bindings(
    .npRmpi_bw_native_selector_plan = function(...) {
      z <- native(...)
      if (!is.null(z)) seen$native[[length(seen$native) + 1L]] <<- z
      z
    },
    .np_bw_selector_formula_call = function(call, generic.name, method.frame,
                                            method.name = "formula") {
      z <- finish(call, generic.name, method.frame, method.name)
      if (identical(method.name, "NULL"))
        seen$null[[length(seen$null) + 1L]] <<- z
      z
    },
    .npRmpi_lease_resolve = function(value, consumer, argument) {
      z <- resolve(value, consumer, argument)
      if (!is.null(z)) seen$valid <<- c(seen$valid, consumer)
      z
    }, .package = "npRmpi")
  i <- seq_len(24L); x <- data.frame(x = sin(i * .7)); y <- cos(i * .3)
  for (family in c("npcdens", "npcdist")) {
    seen <- list(native = list(), null = list(), valid = character())
    fit <- get(family, ns)(txdat = x, tydat = y, bws = c(10, 10),
      bwtype = "generalized_nn", bwscaling = FALSE, bandwidth.compute = FALSE,
      regtype = "ll", se = FALSE)
    expect_true(length(seen$native) > 0L, info = family)
    expect_true(family %in% seen$valid, info = family)
    expect_identical(fit$bws$call$xdat, quote(x))
    expect_identical(fit$bws$call$ydat, quote(y))
  }
  # Ordinary LC primary fits dispatch the whole fit and do not create a
  # reusable constructor lease. Native mode exercises the actual NULL selector
  # chain; its public call must be finalized even though publication is excluded.
  y <- factor(rep(c("a", "b"), 12L))
  seen <- list(native = list(), null = list(), valid = character())
  fit <- npconmode(txdat = x, tydat = y, exdat = x[1:3, , drop = FALSE],
                   nmulti = 1L, itmax = 1L, powell.remin = FALSE)
  expect_true(length(seen$null) > 0L)
  expect_true(any(vapply(seen$null, function(z)
    identical(z$xdat, quote(x)) && identical(z$ydat, quote(y)), logical(1L))))
})

test_that("temporary conditional ownership is positive only on the real chain", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI runtime not available")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  temporary <- getFromNamespace(".npRmpi_bw_publication_temporary", "npRmpi")
  publication <- getFromNamespace(".npRmpi_lease_publication_plan", "npRmpi")
  resolve <- getFromNamespace(".npRmpi_lease_resolve", "npRmpi")
  decisions <- logical(); publications <- character()
  testthat::local_mocked_bindings(
    .npRmpi_bw_publication_temporary = function(...) {
      z <- temporary(...); decisions <<- c(decisions, z); z
    },
    .npRmpi_lease_publication_plan = function(...) {
      z <- publication(...); publications <<- c(publications, z$kind); z
    }, .package = "npRmpi")
  i <- seq_len(24L); d <- data.frame(x = sin(i * .7), y = cos(i * .3))
  for (family in c("npqreg", "npconmode")) {
    dd <- d
    if (family == "npconmode") dd$y <- factor(rep(c("a", "b"), 12L))
    h <- c(.5, if (family == "npconmode") .2 else .6)
    decisions <- logical(); publications <- character()
    fit <- get(family, ns)(y ~ x, data = dd, bws = h,
      bandwidth.compute = FALSE, se = FALSE)
    expect_true(any(decisions), info = family)
    expect_false("leased_public" %in% publications, info = family)
    primary <- if (family == "npqreg") "npcdist" else "npcdens"
    b <- get(paste0(primary, "bw"), ns)(y ~ x, data = dd, bws = h,
                                          bandwidth.compute = FALSE)
    expect_false(is.null(resolve(b, primary, "bws")))
    decisions <- logical(); publications <- character()
    fit <- get(family, ns)(bws = b, exdat = dd[1:3, "x", drop = FALSE], se = FALSE)
    expect_false(any(decisions), info = family)
    expect_false("leased_public" %in% publications, info = family)
    expect_false(is.null(resolve(b, primary, "bws")))
    expect_s3_class(get(primary, ns)(bws = b, se = FALSE),
                    if (primary == "npcdist") "condistribution" else "condensity")
  }
  # Same-name user wrappers and a nested data callback cannot acquire the
  # secondary owner's certificate or overwrite the inner public call.
  npqreg.default <- function() npcdistbw(y ~ x, data = d, bws = c(.5, .6),
                                        bandwidth.compute = FALSE)
  decisions <- logical(); b <- npqreg.default()
  expect_false(any(decisions)); expect_false(is.null(resolve(b, "npcdist", "bws")))
  rm(npqreg.default)
  child.data <- d[1:16, , drop = FALSE]
  child.data$x <- child.data$x / 2
  expected.child <- npcdistbw(xdat = child.data["x"], ydat = child.data["y"],
                              bws = c(.5, .6), bandwidth.compute = FALSE)
  child <- NULL; calls <- 0L
  factory <- function() {
    calls <<- calls + 1L
    child <<- npcdistbw(xdat = child.data["x"], ydat = child.data["y"], bws = c(.5, .6),
                       bandwidth.compute = FALSE)
    d
  }
  fit <- npqreg(y ~ x, data = factory(), bws = c(.5, .6),
                 bandwidth.compute = FALSE, se = FALSE)
  expect_identical(calls, 1L)
  fingerprint <- getFromNamespace(".npRmpi_autodispatch_fingerprint", "npRmpi")
  expect_identical(fingerprint(child$call), fingerprint(expected.child$call))
  expect_identical(child$bw, expected.child$bw)
  expect_false(is.null(resolve(child, "npcdist", "bws")))
})
