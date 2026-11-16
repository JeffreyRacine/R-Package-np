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
