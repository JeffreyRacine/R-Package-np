test_that("leased native replacement preserves effective data and user owners", {
  skip_on_cran()
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  withr::defer(close_mpi_slaves(force = TRUE))
  withr::local_options(np.messages = FALSE, np.tree = FALSE)
  ns <- asNamespace("npRmpi")
  command <- get(".npRmpi_bcast_cmd_expr", ns)
  remote <- get("mpi.remote.exec", ns)
  diagnostics <- get(".npRmpi_lease_diagnostics", ns)
  options(np.test.owner.hits = 0L)
  trace(".npRmpi_bw_restore_public_owner", where = ns, print = FALSE,
    tracer = quote(options(np.test.owner.hits = getOption("np.test.owner.hits", 0L) + length(originals))))
  withr::defer({untrace(".npRmpi_bw_restore_public_owner", where = ns); options(np.test.owner.hits = NULL)})
  command(quote({
    options(np.test.fit.completed = 0L)
    trace(".npreg_fit_impl", where = asNamespace("npRmpi"), print = FALSE,
      exit = quote(options(np.test.fit.completed = getOption("np.test.fit.completed", 0L) + 1L)))
  }), caller.execute = TRUE)
  withr::defer(command(quote({
    untrace(".npreg_fit_impl", where = asNamespace("npRmpi"))
    options(np.test.fit.completed = NULL)
  }), caller.execute = TRUE))
  pids <- remote(Sys.getpid(), simplify = FALSE)
  session <- diagnostics()$session
  set.seed(618)
  d <- data.frame(x = runif(40), z = runif(40))
  d$y <- sin(d$x) + d$z + seq_len(40)/40
  d2 <- d[1:30, ]; d2$y <- d2$y + 1 + .2*d2$x
  for (family in c("npreg", "npplreg", "npscoef", "npindex")) {
    ctor <- get(paste0(family, "bw")); fitfun <- get(family)
    semi <- family %in% c("npplreg", "npscoef")
    formula <- if (semi) y ~ x | z else if (family == "npindex") y ~ x + z else y ~ x
    bandwidth <- if (family == "npplreg") matrix(.4, 2, 1) else if (family == "npindex") c(1, .3, .4) else .4
    b <- ctor(formula, data = d, bws = bandwidth, bandwidth.compute = FALSE)
    if (is.na(session)) session <- diagnostics()$session
    original.frame <- b[[".np.formula.training"]]$frame
    original.owner <- environment(b$formula)
    args <- list(bws = b, txdat = d2[if (family == "npindex") c("x", "z") else "x"], tydat = d2$y, se = FALSE)
    if (semi) args$tzdat <- d2["z"]
    options(np.test.owner.hits = 0L)
    kind <- RNGkind()
    set.seed(791); got <- do.call(fitfun, args); seed <- .Random.seed
    if (family == "npindex") {
      # Index fitting retains native data before dispatch, invalidating the
      # constructor lease on replacement even in the predecessor.
      expect_identical(getOption("np.test.owner.hits"), 0L)
    } else expect_true(getOption("np.test.owner.hits") > 0L, info = family)
    expect_identical(got$bws[[".np.formula.training"]]$frame$y, d2$y, info = family)
    expect_identical(row.names(got$bws[[".np.formula.training"]]$frame), row.names(d2), info = family)
    expect_identical(environment(got$bws$formula), original.owner)
    expect_identical(environment(attr(got$bws[[".np.formula.training"]]$frame, "terms")), original.owner)
    expect_identical(b[[".np.formula.training"]]$frame, original.frame)
    args$bws <- unserialize(serialize(b, NULL))
    options(np.test.owner.hits = 0L)
    set.seed(791); cold <- do.call(fitfun, args)
    expect_identical(.Random.seed, seed)
    expect_identical(RNGkind(), kind)
    expect_equal(fitted(got), fitted(cold), tolerance = 0, info = family)
    expect_identical(getOption("np.test.owner.hits"), 0L)
    expect_equal(fitted(fitfun(got$bws, se = FALSE)), fitted(got), tolerance = 0, info = family)
    expect_equal(predict(got, newdata = d2[1:5, ]), predict(cold, newdata = d2[1:5, ]), tolerance = 0)
    expect_identical(remote(Sys.getpid(), simplify = FALSE), pids)
    expect_identical(diagnostics()$session, session)
    rm(got, cold, b, args)
    invisible(gc())
    get(".npRmpi_lease_drain_pending", ns)()
  }
  workers <- remote(getOption("np.test.fit.completed"), simplify = FALSE)
  expect_true(all(unlist(workers) > 0L))
  expect_identical(diagnostics()$pending_count, 0L)
})
