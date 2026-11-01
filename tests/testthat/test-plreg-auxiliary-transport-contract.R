test_that("cell preparation carries categorical children to participating ranks", {
  skip_on_cran()
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)),
          "Transport ownership requires distributed workers")
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  withr::defer(close_mpi_slaves(force = TRUE))
  withr::local_options(np.messages = FALSE)
  ns <- asNamespace("npRmpi")
  command <- get(".npRmpi_bcast_cmd_expr", ns)
  remote <- get("mpi.remote.exec", ns)
  command(quote({
    options(np.test.plreg.transport = list())
    trace(".npplregbw_build_plbandwidth_bcast_payload",
      where = asNamespace("npRmpi"), print = FALSE, tracer = quote({
        seen <- getOption("np.test.plreg.transport")
        seen[[length(seen) + 1L]] <- list(rank = mpi.comm.rank(1L),
          x = xdat, y = ydat, degree = reg.args$degree)
        options(np.test.plreg.transport = seen)
      }))
  }), caller.execute = TRUE)
  withr::defer(command(quote({
    untrace(".npplregbw_build_plbandwidth_bcast_payload",
            where = asNamespace("npRmpi"))
    options(np.test.plreg.transport = NULL)
  }), caller.execute = TRUE))

  set.seed(691)
  n <- 24L
  z <- data.frame(z = runif(n))
  y <- sin(3 * z$z) + seq_len(n) / n
  cats <- list(factor(rep(c("30", "10", "50"), 8L),
                      levels = c("50", "10", "30", "70")),
               ordered(rep(c("30", "10", "50"), 8L),
                       levels = c("10", "30", "50", "70")))
  scores <- list(rep(c(3, 2, 1), 8L), rep(c(30, 10, 50), 8L))
  for (j in seq_along(cats)) for (upper in 0:1) {
    args <- list(xdat = data.frame(x = cats[[j]]), ydat = y, zdat = z,
      regtype = "lp", degree.select = "coordinate", degree.min = 0L,
      degree.max = upper, search.engine = "cell", nmulti = 1L,
      itmax = 2L, powell.remin = FALSE)
    set.seed(172)
    b <- do.call(npplregbw, args)
    args$xdat <- data.frame(x = scores[[j]])
    set.seed(172)
    numeric <- do.call(npplregbw, args)
    expect_identical(b$bw[[2L]][[".np.native.training"]]$ydat, cats[[j]])
    expect_identical(b$bw[[2L]]$bw, numeric$bw[[2L]]$bw)
    expect_identical(b$fval, numeric$fval)
    expect_identical(b$num.feval, numeric$num.feval)
    expect_identical(fitted(npplreg(b, se = FALSE)),
                     fitted(npplreg(numeric, se = FALSE)))
  }
  seen <- c(list(getOption("np.test.plreg.transport")),
    remote(getOption("np.test.plreg.transport"), simplify = FALSE))
  expect_true(all(vapply(seen, function(rank.events) {
    any(vapply(rank.events, function(event) {
      is.factor(event$x[[1L]]) && is.numeric(event$y) &&
        identical(levels(event$x[[1L]])[4L], "70")
    }, logical(1L)))
  }, logical(1L))))
})

test_that("cell children remain private after serialization and public rejection", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  withr::defer(close_mpi_slaves(force = TRUE))
  withr::local_options(np.messages = FALSE)
  set.seed(171)
  n <- 24L
  z <- data.frame(z = runif(n))
  x <- data.frame(x = factor(rep(c("30", "10", "50"), 8L),
                            levels = c("50", "10", "30", "70")))
  y <- cos(z$z) + seq_len(n) / n
  b <- npplregbw(xdat = x, ydat = y, zdat = z, regtype = "lp",
    degree.select = "coordinate", degree.min = 0L, degree.max = 1L,
    search.engine = "cell", nmulti = 1L, itmax = 2L, powell.remin = FALSE)
  cold <- unserialize(serialize(b, NULL))
  expect_identical(fitted(npplreg(cold, se = FALSE)),
                   fitted(npplreg(b, se = FALSE)))
  child <- cold$bw[[2L]]
  expect_error(npreg(child), "must be a numeric")
  expect_error(npregbw(xdat = z, ydat = x$x, bws = child,
                       bandwidth.compute = FALSE), "must be a numeric")
  expect_identical(fitted(npplreg(b, se = FALSE)),
                   fitted(npplreg(cold, se = FALSE)))
})

test_that("preparation failures preserve participating workers and session reuse", {
  skip_on_cran()
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)),
          "Transport recovery requires distributed workers")
  if (!spawn_mpi_slaves()) skip("MPI pool unavailable")
  withr::defer(close_mpi_slaves(force = TRUE))
  withr::local_options(np.messages = FALSE)
  ns <- asNamespace("npRmpi")
  command <- get(".npRmpi_bcast_cmd_expr", ns)
  remote <- get("mpi.remote.exec", ns)
  withr::defer(command(quote({
    if (inherits(get(".npplregbw_build_plbandwidth",
                     envir = asNamespace("npRmpi")), "functionWithTrace"))
      untrace(".npplregbw_build_plbandwidth", where = asNamespace("npRmpi"))
  }), caller.execute = TRUE))
  pids <- remote(Sys.getpid(), simplify = FALSE)
  set.seed(703)
  n <- 24L
  z <- data.frame(z = runif(n))
  y <- sin(3 * z$z) + seq_len(n) / n
  f <- factor(rep(c("30", "10", "50"), 8L),
              levels = c("50", "10", "30", "70"))
  args <- list(xdat = data.frame(x = f), ydat = y, zdat = z,
    regtype = "lp", degree.select = "coordinate", degree.min = 0L,
    degree.max = 1L, search.engine = "cell", nmulti = 1L,
    itmax = 2L, powell.remin = FALSE)
  for (fault in c("metadata", "rank1")) {
    if (fault == "metadata") {
      args$xdat <- data.frame(x = ordered(f, levels = levels(f)))
      expect_error(do.call(npplregbw, args), "finite and increasing")
    } else {
      command(quote({
        if (mpi.comm.rank(1L) == 1L)
          trace(".npplregbw_build_plbandwidth",
            where = asNamespace("npRmpi"), print = FALSE,
            tracer = quote(stop("injected preparation error", call. = FALSE)))
      }), caller.execute = TRUE)
      expect_error(do.call(npplregbw, args), "injected preparation error")
      command(quote({
        if (mpi.comm.rank(1L) == 1L)
          untrace(".npplregbw_build_plbandwidth", where = asNamespace("npRmpi"))
      }), caller.execute = TRUE)
    }
    expect_identical(remote(Sys.getpid(), simplify = FALSE), pids)
    args$xdat <- data.frame(x = f)
    b <- do.call(npplregbw, args)
    expect_true(all(is.finite(fitted(npplreg(b, se = FALSE)))))
  }
})
