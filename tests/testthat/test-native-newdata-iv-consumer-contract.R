test_that("IV derivative native evaluation keeps positional density consumers", {
  # Stop at the next (unchanged) regression stage: no search/iteration is
  # needed to exercise the four real PDF/CDF prediction consumers.
  if (exists("spawn_mpi_slaves", mode = "function")) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  fun <- npregivderiv.default
  owner <- new.env(parent = environment(fun))
  environment(fun) <- owner
  # Keep this instrumented prefix on the coordinator. Its real density and
  # distribution consumers still exercise their own installed MPI dispatch.
  # Dispatching the whole cloned IV function would discard these local hooks.
  owner$.npRmpi_autodispatch_active <- function() FALSE
  seen <- list()
  owner$predict <- function(object, ...) {
    args <- list(...)
    expect_true("edat" %in% names(args))
    expect_false("newdata" %in% names(args))
    value <- stats::predict(object, ...)
    expect_equal(value, stats::predict(object, edat = unname(as.matrix(args$edat))),
                 tolerance = 0)
    seen[[length(seen) + 1L]] <<- value
    value
  }
  owner$.npreg_complete <- function(...) stop("IV_PREDICTION_CONSUMERS_COMPLETE")
  z <- seq(-1, 1, length.out = 30)
  w <- z + sin(seq_along(z))/10
  y <- sin(z)
  expect_error(fun(y=y, z=z, w=w, zeval=data.frame(different=c(-.6,.1,.7)),
                   iterate.max=2, nmulti=1),
               "IV_PREDICTION_CONSUMERS_COMPLETE", fixed=TRUE)
  expect_identical(length(seen), 4L)
  expect_identical(lengths(seen), c(30L,30L,3L,3L))
})
