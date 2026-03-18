test_that("plot.plregression restores richer bws and infers overlay data from fitted object", {
  plot.plreg <- getFromNamespace(".np_plot_plregression", "npRmpi")

  bws.current <- structure(list(formula = NULL), class = "plbandwidth")
  bws.orig <- structure(list(formula = y ~ x | z), class = "plbandwidth")
  captured <- NULL

  local_mocked_bindings(
    .np_eval_call_arg = function(call, arg, caller_env = parent.frame()) bws.orig,
    .np_plot_call_method = function(method, bws, ..., where = "plot()") {
      captured <<- list(method = method, bws = bws, dots = list(...), where = where)
      invisible("ok")
    },
    .package = "npRmpi"
  )

  object <- structure(
    list(
      bws = bws.current,
      bw = bws.current,
      call = quote(npplreg(bws = saved_bws)),
      trainiseval = TRUE,
      evalx = data.frame(x = c(1, 2, 3)),
      evalz = data.frame(z = c(4, 5, 6)),
      mean = c(0.5, 1.5, 2.5),
      resid = c(0.1, 0.2, 0.3)
    ),
    class = "plregression"
  )

  plot.plreg(object, main = "kept")

  expect_identical(captured$bws, bws.orig)
  expect_identical(captured$dots$xdat, object$evalx)
  expect_identical(captured$dots$zdat, object$evalz)
  expect_equal(captured$dots$ydat, object$mean + object$resid)
  expect_identical(captured$dots$main, "kept")
  expect_identical(captured$where, "plot.plregression()")
})

test_that("plot.smoothcoefficient restores richer bws and infers overlay data from fitted object", {
  plot.scoef <- getFromNamespace(".np_plot_smoothcoefficient", "npRmpi")

  bws.current <- structure(list(formula = NULL), class = "scbandwidth")
  bws.orig <- structure(list(formula = y ~ x | z), class = "scbandwidth")
  captured <- NULL

  local_mocked_bindings(
    .np_eval_call_arg = function(call, arg, caller_env = parent.frame()) bws.orig,
    .np_plot_call_method = function(method, bws, ..., where = "plot()") {
      captured <<- list(method = method, bws = bws, dots = list(...), where = where)
      invisible("ok")
    },
    .package = "npRmpi"
  )

  object <- structure(
    list(
      bws = bws.current,
      call = quote(npscoef(bws = saved_bws)),
      trainiseval = TRUE,
      eval = list(
        exdat = data.frame(x = c(1, 2, 3)),
        ezdat = data.frame(z = c(4, 5, 6))
      ),
      mean = c(0.5, 1.5, 2.5),
      resid = c(0.1, 0.2, 0.3)
    ),
    class = "smoothcoefficient"
  )

  plot.scoef(object, col = "red")

  expect_identical(captured$bws, bws.orig)
  expect_identical(captured$dots$xdat, object$eval$exdat)
  expect_identical(captured$dots$zdat, object$eval$ezdat)
  expect_equal(captured$dots$ydat, object$mean + object$resid)
  expect_identical(captured$dots$col, "red")
})

test_that("plot.plregression recovers training data for direct formula fits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(3001)
  n <- 50
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- 1 + 0.6 * x + sin(z) + rnorm(n, sd = 0.1)

  fit <- npplreg(y ~ x | z, regtype = "ll")

  assign("overlay_calls", 0L, envir = .GlobalEnv)
  trace(".np_plot_overlay_points_1d",
        tracer = quote({
          assign("overlay_calls", get("overlay_calls", envir = .GlobalEnv) + 1L, envir = .GlobalEnv)
        }),
        where = asNamespace("npRmpi"),
        print = FALSE)
  on.exit(untrace(".np_plot_overlay_points_1d", where = asNamespace("npRmpi")), add = TRUE)
  on.exit(rm("overlay_calls", envir = .GlobalEnv), add = TRUE)

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  out <- suppressWarnings(plot(
    fit,
    perspective = FALSE,
    plot.behavior = "plot-data",
    plot.errors.method = "none",
    pch = 4,
    cex = 0.7
  ))

  expect_type(out, "list")
  overlay_calls <- get("overlay_calls", envir = .GlobalEnv)
  expect_gte(overlay_calls, 1L)
})
