test_that("selector call finalization preserves syntax and is idempotent", {
  finish <- getFromNamespace(".np_bw_selector_finish_call", "npRmpi")
  owner <- new.env(parent = baseenv())
  original <- quote(npregbw(formula = quote(y ~ x), data = handoff(),
    ... = pairlist(na.action = handoff(), bws = .5)))
  # Use a real pairlist, as match.call(expand.dots = FALSE) does.
  original[["..."]] <- pairlist(na.action = quote(handoff()), bws = .5)
  environment(original) <- owner
  expressions <- list(data = quote(d), na.action = quote(policy))
  got <- finish(original, expressions, quote(y ~ x))
  expect_identical(got$data, quote(d))
  expect_identical(got[["..."]]$na.action, quote(policy))
  expect_identical(got$formula, quote(y ~ x))
  expect_identical(environment(got), owner)
  expect_identical(finish(got, expressions, quote(y ~ x)), got)
  expect_identical(names(got), names(original))
  expect_null(finish(original, list(data = NULL))$data)
  reader <- getFromNamespace(".np_bw_selector_formula_call", "npRmpi")
  expect_identical(reader(original, "npregbw", environment()), original)
})

test_that("direct formula fits consume their already-finalized worker bandwidth", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  ns <- asNamespace("npRmpi")
  resolve <- getFromNamespace(".npRmpi_lease_resolve", "npRmpi")
  fields <- getFromNamespace(".npRmpi_lease_handle_fields", "npRmpi")
  rows <- list()
  testthat::local_mocked_bindings(.npRmpi_lease_resolve = function(value, consumer, argument) {
    answer <- resolve(value, consumer, argument)
    if (!is.null(fields(value)))
      rows[[length(rows) + 1L]] <<- list(consumer = consumer, valid = !is.null(answer))
    answer
  }, .package = "npRmpi")
  i <- seq_len(40L)
  d <- data.frame(x = sin(i * .7), y = cos(i * .3), z = cos(i * .4))
  count <- 0L
  factory <- function() { count <<- count + 1L; d }
  reg <- npreg(y ~ x, data = factory(), bws = .5, bandwidth.compute = FALSE, se = FALSE)
  expect_identical(count, 1L)
  npplreg(y ~ x | z, data = d, bws = matrix(.5, 2L, 1L), bandwidth.compute = FALSE, se = FALSE)
  npscoef(y ~ x | z, data = d, bws = .5, bandwidth.compute = FALSE, se = FALSE)
  for (family in c("npreg", "npplreg", "npscoef")) {
    hits <- Filter(function(x) identical(x$consumer, family), rows)
    expect_true(length(hits) > 0L, info = family)
    expect_true(all(vapply(hits, `[[`, logical(1L), "valid")), info = family)
  }
  expect_identical(reg$bws$call$data, quote(factory()))
  expect_identical(reg$bws$call$formula, quote(y ~ x))
})

test_that("a callback replay is not the outer selector's own constructor", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  # The callback owns only its observation state, not the test runner or the
  # fitted objects/closures retained by preceding tests in this file.
  observer <- new.env(parent = baseenv())
  observer$child <- NULL
  observer$busy <- FALSE
  policy <- eval(quote(function(frame) {
    if (!busy) {
      busy <<- TRUE
      on.exit({ busy <<- FALSE }, add = TRUE)
      frames <- sys.frames()
      at <- which(vapply(seq_along(frames), function(i)
        identical(sys.function(i), get(".np_eval_bw_call",
          envir = asNamespace("npRmpi"))), logical(1L)))
      testthat::expect_gt(length(at), 0L)
      owner <- frames[[at[[length(at)]]]]
      replay <- get("native", owner)$call
      # The prepared-frame capture belongs to one fitter invocation and is
      # intentionally single-use. The child replays the exact formula/data
      # promises, with its own ordinary model-frame transaction.
      replay[[".np.formula.state"]] <- NULL
      if (".np.formula.state" %in% names(replay[["..."]]))
        replay[["..."]][[".np.formula.state"]] <- NULL
      got <- eval(replay, get("value.env", owner))
      # Observe the ownership contract without retaining the synthetic child's
      # handoff closures in the policy later transported with the outer fit.
      child <<- list(is.bandwidth = inherits(got, "rbandwidth"),
                     formula.head = got$call$formula[[1L]],
                     data.is.call = is.call(got$call$data))
    }
    stats::na.omit(frame)
  }), observer)
  i <- seq_len(24L); d <- data.frame(x = sin(i), y = cos(i))
  fit <- npreg(y ~ x, data = d, na.action = policy, bws = .5,
               bandwidth.compute = FALSE, se = FALSE)
  child <- observer$child
  expect_true(child$is.bandwidth)
  expect_identical(child$formula.head, quote(quote))
  expect_true(child$data.is.call)
  expect_identical(fit$bws$call$formula, quote(y ~ x))
  expect_identical(fit$bws$call$data, quote(d))
})
