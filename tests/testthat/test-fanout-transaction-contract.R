# Pure protocol tests except for the explicitly opted-in real SIGINT subprocess.
fanout_contract_tx <- function(scheduler = "dynamic", n = 3L) {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  tx <- make(101L, 1L, n, scheduler, NULL, "test-session", "test-operation", 18432L)
  tx$phase <- "active"
  tx$rank[] <- "busy"
  tx$assigned[[1L]] <- 1L
  tx$tags[] <- 1L
  tx$stop.tag <- n + 1L
  tx
}

fanout_contract_accept <- function() {
  accept <- getFromNamespace(".npRmpi_fanout_accept", "npRmpi")
  env <- new.env(parent = environment(accept))
  env$.npRmpi_fanout_native <- function(...) NULL
  environment(accept) <- env
  accept
}

test_that("a real fanout SIGINT is not swallowed by ordinary error handlers", {
  skip_on_cran()
  skip_if(Sys.getenv("NP_RMPI_RUN_FANOUT_SIGINT_TESTS") != "TRUE",
          "opt-in real SIGINT subprocess")
  skip_on_os("windows")
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "installed npRmpi unavailable")
  result <- npRmpi_run_rscript_subprocess(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "main <- function() {",
    "  npRmpi.init(nslaves=1, quiet=TRUE)",
    "  on.exit(npRmpi.quit(force=TRUE), add=TRUE)",
    "  ns <- asNamespace('npRmpi')",
    "  original <- get('.npRmpi_fanout_receive', ns)",
    "  state <- new.env(parent=emptyenv())",
    "  replacement <- function(tx, ...) {",
    "    value <- original(tx, ...)",
    "    if (!state$sent && identical(value$kind, 'result')) {",
    "      state$sent <- TRUE",
    "      tools::pskill(Sys.getpid(), 2L)",
    "      Sys.sleep(0.01)",
    "    }",
    "    value",
    "  }",
    "  unlockBinding('.npRmpi_fanout_receive', ns)",
    "  assign('.npRmpi_fanout_receive', replacement, ns)",
    "  lockBinding('.npRmpi_fanout_receive', ns)",
    "  for (wrapper in c('try', 'error')) {",
    "    state$sent <- state$aborted <- state$resumed <- FALSE",
    "    run <- function() npRmpi:::mpi.iapplyLB(1:3, function(x) { Sys.sleep(0.02); x }, sleep=0.001)",
    "    withRestarts({",
    "      if (wrapper == 'try') try(run(), silent=TRUE)",
    "      else tryCatch(run(), error=function(e) NULL)",
    "      state$resumed <- TRUE",
    "    }, abort=function() { state$aborted <- TRUE; NULL })",
    "    stopifnot(state$sent, state$aborted, !state$resumed,",
    "      length(ls(get('.npRmpi_fanout_state', ns))) == 0L)",
    "    stopifnot(identical(npRmpi:::mpi.iapplyLB(1:3, identity, sleep=0.001), as.list(1:3)))",
    "  }",
    "  unlockBinding('.npRmpi_fanout_receive', ns)",
    "  assign('.npRmpi_fanout_receive', original, ns)",
    "  lockBinding('.npRmpi_fanout_receive', ns)",
    "}",
    "main()",
    "cat('REAL_SIGINT_REUSE_QUIT_OK\\n')"
  ), timeout = 20L, env = env, cleanup = FALSE)
  expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(any(grepl("REAL_SIGINT_REUSE_QUIT_OK", result$output, fixed = TRUE)),
              info = paste(result$output, collapse = "\n"))
})

test_that("private native actions expose literal registered calls", {
  bridge <- getFromNamespace(".npRmpi_fanout_native", "npRmpi")
  env <- new.env(parent = environment(bridge))
  env$.Call <- function(name, ...) list(name = name, args = list(...))
  environment(bridge) <- env
  for (action in c("begin", "send", "poll", "finish", "owner")) {
    args <- switch(action, begin = list(new.env()), send = list(1L, 2L), list())
    value <- do.call(bridge, c(list(action, 7), args))
    expect_identical(value$name, paste0("np_mpi_fanout_", action))
    expect_identical(value$args, c(list(7L), args, list(PACKAGE = "npRmpi")))
  }
  expect_error(bridge("unknown", 7L), "unknown private")
})

test_that("an abandoned cooperative cleanup is retained and explicitly closable", {
  skip_on_cran()
  skip_if(Sys.getenv("NP_RMPI_RUN_FANOUT_SIGINT_TESTS") != "TRUE",
          "opt-in cancellation subprocess")
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "installed npRmpi unavailable")
  result <- npRmpi_run_rscript_subprocess(c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "npRmpi.init(nslaves=1, quiet=TRUE)",
    "ns <- asNamespace('npRmpi')",
    "bind <- function(n, x) { unlockBinding(n, ns); assign(n, x, ns); lockBinding(n, ns) }",
    "mpi.bcast.cmd(local({",
    "  ns <- asNamespace('npRmpi')",
    "  original <- get('.npRmpi_fanout_worker_terminal', ns)",
    "  delayed <- function(...) { Sys.sleep(.4); original(...) }",
    "  unlockBinding('.npRmpi_fanout_worker_terminal', ns)",
    "  assign('.npRmpi_fanout_worker_terminal', delayed, ns)",
    "  lockBinding('.npRmpi_fanout_worker_terminal', ns)",
    "}))",
    "receive <- get('.npRmpi_fanout_receive', ns)",
    "bind('.npRmpi_fanout_receive', function(tx, ...) {",
    "  if (tx$phase == 'active') stop('collector witness')",
    "  receive(tx, ...)",
    "})",
    "options(npRmpi.session.recv.timeout=.01, np.messages=FALSE)",
    "e <- tryCatch(npRmpi:::mpi.iapplyLB(1:8, identity), error=identity)",
    "stopifnot(identical(conditionMessage(e), 'collector witness'))",
    "bind('.npRmpi_fanout_receive', receive)",
    "tx <- get('.npRmpi_fanout_retained', ns)(1L)",
    "stopifnot(tx$phase == 'quarantined', !tx$owner.active)",
    "gc(FALSE)",
    "stopifnot(identical(.Call('np_mpi_fanout_owner', 1L, PACKAGE='npRmpi'), tx))",
    "e <- tryCatch(npRmpi:::mpi.iapplyLB(1:8, identity), error=identity)",
    "stopifnot(inherits(e, 'error'))",
    "options(npRmpi.session.recv.timeout=3)",
    "npRmpi.quit(force=TRUE)",
    "stopifnot(is.null(get('.npRmpi_fanout_retained', ns)(1L)))",
    "cat('COOPERATIVE_QUARANTINE_EXPLICIT_CLOSE_OK\\n')"
  ), timeout = 20L, env = env, cleanup = FALSE)
  expect_identical(result$status, 0L, info = paste(result$output, collapse = "\n"))
  expect_true(any(grepl("COOPERATIVE_QUARANTINE_EXPLICIT_CLOSE_OK", result$output, fixed = TRUE)),
              info = paste(result$output, collapse = "\n"))
})

test_that("fanout results and receipts have distinct identity and state duties", {
  accept <- fanout_contract_accept()
  envelope <- getFromNamespace(".npRmpi_fanout_envelope", "npRmpi")
  header <- getFromNamespace(".npRmpi_fanout_header", "npRmpi")
  tx <- fanout_contract_tx()
  result <- envelope(header(tx), "result", 1L, list(42L))
  expect_identical(accept(tx, result, 1L, 1L)$value, list(42L))
  expect_identical(tx$rank, "busy")
  expect_identical(tx$seen, c(TRUE, FALSE, FALSE))
  expect_error(accept(tx, result, 1L, 1L), "unexpected")
  expect_identical(accept(tx, envelope(header(tx), "ready", 1L), 1L, tx$control)$kind, "ready")
  expect_identical(tx$rank, "ready")
  expect_error(accept(tx, envelope(header(tx), "terminal"), 1L, tx$control), "unexpected")
  tx$rank[] <- "stopping"
  tx$assigned[[1L]] <- integer()
  expect_identical(accept(tx, envelope(header(tx), "terminal"), 1L, tx$control)$kind, "terminal")
  expect_identical(tx$rank, "terminal")
  expect_error(accept(tx, envelope(header(tx), "terminal"), 1L, tx$control), "unexpected")
  expect_identical(accept(tx, envelope(header(tx), "closed"), 1L, tx$control)$kind, "closed")
  expect_identical(tx$rank, "closed")
  expect_error(accept(tx, envelope(header(tx), "closed"), 1L, tx$control), "unexpected")
})

test_that("foreign, malformed and task-tag-colliding replies cannot publish", {
  accept <- fanout_contract_accept()
  envelope <- getFromNamespace(".npRmpi_fanout_envelope", "npRmpi")
  header <- getFromNamespace(".npRmpi_fanout_header", "npRmpi")
  tx <- fanout_contract_tx(n = 18432L)
  tx$assigned[[1L]] <- tx$control
  tx$tags[] <- tx$control
  valid <- envelope(header(tx), "result", tx$control, list(9L))
  wrong <- list(valid, valid, valid, valid, valid)
  wrong[[1L]]$session <- "other-session"
  wrong[[2L]]$operation <- "previous-operation"
  wrong[[3L]]$tasks <- 1L
  wrong[[4L]]$tasks <- c(tx$control, tx$control)
  wrong[[5L]]$value <- 9L
  near <- valid
  names(near)[[1L]] <- "session_extra"
  duplicate <- c(valid, list(session = tx$session))
  wrong <- c(wrong, list(near, duplicate))
  for (value in wrong) expect_error(accept(tx, value, 1L, tx$control), "unexpected")
  expect_error(accept(tx, valid, 2L, tx$control), "unexpected")
  expect_error(accept(tx, valid, 1L, tx$control - 1L), "unexpected")
  expect_false(any(tx$seen))
  expect_null(accept(tx, wrong[[2L]], 1L, tx$control, discard = TRUE))
  expect_identical(accept(tx, valid, 1L, tx$control)$value, list(9L))
  expect_identical(accept(tx, envelope(header(tx), "ready", tx$control), 1L, tx$control)$kind, "ready")
  expect_true(tx$seen[[tx$control]])
})

test_that("consumed control receipts remain committed after collector failure", {
  accept <- fanout_contract_accept()
  envelope <- getFromNamespace(".npRmpi_fanout_envelope", "npRmpi")
  header <- getFromNamespace(".npRmpi_fanout_header", "npRmpi")
  tx <- fanout_contract_tx()
  expect_error(accept(tx, envelope(header(tx), "ready", 1L), 1L, tx$control), "unexpected")
  expect_identical(tx$rank, "ready")
  tx <- fanout_contract_tx("scatter")
  expect_error(accept(tx, envelope(header(tx), "terminal"), 1L, tx$control), "unexpected")
  expect_identical(tx$rank, "terminal")
  tx <- fanout_contract_tx("bundle")
  expect_null(accept(tx, envelope(header(tx), "result", 1L, list(raw(4L))), 1L, 1L, discard = TRUE)$value)
  expect_identical(accept(tx, envelope(header(tx), "terminal"), 1L, tx$control, discard = TRUE)$kind, "terminal")
})

test_that("dynamic next-task ownership follows result order not READY order", {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  accept <- fanout_contract_accept()
  envelope <- getFromNamespace(".npRmpi_fanout_envelope", "npRmpi")
  header <- getFromNamespace(".npRmpi_fanout_header", "npRmpi")
  tx <- make(101L, 2L, 6L, "dynamic", NULL, "session", "operation", 18432L)
  tx$rank[] <- "busy"
  tx$assigned <- list(1L, 2L)
  tx$tags <- 1:2
  tx$scheduled <- 2L
  accept(tx, envelope(header(tx), "result", 1L, list(10L)), 1L, 1L)
  accept(tx, envelope(header(tx), "result", 2L, list(20L)), 2L, 2L)
  expect_identical(tx$next.task, 3:4)
  accept(tx, envelope(header(tx), "ready", 2L), 2L, tx$control)
  accept(tx, envelope(header(tx), "ready", 1L), 1L, tx$control)
  expect_identical(tx$next.task, 3:4)
  expect_identical(tx$rank, c("ready", "ready"))
})

test_that("transaction cleanup preserves errors and interrupts and bounds metadata", {
  run <- getFromNamespace(".npRmpi_fanout_run", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  state <- new.env(parent = emptyenv())
  state$drained <- 0L
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) {
    state$drained <- state$drained + 1L
    tx$rank[] <- "terminal"
    tx$phase <- "quiescent"
    invisible(TRUE)
  }, .package = "npRmpi")
  set.seed(441L)
  seed <- .Random.seed
  for (condition in list(simpleError("original collector failure"),
      structure(list(), class = c("interrupt", "condition")),
      structure(list(message = "original interrupt", call = NULL), class = c("interrupt", "condition")))) {
    tx <- fanout_contract_tx()
    caught <- tryCatch(run(tx, stop(condition)), error = identity, interrupt = identity)
    expect_identical(caught, condition)
    expect_false(exists(tx$key, registry, inherits = FALSE))
    expect_identical(tx$phase, "quiescent")
  }
  expect_identical(state$drained, 3L)
  expect_identical(.Random.seed, seed)
  expect_false(any(c("FUN", "results", "payload") %in% ls(tx)))
})

test_that("unhandled interrupts abort after cleanup rather than becoming errors", {
  run <- getFromNamespace(".npRmpi_fanout_run", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  state <- new.env(parent = emptyenv())
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) {
    state$drained <- state$drained + 1L
    tx$rank[] <- "terminal"
    tx$phase <- "quiescent"
    invisible(TRUE)
  }, .package = "npRmpi")
  condition <- structure(list(), class = c("interrupt", "condition"))
  for (wrapper in c("none", "try", "error", "calling")) {
    state$drained <- 0L
    state$resumed <- FALSE
    state$aborted <- FALSE
    state$unwound <- FALSE
    state$caught <- NULL
    tx <- fanout_contract_tx()
    invoke <- function() {
      on.exit({ state$unwound <- TRUE }, add = TRUE)
      run(tx, signalCondition(condition))
    }
    withRestarts({
      switch(wrapper,
        none = invoke(),
        try = try(invoke(), silent = TRUE),
        error = tryCatch(invoke(), error = function(e) { state$caught <- e; NULL }),
        calling = withCallingHandlers(invoke(), interrupt = function(e) {
          state$caught <- e
        }))
      state$resumed <- TRUE
    }, abort = function() { state$aborted <- TRUE; NULL })
    expect_true(state$aborted, info = wrapper)
    expect_false(state$resumed, info = wrapper)
    expect_true(state$unwound, info = wrapper)
    expect_identical(state$drained, 1L)
    expect_identical(state$caught, if (wrapper == "calling") condition else NULL)
    expect_identical(tx$phase, "quiescent")
    expect_false(exists(tx$key, registry, inherits = FALSE))
  }
})

test_that("a failed drain retains the original message-less interrupt under warn two", {
  run <- getFromNamespace(".npRmpi_fanout_run", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  tx <- fanout_contract_tx()
  withr::defer(rm(list = tx$key, envir = registry))
  withr::local_options(warn = 2L)
  state <- new.env(parent = emptyenv())
  state$drained <- 0L
  state$poisoned <- FALSE
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) {
    state$drained <- state$drained + 1L
    stop("transport failure")
  }, .npRmpi_lease_poison = function() { state$poisoned <- TRUE; NULL },
  .package = "npRmpi")
  condition <- structure(list(), class = c("interrupt", "condition"))
  expect_identical(tryCatch(run(tx, signalCondition(condition)),
                            interrupt = identity), condition)
  expect_identical(state$drained, 1L)
  expect_true(state$poisoned)
  expect_identical(tx$phase, "uncertain")
  expect_true(exists(tx$key, registry, inherits = FALSE))
})

test_that("uncertain cleanup cannot replace the first error under warn equals two", {
  run <- getFromNamespace(".npRmpi_fanout_run", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  tx <- fanout_contract_tx()
  withr::defer(rm(list = tx$key, envir = registry))
  withr::local_options(warn = 2L)
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) stop("transport failure"),
                        .npRmpi_lease_poison = function() invisible(NULL), .package = "npRmpi")
  condition <- simpleError("original failure")
  expect_identical(tryCatch(run(tx, stop(condition)), error = identity), condition)
  expect_identical(tx$phase, "uncertain")
  expect_true(exists(tx$key, registry, inherits = FALSE))
})

test_that("scatter preparation retains classed padding and forwarded RNG order", {
  apply <- getFromNamespace(".npRmpi_fanout_apply", "npRmpi")
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(mpi.comm.size = function(comm) 4L,
    .npRmpi_fanout_new = function(...) fanout_contract_tx("scatter"),
    .npRmpi_scatter_prepare = function(obj, comm) { captured$obj <- obj; NULL },
    .npRmpi_fanout_run = function(tx, code) { captured$tag <- tx$tags[[1L]]; NULL },
    .package = "npRmpi")
  x <- factor(c("a", "b"))
  apply(x, identity, list(), 1L)
  expect_identical(captured$obj, c(list("master"), as.list(c(x, as.list(0L)))))
  x <- as.Date(c("2026-09-16", "2026-09-17"))
  old <- tryCatch(c(x, as.list(0L)), error = identity)
  expect_s3_class(old, "error")
  expect_error(apply(x, identity, list(), 1L), conditionMessage(old), fixed = TRUE)
  set.seed(445L)
  forwarded <- runif(1L)
  tag <- as.integer(floor(runif(1L, 1, 1000)))
  seed <- .Random.seed
  set.seed(445L)
  apply(1:2, identity, list(z = runif(1L)), 1L)
  expect_identical(captured$tag, tag)
  expect_identical(.Random.seed, seed)
})

test_that("a repeated interrupt escapes recovery without consuming the receipt", {
  receive <- getFromNamespace(".npRmpi_fanout_receive", "npRmpi")
  # Local transport doubles do not alter the namespace or any global binding.
  env <- new.env(parent = environment(receive))
  environment(receive) <- env
  env$polls <- 0L
  env$sleeps <- 0L
  env$mpi.any.source <- env$mpi.any.tag <- function() -1L
  env$mpi.iprobe <- function(...) { env$polls <- env$polls + 1L; env$polls > 1L }
  env$Sys.sleep <- function(...) {
    env$sleeps <- env$sleeps + 1L
    stop(structure(list(message = "second interrupt", call = NULL),
                   class = c("interrupt", "condition")))
  }
  tx <- fanout_contract_tx()
  tx$rank[] <- "stopping"
  tx$assigned[[1L]] <- integer()
  env$mpi.get.sourcetag <- function() c(1L, tx$control)
  env$.npRmpi_recv_raw_probed <- function(...) serialize(list(session = tx$session,
    operation = tx$id, kind = "terminal", tasks = integer(), value = NULL), NULL)
  env$.npRmpi_fanout_accept <- fanout_contract_accept()
  interrupted <- tryCatch(receive(tx, discard = TRUE), interrupt = identity)
  expect_s3_class(interrupted, "interrupt")
  expect_identical(tx$rank, "stopping")
  expect_identical(tx$raw.messages, 0L)
  expect_identical(receive(tx, discard = TRUE)$kind, "terminal")
  expect_identical(tx$rank, "terminal")
  expect_identical(env$sleeps, 1L)
  expect_identical(env$polls, 2L)
  expect_identical(tx$raw.messages, 1L)
})

test_that("recovery descriptors are restored across nested normal and error exits", {
  scope <- getFromNamespace(".npRmpi_with_fanout_recovery", "npRmpi")
  current <- getFromNamespace(".npRmpi_fanout_recovery_owner", "npRmpi")
  a <- list(ref = new.env(parent = emptyenv()), slot = "state", comm = 1L)
  b <- list(ref = new.env(parent = emptyenv()), slot = "value", comm = 1L)
  expect_null(current(1L))
  scope(a, {
    expect_identical(current(1L), a)
    expect_null(current(2L))
    expect_error(scope(b, { expect_identical(current(1L), b); stop("nested") }), "nested")
    expect_identical(current(1L), a)
    expect_identical(scope(b, current(1L)), b)
  })
  expect_null(current(1L))
})

test_that("quarantined work requires explicit cleanup and cannot reenter close", {
  quiesce <- getFromNamespace(".npRmpi_fanout_quiesce", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  tx <- fanout_contract_tx()
  tx$phase <- "quarantined"
  assign(tx$key, tx, registry)
  withr::defer(if (exists(tx$key, registry, inherits = FALSE)) rm(list = tx$key, envir = registry))
  tx$owner.active <- TRUE
  expect_error(quiesce(tx$comm, resume = TRUE), "active fan-out")
  tx$owner.active <- FALSE
  expect_error(quiesce(tx$comm), "implicit cleanup")
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) {
    tx$phase <- "quiescent"
    tx$rank[] <- "closed"
    invisible(TRUE)
  }, .package = "npRmpi")
  expect_true(quiesce(tx$comm, resume = TRUE))
  expect_false(exists(tx$key, registry, inherits = FALSE))
})

test_that("a second interrupt quarantines once without an on-exit retry", {
  run <- getFromNamespace(".npRmpi_fanout_run", "npRmpi")
  registry <- getFromNamespace(".npRmpi_fanout_state", "npRmpi")
  tx <- fanout_contract_tx()
  withr::defer(rm(list = tx$key, envir = registry))
  state <- new.env(parent = emptyenv())
  state$count <- 0L
  interrupt <- structure(list(), class = c("interrupt", "condition"))
  local_mocked_bindings(.npRmpi_fanout_drain = function(tx, recovery = NULL) {
    state$count <- state$count + 1L
    signalCondition(interrupt)
  }, .npRmpi_fanout_notice = function(...) NULL, .package = "npRmpi")
  original <- simpleError("original")
  expect_identical(tryCatch(run(tx, stop(original)), error = identity), original)
  expect_identical(state$count, 1L)
  expect_identical(tx$phase, "quarantined")
  expect_false(tx$owner.active)
  expect_true(exists(tx$key, registry, inherits = FALSE))
})

test_that("recovery reuses only its visible owner and never advances estimator counts", {
  step <- getFromNamespace(".npRmpi_fanout_recovery_step", "npRmpi")
  registry <- getFromNamespace(".np_progress_registry", "npRmpi")
  previous <- registry$active_id
  withr::defer({ registry$active_id <- previous })
  withr::local_options(np.messages = TRUE)
  ref <- new.env(parent = emptyenv())
  ref$state <- list(id = "fanout-test-owner", enabled = TRUE, visible = TRUE,
    known_total = TRUE, total = 10L, last_done = 3L, last_emitted_done = 3L,
    unknown_total_fields = function(...) stop("must not run"),
    renderer = "legacy", throttle_sec = 0.2, message_muffled = FALSE)
  registry$active_id <- ref$state$id
  owner <- list(ref = ref, slot = "state", comm = 1L)
  calls <- new.env(parent = emptyenv())
  calls$force <- logical()
  local_mocked_bindings(.np_progress_now = function() 2,
    .np_progress_step_at = function(state, now, force = FALSE, ...) {
      calls$force <- c(calls$force, force)
      state
    }, .package = "npRmpi")
  step(owner, first = TRUE)
  step(owner)
  expect_identical(calls$force, c(TRUE, FALSE))
  expect_identical(ref$state$id, "fanout-test-owner")
  expect_identical(registry$active_id, ref$state$id)
  expect_identical(ref$state$label, "Waiting for MPI cleanup")
  expect_false(ref$state$known_total)
  expect_null(ref$state$total)
  expect_null(ref$state$last_done)
  expect_null(ref$state$unknown_total_fields)
  expect_identical(ref$state$renderer, "legacy")
  expect_identical(ref$state$throttle_sec, 0.2)
  options(np.messages = FALSE)
  step(owner)
  expect_length(calls$force, 2L)
  options(np.messages = TRUE)
  registry$active_id <- "different-owner"
  step(owner)
  expect_length(calls$force, 2L)
  registry$active_id <- ref$state$id
  ref$state$message_muffled <- TRUE
  step(owner)
  expect_length(calls$force, 2L)
})

test_that("recovery renderer failure disables display without losing its cleanup owner", {
  step <- getFromNamespace(".npRmpi_fanout_recovery_step", "npRmpi")
  abort <- getFromNamespace(".np_progress_abort", "npRmpi")
  registry <- getFromNamespace(".np_progress_registry", "npRmpi")
  previous <- registry$active_id
  withr::defer({ registry$active_id <- previous })
  withr::local_options(np.messages = TRUE)
  ref <- new.env(parent = emptyenv())
  ref$state <- list(id = "failed-render-owner", enabled = TRUE, visible = TRUE)
  registry$active_id <- ref$state$id
  local_mocked_bindings(.np_progress_now = function() 1,
    .np_progress_step_at = function(...) stop("display failed"), .package = "npRmpi")
  owner <- list(ref = ref, slot = "state", comm = 1L)
  expect_silent(step(owner, first = TRUE))
  expect_false(ref$state$enabled)
  expect_identical(registry$active_id, ref$state$id)
  expect_silent(step(owner))
  expect_silent(abort(ref$state))
  expect_null(registry$active_id)
})

test_that("every private worker requires one exact transaction header", {
  validate <- getFromNamespace(".npRmpi_fanout_worker_header", "npRmpi")
  header <- getFromNamespace(".npRmpi_fanout_header", "npRmpi")
  for (scheduler in c("scatter", "dynamic", "bundle")) {
    tx <- fanout_contract_tx(scheduler)
    tx$comm <- 1L
    shared <- list(FUN = identity, dot.arg = list(), transaction = header(tx))
    expect_identical(validate(shared, scheduler, 1L), shared$transaction)
    expect_error(validate(list(FUN = identity, dot.arg = list()), scheduler, 1L), "transaction header")
    wrong <- shared
    names(wrong)[[3L]] <- "transaction_extra"
    expect_error(validate(wrong, scheduler, 1L), "transaction header")
    wrong <- shared
    wrong$transaction$session <- NA_character_
    expect_error(validate(wrong, scheduler, 1L), "transaction header")
    wrong <- shared
    wrong$transaction <- c(wrong$transaction, list(session = "duplicate"))
    expect_error(validate(wrong, scheduler, 1L), "transaction header")
  }
})
