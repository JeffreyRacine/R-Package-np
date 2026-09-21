test_that("fanout deadlines default to finite positive completion-only budgets", {
  budget <- getFromNamespace(".npRmpi_fanout_budget", "npRmpi")
  old <- options(npRmpi.fanout.terminal.timeout = NULL,
                 npRmpi.fanout.cleanup.timeout = NULL,
                 npRmpi.session.recv.timeout = .01)
  on.exit(options(old), add = TRUE)
  expect_identical(budget("terminal"), 30)
  expect_identical(budget("cleanup"), 30)
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2), "30")) {
    options(npRmpi.fanout.terminal.timeout = bad)
    expect_error(budget("terminal"), "finite positive")
  }
})

test_that("numerical completion accounting never clocks incomplete work", {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  credit <- getFromNamespace(".npRmpi_fanout_result_complete", "npRmpi")
  tx <- make(101L, 3L, 4L, "bundle", NULL, "session", "operation", 18432L)
  state <- new.env(parent = emptyenv()); state$calls <- 0L
  env <- new.env(parent = environment(credit))
  env$proc.time <- function() {
    state$calls <- state$calls + 1L
    c(elapsed = 100)
  }
  environment(credit) <- env
  credit(tx, 2L); credit(tx, 1L)
  expect_identical(state$calls, 0L)
  expect_true(is.na(tx$terminal.started))
  credit(tx, 1L)
  expect_identical(state$calls, 1L)
  expect_identical(tx$terminal.started, 100)
  expect_identical(tx$terminal.deadline, rep(Inf, 3L))
  credit(tx, 0L)
  expect_identical(state$calls, 1L)
})

test_that("terminal expiry quarantines without stacking a cleanup budget", {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  boundary <- getFromNamespace(".npRmpi_fanout_terminal_boundary", "npRmpi")
  cleanup <- getFromNamespace(".npRmpi_fanout_cleanup_attempt", "npRmpi")
  tx <- make(101L, 1L, 1L, "scatter", NULL, "session", "operation", 18432L)
  tx$phase <- "active"; tx$owner.active <- TRUE
  tx$terminal.started <- -100
  tx$terminal.deadline[] <- -70
  expect_error(boundary(tx), class = "npRmpi_cleanup_pending")
  local_mocked_bindings(.npRmpi_fanout_notice = function(...) NULL,
    .npRmpi_fanout_drain = function(...) stop("must not drain"), .package = "npRmpi")
  expect_s3_class(cleanup(tx), "npRmpi_cleanup_pending")
  expect_identical(tx$phase, "quarantined")
  tx$owner.active <- FALSE
  expect_match(conditionMessage(cleanup(tx)), "must not drain")
})

test_that("terminal stage budgets are per rank and do not clock numerical work", {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  stage <- getFromNamespace(".npRmpi_fanout_terminal_stage", "npRmpi")
  boundary <- getFromNamespace(".npRmpi_fanout_terminal_boundary", "npRmpi")
  tx <- make(101L, 2L, 2L, "scatter", NULL, "session", "operation", 18432L)
  env <- new.env(parent = environment(stage))
  env$proc.time <- function() stop("no clock during work")
  environment(stage) <- env
  expect_silent(stage(tx, 1L))
  tx$terminal.started <- 100
  tx$terminal.deadline[] <- 130
  tx$rank[] <- "stopping"
  env$proc.time <- function() c(elapsed = 140)
  environment(boundary) <- env
  stage(tx, 1L)
  expect_identical(tx$terminal.deadline, c(Inf, 130))
  expect_error(boundary(tx), class = "npRmpi_cleanup_pending")
  expect_identical(tx$terminal.deadline, c(170, 130))
  tx$rank[[2L]] <- "closed"
  expect_silent(boundary(tx))
  tx$finish.started <- 105
  expect_error(boundary(tx), class = "npRmpi_cleanup_pending")
})

test_that("queued terminal receipts are consumed before idle deadline checks", {
  make <- getFromNamespace(".npRmpi_fanout_metadata", "npRmpi")
  receive <- getFromNamespace(".npRmpi_fanout_receive", "npRmpi")
  tx <- make(101L, 1L, 1L, "scatter", NULL, "session", "operation", 18432L)
  tx$terminal.started <- -100
  tx$terminal.deadline[] <- -70
  tx$rank[] <- "busy"
  env <- new.env(parent = environment(receive))
  env$mpi.iprobe <- function(...) TRUE
  env$mpi.get.sourcetag <- function() c(1L, 18432L)
  env$.npRmpi_recv_raw_probed <- function(...) serialize("queued", NULL)
  env$.npRmpi_fanout_accept <- function(tx, message, ...) message
  env$.npRmpi_fanout_terminal_boundary <- function(...) stop("expired before probe")
  environment(receive) <- env
  expect_identical(receive(tx, poll = FALSE), "queued")
})

test_that("new terminal polling is lean without changing explicit work polling", {
  receive <- getFromNamespace(".npRmpi_fanout_receive", "npRmpi")
  env <- new.env(parent = environment(receive))
  state <- new.env(parent = emptyenv())
  env$.npRmpi_fanout_terminal_boundary <- function(...) NULL
  env$mpi.iprobe <- function(...) {
    state$probes <- state$probes + 1L
    if (state$probes == 2L) stop("probe sentinel")
    FALSE
  }
  env$mpi.probe <- function(...) stop("blocking sentinel")
  env$Sys.sleep <- function(seconds) { state$slept <- seconds }
  environment(receive) <- env
  for (complete in c(FALSE, TRUE)) for (poll in c(FALSE, TRUE)) {
    state$probes <- 0L; state$slept <- NULL
    tx <- list(comm = 101L, terminal.started = if (complete) 1 else NA_real_)
    expect_error(receive(tx, poll = poll, sleep = .25), "sentinel")
    expect_identical(state$slept,
      if (complete && !poll) .0005 else if (poll) .25 else NULL)
  }
})

test_that("real terminal deadlines spare running tasks and retain recoverable pools", {
  skip_on_cran()
  skip_if(Sys.getenv("NP_RMPI_RUN_FANOUT_SIGINT_TESTS") != "TRUE",
          "opt-in terminal-deadline subprocess")
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "installed npRmpi unavailable")
  for (workers in c(1L, 3L)) {
    result <- npRmpi_run_rscript_subprocess(c(
      "suppressPackageStartupMessages(library(npRmpi))",
      sprintf("npRmpi.init(nslaves=%d, quiet=TRUE)", workers),
      "options(npRmpi.fanout.terminal.timeout=.15, npRmpi.fanout.cleanup.timeout=.15)",
      "ns <- asNamespace('npRmpi')",
      "apply <- get('.npRmpi_fanout_apply', ns)",
      "retained <- get('.npRmpi_fanout_retained', ns)",
      sprintf("x <- seq_len(%dL)", workers),
      "for (poll in c(FALSE, TRUE)) {",
      "  cat('RUNNING_TASK',poll,'\\n'); flush.console()",
      "  v <- apply(x, function(z) { Sys.sleep(.3); z }, list(), 1L, poll=poll)",
      "  stopifnot(identical(v, as.list(x)), is.null(retained(1L)))",
      "}",
      "credit <- get('.npRmpi_fanout_result_complete', ns)",
      "delayed <- function(tx,count) {credit(tx,count); if(tx$completed==tx$n) Sys.sleep(.4)}",
      "unlockBinding('.npRmpi_fanout_result_complete',ns)",
      "assign('.npRmpi_fanout_result_complete',delayed,ns)",
      "lockBinding('.npRmpi_fanout_result_complete',ns)",
      "for(dynamic in c(FALSE,TRUE)) for(poll in c(FALSE,TRUE)) {",
      "  cat('QUEUED_RECEIPTS',dynamic,poll,'\\n'); flush.console()",
      "  z <- if(dynamic) seq_len(2L*length(x)+1L) else x",
      "  v <- apply(z,identity,list(),1L,dynamic=dynamic,poll=poll)",
      "  stopifnot(identical(v,as.list(z)),is.null(retained(1L)))",
      "}",
      "unlockBinding('.npRmpi_fanout_result_complete',ns)",
      "assign('.npRmpi_fanout_result_complete',credit,ns)",
      "lockBinding('.npRmpi_fanout_result_complete',ns)",
      "stopifnot(identical(apply(x,identity,list(),1L),as.list(x)))",
      "cat('MISSING_TERMINAL\\n'); flush.console()",
      "mpi.bcast.cmd(local({",
      "  ns <- asNamespace('npRmpi')",
      "  original <- get('.npRmpi_fanout_worker_terminal', ns)",
      "  delayed <- function(...) { Sys.sleep(1.5); original(...) }",
      "  unlockBinding('.npRmpi_fanout_worker_terminal', ns)",
      "  assign('.npRmpi_fanout_worker_terminal', delayed, ns)",
      "  lockBinding('.npRmpi_fanout_worker_terminal', ns)",
      "}))",
      "started <- proc.time()[['elapsed']]",
      "e <- tryCatch(apply(x, identity, list(), 1L), error=identity)",
      "stopifnot(inherits(e, 'npRmpi_cleanup_pending'),",
      "  proc.time()[['elapsed']] - started < 1,",
      "  retained(1L)$completed == length(x), retained(1L)$phase == 'quarantined')",
      "e <- tryCatch(apply(x, identity, list(), 1L), error=identity)",
      "stopifnot(inherits(e, 'error'))",
      "started <- proc.time()[['elapsed']]",
      "e <- tryCatch(npRmpi.quit(force=TRUE), error=identity)",
      "stopifnot(inherits(e, 'npRmpi_cleanup_pending'),",
      "  proc.time()[['elapsed']] - started < 1)",
      "options(npRmpi.fanout.cleanup.timeout=4)",
      "npRmpi.quit(force=TRUE)",
      "stopifnot(is.null(retained(1L)))",
      "cat('TERMINAL_DEADLINE_RAW_ZERO\\n')"
    ), timeout = 20L, env = env, cleanup = FALSE)
    expect_identical(result$status, 0L,
                     info = paste(result$output, collapse = "\n"))
    expect_true(any(grepl("TERMINAL_DEADLINE_RAW_ZERO", result$output, fixed = TRUE)),
                info = paste(result$output, collapse = "\n"))
  }
})
