local({
  if (!spawn_mpi_slaves(1L)) stop("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  ns <- asNamespace("npRmpi")
  old.opts <- options(np.messages = FALSE, np.tree = FALSE,
                      np.progress.interval.sec = 0)
  on.exit(options(old.opts), add = TRUE)
  replace <- function(name, value) {
    unlockBinding(name, ns)
    assign(name, value, ns)
    lockBinding(name, ns)
  }
  old.interactive <- get(".np_progress_is_interactive", ns)
  old.activity <- get(".np_progress_bandwidth_activity_step", ns)
  on.exit({
    replace(".np_progress_is_interactive", old.interactive)
    replace(".np_progress_bandwidth_activity_step", old.activity)
  }, add = TRUE)
  replace(".np_progress_is_interactive", function() TRUE)
  dispatch <- get(".np_progress_nomad_native_observer_dispatch", ns)
  state <- new.env(parent = emptyenv())
  condition <- structure(list(message = "pooled observer probe", token = 17L),
    class = c("observerProbe", "error", "condition"))
  injection <- function(...) {
    frames <- which(vapply(seq_len(sys.nframe()), function(i)
      identical(sys.function(i), dispatch), logical(1)))
    if (!state$fired && length(frames)) {
      state$fired <- TRUE
      state$propagate <- get("propagate", sys.frame(tail(frames, 1L)))
      state$broadcast <- getOption("npRmpi.manual.bcast.context", FALSE)
      stop(condition)
    }
    old.activity(...)
  }
  set.seed(3)
  dat <- data.frame(x = runif(60), z = runif(60))
  dat$y <- sin(3 * dat$x) + dat$z + rnorm(60, sd = .3)
  search <- function() {
    set.seed(42)
    b <- npregbw(y ~ x + z, data = dat, regtype = "lp",
      degree = c(1L, 1L), bwsolver = "mads", nmulti = 1L,
      nomad.opts = list(MAX_BB_EVAL = 12L))
    c(b$bw, b$fval, b$num.feval)
  }
  for (local.mode in c(FALSE, TRUE)) {
    run <- function() if (local.mode)
      get(".npRmpi_with_local_regression", ns)(search()) else search()
    cold <- run()
    state$fired <- FALSE
    replace(".np_progress_bandwidth_activity_step", injection)
    options(np.messages = TRUE)
    observed <- tryCatch(suppressMessages(run()), error = identity)
    options(np.messages = FALSE)
    replace(".np_progress_bandwidth_activity_step", old.activity)
    if (!state$fired || !identical(state$propagate, local.mode) ||
        !identical(state$broadcast, !local.mode))
      stop("native observer execution scope was not witnessed")
    if (!identical(observed, if (local.mode) condition else cold))
      stop("native observer failure violated its execution-scope contract")
    if (!identical(run(), cold)) stop("search recovery differs from cold")
  }
  cat("NATIVE_POOLED_OBSERVER_PASS\n")
})
