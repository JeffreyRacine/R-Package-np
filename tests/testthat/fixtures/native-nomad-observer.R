# Isolated child; namespace hooks and elapsed limits are confined to this process.
suppressPackageStartupMessages(library(spec$pkg, character.only = TRUE))
ns <- asNamespace(spec$pkg)
options(np.messages = FALSE, np.tree = FALSE)
run <- function(fun) {
  if (spec$pkg == "npRmpi") get(".npRmpi_with_local_regression", ns)(fun())
  else fun()
}
set.seed(3)
dat <- data.frame(x = runif(60), z = runif(60))
dat$y <- sin(3 * dat$x) + dat$z + rnorm(60, sd = .3)
assignInNamespace(".np_progress_is_interactive", function() TRUE, spec$pkg)
options(np.progress.interval.sec = 0)
for (family in c("npregbw", "npudensbw", "npudistbw", "npcdensbw", "npcdistbw")) {
  search <- function() run(function() {
    set.seed(42)
    args <- list(formula = y ~ x + z, data = dat, bwsolver = "mads",
      nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 12L))
    if (family == "npregbw")
      args <- c(args, list(regtype = "lp", degree = c(1L, 1L)))
    if (family %in% c("npudensbw", "npudistbw"))
      args <- list(dat = dat[c("x", "z")], bwsolver = "mads", nmulti = 1L,
                   nomad.opts = list(MAX_BB_EVAL = 12L))
    if (family == "npcdistbw")
      args <- c(args, list(bwtype = "generalized_nn", cxkertype = "beta",
        cykertype = "beta", cxkerbound = "range", cykerbound = "range"))
    b <- do.call(get(family, ns), args)
    list(b$bw, b$xbw, b$ybw, b$fval, b$num.feval)
  })
  before <- search()
  for (kind in c("error", "interrupt", "time")) {
    state <- new.env()
    state$fired <- FALSE
    state$after <- 0L
    state$original <- NULL
    condition <- structure(list(message = "observer lifetime probe",
      call = quote(observer_probe()), token = new.env()),
      class = c("observerProbe", "error", "condition"))
    inject <- function() {
      # Do not fire in pre-native progress: this is the contained observer call.
      if (!any(vapply(sys.frames(), function(e)
        isTRUE(get0("propagate", envir = e, inherits = FALSE)), logical(1))))
        return(invisible(NULL))
      if (state$fired) {
        state$after <- state$after + 1L
        stop("secondary observer entry")
      }
      state$fired <- TRUE
      withCallingHandlers({
        gc()
        if (kind == "error") stop(condition)
        if (kind == "interrupt") stop(structure(list(message = "interrupt", call = NULL),
          class = c("interrupt", "condition")))
        setTimeLimit(elapsed = .01, transient = TRUE)
        for (i in seq_len(1e7)) sqrt(i)
        stop("elapsed-time injection did not fire")
      }, error = function(e) state$original <- e,
         interrupt = function(e) state$original <- e)
    }
    hooks <- c(".np_progress_nomad_native_step_from_c",
               ".np_progress_bandwidth_activity_step")
    for (hook in hooks) suppressMessages(trace(hook, where = ns,
      tracer = quote(inject()), print = FALSE))
    options(np.messages = TRUE)
    caught <- tryCatch(search(), error = identity, interrupt = identity,
                       finally = setTimeLimit())
    for (hook in hooks) suppressMessages(untrace(hook, where = ns))
    options(np.messages = FALSE)
    stopifnot(state$fired, inherits(caught, "condition"),
              identical(caught, state$original), state$after == 0L)
    if (kind == "time") stopifnot(grepl("elapsed time limit", conditionMessage(caught)))
    stopifnot(identical(search(), before), identical(search(), before))
    cat("OBSERVER_RECOVERY", family, kind, "\n")
  }
}
if (spec$pkg == "npRmpi") get("mpi.finalize", ns)()
cat("NATIVE_OBSERVER_PASS\n")
