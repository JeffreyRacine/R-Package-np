# Real R interrupt, deliberately without an interrupt handler in the caller.
suppressPackageStartupMessages(library(spec$pkg, character.only = TRUE))
ns <- asNamespace(spec$pkg)
dll <- dyn.load(spec$dll)
intr <- getNativeSymbolInfo("probe_interrupt", dll)
options(np.messages = TRUE, np.tree = FALSE, np.progress.interval.sec = 0)
assignInNamespace(".np_progress_is_interactive", function() TRUE, spec$pkg)
state <- new.env()
state$fired <- FALSE
inject <- function() {
  if (!any(vapply(sys.frames(), function(e)
    isTRUE(get0("propagate", envir = e, inherits = FALSE)), logical(1))))
    return(invisible(NULL))
  state$fired <- TRUE
  cat("REAL_OBSERVER_INTERRUPT\n", file = stderr())
  .Call(intr)
}
for (hook in c(".np_progress_nomad_native_step_from_c",
               ".np_progress_bandwidth_activity_step"))
  suppressMessages(trace(hook, where = ns, tracer = quote(inject()), print = FALSE))
set.seed(3)
dat <- data.frame(x = runif(60), z = runif(60))
dat$y <- sin(3 * dat$x) + dat$z + rnorm(60, sd = .3)
search <- function() {
  on.exit(cat("REAL_OBSERVER_CALLER_UNWOUND\n"), add = TRUE)
  body <- function() get("npregbw", ns)(y ~ x + z, data = dat,
    regtype = "lp", degree = c(1L, 1L), bwsolver = "mads", nmulti = 1L,
    nomad.opts = list(MAX_BB_EVAL = 12L))
  if (spec$pkg == "npRmpi") get(".npRmpi_with_local_regression", ns)(body())
  else body()
}
switch(spec$wrap,
  top = search(),
  try = try(search(), silent = TRUE),
  error = tryCatch(search(), error = identity))
cat("REAL_OBSERVER_INCORRECTLY_CONTINUED\n")
