# Isolated child: a poisoned native solver must not endanger the test runner.
suppressPackageStartupMessages(library(spec$pkg, character.only = TRUE))
options(np.messages = FALSE, np.tree = FALSE)
dll <- dyn.load(spec$dll)
arm <- getNativeSymbolInfo("probe_arm", dll)
off <- getNativeSymbolInfo("probe_off", dll)
intr <- getNativeSymbolInfo("probe_interrupt", dll)
run <- function(fun) {
  if (spec$pkg == "npRmpi")
    get(".npRmpi_with_local_regression", asNamespace(spec$pkg))(fun())
  else fun()
}
set.seed(62031)
dat <- data.frame(x = runif(60), z = runif(60))
dat$y <- sin(4 * dat$x) + dat$z + rnorm(60, sd = .2)
state <- new.env()
condition <- structure(list(message = "native lifetime probe", call = NULL,
  token = new.env()), class = c("lifetimeProbe", "error", "condition"))
wrappers <- c(npregbw = "npNomadNativeSearchRegression",
  npcdensbw = "npPreparedObjectiveFixedNativeSearchConditionalDensity")
if (spec$pkg == "npRmpi") wrappers[] <- c(
  "npRmpiPreparedObjectiveNativeSearchRegression",
  "npRmpiPreparedObjectiveFixedNativeSearchConditionalDensity")
for (family in names(wrappers)) {
  fun <- get(family, asNamespace(spec$pkg))
  search <- function() run(function() {
    set.seed(42)
    b <- fun(y ~ x + z, data = dat, regtype = "lp", degree = c(1L, 1L),
      bwsolver = "mads", nmulti = 1L, nomad.opts = list(MAX_BB_EVAL = 12L))
    list(b$bw, b$xbw, b$ybw, b$fval, b$num.feval)
  })
  before <- search()
  for (kind in c("error", "interrupt")) {
    state$fired <- FALSE
    injection <- if (kind == "error") quote({state$fired <- TRUE; gc(); stop(condition)})
      else quote({state$fired <- TRUE; gc(); .Call(intr)})
    suppressMessages(trace(wrappers[[family]], where = asNamespace(spec$pkg),
      print = FALSE, tracer = quote(.Call(arm, injection))))
    caught <- tryCatch(search(), error = identity, interrupt = identity)
    suppressMessages(untrace(wrappers[[family]], where = asNamespace(spec$pkg)))
    .Call(off)
    stopifnot(state$fired, inherits(caught, kind))
    if (kind == "error") stopifnot(identical(caught, condition))
    stopifnot(identical(search(), before), identical(search(), before))
    cat("NATIVE_UNWIND_RECOVERY", family, kind, "\n")
  }
}
dyn.unload(spec$dll)
if (spec$pkg == "npRmpi") get("mpi.finalize", asNamespace(spec$pkg))()
cat("NATIVE_UNWIND_PASS\n")
