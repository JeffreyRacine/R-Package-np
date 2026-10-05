# Isolated child: a poisoned native solver must not endanger the test runner.
suppressPackageStartupMessages(library(spec$pkg, character.only = TRUE))
options(np.messages = FALSE, np.tree = FALSE)
dll <- dyn.load(spec$dll)
arm <- getNativeSymbolInfo("probe_arm", dll)
off <- getNativeSymbolInfo("probe_off", dll)
intr <- getNativeSymbolInfo("probe_interrupt", dll)
arm_at <- getNativeSymbolInfo("probe_arm_at", dll)
hits <- getNativeSymbolInfo("probe_native_hits", dll)
# nm/backtrace provide a stronger, test-only owner witness on unstripped
# Darwin/GNU builds. If unavailable, retain the unrelated lifetime cases and
# report this specific coverage gap explicitly; never call it native proof.
poll_range <- NULL
if (nzchar(Sys.which("nm"))) {
  package_dll <- getLoadedDLLs()[[spec$pkg]][["path"]]
  output <- suppressWarnings(system2("nm", c("-n", shQuote(package_dll)),
                                     stdout = TRUE, stderr = TRUE))
  if (is.null(attr(output, "status"))) {
    symbols <- strsplit(trimws(output), "[[:space:]]+")
    symbols <- Filter(function(x) length(x) == 3L && x[[2L]] %in% c("t", "T"), symbols)
    symbol_names <- vapply(symbols, `[`, "", 3L)
    owner <- which(sub("^_", "", symbol_names) ==
      "np_beta_conditional_distribution_bw_objective_ls_ctx")
    if (length(owner) == 1L && owner < length(symbols)) {
      hex <- function(x) sum(strtoi(strsplit(tolower(x), "")[[1L]], 16L) *
                                 16 ^ rev(seq_len(nchar(x)) - 1L))
      candidate <- vapply(symbols[c(owner, owner + 1L)], function(x) hex(x[[1L]]), 0)
      if (all(is.finite(candidate)) && candidate[[2L]] > candidate[[1L]])
        poll_range <- candidate
    }
  }
}
if (!is.null(poll_range)) {
  anchor <- getNativeSymbolInfo("C_np_distribution_conditional_bw", spec$pkg)$address
} else {
  cat("CF234_NATIVE_OWNER_SKIPPED: nm or unstripped beta CV symbol unavailable\n")
}
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
  npcdensbw = "npPreparedObjectiveFixedNativeSearchConditionalDensity",
  npcdistbw = "npNomadNativeSearchConditionalDistribution")
if (spec$pkg == "npRmpi") wrappers[] <- c(
  "npRmpiPreparedObjectiveNativeSearchRegression",
  "npRmpiPreparedObjectiveFixedNativeSearchConditionalDensity",
  "npRmpiNomadNativeSearchConditionalDistribution")
if (is.null(poll_range)) wrappers <- wrappers[names(wrappers) != "npcdistbw"]
for (family in names(wrappers)) {
  fun <- get(family, asNamespace(spec$pkg))
  search <- function() run(function() {
    set.seed(42)
    args <- list(formula = y ~ x + z, data = dat, regtype = "lp",
      degree = c(1L, 1L), bwsolver = "mads", nmulti = 1L,
      nomad.opts = list(MAX_BB_EVAL = 12L))
    if (family == "npcdistbw")
      args <- c(args, list(bwtype = "generalized_nn", cxkertype = "beta",
        cykertype = "beta", cxkerbound = "range", cykerbound = "range"))
    b <- do.call(fun, args)
    list(b$bw, b$xbw, b$ybw, b$fval, b$num.feval)
  })
  before <- search()
  for (kind in c("error", "interrupt")) {
    state$fired <- FALSE
    injection <- if (kind == "error") quote({state$fired <- TRUE; gc(); stop(condition)})
      else quote({state$fired <- TRUE; gc(); .Call(intr)})
    tracer <- if (family == "npcdistbw")
      quote(.Call(arm_at, injection, anchor, poll_range)) else quote(.Call(arm, injection))
    suppressMessages(trace(wrappers[[family]], where = asNamespace(spec$pkg),
      print = FALSE, tracer = tracer))
    caught <- tryCatch(search(), error = identity, interrupt = identity)
    suppressMessages(untrace(wrappers[[family]], where = asNamespace(spec$pkg)))
    .Call(off)
    if (family == "npcdistbw") {
      stopifnot(identical(.Call(hits), 3L))
      cat("CF234_NATIVE_OWNER_WITNESS", kind, "\n")
    }
    stopifnot(state$fired, inherits(caught, kind))
    if (kind == "error") stopifnot(identical(caught, condition))
    stopifnot(identical(search(), before), identical(search(), before))
    cat("NATIVE_UNWIND_RECOVERY", family, kind, "\n")
  }
}
dyn.unload(spec$dll)
if (spec$pkg == "npRmpi") get("mpi.finalize", asNamespace(spec$pkg))()
cat("NATIVE_UNWIND_PASS\n")
