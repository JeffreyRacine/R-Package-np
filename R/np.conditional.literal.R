# Conditional evaluation/replay has a literal categorical X/Y target.  The
# original regression and kernel-sum closures remain in the namespace unchanged.
.np_conditional_literal_holder <- new.env(parent = emptyenv())

.np_conditional_literal_rewrite <- function(expr, old) {
  count <- 0L
  rewrite <- function(expr) {
    if (is.null(expr) || (!is.call(expr) && !is.pairlist(expr)))
      return(expr)
    if (is.call(expr) && length(expr) >= 2L &&
        is.symbol(expr[[1L]]) &&
        as.character(expr[[1L]]) %in% c(".Call", "eval.native") &&
        identical(expr[[2L]], old)) {
      expr[[2L]] <- as.name(paste0(".np_literal_", old))
      count <<- count + 1L
    }
    for (i in seq_along(expr)) {
      # NULL is a pairlist too; assigning it through [[ would delete a node.
      if (is.null(expr[[i]]) || identical(expr[[i]], quote(expr = )))
        next
      if (is.call(expr[[i]]) || is.pairlist(expr[[i]]))
        expr[[i]] <- rewrite(expr[[i]])
    }
    expr
  }
  result <- rewrite(expr)
  list(body = result, count = count)
}

.np_conditional_literal_initialize <- function(pkg) {
  holder <- .np_conditional_literal_holder
  if (environmentIsLocked(holder))
    return(invisible(NULL))
  ns <- asNamespace(pkg)
  lex <- new.env(parent = ns)
  names <- c(
    ".np_regression_direct", ".npreghat_exact_matrix_from_core",
    ".npreghat_exact_lp_apply_from_regression_core",
    ".npreghat_exact_lp_matrix_from_regression_core",
    ".npreghat_exact_lp_matrix_from_kernel_weights", ".np_kernel_weights_direct",
    ".npreghat_exact_lc_matrix_from_kernel_weights",
    ".npreghat_exact_lc_derivative_matrix_from_npksum_chunked",
    ".np_local_operator_ksum", ".np_direct_operator_matrix",
    ".np_direct_operator_apply", ".np_exact_operator_matrix",
    ".np_exact_operator_apply", ".npcdhat_make_kernel_matrix",
    ".npcdhat_make_xhat_matrix", "npksum.default",
    ".np_plot_density_derivative_matrix")
  # MPI uses this closed-value forwarding helper where serial evaluates the
  # named kernel call directly.  It must resolve the same private wrapper.
  if (identical(pkg, "npRmpi"))
    names <- c(names, ".npksum_closed_values")
  substitutions <- c(
    .np_regression_direct = "C_np_regression",
    .npreghat_exact_lp_apply_from_regression_core = "C_np_regression_lp_apply_conditional_ctx",
    .np_kernel_weights_direct = "C_np_kernelsum",
    .np_local_operator_ksum = "C_np_kernelsum",
    npksum.default = "C_np_kernelsum")
  entries <- c(C_np_regression = "C_np_regression_conditional_literal",
    C_np_regression_lp_apply_conditional_ctx = "C_np_regression_lp_apply_conditional_ctx_literal",
    C_np_kernelsum = "C_np_kernelsum_conditional_literal")
  for (old in names(entries))
    assign(paste0(".np_literal_", old),
      getNativeSymbolInfo(entries[[old]], PACKAGE = pkg), envir = lex)
  for (name in names) {
    fun <- get(name, envir = ns, inherits = FALSE)
    if (!is.function(fun) || !identical(environment(fun), ns))
      stop("conditional literal closure source contract changed: ", name)
    environment(fun) <- lex
    if (name %in% names(substitutions)) {
      old <- substitutions[[name]]
      rewritten <- .np_conditional_literal_rewrite(body(fun), old)
      if (rewritten$count != 1L)
        stop("conditional literal native-call contract changed: ", name)
      body(fun) <- rewritten$body
      fun <- compiler::cmpfun(fun)
    }
    assign(name, fun, envir = lex)
  }
  # The engine stays local.  Quoted calls in the beta derivative helper enter
  # the namespace wrapper, which also provides the MPI transport owner frame.
  lex$.np_conditional_npksum_engine <- lex$npksum.default
  lex$npksum.default <- get("npksum.conditional_literal", envir = ns,
                            inherits = FALSE)
  if (identical(pkg, "npRmpi")) {
    # Preserve the ordinary dotted-call classification (including publication
    # bundling), but keep its selected private owner after materialization.
    # Both unchanged dispatch closures must resolve the private materializer.
    transport <- new.env(parent = ns)
    for (name in c(".npRmpi_autodispatch_call", ".npRmpi_distributed_call_impl")) {
      fun <- get(name, envir = ns, inherits = FALSE)
      if (!is.function(fun) || !identical(environment(fun), ns))
        stop("conditional literal transport source contract changed: ", name)
      environment(fun) <- transport
      assign(name, fun, envir = transport)
    }
    transport$.np_conditional_shared_materialize_call <- get(
      ".npRmpi_autodispatch_materialize_call", envir = ns, inherits = FALSE)
    materialize <- function(mc, caller_env, comm = 1L) {
      owner <- as.name("npksum.conditional_literal")
      if (!is.call(mc) || !identical(mc[[1L]], owner) ||
          !isTRUE(attr(mc, ".npRmpi.bound.call.owner", exact = TRUE)))
        stop("unexpected conditional kernel-sum materialization owner")
      prepared <- .np_conditional_shared_materialize_call(
        mc, caller_env = caller_env, comm = comm)
      if (!is.call(prepared$call) ||
          !identical(prepared$call[[1L]], as.name("npksum")) ||
          !isTRUE(attr(prepared$call, ".npRmpi.bound.call.owner", exact = TRUE)))
        stop("conditional kernel-sum materialization contract changed")
      prepared$call[[1L]] <- owner
      prepared
    }
    environment(materialize) <- transport
    transport$.npRmpi_autodispatch_materialize_call <- materialize
    lockEnvironment(transport, bindings = TRUE)
    dispatch <- get(".npRmpi_autodispatch_call", envir = transport,
                    inherits = FALSE)
    lex$.npRmpi_autodispatch_call <- local({
      original <- dispatch
      engine <- lex$.np_conditional_npksum_engine
      function(mc, caller_env = parent.frame(), comm = 1L, owner.name = NULL) {
        if (!identical(owner.name, "npksum.default"))
          stop("unexpected conditional kernel-sum dispatch owner")
        caller <- sys.parent()
        if (!is.call(mc) || caller < 1L ||
            !identical(sys.function(caller), engine))
          stop("unexpected conditional kernel-sum dispatch engine")
        # The wrapper evaluates in its caller's frame, so it need not occur
        # on the engine's dynamic parent chain.  Bind only this trusted engine.
        mc[[1L]] <- as.name("npksum.conditional_literal")
        attr(mc, ".npRmpi.bound.call.owner") <- TRUE
        original(mc, caller_env = caller_env, comm = comm,
                 owner.name = "npksum.conditional_literal")
      }
    })
  }
  lockEnvironment(lex, bindings = TRUE)
  holder$lex <- lex
  lockEnvironment(holder, bindings = TRUE)
  lockBinding(".np_conditional_literal_holder", ns)
  invisible(NULL)
}

# Full formals preserve missingness and the existing MPI promise binder.  No
# NativeSymbolInfo or private lexical environment is sent to another rank.
npksum.conditional_literal <- function(bws,
  txdat = stop("training data 'txdat' missing"), tydat = NULL, exdat = NULL,
  weights = NULL, bandwidth.divide = FALSE, compute.ocg = FALSE,
  compute.score = FALSE, kernel.pow = 1.0, leave.one.out = FALSE,
  operator = names(ALL_OPERATORS),
  permutation.operator = names(PERMUTATION_OPERATORS),
  return.kernel.weights = FALSE, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- .np_conditional_literal_holder$lex$.np_conditional_npksum_engine
  eval(mc, envir = parent.frame())
}

.np_conditional_literal_kernel_matrix <- function(...)
  .np_conditional_literal_holder$lex$.npcdhat_make_kernel_matrix(...)

.np_conditional_literal_xhat_matrix <- function(...)
  .np_conditional_literal_holder$lex$.npcdhat_make_xhat_matrix(...)

.np_conditional_literal_kernel_weights <- function(...)
  .np_conditional_literal_holder$lex$.np_kernel_weights_direct(...)

.np_conditional_literal_derivative_matrix <- function(...)
  .np_conditional_literal_holder$lex$.np_plot_density_derivative_matrix(...)

.np_conditional_literal_density_matrix <- function(kbw, x, ex, operator) {
  .np_conditional_literal_kernel_matrix(kbw, x, ex, operator = operator)
}

.np_conditional_policy_stamp <- function(object) {
  object[["conditional.kernel.policy"]] <- "literal-xy-v1"
  object
}

.np_conditional_require_replay_policy <- function(object, where) {
  if (!identical(object[["conditional.kernel.policy", exact = TRUE]],
                 "literal-xy-v1"))
    stop(where, ": stored conditional values predate the literal categorical policy; refit the object before reconstructing uncertainty", call. = FALSE)
  invisible(NULL)
}
