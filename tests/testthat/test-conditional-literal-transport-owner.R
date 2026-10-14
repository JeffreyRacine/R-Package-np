# The conditional-only dispatch copies retain shared publication policy while
# preventing S3 genericization from replacing the selected literal owner.
test_that("conditional kernel transport preserves owner and publication policy", {
  ns <- asNamespace("npRmpi")
  wrapper <- get("npksum.conditional_literal", ns, inherits = FALSE)
  lex <- get(".np_conditional_literal_holder", ns, inherits = FALSE)$lex
  expect_identical(lex$npksum.default, wrapper)
  expect_identical(formals(wrapper), formals(get("npksum.default", ns)))
  expect_identical(environment(wrapper), ns)
  expect_false("npksum.conditional_literal" %in% getNamespaceExports(ns))

  if (identical(unname(getNamespaceName(ns)), "npRmpi")) {
    closed <- get(".npksum_closed_values", lex, inherits = FALSE)
    shared.closed <- get(".npksum_closed_values", ns, inherits = FALSE)
    expect_identical(environment(closed), lex)
    expect_identical(environment(shared.closed), ns)
    expect_identical(body(closed), body(shared.closed))
    expect_identical(formals(closed), formals(shared.closed))
    expect_true(bindingIsLocked(".npksum_closed_values", lex))
    adapter <- lex$.npRmpi_autodispatch_call
    dispatch <- get("original", envir = environment(adapter), inherits = FALSE)
    transport <- environment(dispatch)
    expect_error(adapter(quote(npksum.default()), caller_env = environment(),
                         owner.name = "npksum.default"),
                 "unexpected conditional kernel-sum dispatch engine")
    # Capture the actual adapter with a test-local, explicitly trusted engine.
    # The real engine/worker execution is a separate focused MPI proof.
    captured.adapter <- adapter
    capture <- new.env(parent = environment(adapter))
    capture$original <- function(mc, caller_env, comm, owner.name)
      list(call = mc, caller = caller_env, comm = comm, owner = owner.name)
    trusted.engine <- function() captured.adapter(
      quote(npksum.default(txdat = owned.data, operator = "normal")),
      caller_env = environment(), comm = 1L, owner.name = "npksum.default")
    capture$engine <- trusted.engine
    environment(captured.adapter) <- capture
    captured <- trusted.engine()
    expect_identical(captured$call[[1L]], as.name("npksum.conditional_literal"))
    expect_true(attr(captured$call, ".npRmpi.bound.call.owner", exact = TRUE))
    expect_identical(as.list(captured$call)[-1L],
                     as.list(quote(npksum.default(txdat = owned.data,
                                                  operator = "normal")))[-1L])
    expect_identical(captured$owner, "npksum.conditional_literal")
    expect_identical(captured$comm, 1L)
    expect_true(environmentIsLocked(transport))
    expect_identical(parent.env(transport), ns)
    for (name in c(".npRmpi_autodispatch_call", ".npRmpi_distributed_call_impl")) {
      original <- get(name, ns, inherits = FALSE)
      copied <- get(name, transport, inherits = FALSE)
      expect_identical(environment(original), ns, info = name)
      expect_identical(environment(copied), transport, info = name)
      expect_identical(body(copied), body(original), info = name)
      expect_identical(formals(copied), formals(original), info = name)
      expect_true(bindingIsLocked(name, transport), info = name)
    }
    shared <- get(".npRmpi_autodispatch_materialize_call", ns)
    private <- get(".npRmpi_autodispatch_materialize_call", transport,
                   inherits = FALSE)
    expect_identical(get(".np_conditional_shared_materialize_call", transport), shared)
    expect_identical(environment(shared), ns)
    expect_identical(environment(private), transport)
    expect_true(bindingIsLocked(".npRmpi_autodispatch_materialize_call", transport))

    # Force multiple ordinary prepared publications without publishing to MPI.
    old <- options(npRmpi.autodispatch.arg.broadcast.threshold = 1L)
    on.exit(options(old), add = TRUE)
    bw.value <- list(type = "fixed")
    training <- data.frame(x = c(.2, .4, .6))
    evaluation <- data.frame(x = c(.3, .5))
    mc <- quote(npksum.conditional_literal(bws = bw.value,
      txdat = training, exdat = evaluation, operator = "normal"))
    bind <- get(".npRmpi_autodispatch_bind_call_owner", ns)
    bound <- bind(mc, "npksum.conditional_literal", mc, wrapper)
    expect_true(attr(bound, ".npRmpi.bound.call.owner", exact = TRUE))
    ordinary <- shared(bound, caller_env = environment())
    conditional <- private(bound, caller_env = environment())
    expected <- ordinary
    expected$call[[1L]] <- as.name("npksum.conditional_literal")
    expect_identical(conditional, expected)
    expect_identical(conditional[setdiff(names(conditional), "call")],
                     ordinary[setdiff(names(ordinary), "call")])
    expect_gt(length(conditional$prepublish), 1L)
    call.name <- get(".npRmpi_autodispatch_call_name", ns)(bound)
    # This is the unchanged distributed body's exact bundle decision.
    expect_identical(sub("\\..*$", "", call.name), "npksum")
    expect_true(identical(sub("\\..*$", "", call.name), "npksum") &&
                  length(conditional$prepublish) > 1L)
    expect_identical(get(".npRmpi_spmd_opcode_from_call", ns)(bound, environment()),
                     "autodispatch.npksum.core")
    # Exercise the unchanged worker guard with an evaluator stub only; a
    # literal npksum allow-list is narrower than its registered core policy.
    guard <- get(".npRmpi_spmd_eval_payload_call_guard", ns)
    guard.capture <- new.env(parent = environment(guard))
    guard.capture$.npRmpi_spmd_eval_payload <- function(payload, envelope)
      list(payload = payload, envelope = envelope)
    environment(guard) <- guard.capture
    envelope <- list(opcode = "autodispatch.npksum.core")
    admitted <- guard(list(call = conditional$call), envelope,
                      allowed_calls = "npksum")
    expect_identical(admitted$payload$call, conditional$call)
    expect_identical(admitted$envelope, envelope)
    expect_error(guard(list(call = quote(npreg())), envelope,
                       allowed_calls = "npksum"), "restricted to")
    expect_identical(get(".npRmpi_lease_publication_plan", ns)(bound),
                     list(kind = "none", producer = "npksum"))
    expect_identical(conditional$lease.bindings, ordinary$lease.bindings)
    expect_identical(get(".npRmpi_autodispatch_large_arg_threshold_for_call", ns)(bound), 1L)
    expect_identical(shared(quote(npksum.default()), environment())$call[[1L]],
                     as.name("npksum"))
    expect_error(private(quote(npksum.default()), environment()),
                 "unexpected conditional kernel-sum materialization owner")
    expect_error(private(mc, environment()),
                 "unexpected conditional kernel-sum materialization owner")
  }
})
