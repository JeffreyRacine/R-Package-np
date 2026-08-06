test_that("npplreg child searches use the canonical prepared regression owner", {
  ns <- asNamespace("npRmpi")
  child.search <- get(".npplregbw_child_specific_nomad_search", ns,
                      inherits = FALSE)
  child.args <- get(".npplregbw_child_nomad_call_args", ns, inherits = FALSE)
  child.search.text <- paste(deparse(body(child.search)), collapse = "\n")
  child.args.text <- paste(deparse(body(child.args)), collapse = "\n")

  expect_match(child.search.text, ".npplregbw_child_nomad_call_args",
               fixed = TRUE)
  expect_match(child.search.text, "do.call(npregbw, child.args)", fixed = TRUE)
  expect_match(child.args.text, "nomad = TRUE", fixed = TRUE)
  expect_match(child.args.text, "search.engine = degree.search$engine",
               fixed = TRUE)

  dead.symbols <- c(
    ".npplregbw_nomad_prepare_state",
    ".npplregbw_nomad_point_to_matrix",
    ".npplregbw_eval_child_payload_subset",
    ".npplregbw_eval_child_payload_serial",
    ".npplregbw_eval_child_payload_collective",
    "npRmpiNomadShadowSearchPlreg",
    "npRmpiNomadSessionServicePlreg"
  )
  expect_false(any(vapply(dead.symbols, exists, logical(1L), envir = ns,
                          inherits = FALSE)))
})
