test_that("master-local fanout captures chunk errors before draining replies", {
  ns <- asNamespace("npRmpi")
  code <- paste(deparse(body(get(".npRmpi_bootstrap_run_fanout", ns))),
                collapse = "\n")
  expect_match(code, "parts.out\\[\\[task.local.idx\\]\\] <- tryCatch\\(")
  expect_match(code, "class = \"try-error\", condition = e")
  condition <- simpleError("required chunk unavailable")
  part <- structure(conditionMessage(condition), class = "try-error",
                    condition = condition)
  expect_identical(attr(part, "condition"), condition)
  collect <- get(".npRmpi_bootstrap_collect_chunks", ns)
  expect_error(collect(list(part, matrix(2)),
                       list(list(bsz = 1L), list(bsz = 1L)), 1L),
               "fan-out worker error detected \\(required chunk unavailable\\)")
})
