test_that("master-local fanout captures chunk errors before draining replies", {
  ns <- asNamespace("npRmpi")
  code <- paste(deparse(body(get(".npRmpi_fanout_bootstrap", ns))),
                collapse = "\n")
  expect_match(code, ".npRmpi_fanout_worker_call(worker", fixed = TRUE)
  condition <- simpleError("required chunk unavailable")
  part <- get(".npRmpi_fanout_worker_call", ns)(function() stop(condition), list())
  expect_s3_class(part, "try-error")
  expect_identical(attr(part, "condition"), condition)
  collect <- get(".npRmpi_bootstrap_collect_chunks", ns)
  expect_error(collect(list(part, matrix(2)),
                       list(list(bsz = 1L), list(bsz = 1L)), 1L),
               "fan-out worker error detected \\(required chunk unavailable\\)")
})
