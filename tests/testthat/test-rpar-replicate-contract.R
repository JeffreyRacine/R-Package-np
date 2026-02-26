test_that(".mpi_make_replicate_fun captures caller environment semantics", {
  env <- new.env(parent = baseenv())
  env$a <- 5L
  helper <- getFromNamespace(".mpi_make_replicate_fun", "npRmpi")
  f <- helper(quote(a), env = env)

  expect_identical(f(), 5L)
  env$a <- 9L
  expect_identical(f(), 9L)
})

test_that("replicate helpers no longer use eval.parent", {
  par_body <- paste(deparse(body(mpi.parReplicate), width.cutoff = 500L), collapse = " ")
  ipar_body <- paste(deparse(body(mpi.iparReplicate), width.cutoff = 500L), collapse = " ")

  expect_false(grepl("eval.parent", par_body, fixed = TRUE))
  expect_false(grepl("eval.parent", ipar_body, fixed = TRUE))
  expect_match(par_body, "\\.mpi_make_replicate_fun")
  expect_match(ipar_body, "\\.mpi_make_replicate_fun")
})
