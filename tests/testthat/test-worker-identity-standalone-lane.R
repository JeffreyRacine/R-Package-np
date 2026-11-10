test_that("worker identity quarantine belongs to the existing standalone lane", {
  path <- test_path("..", "validation", "run_full_mpi_test_shard.R")
  expr <- parse(path)
  # Read only definitions: do not execute the shard runner's command-line entry.
  env <- new.env(parent=baseenv())
  for (node in expr) if(is.call(node) && identical(node[[1L]],quote(`<-`)) &&
    identical(node[[2L]],quote(npRmpi_shard_local_only_files))) eval(node,env)
  expect_true("test-worker-identity-public-admission.R" %in% env$npRmpi_shard_local_only_files())
})
