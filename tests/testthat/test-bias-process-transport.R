test_that("single-chunk corrected process preserves factor payload and RNG", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  F <- matrix(seq_len(27)/27,3L,9L)
  for(law in c("gaussian","mammen","rademacher","inid","fixed","geom")) {
    block <- if(law %in% c("fixed","geom"))2L else NULL
    set.seed(61770)
    expected <- F %*% npRmpi:::.np_bias_process_innovations(9L,1L,31L,31L,law,block)
    rng <- .Random.seed
    set.seed(61770)
    actual <- npRmpi:::.np_bias_process_draws(F,31L,law,block)
    expect_identical(actual,expected)
    expect_identical(.Random.seed,rng)
  }
})
