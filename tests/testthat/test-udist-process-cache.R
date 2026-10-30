test_that("CDF block innovation cache admits only bounded reusable work", {
  budget <- list(bytes = 0, processes = 1L)
  plan <- function(B, m = 26L, L = 2L, bytes = 1e10)
    .np_udist_cache_plan(256L, B, m, 25L, L, budget, bytes, 626L)
  expect_false(is.null(plan(8192L)))
  expect_null(plan(8193L))
  expect_null(plan(8192L, m = 25L))
  expect_null(plan(8192L, L = 1L))
  needed <- plan(33L)$bytes
  expect_null(plan(33L, bytes = needed-1))
  expect_false(is.null(plan(33L, bytes = needed)))
})

test_that("CDF cache preserves exact innovations and chunk RNG states", {
  for (law in c("fixed", "geom")) for (B in c(1L, 2L, 33L)) {
    set.seed(99102906)
    tasks <- .np_udist_cache_prepare(31L, B, law, 2L)
    final <- .Random.seed
    for (task in tasks) {
      assign(".Random.seed", task$rng_state, .GlobalEnv)
      weights <- .np_bias_process_innovations(31L, task$start,
        task$start+task$bsz-1L, B, law, 2L)
      expect_identical(weights, task$weights)
      expect_identical(.Random.seed, task$rng_end)
      F <- matrix(seq_len(93L)/93, 3L, 31L)
      assign(".Random.seed", task$rng_state, .GlobalEnv)
      expect_identical(.np_udist_cache_chunk(F, task), F %*% weights)
      expect_identical(.Random.seed, task$rng_end)
    }
    expect_identical(.Random.seed, final)
    expect_identical(attr(tasks, "rng_final_state"), final)
  }
})
