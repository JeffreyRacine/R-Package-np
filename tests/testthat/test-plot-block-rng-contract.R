test_that("block tasks replay the master's chunked stream in any execution order", {
  for (sim in c("fixed", "geom")) {
    for (blocklen in c(1L, 3L)) {
      drawer <- .np_block_counts_drawer(n = 12L, B = 7L,
                                       blocklen = blocklen, sim = sim)
      set.seed(7613)
      tasks <- .npRmpi_bootstrap_counts_drawer_tasks(7L, 3L, drawer)
      final <- .Random.seed
      set.seed(7613)
      direct <- lapply(tasks, function(task) {
        drawer(task$start, task$start + task$bsz - 1L)
      })
      expect_identical(.Random.seed, final)
      expect_identical(attr(tasks, "rng_final_state", exact = TRUE), final)
      for (i in rev(seq_along(tasks))) {
        set.seed(100L + i)
        got <- .npRmpi_bootstrap_task_counts_drawer(tasks[[i]], drawer)
        expect_identical(got, direct[[i]])
        expect_equal(colSums(got), rep(12, tasks[[i]]$bsz))
        expect_true(all(is.finite(got) & got >= 0 & got == floor(got)))
      }
      expect_identical(vapply(tasks, `[[`, integer(1L), "bsz"), c(3L, 3L, 1L))
      expect_true(all(vapply(tasks, function(task) is.null(task$counts), logical(1L))))
    }
  }
})

test_that("block length one retains the direct multinomial law", {
  drawer <- .np_block_counts_drawer(12L, 7L, blocklen = 1L, sim = "geom")
  set.seed(7614)
  tasks <- .npRmpi_bootstrap_counts_drawer_tasks(7L, 3L, drawer)
  final <- .Random.seed
  got <- do.call(cbind, lapply(tasks, .npRmpi_bootstrap_task_counts_drawer,
                              counts.drawer = drawer))
  set.seed(7614)
  expected <- stats::rmultinom(7L, 12L, rep(1 / 12, 12L))
  expect_identical(got, expected)
  expect_identical(.Random.seed, final)
})
