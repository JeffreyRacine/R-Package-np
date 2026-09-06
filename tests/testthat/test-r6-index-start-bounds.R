test_that("automatic index starts obey bounds without hiding invalid arithmetic", {
  for (h in c(-1, 0, .1, .1000001, 1, 2, 3)) {
    expect_identical(.npindex_bound_automatic_start(h, .1, 2),
                     min(2, max(.1, h)))
  }
  for (h in c(NA_real_, NaN, Inf, -Inf))
    expect_identical(.npindex_bound_automatic_start(h, .1, 2), h)
  x <- cbind(x = seq_len(24), z = -2 * seq_len(24) + sin(seq_len(24)))
  coord <- .npindex_beta_coordinate_setup(x)
  starts <- cbind(beta = coord$to_search(c(0, .5, 20)), h = c(.5, .1, 2),
                  degree = c(0, 1, 2))
  scale <- .npindex_start_bandwidth_scale(x[, 1], nrow(x))
  set.seed(617)
  before <- .Random.seed
  got <- .npindexbw_prepare_fixed_starts(starts, x, coord, scale, .1, 2)
  expect_identical(.Random.seed, before)
  expect_identical(got[1, ], starts[1, ])
  expect_identical(got[, -2], starts[, -2])
  expect_identical(as.double(got[-1, 2]), c(.1, 2))
})
