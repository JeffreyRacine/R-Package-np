test_that("regression percentile inclusive scores match exhaustive tiny families", {
  for (B in c(1L, 2L, 3L, 5L, 8L)) {
    values <- matrix(rep(c(-1, 0, 0, 2, 1), length.out=2L*B), 2L)
    for (truth in list(c(-2, 3), c(0, 0), c(1, -1), c(2, 2))) {
      info <- .np_reg_percentile_rank_info(values)
      expect_identical(.np_reg_percentile_least_content_fast(values, truth),
        as.double(.np_reg_percentile_least_content(info, truth)))
      z <- .np_reg_percentile_e6_scores(values, truth)
      depths <- vapply(seq_len(2L), function(i)
        min(sum(values[i, ] <= truth[i]), sum(values[i, ] >= truth[i])), 0)
      for (i in seq_len(2L)) {
        hits <- vapply(0:B, function(j) {
          k <- .np_reg_percentile_e6_index(B, j)
          info$sorted[i, k] <= truth[i] && truth[i] <= info$sorted[i, B+1L-k]
        }, TRUE)
        expect_identical(z$marginal[i], if (any(hits)) as.double(which(hits)[1L]-1L) else Inf)
      }
      hits <- vapply(0:B, function(j) {
        k <- .np_reg_percentile_e6_family_index(info, j)
        all(info$sorted[, k] <= truth & truth <= info$sorted[, B+1L-k])
      }, TRUE)
      expect_identical(z$family, if(any(hits))as.double(which(hits)[1L]-1L)else Inf)
    }
  }
})

test_that("regression percentile resampling centers and normalizes empirical shape", {
  p <- list(residual=c(-2,-1,1,4),row.variance=c(.5,.6,.7,.8),sigma=c(1,2,3,2))
  s <- .np_reg_percentile_shape(p)
  expect_equal(mean(s$u), 0, tolerance=1e-15)
  expect_equal(mean(s$u^2), 1, tolerance=1e-15)
  set.seed(17); before <- .Random.seed
  .np_reg_percentile_with_seed(91L, .np_reg_percentile_draw_innovations(p, 7L))
  expect_identical(.Random.seed, before)
})
