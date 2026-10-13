test_that("conmode interval meshes retain explicit grid identity and missing endpoints", {
  helper <- getFromNamespace(".np_plot_conmode_surface_intervals", "np")
  ij <- expand.grid(i = 1:3, j = 1:2)
  point <- as.double(10 * ij$i + ij$j) / 40
  dat <- data.frame(grid_index1 = ij$i, grid_index2 = ij$j,
                    probability = point, lower = point - .2, upper = point + .2)
  for (band in c("pointwise", "simultaneous", "bonferroni")) {
    dat[[paste0(band, ".lower")]] <- point - .2
    dat[[paste0(band, ".upper")]] <- point + .2
  }
  dat$simultaneous.lower[c(2L, 5L)] <- NA_real_
  dat$simultaneous.upper[c(2L, 5L)] <- NA_real_
  p <- list(surface = dat, grid = list(dims = c(3L, 2L)))
  actual <- helper(p, "all")
  expected <- matrix(point, nrow = 3L, ncol = 2L)
  expect_identical(actual$z, expected)
  expect_identical(which(is.na(actual$lower.all$simultaneous)), c(2L, 5L))
  permuted <- p
  permuted$surface <- dat[c(6L, 1L, 4L, 3L, 5L, 2L), ]
  expect_identical(helper(permuted, "all"), actual)
  expect_identical(p$surface, dat)
  bad <- p
  bad$surface$grid_index1[2L] <- bad$surface$grid_index1[1L]
  expect_error(helper(bad, "all"), "duplicated cells", fixed = TRUE)
  bad <- p
  bad$surface$grid_index1[2L] <- 1.5
  expect_error(helper(bad, "all"), "indices are invalid", fixed = TRUE)
  bad <- p
  bad$surface$simultaneous.lower <- NULL
  expect_error(helper(bad, "all"), "column is missing", fixed = TRUE)
})
