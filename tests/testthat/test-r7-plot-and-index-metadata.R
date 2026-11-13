test_that("coefficient plot controls retain their value validation", {
  for (bad in list(0, -1, "a", c(1, 2), NA_real_, Inf, NULL)) {
    expect_error(.np_plot_normalize_public_dots(list(coef.index = bad)),
                 "coef.index", fixed = TRUE)
    expect_error(.np_plot_scbandwidth_engine(coef.index = bad),
                 "coef.index", fixed = TRUE)
  }
  for (bad in list(NA, "a", NULL)) {
    expect_error(.np_plot_normalize_public_dots(list(common.scale = bad)),
                 "common.scale", fixed = TRUE)
    expect_error(.np_plot_scbandwidth_engine(common.scale = bad),
                 "common.scale", fixed = TRUE)
  }
  expect_identical(.np_plot_normalize_public_dots(
    list(coef.index = 2L, common.scale = FALSE)),
    list(coef.index = 2L, common.scale = FALSE))
  # Preserve the previous numeric-index coercion; this is not an API redesign.
  expect_identical(.np_plot_match_coef_index(1.5), 1L)
  expect_true(.np_plot_normalize_public_dots(
    list(common.scale = c(TRUE, FALSE)))$common.scale)
})
