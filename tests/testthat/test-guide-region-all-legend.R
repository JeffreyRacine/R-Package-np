test_that("all-band guide labels identify the actual envelope construction", {
  labels <- getFromNamespace(".np_plot_all_band_legend_labels", "npRmpi")
  for (method in c("bootstrap", "asymptotic")) {
    z <- labels(method, "guide-region")
    expect_named(z, c("pointwise", "simultaneous", "bonferroni"))
    expect_true(all(grepl("guide + residual envelope", z, fixed=TRUE)))
    expect_match(z[["simultaneous"]], "finite family", fixed=TRUE)
    expect_false(any(grepl("SE|Gaussian|rank-based", z)))
  }
  expect_match(labels("asymptotic")[["pointwise"]], "SE", fixed=TRUE)
  expect_match(labels("bootstrap")[["simultaneous"]], "rank-based", fixed=TRUE)
})
