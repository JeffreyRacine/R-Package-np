test_that("guide-region annotations describe their operator and common family", {
  descriptor <- getFromNamespace(".np_plot_variability_descriptor", "npRmpi")
  for (method in c("asymptotic", "bootstrap"))
    for (band in c("pointwise", "simultaneous"))
      for (shared in c(FALSE, TRUE)) {
        z <- descriptor(method, band, .05, "bias-corrected", shared=shared,
                        construction="guide-region")$candidates
        expect_true(all(grepl("guide", z, fixed=TRUE)))
        expect_true(all(grepl("residual envelope", z, fixed=TRUE)))
        expect_false(any(grepl("SE|within panel|within each panel|normal quantile|rank-based", z)))
        expect_identical(all(grepl("supplied finite family", z, fixed=TRUE)),
                         identical(band, "simultaneous"))
      }
})
