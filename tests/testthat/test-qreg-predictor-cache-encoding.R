test_that("batched quantile cache keys preserve scalar predictor and ordinate encodings", {
  d <- data.frame(sep = c(0, -0, Inf, -Inf, NA, NaN),
                  g = ordered(c("a", "b", NA, "b", "a", "a")),
                  t = as.Date("2020-01-01") + 0:5,
                  collapse = c("a", "a\rb", "", NA, "x", "x"))
  for (frame in list(d, d[FALSE, ], d[, FALSE, drop = FALSE],
                     transform(d, listcol = I(as.list(seq_len(nrow(d))))))) {
    baseline <- vapply(seq_len(nrow(frame)), function(i)
      .npqreg_selected_cdf_cache_row_key(frame, i), character(1L))
    expect_identical(.npqreg_selected_cdf_cache_row_keys(frame), baseline)
  }
  xy <- .npqreg_selected_cdf_cache_row_keys(d)
  old <- vapply(seq_along(xy), function(i)
    paste(c(xy[i], paste0("y:", sprintf("%a", d$sep[i]))), collapse = "\r"), character(1L))
  expect_identical(.npqreg_selected_cdf_cache_key(xy, d$sep), old)
})
