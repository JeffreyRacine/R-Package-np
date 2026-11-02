test_that("conditional plot schema follows the existing categorical policy", {
  numeric <- data.frame(v = c(1, 2))
  unordered <- data.frame(v = factor(c("a", "b")))
  ordered <- data.frame(v = ordered(c("a", "b")))
  vals <- list(con = numeric, uno = unordered, ord = ordered)
  for (xt in names(vals)) for (yt in names(vals)) {
    b <- list(xbw = .2, ybw = .2)
    for (axis in c("x", "y")) for (kind in names(vals))
      b[[paste0("i", axis, kind)]] <- identical(kind, if (axis == "x") xt else yt)
    expect_null(.np_conditional_training_types(b, vals[[xt]], vals[[yt]]))
    if (xt != "con") expect_error(.np_conditional_training_types(b, numeric, vals[[yt]]), "do not match 'xdat' in type", fixed = TRUE)
    if (yt != "con") expect_error(.np_conditional_training_types(b, vals[[xt]], numeric), "do not match 'ydat' in type", fixed = TRUE)
    expect_error(.np_conditional_training_types(b, cbind(vals[[xt]], vals[[xt]]), vals[[yt]]), "number of columns", fixed = TRUE)
  }
})
