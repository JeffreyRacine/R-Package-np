library(npRmpi)

test_that("conditional plot engines use the local conditional helper", {
  con.body <- paste(deparse(body(getFromNamespace(".np_plot_conbandwidth_engine", "npRmpi"))), collapse = "\n")
  cond.body <- paste(deparse(body(getFromNamespace(".np_plot_condbandwidth_engine", "npRmpi"))), collapse = "\n")

  expect_match(con.body, "\\.np_plot_conditional_eval\\(", perl = TRUE)
  expect_match(cond.body, "\\.np_plot_conditional_eval\\(", perl = TRUE)

  expect_no_match(con.body, "do\\.call\\(method\\.fun", perl = TRUE)
  expect_no_match(cond.body, "do\\.call\\(method\\.fun", perl = TRUE)
})
