library(npRmpi)

gen_bw_sel_str <- getFromNamespace("genBwSelStr", "npRmpi")

test_that("bandwidth summaries map base-run multistart labels to one", {
  expect_match(
    gen_bw_sel_str(list(fval = 1.23, ifval = 0, num.feval = NA, num.feval.fast = NA)),
    "achieved on multistart 1",
    fixed = TRUE
  )

  expect_match(
    gen_bw_sel_str(list(fval = 1.23, ifval = 2, num.feval = NA, num.feval.fast = NA)),
    "achieved on multistart 2",
    fixed = TRUE
  )
})
