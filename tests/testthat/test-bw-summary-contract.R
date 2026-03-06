test_that("genBwSelStr suppresses unknown fast counts but prints known ones", {
  gen_bw_sel <- getFromNamespace("genBwSelStr", "npRmpi")

  s.na <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = NA_integer_))
  s.zero <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 0L))
  s.pos <- gen_bw_sel(list(num.feval = 26L, num.feval.fast = 3L))

  expect_match(s.na, "Number of Function Evaluations: 26")
  expect_false(grepl("fast =", s.na, fixed = TRUE))
  expect_true(grepl("Number of Function Evaluations: 26 (fast = 0)", s.zero, fixed = TRUE))
  expect_true(grepl("Number of Function Evaluations: 26 (fast = 3)", s.pos, fixed = TRUE))
})
